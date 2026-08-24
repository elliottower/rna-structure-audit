"""A container killed mid-run must resume to the same numbers it would have reached.

Checkpointing is only safe if resuming is invisible in the output. These tests
run a stage straight through, run it again against a checkpoint that already
holds part of the work, and require the two results to be equal -- including the
aggregates, which are the part a naive checkpoint gets wrong, because a list
appended inside the loop holds only the families the second container did.
"""

import json

import pytest
import torch

from family_checkpoint import FamilyCheckpoint
from phase6_compensatory_mutation import load_rfam_families, run_phase6
from phases_1_to_5 import run_mutation_sensitivity

VOCAB = {"A": 1, "C": 2, "G": 3, "U": 4, "T": 4, "N": 5}


class StubAdapter:
    """A deterministic stand-in for a language model, with position coupling.

    Every position's embedding mixes every token through a 1/(1+distance)
    kernel, so mutating one nucleotide moves the embedding everywhere and the
    perturbation profiles are non-degenerate. Nothing here is random: the same
    sequence gives the same embeddings in every process, which is what lets the
    interrupted run and the uninterrupted one be compared at all.
    """

    name = "stub"
    d_model = 16
    n_layers = 4

    def tokenize(self, sequence):
        return torch.tensor([[VOCAB.get(c, 5) for c in sequence]])

    def get_all_layer_embeddings(self, tokens):
        ids = tokens[0].to(torch.float64)
        n = len(ids)
        pos = torch.arange(n, dtype=torch.float64)
        kernel = 1.0 / (1.0 + (pos[:, None] - pos[None, :]).abs())
        k = torch.arange(self.d_model, dtype=torch.float64)
        out = []
        for layer in range(self.n_layers):
            site = torch.sin((ids[:, None] + 1.0) * (k[None, :] + 1.0) * 0.37
                             + (layer + 1) * 0.11)
            out.append(kernel @ site)
        return out


class ExplodingAdapter(StubAdapter):
    """Dies on one family, so earlier families are finished and later ones are not.

    Counting calls instead would make the test depend on how many forward passes
    a family happens to need, and would die inside the first one.
    """

    def __init__(self, explode_on):
        self.explode_on = explode_on

    def tokenize(self, sequence):
        if sequence == self.explode_on:
            raise RuntimeError("container reclaimed")
        return super().tokenize(sequence)


@pytest.fixture(scope="module")
def scored_families():
    families = sorted(load_rfam_families(), key=lambda f: len(f["sequence"]))
    kept = []
    for fam in families:
        result = run_phase6(StubAdapter(), [fam], compute_null=False)
        if not result["per_rna"][fam["name"]].get("skipped"):
            kept.append(fam)
        if len(kept) == 3:
            return kept
    pytest.skip("fewer than three families survive the registered filters")


def test_phase6_resume_matches_uninterrupted_run(scored_families, tmp_path):
    straight = run_phase6(StubAdapter(), scored_families, compute_null=True)

    path = tmp_path / "phase6.json"
    run_phase6(StubAdapter(), scored_families[:1], compute_null=True,
               checkpoint=FamilyCheckpoint(path))
    resumed = run_phase6(StubAdapter(), scored_families, compute_null=True,
                         checkpoint=FamilyCheckpoint(path))

    assert json.loads(path.read_text())["done"].keys() == {
        f["name"] for f in scored_families}
    assert resumed == straight


def test_mutation_sensitivity_resume_matches_uninterrupted_run(scored_families, tmp_path):
    straight = run_mutation_sensitivity(StubAdapter(), "rinalmo", scored_families,
                                        device="cpu", n_permutations=200)

    path = tmp_path / "mutation.json"
    run_mutation_sensitivity(StubAdapter(), "rinalmo", scored_families[:2], device="cpu",
                             n_permutations=200, checkpoint=FamilyCheckpoint(path))
    resumed = run_mutation_sensitivity(StubAdapter(), "rinalmo", scored_families,
                                       device="cpu", n_permutations=200,
                                       checkpoint=FamilyCheckpoint(path))

    assert resumed["n_families_scored"] == straight["n_families_scored"]
    assert resumed["mean_best_ratio"] == pytest.approx(straight["mean_best_ratio"])
    assert resumed["median_best_ratio"] == pytest.approx(straight["median_best_ratio"])
    assert resumed == straight


def test_a_crash_leaves_the_finished_families_recoverable(scored_families, tmp_path):
    path = tmp_path / "phase6.json"
    exploding = ExplodingAdapter(explode_on=scored_families[1]["sequence"])
    with pytest.raises(RuntimeError):
        run_phase6(exploding, scored_families, compute_null=True,
                   checkpoint=FamilyCheckpoint(path))

    salvaged = json.loads(path.read_text())["done"]
    assert set(salvaged) == {scored_families[0]["name"]}, (
        "the crash should leave exactly the family that finished before it")

    resumed = run_phase6(StubAdapter(), scored_families, compute_null=True,
                         checkpoint=FamilyCheckpoint(path))
    straight = run_phase6(StubAdapter(), scored_families, compute_null=True)
    assert resumed == straight


def _poison(value):
    """Every number replaced, so a resumed shard cannot match by coincidence."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return 999.0
    if isinstance(value, dict):
        return {k: _poison(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_poison(v) for v in value]
    return value


def _write_shards(path, families, stamp):
    run_phase6(StubAdapter(), families, compute_null=True,
               checkpoint=FamilyCheckpoint(path, stamp=stamp))
    stored = json.loads(path.read_text())
    stored["done"] = _poison(stored["done"])
    path.write_text(json.dumps(stored))
    return stored["done"]


def test_shards_from_a_different_stamp_do_not_reach_the_result(scored_families, tmp_path):
    path = tmp_path / "phase6.json"
    poisoned = _write_shards(path, scored_families[:1],
                             {"commit": "aaa", "libraries": {"torch": "2.1.2"}})
    assert poisoned, "the fixture must leave a shard for the guard to discard"

    resumed = run_phase6(StubAdapter(), scored_families, compute_null=True,
                         checkpoint=FamilyCheckpoint(
                             path, stamp={"commit": "bbb", "libraries": {"torch": "2.1.2"}}))
    straight = run_phase6(StubAdapter(), scored_families, compute_null=True)
    assert resumed == straight


def test_a_library_version_alone_discards_the_shards(scored_families, tmp_path):
    path = tmp_path / "phase6.json"
    _write_shards(path, scored_families[:1],
                  {"commit": "aaa", "libraries": {"torch": "2.1.2"}})

    resumed = run_phase6(StubAdapter(), scored_families, compute_null=True,
                         checkpoint=FamilyCheckpoint(
                             path, stamp={"commit": "aaa", "libraries": {"torch": "2.4.1"}}))
    straight = run_phase6(StubAdapter(), scored_families, compute_null=True)
    assert resumed == straight


def test_shards_written_before_stamps_existed_are_discarded(scored_families, tmp_path):
    path = tmp_path / "phase6.json"
    run_phase6(StubAdapter(), scored_families[:1], compute_null=True,
               checkpoint=FamilyCheckpoint(path))
    assert "stamp" not in json.loads(path.read_text()) or \
        json.loads(path.read_text())["stamp"] is None

    resumed = run_phase6(StubAdapter(), scored_families, compute_null=True,
                         checkpoint=FamilyCheckpoint(path, stamp={"commit": "aaa"}))
    straight = run_phase6(StubAdapter(), scored_families, compute_null=True)
    assert resumed == straight


def test_a_matching_stamp_still_resumes(scored_families, tmp_path):
    path = tmp_path / "phase6.json"
    stamp = {"commit": "aaa", "libraries": {"torch": "2.1.2"}}
    poisoned = _write_shards(path, scored_families[:1], stamp)

    resumed = run_phase6(StubAdapter(), scored_families, compute_null=True,
                         checkpoint=FamilyCheckpoint(path, stamp=stamp))
    first = scored_families[0]["name"]
    assert resumed["per_rna"][first] == poisoned[first], (
        "an unchanged stamp must reuse the stored shard rather than recompute it")
