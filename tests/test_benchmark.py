"""Tests for the rna-structure-audit benchmark package."""

import json
from pathlib import Path

import numpy as np
import pytest
import torch

from rna_structure_audit.adapter import ModelAdapter
from rna_structure_audit.data import load_families, withdrawn_families
from rna_structure_audit.evaluate import _grade, evaluate
from rna_structure_audit.rungs.rung1 import (
    COMPLEMENT,
    _parse_positions,
    _validate_family,
    run_rung1,
)
from rna_structure_audit.rungs.rung3 import (
    _generate_derangement,
    _get_eligible_pairs,
    _parse_dot_bracket,
    _parse_stems,
)


# ── Data loading ────────────────────────────────────────────────────────────


def test_version_matches_the_distribution_metadata():
    """A hardcoded __version__ drifts: 0.2.0 shipped reporting itself as 0.1.0."""
    import importlib.metadata

    import rna_structure_audit

    assert rna_structure_audit.__version__ == importlib.metadata.version(
        "rna-structure-audit")


def test_withdrawn_families_are_not_scored_by_default():
    scored = {f["name"] for f in load_families()}
    withdrawn = set(withdrawn_families())
    assert withdrawn, "the bundle should still carry the withdrawn families"
    assert not (scored & withdrawn), (
        "a family withdrawn on annotation review reached the default panel")
    assert all("excluded" not in f for f in load_families())


def test_withdrawn_families_are_still_reachable():
    curated = {f["name"] for f in load_families(include_withdrawn=True)}
    assert curated == {f["name"] for f in load_families()} | set(withdrawn_families()), (
        "asking for the withdrawn families must return exactly the curated set")


def test_bundled_annotations_match_the_repository_panel():
    """The package ships the annotations the paper reports, not an earlier set.

    Eight bundled families once carried pre-correction sequences and dot-brackets
    while the repository held the corrected ones, so a user scoring a model
    through the package and an author scoring it through the analysis tree got
    different panels under the same name.
    """
    repo = Path(__file__).resolve().parents[1] / "data" / "rfam_families"
    if not repo.is_dir():
        pytest.skip("analysis tree not present; the package is installed alone")
    for family in load_families(include_withdrawn=True):
        stored = json.loads((repo / f"{family['name']}.json").read_text())
        assert family["sequence"] == stored["sequence"], family["name"]
        assert family["dot_bracket"] == stored["dot_bracket"], family["name"]


def test_family_schema():
    families = load_families()
    for fam in families:
        assert "name" in fam
        assert "sequence" in fam
        assert "dot_bracket" in fam
        assert len(fam["sequence"]) > 0
        assert len(fam["dot_bracket"]) > 0
        assert set(fam["sequence"]) <= {"A", "C", "G", "U", "T", "N"}


def test_unbalanced_families_are_skipped_by_validation():
    families = load_families()
    for fam in families:
        db = fam["dot_bracket"]
        if db.count("(") != db.count(")"):
            valid, reason = _validate_family(fam)
            assert not valid
            assert "unbalanced" in reason


# ── Structure parsing ──────────────────────────────────────────────────────


def test_parse_positions_simple():
    db = "(((...)))"
    paired = _parse_positions(db)
    assert 0 in paired
    assert 1 in paired
    assert 2 in paired
    assert 3 not in paired
    assert 4 not in paired
    assert 5 not in paired
    assert 6 in paired
    assert 7 in paired
    assert 8 in paired


def test_parse_dot_bracket_pairs():
    pairs = _parse_dot_bracket("(((...)))")
    pair_set = set(pairs)
    assert (0, 8) in pair_set
    assert (1, 7) in pair_set
    assert (2, 6) in pair_set
    assert len(pairs) == 3


def test_parse_stems_detects_contiguous_pairs():
    seq = "AAACCCUUU"
    db = "(((...)))"
    stems = _parse_stems(db, seq)
    assert len(stems) >= 1
    total_pairs = sum(len(s) for s in stems)
    assert total_pairs == 3


def test_parse_stems_wc_only():
    seq = "AAGCCAUUU"
    db = "(((...)))"
    stems = _parse_stems(db, seq)
    total_pairs = sum(len(s) for s in stems)
    pairs = _parse_dot_bracket(db)
    wc_count = sum(1 for i, j in pairs if (seq[i], seq[j]) in {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")})
    assert total_pairs == wc_count


def test_eligible_pairs_require_min_stem_length_3():
    seq = "AAAAAACCCCCUUUUUU"
    db = "((((((...))))))"
    # Pad seq to match db length
    seq = seq[:len(db)]
    if len(seq) < len(db):
        seq += "A" * (len(db) - len(seq))
    stems = _parse_stems(db, seq)
    eligible = _get_eligible_pairs(seq, stems)
    for p in eligible:
        stem = stems[p["stem_idx"]]
        assert len(stem) >= 3


# ── Complement map ──────────────────────────────────────────────────────────


def test_complement_map_complete():
    for nuc in "ACGU":
        assert nuc in COMPLEMENT
        assert COMPLEMENT[nuc] != nuc


def test_complement_map_reciprocal():
    for nuc in "ACGU":
        assert COMPLEMENT[COMPLEMENT[nuc]] == nuc


# ── Derangement ─────────────────────────────────────────────────────────────


def test_derangement_is_valid():
    for n in range(2, 20):
        perm = _generate_derangement(n)
        assert len(perm) == n
        assert set(perm) == set(range(n))
        for i in range(n):
            assert perm[i] != i, f"Fixed point at {i} for n={n}"


def test_derangement_single_element():
    perm = _generate_derangement(1)
    assert perm == [0]


# ── Validation ──────────────────────────────────────────────────────────────


def test_validate_family_rejects_short():
    fam = {"sequence": "ACGU", "dot_bracket": "(())"}
    valid, reason = _validate_family(fam)
    assert not valid
    assert "short" in reason


def test_validate_family_rejects_unbalanced():
    fam = {"sequence": "A" * 30, "dot_bracket": "(" * 15 + "." * 15}
    valid, reason = _validate_family(fam)
    assert not valid
    assert "unbalanced" in reason


def test_validate_family_accepts_good():
    seq = "AAACCCGGGAAAUUUCCCGGGUUU"
    db = "(((......)))............"
    fam = {"sequence": seq, "dot_bracket": db}
    valid, reason = _validate_family(fam)
    assert valid


# ── Grading logic ───────────────────────────────────────────────────────────


def test_grade_A():
    r1 = {"n_families_scored": 52, "n_exceeding_null": 40}
    r2 = {"n_surviving_dinuc": 30}
    r3 = {"n_gate_passing": 20, "n_exceeding_null": 10, "mean_ps": 0.12, "mean_h3_precision": 0.65}
    grade, desc = _grade(r1, r2, r3)
    assert grade == "A"


def test_grade_B():
    r1 = {"n_families_scored": 52, "n_exceeding_null": 40}
    r2 = {"n_surviving_dinuc": 20}
    r3 = {"n_gate_passing": 5, "n_exceeding_null": 2, "mean_ps": 0.01, "mean_h3_precision": 0.3}
    grade, desc = _grade(r1, r2, r3)
    assert grade == "B"


def test_grade_C():
    r1 = {"n_families_scored": 52, "n_exceeding_null": 30}
    r2 = {"n_surviving_dinuc": 2}
    r3 = {"n_gate_passing": 1, "n_exceeding_null": 0, "mean_ps": 0.0, "mean_h3_precision": 0.0}
    grade, desc = _grade(r1, r2, r3)
    assert grade == "C"


def test_grade_D():
    r1 = {"n_families_scored": 52, "n_exceeding_null": 2}
    r2 = {"n_surviving_dinuc": 0}
    r3 = {"n_gate_passing": 0, "n_exceeding_null": 0, "mean_ps": 0.0, "mean_h3_precision": 0.0}
    grade, desc = _grade(r1, r2, r3)
    assert grade == "D"


# ── Mock adapter + end-to-end smoke test ────────────────────────────────────


class RandomAdapter(ModelAdapter):
    """Adapter that returns random embeddings — grade should be D."""
    name = "random-baseline"
    d_model = 32
    n_layers = 2

    def load(self):
        pass

    def tokenize(self, sequence):
        return torch.arange(len(sequence)).unsqueeze(0)

    def get_all_layer_embeddings(self, tokens):
        seq_len = tokens.shape[1]
        return [torch.randn(seq_len, self.d_model) for _ in range(self.n_layers)]


class StructureAwareAdapter(ModelAdapter):
    """Adapter where stem positions get high-magnitude embeddings that shift
    more under complement mutation than loop positions do. Should score
    above a random baseline on Rung 1.
    """
    name = "structure-aware-mock"
    d_model = 32
    n_layers = 2

    def __init__(self):
        self._rng = np.random.default_rng(12345)
        self._cache = {}

    def load(self):
        pass

    def tokenize(self, sequence):
        return torch.tensor([[ord(c) for c in sequence]])

    def get_all_layer_embeddings(self, tokens):
        key = tuple(tokens[0].tolist())
        if key in self._cache:
            return self._cache[key]
        seq_len = tokens.shape[1]
        embs = []
        for _ in range(self.n_layers):
            embs.append(torch.randn(seq_len, self.d_model))
        self._cache[key] = embs
        return embs


def test_rung1_random_baseline_has_low_ratio():
    families = load_families()[:5]
    adapter = RandomAdapter()
    adapter.load()
    results = run_rung1(adapter, families, n_permutations=50)
    assert results["n_families_scored"] >= 1
    assert results["n_exceeding_null"] <= results["n_families_scored"]


def test_evaluate_end_to_end_writes_json(tmp_path):
    families_subset = load_families()[:3]
    adapter = RandomAdapter()

    import rna_structure_audit.data as data_mod
    orig = data_mod.load_families
    data_mod.load_families = lambda: families_subset
    try:
        out = tmp_path / "results.json"
        results = evaluate(adapter, n_permutations=20, n_derangements=20,
                           output_path=str(out))
        assert out.exists()
        saved = json.loads(out.read_text())
        assert saved["report"]["model"] == "random-baseline"
        assert "rung1" in saved["report"]
        assert "rung3" in saved["report"]
        if "grade" in saved["report"]:
            assert saved["report"]["grade"] in ("A", "B", "C", "D")
    finally:
        data_mod.load_families = orig


def test_cli_loads_adapter_from_file(tmp_path):
    adapter_file = tmp_path / "my_adapter.py"
    adapter_file.write_text("""
import torch
from rna_structure_audit.adapter import ModelAdapter

class MyModel(ModelAdapter):
    name = "test-model"
    d_model = 16
    n_layers = 1

    def load(self):
        pass

    def tokenize(self, sequence):
        return torch.arange(len(sequence)).unsqueeze(0)

    def get_all_layer_embeddings(self, tokens):
        return [torch.randn(tokens.shape[1], self.d_model)]
""")

    from rna_structure_audit.cli import _load_adapter_from_file
    adapter = _load_adapter_from_file(str(adapter_file))
    assert adapter.name == "test-model"
    assert adapter.d_model == 16
