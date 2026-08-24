"""The generated tables must reproduce the registered aggregations, not near ones.

Every quantity here has a registered definition that a plausible reimplementation
gets wrong in a specific way: the Rung 3 mean is taken over non-quarantined
families rather than gate-passing ones, retention is denominated in the families
that exceeded the first-order null rather than in every scored family, and the
count criterion scales with the panel. Each test is built so it fails under the
plausible wrong version.
"""

import json
import math

import numpy as np
import pytest

import generate_results_tables as g


def phase6_entry(ps, gate=True, quarantined=False, exceeds=True,
                 partner_max=6, n_pairs=9):
    return {
        "best_ps": ps,
        "positive_control": {"pass": gate},
        "quarantined": quarantined,
        "exceeds_null_primary": exceeds,
        "exceeds_null_conservative": exceeds,
        "h3_precision": {"fraction": partner_max / n_pairs, "n": n_pairs,
                         "partner_max_count": partner_max, "p_value": 0.01},
    }


def mutation_entry(ratio, exceeds_nuc, exceeds_dinuc=None):
    entry = {"best_ratio": ratio, "best_layer": 3, "nuc_null_95th": 1.0,
             "exceeds_nuc_null": exceeds_nuc, "n_stem": 20, "n_loop": 20}
    if exceeds_nuc:
        entry["dinuc_null_95th"] = 1.0
        entry["exceeds_dinuc_null"] = bool(exceeds_dinuc)
    return entry


def run(phase6=None, mutation=None, attention=None, accuracy=0.6):
    return {
        "phases": {
            "mutation_trained": {"per_rna": mutation or {}},
            "attention_trained": ({"per_rna": attention} if attention is not None
                                  else {"skipped": "no attention (SSM architecture)"}),
            "probing": {"best_accuracy": accuracy, "best_layer": 7},
        },
        "phase6": {"per_rna": phase6 or {}},
        "stamp": {},
    }


def test_registered_mean_ps_is_over_non_quarantined_not_gate_passing():
    entries = {
        "a": phase6_entry(0.20, gate=True),
        "b": phase6_entry(0.10, gate=True),
        # A gate-failing family still counts toward the mean, which is the whole
        # difference between the registered aggregation and the stored one.
        "c": phase6_entry(0.00, gate=False),
        "tRNA_Phe_yeast": phase6_entry(0.90, gate=True, quarantined=True),
    }
    stats = g.rung3_stats(run(phase6=entries))

    assert stats["mean_ps"] == pytest.approx(0.10)
    assert stats["mean_ps_gated"] == pytest.approx(0.15)
    assert stats["eligible"] == 3
    assert stats["gate"] == 2


def test_h3_precision_is_the_mean_of_fractions_not_the_pooled_count():
    entries = {
        "a": phase6_entry(0.2, partner_max=1, n_pairs=1),
        "b": phase6_entry(0.2, partner_max=0, n_pairs=99),
    }
    stats = g.rung3_stats(run(phase6=entries))

    assert stats["h3"] == pytest.approx(0.5)
    assert stats["h3_counts"] == (1, 100)


def test_retention_is_denominated_in_families_exceeding_the_first_order_null():
    mutation = {
        "a": mutation_entry(1.5, True, True),
        "b": mutation_entry(1.4, True, False),
        "c": mutation_entry(1.0, False),
        "d": mutation_entry(0.9, False),
    }
    stats = g.mutation_stats(run(mutation=mutation), np.random.default_rng())

    assert stats["n_scored"] == 4
    assert stats["exceeds_nuc"] == 2
    assert stats["survives_dinuc"] == 1
    assert 0.0 <= stats["retention"]["ci_lower"] <= 0.5 <= stats["retention"]["ci_upper"]
    assert stats["retention"]["point_estimate"] == pytest.approx(0.5)


def test_skipped_families_reach_no_aggregate():
    mutation = {"a": mutation_entry(2.0, True, True), "b": {"skipped": "too short"}}
    stats = g.mutation_stats(run(mutation=mutation), np.random.default_rng())

    assert stats["n_scored"] == 1
    assert stats["mean_ratio"] == pytest.approx(2.0)


def test_gate_threshold_scales_with_the_panel():
    assert g.gate_threshold(32) == 7
    assert g.gate_threshold(36) == 8
    assert all(g.gate_threshold(n) == math.ceil(0.2 * n) for n in range(1, 200))


def test_h1_six_needs_both_the_count_and_a_positive_wilcoxon():
    strong = {f"f{i}": phase6_entry(0.2 + 0.001 * i) for i in range(20)}
    at_panel, at_registered, p_value = g.h1_six(g.rung3_stats(run(phase6=strong)))
    assert p_value < 0.0167
    assert at_panel and at_registered

    # Same count of families exceeding the null, but PS centered on zero, so the
    # Wilcoxon condition fails and the decision must fail with it.
    flat = {f"f{i}": phase6_entry(0.001 if i % 2 else -0.001) for i in range(20)}
    at_panel, at_registered, _p = g.h1_six(g.rung3_stats(run(phase6=flat)))
    assert not at_panel and not at_registered


def test_a_model_failing_every_gate_prints_no_null_column():
    entries = {f"f{i}": phase6_entry(0.01, gate=False) for i in range(5)}
    stats = g.rung3_stats(run(phase6=entries))
    rung3 = {key: stats for key in
             [k for k, *_ in g.MODELS] + [f"{k}_untrained" for k in g.UNTRAINED_KEYS]}
    table = g.rung3_table({}, rung3)

    # The gate column still reports 0 of 5, and both null columns go blank
    # rather than reporting an exceedance out of an empty denominator.
    assert table.count("& 0/5 & --- & --- &") == len(g.MODELS) + len(g.UNTRAINED_KEYS)
    assert "/0" not in table


def _write_run(directory, key, stamp):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "stamp.json").write_text(json.dumps(stamp))
    (directory / f"{key}_phases_1_to_5.json").write_text(json.dumps(
        {"mutation_trained": {"per_rna": {"a": mutation_entry(1.2, True, True)}},
         "attention_trained": {"per_rna": {"a": {"best_corr": 0.1, "best_layer": 2}}},
         "probing": {"best_accuracy": 0.61, "best_layer": 4}}))
    (directory / f"{key}_phase6_ps.json").write_text(json.dumps(
        {"results": {"per_rna": {"a": phase6_entry(0.05)}}}))


def _stamp(**overrides):
    base = {"commit": "abc123ff", "panel_sha256": "f" * 64, "n_families": g.PANEL_N,
            "withdrawn": ["a", "b", "c", "d", "e"], "seeding": "family name",
            "device": "cuda",
            "libraries": {"torch": "2.4.1", "transformers": "4.44.2",
                          "multimolecule": None, "flash-attn": None,
                          "mamba-ssm": None}}
    base.update(overrides)
    return base


def _all_keys():
    return [k for k, *_ in g.MODELS] + [f"{k}_untrained" for k in g.UNTRAINED_KEYS]


@pytest.mark.parametrize("field, value", [
    ("panel_sha256", "e" * 64),
    ("withdrawn", ["a", "b", "c", "d", "z"]),
    ("seeding", "a fixed global seed"),
])
def test_one_run_disagreeing_on_the_data_aborts_the_build(tmp_path, monkeypatch,
                                                          field, value):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    for key in _all_keys():
        _write_run(tmp_path / key, key, _stamp())
    assert len(g.load_all()) == len(_all_keys())

    _write_run(tmp_path / "evo", "evo", _stamp(**{field: value}))
    with pytest.raises(ValueError, match=field):
        g.load_all()


def test_a_differing_commit_is_carried_into_the_provenance_table(tmp_path, monkeypatch):
    """The one stamp field allowed to differ, because no stack loads every model.

    Allowing it is only defensible if the difference reaches the reader, so the
    build must not abort and every distinct commit must appear in the table.
    """
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    for key in _all_keys():
        _write_run(tmp_path / key, key, _stamp())
    _write_run(tmp_path / "evo", "evo", _stamp(commit="def4567"))
    _write_run(tmp_path / "nt_untrained", "nt_untrained",
               _stamp(commit="9876543", libraries={"torch": "2.1.2",
                                                   "transformers": "4.49.0",
                                                   "flash-attn": "2.5.8"}))

    runs = g.load_all()
    table = g.provenance_table(runs)

    assert table.count(r"\\") == len(runs) + 1  # one per run, plus the header
    for short in ("abc123f", "def4567", "9876543"):
        assert short in table
    assert "flash-attn 2.5.8" in table
    assert r"\provenanceCommits}{3}" in g.macros(runs, {
        key: g.rung3_stats(runs[key]) for key in runs})


def test_a_run_on_the_unrepaired_panel_aborts_the_build(tmp_path, monkeypatch):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    for key in [k for k, *_ in g.MODELS] + [f"{k}_untrained" for k in g.UNTRAINED_KEYS]:
        _write_run(tmp_path / key, key, _stamp(n_families=52, withdrawn=[]))

    with pytest.raises(ValueError, match="not the repaired panel"):
        g.load_all()
