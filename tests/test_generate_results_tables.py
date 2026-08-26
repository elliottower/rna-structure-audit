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


def _write_run(directory, key, stamp, mutation=None):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "stamp.json").write_text(json.dumps(stamp))
    (directory / f"{key}_phases_1_to_5.json").write_text(json.dumps(
        {"mutation_trained": {"per_rna": mutation
                              or {"a": mutation_entry(1.2, True, True)}},
         "attention_trained": {"per_rna": {"a": {"best_corr": 0.1, "best_layer": 2}}},
         "probing": {"best_accuracy": 0.61, "best_layer": 4}}))
    (directory / f"{key}_phase6_ps.json").write_text(json.dumps(
        {"results": {"per_rna": {"a": phase6_entry(0.05)}}}))


TRANSVERSION = {"A": "C", "U": "C", "C": "A", "G": "U", "T": "C"}


def _write_transversion(root, key, ratio=1.15, alphabet=None):
    directory = root / f"{key}_transversion"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"{key}_transversion.json").write_text(json.dumps(
        {"experiment": "transversion_control",
         "complement": alphabet or TRANSVERSION,
         "mutation_trained": {"per_rna": {"a": mutation_entry(ratio, True, True)}}}))


def _write_panel(root, mutation=None, transversion_ratio=1.15):
    for key in _all_keys():
        _write_run(root / key, key, _stamp(), mutation=mutation)
    for key, *_ in g.MODELS:
        _write_transversion(root, key, ratio=transversion_ratio)


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
    for key, *_ in g.MODELS:
        _write_transversion(tmp_path, key)
    rung3 = {key: g.rung3_stats(runs[key]) for key in runs}
    written = g.macros(runs, rung3, g.mutation_all(runs), g.load_transversion())
    assert r"\provenanceCommits}{3}" in written


def test_a_run_on_the_unrepaired_panel_aborts_the_build(tmp_path, monkeypatch):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    for key in [k for k, *_ in g.MODELS] + [f"{k}_untrained" for k in g.UNTRAINED_KEYS]:
        _write_run(tmp_path / key, key, _stamp(n_families=52, withdrawn=[]))

    with pytest.raises(ValueError, match="not the repaired panel"):
        g.load_all()


# ---------------------------------------------------------------------------
# One value per quantity
# ---------------------------------------------------------------------------


def test_the_table_and_the_macro_quote_the_same_interval(tmp_path, monkeypatch):
    """v12 drew each interval from wherever the shared generator happened to be.

    The bootstrap needs a stream, and threading one generator through the
    renderers makes an interval depend on which renderer ran first: the table
    prints one draw and the macro quoting that table prints another. The
    families here are spread widely enough that two draws would disagree in the
    second decimal, which is the digit both the table and the macro print.
    """
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    spread = {f"f{i}": mutation_entry(1.0 + 0.2 * i, True, True) for i in range(12)}
    _write_panel(tmp_path, mutation=spread)

    runs = g.load_all()
    mut = g.mutation_all(runs)
    table = g.rung1_table(runs, mut)
    written = g.macros(runs, {key: g.rung3_stats(runs[key]) for key in runs},
                       mut, g.load_transversion())

    for key, _short, _size, _domain, name in g.MODELS:
        interval = g.fmt_ci(mut[key]["ci"])
        assert interval in table
        assert "\\newcommand{\\ci" + name + "}{" + interval + "}" in written


def test_an_intervals_value_does_not_depend_on_what_ran_before_it(tmp_path, monkeypatch):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    spread = {f"f{i}": mutation_entry(1.0 + 0.2 * i, True, True) for i in range(12)}
    _write_panel(tmp_path, mutation=spread)
    runs = g.load_all()

    once = g.mutation_all(runs)
    again = g.mutation_all({key: runs[key] for key in reversed(list(runs))})

    for key in runs:
        assert once[key]["ci"] == again[key]["ci"]
        assert once[key]["retention"] == again[key]["retention"]


# ---------------------------------------------------------------------------
# The transversion control
# ---------------------------------------------------------------------------


def test_an_alphabet_that_leaves_a_nucleotide_alone_aborts_the_build(tmp_path, monkeypatch):
    """A fixed point drops those positions from both means without saying so."""
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    _write_panel(tmp_path)
    _write_transversion(tmp_path, "evo",
                        alphabet={**TRANSVERSION, "G": "G"})

    with pytest.raises(ValueError, match="unchanged"):
        g.load_transversion()


def test_two_alphabets_across_the_panel_abort_the_build(tmp_path, monkeypatch):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    _write_panel(tmp_path)
    _write_transversion(tmp_path, "evo", alphabet={**TRANSVERSION, "G": "C"})

    with pytest.raises(ValueError, match="alphabets"):
        g.load_transversion()


def test_a_missing_control_aborts_rather_than_reporting_the_models_that_ran(
        tmp_path, monkeypatch):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    _write_panel(tmp_path)
    for path in (tmp_path / "hyenadna_transversion").iterdir():
        path.unlink()

    with pytest.raises(FileNotFoundError, match="hyenadna"):
        g.load_transversion()


def test_the_bound_the_prose_quotes_covers_the_widest_model(tmp_path, monkeypatch):
    """v12 printed a model at 7.2% four lines under a sentence claiming 5%.

    The bound is a macro rather than a typed number precisely so that a model
    drifting past it moves the sentence instead of contradicting it.
    """
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    _write_panel(tmp_path, mutation={"a": mutation_entry(1.0, True, True)},
                 transversion_ratio=1.072)
    runs = g.load_all()
    controls = g.load_transversion()
    mut = g.mutation_all(runs)

    written = g.macros(runs, {key: g.rung3_stats(runs[key]) for key in runs},
                       mut, controls)

    assert "{\\transversionBound}{8\\%}" in written
    for key, *_ in g.MODELS:
        deviation = abs(100 * (g.transversion_ratio(controls[key])
                               - mut[key]["mean_ratio"]) / mut[key]["mean_ratio"])
        assert deviation <= 8


# ---------------------------------------------------------------------------
# Claims the numbers can reverse
# ---------------------------------------------------------------------------


def test_a_reversed_ranking_stops_the_build_and_names_the_sentence(tmp_path, monkeypatch):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    _write_panel(tmp_path)
    # Every model identical, so nothing the body calls "highest" is highest.
    runs = g.load_all()
    mut = g.mutation_all(runs)
    rung3 = {key: g.rung3_stats(runs[key]) for key in runs}

    failures = g.claim_failures(runs, rung3, mut, g.load_transversion())

    joined = "\n".join(failures)
    assert "RNA-FM has the highest mean ratio" in joined
    assert "RiNALMo, ERNIE-RNA and Caduceus are the top three" in joined
    with pytest.raises(ValueError, match="contradicts the body"):
        g.check_claims(runs, rung3, mut, g.load_transversion())


def test_every_reversed_claim_is_reported_not_only_the_first(tmp_path, monkeypatch):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    _write_panel(tmp_path)
    runs = g.load_all()

    failures = g.claim_failures(runs, {key: g.rung3_stats(runs[key]) for key in runs},
                                g.mutation_all(runs), g.load_transversion())

    assert len(failures) > 3, "one message per run makes the second re-run pointless"


def test_peak_layer_reports_the_modal_layer_and_how_many_families_reach_it():
    entries = {f"f{i}": {**phase6_entry(0.2), "best_layer": 33 if i else 4,
                         "per_layer_ps": [0.0] * 34}
               for i in range(5)}
    stats = g.rung3_stats(run(phase6=entries))

    assert g.peak_layer(stats) == (33, 4, 5)
    assert stats["n_ps_layers"] == 34


def test_the_transversion_table_leads_with_the_model_that_moved_most(tmp_path, monkeypatch):
    monkeypatch.setattr(g, "RESULTS", tmp_path)
    for key in _all_keys():
        _write_run(tmp_path / key, key, _stamp())
    for index, (key, *_) in enumerate(g.MODELS):
        # The widest sits mid-list, so MODELS order alone cannot pass this.
        _write_transversion(tmp_path, key,
                            ratio=3.0 if index == 6 else 1.0 + 0.02 * index)
    runs = g.load_all()

    table = g.transversion_table(g.mutation_all(runs), g.load_transversion())

    body = table.split(r"\midrule")[1]
    assert all(short in body for _key, short, *_ in g.MODELS), (
        "a model missing from the table is a control silently dropped")
    order = sorted((short for _key, short, *_ in g.MODELS), key=body.index)
    mut, controls = g.mutation_all(runs), g.load_transversion()
    widest = max(g.MODELS, key=lambda m: abs(
        g.transversion_ratio(controls[m[0]]) - mut[m[0]]["mean_ratio"]))
    assert order[0] == widest[1], (
        "the widest deviation belongs at the top, beside the bound the prose quotes")
