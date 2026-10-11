"""The registered decision pipeline, on panels built by hand.

Each family record carries per-layer PS and per-layer decoy rates chosen so the
correct layer, population and decision are known in advance.
"""

import json
import math

import pytest

import decoy_audit
from decoy_audit import PanelMismatch, audit, best_layer, sign_test

N_LAYERS = 3


def cell(rate, pairs=5):
    """One decoy tag: the same rate, pairs and counts at every layer."""
    if isinstance(rate, (int, float)) or rate is None:
        rate = [rate] * N_LAYERS
    if isinstance(pairs, int):
        pairs = [pairs] * N_LAYERS
    return {"rate": list(rate), "pairs": list(pairs),
            "comparisons": [2 * n for n in pairs], "ties": [0] * N_LAYERS}


def family(ps, all_rate=0.9, close_rate=0.9, gate=True, **extra):
    entry = {"per_layer_ps": list(ps), "positive_control": {"pass": gate},
             "decoys": {"family_all": cell(all_rate), "family_close": cell(close_rate),
                        "stem_all": cell(None, 0), "stem_close": cell(None, 0)}}
    entry.update(extra)
    return entry


def write(tmp_path, families, key="rinalmo"):
    folder = tmp_path / key
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{key}_phase6_ps.json").write_text(
        json.dumps({"results": {"per_rna": families}}))
    return tmp_path


def uniform(n, **kwargs):
    return {f"fam{k:02d}": family([0.0, 1.0, 0.5], **kwargs) for k in range(n)}


@pytest.mark.parametrize("above, below, at_half, expected", [
    (19, 7, 0, 971712 / 2 ** 26),
    (13, 4, 0, 3214 / 2 ** 17),
    (8, 1, 0, 10 / 2 ** 9),
    (5, 0, 3, 1 / 2 ** 5),
    (0, 0, 4, 1.0),
])
def test_sign_test_is_the_exact_binomial_tail(above, below, at_half, expected):
    result = sign_test([0.9] * above + [0.1] * below + [0.5] * at_half)
    assert result["p"] == pytest.approx(expected)
    assert (result["above"], result["below"], result["at_half"]) == (above, below, at_half)
    assert result["n_sign"] == above + below


def test_best_layer_takes_the_largest_mean_and_the_earliest_on_a_tie():
    assert best_layer([[0.0, 1.0, 0.2], [0.0, 0.2, 0.9]]) == 1
    assert best_layer([[0.3, 0.3, 0.1]]) == 0


def test_a_family_cannot_choose_its_own_layer(tmp_path):
    # Fifteen families peak at layer 1. One family has a huge PS at layer 2,
    # where its decoy rate is 1; at layer 1, where the others send it, it is 0.
    families = uniform(15)
    families["loud"] = family([0.0, 0.0, 100.0], all_rate=[0.5, 0.0, 1.0])
    result = audit(write(tmp_path, families), "rinalmo")
    assert result["panel_layer"] == 2
    assert result["held_out_layers"]["loud"] == 1
    assert result["family_all"]["family_rates"]["loud"] == 0.0
    # And it drags the other families' layer, which is what leave-one-out means.
    assert result["held_out_layers"]["fam00"] == 2


def test_quarantined_skipped_and_gate_failing_families_are_excluded_with_reasons(tmp_path):
    families = uniform(15)
    families["tRNA_Phe_yeast"] = family([9.0, 0.0, 0.0])
    families["short"] = {"skipped": True, "reason": "< 15 WC pairs (9 found)"}
    families["flat"] = family([9.0, 0.0, 0.0], gate=False)
    result = audit(write(tmp_path, families), "rinalmo")
    assert result["population_families"] == 15
    assert set(result["excluded_families"]) == {"tRNA_Phe_yeast", "short", "flat"}
    assert "positive control" in result["excluded_families"]["flat"]
    # None of the three moved the layer: each would have voted for layer 0.
    assert set(result["held_out_layers"].values()) == {1}


def test_the_pair_minimum_is_read_at_the_selected_layer(tmp_path):
    families = uniform(15)
    # Five pairs everywhere except the selected layer 1, where two contribute.
    thin = family([0.0, 1.0, 0.5])
    thin["decoys"]["family_all"] = cell(0.9, pairs=[5, 2, 5])
    families["thin"] = thin
    result = audit(write(tmp_path, families), "rinalmo")
    assert "thin" not in result["family_all"]["family_rates"]
    assert result["family_all"]["families_below_pair_minimum"] == {"thin": 2}
    assert result["family_all"]["families"] == 15


def test_families_at_one_half_are_counted_and_leave_the_denominator(tmp_path):
    families = uniform(14)
    families.update({f"tie{k}": family([0.0, 1.0, 0.5], all_rate=0.5) for k in range(3)})
    result = audit(write(tmp_path, families), "rinalmo")["family_all"]
    assert (result["families"], result["at_half"], result["n_sign"]) == (17, 3, 14)
    assert result["p"] == pytest.approx(1 / 2 ** 14)


@pytest.mark.parametrize("n, void", [(13, True), (14, False)])
def test_identity_criterion_is_void_below_fourteen_families(tmp_path, n, void):
    result = audit(write(tmp_path, uniform(n)), "rinalmo")
    assert result["identity_criterion"] == {"void": void, "met": not void}


@pytest.mark.parametrize("n_close, void", [(7, True), (8, False)])
def test_distance_guard_is_void_below_eight_families(tmp_path, n_close, void):
    families = uniform(14)
    for k, name in enumerate(sorted(families)):
        if k >= n_close:
            families[name]["decoys"]["family_close"] = cell(None, 0)
    result = audit(write(tmp_path, families), "rinalmo")
    assert result["distance_guard"]["void"] is void
    assert result["distance_guard"]["met"] is (not void)


def test_the_guard_is_confirmatory_only_behind_a_met_identity_criterion(tmp_path):
    # Identity fails (rates below one half); the close set alone would pass.
    families = uniform(14, all_rate=0.1, close_rate=0.9)
    result = audit(write(tmp_path, families), "rinalmo")
    assert result["identity_criterion"]["met"] is False
    assert result["family_close"]["p"] < decoy_audit.ALPHA
    assert result["distance_guard"] == {"confirmatory": False, "void": False, "met": False}
    assert result["decisions"]["H1"]["met"] is False
    assert result["decisions"]["G1"]["met"] is False


def test_decisions_carry_the_registered_names(tmp_path):
    panel = write(tmp_path, uniform(14), "ernierna")
    write(tmp_path, uniform(14), "caduceus")
    assert set(audit(panel, "ernierna")["decisions"]) == {"H2", "G2"}
    assert "decisions" not in audit(panel, "caduceus")


def test_a_missing_decoy_tag_is_refused(tmp_path):
    families = uniform(14)
    del families["fam03"]["decoys"]["family_close"]
    with pytest.raises(PanelMismatch, match="family_close"):
        audit(write(tmp_path, families), "rinalmo")


def test_mismatched_layer_counts_are_refused(tmp_path):
    families = uniform(14)
    families["fam03"]["decoys"]["family_all"]["rate"].append(0.9)
    with pytest.raises(PanelMismatch, match="layers"):
        audit(write(tmp_path, families), "rinalmo")
    families = uniform(14)
    families["fam03"]["per_layer_ps"] = [0.0, 1.0]
    with pytest.raises(PanelMismatch, match="PS layers"):
        audit(write(tmp_path, families), "rinalmo")


def test_a_ps_that_is_not_finite_is_refused(tmp_path):
    families = uniform(14)
    families["fam03"]["per_layer_ps"] = [0.0, math.nan, 0.5]
    with pytest.raises(PanelMismatch, match="not finite"):
        audit(write(tmp_path, families), "rinalmo")


def test_a_panel_without_decoy_rates_is_refused(tmp_path):
    families = uniform(14)
    del families["fam00"]["decoys"]
    with pytest.raises(PanelMismatch, match="v3"):
        audit(write(tmp_path, families), "rinalmo")
