"""The registered test of the matched-decoy comparison.

Reads a panel written by `phase6_compensatory_mutation_v3.py`. For each model:

  - the population is the Phase 6 population: confirmatory families (not
    quarantined, not skipped) that pass the stem-versus-loop positive control;
  - the layer for a family is chosen without that family: the layer maximizing
    mean PS over the model's other population families, whether or not they
    have a matched decoy (leave one family out);
  - the family's decoy rate is read at that layer, and the family enters a test
    when three or more of its pairs contribute there;
  - an exact one-sided sign test asks whether more families sit above 0.5 than
    below it, families at exactly 0.5 counted and dropped from the denominator.

Decisions, under family scope. The identity criterion (decoy set "all") is met
at p < ALPHA. The distance guard (decoy set "close") is read as a confirmation
only when the identity criterion is met, and is then met at p < ALPHA. Either
is void below its MIN_FAMILIES.

    uv run --no-project python scripts/decoy_audit.py --panel results/<panel>
"""

import argparse
import json
import math
import statistics
from pathlib import Path

from decoy_comparison import DECOY_SETS, SCOPES

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "audit" / "decoy_audit.json"

QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
MODELS = ["rinalmo", "ernierna", "caduceus", "rnafm", "utrlm", "splicebert",
          "hyenadna", "evo", "ernierna_untrained", "ernierna_untrained_noattnbias"]
# Hypothesis and guard names of the registration, per confirmatory model.
CONFIRMATORY = {"rinalmo": ("H1", "G1"), "ernierna": ("H2", "G2")}
ALPHA = 0.025
MIN_PAIRS = 3
# Four fifths of the 17 and 9 families the structure makes available.
MIN_FAMILIES = {"all": 14, "close": 8}
TAGS = [f"{scope}_{decoy_set}" for scope in SCOPES for decoy_set in DECOY_SETS]
FIELDS = ("rate", "pairs", "comparisons", "ties")


class PanelMismatch(RuntimeError):
    """The panel does not carry what the registered analysis reads."""


def sign_test(values: list[float]) -> dict:
    """Exact one-sided test that more values exceed 0.5 than fall below it."""
    above = sum(value > 0.5 for value in values)
    below = sum(value < 0.5 for value in values)
    n = above + below
    p = (sum(math.comb(n, k) for k in range(above, n + 1)) / 2 ** n) if n else 1.0
    return {"above": above, "below": below, "at_half": len(values) - n,
            "n_sign": n, "p": p}


def best_layer(ps_rows: list[list[float]]) -> int:
    """The layer with the largest mean PS; ties go to the earliest layer."""
    if not ps_rows:
        raise PanelMismatch("no family left to choose a layer from")
    n_layers = len(ps_rows[0])
    means = [sum(row[layer] for row in ps_rows) / len(ps_rows)
             for layer in range(n_layers)]
    return max(range(n_layers), key=means.__getitem__)


def load(panel: Path, key: str) -> tuple[dict, dict]:
    """(population families, excluded families with the reason for each)."""
    stored = json.loads((panel / key / f"{key}_phase6_ps.json").read_text())
    families, excluded = {}, {}
    for name, entry in sorted(stored["results"]["per_rna"].items()):
        if name in QUARANTINED:
            excluded[name] = "quarantined pilot family"
        elif entry.get("skipped"):
            excluded[name] = f"skipped: {entry.get('reason', 'no reason recorded')}"
        elif not entry.get("positive_control", {}).get("pass", False):
            excluded[name] = "failed the stem-versus-loop positive control"
        else:
            families[name] = entry
    if not families:
        raise PanelMismatch(f"{key}: no family in the population")

    n_layers = len(next(iter(families.values()))["per_layer_ps"])
    for name, entry in families.items():
        ps = entry["per_layer_ps"]
        if len(ps) != n_layers:
            raise PanelMismatch(f"{key}/{name}: {len(ps)} PS layers, expected {n_layers}")
        if any(value is None or not math.isfinite(value) for value in ps):
            raise PanelMismatch(f"{key}/{name}: a per-layer PS is missing or not finite")
        decoys = entry.get("decoys")
        if decoys is None:
            raise PanelMismatch(f"{key}/{name}: no decoy rates; rerun with the v3 script")
        for tag in TAGS:
            if tag not in decoys:
                raise PanelMismatch(f"{key}/{name}: decoy tag {tag!r} is missing")
            for field in FIELDS:
                if len(decoys[tag][field]) != n_layers:
                    raise PanelMismatch(
                        f"{key}/{name}: decoys[{tag!r}][{field!r}] has "
                        f"{len(decoys[tag][field])} layers, expected {n_layers}")
    return families, excluded


def evaluate_tag(families: dict, layers: dict, tag: str) -> dict:
    rates, below_minimum = {}, {}
    pairs = comparisons = ties = 0
    for name, entry in families.items():
        stored = entry["decoys"][tag]
        layer = layers[name]
        if stored["rate"][layer] is None or stored["pairs"][layer] < MIN_PAIRS:
            below_minimum[name] = stored["pairs"][layer]
            continue
        rates[name] = stored["rate"][layer]
        pairs += stored["pairs"][layer]
        comparisons += stored["comparisons"][layer]
        ties += stored["ties"][layer]
    values = list(rates.values())
    return {
        "families": len(values),
        "pairs": pairs,
        "comparisons": comparisons,
        "tie_fraction": ties / comparisons if comparisons else None,
        "mean_rate": statistics.fmean(values) if values else None,
        "median_rate": statistics.median(values) if values else None,
        **sign_test(values),
        "family_rates": rates,
        "families_below_pair_minimum": below_minimum,
    }


def audit(panel: Path, key: str) -> dict:
    families, excluded = load(panel, key)
    names = list(families)
    panel_layer = best_layer([families[name]["per_layer_ps"] for name in names])
    if len(names) > 1:
        held_out = {name: best_layer([families[other]["per_layer_ps"]
                                      for other in names if other != name])
                    for name in names}
    else:
        held_out = {}

    result = {"model": key, "population_families": len(names),
              "excluded_families": excluded,
              "panel_layer": panel_layer, "held_out_layers": held_out}
    for tag in TAGS:
        result[tag] = evaluate_tag(families if held_out else {}, held_out, tag)
        # Description only: the same statistic at the one layer the whole
        # population chooses.
        at_panel_layer = evaluate_tag(families, dict.fromkeys(names, panel_layer), tag)
        result[tag]["panel_layer_mean_rate"] = at_panel_layer["mean_rate"]

    identity, guard = result["family_all"], result["family_close"]
    identity_void = identity["families"] < MIN_FAMILIES["all"]
    identity_met = (not identity_void) and identity["p"] < ALPHA
    guard_void = guard["families"] < MIN_FAMILIES["close"]
    result["identity_criterion"] = {"void": identity_void, "met": identity_met}
    result["distance_guard"] = {
        # Read as a confirmation only behind a met identity criterion.
        "confirmatory": identity_met,
        "void": guard_void,
        "met": identity_met and (not guard_void) and guard["p"] < ALPHA,
    }
    if key in CONFIRMATORY:
        hypothesis, guard_name = CONFIRMATORY[key]
        result["decisions"] = {hypothesis: result["identity_criterion"],
                               guard_name: result["distance_guard"]}
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--models", nargs="+", default=MODELS)
    args = parser.parse_args()

    rows = [audit(args.panel, key) for key in args.models]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "panel": str(args.panel), "alpha": ALPHA, "min_pairs": MIN_PAIRS,
        "min_families": MIN_FAMILIES, "confirmatory_models": CONFIRMATORY,
        "models": rows}, indent=2))
    for row in rows:
        cells = "  ".join(
            f"{s}: {row[f'family_{s}']['above']}/{row[f'family_{s}']['n_sign']} "
            f"p={row[f'family_{s}']['p']:.4f}" for s in DECOY_SETS)
        print(f"  {row['model']:32s} {cells}  identity={row['identity_criterion']}  "
              f"guard={row['distance_guard']}")
    print(f"wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
