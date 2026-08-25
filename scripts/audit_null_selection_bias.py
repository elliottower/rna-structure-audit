"""Why untrained controls exceed a 95th-percentile null about half the time.

`best_ps` is the maximum of mean PS over all layers
(`phase6_compensatory_mutation.py:243`). `null_95th_primary` is the 95th
percentile of the derangement null evaluated at that one selected layer
(line 341, 351). `exceeds_null_primary` compares the two.

A maximum over L layers tested against a threshold calibrated for a single
layer fires, under the null, with probability 1 - 0.95^L rather than 0.05. The
`conservative` variant takes the maximum over layers on the null side too
(line 348), which is the matched comparison.

This reports, per model, the predicted false-positive rate from its layer count
against the rate actually observed in the randomly initialized controls, where
the true PS is zero and every exceedance is a false positive.

    uv run --no-project --python 3.12 python scripts/audit_null_selection_bias.py \
        --results results/repaired_panel
"""

import argparse
import json
from pathlib import Path
from statistics import median

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
ALPHA = 0.05


def read(path):
    per_rna = json.loads(path.read_text())["results"]["per_rna"]
    return {name: entry for name, entry in per_rna.items()
            if name not in QUARANTINED and entry.get("best_ps") is not None}


def layer_count(families):
    """Layers the maximum is taken over, from the stored per-layer profile."""
    lengths = [len(e["per_layer_ps"]) for e in families.values()
               if e.get("per_layer_ps")]
    return int(median(lengths)) if lengths else 0


def rate(families, field):
    decided = [e for e in families.values() if e.get(field) is not None]
    if not decided:
        return None, 0, 0
    hits = sum(1 for e in decided if e[field])
    return hits / len(decided), hits, len(decided)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=str(REPO / "results/repaired_panel"))
    args = parser.parse_args()
    root = Path(args.results)

    print("A maximum over L layers, tested against a threshold calibrated for one\n"
          "layer, fires under the null with probability 1 - 0.95^L.\n")
    print(f"  {'model':22s} {'L':>3s} {'predicted':>10s} {'primary':>14s} "
          f"{'conservative':>14s}")
    for directory in sorted(root.iterdir()):
        path = directory / f"{directory.name}_phase6_ps.json"
        if not path.exists():
            continue
        families = read(path)
        if not families:
            continue
        n_layers = layer_count(families)
        predicted = 1 - (1 - ALPHA) ** n_layers if n_layers else float("nan")
        primary, p_hits, p_n = rate(families, "exceeds_null_primary")
        cons, c_hits, c_n = rate(families, "exceeds_null_conservative")
        primary_s = "--" if primary is None else f"{p_hits}/{p_n} = {primary:.2f}"
        cons_s = "--" if cons is None else f"{c_hits}/{c_n} = {cons:.2f}"
        print(f"  {directory.name:22s} {n_layers:3d} {predicted:10.2f} "
              f"{primary_s:>14s} {cons_s:>14s}")
    print("\nThe randomly initialized controls have no signal, so every exceedance\n"
          "there is a false positive and the observed rate is the test's true size.")


if __name__ == "__main__":
    main()
