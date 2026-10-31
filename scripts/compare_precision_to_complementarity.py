"""Stored Rung 3 precision against what nucleotide identity alone can reach.

A representation that moves every position carrying the partner's nucleotide
equally wins the Rung 3 comparison only where both stem neighbours of the
partner carry another nucleotide. `results/audit/decoy_control_validation.json`
holds that precision for every family. This sets each model's stored per-pair
precision beside it, family by family, at one layer per model: the layer
maximizing mean PS over the model's families, so no family picks its own.

Descriptive and after the fact: every model outcome read here was already
reported. No forward pass is run.

    uv run --no-project python scripts/compare_precision_to_complementarity.py
"""

import json
import math
import statistics
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results" / "repaired_panel_v3"
CONTROL = REPO / "results" / "audit" / "decoy_control_validation.json"
OUT = REPO / "results" / "audit" / "precision_vs_complementarity.json"

QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
MODELS = ["rinalmo", "ernierna", "caduceus", "evo", "splicebert", "utrlm", "hyenadna",
          "ernierna_untrained", "rinalmo_untrained"]


def panel_layer(entries: dict) -> int:
    n_layers = len(next(iter(entries.values()))["per_layer_ps"])
    means = []
    for layer in range(n_layers):
        column = [entry["per_layer_ps"][layer] for entry in entries.values()
                  if not math.isnan(entry["per_layer_ps"][layer])]
        means.append(sum(column) / len(column) if column else float("-inf"))
    return max(range(n_layers), key=means.__getitem__)


def compare(key: str, ceiling: dict) -> dict:
    stored = json.loads((PANEL / key / f"{key}_phase6_ps.json").read_text())
    entries = {name: entry for name, entry in stored["results"]["per_rna"].items()
               if name not in QUARANTINED and not entry.get("skipped")
               and entry.get("precision_per_layer")
               and all(v is not None for v in entry["precision_per_layer"])
               and ceiling.get(name) is not None}
    layer = panel_layer(entries)
    rows = {name: {"model": entry["precision_per_layer"][layer],
                   "complementarity": ceiling[name]}
            for name, entry in sorted(entries.items())}
    above = sum(row["model"] > row["complementarity"] for row in rows.values())
    below = sum(row["model"] < row["complementarity"] for row in rows.values())
    n = above + below
    p = sum(math.comb(n, k) for k in range(above, n + 1)) / 2 ** n if n else 1.0
    return {
        "model": key, "panel_layer": layer, "families": len(rows),
        "mean_model_precision": statistics.fmean(r["model"] for r in rows.values()),
        "mean_complementarity_precision": statistics.fmean(
            r["complementarity"] for r in rows.values()),
        "families_above": above, "families_below": below,
        "families_equal": len(rows) - n, "sign_test_p_one_sided": p,
        "per_family": rows,
    }


def main() -> int:
    control = json.loads(CONTROL.read_text())
    ceiling = {name: row["complementarity"]["rung3_precision"]
               for name, row in control["families"].items()}
    rows = [compare(key, ceiling) for key in MODELS]
    OUT.write_text(json.dumps({"panel": str(PANEL.relative_to(REPO)),
                               "control": str(CONTROL.relative_to(REPO)),
                               "models": rows}, indent=2))
    for row in rows:
        print(f"  {row['model']:22s} layer {row['panel_layer']:2d}  "
              f"precision {row['mean_model_precision']:.3f} vs "
              f"{row['mean_complementarity_precision']:.3f}  "
              f"above/below/equal {row['families_above']}/{row['families_below']}/"
              f"{row['families_equal']}  p={row['sign_test_p_one_sided']:.2g}")
    print(f"wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
