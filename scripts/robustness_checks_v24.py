"""The robustness figures the manuscript quotes, from the stored Phase 6 files.

Four checks, each computed over the confirmatory families from `results/panel_v4`:

  aggregation   mean against median of per-family best-layer PS
  last layer    mean PS read at the final layer against each family's best layer
  no-op floor   mean PS between two forward passes of one unmutated sequence,
                beside the mean PS of the mutated comparison
  controls      mean PS of each randomized-weight control

Checks that need pair-level responses the runs did not store (a wider neighbor
window, excluding wobble pairs) are not here, and the manuscript does not quote
them.

    uv run --no-project --with numpy python scripts/robustness_checks_v24.py
"""

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results" / "panel_v4"
OUT = REPO / "results" / "audit" / "robustness_v24.json"
GENERATED = REPO / "paper" / "generated_v24"
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# (key, macro-safe name), the ten trained arms.
TRAINED = [("rinalmo", "RiNALMo"), ("ernierna", "ErnieRNA"), ("rnafm", "RNAFM"),
           ("utrlm", "UTRLM"), ("splicebert", "SpliceBERT"), ("nt", "NTvTwo"),
           ("dnabert2", "DNABERTTwo"), ("hyenadna", "HyenaDNA"),
           ("caduceus", "Caduceus"), ("evo", "Evo")]
CONTROLS = ["ernierna", "rinalmo", "rnafm", "utrlm", "splicebert", "nt", "dnabert2"]
POSITIVE = ["rinalmo", "ernierna", "rnafm"]
NEAR_ZERO = 1e-8


def scored(arm: str) -> dict:
    """Confirmatory families of an arm that carry a PS."""
    per_rna = json.loads((PANEL / arm / f"{arm}_phase6_ps.json").read_text())["results"]["per_rna"]
    return {name: entry for name, entry in sorted(per_rna.items())
            if name not in QUARANTINED and isinstance(entry.get("best_ps"), (int, float))}


def sci(value: float) -> str:
    mantissa, exponent = f"{value:.1e}".split("e")
    return rf"{mantissa} \times 10^{{{int(exponent)}}}"


def main() -> int:
    out = {"panel": str(PANEL.relative_to(REPO)), "near_zero": NEAR_ZERO,
           "aggregation": {}, "last_layer": {}, "noop": {}, "controls": {}}
    for key in POSITIVE:
        families = scored(key)
        best = np.array([entry["best_ps"] for entry in families.values()])
        last = np.array([entry["per_layer_ps"][-1] for entry in families.values()])
        out["aggregation"][key] = {
            "n": len(best), "mean": float(best.mean()), "median": float(np.median(best)),
            "median_over_mean_change": float(np.median(best) / best.mean() - 1)}
        out["last_layer"][key] = {
            "n": len(best), "mean_best_layer": float(best.mean()),
            "mean_last_layer": float(last.mean()),
            "change": float(last.mean() / best.mean() - 1)}
    for key, _ in TRAINED:
        families = scored(key)
        noop = [entry["noop_floor"]["best_ps_noop"] for entry in families.values()
                if entry.get("noop_floor") and entry["noop_floor"].get("best_ps_noop") is not None]
        out["noop"][key] = {
            "n": len(families), "n_with_floor": len(noop),
            "mean_ps": float(np.mean([entry["best_ps"] for entry in families.values()])),
            "mean_noop_ps": float(np.mean(noop)) if noop else None,
            "max_abs_noop_ps": float(np.max(np.abs(noop))) if noop else None}
    for key in CONTROLS:
        families = scored(f"{key}_untrained")
        out["controls"][key] = {
            "n": len(families),
            "mean_ps": float(np.mean([entry["best_ps"] for entry in families.values()]))}
    near_zero = [key for key, row in out["controls"].items() if abs(row["mean_ps"]) < NEAR_ZERO]
    out["controls_near_zero"] = near_zero
    quiet = [key for key, row in out["noop"].items()
             if row["max_abs_noop_ps"] is not None and row["max_abs_noop_ps"] < 1e-12]
    out["noop_quiet"] = quiet

    # Written before anything is formatted, so a formatting error cannot lose it.
    OUT.write_text(json.dumps(out, indent=2) + "\n")

    names = dict(TRAINED)
    macros = ["% Written by scripts/robustness_checks_v24.py. Do not edit."]
    for key in POSITIVE:
        macros += [
            rf"\newcommand{{\robustMedianChange{names[key]}}}{{{100 * out['aggregation'][key]['median_over_mean_change']:+.1f}\%}}",
            rf"\newcommand{{\robustLastLayerChange{names[key]}}}{{{100 * out['last_layer'][key]['change']:+.1f}\%}}"]
    for key, name in TRAINED:
        row = out["noop"][key]
        if row["mean_noop_ps"] is not None:
            macros += [rf"\newcommand{{\noopMeanPS{name}}}{{{row['mean_noop_ps']:+.4f}}}",
                       rf"\newcommand{{\noopMutatedPS{name}}}{{{row['mean_ps']:+.4f}}}"]
    for key in CONTROLS:
        macros.append(rf"\newcommand{{\controlMeanPS{names[key]}}}{{{sci(out['controls'][key]['mean_ps'])}}}")
    macros += [rf"\newcommand{{\controlsNearZero}}{{{len(near_zero)}}}",
               rf"\newcommand{{\controlsTotal}}{{{len(CONTROLS)}}}",
               rf"\newcommand{{\noopQuietArms}}{{{len(quiet)}}}", ""]
    (GENERATED / "robustness_macros.tex").write_text("\n".join(macros))

    for section in ("aggregation", "last_layer", "noop", "controls"):
        for key, row in out[section].items():
            print(f"  {section:12s} {key:12s} "
                  + "  ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in row.items()))
    print(f"  controls near zero: {near_zero}\n  quiet no-op arms: {quiet}")
    print(f"  wrote {OUT.relative_to(REPO)} and robustness_macros.tex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
