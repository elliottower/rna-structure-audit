"""Where per-pair precision comes from when the layer is chosen from the data.

Precision at Rung 3 is read at the layer maximizing mean PS, and until now the
families that chose the layer were also the families that reported the number.
Two things follow from that, and they are separable.

The layer is chosen per family: 35 independent argmax draws, one per family,
rather than one draw for the model. A model whose partner specificity sits at a
consistent layer loses nothing when the choice is made once for the whole panel.
A model with no effect gains the maximum of a noise profile in every family. The
size of that gain is not ordered by how many layers there are to draw from, and
the correlation is computed here rather than assumed.

Selection also inflates the number even when the layer is right, because the
same families supply both the argmax and the estimate. Splitting the panel
removes that: one half chooses the layer, the other half reports precision at
it, and the reported number is selection-free by construction.

This computes all three -- precision at each family's own best layer, precision
at one layer chosen for the whole panel, and the held-out split -- against each
model's own within-stem derangement chance rate.

Supersedes `layer_selection_audit.py`, which is left as it was run. That script
compared all three precisions with one chance rate per family, measured at the
family's own best layer and over every eligible pair. Here each precision has
the chance rate of the layer it is read at, over the pairs it is read on, which
needs the `chance_per_layer` field written by
`phase6_compensatory_mutation_v3.py`. A panel without that field is refused.

    uv run --no-project --with "numpy<2" --python 3.12 \
        python scripts/layer_selection_audit_v2.py --panel results/<panel>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "audit" / "layer_selection_v2.json"

# The two pilot families are quarantined from every confirmatory statistic.
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
N_SPLITS = 2000
SEED = 42

MODELS = ["ernierna", "rinalmo", "caduceus", "rnafm", "evo", "splicebert",
          "utrlm", "hyenadna", "ernierna_untrained", "rinalmo_untrained",
          "rnafm_untrained", "splicebert_untrained", "utrlm_untrained"]

LABELS = {
    "ernierna": "ERNIE-RNA",
    "rinalmo": "RiNALMo",
    "caduceus": "Caduceus",
    "rnafm": "RNA-FM",
    "evo": "Evo",
    "splicebert": "SpliceBERT",
    "utrlm": "UTR-LM",
    "hyenadna": "HyenaDNA",
    "ernierna_untrained": "ERNIE-RNA untrained",
    "rinalmo_untrained": "RiNALMo untrained",
    "rnafm_untrained": "RNA-FM untrained",
    "splicebert_untrained": "SpliceBERT untrained",
    "utrlm_untrained": "UTR-LM untrained",
}


class PanelMismatch(RuntimeError):
    """A stored number and its recomputation disagree."""


def load(panel: Path, key: str) -> dict:
    """Per-layer PS, per-layer precision and per-layer chance rate for each family.

    Only families carrying a matched derangement chance rate at every layer are
    kept, which is the set every precision figure in the manuscript is taken
    over.
    """
    stored = json.loads((panel / key / f"{key}_phase6_ps.json").read_text())
    per_rna = stored["results"]["per_rna"]
    if not any("chance_per_layer" in entry for entry in per_rna.values()):
        raise PanelMismatch(
            f"{key}: no family carries chance_per_layer. This panel was written "
            "by phase6_compensatory_mutation.py; rerun with the v3 script.")

    names, ps, precision, chance, own, best = [], [], [], [], [], []
    excluded = {}
    # Sorted, so two arms with the same families are split the same way.
    for name in sorted(per_rna):
        entry = per_rna[name]
        if name in QUARANTINED:
            excluded[name] = "quarantined pilot family"
            continue
        if entry.get("skipped"):
            excluded[name] = entry.get("reason", "skipped")
            continue
        layers = entry.get("precision_per_layer")
        if not layers or any(v is None for v in layers):
            excluded[name] = "no deep-interior pair at some layer"
            continue
        per_layer_chance = entry.get("chance_per_layer")
        if not per_layer_chance or any(v is None for v in per_layer_chance):
            # Precision exists and its chance rate does not: the family would
            # vanish from the comparison without a trace.
            raise PanelMismatch(
                f"{key}/{name}: precision_per_layer is present and "
                "chance_per_layer is missing or incomplete")
        if len(per_layer_chance) != len(layers):
            raise PanelMismatch(
                f"{key}/{name}: {len(per_layer_chance)} chance layers against "
                f"{len(layers)} precision layers")
        if len(layers) != len(entry["per_layer_ps"]):
            raise PanelMismatch(
                f"{key}/{name}: {len(layers)} precision layers against "
                f"{len(entry['per_layer_ps'])} PS layers")
        # The stored fraction is precision at this family's own best layer.
        # If the recomputation disagrees, the per-layer array is not the array
        # the reported number came from, and nothing below means anything.
        recomputed = layers[entry["best_layer"]]
        if abs(recomputed - entry["h3_precision"]["fraction"]) > 1e-9:
            raise PanelMismatch(
                f"{key}/{name}: stored precision {entry['h3_precision']['fraction']} "
                f"against {recomputed} at layer {entry['best_layer']}")
        names.append(name)
        ps.append(entry["per_layer_ps"])
        precision.append(layers)
        chance.append(per_layer_chance)
        own.append(recomputed)
        best.append(entry["best_layer"])

    if not names:
        raise PanelMismatch(f"{key}: no family survives the filters")
    ps_array = np.array(ps, dtype=float)
    if np.isnan(ps_array).all(axis=0).any():
        raise PanelMismatch(f"{key}: a layer has no finite PS in any family")
    return {"names": names, "ps": ps_array,
            "precision": np.array(precision, dtype=float),
            "chance": np.array(chance, dtype=float),
            "own": np.array(own, dtype=float),
            "best": np.array(best, dtype=int),
            "excluded": excluded}


def heldout(ps: np.ndarray, precision: np.ndarray, chance: np.ndarray,
            rng) -> tuple[np.ndarray, np.ndarray]:
    """Precision, and its chance rate, on the half that did not choose the layer."""
    n = len(ps)
    draws = np.empty(N_SPLITS)
    chance_draws = np.empty(N_SPLITS)
    for i in range(N_SPLITS):
        order = rng.permutation(n)
        chooser, reporter = order[: n // 2], order[n // 2:]
        layer = int(np.nanmean(ps[chooser], axis=0).argmax())
        draws[i] = precision[reporter, layer].mean()
        chance_draws[i] = chance[reporter, layer].mean()
    return draws, chance_draws


def audit(panel: Path, key: str) -> dict:
    data = load(panel, key)
    # One generator per arm from one seed: arms with the same sorted families
    # receive the same sequence of half-splits.
    rng = np.random.default_rng(SEED)
    ps, precision, chance = data["ps"], data["precision"], data["chance"]
    families = np.arange(len(ps))
    global_layer = int(np.nanmean(ps, axis=0).argmax())
    draws, chance_draws = heldout(ps, precision, chance, rng)
    low, high = (float(v) for v in np.percentile(draws, [2.5, 97.5]))
    return {
        "model": key,
        "label": LABELS[key],
        "n_families": len(data["names"]),
        "families": data["names"],
        "excluded": data["excluded"],
        "n_layers": int(ps.shape[1]),
        # Each chance rate is read at the layer its precision is read at.
        "chance": float(chance[families, data["best"]].mean()),
        "chance_global_layer": float(chance[:, global_layer].mean()),
        "chance_heldout": float(chance_draws.mean()),
        "per_family_layer": float(data["own"].mean()),
        "global_layer": float(precision[:, global_layer].mean()),
        "global_layer_index": global_layer,
        "heldout": float(draws.mean()),
        "heldout_ci": [low, high],
        "n_splits": N_SPLITS,
    }


def depth_correlation(rows: list[dict]) -> dict:
    """Whether the per-family choice buys more in a model with more layers.

    Among the randomly initialized character-level controls alone the excess
    over chance does rise with depth, which is where the expectation came from.
    Across the whole panel it does not, so the gap between the two readings
    cannot be attributed to the number of layers searched.
    """
    layers = [row["n_layers"] for row in rows]
    gain = [row["per_family_layer"] - row["global_layer"] for row in rows]
    rho, p = spearmanr(layers, gain)
    return {"n_arms": len(rows), "rho": float(rho), "p": float(p)}


# LaTeX control sequences take letters only, so each arm gets a spelled-out
# name. Every arm is emitted whether or not the prose currently quotes it: a
# number typed as a digit is one a later edit can leave stale, and several of
# these collide with unrelated quantities elsewhere in the manuscript.
SUFFIX = {
    "ernierna": "ErnieRNA",
    "rinalmo": "RiNALMo",
    "caduceus": "Caduceus",
    "rnafm": "RNAFM",
    "evo": "Evo",
    "splicebert": "SpliceBERT",
    "utrlm": "UTRLM",
    "hyenadna": "HyenaDNA",
}


def macro_name(key: str) -> str:
    base, _, untrained = key.partition("_")
    return SUFFIX[base] + ("Untrained" if untrained else "")


def macros(depth: dict, rows: list[dict]) -> str:
    """Preamble macros, so no number in the prose can be typed as a stale digit."""
    lines = [
        r"% Written by scripts/layer_selection_audit_v2.py. Do not edit.",
        rf"\newcommand{{\selectionDepthRho}}{{{depth['rho']:+.2f}}}",
        rf"\newcommand{{\selectionDepthP}}{{{depth['p']:.2f}}}",
        rf"\newcommand{{\selectionArms}}{{{depth['n_arms']}}}",
    ]
    for row in rows:
        name = macro_name(row["model"])
        low, high = row["heldout_ci"]
        lines += [
            rf"\newcommand{{\selLayers{name}}}{{{row['n_layers']}}}",
            rf"\newcommand{{\selChance{name}}}{{{row['chance']:.3f}}}",
            rf"\newcommand{{\selChancePanel{name}}}"
            rf"{{{row['chance_global_layer']:.3f}}}",
            rf"\newcommand{{\selChanceHeldOut{name}}}"
            rf"{{{row['chance_heldout']:.3f}}}",
            rf"\newcommand{{\selPerFamily{name}}}{{{row['per_family_layer']:.3f}}}",
            rf"\newcommand{{\selPanel{name}}}{{{row['global_layer']:.3f}}}",
            rf"\newcommand{{\selHeldOut{name}}}{{{row['heldout']:.3f}}}",
            rf"\newcommand{{\selHeldOutCI{name}}}{{[{low:.3f}, {high:.3f}]}}",
            # Excess over the arm's own chance rate, which the prose quotes
            # directly rather than leaving the reader to subtract.
            rf"\newcommand{{\selExcess{name}}}"
            rf"{{{row['per_family_layer'] - row['chance']:.3f}}}",
            # Distance from the exchangeability value the manuscript contrasts
            # each model's own chance rate against.
            rf"\newcommand{{\selVsThird{name}}}"
            rf"{{{abs(row['per_family_layer'] - 1 / 3):.3f}}}",
        ]
    # The largest distance the panel-layer repair moves either model that
    # clears Rung 3, quoted wherever the repair is said to cost them nothing.
    positive = [row for row in rows if row["model"] in ("rinalmo", "ernierna")]
    shift = max(abs(row["per_family_layer"] - row["global_layer"])
                for row in positive)
    lines.append(rf"\newcommand{{\selMaxShiftPositive}}{{{shift:.3f}}}")
    lines.append("")
    return "\n".join(lines)


def latex(rows: list[dict]) -> str:
    """Supplementary table: the three ways of choosing a layer, side by side."""
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Per-pair precision under three layer choices, each against "
        r"the model's own within-stem derangement chance rate at the same layer "
        r"and on the same pairs (in parentheses). Per-family: the "
        r"layer maximizing mean PS within each family, which is how "
        r"Table~\ref{tab:rung3} is computed. Panel: one layer maximizing mean "
        r"PS across all families. Held-out: the layer chosen on half the "
        r"families and precision reported on the other half, over "
        rf"{N_SPLITS:,} splits; the bracketed interval is the 2.5th to 97.5th "
        r"percentile over random half-splits and is not a confidence interval "
        r"over families.}",
        r"\label{tab:layer-selection}",
        r"\small",
        r"\begin{tabular}{@{}lrrrr@{}}",
        r"\toprule",
        r"\textbf{Model} & \textbf{Layers} & "
        r"\textbf{Per-family} & \textbf{Panel} & \textbf{Held-out} \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            f"{row['label']} & {row['n_layers']} & "
            f"{row['per_family_layer']:.3f} ({row['chance']:.3f}) & "
            f"{row['global_layer']:.3f} ({row['chance_global_layer']:.3f}) & "
            f"{row['heldout']:.3f} ({row['chance_heldout']:.3f}) "
            f"[{row['heldout_ci'][0]:.3f}, {row['heldout_ci'][1]:.3f}] \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--panel", type=Path, required=True,
                        help="Panel directory written by the v3 phase 6 script")
    args = parser.parse_args()

    rows = [audit(args.panel, key) for key in MODELS]
    depth = depth_correlation(rows)

    # Written before anything is printed, so a formatting error cannot lose it.
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"n_splits": N_SPLITS, "seed": SEED,
                               "panel": str(args.panel),
                               "quarantined": sorted(QUARANTINED),
                               "depth_correlation": depth,
                               "models": rows}, indent=2))

    print(f"  {'model':22s} {'layers':>6s} {'per-family':>11s} {'panel':>7s} "
          f"{'held-out':>9s}   excess over the chance rate at the same layer")
    for row in rows:
        print(f"  {row['label']:22s} {row['n_layers']:6d} "
              f"{row['per_family_layer'] - row['chance']:+11.3f} "
              f"{row['global_layer'] - row['chance_global_layer']:+7.3f} "
              f"{row['heldout'] - row['chance_heldout']:+9.3f}")
    print(f"\n  depth against per-family gain, {depth['n_arms']} arms: "
          f"Spearman rho = {depth['rho']:+.3f}, p = {depth['p']:.3f}")

    generated = REPO / "paper" / "generated"
    table = generated / "layer_selection_table_v2.tex"
    table.write_text(latex(rows))
    macro_file = generated / "layer_selection_macros_v2.tex"
    macro_file.write_text(macros(depth, rows))
    print(f"\n  wrote {OUT.relative_to(REPO)}, {table.relative_to(REPO)} "
          f"and {macro_file.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
