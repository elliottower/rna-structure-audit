"""Where per-pair precision comes from when the layer is chosen from the data.

Precision at Rung 3 is read at the layer maximizing mean PS, and until now the
families that chose the layer were also the families that reported the number.
Two things follow from that, and they are separable.

The layer is chosen per family: 35 independent argmax draws, one per family, not
one draw for the model. A model whose partner specificity sits at a consistent
layer loses nothing when the choice is made once for the whole panel. A model
with no effect gains the maximum of a noise profile in every family, and the
gain grows with the number of layers it can draw from.

Selection also inflates the number even when the layer is right, because the
same families supply both the argmax and the estimate. Splitting the panel
removes that: one half chooses the layer, the other half reports precision at
it, and the reported number is selection-free by construction.

This computes all three -- precision at each family's own best layer, precision
at one layer chosen for the whole panel, and the held-out split -- against each
model's own within-stem derangement chance rate.

    uv run --no-project --with "numpy<2" --python 3.12 \
        python scripts/layer_selection_audit.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results" / "repaired_panel_v3"
OUT = REPO / "results" / "audit" / "layer_selection.json"

# The two pilot families are quarantined from every confirmatory statistic.
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
N_SPLITS = 2000
SEED = 42

MODELS = ["ernierna", "rinalmo", "ernierna_untrained", "rinalmo_untrained",
          "rnafm_untrained", "splicebert_untrained", "utrlm_untrained"]

LABELS = {
    "ernierna": "ERNIE-RNA",
    "rinalmo": "RiNALMo",
    "ernierna_untrained": "ERNIE-RNA untrained",
    "rinalmo_untrained": "RiNALMo untrained",
    "rnafm_untrained": "RNA-FM untrained",
    "splicebert_untrained": "SpliceBERT untrained",
    "utrlm_untrained": "UTR-LM untrained",
}


class PanelMismatch(RuntimeError):
    """A stored number and its recomputation disagree."""


def load(key: str) -> dict:
    """Per-layer PS, per-layer precision and chance rate for each family.

    Only families carrying a derangement chance rate are kept, which is the set
    every precision figure in the manuscript is taken over.
    """
    stored = json.loads((PANEL / key / f"{key}_phase6_ps.json").read_text())
    per_rna = stored["results"]["per_rna"]

    names, ps, precision, chance, own = [], [], [], [], []
    for name, entry in per_rna.items():
        if name in QUARANTINED or entry.get("h3_chance_fraction") is None:
            continue
        layers = entry.get("precision_per_layer")
        if not layers or any(v is None for v in layers):
            continue
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
        chance.append(entry["h3_chance_fraction"])
        own.append(recomputed)

    return {"names": names, "ps": np.array(ps, dtype=float),
            "precision": np.array(precision, dtype=float),
            "chance": np.array(chance, dtype=float),
            "own": np.array(own, dtype=float)}


def heldout(ps: np.ndarray, precision: np.ndarray, rng) -> np.ndarray:
    """Precision on the half of the panel that did not choose the layer."""
    n = len(ps)
    draws = np.empty(N_SPLITS)
    for i in range(N_SPLITS):
        order = rng.permutation(n)
        chooser, reporter = order[: n // 2], order[n // 2:]
        layer = int(np.nanmean(ps[chooser], axis=0).argmax())
        draws[i] = precision[reporter, layer].mean()
    return draws


def audit(key: str, rng) -> dict:
    data = load(key)
    ps, precision, chance = data["ps"], data["precision"], data["chance"]
    global_layer = int(np.nanmean(ps, axis=0).argmax())
    draws = heldout(ps, precision, rng)
    low, high = (float(v) for v in np.percentile(draws, [2.5, 97.5]))
    return {
        "model": key,
        "label": LABELS[key],
        "n_families": len(data["names"]),
        "n_layers": int(ps.shape[1]),
        "chance": float(chance.mean()),
        "per_family_layer": float(data["own"].mean()),
        "global_layer": float(precision[:, global_layer].mean()),
        "global_layer_index": global_layer,
        "heldout": float(draws.mean()),
        "heldout_ci": [low, high],
        "n_splits": N_SPLITS,
    }


def latex(rows: list[dict]) -> str:
    """Supplementary table: the three ways of choosing a layer, side by side."""
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Per-pair precision under three layer choices, each against "
        r"the model's own within-stem derangement chance rate. Per-family: the "
        r"layer maximizing mean PS within each family, which is how "
        r"Table~\ref{tab:rung3} is computed. Panel: one layer maximizing mean "
        r"PS across all families. Held-out: the layer chosen on half the "
        r"families and precision reported on the other half, over "
        rf"{N_SPLITS:,} splits, with a 95\% interval over splits. "
        r"Choosing per family costs the two trained models nothing and buys "
        r"the randomly initialized controls up to $+0.18$.}",
        r"\label{tab:layer-selection}",
        r"\small",
        r"\begin{tabular}{@{}lrrrrr@{}}",
        r"\toprule",
        r"\textbf{Model} & \textbf{Layers} & \textbf{Chance} & "
        r"\textbf{Per-family} & \textbf{Panel} & \textbf{Held-out} \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            f"{row['label']} & {row['n_layers']} & {row['chance']:.3f} & "
            f"{row['per_family_layer']:.3f} & {row['global_layer']:.3f} & "
            f"{row['heldout']:.3f} "
            f"[{row['heldout_ci'][0]:.3f}, {row['heldout_ci'][1]:.3f}] \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def main() -> int:
    rng = np.random.default_rng(SEED)
    rows = [audit(key, rng) for key in MODELS]

    header = (f"  {'model':22s} {'layers':>6s} {'chance':>7s} {'per-family':>11s} "
              f"{'panel':>7s} {'held-out':>19s}")
    print(header)
    for row in rows:
        print(f"  {row['label']:22s} {row['n_layers']:6d} {row['chance']:7.3f} "
              f"{row['per_family_layer']:11.3f} {row['global_layer']:7.3f} "
              f"{row['heldout']:7.3f} "
              f"[{row['heldout_ci'][0]:.3f},{row['heldout_ci'][1]:.3f}]")

    print("\n  excess over each model's own chance rate")
    print(f"  {'model':22s} {'per-family':>11s} {'panel':>7s} {'held-out':>9s}")
    for row in rows:
        print(f"  {row['label']:22s} "
              f"{row['per_family_layer'] - row['chance']:+11.3f} "
              f"{row['global_layer'] - row['chance']:+7.3f} "
              f"{row['heldout'] - row['chance']:+9.3f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"n_splits": N_SPLITS, "seed": SEED,
                               "quarantined": sorted(QUARANTINED),
                               "models": rows}, indent=2))
    table = REPO / "paper" / "generated" / "layer_selection_table.tex"
    table.write_text(latex(rows))
    print(f"\n  wrote {OUT.relative_to(REPO)} and {table.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
