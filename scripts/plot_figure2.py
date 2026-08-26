"""Generate Figure 2: Composition confound and raw vs null-controlled ratios.

Panel A: GC composition in stems vs loops across 52 Rfam families.
Panel B: Raw stem-loop sensitivity ratio vs fraction exceeding null, per model.
"""

import json
import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent

# The panel every figure reads. Each of these named a separate July run
# directory at an unrecorded commit, so the figures and the tables could disagree
# about the same model. `repaired_panel_v3` is one commit per cell.
PANEL = ROOT / "results" / "repaired_panel_v3"

MODEL_FILES = {
    "RNA-FM": PANEL / "rnafm" / "rnafm_phases_1_to_5.json",
    "RiNALMo": PANEL / "rinalmo" / "rinalmo_phases_1_to_5.json",
    "ERNIE-RNA": PANEL / "ernierna" / "ernierna_phases_1_to_5.json",
    "SpliceBERT": PANEL / "splicebert" / "splicebert_phases_1_to_5.json",
    "UTR-LM": PANEL / "utrlm" / "utrlm_phases_1_to_5.json",
    "NT v2": PANEL / "nt" / "nt_phases_1_to_5.json",
    "DNABERT-2": PANEL / "dnabert2" / "dnabert2_phases_1_to_5.json",
    "Caduceus": PANEL / "caduceus" / "caduceus_phases_1_to_5.json",
    "HyenaDNA": PANEL / "hyenadna" / "hyenadna_phases_1_to_5.json",
    "Evo": PANEL / "evo" / "evo_phases_1_to_5.json",
}

RNA_MODELS = {"RNA-FM", "RiNALMo", "ERNIE-RNA", "SpliceBERT", "UTR-LM"}
DNA_MODELS = {"NT v2", "DNABERT-2", "Caduceus", "HyenaDNA", "Evo"}


def load_gc_data():
    """Compute stem and loop GC% for each family."""
    fam_files = sorted(glob.glob(str(ROOT / "data/rfam_families/*.json")))
    stems_gc, loops_gc, names = [], [], []
    for fpath in fam_files:
        with open(fpath) as f:
            d = json.load(f)
        seq = d["sequence"].upper()
        db = d["dot_bracket"]
        stem_pos = [i for i, c in enumerate(db) if c in "()"]
        loop_pos = [i for i, c in enumerate(db) if c == "."]
        if not stem_pos or not loop_pos:
            continue
        sgc = sum(1 for i in stem_pos if seq[i] in "GC") / len(stem_pos) * 100
        lgc = sum(1 for i in loop_pos if seq[i] in "GC") / len(loop_pos) * 100
        stems_gc.append(sgc)
        loops_gc.append(lgc)
        names.append(d["name"])
    return np.array(stems_gc), np.array(loops_gc), names


def load_model_data():
    """Load per-model mean ratio and fraction exceeding null."""
    results = {}
    for model, path in MODEL_FILES.items():
        with open(path) as f:
            d = json.load(f)
        pr = d["mutation_trained"]["per_rna"]
        ratios = [v["best_ratio"] for v in pr.values() if "best_ratio" in v]
        exceeds = [v["exceeds_nuc_null"] for v in pr.values() if "exceeds_nuc_null" in v]
        n_total = len(exceeds)
        n_exceed = sum(1 for e in exceeds if e)
        results[model] = {
            "mean_ratio": np.mean(ratios),
            "n_exceed": n_exceed,
            "n_total": n_total,
            "frac_exceed": n_exceed / n_total if n_total else 0,
        }
    return results


def main():
    stems_gc, loops_gc, fam_names = load_gc_data()
    model_data = load_model_data()

    # Sort models by raw ratio descending
    sorted_models = sorted(model_data.keys(), key=lambda m: model_data[m]["mean_ratio"], reverse=True)

    # Colors
    col_rna = "#2166AC"
    col_dna = "#B2182B"
    col_scatter_lo = "#999999"
    col_scatter_hi = "#D6604D"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 3.0))

    # ── Panel A: GC scatter ──
    diff = stems_gc - loops_gc
    norm_diff = (diff - diff.min()) / (diff.max() - diff.min() + 1e-9)
    colors_a = [plt.cm.RdGy_r(0.3 + 0.5 * nd) for nd in norm_diff]

    ax1.scatter(loops_gc, stems_gc, c=colors_a, s=28, edgecolors="0.3",
                linewidths=0.4, zorder=3)
    lims = [10, 100]
    ax1.plot(lims, lims, ls="--", color="0.55", lw=0.8, zorder=1)
    ax1.set_xlim(15, 85)
    ax1.set_ylim(15, 95)
    ax1.set_xlabel("Loop GC content (%)", fontsize=8)
    ax1.set_ylabel("Stem GC content (%)", fontsize=8)
    ax1.tick_params(labelsize=7)

    mean_stem = stems_gc.mean()
    mean_loop = loops_gc.mean()
    ax1.annotate(
        f"Stem {mean_stem:.1f}% vs Loop {mean_loop:.1f}%\n"
        f"(+{mean_stem - mean_loop:.1f} pp)",
        xy=(0.03, 0.97), xycoords="axes fraction", va="top", ha="left",
        fontsize=6.5, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9),
    )
    n_above = np.sum(stems_gc > loops_gc)
    ax1.annotate(
        f"{n_above}/{len(stems_gc)} families above diagonal",
        xy=(0.97, 0.03), xycoords="axes fraction", va="bottom", ha="right",
        fontsize=6, color="0.4",
    )
    ax1.set_title("(A) GC enrichment in stems", fontsize=8.5, fontweight="bold", loc="left")

    # ── Panel B: Simple bar — families exceeding composition null ──
    sorted_by_exceed = sorted(model_data.keys(),
                              key=lambda m: model_data[m]["n_exceed"], reverse=True)
    y_pos = np.arange(len(sorted_by_exceed))
    n_exceeds = [model_data[m]["n_exceed"] for m in sorted_by_exceed]
    n_totals = [model_data[m]["n_total"] for m in sorted_by_exceed]

    bar_colors = [col_rna if m in RNA_MODELS else col_dna for m in sorted_by_exceed]

    ax2.barh(y_pos, n_exceeds, height=0.6, color=bar_colors, alpha=0.85,
             edgecolor="0.3", linewidth=0.4, zorder=2)

    # Expected false-positive count under independence
    expected_fp = 52 * 0.05  # 2.6
    ax2.axvline(expected_fp, color="0.4", ls="--", lw=0.8, zorder=1)

    # Count labels
    for i, ne in enumerate(n_exceeds):
        ax2.text(ne + 0.4, i, str(ne), va="center", ha="left",
                 fontsize=6.5, color="0.3")

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(sorted_by_exceed, fontsize=7)
    ax2.invert_yaxis()
    ax2.set_xlabel("Families exceeding composition null", fontsize=8)
    ax2.set_xlim(0, 38)
    ax2.tick_params(labelsize=7)

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    leg_handles = [
        Patch(facecolor=col_rna, alpha=0.85, edgecolor="0.3", linewidth=0.4,
              label="RNA-pretrained"),
        Patch(facecolor=col_dna, alpha=0.85, edgecolor="0.3", linewidth=0.4,
              label="DNA-pretrained"),
        Line2D([0], [0], color="0.4", ls="--", lw=0.8,
               label=f"Expected FP ({expected_fp:.1f})"),
    ]
    ax2.legend(handles=leg_handles, loc="lower right", fontsize=6,
               framealpha=0.9, edgecolor="0.7")

    ax2.set_title("(B) Families exceeding composition null",
                  fontsize=8.5, fontweight="bold", loc="left")

    fig.tight_layout(w_pad=2.5)

    out_dir = ROOT / "figures"
    out_dir.mkdir(exist_ok=True)
    fig.savefig(out_dir / "figure2_composition.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / "figure2_composition.png", bbox_inches="tight", dpi=300)
    print(f"Saved to {out_dir / 'figure2_composition.pdf'}")
    print(f"Saved to {out_dir / 'figure2_composition.png'}")
    plt.close(fig)


if __name__ == "__main__":
    main()
