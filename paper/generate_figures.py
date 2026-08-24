"""Generate paper figures from result JSON files."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

RESULTS_CAUSAL = Path("/Users/elliottower/Documents/GitHub/causal-rna/results")
RESULTS_EXPANDED = Path("/Users/elliottower/Documents/GitHub/rna-structure-awareness/data/expanded_rfam")
OUT_DIR = Path("/Users/elliottower/Documents/GitHub/rna-structure-awareness/paper/figures")
OUT_DIR.mkdir(exist_ok=True)

MODEL_FILES = {
    "ERNIE-RNA": RESULTS_CAUSAL / "ernierna_phases15_dinuc.json",
    "RiNALMo": RESULTS_CAUSAL / "rinalmo_phases_1_to_5.json",
    "RNA-FM": RESULTS_EXPANDED / "rnafm_expanded_20260716.json",
    "UTR-LM": RESULTS_EXPANDED / "utrlm_phases_1_to_5.json",
    "SpliceBERT": RESULTS_CAUSAL / "splicebert_phases15_dinuc.json",
    "NT v2": RESULTS_EXPANDED / "nt_expanded_20260716.json",
    "DNABERT-2": RESULTS_CAUSAL / "dnabert2_phases15_dinuc.json",
    "HyenaDNA": RESULTS_EXPANDED / "hyenadna_expanded_20260716.json",
    "Caduceus": RESULTS_EXPANDED / "caduceus_expanded_20260716.json",
    "Evo": RESULTS_EXPANDED / "evo_expanded_20260716.json",
}

MODEL_DOMAIN = {
    "ERNIE-RNA": "RNA", "RiNALMo": "RNA", "RNA-FM": "RNA",
    "UTR-LM": "RNA", "SpliceBERT": "RNA",
    "NT v2": "DNA", "DNABERT-2": "DNA", "HyenaDNA": "DNA",
    "Caduceus": "DNA", "Evo": "DNA",
}

RUNG3_PS = {
    "RiNALMo": 0.2150, "ERNIE-RNA": 0.1204,
    "Caduceus": 0.0034, "Evo": 0.0011,
    "SpliceBERT": 0.0002, "HyenaDNA": 0.0003,
    "RNA-FM": 0.0001, "UTR-LM": 0.00001,
    "NT v2": 0.001, "DNABERT-2": -0.016,
}

RUNG3_H3 = {
    "RiNALMo": 0.882, "ERNIE-RNA": 0.874,
    "Caduceus": 0.416, "RNA-FM": 0.325,
    "UTR-LM": 0.350, "SpliceBERT": 0.276,
    "HyenaDNA": 0.083, "Evo": 0.278,
}

DINUC_DATA = {
    "ERNIE-RNA": (31, 30), "RiNALMo": (18, 16),
    "NT v2": (14, 11), "HyenaDNA": (10, 8),
    "Evo": (10, 7), "Caduceus": (8, 2),
    "DNABERT-2": (7, 2), "SpliceBERT": (4, 2),
    "UTR-LM": (3, 2), "RNA-FM": (2, 2),
}


def load_per_family(path):
    with open(path) as f:
        data = json.load(f)
    return data["mutation_trained"]["per_rna"]


def fig_overview():
    """Figure 2: Multi-panel results overview."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"wspace": 0.4})

    # Panel A: Dinucleotide retention
    models_sorted = sorted(DINUC_DATA.keys(), key=lambda m: DINUC_DATA[m][1] / max(DINUC_DATA[m][0], 1))
    nuc_counts = [DINUC_DATA[m][0] for m in models_sorted]
    dinuc_counts = [DINUC_DATA[m][1] for m in models_sorted]
    retention = [d / max(n, 1) * 100 for n, d in zip(nuc_counts, dinuc_counts)]

    colors_a = []
    for m in models_sorted:
        if MODEL_DOMAIN[m] == "RNA":
            colors_a.append("#2171b5")
        else:
            colors_a.append("#cb181d")

    y_pos = np.arange(len(models_sorted))
    bars_nuc = ax1.barh(y_pos, nuc_counts, height=0.6, color="#d9d9d9", edgecolor="none", label="Nuc null only")
    bars_dinuc = ax1.barh(y_pos, dinuc_counts, height=0.6, color=[c for c in colors_a], edgecolor="none", label="Survive dinuc")

    for i, (m, n, d) in enumerate(zip(models_sorted, nuc_counts, dinuc_counts)):
        pct = d / max(n, 1) * 100
        ax1.text(n + 0.5, i, f"{d}/{n} ({pct:.0f}%)", va="center", fontsize=8, color="#333333")

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(models_sorted, fontsize=9)
    ax1.set_xlabel("Families exceeding null", fontsize=10)
    ax1.set_title("A.  Dinucleotide null retention", fontsize=11, fontweight="bold", loc="left")
    ax1.set_xlim(0, 40)
    ax1.legend(fontsize=8, loc="lower right")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Perturbation specificity (log scale)
    ps_models = sorted([m for m in RUNG3_PS if RUNG3_PS[m] > 0],
                       key=lambda m: RUNG3_PS[m])
    ps_vals = [RUNG3_PS[m] for m in ps_models]

    colors_b = []
    for m in ps_models:
        if MODEL_DOMAIN[m] == "RNA":
            colors_b.append("#2171b5")
        else:
            colors_b.append("#cb181d")

    y_pos2 = np.arange(len(ps_models))
    ax2.barh(y_pos2, ps_vals, height=0.6, color=colors_b, edgecolor="none")

    for i, (m, v) in enumerate(zip(ps_models, ps_vals)):
        if m in RUNG3_H3:
            ax2.text(v * 1.3, i, f"H3={RUNG3_H3[m]:.2f}", va="center", fontsize=7.5, color="#555555")

    ax2.set_xscale("log")
    ax2.set_yticks(y_pos2)
    ax2.set_yticklabels(ps_models, fontsize=9)
    ax2.set_xlabel("Mean perturbation specificity (PS)", fontsize=10)
    ax2.set_title("B.  Partner specificity (PS)", fontsize=11, fontweight="bold", loc="left")
    ax2.axvline(1/3 * 0.001, color="#999999", linestyle=":", linewidth=0.8)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Legend for domain colors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor="#2171b5", label="RNA-pretrained"),
                       Patch(facecolor="#cb181d", label="DNA-pretrained")]
    ax2.legend(handles=legend_elements, fontsize=8, loc="lower right")

    fig.savefig(OUT_DIR / "fig_overview.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "fig_overview.png", bbox_inches="tight", dpi=300)
    print(f"Saved {OUT_DIR / 'fig_overview.png'}")
    plt.close(fig)


def _build_matrix():
    """Load per-family data and build the ratio matrix."""
    all_per_family = {}
    for model_name, path in MODEL_FILES.items():
        all_per_family[model_name] = load_per_family(path)

    all_families = set()
    for pf in all_per_family.values():
        all_families.update(pf.keys())
    families = sorted(all_families)

    model_order = [
        "ERNIE-RNA", "RiNALMo", "RNA-FM", "UTR-LM", "SpliceBERT",
        "NT v2", "DNABERT-2", "HyenaDNA", "Caduceus", "Evo",
    ]

    matrix = np.full((len(families), len(model_order)), np.nan)
    for j, model in enumerate(model_order):
        pf = all_per_family[model]
        for i, fam in enumerate(families):
            if fam in pf and isinstance(pf[fam], dict):
                matrix[i, j] = pf[fam].get("best_ratio", np.nan)

    # Drop families with no data in any model
    has_data = ~np.all(np.isnan(matrix), axis=1)
    matrix = matrix[has_data]
    families = [f for f, h in zip(families, has_data) if h]

    mean_per_fam = np.nanmean(matrix, axis=1)
    sort_idx = np.argsort(mean_per_fam)[::-1]
    matrix = matrix[sort_idx]
    families = [families[i] for i in sort_idx]

    display_names = []
    for f in families:
        name = f.replace("_", " ")
        if len(name) > 25:
            name = name[:22] + "..."
        display_names.append(name)

    return matrix, model_order, families, display_names


def fig_heatmap_tall():
    """Heatmap with families as rows, models as columns (tall)."""
    matrix, model_order, families, display_names = _build_matrix()

    fig, ax = plt.subplots(figsize=(8, 14))
    log_matrix = np.log2(np.clip(matrix, 0.5, 8.0))
    im = ax.imshow(log_matrix, aspect="auto", cmap="RdYlBu_r",
                   vmin=np.log2(0.7), vmax=np.log2(4.0),
                   interpolation="nearest")

    ax.set_xticks(np.arange(len(model_order)))
    ax.set_xticklabels(model_order, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(len(families)))
    ax.set_yticklabels(display_names, fontsize=6.5)

    for j, model in enumerate(model_order):
        color = "#2171b5" if MODEL_DOMAIN[model] == "RNA" else "#cb181d"
        ax.plot(j, -1.5, "s", color=color, markersize=8, clip_on=False)

    cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label("Mutation sensitivity ratio (log₂)", fontsize=9)
    cbar.set_ticks([np.log2(v) for v in [0.7, 1.0, 1.5, 2.0, 3.0, 4.0]])
    cbar.set_ticklabels(["0.7", "1.0", "1.5", "2.0", "3.0", "4.0"])

    ax.set_title("Per-family mutation sensitivity ratio across models",
                 fontsize=12, fontweight="bold", pad=15)

    fig.savefig(OUT_DIR / "fig_heatmap_tall.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "fig_heatmap_tall.png", bbox_inches="tight", dpi=300)
    print(f"Saved {OUT_DIR / 'fig_heatmap_tall.png'}")
    plt.close(fig)


def fig_heatmap_wide():
    """Heatmap with models as rows, families as columns (wide)."""
    matrix, model_order, families, display_names = _build_matrix()
    matrix_t = matrix.T  # now (models, families)

    fig, ax = plt.subplots(figsize=(18, 5))
    log_matrix = np.log2(np.clip(matrix_t, 0.5, 8.0))
    im = ax.imshow(log_matrix, aspect="auto", cmap="RdYlBu_r",
                   vmin=np.log2(0.7), vmax=np.log2(4.0),
                   interpolation="nearest")

    ax.set_yticks(np.arange(len(model_order)))
    ax.set_yticklabels(model_order, fontsize=9)
    ax.set_xticks(np.arange(len(families)))
    ax.set_xticklabels(display_names, rotation=90, ha="center", fontsize=5.5)

    # Domain color markers on left
    for i, model in enumerate(model_order):
        color = "#2171b5" if MODEL_DOMAIN[model] == "RNA" else "#cb181d"
        ax.plot(-1.5, i, "s", color=color, markersize=8, clip_on=False)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor="#2171b5", label="RNA-pretrained"),
                       Patch(facecolor="#cb181d", label="DNA-pretrained")]
    ax.legend(handles=legend_elements, fontsize=8, loc="upper left",
              bbox_to_anchor=(0, 1.08), ncol=2, frameon=False)

    cbar_ax = fig.add_axes([0.25, 1.02, 0.5, 0.02])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Mutation sensitivity ratio (log₂)", fontsize=9)
    cbar.set_ticks([np.log2(v) for v in [0.7, 1.0, 1.5, 2.0, 3.0, 4.0]])
    cbar.set_ticklabels(["0.7", "1.0", "1.5", "2.0", "3.0", "4.0"])

    fig.suptitle("Per-family mutation sensitivity ratio across models",
                 fontsize=12, fontweight="bold", y=1.1)

    fig.savefig(OUT_DIR / "fig_heatmap_wide.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "fig_heatmap_wide.png", bbox_inches="tight", dpi=300)
    print(f"Saved {OUT_DIR / 'fig_heatmap_wide.png'}")
    plt.close(fig)


def fig_htt():
    """HTT case study: cosine distance vs CAG repeat count + ViennaRNA MFE."""
    htt_dir = RESULTS_CAUSAL / "htt"

    # Load ViennaRNA MFE
    vienna_path = list(htt_dir.glob("htt_viennarna_*/viennarna_mfe.json"))
    if not vienna_path:
        print("WARNING: no ViennaRNA MFE data found, skipping fig_htt")
        return
    with open(vienna_path[0]) as f:
        vienna = json.load(f)
    vienna_counts = sorted(int(k) for k in vienna.keys())
    vienna_mfe = [vienna[str(n)]["mfe"] for n in vienna_counts]

    # Load model HTT distances
    model_results = {}
    for d in sorted(htt_dir.iterdir()):
        if not d.is_dir() or "viennarna" in d.name:
            continue
        json_files = list(d.glob("*_htt_distances.json"))
        if not json_files:
            continue
        with open(json_files[0]) as f:
            data = json.load(f)
        model_name = data.get("model", d.name)
        if "error" in data:
            print(f"  HTT {model_name}: SKIPPED ({data['error'][:60]})")
            continue
        best_layer = str(data.get("best_layer", 0))
        dists = data.get("distances_per_layer", {}).get(best_layer, [])
        counts = data.get("repeat_counts", [])
        if dists and counts:
            model_results[model_name] = (counts, dists, data.get("best_correlation", 0))

    DISPLAY = {
        "ernierna": "ERNIE-RNA", "rinalmo": "RiNALMo", "splicebert": "SpliceBERT",
        "utrlm": "UTR-LM", "nt": "NT v2", "hyenadna": "HyenaDNA",
        "evo": "Evo", "caduceus": "Caduceus", "dnabert2": "DNABERT-2",
    }
    COLORS = {
        "ernierna": "#d62728", "rinalmo": "#1f77b4", "splicebert": "#ff7f0e",
        "utrlm": "#2ca02c", "hyenadna": "#9467bd", "dnabert2": "#8c564b",
        "nt": "#e377c2", "evo": "#7f7f7f", "caduceus": "#bcbd22",
    }

    fig, ax1 = plt.subplots(figsize=(9, 5))

    for model_key in sorted(model_results.keys(), key=lambda m: -model_results[m][2]):
        counts, dists, rho = model_results[model_key]
        label = DISPLAY.get(model_key, model_key)
        color = COLORS.get(model_key, "#333333")
        ax1.plot(counts, dists, "o-", color=color, label=f"{label} (ρ={rho:.2f})",
                 markersize=4, linewidth=1.5, alpha=0.85)

    ax1.set_xlabel("CAG repeat count (N)", fontsize=11)
    ax1.set_ylabel("Cosine distance from baseline (N=10)", fontsize=11)
    ax1.set_title("HTT repeat-expansion embedding distances", fontsize=12, fontweight="bold")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.legend(fontsize=7.5, loc="center left", bbox_to_anchor=(0, 0.55), frameon=False)

    ax2 = ax1.twinx()
    ax2.plot(vienna_counts, vienna_mfe, "k--", linewidth=2, alpha=0.5, label="ViennaRNA MFE")
    ax2.set_ylabel("MFE (kcal/mol)", fontsize=11, color="#555555")
    ax2.tick_params(axis="y", labelcolor="#555555")
    ax2.spines["top"].set_visible(False)
    ax2.legend(fontsize=8, loc="lower right", frameon=False)

    fig.savefig(OUT_DIR / "fig_htt.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "fig_htt.png", bbox_inches="tight", dpi=300)
    print(f"Saved {OUT_DIR / 'fig_htt.png'}")
    plt.close(fig)


def fig_transversion():
    """Transversion vs WC-complement mutation sensitivity comparison."""
    WC_RATIOS = {
        "ERNIE-RNA": 1.242, "RiNALMo": 1.128, "SpliceBERT": 1.040,
        "UTR-LM": 1.054, "HyenaDNA": 1.161, "DNABERT-2": 0.992,
    }

    trans_dir = RESULTS_CAUSAL / "transversion"
    trans_ratios = {}
    KEY_MAP = {
        "ernierna": "ERNIE-RNA", "rinalmo": "RiNALMo", "splicebert": "SpliceBERT",
        "utrlm": "UTR-LM", "hyenadna": "HyenaDNA", "dnabert2": "DNABERT-2",
    }
    for d in sorted(trans_dir.iterdir()):
        if not d.is_dir():
            continue
        json_files = list(d.glob("*_transversion.json"))
        if not json_files:
            continue
        with open(json_files[0]) as f:
            data = json.load(f)
        model_key = data.get("model", "")
        display = KEY_MAP.get(model_key, model_key)
        ratio = data.get("mutation_results", {}).get("mean_best_ratio")
        if ratio and display in WC_RATIOS:
            trans_ratios[display] = ratio

    models = sorted(WC_RATIOS.keys(), key=lambda m: WC_RATIOS[m], reverse=True)
    models = [m for m in models if m in trans_ratios]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(models))
    width = 0.35

    wc_vals = [WC_RATIOS[m] for m in models]
    tr_vals = [trans_ratios[m] for m in models]

    bars1 = ax.bar(x - width/2, wc_vals, width, label="WC complement", color="#4292c6", edgecolor="none")
    bars2 = ax.bar(x + width/2, tr_vals, width, label="Transversion", color="#ef6548", edgecolor="none")

    ax.axhline(1.0, color="#999999", linestyle=":", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(models, fontsize=9)
    ax.set_ylabel("Mean mutation sensitivity ratio", fontsize=10)
    ax.set_title("WC-complement vs transversion mutations", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for i, (wc, tr) in enumerate(zip(wc_vals, tr_vals)):
        pct = (tr - wc) / wc * 100
        sign = "+" if pct >= 0 else ""
        ax.text(i, max(wc, tr) + 0.01, f"{sign}{pct:.1f}%", ha="center", fontsize=7.5, color="#555555")

    fig.savefig(OUT_DIR / "fig_transversion.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "fig_transversion.png", bbox_inches="tight", dpi=300)
    print(f"Saved {OUT_DIR / 'fig_transversion.png'}")
    plt.close(fig)


def fig_ablation():
    """ERNIE-RNA attention bias ablation: trained vs no-bias vs untrained."""
    conditions = ["Trained", "No attn bias", "Untrained"]
    ratios = [1.242, 1.235, 1.848]
    n_exceed_nuc = [31, 32, 3]
    dinuc_retain = [30, 30, None]
    probing = [0.666, 0.666, None]

    fig, axes = plt.subplots(1, 3, figsize=(11, 4), gridspec_kw={"wspace": 0.35})
    colors = ["#2171b5", "#6baed6", "#bdbdbd"]

    # Panel A: Mean ratio
    ax = axes[0]
    bars = ax.bar(conditions, ratios, color=colors, edgecolor="none", width=0.6)
    ax.axhline(1.0, color="#999999", linestyle=":", linewidth=0.8)
    ax.set_ylabel("Mean mutation sensitivity ratio", fontsize=9)
    ax.set_title("A.  Mean ratio", fontsize=10, fontweight="bold", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bar, val in zip(bars, ratios):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.02, f"{val:.3f}",
                ha="center", fontsize=8.5)
    ax.set_ylim(0, 2.1)
    ax.tick_params(axis="x", labelsize=8)

    # Panel B: Families exceeding nulls
    ax = axes[1]
    x = np.arange(3)
    w = 0.35
    ax.bar(x - w/2, n_exceed_nuc, w, label="Nuc null", color=["#2171b5", "#6baed6", "#bdbdbd"], edgecolor="none")
    dinuc_vals = [30, 30, 0]
    ax.bar(x + w/2, dinuc_vals, w, label="Dinuc null", color=["#08519c", "#4292c6", "#969696"], edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(conditions, fontsize=8)
    ax.set_ylabel("Families exceeding null", fontsize=9)
    ax.set_title("B.  Exceedance counts", fontsize=10, fontweight="bold", loc="left")
    ax.legend(fontsize=7.5, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for i, (n, d) in enumerate(zip(n_exceed_nuc, dinuc_vals)):
        ax.text(i - w/2, n + 0.5, str(n), ha="center", fontsize=8)
        if d > 0:
            ax.text(i + w/2, d + 0.5, str(d), ha="center", fontsize=8)

    # Panel C: Probing accuracy
    ax = axes[2]
    probe_vals = [0.666, 0.666, 0.554]
    bars = ax.bar(conditions, probe_vals, color=colors, edgecolor="none", width=0.6)
    ax.axhline(0.5, color="#999999", linestyle=":", linewidth=0.8, label="Chance")
    ax.set_ylabel("Probing balanced accuracy", fontsize=9)
    ax.set_title("C.  Structure probing", fontsize=10, fontweight="bold", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bar, val in zip(bars, probe_vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.005, f"{val:.3f}",
                ha="center", fontsize=8.5)
    ax.set_ylim(0.4, 0.72)
    ax.tick_params(axis="x", labelsize=8)
    ax.legend(fontsize=7.5, loc="lower right", frameon=False)

    fig.savefig(OUT_DIR / "fig_ablation.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "fig_ablation.png", bbox_inches="tight", dpi=300)
    print(f"Saved {OUT_DIR / 'fig_ablation.png'}")
    plt.close(fig)


if __name__ == "__main__":
    fig_overview()
    fig_heatmap_tall()
    fig_heatmap_wide()
    fig_htt()
    fig_transversion()
    fig_ablation()
    print("Done.")
