"""Figure 3: partner specificity and the matched-decoy test.

Two panels, every value read from a result file:
  (A) Per-family PS for RiNALMo, ERNIE-RNA and RNA-FM over the confirmatory
      families, with the mean, a bootstrap 95% interval, and the mean of each
      model's randomized-weight control.
  (B) The matched-decoy rate of every family the registered audit scored, per
      model, against all decoys; 0.5 is the value of a response shared by every
      position carrying the partner's nucleotide.

Supersedes `plot_figure3_v24.py`, whose panel B drew ten mean-PS values typed
into the script.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Load per-family PS data for models with it
# ---------------------------------------------------------------------------

QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}


def load_ps(path):
    """Best-layer PS of each confirmatory family that was scored."""
    with open(path) as f:
        d = json.load(f)
    pr = d["results"]["per_rna"]
    return np.array([entry["best_ps"] for name, entry in sorted(pr.items())
                     if name not in QUARANTINED and entry.get("best_ps") is not None])


PANEL = ROOT / "results" / "panel_v4"
AUDIT = ROOT / "results" / "audit" / "decoy_audit.json"

# (label, arm, randomized-weight arm), for the three models that clear Rung 3.
POSITIVE = [
    ("RiNALMo", "rinalmo", "rinalmo_untrained"),
    ("ERNIE-RNA", "ernierna", "ernierna_untrained"),
    ("RNA-FM", "rnafm", "rnafm_untrained"),
]
trained_data = [load_ps(PANEL / arm / f"{arm}_phase6_ps.json") for _, arm, _ in POSITIVE]
rand_data = [load_ps(PANEL / ctrl / f"{ctrl}_phase6_ps.json") for _, _, ctrl in POSITIVE]

# Matched-decoy family rates, as the registered audit stored them. Nothing is
# recomputed here: a rate drawn in panel B is a number in the audit file.
DECOY_ROWS = [
    ("RiNALMo", "rinalmo", "RNA"),
    ("ERNIE-RNA", "ernierna", "RNA"),
    ("RNA-FM", "rnafm", "RNA"),
    ("Evo", "evo", "DNA"),
    ("SpliceBERT", "splicebert", "RNA"),
    ("Caduceus", "caduceus", "DNA"),
    ("UTR-LM", "utrlm", "RNA"),
    ("HyenaDNA", "hyenadna", "DNA"),
    ("ERNIE-RNA\nrandom init.", "ernierna_untrained", "CTRL"),
    ("ERNIE-RNA\nrandom init.,\nbias ablated", "ernierna_untrained_noattnbias", "CTRL"),
]
audit_models = {m["model"]: m for m in json.loads(AUDIT.read_text())["models"]}
decoy_rates = {arm: np.array(sorted(audit_models[arm]["family_all"]["family_rates"].values()))
               for _, arm, _ in DECOY_ROWS}

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------

COL_RNA = "#2166ac"      # blue for RNA-pretrained
COL_DNA = "#b2182b"      # red for DNA-pretrained
COL_RAND = "#999999"     # grey for randomized controls
COL_RINALMO = "#1b7837"  # green accent
COL_ERNIE = "#762a83"    # purple accent
COL_RNAFM = "#b35806"    # brown accent

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

fig = plt.figure(figsize=(7, 6.2))
gs = fig.add_gridspec(2, 1, height_ratios=[1.1, 1], hspace=0.35)
ax_top = fig.add_subplot(gs[0])
ax_bot = fig.add_subplot(gs[1])

# ---- Panel A: Per-family distributions for the three positive models --------

positions_a = [0, 1, 2]
labels_a = [label for label, _, _ in POSITIVE]
cols_a = [COL_RINALMO, COL_ERNIE, COL_RNAFM]

for i, (pos, lab, td, rd, col) in enumerate(
    zip(positions_a, labels_a, trained_data, rand_data, cols_a)
):
    parts = ax_top.violinplot(
        td, positions=[pos], showextrema=False, widths=0.6
    )
    for pc in parts["bodies"]:
        pc.set_facecolor(col)
        pc.set_alpha(0.3)
        pc.set_edgecolor(col)

    # Mean + bootstrap CI
    mean_val = td.mean()
    rng = np.random.default_rng(42)
    boots = np.array([rng.choice(td, size=len(td), replace=True).mean() for _ in range(10000)])
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])

    ax_top.errorbar(
        pos, mean_val, yerr=[[mean_val - ci_lo], [ci_hi - mean_val]],
        fmt="o", color=col, markersize=7, capsize=4, capthick=1.5,
        linewidth=1.5, zorder=5, label=f"{lab}",
    )

    # Randomized control
    ax_top.plot(
        pos + 0.25, rd.mean(), "D", color=COL_RAND, markersize=5,
        markeredgecolor="black", markeredgewidth=0.5, zorder=5,
        label="Randomized weights" if i == 0 else None,
    )

ax_top.axhline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.5)
ax_top.set_xticks(positions_a)
ax_top.set_xticklabels(labels_a, fontsize=10)
ax_top.set_ylabel("Partner specificity (PS)", fontsize=10)
ax_top.legend(fontsize=8, loc="upper right", framealpha=0.9)
ax_top.set_xlim(-0.5, 2.7)
ax_top.text(
    -0.08, 1.05, "A", transform=ax_top.transAxes,
    fontsize=13, fontweight="bold", va="top",
)

# ---- Panel B: matched-decoy rate of every family, per model ---------------

rng_jitter = np.random.default_rng(0)
for i, (label, arm, domain) in enumerate(DECOY_ROWS):
    rates = decoy_rates[arm]
    col = {"RNA": COL_RNA, "DNA": COL_DNA, "CTRL": COL_RAND}[domain]
    if len(rates):
        jitter = rng_jitter.uniform(-0.18, 0.18, size=len(rates))
        ax_bot.plot(i + jitter, rates, "o", color=col, markersize=4, alpha=0.75,
                    markeredgecolor="black", markeredgewidth=0.3, zorder=4)
        ax_bot.plot([i - 0.28, i + 0.28], [rates.mean()] * 2, color="black",
                    linewidth=1.6, zorder=5)
    ax_bot.text(i, 1.06, f"n={len(rates)}", ha="center", va="bottom", fontsize=7)

ax_bot.axhline(0.5, color="black", linewidth=0.8, linestyle="--", alpha=0.6)
ax_bot.set_ylim(-0.04, 1.16)
ax_bot.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_bot.set_xticks(np.arange(len(DECOY_ROWS)))
ax_bot.set_xticklabels([label for label, _, _ in DECOY_ROWS], fontsize=7.5,
                       rotation=35, ha="right")
ax_bot.set_ylabel("Matched-decoy rate\n(per family, all decoys)", fontsize=10)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_RNA,
           markeredgecolor="black", markersize=6, label="RNA-pretrained"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_DNA,
           markeredgecolor="black", markersize=6, label="DNA-pretrained"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_RAND,
           markeredgecolor="black", markersize=6, label="Randomized weights"),
    Line2D([0], [0], color="black", linewidth=1.6, label="Mean over families"),
    Line2D([0], [0], color="black", linewidth=0.8, linestyle="--", label="Chance (0.5)"),
]
ax_bot.legend(handles=legend_elements, fontsize=7, loc="lower left",
              framealpha=0.9, ncol=1)

ax_bot.text(
    -0.08, 1.12, "B", transform=ax_bot.transAxes,
    fontsize=13, fontweight="bold", va="top",
)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

fig.tight_layout()
out_dir = ROOT / "figures_v24"
out_dir.mkdir(exist_ok=True)
fig.savefig(out_dir / "figure3_partner_specificity.pdf", bbox_inches="tight")
fig.savefig(out_dir / "figure3_partner_specificity.png", dpi=300, bbox_inches="tight")
print(f"Saved to {out_dir / 'figure3_partner_specificity.pdf'}")
print(f"Saved to {out_dir / 'figure3_partner_specificity.png'}")
plt.close(fig)
