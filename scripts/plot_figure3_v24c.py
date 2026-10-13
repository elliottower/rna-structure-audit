"""Figure 3: Rung 3 and Rung 4 for the same ten arms.

Two panels on one model axis, every value read from a result file:
  (A) Per-pair precision with the layer chosen on half the families and
      scored on the other half (point, with the 2.5th to 97.5th percentile over
      2,000 half-splits) against the arm's matched within-stem derangement
      chance rate at that layer (bar), for the eight character-tokenized models
      and the two randomly initialized ERNIE-RNA arms.
  (B) The matched-decoy rate of every family the registered audit scored, per
      arm, against all decoys; 0.5 is the value of a response shared by every
      position carrying the partner's nucleotide.

Supersedes `plot_figure3_v24b.py`, whose panel A chose the layer within each
family, the choice under which Caduceus reads 0.340 against 0.002 held out.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent

AUDIT = ROOT / "results" / "audit" / "decoy_audit.json"
LAYERS = ROOT / "results" / "audit" / "layer_selection_v2.json"
ABLATED_ARM = ROOT / "results" / "audit" / "layer_selection_v3_ablated_arm.json"

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
# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(7, 6.4), sharex=True,
    gridspec_kw={"height_ratios": [1, 1], "hspace": 0.12})
DOMAIN_COL = {"RNA": COL_RNA, "DNA": COL_DNA, "CTRL": COL_RAND}

# ---- Panel A: per-pair precision against matched chance, per arm ----------

layer_models = {m["model"]: m for path in (LAYERS, ABLATED_ARM)
                for m in json.loads(path.read_text())["models"]}

for i, (label, arm, domain) in enumerate(DECOY_ROWS):
    row = layer_models[arm]
    precision, chance = row["heldout"], row["chance_heldout"]
    low, high = row["heldout_ci"]
    col = DOMAIN_COL[domain]
    ax_top.plot([i - 0.28, i + 0.28], [chance] * 2, color="black", linewidth=1.6, zorder=4)
    ax_top.errorbar(i, precision, yerr=[[precision - low], [high - precision]],
                    fmt="o", color=col, markersize=8, markeredgecolor="black",
                    markeredgewidth=0.5, ecolor="black", elinewidth=1.0, capsize=3,
                    zorder=5)

ax_top.set_ylim(-0.04, 1.04)
ax_top.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_top.set_ylabel("Held-out per-pair precision\n(Rung 3)", fontsize=10)
ax_top.text(
    -0.11, 1.04, "A", transform=ax_top.transAxes,
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
    ax_bot.text(i, 1.05, f"n={len(rates)}", ha="center", va="bottom", fontsize=7)

ax_bot.axhline(0.5, color="black", linewidth=0.8, linestyle="--", alpha=0.6)
ax_bot.set_ylim(-0.04, 1.16)
ax_bot.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_bot.set_xticks(np.arange(len(DECOY_ROWS)))
ax_bot.set_xticklabels([label for label, _, _ in DECOY_ROWS], fontsize=7.5,
                       rotation=35, ha="right")
ax_bot.set_ylabel("Matched-decoy rate\n(Rung 4)", fontsize=10)

legend_top = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_RNA,
           markeredgecolor="black", markersize=6, label="RNA-pretrained"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_DNA,
           markeredgecolor="black", markersize=6, label="DNA-pretrained"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_RAND,
           markeredgecolor="black", markersize=6, label="Randomly initialized"),
    Line2D([0], [0], color="black", linewidth=1.6, label="Matched chance rate"),
]
ax_top.legend(handles=legend_top, fontsize=7, loc="upper right", framealpha=0.9)
legend_bot = [
    Line2D([0], [0], color="black", linewidth=1.6, label="Mean over families"),
    Line2D([0], [0], color="black", linewidth=0.8, linestyle="--", label="Shared response (0.5)"),
]
ax_bot.legend(handles=legend_bot, fontsize=7, loc="lower left", framealpha=0.9)

ax_bot.text(
    -0.11, 1.04, "B", transform=ax_bot.transAxes,
    fontsize=13, fontweight="bold", va="top",
)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

out_dir = ROOT / "figures_v24"
out_dir.mkdir(exist_ok=True)
fig.savefig(out_dir / "figure3_rungs_v24c.pdf", bbox_inches="tight")
fig.savefig(out_dir / "figure3_rungs_v24c.png", dpi=300, bbox_inches="tight")
print(f"Saved to {out_dir / 'figure3_rungs_v24c.pdf'}")
print(f"Saved to {out_dir / 'figure3_rungs_v24c.png'}")
plt.close(fig)
