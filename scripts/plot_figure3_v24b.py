"""Figure 3: Rung 3 and Rung 4 for the same ten arms.

Two panels on one model axis, every value read from a result file:
  (A) Per-pair precision at the per-family layer (point) against the arm's
      matched within-stem derangement chance rate (bar), for the eight
      character-tokenized models and the two randomly initialized ERNIE-RNA
      arms.
  (B) The matched-decoy rate of every family the registered audit scored, per
      arm, against all decoys; 0.5 is the value of a response shared by every
      position carrying the partner's nucleotide.

Supersedes `plot_figure3_v24a.py`, whose panel A drew PS for three models only.
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
ABLATION = ROOT / "results" / "audit" / "ernierna_ablation_v24.json"

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

layer_models = {m["model"]: m for m in json.loads(LAYERS.read_text())["models"]}
ablated = json.loads(ABLATION.read_text())["arms"]["untrained ablated"]
rung3 = {arm: (layer_models[arm]["per_family_layer"], layer_models[arm]["chance"])
         for _, arm, _ in DECOY_ROWS if arm in layer_models}
rung3["ernierna_untrained_noattnbias"] = (ablated["precision"], ablated["chance_matched"])

for i, (label, arm, domain) in enumerate(DECOY_ROWS):
    precision, chance = rung3[arm]
    col = DOMAIN_COL[domain]
    ax_top.plot([i, i], [chance, precision], color=col, linewidth=1.2, alpha=0.6, zorder=3)
    ax_top.plot([i - 0.28, i + 0.28], [chance] * 2, color="black", linewidth=1.6, zorder=4)
    ax_top.plot(i, precision, "o", color=col, markersize=8,
                markeredgecolor="black", markeredgewidth=0.5, zorder=5)

ax_top.set_ylim(-0.04, 1.04)
ax_top.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_top.set_ylabel("Per-pair precision\n(Rung 3)", fontsize=10)
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
fig.savefig(out_dir / "figure3_rungs_v24b.pdf", bbox_inches="tight")
fig.savefig(out_dir / "figure3_rungs_v24b.png", dpi=300, bbox_inches="tight")
print(f"Saved to {out_dir / 'figure3_rungs_v24b.pdf'}")
print(f"Saved to {out_dir / 'figure3_rungs_v24b.png'}")
plt.close(fig)
