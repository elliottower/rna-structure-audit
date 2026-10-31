"""Figure 2 of the manuscript: Rung 3 and both Rung 4 comparisons for the same ten arms.

Three panels on one model axis, every value read from a result file:
  (A) Held-out per-pair precision (point, with the 2.5th to 97.5th percentile
      over 2,000 half-splits) against the arm's matched chance rate (bar).
  (B) Matched-decoy rate of every family against all same-nucleotide decoys.
  (C) The same against close decoys, those no farther from the mutated
      position than the partner.
Above each column of B and C: families above 0.5 over families not at 0.5,
the count the registered sign test is computed on. An asterisk marks the two
models that carry a registered hypothesis.

Supersedes `plot_figure3_v24d.py`, which drew no close-decoy panel, so RNA-FM
(7 of 9 against close decoys) looked the same as the two registered models.
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
    ("RiNALMo*", "rinalmo", "RNA"),
    ("ERNIE-RNA*", "ernierna", "RNA"),
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
decoy = {kind: {arm: audit_models[arm][f"family_{kind}"] for _, arm, _ in DECOY_ROWS}
         for kind in ("all", "close")}

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------

COL_RNA = "#2166ac"      # blue for RNA-pretrained
COL_DNA = "#b2182b"      # red for DNA-pretrained
COL_RAND = "#999999"     # grey for randomized controls
# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

fig, (ax_top, ax_bot, ax_close) = plt.subplots(
    3, 1, figsize=(7, 8.6), sharex=True,
    gridspec_kw={"height_ratios": [1, 1, 1], "hspace": 0.12})
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
    # The interval is a percentile range over half-splits and need not contain
    # the mean (Caduceus: mean 0.002, range [0, 0]), so it is drawn as its own
    # segment and not as an error bar about the point.
    ax_top.plot([i, i], [low, high], color="black", linewidth=1.0, zorder=4)
    ax_top.plot(i, precision, "o", color=col, markersize=8,
                markeredgecolor="black", markeredgewidth=0.5, zorder=5)

ax_top.set_ylim(-0.04, 1.04)
ax_top.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_top.set_ylabel("Held-out per-pair precision\n(Rung 3)", fontsize=10)
ax_top.text(
    -0.11, 1.04, "A", transform=ax_top.transAxes,
    fontsize=13, fontweight="bold", va="top",
)

# ---- Panels B and C: matched-decoy rate of every family, per arm ---------

rng_jitter = np.random.default_rng(0)
for ax, kind, letter, ylabel in (
        (ax_bot, "all", "B", "Matched-decoy rate\n(Rung 4, all decoys)"),
        (ax_close, "close", "C", "Matched-decoy rate\n(Rung 4, close decoys)")):
    for i, (label, arm, domain) in enumerate(DECOY_ROWS):
        record = decoy[kind][arm]
        rates = np.array(sorted(record["family_rates"].values()))
        col = DOMAIN_COL[domain]
        if len(rates):
            jitter = rng_jitter.uniform(-0.18, 0.18, size=len(rates))
            ax.plot(i + jitter, rates, "o", color=col, markersize=4, alpha=0.75,
                    markeredgecolor="black", markeredgewidth=0.3, zorder=4)
            ax.plot([i - 0.28, i + 0.28], [rates.mean()] * 2, color="black",
                    linewidth=1.6, zorder=5)
        # The sign test's own count: families above 0.5 over families not at 0.5.
        # A column with no family off 0.5 reads "tied" rather than 0/0, which
        # reads as a counting error over points that are present but all at 0.5.
        denom = record["above"] + record["below"]
        count_label = "tied" if denom == 0 else f"{record['above']}/{denom}"
        ax.text(i, 1.05, count_label, ha="center", va="bottom", fontsize=7)
    ax.axhline(0.5, color="black", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.set_ylim(-0.04, 1.16)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel(ylabel, fontsize=10)
    ax.text(-0.11, 1.04, letter, transform=ax.transAxes,
            fontsize=13, fontweight="bold", va="top")

ax_close.set_xticks(np.arange(len(DECOY_ROWS)))
ax_close.set_xticklabels([label for label, _, _ in DECOY_ROWS], fontsize=7.5,
                         rotation=35, ha="right")

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
# No legend box on B/C: the solid bar (mean over families) and the dashed line
# (0.5) are stated in the caption, and any in-panel box covered either the low
# Evo family near 0.17 or the per-column counts along the top.

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

out_dir = ROOT / "figures_v24"
out_dir.mkdir(exist_ok=True)
fig.savefig(out_dir / "figure2_rungs_v24f.pdf", bbox_inches="tight")
fig.savefig(out_dir / "figure2_rungs_v24f.png", dpi=300, bbox_inches="tight")
print(f"Saved to {out_dir / 'figure2_rungs_v24f.pdf'}")
print(f"Saved to {out_dir / 'figure2_rungs_v24f.png'}")
plt.close(fig)
