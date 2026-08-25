"""Figure 3: Partner specificity — trained vs randomized.

Two-panel figure:
  (A) Per-family PS distributions for RiNALMo and ERNIE-RNA (linear scale),
      with trained vs randomized controls and bootstrap 95% CIs.
  (B) Mean PS for all 10 models on symlog scale, showing the two-order-of-
      magnitude gap between the positive pair and the rest.
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

def load_ps(path):
    with open(path) as f:
        d = json.load(f)
    pr = d["results"]["per_rna"]
    return np.array([pr[f]["best_ps"] for f in pr if pr[f].get("best_ps") is not None])


PANEL = ROOT / "results" / "repaired_panel_v3"

rinalmo_trained = load_ps(PANEL / "rinalmo" / "rinalmo_phase6_ps.json")
rinalmo_rand = load_ps(PANEL / "rinalmo_untrained" / "rinalmo_untrained_phase6_ps.json")
ernie_trained = load_ps(PANEL / "ernierna" / "ernierna_phase6_ps.json")
ernie_rand = load_ps(PANEL / "ernierna_untrained" / "ernierna_untrained_phase6_ps.json")

# ---------------------------------------------------------------------------
# Mean PS for all models (from paper Table 2)
# ---------------------------------------------------------------------------

all_models = [
    ("RiNALMo",    0.215,    0.152, 0.280, "RNA"),
    ("ERNIE-RNA",  0.120,    0.095, 0.144, "RNA"),
    ("Caduceus",   4.2e-3,   None,  None,  "DNA"),
    ("Evo",        1.1e-3,   None,  None,  "DNA"),
    ("SpliceBERT", 3.5e-4,   None,  None,  "RNA"),
    ("UTR-LM",     1.4e-5,   None,  None,  "RNA"),
    ("HyenaDNA",   9.7e-7,   None,  None,  "DNA"),
    ("RNA-FM",     7.5e-9,   None,  None,  "RNA"),
    ("NT v2",      1e-12,    None,  None,  "DNA"),   # effectively 0
    ("DNABERT-2", -0.016,    None,  None,  "DNA"),
]

rand_controls = {
    "RiNALMo":   8.0e-9,
    "ERNIE-RNA": 2.5e-8,
}

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------

COL_RNA = "#2166ac"      # blue for RNA-pretrained
COL_DNA = "#b2182b"      # red for DNA-pretrained
COL_RAND = "#999999"     # grey for randomized controls
COL_RINALMO = "#1b7837"  # green accent
COL_ERNIE = "#762a83"    # purple accent

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

fig = plt.figure(figsize=(7, 5.5))
gs = fig.add_gridspec(2, 1, height_ratios=[1.1, 1], hspace=0.35)
ax_top = fig.add_subplot(gs[0])
ax_bot = fig.add_subplot(gs[1])

# ---- Panel A: Per-family distributions for the two positive models --------

positions_a = [0, 1]
labels_a = ["RiNALMo", "ERNIE-RNA"]
trained_data = [rinalmo_trained, ernie_trained]
rand_data = [rinalmo_rand, ernie_rand]
cols_a = [COL_RINALMO, COL_ERNIE]

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
        label="Randomized" if i == 0 else None,
    )

ax_top.axhline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.5)
ax_top.set_xticks(positions_a)
ax_top.set_xticklabels(labels_a, fontsize=10)
ax_top.set_ylabel("Partner specificity (PS)", fontsize=10)
ax_top.legend(fontsize=8, loc="upper right", framealpha=0.9)
ax_top.set_xlim(-0.5, 1.7)
ax_top.text(
    -0.08, 1.05, "A", transform=ax_top.transAxes,
    fontsize=13, fontweight="bold", va="top",
)

# ---- Panel B: All 10 models, symlog scale ---------------------------------

names = [m[0] for m in all_models]
ps_vals = np.array([m[1] for m in all_models])
domains = [m[4] for m in all_models]

x_pos = np.arange(len(names))

# Plot |PS| on log scale; annotate negative values
abs_ps = np.abs(ps_vals)
abs_ps[abs_ps == 0] = 1e-13  # floor for NT v2 (PS=0)

for i, (name, ps, aps, domain) in enumerate(
    zip(names, ps_vals, abs_ps, domains)
):
    col = COL_RNA if domain == "RNA" else COL_DNA
    marker = "v" if ps < 0 else "o"
    ax_bot.plot(
        i, aps, marker, color=col, markersize=7,
        markeredgecolor="black", markeredgewidth=0.5, zorder=5,
    )
    if ps < 0:
        ax_bot.annotate(
            "neg.", (i, aps), textcoords="offset points",
            xytext=(0, -14), ha="center", fontsize=7, color=COL_DNA,
        )

# Randomized controls
for name, rand_ps in rand_controls.items():
    idx = names.index(name)
    ax_bot.plot(
        idx, rand_ps, "D", color=COL_RAND, markersize=5,
        markeredgecolor="black", markeredgewidth=0.5, zorder=4,
    )

ax_bot.set_yscale("log")
ax_bot.set_ylim(1e-14, 2)
ax_bot.axhline(1/3, color="black", linewidth=0.8, linestyle=":", alpha=0.6)

ax_bot.set_xticks(x_pos)
ax_bot.set_xticklabels(names, fontsize=8, rotation=35, ha="right")
ax_bot.set_ylabel("|Mean PS|  (log scale)", fontsize=10)

# Legend for domains
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_RNA,
           markeredgecolor="black", markersize=7, label="RNA-pretrained"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_DNA,
           markeredgecolor="black", markersize=7, label="DNA-pretrained"),
    Line2D([0], [0], marker="D", color="w", markerfacecolor=COL_RAND,
           markeredgecolor="black", markersize=5, label="Randomized ctrl"),
]
ax_bot.legend(handles=legend_elements, fontsize=7.5, loc="lower left",
              framealpha=0.9, ncol=1)

ax_bot.text(
    -0.08, 1.05, "B", transform=ax_bot.transAxes,
    fontsize=13, fontweight="bold", va="top",
)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

fig.tight_layout()
out_dir = ROOT / "figures"
out_dir.mkdir(exist_ok=True)
fig.savefig(out_dir / "figure3_partner_specificity.pdf", bbox_inches="tight")
fig.savefig(out_dir / "figure3_partner_specificity.png", dpi=300, bbox_inches="tight")
print(f"Saved to {out_dir / 'figure3_partner_specificity.pdf'}")
print(f"Saved to {out_dir / 'figure3_partner_specificity.png'}")
plt.close(fig)
