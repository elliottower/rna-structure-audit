"""Why does an untrained model out-score its trained self on the Rung 1 ratio?

Run:  uv run --no-project --with numpy --with scipy --python 3.12 python \
          scripts/audit_trained_untrained_shift.py

Three registered hypotheses predict that training raises the stem-versus-loop
ratio in most families -- H10 for RNA-FM, H15 for RiNALMo, H17 for ERNIE-RNA --
and all three fail in the same direction: the randomly initialized run wins in
most families, and for RiNALMo in every one of them.

Two readings compete. Either the untrained per-family ratios are heavy-tailed, so
a few families carry the comparison and the paired counts are an artifact of
spread; or the ratio is a composition-driven quantity and training moves a model
away from responding to raw nucleotide composition, which puts trained below
untrained systematically. The two predict different things about the paired
values, so the stored per-family ratios settle it: a tail effect leaves the
paired median near 1 and collapses when the tail families are dropped, and a
shift does not.

Every quantity comes from stored per-family values. No model is re-run.
"""

import json
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results/repaired_panel"

PAIRS = ["rinalmo", "ernierna", "rnafm", "splicebert", "utrlm", "nt", "dnabert2"]

DOMAIN = {"rinalmo": "RNA", "ernierna": "RNA", "rnafm": "RNA",
          "splicebert": "RNA", "utrlm": "RNA", "nt": "DNA", "dnabert2": "DNA"}


def ratios(key: str) -> dict[str, float]:
    path = PANEL / key / f"{key}_phases_1_to_5.json"
    per_rna = json.loads(path.read_text())["mutation_trained"]["per_rna"]
    return {name: entry["best_ratio"] for name, entry in per_rna.items()
            if entry.get("best_ratio") is not None}


def main() -> None:
    print("Paired per-family Rung 1 ratios, trained against randomly initialized\n")
    print(f"  {'model':<12}{'pre':<5}{'n':>4}{'trained>':>9}{'median':>9}{'skew':>7}"
          f"{'max/med':>9}{'trimmed':>9}{'p':>9}")
    print(f"  {'':12}{'':5}{'':4}{'untr.':>9}{'t/u':>9}{'untr.':>7}{'untr.':>9}"
          f"{'>untr.':>9}{'':9}")

    for key in PAIRS:
        trained, untrained = ratios(key), ratios(f"{key}_untrained")
        names = sorted(set(trained) & set(untrained))
        t = np.array([trained[n] for n in names])
        u = np.array([untrained[n] for n in names])

        wins = int((t > u).sum())
        live = u > 0
        paired = t[live] / u[live]
        centered = u - u.mean()
        skew = float((centered ** 3).mean() / u.std() ** 3)
        heavy = u >= np.percentile(u, 90)
        trimmed = f"{int((t[~heavy] > u[~heavy]).sum())}/{int((~heavy).sum())}"
        stat = wilcoxon(t, u).pvalue

        zeros = int((~live).sum())
        print(f"  {key:<12}{DOMAIN[key]:<5}{len(names):>4}{wins:>5}/{len(names):<3}"
              f"{np.median(paired):>9.3f}{skew:>7.1f}"
              f"{u.max() / np.median(u):>9.1f}{trimmed:>9}{stat:>9.1e}"
              f"{'  ' + str(zeros) + ' untrained ratios are zero' if zeros else ''}")

    print("\n  Dropping the heaviest decile of untrained ratios moves no paired count\n"
          "  by more than four families, and RiNALMo's stays at zero, so the counts\n"
          "  are not carried by a tail. RiNALMo's untrained ratios have the widest\n"
          "  spread of the five RNA-pretrained models -- a maximum 34 times the\n"
          "  median -- and its trained run still loses in all 47 families.\n")
    print("  The five RNA-pretrained models have paired medians at or below 1 and\n"
          "  the two DNA-pretrained models above it. The split is confounded: the\n"
          "  two models above 1 are the two whose per-position embeddings come out\n"
          "  of the expansion this repository records as D11, so the comparison for\n"
          "  those two waits on the re-run.")


if __name__ == "__main__":
    main()
