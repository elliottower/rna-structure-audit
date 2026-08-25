"""Every H3 verdict as registered, and as corrected, side by side.

`PREREGISTRATION_PHASE6_V2.md` decides H3 by a pooled binomial test of the
partner-is-max fraction against 1/3, "partner is max among three adjacent
positions by chance". That constant assumes the partner and its two stem
neighbours are exchangeable, and they are not: measured on each model's own
derangement, where the pairing is false by construction, the chance rate runs
0.000 to 0.316 (D16).

Printing only the corrected verdict would hide a baseline that was replaced after
the registered one had been seen. Both are printed, with the registered test
computed exactly as registered -- pooled across pairs, against 1/3 -- and the
corrected one against each model's own rate with a family-clustered interval.

    uv run --no-project --with "numpy<2" --with scipy --python 3.12 \\
        python scripts/report_registered_vs_corrected.py
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
REGISTERED_NULL = 1.0 / 3.0
ALPHA = 0.0167          # Bonferroni, three confirmatory hypotheses
N_BOOTSTRAP = 10_000
SEED = 42


def cell(path):
    per_rna = json.loads(path.read_text())["results"]["per_rna"]
    fractions, chances, hits, trials = [], [], 0, 0
    for name, entry in per_rna.items():
        if name in QUARANTINED or entry.get("h3_chance_fraction") is None:
            continue
        precision = entry["h3_precision"]
        fractions.append(precision["fraction"])
        chances.append(entry["h3_chance_fraction"])
        hits += precision["partner_max_count"]
        trials += precision["n"]
    return np.array(fractions), np.array(chances), hits, trials


def excess_interval(fractions, chances, rng):
    """Family-clustered 95% interval on mean precision minus mean chance."""
    draws = np.empty(N_BOOTSTRAP)
    for b in range(N_BOOTSTRAP):
        pick = rng.integers(0, len(fractions), size=len(fractions))
        draws[b] = fractions[pick].mean() - chances[pick].mean()
    return np.percentile(draws, [2.5, 97.5])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=str(REPO / "results/repaired_panel_v3"))
    args = parser.parse_args()
    root = Path(args.results)
    rng = np.random.default_rng(SEED)

    print("  H3 as registered (pooled binomial against 1/3) and as corrected\n"
          "  (against each model's own derangement rate, family-clustered).\n")
    print(f"  {'model':22s} {'precision':>9s} | {'vs 1/3':>8s} {'p':>9s} {'':>4s} | "
          f"{'chance':>7s} {'excess':>7s} {'95% CI':>18s} {'':>4s}")
    for directory in sorted(root.iterdir()):
        path = directory / f"{directory.name}_phase6_ps.json"
        if not path.exists():
            continue
        fractions, chances, hits, trials = cell(path)
        if not len(fractions):
            continue

        registered_p = binomtest(hits, trials, REGISTERED_NULL,
                                 alternative="greater").pvalue
        registered = ("PASS" if (fractions.mean() > REGISTERED_NULL
                                 and registered_p < ALPHA) else "fail")

        low, high = excess_interval(fractions, chances, rng)
        corrected = "PASS" if low > 0 else "fail"

        print(f"  {directory.name:22s} {fractions.mean():9.3f} | "
              f"{fractions.mean() - REGISTERED_NULL:+8.3f} {registered_p:9.2e} "
              f"{registered:>4s} | {chances.mean():7.3f} "
              f"{fractions.mean() - chances.mean():+7.3f} "
              f"[{low:+.3f},{high:+.3f}]{'':>3s} {corrected:>4s}")

    print("\n  A model that changes column is one the registered constant judged\n"
          "  against a rate that is not its own.")


if __name__ == "__main__":
    main()
