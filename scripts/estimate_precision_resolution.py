"""What difference in per-pair precision can this panel actually resolve?

A registered margin should be a number the design can deliver, not a round one.
Per-pair precision is a family-clustered proportion, so the resolution is set by
the between-family variance, which is measurable from any run already on disk.

The registered quantity is not a difference of precisions but a difference of
*excesses over per-arm chance rates*, which is four estimated numbers rather than
two:

    (p_intact - c_intact) - (p_ablated - c_ablated)

Each derangement rate carries its own sampling error over the same families, and
because the deltas feeding the null are the deltas feeding the observed
precision, those errors are correlated by an amount whose sign is not guessable.
So the full statistic is bootstrapped directly -- families are resampled once per
draw and all four quantities recomputed inside that resample -- rather than
combining separately estimated half-widths.

    uv run --no-project --with "numpy<2" --python 3.12 \
        python scripts/estimate_precision_resolution.py
"""

import argparse
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
N_BOOTSTRAP = 10_000
SEED = 42


def per_family(path):
    """Per-family observed precision and derangement chance rate, aligned by name.

    Both come from the same family, so they must be resampled together: a draw
    that takes a family's precision and another family's chance rate would
    break the correlation the bootstrap exists to capture.
    """
    per_rna = json.loads(path.read_text())["results"]["per_rna"]
    names, observed, chance = [], [], []
    for name, entry in per_rna.items():
        if name in QUARANTINED or not isinstance(entry.get("h3_precision"), dict):
            continue
        if entry.get("h3_chance_fraction") is None:
            continue
        names.append(name)
        observed.append(entry["h3_precision"]["fraction"])
        chance.append(entry["h3_chance_fraction"])
    return names, np.array(observed), np.array(chance)


def resample_excess(observed, chance, index):
    """Mean excess over chance on one resample of families."""
    return observed[index].mean() - chance[index].mean()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=str(REPO / "results/repaired_panel_v3"))
    parser.add_argument("--models", default="ernierna,ernierna_untrained")
    args = parser.parse_args()

    rng = np.random.default_rng(SEED)
    root = Path(args.results)
    models = [m.strip() for m in args.models.split(",") if m.strip()]

    arms = {}
    for model in models:
        names, observed, chance = per_family(
            root / model / f"{model}_phase6_ps.json")
        arms[model] = (names, observed, chance)
        print(f"  {model:22s} n={len(names):3d}  precision={observed.mean():.3f}  "
              f"chance={chance.mean():.3f}  excess={observed.mean() - chance.mean():+.3f}")

    if len(models) != 2:
        return

    # The two arms cover the same families, so one index draw serves both. That
    # is also the design of the real comparison: the ablated arm is the same
    # panel, not an independent sample.
    first, second = (arms[m] for m in models)
    shared = [n for n in first[0] if n in set(second[0])]
    take = lambda arm: [arm[0].index(n) for n in shared]
    i_first, i_second = take(first), take(second)
    print(f"\n  families common to both arms: {len(shared)}")

    draws = np.empty(N_BOOTSTRAP)
    for b in range(N_BOOTSTRAP):
        pick = rng.integers(0, len(shared), size=len(shared))
        a = resample_excess(first[1][i_first], first[2][i_first], pick)
        c = resample_excess(second[1][i_second], second[2][i_second], pick)
        draws[b] = a - c

    point = (resample_excess(first[1][i_first], first[2][i_first],
                             np.arange(len(shared)))
             - resample_excess(second[1][i_second], second[2][i_second],
                               np.arange(len(shared))))
    low, high = np.percentile(draws, [2.5, 97.5])
    print(f"  difference of excesses: {point:+.3f}  "
          f"95% CI [{low:+.3f}, {high:+.3f}]  half-width={(high - low) / 2:.3f}")
    print(f"\n  This is the statistic the registration tests, bootstrapped whole.\n"
          f"  A difference smaller than {(high - low) / 2:.2f} is below what this\n"
          f"  design resolves at 95 percent, whatever the point estimate says.")


if __name__ == "__main__":
    main()
