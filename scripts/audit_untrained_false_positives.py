"""How often does the Rung 1 null reject for a model that cannot have learned anything?

Run:  uv run --no-project --with scipy --python 3.12 python \
          scripts/audit_untrained_false_positives.py

Seven models were run with randomly initialized weights. Whatever a random
network's stem-versus-loop ratio is, it is not knowledge of RNA structure, so the
rate at which those runs exceed their own composition-preserving null is an
estimate of the null's false-positive rate. The registered null is a 95th
percentile, so that rate should be 0.05.

The null permutes stem and loop labels within nucleotide strata. A permutation
draw scatters the labels; the real annotation does not, because stems and loops
are contiguous runs. Any embedding property that varies smoothly along the
sequence -- a positional encoding, a convolution's receptive field, an attention
window -- therefore separates the real labels from the permuted ones for reasons
that have nothing to do with pairing, and the null underestimates the spread of
the statistic it is thresholding.

This script measures the rate rather than arguing for it. It also reports which
families carry the exceedances, because a null that fails through run structure
fails on the same sequences for every model.

No model is run.
"""

import json
from collections import Counter
from pathlib import Path

from scipy.stats import binomtest

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results/repaired_panel"

UNTRAINED = ["ernierna", "rinalmo", "rnafm", "utrlm", "splicebert", "nt", "dnabert2"]
TRAINED = UNTRAINED + ["hyenadna", "caduceus", "evo"]

# Models whose tokens are not single nucleotides, so their per-position
# embeddings are produced by the expansion in `phases_1_to_5.py:61-76` rather
# than read off directly. Kept in step with `NON_CHARACTER_TOKENIZERS` in
# `scripts/phase6_compensatory_mutation.py`.
NON_CHARACTER = ["nt", "dnabert2"]

NOMINAL = 0.05


def load_rung1(key: str) -> dict:
    path = PANEL / key / f"{key}_phases_1_to_5.json"
    return json.loads(path.read_text())["mutation_trained"]["per_rna"]


def scored(per_rna: dict) -> dict:
    """Families that reached a ratio. A family with too few stem or loop positions
    at every layer is not scored and is not a trial."""
    return {name: entry for name, entry in per_rna.items()
            if entry.get("best_ratio") is not None and "nuc_null_95th" in entry}


def main() -> None:
    print("Rung 1 exceedances, randomly initialized weights\n")
    print(f"  {'run':<14}{'scored':>8}{'>null':>8}{'rate':>9}"
          f"{'>dinuc':>9}   families exceeding")

    pooled_hits = pooled_trials = 0
    by_family: Counter = Counter()
    for key in UNTRAINED:
        per_rna = scored(load_rung1(f"{key}_untrained"))
        hits = sorted(name for name, e in per_rna.items() if e["exceeds_nuc_null"])
        dinuc = sum(1 for name in hits
                    if per_rna[name].get("exceeds_dinuc_null"))
        pooled_hits += len(hits)
        pooled_trials += len(per_rna)
        by_family.update(hits)
        rate = len(hits) / len(per_rna) if per_rna else 0.0
        print(f"  {key:<14}{len(per_rna):>8}{len(hits):>8}{rate:>9.3f}"
              f"{dinuc:>9}   {', '.join(hits)}")

    rate = pooled_hits / pooled_trials
    test = binomtest(pooled_hits, pooled_trials, NOMINAL, alternative="greater")
    print(f"\n  pooled {pooled_hits}/{pooled_trials} = {rate:.3f} against a nominal "
          f"{NOMINAL:.2f}")
    print(f"  exact binomial, one-sided: p = {test.pvalue:.2e}")
    print(f"  95% CI on the rate: [{test.proportion_ci().low:.3f}, "
          f"{test.proportion_ci().high:.3f}]")
    print(f"  the null is anticonservative by a factor of {rate / NOMINAL:.1f}\n")

    print("  stratified by tokenizer. `nt` and `dnabert2` do not tokenize one\n"
          "  nucleotide per token, so their embeddings reach the per-position\n"
          "  metric through the expansion in `phases_1_to_5.py:61-76`; the other\n"
          "  five do not pass through it at all.\n")
    for label, keys in (("character-level", [k for k in UNTRAINED
                                             if k not in NON_CHARACTER]),
                        ("non-character", NON_CHARACTER)):
        hits = trials = 0
        for key in keys:
            per_rna = scored(load_rung1(f"{key}_untrained"))
            hits += sum(1 for e in per_rna.values() if e["exceeds_nuc_null"])
            trials += len(per_rna)
        stratum = binomtest(hits, trials, NOMINAL, alternative="greater")
        low, high = stratum.proportion_ci(confidence_level=0.95,
                                          method="exact")
        print(f"      {label:<18}{hits:>4}/{trials:<6}= {hits / trials:.3f}   "
              f"p = {stratum.pvalue:.3f}   one-sided 95% lower bound "
              f"{low:.3f}")
    print()

    print("  families carrying the exceedances, of 7 untrained runs each\n")
    for name, count in by_family.most_common():
        if count >= 2:
            print(f"      {name:<26}{count}")
    once = [name for name, count in by_family.items() if count == 1]
    print(f"\n      {len(once)} further families exceed in exactly one run")

    print("\n  trained runs, for scale -- not a false-positive rate, because a "
          "trained\n  model may exceed the null for the reason the null is "
          "testing for\n")
    for key in TRAINED:
        per_rna = scored(load_rung1(key))
        hits = sum(1 for e in per_rna.values() if e["exceeds_nuc_null"])
        print(f"      {key:<14}{hits:>3}/{len(per_rna)}")


if __name__ == "__main__":
    main()
