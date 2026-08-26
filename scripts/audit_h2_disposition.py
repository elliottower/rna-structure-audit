"""The registered H2 statistic, computed from stored results.

H2 compares per-model mean PS between the five RNA-pretrained and five
DNA-pretrained models by rank-biserial correlation, with a Mann-Whitney U
p-value, and is decided on rank-biserial > 0.5
(`PREREGISTRATION_PHASE6_V2.md:99`). The unit is the model, not the family, so
ten numbers decide it.

That makes the question for H2 not whether a given family's PS is resolvable but
whether the ten model means rank stably. This prints the ranking, the statistic,
and the same statistic over the character-tokenized models only, which the
registration calls for as E4 if H2 passes.

    uv run --no-project --python 3.12 python scripts/audit_h2_disposition.py \
        --results results/repaired_panel
"""

import argparse
import json
from pathlib import Path
from statistics import mean

from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

RNA_PRETRAINED = ["rnafm", "rinalmo", "utrlm", "ernierna", "splicebert"]
DNA_PRETRAINED = ["nt", "hyenadna", "caduceus", "evo", "dnabert2"]
# The registration's own note: three DNA models use multi-nucleotide tokens
# against none of the RNA models, which is what E4 exists to separate.
NON_CHARACTER = {"nt", "dnabert2", "evo"}


def model_mean_ps(path):
    per_rna = json.loads(path.read_text())["results"]["per_rna"]
    kept = [e["best_ps"] for name, e in per_rna.items()
            if name not in QUARANTINED and e.get("best_ps") is not None]
    return mean(kept) if kept else None


def rank_biserial(rna, dna):
    """Rank-biserial from Mann-Whitney U, oriented so RNA > DNA is positive."""
    result = mannwhitneyu(rna, dna, alternative="two-sided")
    return 2.0 * result.statistic / (len(rna) * len(dna)) - 1.0, result.pvalue


def report(label, means):
    rna = [means[m] for m in RNA_PRETRAINED if means.get(m) is not None]
    dna = [means[m] for m in DNA_PRETRAINED if means.get(m) is not None]
    if len(rna) < 2 or len(dna) < 2:
        print(f"  {label}: not enough models present ({len(rna)} vs {len(dna)})")
        return
    r, p = rank_biserial(rna, dna)
    verdict = "PASS" if r > 0.5 else "FAIL"
    print(f"  {label:34s} {len(rna)} vs {len(dna)}   "
          f"rank-biserial {r:+.3f}   p = {p:.4f}   {verdict}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=str(REPO / "results/repaired_panel"))
    args = parser.parse_args()
    root = Path(args.results)

    means = {}
    for model in RNA_PRETRAINED + DNA_PRETRAINED:
        path = root / model / f"{model}_phase6_ps.json"
        means[model] = model_mean_ps(path) if path.exists() else None

    print("Per-model mean PS over non-quarantined scored families, ranked:\n")
    ranked = sorted((v, m) for m, v in means.items() if v is not None)
    for rank, (value, model) in enumerate(ranked, start=1):
        group = "RNA" if model in RNA_PRETRAINED else "DNA"
        token = "" if model not in NON_CHARACTER else "  (multi-nucleotide tokens)"
        print(f"  {rank:2d}. {model:12s} {group}  {value:+.6f}{token}")
    missing = [m for m, v in means.items() if v is None]
    if missing:
        print(f"\n  not present in this directory: {', '.join(missing)}")

    print("\nDecision criterion is rank-biserial > 0.5.\n")
    report("H2, all ten models", means)
    character_only = {m: v for m, v in means.items() if m not in NON_CHARACTER}
    report("E4, character tokenizers only", character_only)
    print()


if __name__ == "__main__":
    main()
