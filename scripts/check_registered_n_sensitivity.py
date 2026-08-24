"""Does the registered Phase 6 decision rule survive the annotation repairs?

Run:  uv run --no-project --with scipy --python 3.12 python \
        scripts/check_registered_n_sensitivity.py

Two data corrections move the Rung 3 confirmatory set. Commit 5fc8914 made
mir_21_precursor eligible, taking N from 32 to 33. The repairs scoped in
`repair_annotations.py` make HDV_ribozyme, SAM_riboswitch and TPP_riboswitch
eligible, taking N to 36. The question both raise is whether a registered
decision threshold moves with N, and if so whether the verdict moves with it.

One confirmatory criterion carries N in its formula: H1 condition (a), the
minimum number of families whose PS must exceed their own null 95th percentile,
registered as ceil(4 * 0.05 * N). H1(b), H2 and H3 fix their thresholds
independently of N -- p < 0.0167, rank-biserial > 0.5, and a binomial test
against 1/3 -- and depend on the sample only through the data.

The gate does move, from 7 to 8. H1(a) nonetheless cannot flip, and the reason is
structural rather than empirical: every null in this pipeline is computed within a
family, so a family's exceedance flag does not depend on which other families are
present. No repair removes a scored family from the set -- none of the eleven
unattributed records carries a scored Rung 3 record -- so the exceedance count is
non-decreasing under the repairs while the gate rises by one. A model already
clearing the gate by a margin larger than one clears it afterwards whatever the
entering families do.

H1(b) has no such argument. The Wilcoxon statistic is computed across families,
so the entering families can move it in either direction, and its outcome is not
predictable from the stored values.

Nothing is re-run; every value is read from results/.
"""

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "docs" / "annotation_repair_manifest.json"
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
ALPHA = 0.05
MULTIPLIER = 4

N_REGISTERED = 32
N_AFTER_MIRNA_FIX = 33

MODELS = {
    "RiNALMo": "results/rinalmo_phase6_ps.json",
    "ERNIE-RNA": ("results/audit/phase6_ernierna_20260714_060520/"
                  "ernierna_phase6_ps.json"),
}


def gate(n: int) -> int:
    """H1 condition (a), verbatim from the preregistration."""
    return math.ceil(MULTIPLIER * ALPHA * n)


def exceedances(rel: str) -> tuple[int, int]:
    """(families exceeding their own null 95th percentile, families scored)."""
    body = json.loads((REPO / rel).read_text())["results"]["per_rna"]
    scored = {name: record for name, record in body.items()
              if isinstance(record, dict) and not record.get("skipped")
              and isinstance(record.get("best_ps"), (int, float))
              and name not in QUARANTINE}
    passing = [r for r in scored.values()
               if (r.get("positive_control") or {}).get("pass")]
    over = [r for r in passing
            if isinstance(r.get("null_95th_primary"), (int, float))
            and r["best_ps"] > r["null_95th_primary"]]
    return len(over), len(passing)


def main() -> None:
    n_after = json.loads(MANIFEST.read_text())["n_after"] \
        if MANIFEST.exists() else N_AFTER_MIRNA_FIX
    stages = [(N_REGISTERED, "registered"),
              (N_AFTER_MIRNA_FIX, "after the mir_21 correction"),
              (n_after, "after the annotation repairs")]

    print("H1 condition (a): exceedances >= ceil(4 * 0.05 * N)\n")
    for n, label in stages:
        print(f"  N = {n}   ceil({MULTIPLIER} * {ALPHA} * {n})"
              f" = ceil({MULTIPLIER * ALPHA * n:.1f}) = {gate(n)}    {label}")

    lo, hi = gate(N_REGISTERED), gate(n_after)
    if lo == hi:
        print(f"\n  The gate is {lo} throughout. The registered decision "
              f"threshold does not move.")
    else:
        print(f"\n  The gate moves from {lo} to {hi}. Every registered verdict "
              f"is re-evaluated\n  against {hi}.")

    print("\nStored exceedance counts (pre-repair annotations, not re-run):\n")
    for model, rel in MODELS.items():
        over, passing = exceedances(rel)
        verdicts = ", ".join(
            f"N={n}: {'passes' if over >= gate(n) else 'fails'}"
            for n, _ in stages)
        print(f"  {model:<10} {over} of {passing} exceed their null   {verdicts}")
        margin = over - hi
        print(f"{'':<12}clears the raised gate of {hi} by {margin}; "
              f"{'cannot' if margin > 0 else 'could'} flip when families enter")

    print("\n  Exceedance flags are per-family: the null is a within-family "
          "derangement, so a\n  family's flag does not depend on which other "
          "families are present. No repair\n  removes a scored family, so these "
          "counts can only rise. H1(a) is therefore\n  settled at the raised "
          "gate before the re-run. H1(b), a Wilcoxon across families,\n  is not, "
          "and its outcome is left open here.")


if __name__ == "__main__":
    main()
