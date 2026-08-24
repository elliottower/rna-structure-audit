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

The gate does move, from 7 to 8, and H1(a) holds at the raised gate. Each null is
built inside one family, so the distribution a family is compared against does not
depend on which other families are present -- but the draw did, until
`family_seed.py`: both rungs seeded from a family's position in the loaded panel,
so withdrawing five records would have shifted the stream every later family
received. Reseeding from the family name redraws every null once, so the no-flip
argument rests on how much room each family has rather than on its flag being
fixed. `check_h1a_margin.py` measures that room. No repair removes a scored family
-- none of the eleven unattributed records carries a scored Rung 3 record -- so
entering families can only add to the count.

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


def scored_families(rel: str) -> dict[str, dict]:
    """The families the registered Rung 3 analysis scores, from a results file.

    A family counts when it carries a numeric PS, is not skipped, is not
    quarantined for pilot foreknowledge, and passes its positive control. This is
    the one definition of the confirmatory set; every script that needs it imports
    from here rather than restating the filter, because a second copy drifts.
    """
    body = json.loads((REPO / rel).read_text())["results"]["per_rna"]
    scored = {name: record for name, record in body.items()
              if isinstance(record, dict) and not record.get("skipped")
              and isinstance(record.get("best_ps"), (int, float))
              and name not in QUARANTINE}
    return {name: record for name, record in scored.items()
            if (record.get("positive_control") or {}).get("pass")}


def exceedances(rel: str) -> tuple[int, int]:
    """(families exceeding their own null 95th percentile, families scored)."""
    passing = list(scored_families(rel).values())
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
        print(f"{'':<12}clears the raised gate of {hi} by {margin} against the "
              f"entering families,\n{'':<12}which can only add; the reseed is the "
              f"perturbation this count cannot bound")

    print("\n  No repair removes a scored family, so these counts can only "
          "rise. The nulls are\n  redrawn once, because they are now seeded from "
          "family names rather than from\n  panel positions, so run "
          "check_h1a_margin.py for how much room each family has\n  above its own "
          "threshold. H1(b), a Wilcoxon across families, is left open here\n  and "
          "is reported whichever way it lands.")


if __name__ == "__main__":
    main()
