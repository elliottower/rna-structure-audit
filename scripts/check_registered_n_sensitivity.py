"""Does the registered Phase 6 decision rule change when N moves from 32 to 33?

Run:  uv run --no-project --with scipy --python 3.12 python scripts/check_registered_n_sensitivity.py

Correcting the four family annotations (commit 5fc8914) makes mir_21_precursor
eligible for Rung 3 under the criteria PREREGISTRATION_PHASE6_V2.md froze, so
the confirmatory set moves from 32 families to 33. The question this raises is
whether a registered decision threshold moves with it.

One confirmatory criterion carries N in its formula: H1 condition (a), the
minimum number of families whose PS must exceed their own null 95th percentile,
registered as ceil(4 * 0.05 * N). H1(b), H2 and H3 fix their thresholds
independently of N -- p < 0.0167, rank-biserial > 0.5, and a binomial test
against 1/3 -- and depend on the sample only through the data.

This script recomputes the H1(a) gate at both N and re-evaluates it against the
stored per-family exceedance flags. Nothing is re-run; every value is read from
results/.
"""

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
ALPHA = 0.05
MULTIPLIER = 4

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
    print("H1 condition (a): exceedances >= ceil(4 * 0.05 * N)\n")
    for n in (32, 33):
        print(f"  N = {n}   ceil({MULTIPLIER} * {ALPHA} * {n})"
              f" = ceil({MULTIPLIER * ALPHA * n:.1f}) = {gate(n)}")
    if gate(32) == gate(33):
        print(f"\n  The gate is {gate(32)} at both N. The registered decision "
              f"threshold does not move.")
    else:
        print(f"\n  The gate moves from {gate(32)} to {gate(33)}.")

    print("\nStored exceedance counts (pre-fix annotations, not re-run):\n")
    for model, rel in MODELS.items():
        over, passing = exceedances(rel)
        verdicts = ", ".join(f"N={n}: {'passes' if over >= gate(n) else 'fails'}"
                             for n in (32, 33))
        print(f"  {model:<10} {over} of {passing} exceed their null   {verdicts}")

    print("\n  These counts are computed against the superseded annotations and "
          "will change\n  when the four families are re-run. They are shown to "
          "locate the verdict\n  relative to the gate, which is 7 either way.")


if __name__ == "__main__":
    main()
