"""How far is each family from its own null, and can the H1(a) count move?

Run:  uv run --no-project --with numpy --python 3.12 \
          python scripts/check_h1a_margin.py

The deviation entry registered 2026-08-24 claims H1 condition (a) cannot fail
under the annotation repair. The claim was stated structurally: the null is a
within-family derangement, so a family's exceedance flag does not depend on which
other families are present. That is true of the statistic and false of the
implementation. Both rungs seed from panel position --
`phases_1_to_5.py:265` uses `default_rng(42 + rna_idx * 1000)`, and
`phase6_compensatory_mutation.py:595` seeds one global generator that every family
then draws from in order -- so withdrawing five families shifts the stream every
later family receives. Each family's null is redrawn from the same distribution,
not reused.

So exceedance flags can move, and the question is how many. A family whose PS sits
far above its null 95th percentile will not cross it under a redraw; one sitting on
the threshold may. This script measures the margin for every scored family and
reports how many are close enough to be at risk, which is what the no-flip claim
actually rests on.

No model is run; every number comes from the stored results.
"""

import json
from pathlib import Path

from check_registered_n_sensitivity import scored_families

REPO = Path(__file__).resolve().parents[1]
MODELS = {
    "RiNALMo": "results/rinalmo_phase6_ps.json",
    "ERNIE-RNA": "results/audit/phase6_ernierna_20260714_060520/ernierna_phase6_ps.json",
}


def margins(rel: str) -> list[tuple[str, float, float, float]]:
    """(family, PS, primary null 95th, distance to the conservative null).

    Nulls here are frequently negative -- PS is a signed quantity -- so a margin
    expressed as a fraction of the threshold is meaningless and, worse, silently
    changes sign. The comparison used instead is absolute: how far PS sits above
    the primary null, against how far the two registered null constructions sit
    from each other. Swapping the null construction is a far larger perturbation
    than reseeding one of them, so a family clearing the first distance by more
    than the second will not cross its threshold on a redraw.
    """
    out = []
    for name, record in scored_families(rel).items():
        ps = record["best_ps"]
        primary = record.get("null_95th_primary")
        conservative = record.get("null_95th_conservative")
        if not isinstance(primary, (int, float)):
            continue
        spread = (abs(conservative - primary)
                  if isinstance(conservative, (int, float)) else 0.0)
        out.append((name, ps, primary, spread))
    return out


def main() -> None:
    gate_before, gate_after = 7, 8
    print(f"H1(a) gate: {gate_before} at the registered N = 32, "
          f"{gate_after} at N = 36\n")

    for model, rel in MODELS.items():
        rows = margins(rel)
        exceed = [r for r in rows if r[1] > r[2]]
        risky = [r for r in exceed if (r[1] - r[2]) <= r[3]]
        worst = sorted(exceed, key=lambda r: r[1] - r[2])

        print(f"{model}: {len(exceed)} of {len(rows)} scored families exceed "
              f"their null")
        print(f"  {'family':<28}{'PS':>10}{'null':>10}{'gap':>10}"
              f"{'null spread':>13}")
        for name, ps, null, spread in worst[:4]:
            print(f"      {name:<24}{ps:>10.4f}{null:>10.4f}{ps - null:>10.4f}"
                  f"{spread:>13.4f}")
        print(f"  families whose gap is smaller than the distance between the two "
              f"null\n  constructions, and so could plausibly cross on a redraw: "
              f"{len(risky)}")
        for name, ps, null, spread in sorted(risky):
            print(f"      {name:<24}{ps:>10.4f}{null:>10.4f}{ps - null:>10.4f}"
                  f"{spread:>13.4f}")

        floor = len(exceed) - len(risky)
        verdict = "HOLDS" if floor >= gate_after else "AT RISK"
        print(f"  worst case, every at-risk family crosses the wrong way: {floor}")
        print(f"  against the raised gate of {gate_after}: {verdict} "
              f"(margin {floor - gate_after:+d})\n")

    print("Entering families can only add to the count, so the floor above is the\n"
          "binding case. H1(b) is a Wilcoxon across families and is not bounded by\n"
          "this argument; its outcome is left open.")


if __name__ == "__main__":
    main()
