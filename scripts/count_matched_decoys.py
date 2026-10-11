"""Coverage of the matched-decoy comparison, from sequences and structures only.

Supersedes `count_complementary_decoys.py`, whose decoys could be shallow pairs
and whose pair totals ran over families below the inclusion threshold. No model
output is read. For each scope and decoy set of `decoy_comparison.matched_decoys`
this reports the families that enter the test (three or more deep-interior
pairs with a decoy), the pairs and comparisons inside those families, and how
much farther from the mutated position the decoys are than the partner.

    uv run --no-project python scripts/count_matched_decoys.py
"""

import json
from pathlib import Path

from decoy_comparison import DECOY_SETS, SCOPES, matched_decoys
from stem_pairs import get_eligible_pairs, is_deep_interior, parse_stems

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data" / "rfam_families"
OUT = REPO / "results" / "audit" / "matched_decoy_counts.json"

QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
MIN_WC_PAIRS = 15
MIN_ELIGIBLE = 5
MIN_PAIRS_PER_FAMILY = 3


def confirmatory_families() -> dict:
    """name -> (sequence, eligible pairs) for every family Rung 3 qualifies."""
    families = {}
    for path in sorted(FAMILIES.glob("*.json")):
        family = json.loads(path.read_text())
        if "excluded" in family or family["name"] in QUARANTINED:
            continue
        sequence, structure = family["sequence"], family["dot_bracket"]
        if len(sequence) != len(structure):
            continue
        stems = parse_stems(structure, sequence)
        if sum(len(stem) for stem in stems) < MIN_WC_PAIRS:
            continue
        eligible = get_eligible_pairs(sequence, stems)
        if len(eligible) < MIN_ELIGIBLE:
            continue
        families[family["name"]] = (sequence, eligible)
    return families


def count(sequence, eligible, scope, decoy_set) -> dict:
    pairs = comparisons = farther = 0
    excess = []
    for index, pair in enumerate(eligible):
        if not is_deep_interior(pair):
            continue
        found = matched_decoys(sequence, eligible, index, scope, decoy_set)
        if not found:
            continue
        pairs += 1
        comparisons += len(found)
        reach = abs(pair["i"] - pair["j"])
        for decoy in found:
            distance = abs(pair["i"] - decoy)
            excess.append(distance - reach)
            farther += distance > reach
    return {"pairs": pairs, "comparisons": comparisons, "decoys_farther": farther,
            "sum_excess_distance": sum(excess)}


def main() -> int:
    families = confirmatory_families()
    per_family = {name: {"deep_interior_pairs": sum(is_deep_interior(p) for p in eligible)}
                  for name, (_sequence, eligible) in families.items()}
    summary = {"families_confirmatory": len(families),
               "deep_interior_pairs": sum(row["deep_interior_pairs"]
                                          for row in per_family.values())}
    for scope in SCOPES:
        for decoy_set in DECOY_SETS:
            key = f"{scope}_{decoy_set}"
            analyzed = []
            for name, (sequence, eligible) in families.items():
                row = count(sequence, eligible, scope, decoy_set)
                per_family[name][key] = row
                if row["pairs"] >= MIN_PAIRS_PER_FAMILY:
                    analyzed.append(name)
            rows = [per_family[name][key] for name in analyzed]
            comparisons = sum(row["comparisons"] for row in rows)
            summary[key] = {
                "families_analyzed": len(analyzed),
                "pairs": sum(row["pairs"] for row in rows),
                "comparisons": comparisons,
                "fraction_decoys_farther": (
                    sum(row["decoys_farther"] for row in rows) / comparisons
                    if comparisons else None),
                "mean_excess_distance": (
                    sum(row["sum_excess_distance"] for row in rows) / comparisons
                    if comparisons else None),
                "families": analyzed,
            }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "inputs": str(FAMILIES.relative_to(REPO)),
        "quarantined": sorted(QUARANTINED),
        "thresholds": {"min_wc_pairs": MIN_WC_PAIRS, "min_eligible": MIN_ELIGIBLE,
                       "min_pairs_per_family": MIN_PAIRS_PER_FAMILY},
        "summary": summary,
        "families": per_family,
    }, indent=2))
    for key, value in summary.items():
        if isinstance(value, dict):
            value = {k: v for k, v in value.items() if k != "families"}
        print(key, value)
    print(f"wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
