"""How many pairs and families the complementary-decoy comparison can use.

Computed from sequences and annotated structures only; no model output is read.
For each family that qualifies for Rung 3, counts the deep-interior pairs that
have at least one complementary decoy in their own stem and in the whole family
(`stem_pairs.complementary_decoys`), and the number of (pair, decoy) comparisons.

    uv run --no-project python scripts/count_complementary_decoys.py
"""

import json
from pathlib import Path

from stem_pairs import (
    complementary_decoys,
    get_eligible_pairs,
    is_deep_interior,
    parse_stems,
)

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data" / "rfam_families"
OUT = REPO / "results" / "audit" / "complementary_decoy_counts.json"

QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
MIN_WC_PAIRS = 15
MIN_ELIGIBLE = 5
MIN_PAIRS_PER_FAMILY = 3


def count(family: dict) -> dict | None:
    """Counts for one family, or None where Rung 3 does not qualify it."""
    sequence, structure = family["sequence"], family["dot_bracket"]
    if len(sequence) != len(structure):
        return None
    stems = parse_stems(structure, sequence)
    if sum(len(stem) for stem in stems) < MIN_WC_PAIRS:
        return None
    eligible = get_eligible_pairs(sequence, stems)
    if len(eligible) < MIN_ELIGIBLE:
        return None
    deep = [k for k, pair in enumerate(eligible) if is_deep_interior(pair)]
    row = {"eligible_pairs": len(eligible), "deep_interior_pairs": len(deep)}
    for scope in ("stem", "family"):
        sizes = [len(complementary_decoys(sequence, eligible, k, scope)) for k in deep]
        row[scope] = {"pairs_with_decoy": sum(size > 0 for size in sizes),
                      "comparisons": sum(sizes)}
    return row


def main() -> int:
    families = {}
    for path in sorted(FAMILIES.glob("*.json")):
        family = json.loads(path.read_text())
        if "excluded" in family or family["name"] in QUARANTINED:
            continue
        row = count(family)
        if row is not None:
            families[family["name"]] = row

    summary = {"families_qualifying": len(families)}
    for scope in ("stem", "family"):
        usable = [row for row in families.values()
                  if row[scope]["pairs_with_decoy"] >= MIN_PAIRS_PER_FAMILY]
        summary[scope] = {
            "families_with_min_pairs": len(usable),
            "pairs_with_decoy": sum(row[scope]["pairs_with_decoy"]
                                    for row in families.values()),
            "comparisons": sum(row[scope]["comparisons"] for row in families.values()),
        }
    summary["deep_interior_pairs"] = sum(row["deep_interior_pairs"]
                                         for row in families.values())

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "inputs": str(FAMILIES.relative_to(REPO)),
        "quarantined": sorted(QUARANTINED),
        "thresholds": {"min_wc_pairs": MIN_WC_PAIRS, "min_eligible": MIN_ELIGIBLE,
                       "min_pairs_per_family": MIN_PAIRS_PER_FAMILY},
        "summary": summary,
        "families": families,
    }, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
