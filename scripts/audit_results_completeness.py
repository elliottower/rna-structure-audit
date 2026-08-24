"""Which model-by-family cells does results/ actually contain?

Run:  uv run --no-project --python 3.12 python scripts/audit_results_completeness.py

The manuscript reports ten models over 52 families across three rungs. A missing
cell is invisible in an aggregate: a mean over the families that happen to have a
stored record reads the same as a mean over all of them. This script walks every
JSON under results/, finds the blocks keyed by family name, and reports how many
of the 52 each file covers.

Family coverage is counted against the deposited panel, so a file naming a family
that no longer exists is reported as unknown rather than silently dropped.

No model is run.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PANEL = {json.loads(p.read_text())["name"]
         for p in (REPO / "data/rfam_families").glob("*.json")}


def family_blocks(node, found):
    """Every dict whose keys are mostly family names."""
    if isinstance(node, dict):
        keys = set(node)
        overlap = keys & PANEL
        if len(overlap) >= 3:
            found.append((overlap, keys - PANEL))
        for value in node.values():
            family_blocks(value, found)
    elif isinstance(node, list):
        for value in node:
            family_blocks(value, found)
    return found


def main() -> None:
    files = sorted((REPO / "results").rglob("*.json"))
    print(f"{len(files)} JSON files under results/, panel of {len(PANEL)} families\n")
    print(f"  {'file':<62}{'families':>9}{'unknown':>9}")

    total_unknown, no_block = set(), []
    for path in files:
        try:
            body = json.loads(path.read_text())
        except json.JSONDecodeError:
            print(f"  {str(path.relative_to(REPO)):<62}{'unreadable':>9}")
            continue
        blocks = family_blocks(body, [])
        if not blocks:
            no_block.append(path.relative_to(REPO))
            continue
        covered = set().union(*(b[0] for b in blocks))
        unknown = set().union(*(b[1] for b in blocks))
        unknown = {k for k in unknown if not k.startswith("_")}
        total_unknown |= unknown
        flag = "" if len(covered) == len(PANEL) else "  <"
        print(f"  {str(path.relative_to(REPO)):<62}{len(covered):>9}"
              f"{len(unknown):>9}{flag}")

    print(f"\n  files with no family-keyed block: {len(no_block)}")
    for path in no_block:
        print(f"      {path}")

    if total_unknown:
        print(f"\n  keys appearing in results/ that are not families in the "
              f"deposited panel:\n      {', '.join(sorted(total_unknown))}")


if __name__ == "__main__":
    main()
