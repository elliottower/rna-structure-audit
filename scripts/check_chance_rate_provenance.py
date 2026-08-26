"""Were all 17 derangement chance rates computed by the same code on the same set?

The reporting is about to rest on `h3_chance_fraction` as the only baseline, and
"it is already on disk" is the sentence that made the July ablation comparison
look valid. A chance rate computed under the float32 metric is not
interchangeable with one computed under float64: that defect moved ERNIE-RNA's
untrained precision from 0.161 to 0.476.

This checks three things per cell: the commit its stamp names, that the commit is
one of the two the final pass ran under, and that its analysis set is the same 35
families rather than 35 by coincidence.

    uv run --no-project --python 3.12 python scripts/check_chance_rate_provenance.py
"""

import argparse
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# The commits the final resolution pass ran under. Both contain the float64
# metric, `h3_chance_fraction`, and the token-span repair.
ACCEPTED = {"1350db19dddc", "3de3ba06965e"}
EXPECTED_DIGEST = "10bc3b36b3fe0d08"


def analysis_set(per_rna):
    return sorted(name for name, entry in per_rna.items()
                  if name not in QUARANTINED
                  and entry.get("h3_chance_fraction") is not None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=str(REPO / "results/repaired_panel_v3"))
    args = parser.parse_args()
    root = Path(args.results)

    print(f"  {'cell':22s} {'commit':>14s} {'n':>4s} {'set digest':>12s} {'':>8s}")
    problems = []
    for directory in sorted(root.iterdir()):
        path = directory / f"{directory.name}_phase6_ps.json"
        stamp_path = directory / "stamp.json"
        if not path.exists():
            continue
        stamp = json.loads(stamp_path.read_text())
        commit = stamp["commit"][:12]
        names = analysis_set(json.loads(path.read_text())["results"]["per_rna"])
        digest = hashlib.sha256(",".join(names).encode()).hexdigest()[:16]

        notes = []
        if commit not in ACCEPTED:
            notes.append("STALE COMMIT")
        if digest != EXPECTED_DIGEST:
            notes.append("DIFFERENT SET")
        if notes:
            problems.append((directory.name, notes))
        print(f"  {directory.name:22s} {commit:>14s} {len(names):4d} "
              f"{digest[:12]:>12s} {' '.join(notes):>8s}")

    print()
    if problems:
        for name, notes in problems:
            print(f"  PROBLEM {name}: {', '.join(notes)}")
        raise SystemExit(1)
    print(f"  All cells at an accepted commit, all on the same {len(names)} "
          f"families.\n  The chance rates are interchangeable across cells.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
