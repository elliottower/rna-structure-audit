"""Did each preregistration freeze land before the results it governs?

Run:  uv run --no-project --python 3.12 python scripts/check_freeze_order.py

A preregistration is worth nothing if the data it predicts was already in hand.
The result filenames in this repo carry their own run timestamps, and each
prereg has a commit date, so the ordering is checkable rather than a matter of
recollection.

Reports, per phase, the freeze date and the earliest and latest run it governs.
A run earlier than its freeze is printed as VIOLATION. This measures order only:
it cannot see whether a result was looked at before the freeze, which is a
separate question the preregistration text has to answer for itself.
"""

import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
STAMP = re.compile(r"(20\d{6})_?(\d{6})?")

# prereg file -> the result paths whose runs it governs
PHASES = {
    "Phase 1 (structure metrics)": ("PREREGISTRATION_STRUCTURE_METRICS.md",
                                    ["data/phases_1_to_5", "data/structure_metrics"]),
    "Phase 2 (expanded Rfam)": ("docs/PREREGISTRATION_EXPANDED_RFAM.md",
                                ["data/expanded_rfam"]),
    "Phase 6 V1 (compensatory)": ("PREREGISTRATION_PHASE6_COMPENSATORY_MUTATION.md",
                                  ["data/gpu_results/phase6_compensatory"]),
    # V2 replaces V1 and defines a different metric for the same question. Runs
    # dated before the V2 freeze were registered under V1, so they are not
    # violations; V2 handles them by quarantining the two families whose values
    # had been seen. Its confirmatory runs are the post-freeze results/audit set.
    "Phase 6 V2": ("PREREGISTRATION_PHASE6_V2.md", ["results/audit"]),
    "Phase 6 untrained": ("preregistration/PREREGISTRATION_PHASE6_UNTRAINED_RINALMO.md",
                          ["data/gpu_results/phase6_untrained"]),
}


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args],
                          capture_output=True, text=True).stdout.strip()


def freeze_date(path: str) -> str:
    """Date the file first entered history -- the freeze, not later edits."""
    return git("log", "--all", "--diff-filter=A", "--format=%ad", "--date=short",
               "--", path).split("\n")[-1]


def run_dates(prefixes: list[str]) -> list[tuple[str, str]]:
    found = []
    for prefix in prefixes:
        for path in git("ls-tree", "-r", "--name-only", "HEAD", "--", prefix).split("\n"):
            if not path.strip():
                continue
            match = STAMP.search(path)  # the stamp is often on the directory
            if match:
                day = match.group(1)
                found.append((f"{day[:4]}-{day[4:6]}-{day[6:]}", path))
    return sorted(found)


def main() -> int:
    violations = 0
    for phase, (prereg, prefixes) in PHASES.items():
        frozen = freeze_date(prereg)
        runs = run_dates(prefixes)
        print(f"\n{phase}")
        print(f"  prereg   {prereg}")
        if not frozen:
            print("  frozen   NOT IN THIS HISTORY")
            continue
        print(f"  frozen   {frozen}")
        if not runs:
            print("  runs     none found under " + ", ".join(prefixes))
            continue
        early, late = runs[0], runs[-1]
        print(f"  runs     {len(runs)} files, {early[0]} .. {late[0]}")
        before = [(d, p) for d, p in runs if d < frozen]
        if before:
            violations += len(before)
            print(f"  VIOLATION  {len(before)} run(s) predate the freeze:")
            for day, path in before[:8]:
                print(f"    {day}  {path}")
        else:
            print("  OK       every run is dated on or after the freeze")

    print(f"\n{'-' * 60}")
    print("Order is necessary and not sufficient. PREREGISTRATION_PHASE6_V2.md "
          "discloses\nthat RNA-FM values for tRNA_Phe_yeast and tRNA_Ala_human "
          "were seen before it was\nwritten, and quarantines both families. "
          "Whether that quarantine was applied is\nchecked by "
          "recompute_phase6_quarantined.py, not here.")
    print(f"{violations} run(s) dated before their preregistration."
          if violations else
          "No run is dated before its preregistration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
