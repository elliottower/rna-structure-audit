"""Copy the repaired panel out of the Modal volume into results/repaired_panel/.

One file at a time. `modal volume get` given a directory concatenates what it
finds into a single local file, which looks like a successful fetch and leaves
a directory's worth of runs in one unparseable blob, so every path here is a
file path.

The `_shards_*.json` checkpoints are left in the volume. They exist so a
reclaimed container resumes rather than restarts; every number in them is
reproduced in the stage's final file, and copying them doubles the size of the
deposit with a resume mechanism nobody reads.

A run whose files are not all there is reported and not fetched -- a directory
holding four of five stages is what a stopped app looks like, and half a run
copied into the results tree is indistinguishable from a finished one once the
volume is gone.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DESTINATION = REPO / "results" / "repaired_panel"
VOLUME = "rna-repaired-panel-results"

# What a finished directory holds, by the kind of run it is. `<key>` stands for
# the directory name, which is also the prefix every stage writes under.
EXPECTED = {
    # A transversion directory is named `<model>_transversion`, so `<key>`
    # already carries the suffix the runner puts in the file name.
    "transversion": ["<key>.json", "<key>_positions.json", "stamp.json"],
    "run": ["<key>_phases_1_to_5.json", "<key>_rung1_positions.json",
            "<key>_phase6_ps.json", "stamp.json"],
}


def volume_ls(path: str = "") -> list[str]:
    """Names directly under `path` in the volume, without their parent prefix."""
    listed = subprocess.run(["modal", "volume", "ls", VOLUME, path],
                            capture_output=True, text=True, check=True)
    return [line.split("/")[-1] for line in listed.stdout.split("\n") if line.strip()]


def expected_files(directory: str) -> list[str]:
    kind = "transversion" if directory.endswith("_transversion") else "run"
    return [name.replace("<key>", directory) for name in EXPECTED[kind]]


def fetch(directory: str, name: str, force: bool) -> bool:
    """One file. Returns whether it was fetched rather than already present."""
    local = DESTINATION / directory / name
    if local.exists() and not force:
        return False
    local.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["modal", "volume", "get", "--force", VOLUME,
                    f"{directory}/{name}", str(local)],
                   check=True, capture_output=True, text=True)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="re-fetch files already on disk")
    args = parser.parse_args()

    fetched = skipped = 0
    incomplete: list[tuple[str, list[str]]] = []
    for directory in sorted(volume_ls()):
        present = set(volume_ls(directory))
        wanted = expected_files(directory)
        missing = [name for name in wanted if name not in present]
        if missing:
            incomplete.append((directory, missing))
            continue
        for name in wanted:
            if fetch(directory, name, args.force):
                fetched += 1
                print(f"  {directory}/{name}")
            else:
                skipped += 1

    print(f"\n{fetched} fetched, {skipped} already on disk, "
          f"into {DESTINATION.relative_to(REPO)}")
    for directory, missing in incomplete:
        print(f"INCOMPLETE {directory}: missing {', '.join(missing)}")
    return 1 if incomplete else 0


if __name__ == "__main__":
    raise SystemExit(main())
