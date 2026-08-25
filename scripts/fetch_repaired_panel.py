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
import gzip
import shutil
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


ABLATION_SUFFIX = "_noattnbias"


def expected_files(directory: str) -> list[str]:
    """The files a finished directory holds, by the kind of run it is.

    An ablated run writes into its own directory but names its files after the
    model rather than the directory, because `result_dir` carries the suffix and
    the file writer does not. The transversion path does carry it, so the two
    layouts differ. This reads what is on the volume rather than what would have
    been tidier; renaming would orphan files already written and cost a re-run.
    """
    kind = "transversion" if directory.endswith("_transversion") else "run"
    key = (directory[:-len(ABLATION_SUFFIX)]
           if directory.endswith(ABLATION_SUFFIX) else directory)
    return [name.replace("<key>", key) for name in EXPECTED[kind]]


def fetch(directory: str, name: str, force: bool, destination: Path) -> bool:
    """One file. Returns whether it was fetched rather than already present."""
    local = destination / directory / name
    stored = compressed(local)
    if stored.exists() and not force:
        return False
    local.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["modal", "volume", "get", "--force", VOLUME,
                    f"{directory}/{name}", str(local)],
                   check=True, capture_output=True, text=True)
    if stored != local:
        with local.open("rb") as raw, gzip.open(stored, "wb") as out:
            shutil.copyfileobj(raw, out)
        local.unlink()
    return True


def compressed(local: Path) -> Path:
    """Where a fetched file is kept.

    The per-position sidecars are the bulk of a run -- 47 MB of the 48 -- and
    they are arrays of rounded floats, which gzip to about a sixth. They are
    stored compressed so a repository that may go public does not carry the
    uncompressed form in its history forever. Everything else is small enough
    that compressing it would only make it harder to read.
    """
    return (local.with_suffix(local.suffix + ".gz")
            if local.name.endswith("_positions.json") else local)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="re-fetch files already on disk")
    # A re-run is a new set of numbers, and writing it over the directory the
    # deposited results live in destroys the only copy of what was reported.
    parser.add_argument("--dest", type=Path, default=DESTINATION,
                        help="directory to write into (default: the deposited one)")
    args = parser.parse_args()
    destination = args.dest if args.dest.is_absolute() else REPO / args.dest

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
            if fetch(directory, name, args.force, destination):
                fetched += 1
                print(f"  {directory}/{name}")
            else:
                skipped += 1

    print(f"\n{fetched} fetched, {skipped} already on disk, "
          f"into {destination.relative_to(REPO)}")
    for directory, missing in incomplete:
        print(f"INCOMPLETE {directory}: missing {', '.join(missing)}")
    return 1 if incomplete else 0


if __name__ == "__main__":
    raise SystemExit(main())
