"""Assemble the panel the v24 tables read, from the two runs that hold its parts.

`results/repaired_panel_v3` is the run behind manuscript v23. `results/decoy_panel_v1`
is the registered rerun: Phase 6 for every arm, written by
`phase6_compensatory_mutation_v3.py` with a chance rate on matched pairs at every
layer and the matched-decoy comparisons, and every stage for the RNA-FM arms, whose
adapter was replaced.

The rule, per arm directory:

    an RNA-FM arm              every file from the rerun
    any other arm in the rerun its Phase 6 file from the rerun, the rest from v3
    an arm only in v3          every file from v3

Nothing is recomputed. `MANIFEST.json` names the source run and the sha256 of every
file, so a number in a table can be traced to the run that produced it. Checkpoint
shards are left out.

    uv run --no-project python scripts/build_panel_v4.py
"""

import hashlib
import json
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
V3 = REPO / "results" / "repaired_panel_v3"
RERUN = REPO / "results" / "decoy_panel_v1"
OUT = REPO / "results" / "panel_v4"


class PanelConflict(RuntimeError):
    """The two runs cannot be assembled into one panel."""


def is_rnafm(arm: str) -> bool:
    return arm == "rnafm" or arm.startswith("rnafm_")


def files(directory: Path) -> list[Path]:
    return sorted(p for p in directory.iterdir()
                  if p.is_file() and not p.name.startswith("_shards"))


def plan() -> dict[str, dict[str, Path]]:
    """For each arm, the source path of each file name."""
    arms = sorted({d.name for d in V3.iterdir() if d.is_dir()}
                  | {d.name for d in RERUN.iterdir() if d.is_dir()})
    chosen: dict[str, dict[str, Path]] = {}
    for arm in arms:
        old, new = V3 / arm, RERUN / arm
        if is_rnafm(arm):
            if not new.is_dir():
                raise PanelConflict(
                    f"{arm}: an RNA-FM arm of v3 with no rerun, so its numbers come "
                    "from the superseded adapter")
            chosen[arm] = {p.name: p for p in files(new)}
            continue
        picked = {p.name: p for p in files(old)} if old.is_dir() else {}
        if new.is_dir():
            phase6 = [p for p in files(new) if p.name.endswith("_phase6_ps.json")]
            if len(phase6) != 1:
                raise PanelConflict(f"{arm}: {len(phase6)} Phase 6 files in the rerun")
            picked[phase6[0].name] = phase6[0]
            if not old.is_dir():
                picked.update({p.name: p for p in files(new)})
        chosen[arm] = picked
    return chosen


def main() -> int:
    if OUT.exists():
        raise PanelConflict(f"{OUT} exists; a panel is assembled once and not rebuilt in place")
    chosen = plan()
    manifest = {"rule": __doc__.split("The rule, per arm directory:")[1].split("Nothing is")[0].strip(),
                "arms": {}}
    for arm, picked in chosen.items():
        (OUT / arm).mkdir(parents=True)
        entry = {}
        for name, source in sorted(picked.items()):
            shutil.copyfile(source, OUT / arm / name)
            entry[name] = {
                "from": str(source.relative_to(REPO)),
                "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            }
        manifest["arms"][arm] = entry
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    from_rerun = sum(1 for e in manifest["arms"].values() for f in e.values()
                     if f["from"].startswith("results/decoy_panel_v1"))
    total = sum(len(e) for e in manifest["arms"].values())
    print(f"wrote {OUT.relative_to(REPO)}: {len(chosen)} arms, {total} files, "
          f"{from_rerun} from the rerun")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
