"""Resume a per-family loop that a preempted container left half-finished.

Every expensive stage here is a loop over families, and a stage that writes only
when the loop ends loses all of it when Modal reclaims the container. A stage
given a checkpoint writes after each family and reads back whatever an earlier
container already wrote, so a restarted run pays only for the families it has
not done.

The write is atomic. A container killed midway through `json.dump` leaves a
truncated file, and a truncated checkpoint is worse than none: it loads, it is
missing its tail, and the resumed run silently recomputes fewer families than it
reports. Writing a sibling and renaming means the checkpoint is either the old
complete one or the new complete one.

Shards carry the stamp of the run that wrote them and are discarded when a
later run's stamp differs. A volume outlives the code: shards written under one
commit, one panel and one set of library versions are still sitting there when
the next run starts, and a stage that resumes them produces a row assembled from
two different pipelines while reporting the second one's provenance. Discarding
costs a recomputation; resuming costs a number nobody can attribute.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable


class FamilyCheckpoint:
    """Per-family shards for one stage of one model.

    `on_write` is called after each write. On Modal it is the volume commit,
    which is what actually makes the shard survive the container; without it the
    file exists only in the container's view of the volume.

    `stamp` identifies the run. Any difference at all discards the shards --
    including a library version, which looks harmless and is exactly how two
    numerics end up in one column. Passing no stamp disables the check, which is
    what local runs and tests do.
    """

    def __init__(self, path: Path | str | None,
                 stamp: dict | None = None,
                 on_write: Callable[[], None] | None = None) -> None:
        self.path = Path(path) if path is not None else None
        self.stamp = stamp
        self.on_write = on_write
        self.done: dict[str, dict] = {}
        if self.path is None or not self.path.exists():
            return

        stored = json.loads(self.path.read_text())
        if stamp is None:
            self.done = stored.get("done", {})
            return
        if stored.get("stamp") != stamp:
            print(f"  discarding {len(stored.get('done', {}))} shard(s) in "
                  f"{self.path.name}: {self._mismatch(stored.get('stamp'), stamp)}")
            return
        self.done = stored.get("done", {})
        if self.done:
            print(f"  resuming {len(self.done)} family/families from {self.path.name}")

    @staticmethod
    def _mismatch(stored: dict | None, current: dict) -> str:
        """Which fields differ, so the log says why the work was thrown away."""
        if stored is None:
            return "written before shards carried a stamp"
        differing = sorted(
            key for key in set(stored) | set(current)
            if stored.get(key) != current.get(key)
        )
        return "differs in " + ", ".join(differing)

    def record(self, name: str, entry: dict) -> None:
        self.done[name] = entry
        if self.path is None:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_name(self.path.name + ".partial")
        tmp.write_text(json.dumps({"stamp": self.stamp, "done": self.done},
                                  indent=2, default=str))
        tmp.replace(self.path)
        if self.on_write is not None:
            self.on_write()


def no_checkpoint() -> FamilyCheckpoint:
    """A checkpoint that remembers nothing, for local runs and tests."""
    return FamilyCheckpoint(None)
