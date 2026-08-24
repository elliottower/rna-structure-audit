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
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable


class FamilyCheckpoint:
    """Per-family shards for one stage of one model.

    `commit` is called after each write. On Modal it is the volume commit, which
    is what actually makes the shard survive the container; without it the file
    exists only in the container's view of the volume.
    """

    def __init__(self, path: Path | str | None,
                 commit: Callable[[], None] | None = None) -> None:
        self.path = Path(path) if path is not None else None
        self.commit = commit
        self.done: dict[str, dict] = {}
        if self.path is not None and self.path.exists():
            self.done = json.loads(self.path.read_text())

    def record(self, name: str, entry: dict) -> None:
        self.done[name] = entry
        if self.path is None:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_name(self.path.name + ".partial")
        tmp.write_text(json.dumps(self.done, indent=2, default=str))
        tmp.replace(self.path)
        if self.commit is not None:
            self.commit()


def no_checkpoint() -> FamilyCheckpoint:
    """A checkpoint that remembers nothing, for local runs and tests."""
    return FamilyCheckpoint(None)
