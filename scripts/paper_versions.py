"""Resolve the newest manuscript version.

Versions are never edited in place: paper_vN.tex is generated from paper_vN-1.tex
by a script in paper/patches/, and the older file stays as the record of what it
said. Verification scripts therefore have to find the current version rather than
name one, or they silently keep checking a superseded file.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO / "paper"
VERSIONED = re.compile(r"^paper_v(\d+)")


def newest_paper() -> Path:
    """The highest-numbered paper_vN*.tex in paper/."""
    versions = [(int(match.group(1)), path)
                for path in PAPER_DIR.glob("paper_v*.tex")
                if (match := VERSIONED.match(path.name))]
    if not versions:
        raise FileNotFoundError(f"no paper_vN.tex in {PAPER_DIR}")
    return max(versions)[1]
