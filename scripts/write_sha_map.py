"""Write docs/sha_map.md from the filter-repo commit map.

Run:  uv run --no-project --python 3.12 python scripts/write_sha_map.py

This repository was built from `rna-structure-awareness` with git-filter-repo,
which rewrites every commit id. Preregistrations and result notes written before
the rewrite quote the old ids, and those ids resolve in no repository that now
exists. The frozen documents are left byte-identical -- editing a preregistration
to fix its own self-reference is worse than a stale id -- so the mapping goes
here instead.

.git/filter-repo/commit-map is written by the rewrite and is not part of the
committed tree, so this script regenerates nothing: run it once, keep the output.
"""

import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
COMMIT_MAP = REPO / ".git" / "filter-repo" / "commit-map"
OUT = REPO / "docs" / "sha_map.md"
SHA = re.compile(r"\b[0-9a-f]{7}\b")

HEADER = """# Commit id map

`{n}` commits, rewritten from `rna-structure-awareness` into this repository with
git-filter-repo. Every id in this repository's history is the new one; documents
frozen before the rewrite quote the old one. Both are given at seven characters,
the length used throughout the manuscripts.

Documents quoting old ids are left unedited: a preregistration is a record of what
was frozen, and correcting its self-reference after the fact would make it a record
of something else. Resolve the id here instead.

| old | new | subject |
|---|---|---|
"""

UNCOVERED = """
## Ids this map does not cover

{ids}

These are quoted in superseded drafts and preregistrations under `docs/` and are not
in the rewrite's commit map. Two of them (`0aff46c`, `36da8a1`) resolve in the
pre-rewrite bundle `rna-structure-audit-BACKUP.git`, which is kept outside this
repository; the rest predate that bundle or are not commit ids at all -- `bbae163`
is a *Briefings in Bioinformatics* article number. No manuscript in `paper/` cites
any of them: `paper_v11.tex` pins `a207535`, `ae6712e` and `891d6af`, all of which
resolve here.
"""


def main() -> int:
    assert COMMIT_MAP.exists(), (
        f"{COMMIT_MAP} is missing -- it is written by git-filter-repo and is not "
        "committed. docs/sha_map.md must be kept as-is if the rewrite is not rerun.")

    pairs = []
    for line in COMMIT_MAP.read_text().splitlines():
        parts = line.split()
        if len(parts) != 2 or not all(len(p) == 40 for p in parts):
            continue
        old, new = parts
        if set(new) == {"0"}:  # commit dropped by the rewrite
            continue
        pairs.append((old, new))

    # Ordered oldest first, which is the order the map is written in.
    rows = []
    for old, new in pairs:
        subject = subprocess.run(
            ["git", "-C", str(REPO), "log", "--format=%s", "-1", new],
            capture_output=True, text=True, check=True).stdout.strip()
        rows.append(f"| `{old[:7]}` | `{new[:7]}` | {subject} |")

    OUT.write_text(HEADER.format(n=len(rows)) + "\n".join(rows) + "\n")
    print(f"wrote {OUT.relative_to(REPO)} ({len(rows)} commits)")

    # Every pre-rewrite id quoted anywhere in the tree must be resolvable here.
    quoted = set()
    for path in REPO.rglob("*"):
        if ".git/" in str(path) or not path.is_file():
            continue
        if path.suffix not in {".md", ".tex", ".py"} or path == OUT:
            continue
        quoted |= {m for m in SHA.findall(path.read_text(errors="ignore"))}
    known = {old[:7] for old, _ in pairs} | {new[:7] for _, new in pairs}
    unresolved = sorted(
        s for s in quoted
        if s not in known and not re.fullmatch(r"[0-9]{7}|[a-f]{7}", s))
    if unresolved:
        with OUT.open("a") as handle:
            handle.write(UNCOVERED.format(ids=", ".join(f"`{s}`" for s in unresolved)))
        print(f"appended {len(unresolved)} ids the map does not cover: "
              f"{' '.join(unresolved)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
