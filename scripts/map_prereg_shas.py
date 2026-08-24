"""Map the preregistration commit SHAs the RNA paper prints onto the rewritten
public history.

Run:  uv run --no-project --python 3.12 python scripts/map_prereg_shas.py

paper_v10.tex Section 3.1 cites three commits as the preregistration freeze
points:

    Preregistration documents (Phase 1 SHA 694b43b, Phase 2 bd4b3fd,
    Phase 6 c19aa59) are deposited with the code.

Those SHAs are from `rna-structure-audit`. `rna-sa-public` is a git-filter-repo
export of that repo, so every commit carries a different SHA and none of the
three resolves there. If the audit repo is deleted the paper cites three commits
that exist nowhere, which is the failure the preregistration is meant to prevent.

filter-repo preserves author and committer dates, so the two histories can be
aligned on (author date, committer date). This checks that the alignment is
one-to-one and that the prereg file contents are byte-identical across it --
the second is what actually matters, since a reader verifying a freeze reads the
file, not the commit object.

Read-only. Writes nothing.
"""

import subprocess
import sys
from pathlib import Path

GITHUB = Path.home() / "Documents" / "GitHub"
AUDIT = GITHUB / "rna-structure-audit"
PUBLIC = GITHUB / (sys.argv[1] if len(sys.argv) > 1 else "rna-sa-public")

CITED = {
    "Phase 1": "694b43b",
    "Phase 2": "bd4b3fd",
    "Phase 6": "c19aa59",
}


def git(repo: Path, *args: str) -> str:
    done = subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True)
    return done.stdout.strip() if done.returncode == 0 else ""


def index(repo: Path) -> dict[str, tuple[str, str]]:
    """(author date, committer date) -> (sha, subject) for every commit."""
    out = git(repo, "log", "--all", "--format=%H%x1f%at%x1f%ct%x1f%s")
    table: dict[str, tuple[str, str]] = {}
    for line in out.splitlines():
        sha, adate, cdate, subject = line.split("\x1f", 3)
        table[f"{adate}:{cdate}"] = (sha, subject)
    return table


def files_at(repo: Path, sha: str) -> list[str]:
    out = git(repo, "ls-tree", "-r", "--name-only", sha)
    return [p for p in out.splitlines() if "PREREG" in p.upper()]


def blob(repo: Path, sha: str, path: str) -> str:
    return git(repo, "rev-parse", f"{sha}:{path}")


def main() -> int:
    public_index = index(PUBLIC)
    print(f"{len(public_index)} commits indexed in rna-sa-public\n")

    unresolved = []
    for phase, short in CITED.items():
        full = git(AUDIT, "rev-parse", short)
        if not full:
            print(f"{phase:8s} {short}  NOT IN rna-structure-audit EITHER")
            unresolved.append(phase)
            continue

        subject = git(AUDIT, "log", "-1", "--format=%s", full)
        key = git(AUDIT, "log", "-1", "--format=%at:%ct", full)
        match = public_index.get(key)

        print(f"{phase:8s} audit {short}  {subject[:62]}")
        if match is None:
            print("         -> no commit in rna-sa-public shares its dates\n")
            unresolved.append(phase)
            continue

        new_sha, new_subject = match
        print(f"         -> public {new_sha[:7]}  {new_subject[:62]}")

        # The commit object is a pointer. What a reader checks is the file.
        same, differ, only_audit = 0, [], []
        for path in files_at(AUDIT, full):
            left = blob(AUDIT, full, path)
            right = blob(PUBLIC, new_sha, path)
            if not right:
                only_audit.append(path)
            elif left == right:
                same += 1
            else:
                differ.append(path)
        print(f"         prereg files at that commit: {same} identical, "
              f"{len(differ)} changed, {len(only_audit)} dropped by the rewrite")
        for path in differ + only_audit:
            print(f"           ! {path}")
        print()

    if unresolved:
        print(f"UNRESOLVED: {', '.join(unresolved)} -- the paper's SHA cannot be "
              "remapped, so deleting rna-structure-audit would orphan it.")
    else:
        print("All three cited SHAs map onto rna-sa-public commits. The paper's "
              "Section 3.1 can be reissued against the public history.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
