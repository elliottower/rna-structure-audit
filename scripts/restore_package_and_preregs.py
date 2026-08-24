"""Restore the files the filter-repo rewrite over-filtered out of rna-sa-public.

Run:  uv run --no-project --python 3.12 python scripts/restore_package_and_preregs.py [--write]

`rna-sa-public` was produced by rewriting `rna-structure-audit` to strip the
weight-geometry/DAS material, which belongs to a separate unpublished paper. The
filter removed more than it should have. Three things the RNA paper depends on
went with it:

1. `src/rna_structure_audit/` -- the entire package published to PyPI as
   rna-structure-audit 0.1.0, whose Homepage metadata points at the audit repo.
   `pyproject.toml` survived the rewrite; the code it packages did not.
2. `PREREGISTRATION_STRUCTURE_METRICS.md` -- the Phase 1 preregistration. It is
   the only file commit 694b43b touched, so dropping it pruned that commit to
   empty and the rewrite discarded it. paper_v10.tex Section 3.1 cites
   `694b43b` by SHA; the SHA now resolves in no repository that would survive
   deleting the audit repo.
3. `.github/workflows/publish.yml` -- the trusted-publisher workflow. Without it
   the package cannot be released again from this repo.

The weight-paper material stays out. `tests/test_rnafm.py` is excluded with it:
despite the RNA name it imports `lib.factorized_rnafm`, `lib.factorized_model`
and `lib.distillation`, none of which are in this repo.

This copies from the audit repo's working HEAD into the public tree. It does not
rewrite history -- another working copy has pinned this repo at 1c137616ddc3,
and a second rewrite would make that pin unreachable.
"""

import shutil
import subprocess
import sys
from pathlib import Path

GITHUB = Path.home() / "Documents" / "GitHub"
AUDIT = GITHUB / "rna-structure-audit"
PUBLIC = GITHUB / "rna-sa-public"

# Whole directories to restore, minus anything in EXCLUDE.
TREES = ["src", "tests", "notebooks", ".github/workflows"]

# Individual files: RNA-paper preregistrations and results the rewrite dropped.
FILES = [
    "PREREGISTRATION_STRUCTURE_METRICS.md",
    "PREREGISTRATION_PHASE6_V2_predraft.md",
    "docs/STRATIFIED_NULL_RESULTS.md",
]

# Weight-paper content that must not enter the public repo.
EXCLUDE = {"tests/test_rnafm.py"}


def tracked(repo: Path, *paths: str) -> list[str]:
    done = subprocess.run(
        ["git", "-C", str(repo), "ls-tree", "-r", "--name-only", "HEAD", "--", *paths],
        capture_output=True, text=True, check=True)
    return done.stdout.split()


def main() -> int:
    write = "--write" in sys.argv

    manifest = [p for p in tracked(AUDIT, *TREES) if p not in EXCLUDE]
    manifest += [p for p in FILES if (AUDIT / p).exists()]

    missing_source = [p for p in FILES if not (AUDIT / p).exists()]
    assert not missing_source, f"not in the audit repo: {missing_source}"

    already = [p for p in manifest if (PUBLIC / p).exists()]
    assert not already, f"already in the public repo, refusing to overwrite: {already}"

    excluded_present = [p for p in EXCLUDE if (PUBLIC / p).exists()]
    assert not excluded_present, f"weight-paper file already leaked: {excluded_present}"

    groups: dict[str, int] = {}
    for path in manifest:
        groups[path.split("/")[0]] = groups.get(path.split("/")[0], 0) + 1
    for top, count in sorted(groups.items()):
        print(f"  {count:3d}  {top}")
    print(f"\n{len(manifest)} files"
          f"{' -- copying' if write else ' -- dry run, pass --write to copy'}")

    if not write:
        return 0

    for path in manifest:
        target = PUBLIC / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(AUDIT / path, target)

    # The package must be importable from the restored tree, and the excluded
    # weight-paper test must not have come along with the directory copy.
    assert (PUBLIC / "src" / "rna_structure_audit" / "__init__.py").exists()
    assert (PUBLIC / "pyproject.toml").exists(), "pyproject.toml is missing"
    for path in EXCLUDE:
        assert not (PUBLIC / path).exists(), f"{path} was copied despite exclusion"
    leaked = subprocess.run(
        ["grep", "-rlIE", "factoriz|weight_geom|batch2",
         str(PUBLIC / "src"), str(PUBLIC / "tests"), str(PUBLIC / "notebooks")],
        capture_output=True, text=True).stdout.split()
    assert not leaked, f"weight-paper references in restored files: {leaked}"

    print(f"\nrestored {len(manifest)} files into {PUBLIC.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
