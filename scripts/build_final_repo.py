"""Build the final public RNA repo: filtered audit history + the newer work.

Run:  uv run --no-project --python 3.12 python scripts/build_final_repo.py [--write]

Neither existing repo is complete.

`rna-structure-audit` holds all 36 commits, every preregistration, the package
published to PyPI, the tests and the trusted-publisher workflow -- but its last
commit is 2026-07-18, so it has none of the current manuscript (paper_v10.tex),
the H12-H21 deviations, the nine untrained controls or the multi-sequence
results. It also carries the weight-geometry/DAS material, which belongs to a
separate unpublished paper.

`rna-sa-public` holds that newer work, but it was produced by a filter-repo run
whose path list was wrong: it deleted the RNA Phase 1 preregistration, the whole
src/ package, the tests and the publish workflow along with the weight material.
Deleting the Phase 1 prereg pruned commit 694b43b -- which touched nothing else
-- to empty, so the rewrite discarded it, and the SHA the manuscript prints
resolves nowhere in that history.

This takes the audit repo, filters it with the corrected path list, and lays the
newer work on top. The Phase 1 preregistration survives the filter this time, so
its commit survives with it and all three preregistration SHAs the manuscript
cites have a counterpart in the result.

The source repos are not modified.
"""

import shutil
import subprocess
import sys
from pathlib import Path

GITHUB = Path.home() / "Documents" / "GitHub"
AUDIT = GITHUB / "rna-structure-audit"
PUBLIC = GITHUB / "rna-sa-public"
FINAL = GITHUB / "rna-structure-audit-FINAL"

# Everything belonging to the weight-geometry / DAS / factorization paper.
# Removed from every commit, not just the tip.
REMOVE = [
    "FACTORIZATION_SPEC.md",
    "PREREGISTRATION_DAS.md",
    "PREREGISTRATION_WEIGHT_GEOMETRY.md",
    "batch2/",
    "data/das_results_20260711_153152.json",
    "data/das_results_20260711_160612.json",
    "data/weight_geometry/",
    "direction_instability.py",
    "direction_instability_v1.py",
    "docs/batch2_results_writeup.md",
    "docs/factorization_whitepaper_v1.md",
    "docs/factorization_whitepaper_v2.md",
    "docs/frozen_prereg_weight_geom/",
    "experiment_plan.md",
    "figures/das_counterfactual.png",
    "figures/eap_ig_circuit.png",
    "lib/distillation.py",
    "lib/factorized_model.py",
    "lib/factorized_rnafm.py",
    "lib/rnafm_weights.py",
    "scripts/analyze_factorization.py",
    "scripts/eap_ig_ablation.py",
    "scripts/freeze_prereg_weight_geometry.py",
    "scripts/modal_benchmark_test.py",
    "scripts/modal_train.py",
    "scripts/run_hdas6.py",
    "scripts/weight_geometry_analysis.py",
    "tests/test_rnafm.py",
]

# Kept deliberately, despite matching a weight-ish name. The ERNIE-RNA ablation
# zeroes that model's pairwise bias buffer to ask whether its structure
# awareness comes from the bias or the learned weights -- an RNA result.
KEEP_ANYWAY = ["scripts/modal_ernierna_ablation.py", "results/ernierna_ablation/"]

# Files the manuscript and the package cannot lose.
MUST_SURVIVE = [
    "PREREGISTRATION_STRUCTURE_METRICS.md",
    "PREREGISTRATION_PHASE6_V2_predraft.md",
    "pyproject.toml",
    ".github/workflows/publish.yml",
    "src/rna_structure_audit/__init__.py",
    "src/rna_structure_audit/cli.py",
    "tests/test_benchmark.py",
    "docs/STRATIFIED_NULL_RESULTS.md",
]

CITED_SHAS = {"Phase 1": "694b43b", "Phase 2": "bd4b3fd", "Phase 6": "c19aa59"}
AUDIT_LAST_DAY = "2026-07-19"  # the audit repo's history stops on the 18th


def git(repo: Path, *args: str, check: bool = True) -> str:
    done = subprocess.run(["git", "-C", str(repo), *args],
                          capture_output=True, text=True)
    if check and done.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed:\n{done.stderr}")
    return done.stdout.strip()


def newer_from_public() -> list[str]:
    """Paths rna-sa-public added or changed after the audit history ends."""
    changed = set(git(PUBLIC, "log", f"--since={AUDIT_LAST_DAY}", "--pretty=format:",
                      "--name-only").split("\n"))
    tracked = set(git(PUBLIC, "ls-tree", "-r", "--name-only", "HEAD").split("\n"))
    audit_ever = set(git(AUDIT, "log", "--all", "--pretty=format:", "--name-only",
                         "--diff-filter=A").split("\n"))
    never_in_audit = tracked - audit_ever
    candidates = sorted(p for p in (changed | never_in_audit) & tracked if p.strip())

    # Most of what the "changed since" filter catches is the restore commit, which
    # copied files back out of the audit repo unchanged. Those already arrive
    # through the filtered history; carrying them again would claim as new work
    # what is byte-for-byte the audit's own content.
    fresh = []
    for path in candidates:
        here = git(PUBLIC, "rev-parse", f"HEAD:{path}", check=False)
        there = git(AUDIT, "rev-parse", f"HEAD:{path}", check=False)
        if here != there:
            fresh.append(path)
    return fresh


def main() -> int:
    write = "--write" in sys.argv
    carry = newer_from_public()
    dropped = [p for p in carry if any(p.startswith(r) for r in REMOVE)]
    carry = [p for p in carry if p not in dropped]

    print(f"removing {len(REMOVE)} weight-paper paths from all history")
    print(f"carrying {len(carry)} newer files over from rna-sa-public")
    if dropped:
        print(f"  ({len(dropped)} of them are weight-paper and stay out: {dropped})")
    for path in carry:
        print(f"    {path}")
    if not write:
        print("\ndry run -- pass --write to build")
        return 0

    assert not FINAL.exists(), f"{FINAL} exists; move it aside first"

    subprocess.run(["git", "clone", "--no-local", str(AUDIT), str(FINAL)], check=True)
    listing = FINAL / ".git" / "filter-paths.txt"
    listing.write_text("\n".join(REMOVE) + "\n")
    subprocess.run(["git", "filter-repo", "--invert-paths",
                    "--paths-from-file", str(listing), "--force"],
                   cwd=str(FINAL), check=True)
    listing.unlink()

    for path in carry:
        target = FINAL / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(PUBLIC / path, target)
    git(FINAL, "add", *carry)
    git(FINAL, "commit", "-q", "-m",
        "Add the current manuscript and the results recorded after 2026-07-18\n\n"
        "The audit history this repo is filtered from ends on 2026-07-18. The\n"
        "manuscript (paper_v10.tex), the H12-H21 deviations, the nine untrained\n"
        "controls and the multi-sequence results were recorded after that date in a\n"
        "separate working copy. They are brought over here unchanged.")

    # No weight-paper file may survive anywhere in the rewritten history.
    ever = set(git(FINAL, "log", "--all", "--pretty=format:", "--name-only",
                   "--diff-filter=A").split("\n"))
    leaked = sorted(p for p in ever if any(p.startswith(r) for r in REMOVE))
    assert not leaked, f"weight-paper files survived the filter: {leaked}"

    tip = set(git(FINAL, "ls-tree", "-r", "--name-only", "HEAD").split("\n"))
    for path in MUST_SURVIVE:
        assert path in tip, f"{path} did not survive the filter"
    for path in KEEP_ANYWAY:
        assert any(p.startswith(path) for p in tip), f"{path} was removed by mistake"
    assert "paper/paper_v10.tex" in tip, "the current manuscript is missing"

    print(f"\n{len(tip)} files at HEAD, {git(FINAL, 'rev-list', '--all', '--count')} commits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
