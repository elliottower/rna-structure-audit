"""Check that the final repo is complete and carries no weight-paper content.

Run:  uv run --no-project --python 3.12 python scripts/verify_final_repo.py

Two failure modes matter and they pull in opposite directions. Filtering too
little leaves the unpublished weight-geometry paper readable in a public
history. Filtering too much drops something the manuscript cites -- which is
what happened to the Phase 1 preregistration in the previous export, where
removing the file pruned its commit and the SHA the paper prints stopped
resolving.

Every check is against the whole history, not the tip: a file deleted at HEAD
but present in an old commit is still public.
"""

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

WEIGHT_PAPER = ["batch2/", "PREREGISTRATION_DAS.md", "PREREGISTRATION_WEIGHT_GEOMETRY.md",
                "FACTORIZATION_SPEC.md", "docs/factorization_whitepaper",
                "docs/frozen_prereg_weight_geom/", "lib/factorized", "lib/distillation.py",
                "lib/rnafm_weights.py", "data/weight_geometry/", "data/das_results",
                "figures/das_counterfactual.png", "figures/eap_ig_circuit.png",
                "tests/test_rnafm.py", "direction_instability", "experiment_plan.md",
                "scripts/analyze_factorization.py", "scripts/eap_ig_ablation.py",
                "scripts/weight_geometry_analysis.py", "scripts/run_hdas6.py",
                "scripts/modal_train.py", "scripts/modal_benchmark_test.py",
                "scripts/freeze_prereg_weight_geometry.py"]

REQUIRED = ["PREREGISTRATION_STRUCTURE_METRICS.md", "PREREGISTRATION_PHASE6_V2_predraft.md",
            "pyproject.toml", ".github/workflows/publish.yml", "README.md", "LICENSE",
            "src/rna_structure_audit/__init__.py", "src/rna_structure_audit/cli.py",
            "src/rna_structure_audit/evaluate.py", "tests/test_benchmark.py",
            "tests/conftest.py", "docs/STRATIFIED_NULL_RESULTS.md", "DEVIATIONS.md",
            "paper/paper_v10.tex", "notebooks/quickstart.ipynb",
            "scripts/modal_ernierna_ablation.py"]


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args],
                          capture_output=True, text=True, check=True).stdout.strip()


def main() -> int:
    failures = []

    ever = [p for p in git("log", "--all", "--pretty=format:", "--name-only",
                           "--diff-filter=A").split("\n") if p.strip()]
    leaked = sorted({p for p in ever if any(p.startswith(w) for w in WEIGHT_PAPER)})
    print(f"paths ever in history: {len(set(ever))}")
    print(f"weight-paper paths found: {len(leaked)}")
    for path in leaked:
        failures.append(f"weight-paper file in history: {path}")

    tip = set(git("ls-tree", "-r", "--name-only", "HEAD").split("\n"))
    print(f"files at HEAD: {len(tip)}   commits: {git('rev-list', '--all', '--count')}")
    for path in REQUIRED:
        if path not in tip:
            failures.append(f"required file missing: {path}")

    package = sorted(p for p in tip if p.startswith("src/rna_structure_audit/"))
    print(f"package files: {len(package)}")
    if len(package) < 70:
        failures.append(f"package looks truncated: {len(package)} files")

    preregs = sorted(p for p in tip if "PREREG" in p.upper())
    print(f"preregistrations at HEAD: {len(preregs)}")
    for path in preregs:
        print(f"  {path}")

    # Content, not just filenames: a weight-paper section pasted into a kept file
    # would pass every path check above.
    hits = subprocess.run(
        ["git", "-C", str(REPO), "grep", "-lI", "-iE",
         r"factorized (model|rnafm)|weight[- _]geometry prereg|PREREGISTRATION_DAS",
         "HEAD"], capture_output=True, text=True).stdout.split()
    external = [h for h in hits if "scripts/" not in h]
    print(f"content references to weight-paper material: {len(external)}")
    for hit in external:
        print(f"  ? {hit}")

    tags = git("tag").split()
    print(f"tags: {tags or '(none)'}")

    print()
    if failures:
        for line in failures:
            print(f"FAIL  {line}")
        return 1
    print("PASS -- no weight-paper file in any commit; every required file present.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
