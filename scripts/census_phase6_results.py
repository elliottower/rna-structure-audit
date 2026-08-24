"""Census every stored Phase 6 perturbation-specificity result across the RNA repos.

Run:  uv run python scripts/census_phase6_results.py

Table 5 reports mean perturbation specificity for eleven models. The result files
those numbers came from are spread over four repositories and, within this one,
over three directories, with more than one file per model in some cases. A number
that appears once in the manuscript and twice on disk cannot be reproduced, so
this lists every file that carries a value, which repository it is in, and what
it says.

Read-only. Writes nothing, and does not decide which value is correct.
"""

import json
import re
from pathlib import Path

GITHUB = Path.home() / "Documents" / "GitHub"
REPOS = [
    "rna-sa-public",
    "rna-sa-public-SUPERSEDED",
    "rna-structure-audit",
    "rna-structure-awareness",
]

# Table 5's rows, and the substrings that identify each model's files on disk.
MODELS = {
    "RiNALMo": ("rinalmo",),
    "ERNIE-RNA": ("ernierna", "ernie_rna"),
    "Caduceus": ("caduceus",),
    "Evo": ("evo",),
    "HyenaDNA": ("hyenadna",),
    "SpliceBERT": ("splicebert",),
    "RNA-FM": ("rnafm", "rna_fm"),
    "UTR-LM": ("utrlm", "utr_lm"),
    "NT v2": ("_nt_", "/nt_", "nucleotide"),
    "DNABERT-2": ("dnabert",),
}
VALUE_KEYS = ("mean_best_ps", "mean_ps", "mean_perturbation_specificity")


def value_of(path: Path) -> tuple[str, object] | None:
    try:
        payload = json.loads(path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError, OSError):
        return None
    for scope in (payload, payload.get("results") if isinstance(payload, dict) else None):
        if isinstance(scope, dict):
            for key in VALUE_KEYS:
                if key in scope:
                    return key, scope[key]
    return None


def model_of(relative: str) -> str | None:
    lowered = relative.lower()
    hits = [name for name, needles in MODELS.items()
            if any(needle in lowered for needle in needles)]
    # "evo" matches inside other words; prefer a longer, more specific match.
    if len(hits) > 1:
        hits = [h for h in hits if h != "Evo"] or hits
    return hits[0] if len(hits) == 1 else None


def main() -> int:
    rows: dict[str, list[tuple[str, str, str, object, bool]]] = {n: [] for n in MODELS}
    for repo in REPOS:
        root = GITHUB / repo
        if not root.exists():
            print(f"  (missing repo: {repo})")
            continue
        for path in root.rglob("*.json"):
            if ".git/" in str(path):
                continue
            relative = str(path.relative_to(root))
            if "phase6" not in relative.lower() and "_ps" not in relative.lower():
                continue
            found = value_of(path)
            if found is None:
                continue
            name = model_of(relative)
            if name is None:
                continue
            untrained = "untrained" in relative.lower()
            rows[name].append((repo, relative, found[0], found[1], untrained))

    for name in MODELS:
        entries = sorted(rows[name], key=lambda e: (e[4], e[0], e[1]))
        trained = [e for e in entries if not e[4]]
        distinct = {repr(e[3]) for e in trained}
        flag = "  <-- DISAGREEING FILES" if len(distinct) > 1 else ""
        print(f"\n=== {name}: {len(trained)} trained file(s), "
              f"{len(distinct)} distinct value(s){flag}")
        for repo, relative, key, value, untrained in entries:
            label = "untrained" if untrained else "trained  "
            print(f"    {label}  {value!r:<26} [{repo}] {relative}")
        if not entries:
            print("    (no stored file in any repository)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
