"""Which result file produces each attention rho printed in Table 1?

Run:  uv run --no-project --with numpy --python 3.12 python \
          scripts/audit_attention_rho.py

The attention contact correlation is the one Table 1 column that
scripts/verify_paper_rung12_figures.py does not check, because it is not carried
in results/bootstrap_cis.json in the form the table prints. Two of its cells
carry footnotes saying they were recomputed after the artifact was written:
NT v2 as 0.322 "corrected from 0.230 (SDPA backend); eager attention computation
across 52 families", and DNABERT-2 as 0.257 from "content-only attention from
Wqkv hooks, excluding ALiBi position biases".

This script reads every result file in the repository, computes the mean of
`attention_trained.per_rna[*].best_corr` over scored families the way
scripts/compute_bootstrap_cis.py does, and reports which file, if any, prints as
each value the table gives. A printed value no file produces is a figure with no
record behind it.

Candidates are found by parsing every file rather than by globbing *.json: three
result files under data/gpu_results/expanded_rfam_rerun/ carry no extension, and
they are the 52-family reruns two of the footnotes refer to. Files are attributed
to a model by the `model` field they carry, not by their name.
"""

import json
from pathlib import Path

import numpy as np

from paper_versions import REPO, newest_paper

PAPER = newest_paper()
SEARCH = [REPO / "results", REPO / "data" / "gpu_results"]

# Printed in Table 1. "---" rows carry no attention measurement.
PRINTED = {
    "ERNIE-RNA": "0.100", "RiNALMo": "0.095", "RNA-FM": "0.044",
    "UTR-LM": "0.046", "SpliceBERT": "0.050", "NT v2": "0.322",
    "DNABERT-2": "0.257", "ERNIE-RNA untrained": "0.072",
    "RNA-FM untrained": "0.061",
}

# The `model` field each result file carries. Untrained controls are written
# either as their own model name or as the trained name in a path saying so.
MODEL_FIELD = {
    "ERNIE-RNA": "ernierna", "RiNALMo": "rinalmo", "RNA-FM": "rnafm",
    "UTR-LM": "utrlm", "SpliceBERT": "splicebert", "NT v2": "nt",
    "DNABERT-2": "dnabert2", "ERNIE-RNA untrained": "ernierna_untrained",
    "RNA-FM untrained": "rnafm_untrained",
}

# Two printed correlations have no per-family output in this repository, in the
# old tree it was rewritten from, or under any other aggregation of the files
# that are here. Both are recomputations the table footnotes describe: NT v2's
# with eager attention rather than the SDPA backend, and untrained RNA-FM's. The
# Availability section names them. This script fails on a gap not listed here.
KNOWN_UNRECORDED = {"NT v2", "RNA-FM untrained"}


def mean_rho(data: dict) -> tuple[float, int] | None:
    """Mean best_corr over scored families, as compute_bootstrap_cis computes it."""
    attention = data.get("attention_trained")
    if not isinstance(attention, dict) or attention.get("skipped"):
        return None
    corrs = [entry["best_corr"]
             for entry in attention.get("per_rna", {}).values()
             if isinstance(entry, dict) and not entry.get("skipped")
             and entry.get("best_corr") is not None]
    return (float(np.mean(corrs)), len(corrs)) if corrs else None


def candidates() -> list[tuple[Path, str, float, int]]:
    """Every file in the tree that parses as JSON and carries an attention block."""
    found = []
    for root in SEARCH:
        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue
            try:
                data = json.loads(path.read_text())
            except (json.JSONDecodeError, UnicodeDecodeError, OSError):
                continue
            if not isinstance(data, dict):
                continue
            result = mean_rho(data)
            if result:
                found.append((path, str(data.get("model", "")), *result))
    return found


def belongs(model_field: str, model: str) -> bool:
    return model_field == MODEL_FIELD[model]


def main() -> int:
    print(f"attention rho printed in {PAPER.name}, against every result file\n")
    files = candidates()
    print(f"{len(files)} result files carry an attention measurement\n")

    unmatched = []
    for model, printed in PRINTED.items():
        places = len(printed.split(".")[1])
        own = [(path, value, n) for path, field, value, n in files
               if belongs(field, model)]
        hits = [row for row in own if f"{row[1]:.{places}f}" == printed]
        print(f"{model:<22} {printed}")
        if hits:
            for path, value, n in hits[:3]:
                print(f"    {value:.4f} over {n} families  {path.relative_to(REPO)}")
        else:
            print("    no file in the repository prints as this value")
            for _, path, value, n in sorted(
                    (abs(value - float(printed)), path, value, n)
                    for path, value, n in own)[:3]:
                print(f"    nearest for this model: {value:.4f} over {n} families"
                      f"  {path.relative_to(REPO)}")
            unmatched.append(model)

    unexpected = sorted(set(unmatched) - KNOWN_UNRECORDED)
    recovered = sorted(KNOWN_UNRECORDED - set(unmatched))
    if unexpected:
        print(f"\n{len(unexpected)} printed values have no file behind them and "
              f"are not declared: {', '.join(unexpected)}")
        return 1
    if recovered:
        print(f"\n{', '.join(recovered)} now reproduces; drop it from "
              f"KNOWN_UNRECORDED and from the Availability section.")
        return 1
    print(f"\nEvery attention rho in Table 1 is produced by a file here, except "
          f"the {len(KNOWN_UNRECORDED)} the Availability section names.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
