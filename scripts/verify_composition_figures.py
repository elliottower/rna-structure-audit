"""Recompute the stem/loop GC figures from the deposited family annotations.

Run:  uv run --no-project --with numpy --python 3.12 python \
          scripts/verify_composition_figures.py

Two manuscripts print different composition confounds. paper_v12.tex reports
68.5% GC in stems against 43.9% in loops, a 24.7 percentage-point enrichment,
with no family count. The submitted rna-structure-audit_v14.tex reports 59.7%
against 45.6%, 14.1 points, across the 52 Rfam families, with 42 of 52 above
the diagonal.

The 52 family annotations are in data/rfam_families/, so the enrichment is a
lookup: stem positions are the dot-bracket parentheses, loop positions the
dots. PREREGISTRATION_STRUCTURE_METRICS.md attributes 68.5/43.9 to the
12-family pilot set that preceded the Rfam expansion; that set is not
enumerated anywhere in the repository, so this script measures the 52 and
leaves the pilot figure to the registration that states it.

No model is run. The annotations are the stored inputs.
"""

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data" / "rfam_families"

# (label, stem GC, loop GC) as each manuscript prints them.
PRINTED = [
    ("paper_v12.tex", 68.5, 43.9),
    ("rna-structure-audit_v14.tex", 59.7, 45.6),
]


def per_family_gc() -> tuple[np.ndarray, np.ndarray, list[str]]:
    stems, loops, names = [], [], []
    for path in sorted(FAMILIES.glob("*.json")):
        record = json.loads(path.read_text())
        # Withdrawn records stay on disk and out of the panel; see DEVIATIONS.md.
        if "excluded" in record:
            continue
        sequence = record["sequence"].upper()
        brackets = record["dot_bracket"]
        stem_positions = [i for i, c in enumerate(brackets) if c in "()"]
        loop_positions = [i for i, c in enumerate(brackets) if c == "."]
        if not stem_positions or not loop_positions:
            continue
        stems.append(sum(sequence[i] in "GC" for i in stem_positions)
                     / len(stem_positions) * 100)
        loops.append(sum(sequence[i] in "GC" for i in loop_positions)
                     / len(loop_positions) * 100)
        names.append(record["name"])
    return np.array(stems), np.array(loops), names


def main() -> int:
    stems, loops, names = per_family_gc()
    stem_mean, loop_mean = stems.mean(), loops.mean()
    enrichment = stem_mean - loop_mean
    above = int(np.sum(stems > loops))

    print(f"families with both stem and loop positions: {len(names)}")
    print(f"mean stem GC: {stem_mean:.1f}%")
    print(f"mean loop GC: {loop_mean:.1f}%")
    print(f"enrichment:   {enrichment:+.1f} percentage points")
    print(f"stem GC > loop GC in {above} of {len(names)} families\n")

    for label, printed_stem, printed_loop in PRINTED:
        reproduces = (f"{stem_mean:.1f}" == f"{printed_stem:.1f}"
                      and f"{loop_mean:.1f}" == f"{printed_loop:.1f}")
        verdict = "reproduces" if reproduces else "DOES NOT REPRODUCE"
        print(f"{label:<30} {printed_stem}/{printed_loop} "
              f"({printed_stem - printed_loop:+.1f} pp)  ->  {verdict}")

    matching = [label for label, s, l in PRINTED
                if f"{stem_mean:.1f}" == f"{s:.1f}" and f"{loop_mean:.1f}" == f"{l:.1f}"]
    print()
    if matching:
        print(f"the annotations on this branch are the ones {matching[0]} was "
              f"computed from")
    else:
        print("no manuscript reproduces from these annotations. On main this is "
              "the expected\nstate: commit 5fc8914 corrected four family "
              "annotations and no run has been\nrepeated, so results/ still "
              "carries values computed against the pre-fix data.\nThe "
              "as-submitted branch holds those annotations.")
    assert len(matching) <= 1, f"two manuscripts cannot both match: {matching}"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
