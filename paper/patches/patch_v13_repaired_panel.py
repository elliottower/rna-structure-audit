r"""Build paper_v13_repaired_panel.tex from paper_v12.tex.

Run:  uv run --no-project --python 3.12 python paper/patches/patch_v13_repaired_panel.py

Two groups of edits, and the script refuses to run until both are decidable.

GROUP 1 -- the panel. Determined by the annotation repair alone, and reviewable
without the re-run. Five records are withdrawn and four annotations are repaired
against Rfam seed alignments, so the panel is 47 of 52 curated families. The
manuscript's class breakdown was carried by hand and never described the
deposited panel: `audit_panel_classes.py` finds five of the nine printed counts
wrong in classes whose membership the filenames fix -- tRNAs printed 7 against 2
deposited, rRNAs 4 against 1, riboswitches 8 against 14, CRISPR repeats 3 against
1, miRNA precursors 5 against 4. The printed counts sum to 52, which is how a
census wrong in five of nine entries survived four revisions. Every panel figure
therefore comes from \input{generated/panel_description} rather than from a
number typed into the body, so the next correction to the data corrects the
manuscript.

The panel holds one rRNA family before the repair and none after: 5S_rRNA_ecoli
is withdrawn because RF00001's seed alignment has no E. coli member and
substituting one would falsify the record's own name. The absence is stated in
the dataset description, in Limitations, and in DEVIATIONS.md.

GROUP 2 -- the results. Every rung is recomputed on the repaired panel with nulls
seeded from family names rather than panel positions, so no cell of any results
table survives and none can be carried over. The tables are replaced by
\input{generated/...} files that scripts/generate_results_tables.py writes from
the stamped run, for the same reason the panel description is generated: a
hand-carried number is a defect that reappears at the next correction.

The script asserts that the results it reads were produced at the repaired panel
by checking each stamp's panel_sha256 against a hash recomputed from
data/rfam_families/. A results directory left over from an earlier run has a
different hash and aborts the build rather than silently producing a v13 whose
tables describe a panel its Methods section does not.

v12 is not modified.
"""

import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "paper" / "paper_v12.tex"
DST = REPO / "paper" / "paper_v13_repaired_panel.tex"
GENERATED = REPO / "paper" / "generated"
RESULTS = REPO / "results" / "repaired_panel"
FAMILIES = REPO / "data" / "rfam_families"

# Written by generate_results_tables.py from the stamped run. Named here so a
# missing one aborts the build instead of producing a document that silently
# keeps a v12 table.
REQUIRED_GENERATED = [
    "panel_description.tex",
    "rung1_table.tex",
    "rung2_table.tex",
    "rung3_table.tex",
    "hypotheses_table.tex",
    "results_macros.tex",
]

DATASET = r"""\input{generated/panel_description}
\panelN{} RNA families drawn from Rfam~14.10~\citep{rfam2021}, spanning
\panelClasses{}.
Sequence lengths range from \panelLengthRange{}~nt (median
\panelLengthMedian{}~nt, IQR \panelLengthIQR{}).
This sample covers $<$2\% of the $>$4{,}000 Rfam families but
spans the major structural classes of well-characterized ncRNAs, with
two exceptions. Long-range pairing ($>$150~nt), multi-domain ribozymes,
and thermoswitches are underrepresented. Ribosomal RNA is absent: the
curated set contained a single rRNA record, and it is withdrawn
(\S\ref{sec:limitations}).

Of \panelCurated{} curated records, \panelWithdrawn{} are withdrawn and
the remaining \panelN{} are analyzed. Eleven records named no Rfam seed
member; nine of those carried a dot-bracket that does not describe the
sequence it is attached to, pairing 24--77\% of positions
non-canonically where all 41 records naming a seed member reach at most
18\%. The two distributions do not overlap. Each of the eleven was
resolved by a single rule applied before any model was run: adopt the
structure of the named seed member where one exists, withdraw the record
where none does, and leave records whose annotation is already supported.
Withdrawn records keep their file and carry an exclusion block naming the
reason, so the panel table below accounts for all \panelCurated{}.
Table~\ref{tab:panel} gives every record with its accession and
disposition. The repair is registered in \texttt{DEVIATIONS.md}.
"""


def panel_sha256() -> str:
    """The same hash modal_repaired_panel.py stamps, recomputed from the data."""
    records = [json.loads(p.read_text()) for p in sorted(FAMILIES.glob("*.json"))]
    loaded = sorted((r for r in records if "excluded" not in r), key=lambda r: r["name"])
    payload = json.dumps(
        [[r["name"], r["sequence"], r["dot_bracket"]] for r in loaded],
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


EDITS: list[tuple[str, str, str]] = [
    (
        "Methods: every panel figure comes from the generated file",
        """52 RNA families drawn from Rfam~14.10~\\citep{rfam2021}, spanning
tRNAs (7 families), rRNAs (4), ribozymes (6), riboswitches (8),
snRNAs (5), cis-regulatory elements (9), miRNA precursors (5),
CRISPR repeats (3), and other ncRNAs (5).
Sequence lengths range from 22 to 301~nt (median 89~nt, IQR 56--132).
This sample covers $<$2\\% of the $>$4{,}000 Rfam families but
spans the major structural classes of well-characterized ncRNAs.
Long-range pairing ($>$150~nt), multi-domain ribozymes, and
thermoswitches are underrepresented; extending the evaluation to
larger, more complex structures is a priority for future work.
Three families were excluded post hoc for unbalanced dot-bracket
annotations (ambiguous stem/loop assignments), leaving 49 evaluable
families for Rungs~1--2.
""",
        DATASET,
    ),
    (
        "Abstract: the panel size",
        "We evaluate ten models across 52 Rfam families using a three-rung",
        "We evaluate ten models across \\panelN{} Rfam families using a three-rung",
    ),
    (
        "Introduction: the panel size",
        "models across 52 Rfam families~\\citep{rfam2021}. The",
        "models across \\panelN{} Rfam families~\\citep{rfam2021}. The",
    ),
    (
        "Registered predictions: H1_6 carries both gates, because the repair "
        "raises the registered one",
        """\\item H1$_6$: At least one model has mean PS $> 0$ with $\\geq 7$
  families exceeding derangement null (Wilcoxon)""",
        """\\item H1$_6$: At least one model has mean PS $> 0$ with
  $\\geq \\lceil 4 \\times 0.05 N \\rceil$ families exceeding the derangement
  null (Wilcoxon). The gate was registered as 7 at $N = 32$; the repair
  admits three families that fabricated annotations had held out of the
  rung, and the gate is 8 at the repaired $N$. Both are reported.""",
    ),
    (
        "Limitations: the panel holds no rRNA family",
        """\\paragraph{Structure representation.}""",
        """\\paragraph{No ribosomal RNA.}
The curated set contained one rRNA record, a 5S sequence attributed to
\\emph{E. coli}. Its dot-bracket did not describe its sequence, and RF00001's
seed alignment contains no \\emph{E. coli} member from which to repair it, so
substituting a seed member would have produced a record whose name and contents
disagree. The record is withdrawn and the panel contains no rRNA. Ribosomal RNA
is the structural class with the most experimental support and the longest
range pairing in the ncRNA repertoire, so its absence bears on exactly the
regime the panel is already thin in. A 5S structure taken from a PDB entry
rather than from an Rfam seed is the cleanest single addition to the panel and
is the first extension we intend.

\\paragraph{Structure representation.}""",
    ),
    (
        "Limitations: the pseudoknot count is denominated in the analyzed panel",
        "other non-nested base pairings. Of the 52 families in our panel,",
        "other non-nested base pairings. Of the \\panelN{} families in our panel,",
    ),
    (
        "Limitations: the development set is the analyzed panel, and the class "
        "count comes from the generated census",
        """were developed on the same 52 families used for evaluation. A""",
        """were developed on the same \\panelN{} families used for evaluation. A""",
    ),
    (
        "Availability: what is deposited is every curated record, not only the "
        "analyzed ones",
        "scripts, per-family result files, 52 Rfam structures, and",
        "scripts, per-family result files, all \\panelCurated{} curated Rfam "
        "records (\\panelN{} analyzed, \\panelWithdrawn{} withdrawn with their "
        "exclusion reason), and",
    ),
]


def main() -> int:
    missing = [name for name in REQUIRED_GENERATED if not (GENERATED / name).exists()]
    if missing:
        print("ABORT -- generated tables absent; run "
              "scripts/generate_results_tables.py against the completed "
              f"re-run first:\n  " + "\n  ".join(missing))
        return 1

    expected = panel_sha256()
    stamps = sorted(RESULTS.glob("*/stamp.json"))
    if not stamps:
        print(f"ABORT -- no stamped results under {RESULTS.relative_to(REPO)}.")
        return 1
    wrong = [(p.parent.name, json.loads(p.read_text())["panel_sha256"])
             for p in stamps
             if json.loads(p.read_text())["panel_sha256"] != expected]
    if wrong:
        print("ABORT -- these results were produced on a different panel than "
              f"data/rfam_families/ now holds (expected {expected[:12]}):")
        for model, found in wrong:
            print(f"  {model:12s} {found[:12]}")
        return 1
    print(f"  {len(stamps)} models, all stamped with panel {expected[:12]}")

    text = SRC.read_text()
    for label, old, new in EDITS:
        count = text.count(old)
        if count != 1:
            print(f"ABORT -- [{count} matches] {label}")
            return 1
        text = text.replace(old, new, 1)
        print(f"  applied: {label}")

    # The hand-carried census is gone, in every class it named.
    for stale in ("tRNAs (7 families)", "rRNAs (4)", "riboswitches (8)",
                  "CRISPR repeats (3)", "spanning 8 major ncRNA classes"):
        assert stale not in text, f"the printed census survived: {stale}"
    assert text.count("\\input{generated/panel_description}") == 1
    assert "\\panelN{}" in text
    # Both gates reach the reader, and neither replaces the other.
    assert "registered as 7 at $N = 32$" in text
    assert "the gate is 8 at the repaired $N$" in text

    DST.write_text(text)
    print(f"\n{len(EDITS)} edits applied.\nwrote {DST.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
