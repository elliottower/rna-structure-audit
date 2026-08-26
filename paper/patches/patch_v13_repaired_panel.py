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
table survives and none can be carried over. The five results tables are
replaced whole by \input{generated/...} files that
scripts/generate_results_tables.py writes from the stamped runs, for the same
reason the panel description is generated: a hand-carried number is a defect
that reappears at the next correction. The labels are preserved by the generated
files, so every \ref in the body still resolves.

Two tables are added to the appendix: the curated panel with each record's
accession and disposition, and the provenance of all seventeen runs. The
provenance table exists because no stack loads all ten models, so the runs carry
different commits and different library versions, and a reader given one commit
for the whole results section would be given a commit that does not describe most
of it.

The macros both generators write are \input in the preamble rather than in
Methods, because the abstract prints the panel size and reaches LaTeX before
Methods does.

The script asserts that the results it reads were produced at the repaired panel
by checking each stamp's panel_sha256 against a hash recomputed from
data/rfam_families/. A results directory left over from an earlier run has a
different hash and aborts the build rather than silently producing a v13 whose
tables describe a panel its Methods section does not.

v12 is not modified.
"""

import argparse
import difflib
import hashlib
import json
import re
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
    "panel_table.tex",
    "rung1_table.tex",
    "rung2_table.tex",
    "rung3_table.tex",
    "attention_table.tex",
    "hypotheses_table.tex",
    "provenance_table.tex",
    "transversion_table.tex",
    "results_macros.tex",
]

# Each v12 table environment, and the file that replaces it. The generated files
# carry the same \label, so the swap does not orphan a \ref.
GENERATED_TABLES = [
    ("tab:rung1", "rung1_table"),
    ("tab:attn_trained_untrained", "attention_table"),
    ("tab:rung2", "rung2_table"),
    ("tab:rung3", "rung3_table"),
    ("tab:hypotheses", "hypotheses_table"),
]

DATASET = r"""\panelN{} RNA families drawn from Rfam~14.10~\citep{rfam2021}, spanning
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


def replace_table(text: str, label: str, generated: str) -> str:
    r"""Swap a v12 table environment for \input{generated/<generated>}.

    Anchored on the \label, not on the table's contents: the contents are the
    thing being replaced, and reproducing a caption here to find them would be
    one more copy of a number to keep in step. `\end{tabular}` cannot be
    mistaken for `\end{table}` -- the brace has to follow `table` immediately.
    """
    marker = "\\label{" + label + "}"
    if text.count(marker) != 1:
        raise SystemExit(f"{label}: {text.count(marker)} labels in v12, expected 1")
    at = text.index(marker)
    start = text.rindex("\\begin{table}", 0, at)
    end = text.index("\\end{table}", at) + len("\\end{table}")
    return text[:start] + "\\input{generated/" + generated + "}" + text[end:]


EDITS: list[tuple[str, str, str]] = [
    (
        "Preamble: longtable, and the generated macros before the abstract uses them",
        """\\usetikzlibrary{decorations.pathreplacing}
""",
        """\\usetikzlibrary{decorations.pathreplacing}
\\usepackage{longtable}

% Every figure this manuscript prints about the panel or the runs is defined
% here, computed from data/rfam_families/ and results/repaired_panel/ by
% scripts/generate_panel_description.py and scripts/generate_results_tables.py.
% In the preamble because the abstract prints the panel size.
\\input{generated/panel_description}
\\input{generated/results_macros}
""",
    ),
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
        "Limitations: the class count is printed once, in the generated census",
        """diversity criteria (spanning 8 major ncRNA classes) with no reference""",
        """diversity criteria---the classes and their sizes are given in
Section~\\ref{sec:dataset}---with no reference""",
    ),
    (
        "Methods: the dataset section is referred to from Limitations",
        """\\subsection{RNA structure dataset}""",
        """\\subsection{RNA structure dataset}
\\label{sec:dataset}""",
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
        "Results: the attention section reports the measurement, not its history",
        """An earlier analysis using HuggingFace's default SDPA attention
backend reported $\\rho = 0.000$ for untrained NT~v2. The corrected
analysis (eager attention; see Section~\\ref{sec:methods-attn})
reverses this interpretation entirely.

""",
        "",
    ),
    (
        "Results: the hypothesis summary is a subsection, not a bare float",
        """\\subsection{Preregistered hypothesis summary}

\\begin{table}[htbp]""",
        """\\subsection{Preregistered hypothesis summary}

Table~\\ref{tab:hypotheses} gives every hypothesis in the four
registrations with its criterion, the panel it is decided on and its
outcome.

\\begin{table}[htbp]""",
    ),
    (
        "Availability: the re-run deposits per-family output for every run",
        """ Two attention correlations in
Table~\\ref{tab:rung1} come from runs whose per-family output is not
among the deposited files: the eager-attention recomputation for
NT~v2 and the attention run for untrained RNA-FM.""",
        "",
    ),
    (
        "Appendix: the curated panel and the provenance of every run",
        """\\section*{Declarations}""",
        """\\section{The curated panel}
\\label{app:panel}

\\input{generated/panel_table}

\\section{Provenance}
\\label{app:provenance}

No library stack loads all ten models. DNABERT-2's bundled Triton attention
requires a Torch and Triton pair that the multimolecule build the RNA models
need will not run under, and Evo requires a third combination again. The runs
behind the results tables were therefore launched in groups, under
\\provenanceStacks{} library stacks and \\provenanceCommits{} commits.
Table~\\ref{tab:provenance} gives the commit, the weights, the device and the
library stack for each of the \\provenanceRuns{} runs, so any cell of any table
can be attributed to the container that produced it.

\\input{generated/provenance_table}

\\section*{Declarations}""",
    ),
    (
        "Rung 1: the counts, the intervals and the probing peaks come from the run",
        """Two patterns stand out. First, ERNIE-RNA exceeds the null in 31 of
49 families---the highest consistency of any model---with the
strongest probing signal ($\\Delta = +0.166$, balanced accuracy 0.666
at layer~11). RNA-FM has the highest mean ratio (1.749) but exceeds
the null in only 2 of 49 families: high mean, high variance, low
per-family reliability. At $\\alpha = 0.05$ per family, 2.5 of 49
families are expected to exceed the null by chance. The untrained
controls (3/49 for both ERNIE-RNA and RNA-FM) are consistent with
this false-positive rate, validating the null calibration.
Bootstrap 95\\% CIs (Table~\\ref{tab:rung1}) sharpen the picture:
ERNIE-RNA's CI $[1.19, 1.31]$ and RiNALMo's $[1.10, 1.15]$ exclude
1.0, confirming reliable elevation. DNABERT-2's CI $[0.96, 1.03]$
straddles 1.0, and models with 3--4 families exceeding the null
(UTR-LM, SpliceBERT) sit barely above the 2.5 expected false-positive
count. RiNALMo follows ERNIE-RNA at 18/49, with the second-strongest
probing signal ($\\Delta = +0.123$, balanced accuracy 0.623 at
layer~32). The per-family heatmap
(Figure~\\ref{fig:heatmap}) shows that ERNIE-RNA's signal is
distributed broadly across families rather than driven by a few
outliers.""",
        """ERNIE-RNA exceeds the null in \\exceedErnieRNA{} of \\familiesScored{}
families, more than any other model, and carries the strongest probing
signal ($\\Delta = \\probeDeltaErnieRNA{}$, balanced accuracy
\\probeErnieRNA{} at layer~\\probeLayerErnieRNA{}). RNA-FM has the
highest mean ratio (\\ratioRNAFM{}) and exceeds the null in
\\exceedRNAFM{} families: a high mean over high per-family variance.
At $\\alpha = 0.05$ per family, \\expectedFalsePositives{} of
\\familiesScored{} families are expected to exceed the null by chance,
and the untrained controls exceed it in \\exceedUntrainedErnieRNA{}
(ERNIE-RNA) and \\exceedUntrainedRNAFM{} (RNA-FM) families, which is
what the null calibration predicts. Bootstrap 95\\% CIs on the mean
ratio bound the elevation (Table~\\ref{tab:rung1}): ERNIE-RNA's
$\\ciErnieRNA{}$ and RiNALMo's $\\ciRiNALMo{}$ exclude 1.0, while
DNABERT-2's $\\ciDNABERTTwo{}$ straddles it. UTR-LM and SpliceBERT
exceed the null in \\exceedUTRLM{} and \\exceedSpliceBERT{} families,
against the \\expectedFalsePositives{} expected by chance. RiNALMo
follows ERNIE-RNA at \\exceedRiNALMo{} families, with the
second-strongest probing signal ($\\Delta = \\probeDeltaRiNALMo{}$,
balanced accuracy \\probeRiNALMo{} at layer~\\probeLayerRiNALMo{}).
The per-family heatmap (Figure~\\ref{fig:heatmap}) shows ERNIE-RNA's
signal spread across families rather than carried by a few
outliers.""",
    ),
    (
        "Rung 1: the sign test counts families from the run",
        """Beyond the per-family null exceedances, the trained-versus-untrained
sign test provides a more sensitive measure. Across all 52 Rfam
families, ERNIE-RNA trained exceeds its untrained counterpart in
52/52 families (100\\%). RiNALMo exceeds in 49/52 (94\\%). SpliceBERT
shows 48/52 (92\\%) despite negligible effect size. The directional
effect is consistent across RNA-pretrained models; the effect
magnitude, rather than its existence, separates the leaders.""",
        """The trained-versus-untrained sign test counts families rather than
comparing means, so a shift too small to move a mean registers where
it is consistent across the panel. ERNIE-RNA trained exceeds its
untrained counterpart in \\signErnieRNA{} families
(\\signPctErnieRNA{}), RiNALMo in \\signRiNALMo{}
(\\signPctRiNALMo{}), and SpliceBERT in \\signSpliceBERT{}
(\\signPctSpliceBERT{}) at an effect size near zero. The direction is
shared across the RNA-pretrained models; the magnitude separates
them.""",
    ),
    (
        "Rung 1: NT v2's attention deltas come from the run",
        """NT~v2 produces the highest attention-contact correlation of any model
($\\rho = 0.322$; Table~\\ref{tab:attn_trained_untrained}). Comparing
trained and untrained (randomly initialized) weights reveals that this
signal is entirely architectural: untrained NT~v2 produces
$\\rho = 0.328$---slightly \\emph{higher} than trained.
The per-family delta (trained $-$ untrained) averages $-0.006$
with no systematic direction: 22 families show higher trained
correlation, 30 show higher untrained. NT~v2's 6-mer tokenization
explains the effect.""",
        """NT~v2 produces the highest attention-contact correlation of any model
($\\rho = \\attnNTvTwo{}$; Table~\\ref{tab:attn_trained_untrained}).
Untrained NT~v2 produces $\\rho = \\attnUntrainedNTvTwo{}$, above the
trained value, so the signal is architectural. The per-family delta
(trained $-$ untrained) averages $\\attnFamilyDeltaNTvTwo{}$ with no
systematic direction: \\attnHigherTrainedNTvTwo{} families show higher
trained correlation and \\attnHigherUntrainedNTvTwo{} higher
untrained. NT~v2's 6-mer tokenization accounts for the effect.""",
    ),
    (
        "Rung 1: DNABERT-2's attention and ratio come from the run, and it is not "
        "character-tokenized",
        """DNABERT-2 is the exception among character-tokenized models: its
content-only attention ($\\rho = 0.257$ trained vs $0.000$ untrained)
shows learned structure-aware content patterns. Random Wqkv weights
produce uniform content attention (all Spearman correlations zero)
because random query-key products are unstructured.
DNABERT-2's mean ratio falls below 1.0, reflecting BPE tokenization
that merges nucleotides into multi-character tokens;
single-nucleotide mutations alter token boundaries unpredictably,
confounding the embedding comparison.""",
        """DNABERT-2 is the exception: its content-only attention
($\\rho = \\attnDNABERTTwo{}$ trained against $\\attnUntrainedDNABERTTwo{}$
untrained) shows learned structure-aware content patterns. Random
Wqkv weights produce uniform content attention---every Spearman
correlation zero---because random query-key products are
unstructured. DNABERT-2's mean ratio (\\ratioDNABERTTwo{}) falls
below 1.0, reflecting a BPE tokenizer that merges nucleotides into
multi-character tokens: a single-nucleotide mutation alters token
boundaries unpredictably, which confounds the embedding
comparison.""",
    ),
    (
        "Rung 1: the character-tokenized RNA deltas come from the run",
        """(Table~\\ref{tab:attn_trained_untrained}). RiNALMo and ERNIE-RNA
show the largest deltas ($+0.019$ and $+0.024$ respectively), but
both are small. RNA-FM attention carries no structural information
(0.051 trained vs 0.060 untrained). Attention-contact correlation
and embedding-level structure encoding operate as independent
channels: RNA-FM ranks first on embedding mutation sensitivity but
last on attention, while NT~v2 ranks first on attention but
mid-pack on embeddings.""",
        """(Table~\\ref{tab:attn_trained_untrained}). RiNALMo and ERNIE-RNA show
the largest deltas ($\\attnDeltaRiNALMo{}$ and $\\attnDeltaErnieRNA{}$),
both small. RNA-FM attention carries no structural information
(\\attnRNAFM{} trained against \\attnUntrainedRNAFM{} untrained).
Attention-contact correlation and embedding-level structure encoding
operate as independent channels: RNA-FM ranks first on embedding
mutation sensitivity and last on attention, while NT~v2 ranks first
on attention and mid-pack on embeddings.""",
    ),
    (
        "Rung 1: Evo against Caduceus comes from the run",
        """Evo (7B, byte-level, DNA) achieves a higher mean ratio (1.408) than
Caduceus (14M, character, DNA; 1.211), but both models' probing
accuracy peaks at layer~0--1---consistent with positional encoding
rather than learned computation. Evo's $500\\times$ parameter
advantage yields only marginal gains in mutation sensitivity and no
advantage in partner specificity (Table~\\ref{tab:rung3}).""",
        """Evo (7B, byte-level, DNA) reaches a higher mean ratio (\\ratioEvo{})
than Caduceus (14M, character, DNA; \\ratioCaduceus{}), and both
models' probing accuracy peaks in the first two layers
(layer~\\probeLayerEvo{} for Evo, layer~\\probeLayerCaduceus{} for
Caduceus), consistent with positional encoding rather than learned
computation. Evo's $500\\times$ parameter advantage yields a marginal
gain in mutation sensitivity and none in partner specificity
(Table~\\ref{tab:rung3}).""",
    ),
    (
        "Rung 1: the transversion control becomes a table, and the alphabet is "
        "described by what it does",
        """A reviewer concern is that Watson-Crick complement swaps (A$\\leftrightarrow$U,
C$\\leftrightarrow$G) at stem positions might create reversed but still
valid base pairs (e.g., A-U $\\to$ U-A), underestimating structural
disruption. To test this, we repeated the mutation sensitivity
analysis using transversion mutations (purine $\\to$ pyrimidine),
which guarantee disrupted pairing, on all ten models.
Transversion and WC-complement ratios match closely
across the full panel:
ERNIE-RNA 1.227 vs 1.242 ($-1.2\\%$), RiNALMo 1.129 vs 1.128,
RNA-FM 1.693 vs 1.749 ($-3.2\\%$), Evo 1.434 vs 1.408 ($+1.9\\%$),
Caduceus 1.154 vs 1.211 ($-4.7\\%$), NT~v2 1.120 vs 1.103 ($+1.6\\%$),
SpliceBERT 1.046 vs 1.040, UTR-LM 1.055 vs 1.054,
DNABERT-2 1.008 vs 0.992. HyenaDNA shows the largest divergence
(1.077 vs 1.161, $-7.2\\%$), suggesting that a fraction of its
original signal reflects complement-symmetric features---consistent
with a long-range convolutional model learning A$\\leftrightarrow$U
and G$\\leftrightarrow$C symmetries from DNA pretraining.
All ten models fall within 5\\% of their WC-complement ratios,
confirming that WC-complement mutations are genuine structure
disruptors: reversed base pairs differ from correct pairs enough
to register as embedding perturbations.""",
        """The Rung~1 substitution replaces a nucleotide with its Watson-Crick
complement, which is itself a pairing base, so a model responding to
complement-symmetric sequence features rather than to pairing would
register the substitution for the wrong reason. We repeated the
mutation sensitivity analysis under an alphabet chosen so that no
substituted nucleotide pairs with the partner its position already
had (\\transversionAlphabet{}), on all ten models. Every model's mean
ratio moves by at most \\transversionBound{}
(Table~\\ref{tab:transversion}). \\transversionWidest{} moves most, at
\\transversionWidestDeviation{}, consistent with a fraction of its
signal reflecting complement-symmetric features that a long-range
convolutional model can pick up from DNA pretraining. The
Watson-Crick substitution is a structural disruption rather than a
complement-symmetric artifact.

\\input{generated/transversion_table}""",
    ),
    (
        "Rung 1: the heatmap caption counts the panel from the panel",
        """\\caption{Per-family mutation sensitivity ratio (log$_2$ scale) across
10 models and 49 Rfam families (sorted by mean ratio, descending).""",
        """\\caption{Per-family mutation sensitivity ratio (log$_2$ scale) across
ten models and \\panelN{} Rfam families (sorted by mean ratio,
descending).""",
    ),
    (
        "Rung 2: the retention counts come from the run",
        """The dinucleotide null produces a sharp separation between models.
ERNIE-RNA retains 30 of 31 first-order survivors (97\\%): its
embedding-level structure signal persists after controlling for both
nucleotide and dinucleotide composition. RiNALMo retains 16 of 18
(89\\%). Three DNA models---NT~v2 (11/14), HyenaDNA (8/10), and Evo
(7/10)---retain the majority of their families, suggesting
composition-independent structure encoding despite DNA-only
pretraining.

The remaining models collapse. Caduceus retains 2 of 8 (25\\%),
DNABERT-2 retains 2 of 7 (29\\%), and SpliceBERT retains 2 of 4.
ERNIE-RNA with randomized weights retains 1 of 3---consistent with
the architectural attention bias providing capacity without content.""",
        """The dinucleotide null separates the models. ERNIE-RNA retains
\\retainErnieRNA{} first-order survivors (\\retainPctErnieRNA{}): its
embedding-level structure signal persists after controlling for both
nucleotide and dinucleotide composition. RiNALMo retains
\\retainRiNALMo{} (\\retainPctRiNALMo{}). Three DNA models---NT~v2
(\\retainNTvTwo{}), HyenaDNA (\\retainHyenaDNA{}) and Evo
(\\retainEvo{})---retain the majority of their families, which
suggests composition-independent structure encoding despite DNA-only
pretraining.

The remaining models collapse. Caduceus retains \\retainCaduceus{}
(\\retainPctCaduceus{}), DNABERT-2 \\retainDNABERTTwo{}
(\\retainPctDNABERTTwo{}) and SpliceBERT \\retainSpliceBERT{}.
ERNIE-RNA with randomized weights retains
\\retainUntrainedErnieRNA{}, consistent with the architectural
attention bias providing capacity without content.""",
    ),
    (
        "Rung 2: the overview caption reads its retentions from the run",
        """ERNIE-RNA retains 30/31 (97\\%), RiNALMo 16/18 (89\\%).""",
        """ERNIE-RNA retains \\retainErnieRNA{} (\\retainPctErnieRNA{}), RiNALMo
\\retainRiNALMo{} (\\retainPctRiNALMo{}).""",
    ),
    (
        "Rung 3: the PS values, gate counts and H3 precisions come from the run",
        """RiNALMo achieves the highest perturbation specificity (PS $= 0.210$),
with 28 of 29 gate-passing families exceeding the derangement null
and a partner-is-max precision of 0.874.
ERNIE-RNA follows at PS $= 0.114$, with 28/30 families exceeding the
null and H3 precision of 0.870. Both are more than an order
of magnitude above the third model and far exceed the
$1/3$ chance baseline ($p \\ll 0.001$, binomial test).

Caduceus (14M, DNA, BiMamba SSM) ranks third at PS $= 0.004$, with
13/16 gate-passing families exceeding the null and H3 precision of
0.416---above chance but 32-fold below ERNIE-RNA and 60-fold below
RiNALMo. The remaining models cluster near zero PS. SpliceBERT
exceeds the null in 14 of 21 gate-passing families, but its PS
magnitude (0.0002) is roughly three orders of magnitude below the two
leaders. HyenaDNA's H3 fraction (0.070) falls \\emph{below} chance,
meaning it systematically perturbs partners less than adjacent
positions. The gap between the two leaders and the rest is
qualitative.""",
        """RiNALMo reaches the highest perturbation specificity
(PS $= \\psRiNALMo{}$), with \\nullRiNALMo{} gate-passing families
exceeding the derangement null and a partner-is-max precision of
\\hthreeRiNALMo{}. ERNIE-RNA follows at PS $= \\psErnieRNA{}$, with
\\nullErnieRNA{} families exceeding the null and H3 precision
\\hthreeErnieRNA{}. Both stand more than an order of magnitude above
the third model and far above the $1/3$ chance baseline
($p \\ll 0.001$, binomial test).

Caduceus (14M, DNA, BiMamba SSM) ranks third at PS $= \\psCaduceus{}$,
with \\nullCaduceus{} gate-passing families exceeding the null and H3
precision \\hthreeCaduceus{}---above chance, and more than an order of
magnitude below both leaders. The remaining models cluster near zero
PS. SpliceBERT exceeds the null in \\nullSpliceBERT{} gate-passing
families at PS $= \\psSpliceBERT{}$. HyenaDNA's H3 fraction
(\\hthreeHyenaDNA{}) falls \\emph{below} chance: it perturbs partners
less than adjacent positions.""",
    ),
    (
        "Rung 3: the untrained control's PS and gate counts come from the run",
        """ERNIE-RNA with randomized weights---preserving the architectural
attention bias but destroying all learned parameters---produces
$\\text{PS} = 7.3 \\times 10^{-8}$, six orders of magnitude below
the trained model. Only 10 of 32 families pass the positive control
gate (versus 30/32 trained). Of those 10 gate-passing families,
9 exceed the derangement null---a seemingly high rate that reflects
the gate's role as a filter: families that pass the gate with random
weights happen to have stem/loop distributions where even random
perturbations produce above-chance stem sensitivity, and the
derangement null within those selected families is correspondingly
easy to exceed. The absolute PS magnitude ($10^{-7}$) confirms that
no meaningful partner specificity is present. Architecture without
training produces nothing.""",
        """ERNIE-RNA with randomized weights---preserving the architectural
attention bias and destroying every learned parameter---produces
$\\text{PS} = \\psUntrainedErnieRNA{}$, orders of magnitude below the
trained model. \\gateUntrainedErnieRNA{} families pass the positive
control gate, against \\gateErnieRNA{} trained. Of those,
\\nullUntrainedErnieRNA{} exceed the derangement null---a high rate
that reflects the gate's role as a filter: families passing the gate
with random weights have stem/loop distributions where even random
perturbations produce above-chance stem sensitivity, and the
derangement null within those families is correspondingly easy to
exceed. The absolute PS magnitude confirms that no partner
specificity is present.""",
    ),
    (
        "Rung 3: the bias ablation names the panel it ran on",
        """full Rung~1 pipeline yields a mean ratio of 1.235 (versus 1.242
trained), 32/49 families exceeding the nucleotide null (versus 31),
30/49 surviving the dinucleotide null (identical), attention contact
$\\rho = 0.101$ (versus 0.100), and probing accuracy 0.666
(identical). The ablated model matches or marginally exceeds the trained model on
every metric. ERNIE-RNA's structure awareness resides entirely in learned
weights; the architectural bias provides capacity but is functionally
redundant given sufficient pretraining. This result---combined with the
untrained control showing that the bias alone produces nothing---indicates
that the bias accelerates or stabilizes learning of pairing patterns
during pretraining, but the learned representations carry the signal
at inference time.""",
        """full Rung~1 pipeline on the pre-repair panel of 49 families (the
ablation is not re-run, so its numbers are not comparable cell for
cell with Table~\\ref{tab:rung1}) yields a mean ratio of 1.235 against
1.242 for the trained model, 32 of 49 families exceeding the
nucleotide null against 31, 30 of 49 surviving the dinucleotide null
in both, attention contact $\\rho = 0.101$ against 0.100, and probing
accuracy 0.666 in both. Within that panel the ablated model matches
the trained one on every metric. ERNIE-RNA's structure awareness
resides in learned weights, and the architectural bias provides
capacity that pretraining makes redundant. Read alongside the
untrained control, where the bias alone produces nothing, the bias
accelerates or stabilizes learning of pairing patterns during
pretraining while the learned representations carry the signal at
inference time.""",
    ),
    (
        "Rung 3: the PS peak layers come from the run",
        """RiNALMo's PS peaks at layer~33 (of 34) in 30 of 31 gate-passing
families. ERNIE-RNA peaks at layers~10--12 (of 12). The weak-signal
models peak at layers~0--1, consistent with tokenization and
positional-encoding artifacts.""",
        """RiNALMo's PS peaks at layer~\\psPeakRiNALMo{} of \\psLayersRiNALMo{}
in \\psPeakCountRiNALMo{} gate-passing families. ERNIE-RNA peaks at
layer~\\psPeakErnieRNA{} of \\psLayersErnieRNA{} in
\\psPeakCountErnieRNA{}. The weak-signal models peak at layers~0--1,
consistent with tokenization and positional-encoding artifacts.""",
    ),
    (
        "Rung 3: the tokenization counts come from the run",
        """tokens); only 1 of 4 evaluable families passed the gate. DNABERT-2's
BPE tokenizer limits evaluation to 8 of 52 families; none pass. These
models' near-zero or negative PS reflects tokenization constraints as
much as architectural ones; their Rung~1 attention-contact signal
(above) shows they encode structure through a different channel.""",
        """tokens), and \\gateCountNTvTwo{} of its \\eligibleNTvTwo{} evaluable
families passes the gate. DNABERT-2's BPE tokenizer limits evaluation
to \\eligibleDNABERTTwo{} of \\panelN{} families, of which
\\gateCountDNABERTTwo{} pass. The near-zero and negative PS of these
models reflects tokenization as much as architecture; their Rung~1
attention-contact signal shows they encode structure through a
different channel.""",
    ),
    (
        "Rung 3: H2$_6$ is computed once and quoted from the same test the table uses",
        """This fails ($r_b = -0.28$, $p = 0.265$): Caduceus (14M, DNA) exceeds
three RNA models (SpliceBERT, RNA-FM, UTR-LM) on PS. The two leaders
happen to be RNA-pretrained, but pretraining domain alone does not
explain their success---ERNIE-RNA's bias accelerates learning at
smaller scale (86M) and RiNALMo's 650M parameters reach the same
capability without it.""",
        """It fails ($r_b = \\domainRankBiserial{}$, $p = \\domainP{}$): Caduceus
(14M, DNA) exceeds three RNA models---SpliceBERT, RNA-FM and
UTR-LM---on PS. The two leaders are RNA-pretrained, and pretraining
domain alone does not explain them: ERNIE-RNA's bias accelerates
learning at 86M parameters and RiNALMo reaches the same capability at
650M without one.""",
    ),
    (
        "Discussion: the attention paragraph quotes the same macros as the Results",
        """knowledge, but with an important caveat. NT~v2 produces the highest
attention-contact correlation ($\\rho = 0.322$), yet this signal is
entirely architectural: untrained NT~v2 produces $\\rho = 0.328$,
and the per-family delta averages $-0.006$ with no systematic
direction (Table~\\ref{tab:attn_trained_untrained}). NT~v2's 6-mer
tokenization aggregates base-pairing partners into adjacent tokens,
producing structured attention from positional proximity alone.
DNABERT-2's content-only attention ($\\rho = 0.257$ trained vs $0.000$
untrained) is the sole model showing learned attention-contact
structure; the practical attention including ALiBi position biases
would produce a lower correlation, so $0.257$ is an upper bound.""",
        """knowledge, with one caveat. NT~v2 produces the highest
attention-contact correlation ($\\rho = \\attnNTvTwo{}$), and the signal
is architectural: untrained NT~v2 produces
$\\rho = \\attnUntrainedNTvTwo{}$, and the per-family delta averages
$\\attnFamilyDeltaNTvTwo{}$ with no systematic direction
(Table~\\ref{tab:attn_trained_untrained}). NT~v2's 6-mer tokenization
aggregates base-pairing partners into adjacent tokens, producing
structured attention from positional proximity alone. DNABERT-2's
content-only attention ($\\rho = \\attnDNABERTTwo{}$ trained against
$\\attnUntrainedDNABERTTwo{}$ untrained) is the only learned
attention-contact structure in the panel; the practical attention
including ALiBi position biases would produce a lower correlation, so
\\attnDNABERTTwo{} is an upper bound.""",
    ),
    (
        "Availability: what is deposited is every curated record, not only the "
        "analyzed ones",
        "scripts, per-family result files, 52 Rfam structures, and",
        """scripts, per-family result files, all \\panelCurated{} curated Rfam
records (\\panelN{} analyzed, \\panelWithdrawn{} withdrawn, each with its
exclusion reason), and""",
    ),
]


def macros_used() -> set[str]:
    r"""Every \foo{} this patch writes into the manuscript.

    Only the ones the patch introduces: v12 uses no generated macro, so an
    undefined one can only come from here. LaTeX's failure for a missing macro
    is `Undefined control sequence` against a line number in a file nobody
    edited by hand, which is a bad way to find out.
    """
    return {name for _, _, new in EDITS
            for name in re.findall(r"\\([a-zA-Z]+)\{\}", new)}


def unbraced_macros() -> set[str]:
    r"""Generated macros this patch writes without their `{}`.

    `macros_used` finds `\foo{}`, so `$\psRiNALMo$` -- valid LaTeX, and the
    natural thing to write inside math -- is invisible to the check that every
    macro is defined. A generated name is camelCase and no LaTeX primitive this
    patch writes is, which is enough to tell them apart.
    """
    return {name for _, _, new in EDITS
            for name in re.findall(r"\\([a-z][a-zA-Z]*[A-Z][a-zA-Z]*)(?![a-zA-Z]|\{\})",
                                   new)}


def macros_defined() -> set[str]:
    """Every macro the two generators wrote into paper/generated/."""
    return {name for filename in REQUIRED_GENERATED
            for name in re.findall(r"\\newcommand\{\\([a-zA-Z]+)\}",
                                   (GENERATED / filename).read_text())}


def build(text: str, verbose: bool = False) -> str:
    """paper_v12.tex in, paper_v13's text out.

    Separate from `main` so the edits can be checked against the real v12
    without a completed re-run on disk: every anchor here is a verbatim string
    from a document that four revisions have already moved things around in,
    and an anchor that stopped matching is not visible until the day the
    results land.
    """
    loose = unbraced_macros()
    assert not loose, ("a generated macro written without its braces is invisible "
                       "to the definedness check: " + ", ".join(sorted(loose)))

    for label, old, new in EDITS:
        count = text.count(old)
        if count != 1:
            raise SystemExit(f"ABORT -- [{count} matches] {label}")
        text = text.replace(old, new, 1)
        if verbose:
            print(f"  applied: {label}")

    for tex_label, generated in GENERATED_TABLES:
        text = replace_table(text, tex_label, generated)
        if verbose:
            print(f"  replaced: {tex_label} -> generated/{generated}.tex")

    # The hand-carried census is gone, in every class it named.
    for stale in ("tRNAs (7 families)", "rRNAs (4)", "riboswitches (8)",
                  "CRISPR repeats (3)", "spanning 8 major ncRNA classes"):
        assert stale not in text, f"the printed census survived: {stale}"

    # Every generated file is inputted exactly once, and no results table is
    # left typed into the body. A table swapped out but never inputted, or
    # inputted twice, both compile.
    for name in REQUIRED_GENERATED:
        stem = name[: -len(".tex")]
        found = text.count("\\input{generated/" + stem + "}")
        assert found == 1, f"generated/{stem}.tex is inputted {found} times"
    for tex_label, _ in GENERATED_TABLES:
        assert "\\label{" + tex_label + "}" not in text, \
            f"a v12 table survived the swap: {tex_label}"
        assert "\\ref{" + tex_label + "}" in text, \
            f"nothing refers to {tex_label} any more; the swap orphaned a table"

    assert "\\usepackage{longtable}" in text, "generated/panel_table.tex needs it"
    assert "\\panelN{}" in text
    # Both gates reach the reader, and neither replaces the other.
    assert "registered as 7 at $N = 32$" in text
    assert "the gate is 8 at the repaired $N$" in text
    # A paper reports the measurement, not the reading it replaced.
    for changelog in ("An earlier analysis", "reverses this interpretation",
                      "in a prior run", "Corrected from 0.230",
                      "H11 originally passed"):
        assert changelog not in text, f"changelog prose survived: {changelog}"

    return text


def preview() -> int:
    """The diff the patch would produce, without the run on disk.

    The generated tables and the stamps gate the real build, and neither exists
    until the containers finish; the prose edits are reviewable now.
    """
    v12 = SRC.read_text()
    diff = difflib.unified_diff(v12.splitlines(keepends=True),
                                build(v12).splitlines(keepends=True),
                                fromfile="paper_v12.tex",
                                tofile="paper_v13_repaired_panel.tex")
    print("".join(diff), end="")
    return 0


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

    undefined = sorted(macros_used() - macros_defined())
    if undefined:
        print("ABORT -- the patch writes macros no generator defines: "
              + ", ".join("\\" + name for name in undefined))
        return 1

    text = build(SRC.read_text(), verbose=True)
    DST.write_text(text)
    print(f"\n{len(EDITS)} edits and {len(GENERATED_TABLES)} table replacements "
          f"applied.\nwrote {DST.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diff", action="store_true",
                        help="print the diff against v12 and write nothing")
    raise SystemExit(preview() if parser.parse_args().diff else main())
