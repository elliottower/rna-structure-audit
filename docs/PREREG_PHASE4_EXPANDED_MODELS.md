# Phase 4 Preregistration: Expanded Model Coverage

## Motivation

Phases 1-3 evaluated seven models across five architectures. The
attention-contact finding reversed: NT v2's correlation (rho = 0.322) is
architectural (untrained = 0.328), not learned. All four attention-capable
models show trained attention indistinguishable from untrained.

Three gaps remain:

1. RNA-pretrained BERT coverage is n=3 (RNA-FM, RiNALMo, UTR-LM) — a fourth
   strengthens the null attention claim across BERT-style RNA models.
2. No pre-mRNA model has been tested. SpliceBERT is trained on splice-site
   context, which includes stem-loop structures near splice junctions.
3. No BPE-tokenized DNA model has been tested. NT v2's architectural attention
   correlation comes from 6-mer tokenization; BPE tokenization may produce
   a different pattern.

## Timing attestation

This preregistration is frozen before any Phase 4 computation. No
attention weights, mutation sensitivity ratios, or probing accuracies
have been computed for ERNIE-RNA, SpliceBERT, or DNABERT-2 at the time
of this SHA commit. No Modal images for these models have been built.
No adapter code for these models exists yet. The freeze precedes all
data collection and analysis for Phase 4.

## Models

- **ERNIE-RNA** (86M params, 12-layer BERT, 12 heads, 768 hidden,
  character-level, RNA-pretrained on non-coding RNA).
  HuggingFace: `multimolecule/ernierna`.
  Package: `multimolecule` (same image as RiNALMo/UTR-LM).
  Rationale: Fourth RNA-pretrained BERT model. Same architecture family as
  RNA-FM (86M vs 99M), different pretraining corpus. Tests whether the null
  attention result is specific to RNA-FM's training data or general to
  RNA BERT models.

- **SpliceBERT** (19.4M params, 6-layer BERT, 16 heads, 512 hidden,
  character-level, pre-mRNA from 72 vertebrate species, max 1024 nt).
  HuggingFace: `multimolecule/splicebert`.
  Package: `multimolecule`.
  Rationale: Pre-mRNA pretraining includes intronic stem-loop structures
  near splice sites. Smallest transformer in the study (1.2M UTR-LM is
  smaller but UTR-focused). Tests whether splice-context pretraining
  encodes general secondary structure.

- **DNABERT-2** (117M params, BERT + ALiBi positional encoding, BPE
  tokenization with 4096 vocab, DNA multi-species).
  HuggingFace: `zhihan1996/DNABERT-2-117M`.
  Package: custom (MAGICS-LAB/DNABERT_2 GitHub).
  Rationale: BPE tokenization produces variable-length tokens, unlike
  NT v2's fixed 6-mer. ALiBi attention has no learned positional
  embeddings. Tests whether the architectural attention-contact
  correlation is specific to 6-mer tokenization or general to DNA
  transformers.

## Hypotheses

All hypotheses use N=52 Rfam families. Each specifies: metric, threshold,
test, and pass/fail interpretation.

**H16** (ERNIE-RNA attention null): ERNIE-RNA trained attention-contact
mean Spearman rho < 0.10 across 52 families.

- Metric: mean best-head-best-layer Spearman rho between symmetrized
  attention and binary contact map, across 52 families.
- Threshold: 0.10 (same as H12).
- Test: one-sample comparison of mean against threshold.
- PASS (rho < 0.10): The null attention result holds across 4/4
  RNA-pretrained BERT models.
- FAIL (rho >= 0.10): ERNIE-RNA's pretraining corpus produces attention
  structure that RNA-FM's does not.

**H17** (ERNIE-RNA mutation sign test): ERNIE-RNA trained mutation
sensitivity exceeds untrained in >= 75% of 52 families.

- Metric: per-family indicator (trained ratio > untrained ratio).
- Threshold: 39/52 (75%).
- Test: binomial sign test.
- PASS: Directional structure effect confirmed for ERNIE-RNA.
- FAIL: ERNIE-RNA shows no consistent directional advantage over
  untrained baseline.

**H18** (SpliceBERT attention null): SpliceBERT trained attention-contact
mean Spearman rho < 0.10 across 52 families.

- Metric: same as H16.
- Threshold: 0.10.
- PASS: Splice-context pretraining does not produce structure-aware
  attention.
- FAIL: Splice-site pretraining produces attention structure that
  general RNA/DNA pretraining does not.

**H19** (SpliceBERT mutation sensitivity): SpliceBERT mean
trained/untrained mutation ratio < 2.0 across 52 families.

- Metric: mean across 52 families of (trained stem/loop cosine distance
  ratio) / (untrained stem/loop cosine distance ratio).
- Threshold: 2.0 (same as H6/H14 — uniform across all models).
- Justification for keeping the uniform threshold: applying a different
  threshold per model based on expected performance would introduce
  post-hoc flexibility. The 2.0 threshold was chosen in Phase 1 to
  represent a practically significant effect size regardless of model
  size. SpliceBERT is expected to fall well below 2.0; reporting this
  as a PASS is informative, not a failure of the experiment.
- PASS (ratio < 2.0): SpliceBERT does not show practically significant
  structure sensitivity.
- FAIL (ratio >= 2.0): SpliceBERT's splice-context pretraining produces
  strong structure sensitivity, contradicting the null pattern.

**H20** (DNABERT-2 attention, trained-vs-untrained contrast):
|DNABERT-2 trained attention rho - DNABERT-2 untrained attention rho|
< 0.05 across 52 families.

- Metric: absolute difference between mean trained and mean untrained
  attention-contact Spearman rho across 52 families.
- Threshold: 0.05 (same margin as Phases 2-3 null models: RNA-FM
  |0.051 - 0.060| = 0.009, UTR-LM |0.052 - 0.053| = 0.001, NT v2
  |0.322 - 0.328| = 0.006 — all well within 0.05).
- Test: one-sample comparison of |Delta| against 0.05.
- PASS (|Delta| < 0.05): DNABERT-2's attention-contact correlation
  is present at initialization, confirming the architectural pattern
  extends to BPE-tokenized DNA transformers.
- FAIL (|Delta| >= 0.05): DNABERT-2 shows learned attention structure,
  distinguishing it from NT v2's architectural pattern.

**H21** (DNABERT-2 attention magnitude): DNABERT-2 trained
attention-contact mean Spearman rho < 0.20 across 52 families.

- Metric: mean best-head-best-layer Spearman rho (same as H16/H18).
- Threshold: 0.20 (midpoint between the BERT-style null cluster at
  0.05-0.10 and NT v2's architectural 0.322 — tests whether BPE
  tokenization produces a weaker architectural baseline than 6-mer).
- PASS (rho < 0.20): BPE tokenization produces a weaker architectural
  attention-contact correlation than NT v2's 6-mer tokenization.
  The 6-mer aggregation mechanism specifically drives NT v2's high
  architectural correlation.
- FAIL (rho >= 0.20): The architectural correlation is a general
  property of DNA transformers, not 6-mer-specific.

## Methods

Identical to Phases 2-3:
- Mutation sensitivity ratio under complement swap (A<->U, C<->G)
- Nucleotide-stratified permutation null (100 permutations, max-over-layers)
- Attention-contact Spearman correlation (trained vs untrained)
- 52 Rfam families from Phase 2
- Per-position cosine distance at mutated position (not mean-pooled)
- attn_implementation="eager" for all models (SDPA fix from Phase 3)

Torch version pinned to match Phase 3 (2.13.0).

ERNIE-RNA and SpliceBERT use the multimolecule image (transformers 5.x).
DNABERT-2 uses a separate image with its custom codebase.

## Analysis code

Adapters to write (after SHA freeze):
- `lib/adapters/ernierna.py`
- `lib/adapters/splicebert.py`
- `lib/adapters/dnabert2.py`

Runner functions to add to `scripts/modal_expanded_rfam.py`:
- `run_ernierna_expanded`
- `run_splicebert_expanded`
- `run_dnabert2_expanded`

## Local analysis (no GPU, exploratory)

After Phase 4 GPU runs complete, produce the following from existing +
new data (all 10 models). These analyses are exploratory and labeled
as such in the paper.

1. **Mutation sensitivity forest plot**: Per-model mean mutation
   sensitivity ratio with 95% bootstrap CI (1000 resamples of the 52
   per-family ratios). One row per model, sorted by mean ratio.
   Horizontal reference line at 2.0 practical threshold. Effect size
   metric: the mutation sensitivity ratio itself (trained stem/loop
   cosine distance ratio).

2. **Attention trained-vs-untrained scatter**: Single panel showing
   trained rho (x-axis) vs untrained rho (y-axis) for all
   attention-capable models (up to 10). Identity line y=x; points
   near the line = architectural, points below = learned. One point
   per model (mean across 52 families).

3. **Size-sensitivity scatter (illustrative)**: log10(params) vs mean
   mutation ratio across all 10 models. Reported as illustrative with
   Spearman rank correlation and 95% bootstrap CI. No trend line
   fitted — n=10 is insufficient for a scaling law.

4. **RiNALMo attention heatmap**: Per-family Delta (trained - untrained)
   attention for RiNALMo's borderline result (+0.019 mean). Families
   sorted by structural class (tRNA, riboswitch, ribozyme, etc.).
   Tests whether the weak signal concentrates in specific motifs.

## Preregistration date

2026-07-13

## Status

FROZEN — awaiting git commit for SHA. No Phase 4 data exists at time
of freeze.
