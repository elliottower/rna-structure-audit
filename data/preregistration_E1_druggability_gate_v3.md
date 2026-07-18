# Preregistration — Experiment E1: Composition-Null Gate as a Druggability Discriminator

**Project:** Composition-Controlled Evaluation of RNA Structure Awareness — Phase 6 (Prospective Validation)
**Author:** Elliot Tower
**Date frozen:** 2026-07-13
**Companion papers:** cross_architecture_v12 (Paper C), paper_d_v7 (Paper D), cross-design-evidence-discordance (Paper 3)
**Benchmark file:** druggability_benchmark_v2.tsv (SHA to be recorded at freeze)

---

## 1. Background and motivation

Papers C and D established a *negative* result: no RNA foundation model encodes secondary
structure beyond composition/architecture artifacts, and standard scIB metrics do not predict
transfer. Both are audits. This experiment tests the *positive, decision-relevant* claim that
motivates the TRDNT proposal: **the composition-null gate discriminates druggable RNA targets
from non-druggable ones better than the standard baseline (probing accuracy).**

This converts the work from "the gate says NO to bad models" into "the gate makes a better
go/no-go target-selection decision than current practice," using targets with independent
ground-truth druggability from R-BIND 2.0, HARIBOSS, and R-scape covariation analysis.

**Methodological precedent:** Paper 3 (Tower 2026, cross-design evidence discordance)
demonstrated that pre-registered evaluation gates predict real-world drug outcomes: a
two-criterion classification rule correctly classified 14/15 drug mechanism families across
three disease domains (neuro, cardio, autoimmune), with two live prospective predictions
pending (Lp(a)/pelacarsen, IL-6/ziltivekimab). Note: the 14/15 family tally (Paper 3
prereg) and the 195-pair balanced-accuracy analysis (Paper 7, pQTL MR validity) are
separate results from separate manuscripts. The composition-null gate follows the same
methodological pattern as Paper 3 — freeze the rule, declare the targets, score against
labeled outcomes — in a new domain (RNA foundation model target nomination).

## 2. Design (confirmatory)

- **Units:** 16 RNA targets x 7 models = 112 model-target evaluations.
- **Models:** RNA-FM (99M), RiNALMo (650M), ERNIE-RNA (86M), NTv2 (56M), SpliceBERT (19.4M),
  UTR-LM (1.2M), DNABERT-2 (117M). The first six match Paper C Phase 5; DNABERT-2 adds the
  only model with demonstrated learned attention signal (Paper C Phase 4).
- **Classes:** 8 positive (druggable, structure-validated), 5 negative (no conserved
  targetable structure / synthetic composition controls), 3 graded HTT repeat-length controls.
- **All foundation models are frozen, pretrained, third-party.** The only learned
  component is a lightweight linear probe used as the baseline metric (Metric A).
  This is an evaluation-only protocol by design; the deliverable is a decision gate,
  not a model.

## 3. Hypotheses (frozen before running)

- **H-E1 (primary):** AUC(gate vs druggability label) > AUC(probing accuracy vs label),
  with a >= 0.10 AUC margin. Test: Wilcoxon signed-rank test on 6 paired differences
  (one per gate-eligible model: gate AUC_i minus probing AUC_i), one-sided, alpha = 0.05.
  DNABERT-2 is excluded because BPE tokenization prevents mutation sensitivity (see
  Section 7). DeLong is underpowered at n = 13 binary targets; the paired-model design
  with 6 models yields a more appropriate unit of analysis. **Sensitivity note:** with
  6 pairs, the minimum detectable effect at 80% power (one-sided alpha = 0.05) requires
  all 6 differences to share the same sign — i.e., the gate must outperform probing on
  every gate-eligible model. If any model reverses direction, H-E1 fails. This is
  declared: the test is intentionally conservative for a proof-of-concept at this sample
  size.
- **H-E2:** On graded HTT controls, gate-passing models show Spearman rho >= 0.8 between
  embedding distance and repeat length; gate-failing models show |rho| < 0.3.
  (Predicted: RiNALMo passes/monotonic, RNA-FM fails/flat.)
- **H-E3 (two-gate necessity+sufficiency):** At least one target passes the composition null
  but fails the ViennaRNA base-pair-correspondence gate (degeneracy trap demonstrated).
- **H-E4 (retrospective screen):** On negatives, the gate's FAIL rate exceeds probing's FAIL
  rate (the gate flags hallucinated structure that probing accepts).

## 4. Metrics (exactly as in Paper C)

- **Baseline (Metric A):** logistic-regression probing balanced accuracy, GroupKFold by target.
- **Our gate (Metric B):** mutation-sensitivity ratio R under complement swap with
  nucleotide-stratified permutation null (200 permutations; 95th percentile threshold;
  max-over-layers null selection). Dinucleotide-stratified second-order null on all first-order
  survivors. PASS = exceeds both nulls.
- **Second gate (Metric C):** Spearman correlation between attention (best head, max over
  layers) and ViennaRNA base-pair contact matrix >= 0.10. Threshold derived from Paper C
  Phase 4 results: the Rfam 52-family mean attention-contact correlation across all models
  is 0.05–0.10, with DNABERT-2 the only model showing learned signal above this range.
  Setting the threshold at rho >= 0.10 marks the upper end of the null distribution observed
  in Paper C. Used only for H-E3.

## 5. Analysis plan

1. Compute A and B for all 112 pairs; freeze raw JSON.
2. ROC + AUC for A and B against the binary label; Wilcoxon signed-rank test on 6 paired
   model-level AUC differences for H-E1 (DNABERT-2 excluded from gate).
3. Spearman monotonicity on CTL-01..03 for H-E2.
4. Two-gate crosstab for H-E3.
5. Confusion matrices on negatives for H-E4.

## 6. Stopping rule and multiplicity

Fixed n = 112 evaluations (6 models x 13 targets for gate, 7 x 13 for probing; 7 x 3 for
HTT graded); no optional stopping. Four hypotheses; H-E1 is primary (Wilcoxon signed-rank,
alpha = 0.05); H-E2..E4 are secondary/confirmatory at alpha = 0.05, reported with Bonferroni
correction across H-E2..E4.

## 7. Known limitations (declared in advance)

- Sequence windows and dot-bracket structures in the benchmark are curated from source
  databases and MUST be re-folded with ViennaRNA and verified against reference FASTA before
  freeze. Rows with dot-bracket = '.' are placeholders pending ViennaRNA RNAfold.
- Negative class relies on absence-of-covariation (R-scape null) as ground truth, which is
  evidence of no conserved structure, not proof of non-druggability.
- n = 16 targets is small; result is a proof-of-concept discriminator, not a clinical claim.
- Passing the gate is necessary, not sufficient (this is why H-E3 exists).
- Positive targets come from R-BIND 2.0 and HARIBOSS; selection is based on availability of
  bioactivity data and solved structure, not on expected gate performance.
- NEG-04 and NEG-05 (synthetic controls) are generated at runtime with frozen random seeds
  to ensure reproducibility. The seed is recorded in the SHA-frozen code.
- DNABERT-2's BPE tokenization prevents per-position mutation sensitivity; it is scored on
  attention-contact and probing only. Missing mutation-sensitivity cells are excluded from
  the ROC analysis (not imputed). For H-E1, this means DNABERT-2 contributes a probing AUC
  but no gate AUC, so the Wilcoxon signed-rank test runs on 6 paired differences (not 7).
  The effective target-level ROC for the gate uses n = 13 targets x 6 models = 78 pairs
  (not 91).

## 8. Sequential gating pipeline (TRDNT framing)

This experiment validates one stage of a proposed sequential decision pipeline for
RNA-targeted therapeutic development:

1. **ML screening** (cost: pennies, time: minutes) — RNA foundation model nominates
   structural targets from sequence.
2. **Composition-null gate** (cost: pennies, time: minutes) — filters targets where
   model confidence is composition artifact. **This experiment validates this stage.**
3. **Mendelian randomization check** (cost: free, time: days) — uses existing GWAS
   summary statistics to verify causal link between mechanism and disease (validated
   independently in Paper 3: 14/15 correct, two live predictions pending).
4. **Screening / RCT** (cost: millions, time: years) — only for targets passing all gates.

Each gate is orders of magnitude cheaper than the next. The composition-null gate and MR
discordance gate have been independently validated; this experiment tests whether the
composition-null gate discriminates real druggable targets from undruggable ones.

## 9. Benchmark sequence verification

All sequences in druggability_benchmark_v2.tsv were verified against primary databases:

| Target | Source | Verification |
|--------|--------|-------------|
| POS-01 SMN2 | NM_017411.3 pos 998-1051 | NCBI Entrez |
| POS-02 FMN riboswitch | Rfam RF00050 | Already in Rfam pipeline |
| POS-03 HCV IRES | PDB 1P5O / Rfam RF00061 | Already in Rfam pipeline |
| POS-04 HIV-1 TAR | PDB 1ARJ | PDB crystal structure |
| POS-05 PreQ1 | Rfam RF00522 | Already in Rfam pipeline |
| POS-06 MALAT1 ENE | PDB 4PLX | PDB crystal structure |
| POS-07 EV71 IRES SLII | PDB 6XB7 | PDB NMR structure |
| POS-08 SARS-CoV-2 5'UTR | NC_045512.2 pos 56-101 | NCBI RefSeq genome |
| NEG-01 HOTAIR | NR_047517.1 pos 1-100 | NCBI Entrez efetch |
| NEG-02 SRA | NR_045587.1 pos 1-100 | NCBI Entrez efetch |
| NEG-03 Xist repeat A | NR_001564.2 pos 370-480 | NCBI Entrez efetch |
| NEG-04 Shuffled tRNA | Generated at runtime | uShuffle, seed frozen |
| NEG-05 Random GC-matched | Generated at runtime | Seed frozen |
| CTL-01/02/03 HTT | NM_002111.7 | Full constructs with flanks, verified in Paper C |

Perplexity-generated sequences from v1 were rejected: 13/16 were fabricated or wrong.
All v2 sequences pulled from NCBI Entrez API or existing verified pipeline.

**Pipeline structure provenance note:** The Paper C Rfam pipeline (`modal_expanded_rfam.py`)
had 39/52 families with hardcoded structure strings longer than their sequences (by 1–6
characters). These were corrected: 38 truncated from the 3' end to match sequence length,
1 padded (Histone_3prime, 1 character shorter). Truncation orphaned 5 base pairs across 5
families (7SK_RNA, CRISPR_leader, Bacterial_SRP, IRE_stem_loop, Corona_5UTR), which were
rebalanced by converting the unmatched bracket to `.`. This fix affects the Paper C 52-family
pipeline only. The three E1 benchmark targets sourced from the Rfam pipeline (POS-02 FMN,
POS-03 HCV IRES, POS-05 PreQ1) share *sequences* with the pipeline but use `.` placeholder
structures (folded at runtime by ViennaRNA), so the truncation fix does not alter any E1
benchmark data.

## 10. Infrastructure

Parallelized across Modal containers (one per model-target block), using the same
infrastructure as Paper C Phase 5 (24 containers, ~45 min wall-clock for full sweep).
All Modal images include ViennaRNA 2.7.0 for structure prediction.

**ViennaRNA fold parameters (pinned for reproducibility):** All folds use
`RNA.fold_compound(seq, md).mfe()` with explicitly constructed `RNA.md()` model details:
temperature = 37.0°C, dangles = 2, noLP = 0 (lonely pairs permitted), energy parameters =
Turner 2004 (ViennaRNA 2.7.0 default). Parameters are set per-call via `fold_compound`,
not via `RNA.fold()` which reads mutable global state. Identical sequence + identical
ViennaRNA version + identical md parameters = identical dot-bracket output. The SHA freeze
covers the code that constructs these parameters; the benchmark TSV's `.` placeholders are
replaced at runtime, not pre-computed.

**Permutation null draws are shared across models per target.** The RNG seed for the
composition-null permutations is `42 + rna_idx * 1000`, keyed on target index only (not
model). All 6 gate-eligible models receive identical label shuffles for each target,
ensuring the Wilcoxon signed-rank pairing compares model-level gate scores under the same
null draws.

**NEG-04/05 reproducibility.** Runtime-generated sequences (NEG-04: dinucleotide-shuffled
tRNA-Phe, NEG-05: random GC-matched) use `numpy.random.default_rng(42)` (PCG64 bit
generator, version-stable since NumPy 1.17). Modal images pin numpy==2.4.6 (numpy==1.26.4
for DNABERT-2). The frozen TSV contains `GENERATE_AT_RUNTIME` for these two targets;
their actual sequences are determined by the frozen code + frozen numpy version + seed=42.

## 11. Changes from v2

v3 addresses five issues raised in external review:

1. **H-E1 statistical test reframed (blocker).** DeLong test at n = 13 binary targets is
   underpowered. Replaced with Wilcoxon signed-rank test on paired model-level AUC
   differences (6 pairs, one-sided). Sensitivity note added: the test requires all 6
   gate-eligible models to show the same sign.

2. **Placeholder dot-brackets are by design (non-issue).** Structures are set to '.' because
   they must be folded at runtime with ViennaRNA on GPU. Section 7 already declares this.
   The SHA freeze covers the code that calls `fold_compound()` with pinned parameters.

3. **Paper 3 "14/15" claim verified (non-issue).** The external review cited Paper 7
   (`07_pqtl_mr_validity.pdf`) as its source, not Paper 3. Paper 3's frozen PREREGISTRATION.md
   explicitly reports 14/15 across neuro + cardio + autoimmune, with two live prospective
   predictions. No correction needed. Added disambiguation footnote in Section 1 to prevent
   future confusion between Paper 3 (14/15 families) and Paper 7 (195-pair BA analysis).

4. **H-E3 threshold pre-declared (fix).** Metric C now specifies rho >= 0.10, derived from
   Paper C Phase 4 Rfam 52-family attention-contact correlation distribution (mean 0.05–0.10
   across all models).

5. **DNABERT-2 effective denominator stated (fix).** Section 7 now explicitly states that
   DNABERT-2 contributes probing but not gate scores, yielding 6 (not 7) Wilcoxon pairs and
   78 (not 91) target-level gate evaluations.

6. **ViennaRNA fold enforced in code, not just documented (v3 blocker fix).** Replaced
   `RNA.fold()` (reads mutable global state) with `RNA.fold_compound(seq, md).mfe()` using
   explicit `RNA.md()` per call. Section 10 documents the parameters; the code now enforces
   them. Without this fix, the provenance string would have been false documentation.

10. **NEG-04 shuffle corrected from mononucleotide to dinucleotide-preserved (v3 blocker
    fix).** The benchmark TSV and Section 7 claimed dinucleotide-preserved shuffling (uShuffle),
    but the code used `rng.shuffle(nucs)` — a plain mononucleotide permutation. Replaced with
    an Altschul-Erickson Euler path implementation (`_dinucleotide_shuffle`) that preserves both
    mono- and dinucleotide composition. Verified: output length = input length, mononucleotide
    counts match, dinucleotide counts match, sequence differs from input. The benchmark TSV
    source column ("uShuffle") is updated to "Altschul-Erickson Euler path (equivalent to
    uShuffle)" to match the actual implementation. Without this fix, NEG-04 would have been a
    weaker null than documented, since mononucleotide shuffle destroys dinucleotide composition
    that foundation models may encode.

7. **Pipeline truncation provenance documented (v3 addition).** Section 9 now records that
   39/52 Paper C family structures were length-corrected (38 truncated, 1 padded, 5
   rebalanced). Verified that this fix does not affect any E1 benchmark data — shared
   targets use `.` placeholder structures folded at runtime.

8. **Shared permutation draws documented (v3 addition).** Section 10 now states that
   composition-null permutation seeds are keyed on target index only, ensuring identical
   null draws across models for valid Wilcoxon pairing.

9. **NEG-04/05 NumPy version pinned (v3 addition).** Section 10 now documents that runtime
   sequence generation depends on PCG64 (version-stable since NumPy 1.17) and that Modal
   images pin numpy==2.4.6 / 1.26.4.

## 12. Freeze record

- Code SHA: __________ (to fill at freeze)
- Benchmark TSV SHA: __________
- ViennaRNA version: 2.7.0
- Preregistration timestamp: 2026-07-13
