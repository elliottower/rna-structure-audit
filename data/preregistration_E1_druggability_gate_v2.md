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
pending (Lp(a)/pelacarsen, IL-6/ziltivekimab). The composition-null gate follows the same
methodological pattern — freeze the rule, declare the targets, score against labeled
outcomes — in a new domain (RNA foundation model target nomination).

## 2. Design (confirmatory)

- **Units:** 16 RNA targets x 7 models = 112 model-target evaluations.
- **Models:** RNA-FM (99M), RiNALMo (650M), ERNIE-RNA (86M), NTv2 (56M), SpliceBERT (19.4M),
  UTR-LM (1.2M), DNABERT-2 (117M). The first six match Paper C Phase 5; DNABERT-2 adds the
  only model with demonstrated learned attention signal (Paper C Phase 4).
- **Classes:** 8 positive (druggable, structure-validated), 5 negative (no conserved
  targetable structure / synthetic composition controls), 3 graded HTT repeat-length controls.
- **We train nothing.** All models are frozen, pretrained, third-party. This is an
  evaluation-only protocol by design; the deliverable is a decision gate, not a model.

## 3. Hypotheses (frozen before running)

- **H-E1 (primary):** AUC(gate vs druggability label) > AUC(probing accuracy vs label),
  with a >= 0.10 AUC margin, one-sided DeLong test at alpha = 0.05.
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
- **Second gate (Metric C):** Spearman correlation between attention and ViennaRNA base-pair
  matrix >= preregistered threshold; used only for H-E3.

## 5. Analysis plan

1. Compute A and B for all 112 pairs; freeze raw JSON.
2. ROC + AUC for A and B against the binary label; DeLong test for H-E1.
3. Spearman monotonicity on CTL-01..03 for H-E2.
4. Two-gate crosstab for H-E3.
5. Confusion matrices on negatives for H-E4.

## 6. Stopping rule and multiplicity

Fixed n = 112; no optional stopping. Four hypotheses; H-E1 is primary and Bonferroni-protected
(alpha = 0.05/1); H-E2..E4 are secondary/confirmatory at alpha = 0.05, reported with correction.

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
  the ROC analysis (not imputed).

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

## 10. Infrastructure

Parallelized across Modal containers (one per model-target block), using the same
infrastructure as Paper C Phase 5 (24 containers, ~45 min wall-clock for full sweep).
All Modal images include ViennaRNA 2.7.0 for structure prediction.

## 11. Freeze record

- Code SHA: __________ (to fill at freeze)
- Benchmark TSV SHA: __________
- ViennaRNA version: 2.7.0
- Preregistration timestamp: 2026-07-13
