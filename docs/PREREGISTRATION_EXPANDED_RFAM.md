# Preregistration: Expanded RNA Structure Awareness Evaluation

**Date**: 2026-07-12
**Authors**: Elliot Tower
**Repos**: causal-rna, factorization-unified
**Prior prereg SHA**: `694b43b` (causal-rna), `a7d10f5` (factorization-unified)

## Context

This preregistration covers the second phase of the RNA structure awareness evaluation. Phase 1 (SHA `694b43b`) tested 12 hand-curated RNA families across 5 models. Key results from Phase 1:

- H1 (trained/untrained ratio >= 2.0): FAIL at 1.74x
- H1b (>= 6/12 families exceed nucleotide-stratified null): FAIL at 2/12
- H5 (probing > baseline): FAIL (margins +0.013 to +0.041, within noise at N=12)
- NT v2 attention-contact correlation: 0.278 trained vs 0.000 untrained (observed, not preregistered)
- Exploratory: trained > untrained in 12/12 families (p=0.000244, not preregistered)

Phase 2 addresses three limitations identified in Phase 1:
1. **Statistical power**: N=12 families was underpowered for probing and cross-family tests
2. **Second-order composition confound**: nucleotide-stratified null does not control dinucleotide context
3. **Selection bias in null**: Phase 1 null was computed per-layer; the metric used best-of-layers, inflating the real ratio relative to the null

## Experimental design

### Dataset

- **Original 12 families** from Phase 1 (unchanged)
- **38+ additional Rfam families** spanning ribozymes, riboswitches, snRNAs, cis-regulatory elements, miRNA precursors, CRISPR repeats, and other ncRNAs
- Total: N >= 50 families
- All families have experimentally supported secondary structures (Rfam consensus or PDB crystal)
- Minimum 5 stem and 5 loop positions per family
- Sequence length 30-200 nt

### Models

Same 5 models as Phase 1:
- RNA-FM (99M, BERT, character, RNA-pretrained)
- NT v2 (56M, ESM+GLU, 6-mer, DNA-pretrained)
- HyenaDNA (5.4M, Hyena SSM, character, DNA-pretrained)
- Caduceus (14M, BiMamba, character, DNA-pretrained)
- Evo (7B, StripedHyena, byte-level, DNA-pretrained)

### Metrics

1. **Mutation sensitivity ratio** (complement swap, cosine distance at mutated position)
2. **Structure probing** (logistic regression, GroupKFold by family)
3. **Attention-contact correlation** (Spearman correlation, trained vs randomized weights)

### Null models

1. **Nucleotide-stratified permutation null** (same as Phase 1): shuffle stem/loop labels within each nucleotide type (A, U, G, C). 100 permutations. **Corrected**: null ratio computed as max-over-layers per permutation, matching the real-data selection procedure.

2. **Dinucleotide-stratified permutation null** (new): shuffle stem/loop labels within each dinucleotide type (AA, AC, AG, AU, CA, CC, ...). 100 permutations. Max-over-layers per permutation. Run on all families that exceed the nucleotide-stratified null, plus tRNA-Ala and SAM riboswitch (Phase 1 survivors, regardless of Phase 2 nucleotide-null result). This avoids conditioning the test only on Phase 2 survivors.

### Threshold for null exceedance

95th percentile of the permutation distribution, consistent with Phase 1.

## Hypotheses

### H6: RNA-FM trained/untrained ratio at N >= 50

**Criterion**: RNA-FM mean trained/untrained mutation sensitivity ratio >= 2.0 across N >= 50 families.

**Direction commitment**: We expect H6 to FAIL. Phase 1 showed 1.74x at N=12; we do not expect the ratio to increase with more families. H6 is a powered replication of H1, testing whether the Phase 1 null result holds at scale.

**If H6 passes** (ratio >= 2.0): the Phase 1 failure was a power artifact, and RNA-FM does encode structure at the preregistered magnitude. This would be a meaningful revision of our Phase 1 conclusion.

**If H6 fails** (ratio < 2.0): the Phase 1 null replicates, confirming that the composition-controlled effect does not reach practical significance.

### H6b: Families exceeding nucleotide-stratified null at N >= 50

**Criterion**: RNA-FM exceeds the max-over-layers nucleotide-stratified null in >= ceil(0.17 * N) families.

**Threshold rationale**: 0.17 is the Phase 1 rate (2/12 = 16.7%). This tests whether the same proportion of families exceed the null in the expanded set — a replication threshold, not an escalation. For N=50, the threshold is >= 9 families.

**If H6b passes**: the Phase 1 pattern replicates — a consistent minority of families show structure sensitivity above composition baselines.

**If H6b fails** (fewer than ceil(0.17 * N) exceed null): the Phase 1 survivors were likely noise, and the max-over-layers correction absorbed the remaining signal.

### H7: Dinucleotide null on surviving families

**Criterion**: Among families that exceed the nucleotide-stratified null (plus tRNA-Ala and SAM riboswitch regardless), at least one family also exceeds the dinucleotide-stratified null.

**Scope**: The dinucleotide null is run on:
  (a) All families exceeding the first-order (nucleotide) null in the expanded run
  (b) tRNA-Ala and SAM riboswitch (Phase 1 survivors, included regardless of Phase 2 first-order result)

**If H7 passes**: at least one family's structure sensitivity survives both first-order and second-order composition controls, providing genuine evidence of learned structure encoding.

**If H7 fails** (zero families exceed the dinucleotide null): all apparent structure sensitivity is attributable to dinucleotide-level composition differences between stems and loops.

### H8: Probing at N >= 50

**Criterion**: Best probing accuracy > majority baseline by >= 0.02 (2 percentage points) for at least one model, with GroupKFold cross-validation across N >= 50 families.

**Threshold rationale**: 0.02 is the minimum margin we consider distinguishable from noise. Phase 1 margins ranged +0.013 to +0.041 at N=12 with high variance; N >= 50 should reduce variance enough to resolve whether the signal is real.

**If H8 passes**: at least one model's embeddings encode stem/loop information beyond what input features provide.

**If H8 fails**: probing does not detect cross-family-generalizable structure information in any model's embeddings.

### H10: Sign test — trained > untrained at N >= 50

**Criterion**: RNA-FM trained mutation sensitivity ratio > RNA-FM untrained ratio in >= 0.75 * N families (sign test, one-sided binomial p < 0.01).

**Rationale**: Phase 1 showed trained > untrained in 12/12 families (p=0.000244), reported as exploratory because the sign test was not preregistered. This hypothesis promotes the observation to confirmatory status. The 75% threshold (rather than 100%) allows for noise families while still requiring a clear directional effect.

**If H10 passes**: there is a confirmed, directionally consistent learned effect — RNA-FM representations are more structure-sensitive than untrained weights across a broad set of RNA families. The magnitude may be modest (H6 tests magnitude), but the direction is real.

**If H10 fails**: the Phase 1 unanimity (12/12) was a small-sample artifact.

### H11: NT v2 attention-contact correlation on expanded families

**Criterion**: NT v2 mean trained attention-contact correlation > 0.15 across N >= 50 families, AND mean untrained correlation < 0.05.

**Threshold rationale**: Phase 1 showed trained = 0.278, untrained = 0.000. The 0.15 trained threshold is set at roughly half the Phase 1 value, allowing for weaker families (complex structures, long-range contacts) to pull down the mean. The 0.05 untrained threshold confirms the baseline is near zero.

**If H11 passes**: NT v2 attention encodes base-pairing structure as a genuine, replicable, preregistered finding — the paper's headline positive result.

**If H11 fails** (trained < 0.15 or untrained > 0.05): the Phase 1 attention result was inflated by the 12-family selection or the untrained baseline is not as clean as Phase 1 suggested.

## Methodological correction (not a hypothesis)

The Phase 1 nucleotide-stratified null was computed per-layer, while the metric used best-of-layers. This selection mismatch gave the real data an advantage the null did not receive. Phase 2 corrects this: both the nucleotide and dinucleotide nulls compute the permuted ratio at every layer and take the maximum, matching the real-data selection procedure.

The corrected null will also be applied to the original 12 families as a reanalysis. This reanalysis is reported transparently as a methodological correction, not as a preregistered test (the data have been seen). The reanalysis result may change the number of families exceeding the null from Phase 1's reported 2/12.

## Analysis script

The analysis code is committed in the same SHA as this document:
- `factorization-unified/scripts/modal_expanded_rfam.py` — Modal deployment script
- All adapter code in `factorization-unified/lib/adapters/`

## Decision table

| Hypothesis | Criterion | Pass means | Fail means |
|-----------|-----------|------------|------------|
| H6 | ratio >= 2.0 at N>=50 | Phase 1 null was underpowered | Phase 1 null replicates |
| H6b | >= ceil(0.17*N) exceed null | Replication of Phase 1 rate | Phase 1 survivors were noise |
| H7 | >= 1 family exceeds dinuc null | Genuine structure signal | All signal is composition |
| H8 | probe > baseline + 0.02 | Embeddings encode structure | Probing is noise |
| H10 | trained > untrained in >= 0.75*N | Directional effect confirmed | 12/12 was artifact |
| H11 | NT v2 attn > 0.15, untrained < 0.05 | Attention learns base-pairing | Phase 1 was inflated |

## What is NOT preregistered

- The exact number of expanded families (target N >= 50, actual count depends on validation)
- Which specific new families are included (all that pass validation: >= 20 nt, >= 5 stem + 5 loop positions)
- Evo, HyenaDNA, Caduceus mutation sensitivity thresholds (only RNA-FM and NT v2 have specific hypotheses; the other models are reported descriptively)
- Any analysis of the per-family composition statistics (GC%, dinucleotide distribution) — these are computed for transparency but do not have decision criteria
