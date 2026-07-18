# Phase 6 Results: Perturbation Specificity (Partial — 4 of 7 models)

**Date:** 2026-07-14
**Prereg:** PREREGISTRATION_PHASE6_V2.md (SHA `c19aa59`)
**Status:** 4 models complete (ernierna, splicebert, utrlm, hyenadna). Pending: rnafm, rinalmo, evo.

---

## Summary

Phase 6 tests whether RNA foundation models resolve specific base-pair partners, a stronger form of structure awareness than the stem-vs-loop discrimination tested in Phases 1--5. The metric is perturbation specificity (PS): after swapping one base to its Watson-Crick complement, does the model perturb the paired partner more than adjacent stem positions? The within-stem derangement null controls for positional and compositional confounds by permuting partner labels within each stem.

ERNIE-RNA shows unambiguous partner specificity. The remaining three completed models do not.

---

## Per-model results

### ERNIE-RNA (multimolecule, encoder, 86M)

| Metric | Value |
|--------|-------|
| Active families | 32 |
| Positive control gate pass | 30/32 |
| Mean PS (gate-pass) | 0.1170 |
| Median PS | 0.1388 |
| Max PS | 0.1973 (glycine riboswitch) |
| Families exceeding primary null | 28/30 |
| Families exceeding conservative null | 28/30 |
| H3 precision (partner-is-max fraction) | 151/171 = 0.883 |
| Peak layer range | 10--12 (median 12, final layer) |

ERNIE-RNA resolves specific base-pair partners in 28 of 30 gate-passing families. The effect is large: PS values range from 0.06 to 0.20, meaning the partner position is perturbed 6--20 percentage points more (in cosine distance) than the best adjacent stem position. Both the primary and conservative nulls agree on 28/30. The partner-is-max fraction of 0.883 far exceeds the 1/3 chance baseline, indicating single-nucleotide precision.

All peak layers fall between 10 and 12 (of 12 total), concentrating in the final layer. The signal is in the deepest representations, consistent with learned pairing rather than positional encoding.

Top 5 families by PS:

| Family | PS | Layer | Eligible pairs |
|--------|-----|-------|---------------|
| glycine riboswitch | 0.197 | 12 | 9 |
| purine riboswitch | 0.197 | 12 | 16 |
| U5 snRNA | 0.191 | 12 | 18 |
| mir-122 precursor | 0.191 | 12 | 12 |
| mir-155 precursor | 0.191 | 12 | 12 |

Gate-fail families: 7SK RNA, U1 snRNA (stem perturbation not significantly above loop perturbation).

Quarantined pilot families (excluded from confirmatory counts):
- tRNA Phe yeast: PS = 0.281 at layer 11 (would have been the highest value)
- tRNA Ala human: PS = 0.178 at layer 11

Neither quarantined family would have changed the outcome. H1 passes at 28/30 regardless.

### SpliceBERT (multimolecule, encoder, 19M)

| Metric | Value |
|--------|-------|
| Active families | 32 |
| Gate pass | 21/32 |
| Mean PS (gate-pass) | 0.000350 |
| Median PS | 0.000145 |
| Families exceeding primary null | 14/21 |
| H3 precision | 30/117 = 0.256 |
| Peak layer range | 0--6 (median 1) |

SpliceBERT shows marginal partner specificity. Mean PS is three orders of magnitude smaller than ERNIE-RNA. 14 families exceed null, which exceeds the H1 threshold of 7, but the effect size is negligible. The H3 partner-is-max fraction (0.256) is below the 1/3 chance baseline. Peak layers cluster at 0--1, suggesting the signal is in the tokenizer or positional encoding rather than learned representations.

11 families fail the positive control gate, the highest failure rate of the four models tested so far.

### UTR-LM (multimolecule, encoder, ~2M)

| Metric | Value |
|--------|-------|
| Active families | 32 |
| Gate pass | 16/32 |
| Mean PS (gate-pass) | 0.000014 |
| Median PS | 0.000003 |
| Families exceeding primary null | 13/16 |
| H3 precision | 28/84 = 0.333 |
| Peak layer range | 0--6 (median 2) |

UTR-LM shows no meaningful partner specificity. PS values are in the 10^-5 range. The H3 fraction of exactly 0.333 matches the chance baseline. Half the families fail the positive control gate. Despite 13/16 families exceeding null, the effect size is zero for practical purposes.

### HyenaDNA (SSM, 5.4M, DNA-pretrained)

| Metric | Value |
|--------|-------|
| Active families | 32 |
| Gate pass | 25/32 |
| Mean PS (gate-pass) | 0.000001 |
| Median PS | -0.000000 |
| Families exceeding primary null | 10/25 |
| H3 precision | 11/141 = 0.078 |
| Peak layer range | 0--2 (median 0) |

HyenaDNA shows no partner specificity. Mean PS is effectively zero (10^-6). The H3 fraction (0.078) is below chance, meaning the partner is systematically *less* perturbed than adjacent positions. Peak layers at 0 indicate the variation is in the input embedding, not learned computation. Only 10/25 exceed null, below the H1 threshold of 7 when scaled to gate-pass count.

---

## Preliminary hypothesis decisions (4 of 7 models)

### H1: At least one model has mean PS > 0 exceeding the derangement null

**Provisional: PASS (ERNIE-RNA)**

ERNIE-RNA satisfies both H1 conditions:
- Condition A: 28 families exceed null >= 7 threshold
- Condition B: Mean PS = 0.117 > 0 (Wilcoxon p will be significant given all 30 values are positive and large)

No other model has both conditions. SpliceBERT has 14 families exceeding null but near-zero mean PS; the Wilcoxon test will determine whether condition B holds.

### H2: RNA-pretrained PS > DNA-pretrained PS

**Cannot evaluate.** Requires all 10 models (5 RNA, 5 DNA). Currently have 3 RNA (ernierna, splicebert, utrlm) and 1 DNA (hyenadna). Pending: rnafm, rinalmo (RNA), evo (DNA).

### H3: Partner specificity at single-nucleotide precision

**Provisional: PASS (ERNIE-RNA)**

ERNIE-RNA: 151/171 = 0.883, far exceeding the 0.5 threshold and 1/3 baseline. Binomial test against 1/3 will be significant (p << 0.001).

No other model passes. SpliceBERT (0.256) and HyenaDNA (0.078) are below chance. UTR-LM (0.333) matches chance exactly.

---

## Interpretation

ERNIE-RNA is the only model among the four tested that encodes which specific position pairs with which. The effect is large, consistent across 28/30 families, localized to the final layers (10--12 of 12), and survives both primary and conservative derangement nulls.

The three-order-of-magnitude gap between ERNIE-RNA (PS ~ 0.12) and the next-best model (SpliceBERT, PS ~ 0.0004) suggests a qualitative difference rather than a graded spectrum. ERNIE-RNA appears to have internalized base-pairing topology in its representations; the others have not, at least not at a level detectable by single-nucleotide complement swaps.

The peak-layer distribution reinforces this interpretation. ERNIE-RNA peaks at layer 12 (the final layer), indicating the pairing information is computed through deep processing. SpliceBERT and UTR-LM peak near layer 0--1, and HyenaDNA at layer 0, consistent with positional or tokenization artifacts rather than learned structure.

### Open questions for remaining models

- **RNA-FM** (99M, BERT, RNA-pretrained): The only RNA-pretrained model with comparable architecture to ERNIE-RNA. Will it show similar partner specificity?
- **RiNALMo** (650M, encoder, RNA-pretrained): The largest RNA-pretrained model. Does scale help?
- **Evo** (7B, StripedHyena, DNA-pretrained): The largest model overall. Phases 1--5 showed Evo has the weakest structure sensitivity despite its size. Does that pattern hold for partner specificity?

---

## Data files

- `data/gpu_results/phase6_compensatory/ernierna_phase6_ps.json`
- `data/gpu_results/phase6_compensatory/splicebert_phase6_ps.json`
- `data/gpu_results/phase6_compensatory/utrlm_phase6_ps.json`
- `data/gpu_results/phase6_compensatory/hyenadna_phase6_ps.json`
