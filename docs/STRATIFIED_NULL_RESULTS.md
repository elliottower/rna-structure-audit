# Structure Metrics Results — Nucleotide-Stratified Null

**Date:** 2026-07-12
**Prereg SHA:** `694b43b` (causal-rna), `a7d10f5` (factorization-unified)
**Run:** Modal A10G, all 12 RNA families, all positions mutated

## Summary

The nucleotide-stratified composition control collapsed the mutation sensitivity signal from ~5.5x (confounded) to ~1.8x (composition-controlled). The composition confound was the dominant contributor to the original signal. The residual structure-dependent signal is real but modest, and falls below the preregistered 2.0 trained/untrained threshold.

This is the correct result. The composition control did exactly what it was designed to do.

## Mutation Sensitivity — Per-Family

| Family | RNA-FM trained | RNA-FM untrained | NT v2 | HyenaDNA | Caduceus | Evo |
|--------|---------------|-----------------|-------|----------|----------|-----|
| tRNA_Phe_yeast | 2.704 (✗ null) | 1.050 | 1.180 (✓) | 1.185 | 1.374 (✓) | 0.927 |
| tRNA_Ala_human | 3.210 (✓ null) | 1.091 | 1.207 (✓) | 1.284 | 1.328 | 0.899 |
| 5S_rRNA_ecoli | 1.157 | 1.008 | 0.963 | 1.012 | 1.075 | 1.048 |
| hammerhead_ribozyme | 1.329 | 1.076 | 1.118 | 1.292 | 1.456 | 1.147 |
| SAM_riboswitch | 2.515 (✓ null) | 1.056 | 1.042 (✓) | 1.074 | 1.269 | 1.283 |
| TPP_riboswitch | 1.441 | 1.093 | 1.058 | 1.227 | 1.324 | 1.888 |
| SRP_RNA_helix8 | 1.350 | 1.035 | 1.168 (✓) | 1.029 | 1.136 | 0.925 |
| mir_21_precursor | 1.312 | 1.049 | 1.086 | 1.868 | 1.341 | 1.156 |
| IRES_HCV_domainII | 1.804 | 1.045 | 0.928 | 1.190 | 1.225 | 0.994 |
| HDV_ribozyme | 1.318 | 1.021 | 1.175 (✓) | 1.150 | 1.291 (✓) | 1.316 |
| RNaseP_specificity | 1.510 | 1.034 | 0.998 | 1.217 | 1.086 | 1.217 |
| U2_snRNA_stem | 2.200 | 1.024 | 1.203 (✓) | 1.144 | 1.173 (✓) | 1.234 |
| **Mean** | **1.821** | **1.048** | **1.094** | **1.223** | **1.257** | **1.170** |
| **Median** | **1.475** | **1.047** | **1.102** | **1.188** | **1.280** | **1.151** |
| **Families > null** | **2/12** | **2/12** | **6/12** | **0/12** | **3/12** | **0/12** |

(✓ null) = real ratio exceeds 95th percentile of nucleotide-stratified null for that family. (✗ null) = within null but close.

## Hypothesis Results

### H1: Trained/untrained ratio ≥ 2.0 — FAIL

- Trained mean: 1.821
- Untrained mean: 1.048
- Ratio: **1.74x** (threshold: ≥ 2.0)
- Wilcoxon (trained > untrained): p = 0.000244, **12/12 families**

The preregistered decision criterion was "trained/untrained ratio ≥ 2.0" with "any trained/untrained ratio < 2.0 is a failure, regardless of the absolute value." 1.74x < 2.0. **H1 fails by the preregistered criterion.**

**Exploratory observation (not confirmatory):** Trained > untrained in 12/12 families (Wilcoxon p = 0.000244). This consistent direction across 12 independent families suggests a genuine learned effect, but the per-family unanimity was not named as a decision criterion in the prereg, so this is exploratory evidence.

Per-family trained/untrained ratios range from 1.15x (5S rRNA) to 2.94x (tRNA Ala). The tRNAs and SAM riboswitch show the strongest trained advantage (2.4–2.9x), while 5S rRNA and hammerhead show minimal advantage (1.15–1.24x).

### H1b: Exceeds nucleotide-stratified null — FAIL

- Families exceeding null: **2/12** (tRNA_Ala, SAM_riboswitch)
- Wilcoxon (real - null_95th > 0): p = 0.997
- Most families' null 95th percentiles are *higher* than the real ratios

This is the key finding. The stratified null is high because within-nucleotide variation in stem vs loop positions is substantial. Only 2/12 families show real ratios above the composition-matched null. The signal that survives composition matching is not significant at the population level.

### H2: Layer localization — FAIL (marginal)

- Median best layer: 6 (of 12)
- Wilcoxon (best_layer - 3 > 0): p = 0.036 (does not clear Bonferroni 0.0083)
- Best layers: 1, 1, 1, 1, 4, 5, 6, 7, 7, 11, 11, 12

Four families peak at layer 1 (near embedding), four at layers 4–7 (mid-network), four at layers 11–12 (late). The distribution is bimodal, not clearly shifted to deep layers.

### H3: NT v2 ratio < 2.0 — PASS (trivially)

- NT v2 mean ratio: 1.094
- Well below 2.0 as predicted

### H4: HyenaDNA ratio < 2.0 — PASS (trivially)

- HyenaDNA mean ratio: 1.223
- Well below 2.0 as predicted

### H5: Probing > baseline — FAIL

| Model | Best Accuracy | Baseline | Δ | Best Layer |
|-------|--------------|----------|---|------------|
| RNA-FM | 0.600 | 0.587 | +0.013 | 0 |
| NT v2 | 0.603 | 0.587 | +0.016 | 5 |
| HyenaDNA | 0.621 | 0.587 | +0.034 | 0 |
| Caduceus | 0.628 | 0.587 | +0.041 | 1 |
| Evo | 0.615 | 0.587 | +0.028 | 0 |

All five models barely exceed the 58.7% majority baseline. The GroupKFold (by family) makes this a cross-family generalization test — the probe must predict stem/loop on RNA types it hasn't seen. With only 12 families and 5-fold CV (2–3 families per fold), variance is high and the margin is not significant.

RNA-FM, HyenaDNA, and Evo peak at layer 0 (embedding), Caduceus at layer 1. Four of five models peak at or near the embedding, which means the probe is using input features, not learned representations. Caduceus achieves the highest probing accuracy (0.628) among all models, but the margin over baseline (+0.041) is still within noise given N=12 families. Evo (7B) at 0.615 does not outperform the 14M Caduceus or the 5.4M HyenaDNA despite being 500–1300x larger.

## Attention-Contact Correlation

| Model | Mean Best Corr | Max Best Corr |
|-------|---------------|---------------|
| RNA-FM trained | 0.067 | 0.099 |
| RNA-FM untrained | 0.061 | 0.098 |
| NT v2 trained | 0.278 | 0.509 |
| NT v2 untrained | 0.000 | 0.000 |

NT v2 trained attention correlates moderately with base-pairing contacts (mean 0.278, max 0.509 on hammerhead). NT v2 untrained attention is exactly zero — randomized weights produce uniform attention with no contact correlation. The 0.278 signal is entirely learned, not architectural. This is the strongest positive result in the experiment.

RNA-FM trained attention (0.067) is indistinguishable from RNA-FM untrained (0.061). RNA-FM's attention encodes zero learned structure.

**NT v2 per-family attention correlation:**

| Family | Best Corr | Layer | Head |
|--------|----------|-------|------|
| hammerhead_ribozyme | 0.509 | 4 | 11 |
| HDV_ribozyme | 0.315 | 9 | 8 |
| U2_snRNA_stem | 0.313 | 3 | 2 |
| SRP_RNA_helix8 | 0.309 | 2 | 5 |
| tRNA_Phe_yeast | 0.284 | 9 | 3 |
| tRNA_Ala_human | 0.284 | 3 | 10 |
| TPP_riboswitch | 0.285 | 4 | 12 |
| IRES_HCV_domainII | 0.265 | 11 | 12 |
| mir_21_precursor | 0.218 | 9 | 1 |
| SAM_riboswitch | 0.194 | 7 | 4 |
| RNaseP_specificity | 0.190 | 1 | 10 |
| 5S_rRNA_ecoli | 0.168 | 4 | 1 |

The hammerhead ribozyme (simple stem-loop) has the highest attention-contact correlation. Larger, more complex structures (5S rRNA, RNaseP) have the lowest. This pattern makes sense — simple stem-loops produce clear contact matrices; multi-branch junctions spread base-pairing across many positions.

**RNA-FM per-family:** All families show best correlation at layer 0, head 4, with values between 0.048 and 0.099. This one head captures a weak positional proximity signal, not structure-specific attention.

## What This Means

### The composition confound was real and dominant

The old (unfixed) mutation sensitivity showed RNA-FM at ~5.5x. The stratified null reveals that ~3.7x of that was GC enrichment in stems — complement swaps at G/C positions produce larger embedding perturbations than at A/U positions, and stems are 68.5% GC vs 43.9% in loops. The composition control eliminated this inflated signal, leaving a residual 1.82x that represents genuine (but modest) structure sensitivity.

### H1 fails. The per-family direction is exploratory.

1.74x < 2.0. The preregistered criterion is not met. Trained > untrained in 12/12 families (p < 0.001) is suggestive of a genuine learned effect, but this was not the preregistered statistic. Reporting it as confirmatory evidence would be a post-hoc endpoint switch. The honest conclusion: by the preregistered criterion, RNA-FM's structure sensitivity does not substantially exceed untrained weights. The per-family unanimity is noted as exploratory.

### Probing is underpowered at N=12

With 12 families in GroupKFold (2–3 families per holdout fold), the probe can't reliably learn cross-family generalization. The 1.3–3.4 percentage point margins above baseline are consistent with noise. Expanding to 50+ Rfam families would either confirm or kill the probing signal.

### NT v2 attention: validated learned structure signal

NT v2 trained attention correlates with base-pairing contacts at mean 0.278, max 0.509 (hammerhead ribozyme). NT v2 untrained attention is exactly 0.000 — randomized weights produce uniform attention that correlates with nothing. The entire 0.278 signal is learned during pretraining.

This is the strongest positive result in the experiment: NT v2's 6-mer tokenization and ESM+GLU architecture learn to attend to base-pairing partners, and this signal survives the most conservative comparison (trained vs random weights). RNA-FM, by contrast, shows trained attention (0.067) indistinguishable from untrained (0.061) — zero learned attention structure.

### The paper's contribution is the null, not the effect

The original d_g metric died because a 3-mer baseline matched it. This experiment designed a composition control that would have caught the same failure mode in mutation sensitivity — and it did catch a substantial composition confound. The paper's value is the methodology: preregistered hypotheses, measured composition enrichment (24.7% GC difference), nucleotide-stratified permutation, and an honest report that the confirmatory hypotheses failed.

### Two families surviving the null may still be composition-confounded

tRNA-Ala (3.21x) and SAM riboswitch (2.52x) are the only families that clear the stratified null, and both have high-GC stems with low-GC loops. The per-nucleotide stratification controls for global composition (same number of G's labeled "stem" in the null as in reality), but does not control for second-order effects: dinucleotide context (stem-G in GC pairs vs loop-G in GA contexts), positional proximity to other base-paired residues, or other features that distinguish stem-G from loop-G beyond nucleotide identity. A dinucleotide-stratified null would test whether this residual is genuine.

### Caduceus (BiMamba SSM): consistent with the pattern

Caduceus (14M params, BiMamba architecture) shows mean ratio 1.257, median 1.280 — between NT v2 (1.094) and HyenaDNA (1.223) on the low end and RNA-FM (1.821) on the high end. Only 3/12 families exceed the stratified null, and the margins are razor-thin (tRNA-Phe: 1.374 vs 1.321 null, HDV: 1.291 vs 1.290, U2: 1.173 vs 1.151). Signal concentrates at layer 13 (81% depth). No attention metric (SSM architecture). Probing: 0.628, the highest of any model but still only +0.041 above baseline.

As an SSM without attention, Caduceus provides a clean test of whether the structure signal appears in recurrent state-space representations. The answer: it does, at roughly the same magnitude as the attention-based models, and it similarly fails the composition-stratified null for most families.

### Evo (StripedHyena 7B): scale does not buy structure sensitivity

Evo (7B params, StripedHyena hybrid Hyena+attention) shows the weakest mutation sensitivity of all trained models: mean 1.170, median 1.151, 0/12 families exceed the stratified null. Two families (tRNA-Phe: 0.927, tRNA-Ala: 0.899) have ratios below 1.0 — Evo is *less* sensitive to stem mutations than loop mutations for those structures. TPP riboswitch (1.888) is the only family approaching the other models' range, but it still falls below null.

Despite being 500x larger than Caduceus (14M) and 70x larger than RNA-FM (99M), Evo shows weaker structure sensitivity than both. Probing accuracy (0.615) peaks at layer 0 (embedding) and is lower than Caduceus (0.628) or HyenaDNA (0.621). No attention metric is available because StripedHyena's hybrid architecture does not reliably expose attention weights.

Evo was pretrained on DNA (not RNA) with byte-level tokenization, which likely explains the weak RNA structure signal. The model's 7B parameters are spent on DNA-relevant features (regulatory motifs, coding patterns) rather than RNA secondary structure. This result reinforces the finding that structure sensitivity depends on pretraining data domain, not model scale.

### Next steps

1. ~~Untrained NT v2 attention baseline~~ — **DONE**: untrained = 0.000, trained 0.278 is fully learned
2. ~~Evo~~ — **DONE**: mean 1.170, 0/12 exceed null, weakest of all trained models
3. **Dinucleotide-stratified null** on tRNA-Ala and SAM — robustness check for second-order composition confound
4. **Expand to 50+ Rfam families** — probing and mutation sensitivity are underpowered at N=12
