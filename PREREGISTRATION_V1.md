# Pre-registration: Direction Instability and Cross-View Validation of RNA-FM on HTT mRNA

**Title:** Direction instability of RNA foundation model embeddings across CAG repeat variants, with thermodynamic cross-view validation

**Authors:** Elliot Tower

**Date:** 2026-07-11 (pre SHA freeze)

**Commit SHA:** _TO BE FILLED AFTER COMMIT_

---

## Overview

This document pre-registers the hypotheses, decision criteria, and analysis
pipeline for two analyses extending the causal geometry audit of RNA-FM on
HTT mRNA. The initial audit (bracket norm, Grassmannian transportability,
weight SVD) has already been conducted and is reported separately. This
pre-registration covers experiments not yet run: (1) direction instability
of RNA-FM embeddings across CAG repeat lengths, and (2) cross-view validation
using EternaFold thermodynamic structure prediction.

## Prior results (already observed, NOT pre-registered)

The following results from the initial audit are known and motivate this
pre-registration:

- Grassmannian transportability collapses categorically: geodesic distance
  ~3.7 for ANY non-WT repeat count (including normal-range r=27)
- Layer-wise effective dimensionality: 51 → 1 across layers 0–12
- W_Q effective rank: 265 → 86 across layers (3x compression)
- W_Q–W_K subspace overlap: 0.81 in layer 3
- Per-position embedding norms are constant (LayerNorm effect)

These results are NOT hypotheses — they are facts that inform the design below.

---

## Experiment 1: Direction Instability Across CAG Repeat Variants

### Background

The transportability analysis found that PCA subspaces rotate near-completely
for any repeat-count change. Direction instability (DI) — 1 minus mean pairwise
cosine of unit-normalized signatures — provides a complementary per-position
metric that does not depend on PCA or subspace dimensionality.

### Hypotheses

#### H1: Direction instability increases with repeat-count deviation from WT

DI computed on per-nucleotide RNA-FM embeddings will be higher for CAG repeat
variants further from the WT (r=21) than for variants closer to WT.

**Decision criterion:** Spearman correlation between |r - 21| and mean DI > 0.5
with p < 0.05 (two-sided). Computed across r ∈ {10, 15, 18, 21, 24, 27, 30, 36, 40, 50, 60, 80}.

#### H2: Flanking-sequence DI is lower than CAG-region DI

The CAG repeat region will show higher direction instability than the
50-nt flanking regions, because the repeat itself is the primary structural
variable while flanking sequence is conserved.

**Decision criterion:** Mean DI in CAG positions > mean DI in flanking positions,
with effect size (Cohen's d) > 0.5. Tested per-position with Mann-Whitney U.

#### H3: Direction instability is layer-dependent

Middle layers (4–8) will show more informative DI patterns than early layers
(0–2, dominated by positional encoding) or late layers (9–12, collapsed to
rank-1). Specifically, the Spearman correlation from H1 will be strongest in
layers 4–8.

**Decision criterion:** The maximum Spearman |ρ| across layers occurs in
layers 4–8, not in layers 0–2 or 9–12.

#### H4: DI dissociates from magnitude instability

Some CAG variants will show high direction instability but low magnitude
instability (or vice versa), confirming that DI captures structural change
rather than scale change.

**Decision criterion:** |Spearman(DI, magnitude_CV)| < 0.7, indicating the
two metrics are not redundant.

### Scorer (frozen)

```python
def direction_instability(embeddings_list):
    """
    embeddings_list: list of (n_positions, d_model) tensors,
                     one per CAG variant, aligned to the same positions.
    Returns: per-position DI values (n_positions,)
    """
    n_variants = len(embeddings_list)
    n_pos = embeddings_list[0].shape[0]

    # Unit-normalize each position's embedding
    normed = [e / e.norm(dim=-1, keepdim=True) for e in embeddings_list]

    # Mean pairwise cosine per position
    di = np.zeros(n_pos)
    count = 0
    for i in range(n_variants):
        for j in range(i + 1, n_variants):
            cos = (normed[i] * normed[j]).sum(dim=-1).numpy()
            di += cos
            count += 1
    di = 1 - (di / count)
    return di


def magnitude_instability(embeddings_list):
    """CV of embedding norms per position across variants."""
    norms = np.stack([e.norm(dim=-1).numpy() for e in embeddings_list])
    return norms.std(axis=0) / (norms.mean(axis=0) + 1e-8)
```

### Position alignment strategy

Because different CAG repeat counts produce different numbers of nucleotides,
embeddings are NOT directly position-aligned in the repeat region. We align on:
- Left flank: positions -50 to -1 relative to CAG start (always 50 positions)
- Right flank: positions 0 to +49 relative to CAG end (always 50 positions)
- CAG region: compute DI per-codon (average embedding over each CAG triplet),
  using the first min(r, 21) codons from each variant for alignment

This ensures position-wise comparison is meaningful despite varying repeat counts.

---

## Experiment 2: EternaFold Cross-View Validation

### Background

RNA-FM embeddings provide one "view" of RNA structure (learned from evolutionary
data). EternaFold provides an independent physics-based view (thermodynamic
ensemble prediction). Cross-view validation checks whether structural features
identified by RNA-FM correspond to genuine structural properties.

### Hypotheses

#### H5: EternaFold base-pair probabilities correlate with RNA-FM embedding PCA

The per-nucleotide base-pair probability profile from EternaFold will correlate
with the projection onto the top principal component of middle-layer (layer 6)
RNA-FM embeddings.

**Decision criterion:** Pearson r > 0.3 between EternaFold BP probability
and RNA-FM PC1 projection at layer 6, computed on the WT HTT CAG + flanking region.

#### H6: Cross-view subspace overlap exceeds chance

The Grassmannian distance between the PCA subspace of RNA-FM embeddings and
the PCA subspace of EternaFold structural features (BP probability, accessibility,
loop probability) will be smaller than the distance obtained from a null model
(shuffled sequences).

**Decision criterion:** Geodesic distance (real) < 5th percentile of geodesic
distances from 100 shuffled-sequence null models.

#### H7: EternaFold transportability parallels RNA-FM transportability

If EternaFold structural subspaces also rotate categorically across repeat lengths
(matching the RNA-FM result), this confirms the collapse is a genuine property of
HTT CAG structure, not an artifact of the RNA-FM architecture. If EternaFold
shows gradual degradation where RNA-FM shows categorical collapse, this implies
RNA-FM's representation is more sensitive to repeat count than the underlying physics.

**Decision criterion:** Report and compare. If EternaFold geodesic distances
show Spearman ρ > 0.7 with RNA-FM geodesic distances, the collapse is structural.
If ρ < 0.3, the collapse is architectural (RNA-FM-specific).

### EternaFold features (frozen)

For each HTT variant sequence, compute:
1. `bp_prob[i]`: probability that nucleotide i is base-paired (any partner)
2. `accessibility[i]`: probability that nucleotide i is unpaired
3. `loop_prob[i]`: probability that nucleotide i is in a loop
4. `mfe_structure`: minimum free energy secondary structure (dot-bracket)
5. `ensemble_energy`: partition function free energy

Stack [bp_prob, accessibility, loop_prob] as a (n_positions, 3) feature matrix
for PCA subspace analysis. Use the same position-alignment strategy as Experiment 1.

---

## Statistical corrections

All hypothesis tests use Bonferroni correction across the 7 hypotheses
(α = 0.05/7 ≈ 0.007). Effect sizes are reported regardless of significance.

## Data and code availability

All scorer code, analysis scripts, and results will be committed to the
causal-rna repository. The commit SHA of this pre-registration document
proves that hypotheses and decision criteria were frozen before the
direction instability and EternaFold experiments were run.
