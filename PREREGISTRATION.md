# Frozen Analysis Protocol: Direction Instability and Cross-View Validation of RNA-FM on HTT mRNA

**Title:** Direction instability of RNA foundation model embeddings across CAG repeat variants, with thermodynamic cross-view validation

**Authors:** Elliot Tower

**Date:** 2026-07-11 (pre SHA freeze)

**Commit SHA:** _TO BE FILLED AFTER COMMIT_

**Epistemic status:** This is a **frozen analysis protocol**, not a pre-registration
in the strict sense. The hypotheses below are informed by results from a prior
analysis (the causal geometry audit of RNA-FM, reported separately). We observed
a categorical transportability collapse (geodesic distance ~3.7 for any non-WT
repeat count) and extreme late-layer rank-1 compression before formulating these
hypotheses. The SHA freeze guarantees the **scorer code and decision criteria**
were fixed before the direction instability and cross-view experiments ran, but
the hypotheses themselves encode expectations formed after seeing the geometry.
We state this openly rather than overclaiming pre-registration status.

---

## Prior results (already observed, motivate this protocol)

The following results from the initial audit are known and inform the design below:

- Grassmannian transportability collapses categorically: geodesic distance
  ~3.7 for ANY non-WT repeat count, including normal-range r=27.
  The collapse is binary (0.14 same / ~3.7 different), not gradual.
- Layer-wise effective dimensionality: 51 → 1 across layers 0–12.
  Layers 9–12 have effective dim = 1 (97.9% variance in PC1).
- W_Q effective rank: 265 → 86 across layers (3x compression).
- Per-position embedding norms are constant at 28.667 (LayerNorm effect).

---

## Experiment 1: Direction Instability Across CAG Repeat Variants

### Background

The transportability analysis found that PCA subspaces rotate near-completely
for any repeat-count change. This binary behavior could arise from two sources:
(a) the PCA subspace is dominated by length-dependent structure (different numbers
of positions create different principal components), or (b) RNA-FM genuinely
represents different repeat counts as categorically different structures.

Direction instability (DI) — computed as a per-variant scalar metric against
the WT reference — provides a complementary test. Unlike Grassmannian distance,
DI operates on aligned individual positions rather than global subspaces, so it
can detect gradual per-position changes that the subspace metric averages away.

**Critical design note:** Because Paper 1 already showed the subspace metric
is binary, we specifically test whether DI *also* shows binary behavior (H1a)
or smooth monotonic behavior (H1b). Pre-committing to both possibilities avoids
designing H1 to fail.

### Hypotheses

#### Confirmatory hypotheses (Bonferroni-corrected, α = 0.05/4 = 0.0125)

**H1: DI is structured with respect to repeat-count deviation.**
For each CAG variant r, compute the mean cosine distance of its aligned
embedding (at layer 6) against the WT (r=21) embedding. Test whether
these per-variant DI values correlate with |r − 21|.

- **H1a (smooth):** Spearman ρ > 0.5 with p < 0.0125.
- **H1b (binary, matching Paper 1):** DI is bimodal — WT cluster
  (r ∈ {21}) vs all others — with Mann-Whitney U p < 0.0125
  separating the two groups.
- **Decision:** Report which pattern obtains. If both pass, smooth
  subsumes binary. If neither passes, the metric does not capture
  repeat-count variation at this layer.
- **Power note:** With n=12 variants and only 1 WT observation, the
  Spearman test has limited power. We report the effect size regardless
  of significance.

**H2: CAG-region DI exceeds flanking-region DI.**
The 50-nt flanking regions (conserved across variants) will show lower DI
than the CAG repeat region. Because mean-pooling the CAG to a single
position (as used for H1) would leave n=1 in the CAG group — insufficient
for a group comparison — H2 uses **per-codon alignment**: each CAG codon
(3 nt) is mean-pooled to one vector, and the first 10 codons (the minimum
repeat count in the variant set) are retained. This gives 10 CAG positions
vs 100 flank positions for the Mann-Whitney test.

- **Decision criterion:** Cohen's d > 0.5 and Mann-Whitney U p < 0.0125.
- **Null baseline:** Shuffled-sequence DI. For each variant, shuffle the
  nucleotide order of the full sequence (preserving composition), extract
  embeddings with the same per-codon alignment, compute DI. Repeat 50
  times. H2 passes only if real CAG DI exceeds the 95th percentile of
  shuffled CAG DI.

**H5: ViennaRNA base-pair probabilities correlate with RNA-FM embedding structure.**
Per-nucleotide base-pair probability from ViennaRNA partition function will
correlate with the projection onto PC1 of middle-layer (layer 6) RNA-FM
embeddings on the WT HTT CAG + flanking region.

- **Decision criterion:** Pearson |r| > 0.3 with p < 0.0125.

**H6: Cross-view subspace overlap exceeds chance.**
The Grassmannian distance between the PCA **score** subspaces (k=2) of
RNA-FM layer-6 embeddings and ViennaRNA structural features (BP probability,
loop type) will be smaller than the null distribution from 100
composition-matched shuffled sequences. Both PCA score matrices are
(n_positions, k), living in the shared position-space R^n_positions, so
`subspace_angles` is well-defined. Loadings live in incompatible feature
spaces (R^640 vs R^2) and cannot be compared directly.

- **Decision criterion:** Real geodesic distance < 5th percentile of
  shuffled-sequence geodesic distances.

#### Exploratory analyses (effect sizes reported, no p-threshold)

**H3: DI informativeness is layer-dependent.**
Compute H1's per-variant DI metric at all 13 layers (0–12). Report which
layer yields the strongest effect. Prediction (from Paper 1's layer-wise
compression): middle layers (4–8) will be most informative.

- **Report:** |Spearman ρ| per layer, best layer index. No statistical
  test — this is a descriptive exploration of layer sensitivity.

**H4: DI dissociates from magnitude instability.**
Some CAG variants may show high direction instability but low magnitude
instability (embedding norm CV), confirming DI captures structural change
rather than scale change.

- **Report:** Spearman ρ between per-variant DI and per-variant magnitude
  CV. Dissociation indicated if |ρ| < 0.7.

**H7: ViennaRNA transportability compared to RNA-FM transportability.**
Compute within-modality geodesic distances: for each CAG variant, fit
PCA on its (n_positions, n_features) matrix, QR-orthonormalize the score
matrix, then measure subspace_angles between the variant's orthonormalized
scores and WT's scores in position-space. Both modalities use k=2 so
geodesic magnitudes are on the same scale (ViennaRNA has effective rank 2,
so k=2 is the maximum; RNA-FM uses the same k for comparability).
Compare geodesic distances to WT across the two modalities.

- **Decision:** Report Spearman ρ between ViennaRNA and RNA-FM geodesic
  distances across variants.
  - ρ > 0.7: The collapse is a genuine property of HTT CAG structure.
  - ρ < 0.3: The collapse is architectural (RNA-FM-specific).
  - 0.3 ≤ ρ ≤ 0.7: Partial concordance — both structural and
    architectural factors contribute. Report both distances and note
    which variants diverge.
- **Dimensionality caveat:** ViennaRNA geodesics are computed on a 2D
  feature space ([bp_prob, loop_type]) and are therefore a coarse
  comparator. A weak ρ may reflect ViennaRNA's low feature dimensionality
  rather than genuine architectural divergence between the two views.

---

## Scorer Functions (frozen)

### Per-variant DI against WT (for H1, H3, H4)

```python
def di_vs_wt(embeddings_dict, wt_key=21):
    """
    embeddings_dict: {repeat_count: (n_aligned, d_model) tensor}
    All tensors must be position-aligned (same n_aligned).
    Returns: {repeat_count: scalar DI} — mean cosine distance to WT.
    """
    wt = embeddings_dict[wt_key]
    wt_normed = wt / (wt.norm(dim=-1, keepdim=True) + 1e-8)
    result = {}
    for r, emb in embeddings_dict.items():
        emb_normed = emb / (emb.norm(dim=-1, keepdim=True) + 1e-8)
        cos = (wt_normed * emb_normed).sum(dim=-1)  # (n_aligned,)
        result[r] = float(1 - cos.mean())
    return result
```

### Per-position DI across all variants (for H2)

```python
def direction_instability_profile(embeddings_list):
    """
    embeddings_list: list of (n_aligned, d_model) tensors.
    Returns: per-position DI (n_aligned,) — 1 minus mean pairwise cosine.
    """
    n = len(embeddings_list)
    normed = [e / (e.norm(dim=-1, keepdim=True) + 1e-8) for e in embeddings_list]
    cos_sum = np.zeros(normed[0].shape[0])
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            cos_sum += (normed[i] * normed[j]).sum(dim=-1).numpy()
            count += 1
    return 1 - (cos_sum / count)
```

### Magnitude instability (for H4)

```python
def magnitude_instability(embeddings_dict, wt_key=21):
    """CV of embedding norms for each variant relative to WT."""
    result = {}
    wt_norms = embeddings_dict[wt_key].norm(dim=-1).numpy()
    for r, emb in embeddings_dict.items():
        norms = emb.norm(dim=-1).numpy()
        result[r] = float(np.std(norms - wt_norms) / (np.mean(np.abs(wt_norms)) + 1e-8))
    return result
```

---

## Position Alignment Strategy

Because different CAG repeat counts produce different numbers of nucleotides,
direct positional alignment in the repeat region is not possible. Two alignment
modes are used, depending on the hypothesis:

### Mode A: Mean-pool alignment (H1, H3, H4, H5, H6, H7)

1. **Left flank:** Positions −50 to −1 relative to CAG start. Always 50 positions,
   identical sequence across all variants. Directly aligned.

2. **Right flank:** Positions +1 to +50 relative to CAG end. Always 50 positions,
   identical sequence across all variants. Directly aligned.

3. **CAG region:** Mean-pool all CAG nucleotide embeddings into a single
   (1, d_model) vector per variant. This avoids the per-codon alignment problem
   (codon 5 of an 80-mer is in a different structural context than codon 5 of
   a 21-mer). The mean-pooled vector summarizes the aggregate representation of
   the repeat region regardless of length.

4. **Full aligned representation:** Concatenate [left_flank (50,d), cag_mean (1,d),
   right_flank (50,d)] → (101, d_model) per variant.

**Worked example for r=40:**
- Input sequence: `[50nt flank_L] [CAG×40 = 120nt] [50nt flank_R]` = 220 nt
- RNA-FM embeddings: (220, 640) after stripping CLS/EOS
- Left flank: positions 0–49 → (50, 640)
- CAG mean: mean(positions 50–169) → (1, 640)
- Right flank: positions 170–219 → (50, 640)
- Aligned: (101, 640)

### Mode B: Per-codon alignment (H2 only)

H2 compares the DI *distribution* at CAG positions against flank positions,
so it requires multiple CAG positions (not a single mean-pooled vector).
Each CAG codon (3 nucleotides) is mean-pooled to one (1, d_model) vector.
The first 10 codons (= minimum repeat count in the variant set) are retained,
giving a fixed-size CAG representation across all variants.

**Full aligned representation:** Concatenate [left_flank (50,d), cag_codons (10,d),
right_flank (50,d)] → (110, d_model) per variant.

This provides n=10 CAG DI values vs n=100 flank DI values for the
Mann-Whitney test, avoiding the n=1 problem that mean-pool alignment
would create.

**CLS/EOS handling:** Special tokens (CLS at position 0, EOS at final position)
are stripped before alignment. All indices above refer to post-stripping positions.

---

## Null Model (for H2, H6)

Null sequences are generated by **full-sequence shuffle**: the entire variant
sequence (flanks + CAG region) is randomly permuted, preserving nucleotide
composition and length but destroying all local structure including the
CAG/flank boundary. Embeddings and ViennaRNA features are then extracted
from the shuffled sequence and sliced at the **same positional indices** as
the real sequence (i.e., positions 0–49 as "left flank," the middle region
as "CAG," and the final 50 as "right flank"). After shuffling, these
indices no longer correspond to flanks or CAG repeats — they are arbitrary
positions in a random sequence. This is the intended null: it tests whether
the observed DI or cross-view alignment reflects genuine structural
sensitivity or is an artifact of embedding geometry that would arise for
any composition-matched sequence.

- **H2 null:** 50 shuffles per variant, per-codon alignment (Mode B).
  H2 passes only if real CAG DI exceeds the 95th percentile of shuffled
  CAG-index DI.
- **H6 null:** 100 shuffles of the WT sequence, mean-pool alignment
  (Mode A). H6 passes only if real cross-view geodesic < 5th percentile
  of shuffled geodesics.

---

## Thermodynamic Features (frozen)

### Primary tool: ViennaRNA

ViennaRNA (installed, verified via `RNA.fold()`) provides the baseline
physics-based view using Turner 2004 thermodynamic parameters.

For each HTT variant sequence, compute:

1. `bp_prob[i]`: Sum of base-pair probabilities involving position i.
   **Ensemble quantity** — computed from the full partition function
   via `RNA.fold_compound(seq).pf()` then `fc.bpp()`. Sum upper triangle
   entries for each position i.
2. `loop_type[i]`: Binary indicator — 1 if position i is unpaired in the
   **MFE structure only** (from `RNA.fold(seq)`), 0 if base-paired.
   **This is NOT an ensemble quantity.** It is a single-structure property
   from the minimum free energy fold. We use MFE because ViennaRNA does
   not provide per-position loop-type probabilities from the partition
   function. This limitation is stated explicitly.

Stack [bp_prob, loop_type] as a (n_positions, 2) feature matrix.
Accessibility (1 − bp_prob) is omitted because it is linearly dependent
on bp_prob and would create a rank-deficient feature matrix. PCA on
ViennaRNA features uses k=2 throughout (H6, H7).

Use the same position-alignment strategy as Experiment 1 (mean-pool the CAG
region, keep flanks aligned).

### Future extension: EternaFold

EternaFold (learned parameters from Eterna crowdsourced + experimental probing
data) produces base-pair probabilities that agree better with experimental
chemical probing than ViennaRNA's Turner parameters. EternaFold is not
pip-installable (requires C++ compilation from source) and was not available
in this session. It is deferred to a follow-up analysis.

When EternaFold is added, the same two-feature stack (bp_prob, loop_type)
will be computed from both tools, and an additional hypothesis (H8) will
test cross-tool agreement:
ViennaRNA and EternaFold BP-probability profiles should correlate at Pearson
r > 0.7 on the WT sequence, establishing baseline agreement before interpreting
disagreements on variant sequences. If RNA-FM tracks EternaFold but not ViennaRNA,
this would suggest RNA-FM learned probing-consistent structure beyond Turner
thermodynamics.

---

## Frozen sequence data

All variant sequences are constructed from the HTT mRNA reference NM_002111.7
(13,669 nt, downloaded from NCBI, stored at `data/HTT_NM_002111.7.fasta`).
The native CAG repeat region is at positions 148–211 (21 repeats, identified
by regex `(CAG){5,}`). Flanking sequences are the 50 nt immediately upstream
(positions 98–147) and downstream (positions 212–261) of the native repeat.

Variant set: r ∈ {10, 15, 18, 21, 24, 27, 30, 36, 40, 50, 60, 80}.

---

## Statistical corrections

Confirmatory hypotheses {H1, H2, H5, H6} use Bonferroni correction
(α = 0.05/4 = 0.0125). Exploratory analyses {H3, H4, H7} report effect
sizes without p-value thresholds.

All results are reported regardless of significance. Effect sizes (ρ, d, r)
are primary; p-values are secondary.

---

## Reproducibility

Random seed `SEED = 20260711` is set via `np.random.seed()` and
`torch.manual_seed()` before any stochastic operations (null model shuffles).
Re-running the frozen scorer with the same seed produces identical null
distributions and therefore identical p-values and pass/fail decisions.

---

## Data and code availability

All scorer code, analysis scripts, raw sequences, and results will be committed
to the causal-rna repository. The commit SHA of this document proves that the
scorer and decision criteria were frozen before the direction instability and
ViennaRNA experiments were run.
