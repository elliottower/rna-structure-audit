# Preregistration: Cross-Architecture RNA Structure-Awareness Study (Re-run with Corrected Data)

**Title:** Do RNA foundation models encode secondary structure? Mutation sensitivity, structure probing, and attention-to-contact analysis across architectures

**Author:** Elliot Tower (elliot@elliottower.ai)

**Date:** 2026-07-12 (pre SHA freeze)

**Commit SHA:** `[TO BE FILLED AFTER COMMIT]`

**Epistemic status:** This is a preregistration for a **corrected re-run** of
an initial cross-architecture study. The initial run produced interpretable
mutation sensitivity results for RNA-FM but suffered from data quality issues
(11/12 RNA families had sequence/structure length mismatches due to truncation)
and methodology gaps (NT v2 probing and attention failed due to 6-mer
tokenization, probing was underpowered). The hypotheses below are informed by
the initial results and therefore encode expectations from that run. We state
this openly: the SHA freeze guarantees the scorer code, hypotheses, and kill
criteria were fixed before the corrected data is processed, but the hypotheses
themselves are not blind to the initial (flawed) signal.

---

## Prior results (from initial run, known to inform this protocol)

The following results from the initial (flawed) run motivate the design below.
**All numbers in this section are from truncated/mismatched data, are unreliable,
and must not be cited as results. They appear here solely to explain why the
hypotheses below have the directional expectations they do.**

- **RNA-FM** (BERT architecture, 99M params): Mean mutation sensitivity ratio
  5.47x (stem cosine distance / loop cosine distance) across complement-swap
  mutations. Best layers varied by family but were consistently in middle-to-late
  layers. **Unreliable — 11/12 families had mismatched data.**
- **NT v2** (ESM architecture, 56M params): Mean ratio 1.37x. Near-chance,
  suggesting minimal structure encoding via this metric.
  **Unreliable — same data quality issues, plus 6-mer alignment not implemented.**
- **HyenaDNA** (Hyena SSM, 5.4M params): Mean ratio 1.77x. Slightly above
  chance but far below RNA-FM. **Unreliable — same data quality issues.**
- **Data quality issues:** 11 of 12 RNA structure entries had sequence/structure
  length mismatches (truncated sequences or structures). Only 1 family had clean
  data. All results above are therefore unreliable and require re-running.
- **NT v2 failures:** Structure probing and attention-to-contact analyses
  produced no results for NT v2 because 6-mer tokenization creates tokens that
  span multiple nucleotide positions, making per-nucleotide label assignment
  ill-defined. This requires explicit 6-mer-to-nucleotide alignment.
- **Probing underpowered:** 926 total positions across all families, with
  class balance matching majority baseline. Probe accuracy matched chance,
  but this is likely a power issue rather than a genuine null result.

---

## Hypotheses

### Confirmatory hypotheses (Bonferroni-corrected, alpha = 0.05/6 = 0.0083)

**H1: RNA-FM mutation sensitivity exceeds untrained baseline by a factor of ≥ 2.**

RNA-FM's BERT architecture with self-attention over individual nucleotides
should encode base-pairing relationships that random weights do not. The
trained model's mutation sensitivity ratio should substantially exceed the
untrained baseline's ratio, establishing that the signal is learned rather
than architectural.

- **Metric:** Trained RNA-FM mean ratio / untrained RNA-FM mean ratio
  (both computed as mean across families of stem/loop cosine distance ratio).
- **Decision criterion:** Trained/untrained ratio ≥ 2.0 with the trained
  model's mean ratio also exceeding the nucleotide-stratified null (H1b below).
- **Effect size:** Report the trained/untrained ratio and Cohen's d of
  trained vs untrained per-family ratios.
- **Outcome if trained/untrained < 2.0:** Structure-awareness not supported
  — the signal does not substantially exceed what random weights produce.
  No dead zone: any trained/untrained ratio < 2.0 is a failure, regardless
  of the absolute value.

**H1b: RNA-FM mutation sensitivity exceeds the nucleotide-stratified null.**

The mutation sensitivity metric could be confounded by nucleotide composition
differences between stem and loop positions. Stems are GC-enriched (68.5% GC
vs 43.9% in loops across our 12 families, 12/12 families show stem GC > loop
GC, Wilcoxon p < 0.001). If complement-swap distances differ by nucleotide
identity (G→C producing larger perturbations than A→U), a composition-blind
null would be systematically deflated. To control for this, we use a
nucleotide-stratified permutation: for each nucleotide type (A, U, G, C)
independently, shuffle which positions of that nucleotide are labeled "stem"
vs "loop." This preserves the nucleotide composition of each group exactly.
Recompute the ratio 100 times to build a null distribution. The real ratio
must exceed the 95th percentile of this composition-matched null.

- **Metric:** For each family, compute the real stem/loop ratio and the
  distribution of 100 nucleotide-stratified shuffled-label ratios.
- **Decision criterion:** RNA-FM's real mean ratio exceeds the 95th percentile
  of the stratified null distribution, with a Wilcoxon signed-rank test
  (real ratio - 95th percentile null ratio, paired across families) at
  p < 0.0083.
- **Rationale:** This is the composition control. d_g was falsified because
  a zero-parameter 3-mer baseline reproduced the signal. A global label
  shuffle would not properly control for composition because it destroys the
  GC enrichment correlation with stem labels (see Composition Control section).
  The stratified null tests whether mutation sensitivity exceeds what
  nucleotide identity alone predicts — the shuffled "stem" and "loop" groups
  have identical nucleotide makeup to the real groups, so any excess in the
  real ratio reflects structural context beyond composition.

**H2: RNA-FM shows layer localization of structure encoding.**

If structure-awareness is learned (rather than an artifact of tokenization or
embedding initialization), the best mutation sensitivity ratios should appear
in intermediate or late layers, not at layer 0 (embedding lookup).

- **Metric:** For each RNA family, identify the layer with the highest
  stem/loop ratio. Report the distribution of best-layer indices.
- **Decision criterion:** Median best-layer index > 3 (i.e., past the
  embedding layer and first few transformer blocks). Wilcoxon signed-rank
  test on (best_layer - 3) against 0 at p < 0.0083 (tests whether the
  distribution of best layers is shifted above 3, not just above 0).

**H3: NT v2 mutation sensitivity mean ratio < 2.0.**

NT v2 uses 6-mer tokenization, which obscures individual nucleotide identity.
A single nucleotide mutation changes at most 6 overlapping tokens, distributing
the signal. This architectural choice should limit structure sensitivity.

- **Metric:** Same as H1, applied to NT v2.
- **Decision criterion:** Mean ratio < 2.0. One-sided t-test with
  H_0: mu >= 2.0, reject at p < 0.0083.

**H4: HyenaDNA mutation sensitivity mean ratio < 2.0.**

HyenaDNA uses a Hyena SSM without explicit attention. SSMs process sequences
through convolution-like operations that capture local dependencies but may
not encode the long-range base-pairing relationships that define RNA secondary
structure.

- **Metric:** Same as H1, applied to HyenaDNA.
- **Decision criterion:** Same as H3.

**H5: RNA-FM probing accuracy > majority baseline with corrected data.**

A linear probe trained on RNA-FM embeddings to classify paired vs unpaired
positions should exceed the majority-class baseline when trained on corrected,
properly aligned sequence/structure data.

- **Metric:** 5-fold cross-validated balanced accuracy of a logistic regression
  probe on per-nucleotide embeddings labeled as paired (in a stem) or unpaired
  (in a loop/bulge/junction).
- **Decision criterion:** Best-layer mean balanced accuracy > majority
  baseline + 5 percentage points, with paired t-test across folds at
  p < 0.0083. "Best layer" is selected by highest mean balanced accuracy;
  this selection is part of the metric definition, not a multiplicity issue
  (analogous to reporting the best-layer mutation sensitivity ratio).
- **Baseline:** Majority class frequency in the corrected dataset.

### Exploratory hypotheses (effect sizes reported, no p-threshold)

**H6: Caduceus (bidirectional Mamba, 14M params) mutation sensitivity ratio.**

Caduceus uses a bidirectional SSM (Mamba) without attention. Bidirectionality
allows information from both ends of a stem to meet, which could enable
structure encoding even without attention. Predicted ratio < 2.0 based on
the SSM limitation, but the bidirectional design makes this less certain
than the HyenaDNA prediction.

- **Metric:** Same as H1, applied to Caduceus.
- **Report:** Mean ratio, per-family ratios, best-layer indices.

**H7: Evo (hybrid Hyena+attention, 7B params) mutation sensitivity ratio.**

Evo combines Hyena SSM layers with sparse attention layers and is substantially
larger (7B params) than the other models. The hybrid architecture and scale
could produce intermediate structure-awareness. Predicted ratio between
HyenaDNA (1.77x in initial run) and RNA-FM (5.47x).

- **Metric:** Same as H1, applied to Evo.
- **Report:** Mean ratio, per-family ratios, best-layer indices.
  Also report whether best layers coincide with attention layers or SSM layers.

**H8: Bidirectionality sufficiency test.**

If Caduceus shows a ratio > 3.0, this suggests bidirectional context
propagation is sufficient for structure encoding and explicit attention is
not required. This would be a meaningful finding about architectural
requirements for structure-awareness.

- **Metric:** Caduceus mean ratio from H6.
- **Decision:** If ratio > 3.0, bidirectionality suffices. If ratio < 2.0,
  attention appears necessary. If 2.0--3.0, inconclusive.

**Contingency for H6-H8:** Caduceus requires `mamba_ssm` compilation on a
CUDA devel image, which has failed in previous attempts. Evo requires A100
GPUs. If either model fails to build or run, the corresponding hypotheses
are reported as "not executed" with the specific build/runtime error, not
silently dropped. The confirmatory hypotheses (H1-H5) do not depend on
Caduceus or Evo results.

---

## Kill criteria

These criteria define what would falsify the claim that RNA-FM is
structure-aware via learned representations.

**K1: RNA-FM trained/untrained ratio < 2.0.**

H1 requires the trained model to exceed the untrained baseline by a factor
of ≥ 2. If this fails, the signal is not substantially learned. There is no
dead zone: K1 fires whenever H1 fails. (The absolute trained ratio is also
reported but is not the decision variable — an absolute ratio of 3.0 that
the untrained model also achieves at 2.5 is not evidence of learned
structure-awareness.)

**K4: RNA-FM real ratio does not exceed nucleotide-stratified null.**

If the real stem/loop ratio falls within the 95th percentile of the
nucleotide-stratified null distribution (H1b), the mutation sensitivity metric
is confounded by local sequence composition differences between stem and
loop neighborhoods, not by secondary structure. This is the composition
kill criterion — analogous to the 3-mer baseline that killed d_g.

**K2: RNA-FM best ratios appear at layer 0 (embedding layer).**

If the highest stem/loop ratios consistently appear at the embedding layer
rather than intermediate transformer layers, the signal reflects tokenization
or embedding initialization rather than learned representations. The claim
requires that structure encoding emerges through processing, not input
formatting.

**K3: Untrained RNA-FM baseline produces comparable ratios.**

If an RNA-FM model with randomly initialized weights (same architecture,
no pretraining) produces mutation sensitivity ratios within 80% of the trained
model's ratios, the signal is architectural (e.g., positional encoding
interacting with sequence similarity) rather than learned from RNA data.

- **Metric:** Ratio of untrained RNA-FM mean ratio to trained RNA-FM mean
  ratio.
- **Kill threshold:** Untrained ratio / trained ratio > 0.8.

---

## Analysis plan

### Primary metric: mutation sensitivity ratio

For each RNA family and each model:

1. Identify stem (base-paired) and loop (unpaired) positions from the
   secondary structure annotation.
2. For EVERY position in the sequence, construct a complement-swap mutant:
   - A <-> U, C <-> G (Watson-Crick complements).
   - Replace the nucleotide at that position with its complement.
   - Positions where the complement equals the original (should not occur
     for standard RNA nucleotides) are skipped.
3. Extract embeddings from the original and each mutant sequence at all layers.
4. Compute cosine distance between original and mutant embeddings at the
   mutated position, for each layer.
5. For each layer, compute:
   - Mean stem cosine distance (across ALL stem positions with valid mutations)
   - Mean loop cosine distance (across ALL loop positions with valid mutations)
   - Ratio = stem mean / loop mean
6. Report:
   - Best ratio across layers, and the layer index where it occurs
   - Per-family table of ratios and best layers
   - Number of stem and loop positions mutated per family
   - Mean ratio across families as summary statistic (unweighted)
   - **Median** ratio across families as robustness check (protects against
     a single high-ratio outlier family dominating the mean)
   - Bootstrap 95% CI on the mean ratio (10,000 resamples of families)
7. **Small-N handling:**
   - Families with fewer than 3 stem or 3 loop positions are excluded from
     the per-family ratio and flagged in the results table.
   - If the mean and median ratios disagree in direction (one > 2.0,
     the other < 2.0), report the result as inconclusive for H1.
9. **Statistical power limitation (N=12 families):**
   - With 12 families and Bonferroni alpha=0.0083, hypothesis tests have
     limited power to detect moderate effects. Bootstrap CIs on 12
     observations will be wide. GroupKFold probing (5 folds of ~2-3
     families each) produces high-variance fold-level estimates.
   - Confirmatory claims are therefore restricted to strong effects
     (trained/untrained >= 2.0, clear null exceedance). Any marginal
     result (e.g., p between 0.01 and 0.05 before correction, or
     trained/untrained between 1.5 and 2.0) is treated as exploratory,
     not confirmatory, regardless of nominal significance.
   - This is a hard constraint of the current 12-family dataset, not an
     analysis choice. Expanding to 50+ families from Rfam is the path
     to confirmatory power for moderate effects.

### Composition control: nucleotide-stratified null

For each model and each RNA family:

1. Compute per-position cosine distances for ALL positions (as in the primary
   metric above). This produces a distance vector over all valid positions.
2. For each permutation, shuffle stem/loop labels **within each nucleotide
   type independently** (all A positions are permuted among themselves, all U
   positions among themselves, etc.). This is a nucleotide-stratified
   permutation. Seed: 42 + rna_idx * 1000 + layer_idx.
3. Recompute the stem/loop cosine distance ratio using the shuffled labels
   applied to the same pre-computed distances.
4. Repeat 100 times per family per layer to build a null distribution.
5. Compare the real ratio to the null distribution. Report the percentile
   rank of the real ratio within the null.
6. For the aggregate test (H1b), compute the mean real ratio and the 95th
   percentile of mean shuffled ratios across families.

**Why nucleotide-stratified, not global label shuffle:** Stem positions are
GC-enriched relative to loop positions (pooled: 68.5% GC in stems vs 43.9%
in loops, 12/12 families, Wilcoxon p < 0.001). A global label shuffle would
destroy this composition correlation, producing shuffled "stem" groups with
lower GC than real stems. Since complement-swap distances may differ by
nucleotide identity (G→C vs A→U), the global null would be systematically
deflated — making it too easy for the real ratio to exceed it. Stratified
permutation preserves the per-nucleotide composition of each group: shuffled
"stem" positions contain the same number of A, U, G, C as real stems. Any
remaining signal above the stratified null reflects structure beyond what
nucleotide identity alone predicts.

**Zero-entropy strata:** In some families, a nucleotide type occurs exclusively
in stems or exclusively in loops (e.g., hammerhead: all A's and U's are in
loops; mir-21: all G's are in stems). These strata have C(n, 0) = 1 distinct
relabeling and contribute zero shuffle entropy — they are uninformative by
construction. The test draws its power from the remaining strata that contain
both labels. This is a correct property of the stratification, not a
limitation: if a nucleotide is perfectly structure-correlated, there is no
within-nucleotide variation to test against, and the null should not pretend
otherwise.

### Untrained baseline control

For RNA-FM only:

1. Initialize an RNA-FM model with the same architecture but random weights
   (no pretrained checkpoint).
2. Run the identical mutation sensitivity pipeline on the untrained model.
3. Compare untrained ratios to trained ratios per family and in aggregate.
4. This is the primary control for K3.

### Secondary metric: structure probing

For each model with corrected data:

1. Extract per-nucleotide embeddings at each layer for all RNA families.
2. Label each position as paired (1) or unpaired (0) from the structure
   annotation.
3. Train a logistic regression probe (sklearn, default regularization) using
   5-fold cross-validation, stratified by RNA family.
4. Report balanced accuracy per fold, mean across folds, and comparison to
   majority baseline.
5. For NT v2: map 6-mer token embeddings to nucleotide positions by
   assigning each nucleotide i the embedding of token i // 6 (the
   non-overlapping 6-mer that covers position i).

### Secondary metric: attention-to-contact correlation

For attention-based models (RNA-FM, NT v2) only:

1. Construct binary contact maps from secondary structure annotations
   (1 if positions i,j are base-paired, 0 otherwise).
2. Extract attention matrices from all heads at all layers.
3. For NT v2: align 6-mer attention to nucleotide-level attention by
   distributing each 6-mer token's attention weight uniformly across its
   6 constituent nucleotides (for both the attending and attended dimensions).
4. Compute Spearman correlation between each attention head's weights and
   the contact map (upper triangle only, excluding diagonal).
5. Report: best head per layer, best head overall, distribution of
   correlations across heads.
6. **NT v2 resolution caveat:** The uniform distribution of 6-mer attention
   across constituent nucleotides is an assumption that may create artifacts.
   Near-zero NT v2 correlation should be interpreted as "no signal at this
   resolution," not as definitive evidence that NT v2 lacks structure-awareness.
   This metric is secondary for NT v2; mutation sensitivity (which uses
   per-position embeddings at the mutated nucleotide, with NT v2 6-mer
   tokens expanded to nucleotide resolution) is the primary metric for
   all models.

---

## Data

### RNA families

12 RNA families with secondary structure annotations sourced from PDB
crystal structures or Rfam consensus structures. The initial run used
entries with sequence/structure length mismatches (11/12 truncated).
The corrected dataset will:

1. Re-download or re-extract sequences to ensure sequence length matches
   structure annotation length exactly.
2. Verify that every position in the sequence has a corresponding structure
   character (dot-bracket notation).
3. Validate that annotated base pairs are Watson-Crick or wobble pairs
   (G-U) at the sequence level.

Families will span a range of structure complexity (simple stem-loops,
multi-branch junctions, pseudoknots if available in dot-bracket).

### Reproducibility parameters

- **Random seed:** 42 + rna_idx for per-family reproducibility;
  42 + rna_idx * 1000 + layer_idx for shuffled null permutations.
- **Random seed for untrained baseline:** 42 for weight initialization.
- **Probing CV:** 5-fold, stratified by RNA family (GroupKFold).
- **Positions mutated per family:** ALL positions (no sampling).

---

## Models

| Model | Architecture | Params | Tokenization | Directionality |
|-------|-------------|--------|--------------|----------------|
| RNA-FM | BERT (transformer encoder) | 99M | Per-nucleotide | Bidirectional |
| NT v2 | ESM (transformer encoder) | 56M | 6-mer | Bidirectional |
| HyenaDNA | Hyena SSM | 5.4M | Per-nucleotide | Unidirectional |
| Caduceus | Mamba SSM | 14M | Per-nucleotide | Bidirectional |
| Evo | Hybrid Hyena+attention | 7B | Per-nucleotide | Primarily unidirectional |

---

## Statistical corrections

Confirmatory hypotheses {H1, H1b, H2, H3, H4, H5} use Bonferroni correction
(alpha = 0.05/6 = 0.0083). Exploratory hypotheses {H6, H7, H8} report effect
sizes without p-value thresholds.

All results are reported regardless of significance. Effect sizes are primary;
p-values are secondary.

---

## Reporting plan

Results will be reported in full regardless of outcome. The report will include:

1. **Per-family table:** For each model and each RNA family, the best-layer
   mutation sensitivity ratio and the layer index.
2. **Summary statistics:** Mean ratio across families per model, with
   confidence intervals.
3. **Kill criteria evaluation:** Explicit pass/fail for K1, K2, K3, K4 with
   the numerical values that determined each decision.
4. **Probing results:** Balanced accuracy per model per layer, compared to
   majority baseline.
5. **Attention-to-contact results:** Best Spearman correlation per model,
   with the responsible head identified.
6. **Untrained baseline comparison:** Side-by-side trained vs untrained
   RNA-FM ratios (H1 decision variable).
7. **Nucleotide-stratified null:** Per-family percentile rank of real ratio
   within nucleotide-stratified null distribution, and aggregate pass/fail
   for H1b/K4.

---

## Commit SHA protocol

After this document is committed, the SHA will be recorded above. The SHA
proves that hypotheses, kill criteria, metrics, and analysis code were frozen
before the corrected data experiments were run. Any post-hoc additions will
be clearly marked as such in the results document.

**Commit SHA:** `[TO BE FILLED AFTER COMMIT]`
