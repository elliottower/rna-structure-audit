# Preregistration: Phase 6 — Pairwise Coupling Test

**Title:** Do RNA foundation models resolve specific base-pair partners? Within-stem perturbation specificity across ten models

**Author:** Elliot Tower (elliot@elliottower.ai)

**Date:** 2026-07-13

**Commit SHA:** [TO BE FROZEN BEFORE GPU EXECUTION]

**Replaces:** PREREGISTRATION_PHASE6_COMPENSATORY_MUTATION.md (SHA `e9f2149`), which used a confounded metric (destructive-vs-compensatory ratio conflates number-of-swaps with structure preservation). This preregistration defines a different metric for the same scientific question.

**Pilot data disclosure:** Coupling ratio values (a weaker variant of the metric defined here) were observed for tRNA_Phe_yeast and tRNA_Ala_human using RNA-FM before this preregistration was written. Those families are quarantined from confirmatory analysis.

---

## Motivation

Phases 1-5 tested whether models distinguish stem positions from loop positions in their response to single-nucleotide mutations. No model exceeded the composition-controlled null. That metric tests a weak form of structure awareness: stem-vs-loop discrimination.

This phase tests a stronger form: does the model encode which specific position pairs with which? A model that distinguishes stems from loops could still be ignorant of the specific pairing topology. Conversely, a model that resolves specific partners necessarily encodes pairing.

---

## Experiment design

### RNA families and models

Same 52 Rfam families and 10 models as Phases 2-4. Families with fewer than 15 canonical Watson-Crick pairs are excluded. The exact count of eligible families depends on data quality (see Families Pending below).

**Eligible families:** 34 of 52 families pass the 15-WC-pair threshold. After quarantine (2 families), N = 32 confirmatory families.

**Quarantined pilot families:** tRNA_Phe_yeast and tRNA_Ala_human are excluded from all confirmatory hypothesis tests. Pilot values were observed before this preregistration.

### Stems and eligible pairs

Parse the dot-bracket annotation into stems (maximal runs of stacked base pairs). Keep only stems with >= 3 consecutive canonical WC pairs. Within each stem, exclude the two terminal pairs (they lack a neighbor on one side). The remaining interior pairs are eligible for the primary metric.

This filtering ensures every eligible pair (i, j) has stem neighbors on both sides of j, so the within-stem comparison is symmetric.

### Primary metric: perturbation specificity (PS)

For each eligible base pair (i, j) in a stem:

1. Compute wild-type embeddings at layer L: E_wt = model(sequence).
2. Swap position i to its Watson-Crick complement. Compute mutant embeddings: E_mut = model(mutant_sequence).
3. Compute perturbation at every position k: Delta_k = cosine_distance(E_wt[k], E_mut[k]).
4. **Partner perturbation:** Delta_partner = Delta_j.
5. **Adjacent-stem perturbation:** Delta_adj = max(Delta_{j-1}, Delta_{j+1}), where j-1 and j+1 are the adjacent positions on j's strand within the same stem. Using max (not mean) is conservative: the partner must exceed the best alternative, not just the average.
6. **Per-pair PS:** PS(i,j) = Delta_partner - Delta_adj.

PS > 0 means the model perturbs i's specific partner more than the best adjacent stem position. PS = 0 means the model treats the stem uniformly. PS < 0 means the partner is less perturbed than its neighbors.

Aggregate per-family: mean PS across all eligible pairs. Compute at every layer; report the maximum (max-over-layers), consistent with Phases 1-5.

### Positive control: stem vs loop prerequisite

Before testing partner specificity, verify that the model distinguishes stem from loop at all:

- For each mutation at stem position i, compute mean Delta at other stem positions vs mean Delta at loop positions.
- If stem Delta is not significantly greater than loop Delta (paired t-test, p < 0.05), the model has no detectable structural awareness for that family and the partner-specificity test is moot. Report as "no structural awareness detected" and exclude that model-family pair from confirmatory analysis. This prerequisite is a screening gate applied per model per family, not a confirmatory hypothesis — it is outside the Bonferroni-corrected family.

### Within-stem derangement null

The null tests whether PS > 0 could arise from positional, compositional, or backbone-propagation confounds rather than genuine partner resolution.

For each stem with K eligible pairs:

1. The K "right-side" positions (j_1, ..., j_K) are permuted via derangement (no position assigned to its true partner).
2. Using the same Delta profiles already computed (no model re-runs), relabel which position counts as "partner" and which as "neighbor" according to the derangement.
3. Compute PS on the relabeled pairs.
4. Repeat 1000 derangements per stem.
5. Stems with K < 3 eligible pairs cannot produce a meaningful derangement null (K=1 has no derangement, K=2 has exactly one). These stems are excluded from the null computation. Only stems with K >= 3 eligible pairs participate in the derangement null. Eligible pairs from excluded stems still contribute to the real PS metric but have no per-stem null comparison.
6. The per-stem null distribution is the set of 1000 mean-PS values.
7. The per-family null 95th percentile is the threshold.

This null preserves: nucleotide identity at every position, stem membership, and base-pair type distribution within the stem. Sequence distance from the mutation site is approximately but not exactly preserved (deranged partners are in the same stem but may be at slightly different positions). It destroys only the specific partner assignment.

Max-over-layers is applied to the real metric. For the null, each derangement's PS is computed at the same best layer as the real metric (not independently max'd across layers per derangement), to avoid inflating the null. As a sensitivity analysis, the independently-max'd null (each derangement selects its own best layer) is always reported alongside the primary null. If the two nulls disagree on which families exceed threshold, both counts are reported.

### Secondary null: loop mutation control

Mutate an unpaired (loop) position to its complement. Compute the same Delta profile. Check whether any stem partner shows PS > 0 relative to its stem neighbors. If loop mutations produce positive PS, the metric is picking up non-specific perturbation propagation, not pair disruption. Report as a sanity check, not a formal hypothesis test.

---

## Hypotheses

### Confirmatory (Bonferroni-corrected, alpha = 0.05/3 = 0.0167)

All confirmatory hypotheses exclude quarantined families. N = 32 non-quarantined families passing the stem-filtering criteria.

**H1: At least one model has mean PS > 0 across N families, exceeding the within-stem derangement null.**

- Metric: Per-family PS (max-over-layers), compared to derangement null.
- Decision criterion: H1 is a single confirmatory decision requiring both conditions; the two conditions are not separately corrected. (a) The number of families where PS exceeds the per-family null 95th percentile is >= ceil(4 * 0.05 * 32) = 7. (b) The mean PS across families is > 0 (one-sample Wilcoxon signed-rank, p < 0.0167). Both must hold for at least one model.
- Outcome if no model passes: current RNA foundation models do not resolve specific base-pair partners at detectable levels.

**H2: PS is higher for RNA-pretrained models than DNA-pretrained models.**

- RNA-pretrained: RNA-FM, RiNALMo, UTR-LM, ERNIE-RNA, SpliceBERT.
- DNA-pretrained: NT v2, HyenaDNA, Caduceus, Evo, DNABERT-2.
- Metric: Rank-biserial correlation of per-model mean PS between groups, with Mann-Whitney U p-value.
- Decision criterion: Rank-biserial > 0.5 (large effect). At 5 vs 5, significance is near-binary; effect size is primary.
- Tokenization note: NT v2 (6-mer) and DNABERT-2 (BPE) use multi-nucleotide tokens. Single-nucleotide swaps change sub-token content, potentially diluting PS. 3/5 DNA models use non-character tokenization vs 0/5 RNA models. If H2 passes, E4 reports the character-tokenized-only comparison.

**H3: Partner specificity is resolved at single-nucleotide precision.**

- For eligible pairs in stem interiors (>= 2 pairs from each stem end), compute the fraction where Delta_partner > max(Delta_{j-1}, Delta_{j+1}). The >= 2 filter guarantees that j-1 and j+1 are always within the stem for every pair counted, so the candidate set is always exactly {j, j-1, j+1}. Under no partner specificity, this fraction = 1/3 (partner is max among three adjacent positions by chance).
- Decision criterion: One-sample binomial test against 1/3, p < 0.0167. If the fraction exceeds 0.5 for at least one model, that model resolves partners at single-nucleotide precision.

### Exploratory

**E1: Layer localization.** Report the distribution of peak-PS layers across models and families. Test whether median peak layer > 50% of network depth (Wilcoxon, p < 0.05). If PS peaks at layer 0, the signal is in tokenization/positional encoding, not learned representations.

**E2: GC vs AU pair-type interaction.** Report PS separately for GC and AU pairs. GC pairs have three hydrogen bonds; AU pairs have two. If PS differs by pair type, report the interaction.

**E3: Per-position asymmetry.** For each pair (i, j), compute PS_forward (swap i, measure at j) and PS_reverse (swap j, measure at i). Report mean |PS_forward - PS_reverse| and whether asymmetry correlates with 5'/3' strand assignment.

**E4: Character-tokenized-only comparison.** If H2 passes, repeat the RNA-vs-DNA comparison restricted to character-tokenized models (5 RNA vs 3 DNA: HyenaDNA, Caduceus, Evo).

**E5: Loop mutation sanity check.** Report PS values from loop mutations. If PS > 0 for loop mutations, flag the metric as potentially confounded by non-specific propagation.

---

## Kill criteria

- If no model achieves mean PS > 0.001 (cosine distance difference units — PS is Delta_partner minus Delta_adj, not an absolute distance) across non-quarantined families, report as "no partner specificity detected" and do not test against the null.
- If fewer than 15 canonical WC pairs are available for a family, exclude it.
- If fewer than 5 eligible interior pairs (after stem-length and terminal-pair filtering) remain for a family, exclude it.
- If the positive control (stem > loop) fails for a model, do not test that model for partner specificity.

---

## Known limitations

**Covariation vs pairing.** A model trained on natural RNA sequences may have learned sequence covariations (MSA-like patterns) rather than geometric base pairing. Complement-swapping position i creates a covariation violation that concentrates at j because j's identity is statistically predictable from i's in the training distribution. A positive PS result is consistent with both "the model encodes pairing geometry" and "the model encodes sequence covariation." Distinguishing these would require synthetic sequences with known structure but no evolutionary covariation. This is noted as a limitation, not tested.

**Tokenization confound.** Multi-nucleotide tokenizers (NT v2, DNABERT-2) may dilute per-position perturbation signals. Results for these models are reported but interpreted with this caveat.

**Isostericity.** Only canonical WC pairs (AU, UA, CG, GC) are tested. Wobble pairs (GU) are excluded because complement-swapping a wobble pair changes base-pairing geometry.

---

## Quarantined pilot data

tRNA_Phe_yeast and tRNA_Ala_human were used during metric development. Values under a weaker metric variant (coupling ratio: partner perturbation vs unpaired-position perturbation, not within-stem comparison) were observed before this preregistration:

| Family | Model | Coupling ratio | Best layer | Date |
|--------|-------|---------------|------------|------|
| tRNA_Phe_yeast | RNA-FM | 1.40 | 12 | 2026-07-13 |
| tRNA_Ala_human | RNA-FM | 1.05 | 12 | 2026-07-13 |

These families are excluded from confirmatory analyses. Their PS values under the preregistered metric will be computed and reported as exploratory pilot observations. If either family would have changed the outcome of a confirmatory test, this is disclosed.

**No PS values have been computed for any family.** The pilot observations above used the coupling ratio (partner perturbation vs unpaired-position perturbation), a different and weaker metric. The perturbation specificity metric defined in this document (partner vs adjacent-stem comparison) has never been computed on any data. PS is fully a priori across all families, including the quarantined ones.

---

## Changelog

- **v1** (SHA `e9f2149`): Defined CR = d_dest/d_comp (destructive vs compensatory perturbation at mutated position). Confounded: swapping two tokens produces mechanically more perturbation than swapping one, independent of structure.
- **v2** (this document): Replaced with perturbation specificity (PS), which compares partner perturbation to adjacent-stem perturbation after a single swap. Controls for stem-vs-loop confound by comparing within the same stem. Null uses derangement of partner labels within stems rather than pair-type-stratified shuffling across the full sequence.
