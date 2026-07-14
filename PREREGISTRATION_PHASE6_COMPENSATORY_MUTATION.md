# Preregistration: Phase 6 — Compensatory Double Mutation Test

**Title:** Does pairwise coupling survive the composition null? Compensatory vs destructive mutations across ten RNA foundation models

**Author:** Elliot Tower (elliot@elliottower.ai)

**Date:** 2026-07-13

**Commit SHA:** `36da8a1` (causal-rna repo)

**Epistemic status:** This is a preregistration for a new experiment extending Phases 1-5. The compensatory mutation test was designed after observing that no model exceeded the composition-controlled threshold in Phases 1-5 using single-position complement swaps. The hypothesis that compensatory mutations should produce smaller perturbations than destructive mutations is a prediction from RNA secondary structure theory — a structure-aware model should "know" that swapping both partners of a base pair preserves the pairing, while swapping only one partner breaks it. This test was not available in earlier phases and addresses a limitation of the single-swap metric: single-swap sensitivity measures response to breaking one base pair, but does not test whether the model encodes the pairwise coupling between partners.

---

## Motivation

Phases 1-5 measured structure awareness via single-position complement swap: for each position, swap A↔U or C↔G and measure the cosine distance between wild-type and mutant embeddings. The ratio of mean cosine distance at stem positions to mean cosine distance at loop positions (the "mutation sensitivity ratio") was evaluated against a nucleotide-stratified permutation null.

This metric has a limitation: it tests whether stem positions respond differently to mutations than loop positions, but it does not test whether the model encodes the *coupling* between base-pair partners. A model could show elevated stem sensitivity simply because stem nucleotides are in a different compositional context (GC-enriched), without encoding anything about which positions are paired.

The compensatory mutation test addresses this directly. For each base pair (i, j), we define:

- **Destructive mutation:** Swap only position i (e.g., C→G at position i, leaving position j as the original G). This breaks the C-G pair.
- **Compensatory mutation:** Swap both positions i and j (e.g., C→G at i and G→C at j). This preserves base pairing — the new G-C pair is isosteric with the original C-G pair.

A structure-aware model should produce a **larger** embedding perturbation under destructive mutation (structure broken) than under compensatory mutation (structure preserved). A model that only encodes nucleotide identity or local composition should show **similar** perturbations under both conditions, since both involve the same number and type of nucleotide substitutions.

---

## Experiment design

### RNA families and models

Same 52 Rfam families and 10 models as Phases 2-4, reduced to 44 families after excluding those with fewer than 15 canonical WC pairs (see Kill Criteria). For each family, the dot-bracket annotation defines base pairs. Only canonical Watson-Crick pairs (A-U, U-A, C-G, G-C) at positions where both partners are within the model's context window are included.

### Metrics

**Compensatory ratio (CR).** For each base pair (i, j):

1. Compute wild-type embedding: E_wt = model(sequence)
2. Compute destructive embedding: E_dest = model(sequence with position i swapped)
3. Compute compensatory embedding: E_comp = model(sequence with positions i AND j swapped)
4. d_dest = cosine_distance(E_wt[i], E_dest[i])
5. d_comp = cosine_distance(E_wt[i], E_comp[i])

The compensatory ratio CR = mean(d_dest) / mean(d_comp) across all base pairs in a family. CR > 1 means destructive mutations perturb more than compensatory mutations — evidence for pairwise coupling. CR = 1 means the model does not distinguish structure-preserving from structure-breaking mutations.

We also compute the compensatory ratio at position j (the partner position) and report the average.

We compute CR at every layer and report the maximum (max-over-layers), consistent with Phases 1-5.

**Compensatory null.** To control for composition, we apply a base-pair-type-stratified permutation null. The null must absorb the dominant confound identified in Phases 1-2: GC enrichment in stems. Since compensatory mutations operate on pairs, the confound acts at the pair-type level (GC pairs vs AU pairs produce inherently different perturbation magnitudes). The null therefore stratifies by base-pair type:

For each permutation:
1. Partition all base pairs by type: GC pairs (C-G and G-C) and AU pairs (A-U and U-A).
2. Within each type, shuffle which left-side positions pair with which right-side positions. A GC pair is always reassigned to another GC pair; an AU pair to another AU pair. This preserves: (a) the number of GC vs AU pairs, (b) the nucleotide identity at every position, (c) the stem/loop status of every position. It destroys only the specific pairing topology — which position is coupled to which.
3. Compute CR on the shuffled pairing at every layer; take the max-over-layers value.
4. Repeat 100 times; the 95th percentile of these 100 max-over-layers values defines the null threshold.

This ensures the null matches the real data's base-pair-type composition, preventing the confound that collapsed first-order survivors in Phase 2. The max-over-layers selection in the null matches the real metric (the Phase 1→2 bug where the null was computed per-layer independently while the metric used max-over-layers is avoided by explicit max-over-layers in both).

### Per-position perturbation symmetry

Additionally, for each base pair (i, j), compute:
- d_comp_at_i = cosine_distance(E_wt[i], E_comp[i])
- d_comp_at_j = cosine_distance(E_wt[j], E_comp[j])

Under compensatory mutation, both i and j were swapped but the pair is preserved. If the model encodes pairing, perturbation at the mutated positions should be smaller than under destructive mutation but not zero (the nucleotide identity changed).

---

## Hypotheses

### Confirmatory (Bonferroni-corrected, alpha = 0.05/3 = 0.0167)

**H1: At least one model has mean CR > 1.0 across 44 families (one-sample Wilcoxon signed-rank on per-family CR values, p < 0.0167).**

If any model encodes pairwise coupling, compensatory mutations should produce systematically smaller perturbations than destructive mutations. This is the primary test.

- **Metric:** Per-family CR (max-over-layers), tested against 1.0.
- **Decision criterion:** Wilcoxon signed-rank test of per-family CR - 1.0, p < 0.0167 for at least one model.
- **Outcome if no model passes:** Pairwise coupling is absent from all tested models at detectable levels.

**H2: RiNALMo CR exceeds the compensatory null in >= 9/44 families.**

RiNALMo showed the strongest directional signal in Phase 3 (49/52 sign test, 17/52 exceeding the single-swap null). If any model encodes pairwise coupling beyond composition, RiNALMo is the most likely candidate.

- **Metric:** Number of families where CR exceeds the 95th percentile of the compensatory null.
- **Decision criterion:** >= 9/44 families exceed the null. Rationale: the 95th-percentile null produces a ~5% per-family false-positive rate under the null hypothesis. At 44 families, the expected false-positive count is ~2.2; 9/44 (20.5%) is ~4× the expected false-positive rate, comparable to the observed/expected ratio for Phase 2 single-swap survivors (28/52 for NT v2 vs ~2.6 expected = ~11×, but 1/52 for RNA-FM vs ~2.6 expected = ~0.4×). The 4× threshold is conservative enough to distinguish signal from noise while remaining achievable if coupling is present at the strength suggested by RiNALMo's Phase 3 directional result.
- **Outcome if < 9:** Pairwise coupling not detected for RiNALMo beyond composition control.

**H3: CR is higher for RNA-pretrained models than DNA-pretrained models (effect-size-primary).**

RNA-pretrained models (RNA-FM, RiNALMo, UTR-LM, ERNIE-RNA, SpliceBERT) should show stronger pairwise coupling than DNA-pretrained models (NT v2, HyenaDNA, Caduceus, Evo, DNABERT-2), since RNA secondary structure is domain-relevant.

- **Metric:** Rank-biserial correlation (effect size) of per-model mean CR between pretraining-domain groups, with Mann-Whitney U p-value reported for completeness.
- **Decision criterion:** Rank-biserial correlation > 0.5 (large effect). The p-value is reported but this hypothesis is effect-size-primary: at n = 5 vs n = 5, Mann-Whitney U can reach p < 0.0167 only under near-complete rank separation (minimum achievable two-sided p at 5v5 is ~0.008), so the significance test is near-binary and its failure is uninformative about effect presence. A large effect size with a non-significant p is reported as "consistent with H3 but underpowered"; a small effect size is reported as "H3 not supported."
- **Tokenization note:** NT v2 (6-mer) and DNABERT-2 (BPE) tokenize multiple nucleotides per token. A single-nucleotide complement swap changes sub-token content for these models, diluting the per-position signal relative to character-tokenized models. CR is computed for all 10 models, but the RNA-vs-DNA comparison is confounded by tokenization (3/5 DNA models use non-character tokenization vs 0/5 RNA models). This confound is noted; if H3 passes, a secondary analysis restricted to character-tokenized models (RNA-FM, RiNALMo, UTR-LM, ERNIE-RNA, SpliceBERT vs HyenaDNA, Caduceus; 5 vs 2) is reported as exploratory.

### Exploratory

**E1: Directional sign test for compensatory ratio (trained vs untrained).**

Report the fraction of families where trained CR > untrained CR. This is demoted from confirmatory (cf. Phase 2 H10) because CR is a ratio of two perturbation distances that, for untrained models with random weights, hovers near 1.0 with high variance. A bare sign test on a noisy ratio centered at 1.0 can pass at 75% from directionless noise breaking slightly one way. Without pilot CR data establishing that the trained-vs-untrained CR distribution is directional (as Phase 1's 12/12 established for single-swap sensitivity before Phase 2 promoted it to confirmatory), a sign test here risks a hollow "PASS." Report the sign count, the per-family delta distribution, and the median trained-minus-untrained CR for interpretation.

**E2: Interaction between CR and base-pair type.**

G-C pairs have three hydrogen bonds; A-U pairs have two. Compensatory mutations at G-C pairs may produce different CR patterns than at A-U pairs. Report CR separately for G-C and A-U base pairs.

**E3: Layer localization of compensatory ratio.**

Report the distribution of best-layer indices for CR. Compare to the best-layer distribution for single-swap mutation sensitivity. If CR peaks at different layers than single-swap sensitivity, this suggests pairwise coupling and single-position sensitivity are encoded at different depths.

**E4: Per-position asymmetry.**

For each base pair (i, j) under compensatory mutation, compare d_comp_at_i vs d_comp_at_j. If the model treats the two partners asymmetrically (e.g., 5' partner always perturbed more), this reveals directional biases in the representation.

**E5: Character-tokenized-only RNA vs DNA comparison.**

If H3 shows a pretraining-domain effect, report the same comparison restricted to character-tokenized models only (RNA-FM, RiNALMo, UTR-LM, ERNIE-RNA, SpliceBERT vs HyenaDNA, Caduceus; 5 vs 2) to disentangle tokenization from pretraining domain.

---

## Kill criteria

- If no model achieves mean CR > 1.05 (i.e., all are within 5% of 1.0), report as "no pairwise coupling detected" and do not test against the null.
- If fewer than 15 canonical Watson-Crick base pairs (as returned by `filter_wc_pairs()`, not paired positions) are available for any family, exclude that family from the analysis. Pair counts are estimated from Phase 2 `n_stem` values as `n_pairs = n_stem / 2` (each pair contributes two stem positions). At this threshold, 44 of 52 families are retained (8 excluded: Histone 3' [~5 pairs], hammerhead [~8], IRE stem-loop [~10], CRISPR leader [~12], glycine riboswitch [~13], hatchet ribozyme [~13], preQ1 riboswitch [~13], pistol ribozyme [~14]). The exact pair count will be confirmed at runtime from dot-bracket parsing; if `filter_wc_pairs()` returns a different count (e.g., because some annotated pairs are non-canonical wobbles), the exclusion list will update accordingly and the change will be logged. The 15-pair minimum ensures enough per-pair observations for a stable mean while retaining 85% of the family set.

---

## Analysis code

The experiment reuses the existing multi-model adapter framework (multi_model_audit.py ModelAdapter classes) and the nucleotide-stratified null infrastructure from Phases 1-5. New code:

1. `compute_compensatory_ratio(model, sequence, dot_bracket, layer)` — computes d_dest and d_comp for each base pair at one layer.
2. `compensatory_null(sequence, dot_bracket, n_permutations=100)` — pair-shuffled null.
3. `run_phase6(models, families)` — orchestrates across models and families, saves results to `data/gpu_results/phase6_compensatory/`.

Results format (JSON):
```json
{
  "model": "ernierna",
  "timestamp": "...",
  "compensatory_trained": {
    "per_rna": {
      "tRNA_Phe_yeast": {
        "best_cr": 1.15,
        "best_layer": 8,
        "comp_null_95th": 1.08,
        "exceeds_comp_null": true,
        "n_pairs": 21,
        "cr_gc_pairs": 1.18,
        "cr_au_pairs": 1.10,
        "d_dest_mean": 0.045,
        "d_comp_mean": 0.039
      }
    },
    "mean_best_cr": 1.12,
    "families_exceeding_comp_null": 15,
    "families_total": 44
  },
  "compensatory_untrained": { ... }
}
```

Families with fewer than 15 base pairs are recorded in the JSON with `{"skipped": true, "reason": "< 15 pairs"}` but excluded from all summary statistics and hypothesis tests. The `families_total` field counts only non-skipped families.

**Isostericity note.** Compensatory swaps (A-U → U-A, C-G → G-C) preserve canonical Watson-Crick geometry. Only canonical WC pairs from the dot-bracket annotation are included; G-U wobble pairs (which some Rfam annotations include as paired) are excluded by `filter_wc_pairs()` because a complement swap at a wobble pair would change the base-pairing geometry, violating the isostericity assumption.
```

---

## What this adds to the paper

If H1 passes (any model shows CR > 1.0 systematically): the paper gains its first positive result — models encode pairwise coupling even if they don't pass the single-swap composition null. The paper becomes "models encode coupling but not position-level structure," which is more interesting than "everything is null."

If H1 fails: the null result is substantially stronger because it now spans two independent metrics (single-swap sensitivity and pairwise coupling). The limitation that "complement swap doesn't test coupling" is closed.

Either outcome strengthens the paper.
