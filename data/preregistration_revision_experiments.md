# Preregistration — Revision Experiments: Multi-Sequence RiNALMo and Transversion Completion

**Project:** A Graded Evaluation of RNA Structure Awareness in Foundation Models
**Author:** Elliot Tower
**Date frozen:** 2026-07-17
**Context:** Peer review response — addressing reviewer concerns from Genome Biology (R1), NAR (R2), and ML/benchmarking (R3) reviewers.
**Paper version at freeze:** paper_v9.tex

---

## 1. Background and motivation

Three independent reviewers raised two empirical concerns that require GPU experiments:

1. **Within-family variance validated on one model only** (R2/NAR, major concern #2): The 5.7% within-family CV was measured only for ERNIE-RNA. R2 requested extending the multi-sequence evaluation to RiNALMo to validate that the low variance is not model-specific.

2. **Transversion control incomplete** (R2/NAR, minor concern #5): The original transversion run completed for 6 of 10 models. NT v2, Evo, Caduceus, and RNA-FM were missing due to infrastructure constraints. All reviewers have access to the 6-model results; completing the panel strengthens the control.

Both experiments use the same analysis pipelines as the original runs with no methodological changes.

---

## 2. Experiment A: Multi-sequence RiNALMo

### Design

- **Model:** RiNALMo (650M, RNA, character tokenization)
- **Data:** 3-5 alternative Rfam seed sequences per family, 51 families (253 sequences total). Same multi-sequence files used for ERNIE-RNA.
- **Pipeline:** `run_mutation_sensitivity()` from `phases_1_to_5.py`, identical to the ERNIE-RNA multi-sequence run.
- **Metric:** Within-family coefficient of variation (CV) of the mutation sensitivity ratio.
- **GPU:** 1x A100 (RiNALMo requires >24GB VRAM)
- **Script:** `scripts/modal_multiseq_rinalmo.py`

### Hypotheses (frozen before running)

- **H-R1 (primary):** RiNALMo within-family mean CV < 0.10 (i.e., within-family sequence choice contributes less than 10% noise to the mutation sensitivity ratio). ERNIE-RNA achieved CV = 0.057; we expect RiNALMo to be comparable.
- **H-R2 (secondary):** RiNALMo within-family CV is not significantly different from ERNIE-RNA's CV (two-sample t-test on per-family CVs, alpha = 0.05, two-sided). Failure would indicate model-specific sensitivity to sequence choice.

### Analysis plan

1. Compute per-family CV of the mutation sensitivity ratio across the 3-5 sequences.
2. Report overall mean CV, median CV, and distribution across families.
3. Compare to ERNIE-RNA results (mean CV = 0.057, median CV = 0.055).
4. If H-R1 passes: update paper_v9.tex Section 5.1 (Limitations, "Single representative per family") to include RiNALMo CV alongside ERNIE-RNA.
5. If H-R1 fails: report the higher CV and discuss implications for the single-sequence design.

### Outcomes that would change the paper

- CV > 0.15 would require adding per-family error bars to RiNALMo results and potentially re-evaluating the single-sequence design.
- CV < 0.10 strengthens the existing claim with no structural changes.

---

## 3. Experiment B: Transversion completion (4 models)

### Design

- **Models:** RNA-FM (99M, RNA), NT v2 (56M, DNA), Evo (7B, DNA), Caduceus (14M, DNA)
- **Data:** Same 52 Rfam families as all other analyses.
- **Pipeline:** `run_mutation_sensitivity()` with `COMPLEMENT` dict patched to transversion mutations: `{A: C, U: C, C: A, G: U, T: C}`. Identical to the original 6-model transversion run.
- **Metric:** Mean mutation sensitivity ratio under transversion mutations, compared to the WC-complement ratio from the original analysis.
- **GPU:** 4x A10G (one per model, parallel)
- **Script:** `scripts/modal_transversion_remaining4.py`

### Hypotheses (frozen before running)

- **H-R3 (primary):** For each of the 4 models, the transversion ratio is within 10% of the WC-complement ratio. Specifically:
  - RNA-FM: transversion ratio within 10% of 1.749
  - NT v2: transversion ratio within 10% of 1.103
  - Evo: transversion ratio within 10% of 1.408
  - Caduceus: transversion ratio within 10% of 1.211
- **H-R4 (exploratory):** If any model shows a transversion/WC divergence > 10%, the direction will be the same as HyenaDNA (transversion < WC), consistent with complement-symmetric features.

### Analysis plan

1. For each model, compute the transversion mean ratio across all families.
2. Compare to the WC-complement ratio from Table 1 (paper_v9.tex).
3. Compute percent change: `(transversion - WC) / WC * 100`.
4. If all 4 models pass H-R3: update paper_v9.tex transversion paragraph to report all 10 models instead of 6, remove the "infrastructure constraints" caveat.
5. If any model shows > 10% divergence: report the result and discuss the mechanism.

### Outcomes that would change the paper

- All within 10%: replace "six of ten models" with "all ten models" in the transversion paragraph. Strengthens the control.
- One or more > 10% divergence: add the divergent model(s) alongside HyenaDNA's 7.2% drop and discuss whether complement-symmetric features are architecture-specific.

---

## 4. Infrastructure

- **Platform:** Modal (GPU cloud, --detach mode, 86400s timeout)
- **Images:** multimol (torch 2.6.0, transformers 5.14.1, multimolecule 0.2.0) for RNA-FM and RiNALMo; transformers (torch 2.6.0, transformers 5.14.1) for NT v2, Evo, Caduceus
- **Total GPU time estimate:** ~2-4 hours (RiNALMo A100: ~2h; 4x A10G transversion: ~1h each, parallel)
- **Results saved to:** Modal volumes, downloaded to `results/` directory

## 5. Commitment

These hypotheses and analysis plans are frozen as of this document's creation. The SHA of this file will be recorded before launching any GPU job. No changes to the hypotheses, thresholds, or analysis plan will be made after results are available.
