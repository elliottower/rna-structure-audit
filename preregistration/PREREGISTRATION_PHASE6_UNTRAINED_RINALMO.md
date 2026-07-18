# Preregistration: Phase 6 Untrained Control — RiNALMo

**Title:** Does RiNALMo's partner specificity require learned weights? Randomized-weight control for Phase 6

**Author:** Elliot Tower (elliot@elliottower.ai)

**Commit SHA:** `0aff46c`

**Parent preregistration:** PREREGISTRATION_PHASE6_V2.md (SHA `c19aa59`)

**Pilot data disclosure:** RiNALMo Phase 6 results (trained weights, mean PS = 0.227 over 34 non-skipped families) were observed before this preregistration was written. The trained result motivates this control but does not inform the untrained prediction.

---

## Motivation

Phase 6 showed that RiNALMo (650M, BERT+RoPE+SwiGLU, RNA-pretrained, standard attention) achieves mean PS = 0.227 with 28/34 families exceeding the derangement null. An equivalent control already exists for ERNIE-RNA: randomizing all weights yields PS = 2.5e-8, confirming that ERNIE-RNA's partner specificity is learned rather than architectural.

RiNALMo uses standard multi-head attention with no explicit pairing bias, so there is no architectural mechanism that could produce partner specificity without training. This control confirms that expectation and closes the asymmetry in the paper's evidence.

---

## Experiment design

### Model

RiNALMo (650M parameters, 33 transformer layers, BERT+RoPE+SwiGLU architecture). Loaded via multimolecule (HuggingFace: `multimolecule/rinalmo`).

### Weight randomization

Identical procedure to the ERNIE-RNA untrained control:
- All parameters randomized in-place with `torch.no_grad()`
- Matrices (dim >= 2): `xavier_normal_`
- Vectors (dim < 2): `normal_(std=0.02)`
- Architecture, tokenizer, and embedding dimensions preserved
- Parameter count verified unchanged after randomization

### Evaluation

Run Phase 6 (perturbation specificity) on the same 52 Rfam families using the same protocol as PREREGISTRATION_PHASE6_V2.md:
- Same eligible-pair filtering (stems >= 3 WC pairs, interior pairs only)
- Same PS metric: PS(i,j) = Delta_partner - max(Delta_{j-1}, Delta_{j+1})
- Same layer sweep (all 33 layers), report max-over-layers per family
- Same derangement null (within-stem permutation)

### Predictions

**H_null:** Untrained RiNALMo produces mean PS indistinguishable from zero (|PS| < 0.001).

**Rationale:** Standard attention has no architectural pairing bias. Without learned weights, the model has no mechanism to preferentially perturb paired partners over adjacent stem positions. The ERNIE-RNA untrained control (PS = 2.5e-8) sets the precedent.

### Analysis plan

1. Compute mean PS across all non-skipped families
2. Count families exceeding the derangement null
3. Compare to trained RiNALMo (PS = 0.227) — expect >= 4 orders of magnitude difference
4. Report alongside the ERNIE-RNA untrained control in the paper

---

## Decision rule

If untrained RiNALMo PS < 0.001: confirms partner specificity is learned. Report as confirmatory evidence.

If untrained RiNALMo PS >= 0.001: unexpected. Investigate whether RoPE or SwiGLU introduce a structural bias. Report the finding regardless.
