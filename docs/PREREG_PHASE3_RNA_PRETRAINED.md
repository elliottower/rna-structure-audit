# Phase 3 Preregistration: RNA-Pretrained Model Replication

## Motivation

Phase 2 found that NT v2 (DNA-pretrained, ESM+GLU) learns attention patterns
correlating with base-pairing contacts (Spearman rho = 0.322 trained vs 0.000
untrained), while RNA-FM (RNA-pretrained, BERT) shows null attention-contact
correlation (0.067 trained vs 0.061 untrained). This finding rests on n=1
RNA-pretrained model. Phase 3 adds two additional RNA-pretrained transformers
to test whether the null attention result generalizes beyond RNA-FM.

## Models

- **RiNALMo-giga** (650M params, BERT+RoPE+SwiGLU, character-level, RNA-pretrained
  on 36M ncRNA sequences). HuggingFace: `multimolecule/rinalmo-giga`. This is the
  critical test: same architecture family as RNA-FM (BERT-style encoder), 6.5x
  larger, trained on a broader RNA corpus.

- **UTR-LM** (1.2M params, 6-layer transformer, character-level, semi-supervised
  on 5' UTRs with structure prediction auxiliary objective). HuggingFace:
  `multimolecule/utrlm-te_el`. Weaker test (tiny model, UTR-specific) but
  interesting because it was explicitly trained with a structure objective.

## Hypotheses

**H12** (RiNALMo attention-contact): RiNALMo trained attention-contact
Spearman rho < 0.10 across 52 families (i.e., null like RNA-FM, not learned
like NT v2). This tests whether BERT-style RNA pretraining in general fails
to learn structure attention, or whether RNA-FM is idiosyncratic.

- PASS (rho < 0.10): The null attention result generalizes across RNA-pretrained
  BERT-style models. The "architecture > pretraining data" claim is robust.
- FAIL (rho >= 0.10): RNA-FM's null attention is idiosyncratic. The story
  becomes "RNA-FM specifically is weak," not "BERT-style RNA pretraining."

**H13** (UTR-LM attention-contact): UTR-LM trained attention-contact Spearman
rho > 0.15 across 52 families. UTR-LM was trained with a structure prediction
objective, so it should show stronger structure awareness than models trained
with MLM alone. If this fails (rho < 0.15), even explicit structure supervision
doesn't produce attention-level structure encoding in small transformers.

**H14** (RiNALMo mutation sensitivity): RiNALMo trained/untrained mean ratio
< 2.0 across 52 families with corrected max-over-layers null. Same threshold
as H6. Tests whether a 6.5x-larger RNA model shows stronger embedding-level
structure sensitivity than RNA-FM (1.77x).

**H15** (RiNALMo sign test): RiNALMo trained exceeds untrained in >= 75% of
52 families. Same threshold as H10.

## Methods

Identical to Phase 2 (preregistered at SHA bd4b3fd):
- Mutation sensitivity ratio under complement swap (A<->U, C<->G)
- Nucleotide-stratified permutation null (100 permutations, max-over-layers)
- Attention-contact Spearman correlation (trained vs untrained)
- 52 Rfam families from Phase 2

Torch version pinned to 2.13.0 for reproducibility.

## Analysis code

Adapters: `lib/adapters/rinalmo.py`, `lib/adapters/utrlm.py`
Runner: `scripts/modal_expanded_rfam.py` (run_rinalmo_expanded, run_utrlm_expanded)

## Preregistration date

2026-07-13
