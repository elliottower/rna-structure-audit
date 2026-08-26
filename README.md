# rna-structure-audit

Three-rung benchmark for whether RNA and DNA foundation models encode base-pairing structure or nucleotide composition.

> **State of the panel.** Four of the 52 curated families carried wrong
> annotations and five were withdrawn on review, leaving 47 that the benchmark
> scores. The whole panel was then re-run under pinned numerics, and
> `results/repaired_panel_v3/` is that run: one commit per cell, each file
> stamped with the commit, panel hash and library versions that produced it.
> `docs/OPEN_DEFECTS.md` registers every defect found during the audit, each
> with a detector.
>
> The withdrawn families still ship, so a result computed before the review can
> be reproduced, and `load_families(include_withdrawn=True)` returns them.
>
> The `as-submitted` branch reconstructs the state the submitted manuscript's
> numbers were computed from. `main` is the corrected line.

[![Open Quickstart in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/elliottower/rna-structure-audit/blob/main/notebooks/quickstart.ipynb) [![Open BYOM Tutorial in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/elliottower/rna-structure-audit/blob/main/notebooks/bring_your_own_model.ipynb)

## Install

```bash
pip install rna-structure-audit
```

For pre-built adapters (ERNIE-RNA, RiNALMo, RNA-FM, and 7 more):

```bash
pip install rna-structure-audit[all-models]
```

## Quick start

```python
from rna_structure_audit.adapters import ERNIERNAAdapter
from rna_structure_audit.evaluate import evaluate

results = evaluate(ERNIERNAAdapter(), device="cuda")
print(results["report"]["grade"])  # A
```

## Pre-built adapters

| Adapter | Model | Params | Grade |
|---------|-------|--------|-------|
| `ERNIERNAAdapter` | ERNIE-RNA | 86M | A |
| `RiNALMoAdapter` | RiNALMo | 650M | A |
| `NTAdapter` | Nucleotide Transformer v2 | 56M | B |
| `HyenaDNAAdapter` | HyenaDNA | 0.45M | B |
| `EvoAdapter` | Evo | 7B | B |
| `CaduceusAdapter` | Caduceus | 7.7M | C |
| `SpliceBERTAdapter` | SpliceBERT | 19M | D |
| `UTRLMAdapter` | UTR-LM | 1.2M | D |
| `RNAFMAdapter` | RNA-FM | 99M | D |
| `DNABERT2Adapter` | DNABERT-2 | 117M | D |

Parameter counts are the loaded checkpoint's, except Evo and DNABERT-2 where
they are the authors' designation. Grades are `_grade()` applied to the 47-family
panel; a B means families survive the dinucleotide null, which is a weaker claim
than partner specificity and is not evidence that a model resolves pairing.

## Bring your own model

Write an adapter for your model:

```python
from rna_structure_audit.adapter import ModelAdapter
import torch

class MyModelAdapter(ModelAdapter):
    name = "my-model"
    d_model = 640
    n_layers = 12

    def load(self):
        # Load your model
        ...

    def tokenize(self, sequence: str) -> torch.Tensor:
        # Return input_ids tensor
        ...

    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        # Return list of (seq_len, d_model) tensors, one per layer
        ...
```

Then run from the command line:

```bash
rna-structure-audit --adapter my_adapter.py --device cuda -o results.json
```

Or from Python:

```python
from rna_structure_audit.evaluate import evaluate

adapter = MyModelAdapter()
results = evaluate(adapter, device="cuda")
print(results["report"]["grade"])  # A, B, C, or D
```

See the [Bring Your Own Model tutorial](notebooks/bring_your_own_model.ipynb) for a full walkthrough.

## Grading

| Grade | Meaning |
|-------|---------|
| A | Partner-specific: encodes which position pairs with which |
| B | Structure-aware beyond composition: survives dinucleotide controls |
| C | Composition-sensitive: stem/loop signal absorbed by nucleotide null |
| D | No detectable structure signal |

## The three rungs

1. **Mutation sensitivity** — Do stems respond differently than loops to complement mutations? Controlled by nucleotide-stratified permutation null.
2. **Dinucleotide null** — Of families passing Rung 1, how many survive when the null is stratified by dinucleotide context?
3. **Partner specificity** — When position i is mutated, is perturbation at its base-pairing partner j greater than at j's stem-adjacent neighbors? Controlled by within-stem derangement null.

## The paper and its verification

The manuscript is the highest-numbered `paper/rna-structure-audit_vNN.tex`.
Every table it prints is written from `results/` by a generator in `scripts/`
rather than entered by hand, so the paper can be checked against the data
without re-running a model:

```bash
for s in verify_paper_rung12_figures verify_paper_phase6_figures \
         verify_artifact_regenerates audit_attention_rho check_freeze_order; do
    uv run --no-project --with numpy --with scipy --with tqdm --python 3.12 \
        python "scripts/$s.py"
done
```

| Script | What it checks |
|--------|----------------|
| `verify_paper_rung12_figures.py` | every mean ratio, interval, exceedance count and retention rate in Tables 1 and 2, against `results/bootstrap_cis.json` |
| `verify_paper_phase6_figures.py` | Table 5 and the Rung 3 prose figures, recomputed from the per-family values with the two quarantined families excluded |
| `verify_artifact_regenerates.py` | the deposited artifact against the per-family inputs, leaf by leaf |
| `audit_attention_rho.py` | each attention correlation in Table 1 against the run that produced it |
| `check_freeze_order.py` | that no run is dated before the preregistration it is reported against |

The three scripts that read the manuscript resolve the highest-numbered
`paper/paper_vN.tex`, so they cannot silently keep checking a superseded
version; the two `verify_paper_*` scripts also accept an explicit path.

`verify_artifact_regenerates.py` re-runs the bootstrap over the per-family
result files and compares 269 leaf values against the deposited
`results/bootstrap_cis.json`. `docs/provenance.md` records what was measured,
what changed, and which file settled it. The preregistrations are at the
repository root and in `preregistration/`; `docs/sha_map.md` maps the commit
ids frozen documents quote to their present-day equivalents.

## Citation

Tower, E. (2026). A Graded Evaluation of RNA Structure Awareness in Foundation
Models: From Stem-Loop Discrimination to Partner Specificity. Zenodo.
https://doi.org/10.5281/zenodo.21362717
