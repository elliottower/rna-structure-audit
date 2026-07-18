# rna-structure-audit

Three-rung benchmark for evaluating whether RNA/DNA foundation models encode genuine secondary structure or composition shortcuts.

## Install

```bash
pip install rna-structure-audit
```

## Usage

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

## Citation

Tower, E. (2026). Composition Confounds Inflate Apparent Structure Awareness in RNA Foundation Models.
