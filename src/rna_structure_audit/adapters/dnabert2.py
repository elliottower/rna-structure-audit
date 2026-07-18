"""DNABERT-2 adapter (117M, BERT, DNA-pretrained, BPE tokenization)."""

import torch

from rna_structure_audit.adapter import ModelAdapter


class DNABERT2Adapter(ModelAdapter):
    name = "DNABERT-2"
    d_model = 768
    n_layers = 12

    MODEL_ID = "zhihan1996/DNABERT-2-117M"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from transformers import AutoModel, AutoTokenizer

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model = AutoModel.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model.eval()

    def tokenize(self, sequence: str) -> torch.Tensor:
        seq_dna = sequence.replace("U", "T")
        enc = self.tokenizer(seq_dna, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        if hasattr(out, "hidden_states") and out.hidden_states:
            return [hs[0, 1:-1, :] for hs in out.hidden_states]
        last_hidden = out[0] if isinstance(out, tuple) else out
        if last_hidden.dim() == 2:
            last_hidden = last_hidden.unsqueeze(0)
        return [last_hidden[0, 1:-1, :]]
