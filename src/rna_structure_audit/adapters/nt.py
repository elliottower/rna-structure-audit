"""Nucleotide Transformer v2 adapter (56M, ESM+GLU, DNA-pretrained, 6-mer)."""

import torch

from rna_structure_audit.adapter import ModelAdapter


class NTAdapter(ModelAdapter):
    name = "Nucleotide-Transformer-v2"
    d_model = 512
    n_layers = 12

    MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model = AutoModelForMaskedLM.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model.eval()

    def tokenize(self, sequence: str) -> torch.Tensor:
        seq_dna = sequence.replace("U", "T")
        encoded = self.tokenizer(seq_dna, return_tensors="pt")
        return encoded["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        layers = []
        for hs in out.hidden_states:
            emb = hs[0, 1:-1, :]
            layers.append(emb)
        return layers
