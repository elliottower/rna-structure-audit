"""Caduceus adapter (14M, BiMamba, DNA-pretrained, character-level)."""

import torch

from rna_structure_audit.adapter import ModelAdapter


class CaduceusAdapter(ModelAdapter):
    name = "Caduceus"
    d_model = 256
    n_layers = 16

    MODEL_ID = "kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16"

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
        if hasattr(out, "hidden_states") and out.hidden_states:
            return [hs[0] for hs in out.hidden_states]
        return [out.logits[0]]
