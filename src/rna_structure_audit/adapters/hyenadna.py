"""HyenaDNA adapter (5.4M, Hyena SSM, DNA-pretrained, character-level)."""

import torch

from rna_structure_audit.adapter import ModelAdapter


class HyenaDNAAdapter(ModelAdapter):
    name = "HyenaDNA-tiny"
    d_model = 128
    n_layers = 2

    MODEL_ID = "LongSafari/hyenadna-tiny-1k-seqlen-hf"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from transformers import AutoModelForCausalLM, AutoTokenizer

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model = AutoModelForCausalLM.from_pretrained(
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
