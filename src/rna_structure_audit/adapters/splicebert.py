"""SpliceBERT adapter (19M, BERT, RNA-pretrained, character-level)."""

import torch

from rna_structure_audit.adapter import ModelAdapter


class SpliceBERTAdapter(ModelAdapter):
    name = "SpliceBERT"
    d_model = 512
    n_layers = 6

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from multimolecule import RnaTokenizer, SpliceBertModel

        self.tokenizer = RnaTokenizer.from_pretrained("multimolecule/splicebert")
        self.model = SpliceBertModel.from_pretrained(
            "multimolecule/splicebert", attn_implementation="eager",
        )
        self.model.eval()

    def tokenize(self, sequence: str) -> torch.Tensor:
        enc = self.tokenizer(sequence, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        return [hs[0, 1:-1, :] for hs in out.hidden_states]
