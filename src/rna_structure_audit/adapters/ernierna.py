"""ERNIE-RNA adapter (86M, BERT, RNA-pretrained, character-level)."""

import torch

from rna_structure_audit.adapter import ModelAdapter


class ERNIERNAAdapter(ModelAdapter):
    name = "ERNIE-RNA"
    d_model = 768
    n_layers = 12

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from multimolecule import RnaTokenizer, ErnieRnaModel

        self.tokenizer = RnaTokenizer.from_pretrained("multimolecule/ernierna")
        self.model = ErnieRnaModel.from_pretrained(
            "multimolecule/ernierna", attn_implementation="eager",
        )
        self.model.eval()

    def tokenize(self, sequence: str) -> torch.Tensor:
        enc = self.tokenizer(sequence, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        return [hs[0, 1:-1, :] for hs in out.hidden_states]
