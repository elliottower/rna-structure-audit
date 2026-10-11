"""RNA-FM through its own tokenizer and model class.

Supersedes `multi_model_audit.RNAFMAdapter`, which is left as it was run. That
adapter built token ids by hand (A, C, G, U as 5, 6, 7, 8, bracketed by 2 and 3)
and loaded the weights into `transformers.BertModel`. The checkpoint it pins,
`multimolecule/rnafm`, has A, C, G, U at 6, 7, 8, 9 with `<cls>` at 1 and
`<eos>` at 2, so every nucleotide reached the model as a different token, and
`BertModel` is not the architecture the weights were trained in.

This adapter follows the four other multimolecule adapters in
`multi_model_audit.py`: ids come from the tokenizer and the forward pass from
the library's model class, at the same pinned revision.
"""

import torch

from multi_model_audit import ModelAdapter
from token_spans import content_bounds


class RNAFMAdapterV2(ModelAdapter):
    name = "RNA-FM"
    d_model = 640
    n_layers = 12
    token_resolution = "nucleotide"

    HF_REPO = "multimolecule/rnafm"
    HF_REVISION = "7d6e73ad3b48e042b378f9a788a56ccb4d573a27"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from multimolecule import RnaFmModel, RnaTokenizer

        self.tokenizer = RnaTokenizer.from_pretrained(
            self.HF_REPO, revision=self.HF_REVISION)
        self.model = RnaFmModel.from_pretrained(
            self.HF_REPO, revision=self.HF_REVISION, attn_implementation="eager",
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        enc = self.tokenizer(seq, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        start, stop = content_bounds(self.tokenizer, tokens[0].tolist())
        return [hs[0, start:stop, :] for hs in out.hidden_states]
