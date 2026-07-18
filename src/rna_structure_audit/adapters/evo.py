"""Evo adapter (7B, StripedHyena, DNA-pretrained, byte-level). GPU required."""

import torch

from rna_structure_audit.adapter import ModelAdapter


class EvoAdapter(ModelAdapter):
    name = "Evo-1-8k"
    d_model = 4096
    n_layers = 32

    MODEL_ID = "togethercomputer/evo-1-8k-base"

    def __init__(self):
        self.model = None

    def load(self):
        if not torch.cuda.is_available():
            raise RuntimeError("Evo (7B params) requires CUDA GPU")

        from transformers import AutoModelForCausalLM

        self.model = AutoModelForCausalLM.from_pretrained(
            self.MODEL_ID, trust_remote_code=True, torch_dtype=torch.float16,
        )
        self.model.eval()

    def tokenize(self, sequence: str) -> torch.Tensor:
        seq_dna = sequence.replace("U", "T")
        tokens = [ord(c) for c in seq_dna]
        return torch.tensor([tokens])

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        tokens = tokens.to(next(self.model.parameters()).device)
        out = self.model(tokens, output_hidden_states=True)
        if out.hidden_states is not None:
            return [hs[0].float().cpu() for hs in out.hidden_states]
        hidden_states = []
        hooks = []
        backbone = self.model.backbone if hasattr(self.model, "backbone") else self.model.model
        for layer in backbone.blocks:
            h = layer.register_forward_hook(
                lambda mod, inp, out, hs=hidden_states: hs.append(
                    (out[0] if isinstance(out, tuple) else out)[0].float().cpu()
                )
            )
            hooks.append(h)
        self.model(tokens)
        for h in hooks:
            h.remove()
        return hidden_states
