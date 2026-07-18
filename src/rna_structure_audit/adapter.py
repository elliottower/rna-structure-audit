"""Base adapter for RNA/DNA foundation models.

To evaluate a model, subclass ModelAdapter and implement the three required methods.
See adapters/ for examples.
"""

from abc import ABC, abstractmethod

import torch


class ModelAdapter(ABC):
    """Thin wrapper around any RNA/DNA embedding model.

    Subclass this and implement load(), tokenize(), and
    get_all_layer_embeddings(). Then pass an instance to evaluate().
    """

    name: str
    d_model: int
    n_layers: int

    @abstractmethod
    def load(self) -> None:
        """Load model weights into memory."""
        ...

    @abstractmethod
    def tokenize(self, sequence: str) -> torch.Tensor:
        """Tokenize an RNA sequence (ACGU) and return input_ids tensor."""
        ...

    @abstractmethod
    def get_all_layer_embeddings(
        self, tokens: torch.Tensor
    ) -> list[torch.Tensor]:
        """Return per-layer embeddings as a list of (seq_len, d_model) tensors.

        Must return one tensor per layer (including the embedding layer).
        Embeddings should be at nucleotide resolution where possible.
        Strip special tokens (CLS, EOS, etc.) before returning.
        """
        ...

    def get_final_embeddings(self, tokens: torch.Tensor) -> torch.Tensor:
        """Return final-layer embeddings. Override for efficiency."""
        return self.get_all_layer_embeddings(tokens)[-1]
