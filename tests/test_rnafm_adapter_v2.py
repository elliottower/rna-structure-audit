"""RNA-FM must receive the ids its own vocabulary assigns to each nucleotide.

Needs `multimolecule` and the pinned tokenizer files, so it is skipped where
the library is not installed.
"""

import pytest

from multi_model_audit import RNAFMAdapter
from rnafm_adapter_v2 import RNAFMAdapterV2

multimolecule = pytest.importorskip("multimolecule")

SEQUENCE = "ACGUUGCA"


@pytest.fixture(scope="module")
def tokenizer():
    return multimolecule.RnaTokenizer.from_pretrained(
        RNAFMAdapterV2.HF_REPO, revision=RNAFMAdapterV2.HF_REVISION)


def test_adapter_ids_are_the_vocabulary_ids(tokenizer):
    adapter = RNAFMAdapterV2()
    adapter.tokenizer = tokenizer
    ids = adapter.tokenize(SEQUENCE)[0].tolist()
    assert ids[1:-1] == tokenizer.convert_tokens_to_ids(list(SEQUENCE))
    assert ids[0] == tokenizer.cls_token_id
    assert ids[-1] == tokenizer.eos_token_id


def test_superseded_adapter_ids_are_not_the_vocabulary_ids(tokenizer):
    vocabulary = dict(zip("ACGU", tokenizer.convert_tokens_to_ids(list("ACGU"))))
    assert RNAFMAdapter.NUC_TO_ID != vocabulary
    # Each hand-built id lands on a different token than the one intended.
    for nucleotide, old_id in RNAFMAdapter.NUC_TO_ID.items():
        assert tokenizer.convert_ids_to_tokens(old_id) != nucleotide


@pytest.fixture(scope="module")
def loaded():
    adapter = RNAFMAdapterV2()
    adapter.load()
    return adapter


def test_forward_returns_one_finite_state_per_layer_per_nucleotide(loaded):
    import torch

    states = loaded.get_all_layer_embeddings(loaded.tokenize(SEQUENCE))
    # The embedding output plus one state per block.
    assert len(states) == loaded.n_layers + 1
    for state in states:
        assert tuple(state.shape) == (len(SEQUENCE), loaded.d_model)
        assert torch.isfinite(state).all()


def test_last_state_is_the_models_own_final_hidden_state(loaded):
    import torch

    tokens = loaded.tokenize(SEQUENCE)
    states = loaded.get_all_layer_embeddings(tokens)
    with torch.no_grad():
        out = loaded.model(tokens, output_hidden_states=True)
    # Content rows only: one opening and one closing token around the sequence.
    assert torch.equal(states[-1], out.last_hidden_state[0, 1:-1, :])


def test_two_forward_passes_agree(loaded):
    import torch

    tokens = loaded.tokenize(SEQUENCE)
    first = loaded.get_all_layer_embeddings(tokens)
    second = loaded.get_all_layer_embeddings(tokens)
    for a, b in zip(first, second):
        assert torch.equal(a, b)


def test_complement_changes_the_embedding_of_the_mutated_position(loaded):
    import torch

    wild = loaded.get_all_layer_embeddings(loaded.tokenize(SEQUENCE))[0]
    mutant = loaded.get_all_layer_embeddings(loaded.tokenize("UCGUUGCA"))[0]
    assert not torch.equal(wild[0], mutant[0])
    # At the embedding layer no other position has seen the change.
    assert torch.equal(wild[1:], mutant[1:])
