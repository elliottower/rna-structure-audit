"""Which nucleotides a token covers, and which rows carry them, read off the
tokenizer instead of assumed.

Eight of the ten models emit one token per nucleotide, so a row of the hidden
state is a position of the molecule. NT v2 and DNABERT-2 do not, and the closed
forms that stood in for their tokenizers -- `i // 6` and
`i * n_tokens / n_nucleotides` -- misassign 1.2% and 52.8% of nucleotides, and
put 45.9% of DNABERT-2's mutation comparisons on a pair of rows describing
different positions. `scripts/audit_token_alignment.py` measures all three.

Two facts have to come from the same tokenization or they drift:

  bounds  which rows of the hidden state carry sequence, so an adapter strips
          exactly the special tokens its tokenizer emitted
  spans   what each of those rows covers, so a nucleotide can be mapped to the
          row that holds it

Assuming the first is how NT v2 lost a row. Its tokenizer emits `<cls>` and
nothing after the sequence, so the `hidden[0, 1:-1, :]` slice written for a
symmetric wrapper discarded the last content token of all 47 families -- 92
nucleotides of the panel with no embedding row at all. DNABERT-2, which does
emit `[SEP]`, was unaffected. `scripts/audit_special_token_slices.py` reports
both ends for every tokenizer the panel uses.
"""

from __future__ import annotations

import numpy as np

NUCLEOTIDES = set("ACGTUN")


class TokenizationMismatch(ValueError):
    """The token strings do not account for the sequence one-to-one."""


def _bounds(pieces: list[str], unk_token: str | None = None) -> tuple[int, int]:
    """Half-open row range holding the sequence tokens.

    Special tokens are recognized by not being nucleotide strings rather than
    by name, so a tokenizer that brackets the sequence differently is handled
    without a per-model table.

    An unknown token counts as content. A vocabulary that cannot spell part of
    the sequence emits one -- DNABERT-2's byte-pair vocabulary does it for `N`,
    which Rfam seed alignments carry -- and that token still occupies a row
    holding sequence. Treating it as a bracket instead split the row block and
    stopped six of ten models on the same alignment (D23).
    """
    is_content = [bool(p) and (set(p) <= NUCLEOTIDES or p == unk_token)
                  for p in pieces]
    if not any(is_content):
        raise TokenizationMismatch("no token is a nucleotide string")
    start = is_content.index(True)
    stop = len(pieces) - is_content[::-1].index(True)
    if not all(is_content[start:stop]):
        interior = [p for p, c in zip(pieces[start:stop], is_content[start:stop])
                    if not c]
        raise TokenizationMismatch(
            f"special tokens {interior} sit between content tokens; the "
            "embedding rows are not one contiguous block")
    return start, stop


def content_bounds(tokenizer, token_ids) -> tuple[int, int]:
    """Half-open row range of a hidden state that carries sequence.

    Adapters slice their hidden states with this instead of a hardcoded `1:-1`,
    which is right only for a tokenizer that brackets the sequence
    symmetrically. `token_ids` is the flat id list for one sequence.
    """
    return _bounds(tokenizer.convert_ids_to_tokens(list(token_ids)),
                   getattr(tokenizer, "unk_token", None))


def content_spans(tokenizer, sequence: str) -> list[tuple[int, int]]:
    """Half-open nucleotide range of each content token, in row order.

    `sequence` must be the exact string the adapter hands the tokenizer --
    already transcribed to the alphabet that adapter uses.
    """
    unk = getattr(tokenizer, "unk_token", None)
    if getattr(tokenizer, "is_fast", False):
        # The tokenizer's own character offsets, which are the only thing that
        # can place an unknown token: its string is `[UNK]`, so its width is
        # not its length. Reconstructing spans by joining token strings silently
        # loses the nucleotides an unknown token covers, and every span after it
        # shifts.
        encoded = tokenizer(sequence, return_offsets_mapping=True)
        pieces = tokenizer.convert_ids_to_tokens(encoded["input_ids"])
        start, stop = _bounds(pieces, unk)
        spans = [tuple(span) for span in encoded["offset_mapping"][start:stop]]
        covered = sum(end - begin for begin, end in spans)
        if spans[0][0] != 0 or spans[-1][1] != len(sequence) or covered != len(sequence):
            raise TokenizationMismatch(
                f"offsets cover {covered} of {len(sequence)} nucleotides, "
                f"from {spans[0][0]} to {spans[-1][1]}")
        return spans

    ids = tokenizer(sequence)["input_ids"]
    pieces = tokenizer.convert_ids_to_tokens(ids)
    start, stop = _bounds(pieces, unk)
    content = pieces[start:stop]

    joined = "".join(content)
    if joined != sequence:
        raise TokenizationMismatch(
            f"tokens reconstruct {len(joined)} nucleotides, sequence has "
            f"{len(sequence)}")

    spans, at = [], 0
    for piece in content:
        spans.append((at, at + len(piece)))
        at += len(piece)
    return spans


def nucleotide_rows(spans: list[tuple[int, int]], n_nucleotides: int) -> np.ndarray:
    """Row index of the token covering each nucleotide."""
    if not spans or spans[-1][1] != n_nucleotides:
        covered = spans[-1][1] if spans else 0
        raise TokenizationMismatch(
            f"spans cover {covered} nucleotides, asked for {n_nucleotides}")
    rows = np.empty(n_nucleotides, dtype=np.int64)
    for row, (start, end) in enumerate(spans):
        rows[start:end] = row
    return rows
