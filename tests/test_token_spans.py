"""A token's embedding must land on the nucleotides that token actually covers.

The stub tokenizers here reproduce the two schemes in the panel: fixed 6-mers
with the remainder falling to single nucleotides, and a merge rule whose output
depends on the sequence, so a substitution can retokenize the region around it.
Every embedding is a function of its token string alone, which makes the
distance at a mutated position computable in the test without the pipeline: if
a row is placed on the wrong nucleotide the comparison is between two different
token strings and the number moves.
"""

import numpy as np
import pytest
import torch
from scipy.spatial.distance import cosine

from family_checkpoint import FamilyCheckpoint
from phase6_compensatory_mutation import (
    COMPLEMENT,
    complement_swap,
    compute_delta_profiles,
    get_eligible_pairs,
    parse_stems,
)
from phases_1_to_5 import run_mutation_sensitivity
from token_spans import (
    TokenizationMismatch,
    content_bounds,
    content_spans,
    nucleotide_rows,
)

NUC = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3, "N": 0}
MAX_TOKEN = 8


def embed(token: str) -> np.ndarray:
    """Injective for tokens up to MAX_TOKEN long, so two tokens never collide."""
    v = np.zeros(4 * MAX_TOKEN)
    for j, char in enumerate(token):
        v[4 * j + NUC[char]] += 1.0
    return v


class StubTokenizer:
    """Ids are assigned on first sight, so a token always decodes to itself.

    `closing=False` is NT v2's shape: an opening token and nothing after the
    sequence.
    """

    def __init__(self, split, closing=True):
        self.split = split
        self.closing = closing
        self.vocab = {"[CLS]": 0, "[SEP]": 1}
        self.inverse = {0: "[CLS]", 1: "[SEP]"}

    def _id(self, token):
        if token not in self.vocab:
            self.vocab[token] = len(self.vocab)
            self.inverse[self.vocab[token]] = token
        return self.vocab[token]

    def __call__(self, sequence, **kwargs):
        pieces = (["[CLS]"] + self.split(sequence)
                  + (["[SEP]"] if self.closing else []))
        ids = [self._id(p) for p in pieces]
        if kwargs.get("return_tensors") == "pt":
            return {"input_ids": torch.tensor([ids])}
        return {"input_ids": ids}

    def convert_ids_to_tokens(self, ids):
        return [self.inverse[int(i)] for i in ids]


def six_mers(sequence):
    """Nucleotide Transformer's scheme: 6-mers, then the remainder one by one."""
    whole = len(sequence) // 6
    return ([sequence[i * 6:(i + 1) * 6] for i in range(whole)]
            + list(sequence[whole * 6:]))


def merge_aa(sequence):
    """A merge rule that a substitution can change, as byte-pair encoding is."""
    out, i = [], 0
    while i < len(sequence):
        if sequence[i:i + 2] == "AA":
            out.append("AA")
            i += 2
        else:
            out.append(sequence[i])
            i += 1
    return out


class StubAdapter:
    name = "stub"
    token_resolution = "bpe"

    def __init__(self, split, closing=True):
        self.tokenizer = StubTokenizer(split, closing)

    def tokenize(self, sequence):
        return self.tokenizer(sequence, return_tensors="pt")["input_ids"]

    def get_all_layer_embeddings(self, tokens):
        start, stop = content_bounds(self.tokenizer, tokens[0].tolist())
        pieces = self.tokenizer.convert_ids_to_tokens(tokens[0])[start:stop]
        return [torch.tensor(np.stack([embed(p) for p in pieces]))]


def containing(split, sequence, position):
    """The token string covering `position`, computed without the pipeline."""
    start = 0
    for piece in split(sequence):
        if start <= position < start + len(piece):
            return piece
        start += len(piece)
    raise AssertionError("position past the end of the sequence")


def a_family(sequence):
    stem = len(sequence) // 4
    dot_bracket = ("(" * stem + "." * (len(sequence) - 2 * stem) + ")" * stem)
    return {"name": "stub", "sequence": sequence, "dot_bracket": dot_bracket}


def a_sequence(n, period=7):
    """Not a multiple of six, and no run of A long enough to be all one token."""
    return "".join("ACGTCG"[(i * 5 + i // period) % 6] for i in range(n))


def a_hairpin(arm=20):
    """One stem long enough to clear the registered Rung 3 filters.

    The arm is built out of AA so a complement swap inside it breaks a merge
    and retokenizes everything downstream, which is the case a nucleotide index
    into a token array gets wrong rather than merely off by a constant.
    """
    left = ("AAGC" * ((arm + 3) // 4))[:arm]
    right = "".join(COMPLEMENT[c] for c in reversed(left))
    return {"name": "stub",
            "sequence": left + "CUUG" + right,
            "dot_bracket": "(" * arm + "...." + ")" * arm}


def test_spans_partition_the_sequence_and_carry_their_own_token():
    sequence = a_sequence(100)
    tokenizer = StubTokenizer(six_mers)
    spans = content_spans(tokenizer, sequence)

    assert spans[0][0] == 0 and spans[-1][1] == len(sequence)
    assert all(end == spans[i + 1][0] for i, (_, end) in enumerate(spans[:-1]))
    for (start, end), piece in zip(spans, six_mers(sequence)):
        assert sequence[start:end] == piece


def test_a_tokenization_that_drops_nucleotides_is_an_error():
    tokenizer = StubTokenizer(lambda s: six_mers(s)[:-1])
    with pytest.raises(TokenizationMismatch):
        content_spans(tokenizer, a_sequence(100))


def test_a_tokenizer_that_only_opens_keeps_its_last_token():
    # NT v2's shape. A hardcoded hidden[0, 1:-1, :] drops the final content
    # token here, which is how 92 nucleotides of the panel ended up with no
    # embedding row; see docs/OPEN_DEFECTS.md, D14.
    sequence = a_sequence(100)
    tokenizer = StubTokenizer(six_mers, closing=False)

    assert content_bounds(tokenizer, tokenizer(sequence)["input_ids"]) == (
        1, len(six_mers(sequence)) + 1)
    assert content_spans(tokenizer, sequence)[-1][1] == len(sequence)


def test_bounds_and_spans_agree_on_how_many_rows_carry_sequence():
    for closing in (True, False):
        for n in (96, 101):
            sequence = a_sequence(n)
            tokenizer = StubTokenizer(six_mers, closing=closing)
            start, stop = content_bounds(tokenizer,
                                         tokenizer(sequence)["input_ids"])
            assert stop - start == len(content_spans(tokenizer, sequence))


def test_a_special_token_inside_the_sequence_is_an_error():
    # A contiguous block of rows is what a slice can express; an interior
    # special token means no slice is right, and silently taking one would put
    # every later row on the wrong nucleotides.
    tokenizer = StubTokenizer(six_mers)
    original = tokenizer.convert_ids_to_tokens
    tokenizer.convert_ids_to_tokens = lambda ids: (
        original(ids)[:3] + ["[MASK]"] + original(ids)[3:])
    with pytest.raises(TokenizationMismatch):
        content_spans(tokenizer, a_sequence(100))


def test_every_nucleotide_reads_the_token_that_covers_it():
    for n in (96, 97, 101, 103):
        sequence = a_sequence(n)
        spans = content_spans(StubTokenizer(six_mers), sequence)
        rows = nucleotide_rows(spans, n)
        pieces = six_mers(sequence)
        for position in range(n):
            assert pieces[rows[position]] == containing(six_mers, sequence, position)


def test_the_remainder_of_a_six_mer_sequence_does_not_collapse_onto_one_row():
    n = 101
    spans = content_spans(StubTokenizer(six_mers), a_sequence(n))
    rows = nucleotide_rows(spans, n)
    tail = rows[96:]
    assert len(set(tail.tolist())) == 5, (
        "the five leftover nucleotides are five separate tokens; `i // 6` gives "
        "them all row 16")


@pytest.mark.parametrize("split,name", [(six_mers, "nt"), (merge_aa, "dnabert2")])
def test_the_stored_distance_is_between_the_two_containing_tokens(split, name, tmp_path):
    sequence = a_sequence(101)
    path = tmp_path / "positions.json"
    run_mutation_sensitivity(StubAdapter(split), name, [a_family(sequence)],
                             device="cpu", n_permutations=10,
                             positions=FamilyCheckpoint(path))

    import json
    stored = json.loads(path.read_text())["done"]["stub"]
    distances = stored["layer_distances"][0]
    valid = stored["layer_valid"][0]
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    checked = 0
    for position, nucleotide in enumerate(sequence):
        if valid[position] != "1":
            continue
        mutant = (sequence[:position] + complement[nucleotide]
                  + sequence[position + 1:])
        expected = cosine(embed(containing(split, sequence, position)),
                          embed(containing(split, mutant, position)))
        assert distances[position] == pytest.approx(expected, abs=1e-6)
        checked += 1
    assert checked == len(sequence), "every position should be scored"


def test_a_retokenized_mutation_is_recorded_as_a_changed_span(tmp_path):
    # ...AA... at 40; mutating the second A to T splits the merged token.
    sequence = a_sequence(101)[:40] + "AA" + a_sequence(101)[42:]
    path = tmp_path / "positions.json"
    run_mutation_sensitivity(StubAdapter(merge_aa), "dnabert2", [a_family(sequence)],
                             device="cpu", n_permutations=10,
                             positions=FamilyCheckpoint(path))

    import json
    same_span = json.loads(path.read_text())["done"]["stub"]["same_span"]
    assert same_span[40] == "0" and same_span[41] == "0", (
        "both nucleotides of the merged token change span when either mutates")
    assert same_span.count("0") < len(sequence) / 2, (
        "only the merged positions should move")


def test_the_arithmetic_the_spans_replaced_lands_on_the_wrong_token():
    # 101 nucleotides tokenize to 16 six-mers and 5 singles, 21 rows. `i // 6`
    # clamps positions 96 to 100 onto row 16, so four of the five leftovers read
    # a token that does not contain them. This is the 1.2% that
    # scripts/audit_token_alignment.py measures on the real panel.
    n = 101
    sequence = a_sequence(n)
    spans = content_spans(StubTokenizer(six_mers), sequence)
    rows = nucleotide_rows(spans, n)
    guessed = [min(i // 6, len(spans) - 1) for i in range(n)]
    assert sum(g != r for g, r in zip(guessed, rows)) == 4

    spans = content_spans(StubTokenizer(merge_aa), sequence)
    rows = nucleotide_rows(spans, n)
    guessed = [min(int(i * len(spans) / n), len(spans) - 1) for i in range(n)]
    misplaced = sum(g != r for g, r in zip(guessed, rows))
    assert misplaced > n / 3, (
        f"the byte-pair guess should be wrong for most of the sequence, not "
        f"{misplaced}/{n}")


def test_rung3_reads_the_token_holding_each_partner_position():
    # `compute_delta_profiles` indexed the hidden state by nucleotide for every
    # model, so on a subword tokenizer every partner position in a stem fell
    # past the end of the array and the family was recorded as having no valid
    # PS values; see docs/OPEN_DEFECTS.md, D13.
    family = a_hairpin()
    sequence = family["sequence"]
    stems = parse_stems(family["dot_bracket"], sequence)
    eligible = get_eligible_pairs(sequence, stems)
    assert len(eligible) >= 5

    profiles, _, _, n_layers = compute_delta_profiles(
        StubAdapter(merge_aa), sequence, eligible, stems)
    assert n_layers == 1

    discriminating = 0
    for idx, pair in enumerate(eligible):
        deltas = profiles[idx][0]
        mutant = complement_swap(sequence, pair["i"])
        rows_wt, rows_mut = merge_aa(sequence), merge_aa(mutant)
        for key in ("j", "j_prev", "j_next"):
            position = pair[key]
            assert position in deltas, (
                f"{key} of pair {idx} has no delta; the row it needs is past "
                "the end of the token array")
            expected = cosine(embed(containing(merge_aa, sequence, position)),
                              embed(containing(merge_aa, mutant, position)))
            assert deltas[position] == pytest.approx(expected, abs=1e-6)

            # What indexing the token array by nucleotide would have read here,
            # where it lands on a row at all. The pipeline did that for every
            # model, so the two numbers agreeing everywhere would leave the
            # check above unable to tell the repair from the defect.
            if position < min(len(rows_wt), len(rows_mut)):
                naive = cosine(embed(rows_wt[position]), embed(rows_mut[position]))
                discriminating += abs(naive - expected) > 1e-6

    assert discriminating, ("a nucleotide index reads the same rows as a token "
                            "index on this family, so the test proves nothing")


def test_rung3_leaves_a_nucleotide_tokenizer_indexed_by_nucleotide():
    # The eight character-level models must be untouched by the row map, or
    # repairing NT and DNABERT-2 would move numbers that were already right.
    family = a_hairpin()
    sequence = family["sequence"]
    stems = parse_stems(family["dot_bracket"], sequence)
    eligible = get_eligible_pairs(sequence, stems)

    adapter = StubAdapter(list)
    adapter.token_resolution = "nucleotide"
    profiles, _, _, _ = compute_delta_profiles(adapter, sequence, eligible, stems)

    for idx, pair in enumerate(eligible):
        mutant = complement_swap(sequence, pair["i"])
        for position in (pair["j"], pair["j_prev"], pair["j_next"]):
            expected = cosine(embed(sequence[position]), embed(mutant[position]))
            assert profiles[idx][0][position] == pytest.approx(expected, abs=1e-6)


# ── An unknown token holds sequence ─────────────────────────────────────────
#
# DNABERT-2's byte-pair vocabulary cannot spell `N`, which Rfam seed alignments
# carry, so it emits `[UNK]` for that nucleotide. Reading the row block by
# asking which tokens are nucleotide strings classified that row as a bracket
# and split the block; reconstructing spans by joining token strings gave the
# row the width of `[UNK]` and shifted every span after it (D23).


class UnknownEmittingTokenizer:
    """The parts of a fast tokenizer the helpers use.

    `vocab` is the set of substrings it can spell; anything else becomes the
    unknown token, as a byte-pair vocabulary missing `N` behaves.
    """

    unk_token = "[UNK]"
    is_fast = True

    def __init__(self, vocab, width=2):
        self.vocab = vocab
        self.width = width

    def _pieces(self, sequence):
        out, at = [], 0
        while at < len(sequence):
            chunk = sequence[at:at + self.width]
            if chunk not in self.vocab:
                out.append((self.unk_token, at, at + 1))
                at += 1
                continue
            out.append((chunk, at, at + len(chunk)))
            at += len(chunk)
        return out

    def __call__(self, sequence, return_offsets_mapping=False):
        pieces = self._pieces(sequence)
        self._last = pieces
        out = {"input_ids": [0] + list(range(2, 2 + len(pieces))) + [1]}
        if return_offsets_mapping:
            out["offset_mapping"] = ([(0, 0)] + [(s, e) for _, s, e in pieces]
                                     + [(0, 0)])
        return out

    def convert_ids_to_tokens(self, ids):
        return ["[CLS]"] + [p for p, _, _ in self._last] + ["[SEP]"]


UNKNOWN_VOCAB = {"AC", "GT", "TA", "CG", "GG", "AA", "TT", "CC", "AG", "CA",
                 "GA", "TC"}


def test_unknown_token_does_not_split_the_row_block():
    tok = UnknownEmittingTokenizer(UNKNOWN_VOCAB)
    sequence = "ACGTNACGT"
    ids = tok(sequence)["input_ids"]
    start, stop = content_bounds(tok, ids)
    assert (start, stop) == (1, 1 + len(tok._last))


def test_spans_cover_the_sequence_when_a_token_is_unknown():
    tok = UnknownEmittingTokenizer(UNKNOWN_VOCAB)
    sequence = "ACGTNACGT"
    spans = content_spans(tok, sequence)
    assert spans[0][0] == 0 and spans[-1][1] == len(sequence)
    assert all(spans[i][1] == spans[i + 1][0] for i in range(len(spans) - 1))
    assert sum(end - start for start, end in spans) == len(sequence)


def test_the_unknown_token_covers_one_nucleotide_not_the_width_of_its_name():
    tok = UnknownEmittingTokenizer(UNKNOWN_VOCAB)
    spans = content_spans(tok, "ACGTNACGT")
    assert [s for s in spans if s[0] <= 4 < s[1]] == [(4, 5)]
