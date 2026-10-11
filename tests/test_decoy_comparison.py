"""The decoy comparison must separate knowing the partner from knowing its nucleotide.

One hand-built family of two eight-pair stems, with partner-side nucleotides
chosen so that every filter of `matched_decoys` has a candidate it must reject:
another nucleotide, a shallow pair, a stem neighbour of the partner, another
stem, and a position farther from the mutation than the partner is.
"""

import pytest

from decoy_comparison import (
    complementarity_profile,
    matched_decoys,
    oracle_profile,
    partner_beats_decoy,
    precision_at_every_layer,
)
from stem_pairs import get_eligible_pairs, is_deep_interior

# Partner-side nucleotide of each pair, outermost first.
STEM_A = "CCUCCUCC"   # pairs (k, 60 - k)
STEM_B = "CUCCCUCC"   # pairs (100 + k, 160 - k)
COMPLEMENT = {"C": "G", "U": "A"}


def family():
    first = [(k, 60 - k) for k in range(8)]
    second = [(100 + k, 160 - k) for k in range(8)]
    sequence = ["A"] * 161
    for stem, partners in ((first, STEM_A), (second, STEM_B)):
        for (i, j), nucleotide in zip(stem, partners):
            sequence[i], sequence[j] = COMPLEMENT[nucleotide], nucleotide
    sequence = "".join(sequence)
    return sequence, get_eligible_pairs(sequence, [first, second])


def index_of(pairs, i):
    return next(k for k, pair in enumerate(pairs) if pair["i"] == i)


@pytest.mark.parametrize("i, scope, decoy_set, expected", [
    # Pair (3, 57), partner C. 56 is its stem neighbour; 59 and 54 are shallow;
    # 58 and 55 carry U. Only the other stem supplies decoys, all farther away.
    (3, "stem", "all", []),
    (3, "family", "all", [156, 157, 158]),
    (3, "family", "close", []),
    # Pair (103, 157), partner C. 158 and 156 are its stem neighbours; the
    # decoys in the first stem are nearer the mutation than the partner is.
    (103, "stem", "all", []),
    (103, "family", "all", [56, 57]),
    (103, "family", "close", [56, 57]),
    # Pair (2, 58), partner U. 55 is in its own stem and nearer; 155 is not.
    (2, "stem", "all", [55]),
    (2, "stem", "close", [55]),
    (2, "family", "all", [55, 155]),
    (2, "family", "close", [55]),
])
def test_matched_decoys_are_exactly_the_expected_positions(i, scope, decoy_set, expected):
    sequence, pairs = family()
    assert matched_decoys(sequence, pairs, index_of(pairs, i), scope, decoy_set) == expected


def test_unknown_scope_or_decoy_set_is_refused():
    sequence, pairs = family()
    with pytest.raises(ValueError):
        matched_decoys(sequence, pairs, 0, scope="panel")
    with pytest.raises(ValueError):
        matched_decoys(sequence, pairs, 0, decoy_set="near")


def only(pairs, profiles_by_i):
    """Profiles for the named pairs; every other pair has no deltas at all."""
    by_index = {index_of(pairs, i): profile for i, profile in profiles_by_i.items()}
    return [by_index.get(k, {}) for k in range(len(pairs))]


def test_each_pair_weighs_the_same_whatever_its_number_of_decoys():
    sequence, pairs = family()
    # (3, 57) beats all three of its decoys; (2, 58) loses to both of its two.
    wins = {57: 1.0, 156: 0.0, 157: 0.0, 158: 0.0}
    losses = {58: 0.0, 55: 1.0, 155: 1.0}
    profiles = only(pairs, {3: {0: wins}, 2: {0: losses}})
    result = partner_beats_decoy(sequence, pairs, profiles, 1, "family", "all")
    assert result["pairs"] == [2]
    assert result["comparisons"] == [5]
    # Weighting comparisons equally would give 3/5.
    assert result["rate"] == [pytest.approx(0.5)]


def test_a_tie_counts_one_half_and_is_counted():
    sequence, pairs = family()
    mixed = {57: 0.5, 156: 0.5, 157: 0.2, 158: 0.9}
    profiles = only(pairs, {3: {0: mixed}})
    result = partner_beats_decoy(sequence, pairs, profiles, 1, "family", "all")
    assert result["rate"] == [pytest.approx((0.5 + 1.0 + 0.0) / 3)]
    assert result["ties"] == [1]


def test_a_decoy_without_a_delta_is_skipped_at_that_layer_only():
    sequence, pairs = family()
    layer0 = {57: 0.5, 156: 0.5, 157: 0.2}            # 158 absent
    layer1 = {57: 1.0, 156: 0.0, 157: 0.0, 158: 0.0}
    profiles = only(pairs, {3: {0: layer0, 1: layer1}})
    result = partner_beats_decoy(sequence, pairs, profiles, 2, "family", "all")
    assert result["comparisons"] == [2, 3]
    assert result["rate"] == [pytest.approx(0.75), pytest.approx(1.0)]


def test_no_decoy_gives_no_rate():
    stem = [(k, 40 - k) for k in range(5)]
    sequence = ["A"] * 41
    for i, j in stem:
        sequence[i], sequence[j] = "G", "C"
    sequence = "".join(sequence)
    pairs = get_eligible_pairs(sequence, [stem])
    assert sum(is_deep_interior(pair) for pair in pairs) == 1
    profiles = [{0: oracle_profile(sequence, pair)} for pair in pairs]
    result = partner_beats_decoy(sequence, pairs, profiles, 1, "family", "all")
    assert result["rate"] == [None]
    assert result["pairs"] == [0]


@pytest.mark.parametrize("scope", ["stem", "family"])
@pytest.mark.parametrize("decoy_set", ["all", "close"])
def test_oracle_scores_one_and_complementarity_scores_half(scope, decoy_set):
    sequence, pairs = family()
    oracle = [{0: oracle_profile(sequence, pair)} for pair in pairs]
    blind = [{0: complementarity_profile(sequence, pair)} for pair in pairs]
    best = partner_beats_decoy(sequence, pairs, oracle, 1, scope, decoy_set)
    half = partner_beats_decoy(sequence, pairs, blind, 1, scope, decoy_set)
    assert best["pairs"][0] > 0
    assert best["rate"] == [1.0]
    assert half["rate"] == [0.5]
    assert half["ties"] == half["comparisons"]


def test_complementarity_alone_passes_rung3_where_neighbours_differ():
    stem = [(k, 100 - k) for k in range(12)]
    sequence = ["A"] * 101
    for rank, (i, j) in enumerate(stem):
        sequence[i], sequence[j] = ("G", "C") if rank % 2 == 0 else ("A", "U")
    sequence = "".join(sequence)
    pairs = get_eligible_pairs(sequence, [stem])
    blind = [{0: complementarity_profile(sequence, pair)} for pair in pairs]
    # Every partner's stem neighbours carry the other nucleotide, so the Rung 3
    # comparison is won on every pair with no knowledge of the pairing.
    assert precision_at_every_layer(pairs, blind, 1) == [1.0]
    assert partner_beats_decoy(sequence, pairs, blind, 1, "stem", "all")["rate"] == [0.5]
