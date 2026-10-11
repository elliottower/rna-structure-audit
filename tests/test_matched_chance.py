"""The derangement chance rate must be measured on the pairs precision is read on.

Delta profiles are built by hand, so the chance rate is known without a model.
An "oracle" profile moves only the true partner: under any derangement the
assigned partner is false and has delta zero, so it can never beat its
neighbours and the chance rate is exactly zero. A "noise" profile gives every
position an independent delta, so a false partner beats its two neighbours about
a third of the time.
"""

import numpy as np
import pytest

from phase6_compensatory_mutation_v3 import (
    derangement_null,
    get_eligible_pairs,
    is_deep_interior,
    precision_at_every_layer,
)

STEM_LENGTH = 12
SPAN = 100


def one_stem():
    """A single stem of G-C pairs (k, SPAN - k), and its eligible pairs."""
    stem = [(k, SPAN - k) for k in range(STEM_LENGTH)]
    sequence = ["A"] * (SPAN + 1)
    for i, j in stem:
        sequence[i], sequence[j] = "G", "C"
    return get_eligible_pairs("".join(sequence), [stem])


def oracle(pair):
    return {pos: (1.0 if pos == pair["j"] else 0.0) for pos in range(SPAN + 1)}


def noise(rng):
    return {pos: float(rng.random()) for pos in range(SPAN + 1)}


def test_oracle_has_zero_chance_and_full_precision_at_every_layer():
    pairs = one_stem()
    profiles = [{0: oracle(p), 1: oracle(p)} for p in pairs]
    null = derangement_null(pairs, profiles, n_layers=2, best_layer=0,
                            rng=np.random.default_rng(), n_derangements=200)
    assert null["chance_per_layer"] == [0.0, 0.0]
    assert null["h3_chance_fraction_matched"] == 0.0
    assert precision_at_every_layer(pairs, profiles, n_layers=2) == [1.0, 1.0]


def test_matched_chance_ignores_pairs_precision_does_not_use():
    rng = np.random.default_rng()
    pairs = one_stem()
    assert any(not is_deep_interior(p) for p in pairs)
    # Noise only where the mutated pair is outside the precision set.
    profiles = [{0: oracle(p) if is_deep_interior(p) else noise(rng)}
                for p in pairs]
    null = derangement_null(pairs, profiles, n_layers=1, best_layer=0,
                            rng=rng, n_derangements=500)
    assert null["h3_chance_fraction"] > 0.0
    assert null["h3_chance_fraction_matched"] == 0.0


def test_chance_is_read_at_the_layer_asked_for():
    rng = np.random.default_rng()
    pairs = one_stem()
    profiles = [{0: oracle(p), 1: noise(rng)} for p in pairs]
    at_oracle = derangement_null(pairs, profiles, n_layers=2, best_layer=0,
                                 rng=rng, n_derangements=500)
    at_noise = derangement_null(pairs, profiles, n_layers=2, best_layer=1,
                                rng=rng, n_derangements=500)
    assert at_oracle["h3_chance_fraction_matched"] == 0.0
    assert at_noise["h3_chance_fraction_matched"] == pytest.approx(
        at_noise["chance_per_layer"][1])
    assert at_noise["chance_per_layer"][0] == 0.0
    assert at_noise["chance_per_layer"][1] > 0.1


def two_stems():
    """Two stems far apart, with every eligible pair of both."""
    first = [(k, 100 - k) for k in range(STEM_LENGTH)]
    second = [(200 + k, 300 - k) for k in range(STEM_LENGTH)]
    sequence = ["A"] * 301
    for i, j in first + second:
        sequence[i], sequence[j] = "G", "C"
    return get_eligible_pairs("".join(sequence), [first, second]), first, second


def test_matched_denominator_is_derangements_times_deep_pairs():
    pairs = one_stem()
    deep = sum(is_deep_interior(p) for p in pairs)
    profiles = [{0: oracle(p), 1: oracle(p)} for p in pairs]
    null = derangement_null(pairs, profiles, n_layers=2, best_layer=0,
                            rng=np.random.default_rng(0), n_derangements=200)
    assert null["chance_per_layer_n"] == [200 * deep, 200 * deep]


def test_false_partners_never_come_from_another_stem():
    pairs, first, second = two_stems()
    j_of = {0: [j for _i, j in first], 1: [j for _i, j in second]}
    profiles = []
    for pair in pairs:
        other = j_of[1 - pair["stem_idx"]]
        # Zero inside the pair's own stem; alternating peaks in the other stem,
        # so a false partner drawn from there would beat its neighbours.
        deltas = {pos: 0.0 for pos in range(301)}
        for rank, pos in enumerate(other):
            deltas[pos] = 1.0 if rank % 2 else 0.5
        profiles.append({0: deltas})
    null = derangement_null(pairs, profiles, n_layers=1, best_layer=0,
                            rng=np.random.default_rng(0), n_derangements=300)
    assert null["h3_chance_fraction"] == 0.0
    assert null["chance_per_layer"] == [0.0]


def test_keys_of_the_superseded_script_are_reproduced_exactly():
    from phase6_compensatory_mutation import derangement_null as superseded

    pairs = one_stem()
    profile_rng = np.random.default_rng(1)
    profiles = [{0: noise(profile_rng), 1: noise(profile_rng)} for _ in pairs]
    old = superseded(pairs, profiles, 2, 1, np.random.default_rng(7), 300)
    new = derangement_null(pairs, profiles, 2, 1, np.random.default_rng(7), 300)
    for key, value in old.items():
        assert new[key] == value, key


def test_moved_precision_function_equals_the_superseded_one():
    from phase6_compensatory_mutation import precision_at_every_layer as superseded

    pairs = one_stem()
    rng = np.random.default_rng(3)
    profiles = [{0: noise(rng), 1: noise(rng)} for _ in pairs]
    assert precision_at_every_layer(pairs, profiles, 2) == superseded(pairs, profiles, 2)
