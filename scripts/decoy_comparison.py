"""Does a mutation move its partner more than it moves a matched decoy?

Rung 3 compares the partner of a mutated position with the partner's two stem
neighbours. A representation that responds wherever the complementary
nucleotide occurs, without encoding which occurrence is the partner, wins that
comparison whenever the neighbours carry other nucleotides. Here the partner is
compared with decoys: partner-side positions of other deep-interior pairs that
carry the partner's nucleotide. A response shared by every such position scores
one half.

No model and no tensor library: every function takes sequences, pairs and
delta profiles.
"""

from stem_pairs import is_deep_interior

SCOPES = ("stem", "family")
DECOY_SETS = ("all", "close")


def matched_decoys(sequence, eligible_pairs, index, scope="family", decoy_set="all"):
    """Decoy positions for the pair at `index`, in ascending order.

    A decoy for the pair (i, j) is the partner-side position j' of another
    eligible pair such that
      - that pair is deep-interior, as the pairs precision is read on are;
      - sequence[j'] == sequence[j], so i is as complementary to j' as to j;
      - j' is not j - 1 or j + 1, which PS already compares the partner with;
      - under scope "stem", it belongs to the stem of (i, j);
      - under decoy_set "close", |i - j'| <= |i - j|: the decoy is no farther
        from the mutated position than the partner is.
    """
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {SCOPES}, got {scope!r}")
    if decoy_set not in DECOY_SETS:
        raise ValueError(f"decoy_set must be one of {DECOY_SETS}, got {decoy_set!r}")
    pair = eligible_pairs[index]
    reach = abs(pair["i"] - pair["j"])
    decoys = []
    for other_index, other in enumerate(eligible_pairs):
        if other_index == index or not is_deep_interior(other):
            continue
        if scope == "stem" and other["stem_idx"] != pair["stem_idx"]:
            continue
        if sequence[other["j"]] != sequence[pair["j"]]:
            continue
        if abs(other["j"] - pair["j"]) <= 1:
            continue
        if decoy_set == "close" and abs(pair["i"] - other["j"]) > reach:
            continue
        decoys.append(other["j"])
    return sorted(decoys)


def partner_beats_decoy(sequence, eligible_pairs, delta_profiles, n_layers,
                        scope="family", decoy_set="all"):
    """Per layer, how often the partner's delta exceeds a matched decoy's.

    For each deep-interior pair the score is the mean over its decoys of
    1 (partner larger), 0.5 (tie) or 0 (decoy larger); the family rate is the
    mean of those pair scores, so a pair with many decoys weighs the same as a
    pair with one. A decoy or partner whose delta is absent at a layer is
    skipped at that layer.

    Returns, per layer: `rate`, `pairs` (pairs contributing), `comparisons`
    and `ties`.
    """
    decoys = {}
    for index, pair in enumerate(eligible_pairs):
        if is_deep_interior(pair):
            found = matched_decoys(sequence, eligible_pairs, index, scope, decoy_set)
            if found:
                decoys[index] = found

    rate, pairs, comparisons, ties = [], [], [], []
    for layer in range(n_layers):
        pair_scores = []
        n_comparisons = n_ties = 0
        for index, found in decoys.items():
            deltas = delta_profiles[index].get(layer, {})
            partner = eligible_pairs[index]["j"]
            if partner not in deltas:
                continue
            scores = []
            for decoy in found:
                if decoy not in deltas:
                    continue
                if deltas[partner] > deltas[decoy]:
                    scores.append(1.0)
                elif deltas[partner] == deltas[decoy]:
                    scores.append(0.5)
                    n_ties += 1
                else:
                    scores.append(0.0)
            if scores:
                pair_scores.append(sum(scores) / len(scores))
                n_comparisons += len(scores)
        rate.append(sum(pair_scores) / len(pair_scores) if pair_scores else None)
        pairs.append(len(pair_scores))
        comparisons.append(n_comparisons)
        ties.append(n_ties)
    return {"scope": scope, "decoy_set": decoy_set, "rate": rate, "pairs": pairs,
            "comparisons": comparisons, "ties": ties}


def precision_at_every_layer(eligible_pairs, delta_profiles, n_layers):
    """Rung 3 per-pair precision at each layer, over deep-interior pairs.

    The fraction of pairs whose partner delta exceeds both stem neighbours'.
    Identical to the function of that name in `phase6_compensatory_mutation.py`.
    """
    out = []
    for layer in range(n_layers):
        hits = total = 0
        for idx, pair in enumerate(eligible_pairs):
            if not is_deep_interior(pair):
                continue
            deltas = delta_profiles[idx].get(layer, {})
            if not all(k in deltas for k in (pair["j"], pair["j_prev"], pair["j_next"])):
                continue
            total += 1
            hits += deltas[pair["j"]] > max(deltas[pair["j_prev"]], deltas[pair["j_next"]])
        out.append(hits / total if total else None)
    return out


def oracle_profile(sequence, pair):
    """Deltas of a representation that moves the true partner and nothing else."""
    return {k: (1.0 if k == pair["j"] else 0.0)
            for k in range(len(sequence)) if k != pair["i"]}


def complementarity_profile(sequence, pair):
    """Deltas of a representation that moves every position carrying the
    partner's nucleotide equally, with no knowledge of which is the partner."""
    target = sequence[pair["j"]]
    return {k: (1.0 if sequence[k] == target else 0.0)
            for k in range(len(sequence)) if k != pair["i"]}
