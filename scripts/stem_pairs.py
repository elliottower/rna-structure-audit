"""Stems, eligible pairs and complementary decoys, from a sequence and its structure.

No model and no tensor library: everything here is a function of the annotated
family alone. The first four functions are the definitions of
`phase6_compensatory_mutation.py`, moved so they can be imported without torch;
`tests/test_stem_pairs.py` holds them equal to that script's on every panel
family.
"""

WC_PAIRS = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}


def parse_dot_bracket(db_string):
    stack = []
    pairs = []
    for idx, char in enumerate(db_string):
        if char == "(":
            stack.append(idx)
        elif char == ")":
            if stack:
                partner = stack.pop()
                pairs.append((partner, idx))
    return sorted(pairs)


def parse_stems(dot_bracket, sequence):
    """Group consecutive stacked WC pairs into stems.

    A stem is a maximal run of pairs (i, j), (i+1, j-1), (i+2, j-2)...
    Only canonical WC pairs are included.
    """
    opens = dot_bracket.count("(")
    closes = dot_bracket.count(")")
    if opens != closes:
        return []
    all_pairs = parse_dot_bracket(dot_bracket)
    wc_pairs = [(i, j) for i, j in all_pairs if (sequence[i], sequence[j]) in WC_PAIRS]
    pair_set = set(wc_pairs)

    used = set()
    stems = []
    for i, j in wc_pairs:
        if (i, j) in used:
            continue
        stem = [(i, j)]
        used.add((i, j))
        ci, cj = i + 1, j - 1
        while ci < cj and (ci, cj) in pair_set and (ci, cj) not in used:
            stem.append((ci, cj))
            used.add((ci, cj))
            ci += 1
            cj -= 1
        stems.append(stem)
    return stems


def get_eligible_pairs(sequence, stems):
    """Return interior pairs from stems with >= 3 WC pairs.

    Each eligible pair includes references to its stem-adjacent j positions.
    Returns list of dicts with keys: i, j, j_prev, j_next, stem_idx, pair_type.
    """
    eligible = []
    for stem_idx, stem in enumerate(stems):
        if len(stem) < 3:
            continue
        for k in range(1, len(stem) - 1):
            i, j = stem[k]
            _, j_prev = stem[k - 1]
            _, j_next = stem[k + 1]
            eligible.append({
                "i": i, "j": j,
                "j_prev": j_prev, "j_next": j_next,
                "stem_idx": stem_idx,
                "pos_in_stem": k,
                "stem_length": len(stem),
                "pair_type": f"{sequence[i]}-{sequence[j]}",
            })
    return eligible


def is_deep_interior(pair):
    """Two or more positions from each stem end: the pairs H3 precision uses."""
    return 2 <= pair["pos_in_stem"] <= pair["stem_length"] - 3


def complementary_decoys(sequence, eligible_pairs, index, scope="stem"):
    """Positions that could pair with the mutated position and are not its partner.

    A decoy for pair (i, j) is the partner-side position j' of another eligible
    pair whose nucleotide equals that of j, so i is as complementary to j' as it
    is to j. `scope="stem"` draws decoys from the pair's own stem and
    `scope="family"` from every stem of the family. Decoys at j - 1 or j + 1 are
    dropped: PS already compares the partner with those two positions.
    """
    if scope not in ("stem", "family"):
        raise ValueError(f"scope must be 'stem' or 'family', got {scope!r}")
    pair = eligible_pairs[index]
    decoys = []
    for other_index, other in enumerate(eligible_pairs):
        if other_index == index:
            continue
        if scope == "stem" and other["stem_idx"] != pair["stem_idx"]:
            continue
        if sequence[other["j"]] != sequence[pair["j"]]:
            continue
        if abs(other["j"] - pair["j"]) <= 1 or other["j"] == pair["i"]:
            continue
        decoys.append(other["j"])
    return decoys
