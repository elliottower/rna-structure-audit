"""Do the per-position embeddings for NT v2 and DNABERT-2 sit at the positions they claim?

Run:  uv run --no-project --with "transformers==4.44.2" --python 3.12 python \
          scripts/audit_token_alignment.py

No model weights are loaded, so no deep-learning framework is needed; the pin is
the one the deposited NT v2 and DNABERT-2 runs used.

Rungs 1 and 2 read the cosine distance at the mutated nucleotide. Eight of the
ten models emit one token per nucleotide, so that position is a row of the
hidden state. NT v2 and DNABERT-2 do not, and `phases_1_to_5.py:61-76` expands
their token embeddings back to nucleotides with two closed-form guesses:
`i // 6` for NT v2, which assumes the sequence tiles into non-overlapping 6-mers,
and `i * n_tokens / n_nucleotides` for DNABERT-2, which assumes every byte-pair
token is the same width. Neither assumption is a property of the tokenizer.

This script recovers the true span of every token by decoding it -- both
tokenizers emit literal nucleotide strings, so the spans are exact and the
reconstruction is checked against the sequence -- and reports three quantities:

  misassigned      nucleotides the expansion hands to a token whose span does
                   not contain them
  retokenized      mutations after which the token count changes, so the
                   expansion's own denominator moves
  wrong pair       mutations where the wild-type row and the mutant row that the
                   metric subtracts describe different positions of the molecule

The third is the one that decides whether a number means anything: a distance
between two positions that are not the same position is not a measurement of the
mutation.

The attention-contact rung maps the other way, from the nucleotide contact map
onto token pairs, and `_aggregate_contacts_to_tokens` at `phases_1_to_5.py:79-89`
gives token `t` the window `[6t, 6t+6)` for both models. That is the 6-mer
assumption again, applied where it was never true, so the same comparison is
reported for the contact windows.

"""

import json
from pathlib import Path

from transformers import AutoTokenizer

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data/rfam_families"

MODELS = [
    ("nt", "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species", "6-mer"),
    ("dnabert2", "zhihan1996/DNABERT-2-117M", "bpe"),
]

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def expansion_index(scheme: str, i: int, n_tokens: int, n_nucleotides: int) -> int:
    """The token `phases_1_to_5.py` assigns to nucleotide `i`."""
    if scheme == "6-mer":
        return min(i // 6, n_tokens - 1)
    return min(int(i * n_tokens / n_nucleotides), n_tokens - 1)


def token_spans(tokenizer, sequence: str) -> list[str] | None:
    """The content tokens as literal strings, or None if they do not tile the
    sequence -- an unknown or non-literal token makes every index below a guess."""
    ids = tokenizer(sequence)["input_ids"]
    pieces = tokenizer.convert_ids_to_tokens(ids)
    content = [p for p in pieces if set(p) <= set("ACGTN") and p]
    return content if "".join(content) == sequence else None


def true_index(spans: list[str], i: int) -> int:
    """The token whose span contains nucleotide `i`."""
    start = 0
    for index, piece in enumerate(spans):
        start += len(piece)
        if i < start:
            return index
    return len(spans) - 1


def span_of(spans: list[str], index: int) -> tuple[int, int]:
    """The half-open nucleotide range a token row describes. Two rows are the
    same position of the molecule when their ranges match; comparing row indices
    instead would call two rows equal across tokenizations that differ."""
    start = sum(len(piece) for piece in spans[:index])
    return start, start + len(spans[index])


def contact_window(index: int, n_nucleotides: int) -> tuple[int, int]:
    """The window `_aggregate_contacts_to_tokens` gives a token row."""
    start = index * 6
    return start, min(start + 6, n_nucleotides)


def audit(key: str, model_id: str, scheme: str, families: list[dict]) -> None:
    tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)

    misassigned = positions = 0
    rows = empty_rows = disjoint_rows = 0
    overlap = 0.0
    retokenized = wrong_pair = trials = 0
    unreadable = []
    widths: dict[int, int] = {}

    for family in families:
        sequence = family["sequence"].replace("U", "T")
        spans = token_spans(tokenizer, sequence)
        if spans is None:
            unreadable.append(family["name"])
            continue
        for piece in spans:
            widths[len(piece)] = widths.get(len(piece), 0) + 1

        n_nuc, n_tok = len(sequence), len(spans)
        for index in range(n_tok):
            rows += 1
            true_start, true_end = span_of(spans, index)
            start, end = contact_window(index, n_nuc)
            if start >= n_nuc:
                empty_rows += 1
            shared = max(0, min(true_end, end) - max(true_start, start))
            if shared == 0:
                disjoint_rows += 1
            overlap += shared / (true_end - true_start)

        for i in range(n_nuc):
            positions += 1
            guessed = expansion_index(scheme, i, n_tok, n_nuc)
            if guessed != true_index(spans, i):
                misassigned += 1

            if sequence[i] not in COMPLEMENT:
                continue
            mutant = sequence[:i] + COMPLEMENT[sequence[i]] + sequence[i + 1:]
            mutant_spans = token_spans(tokenizer, mutant)
            if mutant_spans is None:
                continue
            trials += 1
            if len(mutant_spans) != n_tok:
                retokenized += 1
            mutant_guess = expansion_index(scheme, i, len(mutant_spans), n_nuc)
            if span_of(spans, guessed) != span_of(mutant_spans, mutant_guess):
                wrong_pair += 1

    print(f"\n  {key}  ({model_id}, {scheme})")
    if unreadable:
        print(f"    {len(unreadable)} families whose tokens do not tile the "
              f"sequence, not measured: {', '.join(unreadable)}")
    print(f"    token widths: " + ", ".join(
        f"{w} nt x{count}" for w, count in sorted(widths.items())))
    print(f"    misassigned {misassigned}/{positions} nucleotides "
          f"= {misassigned / positions:.3f}")
    print(f"    retokenized {retokenized}/{trials} mutations "
          f"= {retokenized / trials:.3f}")
    print(f"    wrong pair  {wrong_pair}/{trials} mutations "
          f"= {wrong_pair / trials:.3f}")
    print(f"    contact windows: {empty_rows}/{rows} token rows fall past the end "
          f"of the\n                     sequence and score against an empty contact "
          f"block,\n                     {disjoint_rows} share no nucleotide with the "
          f"token they stand for,\n                     mean overlap {overlap / rows:.3f}")
    if wrong_pair == 0 and misassigned:
        print("    every misassignment is the same in both arms, so the metric "
              "compares\n    one position of the molecule against itself at a "
              "coarser resolution")


def main() -> None:
    records = [json.loads(p.read_text()) for p in sorted(FAMILIES.glob("*.json"))]
    families = [r for r in records if "excluded" not in r]
    print(f"{len(families)} analyzed families of {len(records)} curated\n")
    print("  Eight of the ten models tokenize one nucleotide per token and do not\n"
          "  reach this code at all. These two do.")
    for key, model_id, scheme in MODELS:
        audit(key, model_id, scheme, families)


if __name__ == "__main__":
    main()
