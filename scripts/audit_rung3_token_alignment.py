"""How much of the panel did Rung 3 actually score for NT v2 and DNABERT-2?

Run:  uv run --no-project --with "transformers==4.44.2" --with torch \
          --with "numpy<2" --with scipy --with tqdm --python 3.12 \
          python scripts/audit_rung3_token_alignment.py

No model weights are loaded; only the two tokenizers, at the pin the deposited
runs used.

`compute_delta_profiles` in `phase6_compensatory_mutation.py` reads the
perturbation at nucleotide `k` as `emb_wt[k + offset]`, and every entry of
`ADAPTER_OFFSETS` is zero. For the eight character-level models that row is the
nucleotide. For NT v2 and DNABERT-2 the rows are tokens, so a sequence of 106
nucleotides has roughly 18 of them and the guard `k_off < emb_wt.shape[0]`
is false for every position past the eighteenth.

Two things follow, and this script measures both against the deposited results.

  dropped     A pair contributes only when its partner `j` and both stem
              neighbours clear the guard. Stems pair a 5' position against a 3'
              one, so `j` sits in the far half of the molecule and clears
              nothing. A family with no surviving pair yields no PS value at any
              layer and is recorded as `no valid PS values` -- the same field a
              genuine filter writes, which is why the loss is invisible in the
              stored file.

  misplaced   A pair that does survive is scored at token row `j`, which spans
              nucleotides far from `j` itself. The distance between the intended
              nucleotide and the span actually read is reported below.

The registered Rung 3 quantities -- mean PS over non-quarantined families, the
gate denominator, and H3 precision -- are means over whatever survived, so the
count of survivors is the size of the sample those means describe.
"""

import json
from pathlib import Path

from transformers import AutoTokenizer

from phase6_compensatory_mutation import (
    QUARANTINED,
    complement_swap,
    get_eligible_pairs,
    parse_stems,
)
from token_spans import content_spans

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data/rfam_families"
RESULTS = REPO / "results/repaired_panel"

MODELS = [
    ("nt", "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species"),
    ("dnabert2", "zhihan1996/DNABERT-2-117M"),
]
REFERENCE = "rinalmo"


def rows_for(tokenizer, sequence):
    """Number of embedding rows the adapter hands back for this sequence."""
    return len(content_spans(tokenizer, sequence.replace("U", "T")))


def qualifying(family):
    """The eligible pairs, or None where a registered filter already excluded it.

    The filters are `run_phase6`'s own, in its order, so a family this returns
    None for is one no model scored.
    """
    sequence, dot_bracket = family["sequence"], family["dot_bracket"]
    if len(sequence) != len(dot_bracket):
        return None
    stems = parse_stems(dot_bracket, sequence)
    if sum(len(stem) for stem in stems) < 15:
        return None
    eligible = get_eligible_pairs(sequence, stems)
    return eligible if len(eligible) >= 5 else None


def audit(key, model_id, families):
    tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
    stored = json.loads((RESULTS / key / f"{key}_phase6_ps.json").read_text())
    per_rna = stored["results"]["per_rna"]

    predicted_scored, disagreements = [], []
    total_pairs = surviving_pairs = 0
    displacements = []

    for family in families:
        eligible = qualifying(family)
        if eligible is None:
            continue
        sequence = family["sequence"]
        spans = content_spans(tokenizer, sequence.replace("U", "T"))
        n_wt = len(spans)

        mutant_rows = {}
        survivors = 0
        for pair in eligible:
            total_pairs += 1
            pos_i = pair["i"]
            if pos_i not in mutant_rows:
                mutant_rows[pos_i] = rows_for(tokenizer,
                                              complement_swap(sequence, pos_i))
            bound = min(n_wt, mutant_rows[pos_i])
            if max(pair["j"], pair["j_prev"], pair["j_next"]) < bound:
                survivors += 1
                # Row `j` of the hidden state, which is the row the
                # pipeline reads, spans nucleotides `start` to `end`. The
                # distance to the nearer edge is the smallest claim that can be
                # made about how far off the read is.
                start, end = spans[pair["j"]]
                displacements.append(min(abs(pair["j"] - start),
                                         abs(pair["j"] - (end - 1))))
        surviving_pairs += survivors

        was_scored = not per_rna.get(family["name"], {}).get("skipped", False)
        if survivors > 0:
            predicted_scored.append(family["name"])
        if (survivors > 0) != was_scored:
            disagreements.append((family["name"], survivors, was_scored))

    print(f"\n  {key}  ({model_id})")
    print(f"    {len(predicted_scored)} families keep at least one pair, of "
          f"{len(families)} qualifying")
    print(f"    {surviving_pairs}/{total_pairs} eligible pairs clear the bounds "
          f"check = {surviving_pairs / total_pairs:.3f}")
    if displacements:
        print(f"    a surviving pair is read from a token whose nearest edge "
              f"is a median of\n      {sorted(displacements)[len(displacements) // 2]} "
              f"nucleotides from the position it stands for "
              f"(max {max(displacements)})")
    if disagreements:
        print(f"    {len(disagreements)} families where this prediction and the "
              f"deposited file disagree:")
        for name, survivors, was_scored in disagreements:
            print(f"      {name}: {survivors} surviving pairs, "
                  f"stored as {'scored' if was_scored else 'skipped'}")
    else:
        print("    the prediction reproduces the deposited scored/skipped "
              "partition exactly")

    quarantined = [n for n in predicted_scored if n in QUARANTINED]
    print(f"    contributing to the registered means: "
          f"{len(predicted_scored) - len(quarantined)} non-quarantined families")


def main():
    records = [json.loads(p.read_text()) for p in sorted(FAMILIES.glob("*.json"))]
    families = [r for r in records if "excluded" not in r]
    qualified = [f for f in families if qualifying(f) is not None]

    reference = json.loads(
        (RESULTS / REFERENCE / f"{REFERENCE}_phase6_ps.json").read_text())
    scored_ref = sum(1 for e in reference["results"]["per_rna"].values()
                     if not e.get("skipped"))

    print(f"{len(families)} analyzed families, {len(qualified)} pass the "
          f"registered Rung 3 filters")
    print(f"{REFERENCE}, one token per nucleotide, scores {scored_ref} of them")
    for key, model_id in MODELS:
        audit(key, model_id, qualified)


if __name__ == "__main__":
    main()
