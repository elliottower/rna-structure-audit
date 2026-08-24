"""Does `hidden[0, 1:-1, :]` strip the special tokens each tokenizer emits?

Seven of the ten adapters slice one row off each end of the hidden state to drop
a leading and a trailing special token. That is a property of the tokenizer, not
of transformers, and it is wrong for at least one model here: NT v2 emits `<cls>`
and nothing after the sequence, so the slice discards the last content token.

A missing closing token truncates. A missing opening token would be worse -- the
slice would keep a special row and drop a content one, shifting every position by
one -- so both ends are reported.

Two pinned stacks are needed to cover the panel and neither can load the other's
tokenizers, so the script reports what it could load and names what it could not:

  uv run --no-project --with "transformers==4.44.2" --python 3.12 \
      python scripts/audit_special_token_slices.py
  uv run --no-project --with "multimolecule==0.1.0" --with "transformers==4.49.0" \
      --with "numpy<2" --with torch --python 3.12 \
      python scripts/audit_special_token_slices.py

RNA-FM builds its ids by hand as `[2] + nucleotides + [3]` and HyenaDNA, Evo and
Caduceus do not slice at all, so neither reaches this question.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data/rfam_families"
NUCLEOTIDES = set("ACGTUN")

HF = [
    ("nt", "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species", "dna"),
    ("dnabert2", "zhihan1996/DNABERT-2-117M", "dna"),
]
MULTIMOLECULE = [
    ("rinalmo", "multimolecule/rinalmo-giga", "rna"),
    ("utrlm", "multimolecule/utrlm-te_el", "rna"),
    ("ernierna", "multimolecule/ernierna", "rna"),
    ("splicebert", "multimolecule/splicebert", "rna"),
]


def load_tokenizers():
    """Whatever the installed stack can give, and the names of what it cannot."""
    loaded, missing = [], []
    try:
        from transformers import AutoTokenizer
    except ImportError:
        return loaded, [key for key, _, _ in HF + MULTIMOLECULE]
    for key, model_id, alphabet in HF:
        try:
            loaded.append((key, alphabet, AutoTokenizer.from_pretrained(
                model_id, trust_remote_code=True)))
        except Exception:
            missing.append(key)
    try:
        from multimolecule import RnaTokenizer
    except ImportError:
        missing.extend(key for key, _, _ in MULTIMOLECULE)
        return loaded, missing
    for key, model_id, alphabet in MULTIMOLECULE:
        try:
            loaded.append((key, alphabet, RnaTokenizer.from_pretrained(model_id)))
        except Exception:
            missing.append(key)
    return loaded, missing


def report(key, alphabet, tokenizer, families):
    leading, trailing, uncovered, rows = {}, {}, 0, 0
    for family in families:
        sequence = family["sequence"]
        sequence = (sequence.replace("T", "U") if alphabet == "rna"
                    else sequence.replace("U", "T"))
        pieces = tokenizer.convert_ids_to_tokens(tokenizer(sequence)["input_ids"])
        content = [bool(p) and set(p) <= NUCLEOTIDES for p in pieces]
        if not any(content):
            leading["unreadable"] = leading.get("unreadable", 0) + 1
            continue
        head = content.index(True)
        tail = content[::-1].index(True)
        leading[head] = leading.get(head, 0) + 1
        trailing[tail] = trailing.get(tail, 0) + 1
        rows += len(pieces) - 2
        kept = pieces[1:-1]
        uncovered += len(sequence) - sum(
            len(p) for p in kept if p and set(p) <= NUCLEOTIDES)

    def counts(d):
        return ", ".join(f"{k} x{v}" for k, v in sorted(d.items(), key=str))

    verdict = ("correct" if set(leading) == {1} and set(trailing) == {1}
               else "WRONG")
    print(f"  {key:<12} leading specials: {counts(leading):<10} "
          f"trailing: {counts(trailing):<10} slice 1:-1 is {verdict}")
    if uncovered:
        print(f"               the slice leaves {uncovered} nucleotides of the "
              f"panel with no embedding row")


def main():
    records = [json.loads(p.read_text()) for p in sorted(FAMILIES.glob("*.json"))]
    families = [r for r in records if "excluded" not in r]
    loaded, missing = load_tokenizers()
    print(f"{len(families)} analyzed families\n")
    for key, alphabet, tokenizer in loaded:
        report(key, alphabet, tokenizer, families)
    if missing:
        print(f"\n  not loadable in this stack, not measured: "
              f"{', '.join(sorted(missing))}")


if __name__ == "__main__":
    main()
