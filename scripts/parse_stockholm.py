"""Parse Rfam Stockholm seed alignments into per-sequence (sequence, dot_bracket) pairs.

Maps consensus secondary structure (SS_cons) onto each individual sequence
by projecting through alignment columns, then filters for:
  - Valid RNA (ACGU only)
  - Balanced bracket structure
  - Minimum length

Usage:
    uv run python scripts/parse_stockholm.py data/rfam_seeds/RF00005.sto
    uv run python scripts/parse_stockholm.py data/rfam_seeds/ --all --stats
"""

import argparse
import json
import sys
from pathlib import Path

OPEN_BRACKETS = set("(<[{")
CLOSE_BRACKETS = set(")>]}")
BRACKET_MAP = {"(": ")", "<": ">", "[": "]", "{": "}"}
VALID_RNA = set("ACGU")
GAP_CHARS = set(".-")


def parse_stockholm(sto_path):
    """Parse a Stockholm file, returning (family_id, family_name, sequences, ss_cons).

    sequences: dict of {accession: aligned_sequence_string}
    ss_cons: consensus secondary structure string (same length as aligned seqs)
    """
    text = Path(sto_path).read_text()
    lines = text.strip().split("\n")

    family_id = ""
    family_name = ""
    seqs = {}
    ss_cons_parts = []
    seq_order = []

    for line in lines:
        if line.startswith("#=GF AC "):
            family_id = line.split()[2].strip()
        elif line.startswith("#=GF ID "):
            family_name = line.split(None, 2)[2].strip()
        elif line.startswith("#=GC SS_cons"):
            ss_cons_parts.append(line.split(None, 2)[2])
        elif line.startswith("#") or line.startswith("//") or not line.strip():
            continue
        else:
            parts = line.split()
            if len(parts) == 2:
                name, aligned_seq = parts
                if name not in seqs:
                    seqs[name] = []
                    seq_order.append(name)
                seqs[name].append(aligned_seq)

    ss_cons = "".join(ss_cons_parts)
    for name in seqs:
        seqs[name] = "".join(seqs[name])

    return family_id, family_name, seqs, ss_cons, seq_order


def to_standard_bracket(char):
    if char in OPEN_BRACKETS:
        return "("
    if char in CLOSE_BRACKETS:
        return ")"
    return "."


def map_structure_to_sequence(aligned_seq, ss_cons):
    """Map consensus structure onto an individual sequence by removing gap columns.

    Returns (ungapped_sequence, mapped_dot_bracket) or (None, None) if unbalanced.
    """
    if len(aligned_seq) != len(ss_cons):
        return None, None

    ungapped = []
    mapped_ss = []

    for seq_char, ss_char in zip(aligned_seq, ss_cons):
        if seq_char in GAP_CHARS:
            continue
        ungapped.append(seq_char.upper().replace("T", "U"))
        mapped_ss.append(to_standard_bracket(ss_char))

    seq_str = "".join(ungapped)
    ss_str = "".join(mapped_ss)

    opens = ss_str.count("(")
    closes = ss_str.count(")")
    if opens != closes:
        return seq_str, None

    if not set(seq_str).issubset(VALID_RNA):
        return seq_str, None

    return seq_str, ss_str


def extract_sequences(sto_path, min_length=30):
    """Extract all valid (sequence, dot_bracket) pairs from a Stockholm file.

    Returns list of dicts with keys: accession, sequence, dot_bracket, length.
    Filters for: valid RNA, balanced structure, minimum length.
    """
    family_id, family_name, seqs, ss_cons, seq_order = parse_stockholm(sto_path)

    results = []
    stats = {"total": len(seqs), "valid": 0, "unbalanced": 0, "non_rna": 0, "too_short": 0}

    for acc in seq_order:
        aligned_seq = seqs[acc]
        seq_str, ss_str = map_structure_to_sequence(aligned_seq, ss_cons)

        if seq_str is None:
            stats["non_rna"] += 1
            continue

        if ss_str is None:
            if not set(seq_str).issubset(VALID_RNA):
                stats["non_rna"] += 1
            else:
                stats["unbalanced"] += 1
            continue

        if len(seq_str) < min_length:
            stats["too_short"] += 1
            continue

        stats["valid"] += 1
        results.append({
            "accession": acc,
            "sequence": seq_str,
            "dot_bracket": ss_str,
            "length": len(seq_str),
        })

    return family_id, family_name, results, stats


def main():
    parser = argparse.ArgumentParser(description="Parse Stockholm seed alignments")
    parser.add_argument("path", help="Stockholm file or directory of .sto files")
    parser.add_argument("--all", action="store_true", help="Process all .sto files in directory")
    parser.add_argument("--stats", action="store_true", help="Print stats only, no sequences")
    parser.add_argument("--min-length", type=int, default=30, help="Minimum sequence length")
    parser.add_argument("--output", help="Output JSON file")
    args = parser.parse_args()

    path = Path(args.path)
    if path.is_dir() or args.all:
        sto_dir = path if path.is_dir() else path.parent
        files = sorted(sto_dir.glob("*.sto"))
    else:
        files = [path]

    all_families = {}
    grand_total = 0
    grand_valid = 0

    for sto_file in files:
        fam_id, fam_name, seqs, stats = extract_sequences(sto_file, args.min_length)
        all_families[fam_id] = {
            "family_id": fam_id,
            "family_name": fam_name,
            "sequences": seqs,
            "stats": stats,
        }
        grand_total += stats["total"]
        grand_valid += stats["valid"]

        if args.stats:
            print(
                f"{fam_id} ({fam_name}): {stats['valid']}/{stats['total']} valid "
                f"(unbal={stats['unbalanced']}, non_rna={stats['non_rna']}, "
                f"short={stats['too_short']})"
            )
        else:
            print(f"{fam_id} ({fam_name}): {stats['valid']} valid sequences")

    print(f"\nTotal: {grand_valid}/{grand_total} valid sequences across {len(all_families)} families")

    if args.output:
        with open(args.output, "w") as f:
            json.dump(all_families, f, indent=2)
        print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
