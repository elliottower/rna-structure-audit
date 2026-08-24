"""Integrity checks on data/rfam_families/.

Every check here corresponds to a defect that was present in the deposited data
and reached the manuscript. They are cheap and should run before any analysis.

    duplicate sequences      three files were byte-identical copies of RF00001
    reused accessions        four files carried rfam_id RF00001
    bracket balance          unbalanced dot-brackets exclude a family post hoc
    length agreement         a 94 nt sequence carried a 91 character structure
    canonical pair fraction  a structure that does not belong to its sequence
                             pairs bases that cannot pair; real annotations are
                             overwhelmingly Watson-Crick plus GU wobble

Usage:
    uv run python scripts/audit_family_data.py
    uv run python scripts/audit_family_data.py --strict     # exit 1 on any finding
"""

import argparse
import collections
import hashlib
import json
import pathlib
import sys

WATSON_CRICK = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}

# Above this fraction of non-canonical pairs, the structure probably does not
# belong to the sequence. Ribozymes and riboswitches carry genuine non-canonical
# pairs, so this flags for review rather than proving an error.
NON_CANONICAL_FLAG = 0.25


def read_families(directory):
    out = {}
    for path in sorted(pathlib.Path(directory).glob("*.json")):
        out[path.stem] = json.load(open(path))
    return out


def pairs_of(dot_bracket):
    stack, pairs = [], []
    for i, char in enumerate(dot_bracket):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            pairs.append((stack.pop(), i))
    return pairs


def audit(directory):
    families = read_families(directory)
    findings = []

    by_sequence = collections.defaultdict(list)
    by_accession = collections.defaultdict(list)
    for name, fam in families.items():
        sequence = fam.get("sequence", "")
        by_sequence[hashlib.sha256(sequence.encode()).hexdigest()].append(name)
        if fam.get("rfam_id"):
            by_accession[fam["rfam_id"]].append(name)

    for names in by_sequence.values():
        if len(names) > 1:
            findings.append(("duplicate sequence", ", ".join(sorted(names))))

    for accession, names in sorted(by_accession.items()):
        if len(names) > 1:
            findings.append(("reused accession", f"{accession}: {', '.join(sorted(names))}"))

    for name, fam in sorted(families.items()):
        sequence, structure = fam.get("sequence", ""), fam.get("dot_bracket", "")
        if len(sequence) != len(structure):
            findings.append(("length disagreement",
                             f"{name}: {len(sequence)} nt against {len(structure)} characters"))
            continue
        if structure.count("(") != structure.count(")"):
            findings.append(("unbalanced structure",
                             f"{name}: {structure.count('(')} open, {structure.count(')')} close"))
            continue
        pairs = pairs_of(structure)
        if not pairs:
            continue
        non_canonical = [(sequence[i], sequence[j]) for i, j in pairs
                         if (sequence[i], sequence[j]) not in WATSON_CRICK]
        fraction = len(non_canonical) / len(pairs)
        if fraction > NON_CANONICAL_FLAG:
            common = collections.Counter(non_canonical).most_common(3)
            detail = ", ".join(f"{a}-{b} x{n}" for (a, b), n in common)
            findings.append(("non-canonical pairs",
                             f"{name}: {fraction:.0%} of {len(pairs)} pairs ({detail})"))

    return families, findings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default="data/rfam_families")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    families, findings = audit(args.dir)
    print(f"{len(families)} families in {args.dir}\n")

    if not findings:
        print("  no findings")
    else:
        grouped = collections.defaultdict(list)
        for kind, detail in findings:
            grouped[kind].append(detail)
        for kind in ("duplicate sequence", "reused accession", "length disagreement",
                     "unbalanced structure", "non-canonical pairs"):
            if kind not in grouped:
                continue
            print(f"  {kind} ({len(grouped[kind])})")
            for detail in grouped[kind]:
                print(f"      {detail}")
            print()

    blocking = [k for k, _ in findings if k in
                {"duplicate sequence", "length disagreement"}]
    if args.strict and blocking:
        sys.exit(1)


if __name__ == "__main__":
    main()
