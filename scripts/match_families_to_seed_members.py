"""For each family with no accession behind it, is its sequence in the Rfam seed?

Run:  PYTHONPATH=scripts uv run --no-project --python 3.12 python \
        scripts/match_families_to_seed_members.py

`audit_annotation_validity.py` separates the 52 families on whether the record
names the Rfam seed member it came from: the 41 that do reach at most 18%
non-canonical pairs, while 9 of the 11 that do not run 24-77%. That says the
annotation is wrong. It does not say which half of the record is wrong.

Three repairs follow from three different answers. If the stored sequence is
itself a seed member, only the structure was fabricated, and the repair is to
adopt that member's projected consensus structure. If the stored sequence is a
*substring* of a member -- which four of these families' names predict, since
they are named for subdomains rather than for whole molecules -- the repair is to
slice the member's projected structure to the fragment's coordinates, keeping the
sequence exactly as stored. If the sequence appears nowhere in the seed, both
halves are unattributable and no offline repair keeps the sequence.

Containment, not Jaccard, is the ranking used: a Jaccard score between an 86-nt
fragment and the 353-nt member it was cut from is near zero, and says nothing.

This script asks which, by looking every stored sequence up in the cached seed
alignment for its Rfam family. It changes nothing; it reports what a repair would
have to do. The alignments under data/rfam_seeds/ are already present, so no
download happens and no model is run.
"""

import json
from pathlib import Path

import parse_stockholm
from audit_annotation_validity import non_canonical_fraction

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data/rfam_families"
SEEDS = REPO / "data/rfam_seeds"

# IRES_HCV_domainII stores `rfam_id: null`. The candidate is assigned here rather
# than read from the record, so it is stated rather than assumed; the script
# prints the alignment's own #=GF ID line so the assignment can be checked.
ASSIGNED_RFAM_ID = {"IRES_HCV_domainII": "RF00061"}

KMER = 8


def kmers(sequence: str) -> set[str]:
    return {sequence[i:i + KMER] for i in range(len(sequence) - KMER + 1)}


def containment(query: str, target: str) -> float:
    """Fraction of the query's k-mers present in the target. 1.0 for a substring."""
    kq, kt = kmers(query), kmers(target)
    if not kq or not kt:
        return 0.0
    return len(kq & kt) / len(kq)


def sliced_structure(structure: str, start: int, length: int) -> str | None:
    """The member's structure over a fragment, if every pair it names is internal.

    A slice through a helix leaves an unmatched bracket. Those positions become
    unpaired, which is the correct reading: the partner is outside the fragment.
    """
    window = structure[start:start + length]
    if len(window) != length:
        return None
    depth, out = 0, []
    for char in window:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        out.append(char)
    # Drop brackets whose partner falls outside the window.
    stack, keep = [], [True] * length
    for i, char in enumerate(window):
        if char == "(":
            stack.append(i)
        elif char == ")":
            if stack:
                stack.pop()
            else:
                keep[i] = False
    for i in stack:
        keep[i] = False
    return "".join(c if keep[i] else "." for i, c in enumerate(window))


def hamming_identity(a: str, b: str) -> float | None:
    if len(a) != len(b) or not a:
        return None
    return sum(x == y for x, y in zip(a, b)) / len(a)


def seed_members(rfam_id: str) -> tuple[str, list[dict]]:
    """Every member of a cached seed alignment, with the consensus projected on."""
    path = SEEDS / f"{rfam_id}.sto"
    _, family_name, seqs, ss_cons, order = parse_stockholm.parse_stockholm(path)
    members = []
    for accession in order:
        sequence, structure = parse_stockholm.map_structure_to_sequence(
            seqs[accession], ss_cons)
        # A member whose projected consensus comes out unbalanced has no
        # structure to donate, and cannot be the source of one either.
        if sequence is None or structure is None:
            continue
        members.append({"accession": accession, "sequence": sequence,
                        "dot_bracket": structure})
    return family_name, members


def main() -> None:
    records = {p.stem: json.loads(p.read_text())
               for p in sorted(FAMILIES.glob("*.json"))}
    unattributed = [name for name, r in records.items()
                    if "seed alignment member" not in r.get("source", "")]

    print(f"{len(unattributed)} of {len(records)} records name no seed member.\n")

    for name in sorted(unattributed):
        record = records[name]
        stored = record["sequence"].upper().replace("T", "U")
        bad, npairs = non_canonical_fraction(record)
        rfam_id = record.get("rfam_id") or ASSIGNED_RFAM_ID.get(name)
        tag = "" if record.get("rfam_id") else "  (assigned here)"

        print(f"{name}")
        print(f"  stored     {len(stored)} nt, {npairs} pairs, "
              f"{bad:.0%} non-canonical, source {record['source']!r}")

        if rfam_id is None:
            print("  seed       no Rfam family on the record and none assigned\n")
            continue
        if not (SEEDS / f"{rfam_id}.sto").exists():
            print(f"  seed       {rfam_id} not cached under data/rfam_seeds/\n")
            continue

        family_name, members = seed_members(rfam_id)
        print(f"  seed       {rfam_id} {family_name!r}{tag}, "
              f"{len(members)} members with a projectable structure")

        exact = [m for m in members if m["sequence"] == stored]
        if exact:
            m = exact[0]
            mb, mp = non_canonical_fraction(m)
            print(f"  VERDICT    sequence IS seed member {m['accession']}; only the "
                  f"structure is wrong")
            print(f"             its projected structure: {mp} pairs, "
                  f"{mb:.0%} non-canonical\n")
            continue

        contained = [m for m in members if stored in m["sequence"]]
        if contained:
            m = contained[0]
            start = m["sequence"].index(stored)
            structure = sliced_structure(m["dot_bracket"], start, len(stored))
            repaired = {"sequence": stored, "dot_bracket": structure}
            mb, mp = non_canonical_fraction(repaired)
            print(f"  VERDICT    sequence is a SUBSTRING of member {m['accession']} "
                  f"at offset {start}")
            print(f"             sliced structure: {mp} pairs, {mb:.0%} non-canonical; "
                  f"sequence can be kept as stored\n")
            continue

        ranked = sorted(members, key=lambda m: -containment(stored, m["sequence"]))
        best = ranked[0]
        score = containment(stored, best["sequence"])
        identity = hamming_identity(stored, best["sequence"])
        mb, mp = non_canonical_fraction(best)
        shown = f"{identity:.0%} identity" if identity is not None \
            else f"{len(best['sequence'])} nt against {len(stored)}"
        print(f"  closest    {best['accession']}  {KMER}-mer containment "
              f"{score:.2f}, {shown}")
        print(f"             its projected structure: {mp} pairs, {mb:.0%} non-canonical")
        print(f"  VERDICT    sequence is neither a member nor a substring of one; "
              f"no offline repair keeps it\n")


if __name__ == "__main__":
    main()
