"""Does each family's dot-bracket structure belong to its own sequence?

Run:  uv run --no-project --python 3.12 python scripts/audit_annotation_validity.py

A secondary structure annotation pairs positions. If the annotation belongs to
the sequence it is attached to, almost every pair it names is Watson-Crick or a
GU wobble, because those are the pairs that form. A consensus structure copied
onto a different sequence pairs whatever happens to sit at those positions, and
the non-canonical fraction rises accordingly. The fraction is therefore a test
of whether the annotation and the sequence came from the same record.

Commit 5fc8914 corrected four families whose annotations did not belong to their
sequences. This script asks the same question of all 52, and reports the answer
against whether the record names the Rfam seed member it came from, since a
record with no accession behind it has nothing to check the pairing against.

No model is run.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data/rfam_families"
CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# paper_v12.tex:225-227 excludes three families post hoc for unbalanced brackets.
EXCLUDED_POST_HOC = {"U2_snRNA_stem", "hammerhead_ribozyme", "RNaseP_specificity"}


def pairs_of(structure: str) -> list[tuple[int, int]]:
    stack, pairs = [], []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            pairs.append((stack.pop(), i))
    return pairs


def non_canonical_fraction(record: dict) -> tuple[float, int]:
    sequence = record["sequence"].upper().replace("T", "U")
    pairs = pairs_of(record["dot_bracket"])
    usable = [(i, j) for i, j in pairs
              if i < len(sequence) and j < len(sequence)]
    if not usable:
        return 0.0, 0
    bad = sum((sequence[i], sequence[j]) not in CANONICAL for i, j in usable)
    return bad / len(usable), len(usable)


def rung3_membership() -> dict[str, str]:
    """Which families carry a scored Rung 3 record, from the stored results."""
    body = json.loads((REPO / "results/rinalmo_phase6_ps.json")
                      .read_text())["results"]["per_rna"]
    out = {}
    for name, record in body.items():
        if not isinstance(record, dict):
            continue
        if record.get("skipped") or not isinstance(record.get("best_ps"), (int, float)):
            out[name] = "eligible but not scored"
        else:
            out[name] = "IN the Rung 3 confirmatory set"
    return out


def main() -> None:
    rows = []
    for path in sorted(FAMILIES.glob("*.json")):
        record = json.loads(path.read_text())
        fraction, n_pairs = non_canonical_fraction(record)
        traced = "seed alignment member" in (record.get("source") or "")
        rows.append((fraction, record["name"], n_pairs, traced))
    rows.sort(reverse=True)

    print(f"{len(rows)} families, sorted by fraction of annotated pairs that are "
          f"neither Watson-Crick nor GU\n")
    print(f"  {'family':<28}{'non-WC':>8}{'pairs':>7}  {'accession':<11}status")
    for fraction, name, n_pairs, traced in rows:
        status = []
        if name in QUARANTINE:
            status.append("quarantined")
        if name in EXCLUDED_POST_HOC:
            status.append("excluded post hoc")
        print(f"  {name:<28}{fraction:>7.0%}{n_pairs:>7}  "
              f"{'named' if traced else 'ABSENT':<11}{', '.join(status)}")

    traced_rows = [r for r in rows if r[3]]
    untraced_rows = [r for r in rows if not r[3]]
    print(f"\n  worst non-WC fraction among the {len(traced_rows)} records naming "
          f"a seed member: {max(r[0] for r in traced_rows):.0%}")
    print(f"  worst among the {len(untraced_rows)} records with no accession:      "
          f"        {max(r[0] for r in untraced_rows):.0%}")

    scored = rung3_membership()
    live = [(f, n) for f, n, _, _ in rows
            if f > 0.25 and n not in EXCLUDED_POST_HOC and n not in QUARANTINE]
    print(f"\n  families above 25% that are neither excluded nor quarantined, and "
          f"so are\n  inside the reported analysis: {len(live)}")
    for fraction, name in live:
        where = scored.get(name, "not scored in Rung 3")
        print(f"      {name:<28}{fraction:>7.0%}   {where}")


if __name__ == "__main__":
    main()
