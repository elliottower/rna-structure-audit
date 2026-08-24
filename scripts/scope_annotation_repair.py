"""What repairing the unattributed annotations would cost, and what it would move.

Run:  PYTHONPATH=src:scripts uv run --no-project --with torch --with "numpy<2" \
        --with scipy --with tqdm --python 3.12 python scripts/scope_annotation_repair.py

`match_families_to_seed_members.py` establishes that of the 11 records naming no
seed member, one stores a sequence that is an exact seed member and the rest store
sequences that are neither members nor substrings of members. Two repairs follow:

  ADOPT STRUCTURE  the sequence is a seed member, so its projected consensus
                   structure replaces the fabricated one and the sequence is
                   untouched. Nothing about the family's identity changes.
  ADOPT MEMBER     the sequence is unattributable, so the record is replaced by a
                   seed member of the same Rfam family, sequence and structure
                   together. The family keeps its name and class and becomes a
                   different specific molecule.

The member chosen under ADOPT MEMBER is the one closest in length to the stored
sequence among members whose projected structure is at most 5% non-canonical. The
rule is fixed here rather than chosen per family: length is a reported property of
the panel and the models are length-sensitive, so preserving it keeps the repair
from moving a second thing at the same time.

This script proposes and measures. It writes no data file. What it decides is
whether a repair moves the Rung 3 confirmatory set, which
PREREGISTRATION_PHASE6_V2.md fixes by criteria, and how many stored evaluations
it invalidates.

No model is run.
"""

import json
from pathlib import Path

import parse_stockholm
from audit_annotation_validity import non_canonical_fraction
from match_families_to_seed_members import ASSIGNED_RFAM_ID, seed_members
from rna_structure_audit.rungs.rung3 import _get_eligible_pairs, _parse_stems

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data" / "rfam_families"

MIN_WC_PAIRS = 15
MIN_ELIGIBLE = 5
MAX_NON_CANONICAL = 0.05
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
EXCLUDED_POST_HOC = {"U2_snRNA_stem", "hammerhead_ribozyme", "RNaseP_specificity"}

# Above this the annotation cannot be describing its own sequence: the 41 records
# that do name a seed member reach at most 18%.
DEFECT_THRESHOLD = 0.25


def eligibility(sequence: str, structure: str) -> dict:
    sequence = sequence.upper()
    stems = _parse_stems(structure, sequence)
    eligible = _get_eligible_pairs(sequence, stems)
    n_wc = sum(len(stem) for stem in stems)
    return dict(n_wc=n_wc, n_eligible=len(eligible),
                qualifies=n_wc >= MIN_WC_PAIRS and len(eligible) >= MIN_ELIGIBLE)


def propose(name: str, record: dict) -> dict | None:
    stored = record["sequence"].upper().replace("T", "U")
    rfam_id = record.get("rfam_id") or ASSIGNED_RFAM_ID.get(name)
    if rfam_id is None or not (REPO / f"data/rfam_seeds/{rfam_id}.sto").exists():
        return None
    _, members = seed_members(rfam_id)

    for member in members:
        if member["sequence"] == stored:
            return {"repair": "ADOPT STRUCTURE", "accession": member["accession"],
                    "sequence": stored, "dot_bracket": member["dot_bracket"]}

    clean = [m for m in members
             if non_canonical_fraction(m)[0] <= MAX_NON_CANONICAL]
    if not clean:
        return None
    best = min(clean, key=lambda m: (abs(len(m["sequence"]) - len(stored)),
                                     m["accession"]))
    return {"repair": "ADOPT MEMBER", "accession": best["accession"],
            "sequence": best["sequence"], "dot_bracket": best["dot_bracket"]}


def stored_rung3_state() -> dict:
    per = json.loads((REPO / "results/rinalmo_phase6_ps.json")
                     .read_text())["results"]["per_rna"]
    return {name: not (body.get("skipped") or body.get("best_ps") is None)
            for name, body in per.items() if isinstance(body, dict)}


def main() -> int:
    records = {p.stem: json.loads(p.read_text())
               for p in sorted(FAMILIES.glob("*.json"))}
    stored_state = stored_rung3_state()

    unattributed = [n for n, r in records.items()
                    if "seed alignment member" not in r.get("source", "")]
    defective = [n for n in unattributed
                 if non_canonical_fraction(records[n])[0] > DEFECT_THRESHOLD]

    print(f"{len(unattributed)} records name no seed member; "
          f"{len(defective)} of them exceed {DEFECT_THRESHOLD:.0%} non-canonical.\n")

    header = (f"{'family':<20} {'repair':<16} {'len':>9} {'WC':>9} {'elig':>9} "
              f"{'non-canon':>11} {'qualifies':>11}")
    print(header)
    print("-" * len(header))

    proposals, entering, leaving, unrepairable = {}, [], [], []
    for name in sorted(defective):
        record = records[name]
        before = eligibility(record["sequence"], record["dot_bracket"])
        bad_before = non_canonical_fraction(record)[0]
        plan = propose(name, record)
        if plan is None:
            unrepairable.append(name)
            print(f"{name:<20} {'NO CANDIDATE':<16} {record['length']:>9} "
                  f"{before['n_wc']:>9} {before['n_eligible']:>9} "
                  f"{bad_before:>10.0%} {'yes' if before['qualifies'] else 'no':>11}")
            continue

        after = eligibility(plan["sequence"], plan["dot_bracket"])
        bad_after = non_canonical_fraction(plan)[0]
        proposals[name] = plan
        print(f"{name:<20} {plan['repair']:<16} "
              f"{record['length']:>4} -> {len(plan['sequence']):<2} "
              f"{before['n_wc']:>4} -> {after['n_wc']:<2} "
              f"{before['n_eligible']:>4} -> {after['n_eligible']:<2} "
              f"{bad_before:>5.0%} -> {bad_after:<4.0%} "
              f"{('yes' if before['qualifies'] else 'no'):>4} -> "
              f"{'yes' if after['qualifies'] else 'no':<3}")
        print(f"{'':<20} {'':<16} adopts {plan['accession']}")

        was_scored = stored_state.get(name, False)
        in_play = name not in QUARANTINE and name not in EXCLUDED_POST_HOC
        if in_play and after["qualifies"] and not was_scored:
            entering.append(name)
        if in_play and not after["qualifies"] and was_scored:
            leaving.append(name)

    n_now = sum(1 for n, was in stored_state.items()
                if was and n not in QUARANTINE)
    print(f"\nRung 3 confirmatory set")
    print(f"  N = {n_now} as the deposited records stand.")
    print(f"  Of the {len(defective)} defective families, "
          f"{sum(1 for n in defective if stored_state.get(n))} carry a scored "
          f"Rung 3 record.")
    if entering or leaving:
        print(f"  A repair would move N: entering {entering or 'none'}, "
              f"leaving {leaving or 'none'}.")
    else:
        print("  No repair moves N. Every defective family is outside the "
              "confirmatory set\n  both before and after, so Rung 3 is "
              "untouched by this and the registration\n  is not engaged.")
    if unrepairable:
        print(f"  No offline candidate for: {', '.join(unrepairable)}")

    rung12 = [n for n in defective
              if n not in QUARANTINE and n not in EXCLUDED_POST_HOC]
    files = sorted(set(REPO.glob("results/**/*phases*1*5*.json"))
                   | set(REPO.glob("data/gpu_results/**/*phases*1*5*.json")))
    print(f"\nRungs 1-2")
    print(f"  {len(rung12)} of the defective families are inside the reported "
          f"Rungs 1-2 analysis.")
    print(f"  {len(files)} stored files x {len(rung12)} families = "
          f"{len(files) * len(rung12)} family-model evaluations to repeat.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
