"""Apply one stated repair rule to every record that names no Rfam seed member.

Run:  PYTHONPATH=src:scripts uv run --no-project --with torch --with "numpy<2" \
        --with scipy --with tqdm --python 3.12 python scripts/repair_annotations.py
      ... --apply     to write the repaired records

Eleven of the 52 records name no seed member. Six sit inside the reported Rungs
1-2 analysis, three were excluded post hoc for unbalanced brackets, and two are
quarantined. The rule is applied to all eleven and the outcome recorded for each,
including where the outcome is "no repair" -- repairing only the six would leave
the choice of which records got repaired unexplained, and an unexplained
selection is a forking path in a way that re-scoring under a frozen criterion is
not.

The rule, in order:

  1. Annotation at most 25% non-canonical           -> NO REPAIR. It describes its
     own sequence; nothing is wrong with it.
  2. Sequence is an exact seed member               -> ADOPT STRUCTURE. Only the
     structure was fabricated; the sequence is kept byte for byte.
  3. No seed member projects to a clean structure   -> DROP.
  4. The family name asserts a subdomain            -> DROP. Seed members are whole
     molecules, so a substitute would answer to a different name.
  5. The family name asserts an organism            -> DEFER. A substitute would
     leave the name asserting something false; the repair needs an annotation for
     the sequence actually stored, which is not derivable from a seed alignment.
  6. Closest clean member differs in length by      -> DROP.
     more than 30%
  7. Otherwise                                      -> ADOPT MEMBER, the member
     closest in length among those at most 5% non-canonical.

Rule 6's threshold sits in an empty band. The length differences the candidates
actually produce are 0, 1, 2, 18, 71, 210 and 269 percent, so every threshold
between 18 and 71 partitions them identically and the number is not doing work.

Rules 4 and 5 are the only ones reading anything other than the data, so the
claim each name makes is written out below and can be checked against the name.

No model is run. Without --apply nothing is written but the manifest.
"""

import argparse
import json
import math
from pathlib import Path

from audit_annotation_validity import non_canonical_fraction
from match_families_to_seed_members import ASSIGNED_RFAM_ID, seed_members
from rna_structure_audit.rungs.rung3 import _get_eligible_pairs, _parse_stems

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data" / "rfam_families"
MANIFEST = REPO / "docs" / "annotation_repair_manifest.json"

MIN_WC_PAIRS = 15
MIN_ELIGIBLE = 5
MAX_NON_CANONICAL = 0.05
DEFECT_THRESHOLD = 0.25
MAX_LENGTH_DELTA = 0.30
MULTIPLIER, ALPHA = 4, 0.05          # H1 condition (a), PREREGISTRATION_PHASE6_V2.md:96

QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
EXCLUDED_POST_HOC = {"U2_snRNA_stem", "hammerhead_ribozyme", "RNaseP_specificity"}

# What each name asserts beyond its Rfam family. A seed member is a whole molecule
# from one organism, so a substitution can silently falsify either claim.
SUBDOMAIN_NAMES = {
    "IRES_HCV_domainII": "domain II of the HCV IRES, not the whole IRES",
    "SRP_RNA_helix8": "helix 8 of SRP RNA, not the whole SRP RNA",
    "RNaseP_specificity": "the specificity domain of RNase P, not the whole RNA",
    "U2_snRNA_stem": "a stem of U2 snRNA, not the whole snRNA",
}
ORGANISM_NAMES = {
    "5S_rRNA_ecoli": "E. coli; RF00001's seed carries no E. coli member",
    "tRNA_Phe_yeast": "yeast tRNA-Phe",
    "tRNA_Ala_human": "human tRNA-Ala",
}


def eligibility(sequence: str, structure: str) -> dict:
    sequence = sequence.upper()
    stems = _parse_stems(structure, sequence)
    eligible = _get_eligible_pairs(sequence, stems)
    n_wc = sum(len(stem) for stem in stems)
    return dict(n_wc=n_wc, n_eligible=len(eligible),
                qualifies=n_wc >= MIN_WC_PAIRS and len(eligible) >= MIN_ELIGIBLE)


def decide(name: str, record: dict) -> dict:
    stored = record["sequence"].upper().replace("T", "U")
    fraction = non_canonical_fraction(record)[0]
    out = {"name": name, "non_canonical_before": round(fraction, 4),
           "length_before": record["length"],
           "quarantined": name in QUARANTINE,
           "excluded_post_hoc": name in EXCLUDED_POST_HOC}

    if fraction <= DEFECT_THRESHOLD:
        return {**out, "disposition": "NO REPAIR",
                "reason": f"{fraction:.0%} non-canonical; the annotation describes "
                          f"its own sequence"}

    rfam_id = record.get("rfam_id") or ASSIGNED_RFAM_ID.get(name)
    out["rfam_id"] = rfam_id
    if rfam_id is None or not (REPO / f"data/rfam_seeds/{rfam_id}.sto").exists():
        return {**out, "disposition": "DROP",
                "reason": "no cached seed alignment for this family"}

    _, members = seed_members(rfam_id)
    for member in members:
        if member["sequence"] == stored:
            after = eligibility(stored, member["dot_bracket"])
            return {**out, "disposition": "ADOPT STRUCTURE",
                    "accession": member["accession"],
                    "reason": "the stored sequence is this seed member exactly; "
                              "only the structure was fabricated",
                    "sequence": stored, "dot_bracket": member["dot_bracket"],
                    "length_after": len(stored),
                    "non_canonical_after": round(
                        non_canonical_fraction(member)[0], 4), **after}

    clean = [m for m in members
             if non_canonical_fraction(m)[0] <= MAX_NON_CANONICAL]
    if not clean:
        return {**out, "disposition": "DROP",
                "reason": f"no {rfam_id} member projects to a structure under "
                          f"{MAX_NON_CANONICAL:.0%} non-canonical"}
    if name in SUBDOMAIN_NAMES:
        return {**out, "disposition": "DROP",
                "reason": f"the name asserts {SUBDOMAIN_NAMES[name]}; every seed "
                          f"member is a whole molecule"}
    if name in ORGANISM_NAMES:
        return {**out, "disposition": "DEFER",
                "reason": f"the name asserts {ORGANISM_NAMES[name]}; a substitute "
                          f"would leave the name false, and the stored sequence "
                          f"needs an annotation from outside the seed alignment"}

    best = min(clean, key=lambda m: (abs(len(m["sequence"]) - len(stored)),
                                     m["accession"]))
    delta = abs(len(best["sequence"]) - len(stored)) / len(stored)
    if delta > MAX_LENGTH_DELTA:
        return {**out, "disposition": "DROP",
                "reason": f"closest clean member differs in length by {delta:.0%}"}

    after = eligibility(best["sequence"], best["dot_bracket"])
    return {**out, "disposition": "ADOPT MEMBER", "accession": best["accession"],
            "reason": f"closest clean {rfam_id} member by length ({delta:.0%} apart)",
            "sequence": best["sequence"], "dot_bracket": best["dot_bracket"],
            "length_after": len(best["sequence"]),
            "non_canonical_after": round(non_canonical_fraction(best)[0], 4),
            **after}


def stored_rung3_state() -> dict:
    per = json.loads((REPO / "results/rinalmo_phase6_ps.json")
                     .read_text())["results"]["per_rna"]
    return {name: not (body.get("skipped") or body.get("best_ps") is None)
            for name, body in per.items() if isinstance(body, dict)}


def gate(n: int) -> int:
    return math.ceil(MULTIPLIER * ALPHA * n)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true",
                        help="write the repaired records into data/rfam_families/")
    args = parser.parse_args()

    records = {p.stem: json.loads(p.read_text())
               for p in sorted(FAMILIES.glob("*.json"))}
    unattributed = sorted(n for n, r in records.items()
                          if "seed alignment member" not in r.get("source", ""))
    decisions = [decide(n, records[n]) for n in unattributed]

    header = (f"{'family':<20}{'disposition':<17}{'len':>11}{'non-canon':>14}"
              f"{'WC':>6}{'Rung 3':>12}")
    print(f"The rule applied to all {len(decisions)} records naming no seed member\n")
    print(header)
    print("-" * len(header))
    for d in decisions:
        length = (f"{d['length_before']} -> {d['length_after']}"
                  if "length_after" in d else str(d["length_before"]))
        canon = (f"{d['non_canonical_before']:.0%} -> {d['non_canonical_after']:.0%}"
                 if "non_canonical_after" in d else f"{d['non_canonical_before']:.0%}")
        wc = f"{d['n_wc']}" if "n_wc" in d else "-"
        rung3 = ("qualifies" if d.get("qualifies") else
                 "-" if "qualifies" not in d else "no")
        flag = " *" if d["quarantined"] or d["excluded_post_hoc"] else ""
        print(f"{d['name']:<20}{d['disposition']:<17}{length:>11}{canon:>14}"
              f"{wc:>6}{rung3:>12}{flag}")
        print(f"{'':<20}{d['reason']}")
    print("\n  * already outside the confirmatory analysis "
          "(quarantined or excluded post hoc)")

    stored = stored_rung3_state()
    n_now = sum(1 for n, was in stored.items() if was and n not in QUARANTINE)
    entering = [d["name"] for d in decisions
                if d.get("qualifies") and not stored.get(d["name"])
                and d["name"] not in QUARANTINE]
    # mir_21_precursor entered under the four annotations corrected at 5fc8914.
    n_after = n_now + 1 + len(entering)

    print(f"\nRung 3 confirmatory set")
    print(f"  N = {n_now} deposited, + 1 (mir_21_precursor, corrected at 5fc8914)"
          f" + {len(entering)} entering here = {n_after}")
    print(f"  entering: {', '.join(entering) or 'none'}")
    print(f"  H1 condition (a) gate: ceil({MULTIPLIER} * {ALPHA} * {n_now}) = "
          f"{gate(n_now)} registered, ceil({MULTIPLIER} * {ALPHA} * {n_after}) = "
          f"{gate(n_after)} after repair")
    if gate(n_after) != gate(n_now):
        print("  The registered decision threshold MOVES. This must be logged "
              "before the re-run.")

    # Rungs 1-2 report over the panel minus the three excluded post hoc for
    # unbalanced brackets. Repairing an annotation can return a family to that
    # panel as well as remove one from it.
    n_families = len(records)
    rung12_now = n_families - len(EXCLUDED_POST_HOC)
    removed = [d["name"] for d in decisions
               if d["disposition"] in ("DROP", "DEFER")
               and d["name"] not in EXCLUDED_POST_HOC
               and d["name"] not in QUARANTINE]
    returned = [d["name"] for d in decisions
                if d["disposition"] in ("ADOPT MEMBER", "ADOPT STRUCTURE")
                and d["name"] in EXCLUDED_POST_HOC]
    rung12_after = rung12_now - len(removed) + len(returned)
    print(f"\nRungs 1-2 panel")
    print(f"  {rung12_now} families reported ({n_families} minus the "
          f"{len(EXCLUDED_POST_HOC)} excluded post hoc for unbalanced brackets)")
    print(f"  removed: {', '.join(removed) or 'none'}")
    print(f"  returned by repair: {', '.join(returned) or 'none'}")
    print(f"  {rung12_now} - {len(removed)} + {len(returned)} = {rung12_after}")

    counts = {}
    for d in decisions:
        counts[d["disposition"]] = counts.get(d["disposition"], 0) + 1
    MANIFEST.write_text(json.dumps(
        {"rule": {"defect_threshold": DEFECT_THRESHOLD,
                  "max_non_canonical": MAX_NON_CANONICAL,
                  "max_length_delta": MAX_LENGTH_DELTA},
         "n_deposited": n_now, "n_after": n_after,
         "rung12_now": rung12_now, "rung12_after": rung12_after,
         "rung12_removed": removed, "rung12_returned": returned,
         "gate_deposited": gate(n_now), "gate_after": gate(n_after),
         "counts": counts, "decisions": decisions}, indent=2) + "\n")
    print(f"\n  {counts}")
    print(f"  manifest -> {MANIFEST.relative_to(REPO)}")

    if not args.apply:
        print("\n  Dry run. No family record was written. Pass --apply to write.")
        return 0

    written = 0
    for d in decisions:
        if d["disposition"] not in ("ADOPT STRUCTURE", "ADOPT MEMBER"):
            continue
        path = FAMILIES / f"{d['name']}.json"
        record = json.loads(path.read_text())
        record["sequence"] = d["sequence"]
        record["dot_bracket"] = d["dot_bracket"]
        record["length"] = len(d["sequence"])
        record["rfam_id"] = d["rfam_id"]
        record["source"] = (f"Rfam {d['rfam_id']} seed alignment member "
                            f"({d['accession']})")
        path.write_text(json.dumps(record, indent=2) + "\n")
        written += 1
    print(f"\n  wrote {written} repaired records into "
          f"{FAMILIES.relative_to(REPO)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
