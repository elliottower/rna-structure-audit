"""How much has to be re-run after the four family annotations were corrected.

Run:  uv run --no-project --with torch --with numpy --with tqdm --python 3.12 \
          python scripts/scope_rerun.py

Every null in this pipeline is computed within a family -- Rung 1 shuffles
stem/loop labels inside one family's own positions, Rung 3 deranges partners
inside one family's own stems -- and every reported aggregate is a mean or a
count over `per_rna`. Nothing pools across families. Correcting four annotations
therefore invalidates those four families' records and nothing else: the
aggregates recompute from stored values once the four are replaced.

This script uses the package's own stem parser and eligibility filter rather
than restating the criteria, and reports which of the four still qualify for
Rung 3 under their corrected annotations. That decides N, which
PREREGISTRATION_PHASE6_V2.md fixes at 32, so a change here is a registration
matter.

No model is run.
"""

import json
from pathlib import Path

from rna_structure_audit.rungs.rung3 import _get_eligible_pairs, _parse_stems

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data" / "rfam_families"
CORRECTED = ["mir_122_precursor", "mir_155_precursor",
             "mir_21_precursor", "mir_let7_precursor"]
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

MIN_WC_PAIRS = 15      # registration: "Families with fewer than 15 canonical WC pairs are excluded"
MIN_ELIGIBLE = 5       # registration: "If fewer than 5 eligible interior pairs ... exclude it"

# Every stored file that carries per-family records, by rung.
RESULT_GLOBS = {
    "Rungs 1-2 (mutation sensitivity, probing)": "*phases*1*5*.json",
    "Rung 3 (partner specificity)": "*phase6*.json",
}


def eligibility(name: str) -> dict:
    record = json.loads((FAMILIES / f"{name}.json").read_text())
    sequence = record["sequence"].upper()
    stems = _parse_stems(record["dot_bracket"], sequence)
    eligible = _get_eligible_pairs(sequence, stems)
    n_wc = sum(len(stem) for stem in stems)
    return dict(
        rfam_id=record.get("rfam_id"),
        length=record["length"],
        n_wc=n_wc,
        n_stems=sum(1 for stem in stems if len(stem) >= 3),
        n_eligible=len(eligible),
        qualifies=n_wc >= MIN_WC_PAIRS and len(eligible) >= MIN_ELIGIBLE,
    )


def stored_rung3_state() -> dict:
    """Which families the deposited Rung 3 records actually scored."""
    path = REPO / "results" / "rinalmo_phase6_ps.json"
    per = json.loads(path.read_text())["results"]["per_rna"]
    return {name: not (body.get("skipped") or body.get("best_ps") is None)
            for name, body in per.items() if isinstance(body, dict)}


def main() -> int:
    names = sorted(path.stem for path in FAMILIES.glob("*.json"))
    current = {name: eligibility(name) for name in names}
    stored = stored_rung3_state()

    print("The four corrected families, under their new annotations\n")
    header = (f"{'family':<22} {'rfam':>8} {'len':>5} {'WC':>4} {'stems':>6} "
              f"{'elig':>5} {'qualifies':>10} {'scored before':>14}")
    print(header)
    print("-" * len(header))
    for name in CORRECTED:
        row = current[name]
        print(f"{name:<22} {row['rfam_id']:>8} {row['length']:>5} {row['n_wc']:>4} "
              f"{row['n_stems']:>6} {row['n_eligible']:>5} "
              f"{'yes' if row['qualifies'] else 'no':>10} "
              f"{'yes' if stored.get(name) else 'no':>14}")

    entering = [n for n in CORRECTED if current[n]["qualifies"] and not stored.get(n)]
    leaving = [n for n in CORRECTED if not current[n]["qualifies"] and stored.get(n)]
    n_before = sum(1 for n, was in stored.items() if was and n not in QUARANTINE)
    n_after = n_before + len(entering) - len(leaving)

    print(f"\nRung 3 confirmatory set: N = {n_before} before, N = {n_after} after")
    if entering:
        print(f"  entering: {', '.join(entering)}")
    if leaving:
        print(f"  leaving:  {', '.join(leaving)}")
    if n_after != n_before:
        print("  PREREGISTRATION_PHASE6_V2.md fixes N = 32. A change to the "
              "confirmatory set\n  is a registration matter, not an arithmetic one.")

    print("\nStored files carrying per-family records that the four invalidate\n")
    total = 0
    for label, pattern in RESULT_GLOBS.items():
        files = sorted(set(REPO.glob(f"results/**/{pattern}"))
                       | set(REPO.glob(f"data/gpu_results/**/{pattern}")))
        print(f"  {label}: {len(files)} files")
        total += len(files)
    print(f"\n  {total} stored files x 4 families = {total * 4} family-model "
          f"evaluations to repeat, against")
    print(f"  {total} x {len(names)} = {total * len(names)} in a full re-run. The other {len(names) - 4} families are\n  unchanged and their stored records stay valid.")
    print("\n  This is an upper bound: the count includes audit copies of the "
          "same run and the\n  untrained and synthetic controls, several of "
          "which cover the same model.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
