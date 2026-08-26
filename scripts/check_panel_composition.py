"""Which printed dataset-description numbers move under the corrected annotations?

Run:  uv run --no-project --python 3.12 python scripts/check_panel_composition.py

paper_v12.tex:215-234 describes the evaluation panel: 52 families, a length
range and median with IQR, a class breakdown, and the assertion that every
family's sequence is drawn from an Rfam seed alignment with an experimentally
supported structure. Commit 5fc8914 replaced four annotations, three of which
had been one 5S rRNA sequence under three microRNA names, so the printed length
statistics were computed over a panel that is not the deposited one.

This script recomputes the length statistics from data/rfam_families at the
current checkout and at the pre-fix tree, and reports the source and accession
field of every record so the seed-alignment claim can be checked rather than
assumed. No model is run.
"""

import json
import statistics
import subprocess
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PRE_FIX = "891d6af"
CHANGED = ["mir_122_precursor", "mir_155_precursor",
           "mir_21_precursor", "mir_let7_precursor"]

# paper_v12.tex:219
PRINTED = dict(low=22, high=301, median=89, q1=56, q3=132)


def records_at(ref: str | None) -> dict[str, dict]:
    if ref is None:
        return {p.stem: json.loads(p.read_text())
                for p in sorted((REPO / "data/rfam_families").glob("*.json"))}
    names = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", ref, "--", "data/rfam_families"],
        cwd=REPO, capture_output=True, text=True, check=True).stdout.split()
    out = {}
    for name in names:
        blob = subprocess.run(["git", "show", f"{ref}:{name}"],
                              cwd=REPO, capture_output=True, text=True, check=True)
        out[Path(name).stem] = json.loads(blob.stdout)
    return out


def lengths(records: dict[str, dict]) -> list[int]:
    return sorted(len(r["sequence"]) for r in records.values())


def quartiles(values: list[int]) -> tuple[float, float, float]:
    q1, median, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return q1, median, q3


def describe(label: str, records: dict[str, dict]) -> None:
    values = lengths(records)
    q1, median, q3 = quartiles(values)
    print(f"  {label:<12} n={len(values)}  range {values[0]}-{values[-1]}  "
          f"median {median:g}  IQR {q1:g}-{q3:g}")


def main() -> None:
    now = records_at(None)
    before = records_at(PRE_FIX)

    print("Length statistics (paper_v12.tex:219 prints "
          f"{PRINTED['low']}-{PRINTED['high']} nt, median {PRINTED['median']}, "
          f"IQR {PRINTED['q1']}-{PRINTED['q3']})\n")
    describe("pre-fix", before)
    describe("corrected", now)

    values = lengths(now)
    q1, median, q3 = quartiles(values)
    moved = [name for name, (was, is_) in
             (("range", ((before and lengths(before)[0], lengths(before)[-1]),
                         (values[0], values[-1]))),
              ("median", (quartiles(lengths(before))[1], median)),
              ("IQR", (tuple(round(x, 1) for x in quartiles(lengths(before))[::2]),
                       (round(q1, 1), round(q3, 1)))))
             if was != is_]
    print(f"\n  moves under the correction: {', '.join(moved) if moved else 'nothing'}")

    print("\nThe four corrected families\n")
    for name in CHANGED:
        was, is_ = before[name], now[name]
        print(f"  {name}")
        print(f"    before  {was['rfam_id']}  {len(was['sequence'])} nt  "
              f"source: {was.get('source')}")
        print(f"    after   {is_['rfam_id']}  {len(is_['sequence'])} nt  "
              f"source: {is_.get('source')}")

    print("\nsource field across the deposited panel\n")
    for source, count in Counter(
            r.get("source", "(absent)") for r in now.values()).most_common():
        print(f"  {count:>3}  {source}")

    unverifiable = sorted(name for name, r in now.items()
                          if "seed alignment member" not in r.get("source", ""))
    print(f"\n  records naming a seed-alignment member: "
          f"{len(now) - len(unverifiable)} of {len(now)}")
    print(f"  records with no accession behind them: {len(unverifiable)}\n")
    for name in unverifiable:
        print(f"    {name:<28} {str(now[name]["rfam_id"]):<9} "
              f"{len(now[name]['sequence']):>4} nt   {now[name].get('source')}")
    print("\n  paper_v12.tex:229 states every sequence is the highest-bit-score "
          "seed sequence\n  of its family. For these records the claim cannot be "
          "checked from data/ alone.\n  Three of the four families corrected in "
          "5fc8914 named a seed member that was\n  the wrong family's, and the "
          "fourth named no accession at all.")


if __name__ == "__main__":
    main()
