"""Is a family's PS resolvable against its null, or are both at the float32 floor?

Rung 3 scores a pair by cosine distance between two embeddings of sequences
differing at one nucleotide. When a model's representation barely moves, both
distances round to nearly the same float32 and `1 - cos_sim` is computed by
catastrophic cancellation: the difference of two numbers close to 1, carrying
almost no significant digits.

`exceeds_null_primary` then compares one such number against a null percentile
built from the same quantity, and H1's gate counts the families where it holds.
This reports, per model, how many families are decided at a magnitude float32
cannot represent meaningfully, so the count can be read for what it is.

The threshold is float32 epsilon, 1.19e-07, scaled to the operands: a cosine
distance below it carries no significant digits at all, and one within an order
of magnitude of it carries roughly one.

    uv run --no-project --python 3.12 python scripts/audit_ps_noise_floor.py \
        --results results/repaired_panel
"""

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# numpy.finfo(numpy.float32).eps
FLOAT32_EPS = 1.1920929e-07


def audit(per_rna):
    """Families split by whether their PS decision is above the float32 floor."""
    resolvable, at_floor, exceeding_at_floor, exceeding_total = 0, 0, 0, 0
    for name, entry in per_rna.items():
        if name in QUARANTINED or entry.get("best_ps") is None:
            continue
        ps = abs(entry["best_ps"])
        null = abs(entry.get("null_95th_primary") or 0.0)
        decided_at_floor = max(ps, null) < 10 * FLOAT32_EPS
        at_floor += decided_at_floor
        resolvable += not decided_at_floor
        if entry.get("exceeds_null_primary"):
            exceeding_total += 1
            exceeding_at_floor += decided_at_floor
    return {"resolvable": resolvable, "at_floor": at_floor,
            "exceeding_total": exceeding_total,
            "exceeding_at_floor": exceeding_at_floor}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=str(REPO / "results/repaired_panel"))
    args = parser.parse_args()

    root = Path(args.results)
    print(f"Float32 epsilon {FLOAT32_EPS:.2e}; a family is 'at the floor' when "
          f"both its PS and its null percentile fall below {10 * FLOAT32_EPS:.2e}\n")
    print(f"  {'model':14s} {'scored':>7s} {'resolvable':>11s} {'at floor':>9s} "
          f"{'over null':>10s} {'of those, at floor':>19s}")
    for directory in sorted(root.iterdir()):
        path = directory / f"{directory.name}_phase6_ps.json"
        if not path.exists():
            continue
        stats = audit(json.loads(path.read_text())["results"]["per_rna"])
        scored = stats["resolvable"] + stats["at_floor"]
        if not scored:
            continue
        print(f"  {directory.name:14s} {scored:7d} {stats['resolvable']:11d} "
              f"{stats['at_floor']:9d} {stats['exceeding_total']:10d} "
              f"{stats['exceeding_at_floor']:19d}")
    print()


if __name__ == "__main__":
    main()
