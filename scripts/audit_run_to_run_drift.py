"""How far a model's per-layer PS moves between two runs of the same code.

One model, RNA-FM, returned different numbers across two runs whose code path
for it differs by a comment. `per_layer_ps` is stored per family, so the drift
can be measured directly rather than inferred from whichever summary happened to
move.

The distinction that matters is between drift in the values and instability in
`argmax`. Quantities reported "at the best layer" -- the positive control among
them -- are conditioned on a layer chosen by argmax over the per-layer profile.
When layers are near-tied the argmax is decided by noise, and a quantity read at
layer 12 in one run and layer 5 in the next differs by whatever separates those
layers, which can be orders of magnitude while the profile itself barely moved.

    uv run --no-project --python 3.12 python scripts/audit_run_to_run_drift.py \
        --a results/repaired_panel --b results/repaired_panel_v2
"""

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}


def per_rna(path):
    return json.loads(path.read_text())["results"]["per_rna"]


def drift(a, b):
    """Largest per-layer movement, and how often the argmax moved with it."""
    worst, layer_moves, compared = 0.0, 0, 0
    for name, before in a.items():
        after = b.get(name)
        if (name in QUARANTINED or after is None
                or before.get("per_layer_ps") is None
                or after.get("per_layer_ps") is None):
            continue
        compared += 1
        for x, y in zip(before["per_layer_ps"], after["per_layer_ps"]):
            worst = max(worst, abs(x - y))
        layer_moves += before.get("best_layer") != after.get("best_layer")
    return worst, layer_moves, compared


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", default=str(REPO / "results/repaired_panel"))
    parser.add_argument("--b", default=str(REPO / "results/repaired_panel_v2"))
    args = parser.parse_args()
    a_root, b_root = Path(args.a), Path(args.b)

    print(f"  {'model':22s} {'families':>9s} {'largest per-layer move':>23s} "
          f"{'best_layer moved':>17s}")
    for directory in sorted(a_root.iterdir()):
        name = directory.name
        a_path = directory / f"{name}_phase6_ps.json"
        b_path = b_root / name / f"{name}_phase6_ps.json"
        if not (a_path.exists() and b_path.exists()):
            continue
        worst, moves, compared = drift(per_rna(a_path), per_rna(b_path))
        if not compared:
            continue
        print(f"  {name:22s} {compared:9d} {worst:23.3e} "
              f"{f'{moves}/{compared}':>17s}")
    print("\nA model whose forward pass is deterministic moves by exactly zero.")


if __name__ == "__main__":
    main()
