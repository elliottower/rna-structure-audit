"""The measured resolution floor, H3's real chance rate, and the neighbour asymmetry.

Three quantities this pipeline used to assume, estimate from a proxy, or discard:

`noop_floor` runs one unmutated sequence through the model twice and computes the
same statistic as `best_ps`, whose true value is then exactly zero. What comes
back is the floor the measurement can resolve, per model and per layer. A model
whose forward pass is deterministic returns the cosine's own rounding, around
1e-16; anything larger is the model disagreeing with itself.

`h3_chance_fraction` is the fraction of *deranged* pairs whose assigned partner
still beats both its neighbours. The pairing is false by construction, so this is
H3's chance rate for that model on that family -- matched to the same weights,
sequence and layer, which an untrained control is not (D16). The registration
uses 1/3.

`nearer_neighbor_larger` is how often the neighbour nearer the mutated position
carries the larger perturbation. If perturbation decays with sequence distance,
that is systematically above 0.5, and the decay explains why the chance rate sits
below 1/3. If it is at 0.5, the explanation is wrong and the deficit is
something else.

    uv run --no-project --python 3.12 python scripts/report_resolution_and_chance.py \
        --results results/repaired_panel_v3
"""

import argparse
import json
from pathlib import Path
from statistics import mean, median

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}


def collect(per_rna):
    floors, chance, nearer, observed = [], [], [], []
    for name, entry in per_rna.items():
        if name in QUARANTINED or entry.get("best_ps") is None:
            continue
        floor = entry.get("noop_floor") or {}
        if floor.get("best_ps_noop") is not None:
            floors.append(abs(floor["best_ps_noop"]))
        if entry.get("h3_chance_fraction") is not None:
            chance.append(entry["h3_chance_fraction"])
        if entry.get("nearer_neighbor_larger") is not None:
            nearer.append(entry["nearer_neighbor_larger"])
        if isinstance(entry.get("h3_precision"), dict):
            observed.append(entry["h3_precision"]["fraction"])
    return floors, chance, nearer, observed


def show(value, fmt="{:.3f}"):
    return "--" if value is None else fmt.format(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=str(REPO / "results/repaired_panel_v3"))
    args = parser.parse_args()
    root = Path(args.results)

    print(f"  {'model':22s} {'floor (median)':>15s} {'floor (max)':>12s} "
          f"{'H3 chance':>10s} {'H3 observed':>12s} {'nearer larger':>14s}")
    for directory in sorted(root.iterdir()):
        path = directory / f"{directory.name}_phase6_ps.json"
        if not path.exists():
            continue
        floors, chance, nearer, observed = collect(
            json.loads(path.read_text())["results"]["per_rna"])
        if not observed:
            continue
        print(f"  {directory.name:22s} "
              f"{show(median(floors) if floors else None, '{:.2e}'):>15s} "
              f"{show(max(floors) if floors else None, '{:.2e}'):>12s} "
              f"{show(mean(chance) if chance else None):>10s} "
              f"{show(mean(observed)):>12s} "
              f"{show(mean(nearer) if nearer else None):>14s}")
    print("\n  A deterministic forward pass floors near 1e-16. The registration "
          "puts H3's chance rate at 0.333.\n  Distance decay predicts 'nearer "
          "larger' above 0.5; 0.5 means the deficit is something else.")


if __name__ == "__main__":
    main()
