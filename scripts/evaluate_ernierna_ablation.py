"""The registered verdicts for the ERNIE-RNA pairwise-bias ablation.

Implements `preregistration/PREREGISTRATION_ERNIERNA_BIAS_RUNG3.md`, frozen at
d095d09 before either ablated cell was computed. Nothing here is chosen after the
fact: the criteria, the analysis set, the bootstrap and the void conditions are
all in that file.

    uv run --no-project --with "numpy<2" --python 3.12 \\
        python scripts/evaluate_ernierna_ablation.py
"""

import hashlib
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results/repaired_panel_v3"
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
N_BOOTSTRAP = 10_000
SEED = 42
REGISTERED_DIGEST = "10bc3b36b3fe0d08"
VOID_BELOW = 33

ARMS = {
    "trained intact": ("ernierna", "ernierna"),
    "trained ablated": ("ernierna_noattnbias", "ernierna"),
    "untrained intact": ("ernierna_untrained", "ernierna_untrained"),
    "untrained ablated": ("ernierna_untrained_noattnbias", "ernierna_untrained"),
}


def load(directory, prefix):
    """Per-family precision and chance rate. The ablated cells name their files
    after the model rather than the directory."""
    path = RESULTS / directory / f"{prefix}_phase6_ps.json"
    per_rna = json.loads(path.read_text())["results"]["per_rna"]
    out = {}
    for name, entry in per_rna.items():
        if name in QUARANTINED or entry.get("h3_chance_fraction") is None:
            continue
        out[name] = (entry["h3_precision"]["fraction"], entry["h3_chance_fraction"])
    return out


def main():
    cells = {label: load(*where) for label, where in ARMS.items()}

    print("  Registered criteria from PREREGISTRATION_ERNIERNA_BIAS_RUNG3.md,")
    print("  frozen at d095d09 before either ablated cell existed.\n")

    for label, cell in cells.items():
        digest = hashlib.sha256(",".join(sorted(cell)).encode()).hexdigest()[:16]
        note = "" if digest == REGISTERED_DIGEST else "  (differs from registered set)"
        print(f"  {label:20s} n={len(cell):3d}  digest={digest[:12]}{note}")
        if len(cell) < VOID_BELOW:
            print(f"\n  VOID: {label} carries a chance rate for {len(cell)} "
                  f"families, below the registered floor of {VOID_BELOW}.")
            return

    shared = sorted(set.intersection(*(set(c) for c in cells.values())))
    print(f"\n  analysis set: {len(shared)} families common to all four arms\n")

    def arrays(label):
        cell = cells[label]
        return (np.array([cell[n][0] for n in shared]),
                np.array([cell[n][1] for n in shared]))

    stats = {}
    for label in cells:
        precision, chance = arrays(label)
        stats[label] = (precision, chance)
        print(f"  {label:20s} precision={precision.mean():.3f}  "
              f"chance={chance.mean():.3f}  excess={precision.mean() - chance.mean():+.3f}")

    print("\n  Manipulation check: ablated precision below intact, both arms")
    ok = True
    for arm in ("trained", "untrained"):
        intact = stats[f"{arm} intact"][0].mean()
        ablated = stats[f"{arm} ablated"][0].mean()
        held = ablated < intact
        ok &= held
        print(f"    {arm:10s} {intact:.3f} -> {ablated:.3f}   "
              f"{'holds' if held else 'FAILS'}")
    if not ok:
        print("\n  The ablation did not take effect. No hypothesis is evaluated.")
        return

    rng = np.random.default_rng(SEED)
    picks = [rng.integers(0, len(shared), size=len(shared))
             for _ in range(N_BOOTSTRAP)]

    def excess(label, pick):
        precision, chance = stats[label]
        return precision[pick].mean() - chance[pick].mean()

    def interval(fn):
        draws = np.array([fn(pick) for pick in picks])
        low, high = np.percentile(draws, [2.5, 97.5])
        whole = np.arange(len(shared))
        return fn(whole), low, high

    tests = {
        "H1  untrained loses excess":
            lambda p: excess("untrained intact", p) - excess("untrained ablated", p),
        "H2  trained loses less":
            lambda p: ((excess("untrained intact", p) - excess("untrained ablated", p))
                       - (excess("trained intact", p) - excess("trained ablated", p))),
        "H3  trained keeps excess":
            lambda p: excess("trained ablated", p),
    }

    print()
    for label, fn in tests.items():
        point, low, high = interval(fn)
        verdict = "PASS" if low > 0 else "fail"
        print(f"  {label:30s} {point:+.3f}  95% CI [{low:+.3f}, {high:+.3f}]  {verdict}")

    print("\n  Registered resolution is 0.126; an interval covering zero means an\n"
          "  effect smaller than this design resolves, not an absence of effect.")


if __name__ == "__main__":
    main()
