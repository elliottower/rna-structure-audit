"""The registered verdicts for the ERNIE-RNA pairwise-bias ablation.

Implements `preregistration/PREREGISTRATION_ERNIERNA_BIAS_RUNG3.md`, frozen at
d095d09 before either ablated cell was computed. Nothing here is chosen after the
fact: the criteria, the analysis set, the bootstrap and the void conditions are
all in that file.

Supersedes `evaluate_ernierna_ablation.py`, which printed its verdicts and wrote
no file, and read the chance rate over every eligible pair. This reads the four
arms of `results/decoy_panel_v1`, takes each arm's chance rate on the pairs its
precision is read on (the deviation dated 2026-10-11 in `DEVIATIONS.md`), keeps
the registered criteria, bootstrap and seed, and writes the whole result, the
table and the macros the manuscript quotes.

    uv run --no-project --with numpy python scripts/evaluate_ernierna_ablation_v24.py
"""

import hashlib
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results/decoy_panel_v1"
OUT = REPO / "results" / "audit" / "ernierna_ablation_v24.json"
GENERATED = REPO / "paper" / "generated_v24"
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
        if name in QUARANTINED or entry.get("h3_chance_fraction_matched") is None:
            continue
        out[name] = (entry["h3_precision"]["fraction"],
                     entry["h3_chance_fraction_matched"])
    return out


class AblationVoid(RuntimeError):
    """A registered void condition holds, so no hypothesis is evaluated."""


LABELS = {
    "trained intact": ("Trained, bias intact", "TrainedIntact"),
    "trained ablated": ("Trained, bias ablated", "TrainedAblated"),
    "untrained intact": ("Random-init, bias intact", "UntrainedIntact"),
    "untrained ablated": ("Random-init, bias ablated", "UntrainedAblated"),
}


def main():
    cells = {label: load(*where) for label, where in ARMS.items()}
    digests = {}
    for label, cell in cells.items():
        digests[label] = hashlib.sha256(",".join(sorted(cell)).encode()).hexdigest()[:16]
        if len(cell) < VOID_BELOW:
            raise AblationVoid(f"{label} carries a chance rate for {len(cell)} families, "
                               f"below the registered floor of {VOID_BELOW}")

    shared = sorted(set.intersection(*(set(c) for c in cells.values())))

    stats = {label: (np.array([cell[n][0] for n in shared]),
                     np.array([cell[n][1] for n in shared]))
             for label, cell in cells.items()}
    arms = {label: {"n": len(shared),
                    "precision": float(precision.mean()),
                    "chance_matched": float(chance.mean()),
                    "excess": float(precision.mean() - chance.mean()),
                    "family_set_digest": digests[label]}
            for label, (precision, chance) in stats.items()}

    check = {arm: bool(arms[f"{arm} ablated"]["precision"] < arms[f"{arm} intact"]["precision"])
             for arm in ("trained", "untrained")}
    if not all(check.values()):
        raise AblationVoid(f"the ablation did not lower precision in both arms: {check}")

    rng = np.random.default_rng(SEED)
    picks = [rng.integers(0, len(shared), size=len(shared))
             for _ in range(N_BOOTSTRAP)]

    def excess(label, pick):
        precision, chance = stats[label]
        return precision[pick].mean() - chance[pick].mean()

    def interval(fn):
        draws = np.array([fn(pick) for pick in picks])
        low, high = np.percentile(draws, [2.5, 97.5])
        return {"point": float(fn(np.arange(len(shared)))), "ci_low": float(low),
                "ci_high": float(high), "holds": bool(low > 0)}

    tests = {
        "H1": interval(lambda p: excess("untrained intact", p)
                       - excess("untrained ablated", p)),
        "H2": interval(lambda p: ((excess("untrained intact", p)
                                   - excess("untrained ablated", p))
                                  - (excess("trained intact", p)
                                     - excess("trained ablated", p)))),
        "H3": interval(lambda p: excess("trained ablated", p)),
    }

    # Written before anything is formatted, so a formatting error cannot lose it.
    OUT.write_text(json.dumps({
        "registration": "preregistration/PREREGISTRATION_ERNIERNA_BIAS_RUNG3.md",
        "panel": str(RESULTS.relative_to(REPO)),
        "chance": "h3_chance_fraction_matched",
        "n_bootstrap": N_BOOTSTRAP, "seed": SEED,
        "registered_family_digest": REGISTERED_DIGEST,
        "families": shared, "arms": arms,
        "manipulation_check": check, "tests": tests}, indent=2) + "\n")

    rows = [f"{LABELS[label][0]:<28} & {arm['precision']:.3f} & {arm['chance_matched']:.3f} "
            f"& ${arm['excess']:+.3f}$ \\\\" for label, arm in arms.items()]
    (GENERATED / "ablation_table.tex").write_text(f"""% Written by scripts/evaluate_ernierna_ablation_v24.py. Do not edit.
\\begin{{table}}[htbp]
\\centering
\\caption{{ERNIE-RNA with and without its pairwise bias buffer, trained and
randomly initialized. Per-pair precision is the fraction of deep-interior pairs
whose partner is perturbed more than either of its stem neighbors. Chance is the
same fraction under a within-stem derangement on those pairs, where the assigned
partner is false by construction, and excess is their difference. All four arms
cover the same {len(shared)} families. Ablation zeroes
\\texttt{{pairwise\\_bias\\_map}} and \\texttt{{pairwise\\_bias\\_proj}} and sets the
rebuild guard, without which the first forward pass restores the buffer.}}
\\label{{tab:ablation}}
\\begin{{tabular}}{{@{{}}lccc@{{}}}}
\\toprule
\\textbf{{Arm}} & \\textbf{{Precision}} & \\textbf{{Chance}} & \\textbf{{Excess}} \\\\
\\midrule
{chr(10).join(rows)}
\\bottomrule
\\end{{tabular}}
\\end{{table}}
""")
    macros = ["% Written by scripts/evaluate_ernierna_ablation_v24.py. Do not edit.",
              rf"\newcommand{{\ablationN}}{{{len(shared)}}}"]
    for label, arm in arms.items():
        name = LABELS[label][1]
        macros += [rf"\newcommand{{\ablationExcess{name}}}{{{arm['excess']:+.3f}}}",
                   rf"\newcommand{{\ablationPrecision{name}}}{{{arm['precision']:.3f}}}",
                   rf"\newcommand{{\ablationChance{name}}}{{{arm['chance_matched']:.3f}}}"]
    for key, name in (("H1", "ControlLoss"), ("H2", "LossDifference"), ("H3", "TrainedRemaining")):
        test = tests[key]
        macros += [rf"\newcommand{{\ablation{name}}}{{{test['point']:+.3f}}}",
                   rf"\newcommand{{\ablation{name}CI}}{{[{test['ci_low']:+.3f}, {test['ci_high']:+.3f}]}}"]
    trained_cost = arms["trained intact"]["excess"] - arms["trained ablated"]["excess"]
    macros += [rf"\newcommand{{\ablationTrainedCost}}{{{trained_cost:.3f}}}", ""]
    (GENERATED / "ablation_macros.tex").write_text("\n".join(macros))

    for label, arm in arms.items():
        print(f"  {label:20s} precision={arm['precision']:.3f} chance={arm['chance_matched']:.3f} "
              f"excess={arm['excess']:+.3f}")
    for key, test in tests.items():
        print(f"  {key} {test['point']:+.3f} [{test['ci_low']:+.3f}, {test['ci_high']:+.3f}] "
              f"{'holds' if test['holds'] else 'fails'}")
    print(f"  {len(shared)} families; wrote {OUT.relative_to(REPO)}, ablation_table.tex "
          "and ablation_macros.tex")


if __name__ == "__main__":
    main()
