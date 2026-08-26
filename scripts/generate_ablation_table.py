"""The four-arm ERNIE-RNA bias ablation, as a LaTeX table.

Written from `results/repaired_panel_v3` rather than typed, for the same reason
every other table here is generated: a hand-carried number is a defect waiting
for the next re-run.

The arms are trained and randomly initialized, each with the pairwise bias buffer
intact and ablated. Excess is precision minus that arm's own derangement chance
rate, because a chance rate is a property of a model and the four arms are four
models.

    uv run --no-project --with "numpy<2" --python 3.12 \\
        python scripts/generate_ablation_table.py
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results" / "repaired_panel_v3"
OUT = REPO / "paper" / "generated" / "ablation_table.tex"
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

ARMS = [
    ("Trained, bias intact", "ernierna", "ernierna"),
    ("Trained, bias ablated", "ernierna_noattnbias", "ernierna"),
    ("Random-init, bias intact", "ernierna_untrained", "ernierna_untrained"),
    ("Random-init, bias ablated", "ernierna_untrained_noattnbias",
     "ernierna_untrained"),
]


def arm(directory, prefix):
    per_rna = json.loads(
        (PANEL / directory / f"{prefix}_phase6_ps.json").read_text()
    )["results"]["per_rna"]
    kept = {name: e for name, e in per_rna.items()
            if name not in QUARANTINED and e.get("h3_chance_fraction") is not None}
    precision = sum(e["h3_precision"]["fraction"] for e in kept.values()) / len(kept)
    chance = sum(e["h3_chance_fraction"] for e in kept.values()) / len(kept)
    return len(kept), precision, chance


def main():
    rows, n = [], None
    for label, directory, prefix in ARMS:
        count, precision, chance = arm(directory, prefix)
        n = count if n is None else n
        if count != n:
            raise ValueError(f"{label} covers {count} families, not {n}; the arms "
                             "must share an analysis set for the comparison to pair")
        rows.append(f"{label:<28} & {precision:.3f} & {chance:.3f} "
                    f"& ${precision - chance:+.3f}$ \\\\")

    OUT.write_text(f"""% Written by scripts/generate_ablation_table.py. Do not edit.
\\begin{{table}}[htbp]
\\centering
\\caption{{ERNIE-RNA with and without its pairwise bias buffer, trained and
randomly initialized. Per-pair precision is the fraction of eligible pairs whose
partner is perturbed more than either of its stem neighbors. Chance is the same
fraction under a within-stem derangement, where the assigned partner is false by
construction, and excess is their difference. All four arms cover the same {n}
families. Ablation zeroes \\texttt{{pairwise\\_bias\\_map}} and
\\texttt{{pairwise\\_bias\\_proj}} and sets the rebuild guard, without which the
first forward pass restores the buffer.}}
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
    print(f"  wrote {OUT.relative_to(REPO)}, {n} families per arm")


if __name__ == "__main__":
    main()
