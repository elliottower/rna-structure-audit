"""Tables for the arms that vary the sequences rather than the model.

The synthetic arm keeps each family's structure and reassigns its nucleotides,
so the sequences carry no evolutionary covariation. Reading it off mean PS
alone understates what survives: PS is a magnitude and the synthetic sequences
are less structured in composition, while per-pair precision asks the question
the rung is about -- does a mutation perturb its partner more than that
partner's neighbours -- and is measured against each arm's own derangement
chance rate.

    uv run --no-project --with "numpy<2" --python 3.12 \
        python scripts/generate_arm_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results" / "repaired_panel_v3"
GENERATED = REPO / "paper" / "generated"
OUT = REPO / "results" / "audit" / "synthetic_arm.json"

QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# Ordered as Table 3 is, by natural mean PS.
MODELS = [
    ("rinalmo", "RiNALMo", "650M"),
    ("ernierna", "ERNIE-RNA", "86M"),
    ("caduceus", "Caduceus", "14M"),
    ("evo", "Evo", "7B"),
    ("hyenadna", "HyenaDNA", "5.4M"),
    ("splicebert", "SpliceBERT", "19M"),
    ("rnafm", "RNA-FM", "99M"),
    ("utrlm", "UTR-LM", r"$\sim$2M"),
    ("nt", "NT~v2", "56M"),
    ("dnabert2", "DNABERT-2", "117M"),
]


class ArmMismatch(RuntimeError):
    """An arm is absent, or its two arms were not scored over the same panel."""


def read(directory: str, key: str) -> dict:
    path = PANEL / directory / f"{key}_phase6_ps.json"
    if not path.exists():
        raise ArmMismatch(f"missing {path.relative_to(REPO)}")
    stored = json.loads(path.read_text())
    per_rna = stored["results"]["per_rna"]

    scored = {name: body for name, body in per_rna.items()
              if name not in QUARANTINE and isinstance(body, dict)}
    ps = [body["best_ps"] for body in scored.values()
          if isinstance(body.get("best_ps"), (int, float))]
    rated = [body for body in scored.values()
             if body.get("h3_chance_fraction") is not None]
    if not rated:
        raise ArmMismatch(f"{directory}: no family carries a chance rate")
    precision = float(np.mean([b["h3_precision"]["fraction"] for b in rated]))
    chance = float(np.mean([b["h3_chance_fraction"] for b in rated]))
    return {"n_scored": len(ps), "n_rated": len(rated),
            "mean_ps": float(np.mean(ps)), "precision": precision,
            "chance": chance, "excess": precision - chance,
            # The first panel run predates the field; every arm that varies the
            # sequences was run after it was added, so an absent field is a
            # natural-sequence run and a present one is checked below.
            "sequences": stored["stamp"].get("sequences", "natural")}


def render_ps(value: float) -> str:
    """Small magnitudes in scientific notation, as Table 3 renders them."""
    if abs(value) >= 1e-3:
        return f"{value:.4f}"
    exponent = int(np.floor(np.log10(abs(value)))) if value else 0
    mantissa = value / (10 ** exponent)
    return rf"${mantissa:.1f} \times 10^{{{exponent}}}$"


def table(rows: list[dict]) -> str:
    natural_n = {row["natural"]["n_rated"] for row in rows}
    synthetic_n = {row["synthetic"]["n_rated"] for row in rows}
    if len(natural_n) != 1 or len(synthetic_n) != 1:
        raise ArmMismatch(f"arms scored over different panels: "
                          f"natural {natural_n}, synthetic {synthetic_n}")
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{Partner specificity on natural and synthetic sequences. "
        rf"Synthetic sequences keep each family's dot-bracket structure and "
        rf"reassign its nucleotides -- uniform Watson--Crick pairs at paired "
        rf"positions, uniform elsewhere -- five per family, so they carry no "
        rf"evolutionary covariation. Excess is per-pair precision minus that "
        rf"arm's own within-stem derangement chance rate, over "
        rf"{sorted(natural_n)[0]} natural and {sorted(synthetic_n)[0]} "
        rf"synthetic sequences. Mean PS falls further than excess precision "
        rf"does, because PS is a magnitude and the synthetic sequences are "
        rf"less structured in composition.}}",
        r"\label{tab:synthetic}",
        r"\small",
        r"\begin{tabular}{@{}lrrrrr@{}}",
        r"\toprule",
        r" & \multicolumn{2}{c}{\textbf{Mean PS}} & "
        r"\multicolumn{3}{c}{\textbf{Excess precision over own chance}} \\",
        r"\cmidrule(lr){2-3}\cmidrule(lr){4-6}",
        r"\textbf{Model} & \textbf{Natural} & \textbf{Synthetic} & "
        r"\textbf{Natural} & \textbf{Synthetic} & \textbf{Retained} \\",
        r"\midrule",
    ]
    for row in rows:
        natural, synthetic = row["natural"], row["synthetic"]
        # A ratio of two near-zero excesses is noise reported to three figures,
        # so it is printed only where the natural arm has something to retain.
        retained = (rf"{synthetic['excess'] / natural['excess']:.0%}".replace("%", r"\%")
                    if natural["excess"] >= 0.05 else "---")
        lines.append(
            f"{row['label']} ({row['size']}) & {render_ps(natural['mean_ps'])} & "
            f"{render_ps(synthetic['mean_ps'])} & {natural['excess']:+.3f} & "
            f"{synthetic['excess']:+.3f} & {retained} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def macros(rows: list[dict]) -> str:
    lines = ["% Written by scripts/generate_arm_tables.py. Do not edit."]
    by_key = {row["key"]: row for row in rows}
    for key, macro in (("rinalmo", "RiNALMo"), ("ernierna", "ErnieRNA")):
        row = by_key[key]
        natural, synthetic = row["natural"], row["synthetic"]
        drop = 1 - synthetic["mean_ps"] / natural["mean_ps"]
        lines += [
            rf"\newcommand{{\synPS{macro}}}{{{synthetic['mean_ps']:.4f}}}",
            rf"\newcommand{{\synPSDrop{macro}}}{{{drop:.0%}}}".replace("%", r"\%"),
            rf"\newcommand{{\synExcess{macro}}}{{{synthetic['excess']:+.3f}}}",
            rf"\newcommand{{\natExcess{macro}}}{{{natural['excess']:+.3f}}}",
            rf"\newcommand{{\synRetained{macro}}}"
            rf"{{{synthetic['excess'] / natural['excess']:.0%}}}".replace("%", r"\%"),
        ]
    lines.append(rf"\newcommand{{\synN}}{{{rows[0]['synthetic']['n_rated']}}}")
    lines.append(rf"\newcommand{{\synPerFamily}}{{5}}")
    return "\n".join(lines) + "\n"


def main() -> int:
    rows = []
    for key, label, size in MODELS:
        natural = read(key, key)
        synthetic = read(f"{key}_synthetic", key)
        if "synthetic" not in synthetic["sequences"]:
            raise ArmMismatch(
                f"{key}_synthetic ran on {synthetic['sequences']!r} sequences")
        rows.append({"key": key, "label": label, "size": size,
                     "natural": natural, "synthetic": synthetic})

    GENERATED.mkdir(parents=True, exist_ok=True)
    (GENERATED / "synthetic_table.tex").write_text(table(rows))
    (GENERATED / "arm_macros.tex").write_text(macros(rows))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rows, indent=2))

    print(f"  {'model':12s}{'PS nat':>10s}{'PS syn':>10s}"
          f"{'excess nat':>12s}{'excess syn':>12s}{'retained':>10s}")
    for row in rows:
        natural, synthetic = row["natural"], row["synthetic"]
        retained = (f"{synthetic['excess'] / natural['excess']:.0%}"
                    if natural["excess"] >= 0.05 else "--")
        print(f"  {row['label']:12s}{natural['mean_ps']:10.4f}"
              f"{synthetic['mean_ps']:10.4f}{natural['excess']:+12.3f}"
              f"{synthetic['excess']:+12.3f}{retained:>10s}")
    print(f"\n  wrote paper/generated/synthetic_table.tex, arm_macros.tex and "
          f"{OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
