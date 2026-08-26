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
from scipy.stats import spearmanr

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


def htt_macros(rows: list[dict]) -> str:
    """The counts the HTT prose quotes, so no literal is typed twice."""
    high = [row for row in rows if row["probing"] > 0.97]
    clean = [row for row in high if row["exceeding"] <= 1]
    monotone = [row for row in rows if row["rho_embedding"] > 0.99]
    return "\n".join([
        rf"\newcommand{{\httModels}}{{{len(rows)}}}",
        rf"\newcommand{{\httHighProbing}}{{{len(high)}}}",
        rf"\newcommand{{\httHighProbingClean}}{{{len(clean)}}}",
        rf"\newcommand{{\httMonotone}}{{{len(monotone)}}}",
    ]) + "\n"


HTT = REPO / "results" / "repaired_panel_v3"

HTT_MODELS = [
    ("rinalmo", "RiNALMo"), ("ernierna", "ERNIE-RNA"), ("rnafm", "RNA-FM"),
    ("utrlm", "UTR-LM"), ("splicebert", "SpliceBERT"), ("caduceus", "Caduceus"),
    ("hyenadna", "HyenaDNA"), ("nt", "NT~v2"), ("dnabert2", "DNABERT-2"),
]


def htt_row(key: str, label: str) -> dict:
    """One model's five CAG-repeat fragments, from the repaired-panel HTT run."""
    path = HTT / f"{key}_htt" / f"{key}_phases_1_to_5.json"
    if not path.exists():
        raise ArmMismatch(f"missing {path.relative_to(REPO)}")
    stored = json.loads(path.read_text())
    mutation = stored["mutation_trained"]
    attention = stored["attention_trained"].get("per_rna", {})
    rhos = [body["best_corr"] for body in attention.values()
            if isinstance(body, dict) and body.get("best_corr") is not None]
    distances = stored.get("htt_distances")
    if distances is None:
        raise ArmMismatch(f"{key}: the HTT run stored no embedding distances")
    return {
        "key": key, "label": label,
        "ratio": mutation["mean_best_ratio"],
        "exceeding": mutation["n_exceeding_nuc_null"],
        "scored": mutation["n_families_scored"],
        "attention": float(np.mean(rhos)) if rhos else None,
        "probing": stored["probing"]["best_accuracy"],
        # Distance from the CAG17 fragment at the final layer, which is where
        # the measure was originally read. No layer is selected: the embedding
        # layer already tracks repeat count at rho = +1.00 for every model,
        # because the fragments differ in length and composition, so a layer
        # chosen by that correlation is chosen among ties.
        "cag60_last": distances["at_last_layer"][-1],
        "cag60_embedding": distances["per_layer"]["0"][-1],
        "n_layers": distances["n_layers"],
        "rho_embedding": float(spearmanr(distances["repeat_counts"],
                                         distances["per_layer"]["0"]).statistic),
    }


def htt_table(rows: list[dict]) -> str:
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Five HTT exon~1 fragments with CAG repeat counts 17, 21, 36, "
        r"40 and 60, folded with ViennaRNA. Ratio and $>$~null are "
        r"stem--loop mutation sensitivity against the nucleotide-stratified "
        r"null; Attn~$\rho$ is the mean Spearman correlation between "
        r"self-attention and the contact map, absent for the state-space "
        r"models; Probing is linear-probe accuracy for stem-versus-loop "
        r"classification. CAG60 is cosine distance from the CAG17 fragment, at "
        r"the final layer and at the embedding layer. The repeat region contains "
        r"no U, so a linear probe has an extreme composition bias available to "
        r"it, and the fragments differ in length, so embedding distance rises "
        r"with repeat count at the embedding layer for every model "
        r"($\rho = +1.00$) before any block has run.}",
        r"\label{tab:htt}",
        r"\small",
        r"\begin{tabular}{@{}lrcrrrr@{}}",
        r"\toprule",
        r"\textbf{Model} & \textbf{Ratio} & \textbf{$>$ null} & "
        r"\textbf{Attn $\rho$} & \textbf{Probing} & "
        r"\textbf{CAG60 last} & \textbf{CAG60 best} \\",
        r"\midrule",
    ]
    for row in rows:
        attention = "---" if row["attention"] is None else f"{row['attention']:.3f}"
        lines.append(
            f"{row['label']} & {row['ratio']:.3f} & "
            f"{row['exceeding']}/{row['scored']} & {attention} & "
            f"{row['probing']:.3f} & {render_ps(row['cag60_last'])} & "
            f"{render_ps(row['cag60_embedding'])} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


MULTISEQ_MODELS = [
    ("rnafm", "RNA-FM"), ("evo", "Evo"), ("ernierna", "ERNIE-RNA"),
    ("caduceus", "Caduceus"), ("hyenadna", "HyenaDNA"), ("rinalmo", "RiNALMo"),
    ("nt", "NT~v2"), ("utrlm", "UTR-LM"), ("splicebert", "SpliceBERT"),
    ("dnabert2", "DNABERT-2"),
]


def multiseq_row(key: str, label: str) -> dict:
    path = PANEL / f"{key}_multiseq" / f"{key}_multiseq.json"
    if not path.exists():
        raise ArmMismatch(f"missing {path.relative_to(REPO)}")
    stored = json.loads(path.read_text())
    summary = stored["summary"]
    return {
        "key": key, "label": label,
        "mean_ratio": summary["mean_ratio"],
        "mean_cv": summary["mean_cv"],
        "median_cv": summary["median_cv"],
        "n_families": summary["n_families_scored"],
        "n_sequences": summary["n_total_sequences"],
        "excluded": [e["name"] for e in stored.get("excluded_for_ambiguity", [])],
        "commit": stored["stamp"]["commit"],
    }


def multiseq_table(rows: list[dict]) -> str:
    """Within-family variance of the Rung 1 ratio across seed replicates."""
    sequences = {row["n_sequences"] for row in rows}
    families = {row["n_families"] for row in rows}
    if len(sequences) != 1 or len(families) != 1:
        raise ArmMismatch(
            f"models scored different sets: {sorted(sequences)} sequences over "
            f"{sorted(families)} families; within-family variance is only "
            f"comparable across a common set")
    excluded = {tuple(sorted(row["excluded"])) for row in rows}
    if len(excluded) != 1:
        raise ArmMismatch(f"models excluded different sequences: {excluded}")
    ordered = sorted(rows, key=lambda row: -row["mean_ratio"])
    dropped = ", ".join(f"\\texttt{{{name.replace('_', chr(92) + '_')}}}"
                        for name in sorted(next(iter(excluded))))
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{Multi-sequence replication. Each family contributes 3--5 "
        rf"members of its Rfam seed alignment, giving {sorted(sequences)[0]} "
        rf"sequences over {sorted(families)[0]} families, and every model is "
        rf"scored on the same set. CV is the within-family coefficient of "
        rf"variation of the stem--loop sensitivity ratio, averaged over "
        rf"families. Three sequences carrying IUPAC ambiguity codes are "
        rf"excluded for every model ({dropped}). A low CV means the ratio is a "
        rf"property of the family rather than of the one curated "
        rf"representative.}}",
        r"\label{tab:multiseq}",
        r"\begin{tabular}{@{}lrrr@{}}",
        r"\toprule",
        r"\textbf{Model} & \textbf{Mean ratio} & \textbf{Mean CV} & "
        r"\textbf{Median CV} \\",
        r"\midrule",
    ]
    for row in ordered:
        lines.append(f"{row['label']} & {row['mean_ratio']:.3f} & "
                     f"{row['mean_cv']:.3f} & {row['median_cv']:.3f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def multiseq_macros(rows: list[dict]) -> str:
    by = {row["key"]: row for row in rows}
    stable = sorted((r for r in rows if r["mean_cv"] < 0.06),
                    key=lambda r: r["mean_cv"])
    lines = [
        rf"\newcommand{{\multiseqN}}{{{rows[0]['n_sequences']}}}",
        rf"\newcommand{{\multiseqFamilies}}{{{rows[0]['n_families']}}}",
        rf"\newcommand{{\multiseqStable}}{{{len(stable)}}}",
        rf"\newcommand{{\multiseqStableList}}{{"
        + ", ".join(f"{r['label']} ({r['mean_cv']:.3f})" for r in stable) + r"}",
    ]
    for key, macro in (("rinalmo", "RiNALMo"), ("ernierna", "ErnieRNA"),
                       ("rnafm", "RNAFM"), ("evo", "Evo")):
        lines.append(rf"\newcommand{{\cv{macro}}}{{{by[key]['mean_cv']:.3f}}}")
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

    htt_rows = [htt_row(key, label) for key, label in HTT_MODELS]
    multiseq_rows = [multiseq_row(key, label) for key, label in MULTISEQ_MODELS]

    GENERATED.mkdir(parents=True, exist_ok=True)
    (GENERATED / "synthetic_table.tex").write_text(table(rows))
    (GENERATED / "htt_table.tex").write_text(htt_table(htt_rows))
    (GENERATED / "multiseq_table.tex").write_text(multiseq_table(multiseq_rows))
    (GENERATED / "arm_macros.tex").write_text(
        macros(rows) + htt_macros(htt_rows) + multiseq_macros(multiseq_rows))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rows, indent=2))
    (OUT.parent / "htt_arm.json").write_text(json.dumps(htt_rows, indent=2))
    (OUT.parent / "multiseq_arm.json").write_text(json.dumps(multiseq_rows, indent=2))
    print(f"\n  {'model':12s}{'ratio':>8s}{'meanCV':>9s}{'medCV':>8s}")
    for row in sorted(multiseq_rows, key=lambda r: -r["mean_ratio"]):
        print(f"  {row['label']:12s}{row['mean_ratio']:8.3f}{row['mean_cv']:9.3f}"
              f"{row['median_cv']:8.3f}")
    print(f"\n  {'model':12s}{'ratio':>8s}{'>null':>7s}{'probe':>8s}"
          f"{'CAG60 last':>12s}{'CAG60 emb':>12s}{'layers':>8s}{'rho L0':>7s}")
    for row in htt_rows:
        print(f"  {row['label']:12s}{row['ratio']:8.3f}"
              f"{str(row['exceeding']) + '/' + str(row['scored']):>7s}"
              f"{row['probing']:8.3f}{row['cag60_last']:12.3e}"
              f"{row['cag60_embedding']:12.3e}"
              f"{row['n_layers']:8d}{row['rho_embedding']:+7.2f}")

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
