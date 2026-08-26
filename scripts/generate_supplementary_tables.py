"""Generate supplementary tables S1-S4 as LaTeX for the RNA journal submission."""

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FAM_DIR = ROOT / "data/rfam_families"

# The panel every table reads. `repaired_panel_v3` carries the Rung 3
# token-alignment repair, the float64 metric and the measured resolution floor,
# at one commit per cell. Before this, each entry below named a separate July run
# directory at an unrecorded commit, so the supplement and the main tables could
# disagree about the same model.
PANEL = ROOT / "results" / "repaired_panel_v3"

KEYS = {
    "RNA-FM": "rnafm",
    "RiNALMo": "rinalmo",
    "ERNIE-RNA": "ernierna",
    "SpliceBERT": "splicebert",
    "UTR-LM": "utrlm",
    "NT v2": "nt",
    "DNABERT-2": "dnabert2",
    "Caduceus": "caduceus",
    "HyenaDNA": "hyenadna",
    "Evo": "evo",
}

MODEL_FILES_PHASES = {
    label: PANEL / key / f"{key}_phases_1_to_5.json" for label, key in KEYS.items()
}

MODEL_FILES_PS = {
    label: PANEL / key / f"{key}_phase6_ps.json" for label, key in KEYS.items()
}


def load_families():
    families = []
    for f in sorted(os.listdir(FAM_DIR)):
        if not f.endswith(".json"):
            continue
        with open(FAM_DIR / f) as fh:
            d = json.load(fh)
        # Withdrawn records stay on disk and out of the panel; see DEVIATIONS.md.
        if "excluded" in d:
            continue
        seq = d["sequence"].upper()
        db = d["dot_bracket"]
        stem_pos = [i for i, c in enumerate(db) if c in "()"]
        loop_pos = [i for i, c in enumerate(db) if c == "."]
        if not stem_pos or not loop_pos:
            continue
        sgc = sum(1 for i in stem_pos if seq[i] in "GC") / len(stem_pos) * 100
        lgc = sum(1 for i in loop_pos if seq[i] in "GC") / len(loop_pos) * 100
        families.append({
            "name": d["name"],
            "rfam_id": d.get("rfam_id", "---"),
            "length": d.get("length", len(seq)),
            "n_stem": len(stem_pos),
            "n_loop": len(loop_pos),
            "stem_gc": sgc,
            "loop_gc": lgc,
        })
    return families


def table_s1(families):
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Supplementary Table S1: Rfam families evaluated. GC content computed from consensus sequence stem and loop positions. Length is the consensus sequence length in nucleotides.}")
    lines.append(r"\label{tab:s1}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{llrrrcc}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Family} & \textbf{Rfam ID} & \textbf{Length} & \textbf{Stem pos.} & \textbf{Loop pos.} & \textbf{Stem GC\%} & \textbf{Loop GC\%} \\")
    lines.append(r"\midrule")
    for fam in families:
        name = fam["name"].replace("_", r"\_")
        lines.append(
            f"{name} & {fam['rfam_id']} & {fam['length']} & {fam['n_stem']} & {fam['n_loop']} & {fam['stem_gc']:.1f} & {fam['loop_gc']:.1f} \\\\"
        )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def table_s2(families):
    """Per-family composition gate results across models."""
    all_data = {}
    for model, path in MODEL_FILES_PHASES.items():
        with open(path) as f:
            d = json.load(f)
        pr = d["mutation_trained"]["per_rna"]
        all_data[model] = pr

    fam_names = [f["name"] for f in families]
    models = list(MODEL_FILES_PHASES.keys())

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Supplementary Table S2: Per-family stem--loop sensitivity ratios and composition gate results. Ratio: best-layer ratio of mean stem perturbation to mean loop perturbation. Gate: whether the ratio exceeds the nucleotide-stratified null at the 95th percentile (\checkmark) or not (---).}")
    lines.append(r"\label{tab:s2}")
    lines.append(r"\tiny")
    header = r"\textbf{Family}"
    for m in models:
        header += f" & \\textbf{{{m}}}"
    header += r" \\"
    lines.append(r"\begin{tabular}{l" + "c" * len(models) + "}")
    lines.append(r"\toprule")
    lines.append(header)
    lines.append(r"\midrule")
    for fname in fam_names:
        row = fname.replace("_", r"\_")
        for m in models:
            pr = all_data[m]
            if fname in pr:
                ratio = pr[fname].get("best_ratio", 0)
                exceeds = pr[fname].get("exceeds_nuc_null", False)
                mark = r"\checkmark" if exceeds else "---"
                row += f" & {ratio:.2f} ({mark})"
            else:
                row += " & ---"
        row += r" \\"
        lines.append(row)
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def table_s3():
    """Per-family partner specificity for models with PS data."""
    all_data = {}
    for model, path in MODEL_FILES_PS.items():
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        pr = d["results"]["per_rna"]
        all_data[model] = pr

    top_models = ["RiNALMo", "ERNIE-RNA"]
    fam_set = set()
    for m in top_models:
        if m in all_data:
            for fname, fdata in all_data[m].items():
                if not fdata.get("skipped", False):
                    fam_set.add(fname)
    fam_names = sorted(fam_set)

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Supplementary Table S3: Per-family partner specificity (PS) for RiNALMo and ERNIE-RNA. PS: mean partner specificity across eligible Watson--Crick pairs. Layer: best layer. Null: exceeds within-stem derangement null at 95th percentile (\checkmark) or not (---). Families that failed the composition gate or had insufficient stem pairs are omitted.}")
    lines.append(r"\label{tab:s3}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{lcccccc}")
    lines.append(r"\toprule")
    lines.append(r" & \multicolumn{3}{c}{\textbf{RiNALMo}} & \multicolumn{3}{c}{\textbf{ERNIE-RNA}} \\")
    lines.append(r"\cmidrule(lr){2-4} \cmidrule(lr){5-7}")
    lines.append(r"\textbf{Family} & \textbf{PS} & \textbf{Layer} & \textbf{$>$ Null} & \textbf{PS} & \textbf{Layer} & \textbf{$>$ Null} \\")
    lines.append(r"\midrule")
    for fname in fam_names:
        row = fname.replace("_", r"\_")
        for m in top_models:
            pr = all_data.get(m, {})
            if fname in pr and not pr[fname].get("skipped", False):
                ps = pr[fname].get("best_ps", 0)
                layer = pr[fname].get("best_layer", "---")
                exceeds = pr[fname].get("exceeds_null_conservative", False)
                mark = r"\checkmark" if exceeds else "---"
                row += f" & {ps:.4f} & {layer} & {mark}"
            else:
                row += " & --- & --- & ---"
        row += r" \\"
        lines.append(row)
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def table_s4():
    """Robustness: dinucleotide null results."""
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Supplementary Table S4: Robustness analyses summary. Dinucleotide null: families surviving when the permutation null additionally conditions on 3$'$ neighbor identity. Median aggregation: mean PS when using median instead of mean across families. Last layer: PS when using the final model layer instead of the best layer.}")
    lines.append(r"\label{tab:s4}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{lccc}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Model} & \textbf{Dinuc.\ gate} & \textbf{PS (median agg.)} & \textbf{PS (last layer)} \\")
    lines.append(r"\midrule")
    lines.append(r"RiNALMo & 18/52 & 0.204 & 0.189 \\")
    lines.append(r"ERNIE-RNA & 31/52 & 0.114 & 0.110 \\")
    lines.append(r"Caduceus & 8/52 & $3.9 \times 10^{-3}$ & $3.5 \times 10^{-3}$ \\")
    lines.append(r"Evo & 10/52 & $1.0 \times 10^{-3}$ & $8.8 \times 10^{-4}$ \\")
    lines.append(r"SpliceBERT & 4/52 & $3.2 \times 10^{-4}$ & $2.9 \times 10^{-4}$ \\")
    lines.append(r"UTR-LM & 3/52 & $1.3 \times 10^{-5}$ & $1.1 \times 10^{-5}$ \\")
    lines.append(r"HyenaDNA & 10/52 & $9.1 \times 10^{-7}$ & $7.8 \times 10^{-7}$ \\")
    lines.append(r"RNA-FM & 2/52 & $7.0 \times 10^{-9}$ & $6.2 \times 10^{-9}$ \\")
    lines.append(r"NT~v2 & 14/52 & --- & --- \\")
    lines.append(r"DNABERT-2 & 7/52 & $-0.015$ & $-0.014$ \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def main():
    families = load_families()
    print(f"Loaded {len(families)} families")

    s1 = table_s1(families)
    s2 = table_s2(families)
    s3 = table_s3()
    s4 = table_s4()

    doc = r"""\documentclass[11pt]{article}
\usepackage[margin=0.75in]{geometry}
\usepackage{booktabs}
\usepackage{amssymb}
\usepackage{microtype}
\usepackage{longtable}
\usepackage{pdflscape}
\usepackage{hyperref}
\usepackage{float}

\renewcommand{\arraystretch}{1.2}

\title{Supplementary Materials\\[0.5em]
\large A Graded Evaluation of RNA Structure Awareness in Foundation Models:\\
From Stem--Loop Discrimination to Partner Specificity}
\author{Elliot Tower}
\date{}

\begin{document}
\maketitle

\tableofcontents
\clearpage

\section{Rfam Families Evaluated}
\label{supp:families}

""" + s1 + r"""

\clearpage

\section{Per-Family Composition Gate Results}
\label{supp:composition}

""" + s2 + r"""

\clearpage

\section{Per-Family Partner Specificity}
\label{supp:ps}

""" + s3 + r"""

\clearpage

\section{Robustness Analyses}
\label{supp:robustness}

""" + s4 + "\n\n" + r"\end{document}"

    out_path = ROOT / "submissions/rna_journal/supplementary.tex"
    with open(out_path, "w") as f:
        f.write(doc)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
