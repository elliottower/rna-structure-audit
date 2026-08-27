"""The per-family tables, which are too long for the manuscript body.

The manuscript reports panel-level figures; a reader who wants to know which
family carried them needs the breakdown. The previous supplement was built on
the 52-family panel and lists families withdrawn on annotation review, so it
cannot ship beside a manuscript reporting 47.

    uv run --no-project --with "numpy<2" --python 3.12 \
        python scripts/generate_supplement.py
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results" / "repaired_panel_v3"
OUT = REPO / "paper" / "supplement_v17.tex"

QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
GATE_MODELS = [("ernierna", "ERNIE-RNA"), ("rinalmo", "RiNALMo"),
               ("rnafm", "RNA-FM"), ("nt", "NT~v2"), ("caduceus", "Caduceus"),
               ("evo", "Evo")]
PS_MODELS = [("rinalmo", "RiNALMo"), ("ernierna", "ERNIE-RNA")]


class MissingPanel(FileNotFoundError):
    """A model directory the supplement needs is absent."""


def load(key: str, stage: str) -> dict:
    path = PANEL / key / f"{key}_{stage}.json"
    if not path.exists():
        raise MissingPanel(f"missing {path.relative_to(REPO)}")
    return json.loads(path.read_text())


def escape(name: str) -> str:
    return name.replace("_", r"\_")


def families() -> list[str]:
    per = load("ernierna", "phases_1_to_5")["mutation_trained"]["per_rna"]
    return sorted(per)


def gate_table() -> str:
    runs = {key: load(key, "phases_1_to_5")["mutation_trained"]["per_rna"]
            for key, _ in GATE_MODELS}
    header = " & ".join(rf"\textbf{{{label}}}" for _, label in GATE_MODELS)
    lines = [
        r"\begin{longtable}{@{}l" + "r" * len(GATE_MODELS) + r"@{}}",
        r"\caption{Supplementary Table S2: per-family stem--loop sensitivity "
        r"ratio at each model's best layer. A dagger marks a family whose ratio "
        r"exceeds that model's nucleotide-stratified null at the 95th "
        r"percentile. Six of the ten models are shown; the remaining four never "
        r"exceed the null in more than three families.}\\",
        r"\label{tab:s2}\\",
        r"\toprule",
        rf"\textbf{{Family}} & {header} \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        rf"\textbf{{Family}} & {header} \\",
        r"\midrule",
        r"\endhead",
    ]
    for name in families():
        cells = []
        for key, _ in GATE_MODELS:
            body = runs[key].get(name, {})
            if not isinstance(body, dict) or "best_ratio" not in body:
                cells.append("---")
                continue
            mark = r"$^\dagger$" if body.get("exceeds_nuc_null") else ""
            cells.append(f"{body['best_ratio']:.3f}{mark}")
        lines.append(f"{escape(name)} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{longtable}", ""]
    return "\n".join(lines)


def ps_table() -> str:
    runs = {key: load(key, "phase6_ps")["results"]["per_rna"] for key, _ in PS_MODELS}
    lines = [
        r"\begin{longtable}{@{}lrrrrrr@{}}",
        r"\caption{Supplementary Table S3: per-family partner specificity for "
        r"the two models that clear Rung 3. PS is the mean over eligible pairs; "
        r"layer is where mean PS peaks in that family; precision is the fraction "
        r"of pairs whose partner moves more than either neighbour, against that "
        r"family's own within-stem derangement chance rate. The two pilot "
        r"families are quarantined from every confirmatory statistic and are "
        r"omitted.}\\",
        r"\label{tab:s3}\\",
        r"\toprule",
        r" & \multicolumn{3}{c}{\textbf{RiNALMo}} & "
        r"\multicolumn{3}{c}{\textbf{ERNIE-RNA}} \\",
        r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
        r"\textbf{Family} & \textbf{PS} & \textbf{Layer} & \textbf{Prec.} & "
        r"\textbf{PS} & \textbf{Layer} & \textbf{Prec.} \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"\textbf{Family} & \textbf{PS} & \textbf{Layer} & \textbf{Prec.} & "
        r"\textbf{PS} & \textbf{Layer} & \textbf{Prec.} \\",
        r"\midrule",
        r"\endhead",
    ]
    for name in families():
        if name in QUARANTINE:
            continue
        cells = []
        present = False
        for key, _ in PS_MODELS:
            body = runs[key].get(name, {})
            if (not isinstance(body, dict) or "skipped" in body
                    or not isinstance(body.get("best_ps"), (int, float))):
                cells += ["---", "---", "---"]
                continue
            present = True
            precision = body.get("h3_precision", {}).get("fraction")
            cells += [f"{body['best_ps']:+.4f}", str(body.get("best_layer", "---")),
                      "---" if precision is None else f"{precision:.3f}"]
        if present:
            lines.append(f"{escape(name)} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{longtable}", ""]
    return "\n".join(lines)


def main() -> int:
    document = r"""\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{microtype}

\title{Supplementary Material\\Only Two of Ten RNA and DNA Foundation Models
Retain Base-Pairing Signal Under Composition-Matched Nulls}
\author{Elliot Tower}
\date{}

\begin{document}
\maketitle

The manuscript reports panel-level figures. These are the per-family values they
aggregate, over the """ + str(len(families())) + r""" families of the repaired
panel. The family list, the attention--contact channel, the transversion
control, the registered hypotheses and the run provenance are in the
manuscript's own appendix.

\section*{Per-family composition gate}

""" + gate_table() + r"""

\section*{Per-family partner specificity}

""" + ps_table() + r"""

\end{document}
"""
    OUT.write_text(document)
    print(f"  wrote {OUT.relative_to(REPO)}")
    print(f"  {len(families())} families, {len(GATE_MODELS)} models in S2, "
          f"{len(PS_MODELS)} in S3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
