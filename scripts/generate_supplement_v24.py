"""The per-family tables, which are too long for the manuscript body.

The manuscript reports panel-level figures; a reader who wants to know which
family carried them needs the breakdown.

Supersedes `generate_supplement.py`. This reads `results/panel_v4`, writes the two
tables as files the supplement inputs instead of a whole document, shows the three
models that clear Rung 3 in Table S3, and counts what its caption says about the
models Table S2 leaves out instead of stating it.

    uv run --no-project python scripts/generate_supplement_v24.py
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results" / "panel_v4"
GENERATED = REPO / "paper" / "generated_v24"

QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
GATE_MODELS = [("ernierna", "ERNIE-RNA"), ("rinalmo", "RiNALMo"),
               ("rnafm", "RNA-FM"), ("nt", "NT~v2"), ("caduceus", "Caduceus"),
               ("evo", "Evo")]
PS_MODELS = [("rinalmo", "RiNALMo"), ("ernierna", "ERNIE-RNA"), ("rnafm", "RNA-FM")]
OTHER_MODELS = ["utrlm", "splicebert", "dnabert2", "hyenadna"]


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
    # What the caption says about the models not shown, counted from their runs.
    left_out = max(
        sum(1 for body in load(key, "phases_1_to_5")["mutation_trained"]["per_rna"].values()
            if isinstance(body, dict) and body.get("exceeds_nuc_null"))
        for key in OTHER_MODELS)
    lines = [
        r"\begin{longtable}{@{}l" + "r" * len(GATE_MODELS) + r"@{}}",
        r"\caption{Per-family stem--loop sensitivity "
        r"ratio at each model's best layer. A dagger marks a family whose ratio "
        r"exceeds that model's nucleotide-stratified null at the 95th "
        r"percentile. Six of the ten models are shown; none of the remaining four "
        rf"exceeds the null in more than {left_out} families.}}\\",
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
    column_heads = (r"\textbf{Family} & " + " & ".join(
        [r"\textbf{PS} & \textbf{Layer} & \textbf{Prec.}"] * len(PS_MODELS)) + r" \\")
    lines = [
        r"\small",
        r"\begin{longtable}{@{}l" + "rrr" * len(PS_MODELS) + r"@{}}",
        r"\caption{Per-family partner specificity for "
        r"the three models that clear Rung 3. PS is the mean over eligible pairs; "
        r"layer is where mean PS peaks in that family; precision is the fraction "
        r"of pairs whose partner moves more than either neighbour, against that "
        r"family's own within-stem derangement chance rate. The two pilot "
        r"families are quarantined from every confirmatory statistic and are "
        r"omitted.}\\",
        r"\label{tab:s3}\\",
        r"\toprule",
        " & " + " & ".join(rf"\multicolumn{{3}}{{c}}{{\textbf{{{label}}}}}"
                           for _, label in PS_MODELS) + r" \\",
        "".join(rf"\cmidrule(lr){{{2 + 3 * i}-{4 + 3 * i}}}" for i in range(len(PS_MODELS))),
        column_heads,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        column_heads,
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
    lines += [r"\bottomrule", r"\end{longtable}", r"\normalsize", ""]
    return "\n".join(lines)


def main() -> int:
    GENERATED.mkdir(parents=True, exist_ok=True)
    (GENERATED / "s2_gate_table.tex").write_text(gate_table())
    (GENERATED / "s3_ps_table.tex").write_text(ps_table())
    print(f"  wrote s2_gate_table.tex and s3_ps_table.tex to {GENERATED.relative_to(REPO)}")
    print(f"  {len(families())} families, {len(GATE_MODELS)} models in S2, "
          f"{len(PS_MODELS)} in S3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
