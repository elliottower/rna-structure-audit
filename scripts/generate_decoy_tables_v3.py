"""Write the matched-decoy table and its macros from the registered audit's output.

Reads `results/audit/decoy_audit.json`, which `scripts/decoy_audit_v2.py` wrote, and
computes nothing: every cell is a stored field. A hypothesis row is refused unless the
audit recorded a decision for it.

Supersedes `generate_decoy_tables.py`, whose table ran 104pt past the text block: the
first column now wraps and the two control rows carry shorter labels. The macros are
unchanged.

    uv run --no-project python scripts/generate_decoy_tables_v3.py
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
AUDIT = REPO / "results" / "audit" / "decoy_audit.json"
GENERATED = REPO / "paper" / "generated_v24"

# (key, label, macro-safe name), in the order the table prints.
ROWS = [
    ("rinalmo", "RiNALMo", "RiNALMo"),
    ("ernierna", "ERNIE-RNA", "ErnieRNA"),
    ("rnafm", "RNA-FM", "RNAFM"),
    ("evo", "Evo", "Evo"),
    ("splicebert", "SpliceBERT", "SpliceBERT"),
    ("caduceus", "Caduceus", "Caduceus"),
    ("utrlm", "UTR-LM", "UTRLM"),
    ("hyenadna", "HyenaDNA", "HyenaDNA"),
    ("ernierna_untrained", "ERNIE-RNA, random init.", "ErnieRNAUntrained"),
    ("ernierna_untrained_noattnbias",
     "ERNIE-RNA, random init., pairing bias ablated", "ErnieRNAUntrainedAblated"),
]
REGISTERED = {"rinalmo": ("H1", "G1"), "ernierna": ("H2", "G2")}


class AuditMismatch(RuntimeError):
    """The audit file does not hold what the table is about to print."""


def rate(value) -> str:
    return "--" if value is None else f"{value:.3f}"


def p_value(value: float) -> str:
    return f"{value:.1e}".replace("e-0", r" \times 10^{-").replace("e-", r" \times 10^{-") + "}" \
        if value < 0.001 else f"{value:.3f}"


def cell(tag: dict) -> str:
    """Families above, below and at 0.5, and the mean rate."""
    if not tag["families"]:
        return "-- & --"
    return (f"{tag['above']}/{tag['below']}/{tag['at_half']} & {rate(tag['mean_rate'])}")


def main() -> int:
    audit = json.loads(AUDIT.read_text())
    models = {m["model"]: m for m in audit["models"]}
    missing = [key for key, *_ in ROWS if key not in models]
    if missing:
        raise AuditMismatch(f"the audit has no entry for {missing}")

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Matched-decoy test. For each mutated pair, the response at the annotated "
        r"partner is compared with the response at every other deep-interior paired position "
        r"of the family that carries the partner's nucleotide (all decoys), and separately "
        r"with only the decoys no farther from the mutated position than the partner "
        r"(close decoys). A family's rate is the mean over its pairs of the fraction of "
        r"decoys the partner exceeds, with ties counted as one half, so $0.5$ is the value "
        r"of a response shared equally by every such position. Families: the number above, "
        r"below and at exactly $0.5$, among confirmatory families passing the model's "
        r"positive control with at least three contributing pairs. The layer is chosen "
        r"with the family left out. Decision: the registered outcome, from an exact "
        r"one-sided sign test across families at $\alpha = 0.025$; rows with no decision "
        r"carry no registered hypothesis. Random init.: every parameter re-initialized. "
        r"NT~v2 and DNABERT-2 are outside the registered test because their tokens span "
        r"several nucleotides.}",
        r"\label{tab:decoy}",
        r"\small",
        r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{3.6cm}rrrrl@{}}",
        r"\toprule",
        r" & \multicolumn{2}{c}{\textbf{All decoys}} & \multicolumn{2}{c}{\textbf{Close decoys}} & \\",
        r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
        r"\textbf{Model} & \textbf{Families} & \textbf{Mean rate} & \textbf{Families} & "
        r"\textbf{Mean rate} & \textbf{Decision} \\",
        r"\midrule",
    ]
    macros = [r"% Written by scripts/generate_decoy_tables_v3.py. Do not edit."]
    for key, label, name in ROWS:
        model = models[key]
        every, close = model["family_all"], model["family_close"]
        decision = ""
        if key in REGISTERED:
            identity, guard = REGISTERED[key]
            decisions = model.get("decisions") or {}
            if identity not in decisions or guard not in decisions:
                raise AuditMismatch(f"{key}: the audit recorded no decision for "
                                    f"{identity} or {guard}")
            parts = []
            for hypothesis in (identity, guard):
                d = decisions[hypothesis]
                parts.append(f"{hypothesis} " + ("void" if d["void"] else
                                                 "holds" if d["met"] else "fails"))
            decision = ", ".join(parts)
        lines.append(f"{label} & {cell(every)} & {cell(close)} & {decision} \\\\")
        if key == "hyenadna":
            lines.append(r"\midrule")
        for tag_name, tag in (("All", every), ("Close", close)):
            macros += [
                rf"\newcommand{{\decoy{tag_name}Families{name}}}{{{tag['families']}}}",
                rf"\newcommand{{\decoy{tag_name}Above{name}}}{{{tag['above']}}}",
                rf"\newcommand{{\decoy{tag_name}Below{name}}}{{{tag['below']}}}",
                rf"\newcommand{{\decoy{tag_name}AtHalf{name}}}{{{tag['at_half']}}}",
                rf"\newcommand{{\decoy{tag_name}Rate{name}}}{{{rate(tag['mean_rate'])}}}",
                rf"\newcommand{{\decoy{tag_name}Pairs{name}}}{{{tag['pairs']}}}",
                rf"\newcommand{{\decoy{tag_name}P{name}}}{{{p_value(tag['p'])}}}",
            ]
        macros.append(rf"\newcommand{{\decoyPopulation{name}}}{{{model['population_families']}}}")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    macros.append("")

    GENERATED.mkdir(parents=True, exist_ok=True)
    (GENERATED / "decoy_table_v3.tex").write_text("\n".join(lines))
    (GENERATED / "decoy_macros_v3.tex").write_text("\n".join(macros))
    print(f"wrote decoy_table_v3.tex and decoy_macros_v3.tex to {GENERATED.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
