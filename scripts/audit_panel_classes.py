"""Does the class breakdown printed in the manuscript describe the deposited panel?

Run:  uv run --no-project --python 3.12 python scripts/audit_panel_classes.py

paper_v12.tex:215-218 describes the panel as "tRNAs (7 families), rRNAs (4),
ribozymes (6), riboswitches (8), snRNAs (5), cis-regulatory elements (9), miRNA
precursors (5), CRISPR repeats (3), and other ncRNAs (5)."

Classification is a judgment for some families and not for others. A family named
tRNA_* is a tRNA and a family named *_riboswitch is a riboswitch; whether
Histone_3prime is cis-regulatory or other is arguable. This script assigns every
family by name, reports the counts against the printed ones, and separates the
classes whose membership is fixed by the filenames from the classes where the
assignment is a choice.

No model is run.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FAMILIES = REPO / "data/rfam_families"

PRINTED = {
    "tRNAs": 7, "rRNAs": 4, "ribozymes": 6, "riboswitches": 8, "snRNAs": 5,
    "cis-regulatory elements": 9, "miRNA precursors": 5, "CRISPR repeats": 3,
    "other ncRNAs": 5,
}

# Classes whose membership is determined by the family name, not by judgment.
UNAMBIGUOUS = {"tRNAs", "rRNAs", "ribozymes", "riboswitches", "snRNAs",
               "miRNA precursors", "CRISPR repeats"}

CIS_REGULATORY = {
    "Corona_5UTR", "Corona_s2m", "CrPV_IRES", "Hepatitis_C_IRES_III",
    "IRE_stem_loop", "IRES_HCV_domainII", "SECIS_element", "T_box_leader",
    "Histone_3prime",
}


def classify(name: str) -> str:
    if name.startswith("tRNA_"):
        return "tRNAs"
    if "rRNA" in name:
        return "rRNAs"
    if name.endswith("_ribozyme"):
        return "ribozymes"
    if name.endswith("_riboswitch"):
        return "riboswitches"
    if "snRNA" in name:
        return "snRNAs"
    if name.startswith("mir_"):
        return "miRNA precursors"
    if name.startswith("CRISPR"):
        return "CRISPR repeats"
    if name in CIS_REGULATORY:
        return "cis-regulatory elements"
    return "other ncRNAs"


def main() -> None:
    names = sorted(json.loads(p.read_text())["name"]
                   for p in FAMILIES.glob("*.json"))
    assigned: dict[str, list[str]] = {key: [] for key in PRINTED}
    for name in names:
        assigned[classify(name)].append(name)

    print(f"{len(names)} families in the deposited panel\n")
    print(f"  {'class':<26}{'printed':>8}{'deposited':>11}   basis")
    for key, printed in PRINTED.items():
        found = len(assigned[key])
        basis = "filename" if key in UNAMBIGUOUS else "assignment is a choice"
        mark = " " if found == printed else "*"
        print(f" {mark}{key:<26}{printed:>8}{found:>11}   {basis}")

    print(f"\n  printed total {sum(PRINTED.values())}, "
          f"deposited total {sum(len(v) for v in assigned.values())}")

    wrong = [(k, PRINTED[k], len(assigned[k])) for k in UNAMBIGUOUS
             if len(assigned[k]) != PRINTED[k]]
    print(f"\n  classes whose printed count is wrong and cannot be defended as a "
          f"classification\n  choice, because the filenames fix the membership: "
          f"{len(wrong)}\n")
    for key, printed, found in sorted(wrong):
        print(f"      {key:<26} printed {printed}, deposited {found}")
        print(f"          {', '.join(assigned[key])}")


if __name__ == "__main__":
    main()
