"""Check the Rung 1 and Rung 2 tables against results/bootstrap_cis.json.

Run:  uv run --no-project --python 3.12 python \
          scripts/verify_paper_rung12_figures.py [paper/paper_vNN.tex]

With no argument it checks the highest-numbered manuscript in paper/.

Tables 1 and 2 are parsed out of the .tex and every mean ratio, bootstrap
interval, exceedance count and retention rate is compared to the artifact that
`scripts/compute_bootstrap_cis.py` writes. Intervals are printed to two decimals
in Table 1 and Table 2 and to three in the artifact, so each bound is compared
at the precision it is printed at.

Attention rho and probing delta come from the per-model result files rather than
this artifact and are not checked here.
"""

import json
import re
from pathlib import Path

from paper_versions import REPO, newest_paper

DEFAULT_PAPER = newest_paper()
ARTIFACT = REPO / "results" / "bootstrap_cis.json"

# Table label -> artifact key. The tables spell the untrained control two ways.
NAMES = {
    "ERNIE-RNA": "ERNIE-RNA", "RiNALMo": "RiNALMo", "RNA-FM": "RNA-FM",
    "UTR-LM": "UTR-LM", "SpliceBERT": "SpliceBERT", "NT~v2": "NT v2",
    "DNABERT-2": "DNABERT-2", "HyenaDNA": "HyenaDNA", "Caduceus": "Caduceus",
    "Evo": "Evo", "ERNIE-RNA untrained": "ERNIE-RNA (untrained)",
}

RUNG1 = re.compile(
    r"^(?P<label>[^&]+?)\s*&\s*(?P<ratio>[\d.]+)\s*&\s*(?P<ci>\[[^\]]*\]|---)\s*&"
    r"\s*(?P<fam>[\d]+)/(?P<n>[\d]+)\s*&")
RUNG2 = re.compile(
    r"^(?P<label>[^&]+?)\s*&\s*(?P<nuc>\d+)\s*&\s*(?P<dinuc>\d+)\s*&"
    r"\s*(?P<ret>\d+)\\%\s*&\s*(?P<ci>\[[^\]]*\]|---)\s*\\\\$")


def clean(label: str) -> str:
    label = label.replace(r"($\sim$2M, RNA)", "").replace("$^\\dagger$", "")
    return re.sub(r"\s*\([^)]*\)\s*", " ", label).strip()


def section(text: str, label: str) -> str:
    start = text.index(rf"\label{{{label}}}")
    return text[start:text.index(r"\end{table}", start)]


def bounds(cell: str) -> tuple | None:
    if cell.strip() == "---":
        return None
    lo, hi = cell.strip("[]").split(",")
    return lo.strip(), hi.strip()


def matches(printed: str, value: float) -> bool:
    """Does the stored value print as this literal at the literal's precision?"""
    places = len(printed.split(".")[1]) if "." in printed else 0
    return f"{value:.{places}f}" == printed


def main(argv: list) -> int:
    paper = Path(argv[1]) if len(argv) > 1 else DEFAULT_PAPER
    text = paper.read_text()
    models = json.loads(ARTIFACT.read_text())["models"]
    print(f"checking {paper.relative_to(REPO)} against "
          f"{ARTIFACT.relative_to(REPO)}\n")

    failures, checked = [], 0

    for line in section(text, "tab:rung1").splitlines():
        match = RUNG1.match(line.strip())
        if not match or clean(match.group("label")) not in NAMES:
            continue
        model = NAMES[clean(match.group("label"))]
        stored = models[model]
        checked += 1
        sens = stored["mutation_sensitivity"]
        if not matches(match.group("ratio"), sens["point_estimate"]):
            failures.append(f"Table 1 {model} mean ratio: paper "
                            f"{match.group('ratio')}, artifact "
                            f"{sens['point_estimate']:.4f}")
        printed_ci = bounds(match.group("ci"))
        if printed_ci:
            for printed, value in zip(printed_ci, (sens["ci_lower"], sens["ci_upper"])):
                if not matches(printed, value):
                    failures.append(f"Table 1 {model} CI bound: paper {printed}, "
                                    f"artifact {value:.4f}")
        frac = stored["frac_exceeds_nuc_null"]
        fam, n = int(match.group("fam")), int(match.group("n"))
        if n != frac["n"] or fam != round(frac["point_estimate"] * frac["n"]):
            failures.append(f"Table 1 {model} exceedances: paper {fam}/{n}, "
                            f"artifact {round(frac['point_estimate'] * frac['n'])}"
                            f"/{frac['n']}")

    for line in section(text, "tab:rung2").splitlines():
        match = RUNG2.match(line.strip())
        if not match or clean(match.group("label")) not in NAMES:
            continue
        model = NAMES[clean(match.group("label"))]
        stored = models[model]
        checked += 1
        nuc = stored["frac_exceeds_nuc_null"]
        din = stored["frac_exceeds_dinuc_null"]
        if int(match.group("nuc")) != round(nuc["point_estimate"] * nuc["n"]):
            failures.append(f"Table 2 {model} fam > nuc: paper "
                            f"{match.group('nuc')}, artifact "
                            f"{round(nuc['point_estimate'] * nuc['n'])}")
        if int(match.group("dinuc")) != round(din["point_estimate"] * din["n"]):
            failures.append(f"Table 2 {model} survive dinuc: paper "
                            f"{match.group('dinuc')}, artifact "
                            f"{round(din['point_estimate'] * din['n'])}")
        if int(match.group("ret")) != round(din["point_estimate"] * 100):
            failures.append(f"Table 2 {model} retention: paper "
                            f"{match.group('ret')}%, artifact "
                            f"{din['point_estimate'] * 100:.0f}%")
        printed_ci = bounds(match.group("ci"))
        if printed_ci:
            for printed, value in zip(printed_ci, (din["ci_lower"], din["ci_upper"])):
                if not matches(printed, value):
                    failures.append(f"Table 2 {model} retention CI: paper "
                                    f"{printed}, artifact {value:.4f}")

    print(f"{checked} table rows checked")
    if failures:
        print(f"\n{len(failures)} FAILURES")
        for line in failures:
            print(f"  {line}")
        return 1
    print("Every mean ratio, interval, exceedance count and retention rate in\n"
          "Tables 1 and 2 matches the bootstrap artifact. Attention rho and\n"
          "probing delta are not checked here.")
    return 0


if __name__ == "__main__":
    from sys import argv

    raise SystemExit(main(argv))
