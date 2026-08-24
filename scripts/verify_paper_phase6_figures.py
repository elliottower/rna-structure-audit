"""Check every Phase 6 figure in the manuscript against the stored results.

Run:  uv run --no-project --with scipy --python 3.12 python \
          scripts/verify_paper_phase6_figures.py [paper/paper_v11.tex]

Table 5 is parsed out of the .tex and each cell is recomputed from the
per-family values in results/, with the two quarantined families excluded as
PREREGISTRATION_PHASE6_V2.md requires. The prose figures scattered through
Rung 3, the abstract, the ladder figure, the hypothesis table and the
conclusion are checked the same way, by string.

A cell that no longer matches its stored value is an error, not a warning: the
manuscript and the results directory are the same claim written twice.
"""

import importlib.util
import re
from pathlib import Path

from paper_versions import REPO, newest_paper

DEFAULT_PAPER = newest_paper()

_spec = importlib.util.spec_from_file_location(
    "generate_table5_registered", REPO / "scripts" / "generate_table5_registered.py")
_gen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_gen)

# Table 5 row label -> results file, taken from the generator so the two
# scripts cannot drift apart.
SOURCES = {label.replace(r"\textbf{", "").replace("}", "").split(" (")[0]: rel
           for label, _domain, rel, _bold in _gen.ROWS}

ROW = re.compile(r"^(?P<label>[^&]+?)\s*&\s*(?:RNA|DNA)\s*&\s*(?P<ps>[^&]+?)\s*&"
                 r"\s*(?P<gate>[^&]+?)\s*&\s*(?P<null>[^&]+?)\s*&\s*(?P<h3>[^&]+?)"
                 r"\s*\\\\$")


def table_rows(text: str) -> dict:
    """The body rows of the table labelled tab:rung3, keyed by model."""
    start = text.index(r"\label{tab:rung3}")
    body = text[start:text.index(r"\end{table}", start)]
    rows = {}
    for line in body.splitlines():
        match = ROW.match(line.strip())
        if not match:
            continue
        key = (match.group("label").replace(r"\textbf{", "").replace("}", "")
               .split(" (")[0].strip())
        rows[key] = {name: match.group(name).replace(r"\textbf{", "")
                     .replace("}", "").strip()
                     for name in ("ps", "gate", "null", "h3")}
    return rows


def check_table(text: str) -> list:
    failures = []
    rows = table_rows(text)
    missing = set(SOURCES) - set(rows)
    if missing:
        failures.append(f"Table 5 is missing rows for {sorted(missing)}")
    for model, printed in rows.items():
        stored = _gen.summarize(SOURCES[model])
        expect = {
            "ps": ("---" if stored["mean_ps"] is None
                   else _gen.render_ps(stored["mean_ps"]).replace(r"\textbf{", "")
                   .replace("}", "").strip()),
            "gate": f"{stored['gate']}/{stored['eligible']}",
            "null": "---" if not stored["gate"] else f"{stored['exceed']}/{stored['gate']}",
            "h3": "---" if stored["h3"] is None else f"{stored['h3']:.3f}",
        }
        for column, want in expect.items():
            got = printed[column].replace("$", "").replace(" ", "")
            want_bare = want.replace("$", "").replace(" ", "")
            if got != want_bare:
                failures.append(
                    f"Table 5 {model} {column}: paper says {printed[column]!r}, "
                    f"results give {want!r}")
    return failures


def prose_checks() -> list:
    """(substring that must appear, what it asserts) for the figures in prose."""
    leaders = {name: _gen.summarize(SOURCES[name]) for name in ("RiNALMo", "ERNIE-RNA")}
    rinalmo, ernie = leaders["RiNALMo"], leaders["ERNIE-RNA"]
    untrained = _gen.summarize(SOURCES["ERNIE-RNA untrained"])
    checks = [
        (f"PS $= {rinalmo['mean_ps']:.3f}$", "Rung 3 prose, RiNALMo mean"),
        (f"PS $= {ernie['mean_ps']:.3f}$", "Rung 3 prose, ERNIE-RNA mean"),
        (f"{rinalmo['exceed']} of {rinalmo['gate']} gate-passing families",
         "Rung 3 prose, RiNALMo null count"),
        (f"{ernie['exceed']}/{ernie['gate']} families exceeding the",
         "Rung 3 prose, ERNIE-RNA null count"),
        (f"precision of {rinalmo['h3']:.3f}", "Rung 3 prose, RiNALMo H3"),
        (f"H3 precision of {ernie['h3']:.3f}", "Rung 3 prose, ERNIE-RNA H3"),
        (rf"PS\,=\,{rinalmo['mean_ps']:.3f}", "ladder figure, RiNALMo"),
        (rf"PS\,=\,{ernie['mean_ps']:.3f}", "ladder figure, ERNIE-RNA"),
        (rf"$\text{{PS}} = {rinalmo['mean_ps']:.3f}$", "conclusion, RiNALMo"),
        (rf"$\text{{PS}} = {ernie['mean_ps']:.3f}$", "conclusion, ERNIE-RNA"),
        (f"{rinalmo['h3'] * 100:.1f}\\% precision", "conclusion, RiNALMo H3"),
        (f"{ernie['h3'] * 100:.1f}\\% precision", "conclusion, ERNIE-RNA H3"),
        (f"RiNALMo: {rinalmo['h3']:.3f}; ERNIE-RNA: {ernie['h3']:.3f}",
         "hypothesis table, H3"),
        (f"RiNALMo: {rinalmo['exceed']}/{rinalmo['gate']}; "
         f"ERNIE-RNA: {ernie['exceed']}/{ernie['gate']}",
         "hypothesis table, H1"),
        (f"{untrained['gate']} of {untrained['eligible']}",
         "untrained control, gate count"),
    ]
    return checks


def check_no_stale_denominator(text: str) -> list:
    """Every Rung 3 gate is out of 32, 4 or 8. A /34 there is the old table."""
    start = text.index(r"\label{tab:rung3}")
    window = text[start - 3000:text.index(r"\section{Discussion}")]
    hits = sorted(set(re.findall(r"\b\d+/34\b", window)))
    return [f"a /34 gate survives in the Rung 3 section: {hits}"] if hits else []


def main(argv: list) -> int:
    paper = Path(argv[1]) if len(argv) > 1 else DEFAULT_PAPER
    text = paper.read_text()
    print(f"checking {paper.relative_to(REPO)} against results/\n")

    failures = check_table(text)
    print(f"Table 5: {len(table_rows(text))} rows, "
          f"{'OK' if not failures else str(len(failures)) + ' mismatches'}")

    for needle, what in prose_checks():
        if needle in text:
            print(f"  OK   {what}: {needle}")
        else:
            failures.append(f"{what}: expected {needle!r}, not found")

    failures += check_no_stale_denominator(text)

    if failures:
        print(f"\n{len(failures)} FAILURES")
        for line in failures:
            print(f"  {line}")
        return 1
    print("\nEvery Phase 6 figure in the manuscript matches the stored results.")
    return 0


if __name__ == "__main__":
    from sys import argv

    raise SystemExit(main(argv))
