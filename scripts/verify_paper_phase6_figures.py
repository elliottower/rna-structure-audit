"""Check every Phase 6 figure in the manuscript against the stored results.

Run:  uv run --no-project --with scipy --python 3.12 python \
          scripts/verify_paper_phase6_figures.py [paper/paper_vN.tex]

With no argument it checks the highest-numbered manuscript in paper/.

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

from paper_versions import REPO, expanded_text, newest_paper

DEFAULT_PAPER = newest_paper()

_spec = importlib.util.spec_from_file_location(
    "generate_table5_registered", REPO / "scripts" / "generate_table5_registered.py")
_gen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_gen)

# Table 5 row label -> results file, taken from the generator so the two
# scripts cannot drift apart.
# Built from the panel rather than from a hand-written row list. The list this
# read before (`generate_table5_registered.ROWS`) named one randomly initialized
# arm where the table prints seven, and carried the parameter counts D22
# corrected, so the check raised KeyError on a row it had no source for instead
# of reporting one.
PANEL = REPO / "results" / "repaired_panel_v3"

_KEYS = ["rinalmo", "ernierna", "caduceus", "evo", "splicebert", "hyenadna",
         "rnafm", "utrlm", "nt", "dnabert2"]
_LABELS = {"rinalmo": "RiNALMo", "ernierna": "ERNIE-RNA", "caduceus": "Caduceus",
           "evo": "Evo", "splicebert": "SpliceBERT", "hyenadna": "HyenaDNA",
           "rnafm": "RNA-FM", "utrlm": "UTR-LM", "nt": "NT~v2",
           "dnabert2": "DNABERT-2"}


def _source(key: str) -> str:
    return str((PANEL / key / f"{key}_phase6_ps.json").relative_to(REPO))


SOURCES = {_LABELS[key]: _source(key) for key in _KEYS}
SOURCES.update({f"{_LABELS[key]} untrained": _source(f"{key}_untrained")
                for key in _KEYS
                if (PANEL / f"{key}_untrained").is_dir()})

# Seven columns since the conservative null was added beside the primary one:
# model, domain, mean PS, gate, > null, > conservative, H3 precision. The
# six-column form matched nothing and every row read as missing.
ROW = re.compile(r"^(?P<label>[^&]+?)\s*&\s*(?:RNA|DNA)\s*&\s*(?P<ps>[^&]+?)\s*&"
                 r"\s*(?P<gate>[^&]+?)\s*&\s*(?P<null>[^&]+?)\s*&"
                 r"\s*(?P<cons>[^&]+?)\s*&\s*(?P<h3>[^&]+?)"
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
    unknown = set(rows) - set(SOURCES)
    if unknown:
        failures.append(
            f"Table 5 prints rows with no source in the panel: {sorted(unknown)}")
    for model, printed in rows.items():
        if model not in SOURCES:
            continue
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
    """(value that must appear, what it asserts) for the Rung 3 figures.

    These named the sentence as well as the number -- `PS $= 0.206$`, `85.9\\%
    precision`, a `ladder figure` -- and a rewrite that kept every number right
    failed all fourteen of them, because the manuscript no longer has a ladder
    figure and writes `\\psRiNALMo{}` where it once typed `0.206`. A check that
    fails on a rewrite it should be indifferent to stops being read.

    What survives is the weaker claim worth making: the number is somewhere in
    the manuscript. Its phrasing is not this script's business, and the specific
    drift the old form guarded -- prose disagreeing with the table -- is now
    caught by `literal_report`, which is general over every macro rather than
    over fourteen sentences.
    """
    leaders = {name: _gen.summarize(SOURCES[name]) for name in ("RiNALMo", "ERNIE-RNA")}
    rinalmo, ernie = leaders["RiNALMo"], leaders["ERNIE-RNA"]
    untrained = _gen.summarize(SOURCES["ERNIE-RNA untrained"])
    return [
        (f"{rinalmo['mean_ps']:.4f}", "RiNALMo mean PS"),
        (f"{ernie['mean_ps']:.4f}", "ERNIE-RNA mean PS"),
        (f"{rinalmo['exceed']}/{rinalmo['gate']}", "RiNALMo exceeding the null"),
        (f"{ernie['exceed']}/{ernie['gate']}", "ERNIE-RNA exceeding the null"),
        (f"{rinalmo['h3']:.3f}", "RiNALMo H3 precision"),
        (f"{ernie['h3']:.3f}", "ERNIE-RNA H3 precision"),
        (f"{untrained['gate']}/{untrained['eligible']}",
         "untrained control gate count"),
    ]


def literal_report(paper: Path) -> list:
    """Generated values typed into the manuscript instead of used as macros.

    A macro cannot disagree with the table it is generated from; a literal can,
    and does. The parameter ratio `75` was typed, written first as `80` from a
    nominal parameter count, and then left behind when the counts underneath it
    were corrected -- no table check would ever have looked at it.

    Three or more decimals, and the value unique across macros. At two decimals
    the report was eight hits and eight false positives: 0.05 was a significance
    level rather than the attention bound, 0.241 was ERNIE-RNA's chance rate
    rather than NT v2's attention, 1.00 was a Spearman correlation. A section
    that is all false positives trains the reader to skip it.

    Hits are candidates rather than failures. This cannot know which quantity a
    literal refers to, only that some macro shares its value.
    """
    source = paper.read_text()
    macros = {}
    for generated in sorted((PAPER_DIR := paper.parent).glob("generated/*.tex")):
        macros.update(_macro_defs(generated.read_text()))
    hits = []
    for name, value in sorted(macros.items()):
        bare = value.strip().rstrip("\\%")
        if not re.fullmatch(r"-?\d+\.\d{3,}", bare):
            continue
        if sum(1 for other in macros.values()
               if other.strip().rstrip("\\%") == bare) > 1:
            continue
        for number, line in enumerate(source.splitlines(), start=1):
            if line.lstrip().startswith("%"):
                continue
            if re.search(rf"(?<![\d.]){re.escape(bare)}(?![\d.])", line):
                hits.append(f"line {number}: {bare} is \\{name}; "
                            f"use the macro so it cannot drift")
    return hits


def _macro_defs(text: str) -> dict:
    return {m.group(1): m.group(2)
            for m in re.finditer(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}\s*$",
                                 text, re.M)}


def check_no_stale_denominator(text: str) -> list:
    """A gate denominator the current panel does not produce is an old table.

    This named `/34` as the stale one, from a panel where no model gated 34
    families. ERNIE-RNA and NT v2 both gate 34 of 36 now, so the check failed on
    two correct figures. The denominators come from the panel instead.
    """
    legitimate = {str(_gen.summarize(path)["gate"]) for path in SOURCES.values()}
    legitimate |= {str(_gen.summarize(path)["eligible"]) for path in SOURCES.values()}
    # The table body only. Widened to the surrounding prose it caught 31/47 from
    # Rung 1 and 0/5 from the HTT fragments, neither of which is a gate.
    start = text.index(r"\label{tab:rung3}")
    window = text[start:text.index(r"\end{tabular}", start)]
    hits = sorted({m for m in re.findall(r"\b\d+/(\d+)\b", window)
                   if m not in legitimate})
    if not hits:
        return []
    return [f"gate denominators no model in the panel produces: "
            f"{hits}; the panel gives {sorted(legitimate, key=int)}"]


def main(argv: list) -> int:
    paper = Path(argv[1]).resolve() if len(argv) > 1 else DEFAULT_PAPER
    # The expansion, because the tables live in paper/generated/ now.
    text = expanded_text(paper)
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

    typed = literal_report(paper)
    if typed:
        print(f"\n{len(typed)} literal(s) sharing a generated value, to account "
              f"for rather than to fix:")
        for line in typed:
            print(f"  {line}")

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
