"""Which numbers printed in the manuscript have a stored value behind them?

Run:  uv run --no-project --python 3.12 python scripts/audit_manuscript_numbers.py

The existing checks cover the Rung 1-2 tables, the Phase 6 figures and the
attention correlations. Numbers printed in prose have no systematic check, and
two defects were found by hand that this would have caught: a stem/loop GC pair
carried over from a 12-family pilot, and a length range matching no version of
the panel.

For every decimal literal in the manuscript, this script asks whether any value
stored under results/ rounds to it at the literal's own precision. That is the
same attribution technique as audit_table12_sources.py, applied to the whole
document rather than one table.

A literal with no match is not necessarily wrong. Thresholds, model sizes and
values computed at write-up time all legitimately have no stored source. The
output is a list of numbers to account for, not a list of errors.

No model is run.
"""

import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PAPER = REPO / "paper" / "paper_v12.tex"

# Literals that are never data: significance thresholds, chance baselines,
# corrected alphas, and the arithmetic of the design itself.
DESIGN_CONSTANTS = {
    "0.05", "0.01", "0.001", "0.0167", "0.333", "0.33", "1.96", "0.95", "0.5",
    "0.25", "0.75", "1.0", "0.0", "2.0", "3.0", "1.5",
}
LITERAL = re.compile(r"(?<![\d.])(\d+\.\d+)(?![\d])")
# A mantissa is not a value: "$7.3 \\times 10^{-8}$" prints 7.3e-8, and 7.3
# itself is never stored. Matching the mantissa alone would report every
# scientific-notation figure in the paper as unsourced.
MANTISSA = re.compile(r"(?<![\d.])\d+\.\d+(?=\s*(\\times|\$?\s*\\times))")


def stored_values() -> list[tuple[float, str]]:
    """Every number under results/, with the file it came from."""
    pool = []

    def walk(node, source):
        if isinstance(node, bool):
            return
        if isinstance(node, (int, float)):
            pool.append((float(node), source))
        elif isinstance(node, dict):
            for value in node.values():
                walk(value, source)
        elif isinstance(node, list):
            for value in node:
                walk(value, source)

    # Table 5's Caduceus row is sourced from data/gpu_results, not results/,
    # so a pool built from results/ alone reports sourced figures as unsourced.
    for path in sorted([*(REPO / "results").rglob("*.json"),
                        *(REPO / "data").rglob("*.json")]):
        try:
            walk(json.loads(path.read_text()), str(path.relative_to(REPO)))
        except json.JSONDecodeError:
            continue
    return pool


def body_lines() -> list[tuple[int, str]]:
    """Manuscript lines before the bibliography, comments stripped."""
    out = []
    for number, line in enumerate(PAPER.read_text().splitlines(), start=1):
        if "\\bibitem" in line or "\\begin{thebibliography}" in line:
            break
        stripped = re.sub(r"(?<!\\)%.*$", "", line)
        out.append((number, stripped))
    return out


def main() -> None:
    pool = stored_values()
    print(f"{len(pool)} numeric values stored under results/ and data/\n")

    by_precision: dict[int, dict[str, str]] = {}
    for value, source in pool:
        for places in range(1, 7):
            by_precision.setdefault(places, {}).setdefault(
                f"{value:.{places}f}", source)

    unmatched, matched = [], 0
    for number, line in body_lines():
        mantissas = set(MANTISSA.findall(line)) | {
            m.group(0) for m in MANTISSA.finditer(line)}
        for literal in LITERAL.findall(line):
            if literal in DESIGN_CONSTANTS or literal in mantissas:
                continue
            places = len(literal.split(".")[1])
            if places > 6:
                continue
            if by_precision.get(places, {}).get(literal):
                matched += 1
            else:
                unmatched.append((number, literal, line.strip()))

    total = matched + len(unmatched)
    print(f"{total} decimal literals in the body, excluding design constants")
    print(f"  {matched} have a stored value that rounds to them")
    print(f"  {len(unmatched)} do not\n")

    seen = set()
    for number, literal, context in unmatched:
        if literal in seen:
            continue
        seen.add(literal)
        context = re.sub(r"\s+", " ", context)[:96]
        print(f"  paper_v12.tex:{number:<5} {literal:<12} {context}")

    print(f"\n  {len(seen)} distinct unmatched literals. A literal with no stored "
          f"match is a number\n  to account for, not proof of an error.")


if __name__ == "__main__":
    main()
