"""Which registered hypotheses does the manuscript actually report?

Ten registered hypotheses were absent from the manuscript's summary table, three
of which fail against the paper's own reported values. That was found by reading
every preregistration by hand. This finds it in a second, and will keep finding
it as registrations and drafts change.

The check is deliberately dumb: it scans every preregistration for hypothesis
identifiers, scans the manuscript for the same identifiers, and reports the
difference. It does not try to decide whether a hypothesis passed -- it reports
the registered criterion so a human can.

Usage:
    uv run python scripts/check_hypothesis_coverage.py
    uv run python scripts/check_hypothesis_coverage.py --strict
"""

import argparse
import pathlib
import re
import sys

# Hypothesis identifiers: H1, H6b, H1_6, H_null, and LaTeX H1$_6$.
IDENT = re.compile(r"\bH(?:_null|[0-9]{1,2}[a-z]?(?:_6|\$_6\$)?)\b")

# A registration states a hypothesis on a line beginning with the identifier,
# usually bolded. Capture the sentence that follows for context.
INPUT = re.compile(r"\\(?:input|include)\{([^}]+)\}")

DECLARATION = re.compile(
    r"^\*{0,2}(H(?:_null|[0-9]{1,2}[a-z]?(?:_6)?))\*{0,2}\s*(?:\(([^)]*)\))?\s*[:.]?\s*(.*)",
    re.MULTILINE,
)


def normalise(ident):
    return ident.replace("$_6$", "_6").replace("\\", "").strip()


def scan_registrations(paths):
    """Return {hypothesis: [(file, criterion text), ...]}."""
    found = {}
    for path in paths:
        text = path.read_text(errors="ignore")
        for match in DECLARATION.finditer(text):
            ident = normalise(match.group(1))
            label = (match.group(2) or "").strip()
            rest = " ".join(match.group(3).split())[:150]
            criterion = f"{label}: {rest}" if label else rest
            found.setdefault(ident, []).append((path.name, criterion))
    return found


def read_with_inputs(path, seen=None):
    """Manuscript text with \\input and \\include expanded one file deep.

    The hypothesis table is generated into paper/generated/ and pulled in with
    \\input, so scanning the top-level .tex alone reports every hypothesis as
    absent.
    """
    seen = seen if seen is not None else set()
    path = path.resolve()
    if path in seen or not path.exists():
        return ""
    seen.add(path)
    text = path.read_text(errors="ignore")
    for name in INPUT.findall(text):
        target = path.parent / (name if name.endswith(".tex") else name + ".tex")
        text += "\n" + read_with_inputs(target, seen)
    return text


def scan_manuscript(path):
    return {normalise(m.group(0)) for m in IDENT.finditer(read_with_inputs(path))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper", default="paper/rna-structure-audit_v19.tex")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    root = pathlib.Path(".")
    registrations = sorted(
        p for p in list(root.glob("PREREG*.md"))
        + list(root.glob("preregistration/*.md"))
        + list(root.glob("docs/PREREG*.md"))
    )
    if not registrations:
        sys.exit("no preregistration files found")

    registered = scan_registrations(registrations)
    reported = scan_manuscript(pathlib.Path(args.paper))

    print(f"{len(registrations)} preregistrations, {len(registered)} hypotheses registered")
    for path in registrations:
        idents = sorted(k for k, v in registered.items() if any(f == path.name for f, _ in v))
        print(f"  {path.name}: {', '.join(idents) if idents else '(none parsed)'}")

    missing = sorted(set(registered) - reported)
    print(f"\nregistered but absent from {args.paper}: {len(missing)}")
    for ident in missing:
        source, criterion = registered[ident][0]
        print(f"  {ident:<8} [{source}]")
        if criterion:
            print(f"           {criterion}")

    orphans = sorted(reported - set(registered))
    if orphans:
        print(f"\nappear in the manuscript with no registration parsed: {', '.join(orphans)}")
        print("  (check by hand -- may be a parsing miss, or genuinely unregistered)")

    if args.strict and missing:
        sys.exit(1)


if __name__ == "__main__":
    main()
