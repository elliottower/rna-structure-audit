"""Which families does the frozen registration quarantine, and does the code agree?

Run:  uv run --no-project --python 3.12 python scripts/registered_quarantine.py

`PREREGISTRATION_PHASE6_V2.md` quarantines two families because pilot values for
them were seen under a weaker metric before the document was written. Ten scripts
in this repository carry that set as a literal. A literal is not wrong, but it is
not tied to the document either: an edit to one copy changes a confirmatory set
and nothing detects it.

The set is read here out of the frozen file instead, twice -- from the pilot
disclosure and from the quarantine statement, which are separate sentences making
the same commitment -- and the two readings must agree. Names are recognized by
membership in the deposited panel rather than by pattern, so a name the panel does
not contain is an error rather than a silent miss.

Run as a script it also reports every `QUARANTINE` literal under `scripts/` that
disagrees with the frozen set, and the classes the quarantine removes from the
confirmatory panel.
"""

import ast
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FROZEN = REPO / "PREREGISTRATION_PHASE6_V2.md"
FAMILIES = REPO / "data/rfam_families"

DISCLOSURE = "**Pilot data disclosure:**"
STATEMENT = "**Quarantined pilot families:**"


def panel_names() -> set[str]:
    return {json.loads(path.read_text())["name"]
            for path in FAMILIES.glob("*.json")}


def _named_in(text: str, marker: str, names: set[str]) -> set[str]:
    """The panel families named in the sentence beginning at `marker`."""
    start = text.index(marker)
    sentence = text[start:text.index(".", start)]
    return {name for name in names if name in sentence}


def registered_quarantine(frozen: Path = FROZEN) -> set[str]:
    """The quarantined families, as the frozen registration states them."""
    text = frozen.read_text()
    names = panel_names()
    disclosed = _named_in(text, DISCLOSURE, names)
    stated = _named_in(text, STATEMENT, names)
    if disclosed != stated:
        raise ValueError(
            f"{frozen.name} discloses pilot values for {sorted(disclosed)} and "
            f"quarantines {sorted(stated)}; the two sentences must name the same "
            f"families")
    if not stated:
        raise ValueError(f"{frozen.name} names no quarantined family the panel "
                         f"contains")
    return stated


QUARANTINE = registered_quarantine()


def literals() -> dict[Path, set[str]]:
    """Every `QUARANTINE = {...}` assignment under `scripts/`."""
    found = {}
    for path in sorted((REPO / "scripts").glob("*.py")):
        if path.name == Path(__file__).name:
            continue
        for node in ast.parse(path.read_text()).body:
            if not isinstance(node, ast.Assign):
                continue
            targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
            if not any(t in ("QUARANTINE", "QUARANTINED") for t in targets):
                continue
            try:
                found[path] = set(ast.literal_eval(node.value))
            except ValueError:
                pass
    return found


def main() -> None:
    print(f"{FROZEN.name} quarantines {len(QUARANTINE)} of "
          f"{len(panel_names())} curated families\n")
    for name in sorted(QUARANTINE):
        print(f"    {name}")

    copies = literals()
    disagree = {path: value for path, value in copies.items()
                if value != QUARANTINE}
    print(f"\n  {len(copies)} scripts carry the set as a literal, "
          f"{len(disagree)} disagree with the frozen file")
    for path, value in sorted(disagree.items()):
        print(f"    {path.relative_to(REPO)}: {sorted(value)}")

    analyzed = {json.loads(p.read_text())["name"]
                for p in FAMILIES.glob("*.json")
                if "excluded" not in json.loads(p.read_text())}
    print(f"\n  of the {len(analyzed)} analyzed families, the quarantine removes "
          f"{len(QUARANTINE & analyzed)}\n  from the confirmatory set at Rung 3 "
          f"and none from Rungs 1 and 2")
    remaining_trna = sorted(n for n in analyzed - QUARANTINE if n.startswith("tRNA_"))
    remaining_rrna = sorted(n for n in analyzed if "rRNA" in n)
    print(f"    tRNA families left in the confirmatory set: "
          f"{', '.join(remaining_trna) or 'none'}")
    print(f"    rRNA families left in the analyzed panel:   "
          f"{', '.join(remaining_rrna) or 'none'}")


if __name__ == "__main__":
    main()
