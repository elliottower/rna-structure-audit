"""Resolve the newest manuscript version.

Versions are never edited in place: `<stem>_vN.tex` is copied from
`<stem>_vN-1.tex` and edited, and the older file stays as the record of what it
said. Verification scripts therefore have to find the current version rather than
name one, or they silently keep checking a superseded file.

Globbing for a version number is not enough on its own. This directory holds
several distinct documents that each carry versions -- `cross_architecture_v13`,
`unified_v2`, `main_v4` -- so a resolver that takes the highest N anywhere would
compare the manuscript against a different paper. It also has to survive the
manuscript being renamed, which is what broke it: the stem changed from `paper`
to `rna-structure-audit` at v13, the glob kept matching `paper_v12.tex`, and
every verification script went on reporting success against a manuscript five
versions behind the repaired panel.

The stem is therefore named once, here, and an absent stem raises.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO / "paper"

# The manuscript. Renaming it means changing this line and nothing else; leaving
# it stale means every check silently passes against the wrong file.
MANUSCRIPT_STEM = "rna-structure-audit"


class NoManuscript(FileNotFoundError):
    """No versioned file for the named stem, so there is nothing to check."""


def paper_versions(stem: str = MANUSCRIPT_STEM) -> list[tuple[int, Path]]:
    """Every `<stem>_vN.tex`, ascending by N."""
    pattern = re.compile(rf"^{re.escape(stem)}_v(\d+)\.tex$")
    found = [(int(match.group(1)), path)
             for path in PAPER_DIR.glob(f"{stem}_v*.tex")
             if (match := pattern.match(path.name))]
    return sorted(found)


def newest_paper(stem: str = MANUSCRIPT_STEM) -> Path:
    """The highest-numbered version of the manuscript."""
    versions = paper_versions(stem)
    if not versions:
        raise NoManuscript(
            f"no {stem}_vN.tex in {PAPER_DIR}. If the manuscript was renamed, "
            f"update MANUSCRIPT_STEM in {Path(__file__).name} -- a stale stem "
            f"makes every verification script pass against the wrong file.")
    return versions[-1][1]


def _macro_definitions(text: str) -> dict[str, str]:
    r"""Every `\newcommand{\name}{value}` in the text, by name."""
    found = {}
    for match in re.finditer(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}\s*$", text, re.M):
        found[match.group(1)] = match.group(2)
    return found


def expanded_text(paper: Path | None = None, resolve_macros: bool = True) -> str:
    r"""The manuscript as a reader sees it: `\input`s inlined, macros resolved.

    Two things moved out of the manuscript. Tables became
    `\input{generated/...}`, so a check reading the `.tex` for
    `\label{tab:rung1}` finds nothing. Numbers became macros written by the
    generators, so a check looking for `PS $= 0.206$` finds `\psRiNALMo{}`.
    Either way the check reports on a document nobody reads.

    Expanding both is what makes a verification script test the manuscript
    rather than its source form. Macros are substituted longest-name-first, so
    `\psRiNALMo` is not eaten by a prefix of it.
    """
    paper = paper or newest_paper()
    lines = []
    for line in paper.read_text().splitlines(keepends=True):
        stripped = line.strip()
        if stripped.startswith(r"\input{") and stripped.endswith("}"):
            target = PAPER_DIR / (stripped[len(r"\input{"):-1] + ".tex")
            if target.exists():
                lines.append(target.read_text())
                continue
        lines.append(line)
    text = "".join(lines)

    if not resolve_macros:
        return text
    macros = _macro_definitions(text)
    for name in sorted(macros, key=len, reverse=True):
        value = macros[name]
        text = text.replace(f"\\{name}{{}}", value).replace(f"\\{name} ", value + " ")
    return text
