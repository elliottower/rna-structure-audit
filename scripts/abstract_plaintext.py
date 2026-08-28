"""Render a paper's abstract as plain text, with generated macros resolved.

Submission portals (OpenReview, journal forms) take a plain-text abstract. The
abstract in the .tex is written against macros in paper/generated/, so the text
pasted into a portal has to be built from the same source the PDF is, or the two
drift. This prints it; pipe to pbcopy.

Usage:
    uv run python scripts/abstract_plaintext.py paper/xai4science/main.tex
    uv run python scripts/abstract_plaintext.py paper/xai4science/main.tex | pbcopy
"""

import argparse
import pathlib
import re

NEWCOMMAND = re.compile(r"\\newcommand\{\\(\w+)\}\{(.*?)\}\s*$", re.MULTILINE)
ABSTRACT = re.compile(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", re.DOTALL)

# LaTeX spellings that have a plain-text equivalent, applied after macros resolve.
LITERAL = [
    (r"\\%", "%"), (r"\\&", "&"), (r"\\_", "_"),
    (r"---", "\u2014"), (r"--", "\u2013"),
    (r"``", "\u201c"), (r"''", "\u201d"),
    (r"\\nicefrac\{(\w+)\}\{(\w+)\}", r"\1/\2"),
    (r"\\(?:text(?:it|bf|rm|sc)|emph)\{([^{}]*)\}", r"\1"),
    (r"~", " "), (r"\$", ""), (r"\\,", " "), (r"\{\}", ""),
]


def macros(directory):
    out = {}
    for path in sorted(pathlib.Path(directory).glob("*.tex")):
        for name, value in NEWCOMMAND.findall(path.read_text()):
            out[name] = value
    return out


def render(tex_path, generated_dir):
    tex = pathlib.Path(tex_path).read_text()
    match = ABSTRACT.search(tex)
    if not match:
        raise SystemExit(f"no abstract environment in {tex_path}")

    body = match.group(1)
    table = macros(generated_dir)
    unresolved = set()

    def resolve(m):
        name = m.group(1)
        if name in table:
            return table[name]
        unresolved.add(name)
        return m.group(0)

    for _ in range(5):  # macros may expand into macros
        expanded = re.sub(r"\\(\w+)\{\}", resolve, body)
        if expanded == body:
            break
        body = expanded

    for pattern, replacement in LITERAL:
        body = re.sub(pattern, replacement, body)

    body = re.sub(r"%.*$", "", body, flags=re.MULTILINE)
    return " ".join(body.split()), sorted(unresolved)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tex")
    parser.add_argument("--generated", default="paper/generated")
    args = parser.parse_args()

    text, unresolved = render(args.tex, args.generated)
    print(text)
    if unresolved:
        raise SystemExit(f"\nunresolved macros: {', '.join(unresolved)}")


if __name__ == "__main__":
    main()
