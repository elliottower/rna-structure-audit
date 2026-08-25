"""Resolve the paper's references against OpenAlex and write a BibTeX file.

The manuscript carried sixteen references as inline `\\bibitem` entries with
author lists truncated to "et al.", which is unusable as a source: three venues
need three bibliography styles, and a truncated author list cannot be expanded
by formatting. Rather than invent the missing names, each entry is matched
against OpenAlex by title and the full author list read off the record.

Every match is printed with the year and first author from the manuscript beside
the year and first author OpenAlex returned, so a wrong match is visible rather
than silent. Nothing is written for an entry that does not match.

    uv run --no-project --python 3.12 python scripts/resolve_bibliography.py
"""

import json
import re
import unicodedata
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "paper" / "references.bib"
API = "https://api.openalex.org/works"

# (bibtex key, title to match, expected year, expected first-author surname).
# Titles and years are read off the manuscript's own bibliography; they are the
# search key and the check, not the output.
ENTRIES = [
    ("rnafm", "Interpretable RNA foundation model from unannotated data for highly accurate RNA structure and function predictions", 2022, "Chen"),
    ("chen2024splicebert", "Self-supervised learning on millions of primary RNA sequences from 72 vertebrates improves sequence-based RNA splicing prediction", 2024, "Chen"),
    ("chu2024utrlm", "A 5' UTR language model for decoding untranslated regions of mRNA and function predictions", 2024, "Chu"),
    ("dalla2023nucleotide", "Nucleotide Transformer: building and evaluating robust foundation models for human genomics", 2024, "Dalla-Torre"),
    ("demezer2011mutant", "Mutant CAG repeats of huntingtin transcript fold into hairpins, form nuclear foci and are targets for RNA interference", 2011, "de Mezer"),
    ("kalvari2021rfam", "Rfam 14: expanded coverage of metagenomic, viral and microRNA families", 2021, "Kalvari"),
    ("krzyzosiak2012triplet", "Triplet repeat RNA structure and its role as pathogenic agent and therapeutic target", 2012, "Krzyzosiak"),
    ("lorenz2011viennarna", "ViennaRNA Package 2.0", 2011, "Lorenz"),
    ("nguyen2024evo", "Sequence modeling and design from molecular to genome scale with Evo", 2024, "Nguyen"),
    ("nguyen2024hyenadna", "HyenaDNA: long-range genomic sequence modeling at single nucleotide resolution", 2023, "Nguyen"),
    ("rinalmo", "RiNALMo: general-purpose RNA language models can generalize well on structure prediction tasks", 2025, "Peni"),
    ("schiff2024caduceus", "Caduceus: bi-directional equivariant long-range DNA sequence modeling", 2024, "Schiff"),
    ("szikszai2022generalization", "Deep learning models for RNA secondary structure prediction (probably) do not generalize across families", 2022, "Szikszai"),
    ("yin2024ernierna", "ERNIE-RNA: an RNA language model with structure-enhanced representations", 2025, "Yin"),
    ("zhou2024dnabert2", "DNABERT-2: efficient foundation model and benchmark for multi-species genome", 2024, "Zhou"),
]

# Not scholarly works, so OpenAlex has nothing to match. Written out here rather
# than added to the .bib by hand, which regeneration would silently drop.
SOFTWARE = '''@software{tower2026reproducible,
  author = {Tower, Elliot},
  title = {reproducible-science: prereg, results and citations},
  year = {2026},
  doi = {10.5281/zenodo.22100272},
  url = {https://github.com/elliottower/reproducible-science}
}

@misc{anthropic2025claude,
  author = {{Anthropic}},
  title = {Claude Code},
  year = {2025},
  howpublished = {\\url{https://claude.ai/code}}
}'''


def fetch(title):
    url = f"{API}?search={urllib.parse.quote(title)}&per_page=3"
    request = urllib.request.Request(url, headers={"User-Agent": "bib-resolver"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.loads(response.read()).get("results", [])


def authors(work):
    """Every author, surname first, in the order OpenAlex lists them."""
    names = []
    for authorship in work.get("authorships", []):
        name = (authorship.get("author") or {}).get("display_name")
        if not name:
            continue
        parts = name.split()
        surname_first = (f"{parts[-1]}, {' '.join(parts[:-1])}"
                         if len(parts) > 1 else name)
        names.append(NAME_OVERRIDES.get(surname_first, surname_first))
    return names


def comparable(surname):
    """A surname stripped of diacritics and particles, for matching only.

    OpenAlex writes Krzy\u017cosiak where the manuscript writes Krzyzosiak, and
    indexes "de Mezer" under Mezer. Neither is a different person, and neither
    should be written into the bibliography -- this normalization exists to
    compare two spellings, not to choose one.
    """
    folded = unicodedata.normalize("NFKD", surname)
    folded = "".join(c for c in folded if not unicodedata.combining(c))
    return folded.lower().replace("-", "").replace(" ", "")


# OpenAlex holds one author of the 2011 Nucleic Acids Research paper under a
# Cyrillic transliteration of a Polish name. The manuscript's own bibliography
# has it in Latin script, which is how the paper prints it. Overrides are listed
# by key and reason rather than applied silently, because an index returning the
# wrong script for one author is a reason to check the others, not to patch and
# move on.
NAME_OVERRIDES = {
    "\u041d\u0430\u043f\u0438\u0435\u0440\u0430\u043b\u0430, \u041c\u0430\u0440\u0435\u043a": "Napierala, Marek",
}

# pdflatex has no default mapping for a prime or a Polish crossed l, and a
# journal will not run xelatex. Names and titles are written as LaTeX escapes so
# the .bib compiles anywhere.
LATEX_ESCAPES = {
    "\u2032": "$'$", "\u2010": "-", "\u2013": "--", "\u2014": "---",
    "\u00e1": r"\'a", "\u00e9": r"\'e", "\u00f3": r"\'o", "\u00c9": r"\'E",
    "\u0107": r"\'c", "\u00f6": r'\"o', "\u0142": r"{\l}", "\u017c": r"\.z",
    "\u0160": r"\v{S}", "\u0161": r"\v{s}",
}


def to_latex(text):
    """Every non-ASCII character as a LaTeX escape, or an error naming it."""
    for character, escape in LATEX_ESCAPES.items():
        text = text.replace(character, escape)
    unmapped = sorted({c for c in text if ord(c) > 127})
    if unmapped:
        raise ValueError(
            f"no LaTeX escape for {[(c, hex(ord(c))) for c in unmapped]} in {text!r}; "
            "add it to LATEX_ESCAPES rather than dropping the character")
    return text


def brace_acronyms(title):
    """Protect capitalized tokens so BibTeX styles cannot lowercase them."""
    return re.sub(r"\b([A-Z][A-Za-z0-9-]*[A-Z0-9][A-Za-z0-9-]*)\b", r"{\1}", title)


def entry(key, work):
    venue = ((work.get("primary_location") or {}).get("source") or {}).get("display_name")
    biblio = work.get("biblio") or {}
    fields = [
        ("author", " and ".join(authors(work))),
        ("title", brace_acronyms(work.get("title") or "")),
        ("journal", venue or ""),
        ("year", str(work.get("publication_year") or "")),
        ("volume", biblio.get("volume") or ""),
        ("pages", "--".join(p for p in (biblio.get("first_page"),
                                        biblio.get("last_page")) if p)),
        ("doi", (work.get("doi") or "").replace("https://doi.org/", "")),
    ]
    body = ",\n".join(f"  {name} = {{{to_latex(value)}}}"
                      for name, value in fields if value)
    return f"@article{{{key},\n{body}\n}}"


def main():
    written, missed = [], []
    for key, title, year, first in ENTRIES:
        results = fetch(title)
        time.sleep(0.3)
        if not results:
            missed.append((key, "no result"))
            continue
        work = results[0]
        got_year = work.get("publication_year")
        got_first = (authors(work) or ["?"])[0].split(",")[0]
        # A one-year gap is the usual distance between an online and a print
        # date, and either surname may contain the other once particles and
        # diacritics are set aside.
        a, b = comparable(first), comparable(got_first)
        agrees = abs((got_year or 0) - year) <= 1 and (a in b or b in a)
        print(f"  {key:28s} {'ok ' if agrees else 'CHECK'} "
              f"{first} {year} -> {got_first} {got_year}, "
              f"{len(authors(work))} authors")
        if not agrees:
            missed.append((key, f"{got_first} {got_year}"))
            continue
        written.append(entry(key, work))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n\n".join([SOFTWARE] + written) + "\n")
    print(f"\n  wrote {OUT.relative_to(REPO)} with {len(written) + 1} entries")
    for key, why in missed:
        print(f"  UNRESOLVED {key}: {why}")


if __name__ == "__main__":
    main()
