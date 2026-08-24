"""Assemble the mobile reading packet for the newest manuscript in iCloud Drive.

Run:  uv run --no-project --with numpy --with scipy --with tqdm --python 3.12 \
          python scripts/build_icloud_packet.py

The packet is a folder per manuscript version holding the PDF, the LaTeX
source, EPUB and HTML conversions, the figures the manuscript includes, the
preregistration it reports against, the provenance record, and the output of
the checks run to build it. An existing packet is never overwritten: each
version keeps its own folder, so the packet for a superseded manuscript stays
readable beside the current one.

The checks run in this process rather than as subprocesses, so a packet cannot
be written for a manuscript whose figures do not reconcile: any check returning
non-zero raises before anything is copied.
"""

import contextlib
import io
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import audit_attention_rho
import check_freeze_order
import generate_table5_registered
import verify_artifact_regenerates
import verify_paper_phase6_figures
import verify_paper_rung12_figures
from paper_versions import REPO, VERSIONED, newest_paper

ICLOUD = Path.home() / "Library" / "Mobile Documents" / "com~apple~CloudDocs"
PAPER_DIR = REPO / "paper"
INCLUDED = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}")

# Each check is run with no manuscript argument, so each resolves the newest
# manuscript itself. The two that take argv get a bare argv[0].
CHECKS = [
    ("Tables 1 and 2 against results/bootstrap_cis.json",
     lambda: verify_paper_rung12_figures.main([""])),
    ("Table 5 and the Rung 3 prose against the per-family values",
     lambda: verify_paper_phase6_figures.main([""])),
    ("each attention correlation against the run that produced it",
     audit_attention_rho.main),
    ("the deposited artifact against its per-family inputs",
     verify_artifact_regenerates.main),
    ("preregistration freeze order",
     check_freeze_order.main),
]


def capture(fn) -> str:
    """Run a check, returning its stdout. Raises if the check fails."""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        status = fn()
    if status != 0:
        raise RuntimeError(f"check failed with status {status}:\n{buffer.getvalue()}")
    return buffer.getvalue()


def convert(source: Path, stem: str, dest: Path) -> None:
    """Write EPUB and HTML beside the PDF, with images embedded in the HTML."""
    for suffix, extra in ((".html", ["--embed-resources"]), (".epub", [])):
        subprocess.run(
            ["pandoc", source.name, "--standalone", "--mathml",
             f"--resource-path={PAPER_DIR}", *extra,
             "-o", str(dest / f"{stem}{suffix}")],
            cwd=source.parent, check=True, capture_output=True, text=True,
        )


def paper_as_png(paper: Path) -> Path:
    """A copy of the manuscript pointing at PNG figures, for pandoc.

    The manuscript includes PDF figures, which pandoc cannot embed in HTML or
    EPUB. Every figure exists in both formats.
    """
    handle = tempfile.NamedTemporaryFile(
        dir=PAPER_DIR, suffix=".tex", delete=False, mode="w"
    )
    with handle:
        handle.write(paper.read_text().replace(".pdf}", ".png}"))
    return Path(handle.name)


def main() -> int:
    paper = newest_paper()
    version = VERSIONED.match(paper.name).group(1)
    dest = ICLOUD / f"rna-structure-audit-v{version}"
    if dest.exists():
        raise FileExistsError(
            f"{dest} exists; a packet is never rewritten. Remove it by hand if "
            "the manuscript was rebuilt at the same version."
        )

    print(f"packet for {paper.name}\n")
    checks = [f"Checks run against {paper.name} when this packet was built.\n"]
    for description, fn in CHECKS:
        print(f"  running: {description}")
        checks.append(f"\n{'=' * 72}\n{description}\n{'=' * 72}\n")
        checks.append(capture(fn))
    table5 = capture(generate_table5_registered.main)

    dest.mkdir(parents=True)
    (dest / "figures").mkdir()
    shutil.copy2(paper, dest / paper.name)
    shutil.copy2(paper.with_suffix(".pdf"), dest / f"{paper.stem}.pdf")
    shutil.copy2(REPO / "PREREGISTRATION_PHASE6_V2.md", dest)
    shutil.copy2(REPO / "docs" / "provenance.md", dest)
    (dest / "verification.txt").write_text("".join(checks))
    (dest / "registered_phase6_figures.txt").write_text(table5)

    for reference in INCLUDED.findall(paper.read_text()):
        figure = (PAPER_DIR / reference).with_suffix(".png")
        shutil.copy2(figure, dest / "figures" / figure.name)

    png_paper = paper_as_png(paper)
    try:
        convert(png_paper, paper.stem, dest)
    finally:
        png_paper.unlink()
    convert(REPO / "docs" / "provenance.md", "provenance", dest)

    print(f"\nwrote {dest}")
    for item in sorted(dest.iterdir()):
        print(f"  {item.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
