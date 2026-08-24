"""Which script produced each bound printed in Table 1 and Table 2?

Run:  uv run --no-project --with numpy --python 3.12 python \
          scripts/audit_table12_sources.py

Two scripts in this repository bootstrap the same per-family values and disagree
in the third decimal, because they draw from one seeded generator in different
orders:

  scripts/compute_bootstrap_cis.py  writes results/bootstrap_cis.json, the
                                    committed artifact, and reports exceedance
                                    counts against expected false positives.
  scripts/superseded/compute_cis.py writes the same filename with an
                                    incompatible flat schema, and assigns each
                                    family a p-value from a lookup table keyed on
                                    effect size before running Benjamini-Hochberg
                                    on those values.

This script recomputes both and reports, for every interval bound in Tables 1
and 2, which of them the printed literal came from. Bounds are compared at the
precision they are printed at, and separately against a bound widened outward --
a lower bound rounded down, an upper bound rounded up -- which is the rounding a
conservative author applies to an interval.

It also recomputes the Benjamini-Hochberg counts the Limitations section quotes,
since those come from compute_cis.py's fabricated p-values if they come from
anywhere.

compute_cis.py is imported for its pure functions only.
"""

import importlib.util
import json
import math
import re
from paper_versions import REPO, newest_paper

PAPER = newest_paper()
ARTIFACT = REPO / "results" / "bootstrap_cis.json"

_spec = importlib.util.spec_from_file_location(
    "compute_cis", REPO / "scripts" / "superseded" / "compute_cis.py")
_cis = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cis)

# printed label -> (key in bootstrap_cis.json, key in compute_cis.py output)
NAMES = {
    "ERNIE-RNA": ("ERNIE-RNA", "ERNIE-RNA"), "RiNALMo": ("RiNALMo", "RiNALMo"),
    "RNA-FM": ("RNA-FM", "RNA-FM"), "UTR-LM": ("UTR-LM", "UTR-LM"),
    "SpliceBERT": ("SpliceBERT", "SpliceBERT"), "NT~v2": ("NT v2", "NT v2"),
    "DNABERT-2": ("DNABERT-2", "DNABERT-2"), "HyenaDNA": ("HyenaDNA", "HyenaDNA"),
    "Caduceus": ("Caduceus", "Caduceus"), "Evo": ("Evo", "Evo"),
    "ERNIE-RNA untrained": ("ERNIE-RNA (untrained)", "ERNIE-RNA untrained"),
}

RUNG1 = re.compile(r"^(?P<label>[^&]+?)\s*&\s*[\d.]+\s*&\s*(?P<ci>\[[^\]]*\])\s*&")
RUNG2 = re.compile(r"^(?P<label>[^&]+?)\s*&\s*\d+\s*&\s*\d+\s*&\s*\d+\\%\s*&"
                   r"\s*(?P<ci>\[[^\]]*\])\s*\\\\$")


def clean(label: str) -> str:
    label = label.replace(r"($\sim$2M, RNA)", "").replace("$^\\dagger$", "")
    return re.sub(r"\s*\([^)]*\)\s*", " ", label).strip()


def section(text: str, label: str) -> str:
    start = text.index(rf"\label{{{label}}}")
    return text[start:text.index(r"\end{table}", start)]


def renderings(value: float, places: int, side: str) -> dict[str, str]:
    """How this value could legitimately be printed at this precision."""
    scale = 10 ** places
    widened = (math.floor(value * scale) if side == "lo"
               else math.ceil(value * scale)) / scale
    return {"": f"{value:.{places}f}", " widened": f"{widened:.{places}f}"}


def attribute(printed: str, side: str, candidates: dict[str, float]) -> list[str]:
    places = len(printed.split(".")[1]) if "." in printed else 0
    return [f"{name}{how}"
            for name, value in candidates.items()
            for how, rendered in renderings(value, places, side).items()
            if rendered == printed]


def main() -> int:
    text = PAPER.read_text()
    stored = json.loads(ARTIFACT.read_text())["models"]
    other = {model: _cis.compute_model_stats(model, _cis.load_per_family(path))
             for model, path in _cis.MODEL_FILES.items()}

    print("Table 1 and Table 2 interval bounds\n")
    print(f"{'row':<32} {'printed':>7} {'bootstrap_cis':>14} {'compute_cis':>12}"
          f"   produced by")
    print("-" * 100)
    only_one = {"compute_bootstrap_cis": [], "compute_cis": []}
    neither = []

    for table, pattern, keys in [
            ("tab:rung1", RUNG1, ("ci_95_lo", "ci_95_hi")),
            ("tab:rung2", RUNG2, ("dinuc_ci_lo", "dinuc_ci_hi"))]:
        for line in section(text, table).splitlines():
            match = pattern.match(line.strip())
            if not match or clean(match.group("label")) not in NAMES:
                continue
            artifact_key, other_key = NAMES[clean(match.group("label"))]
            block = (stored[artifact_key]["mutation_sensitivity"]
                     if table == "tab:rung1"
                     else stored[artifact_key]["frac_exceeds_dinuc_null"])
            row = other[other_key]
            lo, hi = (b.strip() for b in match.group("ci").strip("[]").split(","))
            for side, printed, art, oth in [
                    ("lo", lo, block["ci_lower"], row[keys[0]]),
                    ("hi", hi, block["ci_upper"], row[keys[1]])]:
                who = attribute(printed, side,
                                {"compute_bootstrap_cis": art, "compute_cis": oth})
                tag = f"{table.split(':')[1]} {artifact_key} {side}"
                print(f"{tag:<32} {printed:>7} {art:>14.4f} {oth:>12.4f}"
                      f"   {', '.join(who) if who else 'NEITHER'}")
                sources = {name.split(" widened")[0] for name in who}
                if not who:
                    neither.append((tag, printed, art, oth))
                elif len(sources) == 1:
                    only_one[sources.pop()].append(tag)

    print("\nBounds only one script can produce:")
    for name, rows in only_one.items():
        print(f"  {name}: {len(rows)}" + (f" -- {', '.join(rows)}" if rows else ""))
    if neither:
        print(f"\n{len(neither)} bounds match neither script:")
        for tag, printed, art, oth in neither:
            print(f"  {tag} = {printed} (bootstrap_cis {art:.4f}, "
                  f"compute_cis {oth:.4f})")

    print("\n\nBenjamini-Hochberg counts in the Limitations section\n")
    quoted = re.search(
        r"ERNIE-RNA: (\d+) \$\\to\$ (\d+); RiNALMo: (\d+) \$\\to\$ (\d+)", text)
    if quoted:
        print(f"  paper: ERNIE-RNA {quoted.group(1)} -> {quoted.group(2)}; "
              f"RiNALMo {quoted.group(3)} -> {quoted.group(4)}")
    else:
        print(f"  {PAPER.name} quotes no corrected counts.")
    for model in ("ERNIE-RNA", "RiNALMo"):
        row = other[model]
        print(f"  compute_cis.py would give: {model} {row['nuc_exceed_raw']} -> "
              f"{row['nuc_exceed_bh']}, from p-values it assigns by lookup")
    print("  compute_bootstrap_cis.py computes no p-values, so it can produce"
          "\n  no Benjamini-Hochberg count at all.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
