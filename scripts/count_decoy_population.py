"""Families each model keeps for the decoy test once the positive control is applied.

Reads the structure-only family lists of `matched_decoy_counts.json` and the
stored positive-control outcomes of the repaired panel. No decoy rate exists in
those files; the only model outcome read is whether a family passed the
stem-versus-loop positive control, which was already reported.

    uv run --no-project python scripts/count_decoy_population.py
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PANEL = REPO / "results" / "repaired_panel_v3"
COUNTS = REPO / "results" / "audit" / "matched_decoy_counts.json"
OUT = REPO / "results" / "audit" / "decoy_population_counts.json"

QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
MODELS = ["rinalmo", "ernierna", "caduceus", "rnafm", "utrlm", "splicebert",
          "hyenadna", "evo", "ernierna_untrained"]


def main() -> int:
    summary = json.loads(COUNTS.read_text())["summary"]
    structural = {tag: set(summary[tag]["families"])
                  for tag in ("family_all", "family_close")}
    rows = {}
    for key in MODELS:
        per_rna = json.loads(
            (PANEL / key / f"{key}_phase6_ps.json").read_text())["results"]["per_rna"]
        confirmatory = {name for name, entry in per_rna.items()
                        if name not in QUARANTINED and not entry.get("skipped")}
        passing = {name for name in confirmatory
                   if per_rna[name].get("positive_control", {}).get("pass", False)}
        rows[key] = {"confirmatory": len(confirmatory), "positive_control": len(passing),
                     **{tag: len(families & passing)
                        for tag, families in structural.items()}}
    OUT.write_text(json.dumps({
        "panel": str(PANEL.relative_to(REPO)),
        "structural": {tag: len(families) for tag, families in structural.items()},
        "models": rows}, indent=2))
    for key, row in rows.items():
        print(f"  {key:22s} {row}")
    print(f"wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
