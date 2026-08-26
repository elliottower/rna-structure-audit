"""What the four re-annotated families did to the figures computed before them.

Run:  uv run --no-project --python 3.12 python scripts/audit_duplicate_families.py

Commit 5fc8914 replaced the annotations for mir_122_precursor,
mir_155_precursor, mir_21_precursor and mir_let7_precursor. Before it, three of
the four carried one byte-identical sequence and dot-bracket -- K02350.1/1-119,
an RF00001 (5S ribosomal RNA) seed member -- under three microRNA names, and
the fourth carried a 72-nt hairpin also labelled RF00001. Every stored result
in results/ was computed against those annotations and none has been re-run, so
the deposited per-family values still treat one 5S rRNA sequence as three
independent families.

This script reads the pre-fix annotations out of git, identifies the duplicate
groups by content, and reports what the Rung 3 aggregates become when the
duplicates are collapsed to one representative and when they are dropped. Both
are sensitivities. PREREGISTRATION_PHASE6_V2.md fixes N = 32, and changing N is
a registration matter, not an arithmetic one.

No model is run. Every value comes from the stored per-family records.
"""

import json
import subprocess
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PRE_FIX = "891d6af"
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

MODELS = {
    "RiNALMo": "results/rinalmo_phase6_ps.json",
    "ERNIE-RNA": ("results/audit/phase6_ernierna_20260714_060520/"
                  "ernierna_phase6_ps.json"),
}

# What paper_v12.tex prints for these two models over N = 32.
PRINTED = {
    "RiNALMo": dict(mean_ps=0.2098, gate="29/32", exceed="28/29", h3=0.874),
    "ERNIE-RNA": dict(mean_ps=0.1136, gate="30/32", exceed="28/30", h3=0.870),
}


def prefix_annotations() -> dict[str, dict]:
    names = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", PRE_FIX, "--", "data/rfam_families"],
        cwd=REPO, capture_output=True, text=True, check=True).stdout.split()
    records = {}
    for name in names:
        blob = subprocess.run(["git", "show", f"{PRE_FIX}:{name}"],
                              cwd=REPO, capture_output=True, text=True, check=True)
        record = json.loads(blob.stdout)
        records[record["name"]] = record
    return records


def duplicate_groups(records: dict[str, dict]) -> list[list[str]]:
    by_content = defaultdict(list)
    for name, record in records.items():
        by_content[(record["sequence"].upper(), record["dot_bracket"])].append(name)
    return [sorted(group) for group in by_content.values() if len(group) > 1]


def scored(rel: str) -> dict:
    body = json.loads((REPO / rel).read_text())["results"]["per_rna"]
    return {name: record for name, record in body.items()
            if isinstance(record, dict) and not record.get("skipped")
            and isinstance(record.get("best_ps"), (int, float))
            and name not in QUARANTINE}


def aggregate(entries: dict) -> dict:
    passing = [r for r in entries.values()
               if (r.get("positive_control") or {}).get("pass")]
    fractions = [r["h3_precision"]["fraction"] for r in passing
                 if isinstance(r.get("h3_precision"), dict)
                 and isinstance(r["h3_precision"].get("fraction"), (int, float))]
    values = [r["best_ps"] for r in entries.values()]
    return dict(
        n=len(entries),
        mean_ps=sum(values) / len(values) if values else float("nan"),
        gate=len(passing),
        exceed=sum(1 for r in passing if r.get("exceeds_null_primary")),
        h3=sum(fractions) / len(fractions) if fractions else float("nan"),
    )


def line(label: str, row: dict) -> str:
    return (f"  {label:<32} N={row['n']:2d}  mean PS {row['mean_ps']:.4f}  "
            f"gate {row['gate']:2d}/{row['n']:2d}  >null {row['exceed']:2d}/"
            f"{row['gate']:2d}  H3 {row['h3']:.3f}")


def main() -> int:
    records = prefix_annotations()
    groups = duplicate_groups(records)

    print(f"annotations as of {PRE_FIX} (the state every stored result was "
          f"computed against)\n")
    print(f"{len(records)} family files, "
          f"{len({(r['sequence'].upper(), r['dot_bracket']) for r in records.values()})} "
          f"distinct (sequence, dot-bracket) pairs\n")
    for group in groups:
        source = records[group[0]]["source"]
        print(f"  duplicated: {', '.join(group)}")
        print(f"    all carry {source}")
        print(f"    length {records[group[0]]['length']}, "
              f"rfam_id {records[group[0]]['rfam_id']}\n")

    dropped = {name for group in groups for name in group[1:]}
    all_dupes = {name for group in groups for name in group}

    for model, rel in MODELS.items():
        entries = scored(rel)
        printed = PRINTED[model]
        print(f"\n{model}  (paper_v12.tex prints mean PS {printed['mean_ps']}, "
              f"gate {printed['gate']}, >null {printed['exceed']}, "
              f"H3 {printed['h3']})")
        as_stored = aggregate(entries)
        print(line("as stored (registered N = 32)", as_stored))
        assert f"{as_stored['mean_ps']:.4f}" == f"{printed['mean_ps']:.4f}", (
            f"{model}: stored mean {as_stored['mean_ps']:.4f} does not match "
            f"the printed {printed['mean_ps']}")

        collapsed = aggregate({k: v for k, v in entries.items() if k not in dropped})
        print(line("duplicates collapsed to one", collapsed))
        removed = aggregate({k: v for k, v in entries.items() if k not in all_dupes})
        print(line("duplicate group removed", removed))
        for name in sorted(all_dupes & set(entries)):
            record = entries[name]
            print(f"    {name:<24} PS {record['best_ps']:.4f}  "
                  f"gate {'pass' if (record.get('positive_control') or {}).get('pass') else 'fail'}  "
                  f">null {record.get('exceeds_null_primary')}  "
                  f"H3 {record['h3_precision']['fraction']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
