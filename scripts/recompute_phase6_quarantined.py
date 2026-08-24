"""Recompute the Phase 6 table with the preregistered quarantine applied.

Run:  uv run --no-project --python 3.12 python scripts/recompute_phase6_quarantined.py

PREREGISTRATION_PHASE6_V2.md discloses that coupling-ratio values were observed
for tRNA_Phe_yeast and tRNA_Ala_human using RNA-FM before the registration was
written, and quarantines both families:

    Quarantined pilot families: tRNA_Phe_yeast and tRNA_Ala_human are excluded
    from all confirmatory hypothesis tests. Pilot values were observed before
    this preregistration.
    Eligible families: 34 of 52 ... After quarantine (2 families),
    N = 32 confirmatory families.

paper_v10.tex Table 5 reports gate denominators of 34 and never mentions the
quarantine, so the reported analysis appears to include the two families the
registration excludes. The stored results carry per-family values, so the
registered analysis can be recovered without re-running anything.

Prints the reported figure beside the quarantined figure for every model.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# Mean PS as printed in paper_v10.tex Table 5, for comparison.
REPORTED = {
    "rinalmo": 0.2150, "ernierna": 0.1204, "caduceus": 0.0034, "evo": 0.0011,
    "splicebert": 0.0002, "hyenadna": 0.0003, "rnafm": 0.0001, "utrlm": 0.00001,
    "nt": 0.001, "dnabert2": -0.016,
}


# One file per model, each the run that postdates the V2 freeze (2026-07-13).
# The rnafm_phase6_ps_20260713_*.json files under data/gpu_results/ are the
# pre-registration pilot runs and are not used here.
CANONICAL = {
    "rinalmo": "results/rinalmo_phase6_ps.json",
    "ernierna": "results/audit/phase6_ernierna_20260714_060520/ernierna_phase6_ps.json",
    "evo": "results/audit/phase6_evo_20260715_073218/evo_phase6_ps.json",
    "hyenadna": "results/audit/phase6_hyenadna_20260714_060542/hyenadna_phase6_ps.json",
    "nt": "results/audit/phase6_nt_20260714_082334/nt_phase6_ps.json",
    "rnafm": "results/audit/phase6_rnafm_20260714_065238/rnafm_phase6_ps.json",
    "splicebert": "results/audit/phase6_splicebert_20260714_060520/splicebert_phase6_ps.json",
    "utrlm": "results/audit/phase6_utrlm_20260714_060516/utrlm_phase6_ps.json",
    "caduceus": "data/gpu_results/expanded_rfam/caduceus_phase6_ps.json",
}


def per_family(payload: dict) -> dict:
    results = payload.get("results", payload)
    return results.get("per_rna", {})


def score(entries: dict, drop: set) -> tuple[float | None, int, int]:
    """Mean best PS over families that were actually scored, and the gate."""
    kept = {name: body for name, body in entries.items() if name not in drop}
    scored = [body for body in kept.values()
              if isinstance(body, dict) and not body.get("skipped")
              and isinstance(body.get("best_ps"), (int, float))]
    if not scored:
        return None, 0, len(kept)
    values = [body["best_ps"] for body in scored]
    return sum(values) / len(values), len(scored), len(kept)


def main() -> int:
    files = [REPO / p for p in CANONICAL.values()]
    missing = [f for f in files if not f.exists()]
    assert not missing, f"missing result files: {missing}"

    print(f"{'model':<12} {'paper':>10} {'all fams':>11} {'quarantined':>13} "
          f"{'gate':>9} {'shift':>9}")
    print("-" * 70)
    for path in files:
        model = [m for m, v in CANONICAL.items() if str(path).endswith(v)][0]
        entries = per_family(json.loads(path.read_text()))
        if not entries:
            print(f"{model:<12}  no per-family values stored")
            continue
        present = QUARANTINE & set(entries)
        full, gate_full, _ = score(entries, set())
        cut, gate_cut, n_cut = score(entries, QUARANTINE)
        paper = REPORTED.get(model)
        shift = ""
        if full and cut:
            shift = f"{(cut - full) / full * 100:+.1f}%"
        print(f"{model:<12} {paper if paper is not None else '--':>10} "
              f"{full if full is None else round(full, 6):>11} "
              f"{cut if cut is None else round(cut, 6):>13} "
              f"{f'{gate_cut}/{n_cut}':>9} {shift:>9}"
              + ("" if present else "   (neither quarantined family present)"))

    print("\nThe paper prints gate denominators of 34; the registration specifies "
          "N = 32 after quarantine.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
