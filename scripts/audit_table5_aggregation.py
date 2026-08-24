"""Determine which aggregation produced each figure in paper_v10.tex Table 5.

Run:  uv run --no-project --python 3.12 python scripts/audit_table5_aggregation.py

PREREGISTRATION_PHASE6_V2.md fixes the confirmatory set before the runs:

    Eligible families: 34 of 52 families pass the 15-WC-pair threshold. After
    quarantine (2 families), N = 32 confirmatory families.

    [positive control] ... exclude that model-family pair from confirmatory
    analysis. This prerequisite is a screening gate applied per model per
    family, not a confirmatory hypothesis.

So the registered mean is over non-quarantined families that pass the gate, and
the registered gate denominator is 32. Table 5 prints denominators of 34.

Each results file already stores `families_quarantined` and a per-family
`quarantined` flag, so the pipeline knew about the quarantine. This script
recomputes every candidate aggregation from the stored per-family values and
reports which one reproduces each printed figure, so the correction is made
against a known cause rather than a guess. Nothing is re-run: `best_ps`,
`positive_control`, `exceeds_null_primary` and `h3_precision` are all stored.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# One file per model, each the run that postdates the V2 freeze (2026-07-13).
CANONICAL = {
    "rinalmo": "results/rinalmo_phase6_ps.json",
    "ernierna": "results/audit/phase6_ernierna_20260714_060520/ernierna_phase6_ps.json",
    "caduceus": "data/gpu_results/expanded_rfam/caduceus_phase6_ps.json",
    "evo": "results/audit/phase6_evo_20260715_073218/evo_phase6_ps.json",
    "splicebert": "results/audit/phase6_splicebert_20260714_060520/splicebert_phase6_ps.json",
    "hyenadna": "results/audit/phase6_hyenadna_20260714_060542/hyenadna_phase6_ps.json",
    "rnafm": "results/audit/phase6_rnafm_20260714_065238/rnafm_phase6_ps.json",
    "utrlm": "results/audit/phase6_utrlm_20260714_060516/utrlm_phase6_ps.json",
    "nt": "results/audit/phase6_nt_20260714_082334/nt_phase6_ps.json",
}

# Mean PS, gate, and H3 precision exactly as printed in Table 5 of paper_v10.tex.
# Mean PS is kept as the printed literal: the table rounds to as few as one
# significant figure, so a candidate is only comparable at that same precision.
PRINTED = {
    "rinalmo": ("0.2150", 31, 34, 0.882),
    "ernierna": ("0.1204", 31, 34, 0.874),
    "caduceus": ("0.0034", 16, 34, 0.416),
    "evo": ("0.0011", 10, 34, 0.278),
    "splicebert": ("0.0002", 22, 34, 0.276),
    "hyenadna": ("0.0003", 27, 34, 0.083),
    "rnafm": ("0.0001", 8, 34, 0.325),
    "utrlm": ("0.00001", 16, 34, 0.350),
    "nt": ("0.001", 1, 4, None),
}


def sig_figs(literal: str) -> int:
    """Significant figures in a printed decimal literal."""
    return len(literal.replace("0.", "").lstrip("0")) or 1


def scored(entries: dict) -> dict:
    """Families with a real PS value: not skipped for too few WC pairs."""
    return {name: body for name, body in entries.items()
            if isinstance(body, dict) and not body.get("skipped")
            and isinstance(body.get("best_ps"), (int, float))}


def gate_pass(body: dict) -> bool:
    control = body.get("positive_control")
    return bool(control and control.get("pass"))


def mean(values: list) -> float | None:
    return sum(values) / len(values) if values else None


def h3_mean_of_fractions(bodies: list) -> float | None:
    fractions = [b["h3_precision"]["fraction"] for b in bodies
                 if isinstance(b.get("h3_precision"), dict)
                 and isinstance(b["h3_precision"].get("fraction"), (int, float))]
    return mean(fractions)


def h3_pooled(bodies: list) -> float | None:
    hits = sum(b["h3_precision"].get("partner_max_count", 0) for b in bodies
               if isinstance(b.get("h3_precision"), dict))
    total = sum(b["h3_precision"].get("n", 0) for b in bodies
                if isinstance(b.get("h3_precision"), dict))
    return hits / total if total else None


def rounds_to(candidate: float | None, literal: str) -> bool:
    """Does the candidate print as this literal at the literal's precision?"""
    if candidate is None:
        return False
    figs = sig_figs(literal)
    return f"{candidate:.{figs}g}" == f"{float(literal):.{figs}g}"


def main() -> int:
    print(f"{'model':<11} {'printed':>9} {'all34':>9} {'nonq32':>9} "
          f"{'gate/all':>9} {'gate/nonq':>10} {'stored':>9}   which")
    print("-" * 88)
    verdicts = {}
    for model, rel in CANONICAL.items():
        path = REPO / rel
        assert path.exists(), f"missing {rel}"
        payload = json.loads(path.read_text())
        results = payload.get("results", payload)
        entries = scored(results.get("per_rna", {}))
        stored = results.get("mean_best_ps")

        every = list(entries.values())
        nonq = [b for n, b in entries.items() if n not in QUARANTINE]
        gate_all = [b for b in every if gate_pass(b)]
        gate_nonq = [b for b in nonq if gate_pass(b)]

        cand = {
            "all34": mean([b["best_ps"] for b in every]),
            "nonq32": mean([b["best_ps"] for b in nonq]),
            "gate/all": mean([b["best_ps"] for b in gate_all]),
            "gate/nonq": mean([b["best_ps"] for b in gate_nonq]),
        }
        printed_ps, printed_gate_n, printed_gate_d, printed_h3 = PRINTED[model]
        match = [k for k, v in cand.items() if rounds_to(v, printed_ps)] or ["none"]
        verdicts[model] = match

        def show(x):
            return "--" if x is None else f"{x:.6f}"
        print(f"{model:<11} {printed_ps:>9} {show(cand['all34']):>9} "
              f"{show(cand['nonq32']):>9} {show(cand['gate/all']):>9} "
              f"{show(cand['gate/nonq']):>10} {show(stored):>9}   {'+'.join(match)}")

    print("\nGate column: printed numerator/denominator vs recomputed")
    print(f"{'model':<11} {'printed':>9} {'pass/34':>9} {'pass/32':>9}")
    print("-" * 42)
    for model, rel in CANONICAL.items():
        payload = json.loads((REPO / rel).read_text())
        results = payload.get("results", payload)
        entries = scored(results.get("per_rna", {}))
        nonq = {n: b for n, b in entries.items() if n not in QUARANTINE}
        _, gn, gd, _ = PRINTED[model]
        print(f"{model:<11} {f'{gn}/{gd}':>9} "
              f"{f'{sum(gate_pass(b) for b in entries.values())}/{len(entries)}':>9} "
              f"{f'{sum(gate_pass(b) for b in nonq.values())}/{len(nonq)}':>9}")

    print("\nH3 precision: printed vs mean-of-fractions and pooled, over "
          "gate-passing families with and without the quarantine")
    print(f"{'model':<11} {'printed':>9} {'frac/34':>9} {'pool/34':>9} "
          f"{'frac/32':>9} {'pool/32':>9}   which")
    print("-" * 70)
    for model, rel in CANONICAL.items():
        payload = json.loads((REPO / rel).read_text())
        results = payload.get("results", payload)
        entries = scored(results.get("per_rna", {}))
        every = [b for b in entries.values() if gate_pass(b)]
        nonq = [b for n, b in entries.items()
                if n not in QUARANTINE and gate_pass(b)]
        _, _, _, printed_h3 = PRINTED[model]
        cand = {"frac/34": h3_mean_of_fractions(every), "pool/34": h3_pooled(every),
                "frac/32": h3_mean_of_fractions(nonq), "pool/32": h3_pooled(nonq)}
        match = ["none"] if printed_h3 is None else (
            [k for k, v in cand.items()
             if v is not None and abs(v - printed_h3) < 5e-4] or ["none"])
        cells = "".join(f"{'--' if v is None else f'{v:.3f}':>10}"
                        for v in cand.values())
        print(f"{model:<11} {str(printed_h3):>9}{cells}   {'+'.join(match)}")

    print("\nRegistered analysis: mean over non-quarantined gate-passing families, "
          "gate denominator 32.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
