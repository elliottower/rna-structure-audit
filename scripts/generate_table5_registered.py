"""Emit Table 5 of the manuscript as PREREGISTRATION_PHASE6_V2.md specifies it.

Run:  uv run --no-project --python 3.12 python scripts/generate_table5_registered.py

The registration fixes the confirmatory set before the runs: 34 families pass
the 15-WC-pair threshold, two are quarantined as pilot data, and "all
confirmatory hypotheses exclude quarantined families. N = 32". The positive
control is "a screening gate applied per model per family", and a model-family
pair that fails it is excluded from confirmatory analysis.

audit_table5_aggregation.py establishes that paper_v10.tex Table 5 was computed
over all 34: the printed gate numerators reproduce pass/34 for all nine models,
and the printed H3 precisions reproduce the mean of per-family fractions over
gate-passing families out of 34 for all eight that report one.

Each results file already carries `families_quarantined` and a per-family
`quarantined` flag, and its stored `mean_best_ps` already equals the
non-quarantined gate-passing mean, so the registered analysis is recovered from
stored per-family values. No model is re-run.

Prints the corrected LaTeX rows and the figures the surrounding prose quotes.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# (display name, domain, results file). Order follows Table 5.
ROWS = [
    (r"\textbf{RiNALMo} (650M)", "RNA",
     "results/rinalmo_phase6_ps.json", True),
    (r"\textbf{ERNIE-RNA} (86M)", "RNA",
     "results/audit/phase6_ernierna_20260714_060520/ernierna_phase6_ps.json", True),
    ("ERNIE-RNA untrained", "RNA",
     "results/audit/phase6_ernierna_untrained_20260714_172412/"
     "ernierna_untrained_phase6_ps.json", False),
    ("Caduceus (14M)", "DNA",
     "data/gpu_results/expanded_rfam/caduceus_phase6_ps.json", False),
    ("Evo (7B)", "DNA",
     "results/audit/phase6_evo_20260715_073218/evo_phase6_ps.json", False),
    ("SpliceBERT (19M)", "RNA",
     "results/audit/phase6_splicebert_20260714_060520/splicebert_phase6_ps.json", False),
    ("HyenaDNA (5.4M)", "DNA",
     "results/audit/phase6_hyenadna_20260714_060542/hyenadna_phase6_ps.json", False),
    ("RNA-FM (99M)", "RNA",
     "results/audit/phase6_rnafm_20260714_065238/rnafm_phase6_ps.json", False),
    (r"UTR-LM ($\sim$2M)", "RNA",
     "results/audit/phase6_utrlm_20260714_060516/utrlm_phase6_ps.json", False),
    ("NT~v2 (56M)", "DNA",
     "results/audit/phase6_nt_20260714_082334/nt_phase6_ps.json", False),
    ("DNABERT-2 (117M)", "DNA",
     "results/phase6_dnabert2_natural.json", False),
]


def load(rel: str) -> dict:
    path = REPO / rel
    assert path.exists(), f"missing {rel}"
    results = json.loads(path.read_text()).get("results", {})
    return {name: body for name, body in results.get("per_rna", {}).items()
            if isinstance(body, dict) and not body.get("skipped")
            and isinstance(body.get("best_ps"), (int, float))}


def gate_pass(body: dict) -> bool:
    control = body.get("positive_control")
    return bool(control and control.get("pass"))


def render_ps(value: float) -> str:
    """Match the table's existing style: fixed decimals, exponent when tiny."""
    if abs(value) < 1e-4:
        mantissa = f"{value:.1e}".split("e")[0]
        exponent = int(f"{value:.1e}".split("e")[1])
        return f"${mantissa} \\times 10^{{{exponent}}}$"
    return f"{value:.4f}"


def main() -> int:
    print("Registered analysis: non-quarantined families, positive-control gate "
          "applied.\n")
    print(f"{'model':<26} {'mean PS':>12} {'gate':>8} {'>null':>8} {'H3':>7}")
    print("-" * 66)
    latex, prose = [], {}
    for label, domain, rel, bold in ROWS:
        entries = load(rel)
        nonq = {n: b for n, b in entries.items() if n not in QUARANTINE}
        passing = [b for b in nonq.values() if gate_pass(b)]

        n_elig, n_gate = len(nonq), len(passing)
        mean_ps = sum(b["best_ps"] for b in passing) / n_gate if n_gate else None
        exceed = sum(1 for b in passing if b.get("exceeds_null_primary"))
        fractions = [b["h3_precision"]["fraction"] for b in passing
                     if isinstance(b.get("h3_precision"), dict)
                     and isinstance(b["h3_precision"].get("fraction"), (int, float))]
        h3 = sum(fractions) / len(fractions) if fractions else None

        key = label.replace("\\textbf{", "").replace("}", "").split(" (")[0]
        prose[key] = dict(mean_ps=mean_ps, gate=n_gate, elig=n_elig,
                          exceed=exceed, h3=h3)

        ps_cell = "---" if mean_ps is None else render_ps(mean_ps)
        gate_cell = f"{n_gate}/{n_elig}"
        null_cell = "---" if not n_gate else f"{exceed}/{n_gate}"
        h3_cell = "---" if h3 is None else f"{h3:.3f}"
        if bold:
            ps_cell, null_cell = f"\\textbf{{{ps_cell}}}", f"\\textbf{{{null_cell}}}"
        print(f"{key:<26} {ps_cell:>12} {gate_cell:>8} {null_cell:>8} {h3_cell:>7}")
        latex.append(f"{label:<26} & {domain} & {ps_cell} & {gate_cell} & "
                     f"{null_cell} & {h3_cell} \\\\")

    print("\n--- LaTeX rows ---")
    for line in latex:
        print(line)

    print("\n--- figures the prose quotes ---")
    r, e = prose["RiNALMo"], prose["ERNIE-RNA"]
    u, c, s, h = prose["ERNIE-RNA untrained"], prose["Caduceus"], \
        prose["SpliceBERT"], prose["HyenaDNA"]
    print(f"RiNALMo   PS = {r['mean_ps']:.3f}, {r['exceed']}/{r['gate']} exceed null, "
          f"H3 = {r['h3']:.3f}")
    print(f"ERNIE-RNA PS = {e['mean_ps']:.3f}, {e['exceed']}/{e['gate']} exceed null, "
          f"H3 = {e['h3']:.3f}")
    print(f"untrained PS = {u['mean_ps']:.2e}, gate {u['gate']}/{u['elig']} "
          f"(vs {e['gate']}/{e['elig']} trained), {u['exceed']}/{u['gate']} exceed null")
    print(f"Caduceus  PS = {c['mean_ps']:.4f}, {c['exceed']}/{c['gate']} exceed null, "
          f"H3 = {c['h3']:.3f}")
    print(f"SpliceBERT PS = {s['mean_ps']:.2e}, {s['exceed']}/{s['gate']} exceed null")
    print(f"HyenaDNA  H3 = {h['h3']:.3f}")
    print(f"\nleader/third ratio: RiNALMo/Caduceus = "
          f"{r['mean_ps'] / c['mean_ps']:.0f}x, "
          f"ERNIE-RNA/Caduceus = {e['mean_ps'] / c['mean_ps']:.0f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
