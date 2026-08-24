"""Emit the Phase 6 figures as PREREGISTRATION_PHASE6_V2.md specifies them.

Run:  uv run --no-project --with scipy --python 3.12 python \
          scripts/generate_table5_registered.py

The registration fixes the confirmatory set before the runs: 34 families pass
the 15-WC-pair threshold, two are quarantined as pilot data, and "all
confirmatory hypotheses exclude quarantined families. N = 32 non-quarantined
families passing the stem-filtering criteria."

Mean PS is therefore taken over the 32, not over the 32 that also clear the
positive control. The gate is a separate screening criterion with its own
column, and H1's count criterion is already evaluated among gate-passing
families. Taking the mean over gate-passing families as well would be a second
methodological change the manuscript never made; that variant is printed
alongside as a sensitivity.

audit_table5_aggregation.py establishes that paper_v10.tex computed all of this
over 34: the printed gate numerators reproduce pass/34 for all nine models and
the printed H3 precisions reproduce mean-of-fractions over gate-passing
families out of 34 for all eight reporting one.

Every quantity comes from stored per-family values. No model is re-run.
"""

import json
from pathlib import Path

from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parents[1]
QUARANTINE = {"tRNA_Phe_yeast", "tRNA_Ala_human"}

# (display label, domain, results file, bold in the table).
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

# H2 groups as registered.
RNA_MODELS = ["RNA-FM", "RiNALMo", "UTR-LM", "ERNIE-RNA", "SpliceBERT"]
DNA_MODELS = ["NT~v2", "HyenaDNA", "Caduceus", "Evo", "DNABERT-2"]

SYNTHETIC = {
    "ERNIE-RNA": ("results/audit/phase6_ernierna_20260714_060520/"
                  "ernierna_phase6_ps.json",
                  "results/audit/phase6_synthetic_ernierna_20260715_003418/"
                  "ernierna_synthetic_covariation.json"),
    "RiNALMo": ("results/rinalmo_phase6_ps.json",
                "results/audit/phase6_synthetic_rinalmo_20260715_002812/"
                "rinalmo_synthetic_covariation.json"),
}


def scored(rel: str) -> dict:
    path = REPO / rel
    assert path.exists(), f"missing {rel}"
    results = json.loads(path.read_text()).get("results", {})
    return {name: body for name, body in results.get("per_rna", {}).items()
            if isinstance(body, dict) and not body.get("skipped")
            and isinstance(body.get("best_ps"), (int, float))}


def gate_pass(body: dict) -> bool:
    control = body.get("positive_control")
    return bool(control and control.get("pass"))


def mean(values: list) -> float | None:
    return sum(values) / len(values) if values else None


def render_ps(value: float) -> str:
    """Match the table's style: fixed decimals, exponent when below 1e-4."""
    if value == 0:
        return "0.0000"
    if abs(value) < 1e-4:
        mantissa, exponent = f"{value:.1e}".split("e")
        return f"${mantissa} \\times 10^{{{int(exponent)}}}$"
    return f"{value:.4f}"


def summarize(rel: str) -> dict:
    """Every Table 5 quantity for one model, with the quarantine applied."""
    return summarize_entries(scored(rel))


def summarize_entries(entries: dict) -> dict:
    """The registered Table 5 aggregation, over already-loaded per-family entries.

    Split out so generate_results_tables.py computes the repaired-panel table
    through the same code audit_table5_aggregation.py checked against the
    manuscript, rather than through a second implementation that agrees with it
    until one of them is edited.
    """
    nonq = {name: body for name, body in entries.items() if name not in QUARANTINE}
    passing = [body for body in nonq.values() if gate_pass(body)]
    fractions = [body["h3_precision"]["fraction"] for body in passing
                 if isinstance(body.get("h3_precision"), dict)
                 and isinstance(body["h3_precision"].get("fraction"), (int, float))]
    return dict(
        mean_ps=mean([body["best_ps"] for body in nonq.values()]),
        mean_ps_gated=mean([body["best_ps"] for body in passing]),
        eligible=len(nonq),
        gate=len(passing),
        exceed=sum(1 for body in passing if body.get("exceeds_null_primary")),
        h3=mean(fractions),
    )


def synthetic_pair(natural_rel: str, synthetic_rel: str) -> tuple:
    """Natural and synthetic mean PS on matched, non-quarantined families."""
    natural = mean([body["best_ps"] for name, body in scored(natural_rel).items()
                    if name not in QUARANTINE])
    # Synthetic sequences are named <family>_syn<k>; drop those derived from a
    # quarantined family so the two sides cover the same families.
    synthetic = mean([body["best_ps"]
                      for name, body in scored(synthetic_rel).items()
                      if name.rsplit("_syn", 1)[0] not in QUARANTINE
                      and gate_pass(body)])
    return natural, synthetic


def main() -> int:
    stats = {}
    print("Registered analysis: N = 32 non-quarantined families.\n")
    print(f"{'model':<22} {'mean PS':>22} {'gate':>7} {'>null':>7} {'H3':>7} "
          f"{'gated mean':>12}")
    print("-" * 82)
    latex = []
    for label, domain, rel, bold in ROWS:
        row = summarize(rel)
        key = label.replace(r"\textbf{", "").replace("}", "").split(" (")[0]
        stats[key] = row

        ps_cell = "---" if row["mean_ps"] is None else render_ps(row["mean_ps"])
        gate_cell = f"{row['gate']}/{row['eligible']}"
        null_cell = "---" if not row["gate"] else f"{row['exceed']}/{row['gate']}"
        h3_cell = "---" if row["h3"] is None else f"{row['h3']:.3f}"
        if bold:
            ps_cell = f"\\textbf{{{ps_cell}}}"
            null_cell = f"\\textbf{{{null_cell}}}"
        gated = ("--" if row["mean_ps_gated"] is None
                 else f"{row['mean_ps_gated']:.6f}")
        print(f"{key:<22} {ps_cell:>22} {gate_cell:>7} {null_cell:>7} "
              f"{h3_cell:>7} {gated:>12}")
        latex.append(f"{label:<26} & {domain} & {ps_cell} & {gate_cell} & "
                     f"{null_cell} & {h3_cell} \\\\")

    print("\n--- LaTeX rows ---")
    for line in latex:
        print(line)

    rinalmo, ernie = stats["RiNALMo"], stats["ERNIE-RNA"]
    caduceus, evo = stats["Caduceus"], stats["Evo"]
    splice, hyena = stats["SpliceBERT"], stats["HyenaDNA"]
    untrained = stats["ERNIE-RNA untrained"]

    print("\n--- separation from the leaders (lower leader = ERNIE-RNA) ---")
    for name in ["Caduceus", "Evo", "SpliceBERT", "HyenaDNA", "RNA-FM",
                 "UTR-LM", "NT~v2"]:
        value = stats[name]["mean_ps"]
        if value and value > 0:
            print(f"  {name:<12} {ernie['mean_ps'] / value:10.0f}x below "
                  f"ERNIE-RNA, {rinalmo['mean_ps'] / value:10.0f}x below RiNALMo")

    print("\n--- H2: RNA vs DNA rank-biserial ---")
    rna = [stats[m]["mean_ps"] for m in RNA_MODELS]
    dna = [stats[m]["mean_ps"] for m in DNA_MODELS]
    u_stat, p_value = mannwhitneyu(rna, dna, alternative="two-sided")
    rank_biserial = 2 * u_stat / (len(rna) * len(dna)) - 1
    print(f"  rb = {rank_biserial:+.2f}, p = {p_value:.3f}")

    print("\n--- synthetic covariation control ---")
    for model, (natural_rel, synthetic_rel) in SYNTHETIC.items():
        natural, synth = synthetic_pair(natural_rel, synthetic_rel)
        print(f"  {model:<10} natural {natural:.3f} -> synthetic {synth:.3f} "
              f"({(1 - synth / natural) * 100:.0f}% reduction)")

    print("\n--- quarantined families: per-family PS and both means ---")
    print("LaTeX rows: model & tRNA_Phe & tRNA_Ala & PS(N=32) & PS(N=34)")
    for label, _domain, rel, _bold in ROWS:
        key = label.replace(r"\textbf{", "").replace("}", "").split(" (")[0]
        entries = scored(rel)
        held = {name: entries[name]["best_ps"]
                for name in QUARANTINE if name in entries}
        if not held:
            continue
        n32 = mean([b["best_ps"] for n, b in entries.items()
                    if n not in QUARANTINE])
        n34 = mean([b["best_ps"] for b in entries.values()])
        if key == "ERNIE-RNA untrained":
            continue
        print(f"{key.replace('~', ' '):<26} & {render_ps(held['tRNA_Phe_yeast'])} & "
              f"{render_ps(held['tRNA_Ala_human'])} & {render_ps(n32)} & "
              f"{render_ps(n34)} \\\\")

    print("\n--- confirmatory verdicts with the quarantined families restored ---")
    for key, rel in [("RiNALMo", ROWS[0][2]), ("ERNIE-RNA", ROWS[1][2])]:
        entries = scored(rel)
        passing = [b for b in entries.values() if gate_pass(b)]
        exceed = sum(1 for b in passing if b.get("exceeds_null_primary"))
        fractions = [b["h3_precision"]["fraction"] for b in passing
                     if isinstance(b.get("h3_precision"), dict)]
        print(f"  {key:<10} N=34: {exceed}/{len(passing)} exceed null "
              f"(criterion >= 7), H3 = {mean(fractions):.3f} (criterion > 1/3)")
    rna34, dna34 = [], []
    for label, _domain, rel, _bold in ROWS:
        key = label.replace(r"\textbf{", "").replace("}", "").split(" (")[0]
        if key not in RNA_MODELS + DNA_MODELS:
            continue
        value = mean([b["best_ps"] for b in scored(rel).values()])
        (rna34 if key in RNA_MODELS else dna34).append(value)
    u34, p34 = mannwhitneyu(rna34, dna34, alternative="two-sided")
    print(f"  H2 N=34: rb = {2 * u34 / (len(rna34) * len(dna34)) - 1:+.2f}, "
          f"p = {p34:.3f} (criterion rb > 0.5)")

    print("\n--- prose figures ---")
    print(f"  RiNALMo   PS = {rinalmo['mean_ps']:.4f}, "
          f"{rinalmo['exceed']}/{rinalmo['gate']} exceed null, "
          f"H3 = {rinalmo['h3']:.3f}")
    print(f"  ERNIE-RNA PS = {ernie['mean_ps']:.4f}, "
          f"{ernie['exceed']}/{ernie['gate']} exceed null, H3 = {ernie['h3']:.3f}")
    print(f"  Caduceus  PS = {caduceus['mean_ps']:.4f}, "
          f"{caduceus['exceed']}/{caduceus['gate']} exceed null, "
          f"H3 = {caduceus['h3']:.3f}")
    print(f"  SpliceBERT {splice['exceed']}/{splice['gate']} exceed null, "
          f"PS = {splice['mean_ps']:.1e}")
    print(f"  HyenaDNA  H3 = {hyena['h3']:.3f}")
    print(f"  untrained PS = {untrained['mean_ps']:.1e}, "
          f"gate {untrained['gate']}/{untrained['eligible']} "
          f"(vs {ernie['gate']}/{ernie['eligible']} trained), "
          f"{untrained['exceed']}/{untrained['gate']} exceed null, "
          f"{ernie['mean_ps'] / untrained['mean_ps']:.2e}x below trained")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
