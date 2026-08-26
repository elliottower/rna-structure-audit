"""Aggregate Phase 6 PS results into confirmatory H1/H2/H3 decisions.

Reads per-model JSON files from data/gpu_results/phase6_compensatory/
and applies the preregistered decision criteria from PREREGISTRATION_PHASE6_V2.md.

Usage:
    uv run python scripts/analyze_phase6_results.py
    uv run python scripts/analyze_phase6_results.py --results-dir data/gpu_results/phase6_compensatory
"""

import argparse
import json
import math
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DIR = ROOT / "data" / "gpu_results" / "phase6_compensatory"
ALPHA = 0.05 / 3  # Bonferroni
N_CONFIRMATORY = 32
H1_FAMILY_THRESHOLD = math.ceil(4 * 0.05 * N_CONFIRMATORY)  # = 7

RNA_MODELS = {"rnafm", "rinalmo", "utrlm", "ernierna", "splicebert"}
DNA_MODELS = {"nt", "hyenadna", "caduceus", "evo", "dnabert2"}
CHARACTER_DNA = {"hyenadna", "caduceus", "evo"}


def load_results(results_dir):
    """Load most recent result file per model."""
    model_results = {}
    for f in sorted(results_dir.rglob("*_phase6_ps*.json"), reverse=True):
        with open(f) as fh:
            data = json.load(fh)
        model = data["model"]
        if model not in model_results:
            model_results[model] = data
    return model_results


def get_confirmatory_families(results):
    """Extract confirmatory (non-skipped, non-quarantined, gate-passing) families."""
    per_rna = results.get("results", {}).get("per_rna", {})
    active = {}
    for name, v in per_rna.items():
        if v.get("skipped") or v.get("quarantined"):
            continue
        if not v.get("positive_control", {}).get("pass", False):
            continue
        active[name] = v
    return active


def test_h1(model_name, active_families):
    """H1: mean PS > 0 AND >= 7 families exceed null."""
    ps_values = [v["best_ps"] for v in active_families.values()]
    n = len(ps_values)

    if n < 2:
        return {"pass": False, "reason": f"only {n} active families", "model": model_name}

    families_exceeding = sum(
        1 for v in active_families.values()
        if v.get("exceeds_null_conservative") is True
    )
    families_no_null = sum(
        1 for v in active_families.values()
        if v.get("null_available") is False
    )

    mean_ps = float(np.mean(ps_values))

    if all(v == 0.0 for v in ps_values):
        return {
            "pass": False, "model": model_name, "mean_ps": 0.0,
            "n_families": n, "reason": "all PS values are zero",
        }

    try:
        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(ps_values, alternative="greater")
    except ValueError:
        return {
            "pass": False, "model": model_name, "mean_ps": mean_ps,
            "n_families": n, "reason": "Wilcoxon test failed (insufficient variation)",
        }

    condition_a = families_exceeding >= H1_FAMILY_THRESHOLD
    condition_b = wilcoxon_p < ALPHA and mean_ps > 0

    return {
        "pass": bool(condition_a and condition_b),
        "model": model_name,
        "mean_ps": mean_ps,
        "n_families": n,
        "families_exceeding_null": families_exceeding,
        "families_no_null": families_no_null,
        "threshold": H1_FAMILY_THRESHOLD,
        "wilcoxon_p": float(wilcoxon_p),
        "condition_a": bool(condition_a),
        "condition_b": bool(condition_b),
    }


def test_h2(model_ps_means):
    """H2: RNA-pretrained PS > DNA-pretrained PS."""
    rna_ps = [v for k, v in model_ps_means.items() if k in RNA_MODELS]
    dna_ps = [v for k, v in model_ps_means.items() if k in DNA_MODELS]

    if len(rna_ps) < 2 or len(dna_ps) < 2:
        return {"pass": False, "reason": "insufficient models in one group"}

    u_stat, u_p = stats.mannwhitneyu(rna_ps, dna_ps, alternative="greater")
    n1, n2 = len(rna_ps), len(dna_ps)
    rank_biserial = (2 * u_stat) / (n1 * n2) - 1

    return {
        "pass": bool(rank_biserial > 0.5),
        "rank_biserial": float(rank_biserial),
        "mann_whitney_p": float(u_p),
        "rna_mean": float(np.mean(rna_ps)),
        "dna_mean": float(np.mean(dna_ps)),
        "rna_models": {k: v for k, v in model_ps_means.items() if k in RNA_MODELS},
        "dna_models": {k: v for k, v in model_ps_means.items() if k in DNA_MODELS},
    }


def test_h3(model_name, active_families, eligible_pairs_data=None):
    """H3: partner-is-max fraction > 1/3 at single-nucleotide precision."""
    h3_data = [v.get("h3_precision", {}) for v in active_families.values()]
    total_max = sum(d.get("partner_max_count", 0) for d in h3_data)
    total_n = sum(d.get("n", 0) for d in h3_data)

    if total_n == 0:
        return {"pass": False, "reason": "no H3-eligible pairs", "model": model_name}

    fraction = total_max / total_n
    p_val = stats.binomtest(total_max, total_n, 1.0 / 3.0, alternative="greater").pvalue

    return {
        "pass": bool(p_val < ALPHA and fraction > 0.5),
        "model": model_name,
        "fraction": float(fraction),
        "partner_max_count": total_max,
        "total_pairs": total_n,
        "p_value": float(p_val),
        "exceeds_half": bool(fraction > 0.5),
        "significant": bool(p_val < ALPHA),
    }


def main():
    parser = argparse.ArgumentParser(description="Analyze Phase 6 PS results")
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_DIR)
    args = parser.parse_args()

    print(f"Loading results from {args.results_dir}")
    all_results = load_results(args.results_dir)

    if not all_results:
        print("No results found.")
        return

    print(f"Found {len(all_results)} models: {list(all_results.keys())}")
    print(f"Confirmatory N = {N_CONFIRMATORY}, H1 threshold = {H1_FAMILY_THRESHOLD}")
    print(f"Bonferroni alpha = {ALPHA:.4f}")
    print()

    model_ps_means = {}
    h1_results = {}
    h3_results = {}

    for model_name, data in sorted(all_results.items()):
        active = get_confirmatory_families(data)
        print(f"--- {model_name} ---")
        print(f"  Active families: {len(active)}")

        if not active:
            print("  No active families (all failed gate or were excluded)")
            continue

        ps_values = [v["best_ps"] for v in active.values()]
        mean_ps = float(np.mean(ps_values))
        model_ps_means[model_name] = mean_ps
        print(f"  Mean PS: {mean_ps:.6f}")

        h1 = test_h1(model_name, active)
        h1_results[model_name] = h1
        print(f"  H1: {'PASS' if h1['pass'] else 'FAIL'} "
              f"(families>{H1_FAMILY_THRESHOLD}: {h1.get('families_exceeding_null', 0)}, "
              f"Wilcoxon p={h1.get('wilcoxon_p', 'N/A'):.4g})")

        h3 = test_h3(model_name, active)
        h3_results[model_name] = h3
        print(f"  H3: {'PASS' if h3['pass'] else 'FAIL'} "
              f"(fraction={h3.get('fraction', 0):.3f}, "
              f"p={h3.get('p_value', 'N/A'):.4g})")
        print()

    print("=" * 60)
    print("CONFIRMATORY DECISIONS")
    print("=" * 60)

    any_h1 = any(h1.get("pass") for h1 in h1_results.values())
    print(f"\nH1 (any model PS > 0 exceeding null): {'PASS' if any_h1 else 'FAIL'}")
    for model_name, h1 in sorted(h1_results.items()):
        if h1.get("pass"):
            print(f"  {model_name}: PASS")

    h2 = test_h2(model_ps_means)
    print(f"\nH2 (RNA > DNA): {'PASS' if h2['pass'] else 'FAIL'}")
    print(f"  Rank-biserial: {h2.get('rank_biserial', 'N/A')}")
    print(f"  Mann-Whitney p: {h2.get('mann_whitney_p', 'N/A')}")

    any_h3 = any(h3.get("pass") for h3 in h3_results.values())
    print(f"\nH3 (single-nucleotide precision): {'PASS' if any_h3 else 'FAIL'}")
    for model_name, h3 in sorted(h3_results.items()):
        if h3.get("pass"):
            print(f"  {model_name}: PASS (fraction={h3['fraction']:.3f})")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary = {
        "timestamp": timestamp,
        "n_models": len(all_results),
        "n_confirmatory_families": N_CONFIRMATORY,
        "h1_threshold": H1_FAMILY_THRESHOLD,
        "alpha": ALPHA,
        "h1_any_pass": any_h1,
        "h1_per_model": h1_results,
        "h2": h2,
        "h3_any_pass": any_h3,
        "h3_per_model": h3_results,
        "model_ps_means": model_ps_means,
    }

    out_path = args.results_dir / f"phase6_confirmatory_decisions_{timestamp}.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSaved decisions to {out_path}")


if __name__ == "__main__":
    main()
