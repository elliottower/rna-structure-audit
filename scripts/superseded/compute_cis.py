"""Superseded. Bootstrap confidence intervals, with p-values assigned by lookup.

This script produced the intervals printed in Tables 1 and 2 up to paper_v11.
It is retained because scripts/audit_table12_sources.py needs it to attribute
those printed bounds, and because the Benjamini-Hochberg counts v11 quoted came
from nowhere else. It is not run to produce results.

Two defects retire it. Its exceedance p-values are not measured but assigned
from a lookup table keyed on effect size (0.001 above 2x the null threshold,
0.005 above 1.5x, 0.01 above 1.2x, 0.03 otherwise), and Benjamini-Hochberg
applied to those values is arithmetic on invented inputs. It also wrote
results/bootstrap_cis.json, the same filename scripts/compute_bootstrap_cis.py
writes, in an incompatible flat schema, so running it overwrote the shipped
artifact with a file no other script could read.

scripts/compute_bootstrap_cis.py is the source of record. The output path below
is changed so this script can no longer overwrite it.

Usage:
    uv run --no-project --with numpy --python 3.12 python \
        scripts/superseded/compute_cis.py
"""

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
EXPANDED = REPO / "data" / "gpu_results" / "expanded_rfam"

MODEL_FILES = {
    "ERNIE-RNA": REPO / "results/ernierna_phases15_dinuc.json",
    "RiNALMo": REPO / "results/rinalmo_phases_1_to_5.json",
    "RNA-FM": EXPANDED / "rnafm_expanded_20260716.json",
    "UTR-LM": EXPANDED / "utrlm_phases_1_to_5.json",
    "SpliceBERT": REPO / "results/splicebert_phases15_dinuc.json",
    "NT v2": EXPANDED / "nt_expanded_20260716.json",
    "DNABERT-2": REPO / "results/dnabert2_phases15_dinuc.json",
    "HyenaDNA": EXPANDED / "hyenadna_expanded_20260716.json",
    "Caduceus": EXPANDED / "caduceus_expanded_20260716.json",
    "Evo": EXPANDED / "evo_expanded_20260716.json",
    "ERNIE-RNA untrained": REPO / "results/ernierna_untrained_phases15.json",
}

N_BOOTSTRAP = 10_000
ALPHA = 0.05
RNG = np.random.default_rng(42)


def load_per_family(path):
    with open(path) as f:
        data = json.load(f)
    per_rna = data["mutation_trained"]["per_rna"]
    families = {}
    for fam, vals in per_rna.items():
        if isinstance(vals, dict) and "best_ratio" in vals:
            families[fam] = vals
    return families


def bootstrap_ci(values, n_boot=N_BOOTSTRAP, ci=0.95):
    values = np.array(values, dtype=float)
    n = len(values)
    boot_means = np.empty(n_boot)
    for i in range(n_boot):
        idx = RNG.integers(0, n, size=n)
        boot_means[i] = np.mean(values[idx])
    lo = np.percentile(boot_means, (1 - ci) / 2 * 100)
    hi = np.percentile(boot_means, (1 + ci) / 2 * 100)
    return float(lo), float(hi)


def bh_correct(p_values):
    """Benjamini-Hochberg correction. Returns adjusted p-values."""
    p = np.array(p_values, dtype=float)
    n = len(p)
    sorted_idx = np.argsort(p)
    sorted_p = p[sorted_idx]
    adjusted = np.empty(n)
    adjusted[sorted_idx[-1]] = sorted_p[-1]
    for i in range(n - 2, -1, -1):
        adjusted[sorted_idx[i]] = min(
            adjusted[sorted_idx[i + 1]],
            sorted_p[i] * n / (i + 1),
        )
    return adjusted


def compute_model_stats(model_name, families):
    ratios = []
    nuc_exceed = 0
    nuc_tested = 0
    dinuc_exceed = 0
    dinuc_tested = 0

    raw_p_nuc = []
    raw_p_dinuc = []

    for fam, vals in families.items():
        ratio = vals["best_ratio"]
        ratios.append(ratio)
        nuc_tested += 1

        threshold = vals["nuc_null_95th"]
        if vals.get("exceeds_nuc_null", False):
            nuc_exceed += 1
            excess = ratio / threshold if threshold > 0 else 2.0
            if excess > 2.0:
                raw_p_nuc.append(0.001)
            elif excess > 1.5:
                raw_p_nuc.append(0.005)
            elif excess > 1.2:
                raw_p_nuc.append(0.01)
            else:
                raw_p_nuc.append(0.03)
        else:
            shortfall = threshold / ratio if ratio > 0 else 10.0
            raw_p_nuc.append(min(0.5, 0.05 * shortfall))

        if "dinuc_null_95th" in vals:
            dinuc_tested += 1
            dinuc_thresh = vals["dinuc_null_95th"]
            if vals.get("exceeds_dinuc_null", False):
                dinuc_exceed += 1
                excess_d = ratio / dinuc_thresh if dinuc_thresh > 0 else 2.0
                if excess_d > 2.0:
                    raw_p_dinuc.append(0.001)
                elif excess_d > 1.5:
                    raw_p_dinuc.append(0.005)
                elif excess_d > 1.2:
                    raw_p_dinuc.append(0.01)
                else:
                    raw_p_dinuc.append(0.03)
            else:
                raw_p_dinuc.append(0.5)

    ratios = np.array(ratios)
    mean_ratio = float(np.mean(ratios))
    ci_lo, ci_hi = bootstrap_ci(ratios)

    expected_fp_nuc = ALPHA * nuc_tested

    # BH correction on nuc null p-values
    adj_p_nuc = bh_correct(raw_p_nuc)
    n_survive_bh_nuc = int(np.sum(adj_p_nuc < ALPHA))

    # Dinuc retention CI (bootstrap on binary outcomes)
    dinuc_retention = None
    dinuc_ci = None
    n_survive_bh_dinuc = None
    expected_fp_dinuc = None
    if dinuc_tested > 0:
        dinuc_retention = dinuc_exceed / dinuc_tested
        expected_fp_dinuc = ALPHA * dinuc_tested
        if dinuc_tested >= 3:
            outcomes = np.array([1] * dinuc_exceed + [0] * (dinuc_tested - dinuc_exceed))
            boot_rates = np.empty(N_BOOTSTRAP)
            for i in range(N_BOOTSTRAP):
                idx = RNG.integers(0, len(outcomes), size=len(outcomes))
                boot_rates[i] = np.mean(outcomes[idx])
            dinuc_ci = (float(np.percentile(boot_rates, 2.5)),
                        float(np.percentile(boot_rates, 97.5)))

        if len(raw_p_dinuc) > 1:
            adj_p_dinuc = bh_correct(raw_p_dinuc)
            n_survive_bh_dinuc = int(np.sum(adj_p_dinuc < ALPHA))
        elif len(raw_p_dinuc) == 1:
            n_survive_bh_dinuc = 1 if raw_p_dinuc[0] < ALPHA else 0

    return {
        "model": model_name,
        "n_families": nuc_tested,
        "mean_ratio": round(mean_ratio, 4),
        "ci_95_lo": round(ci_lo, 4),
        "ci_95_hi": round(ci_hi, 4),
        "nuc_exceed_raw": nuc_exceed,
        "nuc_exceed_bh": n_survive_bh_nuc,
        "expected_fp_nuc": round(expected_fp_nuc, 1),
        "dinuc_tested": dinuc_tested,
        "dinuc_exceed_raw": dinuc_exceed,
        "dinuc_retention": round(dinuc_retention, 4) if dinuc_retention is not None else None,
        "dinuc_ci_lo": round(dinuc_ci[0], 4) if dinuc_ci else None,
        "dinuc_ci_hi": round(dinuc_ci[1], 4) if dinuc_ci else None,
        "dinuc_exceed_bh": n_survive_bh_dinuc,
        "expected_fp_dinuc": round(expected_fp_dinuc, 1) if expected_fp_dinuc is not None else None,
    }


def main():
    results = []
    for model_name, path in MODEL_FILES.items():
        families = load_per_family(path)
        stats = compute_model_stats(model_name, families)
        results.append(stats)

    out_path = REPO / "results" / "superseded_compute_cis_flat.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"{'Model':<22} {'Mean ratio':>11} {'95% CI':>16} {'Exc nuc':>8} {'BH nuc':>7} {'E[FP]':>6} {'Dinuc':>10} {'Dinuc CI':>16} {'BH din':>7}")
    print("-" * 120)
    for r in results:
        ci_str = f"[{r['ci_95_lo']:.3f}, {r['ci_95_hi']:.3f}]"
        exc_str = f"{r['nuc_exceed_raw']}/{r['n_families']}"
        bh_str = f"{r['nuc_exceed_bh']}/{r['n_families']}"

        if r["dinuc_tested"] and r["dinuc_tested"] > 0:
            din_str = f"{r['dinuc_exceed_raw']}/{r['dinuc_tested']}"
            if r["dinuc_ci_lo"] is not None:
                din_ci = f"[{r['dinuc_ci_lo']:.2f}, {r['dinuc_ci_hi']:.2f}]"
            else:
                din_ci = "---"
            bh_din = f"{r['dinuc_exceed_bh']}/{r['dinuc_tested']}" if r["dinuc_exceed_bh"] is not None else "---"
        else:
            din_str = "---"
            din_ci = "---"
            bh_din = "---"

        print(f"{r['model']:<22} {r['mean_ratio']:>11.4f} {ci_str:>16} {exc_str:>8} {bh_str:>7} {r['expected_fp_nuc']:>6.1f} {din_str:>10} {din_ci:>16} {bh_din:>7}")

    print(f"\nSaved to {out_path}")
    print(f"\nNote: p-values are approximate. Families exceeding the 95th percentile null")
    print(f"are assigned p=0.04; non-exceeding families are assigned p>0.05 based on")
    print(f"the margin between observed ratio and null threshold. Raw permutation")
    print(f"distributions would give exact p-values but are not stored in result files.")


if __name__ == "__main__":
    main()
