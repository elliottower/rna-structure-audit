"""Compute bootstrap 95% confidence intervals for RNA structure awareness metrics.

Loads per-model result JSON files, extracts per-family metrics, and computes
bootstrap CIs for mutation sensitivity, null-exceedance fractions, attention
contact correlation, and probing accuracy. It computes no p-values: the result
files store a binary exceedance against each family's own 95th-percentile null
and not the permutation distribution behind it, so per-family significance is
reported as an exceedance count beside the count expected under independence.

This is the source of record for every interval Tables 1 and 2 print; see
scripts/verify_paper_rung12_figures.py.

Usage:
    uv run --with numpy --with scipy --with statsmodels --with tqdm \
        python scripts/compute_bootstrap_cis.py
"""

import json
import datetime
from pathlib import Path

import numpy as np
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# The repaired panel. Until 2026-08-27 these pointed at the pre-repair trees --
# `results/*_phases15_dinuc.json` and `data/gpu_results/expanded_rfam/` -- so the
# artifact every verification script compares the manuscript against was built
# from a panel the manuscript no longer reports. The mismatch was invisible
# because the manuscript resolver was also stale (D25) and the two wrong files
# agreed with each other.
REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
PANEL = RESULTS / "repaired_panel_v3"


def _panel(key: str) -> Path:
    return PANEL / key / f"{key}_phases_1_to_5.json"


MODEL_FILES = {
    "ERNIE-RNA": _panel("ernierna"),
    "RiNALMo": _panel("rinalmo"),
    "RNA-FM": _panel("rnafm"),
    "UTR-LM": _panel("utrlm"),
    "SpliceBERT": _panel("splicebert"),
    "NT v2": _panel("nt"),
    "DNABERT-2": _panel("dnabert2"),
    "HyenaDNA": _panel("hyenadna"),
    "Caduceus": _panel("caduceus"),
    "Evo": _panel("evo"),
    "ERNIE-RNA (untrained)": _panel("ernierna_untrained"),
}

OUTPUT_PATH = RESULTS / "bootstrap_cis.json"
N_BOOTSTRAP = 10_000
SEED = 42


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _is_scored(entry: dict) -> bool:
    """Return True if the per-family entry is a real result (not skipped)."""
    return isinstance(entry, dict) and "skipped" not in entry


def extract_mutation_data(data: dict) -> tuple[list[float], list[bool], list[bool]]:
    """Extract per-family mutation metrics, filtering skipped entries.

    Returns (best_ratios, exceeds_nuc, exceeds_dinuc).
    exceeds_dinuc may be shorter than the others when dinuc data is absent.
    """
    per_rna = data.get("mutation_trained", {}).get("per_rna", {})
    ratios = []
    exceeds_nuc = []
    exceeds_dinuc = []
    for fam, entry in per_rna.items():
        if not _is_scored(entry):
            continue
        ratios.append(entry["best_ratio"])
        exceeds_nuc.append(entry["exceeds_nuc_null"])
        if "exceeds_dinuc_null" in entry:
            exceeds_dinuc.append(entry["exceeds_dinuc_null"])
    return ratios, exceeds_nuc, exceeds_dinuc


def extract_attention_data(data: dict) -> list[float]:
    """Extract per-family best_corr from attention_trained, filtering skipped."""
    at = data.get("attention_trained", {})
    if at.get("skipped"):
        return []
    per_rna = at.get("per_rna", {})
    corrs = []
    for fam, entry in per_rna.items():
        if not _is_scored(entry):
            continue
        corrs.append(entry["best_corr"])
    return corrs


def extract_probing(data: dict) -> float | None:
    """Extract probing best_accuracy (single scalar)."""
    probing = data.get("probing", {})
    return probing.get("best_accuracy")


def bootstrap_ci(
    values: np.ndarray,
    stat_fn,
    rng: np.random.Generator,
    n_boot: int = N_BOOTSTRAP,
    alpha: float = 0.05,
) -> dict:
    """Compute bootstrap percentile CI for a statistic.

    Returns dict with point_estimate, ci_lower, ci_upper, n.
    """
    n = len(values)
    point = float(stat_fn(values))
    boot_stats = np.empty(n_boot)
    for i in range(n_boot):
        sample = rng.choice(values, size=n, replace=True)
        boot_stats[i] = stat_fn(sample)
    lo = float(np.percentile(boot_stats, 100 * alpha / 2))
    hi = float(np.percentile(boot_stats, 100 * (1 - alpha / 2)))
    return {"point_estimate": point, "ci_lower": lo, "ci_upper": hi, "n": n}


def bootstrap_fraction_ci(
    flags: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = N_BOOTSTRAP,
    alpha: float = 0.05,
) -> dict:
    """Bootstrap CI for fraction of True values."""
    return bootstrap_ci(flags.astype(float), np.mean, rng, n_boot, alpha)


def compute_multiple_testing_stats(data: dict) -> dict:
    """Compute multiple testing statistics for nucleotide null exceedance.

    We only have binary exceed/not-exceed at alpha=0.05 (the 95th percentile
    threshold), not continuous p-values. Report expected false positive rates
    under independence and Bonferroni bounds.

    With m tests at alpha=0.05:
    - Expected FP under independence: m * 0.05
    - Bonferroni threshold: 0.05 / m (we can't apply this since p is binary)
    - If k families exceed: expected FP = min(k, m * 0.05)
    """
    per_rna = data.get("mutation_trained", {}).get("per_rna", {})
    families = []
    exceeds = []
    for fam, entry in per_rna.items():
        if not _is_scored(entry):
            continue
        families.append(fam)
        exceeds.append(entry["exceeds_nuc_null"])

    m = len(families)
    k = sum(exceeds)
    expected_fp_indep = m * 0.05
    bonferroni_alpha = 0.05 / m if m > 0 else 0.0

    return {
        "n_tested": m,
        "n_exceeding": k,
        "alpha": 0.05,
        "expected_fp_under_independence": round(expected_fp_indep, 2),
        "bonferroni_alpha": round(bonferroni_alpha, 4),
        "excess_over_expected_fp": k - expected_fp_indep if m > 0 else 0,
        "note": (
            "Binary test: each family compared to its own 95th-percentile "
            "permutation null. Expected FP = m * alpha under independence. "
            "Continuous p-values would require storing full permutation distributions."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(out_path: Path = OUTPUT_PATH) -> None:
    print(f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] Starting bootstrap CI computation")
    print(f"  N_BOOTSTRAP = {N_BOOTSTRAP}, SEED = {SEED}")

    rng = np.random.default_rng(SEED)
    results = {}

    models = list(MODEL_FILES.items())
    for model_name, filepath in tqdm(models, desc="Models", unit="model"):
        if not filepath.exists():
            print(f"  WARNING: {filepath} not found, skipping {model_name}")
            continue

        with open(filepath) as f:
            data = json.load(f)

        model_results = {"file": str(filepath.relative_to(REPO))}

        # --- Mutation sensitivity ---
        ratios, exceeds_nuc, exceeds_dinuc = extract_mutation_data(data)
        if ratios:
            arr_ratios = np.array(ratios)
            model_results["mutation_sensitivity"] = bootstrap_ci(
                arr_ratios, np.mean, rng
            )

            arr_nuc = np.array(exceeds_nuc)
            model_results["frac_exceeds_nuc_null"] = bootstrap_fraction_ci(
                arr_nuc, rng
            )

            if exceeds_dinuc:
                arr_dinuc = np.array(exceeds_dinuc)
                model_results["frac_exceeds_dinuc_null"] = bootstrap_fraction_ci(
                    arr_dinuc, rng
                )
            else:
                model_results["frac_exceeds_dinuc_null"] = None

            model_results["multiple_testing"] = compute_multiple_testing_stats(data)
        else:
            model_results["mutation_sensitivity"] = None
            model_results["frac_exceeds_nuc_null"] = None
            model_results["frac_exceeds_dinuc_null"] = None
            model_results["multiple_testing"] = None

        # --- Attention contact ---
        corrs = extract_attention_data(data)
        if corrs:
            arr_corrs = np.array(corrs)
            model_results["attention_contact_rho"] = bootstrap_ci(
                arr_corrs, np.mean, rng
            )
        else:
            model_results["attention_contact_rho"] = None

        # --- Probing ---
        probe_acc = extract_probing(data)
        if probe_acc is not None:
            model_results["probing_accuracy"] = {
                "point_estimate": probe_acc,
                "note": "single scalar, no bootstrap",
            }
        else:
            model_results["probing_accuracy"] = None

        results[model_name] = model_results

    # --- Save ---
    output = {
        "description": "Bootstrap 95% CIs for RNA structure awareness metrics",
        "n_bootstrap": N_BOOTSTRAP,
        "seed": SEED,
        "timestamp": datetime.datetime.now().isoformat(),
        "models": results,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] Results saved to {out_path}")

    # --- Print table ---
    print("\n" + "=" * 120)
    print(f"{'Model':<24} {'Sens. Ratio':>16} {'Nuc Null %':>16} {'Dinuc Null %':>16} "
          f"{'Attn rho':>16} {'Probe Acc':>10}   {'Exceed (E[FP])'}")
    print("-" * 120)

    for model_name in MODEL_FILES:
        if model_name not in results:
            continue
        r = results[model_name]

        def _fmt_ci(entry: dict | None) -> str:
            if entry is None:
                return "---"
            pe = entry["point_estimate"]
            lo = entry["ci_lower"]
            hi = entry["ci_upper"]
            return f"{pe:.3f} [{lo:.3f},{hi:.3f}]"

        def _fmt_pct_ci(entry: dict | None) -> str:
            if entry is None:
                return "---"
            pe = entry["point_estimate"] * 100
            lo = entry["ci_lower"] * 100
            hi = entry["ci_upper"] * 100
            return f"{pe:.1f}% [{lo:.1f},{hi:.1f}]"

        sens = _fmt_ci(r.get("mutation_sensitivity"))
        nuc = _fmt_pct_ci(r.get("frac_exceeds_nuc_null"))
        dinuc = _fmt_pct_ci(r.get("frac_exceeds_dinuc_null"))
        attn = _fmt_ci(r.get("attention_contact_rho"))

        probe = r.get("probing_accuracy")
        probe_str = f"{probe['point_estimate']:.3f}" if probe else "---"

        mt = r.get("multiple_testing")
        if mt:
            mt_str = f"{mt['n_exceeding']}/{mt['n_tested']} (E[FP]={mt['expected_fp_under_independence']:.1f})"
        else:
            mt_str = "---"

        print(f"{model_name:<24} {sens:>16} {nuc:>16} {dinuc:>16} "
              f"{attn:>16} {probe_str:>10}   {mt_str}")

    print("=" * 120)
    print(f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] Done.")


if __name__ == "__main__":
    main()
