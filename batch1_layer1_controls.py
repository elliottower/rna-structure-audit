"""
Two controls that determine whether batch-1's candidate findings survive.

Control A — Full DI at layer 1 (effective dim ≈ 51, not rank-1 collapsed).
    If layer-1 DI shows CAG-vs-flank structure, then layer-6 flatness is the
    rank-1 collapse eating signal, and the story becomes "DI is informative only
    before the collapse." If layer-1 DI is also flat, per-position stability
    is genuine and the claim holds.

Control B — Length-matched non-CAG insert (CAA*r instead of CAG*r).
    Tests whether left-flank ρ=0.96 at layer 1 is repeat-count semantics or
    just "model knows the sequence got longer." If left-flank DI tracks insert
    length equally for CAA, it's length-sensitivity and the finding evaporates.
"""
import json
import numpy as np
import torch
from pathlib import Path
from scipy import stats
from tqdm import tqdm
from datetime import datetime
from direction_instability import (
    load_model, get_full_sequence, find_cag_region, get_all_layer_embeddings,
    direction_instability_profile,
    REPEAT_COUNTS, WT_REPEATS, FLANK_SIZE, MIN_CODONS, SEED,
)

ROOT = Path(__file__).parent
DATA = ROOT / "data"
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


def _build_raw_aligned(model, flank_left_seq, flank_right_seq, repeat_counts,
                        insert_unit, layer_idx):
    """Build per-nucleotide aligned embeddings (no codon mean-pool).
    insert_unit: e.g. "CAG" or "CAA" — repeated r times."""
    n_nuc = MIN_CODONS * 3
    aligned = {}
    label = f"L{layer_idx} {insert_unit}"
    for r in tqdm(repeat_counts, desc=label):
        variant_seq = flank_left_seq + insert_unit * r + flank_right_seq
        all_layers = get_all_layer_embeddings(model, variant_seq)
        emb = all_layers[layer_idx]
        left = emb[:FLANK_SIZE]
        cag_raw = emb[FLANK_SIZE:FLANK_SIZE + n_nuc]
        right = emb[FLANK_SIZE + r * 3:FLANK_SIZE + r * 3 + FLANK_SIZE]
        aligned[r] = torch.cat([left, cag_raw, right], dim=0)
    return aligned


def _build_meanpool_aligned(model, flank_left_seq, flank_right_seq,
                             repeat_counts, insert_unit, layer_idx):
    """Build mean-pooled aligned embeddings (101 positions)."""
    aligned = {}
    label = f"L{layer_idx} {insert_unit} meanpool"
    for r in tqdm(repeat_counts, desc=label):
        variant_seq = flank_left_seq + insert_unit * r + flank_right_seq
        all_layers = get_all_layer_embeddings(model, variant_seq)
        emb = all_layers[layer_idx]
        left = emb[:FLANK_SIZE]
        insert_mean = emb[FLANK_SIZE:FLANK_SIZE + r * 3].mean(dim=0, keepdim=True)
        right = emb[FLANK_SIZE + r * 3:FLANK_SIZE + r * 3 + FLANK_SIZE]
        aligned[r] = torch.cat([left, insert_mean, right], dim=0)
    return aligned


def _di_vs_wt_by_region(aligned_dict, n_middle_positions):
    """Per-variant DI vs WT, split by left/middle/right."""
    wt = aligned_dict[WT_REPEATS]
    wt_normed = wt / (wt.norm(dim=-1, keepdim=True) + 1e-8)
    repeats = [r for r in REPEAT_COUNTS if r != WT_REPEATS]
    left, mid, right, full = {}, {}, {}, {}
    for r in repeats:
        emb = aligned_dict[r]
        emb_normed = emb / (emb.norm(dim=-1, keepdim=True) + 1e-8)
        cos = (wt_normed * emb_normed).sum(dim=-1)
        left[r] = float(1 - cos[:FLANK_SIZE].mean())
        mid[r] = float(1 - cos[FLANK_SIZE:FLANK_SIZE + n_middle_positions].mean())
        right[r] = float(1 - cos[FLANK_SIZE + n_middle_positions:].mean())
        full[r] = float(1 - cos.mean())
    return left, mid, right, full


def _rho_table(left_di, mid_di, right_di, full_di, label_mid="insert"):
    """Print and return ρ(DI vs |r-21|) for each region."""
    repeats = sorted(left_di.keys())
    deviations = [abs(r - WT_REPEATS) for r in repeats]
    rho_l, p_l = stats.spearmanr(deviations, [left_di[r] for r in repeats])
    rho_m, p_m = stats.spearmanr(deviations, [mid_di[r] for r in repeats])
    rho_r, p_r = stats.spearmanr(deviations, [right_di[r] for r in repeats])
    rho_f, p_f = stats.spearmanr(deviations, [full_di[r] for r in repeats])
    mean_l = np.mean([left_di[r] for r in repeats])
    mean_m = np.mean([mid_di[r] for r in repeats])
    mean_r = np.mean([right_di[r] for r in repeats])

    print(f"    Left flank  (same abs pos):  ρ={rho_l:.4f}  p={p_l:.4f}  mean DI={mean_l:.6f}")
    print(f"    {label_mid:12s} region:        ρ={rho_m:.4f}  p={p_m:.4f}  mean DI={mean_m:.6f}")
    print(f"    Right flank (shifted pos):   ρ={rho_r:.4f}  p={p_r:.4f}  mean DI={mean_r:.6f}")
    print(f"    Full aligned:                ρ={rho_f:.4f}  p={p_f:.4f}")
    print(f"    Right/left DI ratio:         {mean_r / (mean_l + 1e-10):.2f}x")
    print(f"    {label_mid}/left DI ratio:       {mean_m / (mean_l + 1e-10):.2f}x")

    return {
        "rho_left": float(rho_l), "p_left": float(p_l), "mean_left": float(mean_l),
        "rho_mid": float(rho_m), "p_mid": float(p_m), "mean_mid": float(mean_m),
        "rho_right": float(rho_r), "p_right": float(p_r), "mean_right": float(mean_r),
        "rho_full": float(rho_f), "p_full": float(p_f),
        "right_over_left": float(mean_r / (mean_l + 1e-10)),
        "mid_over_left": float(mean_m / (mean_l + 1e-10)),
    }


def control_a_layer1_di(model, flank_left_seq, flank_right_seq):
    """Full DI analysis at layer 1 — per-nucleotide and mean-pooled."""
    print("=" * 60)
    print("CONTROL A: Full DI at layer 1 (effective dim ≈ 51)")
    print("=" * 60)

    n_nuc = MIN_CODONS * 3

    # --- Per-nucleotide (no codon mean-pool) ---
    print("\n  A1: Per-nucleotide CAG DI at layer 1 (no mean-pool)")
    aligned_raw = _build_raw_aligned(
        model, flank_left_seq, flank_right_seq, REPEAT_COUNTS, "CAG", 1)

    emb_list = [aligned_raw[r] for r in REPEAT_COUNTS]
    di_profile = direction_instability_profile(emb_list)
    di_left = di_profile[:FLANK_SIZE]
    di_cag = di_profile[FLANK_SIZE:FLANK_SIZE + n_nuc]
    di_right = di_profile[FLANK_SIZE + n_nuc:]

    print(f"\n    Pairwise DI profile (layer 1):")
    print(f"      Left flank:       {di_left.mean():.6f}")
    print(f"      CAG raw (30 nuc): {di_cag.mean():.6f}")
    print(f"      Right flank:      {di_right.mean():.6f}")
    print(f"      CAG/left ratio:   {di_cag.mean() / (di_left.mean() + 1e-10):.2f}x")
    print(f"      Right/left ratio: {di_right.mean() / (di_left.mean() + 1e-10):.2f}x")

    # Mann-Whitney: CAG vs flank (the actual H2 test at layer 1)
    di_flank = np.concatenate([di_left, di_right])
    u_stat, p_mw = stats.mannwhitneyu(di_cag, di_flank, alternative="greater")
    pooled_std = float(np.sqrt(
        (di_cag.var() * len(di_cag) + di_flank.var() * len(di_flank)) /
        (len(di_cag) + len(di_flank))))
    cohens_d = (di_cag.mean() - di_flank.mean()) / (pooled_std + 1e-8)

    print(f"      Mann-Whitney CAG > flank: U={u_stat:.0f}, p={p_mw:.6f}")
    print(f"      Cohen's d: {cohens_d:.4f}")

    # --- Per-variant DI vs WT by region ---
    print(f"\n  A2: Per-variant DI vs WT at layer 1 (per-nucleotide)")
    left_di, cag_di, right_di, full_di = _di_vs_wt_by_region(aligned_raw, n_nuc)
    raw_rhos = _rho_table(left_di, cag_di, right_di, full_di, label_mid="CAG raw")

    print(f"\n    Per-variant breakdown:")
    repeats = sorted(left_di.keys())
    print(f"    {'r':>5s}  {'|Δ|':>3s}  {'left':>10s}  {'CAG_raw':>10s}  {'right':>10s}")
    for r in repeats:
        print(f"    {r:5d}  {abs(r-WT_REPEATS):3d}  {left_di[r]:10.6f}  {cag_di[r]:10.6f}  {right_di[r]:10.6f}")

    # --- Also mean-pooled at layer 1 (for comparison with layer 6) ---
    print(f"\n  A3: Mean-pooled DI at layer 1 (1 CAG position)")
    aligned_mp = _build_meanpool_aligned(
        model, flank_left_seq, flank_right_seq, REPEAT_COUNTS, "CAG", 1)
    left_mp, cag_mp, right_mp, full_mp = _di_vs_wt_by_region(aligned_mp, 1)
    mp_rhos = _rho_table(left_mp, cag_mp, right_mp, full_mp, label_mid="CAG mean")

    # Verdict
    has_structure = di_cag.mean() / (di_left.mean() + 1e-10) > 1.5 or cohens_d > 0.5
    if has_structure:
        verdict = "LAYER-6 WAS BLIND: layer-1 DI shows CAG-vs-flank structure that rank-1 collapse erased"
    else:
        verdict = "GENUINELY FLAT: layer-1 DI also shows no CAG-vs-flank structure — per-position stability is real"

    print(f"\n  VERDICT: {verdict}")

    return {
        "profile_left": float(di_left.mean()),
        "profile_cag_raw": float(di_cag.mean()),
        "profile_right": float(di_right.mean()),
        "cag_over_left_ratio": float(di_cag.mean() / (di_left.mean() + 1e-10)),
        "right_over_left_ratio": float(di_right.mean() / (di_left.mean() + 1e-10)),
        "h2_mannwhitney_u": float(u_stat),
        "h2_mannwhitney_p": float(p_mw),
        "h2_cohens_d": float(cohens_d),
        "raw_rhos": raw_rhos,
        "meanpool_rhos": mp_rhos,
        "per_variant_left": {str(r): v for r, v in left_di.items()},
        "per_variant_cag": {str(r): v for r, v in cag_di.items()},
        "per_variant_right": {str(r): v for r, v in right_di.items()},
        "verdict": verdict,
    }


def control_b_length_insert(model, flank_left_seq, flank_right_seq):
    """Length-matched non-CAG insert (CAA*r) at layer 1.
    Tests whether left-flank ρ≈0.96 is repeat-count semantics or length sensitivity."""
    print("\n" + "=" * 60)
    print("CONTROL B: Length-matched insert (CAA*r vs CAG*r) at layer 1")
    print("=" * 60)

    results = {}

    for insert_unit, label in [("CAG", "CAG"), ("CAA", "CAA")]:
        print(f"\n  {label} inserts:")
        aligned = _build_meanpool_aligned(
            model, flank_left_seq, flank_right_seq, REPEAT_COUNTS, insert_unit, 1)

        left_di, mid_di, right_di, full_di = _di_vs_wt_by_region(aligned, 1)
        rhos = _rho_table(left_di, mid_di, right_di, full_di, label_mid=label)

        print(f"\n    Per-variant left-flank DI vs WT:")
        repeats = sorted(left_di.keys())
        for r in repeats:
            print(f"      r={r:3d} (|Δ|={abs(r-WT_REPEATS):2d}): left={left_di[r]:.6f}")

        results[label] = {
            "rhos": rhos,
            "per_variant_left": {str(r): v for r, v in left_di.items()},
            "per_variant_mid": {str(r): v for r, v in mid_di.items()},
            "per_variant_right": {str(r): v for r, v in right_di.items()},
        }

    # Compare: if CAA left-flank ρ ≈ CAG left-flank ρ, it's length sensitivity
    cag_rho = results["CAG"]["rhos"]["rho_left"]
    caa_rho = results["CAA"]["rhos"]["rho_left"]
    cag_mean_left = results["CAG"]["rhos"]["mean_left"]
    caa_mean_left = results["CAA"]["rhos"]["mean_left"]

    print(f"\n  Comparison:")
    print(f"    CAG left-flank ρ: {cag_rho:.4f}  (mean DI: {cag_mean_left:.6f})")
    print(f"    CAA left-flank ρ: {caa_rho:.4f}  (mean DI: {caa_mean_left:.6f})")
    print(f"    Δρ:               {abs(cag_rho - caa_rho):.4f}")

    if abs(caa_rho) > 0.8 and abs(cag_rho - caa_rho) < 0.15:
        verdict = "LENGTH SENSITIVITY: CAA produces same left-flank ρ — the model responds to sequence length, not CAG content"
    elif abs(cag_rho) > abs(caa_rho) + 0.3:
        verdict = "REPEAT-SPECIFIC: CAG produces stronger left-flank ρ than CAA — genuine repeat-count semantics"
    elif abs(caa_rho) > abs(cag_rho) + 0.3:
        verdict = "UNEXPECTED: CAA produces stronger ρ than CAG — needs investigation"
    else:
        verdict = "AMBIGUOUS: both inserts produce moderate ρ — can't cleanly separate length from content"

    print(f"\n  VERDICT: {verdict}")
    results["verdict"] = verdict
    results["cag_rho_left"] = float(cag_rho)
    results["caa_rho_left"] = float(caa_rho)
    results["delta_rho"] = float(abs(cag_rho - caa_rho))

    return results


def main():
    print(f"Layer-1 controls — {TIMESTAMP}")
    print("=" * 60)

    np.random.seed(SEED)
    torch.manual_seed(SEED)

    model = load_model()
    seq = get_full_sequence()
    cag_start, cag_end = find_cag_region(seq)
    flank_left_seq = seq[cag_start - FLANK_SIZE:cag_start]
    flank_right_seq = seq[cag_end:cag_end + FLANK_SIZE]

    ctrl_a = control_a_layer1_di(model, flank_left_seq, flank_right_seq)
    ctrl_b = control_b_length_insert(model, flank_left_seq, flank_right_seq)

    all_results = {
        "timestamp": TIMESTAMP,
        "control_a_layer1_di": ctrl_a,
        "control_b_length_insert": ctrl_b,
    }

    out_path = DATA / f"layer1_controls_{TIMESTAMP}.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    print("\n" + "=" * 60)
    print("CONTROL SUMMARY")
    print("=" * 60)
    print(f"  A (layer-1 DI):       {ctrl_a['verdict']}")
    print(f"  B (length insert):    {ctrl_b['verdict']}")

    print("\n  Implications for Paper 2 thesis:")
    if "BLIND" in ctrl_a["verdict"]:
        print("    → Layer-6 DI flatness is metric failure, not per-position stability")
        print("    → Reframe as: 'DI resolution degrades with rank-1 collapse'")
    else:
        print("    → Per-position stability claim survives across layers")

    if "LENGTH" in ctrl_b["verdict"]:
        print("    → Left-flank ρ=0.96 is length-sensitivity, not repeat-count semantics")
        print("    → Drop the 'contextual propagation' headline")
    elif "SPECIFIC" in ctrl_b["verdict"]:
        print("    → Left-flank ρ is genuinely repeat-count-specific — finding survives")


if __name__ == "__main__":
    main()
