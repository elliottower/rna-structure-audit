"""
Sanity checks for batch-1 results before committing to Paper 2 framing.

Check 1 (H2 artifact): Per-nucleotide CAG DI without codon mean-pooling.
    If DI ≈ 0.0035 at both CAG and flank, the flatness is real.
    If per-nucleotide CAG DI > codon-pooled CAG DI, mean-pooling washed out signal.

Check 2 (H3 positional encoding): Left vs right flank DI at layer 1.
    Left flank has identical absolute positions across variants (0-49).
    Right flank has shifted absolute positions (50+r*3 onward).
    If layer-1 ρ comes from right flank, it's positional encoding.

Check 3 (H5 ceiling): bp_prob vs PC1-PC10, not just PC1.
    PC1 explains ~98% variance but may encode position/CLS-adjacency,
    not thermodynamic structure. A higher PC might correlate.
"""
import json
import numpy as np
import torch
from pathlib import Path
from scipy import stats
from sklearn.decomposition import PCA
from tqdm import tqdm
from datetime import datetime
from direction_instability import (
    load_model, get_full_sequence, find_cag_region, get_all_layer_embeddings,
    build_aligned_embeddings, build_aligned_vienna,
    direction_instability_profile,
    REPEAT_COUNTS, WT_REPEATS, FLANK_SIZE, MIN_CODONS, SEED,
)

ROOT = Path(__file__).parent
DATA = ROOT / "data"
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


def _per_variant_di_by_region(aligned_dict, n_cag_positions):
    """DI vs WT split by left flank, CAG, right flank."""
    wt = aligned_dict[WT_REPEATS]
    wt_normed = wt / (wt.norm(dim=-1, keepdim=True) + 1e-8)
    repeats = [r for r in REPEAT_COUNTS if r != WT_REPEATS]
    left, cag, right = {}, {}, {}
    for r in repeats:
        emb = aligned_dict[r]
        emb_normed = emb / (emb.norm(dim=-1, keepdim=True) + 1e-8)
        cos = (wt_normed * emb_normed).sum(dim=-1)
        left[r] = float(1 - cos[:FLANK_SIZE].mean())
        cag[r] = float(1 - cos[FLANK_SIZE:FLANK_SIZE + n_cag_positions].mean())
        right[r] = float(1 - cos[FLANK_SIZE + n_cag_positions:].mean())
    return left, cag, right


def check_h2_artifact(model, seq):
    """Per-nucleotide CAG DI (no codon mean-pool) vs flank DI."""
    print("=" * 60)
    print("CHECK 1: H2 artifact — per-nucleotide CAG DI (no mean-pool)")
    print("=" * 60)

    cag_start, cag_end = find_cag_region(seq)
    flank_left_seq = seq[cag_start - FLANK_SIZE:cag_start]
    flank_right_seq = seq[cag_end:cag_end + FLANK_SIZE]

    n_nuc = MIN_CODONS * 3  # 30 raw nucleotides

    aligned_raw = {}
    for r in tqdm(REPEAT_COUNTS, desc="Per-nucleotide embeddings (L6)"):
        variant_seq = flank_left_seq + "CAG" * r + flank_right_seq
        all_layers = get_all_layer_embeddings(model, variant_seq)
        emb = all_layers[6]
        left = emb[:FLANK_SIZE]
        cag_raw = emb[FLANK_SIZE:FLANK_SIZE + n_nuc]
        right = emb[FLANK_SIZE + r * 3:FLANK_SIZE + r * 3 + FLANK_SIZE]
        aligned_raw[r] = torch.cat([left, cag_raw, right], dim=0)

    # Pairwise DI profile (all variants, per position)
    emb_list = [aligned_raw[r] for r in REPEAT_COUNTS]
    di_profile = direction_instability_profile(emb_list)

    di_left = di_profile[:FLANK_SIZE]
    di_cag = di_profile[FLANK_SIZE:FLANK_SIZE + n_nuc]
    di_right = di_profile[FLANK_SIZE + n_nuc:]
    di_flank = np.concatenate([di_left, di_right])

    # Per-variant DI vs WT, split by region
    left_di, cag_di, right_di = _per_variant_di_by_region(aligned_raw, n_nuc)

    # Spearman ρ(DI vs |r-21|) for each region
    repeats = sorted(left_di.keys())
    deviations = [abs(r - WT_REPEATS) for r in repeats]
    rho_left, _ = stats.spearmanr(deviations, [left_di[r] for r in repeats])
    rho_cag, _ = stats.spearmanr(deviations, [cag_di[r] for r in repeats])
    rho_right, _ = stats.spearmanr(deviations, [right_di[r] for r in repeats])

    results = {
        "profile_mean_cag_raw": float(di_cag.mean()),
        "profile_mean_left": float(di_left.mean()),
        "profile_mean_right": float(di_right.mean()),
        "profile_mean_flank": float(di_flank.mean()),
        "cag_over_flank_ratio": float(di_cag.mean() / (di_flank.mean() + 1e-10)),
        "rho_left_vs_deviation": float(rho_left),
        "rho_cag_vs_deviation": float(rho_cag),
        "rho_right_vs_deviation": float(rho_right),
        "per_variant_left_di": {str(r): v for r, v in left_di.items()},
        "per_variant_cag_di": {str(r): v for r, v in cag_di.items()},
        "per_variant_right_di": {str(r): v for r, v in right_di.items()},
        "n_cag_nucleotides": n_nuc,
    }

    print(f"\n  Pairwise DI profile (all 12 variants, per position):")
    print(f"    Left flank:       {di_left.mean():.6f}  (same tokens, same abs positions)")
    print(f"    CAG raw (30 nuc): {di_cag.mean():.6f}  (same tokens, same abs positions)")
    print(f"    Right flank:      {di_right.mean():.6f}  (same tokens, SHIFTED abs positions)")
    print(f"    CAG/flank ratio:  {di_cag.mean() / (di_flank.mean() + 1e-10):.3f}")

    print(f"\n  ρ(per-variant DI vs |r-21|) by region:")
    print(f"    Left flank:  ρ={rho_left:.4f}")
    print(f"    CAG raw:     ρ={rho_cag:.4f}")
    print(f"    Right flank: ρ={rho_right:.4f}")

    print(f"\n  Per-variant DI vs WT:")
    print(f"  {'r':>5s}  {'|Δ|':>3s}  {'left':>10s}  {'CAG_raw':>10s}  {'right':>10s}")
    for r in repeats:
        print(f"  {r:5d}  {abs(r - WT_REPEATS):3d}  {left_di[r]:10.6f}  {cag_di[r]:10.6f}  {right_di[r]:10.6f}")

    if di_cag.mean() / (di_flank.mean() + 1e-10) > 2.0:
        verdict = "SIGNAL_WASHED: CAG >> flank — codon pooling destroyed real variation"
    elif abs(rho_right) > abs(rho_left) + 0.3 and abs(rho_right) > abs(rho_cag) + 0.3:
        verdict = "POSITIONAL: right flank drives correlation — position encoding artifact"
    elif abs(di_right.mean() - di_left.mean()) / (di_left.mean() + 1e-10) > 3.0:
        verdict = "POSITIONAL: right flank DI >> left flank DI — position shift dominates"
    else:
        verdict = "REAL: DI is uniformly low across all regions — per-position stability is genuine"

    print(f"\n  VERDICT: {verdict}")
    results["verdict"] = verdict
    return results


def check_h3_positional(model, seq):
    """Left vs right flank DI at layer 1 — positional encoding check."""
    print("\n" + "=" * 60)
    print("CHECK 2: H3 positional encoding — left vs right flank at layers 1 and 6")
    print("=" * 60)

    cag_start, cag_end = find_cag_region(seq)
    flank_left_seq = seq[cag_start - FLANK_SIZE:cag_start]
    flank_right_seq = seq[cag_end:cag_end + FLANK_SIZE]

    results = {}

    for layer_idx in [1, 6]:
        aligned = {}
        for r in tqdm(REPEAT_COUNTS, desc=f"Layer {layer_idx} embeddings"):
            variant_seq = flank_left_seq + "CAG" * r + flank_right_seq
            all_layers = get_all_layer_embeddings(model, variant_seq)
            emb = all_layers[layer_idx]
            left = emb[:FLANK_SIZE]
            cag_mean = emb[FLANK_SIZE:FLANK_SIZE + r * 3].mean(dim=0, keepdim=True)
            right = emb[FLANK_SIZE + r * 3:FLANK_SIZE + r * 3 + FLANK_SIZE]
            aligned[r] = torch.cat([left, cag_mean, right], dim=0)

        wt = aligned[WT_REPEATS]
        wt_normed = wt / (wt.norm(dim=-1, keepdim=True) + 1e-8)

        repeats = [r for r in REPEAT_COUNTS if r != WT_REPEATS]
        deviations = [abs(r - WT_REPEATS) for r in repeats]
        left_di, cag_di, right_di, full_di = [], [], [], []

        for r in repeats:
            emb = aligned[r]
            emb_normed = emb / (emb.norm(dim=-1, keepdim=True) + 1e-8)
            cos = (wt_normed * emb_normed).sum(dim=-1)
            left_di.append(float(1 - cos[:FLANK_SIZE].mean()))
            cag_di.append(float(1 - cos[FLANK_SIZE].item()))
            right_di.append(float(1 - cos[FLANK_SIZE + 1:].mean()))
            full_di.append(float(1 - cos.mean()))

        rho_left, p_left = stats.spearmanr(deviations, left_di)
        rho_cag, p_cag = stats.spearmanr(deviations, cag_di)
        rho_right, p_right = stats.spearmanr(deviations, right_di)
        rho_full, p_full = stats.spearmanr(deviations, full_di)

        layer_key = f"layer_{layer_idx}"
        results[layer_key] = {
            "rho_left": float(rho_left), "p_left": float(p_left),
            "rho_cag": float(rho_cag), "p_cag": float(p_cag),
            "rho_right": float(rho_right), "p_right": float(p_right),
            "rho_full": float(rho_full), "p_full": float(p_full),
            "mean_left_di": float(np.mean(left_di)),
            "mean_right_di": float(np.mean(right_di)),
            "right_over_left_di": float(np.mean(right_di) / (np.mean(left_di) + 1e-10)),
        }

        print(f"\n  Layer {layer_idx} — ρ(DI vs |r-21|) by region:")
        print(f"    Left flank  (same abs pos):  ρ={rho_left:.4f}  (mean DI={np.mean(left_di):.6f})")
        print(f"    CAG only    (mean-pooled):   ρ={rho_cag:.4f}")
        print(f"    Right flank (shifted pos):   ρ={rho_right:.4f}  (mean DI={np.mean(right_di):.6f})")
        print(f"    Full aligned:                ρ={rho_full:.4f}")
        print(f"    Right/left DI ratio:         {np.mean(right_di) / (np.mean(left_di) + 1e-10):.2f}x")

    # Verdict based on layer 1
    l1 = results["layer_1"]
    if l1["right_over_left_di"] > 5.0:
        verdict = "POSITIONAL: right-flank DI >> left-flank DI at layer 1 — position encoding dominates"
    elif abs(l1["rho_right"]) > abs(l1["rho_left"]) + 0.3:
        verdict = "POSITIONAL: right-flank ρ >> left-flank ρ — layer-1 sensitivity is position-driven"
    elif abs(l1["rho_left"]) > 0.5:
        verdict = "REAL: left-flank ρ is substantial — layer-1 captures contextual effects beyond position"
    else:
        verdict = "MIXED: neither flank dominates — partial positional encoding effect"

    print(f"\n  VERDICT: {verdict}")
    results["verdict"] = verdict
    return results


def check_h5_multi_pc(aligned_embs, aligned_vienna):
    """bp_prob vs PC1-PC10 — is thermodynamic structure encoded in a non-dominant PC?"""
    print("\n" + "=" * 60)
    print("CHECK 3: H5 ceiling — bp_prob vs PC1-PC10")
    print("=" * 60)

    wt_emb = aligned_embs[WT_REPEATS].numpy()
    wt_vienna = aligned_vienna[WT_REPEATS]
    bp_prob = wt_vienna[:, 0]
    loop_type = wt_vienna[:, 1]

    n_components = min(10, wt_emb.shape[0])
    pca = PCA(n_components=n_components, svd_solver="full")
    scores = pca.fit_transform(wt_emb)

    results = {"bp_prob": {}, "loop_type": {}, "explained_variance": {}}

    print(f"\n  {'PC':>4s} | {'var%':>6s} | {'|r| bp_prob':>11s} {'p':>9s} | {'|r| loop_type':>13s} {'p':>9s}")
    print(f"  {'----':>4s}-+-{'------':>6s}-+-{'-----------':>11s}-{'---------':>9s}-+-{'-------------':>13s}-{'---------':>9s}")

    best_pc_bp, best_r_bp = None, 0.0
    best_pc_lt, best_r_lt = None, 0.0

    for i in range(n_components):
        var = pca.explained_variance_ratio_[i]
        results["explained_variance"][f"PC{i+1}"] = float(var)

        r_bp, p_bp = stats.pearsonr(scores[:, i], bp_prob)
        results["bp_prob"][f"PC{i+1}"] = {"abs_r": float(abs(r_bp)), "p": float(p_bp)}

        r_lt, p_lt = stats.pearsonr(scores[:, i], loop_type)
        results["loop_type"][f"PC{i+1}"] = {"abs_r": float(abs(r_lt)), "p": float(p_lt)}

        mark_bp = " **" if abs(r_bp) > 0.3 and p_bp < 0.01 else " *" if abs(r_bp) > 0.3 else ""
        mark_lt = " **" if abs(r_lt) > 0.3 and p_lt < 0.01 else " *" if abs(r_lt) > 0.3 else ""

        print(f"  PC{i+1:2d} | {var:5.3f} | {abs(r_bp):11.4f} {p_bp:9.6f}{mark_bp} | {abs(r_lt):13.4f} {p_lt:9.6f}{mark_lt}")

        if abs(r_bp) > best_r_bp:
            best_r_bp, best_pc_bp = abs(r_bp), i + 1
        if abs(r_lt) > best_r_lt:
            best_r_lt, best_pc_lt = abs(r_lt), i + 1

    results["best_bp_pc"] = best_pc_bp
    results["best_bp_r"] = float(best_r_bp)
    results["best_lt_pc"] = best_pc_lt
    results["best_lt_r"] = float(best_r_lt)

    if best_r_bp > 0.3:
        verdict = f"CEILING: PC{best_pc_bp} correlates with bp_prob (|r|={best_r_bp:.3f}) — H5 tested the wrong PC"
    else:
        verdict = f"NO CEILING: best bp_prob correlation is |r|={best_r_bp:.3f} at PC{best_pc_bp} — RNA-FM genuinely doesn't linearly encode bp_prob on this window"

    print(f"\n  VERDICT: {verdict}")
    results["verdict"] = verdict
    return results


def main():
    print(f"Batch-1 sanity checks — {TIMESTAMP}")
    print("=" * 60)

    np.random.seed(SEED)
    torch.manual_seed(SEED)

    model = load_model()
    seq = get_full_sequence()

    check1 = check_h2_artifact(model, seq)
    check2 = check_h3_positional(model, seq)

    print("\nBuilding aligned embeddings for H5 multi-PC check...")
    aligned_embs = build_aligned_embeddings(model, seq, REPEAT_COUNTS, layer_idx=6)
    aligned_vienna = build_aligned_vienna(seq, REPEAT_COUNTS)
    check3 = check_h5_multi_pc(aligned_embs, aligned_vienna)

    all_results = {
        "timestamp": TIMESTAMP,
        "check_h2_artifact": check1,
        "check_h3_positional": check2,
        "check_h5_multi_pc": check3,
    }

    out_path = DATA / f"sanity_checks_{TIMESTAMP}.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    print("\n" + "=" * 60)
    print("SANITY CHECK SUMMARY")
    print("=" * 60)
    print(f"  H2 artifact:   {check1['verdict']}")
    print(f"  H3 positional: {check2['verdict']}")
    print(f"  H5 ceiling:    {check3['verdict']}")


if __name__ == "__main__":
    main()
