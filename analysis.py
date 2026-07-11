"""
Full analysis pipeline for RNA-FM on HTT mRNA.
Runs: bracket norm, transportability, layer-wise PCA, weight factorization.
Saves all results to data/ and figures/.
"""
import os
import re
import json
import numpy as np
import torch
from pathlib import Path
from scipy.linalg import subspace_angles
from sklearn.decomposition import PCA
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
from transformers import BertConfig, BertModel
from tqdm import tqdm

ROOT = Path(__file__).parent
DATA = ROOT / "data"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

NUC_TO_ID = {"A": 5, "C": 6, "G": 7, "U": 8}
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


def load_model():
    state = torch.load(ROOT / "pretrained" / "pytorch_model.bin", map_location="cpu")
    config = BertConfig(
        vocab_size=28, hidden_size=640, num_hidden_layers=12,
        num_attention_heads=20, intermediate_size=5120,
        max_position_embeddings=state["model.embeddings.position_embeddings.weight"].shape[0],
    )
    clean_state = {}
    for k, v in state.items():
        if k.startswith("model."):
            nk = k[len("model."):]
            nk = nk.replace("attention.layer_norm", "attention.output.LayerNorm")
            if ".layer_norm." in nk and "attention" not in nk:
                nk = nk.replace("layer_norm", "output.LayerNorm")
            clean_state[nk] = v
    model = BertModel(config)
    model.load_state_dict(clean_state, strict=False)
    model.eval()
    return model, config


def tokenize(seq, max_len=510):
    seq = seq[:max_len]
    tokens = [2] + [NUC_TO_ID.get(c, 1) for c in seq] + [3]
    return torch.tensor([tokens])


def get_embeddings(model, seq, max_len=510):
    input_ids = tokenize(seq, max_len)
    with torch.no_grad():
        out = model(input_ids, output_hidden_states=True)
    return out.last_hidden_state[0, 1:-1, :], [h[0, 1:-1, :] for h in out.hidden_states]


def get_full_sequence():
    record = SeqIO.read(DATA / "HTT_NM_002111.7.fasta", "fasta")
    return str(record.seq).upper().replace("T", "U")


def find_cag_region(seq):
    matches = [(m.start(), m.end()) for m in re.finditer(r"(CAG){5,}", seq)]
    return max(matches, key=lambda x: x[1] - x[0])


# ── 1. BRACKET NORM ──
def run_bracket_norm(seq, model):
    print("\n=== BRACKET NORM ANALYSIS ===")
    cag_start, cag_end = find_cag_region(seq)
    results = {}

    window_sizes = [100, 200, 300, 400, 500]
    center = (cag_start + cag_end) // 2

    norms_raw = []
    norms_corrected = []
    for ws in tqdm(window_sizes, desc="Bracket norm windows"):
        start = max(0, center - ws // 2)
        end = min(len(seq), center + ws // 2)
        window = seq[start:end]
        emb, _ = get_embeddings(model, window)
        n = emb.shape[0]

        raw_norm = torch.norm(emb, dim=-1).mean().item()
        corrected_norm = raw_norm / np.sqrt(n)

        norms_raw.append(raw_norm)
        norms_corrected.append(corrected_norm)
        print(f"  window={n}: raw_norm={raw_norm:.4f}, BN_corrected={corrected_norm:.4f}")

    results["window_sizes"] = window_sizes
    results["norms_raw"] = norms_raw
    results["norms_corrected"] = norms_corrected

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(window_sizes, norms_raw, "o-", color="#2563eb", linewidth=2)
    axes[0].set_xlabel("Window size (nt)")
    axes[0].set_ylabel("Mean embedding norm")
    axes[0].set_title("Raw embedding norm vs sequence length")

    axes[1].plot(window_sizes, norms_corrected, "s-", color="#dc2626", linewidth=2)
    axes[1].set_xlabel("Window size (nt)")
    axes[1].set_ylabel("BN-corrected norm (÷ √n)")
    axes[1].set_title("Bracket-norm corrected")
    plt.tight_layout()
    plt.savefig(FIGS / "bracket_norm.png", dpi=200)
    plt.close()
    print("  Saved figures/bracket_norm.png")

    # Per-dimension bracket norm
    emb_full, _ = get_embeddings(model, seq[cag_start - 50:cag_end + 50])
    dim_norms = torch.norm(emb_full, dim=0).numpy()
    dim_norms_corrected = dim_norms / np.sqrt(emb_full.shape[0])

    top_dims = np.argsort(dim_norms_corrected)[-20:][::-1]
    results["top_intensive_dims"] = top_dims.tolist()
    results["top_intensive_values"] = dim_norms_corrected[top_dims].tolist()

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(range(640), np.sort(dim_norms_corrected)[::-1], color="#2563eb", alpha=0.7)
    ax.set_xlabel("Dimension (sorted)")
    ax.set_ylabel("BN-corrected norm per dimension")
    ax.set_title("RNA-FM dimension-wise bracket norm on HTT CAG region")
    plt.tight_layout()
    plt.savefig(FIGS / "bracket_norm_per_dim.png", dpi=200)
    plt.close()
    print("  Saved figures/bracket_norm_per_dim.png")
    return results


# ── 2. TRANSPORTABILITY ──
def run_transportability(seq, model):
    print("\n=== GRASSMANNIAN TRANSPORTABILITY ===")
    cag_start, cag_end = find_cag_region(seq)
    n_wt_repeats = (cag_end - cag_start) // 3
    print(f"  WT CAG repeats: {n_wt_repeats}")

    flank_left = seq[max(0, cag_start - 50):cag_start]
    flank_right = seq[cag_end:cag_end + 50]
    cag_unit = "CAG"

    repeat_counts = [15, 21, 27, 36, 40, 50, 60, 80]
    subspace_dim = 10

    # Get WT subspace
    wt_seq = flank_left + cag_unit * n_wt_repeats + flank_right
    wt_emb, _ = get_embeddings(model, wt_seq)
    pca_wt = PCA(n_components=subspace_dim)
    pca_wt.fit(wt_emb.numpy())
    U_wt = pca_wt.components_.T  # (640, k)

    results = {"wt_repeats": n_wt_repeats, "repeat_counts": [], "geodesic_distances": [],
               "principal_angles_deg": [], "explained_var_wt": pca_wt.explained_variance_ratio_.tolist()}

    for n_rep in tqdm(repeat_counts, desc="Transportability"):
        mut_seq = flank_left + cag_unit * n_rep + flank_right
        mut_emb, _ = get_embeddings(model, mut_seq)
        pca_mut = PCA(n_components=subspace_dim)
        pca_mut.fit(mut_emb.numpy())
        U_mut = pca_mut.components_.T

        angles = subspace_angles(U_wt, U_mut)
        geodesic = np.sqrt(np.sum(angles ** 2))
        angles_deg = np.degrees(angles)

        results["repeat_counts"].append(n_rep)
        results["geodesic_distances"].append(float(geodesic))
        results["principal_angles_deg"].append(angles_deg.tolist())
        print(f"  repeats={n_rep}: geodesic={geodesic:.4f}, max_angle={angles_deg[0]:.1f}°")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].plot(results["repeat_counts"], results["geodesic_distances"], "o-",
                 color="#2563eb", linewidth=2, markersize=8)
    axes[0].axvline(x=n_wt_repeats, color="#dc2626", linestyle="--", label=f"WT ({n_wt_repeats})")
    axes[0].axvline(x=36, color="#f59e0b", linestyle="--", alpha=0.7, label="Disease threshold (~36)")
    axes[0].set_xlabel("CAG repeat count")
    axes[0].set_ylabel("Geodesic distance from WT subspace")
    axes[0].set_title("Grassmannian transportability: WT → mutant HTT")
    axes[0].legend()

    angles_matrix = np.array(results["principal_angles_deg"])
    im = axes[1].imshow(angles_matrix.T, aspect="auto", cmap="RdYlBu_r",
                        extent=[repeat_counts[0], repeat_counts[-1], subspace_dim, 1])
    axes[1].set_xlabel("CAG repeat count")
    axes[1].set_ylabel("Principal angle index")
    axes[1].set_title("Principal angles (degrees) between WT and mutant subspaces")
    plt.colorbar(im, ax=axes[1], label="Angle (°)")

    plt.tight_layout()
    plt.savefig(FIGS / "transportability.png", dpi=200)
    plt.close()
    print("  Saved figures/transportability.png")
    return results


# ── 3. LAYER-WISE PCA ──
def run_layerwise_analysis(seq, model):
    print("\n=== LAYER-WISE ANALYSIS ===")
    cag_start, cag_end = find_cag_region(seq)
    center = (cag_start + cag_end) // 2
    window = seq[max(0, center - 250):min(len(seq), center + 250)]

    _, all_hidden = get_embeddings(model, window)
    n_layers = len(all_hidden)
    print(f"  {n_layers} layers, {all_hidden[0].shape[0]} positions")

    results = {"n_layers": n_layers, "explained_variance": [], "effective_dims": []}

    fig, axes = plt.subplots(3, 5, figsize=(20, 12))
    axes = axes.flatten()

    for i, h in enumerate(tqdm(all_hidden, desc="Layer PCA")):
        h_np = h.numpy()
        pca = PCA(n_components=min(50, h_np.shape[0]))
        pca.fit(h_np)
        ev = pca.explained_variance_ratio_
        cumvar = np.cumsum(ev)
        eff_dim = np.searchsorted(cumvar, 0.95) + 1

        results["explained_variance"].append(ev[:20].tolist())
        results["effective_dims"].append(int(eff_dim))

        if i < len(axes):
            axes[i].bar(range(len(ev[:20])), ev[:20], color="#2563eb", alpha=0.7)
            axes[i].set_title(f"Layer {i} (eff_dim={eff_dim})")
            axes[i].set_ylim(0, max(ev[:3]) * 1.2)

    for j in range(len(all_hidden), len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("PCA explained variance by layer (RNA-FM on HTT CAG region)", fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGS / "layerwise_pca.png", dpi=200)
    plt.close()
    print("  Saved figures/layerwise_pca.png")

    # Cross-layer subspace similarity
    subspace_dim = 10
    subspaces = []
    for h in all_hidden:
        pca = PCA(n_components=subspace_dim)
        pca.fit(h.numpy())
        subspaces.append(pca.components_.T)

    sim_matrix = np.zeros((n_layers, n_layers))
    for i in range(n_layers):
        for j in range(n_layers):
            angles = subspace_angles(subspaces[i], subspaces[j])
            sim_matrix[i, j] = np.mean(np.cos(angles))

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(sim_matrix, cmap="RdYlBu", vmin=0, vmax=1)
    ax.set_xlabel("Layer")
    ax.set_ylabel("Layer")
    ax.set_title("Cross-layer subspace similarity (mean cos principal angle)")
    ax.set_xticks(range(n_layers))
    ax.set_yticks(range(n_layers))
    plt.colorbar(im, label="Similarity")
    plt.tight_layout()
    plt.savefig(FIGS / "cross_layer_similarity.png", dpi=200)
    plt.close()
    print("  Saved figures/cross_layer_similarity.png")

    results["cross_layer_similarity"] = sim_matrix.tolist()
    return results


# ── 4. WEIGHT FACTORIZATION (SVD) ──
def run_weight_factorization(model, config):
    print("\n=== WEIGHT FACTORIZATION (SVD) ===")
    results = {"layers": []}

    all_singular_values = []
    weight_names = []

    for layer_idx in range(config.num_hidden_layers):
        layer = model.encoder.layer[layer_idx]
        layer_result = {"layer": layer_idx, "projections": {}}

        weights = {
            "W_Q": layer.attention.self.query.weight.data,
            "W_K": layer.attention.self.key.weight.data,
            "W_V": layer.attention.self.value.weight.data,
            "W_O": layer.attention.output.dense.weight.data,
            "W_in": layer.intermediate.dense.weight.data,
            "W_out": layer.output.dense.weight.data,
        }

        for name, w in weights.items():
            u, s, v = torch.linalg.svd(w, full_matrices=False)
            s_np = s.numpy()
            cumvar = np.cumsum(s_np ** 2) / np.sum(s_np ** 2)
            eff_rank = np.searchsorted(cumvar, 0.95) + 1

            layer_result["projections"][name] = {
                "shape": list(w.shape),
                "top_20_sv": s_np[:20].tolist(),
                "effective_rank_95": int(eff_rank),
                "spectral_decay": float(s_np[0] / s_np[min(9, len(s_np) - 1)]),
            }
            all_singular_values.append(s_np[:50])
            weight_names.append(f"L{layer_idx}.{name}")

        results["layers"].append(layer_result)
        print(f"  Layer {layer_idx}: " + ", ".join(
            f"{n}={p['effective_rank_95']}" for n, p in layer_result["projections"].items()))

    # Plot: effective rank by layer and projection
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    proj_names = ["W_Q", "W_K", "W_V", "W_O", "W_in", "W_out"]
    colors = ["#2563eb", "#3b82f6", "#60a5fa", "#dc2626", "#16a34a", "#65a30d"]
    for pidx, pn in enumerate(proj_names):
        ranks = [r["projections"][pn]["effective_rank_95"] for r in results["layers"]]
        axes[0].plot(range(config.num_hidden_layers), ranks, "o-",
                     color=colors[pidx], label=pn, linewidth=2)
    axes[0].set_xlabel("Layer")
    axes[0].set_ylabel("Effective rank (95% variance)")
    axes[0].set_title("Weight matrix effective rank by layer")
    axes[0].legend(ncol=2)

    # Heatmap of spectral decay
    decay_matrix = np.zeros((config.num_hidden_layers, len(proj_names)))
    for lidx, lr in enumerate(results["layers"]):
        for pidx, pn in enumerate(proj_names):
            decay_matrix[lidx, pidx] = lr["projections"][pn]["spectral_decay"]

    im = axes[1].imshow(decay_matrix, aspect="auto", cmap="viridis")
    axes[1].set_xticks(range(len(proj_names)))
    axes[1].set_xticklabels(proj_names)
    axes[1].set_yticks(range(config.num_hidden_layers))
    axes[1].set_ylabel("Layer")
    axes[1].set_title("Spectral decay ratio (σ₁/σ₁₀)")
    plt.colorbar(im, ax=axes[1])

    plt.tight_layout()
    plt.savefig(FIGS / "weight_factorization.png", dpi=200)
    plt.close()
    print("  Saved figures/weight_factorization.png")

    # Cross-projection subspace overlap (do W_Q and W_K share singular directions?)
    print("\n  Cross-projection subspace overlap (top-k singular vectors)...")
    k = 20
    overlap_results = {}
    for layer_idx in range(config.num_hidden_layers):
        layer = model.encoder.layer[layer_idx]
        weights = {
            "W_Q": layer.attention.self.query.weight.data,
            "W_K": layer.attention.self.key.weight.data,
            "W_V": layer.attention.self.value.weight.data,
            "W_O": layer.attention.output.dense.weight.data,
        }
        svd_results = {}
        for name, w in weights.items():
            u, s, v = torch.linalg.svd(w, full_matrices=False)
            svd_results[name] = v[:k].T.numpy()  # (640, k) right singular vectors

        pairs = [("W_Q", "W_K"), ("W_Q", "W_V"), ("W_K", "W_V"), ("W_Q", "W_O")]
        for a, b in pairs:
            angles = subspace_angles(svd_results[a], svd_results[b])
            mean_cos = float(np.mean(np.cos(angles)))
            key = f"L{layer_idx}.{a}-{b}"
            overlap_results[key] = mean_cos

    results["cross_projection_overlap"] = overlap_results
    print("  Done.")
    return results


# ── 5. POSITION-WISE ANALYSIS ──
def run_position_analysis(seq, model):
    print("\n=== POSITION-WISE ANALYSIS ===")
    cag_start, cag_end = find_cag_region(seq)
    flank = 100
    start = max(0, cag_start - flank)
    end = min(len(seq), cag_end + flank)
    window = seq[start:end]

    emb, _ = get_embeddings(model, window)
    emb_np = emb.numpy()

    cag_local_start = cag_start - start
    cag_local_end = cag_end - start

    # PCA colored by region
    pca = PCA(n_components=3)
    coords = pca.fit_transform(emb_np)

    regions = np.zeros(emb_np.shape[0])
    regions[cag_local_start:cag_local_end] = 1  # CAG = 1, flanking = 0

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for r, label, color in [(0, "Flanking", "#2563eb"), (1, "CAG repeat", "#dc2626")]:
        mask = regions == r
        axes[0].scatter(coords[mask, 0], coords[mask, 1], c=color, label=label, alpha=0.6, s=15)
    axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    axes[0].set_title("RNA-FM embedding PCA: CAG vs flanking positions")
    axes[0].legend()

    # Norm profile along sequence
    norms = np.linalg.norm(emb_np, axis=1)
    axes[1].plot(range(len(norms)), norms, color="#2563eb", linewidth=0.8, alpha=0.7)
    axes[1].axvspan(cag_local_start, cag_local_end, color="#dc2626", alpha=0.15, label="CAG region")
    axes[1].set_xlabel("Position (nt)")
    axes[1].set_ylabel("Embedding norm")
    axes[1].set_title("Position-wise embedding norm")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(FIGS / "position_analysis.png", dpi=200)
    plt.close()
    print("  Saved figures/position_analysis.png")

    return {
        "pca_explained_variance": pca.explained_variance_ratio_.tolist(),
        "cag_region": [int(cag_local_start), int(cag_local_end)],
        "mean_norm_cag": float(norms[cag_local_start:cag_local_end].mean()),
        "mean_norm_flank": float(np.concatenate([norms[:cag_local_start], norms[cag_local_end:]]).mean()),
    }


def main():
    print(f"RNA-FM analysis pipeline — {TIMESTAMP}")
    print("=" * 60)

    print("Loading model...")
    model, config = load_model()
    seq = get_full_sequence()
    print(f"HTT mRNA: {len(seq)} nt")

    all_results = {"timestamp": TIMESTAMP, "sequence_length": len(seq)}

    all_results["bracket_norm"] = run_bracket_norm(seq, model)
    all_results["transportability"] = run_transportability(seq, model)
    all_results["layerwise"] = run_layerwise_analysis(seq, model)
    all_results["weight_factorization"] = run_weight_factorization(model, config)
    all_results["position_analysis"] = run_position_analysis(seq, model)

    # Save all results
    out_path = DATA / f"analysis_results_{TIMESTAMP}.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nAll results saved to {out_path}")

    # Also save a latest symlink
    latest = DATA / "analysis_results_latest.json"
    if latest.exists():
        latest.unlink()
    os.symlink(out_path.name, latest)
    print(f"Symlinked {latest} → {out_path.name}")

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    bn = all_results["bracket_norm"]
    print(f"Bracket norm: raw norms scale {bn['norms_raw'][0]:.2f}→{bn['norms_raw'][-1]:.2f}, "
          f"corrected {bn['norms_corrected'][0]:.4f}→{bn['norms_corrected'][-1]:.4f}")

    tr = all_results["transportability"]
    print(f"Transportability: geodesic distances {[f'{d:.3f}' for d in tr['geodesic_distances']]}")

    lw = all_results["layerwise"]
    print(f"Effective dims by layer: {lw['effective_dims']}")

    pa = all_results["position_analysis"]
    print(f"CAG norm={pa['mean_norm_cag']:.3f}, flanking norm={pa['mean_norm_flank']:.3f}")

    print("\nDONE.")


if __name__ == "__main__":
    main()
