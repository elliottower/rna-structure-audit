"""Multi-model geometric audit of RNA foundation model representations.

Runs the same analysis pipeline from Paper 1 (main_v3.tex) on multiple RNA/DNA
models: RNA-FM, Nucleotide Transformer v2, HyenaDNA, and optionally Evo.

Analyses (all model-agnostic once embeddings are extracted):
  1. Bracket-norm correction
  2. Grassmannian transportability
  3. Layer-wise compression (effective dimensionality per layer)
  4. Direction instability (DI) vs wild-type
  5. Cross-view ViennaRNA independence
  6. Preflight-bio composite score (M1-M4)

Usage:
    uv run python multi_model_audit.py --models rnafm nt hyenadna
    uv run python multi_model_audit.py --models rnafm --skip-preflight
"""

import argparse
import json
import subprocess
import sys
import warnings
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path

import numpy as np
import torch
from scipy import stats
from scipy.linalg import subspace_angles
from sklearn.decomposition import PCA
from tqdm import tqdm

REPEAT_COUNTS = [10, 15, 18, 21, 24, 27, 30, 36, 40, 50, 60, 80]
WT_REPEATS = 21
FLANK_SIZE = 50
SEED = 20260711

HTT_LEFT_FLANK = (
    "ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC"
)
HTT_RIGHT_FLANK = (
    "CAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAG"
)


# ── Model adapter base ──────────────────────────────────────────────────────

class ModelAdapter(ABC):
    name: str
    d_model: int
    n_layers: int
    token_resolution: str  # "nucleotide", "6mer", "byte"

    @abstractmethod
    def load(self) -> None:
        ...

    @abstractmethod
    def tokenize(self, seq: str) -> torch.Tensor:
        ...

    @abstractmethod
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        """Return list of (seq_len, d_model) tensors, one per layer (including embedding layer)."""
        ...

    def get_final_embeddings(self, tokens: torch.Tensor) -> torch.Tensor:
        return self.get_all_layer_embeddings(tokens)[-1]


# ── RNA-FM adapter ───────────────────────────────────────────────────────────

class RNAFMAdapter(ModelAdapter):
    name = "RNA-FM"
    d_model = 640
    n_layers = 12
    token_resolution = "nucleotide"

    NUC_TO_ID = {"A": 5, "C": 6, "G": 7, "U": 8}

    def __init__(self):
        self.model = None

    def load(self):
        from transformers import BertConfig, BertModel

        weight_path = Path("pretrained/pytorch_model.bin")
        if not weight_path.exists():
            raise FileNotFoundError(f"RNA-FM weights not found at {weight_path}")

        state = torch.load(weight_path, map_location="cpu", weights_only=False)
        max_pos = state["model.embeddings.position_embeddings.weight"].shape[0]

        config = BertConfig(
            vocab_size=28, hidden_size=640, num_hidden_layers=12,
            num_attention_heads=20, intermediate_size=5120,
            max_position_embeddings=max_pos,
        )
        self.model = BertModel(config)

        clean_state = {}
        for k, v in state.items():
            if not k.startswith("model."):
                continue
            nk = k[len("model."):]
            nk = nk.replace("attention.layer_norm", "attention.output.LayerNorm")
            if "encoder.layer." in nk and ".layer_norm." in nk and "attention" not in nk:
                nk = nk.replace("layer_norm", "output.LayerNorm")
            elif nk.startswith("embeddings.layer_norm"):
                nk = nk.replace("embeddings.layer_norm", "embeddings.LayerNorm")
            clean_state[nk] = v

        self.model.load_state_dict(clean_state, strict=False)
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        tokens = [2] + [self.NUC_TO_ID.get(c, 1) for c in seq] + [3]
        return torch.tensor([tokens])

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        return [hs[0, 1:-1, :] for hs in out.hidden_states]


# ── Nucleotide Transformer v2 adapter ────────────────────────────────────────

class NTAdapter(ModelAdapter):
    name = "Nucleotide-Transformer-v2"
    d_model = 512
    n_layers = 12
    token_resolution = "6mer"

    MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model = AutoModelForMaskedLM.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        seq_dna = seq.replace("U", "T")
        encoded = self.tokenizer(seq_dna, return_tensors="pt")
        return encoded["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        # Strip CLS/EOS, expand 6-mer tokens to per-nucleotide by repeating
        layers = []
        for hs in out.hidden_states:
            emb = hs[0, 1:-1, :]  # strip special tokens
            layers.append(emb)
        return layers


# ── HyenaDNA adapter ────────────────────────────────────────────────────────

class HyenaDNAAdapter(ModelAdapter):
    name = "HyenaDNA-tiny"
    d_model = 128
    n_layers = 2
    token_resolution = "nucleotide"

    MODEL_ID = "LongSafari/hyenadna-tiny-1k-seqlen-hf"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from transformers import AutoModelForCausalLM, AutoTokenizer

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model = AutoModelForCausalLM.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        seq_dna = seq.replace("U", "T")
        encoded = self.tokenizer(seq_dna, return_tensors="pt")
        return encoded["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        if hasattr(out, "hidden_states") and out.hidden_states:
            return [hs[0] for hs in out.hidden_states]
        # Fallback: only final layer
        return [out.logits[0]]


# ── Evo adapter (GPU only, placeholder) ─────────────────────────────────────

class EvoAdapter(ModelAdapter):
    name = "Evo-1-8k"
    d_model = 4096
    n_layers = 32
    token_resolution = "byte"

    MODEL_ID = "togethercomputer/evo-1-8k-base"

    def __init__(self):
        self.model = None

    def load(self):
        if not torch.cuda.is_available():
            raise RuntimeError("Evo (7B params) requires CUDA GPU")

        from transformers import AutoModelForCausalLM

        self.model = AutoModelForCausalLM.from_pretrained(
            self.MODEL_ID, trust_remote_code=True, torch_dtype=torch.float16,
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        seq_dna = seq.replace("U", "T")
        tokens = [ord(c) for c in seq_dna]
        return torch.tensor([tokens])

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        tokens = tokens.to(next(self.model.parameters()).device)
        out = self.model(tokens, output_hidden_states=True)
        return [hs[0].float().cpu() for hs in out.hidden_states]


# ── Caduceus adapter ───────────────────────────────────────────────────────

class CaduceusAdapter(ModelAdapter):
    name = "Caduceus"
    d_model = 256
    n_layers = 16
    token_resolution = "nucleotide"

    MODEL_ID = "kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model = AutoModelForMaskedLM.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        seq_dna = seq.replace("U", "T")
        encoded = self.tokenizer(seq_dna, return_tensors="pt")
        return encoded["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        if hasattr(out, "hidden_states") and out.hidden_states:
            return [hs[0] for hs in out.hidden_states]
        return [out.logits[0]]


# ── RiNALMo adapter ──────────────────────────────────────────────────────────

class RiNALMoAdapter(ModelAdapter):
    name = "RiNALMo"
    d_model = 1280
    n_layers = 33
    token_resolution = "nucleotide"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from multimolecule import RnaTokenizer, RiNALMoModel

        self.tokenizer = RnaTokenizer.from_pretrained("multimolecule/rinalmo-giga")
        self.model = RiNALMoModel.from_pretrained(
            "multimolecule/rinalmo-giga", attn_implementation="eager",
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        enc = self.tokenizer(seq, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        return [hs[0, 1:-1, :] for hs in out.hidden_states]


# ── UTR-LM adapter ───────────────────────────────────────────────────────────

class UTRLMAdapter(ModelAdapter):
    name = "UTR-LM"
    d_model = 128
    n_layers = 6
    token_resolution = "nucleotide"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from multimolecule import RnaTokenizer, UtrLmModel

        self.tokenizer = RnaTokenizer.from_pretrained("multimolecule/utrlm-te_el")
        self.model = UtrLmModel.from_pretrained(
            "multimolecule/utrlm-te_el", attn_implementation="eager",
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        enc = self.tokenizer(seq, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        return [hs[0, 1:-1, :] for hs in out.hidden_states]


# ── ERNIE-RNA adapter ────────────────────────────────────────────────────────

class ERNIERNAAdapter(ModelAdapter):
    name = "ERNIE-RNA"
    d_model = 768
    n_layers = 12
    token_resolution = "nucleotide"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from multimolecule import RnaTokenizer, ErnieRnaModel

        self.tokenizer = RnaTokenizer.from_pretrained("multimolecule/ernierna")
        self.model = ErnieRnaModel.from_pretrained(
            "multimolecule/ernierna", attn_implementation="eager",
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        enc = self.tokenizer(seq, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        return [hs[0, 1:-1, :] for hs in out.hidden_states]


# ── SpliceBERT adapter ───────────────────────────────────────────────────────

class SpliceBERTAdapter(ModelAdapter):
    name = "SpliceBERT"
    d_model = 512
    n_layers = 6
    token_resolution = "nucleotide"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from multimolecule import RnaTokenizer, SpliceBertModel

        self.tokenizer = RnaTokenizer.from_pretrained("multimolecule/splicebert")
        self.model = SpliceBertModel.from_pretrained(
            "multimolecule/splicebert", attn_implementation="eager",
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        enc = self.tokenizer(seq, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens, output_hidden_states=True)
        return [hs[0, 1:-1, :] for hs in out.hidden_states]


# ── DNABERT-2 adapter ────────────────────────────────────────────────────────

class DNABERT2Adapter(ModelAdapter):
    name = "DNABERT-2"
    d_model = 768
    n_layers = 12
    token_resolution = "bpe"

    MODEL_ID = "zhihan1996/DNABERT-2-117M"

    def __init__(self):
        self.model = None
        self.tokenizer = None

    def load(self):
        from transformers import AutoModel, AutoTokenizer

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model = AutoModel.from_pretrained(
            self.MODEL_ID, trust_remote_code=True,
        )
        self.model.eval()

    def tokenize(self, seq: str) -> torch.Tensor:
        seq_dna = seq.replace("U", "T")
        enc = self.tokenizer(seq_dna, return_tensors="pt")
        return enc["input_ids"]

    @torch.no_grad()
    def get_all_layer_embeddings(self, tokens: torch.Tensor) -> list[torch.Tensor]:
        out = self.model(tokens)
        last_hidden = out[0] if isinstance(out, tuple) else out
        if last_hidden.dim() == 2:
            last_hidden = last_hidden.unsqueeze(0)
        return [last_hidden[0, 1:-1, :]]


# ── HTT mRNA variant construction ───────────────────────────────────────────

def get_htt_variant(repeat_count: int) -> str:
    cag = "CAG" * repeat_count
    return HTT_LEFT_FLANK + cag + HTT_RIGHT_FLANK


def build_aligned_embeddings(
    adapter: ModelAdapter, layer_idx: int, repeat_counts: list[int],
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Build position-aligned embedding matrix for DI analysis.

    Returns:
        wt_embs: (n_positions, d_model) for WT
        variant_embs: dict mapping r -> (n_positions, d_model)
        position_labels: list of position labels
    """
    min_codons = min(repeat_counts) // 3
    n_codon_slots = min(10, min_codons)
    n_positions = FLANK_SIZE + n_codon_slots + FLANK_SIZE

    variant_embs = {}
    for r in tqdm(repeat_counts, desc=f"  Extracting {adapter.name} L{layer_idx}"):
        seq = get_htt_variant(r)
        tokens = adapter.tokenize(seq)
        all_layers = adapter.get_all_layer_embeddings(tokens)

        if layer_idx >= len(all_layers):
            layer_idx = len(all_layers) - 1
        emb = all_layers[layer_idx].numpy()

        # For 6-mer tokenizers, positions don't align 1:1 with nucleotides.
        # Use whatever positions are available and truncate to aligned size.
        n_avail = emb.shape[0]

        # Align: left flank | mean-pooled codons | right flank
        left = emb[:FLANK_SIZE] if n_avail >= FLANK_SIZE else emb[:n_avail]
        if left.shape[0] < FLANK_SIZE:
            left = np.pad(left, ((0, FLANK_SIZE - left.shape[0]), (0, 0)))

        cag_start = FLANK_SIZE
        cag_end = cag_start + r * 3
        if adapter.token_resolution == "6mer":
            # For 6-mer models, tokens cover multiple nt.
            # Use available token positions in the CAG region.
            n_tok = n_avail
            tok_per_nt = n_tok / (FLANK_SIZE * 2 + r * 3)
            cag_tok_start = int(FLANK_SIZE * tok_per_nt)
            cag_tok_end = int(cag_end * tok_per_nt)
            cag_embs = emb[cag_tok_start:cag_tok_end]
        else:
            cag_embs = emb[cag_start:min(cag_end, n_avail)]

        # Mean-pool CAG codons into n_codon_slots
        if cag_embs.shape[0] >= n_codon_slots:
            codon_size = cag_embs.shape[0] // n_codon_slots
            pooled = np.stack([
                cag_embs[i * codon_size:(i + 1) * codon_size].mean(axis=0)
                for i in range(n_codon_slots)
            ])
        else:
            pooled = np.zeros((n_codon_slots, emb.shape[1]))
            pooled[:cag_embs.shape[0]] = cag_embs

        right_start = cag_end if adapter.token_resolution != "6mer" else int(cag_end * tok_per_nt) if adapter.token_resolution == "6mer" else cag_end
        right = emb[right_start:right_start + FLANK_SIZE] if n_avail > right_start else np.zeros((FLANK_SIZE, emb.shape[1]))
        if right.shape[0] < FLANK_SIZE:
            right = np.pad(right, ((0, FLANK_SIZE - right.shape[0]), (0, 0)))

        aligned = np.concatenate([left, pooled, right], axis=0)
        variant_embs[r] = aligned

    wt_embs = variant_embs[WT_REPEATS]
    return wt_embs, variant_embs, (
        [f"L{i}" for i in range(FLANK_SIZE)]
        + [f"CAG{i}" for i in range(n_codon_slots)]
        + [f"R{i}" for i in range(FLANK_SIZE)]
    )


# ── Analysis functions (model-agnostic) ──────────────────────────────────────

def compute_di_vs_wt(wt_emb: np.ndarray, var_emb: np.ndarray) -> np.ndarray:
    """DI = 1 - cosine_similarity per position."""
    wt_norm = wt_emb / (np.linalg.norm(wt_emb, axis=1, keepdims=True) + 1e-10)
    var_norm = var_emb / (np.linalg.norm(var_emb, axis=1, keepdims=True) + 1e-10)
    cos = np.sum(wt_norm * var_norm, axis=1)
    return 1.0 - cos


def analysis_bracket_norm(
    adapter: ModelAdapter, window_sizes: list[int] = [100, 200, 300, 400, 500],
) -> dict:
    seq = get_htt_variant(WT_REPEATS)
    tokens = adapter.tokenize(seq)
    emb = adapter.get_final_embeddings(tokens).numpy()
    n_total = emb.shape[0]

    center = n_total // 2
    results = {}
    for w in window_sizes:
        half = w // 2
        start = max(0, center - half)
        end = min(n_total, center + half)
        window = emb[start:end]
        n = window.shape[0]
        mean_norm = float(np.mean(np.linalg.norm(window, axis=1)))
        bn_corrected = mean_norm / np.sqrt(n)
        results[w] = {"mean_norm": mean_norm, "bn_corrected": bn_corrected, "n": n}

    return results


def analysis_grassmannian_transportability(
    adapter: ModelAdapter, k: int = 10, layer_idx: int = -1,
) -> dict:
    subspaces = {}
    for r in tqdm(REPEAT_COUNTS, desc=f"  Grassmannian {adapter.name}"):
        seq = get_htt_variant(r)
        tokens = adapter.tokenize(seq)
        all_layers = adapter.get_all_layer_embeddings(tokens)
        emb = all_layers[layer_idx].numpy()
        n_components = min(k, emb.shape[0], emb.shape[1])
        pca = PCA(n_components=n_components)
        pca.fit(emb)
        subspaces[r] = pca.components_.T  # (d_model, k)

    wt_sub = subspaces[WT_REPEATS]
    distances = {}
    for r, sub in subspaces.items():
        if r == WT_REPEATS:
            continue
        k_shared = min(wt_sub.shape[1], sub.shape[1])
        angles = subspace_angles(wt_sub[:, :k_shared], sub[:, :k_shared])
        distances[r] = {
            "geodesic": float(np.sqrt(np.sum(angles ** 2))),
            "max_angle_deg": float(np.degrees(np.max(angles))),
        }

    return distances


def analysis_layer_compression(adapter: ModelAdapter) -> dict:
    seq = get_htt_variant(WT_REPEATS)
    tokens = adapter.tokenize(seq)
    all_layers = adapter.get_all_layer_embeddings(tokens)

    results = {}
    for i, emb in enumerate(all_layers):
        emb_np = emb.numpy()
        n_components = min(emb_np.shape[0], emb_np.shape[1])
        pca = PCA(n_components=n_components)
        pca.fit(emb_np)
        cumvar = np.cumsum(pca.explained_variance_ratio_)
        eff_dim = int(np.searchsorted(cumvar, 0.95) + 1)
        pc1_var = float(pca.explained_variance_ratio_[0])
        results[i] = {
            "effective_dim_95": eff_dim,
            "pc1_variance": pc1_var,
            "n_components": n_components,
        }

    return results


def analysis_direction_instability(
    adapter: ModelAdapter, layer_idx: int,
) -> dict:
    wt_embs, variant_embs, pos_labels = build_aligned_embeddings(
        adapter, layer_idx, REPEAT_COUNTS,
    )

    # Per-variant DI
    di_per_variant = {}
    for r, var_emb in variant_embs.items():
        if r == WT_REPEATS:
            continue
        di = compute_di_vs_wt(wt_embs, var_emb)
        di_per_variant[r] = float(np.mean(di))

    # Mean DI profile across all non-WT variants
    di_profiles = []
    for r, var_emb in variant_embs.items():
        if r == WT_REPEATS:
            continue
        di_profiles.append(compute_di_vs_wt(wt_embs, var_emb))
    mean_profile = np.mean(di_profiles, axis=0)

    # Region-wise DI
    left_di = float(np.mean(mean_profile[:FLANK_SIZE]))
    cag_di = float(np.mean(mean_profile[FLANK_SIZE:FLANK_SIZE + 10]))
    right_di = float(np.mean(mean_profile[FLANK_SIZE + 10:]))
    cag_elevation = cag_di / left_di if left_di > 0 else float("inf")

    # DI vs repeat deviation correlation
    deviations = []
    di_scores = []
    for r in sorted(di_per_variant.keys()):
        deviations.append(abs(r - WT_REPEATS))
        di_scores.append(di_per_variant[r])
    rho, p = stats.spearmanr(deviations, di_scores)

    # DI uniformity (how identical are values across positions)
    di_range = float(np.max(mean_profile) - np.min(mean_profile))

    return {
        "di_per_variant": {str(r): v for r, v in di_per_variant.items()},
        "left_flank_di": left_di,
        "cag_di": cag_di,
        "right_flank_di": right_di,
        "cag_elevation_vs_left": cag_elevation,
        "di_range": di_range,
        "rho_di_vs_deviation": float(rho),
        "rho_p": float(p),
    }


def analysis_cross_view_vienna(adapter: ModelAdapter, layer_idx: int) -> dict:
    try:
        import RNA
    except ImportError:
        return {"error": "ViennaRNA not installed (pip install ViennaRNA)"}

    seq = get_htt_variant(WT_REPEATS)
    seq_dna = seq.replace("U", "T") if "U" not in seq else seq

    # ViennaRNA base-pair probabilities
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe()
    fc.pf()
    bpp = fc.bpp()

    n = len(seq)
    bp_prob = np.zeros(n)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if i != j:
                bp_prob[i - 1] += bpp[i][j] if i < j else bpp[j][i]

    # RNA-FM embeddings
    tokens = adapter.tokenize(seq)
    all_layers = adapter.get_all_layer_embeddings(tokens)
    emb = all_layers[layer_idx].numpy()

    # Align lengths (RNA-FM strips CLS/EOS, but bp_prob is full length)
    min_len = min(emb.shape[0], len(bp_prob))
    emb_aligned = emb[:min_len]
    bp_aligned = bp_prob[:min_len]

    # PC-level correlation
    n_pcs = min(10, emb_aligned.shape[0], emb_aligned.shape[1])
    pca = PCA(n_components=n_pcs)
    scores = pca.fit_transform(emb_aligned)

    best_r = 0.0
    best_p = 1.0
    for pc_idx in range(n_pcs):
        r, p = stats.pearsonr(scores[:, pc_idx], bp_aligned)
        if abs(r) > abs(best_r):
            best_r = float(r)
            best_p = float(p)

    # Subspace geodesic
    bp_2d = bp_aligned.reshape(-1, 1)
    from scipy.linalg import qr
    emb_q, _ = qr(scores, mode="economic")
    bp_q, _ = qr(bp_2d, mode="economic")
    k_shared = min(emb_q.shape[1], bp_q.shape[1])
    angles = subspace_angles(emb_q[:, :k_shared], bp_q[:, :k_shared])
    geodesic = float(np.sqrt(np.sum(angles ** 2)))

    return {
        "best_pc_bp_correlation": best_r,
        "best_pc_bp_p": best_p,
        "subspace_geodesic": geodesic,
        "mfe": float(mfe) if isinstance(mfe, (int, float)) else mfe,
    }


def analysis_preflight(adapter: ModelAdapter, layer_idx: int) -> dict:
    """Run preflight-bio M1-M4 on WT vs each variant."""
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent / "preflight-bio" / "src"))
        from preflight.core.runner import run as preflight_run
    except ImportError:
        return {"error": "preflight-bio not importable"}

    # Extract WT embeddings
    wt_seq = get_htt_variant(WT_REPEATS)
    wt_tokens = adapter.tokenize(wt_seq)
    wt_emb = adapter.get_all_layer_embeddings(wt_tokens)[layer_idx].numpy()

    results = {}
    for r in tqdm(REPEAT_COUNTS, desc=f"  Preflight {adapter.name}"):
        if r == WT_REPEATS:
            continue
        var_seq = get_htt_variant(r)
        var_tokens = adapter.tokenize(var_seq)
        var_emb = adapter.get_all_layer_embeddings(var_tokens)[layer_idx].numpy()

        try:
            report = preflight_run(
                source=wt_emb,
                target=var_emb,
                model_name=f"{adapter.name}_r{r}",
                modules=["m1_grassmannian", "m2_domain_shift",
                         "m3_direction_stability", "m4_domain_validity"],
            )
            results[r] = {
                "overall_score": float(report.overall_score),
                "overall_tier": report.overall_tier,
                "tier_label": report.overall_tier_label,
                "modules": {
                    name: {"score": float(m.score), "flag": m.flag}
                    for name, m in report.modules.items()
                },
            }
        except Exception as e:
            results[r] = {"error": str(e)}

    return results


# ── Main pipeline ────────────────────────────────────────────────────────────

MODEL_REGISTRY = {
    "rnafm": RNAFMAdapter,
    "nt": NTAdapter,
    "hyenadna": HyenaDNAAdapter,
    "evo": EvoAdapter,
}


def run_audit(adapter: ModelAdapter, skip_preflight: bool = False, skip_vienna: bool = False) -> dict:
    print(f"\n{'='*60}")
    print(f"  {adapter.name}  (d={adapter.d_model}, L={adapter.n_layers}, res={adapter.token_resolution})")
    print(f"{'='*60}")

    print(f"\n[{datetime.now():%H:%M:%S}] Loading model...")
    adapter.load()

    results = {
        "model": adapter.name,
        "d_model": adapter.d_model,
        "n_layers": adapter.n_layers,
        "token_resolution": adapter.token_resolution,
        "timestamp": datetime.now().isoformat(),
    }

    # 1. Bracket norm
    print(f"[{datetime.now():%H:%M:%S}] Bracket-norm analysis...")
    results["bracket_norm"] = analysis_bracket_norm(adapter)

    # 2. Grassmannian transportability (final layer)
    print(f"[{datetime.now():%H:%M:%S}] Grassmannian transportability (final layer)...")
    results["grassmannian_final"] = analysis_grassmannian_transportability(adapter, layer_idx=-1)

    # 3. Layer compression
    print(f"[{datetime.now():%H:%M:%S}] Layer-wise compression...")
    results["layer_compression"] = analysis_layer_compression(adapter)

    # Find best analysis layer (highest effective dim that isn't the embedding layer)
    compressions = results["layer_compression"]
    best_layer = 0
    best_dim = 0
    for layer_str, info in compressions.items():
        layer_i = int(layer_str)
        if layer_i == 0:
            continue
        if info["effective_dim_95"] > best_dim:
            best_dim = info["effective_dim_95"]
            best_layer = layer_i

    results["best_analysis_layer"] = best_layer
    results["best_analysis_layer_dim"] = best_dim

    # 4. Direction instability (final layer + best layer)
    print(f"[{datetime.now():%H:%M:%S}] Direction instability (final layer)...")
    results["di_final"] = analysis_direction_instability(adapter, layer_idx=-1)

    if best_layer != adapter.n_layers:
        print(f"[{datetime.now():%H:%M:%S}] Direction instability (layer {best_layer}, eff_dim={best_dim})...")
        results["di_best_layer"] = analysis_direction_instability(adapter, layer_idx=best_layer)
        results["di_best_layer_idx"] = best_layer

    # 5. Cross-view ViennaRNA
    if not skip_vienna:
        print(f"[{datetime.now():%H:%M:%S}] Cross-view ViennaRNA (final layer)...")
        results["cross_view_final"] = analysis_cross_view_vienna(adapter, layer_idx=-1)

        if best_layer != adapter.n_layers:
            print(f"[{datetime.now():%H:%M:%S}] Cross-view ViennaRNA (layer {best_layer})...")
            results["cross_view_best"] = analysis_cross_view_vienna(adapter, layer_idx=best_layer)

    # 6. Preflight-bio
    if not skip_preflight:
        print(f"[{datetime.now():%H:%M:%S}] Preflight-bio composite (final layer)...")
        results["preflight_final"] = analysis_preflight(adapter, layer_idx=-1)

        if best_layer != adapter.n_layers:
            print(f"[{datetime.now():%H:%M:%S}] Preflight-bio composite (layer {best_layer})...")
            results["preflight_best"] = analysis_preflight(adapter, layer_idx=best_layer)

    print(f"[{datetime.now():%H:%M:%S}] Done with {adapter.name}.")
    return results


def print_comparison_table(all_results: list[dict]) -> None:
    print(f"\n{'='*80}")
    print("  CROSS-MODEL COMPARISON TABLE")
    print(f"{'='*80}\n")

    headers = ["Metric"] + [r["model"] for r in all_results]
    rows = []

    # Basic info
    rows.append(["d_model"] + [str(r["d_model"]) for r in all_results])
    rows.append(["n_layers"] + [str(r["n_layers"]) for r in all_results])
    rows.append(["token_resolution"] + [r["token_resolution"] for r in all_results])

    # Layer compression
    rows.append(["---"] * len(headers))
    for r in all_results:
        comp = r.get("layer_compression", {})
        if comp:
            last_key = str(max(int(k) for k in comp.keys()))
            first_key = "0"
            rows_data = [
                f"L0={comp.get(first_key, {}).get('effective_dim_95', '?')}",
            ]
    rows.append(["Final layer eff_dim"] + [
        str(r.get("layer_compression", {}).get(
            str(max(int(k) for k in r["layer_compression"].keys())), {}
        ).get("effective_dim_95", "?"))
        if r.get("layer_compression") else "?"
        for r in all_results
    ])
    rows.append(["Best analysis layer"] + [
        f"L{r.get('best_analysis_layer', '?')} (dim={r.get('best_analysis_layer_dim', '?')})"
        for r in all_results
    ])

    # Grassmannian
    rows.append(["---"] * len(headers))
    rows.append(["Grassmannian (mean d_g)"] + [
        f"{np.mean([v['geodesic'] for v in r.get('grassmannian_final', {}).values()]):.2f}"
        if r.get("grassmannian_final") else "?"
        for r in all_results
    ])

    # DI
    rows.append(["---"] * len(headers))
    rows.append(["DI final: CAG elevation"] + [
        f"{r.get('di_final', {}).get('cag_elevation_vs_left', '?'):.2f}x"
        if isinstance(r.get("di_final", {}).get("cag_elevation_vs_left"), (int, float)) else "?"
        for r in all_results
    ])
    rows.append(["DI final: range"] + [
        f"{r.get('di_final', {}).get('di_range', '?'):.6f}"
        if isinstance(r.get("di_final", {}).get("di_range"), (int, float)) else "?"
        for r in all_results
    ])
    rows.append(["DI final: rho(dev)"] + [
        f"{r.get('di_final', {}).get('rho_di_vs_deviation', '?'):.3f}"
        if isinstance(r.get("di_final", {}).get("rho_di_vs_deviation"), (int, float)) else "?"
        for r in all_results
    ])

    if any(r.get("di_best_layer") for r in all_results):
        rows.append(["DI best: CAG elevation"] + [
            f"{r.get('di_best_layer', {}).get('cag_elevation_vs_left', '?'):.2f}x"
            if isinstance(r.get("di_best_layer", {}).get("cag_elevation_vs_left"), (int, float)) else "?"
            for r in all_results
        ])

    # Cross-view
    rows.append(["---"] * len(headers))
    rows.append(["Vienna |r| (final)"] + [
        f"{abs(r.get('cross_view_final', {}).get('best_pc_bp_correlation', 0)):.3f}"
        if isinstance(r.get("cross_view_final", {}).get("best_pc_bp_correlation"), (int, float)) else "?"
        for r in all_results
    ])

    # Print table
    col_widths = [max(len(str(row[i])) for row in [headers] + rows) for i in range(len(headers))]
    fmt = " | ".join(f"{{:<{w}}}" for w in col_widths)
    print(fmt.format(*headers))
    print("-+-".join("-" * w for w in col_widths))
    for row in rows:
        if row[0] == "---":
            print("-+-".join("-" * w for w in col_widths))
        else:
            print(fmt.format(*[str(x) for x in row]))


def main():
    parser = argparse.ArgumentParser(description="Multi-model RNA geometric audit")
    parser.add_argument(
        "--models", nargs="+", default=["rnafm"],
        choices=list(MODEL_REGISTRY.keys()),
        help="Models to audit",
    )
    parser.add_argument("--skip-preflight", action="store_true")
    parser.add_argument("--skip-vienna", action="store_true")
    parser.add_argument("--output-dir", default="data/multi_model", type=str)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for model_key in args.models:
        adapter_cls = MODEL_REGISTRY[model_key]
        adapter = adapter_cls()

        try:
            results = run_audit(
                adapter,
                skip_preflight=args.skip_preflight,
                skip_vienna=args.skip_vienna,
            )
            all_results.append(results)

            out_path = output_dir / f"audit_{model_key}_{datetime.now():%Y%m%d_%H%M%S}.json"
            with open(out_path, "w") as f:
                json.dump(results, f, indent=2, default=str)
            print(f"\nSaved: {out_path}")

        except Exception as e:
            print(f"\nFailed on {adapter.name}: {e}")
            import traceback
            traceback.print_exc()

    if len(all_results) > 1:
        print_comparison_table(all_results)

    # Save combined results
    combined_path = output_dir / f"audit_combined_{datetime.now():%Y%m%d_%H%M%S}.json"
    with open(combined_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nCombined results: {combined_path}")


if __name__ == "__main__":
    main()
