"""Modal wrapper: HTT trinucleotide repeat case study on all RNA/DNA foundation models.

Measures how embedding distance scales with CAG repeat length in the
Huntington's disease gene (HTT exon 1). For each model, constructs HTT mRNA
variants with different CAG repeat counts, extracts embeddings at all layers,
mean-pools across positions, and computes cosine distance from baseline (N=10).

Picks the layer with strongest Spearman correlation between repeat count and
cosine distance. Also runs ViennaRNA MFE predictions (CPU-only).

Models (9 total, RNA-FM excluded -- needs fm package):
  multimol image:      ERNIE-RNA, SpliceBERT, RiNALMo, UTR-LM
  dnabert2 image:      DNABERT-2
  transformers image:  Nucleotide Transformer v2, HyenaDNA, Evo, Caduceus

Usage:
    cd /path/to/causal-rna
    modal run --detach scripts/modal_htt_all_models.py
"""

import modal

app = modal.App("causal-rna-htt-all-models")

# -- Constants ---------------------------------------------------------------

REPEAT_COUNTS = [10, 15, 18, 21, 24, 27, 30, 36, 40, 50, 60, 80]

HTT_LEFT_FLANK = "CCAGCAGCAGCAGCAGCAG"
HTT_RIGHT_FLANK = "CCGCCGCCGCCGCCGCCG"

RNA_MODELS = {"ernierna", "splicebert", "rinalmo", "utrlm"}

MULTIMOL_MODELS = ["ernierna", "splicebert", "rinalmo", "utrlm"]
TRANSFORMERS_MODELS = ["nt", "hyenadna", "evo", "caduceus"]

# -- Images ------------------------------------------------------------------

multimol_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.6.0",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.66.5",
        "transformers==5.14.1",
        "multimolecule==0.2.0",
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
)

dnabert2_image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "torch==2.4.0",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.66.5",
        "transformers==4.28.0",
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
    )
    .add_local_file(
        "scripts/patch_dnabert2_flash_attn.py",
        "/root/patch_dnabert2_flash_attn.py",
        copy=True,
    )
    .run_commands("python /root/patch_dnabert2_flash_attn.py")
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
)

# NOTE: mamba-ssm and causal-conv1d are needed for Caduceus. If image build
# fails (missing CUDA headers for compilation), remove them -- Caduceus will
# then fail gracefully at runtime with an ImportError.
transformers_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.6.0",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.66.5",
        "transformers==5.14.1",
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
)

vienna_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "ViennaRNA==2.7.0",
        "tqdm==4.66.5",
    )
)

vol = modal.Volume.from_name("causal-rna-htt-all-models", create_if_missing=True)

# -- Shared experiment logic -------------------------------------------------


def _run_htt_experiment(adapter, model_name: str, device: str) -> dict:
    """Run HTT repeat-expansion embedding distance experiment for one model.

    Called from within Modal functions where numpy/scipy/tqdm are installed.
    Returns result dict ready for JSON serialization.
    """
    from datetime import datetime

    import numpy as np
    from scipy.stats import spearmanr
    from tqdm import tqdm

    use_rna = model_name in RNA_MODELS

    # {layer_idx: {repeat_count: mean_pooled_embedding}}
    embeddings = {}
    succeeded_counts = []
    failed_counts = {}

    for n in tqdm(REPEAT_COUNTS, desc=f"{model_name} HTT repeats"):
        seq = HTT_LEFT_FLANK + "CAG" * n + HTT_RIGHT_FLANK
        if use_rna:
            seq = seq.replace("T", "U")

        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] {model_name}: N={n}, seq_len={len(seq)}")

        tokens = adapter.tokenize(seq)
        if device == "cuda":
            tokens = tokens.to(device)

        try:
            all_layers = adapter.get_all_layer_embeddings(tokens)
            succeeded_counts.append(n)
        except Exception as e:
            msg = str(e)
            print(f"[{ts}] {model_name}: N={n} FAILED: {msg}")
            failed_counts[str(n)] = msg
            continue

        for layer_idx, layer_emb in enumerate(all_layers):
            pooled = layer_emb.float().mean(dim=0).cpu().numpy()
            if layer_idx not in embeddings:
                embeddings[layer_idx] = {}
            embeddings[layer_idx][n] = pooled

    if 10 not in succeeded_counts:
        return {
            "model": model_name,
            "error": "Baseline N=10 failed -- cannot compute distances",
            "failed_counts": failed_counts,
        }

    # Cosine distances from baseline (N=10) at each layer
    distances_per_layer = {}
    for layer_idx in sorted(embeddings.keys()):
        baseline = embeddings[layer_idx][10]
        baseline_norm = np.linalg.norm(baseline)
        distances = []
        for n in succeeded_counts:
            emb = embeddings[layer_idx][n]
            emb_norm = np.linalg.norm(emb)
            cos_sim = np.dot(baseline, emb) / (baseline_norm * emb_norm + 1e-10)
            distances.append(float(1.0 - cos_sim))
        distances_per_layer[str(layer_idx)] = distances

    # Best layer: highest |Spearman rho| between repeat count and distance
    best_layer = 0
    best_corr = 0.0
    for layer_idx_str, dists in distances_per_layer.items():
        if len(dists) < 3:
            continue
        rho, _ = spearmanr(succeeded_counts, dists)
        if not np.isnan(rho) and abs(rho) > abs(best_corr):
            best_corr = float(rho)
            best_layer = int(layer_idx_str)

    result = {
        "model": model_name,
        "repeat_counts": succeeded_counts,
        "distances_per_layer": distances_per_layer,
        "best_layer": best_layer,
        "best_correlation": best_corr,
        "distances_at_best_layer": distances_per_layer.get(str(best_layer), []),
        "n_layers_returned": len(embeddings),
        "failed_counts": failed_counts,
    }

    ts = datetime.now().strftime("%H:%M:%S")
    print(
        f"[{ts}] {model_name}: best_layer={best_layer}, "
        f"best_corr={best_corr:.4f}, "
        f"{len(succeeded_counts)}/{len(REPEAT_COUNTS)} counts succeeded"
    )

    return result


# -- Modal functions ---------------------------------------------------------


@app.function(
    image=multimol_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_multimol_htt(model_name: str):
    import os
    import sys

    os.chdir("/root/project")
    sys.path.insert(0, "/root/project")

    import json
    from datetime import datetime
    from pathlib import Path

    import torch

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] HTT experiment: {model_name} (multimol image)")

    from multi_model_audit import (
        ERNIERNAAdapter,
        RiNALMoAdapter,
        SpliceBERTAdapter,
        UTRLMAdapter,
    )

    adapter_map = {
        "ernierna": ERNIERNAAdapter,
        "splicebert": SpliceBERTAdapter,
        "rinalmo": RiNALMoAdapter,
        "utrlm": UTRLMAdapter,
    }

    device = "cuda" if torch.cuda.is_available() else "cpu"

    adapter = adapter_map[model_name]()
    try:
        adapter.load()
    except Exception as e:
        print(f"[{datetime.now().strftime('%H:%M:%S')}] FAILED to load {model_name}: {e}")
        result = {"model": model_name, "error": f"Load failed: {e}"}
        out_dir = Path(f"/results/htt_{model_name}_{timestamp}")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{model_name}_htt_distances.json"
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2, default=str)
        vol.commit()
        return {"model": model_name, "path": str(out_path), "error": str(e)}

    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    result = _run_htt_experiment(adapter, model_name, device)
    result["timestamp"] = timestamp
    result["torch_version"] = torch.__version__
    result["device"] = device

    out_dir = Path(f"/results/htt_{model_name}_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_name}_htt_distances.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    vol.commit()

    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} HTT COMPLETE. Saved: {out_path}")
    return {
        "model": model_name,
        "path": str(out_path),
        "best_corr": result.get("best_correlation"),
    }


@app.function(
    image=dnabert2_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_dnabert2_htt():
    import os
    import sys

    os.chdir("/root/project")
    sys.path.insert(0, "/root/project")

    import json
    from datetime import datetime
    from pathlib import Path

    import torch

    model_name = "dnabert2"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] HTT experiment: {model_name}")

    from multi_model_audit import DNABERT2Adapter

    device = "cuda" if torch.cuda.is_available() else "cpu"

    adapter = DNABERT2Adapter()
    try:
        adapter.load()
    except Exception as e:
        print(f"[{datetime.now().strftime('%H:%M:%S')}] FAILED to load {model_name}: {e}")
        result = {"model": model_name, "error": f"Load failed: {e}"}
        out_dir = Path(f"/results/htt_{model_name}_{timestamp}")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{model_name}_htt_distances.json"
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2, default=str)
        vol.commit()
        return {"model": model_name, "path": str(out_path), "error": str(e)}

    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    result = _run_htt_experiment(adapter, model_name, device)
    result["timestamp"] = timestamp
    result["torch_version"] = torch.__version__
    result["device"] = device
    result["note"] = "DNABERT-2 adapter returns only final hidden state (1 layer)"

    out_dir = Path(f"/results/htt_{model_name}_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_name}_htt_distances.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    vol.commit()

    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} HTT COMPLETE. Saved: {out_path}")
    return {
        "model": model_name,
        "path": str(out_path),
        "best_corr": result.get("best_correlation"),
    }


@app.function(
    image=transformers_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_transformers_htt(model_name: str):
    import os
    import sys

    os.chdir("/root/project")
    sys.path.insert(0, "/root/project")

    import json
    from datetime import datetime
    from pathlib import Path

    import torch

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] HTT experiment: {model_name} (transformers image)")

    from multi_model_audit import (
        CaduceusAdapter,
        EvoAdapter,
        HyenaDNAAdapter,
        NTAdapter,
    )

    adapter_map = {
        "nt": NTAdapter,
        "hyenadna": HyenaDNAAdapter,
        "evo": EvoAdapter,
        "caduceus": CaduceusAdapter,
    }

    device = "cuda" if torch.cuda.is_available() else "cpu"

    adapter = adapter_map[model_name]()
    try:
        adapter.load()
    except Exception as e:
        print(f"[{datetime.now().strftime('%H:%M:%S')}] FAILED to load {model_name}: {e}")
        result = {"model": model_name, "error": f"Load failed: {e}"}
        out_dir = Path(f"/results/htt_{model_name}_{timestamp}")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{model_name}_htt_distances.json"
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2, default=str)
        vol.commit()
        return {"model": model_name, "path": str(out_path), "error": str(e)}

    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    result = _run_htt_experiment(adapter, model_name, device)
    result["timestamp"] = timestamp
    result["torch_version"] = torch.__version__
    result["device"] = device

    out_dir = Path(f"/results/htt_{model_name}_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_name}_htt_distances.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    vol.commit()

    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} HTT COMPLETE. Saved: {out_path}")
    return {
        "model": model_name,
        "path": str(out_path),
        "best_corr": result.get("best_correlation"),
    }


@app.function(
    image=vienna_image,
    timeout=3600,
    volumes={"/results": vol},
)
def run_viennarna():
    import json
    from datetime import datetime
    from pathlib import Path

    import RNA
    from tqdm import tqdm

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] ViennaRNA MFE predictions for HTT variants")

    results = {}
    for n in tqdm(REPEAT_COUNTS, desc="ViennaRNA MFE"):
        seq = HTT_LEFT_FLANK + "CAG" * n + HTT_RIGHT_FLANK
        structure, mfe = RNA.fold(seq)
        results[str(n)] = {
            "structure": structure,
            "mfe": float(mfe),
            "length": len(seq),
        }
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] N={n}: len={len(seq)}, MFE={mfe:.2f}")

    out_dir = Path(f"/results/htt_viennarna_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "viennarna_mfe.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    vol.commit()

    print(f"[{datetime.now().strftime('%H:%M:%S')}] ViennaRNA COMPLETE. Saved: {out_path}")
    return {"path": str(out_path)}


# -- Entrypoint --------------------------------------------------------------


@app.local_entrypoint()
def main():
    from datetime import datetime

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] HTT trinucleotide repeat case study - all models")
    print(f"Repeat counts: {REPEAT_COUNTS}")
    print(f"Left flank:  {HTT_LEFT_FLANK}")
    print(f"Right flank: {HTT_RIGHT_FLANK}")
    print()

    handles = []

    # Multimol models (ERNIE-RNA, SpliceBERT, RiNALMo, UTR-LM)
    for model_name in MULTIMOL_MODELS:
        fc = run_multimol_htt.spawn(model_name=model_name)
        print(f"  Spawned {model_name} (multimol image) on A10G: {fc.object_id}")
        handles.append((model_name, fc))

    # DNABERT-2 (separate image with patched flash_attn)
    fc_dnabert2 = run_dnabert2_htt.spawn()
    print(f"  Spawned dnabert2 (dnabert2 image) on A10G: {fc_dnabert2.object_id}")
    handles.append(("dnabert2", fc_dnabert2))

    # Transformers models (NT v2, HyenaDNA, Evo, Caduceus)
    for model_name in TRANSFORMERS_MODELS:
        fc = run_transformers_htt.spawn(model_name=model_name)
        print(f"  Spawned {model_name} (transformers image) on A10G: {fc.object_id}")
        handles.append((model_name, fc))

    # ViennaRNA (CPU only)
    fc_vienna = run_viennarna.spawn()
    print(f"  Spawned ViennaRNA (CPU only): {fc_vienna.object_id}")
    handles.append(("viennarna", fc_vienna))

    print(f"\n{len(handles)} containers launched.")
    print("Use 'modal app logs causal-rna-htt-all-models' to monitor.")
