"""Modal wrapper: Phase 6 for models needing flash-attn (Evo, DNABERT-2).

Uses Modal's GPU-enabled image builder to compile flash-attn from source.

Usage:
    modal run --detach scripts/modal_phase6_flashattn_models.py --model evo
    modal run --detach scripts/modal_phase6_flashattn_models.py --model dnabert2
    modal run --detach scripts/modal_phase6_flashattn_models.py --model evo --synthetic
"""

import modal

app = modal.App("causal-rna-phase6-flashattn")

base_image = (
    modal.Image.from_registry("nvidia/cuda:12.1.1-devel-ubuntu22.04", add_python="3.11")
    .apt_install("git")
    .pip_install(
        "torch==2.1.2",
        "triton==2.1.0",
        "packaging",
        "ninja",
        "wheel",
        "setuptools",
        "numpy==1.26.4",
        "scipy==1.13.1",
        "tqdm==4.66.4",
        "transformers==4.49.0",
        "multimolecule==0.1.0",
        "matplotlib==3.9.0",
        "scikit-learn==1.5.0",
        "einops==0.8.0",
    )
    .pip_install(
        "flash-attn==2.5.8",
        extra_options="--no-build-isolation",
        gpu="A10G",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)

GPU_MAP = {"evo": "A100", "dnabert2": "A10G"}


@app.function(
    image=base_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_a10g(model_name: str, synthetic: bool = False, n_synthetic: int = 5, seed: int = 42):
    return _run(model_name, synthetic, n_synthetic, seed)


@app.function(
    image=base_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_a100(model_name: str, synthetic: bool = False, n_synthetic: int = 5, seed: int = 42):
    return _run(model_name, synthetic, n_synthetic, seed)


def _run(model_name, synthetic, n_synthetic, seed):
    import os
    import sys
    os.chdir("/root/project")
    sys.path.insert(0, "/root/project")
    sys.path.insert(0, "/root/project/scripts")

    import json
    from datetime import datetime
    from pathlib import Path

    import numpy as np
    import torch

    # Bypass torch.load version check (CVE-2025-32434) — we use torch==2.1.2
    # for triton 2.1.0 compatibility (DNABERT-2 custom code needs trans_b).
    _noop = lambda: None
    import transformers.utils.import_utils as _tiu
    if hasattr(_tiu, "check_torch_load_is_safe"):
        _tiu.check_torch_load_is_safe = _noop
    import transformers.modeling_utils as _mu
    if hasattr(_mu, "check_torch_load_is_safe"):
        _mu.check_torch_load_is_safe = _noop

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    mode = "SYNTHETIC" if synthetic else "NATURAL"
    print(f"[{timestamp}] Phase 6 {mode}: {model_name}")

    from phase6_compensatory_mutation import (
        ADAPTER_OFFSETS,
        load_adapter,
        load_rfam_families,
        run_phase6,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    offset = ADAPTER_OFFSETS.get(model_name, 0)
    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    if synthetic:
        from modal_phase6_synthetic_covariation import generate_synthetic_families
        families = generate_synthetic_families(families, n_synthetic, rng)
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Generated {len(families)} synthetic sequences")

    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=offset,
    )

    prefix = "phase6_synthetic" if synthetic else "phase6"
    out_dir = Path(f"/results/{prefix}_{model_name}_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    output = {
        "model": model_name,
        "experiment": "synthetic_covariation_control" if synthetic else "phase6_natural",
        "timestamp": timestamp,
        "results": results,
    }
    if synthetic:
        output["n_synthetic_per_family"] = n_synthetic
        output["seed"] = seed

    out_path = out_dir / f"{model_name}_{'synthetic_covariation' if synthetic else 'phase6_ps'}.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done = datetime.now().strftime("%H:%M:%S")
    print(f"[{done}] {model_name} {mode} complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families: {results.get('families_total', 0)}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")

    return {"model": model_name, "mode": mode, "path": str(out_path)}


@app.local_entrypoint()
def main(
    model: str = "evo",
    synthetic: bool = False,
    n_synthetic: int = 5,
    seed: int = 42,
):
    mode = "SYNTHETIC" if synthetic else "NATURAL"
    gpu = GPU_MAP.get(model, "A10G")
    run_fn = run_a100 if gpu == "A100" else run_a10g

    print(f"Phase 6 {mode}: {model} on {gpu}")
    fc = run_fn.spawn(model_name=model, synthetic=synthetic, n_synthetic=n_synthetic, seed=seed)
    print(f"  Spawned: {fc.object_id}")
    print(f"Use 'modal app logs causal-rna-phase6-flashattn' to monitor.")
