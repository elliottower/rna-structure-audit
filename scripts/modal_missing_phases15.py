"""Modal wrapper: re-run Phase 1-5 for 4 models missing expanded-Rfam results.

Missing models: rnafm, nt, hyenadna, evo.
Caduceus excluded (mamba-ssm CUDA build requirement — deferred).

Saves JSON results to Modal volume, then downloads locally.

Usage:
    cd /path/to/causal-rna
    modal run --detach scripts/modal_missing_phases15.py
"""

import modal

app = modal.App("causal-rna-missing-phases15")

base_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.4.1",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.66.5",
        "transformers==4.44.2",
        "multimolecule==0.2.0",
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)

vol = modal.Volume.from_name("causal-rna-missing-phases15", create_if_missing=True)

MODELS_TO_RUN = ["evo"]

evo_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.1.0-devel-ubuntu22.04",
        add_python="3.11",
    )
    .apt_install("git")
    .pip_install(
        "torch==2.4.1",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.66.5",
        "transformers==4.44.2",
        "multimolecule==0.2.0",
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
        "packaging==24.1",
        "ninja==1.11.1.1",
        "wheel==0.44.0",
        "setuptools==75.1.0",
    )
    .run_commands(
        "pip install --no-build-isolation flash-attn==2.6.3",
        gpu="A100",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)
GPU_OVERRIDES = {"evo": "A100"}


@app.function(
    image=base_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_phases15_a10g(model_name: str, seed: int = 42):
    return _run_phases15(model_name, seed)


@app.function(
    image=evo_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_phases15_a100(model_name: str, seed: int = 42):
    return _run_phases15(model_name, seed)


def _run_phases15(model_name, seed):
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

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] Phase 1-5 only: starting {model_name}")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Loaded {len(families)} families")

    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] {model_name} loaded on {device}")

    out_dir = Path(f"/results/{model_name}_phases15_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running mutation sensitivity...")
    r_mutation = run_mutation_sensitivity(adapter, model_name, families, device=device)

    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running attention contact...")
    r_attention = run_attention_contact(adapter, model_name, families, device=device)

    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running structure probing...")
    r_probing = run_structure_probing(adapter, model_name, families, device=device)

    output = {
        "model": model_name,
        "experiment": "phases_1_to_5",
        "timestamp": timestamp,
        "seed": seed,
        "torch_version": torch.__version__,
        "n_families": len(families),
        "mutation_trained": r_mutation,
        "attention_trained": r_attention,
        "probing": r_probing,
    }

    out_path = out_dir / f"{model_name}_phases_1_to_5.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    vol.commit()

    done_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{done_ts}] {model_name} Phase 1-5 COMPLETE.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}, "
          f"exceed_null={r_mutation.get('n_exceeding_nuc_null', r_mutation.get('families_exceeding_nuc_null', 0))}")
    if r_attention.get("skipped"):
        print(f"  Attention: {r_attention['skipped']}")
    else:
        per_rna = r_attention.get("per_rna", {})
        corrs = [v.get("best_corr", 0) for v in per_rna.values() if isinstance(v, dict)]
        if corrs:
            print(f"  Attention: mean rho = {sum(corrs)/len(corrs):.4f} over {len(corrs)} families")
    if "best_accuracy" in r_probing:
        print(f"  Probing: best_accuracy={r_probing['best_accuracy']:.3f} at layer {r_probing.get('best_layer')}")
    print(f"  Saved: {out_path}")

    return {
        "model": model_name,
        "path": str(out_path),
        "mean_ratio": r_mutation.get("mean_best_ratio"),
    }


@app.local_entrypoint()
def main(seed: int = 42):
    print(f"Launching Phase 1-5 for {len(MODELS_TO_RUN)} missing models: {MODELS_TO_RUN}")
    print(f"Seed: {seed}")

    handles = []
    for model_name in MODELS_TO_RUN:
        if GPU_OVERRIDES.get(model_name) == "A100":
            run_fn = run_phases15_a100
            gpu = "A100"
        else:
            run_fn = run_phases15_a10g
            gpu = "A10G"

        fc = run_fn.spawn(model_name=model_name, seed=seed)
        print(f"  Spawned {model_name} on {gpu}: {fc.object_id}")
        handles.append((model_name, fc))

    print(f"\n{len(handles)} containers launched. Use 'modal app logs causal-rna-missing-phases15' to monitor.")
