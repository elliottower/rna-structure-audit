"""Modal wrapper: run Phase 1-5 + Phase 6 for Caduceus.

Caduceus requires mamba-ssm which needs CUDA compilation (nvcc).
Uses nvidia/cuda devel base image instead of debian_slim.

Usage:
    cd /path/to/causal-rna
    modal run --detach scripts/modal_caduceus_phases.py
"""

import modal

app = modal.App("causal-rna-caduceus")

caduceus_image = (
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
        "pip install --no-build-isolation causal-conv1d==1.4.0",
        gpu="A10G",
    )
    .run_commands(
        "pip install --no-build-isolation mamba-ssm==2.2.4",
        gpu="A10G",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)

vol = modal.Volume.from_name("causal-rna-caduceus-results", create_if_missing=True)


@app.function(
    image=caduceus_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_caduceus(seed: int = 42):
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

    model_name = "caduceus"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] Starting Caduceus (Phase 1-5 + Phase 6)")

    from phase6_compensatory_mutation import (
        ADAPTER_OFFSETS,
        load_adapter,
        load_rfam_families,
        run_phase6,
    )
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    offset = ADAPTER_OFFSETS.get(model_name, 0)
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Loaded {len(families)} families")

    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Caduceus loaded on {device}")

    out_dir = Path(f"/results/caduceus_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1-5
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running mutation sensitivity...")
    r_mutation = run_mutation_sensitivity(adapter, model_name, families, device=device)

    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running attention contact...")
    r_attention = run_attention_contact(adapter, model_name, families, device=device)

    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running structure probing...")
    r_probing = run_structure_probing(adapter, model_name, families, device=device)

    p15_output = {
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

    p15_path = out_dir / f"{model_name}_phases_1_to_5.json"
    with open(p15_path, "w") as f:
        json.dump(p15_output, f, indent=2, default=str)
    vol.commit()
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Phase 1-5 saved: {p15_path}")

    # Phase 6
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running Phase 6: perturbation specificity...")
    p6_results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=offset,
    )

    p6_output = {
        "model": model_name,
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "offset": offset,
        "seed": seed,
        "results": p6_results,
    }

    p6_path = out_dir / f"{model_name}_phase6_ps.json"
    with open(p6_path, "w") as f:
        json.dump(p6_output, f, indent=2)
    vol.commit()

    done_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{done_ts}] Caduceus ALL PHASES complete.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}")
    print(f"  Phase 6 mean PS: {p6_results.get('mean_best_ps', 'N/A')}")

    return {
        "model": model_name,
        "phase15_path": str(p15_path),
        "phase6_path": str(p6_path),
        "mean_ratio": r_mutation.get("mean_best_ratio"),
        "mean_ps": p6_results.get("mean_best_ps"),
    }


@app.local_entrypoint()
def main(seed: int = 42):
    print("Launching Caduceus (Phase 1-5 + Phase 6) on A10G with CUDA devel image")
    fc = run_caduceus.spawn(seed=seed)
    print(f"  Spawned: {fc.object_id}")
    print("Use 'modal app logs causal-rna-caduceus' to monitor.")
