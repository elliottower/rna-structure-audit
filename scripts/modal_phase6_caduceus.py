"""Modal wrapper: Phase 6 for Caduceus using CUDA devel image.

Caduceus requires mamba-ssm which compiles CUDA kernels at install time.
Standard Modal images lack nvcc; this script uses nvidia/cuda devel base.

Usage:
    modal run --detach scripts/modal_phase6_caduceus.py
"""

import modal

app = modal.App("causal-rna-phase6-caduceus")

base_image = (
    modal.Image.from_registry(
        "pytorch/pytorch:2.4.0-cuda12.4-cudnn9-devel",
    )
    .run_commands("pip install packaging ninja")
    .pip_install(
        "causal-conv1d>=1.2.0",
        "mamba-ssm>=1.2.0",
    )
    .pip_install(
        "numpy",
        "scipy",
        "tqdm",
        "transformers>=4.40.0,<5.0",
        "matplotlib",
        "scikit-learn",
        "einops",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)


@app.function(
    image=base_image,
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

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] Phase 6: Caduceus")

    from phase6_compensatory_mutation import load_rfam_families, run_phase6
    from multi_model_audit import CaduceusAdapter

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = CaduceusAdapter()
    adapter.load()
    if device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Caduceus loaded on {device}")

    results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=0,
    )

    output = {
        "model": "caduceus",
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "preregistration": "PREREGISTRATION_PHASE6_V2.md",
        "seed": seed,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_caduceus_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "caduceus_phase6_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] Caduceus complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families total: {results.get('families_total', 0)}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")

    return {"path": str(out_path), "mean_ps": results.get("mean_best_ps")}


@app.local_entrypoint()
def main(seed: int = 42):
    print("Launching Phase 6 for Caduceus (CUDA devel image for mamba-ssm)")
    fc = run_caduceus.spawn(seed=seed)
    print(f"Spawned on A10G: {fc.object_id}")
    print("Use 'modal app logs causal-rna-phase6-caduceus' to monitor.")
