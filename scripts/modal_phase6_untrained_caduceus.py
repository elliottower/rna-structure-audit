"""Modal wrapper: Phase 6 randomized-weight control for Caduceus.

Caduceus requires mamba-ssm which compiles CUDA kernels at install time,
so it needs a CUDA devel base image instead of debian_slim.

Usage:
    modal run --detach scripts/modal_phase6_untrained_caduceus.py
"""

import modal

app = modal.App("causal-rna-phase6-untrained-caduceus")

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
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
        "packaging",
        "ninja",
        "wheel",
        "setuptools",
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
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)


@app.function(
    image=caduceus_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_untrained_caduceus(seed: int = 42):
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
    print(f"[{timestamp}] Phase 6 UNTRAINED Caduceus (CUDA devel + mamba-ssm)")

    from multi_model_audit import CaduceusAdapter
    from phase6_compensatory_mutation import load_rfam_families, run_phase6

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = CaduceusAdapter()
    adapter.load()

    n_params = sum(p.numel() for p in adapter.model.parameters())
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded trained Caduceus ({n_params:,} params)")

    with torch.no_grad():
        for name, param in adapter.model.named_parameters():
            if param.dim() >= 2:
                torch.nn.init.xavier_normal_(param)
            else:
                torch.nn.init.normal_(param, std=0.02)

    n_params_after = sum(p.numel() for p in adapter.model.parameters())
    assert n_params == n_params_after, "Parameter count changed during randomization"
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Randomized all {n_params_after:,} parameters")

    if device == "cuda":
        adapter.model = adapter.model.to(device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running Phase 6 on untrained Caduceus...")
    results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=0,
    )

    output = {
        "model": "caduceus_untrained",
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "control_type": "randomized_weights",
        "description": (
            "CaduceusAdapter with all weights randomized "
            "(xavier_normal_ for matrices, normal_(std=0.02) for vectors). "
            "Architecture preserved, learned weights destroyed."
        ),
        "seed": seed,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_caduceus_untrained_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "caduceus_untrained_phase6_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] UNTRAINED Caduceus complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families total: {results.get('families_total', 0)}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")
    print(f"  Saved to: {out_path}")

    return {"path": str(out_path), "mean_ps": results.get("mean_best_ps")}


@app.local_entrypoint()
def main(seed: int = 42):
    print("Launching Phase 6 on UNTRAINED Caduceus (CUDA devel image for mamba-ssm)")
    fc = run_untrained_caduceus.spawn(seed=seed)
    print(f"Spawned on A10G: {fc.object_id}")
    print("Use 'modal app logs causal-rna-phase6-untrained-caduceus' to monitor.")
