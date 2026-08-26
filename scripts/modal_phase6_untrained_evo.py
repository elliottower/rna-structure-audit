"""Modal wrapper: Phase 6 randomized-weight control for Evo on A100.

Uses the same CUDA devel image that worked for transversion_evo_a100.

Usage:
    modal run --detach scripts/modal_phase6_untrained_evo.py
"""

import modal

app = modal.App("causal-rna-phase6-untrained-evo")

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
    .run_commands(
        "pip install --no-build-isolation flash-attn==2.6.3",
        gpu="A10G",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)


@app.function(
    image=evo_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_untrained_evo(seed: int = 42):
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
    print(f"[{timestamp}] Phase 6 UNTRAINED Evo (A100)")

    from multi_model_audit import EvoAdapter
    from phase6_compensatory_mutation import load_rfam_families, run_phase6

    device = "cuda"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = EvoAdapter()
    adapter.load()

    n_params = sum(p.numel() for p in adapter.model.parameters())
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded Evo ({n_params:,} params)")

    with torch.no_grad():
        for name, param in adapter.model.named_parameters():
            if param.dim() >= 2:
                torch.nn.init.xavier_normal_(param)
            else:
                torch.nn.init.normal_(param, std=0.02)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Randomized all {n_params:,} parameters")

    adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Moved to {device}")

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running Phase 6 on untrained Evo...")
    results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=0,
    )

    output = {
        "model": "evo_untrained",
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "control_type": "randomized_weights",
        "description": (
            "EvoAdapter with all weights randomized "
            "(xavier_normal_ for matrices, normal_(std=0.02) for vectors). "
            "Architecture preserved, learned weights destroyed."
        ),
        "seed": seed,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_evo_untrained_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "evo_untrained_phase6_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] UNTRAINED Evo complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families total: {results.get('families_total', 0)}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")
    print(f"  Saved to: {out_path}")

    return {"path": str(out_path), "mean_ps": results.get("mean_best_ps")}


@app.local_entrypoint()
def main(seed: int = 42):
    print("Launching Phase 6 on UNTRAINED Evo (A100, CUDA devel image)")
    result = run_untrained_evo.remote(seed=seed)
    print(f"Done! Result: {result}")
