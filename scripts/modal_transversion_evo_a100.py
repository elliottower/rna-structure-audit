"""Modal wrapper: transversion control for Evo on A100.

Evo (7B) OOM'd on A10G (24GB). Needs A100 (80GB).

Preregistered prediction:
  - Evo transversion ratio within 10% of WC ratio (1.408)

Usage:
    cd /path/to/causal-rna
    modal run --detach scripts/modal_transversion_evo_a100.py
"""

import modal

app = modal.App("causal-rna-transversion-evo-a100")

TRANSVERSION_COMPLEMENT = {"A": "C", "U": "C", "C": "A", "G": "U", "T": "C"}

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

vol = modal.Volume.from_name("causal-rna-transversion-evo-a100", create_if_missing=True)


@app.function(
    image=evo_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_evo_transversion(seed: int = 42):
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
    print(f"[{timestamp}] Transversion control: evo (A100)")
    print(f"  Patching COMPLEMENT -> {TRANSVERSION_COMPLEMENT}")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families

    import phases_1_to_5
    phases_1_to_5.COMPLEMENT = TRANSVERSION_COMPLEMENT
    from phases_1_to_5 import run_mutation_sensitivity

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = load_adapter("evo")
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] evo loaded on {device}")

    out_dir = Path(f"/results/evo_transversion_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity with transversions...")
    r_mutation = run_mutation_sensitivity(adapter, "evo", families, device=device)

    output = {
        "model": "evo",
        "experiment": "transversion_control_revision",
        "description": (
            "Transversion mutations (purine->pyrimidine) on A100. "
            f"COMPLEMENT dict patched to: {TRANSVERSION_COMPLEMENT}"
        ),
        "complement_used": TRANSVERSION_COMPLEMENT,
        "timestamp": timestamp,
        "seed": seed,
        "torch_version": torch.__version__,
        "n_families": len(families),
        "mutation_results": r_mutation,
    }

    out_path = out_dir / f"evo_transversion.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    vol.commit()

    n_exc_nuc = r_mutation.get("n_exceeding_nuc_null", 0)

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] evo transversion control COMPLETE.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}")
    print(f"  Nuc null: {n_exc_nuc} families exceed")
    print(f"  Saved: {out_path}")

    return {
        "model": "evo",
        "path": str(out_path),
        "mean_ratio": r_mutation.get("mean_best_ratio"),
        "n_exc_nuc": n_exc_nuc,
    }


@app.local_entrypoint()
def main(seed: int = 42):
    print("Transversion control: Evo on A100")
    print(f"  Transversion COMPLEMENT: {TRANSVERSION_COMPLEMENT}")
    print(f"  Seed: {seed}")
    print()

    fc = run_evo_transversion.spawn(seed=seed)
    print(f"  Spawned evo (A100): {fc.object_id}")
    print("Use 'modal app logs causal-rna-transversion-evo-a100' to monitor.")
