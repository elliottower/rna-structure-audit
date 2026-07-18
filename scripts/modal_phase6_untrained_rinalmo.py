"""Modal wrapper: run Phase 6 on RiNALMo with randomized weights.

Tests whether RiNALMo's partner specificity (PS = 0.227) requires
learned weights. Mirrors the ERNIE-RNA untrained control (PS = 2.5e-8).

Usage:
    modal run --detach scripts/modal_phase6_untrained_rinalmo.py
"""

import modal

app = modal.App("causal-rna-phase6-untrained-rinalmo")

base_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.2.0",
        "numpy==1.26.4",
        "scipy==1.13.0",
        "tqdm==4.66.4",
        "transformers==4.44.0",
        "multimolecule==0.0.30",
        "matplotlib==3.9.0",
        "scikit-learn==1.5.0",
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
def run_untrained_rinalmo(seed: int = 42):
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
    print(f"[{timestamp}] Phase 6 UNTRAINED RiNALMo control")

    from multi_model_audit import RiNALMoAdapter
    from phase6_compensatory_mutation import load_rfam_families, run_phase6

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = RiNALMoAdapter()
    adapter.load()

    n_params_before = sum(p.numel() for p in adapter.model.parameters())
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded trained RiNALMo ({n_params_before:,} params)")

    with torch.no_grad():
        for name, param in adapter.model.named_parameters():
            if param.dim() >= 2:
                torch.nn.init.xavier_normal_(param)
            else:
                torch.nn.init.normal_(param, std=0.02)

    n_params_after = sum(p.numel() for p in adapter.model.parameters())
    assert n_params_before == n_params_after, "Parameter count changed during randomization"
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Randomized all {n_params_after:,} parameters")

    if device == "cuda":
        adapter.model = adapter.model.to(device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running Phase 6 on untrained RiNALMo...")
    results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=0,
    )

    output = {
        "model": "rinalmo_untrained",
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "control_type": "randomized_weights",
        "description": "RiNALMo with all weights randomized (xavier_normal_ for matrices, normal_(std=0.02) for vectors). Architecture preserved, learned weights destroyed.",
        "preregistration": "PREREGISTRATION_PHASE6_UNTRAINED_RINALMO.md",
        "seed": seed,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_rinalmo_untrained_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "rinalmo_untrained_phase6_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] UNTRAINED RiNALMo complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families total: {results.get('families_total', 0)}")
    print(f"  Gate pass: {results.get('families_gate_pass', 'N/A')}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")

    return {"path": str(out_path), "mean_ps": results.get("mean_best_ps")}


@app.local_entrypoint()
def main(seed: int = 42):
    print("Launching Phase 6 on UNTRAINED RiNALMo (randomized weights)")
    print("This tests whether partner specificity requires learned weights")
    fc = run_untrained_rinalmo.spawn(seed=seed)
    print(f"Spawned on A10G: {fc.object_id}")
    print("Use 'modal app logs causal-rna-phase6-untrained-rinalmo' to monitor.")
