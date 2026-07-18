"""Modal wrapper: run Phase 6 on ERNIE-RNA with randomized weights.

Tests whether the architectural base-pairing attention bias produces
partner specificity WITHOUT training — the critical control for
attributing ERNIE-RNA's Phase 6 result to learned representations
vs architectural inductive bias.

Usage:
    modal run --detach scripts/modal_phase6_untrained_ernierna.py
"""

import modal

app = modal.App("causal-rna-phase6-untrained-ernierna")

base_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch>=2.2.0",
        "numpy",
        "scipy",
        "tqdm",
        "transformers>=4.40.0",
        "multimolecule",
        "matplotlib",
        "scikit-learn",
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
def run_untrained_ernierna(seed: int = 42):
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
    print(f"[{timestamp}] Phase 6 UNTRAINED ERNIE-RNA control")

    from multi_model_audit import ERNIERNAAdapter
    from phase6_compensatory_mutation import load_rfam_families, run_phase6

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = ERNIERNAAdapter()
    adapter.load()

    n_params_before = sum(p.numel() for p in adapter.model.parameters())
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded trained ERNIE-RNA ({n_params_before:,} params)")

    # Randomize all weights while preserving architecture
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

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running Phase 6 on untrained ERNIE-RNA...")
    results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=0,
    )

    output = {
        "model": "ernierna_untrained",
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "control_type": "randomized_weights",
        "description": "ERNIE-RNA with all weights randomized (xavier_normal_ for matrices, normal_(std=0.02) for vectors). Architecture and attention bias structure preserved, learned weights destroyed.",
        "preregistration": "PREREGISTRATION_PHASE6_V2.md",
        "seed": seed,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_ernierna_untrained_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "ernierna_untrained_phase6_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] UNTRAINED ERNIE-RNA complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families total: {results.get('families_total', 0)}")
    print(f"  Gate pass: {results.get('families_gate_pass', 'N/A')}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")

    return {"path": str(out_path), "mean_ps": results.get("mean_best_ps")}


@app.local_entrypoint()
def main(seed: int = 42):
    print("Launching Phase 6 on UNTRAINED ERNIE-RNA (randomized weights)")
    print("This tests whether architectural attention bias alone produces partner specificity")
    fc = run_untrained_ernierna.spawn(seed=seed)
    print(f"Spawned on A10G: {fc.object_id}")
    print("Use 'modal app logs causal-rna-phase6-untrained-ernierna' to monitor.")
