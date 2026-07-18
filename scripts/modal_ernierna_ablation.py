"""Modal wrapper: ERNIE-RNA attention bias ablation.

ERNIE-RNA has a pairwise_bias_map buffer that encodes structural priors
(WC pairing preferences) into attention. This ablation zeros that buffer
to test whether structure awareness comes from the bias or learned weights.

Three conditions:
1. Full model (existing results — not re-run here)
2. Zeroed pairwise_bias_map (attention bias ablated)
3. Untrained (existing results — not re-run here)

Usage:
    cd /path/to/causal-rna
    modal run --detach scripts/modal_ernierna_ablation.py
"""

import modal

app = modal.App("causal-rna-ernierna-ablation")

image = (
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
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)

vol = modal.Volume.from_name("causal-rna-ernierna-ablation", create_if_missing=True)


@app.function(
    image=image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_ablation(seed: int = 42):
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

    from multi_model_audit import ERNIERNAAdapter
    from phase6_compensatory_mutation import load_rfam_families
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] ERNIE-RNA attention bias ablation")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = ERNIERNAAdapter()
    adapter.load()

    n_params = sum(p.numel() for p in adapter.model.parameters())
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded ERNIE-RNA ({n_params:,} params)")

    # Zero the pairwise bias map — this removes structural attention bias
    # while keeping all learned transformer weights intact
    original_bias = adapter.model.pairwise_bias_map.clone()
    bias_norm = original_bias.norm().item()
    adapter.model.pairwise_bias_map.zero_()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Zeroed pairwise_bias_map "
          f"(was norm={bias_norm:.4f}, shape={list(original_bias.shape)})")

    # Also zero the pairwise_bias_proj to be thorough — even with zero input,
    # the bias term in the projection could contribute
    for name, param in adapter.model.pairwise_bias_proj.named_parameters():
        param_norm = param.norm().item()
        param.data.zero_()
        print(f"  Zeroed pairwise_bias_proj.{name} (was norm={param_norm:.4f})")

    if device == "cuda":
        adapter.model = adapter.model.to(device)

    out_dir = Path(f"/results/ernierna_no_attn_bias_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    model_name = "ernierna_no_attn_bias"

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity...")
    r_mutation = run_mutation_sensitivity(adapter, "ernierna", families, device=device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running attention contact...")
    r_attention = run_attention_contact(adapter, "ernierna", families, device=device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running structure probing...")
    r_probing = run_structure_probing(adapter, "ernierna", families, device=device)

    output = {
        "model": model_name,
        "experiment": "ernierna_attention_bias_ablation",
        "ablation": "pairwise_bias_map and pairwise_bias_proj zeroed",
        "description": (
            "ERNIE-RNA with pairwise_bias_map buffer zeroed and "
            "pairwise_bias_proj weights zeroed. Architecture and all other "
            "learned weights preserved. Tests whether structure awareness "
            "comes from the attention bias or from learned representations."
        ),
        "original_bias_norm": bias_norm,
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

    n_exc_nuc = r_mutation.get("n_exceeding_nuc_null", 0)
    n_scored = r_mutation.get("n_families_scored", 0)
    n_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and "dinuc_null_95th" in v
    )
    n_exc_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and v.get("exceeds_dinuc_null")
    )

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] ERNIE-RNA (no attn bias) Phase 1-5 COMPLETE.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}")
    print(f"  Nuc null: {n_exc_nuc}/{n_scored} families exceed")
    print(f"  Dinuc null: {n_dinuc} tested, {n_exc_dinuc} exceed")
    if "best_accuracy" in r_probing:
        print(f"  Probing: acc={r_probing['best_accuracy']:.3f} at layer {r_probing.get('best_layer')}")
    print(f"  Saved: {out_path}")

    return {
        "model": model_name,
        "path": str(out_path),
        "mean_ratio": r_mutation.get("mean_best_ratio"),
        "n_exc_nuc": n_exc_nuc,
        "n_exc_dinuc": n_exc_dinuc,
    }


@app.local_entrypoint()
def main(seed: int = 42):
    print("Launching ERNIE-RNA attention bias ablation")
    print("  Condition: pairwise_bias_map + pairwise_bias_proj zeroed")
    fc = run_ablation.spawn(seed=seed)
    print(f"  Spawned on A10G: {fc.object_id}")
    print("Use 'modal app logs causal-rna-ernierna-ablation' to monitor.")
