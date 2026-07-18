"""Modal wrapper: Phase 1-5 for RiNALMo (650M).

Computes Rung 1-2 data:
  - Mutation sensitivity (mean ratio, families > null)
  - Attention-contact correlation (Spearman rho)
  - Structure probing (logistic regression)

Usage:
    modal run --detach scripts/modal_rinalmo_phases15.py
"""

import modal

app = modal.App("causal-rna-rinalmo-phases15")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.6.0",
        "numpy==2.1.3",
        "scipy==1.14.1",
        "tqdm==4.67.1",
        "transformers==5.13.1",
        "matplotlib==3.9.3",
        "scikit-learn==1.6.0",
        "einops==0.8.0",
        "multimolecule==0.1.0",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py", copy=True)
    .add_local_dir("scripts", "/root/project/scripts", copy=True)
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families", copy=True)
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)


@app.function(image=image, gpu="A100", timeout=86400, volumes={"/results": vol})
def run_rinalmo_phases15(seed: int = 42):
    import json
    import os
    import sys
    from datetime import datetime
    from pathlib import Path

    import numpy as np
    import torch

    os.chdir("/root/project")
    sys.path.insert(0, "/root/project")
    sys.path.insert(0, "/root/project/scripts")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] Phase 1-5: rinalmo (650M)")
    print(f"  torch={torch.__version__}, cuda={torch.cuda.is_available()}")

    from phase6_compensatory_mutation import load_rfam_families
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    from multi_model_audit import RiNALMoAdapter
    adapter = RiNALMoAdapter()
    adapter.load()
    if device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] RiNALMo loaded on {device}")

    # Sanity check
    test_tokens = adapter.tokenize("ACGUACGU").to(device)
    test_out = adapter.get_all_layer_embeddings(test_tokens)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Sanity check: {len(test_out)} layers, shape {test_out[0].shape}")

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity...")
    r_mutation = run_mutation_sensitivity(adapter, "rinalmo", families, device=device)
    print(f"  mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}, "
          f"exceed_null={r_mutation.get('n_exceeding_nuc_null', 0)}/{r_mutation.get('n_families_scored', 0)}")

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running attention-contact correlation...")
    r_attention = run_attention_contact(adapter, "rinalmo", families, device=device)
    if r_attention.get("skipped"):
        print(f"  Attention: {r_attention['skipped']}")
    else:
        corrs = [v.get("best_corr", 0) for v in r_attention.get("per_rna", {}).values()
                 if not v.get("skipped")]
        if corrs:
            print(f"  Best mean corr: {np.mean(corrs):.4f}, max: {max(corrs):.4f}")

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running structure probing...")
    r_probing = run_structure_probing(adapter, "rinalmo", families, device=device)
    if "best_accuracy" in r_probing:
        print(f"  best_accuracy={r_probing['best_accuracy']:.3f} at layer {r_probing.get('best_layer')}")

    out_dir = Path(f"/results/rinalmo_phases15_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    output = {
        "model": "rinalmo",
        "experiment": "phases_1_to_5",
        "timestamp": timestamp,
        "seed": seed,
        "mutation_trained": r_mutation,
        "attention_trained": r_attention,
        "probing": r_probing,
    }

    out_path = out_dir / "rinalmo_phases_1_to_5.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)

    vol.commit()

    done = datetime.now().strftime("%H:%M:%S")
    print(f"[{done}] rinalmo Phase 1-5 complete.")

    return {"model": "rinalmo", "path": str(out_path)}


@app.local_entrypoint()
def main(seed: int = 42):
    print("Phase 1-5: rinalmo on A100")
    fc = run_rinalmo_phases15.spawn(seed=seed)
    print(f"  Spawned: {fc.object_id}")
    print("Use 'modal app logs causal-rna-rinalmo-phases15' to monitor.")
