"""Modal wrapper: Phase 6 for models with version-sensitive deps.

RiNALMo needs multimolecule==0.1.0 (0.2.0 breaks with newer transformers).
DNABERT-2 needs torch>=2.6 for torch.load security fix (CVE-2025-32434).

Usage:
    modal run --detach scripts/modal_phase6_compat.py --model rinalmo
    modal run --detach scripts/modal_phase6_compat.py --model dnabert2
"""

import modal

app = modal.App("causal-rna-phase6-compat")

rinalmo_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch>=2.2.0",
        "numpy",
        "scipy",
        "tqdm",
        "transformers>=4.40.0",
        "multimolecule==0.1.0",
        "matplotlib",
        "scikit-learn",
        "einops",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
)

dnabert2_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.6.0",
        "numpy",
        "scipy",
        "tqdm",
        "transformers>=4.40.0,<4.50",
        "matplotlib",
        "scikit-learn",
        "einops",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)


@app.function(image=rinalmo_image, gpu="A100", timeout=86400, volumes={"/results": vol})
def run_rinalmo(seed: int = 42):
    return _run("rinalmo", seed)


@app.function(image=dnabert2_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_dnabert2(seed: int = 42):
    return _run("dnabert2", seed)


def _run(model_name, seed):
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
    print(f"[{timestamp}] Phase 6 PS: {model_name}")

    from phase6_compensatory_mutation import (
        ADAPTER_OFFSETS,
        load_adapter,
        load_rfam_families,
        run_phase6,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    offset = ADAPTER_OFFSETS.get(model_name, 0)
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = load_adapter(model_name)
    if model_name == "dnabert2":
        from transformers import AutoModel, AutoTokenizer
        adapter.tokenizer = AutoTokenizer.from_pretrained(
            adapter.MODEL_ID, trust_remote_code=True,
        )
        adapter.model = AutoModel.from_pretrained(
            adapter.MODEL_ID, trust_remote_code=True,
            attn_implementation="eager",
        )
        adapter.model.eval()
    else:
        adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    results = run_phase6(adapter, families, device=device, compute_null=True, offset=offset)

    output = {
        "model": model_name,
        "timestamp": timestamp,
        "phase": 6,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_{model_name}_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_name}_phase6_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done = datetime.now().strftime("%H:%M:%S")
    print(f"[{done}] {model_name} Phase 6 complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families: {results.get('families_total', 0)}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")

    return {"model": model_name, "path": str(out_path)}


@app.local_entrypoint()
def main(model: str = "rinalmo", seed: int = 42):
    if model == "rinalmo":
        fc = run_rinalmo.spawn(seed=seed)
        gpu = "A100"
    elif model == "dnabert2":
        fc = run_dnabert2.spawn(seed=seed)
        gpu = "A10G"
    else:
        raise ValueError(f"Use --model rinalmo or --model dnabert2, got {model}")

    print(f"Phase 6: {model} on {gpu}")
    print(f"  Spawned: {fc.object_id}")
    print("Use 'modal app logs causal-rna-phase6-compat' to monitor.")
