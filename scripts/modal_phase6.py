"""Modal wrapper for Phase 6 perturbation specificity test.

Runs one model per container. Results saved to Modal volume and
downloaded locally after completion.

Usage:
    modal run --detach scripts/modal_phase6.py
    modal run --detach scripts/modal_phase6.py --models rnafm,hyenadna
    modal run --detach scripts/modal_phase6.py --models rnafm --no-null
"""

import modal

app = modal.App("causal-rna-phase6-ps")

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
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)

caduceus_image = (
    modal.Image.from_registry("nvidia/cuda:12.1.0-devel-ubuntu22.04", add_python="3.11")
    .pip_install(
        "torch>=2.2.0",
        "numpy",
        "scipy",
        "tqdm",
        "transformers>=4.40.0",
        "matplotlib",
        "packaging",
        "ninja",
    )
    .pip_install("causal-conv1d", "mamba-ssm")
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)

AVAILABLE_MODELS = [
    "rnafm", "rinalmo", "utrlm", "ernierna", "splicebert",
    "nt", "hyenadna", "caduceus", "evo", "dnabert2",
]
CHARACTER_MODELS = [m for m in AVAILABLE_MODELS if m not in {"nt", "dnabert2"}]

GPU_OVERRIDES = {
    "evo": "A100",
    "rinalmo": "A100",
}
CADUCEUS_MODELS = {"caduceus"}


@app.function(
    image=base_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_a10g(model_name: str, seed: int = 42, compute_null: bool = True):
    return _run_model(model_name, seed, compute_null)


@app.function(
    image=base_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_a100(model_name: str, seed: int = 42, compute_null: bool = True):
    return _run_model(model_name, seed, compute_null)


@app.function(
    image=caduceus_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_caduceus(model_name: str, seed: int = 42, compute_null: bool = True):
    return _run_model(model_name, seed, compute_null)


def _run_model(model_name, seed, compute_null):
    import sys
    sys.path.insert(0, "/root/project")

    import json
    from datetime import datetime
    from pathlib import Path

    import numpy as np
    import torch

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] Phase 6 PS: starting {model_name}")

    sys.path.insert(0, "/root/project/scripts")
    from phase6_compensatory_mutation import (
        ADAPTER_OFFSETS,
        NON_CHARACTER_TOKENIZERS,
        load_adapter,
        load_rfam_families,
        run_phase6,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    offset = ADAPTER_OFFSETS.get(model_name, 0)
    tokenizer_caveated = model_name in NON_CHARACTER_TOKENIZERS

    np.random.seed(seed)

    print(f"[{timestamp}] Loading families...")
    families = load_rfam_families()
    print(f"[{timestamp}] Loaded {len(families)} families")

    print(f"[{timestamp}] Loading {model_name} (device={device}, offset={offset})...")
    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, 'model') and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)

    print(f"[{timestamp}] Running Phase 6...")
    results = run_phase6(
        adapter, families, device=device,
        compute_null=compute_null, offset=offset,
    )

    output = {
        "model": model_name,
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "preregistration": "PREREGISTRATION_PHASE6_V2.md",
        "offset": offset,
        "seed": seed,
        "tokenizer_caveated": tokenizer_caveated,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_{model_name}_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_name}_phase6_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{done_ts}] {model_name} complete. Saved to {out_path}")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families (confirmatory): {results.get('families_total', 0)}")
    print(f"  Exceeding null (primary): {results.get('families_exceeding_null_primary', 'N/A')}")

    return {"model": model_name, "path": str(out_path), "mean_ps": results.get("mean_best_ps")}


@app.local_entrypoint()
def main(
    models: str = ",".join(CHARACTER_MODELS),
    seed: int = 42,
    no_null: bool = False,
    allow_non_character: bool = False,
):
    model_list = [m.strip() for m in models.split(",")]

    print(f"Launching Phase 6 PS for {len(model_list)} models: {model_list}")
    print(f"Seed: {seed}, compute_null: {not no_null}")

    handles = []
    for model_name in model_list:
        if model_name not in AVAILABLE_MODELS:
            print(f"  SKIP {model_name}: unknown model")
            continue
        if model_name in {"nt", "dnabert2"} and not allow_non_character:
            print(f"  SKIP {model_name}: non-character tokenizer (use --allow-non-character)")
            continue

        if model_name in CADUCEUS_MODELS:
            run_fn = run_model_caduceus
            gpu = "A10G (caduceus image)"
        elif GPU_OVERRIDES.get(model_name) == "A100":
            run_fn = run_model_a100
            gpu = "A100"
        else:
            run_fn = run_model_a10g
            gpu = "A10G"
        fc = run_fn.spawn(model_name=model_name, seed=seed, compute_null=not no_null)
        print(f"  Spawned {model_name} on {gpu}: {fc.object_id}")
        handles.append((model_name, fc))

    print(f"\n{len(handles)} containers launched. Use 'modal app logs causal-rna-phase6-ps' to monitor.")
