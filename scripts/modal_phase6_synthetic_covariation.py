"""Modal wrapper: Phase 6 on synthetic sequences (covariation control).

Tests whether partner specificity reflects evolutionary covariation or
geometric pairing knowledge. For each Rfam family, generates synthetic
sequences with the same secondary structure but random nucleotide
assignment (no evolutionary history). Runs all character-tokenized models.

Usage:
    modal run --detach scripts/modal_phase6_synthetic_covariation.py
    modal run --detach scripts/modal_phase6_synthetic_covariation.py --models ernierna,rnafm
"""

import modal

app = modal.App("causal-rna-phase6-synthetic-covariation")

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
        "einops",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)

CHARACTER_MODELS = ["rnafm", "rinalmo", "utrlm", "ernierna", "splicebert", "hyenadna", "evo"]
GPU_OVERRIDES = {"rinalmo": "A100", "evo": "A100"}
SKIP_MODELS = {"caduceus"}


def generate_synthetic_families(families, n_synthetic, rng):
    """For each family, generate n_synthetic random sequences preserving structure."""
    WC_PAIRS = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")]
    NUCS = ["A", "U", "G", "C"]

    synthetic = []
    for fam in families:
        structure = fam.get("dot_bracket", "")
        if not structure:
            continue

        for syn_idx in range(n_synthetic):
            seq = ["N"] * len(structure)
            stack = []
            for i, char in enumerate(structure):
                if char == "(":
                    stack.append(i)
                elif char == ")":
                    j = stack.pop()
                    pair = WC_PAIRS[rng.integers(len(WC_PAIRS))]
                    seq[j] = pair[0]
                    seq[i] = pair[1]
                else:
                    seq[i] = NUCS[rng.integers(len(NUCS))]

            synthetic.append({
                "name": f"{fam['name']}_syn{syn_idx}",
                "sequence": "".join(seq),
                "dot_bracket": structure,
                "original_family": fam["name"],
                "synthetic_index": syn_idx,
            })
    return synthetic


@app.function(
    image=base_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_a10g(model_name: str, n_synthetic: int = 5, seed: int = 42):
    return _run_model(model_name, n_synthetic, seed)


@app.function(
    image=base_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_a100(model_name: str, n_synthetic: int = 5, seed: int = 42):
    return _run_model(model_name, n_synthetic, seed)


def _run_model(model_name, n_synthetic, seed):
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
    print(f"[{timestamp}] Phase 6 SYNTHETIC COVARIATION: {model_name}")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families, run_phase6

    device = "cuda" if torch.cuda.is_available() else "cpu"
    rng = np.random.default_rng(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    synthetic_families = generate_synthetic_families(families, n_synthetic, rng)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Generated {len(synthetic_families)} synthetic sequences")

    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    results = run_phase6(
        adapter, synthetic_families, device=device,
        compute_null=True, offset=0,
    )

    output = {
        "model": model_name,
        "experiment": "synthetic_covariation_control",
        "timestamp": timestamp,
        "description": (
            "Phase 6 on synthetic sequences. Each Rfam family's dot-bracket "
            "structure preserved, nucleotides randomly assigned (uniform WC "
            "pairs at paired positions, uniform random at unpaired). Destroys "
            "evolutionary covariation while preserving exact structure."
        ),
        "n_synthetic_per_family": n_synthetic,
        "seed": seed,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_synthetic_{model_name}_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_name}_synthetic_covariation.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] {model_name} SYNTHETIC complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families total: {results.get('families_total', 0)}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")

    return {"model": model_name, "path": str(out_path), "mean_ps": results.get("mean_best_ps")}


@app.local_entrypoint()
def main(
    models: str = ",".join(CHARACTER_MODELS),
    n_synthetic: int = 5,
    seed: int = 42,
):
    model_list = [m.strip() for m in models.split(",")]

    print("Phase 6 SYNTHETIC COVARIATION CONTROL")
    print(f"  Models: {model_list}")
    print(f"  {n_synthetic} synthetic sequences per family, seed={seed}")

    handles = []
    for model_name in model_list:
        if model_name in SKIP_MODELS:
            print(f"  SKIP {model_name}")
            continue

        if GPU_OVERRIDES.get(model_name) == "A100":
            run_fn = run_model_a100
            gpu = "A100"
        else:
            run_fn = run_model_a10g
            gpu = "A10G"

        fc = run_fn.spawn(model_name=model_name, n_synthetic=n_synthetic, seed=seed)
        print(f"  Spawned {model_name} on {gpu}: {fc.object_id}")
        handles.append((model_name, fc))

    print(f"\n{len(handles)} containers launched.")
    print("Use 'modal app logs causal-rna-phase6-synthetic-covariation' to monitor.")
