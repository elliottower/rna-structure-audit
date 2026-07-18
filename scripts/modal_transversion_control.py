"""Modal wrapper: transversion mutation control for all models.

A reviewer noted that Watson-Crick complement mutations (A<->U, C<->G)
at stem positions still create valid base pairs (e.g., A-U -> U-A).
True structure disruptions require transversions: purine -> pyrimidine
or vice versa (A->C, G->U, etc.).

This script monkey-patches the COMPLEMENT dict in phases_1_to_5.py
to use transversion mutations, then runs run_mutation_sensitivity()
as normal. If models are genuinely structure-aware, transversions
at stem positions should produce larger embedding shifts than WC swaps.

Usage:
    cd /path/to/causal-rna
    modal run --detach scripts/modal_transversion_control.py
"""

import modal

app = modal.App("causal-rna-transversion-control")

TRANSVERSION_COMPLEMENT = {"A": "C", "U": "C", "C": "A", "G": "U", "T": "C"}

multimol_image = (
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

dnabert2_image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "torch==2.4.0",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.66.5",
        "transformers==4.28.0",
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
    )
    .add_local_file("scripts/patch_dnabert2_flash_attn.py", "/root/patch_dnabert2_flash_attn.py", copy=True)
    .run_commands("python /root/patch_dnabert2_flash_attn.py")
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
)

transformers_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.6.0",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.66.5",
        "transformers==5.14.1",
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
)

vol = modal.Volume.from_name("causal-rna-transversion-control", create_if_missing=True)

MULTIMOL_MODELS = ["ernierna", "splicebert", "rinalmo", "utrlm"]
TRANSFORMERS_MODELS = ["nt", "hyenadna", "evo", "caduceus"]


@app.function(
    image=multimol_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_multimol_transversion(model_name: str, seed: int = 42):
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
    print(f"[{timestamp}] Transversion control: {model_name}")
    print(f"  Patching COMPLEMENT -> {TRANSVERSION_COMPLEMENT}")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families

    import phases_1_to_5
    phases_1_to_5.COMPLEMENT = TRANSVERSION_COMPLEMENT
    from phases_1_to_5 import run_mutation_sensitivity

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    out_dir = Path(f"/results/{model_name}_transversion_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity with transversions...")
    r_mutation = run_mutation_sensitivity(adapter, model_name, families, device=device)

    output = {
        "model": model_name,
        "experiment": "transversion_control",
        "description": (
            "Mutation sensitivity with transversion mutations (purine->pyrimidine) "
            "instead of Watson-Crick complements. COMPLEMENT dict patched to: "
            f"{TRANSVERSION_COMPLEMENT}"
        ),
        "complement_used": TRANSVERSION_COMPLEMENT,
        "timestamp": timestamp,
        "seed": seed,
        "torch_version": torch.__version__,
        "n_families": len(families),
        "mutation_results": r_mutation,
    }

    out_path = out_dir / f"{model_name}_transversion.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    vol.commit()

    n_exc_nuc = r_mutation.get("n_exceeding_nuc_null", 0)
    n_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and "dinuc_null_95th" in v
    )
    n_exc_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and v.get("exceeds_dinuc_null")
    )

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] {model_name} transversion control COMPLETE.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}")
    print(f"  Nuc null: {n_exc_nuc} families exceed")
    print(f"  Dinuc null: {n_dinuc} tested, {n_exc_dinuc} exceed")
    print(f"  Saved: {out_path}")

    return {
        "model": model_name,
        "path": str(out_path),
        "mean_ratio": r_mutation.get("mean_best_ratio"),
        "n_exc_nuc": n_exc_nuc,
        "n_exc_dinuc": n_exc_dinuc,
    }


@app.function(
    image=dnabert2_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_dnabert2_transversion(seed: int = 42):
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
    print(f"[{timestamp}] Transversion control: dnabert2")
    print(f"  Patching COMPLEMENT -> {TRANSVERSION_COMPLEMENT}")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families

    import phases_1_to_5
    phases_1_to_5.COMPLEMENT = TRANSVERSION_COMPLEMENT
    from phases_1_to_5 import run_mutation_sensitivity

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = load_adapter("dnabert2")
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] dnabert2 loaded on {device}")

    out_dir = Path(f"/results/dnabert2_transversion_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity with transversions...")
    r_mutation = run_mutation_sensitivity(adapter, "dnabert2", families, device=device)

    output = {
        "model": "dnabert2",
        "experiment": "transversion_control",
        "description": (
            "Mutation sensitivity with transversion mutations (purine->pyrimidine) "
            "instead of Watson-Crick complements. COMPLEMENT dict patched to: "
            f"{TRANSVERSION_COMPLEMENT}"
        ),
        "complement_used": TRANSVERSION_COMPLEMENT,
        "timestamp": timestamp,
        "seed": seed,
        "torch_version": torch.__version__,
        "n_families": len(families),
        "mutation_results": r_mutation,
    }

    out_path = out_dir / "dnabert2_transversion.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    vol.commit()

    n_exc_nuc = r_mutation.get("n_exceeding_nuc_null", 0)
    n_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and "dinuc_null_95th" in v
    )
    n_exc_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and v.get("exceeds_dinuc_null")
    )

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] dnabert2 transversion control COMPLETE.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}")
    print(f"  Nuc null: {n_exc_nuc} families exceed")
    print(f"  Dinuc null: {n_dinuc} tested, {n_exc_dinuc} exceed")
    print(f"  Saved: {out_path}")

    return {
        "model": "dnabert2",
        "path": str(out_path),
        "mean_ratio": r_mutation.get("mean_best_ratio"),
        "n_exc_nuc": n_exc_nuc,
        "n_exc_dinuc": n_exc_dinuc,
    }


@app.function(
    image=transformers_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_transformers_transversion(model_name: str, seed: int = 42):
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
    print(f"[{timestamp}] Transversion control: {model_name}")
    print(f"  Patching COMPLEMENT -> {TRANSVERSION_COMPLEMENT}")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families

    import phases_1_to_5
    phases_1_to_5.COMPLEMENT = TRANSVERSION_COMPLEMENT
    from phases_1_to_5 import run_mutation_sensitivity

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    out_dir = Path(f"/results/{model_name}_transversion_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity with transversions...")
    r_mutation = run_mutation_sensitivity(adapter, model_name, families, device=device)

    output = {
        "model": model_name,
        "experiment": "transversion_control",
        "description": (
            "Mutation sensitivity with transversion mutations (purine->pyrimidine) "
            "instead of Watson-Crick complements. COMPLEMENT dict patched to: "
            f"{TRANSVERSION_COMPLEMENT}"
        ),
        "complement_used": TRANSVERSION_COMPLEMENT,
        "timestamp": timestamp,
        "seed": seed,
        "torch_version": torch.__version__,
        "n_families": len(families),
        "mutation_results": r_mutation,
    }

    out_path = out_dir / f"{model_name}_transversion.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    vol.commit()

    n_exc_nuc = r_mutation.get("n_exceeding_nuc_null", 0)
    n_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and "dinuc_null_95th" in v
    )
    n_exc_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and v.get("exceeds_dinuc_null")
    )

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] {model_name} transversion control COMPLETE.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}")
    print(f"  Nuc null: {n_exc_nuc} families exceed")
    print(f"  Dinuc null: {n_dinuc} tested, {n_exc_dinuc} exceed")
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
    all_models = MULTIMOL_MODELS + ["dnabert2"] + TRANSFORMERS_MODELS
    print(f"Launching transversion control for {len(all_models)} models")
    print(f"  Multimol (multimolecule image): {MULTIMOL_MODELS}")
    print(f"  DNABERT-2 (flash attn patched image)")
    print(f"  Transformers (plain transformers image): {TRANSFORMERS_MODELS}")
    print(f"  Transversion COMPLEMENT: {TRANSVERSION_COMPLEMENT}")
    print(f"  Seed: {seed}")
    print()

    handles = []

    for model_name in MULTIMOL_MODELS:
        fc = run_multimol_transversion.spawn(model_name=model_name, seed=seed)
        print(f"  Spawned {model_name} (multimol image) on A10G: {fc.object_id}")
        handles.append((model_name, fc))

    fc_dnabert2 = run_dnabert2_transversion.spawn(seed=seed)
    print(f"  Spawned dnabert2 (dnabert2 image) on A10G: {fc_dnabert2.object_id}")
    handles.append(("dnabert2", fc_dnabert2))

    for model_name in TRANSFORMERS_MODELS:
        fc = run_transformers_transversion.spawn(model_name=model_name, seed=seed)
        print(f"  Spawned {model_name} (transformers image) on A10G: {fc.object_id}")
        handles.append((model_name, fc))

    print(f"\n{len(handles)} containers launched. All models running in parallel.")
    print("Use 'modal app logs causal-rna-transversion-control' to monitor.")
