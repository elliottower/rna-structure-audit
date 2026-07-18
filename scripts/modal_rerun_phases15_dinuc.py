"""Modal wrapper: re-run Phase 1-5 for models missing dinucleotide null data.

Three models (ERNIE-RNA, SpliceBERT, DNABERT-2) were run on July 13 before
the dinucleotide null was added to the Phase 1-5 pipeline. This re-runs them
to get the dinucleotide null results for Rung 2.

Also runs ERNIE-RNA with randomized weights (untrained control) through
Phase 1-5 for Table 1 completeness.

Usage:
    cd /path/to/causal-rna
    modal run --detach scripts/modal_rerun_phases15_dinuc.py
"""

import modal

app = modal.App("causal-rna-rerun-phases15-dinuc")

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

vol = modal.Volume.from_name("causal-rna-rerun-dinuc", create_if_missing=True)

MULTIMOL_MODELS = ["ernierna", "splicebert"]


@app.function(
    image=multimol_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_phases15(model_name: str, seed: int = 42):
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
    print(f"[{timestamp}] Phase 1-5 re-run (dinuc null): {model_name}")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    out_dir = Path(f"/results/{model_name}_phases15_dinuc_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity...")
    r_mutation = run_mutation_sensitivity(adapter, model_name, families, device=device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running attention contact...")
    r_attention = run_attention_contact(adapter, model_name, families, device=device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running structure probing...")
    r_probing = run_structure_probing(adapter, model_name, families, device=device)

    output = {
        "model": model_name,
        "experiment": "phases_1_to_5_with_dinuc_null",
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
    n_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and "dinuc_null_95th" in v
    )
    n_exc_dinuc = sum(
        1 for v in r_mutation.get("per_rna", {}).values()
        if isinstance(v, dict) and v.get("exceeds_dinuc_null")
    )

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] {model_name} Phase 1-5 COMPLETE.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}")
    print(f"  Nuc null: {n_exc_nuc} families exceed")
    print(f"  Dinuc null: {n_dinuc} tested, {n_exc_dinuc} exceed")
    if r_attention.get("skipped"):
        print(f"  Attention: {r_attention['skipped']}")
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


@app.function(
    image=dnabert2_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_dnabert2_phases15(seed: int = 42):
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
    print(f"[{timestamp}] Phase 1-5 re-run (dinuc null): dnabert2")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = load_adapter("dnabert2")
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] dnabert2 loaded on {device}")

    out_dir = Path(f"/results/dnabert2_phases15_dinuc_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity...")
    r_mutation = run_mutation_sensitivity(adapter, "dnabert2", families, device=device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running attention contact...")
    r_attention = run_attention_contact(adapter, "dnabert2", families, device=device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running structure probing...")
    r_probing = run_structure_probing(adapter, "dnabert2", families, device=device)

    output = {
        "model": "dnabert2",
        "experiment": "phases_1_to_5_with_dinuc_null",
        "timestamp": timestamp,
        "seed": seed,
        "torch_version": torch.__version__,
        "n_families": len(families),
        "mutation_trained": r_mutation,
        "attention_trained": r_attention,
        "probing": r_probing,
    }

    out_path = out_dir / "dnabert2_phases_1_to_5.json"
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
    print(f"[{done_ts}] dnabert2 Phase 1-5 COMPLETE.")
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
    image=multimol_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_untrained_ernierna_phases15(seed: int = 42):
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
    print(f"[{timestamp}] Phase 1-5 UNTRAINED ERNIE-RNA control")

    from multi_model_audit import ERNIERNAAdapter
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

    adapter = ERNIERNAAdapter()
    adapter.load()

    n_params = sum(p.numel() for p in adapter.model.parameters())
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded trained ERNIE-RNA ({n_params:,} params)")

    with torch.no_grad():
        for name, param in adapter.model.named_parameters():
            if param.dim() >= 2:
                torch.nn.init.xavier_normal_(param)
            else:
                torch.nn.init.normal_(param, std=0.02)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Randomized all {n_params:,} parameters")

    if device == "cuda":
        adapter.model = adapter.model.to(device)

    out_dir = Path(f"/results/ernierna_untrained_phases15_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    model_name = "ernierna_untrained"

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity...")
    r_mutation = run_mutation_sensitivity(adapter, "ernierna", families, device=device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running attention contact...")
    r_attention = run_attention_contact(adapter, "ernierna", families, device=device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running structure probing...")
    r_probing = run_structure_probing(adapter, "ernierna", families, device=device)

    output = {
        "model": model_name,
        "experiment": "phases_1_to_5",
        "control_type": "randomized_weights",
        "description": "ERNIE-RNA with all weights randomized (xavier_normal_ for matrices, normal_(std=0.02) for vectors). Architecture preserved, learned weights destroyed.",
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

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] UNTRAINED ERNIE-RNA Phase 1-5 COMPLETE.")
    print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}")
    print(f"  Saved: {out_path}")

    return {
        "model": model_name,
        "path": str(out_path),
        "mean_ratio": r_mutation.get("mean_best_ratio"),
    }


@app.local_entrypoint()
def main(seed: int = 42, dnabert2_only: bool = False):
    if dnabert2_only:
        print(f"Launching DNABERT-2 only (flash attn patched at build time)")
        fc = run_dnabert2_phases15.spawn(seed=seed)
        print(f"  Spawned dnabert2 on A10G: {fc.object_id}")
        print("1 container launched.")
    else:
        print(f"Launching Phase 1-5 re-runs with corrected family data + dinuc null")
        print(f"Models: {MULTIMOL_MODELS} + dnabert2 + ERNIE-RNA untrained")
        print(f"Seed: {seed}")

        handles = []
        for model_name in MULTIMOL_MODELS:
            fc = run_phases15.spawn(model_name=model_name, seed=seed)
            print(f"  Spawned {model_name} (multimol image) on A10G: {fc.object_id}")
            handles.append((model_name, fc))

        fc_dnabert2 = run_dnabert2_phases15.spawn(seed=seed)
        print(f"  Spawned dnabert2 (separate image) on A10G: {fc_dnabert2.object_id}")
        handles.append(("dnabert2", fc_dnabert2))

        fc_untrained = run_untrained_ernierna_phases15.spawn(seed=seed)
        print(f"  Spawned ernierna_untrained on A10G: {fc_untrained.object_id}")
        handles.append(("ernierna_untrained", fc_untrained))

        print(f"\n{len(handles)} containers launched.")
    print("Use 'modal app logs causal-rna-rerun-phases15-dinuc' to monitor.")
