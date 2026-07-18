"""Modal wrapper: run Phase 1-5 + Phase 6 per model in one container.

One container per model. Loads model once, runs all phases, saves results.

Usage:
    modal run --detach scripts/modal_all_phases.py
    modal run --detach scripts/modal_all_phases.py --models ernierna,rinalmo
    modal run --detach scripts/modal_all_phases.py --models nt,dnabert2 --allow-non-character
    modal run --detach scripts/modal_all_phases.py --phase6-only --models rnafm
"""

import modal

app = modal.App("causal-rna-all-phases")

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

vol = modal.Volume.from_name("causal-rna-all-phases-results", create_if_missing=True)

AVAILABLE_MODELS = [
    "rnafm", "rinalmo", "utrlm", "ernierna", "splicebert",
    "nt", "hyenadna", "caduceus", "evo", "dnabert2",
]
CHARACTER_MODELS = [m for m in AVAILABLE_MODELS if m not in {"nt", "dnabert2"}]

GPU_OVERRIDES = {
    "evo": "A100",
    "rinalmo": "A100",
}
SKIP_MODELS = {"caduceus"}


@app.function(
    image=base_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_a10g(model_name: str, seed: int = 42, phase6_only: bool = False):
    return _run_all_phases(model_name, seed, phase6_only)


@app.function(
    image=base_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_a100(model_name: str, seed: int = 42, phase6_only: bool = False):
    return _run_all_phases(model_name, seed, phase6_only)


def _run_all_phases(model_name, seed, phase6_only):
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
    print(f"[{timestamp}] All phases: starting {model_name}")

    from phase6_compensatory_mutation import (
        ADAPTER_OFFSETS,
        NON_CHARACTER_TOKENIZERS,
        load_adapter,
        load_rfam_families,
        run_phase6,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    offset = ADAPTER_OFFSETS.get(model_name, 0)
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{timestamp}] Loaded {len(families)} families")

    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] {model_name} loaded on {device}")

    out_dir = Path(f"/results/all_phases_{model_name}_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1-5: mutation sensitivity, attention-contact, structure probing
    if not phase6_only:
        from phases_1_to_5 import (
            run_mutation_sensitivity,
            run_attention_contact,
            run_structure_probing,
        )

        print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running Phase 1-5: mutation sensitivity...")
        r_mutation = run_mutation_sensitivity(adapter, model_name, families, device=device)

        print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running Phase 1-5: attention contact...")
        r_attention = run_attention_contact(adapter, model_name, families, device=device)

        print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running Phase 1-5: structure probing...")
        r_probing = run_structure_probing(adapter, model_name, families, device=device)

        phase15_output = {
            "model": model_name,
            "timestamp": timestamp,
            "torch_version": torch.__version__,
            "mutation_trained": r_mutation,
            "attention_trained": r_attention,
            "probing": r_probing,
        }

        p15_path = out_dir / f"{model_name}_phases_1_to_5.json"
        with open(p15_path, "w") as f:
            json.dump(phase15_output, f, indent=2, default=str)
        vol.commit()
        print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Phase 1-5 saved: {p15_path}")
        print(f"  Mutation: mean_ratio={r_mutation.get('mean_best_ratio', 0):.4f}, "
              f"exceed_null={r_mutation.get('n_exceeding_nuc_null', 0)}/{r_mutation.get('n_families_scored', 0)}")
        if r_attention.get("skipped"):
            print(f"  Attention: {r_attention['skipped']}")
        if "best_accuracy" in r_probing:
            print(f"  Probing: best_accuracy={r_probing['best_accuracy']:.3f} at layer {r_probing.get('best_layer')}")

    # Phase 6: perturbation specificity
    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running Phase 6: perturbation specificity...")
    p6_results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=offset,
    )

    p6_output = {
        "model": model_name,
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "preregistration": "PREREGISTRATION_PHASE6_V2.md",
        "offset": offset,
        "seed": seed,
        "tokenizer_caveated": model_name in NON_CHARACTER_TOKENIZERS,
        "results": p6_results,
    }

    p6_path = out_dir / f"{model_name}_phase6_ps.json"
    with open(p6_path, "w") as f:
        json.dump(p6_output, f, indent=2)
    vol.commit()

    done_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{done_ts}] {model_name} ALL PHASES complete.")
    print(f"  Phase 6 mean PS: {p6_results.get('mean_best_ps', 'N/A')}")
    print(f"  Families (confirmatory): {p6_results.get('families_total', 0)}")
    print(f"  Exceeding null (primary): {p6_results.get('families_exceeding_null_primary', 'N/A')}")

    return {
        "model": model_name,
        "phase6_path": str(p6_path),
        "phase15_path": str(out_dir / f"{model_name}_phases_1_to_5.json") if not phase6_only else None,
        "mean_ps": p6_results.get("mean_best_ps"),
    }


@app.local_entrypoint()
def main(
    models: str = "",
    seed: int = 42,
    allow_non_character: bool = False,
    phase6_only: bool = False,
):
    if not models:
        model_list = list(CHARACTER_MODELS)
    else:
        model_list = [m.strip() for m in models.split(",")]

    print(f"Launching ALL PHASES for {len(model_list)} models: {model_list}")
    print(f"Seed: {seed}, phase6_only: {phase6_only}")

    handles = []
    for model_name in model_list:
        if model_name not in AVAILABLE_MODELS:
            print(f"  SKIP {model_name}: unknown model")
            continue
        if model_name in {"nt", "dnabert2"} and not allow_non_character:
            print(f"  SKIP {model_name}: non-character tokenizer (use --allow-non-character)")
            continue
        if model_name in SKIP_MODELS:
            print(f"  SKIP {model_name}: requires mamba-ssm CUDA build")
            continue

        if GPU_OVERRIDES.get(model_name) == "A100":
            run_fn = run_model_a100
            gpu = "A100"
        else:
            run_fn = run_model_a10g
            gpu = "A10G"

        fc = run_fn.spawn(model_name=model_name, seed=seed, phase6_only=phase6_only)
        print(f"  Spawned {model_name} on {gpu}: {fc.object_id}")
        handles.append((model_name, fc))

    print(f"\n{len(handles)} containers launched. Use 'modal app logs causal-rna-all-phases' to monitor.")
