"""Modal wrapper: Phase 6 randomized-weight controls for DNABERT-2 and Evo.

These models don't use multimolecule, so they need transformers==4.47.1
(the version used in all prior successful runs).

Usage:
    modal run --detach scripts/modal_phase6_untrained_dnabert2_evo.py --model dnabert2
    modal run --detach scripts/modal_phase6_untrained_dnabert2_evo.py --model evo
"""

import modal

app = modal.App("causal-rna-phase6-untrained-dnabert2-evo")

base_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.5.1",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.67.1",
        "transformers==4.47.1",
        "matplotlib==3.9.3",
        "scikit-learn==1.6.0",
        "einops==0.8.0",
    )
    .add_local_dir("pretrained/evo_code/evo-1-8k-base", "/root/project/pretrained/evo_code/evo-1-8k-base")
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)

ADAPTER_MAP = {
    "dnabert2": "DNABERT2Adapter",
    "evo": "EvoAdapter",
}

LARGE_MODELS = {"evo"}


@app.function(
    image=base_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_untrained_model_a100(model_key: str, seed: int = 42):
    return _run_untrained(model_key, seed)


@app.function(
    image=base_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_untrained_model(model_key: str, seed: int = 42):
    return _run_untrained(model_key, seed)


def _run_untrained(model_key: str, seed: int = 42):
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
    adapter_class_name = ADAPTER_MAP[model_key]
    print(f"[{timestamp}] Phase 6 UNTRAINED {model_key} ({adapter_class_name})")

    import multi_model_audit
    from phase6_compensatory_mutation import load_rfam_families, run_phase6

    AdapterClass = getattr(multi_model_audit, adapter_class_name)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    adapter = AdapterClass()

    if model_key == "dnabert2":
        from transformers import AutoTokenizer, BertModel, BertConfig
        adapter.tokenizer = AutoTokenizer.from_pretrained(
            "zhihan1996/DNABERT-2-117M", trust_remote_code=True,
        )
        config = BertConfig.from_pretrained("zhihan1996/DNABERT-2-117M")
        adapter.model = BertModel(config)
        adapter.model.eval()
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Created DNABERT-2 from config (random init)")
    elif model_key == "evo":
        from transformers import AutoConfig, AutoModelForCausalLM
        config = AutoConfig.from_pretrained(
            "/root/project/pretrained/evo_code/evo-1-8k-base",
            trust_remote_code=True,
            local_files_only=True,
        )
        adapter.model = AutoModelForCausalLM.from_config(
            config, trust_remote_code=True, torch_dtype=torch.float16,
        )
        adapter.model.eval()
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Created Evo from config (random init)")
    else:
        adapter.load()
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {model_key}")
        with torch.no_grad():
            for name, param in adapter.model.named_parameters():
                if param.dim() >= 2:
                    torch.nn.init.xavier_normal_(param)
                else:
                    torch.nn.init.normal_(param, std=0.02)

    n_params = sum(p.numel() for p in adapter.model.parameters())
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Randomized model has {n_params:,} parameters")

    if device == "cuda":
        adapter.model = adapter.model.to(device)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running Phase 6 on untrained {model_key}...")
    results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=0,
    )

    output = {
        "model": f"{model_key}_untrained",
        "timestamp": timestamp,
        "phase": 6,
        "metric": "perturbation_specificity",
        "control_type": "randomized_weights",
        "description": (
            f"{adapter_class_name} with all weights randomized "
            f"(xavier_normal_ for matrices, normal_(std=0.02) for vectors). "
            f"Architecture preserved, learned weights destroyed."
        ),
        "seed": seed,
        "results": results,
    }

    out_dir = Path(f"/results/phase6_{model_key}_untrained_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_key}_untrained_phase6_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] UNTRAINED {model_key} complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families total: {results.get('families_total', 0)}")
    print(f"  Gate pass: {results.get('families_gate_pass', 'N/A')}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")
    print(f"  Saved to: {out_path}")

    return {"model": model_key, "path": str(out_path), "mean_ps": results.get("mean_best_ps")}


@app.local_entrypoint()
def main(model: str, seed: int = 42):
    if model not in ADAPTER_MAP:
        raise ValueError(f"Unknown model '{model}'. Choose from: {list(ADAPTER_MAP.keys())}")

    gpu_label = "A100" if model in LARGE_MODELS else "A10G"
    print(f"Launching Phase 6 randomized-weight control for {model} on {gpu_label}")

    if model in LARGE_MODELS:
        result = run_untrained_model_a100.remote(model_key=model, seed=seed)
    else:
        result = run_untrained_model.remote(model_key=model, seed=seed)

    print(f"Done! Result: {result}")
