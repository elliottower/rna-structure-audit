"""Modal wrapper: Phase 6 for DNABERT-2.

DNABERT-2's custom flash_attn_triton.py uses tl.dot(q, k, trans_b=True),
which was removed in all pip-installable triton versions. The model code
(bert_layers.py) has a fallback: if flash_attn_qkvpacked_func is None,
it uses standard PyTorch attention. We trigger this fallback by blocking
triton import before loading the model.

Usage:
    modal run --detach scripts/modal_phase6_dnabert2.py
    modal run --detach scripts/modal_phase6_dnabert2.py --synthetic
"""

import modal

app = modal.App("causal-rna-phase6-dnabert2")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.6.0",
        "numpy==2.1.3",
        "scipy==1.14.1",
        "tqdm==4.67.1",
        "transformers==4.49.0",
        "tokenizers==0.21.0",
        "huggingface-hub==0.27.1",
        "safetensors==0.5.3",
        "matplotlib==3.9.3",
        "scikit-learn==1.6.0",
        "einops==0.8.0",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py", copy=True)
    .add_local_dir("scripts", "/root/project/scripts", copy=True)
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families", copy=True)
)

vol = modal.Volume.from_name("causal-rna-phase6-results", create_if_missing=True)


@app.function(image=image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_dnabert2(synthetic: bool = False, n_synthetic: int = 5, seed: int = 42):
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
    mode = "SYNTHETIC" if synthetic else "NATURAL"
    print(f"[{timestamp}] Phase 6 {mode}: dnabert2")
    print(f"  torch={torch.__version__}, cuda={torch.cuda.is_available()}")

    # Block triton AFTER torch is loaded but BEFORE DNABERT-2 model init.
    # DNABERT-2's bert_layers.py: `from .flash_attn_triton import flash_attn_qkvpacked_func`
    # flash_attn_triton.py does `import triton` which uses removed API (trans_b).
    # Setting sys.modules["triton"] = None makes `import triton` raise ImportError.
    # bert_layers.py catches this with bare `except:`, sets flash_attn_qkvpacked_func = None.
    # Then: `if self.p_dropout or flash_attn_qkvpacked_func is None:` → True → PyTorch fallback.
    _triton_backup = sys.modules.pop("triton", None)
    _triton_lang_backup = sys.modules.pop("triton.language", None)
    sys.modules["triton"] = None
    sys.modules["triton.language"] = None
    print(f"  triton blocked to force PyTorch attention fallback")

    from phase6_compensatory_mutation import (
        ADAPTER_OFFSETS,
        load_rfam_families,
        run_phase6,
    )
    from transformers import AutoModel, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    offset = ADAPTER_OFFSETS.get("dnabert2", 0)
    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    families = load_rfam_families()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded {len(families)} families")

    if synthetic:
        from modal_phase6_synthetic_covariation import generate_synthetic_families
        families = generate_synthetic_families(families, n_synthetic, rng)
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Generated {len(families)} synthetic sequences")

    from multi_model_audit import DNABERT2Adapter
    adapter = DNABERT2Adapter()
    adapter.tokenizer = AutoTokenizer.from_pretrained(
        adapter.MODEL_ID, trust_remote_code=True,
    )
    adapter.model = AutoModel.from_pretrained(
        adapter.MODEL_ID, trust_remote_code=True,
    )
    adapter.model.eval()

    # Restore triton now that the model is loaded (torch.compile may need it)
    sys.modules.pop("triton", None)
    sys.modules.pop("triton.language", None)
    if _triton_backup is not None:
        sys.modules["triton"] = _triton_backup
    if _triton_lang_backup is not None:
        sys.modules["triton.language"] = _triton_lang_backup

    if device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] dnabert2 loaded on {device}")

    # Sanity check: run one forward pass to verify standard attention works
    test_tokens = adapter.tokenize("ACGTACGT").to(device)
    test_out = adapter.get_all_layer_embeddings(test_tokens)
    assert len(test_out) > 0 and test_out[0].shape[-1] == 768, f"Bad output shape: {[t.shape for t in test_out]}"
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Sanity check passed: output shape {test_out[0].shape}")

    results = run_phase6(
        adapter, families, device=device,
        compute_null=True, offset=offset,
    )

    prefix = "phase6_synthetic" if synthetic else "phase6"
    out_dir = Path(f"/results/{prefix}_dnabert2_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)

    output = {
        "model": "dnabert2",
        "experiment": "synthetic_covariation_control" if synthetic else "phase6_natural",
        "timestamp": timestamp,
        "results": results,
    }
    if synthetic:
        output["n_synthetic_per_family"] = n_synthetic
        output["seed"] = seed

    out_path = out_dir / f"dnabert2_{'synthetic_covariation' if synthetic else 'phase6_ps'}.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    vol.commit()

    done = datetime.now().strftime("%H:%M:%S")
    print(f"[{done}] dnabert2 {mode} complete.")
    print(f"  Mean PS: {results.get('mean_best_ps', 'N/A')}")
    print(f"  Families: {results.get('families_total', 0)}")
    print(f"  Exceeding null: {results.get('families_exceeding_null_primary', 'N/A')}")

    return {"model": "dnabert2", "mode": mode, "path": str(out_path)}


@app.local_entrypoint()
def main(
    synthetic: bool = False,
    n_synthetic: int = 5,
    seed: int = 42,
):
    mode = "SYNTHETIC" if synthetic else "NATURAL"
    print(f"Phase 6 {mode}: dnabert2 on A10G")
    fc = run_dnabert2.spawn(synthetic=synthetic, n_synthetic=n_synthetic, seed=seed)
    print(f"  Spawned: {fc.object_id}")
    print("Use 'modal app logs causal-rna-phase6-dnabert2' to monitor.")
