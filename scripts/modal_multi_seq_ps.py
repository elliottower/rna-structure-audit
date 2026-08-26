"""Modal wrapper: multi-sequence perturbation specificity from Rfam seed alignments.

Runs PS on all valid sequences from all 46 Rfam seed alignments for a given
model. One container per model. Results saved to Modal volume.

Usage:
    modal run --detach scripts/modal_multi_seq_ps.py --model rinalmo
    modal run --detach scripts/modal_multi_seq_ps.py --model ernierna
    modal run --detach scripts/modal_multi_seq_ps.py --model rinalmo --model ernierna
"""

import modal

app = modal.App("rna-multi-seq-ps")

base_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch>=2.2.0",
        "numpy",
        "scipy",
        "tqdm",
        "transformers>=4.40.0,<5.0",
        "multimolecule",
        "matplotlib",
        "scikit-learn",
        "einops",
    )
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
    .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    .add_local_dir("data/rfam_seeds", "/root/project/data/rfam_seeds")
    .add_local_file("pretrained/pytorch_model.bin", "/root/project/pretrained/pytorch_model.bin")
)

vol = modal.Volume.from_name("rna-multi-seq-ps-results", create_if_missing=True)

GPU_MAP = {
    "rinalmo": "A100",
    "ernierna": "A10G",
    "caduceus": "A10G",
    "evo": "A100",
}


@app.function(
    image=base_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_a100(model_name: str, seed: int = 42, compute_null: bool = True,
                   max_seqs_per_family: int = 0):
    return _run(model_name, seed, compute_null, max_seqs_per_family)


@app.function(
    image=base_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_model_a10g(model_name: str, seed: int = 42, compute_null: bool = True,
                   max_seqs_per_family: int = 0):
    return _run(model_name, seed, compute_null, max_seqs_per_family)


def _run(model_name, seed, compute_null, max_seqs_per_family):
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

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{ts}] Multi-seq PS: starting {model_name}")
    print(f"  seed={seed}, compute_null={compute_null}, max_seqs_per_family={max_seqs_per_family or 'all'}")

    np.random.seed(seed)

    from multi_seq_ps import run_multi_seq_ps
    from phase6_compensatory_mutation import ADAPTER_OFFSETS, load_adapter

    device = "cuda" if torch.cuda.is_available() else "cpu"
    offset = ADAPTER_OFFSETS.get(model_name, 0)

    print(f"[{ts}] Loading {model_name} (device={device})...")
    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)

    print(f"[{datetime.now().strftime('%Y%m%d_%H%M%S')}] Running multi-seq PS...")
    results = run_multi_seq_ps(
        adapter,
        seed_dir="/root/project/data/rfam_seeds",
        device=device,
        offset=offset,
        compute_null=compute_null,
        max_seqs_per_family=max_seqs_per_family or None,
    )

    output = {
        "model": model_name,
        "timestamp": ts,
        "experiment": "multi_sequence_perturbation_specificity",
        "source": "rfam_14.10_seed_alignments",
        "seed": seed,
        "offset": offset,
        "compute_null": compute_null,
        "max_seqs_per_family": max_seqs_per_family or "all",
        "summary": results["summary"],
        "per_family": results["per_family"],
    }

    out_dir = Path(f"/results/multi_seq_ps_{model_name}_{ts}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_name}_multi_seq_ps.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    vol.commit()

    done_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary = results["summary"]
    print(f"\n[{done_ts}] {model_name} COMPLETE")
    print(f"  Total sequences processed: {summary['n_total_sequences']}")
    print(f"  Families with data: {summary['n_families_with_data']}/{summary['n_families']}")
    if "global_mean_ps" in summary:
        print(f"  Global mean PS: {summary['global_mean_ps']:.6f}")
        print(f"  Global SE: {summary['global_se_ps']:.6f}")
    if "family_mean_of_means" in summary:
        print(f"  Family mean-of-means PS: {summary['family_mean_of_means']:.6f}")
        print(f"  Family SE-of-means: {summary['family_se_of_means']:.6f}")
    print(f"  Saved to {out_path}")

    return {"model": model_name, "path": str(out_path), "summary": summary}


@app.local_entrypoint()
def main(
    models: str = "rinalmo",
    seed: int = 42,
    no_null: bool = False,
    max_seqs: int = 0,
):
    model_list = [m.strip() for m in models.split(",")]
    print(f"Launching multi-seq PS for models: {model_list}")
    print(f"  seed={seed}, compute_null={not no_null}, max_seqs_per_family={max_seqs or 'all'}")

    handles = []
    for m in model_list:
        gpu = GPU_MAP.get(m, "A10G")
        run_fn = run_model_a100 if gpu == "A100" else run_model_a10g
        fc = run_fn.spawn(
            model_name=m,
            seed=seed,
            compute_null=not no_null,
            max_seqs_per_family=max_seqs,
        )
        print(f"  Spawned {m} on {gpu}: {fc.object_id}")
        handles.append((m, fc))

    print(f"\n{len(handles)} containers launched.")
    print("Monitor: modal app logs rna-multi-seq-ps")
