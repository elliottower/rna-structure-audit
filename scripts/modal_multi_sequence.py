"""Modal wrapper: multi-sequence mutation sensitivity for within-family variance.

Runs mutation sensitivity on 3-5 alternative Rfam seed sequences per family
for three models (ERNIE-RNA, RiNALMo, NT v2). Measures within-family variance
of the stem/loop sensitivity ratio across different sequences from the same
RNA family.

Data: data/multi_sequence/{family}_multi.json files, each containing 3-5
sequences extracted from Rfam seed alignments.

Usage:
    cd /path/to/causal-rna
    modal run --detach scripts/modal_multi_sequence.py
"""

import modal

app = modal.App("causal-rna-multi-sequence")

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
    .add_local_dir("data/multi_sequence", "/root/project/data/multi_sequence")
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
    .add_local_dir("data/multi_sequence", "/root/project/data/multi_sequence")
)

vol = modal.Volume.from_name("causal-rna-multi-sequence", create_if_missing=True)


def _run_multiseq(model_name, seed):
    """Shared logic: load multi-sequence data, run mutation sensitivity, aggregate."""
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
    from tqdm import tqdm

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] Multi-sequence mutation sensitivity: {model_name}")

    from phase6_compensatory_mutation import load_adapter
    from phases_1_to_5 import run_mutation_sensitivity

    device = "cuda" if torch.cuda.is_available() else "cpu"
    np.random.seed(seed)

    # Load all multi-sequence family files
    multi_dir = Path("/root/project/data/multi_sequence")
    multi_files = sorted(multi_dir.glob("*_multi.json"))
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Found {len(multi_files)} multi-sequence family files")

    # Build mapping: family_name -> list of sequences
    family_data = {}
    for mf in multi_files:
        with open(mf) as fh:
            data = json.load(fh)
        family_data[data["name"]] = data

    # Flatten all sequences into synthetic families for the pipeline
    synthetic_families = []
    seq_to_family = {}  # maps synthetic name -> (original_family, seq_index)
    for fam_name, data in family_data.items():
        for i, seq_entry in enumerate(data["sequences"]):
            synth_name = f"{fam_name}_seq{i}"
            synthetic_families.append({
                "name": synth_name,
                "sequence": seq_entry["sequence"],
                "dot_bracket": seq_entry["dot_bracket"],
            })
            seq_to_family[synth_name] = (fam_name, i)

    total_seqs = len(synthetic_families)
    total_families = len(family_data)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {total_seqs} sequences across {total_families} families")

    # Load model
    adapter = load_adapter(model_name)
    adapter.load()
    if hasattr(adapter, "model") and adapter.model is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {model_name} loaded on {device}")

    # Run mutation sensitivity on all synthetic families at once
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Running mutation sensitivity on {total_seqs} sequences...")
    r_mutation = run_mutation_sensitivity(adapter, model_name, synthetic_families, device=device)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Mutation sensitivity complete")

    # Group results back by original family
    per_family = {}
    for fam_name, data in family_data.items():
        n_sequences = len(data["sequences"])
        ratios = []
        exceeds_nuc = []

        for i in range(n_sequences):
            synth_name = f"{fam_name}_seq{i}"
            seq_result = r_mutation.get("per_rna", {}).get(synth_name, {})
            if "best_ratio" in seq_result:
                ratios.append(seq_result["best_ratio"])
                exceeds_nuc.append(seq_result.get("exceeds_nuc_null", False))

        if ratios:
            mean_ratio = float(np.mean(ratios))
            sd_ratio = float(np.std(ratios, ddof=1)) if len(ratios) > 1 else 0.0
            cv_ratio = sd_ratio / mean_ratio if mean_ratio > 1e-10 else 0.0

            per_family[fam_name] = {
                "n_sequences": n_sequences,
                "n_scored": len(ratios),
                "per_sequence_ratios": ratios,
                "per_sequence_exceeds_nuc": exceeds_nuc,
                "mean_ratio": mean_ratio,
                "sd_ratio": sd_ratio,
                "cv_ratio": cv_ratio,
                "n_exceeds_nuc": sum(1 for e in exceeds_nuc if e),
            }
        else:
            per_family[fam_name] = {
                "n_sequences": n_sequences,
                "n_scored": 0,
                "skipped": "all sequences failed validation",
            }

    # Overall statistics across families
    all_cvs = [v["cv_ratio"] for v in per_family.values() if "cv_ratio" in v]
    all_sds = [v["sd_ratio"] for v in per_family.values() if "sd_ratio" in v]
    all_means = [v["mean_ratio"] for v in per_family.values() if "mean_ratio" in v]

    output = {
        "model": model_name,
        "experiment": "multi_sequence_mutation_sensitivity",
        "timestamp": timestamp,
        "seed": seed,
        "torch_version": torch.__version__,
        "n_families": total_families,
        "n_total_sequences": total_seqs,
        "per_family": per_family,
        "overall_mean_cv": float(np.mean(all_cvs)) if all_cvs else 0.0,
        "overall_median_cv": float(np.median(all_cvs)) if all_cvs else 0.0,
        "overall_mean_sd": float(np.mean(all_sds)) if all_sds else 0.0,
        "overall_mean_ratio": float(np.mean(all_means)) if all_means else 0.0,
        "n_families_scored": len(all_cvs),
    }

    # Save results
    out_dir = Path(f"/results/{model_name}_multiseq_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{model_name}_multiseq.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    vol.commit()

    # Print summary
    done_ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{done_ts}] {model_name} multi-sequence COMPLETE.")
    print(f"  Families scored: {len(all_cvs)}/{total_families}")
    print(f"  Overall mean ratio: {output['overall_mean_ratio']:.4f}")
    print(f"  Overall mean CV: {output['overall_mean_cv']:.4f}")
    print(f"  Overall mean SD: {output['overall_mean_sd']:.4f}")
    for fam_name, fam_result in sorted(per_family.items()):
        if "mean_ratio" in fam_result:
            print(f"    {fam_name}: mean={fam_result['mean_ratio']:.3f} "
                  f"SD={fam_result['sd_ratio']:.3f} "
                  f"CV={fam_result['cv_ratio']:.3f} "
                  f"({fam_result['n_exceeds_nuc']}/{fam_result['n_scored']} exceed nuc null)")
        else:
            print(f"    {fam_name}: {fam_result.get('skipped', 'no data')}")
    print(f"  Saved: {out_path}")

    return {
        "model": model_name,
        "path": str(out_path),
        "overall_mean_ratio": output["overall_mean_ratio"],
        "overall_mean_cv": output["overall_mean_cv"],
        "n_families_scored": len(all_cvs),
    }


@app.function(
    image=multimol_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_multiseq_ernierna(seed: int = 42):
    return _run_multiseq("ernierna", seed)


@app.function(
    image=multimol_image,
    gpu="A100",
    timeout=86400,
    volumes={"/results": vol},
)
def run_multiseq_rinalmo(seed: int = 42):
    return _run_multiseq("rinalmo", seed)


@app.function(
    image=transformers_image,
    gpu="A10G",
    timeout=86400,
    volumes={"/results": vol},
)
def run_multiseq_nt(seed: int = 42):
    return _run_multiseq("nt", seed)


@app.local_entrypoint()
def main(seed: int = 42):
    print("Multi-sequence mutation sensitivity: within-family variance")
    print("Models: ERNIE-RNA (A10G), RiNALMo (A100), NT v2 (A10G)")
    print(f"Seed: {seed}")

    fc_ernierna = run_multiseq_ernierna.spawn(seed=seed)
    print(f"  Spawned ernierna (multimol, A10G): {fc_ernierna.object_id}")

    fc_rinalmo = run_multiseq_rinalmo.spawn(seed=seed)
    print(f"  Spawned rinalmo (multimol, A100): {fc_rinalmo.object_id}")

    fc_nt = run_multiseq_nt.spawn(seed=seed)
    print(f"  Spawned nt (transformers, A10G): {fc_nt.object_id}")

    print(f"\n3 containers launched.")
    print("Use 'modal app logs causal-rna-multi-sequence' to monitor.")
