"""Re-run every rung on the repaired 47-family panel, one container per model.

Why this exists rather than a re-invocation of `modal_all_phases.py`: the
deposited results were computed with nulls seeded from a family's position in
the loaded panel, and `family_seed.family_rng` now seeds from the family name.
Every stored null therefore changes once, in every cell, so no cell of the
results table can be carried over from the old run and no cell can be compared
to it digit for digit. The panel itself also changed -- four annotations
repaired against Rfam seed alignments, five records withdrawn -- so this run
replaces the table rather than patching it. See DEVIATIONS.md, 2026-08-24.

Three images, because no single transformers version loads all ten models:
DNABERT-2's remote code predates the 4.29 attention refactor and needs its
flash-attn import patched out, and Caduceus dispatches to a mamba-ssm CUDA
kernel that has to be compiled with nvcc. The packages that *compute* the
statistics -- numpy, scipy, scikit-learn -- are pinned identically across all
three, so the split is in model loading only, and each cell records the stack it
was produced under.

Every output carries the commit it was produced at, the number of families
loaded, the names withdrawn, and a hash over the panel records themselves, so a
reader can tell which panel a number came from without trusting a filename.

The output directory is keyed on the model and not on a timestamp. A timestamped
directory means a restarted container resumes nothing, having written its
shards where the new container will not look.

Usage:
    modal run scripts/modal_repaired_panel.py --smoke-only
    modal run --detach scripts/modal_repaired_panel.py
    modal run --detach scripts/modal_repaired_panel.py --models rinalmo,ernierna
"""

import hashlib
import json
import subprocess
from pathlib import Path

import modal

app = modal.App("rna-repaired-panel")

REPO = Path(__file__).resolve().parent.parent

# Identical in all three images: these compute the statistics, and a table whose
# cells were produced under different numerics is not one table.
ANALYSIS_PACKAGES = (
    "numpy==1.26.4",
    "scipy==1.14.1",
    "scikit-learn==1.5.2",
    "matplotlib==3.9.2",
    "tqdm==4.66.5",
    "einops==0.8.0",
)

# huggingface-hub is deliberately unpinned: transformers 5.x and transformers
# 4.28 require incompatible ranges, and guessing a version is what cost the
# first launch of this script. The resolved version is recorded in every stamp.


def with_project(image):
    """Mount the repo and put it on the path.

    PYTHONPATH rather than a sys.path insert inside the container function, and
    no local checkpoint: RNA-FM falls back to a pinned HF revision when
    `pretrained/pytorch_model.bin` is absent, and it is absent from this repo.
    Local mounts are the last build step Modal permits, so this runs last.
    """
    return (
        image
        .env({"PYTHONPATH": "/root/project:/root/project/scripts"})
        .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
        .add_local_dir("scripts", "/root/project/scripts")
        .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
    )


# multimolecule is installed here but used by only four of the eight; it is
# harmless to the rest and keeps them on one image instead of two that differ by
# a single package.
main_image = with_project(
    modal.Image.debian_slim(python_version="3.11").pip_install(
        "torch==2.6.0",
        "transformers==5.14.1",
        "multimolecule==0.2.0",
        # multimolecule 0.2.0 imports torchmetrics through danling without
        # depending on it, so the container dies at adapter load without this.
        "torchmetrics==1.4.1",
        *ANALYSIS_PACKAGES,
    )
)

dnabert2_image = with_project(
    modal.Image.debian_slim(python_version="3.10")
    .pip_install("torch==2.4.0", "transformers==4.28.0", *ANALYSIS_PACKAGES)
    .add_local_file("scripts/patch_dnabert2_flash_attn.py",
                    "/root/patch_dnabert2_flash_attn.py", copy=True)
    .run_commands("python /root/patch_dnabert2_flash_attn.py")
)

caduceus_image = with_project(
    modal.Image.from_registry("nvidia/cuda:12.1.0-devel-ubuntu22.04", add_python="3.11")
    .apt_install("git")
    .pip_install(
        "torch==2.4.1",
        "transformers==4.44.2",
        *ANALYSIS_PACKAGES,
        "packaging==24.1", "ninja==1.11.1.1", "wheel==0.44.0", "setuptools==75.1.0",
    )
    # Both wheels compile against the installed torch, so they are built with a
    # GPU attached and without build isolation.
    .run_commands("pip install --no-build-isolation causal-conv1d==1.4.0", gpu="A10G")
    .run_commands("pip install --no-build-isolation mamba-ssm==2.2.4", gpu="A10G")
)

vol = modal.Volume.from_name("rna-repaired-panel-results", create_if_missing=True)

# Which image each model loads under, from the partition the existing scripts
# established: multimolecule for the four RNA language models, plain
# transformers for the four that load through trust_remote_code, and one image
# each for the two that need a build.
MAIN_MODELS = ["rnafm", "rinalmo", "utrlm", "ernierna",
               "splicebert", "nt", "hyenadna", "evo"]
AVAILABLE_MODELS = MAIN_MODELS + ["dnabert2", "caduceus"]
A100_MODELS = {"evo", "rinalmo"}
# Two models have no CPU path at all: the Evo adapter raises rather than load 7B
# parameters onto a CPU, and mamba-ssm has no CPU kernel. They smoke on a GPU.
GPU_SMOKE_MODELS = {"evo", "caduceus"}


def panel_stamp(families):
    """A hash over the records actually loaded, plus what was left out.

    Keyed on the fields the analysis reads. A filename can be copied onto the
    wrong file; this cannot.
    """
    payload = json.dumps(
        [[f["name"], f["sequence"], f["dot_bracket"]] for f in sorted(families, key=lambda f: f["name"])],
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _library_versions():
    """Every package that could move a number, as resolved in this container."""
    from importlib.metadata import PackageNotFoundError, version

    names = ["torch", "numpy", "scipy", "scikit-learn", "transformers",
             "multimolecule", "huggingface-hub", "mamba-ssm"]
    out = {}
    for name in names:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = None
    return out


def _run_model(model_name, commit, phase6_only):
    import os
    os.chdir("/root/project")

    from datetime import datetime, timezone

    import torch

    def now():
        return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    from family_checkpoint import FamilyCheckpoint
    from phase6_compensatory_mutation import (
        ADAPTER_OFFSETS,
        NON_CHARACTER_TOKENIZERS,
        load_adapter,
        load_rfam_families,
        run_phase6,
    )
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    families = load_rfam_families()
    withdrawn = sorted(
        json.loads(p.read_text())["name"]
        for p in Path("data/rfam_families").glob("*.json")
        if "excluded" in json.loads(p.read_text())
    )
    stamp = {
        "commit": commit,
        "panel_sha256": panel_stamp(families),
        "n_families": len(families),
        "withdrawn": withdrawn,
        "seeding": "family_seed.family_rng, derived from the family name",
        "deviation": "DEVIATIONS.md, 2026-08-24",
        "libraries": _library_versions(),
        "device": device,
    }
    print(f"[{now()}] {model_name}: {len(families)} families, panel "
          f"{stamp['panel_sha256'][:12]}, withdrawn {withdrawn}")

    out_dir = Path("/results") / model_name
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "stamp.json").write_text(json.dumps(stamp, indent=2) + "\n")
    vol.commit()

    adapter = load_adapter(model_name)
    adapter.load()
    if getattr(adapter, "model", None) is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    print(f"[{now()}] {model_name} loaded on {device}")

    def checkpoint(stage):
        return FamilyCheckpoint(out_dir / f"_shards_{stage}.json", commit=vol.commit)

    def save(name, payload):
        payload = {"model": model_name, "stamp": stamp,
                   "completed": now(), **payload}
        (out_dir / name).write_text(json.dumps(payload, indent=2, default=str))
        vol.commit()
        print(f"[{now()}] wrote {name}")

    if not phase6_only:
        print(f"[{now()}] {model_name}: mutation sensitivity")
        mutation = run_mutation_sensitivity(adapter, model_name, families,
                                            device=device,
                                            checkpoint=checkpoint("mutation"))
        print(f"[{now()}] {model_name}: attention-contact")
        attention = run_attention_contact(adapter, model_name, families, device=device,
                                          checkpoint=checkpoint("attention"))
        # Probing pools every family's embeddings before it fits, so there is no
        # per-family unit to checkpoint; it is saved the moment it returns.
        print(f"[{now()}] {model_name}: structure probing")
        probing = run_structure_probing(adapter, model_name, families, device=device)
        save(f"{model_name}_phases_1_to_5.json",
             {"experiment": "phases_1_to_5", "mutation_trained": mutation,
              "attention_trained": attention, "probing": probing})
        print(f"  mutation: mean_ratio={mutation.get('mean_best_ratio', 0):.4f}, "
              f"exceeding null {mutation.get('n_exceeding_nuc_null', 0)}"
              f"/{mutation.get('n_families_scored', 0)}")
        if probing.get("best_accuracy") is not None:
            print(f"  probing: {probing['best_accuracy']:.3f} at layer {probing.get('best_layer')}")

    print(f"[{now()}] {model_name}: phase 6 perturbation specificity")
    phase6 = run_phase6(adapter, families, device=device, compute_null=True,
                        offset=ADAPTER_OFFSETS.get(model_name, 0),
                        checkpoint=checkpoint("phase6"))
    save(f"{model_name}_phase6_ps.json",
         {"phase": 6, "metric": "perturbation_specificity",
          "preregistration": "PREREGISTRATION_PHASE6_V2.md",
          "offset": ADAPTER_OFFSETS.get(model_name, 0),
          "tokenizer_caveated": model_name in NON_CHARACTER_TOKENIZERS,
          "results": phase6})
    print(f"[{now()}] {model_name} COMPLETE. mean PS {phase6.get('mean_best_ps')}, "
          f"confirmatory {phase6.get('families_total')}, "
          f"exceeding null {phase6.get('families_exceeding_null_primary')}")
    return {"model": model_name, "mean_ps": phase6.get("mean_best_ps"),
            "panel_sha256": stamp["panel_sha256"]}


def _smoke(model_name):
    """Import, load the panel, load the model, score the six shortest families.

    Ten GPU containers that all die on the same missing package cost far more
    than one CPU container that finds it. This runs the same import path and the
    same adapter load as the real thing.
    """
    import os
    os.chdir("/root/project")

    from phase6_compensatory_mutation import load_adapter, load_rfam_families, run_phase6

    families = load_rfam_families()
    adapter = load_adapter(model_name)
    adapter.load()
    shortest = sorted(families, key=lambda f: len(f["sequence"]))[:6]
    scored = run_phase6(adapter, shortest, device="cpu", compute_null=False)
    ok = [n for n, r in scored["per_rna"].items() if not r.get("skipped")]
    print(f"{model_name}: {len(families)} families, adapter loaded, scored {ok}")
    print(f"  libraries: {_library_versions()}")
    if not ok:
        raise RuntimeError(f"{model_name} scored no family among the six shortest")
    return model_name


@app.function(image=main_image, timeout=3600)
def smoke_main(model_name: str):
    return _smoke(model_name)


@app.function(image=dnabert2_image, timeout=3600)
def smoke_dnabert2(model_name: str = "dnabert2"):
    return _smoke(model_name)


@app.function(image=main_image, gpu="A10G", timeout=3600)
def smoke_main_gpu(model_name: str):
    """For the models on the main image that refuse to load without a GPU."""
    return _smoke(model_name)


@app.function(image=caduceus_image, gpu="A10G", timeout=3600)
def smoke_caduceus(model_name: str = "caduceus"):
    """On a GPU, because mamba-ssm has no CPU kernel to fall back to."""
    return _smoke(model_name)


@app.function(image=main_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_main_a10g(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


@app.function(image=main_image, gpu="A100", timeout=86400, volumes={"/results": vol})
def run_main_a100(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


@app.function(image=dnabert2_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_dnabert2(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


@app.function(image=caduceus_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_caduceus(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


def _route(model_name):
    """The (function, label) each model runs under."""
    if model_name == "caduceus":
        return run_caduceus, "caduceus/A10G"
    if model_name == "dnabert2":
        return run_dnabert2, "dnabert2/A10G"
    if model_name in A100_MODELS:
        return run_main_a100, "main/A100"
    return run_main_a10g, "main/A10G"


@app.local_entrypoint()
def main(models: str = "", phase6_only: bool = False, smoke_only: bool = False):
    commit = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                            capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain"],
                           capture_output=True, text=True, check=True).stdout.strip()
    if dirty:
        raise SystemExit(
            "Refusing to launch from a dirty tree: the stamp would name a commit "
            "that does not contain the code being run.\n" + dirty)

    requested = [m.strip() for m in models.split(",") if m.strip()] or AVAILABLE_MODELS
    unknown = [m for m in requested if m not in AVAILABLE_MODELS]
    if unknown:
        raise SystemExit(f"Unknown models: {unknown}")

    if smoke_only:
        print(f"Smoke check, {len(requested)} models: {requested}")
        # return_exceptions, so one broken adapter reports itself alongside the
        # nine that work instead of hiding them behind its own traceback.
        on_cpu = [m for m in requested
                  if m in MAIN_MODELS and m not in GPU_SMOKE_MODELS]
        outcomes = dict(zip(on_cpu, smoke_main.map(on_cpu, return_exceptions=True)))
        for model_name in requested:
            if model_name in on_cpu:
                continue
            if model_name == "caduceus":
                fn = smoke_caduceus
            elif model_name in GPU_SMOKE_MODELS:
                fn = smoke_main_gpu
            else:
                fn = smoke_dnabert2
            try:
                outcomes[model_name] = fn.remote(model_name)
            except Exception as exc:  # noqa: BLE001 -- reported, not swallowed
                outcomes[model_name] = exc

        print()
        failed = []
        for model_name in requested:
            result = outcomes.get(model_name)
            if isinstance(result, BaseException):
                failed.append(model_name)
                print(f"  FAIL {model_name:12s} {type(result).__name__}: {result}")
            else:
                print(f"  ok   {model_name}")
        if failed:
            raise SystemExit(f"\n{len(failed)} of {len(requested)} failed: {failed}")
        print(f"\nAll {len(requested)} models load and score.")
        return

    print(f"Repaired-panel re-run at {commit[:12]}, {len(requested)} models: {requested}")
    for model_name in requested:
        fn, label = _route(model_name)
        handle = fn.spawn(model_name=model_name, commit=commit, phase6_only=phase6_only)
        print(f"  {model_name:12s} {label:14s} {handle.object_id}")
    print("\nmodal app logs rna-repaired-panel")
