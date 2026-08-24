"""Re-run every rung on the repaired 47-family panel, one container per model.

Why this exists rather than a re-invocation of `modal_all_phases.py`: the
deposited results were computed with nulls seeded from a family's position in
the loaded panel, and `family_seed.family_rng` now seeds from the family name.
Every stored null therefore changes once, in every cell, so no cell of the
results table can be carried over from the old run and no cell can be compared
to it digit for digit. The panel itself also changed -- four annotations
repaired against Rfam seed alignments, five records withdrawn -- so this run
replaces the table rather than patching it. See DEVIATIONS.md, 2026-08-24.

Every image is copied from the script that produced the deposited results for
the models it serves, pin for pin. Nothing here was derived by trying versions
until one imported: `modal_all_phases.py` left transformers and multimolecule
unpinned, so what it resolved to in July is unrecoverable, and every attempt to
re-derive that resolution has failed on a different package.

Four stacks, because no single one loads all ten models, and each is the newest
that works for its models rather than the oldest that works for any:

    torch 2.6.0  / transformers 5.14.1 / multimolecule 0.2.0
        RiNALMo, UTR-LM, ERNIE-RNA, SpliceBERT
        from modal_ernierna_ablation.py and modal_multiseq_rinalmo.py
    torch 2.4.1  / transformers 4.44.2
        RNA-FM, NT, HyenaDNA -- NT's remote code calls
        find_pruneable_heads_and_indices, dropped in transformers 5
    torch 2.4.0  / transformers 4.28.0 (python 3.10)
        DNABERT-2, whose remote code predates the 4.29 attention refactor
        from modal_rerun_phases15_dinuc.py
    torch 2.1.2  / transformers 4.49.0 / flash-attn 2.5.8
        Evo alone. Its remote code hard-requires flash-attn, which builds only
        against this torch and Triton 2.1.0. Confined to one container: it is
        old enough that transformers refuses torch.load under CVE-2025-32434,
        and too old for the aten.rms_norm that multimolecule 0.2.0 calls.

Caduceus shares the 2.4.1 pins but needs its own image, because mamba-ssm
compiles against the installed torch with a GPU attached.

The output directory is keyed on the model and not on a timestamp. A timestamped
directory means a restarted container resumes nothing, having written its
shards where the new container will not look.

Every output carries the commit it was produced at, the number of families
loaded, the names withdrawn, and a hash over the panel records themselves, so a
reader can tell which panel a number came from without trusting a filename.

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


MULTIMOL_PACKAGES = ("numpy==1.26.4", "scipy==1.14.1", "tqdm==4.66.5",
                     "matplotlib==3.9.2", "scikit-learn==1.5.2", "einops==0.8.0")

# From `modal_ernierna_ablation.py` and `modal_multiseq_rinalmo.py`, plus
# torchmetrics: multimolecule 0.2.0 imports it through danling without depending
# on it, so the container dies at adapter load without it.
multimol_image = with_project(
    modal.Image.debian_slim(python_version="3.11").pip_install(
        "torch==2.6.0",
        "transformers==5.14.1",
        "multimolecule==0.2.0",
        "torchmetrics==1.4.1",
        *MULTIMOL_PACKAGES,
    )
)

# The three models that load through plain transformers with trust_remote_code.
# NT is the binding constraint: its remote code imports
# find_pruneable_heads_and_indices, which transformers 5 removed.
legacy_image = with_project(
    modal.Image.debian_slim(python_version="3.11").pip_install(
        "torch==2.4.1",
        "transformers==4.44.2",
        *MULTIMOL_PACKAGES,
    )
)

# From `modal_rerun_phases15_dinuc.py`, which produced
# results/dnabert2_phases15_dinuc.json. The patch runs at build time and
# pre-downloads the model, so the cached remote code is rewritten before any
# container imports it.
dnabert2_image = with_project(
    modal.Image.debian_slim(python_version="3.10")
    .pip_install("torch==2.4.0", "transformers==4.28.0", *MULTIMOL_PACKAGES)
    .add_local_file("scripts/patch_dnabert2_flash_attn.py",
                    "/root/patch_dnabert2_flash_attn.py", copy=True)
    .run_commands("python /root/patch_dnabert2_flash_attn.py")
)

# From `modal_phase6_flashattn_models.py`, which produced
# results/phase6_all/phase6_evo_20260715_073218. Evo only.
evo_image = with_project(
    modal.Image.from_registry("nvidia/cuda:12.1.1-devel-ubuntu22.04", add_python="3.11")
    .apt_install("git")
    .pip_install(
        "torch==2.1.2", "triton==2.1.0",
        "packaging", "ninja", "wheel", "setuptools",
        "numpy==1.26.4", "scipy==1.13.1", "tqdm==4.66.4",
        "transformers==4.49.0", "multimolecule==0.1.0",
        "matplotlib==3.9.0", "scikit-learn==1.5.0", "einops==0.8.0",
    )
    .pip_install("flash-attn==2.5.8", extra_options="--no-build-isolation", gpu="A10G")
)

# From `modal_caduceus_phases.py`. Both wheels compile against the installed
# torch, so they are built with a GPU attached and without build isolation.
caduceus_image = with_project(
    modal.Image.from_registry("nvidia/cuda:12.1.0-devel-ubuntu22.04", add_python="3.11")
    .apt_install("git")
    .pip_install(
        "torch==2.4.1", "transformers==4.44.2", "multimolecule==0.2.0",
        *MULTIMOL_PACKAGES,
        "packaging==24.1", "ninja==1.11.1.1", "wheel==0.44.0", "setuptools==75.1.0",
    )
    .run_commands("pip install --no-build-isolation causal-conv1d==1.4.0", gpu="A10G")
    .run_commands("pip install --no-build-isolation mamba-ssm==2.2.4", gpu="A10G")
)

vol = modal.Volume.from_name("rna-repaired-panel-results", create_if_missing=True)

# Which image each model runs under. Evo and DNABERT-2 and Caduceus are one
# model each, so they are named by the routing rather than listed.
MULTIMOL_MODELS = ["rinalmo", "utrlm", "ernierna", "splicebert"]
LEGACY_MODELS = ["rnafm", "nt", "hyenadna"]
AVAILABLE_MODELS = MULTIMOL_MODELS + LEGACY_MODELS + ["evo", "dnabert2", "caduceus"]

A100_MODELS = {"evo", "rinalmo"}
# Three models have no CPU path. The Evo adapter refuses to put 7B parameters on
# a CPU, mamba-ssm has no CPU kernel, and DNABERT-2's Triton attention asserts
# `q.is_cuda`. They smoke on a GPU; the rest smoke on a CPU for a tenth of the
# cost.
GPU_SMOKE_MODELS = {"evo", "dnabert2", "caduceus"}


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


def _allow_torch_load():
    """Let transformers 4.49 load a checkpoint under torch 2.1.2.

    Copied from `modal_phase6_flashattn_models.py`, and reached only on the Evo
    image: transformers refuses `torch.load` below torch 2.6 over
    CVE-2025-32434, and Evo's torch is held at 2.1.2 by flash-attn. The other
    three images are above the threshold and the attribute is absent, so this is
    a no-op there. Evo is loaded from the same pinned revision this project has
    loaded it from since July.
    """
    import transformers.modeling_utils as modeling_utils
    import transformers.utils.import_utils as import_utils

    def permitted():
        return None

    for module in (import_utils, modeling_utils):
        if hasattr(module, "check_torch_load_is_safe"):
            module.check_torch_load_is_safe = permitted


def _library_versions():
    """Every package that could move a number, as resolved in this container."""
    from importlib.metadata import PackageNotFoundError, version

    names = ["torch", "triton", "numpy", "scipy", "scikit-learn", "transformers",
             "multimolecule", "huggingface-hub", "flash-attn", "mamba-ssm"]
    out = {}
    for name in names:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = None
    return out


def _load_on_device(model_name):
    """Adapter, loaded and moved, with the device it landed on."""
    import torch

    from phase6_compensatory_mutation import load_adapter

    _allow_torch_load()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    adapter = load_adapter(model_name)
    adapter.load()
    if getattr(adapter, "model", None) is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    return adapter, device


def _run_model(model_name, commit, phase6_only):
    import os
    os.chdir("/root/project")

    from datetime import datetime, timezone

    def now():
        return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    from family_checkpoint import FamilyCheckpoint
    from phase6_compensatory_mutation import (
        ADAPTER_OFFSETS,
        NON_CHARACTER_TOKENIZERS,
        load_rfam_families,
        run_phase6,
    )
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    families = load_rfam_families()
    withdrawn = sorted(
        json.loads(p.read_text())["name"]
        for p in Path("data/rfam_families").glob("*.json")
        if "excluded" in json.loads(p.read_text())
    )
    adapter, device = _load_on_device(model_name)
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
    print(f"[{now()}] {model_name} loaded on {device}")

    out_dir = Path("/results") / model_name
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "stamp.json").write_text(json.dumps(stamp, indent=2) + "\n")
    vol.commit()

    def checkpoint(stage):
        return FamilyCheckpoint(out_dir / f"_shards_{stage}.json", stamp=stamp,
                                on_write=vol.commit)

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
    """Import, load the panel, load the model, score six families.

    Ten GPU containers that all die on the same missing package cost far more
    than one container that finds it. This runs the same import path, the same
    adapter load and the same device placement as the real thing.

    Families are tried shortest-first until six score, rather than the six
    shortest being taken outright. A 6-mer tokenizer resolves a 30-nucleotide
    family into too few tokens to perturb, so NT scored none of the six shortest
    and HyenaDNA one -- a property of the selection, which would have read here
    as a broken container.
    """
    import os
    os.chdir("/root/project")

    from phase6_compensatory_mutation import load_rfam_families, run_phase6

    families = sorted(load_rfam_families(), key=lambda f: len(f["sequence"]))
    adapter, device = _load_on_device(model_name)
    ok = []
    for family in families:
        scored = run_phase6(adapter, [family], device=device, compute_null=False)
        if not scored["per_rna"][family["name"]].get("skipped"):
            ok.append(family["name"])
        if len(ok) == 6:
            break
    print(f"{model_name}: {len(families)} families, loaded on {device}, scored {ok}")
    print(f"  libraries: {_library_versions()}")
    if not ok:
        raise RuntimeError(f"{model_name} scored no family in the whole panel")
    return model_name


@app.function(image=multimol_image, timeout=3600)
def smoke_multimol(model_name: str):
    return _smoke(model_name)


@app.function(image=legacy_image, timeout=3600)
def smoke_legacy(model_name: str):
    return _smoke(model_name)


@app.function(image=evo_image, gpu="A10G", timeout=3600)
def smoke_evo(model_name: str = "evo"):
    """On a GPU, because the Evo adapter refuses to place 7B parameters on a CPU."""
    return _smoke(model_name)


@app.function(image=dnabert2_image, gpu="A10G", timeout=3600)
def smoke_dnabert2(model_name: str = "dnabert2"):
    """On a GPU, because the bundled Triton attention asserts `q.is_cuda`."""
    return _smoke(model_name)


@app.function(image=caduceus_image, gpu="A10G", timeout=3600)
def smoke_caduceus(model_name: str = "caduceus"):
    """On a GPU, because mamba-ssm has no CPU kernel to fall back to."""
    return _smoke(model_name)


@app.function(image=multimol_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_multimol_a10g(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


@app.function(image=multimol_image, gpu="A100", timeout=86400, volumes={"/results": vol})
def run_multimol_a100(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


@app.function(image=legacy_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_legacy(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


@app.function(image=evo_image, gpu="A100", timeout=86400, volumes={"/results": vol})
def run_evo(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


@app.function(image=dnabert2_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_dnabert2(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


@app.function(image=caduceus_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_caduceus(model_name: str, commit: str, phase6_only: bool = False):
    return _run_model(model_name, commit, phase6_only)


def _route(model_name):
    """The (function, label) each model runs under."""
    if model_name == "evo":
        return run_evo, "evo/A100"
    if model_name == "dnabert2":
        return run_dnabert2, "dnabert2/A10G"
    if model_name == "caduceus":
        return run_caduceus, "caduceus/A10G"
    if model_name in LEGACY_MODELS:
        return run_legacy, "legacy/A10G"
    if model_name in A100_MODELS:
        return run_multimol_a100, "multimol/A100"
    return run_multimol_a10g, "multimol/A10G"


def _smoke_route(model_name):
    """The smoke function for a model, on its own image."""
    if model_name == "evo":
        return smoke_evo
    if model_name == "dnabert2":
        return smoke_dnabert2
    if model_name == "caduceus":
        return smoke_caduceus
    if model_name in LEGACY_MODELS:
        return smoke_legacy
    return smoke_multimol


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
        # Spawned rather than called, so one broken adapter reports itself
        # alongside the nine that work instead of hiding them behind its own
        # traceback, and so the four images build concurrently.
        handles = {m: _smoke_route(m).spawn(m) for m in requested}
        outcomes = {}
        for model_name, handle in handles.items():
            try:
                outcomes[model_name] = handle.get()
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
