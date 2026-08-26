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
        .env({"PYTHONPATH": "/root/project:/root/project/scripts",
              # cuBLAS picks a workspace per call above CUDA 10.2 and the choice
              # is not reproducible, so `use_deterministic_algorithms(True)`
              # refuses without this. Unset, it is the reason one model's no-op
              # floor exceeded its own signal (D18); it must be set before torch
              # initializes CUDA, which is why it is an image variable rather
              # than something the run sets.
              "CUBLAS_WORKSPACE_CONFIG": ":4096:8"})
        .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
        .add_local_dir("scripts", "/root/project/scripts")
        .add_local_dir("data/rfam_families", "/root/project/data/rfam_families")
        # 2.9 MB of Stockholm alignments, for the multi-sequence replication.
        .add_local_dir("data/rfam_seeds", "/root/project/data/rfam_seeds")
        # The 3--5 seed replicates per family that the deposited within-family
        # variance was computed on, already extracted and folded.
        .add_local_dir("data/multi_sequence", "/root/project/data/multi_sequence")
        # Five folded CAG-repeat fragments, in the same record format as the
        # Rfam families, so the case study runs through the panel's own stages.
        .add_local_dir("data/htt_fragments", "/root/project/data/htt_fragments")
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
TRAINED_MODELS = MULTIMOL_MODELS + LEGACY_MODELS + ["evo", "dnabert2", "caduceus"]

# Randomized-weight controls, by the procedure in
# `modal_phase6_untrained_all.py`: the trained adapter is loaded and every
# parameter re-initialized, so the architecture and the tokenizer are held and
# only the learned weights are destroyed.
#
# Seven, not ten, and the seven are chosen by what a registered hypothesis or a
# reported row needs rather than by symmetry. H10 (RNA-FM sign test), H11 (NT
# attention, trained versus untrained), H14 and H15 (RiNALMo ratio and sign
# test), H17 (ERNIE-RNA sign test), H19 (SpliceBERT ratio) and H20 (DNABERT-2
# attention contrast) are all defined as a trained-versus-untrained difference,
# so without the control the hypothesis has no result at all. UTR-LM carries no
# such hypothesis but does carry an untrained row in the attention table.
# HyenaDNA, Caduceus and Evo appear in no trained-versus-untrained comparison
# the manuscript reports, and Evo's control would cost an A100 to produce a
# number nothing reads.
UNTRAINED_MODELS = ["ernierna_untrained", "rnafm_untrained", "nt_untrained",
                    "rinalmo_untrained", "splicebert_untrained",
                    "dnabert2_untrained", "utrlm_untrained"]
RANDOM_INIT_SEED = 42

# Verbatim from `modal_transversion_control.py`, which produced the deposited
# transversion numbers. A Watson-Crick swap at a stem position leaves a valid
# pair (A-U becomes U-A), so the control substitutes a purine for a pyrimidine
# and back, which cannot pair with the original partner. `phases_1_to_5` reads
# COMPLEMENT as a module global, so replacing the dict is the whole change.
TRANSVERSION_COMPLEMENT = {"A": "C", "U": "C", "C": "A", "G": "U", "T": "C"}

AVAILABLE_MODELS = TRAINED_MODELS + UNTRAINED_MODELS

A100_MODELS = {"evo", "rinalmo"}
# Three models have no CPU path. The Evo adapter refuses to put 7B parameters on
# a CPU, mamba-ssm has no CPU kernel, and DNABERT-2's Triton attention asserts
# `q.is_cuda`. They smoke on a GPU; the rest smoke on a CPU for a tenth of the
# cost.
# Read by `_smoke_route`. Evo, DNABERT-2 and Caduceus are here because they have
# no CPU path at all. RNA-FM is here for a different reason: its adapter attaches
# a layer norm to the model by hand, and a module that does not travel with
# `.to(device)` fails only where there is more than one device. A CPU smoke
# cannot see that, because everything agrees on one device.
GPU_SMOKE_MODELS = {"evo", "dnabert2", "caduceus", "rnafm"}


def base_model(model_name):
    """The adapter a run loads, with any control suffix removed."""
    return model_name[:-len("_untrained")] if model_name.endswith("_untrained") else model_name


def _ablate_pairwise_bias(model):
    """Remove ERNIE-RNA's hardcoded pairing prior, and make the removal stick.

    `pairwise_bias_map` is a non-persistent buffer holding a Watson-Crick and
    wobble table -- A-U 2.0, C-G 3.0, G-U 0.8 -- so `_randomize` never reaches it
    and an untrained control keeps it (D20). Zeroing it and the projection above
    it is the ablation registered in
    `preregistration/PREREGISTRATION_ERNIERNA_BIAS_RUNG3.md`.

    `_inited` is the part that matters. `ErnieRnaEmbeddings` rebuilds the buffer
    on the first forward pass while that flag is false, as a workaround for
    transformers v5 leaving non-persistent buffers on the meta device. An
    ablation that does not set it is silently undone by the first sequence, which
    is what `modal_ernierna_ablation.py` did in July (D19).
    """
    import torch

    target = None
    for module in model.modules():
        if hasattr(module, "pairwise_bias_map"):
            target = module
            break
    if target is None:
        raise RuntimeError("no module holds pairwise_bias_map; nothing to ablate")

    with torch.no_grad():
        target.pairwise_bias_map.zero_()
        zeroed = int(target.pairwise_bias_map.numel())
        for name, param in model.named_parameters():
            if "pairwise_bias_proj" in name:
                param.zero_()
                zeroed += param.numel()
    if not hasattr(target, "_inited"):
        raise RuntimeError(
            "the module holding pairwise_bias_map has no _inited flag, so the "
            "rebuild guard this ablation depends on is not the one documented")
    target._inited = True
    print(f"  ablated the pairwise bias: {zeroed:,} values zeroed, "
          f"_inited set so the first forward cannot rebuild it")
    return zeroed


def _reseed_buffers(adapter, seed, device="cpu"):
    """Redraw every non-persistent buffer, so a control keeps no fixed table.

    `_randomize` iterates `named_parameters()`, which is what a random-init
    control is conventionally built from, and a buffer is not a parameter. That
    left ERNIE-RNA's Watson-Crick table intact (D20). Two other controls sit
    above their own chance rates with no buffer identified -- RiNALMo by 0.095
    and RNA-FM by 0.084 -- and RiNALMo holds 34 non-persistent buffers, mostly
    rotary inverse frequencies.

    Replacing every such buffer with noise of the same shape tests whether those
    buffers carry the excess. A control that still exceeds its chance rate after
    this holds nothing outside its parameters that the probe reads.
    """
    import torch

    # A non-persistent buffer materializes from the meta device and every model
    # here repairs its own on the first forward pass (D19). Reseeding before that
    # pass is silently undone, which is what voided the July ablation and what an
    # earlier version of this function did: all three controls came back with
    # numbers identical to their unreseeded runs, including the one whose buffer
    # is known to carry the effect. One forward pass first, so the values being
    # replaced are the ones the model actually computes with.
    model = adapter.model
    with torch.no_grad():
        adapter.get_all_layer_embeddings(
            adapter.tokenize("ACGUACGUACGUACGUACGU").to(device))
    for module in model.modules():
        if hasattr(module, "_inited"):
            module._inited = True

    generator = torch.Generator(device="cpu").manual_seed(seed)
    persistent = set(model.state_dict())
    touched = 0
    with torch.no_grad():
        for name, tensor in model.named_buffers():
            if name in persistent or tensor is None or tensor.numel() == 0:
                continue
            if not tensor.is_floating_point():
                continue
            noise = torch.randn(tensor.shape, generator=generator,
                                dtype=tensor.dtype).to(tensor.device)
            tensor.copy_(noise * float(tensor.std()) if tensor.numel() > 1
                         else noise)
            touched += 1
    # Whether the redraw survives is the whole question. ERNIE-RNA rebuilds its
    # buffer behind an `_inited` flag, which the loop above sets; RiNALMo's 33
    # rotary tables are repaired by some other mechanism and no flag here
    # reaches it. A reseed that a later forward pass overwrites produces numbers
    # identical to the unreseeded run, which is indistinguishable from the
    # buffers not mattering -- the false negative this exists to avoid.
    written = {name: tensor.detach().clone()
               for name, tensor in model.named_buffers()
               if name not in persistent and tensor is not None
               and tensor.numel() and tensor.is_floating_point()}
    with torch.no_grad():
        adapter.get_all_layer_embeddings(
            adapter.tokenize("ACGUACGUACGUACGUACGUACGUACGU").to(device))
    after = dict(model.named_buffers())
    reverted = sorted(name for name, value in written.items()
                      if not torch.equal(after[name], value))
    if reverted:
        raise RuntimeError(
            f"{len(reverted)} of {len(written)} reseeded buffer(s) were restored "
            f"by a later forward pass, so the intervention did not take: "
            f"{reverted[:4]}. Reseeding cannot be verified for this model.")
    print(f"  reseeded {touched} non-persistent float buffer(s) at seed {seed}, "
          f"{len(written)} verified to survive a further forward pass")
    return touched


def result_dir(model_name, transversion, ablate_bias=False, synthetic=False,
               multi_seq=False, htt=False, reseed_buffers=False):
    """Where a run writes, on the volume.

    The transversion control gets its own directory. `_run_model` empties a
    directory whose stored stamp disagrees with the one it is about to write,
    and the control's stamp disagrees on the commit -- it is launched from a
    later one than the Watson-Crick runs, which no stack can re-run cheaply --
    so writing the two into one directory would have the control delete the run
    it exists to be compared against.
    """
    suffix = "_transversion" if transversion else ""
    if synthetic:
        suffix += "_synthetic"
    if multi_seq:
        suffix += "_multiseq"
    if htt:
        suffix += "_htt"
    if reseed_buffers:
        suffix += "_reseeded"
    # The ablation is a different model, not a different metric, so it cannot
    # share a directory with the intact run it is compared against.
    if ablate_bias:
        suffix += "_noattnbias"
    return f"{model_name}{suffix}"


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


# Verbatim from `modal_phase6_synthetic_covariation.py`, which produced the
# deposited synthetic numbers. Structure is held and nucleotides are reassigned:
# a uniform Watson-Crick pair at every paired position, a uniform nucleotide at
# every unpaired one. Evolutionary covariation is destroyed and the geometry is
# not, so a model scoring on covariation loses its score and one scoring on
# pairing geometry keeps it.
WC_PAIRS = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")]
SYNTHETIC_NUCLEOTIDES = ["A", "U", "G", "C"]
# Five per family, as in the run that produced the deposited synthetic numbers.
SYNTHETIC_PER_FAMILY = 5


def synthetic_families(families, n_synthetic, rng):
    """`n_synthetic` sequences per family, each keeping that family's structure."""
    out = []
    for family in families:
        structure = family.get("dot_bracket", "")
        if not structure:
            continue
        for index in range(n_synthetic):
            sequence = ["N"] * len(structure)
            stack = []
            for i, char in enumerate(structure):
                if char == "(":
                    stack.append(i)
                elif char == ")":
                    j = stack.pop()
                    left, right = WC_PAIRS[rng.integers(len(WC_PAIRS))]
                    sequence[j] = left
                    sequence[i] = right
                else:
                    sequence[i] = SYNTHETIC_NUCLEOTIDES[
                        rng.integers(len(SYNTHETIC_NUCLEOTIDES))]
            out.append({
                "name": f"{family['name']}_syn{index}",
                "sequence": "".join(sequence),
                "dot_bracket": structure,
                "original_family": family["name"],
                "synthetic_index": index,
            })
    return out


def _pin_numerics():
    """Turn off reduced-precision matmul and record what the flags actually were.

    `torch.backends.cuda.matmul.allow_tf32` has defaulted to False since torch
    1.12, but `torch.backends.cudnn.allow_tf32` has always defaulted to True, so
    a repository that sets neither is not thereby in full float32. TF32 carries a
    10-bit mantissa, which would put the noise floor near 1e-3 rather than near
    1e-6 and would make most of the Rung 3 panel uninterpretable rather than
    merely imprecise.

    The values are read before they are set, so the returned dict says what the
    previous runs were computed under and not merely what this one asks for.
    """
    import torch

    before = {
        "cuda_matmul_allow_tf32": bool(torch.backends.cuda.matmul.allow_tf32),
        "cudnn_allow_tf32": bool(torch.backends.cudnn.allow_tf32),
        "cudnn_benchmark": bool(torch.backends.cudnn.benchmark),
    }
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    # Kernel selection by timing makes the same input take different code paths
    # on different runs, which is a candidate for the one model whose results
    # move across otherwise identical runs.
    torch.backends.cudnn.benchmark = False
    # `warn_only` was the reason D18's cause stayed a guess: a kernel with no
    # deterministic implementation warned and carried on, so nothing named it.
    # Forcing it either raises on that kernel, which identifies it, or succeeds,
    # which rules determinism out as the explanation. The failure is caught and
    # recorded rather than killing the run, because the point is the diagnosis.
    determinism = "forced"
    try:
        torch.use_deterministic_algorithms(True)
    except Exception as exc:  # noqa: BLE001 -- recorded, not swallowed
        determinism = f"refused: {type(exc).__name__}: {exc}"
        torch.use_deterministic_algorithms(True, warn_only=True)
        print(f"  determinism {determinism}")
    print(f"  numerics: defaults were {before}, now pinned to full float32")
    return {"defaults_observed": before,
            "determinism": determinism,
            "pinned": "tf32 off, cudnn.benchmark off, deterministic algorithms"}


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


def _randomize(model, seed):
    """Destroy the learned weights and keep the architecture.

    Verbatim from `modal_phase6_untrained_all.py`, with the torch generator
    seeded as well: that script seeded numpy, and every initializer it calls
    draws from torch's generator, so its controls were not reproducible. Nothing
    carries over from those runs, so seeding here costs no comparability.
    """
    import torch

    torch.manual_seed(seed)
    before = sum(p.numel() for p in model.parameters())
    with torch.no_grad():
        for _, param in model.named_parameters():
            if param.dim() >= 2:
                torch.nn.init.xavier_normal_(param)
            else:
                torch.nn.init.normal_(param, std=0.02)
    after = sum(p.numel() for p in model.parameters())
    assert before == after, "parameter count changed during randomization"
    return after


def _load_on_device(model_name, ablate_bias=False, reseed_buffers=False):
    """Adapter, loaded and moved, with the device it landed on.

    A `_untrained` suffix loads the trained adapter and then re-initializes it,
    which is what makes the control a control: same tokenizer, same shapes, same
    forward, no learned weights.
    """
    import torch

    from phase6_compensatory_mutation import load_adapter

    _allow_torch_load()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    adapter = load_adapter(base_model(model_name))
    adapter.load()
    if model_name.endswith("_untrained"):
        n = _randomize(adapter.model, RANDOM_INIT_SEED)
        print(f"  randomized {n:,} parameters at seed {RANDOM_INIT_SEED}")
    if ablate_bias:
        _ablate_pairwise_bias(adapter.model)
    if getattr(adapter, "model", None) is not None and device == "cuda":
        adapter.model = adapter.model.to(device)
    # After the move, because reseeding runs a forward pass and a probe built on
    # the CPU against a model on the GPU is the same mismatch that has already
    # cost three runs in this file.
    if reseed_buffers and model_name.endswith("_untrained"):
        _reseed_buffers(adapter, RANDOM_INIT_SEED, device)
    return adapter, device


def _run_model(model_name, commit, phase6_only, transversion=False,
               ablate_bias=False, synthetic=False, multi_seq=False, htt=False,
               reseed_buffers=False):
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
    import phases_1_to_5
    from phases_1_to_5 import (
        run_attention_contact,
        run_mutation_sensitivity,
        run_structure_probing,
    )

    if transversion:
        phases_1_to_5.COMPLEMENT = TRANSVERSION_COMPLEMENT
        print(f"  COMPLEMENT -> {TRANSVERSION_COMPLEMENT}")

    families = load_rfam_families()
    if htt:
        # The case study asks whether probing accuracy tracks structure or the
        # composition a CAG repeat forces, so it runs the panel's own stages on
        # fragments whose structures are folded and versioned rather than
        # predicted inside the container.
        htt_dir = Path("data/htt_fragments")
        families = [json.loads(p.read_text())
                    for p in sorted(htt_dir.glob("*.json"))]
        print(f"  HTT: {len(families)} CAG-repeat fragments, "
              f"{[f['n_cag'] for f in families]} repeats")
    if synthetic:
        import numpy as _np
        n_before = len(families)
        families = synthetic_families(families, SYNTHETIC_PER_FAMILY,
                                      _np.random.default_rng(RANDOM_INIT_SEED))
        print(f"  synthetic: {n_before} families -> {len(families)} sequences, "
              f"{SYNTHETIC_PER_FAMILY} per family at seed {RANDOM_INIT_SEED}")
    multiseq_members: dict[str, list[str]] = {}
    if multi_seq:
        # Replication across the seed replicates of each family, which is what
        # the deposited within-family variance was measured on: the same
        # `run_mutation_sensitivity` over `<family>_seq<i>` pseudo-families,
        # regrouped afterwards. Running it this way rather than through its own
        # loop is what gives it the checkpoint -- a reclaimed container took the
        # whole of the first attempt, which wrote only at the end.
        multi_dir = Path("data/multi_sequence")
        panel = {f["name"] for f in families}
        replicates = []
        for path in sorted(multi_dir.glob("*_multi.json")):
            record = json.loads(path.read_text())
            # A family withdrawn from the panel is withdrawn here too, or the
            # replication covers a set the rest of the manuscript does not.
            if record["name"] not in panel:
                continue
            members = []
            for index, entry in enumerate(record["sequences"]):
                name = f"{record['name']}_seq{index}"
                replicates.append({"name": name, "sequence": entry["sequence"],
                                   "dot_bracket": entry["dot_bracket"]})
                members.append(name)
            multiseq_members[record["name"]] = members
        if not replicates:
            raise RuntimeError(
                f"no seed replicates matched the panel in {multi_dir}")
        print(f"  multi-sequence: {len(replicates)} sequences across "
              f"{len(multiseq_members)} of {len(panel)} panel families")
        families = replicates
    withdrawn = sorted(
        json.loads(p.read_text())["name"]
        for p in Path("data/rfam_families").glob("*.json")
        if "excluded" in json.loads(p.read_text())
    )
    numerics = _pin_numerics()
    adapter, device = _load_on_device(model_name, ablate_bias, reseed_buffers)
    stamp = {
        "commit": commit,
        "panel_sha256": panel_stamp(families),
        "n_families": len(families),
        "withdrawn": withdrawn,
        "seeding": "family_seed.family_rng, derived from the family name",
        "deviation": "DEVIATIONS.md, 2026-08-24",
        "libraries": _library_versions(),
        "numerics": numerics,
        "device": device,
        "ablation": ("pairwise_bias_map and pairwise_bias_proj zeroed, _inited set"
                     if ablate_bias else None),
        "sequences": ("synthetic, structure preserved and nucleotides reassigned"
                      if synthetic else
                      "Rfam seed replicates, 3--5 per family" if multi_seq else
                      "natural"),
        "weights": ("randomized, xavier_normal_ on matrices and normal_(0, 0.02) "
                    f"on vectors, seed {RANDOM_INIT_SEED}"
                    if model_name.endswith("_untrained") else "pretrained"),
    }
    print(f"[{now()}] {model_name}: {len(families)} families, panel "
          f"{stamp['panel_sha256'][:12]}, withdrawn {withdrawn}")
    print(f"[{now()}] {model_name} loaded on {device}")

    out_dir = Path("/results") / result_dir(model_name, transversion, ablate_bias,
                                            synthetic, multi_seq, htt,
                                            reseed_buffers)
    out_dir.mkdir(parents=True, exist_ok=True)

    # The directory is emptied when its stamp changes, so it never holds two
    # runs at once. FamilyCheckpoint already refuses a stale shard, which is
    # enough for a run that finishes; a run that dies partway leaves its result
    # files behind, and the next run would write a new stamp.json beside them.
    # A reader who trusts the stamp then reads last month's numbers.
    previous = out_dir / "stamp.json"
    if previous.exists() and json.loads(previous.read_text()) != stamp:
        stale = sorted(p for p in out_dir.iterdir() if p.is_file())
        print(f"[{now()}] {model_name}: clearing {len(stale)} file(s) from a "
              f"different run: {[p.name for p in stale]}")
        for path in stale:
            path.unlink()
    previous.write_text(json.dumps(stamp, indent=2) + "\n")
    vol.commit()

    def checkpoint(stage):
        return FamilyCheckpoint(out_dir / f"_shards_{stage}.json", stamp=stamp,
                                on_write=vol.commit)

    def positions(stage):
        """Per-position distances, written as a shard file that is kept.

        Same machinery as a checkpoint -- atomic write, volume commit, stamp
        guard -- because the requirement is the same: a container reclaimed
        after thirty families must not take the thirty with it. The file is a
        deliverable rather than scratch, so it is named for the run.
        """
        return FamilyCheckpoint(out_dir / f"{model_name}_{stage}_positions.json",
                                stamp=stamp, on_write=vol.commit)

    def save(name, payload):
        payload = {"model": model_name, "stamp": stamp,
                   "completed": now(), **payload}
        (out_dir / name).write_text(json.dumps(payload, indent=2, default=str))
        vol.commit()
        print(f"[{now()}] wrote {name}")

    if multi_seq:
        # Within-family variance of the stem-loop sensitivity ratio, over the
        # seed replicates built above. Each replicate is scored as its own
        # family, so the checkpoint is per sequence and a reclaimed container
        # resumes at the sequence it died on.
        import numpy as _np

        print(f"[{now()}] {model_name}: multi-sequence replication over "
              f"{len(families)} sequences")
        mutation = run_mutation_sensitivity(
            adapter, base_model(model_name), families, device=device,
            checkpoint=checkpoint("multiseq"))

        per_family = {}
        for family, members in multiseq_members.items():
            scored = [mutation["per_rna"][name] for name in members
                      if "best_ratio" in mutation["per_rna"].get(name, {})]
            entry = {"n_sequences": len(members), "n_scored": len(scored)}
            if scored:
                ratios = [r["best_ratio"] for r in scored]
                exceeds = [bool(r.get("exceeds_nuc_null", False)) for r in scored]
                mean = float(_np.mean(ratios))
                sd = float(_np.std(ratios, ddof=1)) if len(ratios) > 1 else 0.0
                entry.update({
                    "per_sequence_ratios": ratios,
                    "per_sequence_exceeds_nuc": exceeds,
                    "mean_ratio": mean,
                    "sd_ratio": sd,
                    # A ratio near zero would make this meaningless, and no
                    # model in the panel produces one; guarding it silently
                    # would hide the case where one did.
                    "cv_ratio": sd / mean if mean > 1e-10 else None,
                    "n_exceeds_nuc": sum(exceeds),
                })
            per_family[family] = entry

        scored_families = [v for v in per_family.values()
                           if v.get("cv_ratio") is not None]
        summary = {
            "n_families": len(per_family),
            "n_families_scored": len(scored_families),
            "n_total_sequences": sum(v["n_scored"] for v in per_family.values()),
            "mean_cv": float(_np.mean([v["cv_ratio"] for v in scored_families])),
            "median_cv": float(_np.median([v["cv_ratio"] for v in scored_families])),
            "mean_ratio": float(_np.mean([v["mean_ratio"] for v in scored_families])),
        }
        save(f"{model_name}_multiseq.json",
             {"metric": "multi_sequence_mutation_sensitivity",
              "summary": summary, "per_family": per_family})
        print(f"[{now()}] {model_name} MULTI-SEQ COMPLETE: "
              f"mean CV {summary['mean_cv']:.3f} over "
              f"{summary['n_families_scored']} families")
        return model_name

    # Every stage in phases_1_to_5 dispatches on the key it is given -- which
    # tokenizer the model has, whether it has attention at all, and NT's and
    # DNABERT-2's token-to-nucleotide mapping. A `_untrained` suffix matches
    # none of those sets, so a control passed its own name silently takes the
    # character-tokenizer path, reports "no attention (SSM architecture)" for a
    # transformer, and skips the 6-mer offset. `modal_phase6_untrained_all.py`
    # avoids this by keying on the base name and suffixing only the output; the
    # same split is kept here, with `model_name` for paths and the stamp.
    key = base_model(model_name)
    if transversion:
        # The alphabet is recorded in the payload as well as in the directory
        # name, so a file read on its own still says which substitution produced
        # it.
        mutation = run_mutation_sensitivity(adapter, key, families, device=device,
                                            checkpoint=checkpoint("transversion"),
                                            positions=positions("transversion"))
        save(f"{model_name}_transversion.json",
             {"experiment": "transversion_control",
              "complement": TRANSVERSION_COMPLEMENT,
              "mutation_trained": mutation})
        print(f"[{now()}] {model_name} TRANSVERSION COMPLETE. mean ratio "
              f"{mutation.get('mean_best_ratio')}")
        return {"model": model_name, "mean_ratio": mutation.get("mean_best_ratio"),
                "panel_sha256": stamp["panel_sha256"]}

    if not phase6_only:
        print(f"[{now()}] {model_name}: mutation sensitivity")
        mutation = run_mutation_sensitivity(adapter, key, families,
                                            device=device,
                                            checkpoint=checkpoint("mutation"),
                                            positions=positions("rung1"))
        print(f"[{now()}] {model_name}: attention-contact")
        attention = run_attention_contact(adapter, key, families, device=device,
                                          checkpoint=checkpoint("attention"))
        # Probing pools every family's embeddings before it fits, so there is no
        # per-family unit to checkpoint; it is saved the moment it returns.
        print(f"[{now()}] {model_name}: structure probing")
        probing = run_structure_probing(adapter, key, families, device=device)
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
                        offset=ADAPTER_OFFSETS.get(key, 0),
                        checkpoint=checkpoint("phase6"))
    save(f"{model_name}_phase6_ps.json",
         {"phase": 6, "metric": "perturbation_specificity",
          "preregistration": "PREREGISTRATION_PHASE6_V2.md",
          "offset": ADAPTER_OFFSETS.get(key, 0),
          "tokenizer_caveated": key in NON_CHARACTER_TOKENIZERS,
          "results": phase6})
    print(f"[{now()}] {model_name} COMPLETE. mean PS {phase6.get('mean_best_ps')}, "
          f"confirmatory {phase6.get('families_total')}, "
          f"exceeding null {phase6.get('families_exceeding_null_conservative')} "
          f"(primary null: {phase6.get('families_exceeding_null_primary')})")
    return {"model": model_name, "mean_ps": phase6.get("mean_best_ps"),
            "panel_sha256": stamp["panel_sha256"]}


def _inspect_buffers(adapter, model_name, device="cpu"):
    """What each model holds outside `named_parameters()`, and whether it survives.

    `_randomize` builds the random-init control by iterating `named_parameters()`.
    A buffer is not a parameter, so anything held in one survives untouched, and
    the control keeps whatever that buffer encodes. For ERNIE-RNA that is a
    Watson-Crick table and the control climbs all three rungs of the ladder (D20).

    This asks the same question of every model in the panel: what is outside the
    parameters, does it hold real values once a forward pass has run, and is any
    of it plausibly about the task rather than about position or dtype.
    """
    import torch

    model = getattr(adapter, "model", None)
    if model is None:
        print(f"{model_name}: adapter exposes no model")
        return
    persistent = set(model.state_dict())
    rows = [(name, tensor) for name, tensor in model.named_buffers()
            if name not in persistent]

    # Buffers materialize from the meta device after transformers v5 loads a
    # checkpoint, so they are only meaningful once a forward pass has run.
    with torch.no_grad():
        adapter.get_all_layer_embeddings(
            adapter.tokenize("ACGUACGUACGUACGUACGU").to(device))

    after = dict(model.named_buffers())
    total = sum(int(t.numel()) for _, t in rows)
    print(f"{model_name}: {len(rows)} non-persistent buffer(s), "
          f"{total:,} values, {len(persistent)} persistent tensors")
    shapes = {}
    for name, _ in rows:
        tensor = after.get(name)
        if tensor is None:
            continue
        stem = name.split(".")[-1]
        entry = shapes.setdefault(stem, {"count": 0, "numel": 0, "nonzero": 0,
                                         "shape": list(tensor.shape)})
        entry["count"] += 1
        entry["numel"] += int(tensor.numel())
        entry["nonzero"] += int((tensor != 0).sum())
    for stem, entry in sorted(shapes.items(), key=lambda kv: -kv[1]["numel"]):
        print(f"    {stem:34s} x{entry['count']:<3d} shape={entry['shape']} "
              f"nonzero={entry['nonzero']}/{entry['numel']}")


def _smoke(model_name, ablate_bias=False, buffers_only=False,
           reseed_buffers=False):
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

    # The real run pins these before loading anything, and the calls differ by
    # torch version across the four images. A smoke that skips them leaves the
    # one piece of startup code that could kill all 27 containers untested.
    numerics = _pin_numerics()
    families = sorted(load_rfam_families(), key=lambda f: len(f["sequence"]))
    adapter, device = _load_on_device(model_name, ablate_bias, reseed_buffers)
    if buffers_only:
        _inspect_buffers(adapter, model_name, device)
        return model_name
    ok = []
    for family in families:
        scored = run_phase6(adapter, [family], device=device, compute_null=False)
        if not scored["per_rna"][family["name"]].get("skipped"):
            ok.append(family["name"])
        if len(ok) == 6:
            break
    print(f"{model_name}: {len(families)} families, loaded on {device}, scored {ok}")
    print(f"  libraries: {_library_versions()}")
    print(f"  numerics: {numerics['defaults_observed']}")
    if not ok:
        raise RuntimeError(f"{model_name} scored no family in the whole panel")
    return model_name


@app.function(image=multimol_image, timeout=3600)
def smoke_multimol(model_name: str, ablate_bias: bool = False, buffers_only: bool = False,
                   reseed_buffers: bool = False):
    return _smoke(model_name, ablate_bias, buffers_only, reseed_buffers)


@app.function(image=legacy_image, timeout=3600)
def smoke_legacy(model_name: str, ablate_bias: bool = False, buffers_only: bool = False,
                   reseed_buffers: bool = False):
    return _smoke(model_name, ablate_bias, buffers_only, reseed_buffers)


@app.function(image=legacy_image, gpu="A10G", timeout=3600)
def smoke_legacy_gpu(model_name: str, ablate_bias: bool = False, buffers_only: bool = False,
                   reseed_buffers: bool = False):
    """The legacy image with a GPU, for a model whose failure needs two devices."""
    return _smoke(model_name, ablate_bias, buffers_only, reseed_buffers)


@app.function(image=evo_image, gpu="A10G", timeout=3600)
def smoke_evo(model_name: str = "evo", ablate_bias: bool = False, buffers_only: bool = False,
                   reseed_buffers: bool = False):
    """On a GPU, because the Evo adapter refuses to place 7B parameters on a CPU."""
    return _smoke(model_name, ablate_bias, buffers_only, reseed_buffers)


@app.function(image=dnabert2_image, gpu="A10G", timeout=3600)
def smoke_dnabert2(model_name: str = "dnabert2", ablate_bias: bool = False, buffers_only: bool = False,
                   reseed_buffers: bool = False):
    """On a GPU, because the bundled Triton attention asserts `q.is_cuda`."""
    return _smoke(model_name, ablate_bias, buffers_only, reseed_buffers)


@app.function(image=caduceus_image, gpu="A10G", timeout=3600)
def smoke_caduceus(model_name: str = "caduceus", ablate_bias: bool = False, buffers_only: bool = False,
                   reseed_buffers: bool = False):
    """On a GPU, because mamba-ssm has no CPU kernel to fall back to."""
    return _smoke(model_name, ablate_bias, buffers_only, reseed_buffers)


@app.function(image=multimol_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_multimol_a10g(model_name: str, commit: str, phase6_only: bool = False,
                     transversion: bool = False, ablate_bias: bool = False,
                     synthetic: bool = False, multi_seq: bool = False,
                     htt: bool = False,
                     reseed_buffers: bool = False):
    return _run_model(model_name, commit, phase6_only, transversion,
                      ablate_bias, synthetic, multi_seq, htt, reseed_buffers)


@app.function(image=multimol_image, gpu="A100", timeout=86400, volumes={"/results": vol})
def run_multimol_a100(model_name: str, commit: str, phase6_only: bool = False,
                     transversion: bool = False, ablate_bias: bool = False,
                     synthetic: bool = False, multi_seq: bool = False,
                     htt: bool = False,
                     reseed_buffers: bool = False):
    return _run_model(model_name, commit, phase6_only, transversion,
                      ablate_bias, synthetic, multi_seq, htt, reseed_buffers)


@app.function(image=legacy_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_legacy(model_name: str, commit: str, phase6_only: bool = False,
              transversion: bool = False, ablate_bias: bool = False,
                     synthetic: bool = False, multi_seq: bool = False,
                     htt: bool = False,
                     reseed_buffers: bool = False):
    return _run_model(model_name, commit, phase6_only, transversion,
                      ablate_bias, synthetic, multi_seq, htt, reseed_buffers)


@app.function(image=evo_image, gpu="A100", timeout=86400, volumes={"/results": vol})
def run_evo(model_name: str, commit: str, phase6_only: bool = False,
           transversion: bool = False, ablate_bias: bool = False,
                     synthetic: bool = False, multi_seq: bool = False,
                     htt: bool = False,
                     reseed_buffers: bool = False):
    return _run_model(model_name, commit, phase6_only, transversion,
                      ablate_bias, synthetic, multi_seq, htt, reseed_buffers)


@app.function(image=dnabert2_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_dnabert2(model_name: str, commit: str, phase6_only: bool = False,
                transversion: bool = False, ablate_bias: bool = False,
                     synthetic: bool = False, multi_seq: bool = False,
                     htt: bool = False,
                     reseed_buffers: bool = False):
    return _run_model(model_name, commit, phase6_only, transversion,
                      ablate_bias, synthetic, multi_seq, htt, reseed_buffers)


@app.function(image=caduceus_image, gpu="A10G", timeout=86400, volumes={"/results": vol})
def run_caduceus(model_name: str, commit: str, phase6_only: bool = False,
                transversion: bool = False, ablate_bias: bool = False,
                     synthetic: bool = False, multi_seq: bool = False,
                     htt: bool = False,
                     reseed_buffers: bool = False):
    return _run_model(model_name, commit, phase6_only, transversion,
                      ablate_bias, synthetic, multi_seq, htt, reseed_buffers)


def _route(model_name):
    """The (function, label) each model runs under, chosen by its adapter."""
    model_name = base_model(model_name)
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
    model_name = base_model(model_name)
    if model_name == "evo":
        return smoke_evo
    if model_name == "dnabert2":
        return smoke_dnabert2
    if model_name == "caduceus":
        return smoke_caduceus
    if model_name in LEGACY_MODELS:
        return (smoke_legacy_gpu if model_name in GPU_SMOKE_MODELS
                else smoke_legacy)
    return smoke_multimol


@app.local_entrypoint()
def main(models: str = "", phase6_only: bool = False, smoke_only: bool = False,
         transversion: bool = False, ablate_bias: bool = False,
         buffers_only: bool = False, synthetic: bool = False,
         multi_seq: bool = False, htt: bool = False,
         reseed_buffers: bool = False):
    commit = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                            capture_output=True, text=True, check=True).stdout.strip()
    porcelain = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain"],
                               capture_output=True, text=True, check=True).stdout
    # Outputs, not code. A fetched result directory says nothing about whether
    # the commit in the stamp contains the code the containers are about to run,
    # and leaving it in would mean no launch is possible once results land.
    dirty = "\n".join(line for line in porcelain.splitlines()
                       if not line[3:].startswith(("results/", "logs/")))
    if dirty.strip():
        raise SystemExit(
            "Refusing to launch from a dirty tree: the stamp would name a commit "
            "that does not contain the code being run.\n" + dirty)

    requested = [m.strip() for m in models.split(",") if m.strip()] or AVAILABLE_MODELS
    unknown = [m for m in requested if m not in AVAILABLE_MODELS]
    if unknown:
        raise SystemExit(f"Unknown models: {unknown}")

    if buffers_only:
        smoke_only = True
    if smoke_only:
        # A control loads the same adapter on the same image as its trained
        # counterpart, so smoking both twice buys nothing -- unless the thing
        # being checked only happens to a control. Collapsing the suffix then
        # removes the intervention from the check, which is how a reseeding run
        # reached the volume having reseeded nothing.
        if not reseed_buffers:
            requested = sorted({base_model(m) for m in requested},
                               key=lambda m: AVAILABLE_MODELS.index(m))
        print(f"Smoke check, {len(requested)} adapters: {requested}"
              + (", pairwise bias ablated" if ablate_bias else "")
              + (", buffers reseeded" if reseed_buffers else ""))
        # Spawned rather than called, so one broken adapter reports itself
        # alongside the nine that work instead of hiding them behind its own
        # traceback, and so the four images build concurrently.
        handles = {m: _smoke_route(m).spawn(m, ablate_bias=ablate_bias,
                                           buffers_only=buffers_only,
                                           reseed_buffers=reseed_buffers)
                   for m in requested}
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

    if transversion and phase6_only:
        raise SystemExit("--transversion runs the mutation stage; "
                         "--phase6-only skips it")

    if ablate_bias and any(base_model(m) != "ernierna" for m in requested):
        raise SystemExit(
            "--ablate-bias zeroes ERNIE-RNA's pairwise bias buffer and no other "
            f"model has one: {[m for m in requested if base_model(m) != 'ernierna']}")

    stage = "transversion control" if transversion else "re-run"
    if synthetic:
        stage += ", synthetic sequences"
    if multi_seq:
        stage += ", multi-sequence replication"
    if ablate_bias:
        stage += ", pairwise bias ablated"
    print(f"Repaired-panel {stage} at {commit[:12]}, "
          f"{len(requested)} models: {requested}")
    for model_name in requested:
        fn, label = _route(model_name)
        handle = fn.spawn(model_name=model_name, commit=commit,
                          phase6_only=phase6_only, transversion=transversion,
                          ablate_bias=ablate_bias, synthetic=synthetic,
                          multi_seq=multi_seq, htt=htt,
                          reseed_buffers=reseed_buffers)
        print(f"  {model_name:12s} {label:14s} {handle.object_id}")
    print("\nmodal app logs rna-repaired-panel")
