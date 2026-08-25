"""Is ERNIE-RNA's pairwise bias a parameter or a buffer?

`_randomize` in `modal_repaired_panel.py` iterates `named_parameters()`, so a
buffer survives it untouched. If `pairwise_bias_map` is a buffer, the randomly
initialized ERNIE-RNA control keeps the Watson-Crick pairing prior it encodes,
and is therefore not a null model for a partner-specificity test -- which is what
its per-pair precision of 0.476 against a chance rate of 0.241 suggests.

If it is a parameter, randomization destroys it and that explanation is wrong.

The image is copied from `modal_repaired_panel.py`, pin for pin. CPU only: this
loads a model and reads its attribute table, and never runs a forward pass.

    modal run scripts/modal_inspect_ernierna_bias.py
"""

import modal

app = modal.App("rna-inspect-ernierna-bias")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.6.0",
        "transformers==5.14.1",
        "multimolecule==0.2.0",
        "torchmetrics==1.4.1",
        "numpy==1.26.4",
        "scipy==1.14.1",
        "tqdm==4.66.5",
        "matplotlib==3.9.2",
        "scikit-learn==1.5.2",
        "einops==0.8.0",
    )
    .env({"PYTHONPATH": "/root/project:/root/project/scripts"})
    .add_local_file("multi_model_audit.py", "/root/project/multi_model_audit.py")
    .add_local_dir("scripts", "/root/project/scripts")
)


@app.function(image=image, timeout=1800)
def inspect():
    import os
    os.chdir("/root/project")

    import torch
    from multi_model_audit import ERNIERNAAdapter

    adapter = ERNIERNAAdapter()
    adapter.load()
    model = adapter.model

    params = dict(model.named_parameters())
    buffers = dict(model.named_buffers())

    print(f"parameters: {len(params)}   buffers: {len(buffers)}")
    for label, table in (("PARAMETER", params), ("BUFFER", buffers)):
        for name, tensor in table.items():
            if "pairwise" in name or "bias_map" in name:
                finite = torch.isfinite(tensor)
                print(f"  {label:9s} {name}  shape={list(tensor.shape)}  "
                      f"finite={int(finite.sum())}/{tensor.numel()}  "
                      f"nan={int(torch.isnan(tensor).sum())}  "
                      f"inf={int(torch.isinf(tensor).sum())}")
                if int(finite.sum()):
                    print(f"            finite range "
                          f"[{tensor[finite].min():.4f}, {tensor[finite].max():.4f}]")

    # Which cells hold the damage, and whether a nucleotide pair can reach them.
    # The panel is pure ACGU, so a NaN in a row for an ambiguity code or a
    # special token is never indexed and never matters; one in an ACGU row
    # saturates or poisons attention on real sequences.
    bias = model.pairwise_bias_map
    vocab = adapter.tokenizer.get_vocab()
    by_id = {index: token for token, index in vocab.items()}
    nucleotides = {token: index for token, index in vocab.items()
                   if token in ("A", "C", "G", "U")}
    print(f"\n  nucleotide token ids: {nucleotides}")

    bad = (~torch.isfinite(bias)) | (bias.abs() > 1e6)
    rows, cols = torch.nonzero(bad, as_tuple=True)
    print(f"  {len(rows)} cells are NaN or above 1e6:")
    reachable = 0
    for r, c in zip(rows.tolist(), cols.tolist()):
        a, b = by_id.get(r, f"<{r}>"), by_id.get(c, f"<{c}>")
        hit = a in nucleotides and b in nucleotides
        reachable += hit
        print(f"    [{r:2d},{c:2d}] {a!r} x {b!r}  value={bias[r, c].item()}"
              f"{'   <-- REACHABLE from an ACGU sequence' if hit else ''}")
    print(f"\n  cells reachable from an ACGU pair: {reachable}")

    acgu = sorted(nucleotides.values())
    block = bias[acgu][:, acgu]
    print(f"  the 4x4 ACGU block: finite={int(torch.isfinite(block).sum())}/16, "
          f"range [{block[torch.isfinite(block)].min():.4f}, "
          f"{block[torch.isfinite(block)].max():.4f}]")
    for i, ra in zip(acgu, [by_id[i] for i in acgu]):
        print("    " + ra + ": " + "  ".join(
            f"{by_id[j]}={bias[i, j].item():+.3f}" for j in acgu))

    # --- Is the buffer in the checkpoint, or is it whatever memory held? ---
    import hashlib
    import inspect as inspect_mod

    print("\n  --- all non-persistent buffers, full names ---")
    persistent_names = set(model.state_dict())
    for name, tensor in model.named_buffers():
        if name in persistent_names:
            continue
        print(f"    {name:52s} shape={list(tensor.shape)} "
              f"nonzero={int((tensor != 0).sum())}/{tensor.numel()}")

    # A non-persistent buffer is never saved to or loaded from a checkpoint, so
    # whatever it holds was computed in `__init__` rather than trained. A table
    # of sines lies in [-1, 1] and has a mean near zero; a learned embedding
    # table would not be confined that way.
    for name, tensor in list(model.named_buffers()):
        if name in persistent_names or tensor.dim() != 2 or tensor.shape[0] < 100:
            continue
        print(f"\n    {name}: row 5, first 8 values")
        print(f"      {[round(float(v), 4) for v in tensor[5].detach().cpu()[:8]]}")
        print(f"      min={float(tensor.min()):.4f} max={float(tensor.max()):.4f} "
              f"mean={float(tensor.mean()):.4f} std={float(tensor.std()):.4f}")
        print(f"      inside [-1, 1]: "
              f"{int(((tensor >= -1) & (tensor <= 1)).sum())}/{tensor.numel()}")

    print("\n  --- provenance of the buffer ---")
    in_state = "pairwise_bias_map" in model.state_dict()
    print(f"  in model.state_dict() (i.e. registered persistent): {in_state}")

    from huggingface_hub import hf_hub_download
    stored_keys = None
    for filename in ("model.safetensors", "pytorch_model.bin"):
        try:
            path = hf_hub_download("multimolecule/ernierna", filename)
        except Exception:
            continue
        if filename.endswith(".safetensors"):
            from safetensors.torch import load_file
            stored_keys = set(load_file(path).keys())
        else:
            stored_keys = set(torch.load(path, map_location="cpu",
                                         weights_only=False).keys())
        print(f"  checkpoint read: {filename}, {len(stored_keys)} tensors")
        break
    if stored_keys is not None:
        hits = sorted(k for k in stored_keys if "pairwise" in k)
        print(f"  checkpoint keys matching 'pairwise': {hits or 'NONE'}")

    # --- Is it read in the forward path, and filled in __init__? ---
    for label, fn in (("__init__", type(model).__init__),
                      ("forward", type(model).forward)):
        try:
            src = inspect_mod.getsource(fn)
        except Exception as exc:
            print(f"  {label}: source unavailable ({exc})")
            continue
        lines = [line.strip() for line in src.splitlines()
                 if "pairwise_bias_map" in line]
        print(f"  {label}: {len(lines)} line(s) mention pairwise_bias_map")
        for line in lines[:6]:
            print(f"      {line}")

    # --- Is it the same tensor on a second instantiation? ---
    second = ERNIERNAAdapter()
    second.load()
    def digest(tensor):
        return hashlib.sha256(
            tensor.detach().cpu().float().nan_to_num(0.0).numpy().tobytes()
        ).hexdigest()[:16]
    a, b = digest(bias), digest(second.model.pairwise_bias_map)
    print(f"\n  buffer sha256 (first load):  {a}")
    print(f"  buffer sha256 (second load): {b}")
    print(f"  identical across instantiations: {a == b}")

    # --- Is it read anywhere at all, and how is it built? ---
    import multimolecule.models.ernierna as ernierna_mod
    from pathlib import Path as _Path

    print("\n  --- every mention in the ernierna package ---")
    root = _Path(inspect_mod.getfile(ernierna_mod)).parent
    for source_file in sorted(root.rglob("*.py")):
        for number, line in enumerate(source_file.read_text().splitlines(), 1):
            if "pairwise_bias_map" in line:
                print(f"    {source_file.name}:{number}  {line.strip()}")

    # --- Which other models rebuild a non-persistent buffer on first forward? ---
    print("\n  --- the _inited re-registration pattern across multimolecule ---")
    package = root.parent.parent
    hits = {}
    for source_file in sorted(package.rglob("*.py")):
        lines = source_file.read_text().splitlines()
        for number, line in enumerate(lines, 1):
            if "_inited" in line or ("register_buffer" in line
                                     and "persistent=False" in line):
                hits.setdefault(source_file.relative_to(package).as_posix(),
                                []).append((number, line.strip()))
    for name, found in hits.items():
        print(f"    {name}")
        for number, line in found[:4]:
            print(f"      {number}: {line}")
    if not hits:
        print("    none")

    print("\n  --- _build_pairwise_bias_map, and its call site ---")
    for source_file in sorted(root.rglob("modeling_ernierna.py")):
        lines = source_file.read_text().splitlines()
        for start, end, label in ((129, 155, "call site"), (1414, 1450, "builder")):
            print(f"\n    [{label}] {source_file.name}:{start + 1}-{end}")
            for number in range(start, min(end, len(lines))):
                print(f"    {number + 1:5d}  {lines[number]}")

    # --- The same buffer, after one forward pass has run ---
    print("\n  --- after a forward pass ---")
    tokens = adapter.tokenize("ACGUACGUACGUACGU")
    with torch.no_grad():
        model(tokens)
    after = model.pairwise_bias_map
    print(f"  nan={int(torch.isnan(after).sum())}  "
          f"above 1e6={int((after.abs() > 1e6).sum())}  "
          f"nonzero={int((after != 0).sum())}/{after.numel()}")
    for i in sorted(nucleotides.values()):
        print("    " + by_id[i] + ": " + "  ".join(
            f"{by_id[j]}={after[i, j].item():+.3f}"
            for j in sorted(nucleotides.values())))

    second.model(second.adapter_tokens if False else tokens)
    print(f"  identical to a second model's, post-forward: "
          f"{digest(after) == digest(second.model.pairwise_bias_map)}")

    attr = getattr(model, "pairwise_bias_map", None)
    if attr is not None:
        kind = ("Parameter" if isinstance(attr, torch.nn.Parameter)
                else type(attr).__name__)
        survives = kind != "Parameter"
        print(f"\n  model.pairwise_bias_map is a {kind}")
        print(f"  survives _randomize (which iterates named_parameters): {survives}")
    else:
        print("\n  model has no attribute pairwise_bias_map")
    return "done"


@app.local_entrypoint()
def main():
    print(inspect.remote())
