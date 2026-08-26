"""Do the multimolecule models' non-persistent buffers hold real values?

`transformers` v5 leaves a non-persistent buffer on the meta device after
`from_pretrained`, so it materializes as whatever memory held. ERNIE-RNA works
around this by re-registering the buffer on the first forward pass, guarded by an
`_inited` flag. No other model in this panel has that guard.

A buffer that is never rebuilt is either harmless -- a zero-element dtype
reference, or one the forward pass never reads -- or it is garbage the model
computes with on every sequence. This reports, per model, every buffer with its
finiteness before and after a forward pass, so the difference between the two
cases is visible rather than assumed.

    modal run scripts/modal_inspect_multimol_buffers.py
"""

import modal

app = modal.App("rna-inspect-multimol-buffers")

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

MODELS = ["rinalmo", "utrlm", "splicebert", "ernierna"]


@app.function(image=image, timeout=3600)
def inspect():
    import os
    os.chdir("/root/project")

    import torch
    from phase6_compensatory_mutation import load_adapter

    for key in MODELS:
        adapter = load_adapter(key)
        adapter.load()
        model = adapter.model
        persistent = set(model.state_dict())

        print(f"\n=== {key} ===")
        rows = []
        for name, tensor in model.named_buffers():
            if name in persistent:
                continue
            rows.append((name, tensor))
        if not rows:
            print("  no non-persistent buffers")
            continue

        def state(tensor):
            if tensor is None or tensor.numel() == 0:
                return "empty"
            nan = int(torch.isnan(tensor).sum())
            big = int((tensor.abs() > 1e6).sum())
            return (f"n={tensor.numel()} nan={nan} >1e6={big} "
                    f"{'GARBAGE' if (nan or big) else 'clean'}")

        before = {name: state(tensor) for name, tensor in rows}
        with torch.no_grad():
            model(adapter.tokenize("ACGUACGUACGUACGUACGU"))

        after = dict(model.named_buffers())
        for name, _ in rows:
            tensor = after.get(name)
            verdict = state(tensor)
            flag = ""
            if "GARBAGE" in before[name]:
                flag = ("  <-- REPAIRED on first forward"
                        if "GARBAGE" not in verdict
                        else "  <-- STILL GARBAGE after forward")
            print(f"  {name:44s} before[{before[name]}] after[{verdict}]{flag}")
    return "done"


@app.local_entrypoint()
def main():
    print(inspect.remote())
