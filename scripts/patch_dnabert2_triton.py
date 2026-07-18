"""Patch DNABERT-2's flash_attn_triton.py to work with triton >= 2.0.0.

The trans_b and trans_a kwargs to tl.dot() were removed in triton 2.0.0.
The replacement is tl.trans() wrapping the relevant argument.

Run after snapshot_download("zhihan1996/DNABERT-2-117M").
"""

import glob
import re

files = glob.glob(
    "/root/.cache/huggingface/hub/models--zhihan1996--DNABERT-2-117M/snapshots/*/flash_attn_triton.py"
)
print(f"Found {len(files)} flash_attn_triton.py files to patch")

for f in files:
    code = open(f).read()
    orig = code

    # tl.dot(X, Y, trans_b=True) -> tl.dot(X, tl.trans(Y))
    code = re.sub(
        r"tl\.dot\(([^,]+),\s*([^,)]+),\s*trans_b\s*=\s*True\)",
        lambda m: f"tl.dot({m.group(1)}, tl.trans({m.group(2)}))",
        code,
    )

    # tl.dot(X, Y, trans_a=True) -> tl.dot(tl.trans(X), Y)
    code = re.sub(
        r"tl\.dot\(([^,]+),\s*([^,)]+),\s*trans_a\s*=\s*True\)",
        lambda m: f"tl.dot(tl.trans({m.group(1)}), {m.group(2)})",
        code,
    )

    # Remove remaining false flags
    code = code.replace(", trans_a=False", "").replace(", trans_b=False", "")

    n_changed = sum(1 for a, b in zip(orig, code) if a != b)
    open(f, "w").write(code)
    print(f"  Patched {f} ({n_changed} chars changed)")
