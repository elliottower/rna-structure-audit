"""Patch DNABERT-2's flash_attn_triton.py to fix trans_b deprecation.

MosaicBERT's Triton flash attention uses tl.dot(q, k, trans_b=True),
which was removed in Triton 2.0+. This script replaces it with
tl.dot(q, tl.trans(k)).

Run AFTER downloading the DNABERT-2 model to patch the cached files.
"""

import glob
import re
import subprocess
import sys


def main():
    from huggingface_hub import snapshot_download
    snapshot_download("zhihan1996/DNABERT-2-117M")

    from transformers import AutoConfig
    AutoConfig.from_pretrained(
        "zhihan1996/DNABERT-2-117M", trust_remote_code=True,
    )

    patterns = [
        "/root/.cache/huggingface/modules/transformers_modules/"
        "zhihan1996/DNABERT-2-117M/*/flash_attn_triton.py",
        "/root/.cache/huggingface/hub/models--zhihan1996--DNABERT-2-117M/"
        "snapshots/*/flash_attn_triton.py",
    ]
    files = []
    for p in patterns:
        files.extend(glob.glob(p))

    if not files:
        result = subprocess.run(
            ["find", "/root/.cache/huggingface", "-name", "flash_attn_triton.py"],
            capture_output=True, text=True,
        )
        files = [f.strip() for f in result.stdout.strip().split("\n") if f.strip()]

    print(f"Found {len(files)} flash_attn_triton.py files: {files}")

    regex = re.compile(r"tl\.dot\((\w+),\s*(\w+),\s*trans_b\s*=\s*True\)")
    patched_count = 0
    for fpath in files:
        src = open(fpath).read()
        patched = regex.sub(r"tl.dot(\1, tl.trans(\2))", src)
        if patched != src:
            open(fpath, "w").write(patched)
            patched_count += 1
            print(f"Patched: {fpath}")
        else:
            print(f"No match: {fpath}")

    print(f"Patched {patched_count}/{len(files)} files")
    if patched_count == 0 and files:
        print("WARNING: files found but no trans_b pattern matched")
        for f in files:
            content = open(f).read()
            if "trans_b" in content:
                print(f"  trans_b IS in {f} but regex didn't match")
                for i, line in enumerate(content.split("\n")):
                    if "trans_b" in line:
                        print(f"    Line {i+1}: {line.strip()}")
            else:
                print(f"  trans_b NOT in {f}")
    sys.exit(0 if patched_count > 0 or not files else 1)


if __name__ == "__main__":
    main()
