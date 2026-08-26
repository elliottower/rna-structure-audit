"""Modal wrapper: verify the nnsight envoy patch on Linux, on real RNA models.

Local verification on an Intel Mac is limited to torch 2.2.2 and cannot build
nnsight's C extension, so the patch is checked here against a proper install.

Three things are verified, in order of what would block the PR:

  1. nnsight's own test suite still passes (regression)
  2. the repro suite flips from 3 failed to 0 failed
  3. a real multimolecule RNA model produces a picklable envoy class and traces
     correctly, including the single-nucleotide perturbation the paper needs

Runs on CPU: every model here is tiny and the question is serialization, not
throughput.

Usage:
    modal run --detach scripts/ndif_feasibility/modal_verify_patch.py
    modal run --detach scripts/ndif_feasibility/modal_verify_patch.py --branch main
"""

import modal

app = modal.App("nnsight-envoy-patch-verify")

# torch/transformers/multimolecule pinned to the versions the RNA panel already
# runs under (see scripts/modal_rinalmo_phases15.py); git is needed for the
# editable install of the patched checkout.
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential")
    .pip_install(
        "torch==2.6.0",
        "numpy==2.1.3",
        "transformers==5.13.1",
        "multimolecule==0.1.0",
        "pytest==8.3.4",
        "tqdm==4.67.1",
        "einops==0.8.0",
        "setuptools>=61",
        "setuptools-scm>=8",
        # nnsight's own dependencies, listed here rather than resolved by the
        # editable install, so --no-deps can protect the torch/transformers pins
        "astor",
        "cloudpickle",
        "httpx",
        "python-socketio[client]",
        "python-engineio>=4.13.0",
        "pydantic>=2.9.0",
        "accelerate",
        "toml",
        "ipython",
        "rich",
        "zstandard",
        "pytest-cov",
    )
    .add_local_dir(
        "../nnsight", "/root/nnsight", copy=True,
        ignore=["**/__pycache__", "**/*.pyc"],
    )
    .add_local_file(
        "scripts/ndif_feasibility/nnsight_repro/test_envoy_overloaded_mount.py",
        "/root/repro/test_envoy_overloaded_mount.py",
        copy=True,
    )
    .add_local_file(
        "scripts/ndif_feasibility/_rna_envoy_check.py",
        "/root/repro/_rna_envoy_check.py",
        copy=True,
    )
)

vol = modal.Volume.from_name("nnsight-patch-verify", create_if_missing=True)

RNA_MODELS = [
    # (repo, multimolecule model class, tokenizer repo)
    ("multimolecule/rnabert", "RnaBertModel"),
    ("multimolecule/rnafm", "RnaFmModel"),
    ("multimolecule/splicebert", "SpliceBertModel"),
]

WT = "GGCUAGCUAAGGCUAGCC"
MUT = "GGCUAGCUAUGGCUAGCC"  # single substitution at position 10


@app.function(image=image, timeout=86400, volumes={"/results": vol}, cpu=4.0)
def verify(branches: str = "main,verify/both-fixes"):
    """Install once, then run the same checks on each branch and compare.

    `main` is the negative control: it must fail. A patch that passes on the
    fix branch tells you nothing unless the same suite fails without it.
    """
    import json
    import subprocess
    from datetime import datetime, timezone
    from pathlib import Path

    def stamp() -> str:
        return datetime.now(timezone.utc).strftime("%H:%M:%S")

    def run(cmd, cwd=None):
        print(f"[{stamp()}] $ {' '.join(cmd)}")
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
        tail = (p.stdout + p.stderr).strip().splitlines()
        for line in tail[-12:]:
            print(f"    {line}")
        return {"returncode": p.returncode, "tail": tail[-30:]}

    out = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "branches": {},
    }
    out_path = Path("/results") / "verify_branches.json"

    def save():
        out_path.write_text(json.dumps(out, indent=2))
        vol.commit()

    # editable install once; later branch switches change the .py files in place
    out["install"] = run(
        ["pip", "install", "-e", ".", "--no-build-isolation", "--no-deps"],
        cwd="/root/nnsight",
    )
    save()

    for branch in [b.strip() for b in branches.split(",") if b.strip()]:
        print(f"\n[{stamp()}] ===== {branch} =====")
        entry = {}
        entry["checkout"] = run(["git", "checkout", "-f", branch], cwd="/root/nnsight")
        entry["sha"] = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd="/root/nnsight", capture_output=True, text=True,
        ).stdout.strip()

        entry["nnsight_tests"] = run(
            ["python", "-m", "pytest", "tests/test_tiny.py", "--device", "cpu", "-q",
             "-p", "no:cacheprovider"],
            cwd="/root/nnsight",
        )
        entry["repro_tests"] = run(
            ["python", "-m", "pytest", "test_envoy_overloaded_mount.py", "-q",
             "-p", "no:cacheprovider"],
            cwd="/root/repro",
        )

        rna_out = f"/results/rna_{branch.replace('/', '_')}.json"
        entry["rna_check"] = run(
            ["python", "_rna_envoy_check.py", "--out", rna_out], cwd="/root/repro"
        )
        try:
            entry["rna"] = json.loads(Path(rna_out).read_text())
        except Exception as exc:
            entry["rna"] = {"error": str(exc)}

        out["branches"][branch] = entry
        save()

    out["finished_utc"] = datetime.now(timezone.utc).isoformat()
    save()

    print(f"\n[{stamp()}] ===== summary =====")
    for branch, entry in out["branches"].items():
        print(f"  {branch} ({entry.get('sha')})")
        for k in ("nnsight_tests", "repro_tests", "rna_check"):
            print(f"    {k:<16} rc={entry.get(k, {}).get('returncode')}")
        for repo, m in (entry.get("rna", {}).get("models", {}) or {}).items():
            print(f"      {repo:<32} class={m.get('envoy_class')} "
                  f"picklable={m.get('class_picklable')}")
    return out


@app.local_entrypoint()
def main(branches: str = "main,verify/both-fixes"):
    result = verify.remote(branches=branches)
    print(json.dumps(result, indent=2)[:4000])


import json  # noqa: E402  -- used by the local entrypoint
