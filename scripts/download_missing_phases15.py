"""Download Phase 1-5 results from Modal volume to local disk.

Run after modal_missing_phases15.py completes.

Usage:
    cd /path/to/causal-rna
    uv run python scripts/download_missing_phases15.py
"""

import json
import shutil
import subprocess
from pathlib import Path

VOLUME_NAME = "causal-rna-missing-phases15"
LOCAL_DEST = Path("data/gpu_results/expanded_rfam")
RNA_AWARENESS_DEST = Path(__file__).resolve().parent.parent.parent / "rna-structure-awareness" / "data" / "expanded_rfam"

EXPECTED_MODELS = ["rnafm", "hyenadna", "evo", "nt"]


def main():
    LOCAL_DEST.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        ["modal", "volume", "ls", VOLUME_NAME],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"Error listing volume: {result.stderr}")
        return

    print(f"Volume contents:\n{result.stdout}")

    dirs = [line.strip().rstrip("/") for line in result.stdout.strip().split("\n")
            if line.strip() and "phases15" in line]

    if not dirs:
        print("No Phase 1-5 result directories found in volume.")
        return

    for d in dirs:
        print(f"\nDownloading {d}...")
        dl_result = subprocess.run(
            ["modal", "volume", "get", VOLUME_NAME, d, str(LOCAL_DEST / d)],
            capture_output=True, text=True,
        )
        if dl_result.returncode != 0:
            print(f"  Error: {dl_result.stderr}")
            continue

        json_files = list((LOCAL_DEST / d).glob("*phases_1_to_5.json"))
        for jf in json_files:
            with open(jf) as f:
                data = json.load(f)
            model = data.get("model", "unknown")
            ts = data.get("timestamp", "unknown")
            dest_name = f"{model}_expanded_{ts}.json"

            dest_causal = LOCAL_DEST / dest_name
            shutil.copy2(jf, dest_causal)
            print(f"  Copied to {dest_causal}")

            if RNA_AWARENESS_DEST.exists():
                dest_rna = RNA_AWARENESS_DEST / dest_name
                shutil.copy2(jf, dest_rna)
                print(f"  Copied to {dest_rna}")

    print("\nDone. Verify values match paper Table 1:")
    for model in EXPECTED_MODELS:
        matches = list(LOCAL_DEST.glob(f"{model}_expanded_*.json"))
        if matches:
            with open(matches[-1]) as f:
                data = json.load(f)
            mt = data.get("mutation_trained", {})
            ratio = mt.get("mean_best_ratio", "N/A")
            exceed = mt.get("families_exceeding_nuc_null", mt.get("n_exceeding_nuc_null", "N/A"))
            total = mt.get("families_total", mt.get("n_families_scored", "N/A"))
            print(f"  {model}: mean_ratio={ratio}, exceed={exceed}/{total}")
        else:
            print(f"  {model}: NO FILE FOUND")


if __name__ == "__main__":
    main()
