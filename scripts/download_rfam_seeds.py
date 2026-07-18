"""Download Rfam seed alignments and extract multiple sequences per family.

For each of the 52 Rfam families in data/rfam_families/, download the seed
alignment from Rfam's REST API, parse Stockholm format, and extract up to
5 sequences with consensus secondary structure.

Usage:
    cd /path/to/causal-rna
    uv run --with requests --with tqdm python scripts/download_rfam_seeds.py
"""

import json
import os
import re
import time
from pathlib import Path

import requests
from tqdm import tqdm

DATA_DIR = Path("/Users/elliottower/Documents/GitHub/causal-rna/data/rfam_families")
SEED_DIR = Path("/Users/elliottower/Documents/GitHub/causal-rna/data/rfam_seeds")
OUTPUT_DIR = Path("/Users/elliottower/Documents/GitHub/causal-rna/data/multi_sequence")
MAX_SEQS = 5
MAX_SEQ_LEN = 512


def load_families() -> list[dict]:
    families = []
    for f in sorted(DATA_DIR.iterdir()):
        if f.suffix == ".json":
            with open(f) as fh:
                families.append(json.load(fh))
    return families


def download_seed_alignment(rfam_id: str, out_path: Path) -> bool:
    """Download seed alignment in Stockholm format from Rfam."""
    url = f"https://rfam.org/family/{rfam_id}/alignment?acc={rfam_id}&format=stockholm&download=1"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            with open(out_path, "w") as f:
                f.write(resp.text)
            return True
        print(f"  HTTP {resp.status_code} for {rfam_id}")
        return False
    except requests.RequestException as e:
        print(f"  Error downloading {rfam_id}: {e}")
        return False


def parse_stockholm(sto_path: Path) -> list[dict]:
    """Parse Stockholm alignment file, extract sequences and consensus structure.

    Returns list of dicts with keys: accession, sequence, dot_bracket.
    """
    sequences = {}
    gc_ss_cons = ""

    with open(sto_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#=GC SS_cons"):
                gc_ss_cons += line.split(None, 2)[2]
            elif line.startswith("#") or line.startswith("//") or not line.strip():
                continue
            else:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    name, seq = parts
                    if name in sequences:
                        sequences[name] += seq
                    else:
                        sequences[name] = seq

    if not gc_ss_cons:
        return []

    ss_clean = gc_ss_cons.replace(":", ".").replace(",", ".").replace("_", ".")
    ss_clean = re.sub(r"[<{(\[]", "(", ss_clean)
    ss_clean = re.sub(r"[>})\]]", ")", ss_clean)
    ss_clean = re.sub(r"[^().]", ".", ss_clean)

    results = []
    for name, aligned_seq in sequences.items():
        ungapped_seq = ""
        ungapped_ss = ""
        for i, (s, ss) in enumerate(zip(aligned_seq, ss_clean)):
            if s not in ".-~":
                ungapped_seq += s.upper().replace("T", "U")
                ungapped_ss += ss
        if len(ungapped_seq) != len(ungapped_ss):
            continue
        if len(ungapped_seq) < 20 or len(ungapped_seq) > MAX_SEQ_LEN:
            continue
        if not _valid_structure(ungapped_ss):
            continue
        results.append({
            "accession": name,
            "sequence": ungapped_seq.replace("U", "T"),
            "dot_bracket": ungapped_ss,
            "length": len(ungapped_seq),
        })
    return results


def _valid_structure(ss: str) -> bool:
    """Check that dot-bracket has balanced parentheses and at least 3 stems."""
    depth = 0
    n_pairs = 0
    for c in ss:
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            n_pairs += 1
        if depth < 0:
            return False
    return depth == 0 and n_pairs >= 3


def select_diverse_sequences(seqs: list[dict], original_seq: str, n: int = MAX_SEQS) -> list[dict]:
    """Select up to n sequences, preferring those most different from original."""
    if len(seqs) <= n:
        return seqs

    def seq_distance(s):
        a = s["sequence"]
        min_len = min(len(a), len(original_seq))
        if min_len == 0:
            return 0
        mismatches = sum(1 for i in range(min_len) if a[i] != original_seq[i])
        return mismatches / min_len

    scored = [(seq_distance(s), s) for s in seqs]
    scored.sort(key=lambda x: -x[0])
    return [s for _, s in scored[:n]]


def main():
    SEED_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    families = load_families()
    print(f"Loaded {len(families)} families")

    seen_rfam_ids = {}
    for fam in families:
        rfam_id = fam.get("rfam_id")
        if rfam_id and rfam_id != "None":
            if rfam_id not in seen_rfam_ids:
                seen_rfam_ids[rfam_id] = []
            seen_rfam_ids[rfam_id].append(fam["name"])

    unique_ids = set(seen_rfam_ids.keys())
    print(f"Unique Rfam IDs: {len(unique_ids)}")
    for rfam_id, names in seen_rfam_ids.items():
        if len(names) > 1:
            print(f"  {rfam_id} shared by: {', '.join(names)}")

    downloaded = 0
    for rfam_id in tqdm(sorted(unique_ids), desc="Downloading seeds"):
        sto_path = SEED_DIR / f"{rfam_id}.sto"
        if sto_path.exists():
            continue
        if download_seed_alignment(rfam_id, sto_path):
            downloaded += 1
        time.sleep(1)

    print(f"Downloaded {downloaded} new alignments")

    summary = {}
    for fam in tqdm(families, desc="Extracting sequences"):
        rfam_id = fam.get("rfam_id")
        if not rfam_id or rfam_id == "None":
            print(f"  {fam['name']}: no rfam_id, skipping")
            continue

        sto_path = SEED_DIR / f"{rfam_id}.sto"
        if not sto_path.exists():
            print(f"  {fam['name']}: no alignment file for {rfam_id}")
            continue

        all_seqs = parse_stockholm(sto_path)
        if not all_seqs:
            print(f"  {fam['name']}: no valid sequences from {rfam_id}")
            continue

        selected = select_diverse_sequences(all_seqs, fam["sequence"], MAX_SEQS)

        out_data = {
            "name": fam["name"],
            "rfam_id": rfam_id,
            "original_sequence": fam["sequence"],
            "original_dot_bracket": fam["dot_bracket"],
            "n_available": len(all_seqs),
            "n_selected": len(selected),
            "sequences": selected,
        }

        out_path = OUTPUT_DIR / f"{fam['name']}_multi.json"
        with open(out_path, "w") as f:
            json.dump(out_data, f, indent=2)

        summary[fam["name"]] = {
            "rfam_id": rfam_id,
            "n_available": len(all_seqs),
            "n_selected": len(selected),
        }

    with open(OUTPUT_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    n_with_multi = sum(1 for v in summary.values() if v["n_selected"] >= 2)
    print(f"\nSummary: {len(summary)} families processed, {n_with_multi} with 2+ sequences")
    print(f"Results saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
