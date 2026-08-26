"""Multi-sequence perturbation specificity computation.

Runs PS on all valid sequences from Rfam seed alignments for a given model.
Pure Python logic — no Modal dependency. Called by modal_multi_seq_ps.py.

Each family contributes many sequences (from its seed alignment). PS is
computed per-sequence, then aggregated per-family (mean +/- SE) and globally.
"""

import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from parse_stockholm import extract_sequences
from phase6_compensatory_mutation import (
    compute_delta_profiles,
    compute_ps_from_deltas,
    derangement_null,
    get_eligible_pairs,
    h3_precision_test,
    parse_stems,
    positive_control,
)

MIN_WC_PAIRS = 15
MIN_ELIGIBLE = 5
MIN_SEQ_LENGTH = 30


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def process_single_sequence(adapter, seq_str, db_str, device="cpu", offset=0,
                            compute_null=True):
    """Run PS pipeline on a single (sequence, dot_bracket) pair.

    Returns result dict or None if sequence doesn't meet eligibility criteria.
    """
    if len(seq_str) != len(db_str):
        return None

    stems = parse_stems(db_str, seq_str)
    n_wc_total = sum(len(s) for s in stems)

    if n_wc_total < MIN_WC_PAIRS:
        return {"skipped": True, "reason": f"< {MIN_WC_PAIRS} WC pairs ({n_wc_total})"}

    eligible = get_eligible_pairs(seq_str, stems)
    if len(eligible) < MIN_ELIGIBLE:
        return {"skipped": True, "reason": f"< {MIN_ELIGIBLE} eligible pairs ({len(eligible)})"}

    delta_profiles, stem_pos, loop_pos, n_layers = compute_delta_profiles(
        adapter, seq_str, eligible, stems, device=device, offset=offset,
    )

    ps_result = compute_ps_from_deltas(eligible, delta_profiles, n_layers)
    if ps_result is None:
        return {"skipped": True, "reason": "no valid PS"}

    best_layer = ps_result["best_layer"]
    pc = positive_control(eligible, delta_profiles, best_layer, stem_pos, loop_pos)

    null_result = None
    if compute_null:
        null_result = derangement_null(eligible, delta_profiles, n_layers, best_layer,
                                       n_derangements=500)

    h3 = h3_precision_test(ps_result["pair_details"], eligible)

    result = {
        "best_ps": ps_result["best_ps"],
        "best_layer": ps_result["best_layer"],
        "n_eligible_pairs": len(eligible),
        "n_wc_pairs_total": n_wc_total,
        "n_stems": len([s for s in stems if len(s) >= 3]),
        "positive_control": pc,
        "h3_precision": h3,
    }

    if null_result:
        result["null_95th_primary"] = null_result["null_95th_primary"]
        result["null_available"] = null_result.get("null_available", True)
        if null_result["null_95th_primary"] is not None:
            result["exceeds_null_primary"] = ps_result["best_ps"] > null_result["null_95th_primary"]
        else:
            result["exceeds_null_primary"] = None

    return result


def run_multi_seq_ps(adapter, seed_dir, device="cpu", offset=0, compute_null=True,
                     max_seqs_per_family=None):
    """Run PS on all valid sequences from all seed alignments.

    Returns structured results with per-family and global aggregation.
    """
    seed_path = Path(seed_dir)
    sto_files = sorted(seed_path.glob("*.sto"))

    all_results = {}
    global_ps_values = []

    for sto_file in tqdm(sto_files, desc=f"Families [{adapter.name}]"):
        fam_id, fam_name, sequences, parse_stats = extract_sequences(
            sto_file, min_length=MIN_SEQ_LENGTH
        )

        if max_seqs_per_family and len(sequences) > max_seqs_per_family:
            np.random.shuffle(sequences)
            sequences = sequences[:max_seqs_per_family]

        print(f"[{timestamp()}] {fam_id} ({fam_name}): {len(sequences)} valid sequences")

        family_ps_values = []
        family_details = []
        n_skipped = 0
        n_failed_gate = 0

        for seq_entry in tqdm(sequences, desc=f"  {fam_id}", leave=False):
            result = process_single_sequence(
                adapter, seq_entry["sequence"], seq_entry["dot_bracket"],
                device=device, offset=offset, compute_null=compute_null,
            )

            if result is None or result.get("skipped"):
                n_skipped += 1
                continue

            if not result.get("positive_control", {}).get("pass", False):
                n_failed_gate += 1
                continue

            family_ps_values.append(result["best_ps"])
            family_details.append({
                "accession": seq_entry["accession"],
                "length": seq_entry["length"],
                "best_ps": result["best_ps"],
                "best_layer": result["best_layer"],
                "n_eligible_pairs": result["n_eligible_pairs"],
                "exceeds_null": result.get("exceeds_null_primary"),
                "h3_fraction": result.get("h3_precision", {}).get("fraction"),
            })

        family_entry = {
            "family_id": fam_id,
            "family_name": fam_name,
            "n_total_valid": parse_stats["valid"],
            "n_processed": len(sequences),
            "n_skipped": n_skipped,
            "n_failed_gate": n_failed_gate,
            "n_usable": len(family_ps_values),
        }

        if family_ps_values:
            ps_arr = np.array(family_ps_values)
            family_entry["mean_ps"] = float(np.mean(ps_arr))
            family_entry["median_ps"] = float(np.median(ps_arr))
            family_entry["std_ps"] = float(np.std(ps_arr, ddof=1)) if len(ps_arr) > 1 else 0.0
            family_entry["se_ps"] = family_entry["std_ps"] / np.sqrt(len(ps_arr))
            family_entry["min_ps"] = float(np.min(ps_arr))
            family_entry["max_ps"] = float(np.max(ps_arr))
            family_entry["per_sequence"] = family_details
            global_ps_values.extend(family_ps_values)

        all_results[fam_id] = family_entry
        if family_ps_values:
            print(f"  -> {len(family_ps_values)} usable, mean PS = {family_entry['mean_ps']:.4f}")
        else:
            print(f"  -> 0 usable")

    summary = {
        "n_families": len(all_results),
        "n_families_with_data": sum(1 for v in all_results.values() if v.get("mean_ps") is not None),
        "n_total_sequences": sum(v["n_usable"] for v in all_results.values()),
    }

    if global_ps_values:
        ps_global = np.array(global_ps_values)
        summary["global_mean_ps"] = float(np.mean(ps_global))
        summary["global_median_ps"] = float(np.median(ps_global))
        summary["global_std_ps"] = float(np.std(ps_global, ddof=1))
        summary["global_se_ps"] = summary["global_std_ps"] / np.sqrt(len(ps_global))

        family_means = [v["mean_ps"] for v in all_results.values() if v.get("mean_ps") is not None]
        if family_means:
            fm = np.array(family_means)
            summary["family_mean_of_means"] = float(np.mean(fm))
            summary["family_se_of_means"] = float(np.std(fm, ddof=1) / np.sqrt(len(fm))) if len(fm) > 1 else 0.0

    return {"summary": summary, "per_family": all_results}
