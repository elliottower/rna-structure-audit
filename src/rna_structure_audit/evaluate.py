"""Main evaluation entry point."""

import json
from pathlib import Path

from rna_structure_audit.data import load_families
from rna_structure_audit.rungs import run_rung1, run_rung2, run_rung3


def _grade(rung1, rung2, rung3):
    n_fam = rung1["n_families_scored"]
    n_r1 = rung1["n_exceeding_null"]
    n_r2 = rung2["n_surviving_dinuc"]
    n_r3_gate = rung3["n_gate_passing"]
    n_r3_null = rung3["n_exceeding_null"]
    mean_ps = rung3["mean_ps"]
    h3 = rung3["mean_h3_precision"]

    if mean_ps > 0.05 and h3 > 0.5 and n_r3_null > 5:
        return "A", "Partner-specific: encodes which position pairs with which"
    if n_r2 > 5:
        return "B", "Structure-aware beyond composition: survives dinucleotide controls"
    if n_r1 > 5:
        return "C", "Composition-sensitive: stem/loop signal absorbed by nucleotide null"
    return "D", "No detectable structure signal"


def evaluate(adapter, device="cpu", n_permutations=1000, n_derangements=1000,
             output_path=None, rungs=(1, 2, 3)):
    """Run the full RNA structure audit on a model.

    Args:
        adapter: ModelAdapter instance (must implement load/tokenize/get_all_layer_embeddings)
        device: torch device string
        n_permutations: number of permutations for Rung 1/2 null
        n_derangements: number of derangements for Rung 3 null
        output_path: optional path to save JSON results
        rungs: which rungs to run (default all three)

    Returns:
        dict with per-rung results and overall grade
    """
    families = load_families()
    adapter.load()

    report = {"model": adapter.name, "n_families": len(families)}

    rung1_results = None
    if 1 in rungs:
        rung1_results = run_rung1(adapter, families, device=device,
                                  n_permutations=n_permutations)
        report["rung1"] = {
            "mean_ratio": rung1_results["mean_ratio"],
            "families_exceeding_null": rung1_results["n_exceeding_null"],
            "families_scored": rung1_results["n_families_scored"],
        }

    rung2_results = None
    if 2 in rungs and rung1_results is not None:
        rung2_results = run_rung2(adapter, families, rung1_results, device=device,
                                  n_permutations=n_permutations)
        report["rung2"] = {
            "rung1_passing": rung2_results["n_rung1_passing"],
            "surviving_dinuc": rung2_results["n_surviving_dinuc"],
        }

    rung3_results = None
    if 3 in rungs:
        rung3_results = run_rung3(adapter, families, device=device,
                                  n_derangements=n_derangements)
        report["rung3"] = {
            "mean_ps": rung3_results["mean_ps"],
            "mean_h3_precision": rung3_results["mean_h3_precision"],
            "gate_passing": rung3_results["n_gate_passing"],
            "exceeding_null": rung3_results["n_exceeding_null"],
        }

    if rung1_results and rung2_results and rung3_results:
        grade, description = _grade(rung1_results, rung2_results, rung3_results)
        report["grade"] = grade
        report["grade_description"] = description

    full_results = {
        "report": report,
        "rung1": rung1_results,
        "rung2": rung2_results,
        "rung3": rung3_results,
    }

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(full_results, f, indent=2, default=str)

    return full_results
