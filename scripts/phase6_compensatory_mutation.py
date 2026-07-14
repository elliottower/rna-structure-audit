"""Phase 6: Compensatory double mutation experiment.

For each base pair (i, j) in each RNA family, compare:
  - Destructive: swap only position i → breaks the pair
  - Compensatory: swap both i and j → preserves the pair

A structure-aware model should show d_dest > d_comp (CR > 1).

Usage (local test, 1 family, 1 model):
    uv run python scripts/phase6_compensatory_mutation.py --models rnafm --families tRNA_Phe_yeast --device cpu

Usage (Modal, all models, all families):
    modal run scripts/modal_phase6.py
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import torch
from scipy import stats
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

COMPLEMENT = {"A": "U", "U": "A", "C": "G", "G": "C"}
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
ROOT = Path(__file__).resolve().parent.parent
DATA_OUT = ROOT / "data" / "gpu_results" / "phase6_compensatory"


def parse_dot_bracket(db_string):
    """Extract base pairs from dot-bracket notation.

    Returns list of (i, j) tuples where i < j and both positions
    form a Watson-Crick pair in the sequence.
    """
    stack = []
    pairs = []
    for idx, char in enumerate(db_string):
        if char == "(":
            stack.append(idx)
        elif char == ")":
            if stack:
                partner = stack.pop()
                pairs.append((partner, idx))
    return sorted(pairs)


def filter_wc_pairs(sequence, pairs):
    """Keep only canonical Watson-Crick pairs (AU, UA, CG, GC)."""
    wc = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}
    return [(i, j) for i, j in pairs if (sequence[i], sequence[j]) in wc]


def complement_swap(sequence, positions):
    """Swap nucleotides at given positions to their complement."""
    seq_list = list(sequence)
    for pos in positions:
        seq_list[pos] = COMPLEMENT[seq_list[pos]]
    return "".join(seq_list)


def cosine_distance(a, b):
    """Cosine distance between two vectors."""
    sim = torch.nn.functional.cosine_similarity(a.unsqueeze(0), b.unsqueeze(0))
    return (1.0 - sim).item()


def compute_compensatory_ratio(adapter, sequence, dot_bracket, device="cpu"):
    """Compute compensatory ratio across all layers for one RNA family.

    Returns dict with per-layer CR, best layer, per-pair details.
    """
    pairs = parse_dot_bracket(dot_bracket)
    pairs = filter_wc_pairs(sequence, pairs)

    if len(pairs) < 15:
        return None

    tokens_wt = adapter.tokenize(sequence).to(device)
    layers_wt = adapter.get_all_layer_embeddings(tokens_wt)
    n_layers = len(layers_wt)

    per_layer_cr = []

    for layer_idx in range(n_layers):
        emb_wt = layers_wt[layer_idx]

        d_partner_list = []
        d_nonpartner_list = []
        pair_details = []

        all_positions = set(range(len(sequence)))
        paired_positions = set()
        for i, j in pairs:
            paired_positions.add(i)
            paired_positions.add(j)
        unpaired_positions = sorted(all_positions - paired_positions)

        for i, j in pairs:
            seq_mut = complement_swap(sequence, [i])
            tokens_mut = adapter.tokenize(seq_mut).to(device)
            emb_mut = adapter.get_all_layer_embeddings(tokens_mut)[layer_idx]

            # Perturbation at partner j when i is swapped (coupling signal)
            d_partner = cosine_distance(emb_wt[j], emb_mut[j])
            d_partner_list.append(d_partner)

            # Perturbation at non-partner positions when i is swapped (background)
            if unpaired_positions:
                bg_dists = [
                    cosine_distance(emb_wt[k], emb_mut[k])
                    for k in unpaired_positions[:10]
                ]
                d_nonpartner_list.extend(bg_dists)

            pair_type = f"{sequence[i]}-{sequence[j]}"
            pair_details.append({
                "i": i, "j": j, "type": pair_type,
                "d_partner": d_partner,
                "d_direct_i": cosine_distance(emb_wt[i], emb_mut[i]),
            })

        mean_partner = float(np.mean(d_partner_list))
        mean_nonpartner = float(np.mean(d_nonpartner_list)) if d_nonpartner_list else 0.0
        cr = mean_partner / mean_nonpartner if mean_nonpartner > 1e-10 else float("nan")
        per_layer_cr.append(cr)

    valid_crs = [(i, cr) for i, cr in enumerate(per_layer_cr) if not np.isnan(cr)]
    if not valid_crs:
        return None
    best_layer = max(valid_crs, key=lambda x: x[1])[0]
    best_cr = float(per_layer_cr[best_layer])

    gc_pairs = [(i, j) for i, j in pairs if sequence[i] in "GC"]
    au_pairs = [(i, j) for i, j in pairs if sequence[i] in "AU"]

    return {
        "best_cr": best_cr,
        "best_layer": best_layer,
        "per_layer_cr": [float(x) for x in per_layer_cr],
        "n_pairs": len(pairs),
        "n_gc_pairs": len(gc_pairs),
        "n_au_pairs": len(au_pairs),
    }


def compensatory_null(adapter, sequence, dot_bracket, n_permutations=100, device="cpu"):
    """Base-pair-type-stratified null for compensatory ratio.

    Shuffles which positions pair with which, stratified by base-pair
    type (GC vs AU). A GC pair is only reassigned to another GC pair;
    an AU pair only to another AU pair. This preserves the GC/AU pair
    distribution that is the dominant composition confound (Phase 1-2).

    Max-over-layers is applied to each permutation, matching the real
    metric's selection procedure (Phase 1→2 bug fix).
    """
    pairs = parse_dot_bracket(dot_bracket)
    pairs = filter_wc_pairs(sequence, pairs)

    if len(pairs) < 15:
        return None

    gc_nucs = {"G", "C"}
    gc_pairs = [(i, j) for i, j in pairs if sequence[i] in gc_nucs]
    au_pairs = [(i, j) for i, j in pairs if sequence[i] not in gc_nucs]

    tokens_wt = adapter.tokenize(sequence).to(device)
    layers_wt = adapter.get_all_layer_embeddings(tokens_wt)
    n_layers = len(layers_wt)

    null_crs = []

    for perm in tqdm(range(n_permutations), desc="Null permutations", leave=False):
        shuffled_pairs = []
        for type_pairs in [gc_pairs, au_pairs]:
            if not type_pairs:
                continue
            lefts = [i for i, j in type_pairs]
            rights = [j for i, j in type_pairs]
            np.random.shuffle(lefts)
            np.random.shuffle(rights)
            shuffled_pairs.extend(zip(lefts, rights))

        best_cr_this_perm = float("-inf")

        for layer_idx in range(n_layers):
            emb_wt = layers_wt[layer_idx]
            d_dest_list = []
            d_comp_list = []

            for i, j in shuffled_pairs:
                seq_dest = complement_swap(sequence, [i])
                tokens_dest = adapter.tokenize(seq_dest).to(device)
                emb_dest = adapter.get_all_layer_embeddings(tokens_dest)[layer_idx]

                seq_comp = complement_swap(sequence, [i, j])
                tokens_comp = adapter.tokenize(seq_comp).to(device)
                emb_comp = adapter.get_all_layer_embeddings(tokens_comp)[layer_idx]

                d_dest_i = cosine_distance(emb_wt[i], emb_dest[i])
                d_comp_i = cosine_distance(emb_wt[i], emb_comp[i])

                d_dest_list.append(d_dest_i)
                d_comp_list.append(d_comp_i)

            mean_dest = float(np.mean(d_dest_list))
            mean_comp = float(np.mean(d_comp_list))
            cr = mean_dest / mean_comp if mean_comp > 1e-10 else float("nan")
            if not np.isnan(cr):
                best_cr_this_perm = max(best_cr_this_perm, cr)

        null_crs.append(best_cr_this_perm)

    return {
        "null_95th": float(np.percentile(null_crs, 95)),
        "null_mean": float(np.mean(null_crs)),
        "null_std": float(np.std(null_crs)),
        "n_permutations": n_permutations,
    }


def load_rfam_families(family_names=None):
    """Load RNA families from data directory.

    Returns list of dicts with keys: name, sequence, dot_bracket.
    Looks for families in data/rfam_families/ or falls back to
    hardcoded Phase 1 families.
    """
    rfam_dir = ROOT / "data" / "rfam_families"
    if rfam_dir.exists():
        families = []
        for f in sorted(rfam_dir.glob("*.json")):
            with open(f) as fh:
                fam = json.load(fh)
            if family_names and fam["name"] not in family_names:
                continue
            families.append(fam)
        return families

    raise FileNotFoundError(
        f"No rfam_families directory found at {rfam_dir}. "
        "Run the family preparation script first."
    )


def run_phase6(adapter, families, device="cpu", compute_null=True):
    """Run Phase 6 for one model on all families."""
    results = {"per_rna": {}}

    for fam in tqdm(families, desc=f"Phase 6 [{adapter.name}]"):
        name = fam["name"]
        seq = fam["sequence"]
        db = fam["dot_bracket"]

        cr_result = compute_compensatory_ratio(adapter, seq, db, device=device)
        if cr_result is None:
            results["per_rna"][name] = {"skipped": True, "reason": "< 15 pairs"}
            continue

        if compute_null:
            null_result = compensatory_null(adapter, seq, db, device=device)
            cr_result["comp_null_95th"] = null_result["null_95th"]
            cr_result["exceeds_comp_null"] = cr_result["best_cr"] > null_result["null_95th"]

        results["per_rna"][name] = cr_result

    active = {k: v for k, v in results["per_rna"].items() if not v.get("skipped")}
    if active:
        crs = [v["best_cr"] for v in active.values()]
        results["mean_best_cr"] = float(np.mean(crs))
        results["median_best_cr"] = float(np.median(crs))
        results["families_exceeding_comp_null"] = sum(
            1 for v in active.values() if v.get("exceeds_comp_null", False)
        )
    results["families_total"] = len(active)

    return results


def main():
    parser = argparse.ArgumentParser(description="Phase 6: Compensatory mutation test")
    parser.add_argument("--models", nargs="+", default=["rnafm"],
                        help="Models to test")
    parser.add_argument("--families", nargs="+", default=None,
                        help="Specific families to test (default: all)")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--no-null", action="store_true",
                        help="Skip null computation (for quick testing)")
    args = parser.parse_args()

    DATA_OUT.mkdir(parents=True, exist_ok=True)
    families = load_rfam_families(args.families)
    print(f"Loaded {len(families)} families")

    for model_name in args.models:
        print(f"\n{'='*60}")
        print(f"Running Phase 6 for {model_name}")
        print(f"{'='*60}")

        adapter = load_adapter(model_name)
        adapter.load()

        trained_results = run_phase6(
            adapter, families, device=args.device,
            compute_null=not args.no_null,
        )

        output = {
            "model": model_name,
            "timestamp": TIMESTAMP,
            "phase": 6,
            "experiment": "compensatory_double_mutation",
            "compensatory_trained": trained_results,
        }

        out_path = DATA_OUT / f"{model_name}_phase6_{TIMESTAMP}.json"
        with open(out_path, "w") as f:
            json.dump(output, f, indent=2)
        print(f"Saved to {out_path}")

        print(f"\n{model_name} summary:")
        print(f"  Mean CR: {trained_results.get('mean_best_cr', 'N/A'):.4f}")
        print(f"  Families exceeding null: {trained_results.get('families_exceeding_comp_null', 'N/A')}/{trained_results.get('families_total', 0)}")


def load_adapter(model_name):
    """Load a model adapter by name. Uses multi_model_audit adapters."""
    adapters = {}

    try:
        from multi_model_audit import RNAFMAdapter
        adapters["rnafm"] = RNAFMAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import NTAdapter
        adapters["nt"] = NTAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import HyenaDNAAdapter
        adapters["hyenadna"] = HyenaDNAAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import EvoAdapter
        adapters["evo"] = EvoAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import CaduceusAdapter
        adapters["caduceus"] = CaduceusAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import RiNALMoAdapter
        adapters["rinalmo"] = RiNALMoAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import UTRLMAdapter
        adapters["utrlm"] = UTRLMAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import ERNIERNAAdapter
        adapters["ernierna"] = ERNIERNAAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import SpliceBERTAdapter
        adapters["splicebert"] = SpliceBERTAdapter
    except ImportError:
        pass

    try:
        from multi_model_audit import DNABERT2Adapter
        adapters["dnabert2"] = DNABERT2Adapter
    except ImportError:
        pass

    if model_name not in adapters:
        raise ValueError(f"Unknown model: {model_name}. Available: {list(adapters.keys())}")

    return adapters[model_name]()


if __name__ == "__main__":
    main()
