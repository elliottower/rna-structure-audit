"""Phase 6: Perturbation specificity test.

For each interior base pair (i, j) in a stem, swap position i and measure
whether perturbation at partner j exceeds perturbation at j's adjacent
stem neighbors. PS(i,j) = Delta_j - max(Delta_{j-1}, Delta_{j+1}).

Metric defined in PREREGISTRATION_PHASE6_V2.md.

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

from family_checkpoint import FamilyCheckpoint, no_checkpoint
from family_seed import family_rng
from token_spans import content_spans, nucleotide_rows

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

COMPLEMENT = {"A": "U", "U": "A", "C": "G", "G": "C"}
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
ROOT = Path(__file__).resolve().parent.parent
DATA_OUT = ROOT / "data" / "gpu_results" / "phase6_compensatory"
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
WC_PAIRS = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}

ADAPTER_OFFSETS = {
    "rnafm": 0,
    "rinalmo": 0,
    "utrlm": 0,
    "ernierna": 0,
    "splicebert": 0,
    "hyenadna": 0,
    "caduceus": 0,
    "evo": 0,
    "nt": 0,
    "dnabert2": 0,
}
NON_CHARACTER_TOKENIZERS = {"nt", "dnabert2"}
SUBWORD_RESOLUTIONS = {"6mer", "bpe"}


def _row_map(adapter, sequence):
    """Hidden-state row holding each nucleotide, or None when rows are nucleotides.

    Reading `emb[k]` for nucleotide `k` is right only for the eight models
    that emit one token per nucleotide. NT v2 puts six nucleotides in a token,
    so a 106-nucleotide family has about eighteen rows and every partner
    position in a stem falls past the end of the array; the bounds check below
    then drops the pair, and a family that loses all of its pairs is recorded
    as having no valid PS values. That discarded 34 of NT v2's 38 qualifying
    families and 30 of DNABERT-2's; see docs/OPEN_DEFECTS.md, D13.
    """
    if getattr(adapter, "token_resolution", "nucleotide") not in SUBWORD_RESOLUTIONS:
        return None
    spans = content_spans(adapter.tokenizer, sequence.replace("U", "T"))
    return nucleotide_rows(spans, len(sequence))


def parse_dot_bracket(db_string):
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


def parse_stems(dot_bracket, sequence):
    """Group consecutive stacked WC pairs into stems.

    A stem is a maximal run of pairs (i, j), (i+1, j-1), (i+2, j-2)...
    Only canonical WC pairs are included.
    """
    opens = dot_bracket.count("(")
    closes = dot_bracket.count(")")
    if opens != closes:
        return []
    all_pairs = parse_dot_bracket(dot_bracket)
    wc_pairs = [(i, j) for i, j in all_pairs if (sequence[i], sequence[j]) in WC_PAIRS]
    pair_set = set(wc_pairs)

    used = set()
    stems = []
    for i, j in wc_pairs:
        if (i, j) in used:
            continue
        stem = [(i, j)]
        used.add((i, j))
        ci, cj = i + 1, j - 1
        while ci < cj and (ci, cj) in pair_set and (ci, cj) not in used:
            stem.append((ci, cj))
            used.add((ci, cj))
            ci += 1
            cj -= 1
        stems.append(stem)
    return stems


def get_eligible_pairs(sequence, stems):
    """Return interior pairs from stems with >= 3 WC pairs.

    Each eligible pair includes references to its stem-adjacent j positions.
    Returns list of dicts with keys: i, j, j_prev, j_next, stem_idx, pair_type.
    """
    eligible = []
    for stem_idx, stem in enumerate(stems):
        if len(stem) < 3:
            continue
        for k in range(1, len(stem) - 1):
            i, j = stem[k]
            _, j_prev = stem[k - 1]
            _, j_next = stem[k + 1]
            eligible.append({
                "i": i, "j": j,
                "j_prev": j_prev, "j_next": j_next,
                "stem_idx": stem_idx,
                "pos_in_stem": k,
                "stem_length": len(stem),
                "pair_type": f"{sequence[i]}-{sequence[j]}",
            })
    return eligible


def complement_swap(sequence, position):
    seq_list = list(sequence)
    seq_list[position] = COMPLEMENT[seq_list[position]]
    return "".join(seq_list)


def cosine_distance(a, b):
    sim = torch.nn.functional.cosine_similarity(a.unsqueeze(0), b.unsqueeze(0))
    return (1.0 - sim).item()


def compute_delta_profiles(adapter, sequence, eligible_pairs, all_stems, device="cpu", offset=0):
    """Compute perturbation profiles for all eligible pairs at all layers.

    Returns:
        delta_profiles: dict mapping pair_index -> {layer -> {pos -> Delta}}
        stem_positions: set of all positions in any stem
        loop_positions: set of all unpaired positions
    """
    tokens_wt = adapter.tokenize(sequence).to(device)
    layers_wt = adapter.get_all_layer_embeddings(tokens_wt)
    n_layers = len(layers_wt)
    rows_wt = _row_map(adapter, sequence)

    all_positions = set(range(len(sequence)))
    stem_positions = set()
    for stem in all_stems:
        for i, j in stem:
            stem_positions.add(i)
            stem_positions.add(j)
    loop_positions = all_positions - stem_positions

    # Cache mutations: many pairs may share the same mutated position i
    unique_mutations = {}
    for idx, p in enumerate(eligible_pairs):
        pos_i = p["i"]
        if pos_i not in unique_mutations:
            unique_mutations[pos_i] = []
        unique_mutations[pos_i].append(idx)

    delta_profiles = {idx: {} for idx in range(len(eligible_pairs))}

    for pos_i, pair_indices in tqdm(unique_mutations.items(), desc="Computing deltas", leave=False):
        seq_mut = complement_swap(sequence, pos_i)
        tokens_mut = adapter.tokenize(seq_mut).to(device)
        layers_mut = adapter.get_all_layer_embeddings(tokens_mut)
        rows_mut = None if rows_wt is None else _row_map(adapter, seq_mut)

        for layer_idx in range(n_layers):
            emb_wt = layers_wt[layer_idx]
            emb_mut = layers_mut[layer_idx]

            # Compute Delta at all relevant positions for this mutation
            positions_needed = set()
            for pidx in pair_indices:
                p = eligible_pairs[pidx]
                positions_needed.update([p["j"], p["j_prev"], p["j_next"]])
            positions_needed.update(stem_positions)
            positions_needed.update(loop_positions)
            positions_needed.discard(pos_i)

            deltas_this_layer = {}
            for k in positions_needed:
                if rows_wt is None:
                    k_wt = k_mut = k + offset
                else:
                    k_wt, k_mut = int(rows_wt[k]), int(rows_mut[k])
                if k_wt < emb_wt.shape[0] and k_mut < emb_mut.shape[0]:
                    deltas_this_layer[k] = cosine_distance(emb_wt[k_wt], emb_mut[k_mut])

            for pidx in pair_indices:
                delta_profiles[pidx][layer_idx] = deltas_this_layer

    return delta_profiles, stem_positions, loop_positions, n_layers


def compute_ps_from_deltas(eligible_pairs, delta_profiles, n_layers):
    """Compute PS at every layer from precomputed Delta profiles.

    Returns per-layer mean PS and per-pair details at best layer.
    """
    per_layer_ps = []

    for layer_idx in range(n_layers):
        ps_values = []
        for idx, p in enumerate(eligible_pairs):
            deltas = delta_profiles[idx].get(layer_idx, {})
            if p["j"] not in deltas or p["j_prev"] not in deltas or p["j_next"] not in deltas:
                continue
            d_partner = deltas[p["j"]]
            d_prev = deltas[p["j_prev"]]
            d_next = deltas[p["j_next"]]
            d_adj = max(d_prev, d_next)
            ps = d_partner - d_adj
            ps_values.append(ps)
        mean_ps = float(np.mean(ps_values)) if ps_values else float("nan")
        per_layer_ps.append(mean_ps)

    valid = [(i, v) for i, v in enumerate(per_layer_ps) if not np.isnan(v)]
    if not valid:
        return None
    best_layer = max(valid, key=lambda x: x[1])[0]
    best_ps = per_layer_ps[best_layer]

    pair_details_at_best = []
    for idx, p in enumerate(eligible_pairs):
        deltas = delta_profiles[idx].get(best_layer, {})
        if p["j"] not in deltas or p["j_prev"] not in deltas or p["j_next"] not in deltas:
            continue
        d_partner = deltas[p["j"]]
        d_prev = deltas[p["j_prev"]]
        d_next = deltas[p["j_next"]]
        d_adj = max(d_prev, d_next)
        pair_details_at_best.append({
            "i": p["i"], "j": p["j"],
            "pair_type": p["pair_type"],
            "stem_idx": p["stem_idx"],
            "d_partner": d_partner,
            "d_adj": d_adj,
            "ps": d_partner - d_adj,
            "partner_is_max": d_partner > d_adj,
        })

    return {
        "per_layer_ps": [float(x) for x in per_layer_ps],
        "best_ps": float(best_ps),
        "best_layer": best_layer,
        "pair_details": pair_details_at_best,
    }


def generate_derangement(n, rng):
    """Generate a random derangement of range(n). No element maps to itself."""
    if n < 2:
        return list(range(n))
    while True:
        perm = list(range(n))
        rng.shuffle(perm)
        if all(perm[i] != i for i in range(n)):
            return perm


def derangement_null(eligible_pairs, delta_profiles, n_layers, best_layer, rng,
                     n_derangements=1000):
    """Within-stem derangement null.

    Shuffles partner assignments within each stem, recomputes PS from
    existing Delta profiles. Returns null thresholds for primary (at
    best_layer) and conservative (independently max'd) variants.
    """
    stems_map = {}
    for idx, p in enumerate(eligible_pairs):
        sid = p["stem_idx"]
        if sid not in stems_map:
            stems_map[sid] = []
        stems_map[sid].append(idx)

    derangeable = {sid: idxs for sid, idxs in stems_map.items() if len(idxs) >= 3}

    if not derangeable:
        return {
            "null_95th_primary": None,
            "null_95th_conservative": None,
            "null_mean_primary": None,
            "n_derangements": 0,
            "n_stems_in_null": 0,
            "null_available": False,
        }

    null_ps_primary = []
    null_ps_conservative = []

    for _ in range(n_derangements):
        deranged_ps_at_best = []
        deranged_ps_per_layer = [[] for _ in range(n_layers)]

        for sid, pair_idxs in derangeable.items():
            k = len(pair_idxs)
            derangement = generate_derangement(k, rng)
            j_positions = [eligible_pairs[pair_idxs[orig]]["j"] for orig in range(k)]
            j_prev_positions = [eligible_pairs[pair_idxs[orig]]["j_prev"] for orig in range(k)]
            j_next_positions = [eligible_pairs[pair_idxs[orig]]["j_next"] for orig in range(k)]

            for orig_idx_in_stem, deranged_partner_idx in enumerate(derangement):
                pair_idx = pair_idxs[orig_idx_in_stem]
                fake_j = j_positions[deranged_partner_idx]
                fake_j_prev = j_prev_positions[deranged_partner_idx]
                fake_j_next = j_next_positions[deranged_partner_idx]

                for layer_idx in range(n_layers):
                    deltas = delta_profiles[pair_idx].get(layer_idx, {})
                    if fake_j not in deltas or fake_j_prev not in deltas or fake_j_next not in deltas:
                        continue
                    d_partner = deltas[fake_j]
                    d_prev = deltas[fake_j_prev]
                    d_next = deltas[fake_j_next]
                    d_adj = max(d_prev, d_next)
                    ps = d_partner - d_adj
                    deranged_ps_per_layer[layer_idx].append(ps)
                    if layer_idx == best_layer:
                        deranged_ps_at_best.append(ps)

        if deranged_ps_at_best:
            null_ps_primary.append(float(np.mean(deranged_ps_at_best)))

        layer_means = [float(np.mean(lps)) if lps else float("-inf") for lps in deranged_ps_per_layer]
        null_ps_conservative.append(max(layer_means))

    return {
        "null_95th_primary": float(np.percentile(null_ps_primary, 95)) if null_ps_primary else 0.0,
        "null_95th_conservative": float(np.percentile(null_ps_conservative, 95)) if null_ps_conservative else 0.0,
        "null_mean_primary": float(np.mean(null_ps_primary)) if null_ps_primary else 0.0,
        "n_derangements": n_derangements,
        "n_stems_in_null": len(derangeable),
        "null_available": True,
    }


def positive_control(eligible_pairs, delta_profiles, best_layer, stem_positions, loop_positions):
    """Stem > loop screening gate (per prereg: outside Bonferroni family)."""
    paired_stem = []
    paired_loop = []

    for idx, p in enumerate(eligible_pairs):
        deltas = delta_profiles[idx].get(best_layer, {})
        pos_i = p["i"]
        stem_d = [v for k, v in deltas.items() if k in stem_positions and k != pos_i]
        loop_d = [v for k, v in deltas.items() if k in loop_positions]
        if stem_d and loop_d:
            paired_stem.append(float(np.mean(stem_d)))
            paired_loop.append(float(np.mean(loop_d)))

    if len(paired_stem) < 2:
        return {"pass": False, "p_value": 1.0, "reason": "insufficient paired data"}

    t_stat, p_val = stats.ttest_rel(paired_stem, paired_loop)
    passes = p_val < 0.05 and t_stat > 0
    return {
        "pass": bool(passes),
        "p_value": float(p_val),
        "mean_stem_delta": float(np.mean(paired_stem)),
        "mean_loop_delta": float(np.mean(paired_loop)),
    }


def h3_precision_test(pair_details, eligible_pairs):
    """Fraction of deep-interior pairs where partner is max among {j, j-1, j+1}.

    H3 uses pairs >= 2 from each stem end (stricter than the primary
    eligible filter which only drops terminals). Uses pos_in_stem stored
    at creation time in get_eligible_pairs.
    """
    eligible_lookup = {(p["i"], p["j"]): p for p in eligible_pairs}

    h3_details = []
    for p in pair_details:
        ep = eligible_lookup.get((p["i"], p["j"]))
        if ep is None:
            continue
        pos = ep["pos_in_stem"]
        slen = ep["stem_length"]
        if pos >= 2 and pos <= slen - 3:
            h3_details.append(p)

    partner_is_max_count = sum(1 for p in h3_details if p["partner_is_max"])
    n = len(h3_details)
    if n == 0:
        return {"fraction": 0.0, "n": 0, "p_value": 1.0}
    fraction = partner_is_max_count / n
    p_val = stats.binomtest(partner_is_max_count, n, 1.0 / 3.0, alternative="greater").pvalue
    return {
        "fraction": float(fraction),
        "n": n,
        "partner_max_count": partner_is_max_count,
        "p_value": float(p_val),
    }


def load_rfam_families(family_names=None):
    rfam_dir = ROOT / "data" / "rfam_families"
    if rfam_dir.exists():
        families = []
        for f in sorted(rfam_dir.glob("*.json")):
            with open(f) as fh:
                fam = json.load(fh)
            if family_names and fam["name"] not in family_names:
                continue
            # Records whose annotation could not be repaired against the Rfam
            # seed alignment carry an `excluded` block; see DEVIATIONS.md.
            if "excluded" in fam:
                continue
            families.append(fam)
        return families
    raise FileNotFoundError(f"No rfam_families directory at {rfam_dir}.")


def load_adapter(model_name):
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


def run_phase6(adapter, families, device="cpu", compute_null=True, offset=0,
               checkpoint=None):
    results = {"per_rna": {}}
    ckpt = checkpoint if checkpoint is not None else no_checkpoint()

    def record(name, entry):
        results["per_rna"][name] = entry
        ckpt.record(name, entry)

    for fam in tqdm(families, desc=f"Phase 6 [{adapter.name}]"):
        name = fam["name"]
        if name in ckpt.done:
            results["per_rna"][name] = ckpt.done[name]
            continue
        seq = fam["sequence"]
        db = fam["dot_bracket"]
        quarantined = name in QUARANTINED

        if len(seq) != len(db):
            record(name, {
                "skipped": True,
                "reason": f"sequence/dot_bracket length mismatch ({len(seq)} vs {len(db)})",
            })
            continue

        stems = parse_stems(db, seq)
        n_wc_total = sum(len(s) for s in stems)

        if n_wc_total < 15:
            record(name, {
                "skipped": True,
                "reason": f"< 15 WC pairs ({n_wc_total} found)",
                "n_wc_pairs_total": n_wc_total,
            })
            continue

        eligible = get_eligible_pairs(seq, stems)

        if len(eligible) < 5:
            record(name, {
                "skipped": True,
                "reason": f"< 5 eligible interior pairs ({len(eligible)} found)",
                "n_wc_pairs_total": n_wc_total,
                "n_stems": len(stems),
            })
            continue

        delta_profiles, stem_pos, loop_pos, n_layers = compute_delta_profiles(
            adapter, seq, eligible, stems, device=device, offset=offset,
        )

        ps_result = compute_ps_from_deltas(eligible, delta_profiles, n_layers)
        if ps_result is None:
            record(name, {"skipped": True, "reason": "no valid PS values"})
            continue

        best_layer = ps_result["best_layer"]

        pc = positive_control(eligible, delta_profiles, best_layer, stem_pos, loop_pos)

        null_result = None
        if compute_null:
            null_result = derangement_null(eligible, delta_profiles, n_layers,
                                           best_layer, family_rng(name))

        h3 = h3_precision_test(ps_result["pair_details"], eligible)

        gc_ps = [p["ps"] for p in ps_result["pair_details"] if p["pair_type"] in ("C-G", "G-C")]
        au_ps = [p["ps"] for p in ps_result["pair_details"] if p["pair_type"] in ("A-U", "U-A")]

        entry = {
            "best_ps": ps_result["best_ps"],
            "best_layer": ps_result["best_layer"],
            "per_layer_ps": ps_result["per_layer_ps"],
            "n_eligible_pairs": len(eligible),
            "n_stems": len([s for s in stems if len(s) >= 3]),
            "n_wc_pairs_total": sum(len(s) for s in stems),
            "positive_control": pc,
            "h3_precision": h3,
            "ps_gc_pairs": float(np.mean(gc_ps)) if gc_ps else None,
            "ps_au_pairs": float(np.mean(au_ps)) if au_ps else None,
            "quarantined": quarantined,
        }

        if null_result:
            entry["null_95th_primary"] = null_result["null_95th_primary"]
            entry["null_95th_conservative"] = null_result["null_95th_conservative"]
            entry["null_available"] = null_result.get("null_available", True)
            if null_result["null_95th_primary"] is not None:
                entry["exceeds_null_primary"] = ps_result["best_ps"] > null_result["null_95th_primary"]
                entry["exceeds_null_conservative"] = ps_result["best_ps"] > null_result["null_95th_conservative"]
            else:
                entry["exceeds_null_primary"] = None
                entry["exceeds_null_conservative"] = None

        record(name, entry)

    active = {k: v for k, v in results["per_rna"].items()
              if not v.get("skipped") and not v.get("quarantined")
              and v.get("positive_control", {}).get("pass", False)}
    if active:
        ps_values = [v["best_ps"] for v in active.values()]
        results["mean_best_ps"] = float(np.mean(ps_values))
        results["median_best_ps"] = float(np.median(ps_values))
        if compute_null:
            results["families_exceeding_null_primary"] = sum(
                1 for v in active.values() if v.get("exceeds_null_primary") is True)
            results["families_exceeding_null_conservative"] = sum(
                1 for v in active.values() if v.get("exceeds_null_conservative") is True)
            results["families_no_null"] = [
                k for k, v in active.items() if v.get("null_available") is False]
    results["families_total"] = len(active)
    results["families_quarantined"] = [k for k, v in results["per_rna"].items() if v.get("quarantined")]
    results["families_skipped"] = [k for k, v in results["per_rna"].items() if v.get("skipped")]
    results["families_failed_gate"] = [
        k for k, v in results["per_rna"].items()
        if not v.get("skipped") and not v.get("quarantined")
        and not v.get("positive_control", {}).get("pass", False)
    ]

    return results


def main():
    parser = argparse.ArgumentParser(description="Phase 6: Perturbation specificity test")
    parser.add_argument("--models", nargs="+", default=["rnafm"])
    parser.add_argument("--families", nargs="+", default=None)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--no-null", action="store_true",
                        help="Skip null computation (for quick testing)")
    parser.add_argument("--offset", type=int, default=None,
                        help="Override token offset (default: per-adapter lookup)")
    parser.add_argument("--allow-non-character", action="store_true",
                        help="Run non-character tokenizers (NT, DNABERT-2) with caveated results")
    args = parser.parse_args()

    DATA_OUT.mkdir(parents=True, exist_ok=True)
    families = load_rfam_families(args.families)
    print(f"Loaded {len(families)} families")

    for model_name in args.models:
        if model_name in NON_CHARACTER_TOKENIZERS and not args.allow_non_character:
            print(f"\nSKIPPING {model_name}: non-character tokenizer. "
                  f"Per-nucleotide position mapping is undefined for k-mer/BPE models. "
                  f"Use --allow-non-character to run with caveated results.")
            continue

        offset = args.offset if args.offset is not None else ADAPTER_OFFSETS.get(model_name, 0)
        tokenizer_caveated = model_name in NON_CHARACTER_TOKENIZERS

        print(f"\n{'='*60}")
        print(f"Phase 6 PS metric for {model_name} (offset={offset})")
        print(f"{'='*60}")

        adapter = load_adapter(model_name)
        adapter.load()

        trained_results = run_phase6(
            adapter, families, device=args.device,
            compute_null=not args.no_null, offset=offset,
        )

        output = {
            "model": model_name,
            "timestamp": TIMESTAMP,
            "phase": 6,
            "metric": "perturbation_specificity",
            "preregistration": "PREREGISTRATION_PHASE6_V2.md",
            "offset": offset,
            "seeding": "family_seed.family_rng, derived from the family name",
            "tokenizer_caveated": tokenizer_caveated,
            "results": trained_results,
        }

        out_path = DATA_OUT / f"{model_name}_phase6_ps_{TIMESTAMP}.json"
        with open(out_path, "w") as f:
            json.dump(output, f, indent=2)
        print(f"Saved to {out_path}")

        print(f"\n{model_name} summary:")
        print(f"  Mean PS: {trained_results.get('mean_best_ps', 'N/A')}")
        print(f"  Families (confirmatory): {trained_results.get('families_total', 0)}")
        print(f"  Families (quarantined): {trained_results.get('families_quarantined', [])}")
        print(f"  Families (skipped): {trained_results.get('families_skipped', [])}")
        print(f"  Families (failed gate): {trained_results.get('families_failed_gate', [])}")
        if not args.no_null:
            print(f"  Exceeding null (primary): {trained_results.get('families_exceeding_null_primary', 'N/A')}")
            print(f"  Exceeding null (conservative): {trained_results.get('families_exceeding_null_conservative', 'N/A')}")


if __name__ == "__main__":
    main()
