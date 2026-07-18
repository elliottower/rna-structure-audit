"""Rung 3: Partner specificity with within-stem derangement null.

Does the model know WHICH position pairs with which? When position i is
mutated, is the perturbation at its base-pairing partner j greater than
at j's stem-adjacent neighbors?

PS(i,j) = Delta_j - max(Delta_{j-1}, Delta_{j+1})
"""

import numpy as np
import torch
from tqdm import tqdm

from rna_structure_audit.rungs.rung1 import COMPLEMENT

WC_PAIRS = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}


def _parse_dot_bracket(db):
    stack = []
    pairs = []
    for idx, char in enumerate(db):
        if char == "(":
            stack.append(idx)
        elif char == ")":
            if stack:
                pairs.append((stack.pop(), idx))
    return sorted(pairs)


def _parse_stems(dot_bracket, sequence):
    opens = dot_bracket.count("(")
    closes = dot_bracket.count(")")
    if opens != closes:
        return []
    all_pairs = _parse_dot_bracket(dot_bracket)
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


def _get_eligible_pairs(sequence, stems):
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
                "pair_type": f"{sequence[i]}-{sequence[j]}",
            })
    return eligible


def _cosine_dist(a, b):
    sim = torch.nn.functional.cosine_similarity(a.unsqueeze(0), b.unsqueeze(0))
    return (1.0 - sim).item()


def _generate_derangement(n):
    if n < 2:
        return list(range(n))
    while True:
        perm = list(range(n))
        np.random.shuffle(perm)
        if all(perm[i] != i for i in range(n)):
            return perm


def run_rung3(adapter, families, device="cpu", n_derangements=1000,
              min_pairs_per_stem=3):
    """Run Rung 3: partner specificity with derangement null.

    Returns per-family PS scores, null exceedances, and H3 precision.
    """
    results = {"rung": 3, "model": adapter.name, "per_family": {}}
    all_ps = []
    all_h3_fracs = []

    for family in tqdm(families, desc=f"Rung 3 [{adapter.name}]"):
        seq = family["sequence"]
        db = family["dot_bracket"]
        n_pos = min(len(seq), len(db))

        stems = _parse_stems(db[:n_pos], seq[:n_pos])
        eligible = _get_eligible_pairs(seq[:n_pos], stems)

        if len(eligible) < min_pairs_per_stem:
            results["per_family"][family["name"]] = {
                "skipped": f"insufficient eligible pairs ({len(eligible)})"
            }
            continue

        tokens_wt = adapter.tokenize(seq[:n_pos]).to(device)
        layers_wt = adapter.get_all_layer_embeddings(tokens_wt)
        n_layers = len(layers_wt)

        unique_muts = {}
        for idx, p in enumerate(eligible):
            pos_i = p["i"]
            if pos_i not in unique_muts:
                unique_muts[pos_i] = []
            unique_muts[pos_i].append(idx)

        delta_profiles = {idx: {} for idx in range(len(eligible))}

        for pos_i, pair_indices in unique_muts.items():
            mut_seq = list(seq[:n_pos])
            mut_seq[pos_i] = COMPLEMENT.get(mut_seq[pos_i], mut_seq[pos_i])
            tokens_mut = adapter.tokenize("".join(mut_seq)).to(device)
            layers_mut = adapter.get_all_layer_embeddings(tokens_mut)

            for layer_idx in range(n_layers):
                emb_wt = layers_wt[layer_idx]
                emb_mut = layers_mut[layer_idx]

                positions_needed = set()
                for pidx in pair_indices:
                    p = eligible[pidx]
                    positions_needed.update([p["j"], p["j_prev"], p["j_next"]])

                deltas = {}
                for k in positions_needed:
                    if k < emb_wt.shape[0] and k < emb_mut.shape[0]:
                        deltas[k] = _cosine_dist(emb_wt[k], emb_mut[k])

                for pidx in pair_indices:
                    delta_profiles[pidx][layer_idx] = deltas

        per_layer_ps = []
        for layer_idx in range(n_layers):
            ps_vals = []
            for idx, p in enumerate(eligible):
                d = delta_profiles[idx].get(layer_idx, {})
                if p["j"] not in d or p["j_prev"] not in d or p["j_next"] not in d:
                    continue
                ps_vals.append(d[p["j"]] - max(d[p["j_prev"]], d[p["j_next"]]))
            per_layer_ps.append(float(np.mean(ps_vals)) if ps_vals else float("nan"))

        valid = [(i, v) for i, v in enumerate(per_layer_ps) if not np.isnan(v)]
        if not valid:
            results["per_family"][family["name"]] = {"skipped": "no valid PS"}
            continue

        best_layer = max(valid, key=lambda x: x[1])[0]
        best_ps = per_layer_ps[best_layer]

        partner_is_max_count = 0
        total_pairs = 0
        for idx, p in enumerate(eligible):
            d = delta_profiles[idx].get(best_layer, {})
            if p["j"] not in d or p["j_prev"] not in d or p["j_next"] not in d:
                continue
            total_pairs += 1
            if d[p["j"]] > max(d[p["j_prev"]], d[p["j_next"]]):
                partner_is_max_count += 1
        h3_frac = partner_is_max_count / total_pairs if total_pairs > 0 else 0.0

        gate_pass = best_ps > 0

        exceeds_null = False
        if gate_pass:
            stems_map = {}
            for idx, p in enumerate(eligible):
                sid = p["stem_idx"]
                if sid not in stems_map:
                    stems_map[sid] = []
                stems_map[sid].append(idx)
            derangeable = {s: idxs for s, idxs in stems_map.items() if len(idxs) >= 3}

            if derangeable:
                null_ps_values = []
                for _ in range(n_derangements):
                    perm_ps = []
                    for sid, idxs in derangeable.items():
                        deranged = _generate_derangement(len(idxs))
                        for new_k, orig_k in zip(deranged, range(len(idxs))):
                            orig_idx = idxs[orig_k]
                            new_idx = idxs[new_k]
                            p_orig = eligible[orig_idx]
                            p_new = eligible[new_idx]
                            d = delta_profiles[orig_idx].get(best_layer, {})
                            if p_new["j"] not in d or p_new["j_prev"] not in d or p_new["j_next"] not in d:
                                continue
                            perm_ps.append(
                                d[p_new["j"]] - max(d[p_new["j_prev"]], d[p_new["j_next"]])
                            )
                    null_ps_values.append(float(np.mean(perm_ps)) if perm_ps else 0.0)

                null_95 = float(np.percentile(null_ps_values, 95))
                exceeds_null = best_ps > null_95

        results["per_family"][family["name"]] = {
            "best_ps": best_ps,
            "best_layer": int(best_layer),
            "gate_pass": gate_pass,
            "exceeds_null": exceeds_null,
            "h3_precision": h3_frac,
            "n_eligible_pairs": len(eligible),
        }

        if gate_pass:
            all_ps.append(best_ps)
            all_h3_fracs.append(h3_frac)

    results["mean_ps"] = float(np.mean(all_ps)) if all_ps else 0.0
    results["mean_h3_precision"] = float(np.mean(all_h3_fracs)) if all_h3_fracs else 0.0
    n_gate = sum(1 for r in results["per_family"].values() if r.get("gate_pass"))
    n_exceed = sum(1 for r in results["per_family"].values() if r.get("exceeds_null"))
    results["n_gate_passing"] = n_gate
    results["n_exceeding_null"] = n_exceed
    return results
