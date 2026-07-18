"""Rung 2: Dinucleotide-stratified null.

Of families that exceed the nucleotide null (Rung 1), how many survive
when the permutation null is stratified by dinucleotide context?
This absorbs composition artifacts that operate at the dinucleotide level.
"""

import numpy as np
import torch
from scipy.spatial.distance import cosine
from tqdm import tqdm

from rna_structure_audit.rungs.rung1 import COMPLEMENT, _parse_positions, _validate_family


def run_rung2(adapter, families, rung1_results, device="cpu", n_permutations=1000):
    """Run Rung 2: dinucleotide-stratified null on families that passed Rung 1.

    Takes rung1_results to identify which families exceeded the nucleotide null.
    Only evaluates families that passed Rung 1.
    """
    results = {"rung": 2, "model": adapter.name, "n_permutations": n_permutations,
               "per_family": {}}
    n_survive = 0

    passing_families = {
        name for name, r in rung1_results["per_family"].items()
        if r.get("exceeds_nuc_null")
    }

    for fam_idx, family in enumerate(tqdm(families, desc=f"Rung 2 [{adapter.name}]")):
        if family["name"] not in passing_families:
            continue

        valid, reason = _validate_family(family)
        if not valid:
            continue

        seq = family["sequence"]
        db = family["dot_bracket"]
        n_pos = min(len(seq), len(db))

        tokens = adapter.tokenize(seq[:n_pos]).to(device)
        wt_embs = adapter.get_all_layer_embeddings(tokens)
        n_layers = len(wt_embs)

        paired = _parse_positions(db, max_pos=n_pos)
        labels = np.array(["stem" if i in paired else "loop" for i in range(n_pos)])

        per_layer_dists = {l: np.zeros(n_pos) for l in range(n_layers)}
        per_layer_valid = {l: np.zeros(n_pos, dtype=bool) for l in range(n_layers)}

        for pos in range(n_pos):
            comp = COMPLEMENT.get(seq[pos])
            if comp is None or comp == seq[pos]:
                continue
            mut_seq = list(seq[:n_pos])
            mut_seq[pos] = comp
            mut_tokens = adapter.tokenize("".join(mut_seq)).to(device)
            mut_embs = adapter.get_all_layer_embeddings(mut_tokens)

            for layer_idx in range(n_layers):
                wt_layer = wt_embs[layer_idx].cpu().numpy()
                mut_layer = mut_embs[layer_idx].cpu().numpy()
                if pos >= min(wt_layer.shape[0], mut_layer.shape[0]):
                    continue
                per_layer_dists[layer_idx][pos] = float(cosine(wt_layer[pos], mut_layer[pos]))
                per_layer_valid[layer_idx][pos] = True

        real_ratios = {}
        for layer_idx in range(n_layers):
            dists = per_layer_dists[layer_idx]
            vmask = per_layer_valid[layer_idx]
            stem_mask = np.array([labels[i] == "stem" and vmask[i] for i in range(n_pos)])
            loop_mask = np.array([labels[i] == "loop" and vmask[i] for i in range(n_pos)])
            if stem_mask.sum() < 3 or loop_mask.sum() < 3:
                continue
            mean_loop = float(np.mean(dists[loop_mask]))
            ratio = float(np.mean(dists[stem_mask])) / mean_loop if mean_loop > 1e-10 else 0.0
            real_ratios[layer_idx] = ratio

        if not real_ratios:
            continue

        best_layer = max(real_ratios, key=real_ratios.get)
        best_ratio = real_ratios[best_layer]
        valid_layers = sorted(real_ratios.keys())

        rng = np.random.default_rng(42 + fam_idx * 1000)
        dinucs = np.array([
            ("N" if i == 0 else seq[i - 1]) + seq[i] for i in range(n_pos)
        ])

        dinuc_null_ratios = []
        for _ in range(n_permutations):
            shuffled = labels.copy()
            for dinuc in np.unique(dinucs):
                idx = np.where(dinucs == dinuc)[0]
                if len(idx) >= 2:
                    shuffled[idx] = rng.permutation(shuffled[idx])
            max_ratio = 0.0
            for li in valid_layers:
                dists = per_layer_dists[li]
                vmask = per_layer_valid[li]
                vi = np.where(vmask)[0]
                vd, vs = dists[vi], shuffled[vi]
                s, l = vd[vs == "stem"], vd[vs == "loop"]
                if len(s) == 0 or len(l) == 0:
                    continue
                ml = np.mean(l)
                max_ratio = max(max_ratio, np.mean(s) / ml if ml > 1e-10 else 0.0)
            dinuc_null_ratios.append(max_ratio)

        dinuc_null_95 = float(np.percentile(dinuc_null_ratios, 95))
        survives = best_ratio > dinuc_null_95

        results["per_family"][family["name"]] = {
            "best_ratio": best_ratio,
            "best_layer": int(best_layer),
            "dinuc_null_95th": dinuc_null_95,
            "survives_dinuc_null": survives,
        }
        if survives:
            n_survive += 1

    results["n_rung1_passing"] = len(passing_families)
    results["n_surviving_dinuc"] = n_survive
    return results
