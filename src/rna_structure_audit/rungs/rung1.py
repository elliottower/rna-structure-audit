"""Rung 1: Mutation sensitivity with nucleotide-stratified null.

Do stems respond differently than loops to complement mutations?
A nucleotide-stratified permutation null controls for composition bias.
"""

import numpy as np
import torch
from scipy.spatial.distance import cosine
from tqdm import tqdm

COMPLEMENT = {"A": "U", "U": "A", "C": "G", "G": "C", "T": "A"}


def _parse_positions(dot_bracket, max_pos=None):
    paired = set()
    stack = []
    n = len(dot_bracket) if max_pos is None else min(len(dot_bracket), max_pos)
    for i in range(n):
        if dot_bracket[i] == "(":
            stack.append(i)
            paired.add(i)
        elif dot_bracket[i] == ")":
            if stack:
                stack.pop()
            paired.add(i)
    return paired


def _validate_family(family):
    seq = family["sequence"]
    db = family["dot_bracket"]
    n = min(len(seq), len(db))
    if n < 20:
        return False, "too short"
    opens = db.count("(")
    closes = db.count(")")
    if opens != closes:
        return False, f"unbalanced dot-bracket ({opens} vs {closes})"
    paired = _parse_positions(db, max_pos=n)
    unpaired_count = n - len(paired)
    if len(paired) < 5 or unpaired_count < 5:
        return False, "insufficient positions"
    return True, "ok"


def run_rung1(adapter, families, device="cpu", n_permutations=1000):
    """Run Rung 1: mutation sensitivity ratio with nucleotide-stratified null.

    Returns dict with per-family ratios, null exceedances, and summary stats.
    """
    results = {"rung": 1, "model": adapter.name, "n_permutations": n_permutations,
               "per_family": {}}
    all_ratios = []

    for fam_idx, family in enumerate(tqdm(families, desc=f"Rung 1 [{adapter.name}]")):
        valid, reason = _validate_family(family)
        if not valid:
            results["per_family"][family["name"]] = {"skipped": reason}
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
            results["per_family"][family["name"]] = {"skipped": "no valid layers"}
            continue

        best_layer = max(real_ratios, key=real_ratios.get)
        best_ratio = real_ratios[best_layer]

        rng = np.random.default_rng(42 + fam_idx * 1000)
        nucs = np.array(list(seq[:n_pos]))
        valid_layers = sorted(real_ratios.keys())

        null_max_ratios = []
        for _ in range(n_permutations):
            shuffled = labels.copy()
            for nuc in ("A", "U", "G", "C"):
                idx = np.where(nucs == nuc)[0]
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
            null_max_ratios.append(max_ratio)

        null_95 = float(np.percentile(null_max_ratios, 95))
        exceeds = best_ratio > null_95

        results["per_family"][family["name"]] = {
            "best_ratio": best_ratio,
            "best_layer": int(best_layer),
            "nuc_null_95th": null_95,
            "exceeds_nuc_null": exceeds,
        }
        all_ratios.append(best_ratio)

    results["mean_ratio"] = float(np.mean(all_ratios)) if all_ratios else 0.0
    results["n_exceeding_null"] = sum(
        1 for r in results["per_family"].values() if r.get("exceeds_nuc_null")
    )
    results["n_families_scored"] = len(all_ratios)
    return results
