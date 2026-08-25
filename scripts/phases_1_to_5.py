"""Phase 1-5 analysis: mutation sensitivity, attention-contact, structure probing.

Ported from rna-target-pipeline/batch1/code/modal_e1_druggability.py.
Adapted to use Rfam family JSON format (dot_bracket instead of structure).

These are pure analysis functions — no Modal, no model loading.
"""

import math

import numpy as np

from family_checkpoint import FamilyCheckpoint, no_checkpoint
from family_seed import family_rng
from scipy import stats
from scipy.spatial.distance import cosine
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_score
from token_spans import content_bounds, content_spans, nucleotide_rows
from tqdm import tqdm

COMPLEMENT = {"A": "U", "U": "A", "C": "G", "G": "C", "T": "A"}

RNA_TOKENIZER_MODELS = {"rnafm", "rinalmo", "utrlm", "ernierna", "splicebert"}
ATTENTION_MODELS = {"rnafm", "nt", "rinalmo", "utrlm", "ernierna", "splicebert", "dnabert2"}
SUBWORD_MODELS = {"nt", "dnabert2"}


def _get_paired_and_unpaired_positions(structure, max_pos=None):
    paired = []
    unpaired = []
    stack = []
    n = len(structure) if max_pos is None else min(len(structure), max_pos)
    for i in range(n):
        c = structure[i]
        if c == "(":
            stack.append(i)
            paired.append(i)
        elif c == ")":
            if stack:
                stack.pop()
            paired.append(i)
        else:
            unpaired.append(i)
    return paired, unpaired


def _parse_structure_to_contacts(structure):
    n = len(structure)
    contacts = np.zeros((n, n), dtype=np.float32)
    stack = []
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if stack:
                j = stack.pop()
                contacts[i, j] = 1.0
                contacts[j, i] = 1.0
    return contacts


def _tokenizer_sequence(model_key, seq):
    """The exact string the adapter hands its tokenizer.

    Spans have to be computed on that string and not on the stored sequence:
    the RNA tokenizers see U and the DNA tokenizers see T, and DNABERT-2's
    merges depend on which.
    """
    if model_key in RNA_TOKENIZER_MODELS:
        return seq.replace("T", "U")
    return seq.replace("U", "T")


def _family_spans(adapter, model_key, seq):
    """Token spans for the models that put more than one nucleotide in a token."""
    if model_key not in SUBWORD_MODELS:
        return None
    return content_spans(adapter.tokenizer, _tokenizer_sequence(model_key, seq))


def _aggregate_contacts_to_tokens(contact_map, spans):
    """Two tokens are in contact when any nucleotide they cover is paired.

    The spans come from the tokenizer. Taking the window as `ti * 6` instead
    put 88.7% of DNABERT-2's rows on a stretch of sequence sharing no
    nucleotide with the token whose attention that row carried; see
    docs/OPEN_DEFECTS.md, D12.
    """
    n_tokens = len(spans)
    token_contacts = np.zeros((n_tokens, n_tokens), dtype=np.float32)
    for ti, (si, ei) in enumerate(spans):
        for tj, (sj, ej) in enumerate(spans):
            if contact_map[si:ei, sj:ej].any():
                token_contacts[ti, tj] = 1.0
    return token_contacts


def _validate_rna(rna):
    seq = rna["sequence"]
    db = rna["dot_bracket"]
    n = min(len(seq), len(db))
    if n < 20:
        return False, "too short"
    opens = db.count("(")
    closes = db.count(")")
    if opens != closes:
        return False, f"unbalanced dot-bracket ({opens} vs {closes})"
    paired, unpaired = _get_paired_and_unpaired_positions(db, max_pos=n)
    if len(paired) < 5 or len(unpaired) < 5:
        return False, f"insufficient positions (stem={len(paired)}, loop={len(unpaired)})"
    return True, "ok"


def _get_dnabert2_all_hidden_states(model, tokens, bounds=None):
    import torch

    hidden_states = []
    hooks = []

    for layer in model.encoder.layer:
        def _hook(module, input, output, _list=hidden_states):
            out = output[0] if isinstance(output, tuple) else output
            _list.append(out.detach())
        hooks.append(layer.register_forward_hook(_hook))

    with torch.no_grad():
        model(tokens)

    for h in hooks:
        h.remove()

    start, stop = bounds if bounds is not None else (1, -1)
    result = []
    for hs in hidden_states:
        if hs.dim() == 2:
            hs = hs.unsqueeze(0)
        result.append(hs[0, start:stop, :])
    return result


def _extract_dnabert2_attention(model, tokens):
    import torch

    qkv_outputs = []
    hooks = []

    for layer in model.encoder.layer:
        self_attn = layer.attention.self
        if hasattr(self_attn, "Wqkv"):
            def _hook(module, input, output, _list=qkv_outputs):
                _list.append(output.detach())
            hooks.append(self_attn.Wqkv.register_forward_hook(_hook))
        else:
            for h in hooks:
                h.remove()
            return None

    with torch.no_grad():
        model(tokens)

    for h in hooks:
        h.remove()

    if not qkv_outputs:
        return None

    attentions = []
    n_heads = 12
    for qkv in qkv_outputs:
        if qkv.dim() == 2:
            qkv = qkv.unsqueeze(0)
        batch, seq_len, three_d = qkv.shape
        d_model = three_d // 3
        d_head = d_model // n_heads
        q, k, _ = qkv.chunk(3, dim=-1)
        q = q.view(batch, seq_len, n_heads, d_head).transpose(1, 2)
        k = k.view(batch, seq_len, n_heads, d_head).transpose(1, 2)
        scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(d_head)
        attn_weights = torch.softmax(scores, dim=-1)
        attentions.append(attn_weights)

    return attentions


def run_mutation_sensitivity(adapter, model_key, families, device="cuda",
                             n_permutations=1000, checkpoint=None, positions=None):
    """Mutation sensitivity ratio with composition-controlled null.

    Returns per-family best_ratio, nucleotide-stratified null, dinucleotide null.

    `positions` is a second checkpoint. Given one, every per-layer per-position
    cosine distance is written to it with the labels, the nucleotides and, for
    the subword models, whether the mutated position kept its token span. The
    null is a function of those distances and those labels alone, so a null that
    turns out to be wrong -- the registered one has a measured 9.1% false
    positive rate on randomly initialized models, docs/OPEN_DEFECTS.md D10 --
    can be replaced from the stored values without loading a model again.
    """
    import torch

    results = {"model": model_key, "metric": "mutation_sensitivity",
               "n_permutations": n_permutations, "per_rna": {}}
    ckpt = checkpoint if checkpoint is not None else no_checkpoint()

    for rna_idx, rna in enumerate(tqdm(families, desc=f"Phase 1-5 [{model_key}]")):
        if rna["name"] in ckpt.done:
            results["per_rna"][rna["name"]] = ckpt.done[rna["name"]]
            continue
        valid, reason = _validate_rna(rna)
        if not valid:
            results["per_rna"][rna["name"]] = {"skipped": reason}
            ckpt.record(rna["name"], {"skipped": reason})
            continue

        seq = rna["sequence"]
        db = rna["dot_bracket"]
        n_pos = min(len(seq), len(db))

        tokens = adapter.tokenize(_tokenizer_sequence(model_key, seq[:n_pos])).to(device)
        wt_embs = adapter.get_all_layer_embeddings(tokens)
        n_layers = len(wt_embs)
        wt_spans = _family_spans(adapter, model_key, seq[:n_pos])
        wt_rows = (nucleotide_rows(wt_spans, n_pos) if wt_spans is not None
                   else np.arange(n_pos))
        same_span = np.ones(n_pos, dtype=bool)

        paired, unpaired = _get_paired_and_unpaired_positions(db, max_pos=n_pos)
        labels = np.array(["stem" if i in set(paired) else "loop" for i in range(n_pos)])

        per_layer_dists = {l: np.zeros(n_pos) for l in range(n_layers)}
        per_layer_valid = {l: np.zeros(n_pos, dtype=bool) for l in range(n_layers)}

        for pos in range(n_pos):
            nuc = seq[pos]
            comp = COMPLEMENT.get(nuc)
            if comp is None or comp == nuc:
                continue

            mut_seq = list(seq[:n_pos])
            mut_seq[pos] = comp
            mut_seq_str = "".join(mut_seq)

            mut_tokens = adapter.tokenize(
                _tokenizer_sequence(model_key, mut_seq_str)).to(device)
            mut_embs = adapter.get_all_layer_embeddings(mut_tokens)

            if wt_spans is None:
                mut_rows = wt_rows
            else:
                mut_spans = _family_spans(adapter, model_key, mut_seq_str)
                mut_rows = nucleotide_rows(mut_spans, n_pos)
                same_span[pos] = wt_spans[wt_rows[pos]] == mut_spans[mut_rows[pos]]

            wt_row, mut_row = int(wt_rows[pos]), int(mut_rows[pos])
            for layer_idx in range(n_layers):
                wt_layer = wt_embs[layer_idx].cpu().numpy()
                mut_layer = mut_embs[layer_idx].cpu().numpy()
                if wt_row >= wt_layer.shape[0] or mut_row >= mut_layer.shape[0]:
                    continue
                # float64 for the same reason as `phase6.cosine_distance`: the
                # embeddings are float32 and `1 - cos_sim` on two rows that
                # barely differ loses its significant digits to cancellation.
                dist = cosine(wt_layer[wt_row].astype(np.float64),
                              mut_layer[mut_row].astype(np.float64))
                per_layer_dists[layer_idx][pos] = float(dist)
                per_layer_valid[layer_idx][pos] = True

        real_ratios_per_layer = {}
        for layer_idx in range(n_layers):
            dists = per_layer_dists[layer_idx]
            valid_mask = per_layer_valid[layer_idx]
            stem_mask = np.array([labels[i] == "stem" and valid_mask[i] for i in range(n_pos)])
            loop_mask = np.array([labels[i] == "loop" and valid_mask[i] for i in range(n_pos)])
            if int(stem_mask.sum()) < 3 or int(loop_mask.sum()) < 3:
                continue
            mean_stem = float(np.mean(dists[stem_mask]))
            mean_loop = float(np.mean(dists[loop_mask]))
            ratio = mean_stem / mean_loop if mean_loop > 1e-10 else 0.0
            real_ratios_per_layer[layer_idx] = {
                "ratio": ratio, "mean_stem": mean_stem, "mean_loop": mean_loop,
                "n_stem": int(stem_mask.sum()), "n_loop": int(loop_mask.sum()),
            }

        if not real_ratios_per_layer:
            results["per_rna"][rna["name"]] = {"skipped": "no valid layers"}
            continue

        best_layer_idx = max(real_ratios_per_layer, key=lambda l: real_ratios_per_layer[l]["ratio"])
        best_ratio = real_ratios_per_layer[best_layer_idx]["ratio"]

        valid_layers = sorted(real_ratios_per_layer.keys())
        # Seeded from the family name, not its position in the panel: the null
        # is defined within the family, so withdrawing another record must not
        # change this family's draws. See DEVIATIONS.md, 2026-08-24.
        perm_rng = family_rng(rna["name"])
        nucs_full = np.array(list(seq[:n_pos]))

        null_max_ratios = []
        for _ in range(n_permutations):
            shuffled_labels = labels.copy()
            for nuc in ("A", "U", "G", "C"):
                nuc_idx = np.where(nucs_full == nuc)[0]
                if len(nuc_idx) < 2:
                    continue
                shuffled_labels[nuc_idx] = perm_rng.permutation(shuffled_labels[nuc_idx])

            max_perm_ratio = 0.0
            for layer_idx in valid_layers:
                dists = per_layer_dists[layer_idx]
                valid_mask = per_layer_valid[layer_idx]
                vi = np.where(valid_mask)[0]
                vd = dists[vi]
                vs = shuffled_labels[vi]
                perm_stem = vd[vs == "stem"]
                perm_loop = vd[vs == "loop"]
                if len(perm_stem) == 0 or len(perm_loop) == 0:
                    continue
                ml = np.mean(perm_loop)
                max_perm_ratio = max(max_perm_ratio, np.mean(perm_stem) / ml if ml > 1e-10 else 0.0)
            null_max_ratios.append(max_perm_ratio)

        nuc_null_95th = float(np.percentile(null_max_ratios, 95)) if null_max_ratios else 0.0
        exceeds_nuc = best_ratio > nuc_null_95th

        rna_result = {
            "best_ratio": best_ratio, "best_layer": int(best_layer_idx),
            "nuc_null_95th": nuc_null_95th, "exceeds_nuc_null": exceeds_nuc,
            "n_stem": real_ratios_per_layer[best_layer_idx]["n_stem"],
            "n_loop": real_ratios_per_layer[best_layer_idx]["n_loop"],
        }

        if exceeds_nuc:
            dinucs_full = np.array([
                ("N" if i == 0 else seq[i - 1]) + seq[i] for i in range(n_pos)
            ])
            dinuc_null_max_ratios = []
            for _ in range(n_permutations):
                shuffled_labels = labels.copy()
                for dinuc in np.unique(dinucs_full):
                    dinuc_idx = np.where(dinucs_full == dinuc)[0]
                    if len(dinuc_idx) < 2:
                        continue
                    shuffled_labels[dinuc_idx] = perm_rng.permutation(shuffled_labels[dinuc_idx])
                max_perm_ratio = 0.0
                for layer_idx in valid_layers:
                    dists = per_layer_dists[layer_idx]
                    valid_mask = per_layer_valid[layer_idx]
                    vi = np.where(valid_mask)[0]
                    vd = dists[vi]
                    vs = shuffled_labels[vi]
                    perm_stem = vd[vs == "stem"]
                    perm_loop = vd[vs == "loop"]
                    if len(perm_stem) == 0 or len(perm_loop) == 0:
                        continue
                    ml = np.mean(perm_loop)
                    max_perm_ratio = max(max_perm_ratio, np.mean(perm_stem) / ml if ml > 1e-10 else 0.0)
                dinuc_null_max_ratios.append(max_perm_ratio)

            dinuc_null_95th = float(np.percentile(dinuc_null_max_ratios, 95))
            rna_result["dinuc_null_95th"] = dinuc_null_95th
            rna_result["exceeds_dinuc_null"] = best_ratio > dinuc_null_95th

        results["per_rna"][rna["name"]] = rna_result
        ckpt.record(rna["name"], rna_result)
        if positions is not None:
            positions.record(rna["name"], {
                "nucleotides": seq[:n_pos],
                "labels": "".join("S" if lab == "stem" else "L" for lab in labels),
                "same_span": "".join("1" if s else "0" for s in same_span),
                "layer_distances": [[round(float(d), 6) for d in per_layer_dists[l]]
                                    for l in range(n_layers)],
                "layer_valid": ["".join("1" if v else "0" for v in per_layer_valid[l])
                                for l in range(n_layers)],
            })

    # Read back from per_rna rather than accumulating in the loop: a resumed
    # run never executes the loop body for a checkpointed family, so a list
    # appended in the loop would hold only the families this container did.
    all_ratios = [r["best_ratio"] for r in results["per_rna"].values()
                  if "best_ratio" in r]
    results["mean_best_ratio"] = float(np.mean(all_ratios)) if all_ratios else 0.0
    results["median_best_ratio"] = float(np.median(all_ratios)) if all_ratios else 0.0
    n_exceed = sum(1 for r in results["per_rna"].values() if r.get("exceeds_nuc_null"))
    results["n_exceeding_nuc_null"] = n_exceed
    results["n_families_scored"] = len(all_ratios)
    return results


def run_attention_contact(adapter, model_key, families, device="cuda", checkpoint=None):
    """Attention-contact Spearman correlation."""
    import torch

    results = {"model": model_key, "metric": "attention_contact", "per_rna": {}}
    ckpt = checkpoint if checkpoint is not None else no_checkpoint()

    if model_key not in ATTENTION_MODELS:
        results["skipped"] = "no attention (SSM architecture)"
        return results

    for rna in tqdm(families, desc=f"Attention [{model_key}]"):
        if rna["name"] in ckpt.done:
            results["per_rna"][rna["name"]] = ckpt.done[rna["name"]]
            continue
        valid, reason = _validate_rna(rna)
        if not valid:
            results["per_rna"][rna["name"]] = {"skipped": reason}
            ckpt.record(rna["name"], {"skipped": reason})
            continue

        seq = rna["sequence"]
        db = rna["dot_bracket"]
        n_pos = min(len(seq), len(db))

        tokens = adapter.tokenize(_tokenizer_sequence(model_key, seq[:n_pos])).to(device)
        spans = _family_spans(adapter, model_key, seq[:n_pos])
        # RNA-FM builds its own ids and exposes no tokenizer, so the bounds are
        # asked for only where a tokenizer decides them.
        bounds = (content_bounds(adapter.tokenizer, tokens[0].tolist())
                  if spans is not None else None)

        with torch.no_grad():
            if model_key == "dnabert2":
                attentions = _extract_dnabert2_attention(adapter.model, tokens)
            else:
                out = adapter.model(tokens, output_attentions=True, return_dict=True)
                attentions = getattr(out, "attentions", None) if not isinstance(out, tuple) else (out[3] if len(out) > 3 else None)

        if attentions is None:
            results["per_rna"][rna["name"]] = {"skipped": "no attention returned"}
            ckpt.record(rna["name"], {"skipped": "no attention returned"})
            continue

        contact_map = _parse_structure_to_contacts(db)
        contact_flat = contact_map[np.triu_indices(n_pos, k=1)]

        best_corr = 0.0
        best_layer = 0

        for layer_idx, attn_tensor in enumerate(attentions):
            attn = attn_tensor[0].cpu().numpy()
            n_heads = attn.shape[0]

            eff_n = n_pos
            eff_contact_flat = contact_flat
            if model_key in RNA_TOKENIZER_MODELS:
                attn = attn[:, 1:n_pos+1, 1:n_pos+1]
            elif model_key in SUBWORD_MODELS:
                start, stop = bounds
                attn = attn[:, start:stop, start:stop]
                if attn.shape[1] != len(spans):
                    continue
                eff_n = len(spans)
                tok_contact = _aggregate_contacts_to_tokens(contact_map, spans)
                eff_contact_flat = tok_contact[np.triu_indices(eff_n, k=1)]

            if attn.shape[1] != eff_n or attn.shape[2] != eff_n:
                continue

            for head in range(n_heads):
                attn_sym = (attn[head] + attn[head].T) / 2
                attn_flat = attn_sym[np.triu_indices(eff_n, k=1)]
                if np.std(attn_flat) < 1e-10:
                    continue
                rho, _ = stats.spearmanr(attn_flat, eff_contact_flat)
                if np.isnan(rho):
                    rho = 0.0
                if rho > best_corr:
                    best_corr = float(rho)
                    best_layer = layer_idx

        entry = {"best_corr": best_corr, "best_layer": best_layer}
        results["per_rna"][rna["name"]] = entry
        ckpt.record(rna["name"], entry)

    return results


def run_structure_probing(adapter, model_key, families, device="cuda"):
    """Logistic regression probing for stem/loop classification."""
    import torch

    results = {"model": model_key, "metric": "structure_probing", "per_layer": {}}

    all_embeddings = []
    all_labels = []
    all_groups = []

    for rna_idx, rna in enumerate(tqdm(families, desc=f"Probing [{model_key}]")):
        valid, reason = _validate_rna(rna)
        if not valid:
            continue

        seq = rna["sequence"]
        db = rna["dot_bracket"]
        n_pos = min(len(seq), len(db))

        tokens = adapter.tokenize(_tokenizer_sequence(model_key, seq[:n_pos])).to(device)
        spans = _family_spans(adapter, model_key, seq[:n_pos])
        rows = nucleotide_rows(spans, n_pos) if spans is not None else None

        if model_key == "dnabert2":
            layer_embs = _get_dnabert2_all_hidden_states(
                adapter.model, tokens,
                bounds=content_bounds(adapter.tokenizer, tokens[0].tolist()))
        else:
            layer_embs = adapter.get_all_layer_embeddings(tokens)
        paired, unpaired = _get_paired_and_unpaired_positions(db, max_pos=n_pos)
        lab = np.zeros(n_pos, dtype=np.int32)
        for p in paired:
            lab[p] = 1
        groups = np.full(n_pos, rna_idx, dtype=np.int32)

        for layer_idx in range(len(layer_embs)):
            emb = layer_embs[layer_idx].cpu().numpy()
            if rows is not None:
                if emb.shape[0] != len(spans):
                    continue
                emb = emb[rows]
            elif emb.shape[0] < n_pos:
                continue
            emb_trimmed = emb[:n_pos]
            if layer_idx >= len(all_embeddings):
                all_embeddings.append([])
                all_labels.append([])
                all_groups.append([])
            all_embeddings[layer_idx].append(emb_trimmed)
            all_labels[layer_idx].append(lab)
            all_groups[layer_idx].append(groups)

    for layer_idx in range(len(all_embeddings)):
        X = np.concatenate(all_embeddings[layer_idx], axis=0)
        y = np.concatenate(all_labels[layer_idx], axis=0)
        g = np.concatenate(all_groups[layer_idx], axis=0)

        if len(np.unique(y)) < 2:
            continue
        n_groups = len(np.unique(g))
        n_splits = min(5, n_groups)
        if n_splits < 2:
            continue

        clf = LogisticRegression(max_iter=1000, C=1.0)
        gkf = GroupKFold(n_splits=n_splits)
        scores = cross_val_score(clf, X, y, cv=gkf, groups=g, scoring="balanced_accuracy")
        results["per_layer"][str(layer_idx)] = {
            "mean_balanced_accuracy": float(np.mean(scores)),
            "std_balanced_accuracy": float(np.std(scores)),
            "n_samples": int(len(y)),
            "n_targets": n_groups,
        }

    if results["per_layer"]:
        best = max(results["per_layer"].items(), key=lambda x: x[1]["mean_balanced_accuracy"])
        results["best_layer"] = int(best[0])
        results["best_accuracy"] = best[1]["mean_balanced_accuracy"]

    return results
