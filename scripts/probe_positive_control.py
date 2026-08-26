"""Can the Rung 3 probe detect partner specificity in a system known to have it?

Eight of ten models score at their own chance rate. That is consistent with two
readings a reader cannot separate: the models do not resolve base-pair partners,
or the probe cannot detect it. Nothing in the panel distinguishes them, because
every arm is a model under test.

This runs the unmodified `run_phase6` on two constructed adapters whose answer is
known before the pipeline sees them:

  oracle    each position's embedding is a function of its own nucleotide and its
            partner's, so mutating i moves j and moves nothing adjacent to j.
            Partner specificity is present by construction.
  local     each position's embedding is a function of its own nucleotide and its
            two sequence neighbours, and knows nothing about pairing. Partner
            specificity is absent by construction.

A probe that works reports near 1.0 for the oracle and its own chance rate for
the local model. Anything else bounds what the panel's numbers can mean.

    uv run --no-project --with "numpy<2" --with scipy --with tqdm --with torch \\
        --python 3.12 python scripts/probe_positive_control.py
"""

import json
from pathlib import Path

import numpy as np
import torch

from phase6_compensatory_mutation import load_rfam_families, parse_stems, run_phase6

REPO = Path(__file__).resolve().parents[1]
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}
D_MODEL = 32
VOCAB = {"A": 0, "C": 1, "G": 2, "U": 3, "T": 3, "N": 4}


def _signature(index, salt):
    """A fixed pseudo-random vector for an integer, with no global RNG."""
    k = np.arange(D_MODEL, dtype=np.float64)
    return np.sin((index + 1) * (k + 1) * 0.7301 + salt)


class _Base:
    n_layers = 1
    d_model = D_MODEL
    token_resolution = "nucleotide"

    def __init__(self):
        self.structure = None      # set per family by the caller

    def tokenize(self, sequence):
        return torch.tensor([[VOCAB.get(c, 4) for c in sequence]])


class Oracle(_Base):
    """Position k encodes (nucleotide at k, nucleotide at k's partner)."""

    name = "oracle"

    def get_all_layer_embeddings(self, tokens):
        ids = tokens[0].tolist()
        partner = self.structure
        rows = []
        for k, token in enumerate(ids):
            mate = partner.get(k)
            mate_token = ids[mate] if mate is not None and mate < len(ids) else 4
            rows.append(_signature(token, 0.0) + _signature(mate_token, 1.7))
        return [torch.tensor(np.stack(rows))]


class Local(_Base):
    """Position k encodes its own nucleotide and its two sequence neighbours."""

    name = "local"

    def get_all_layer_embeddings(self, tokens):
        ids = tokens[0].tolist()
        rows = []
        for k, token in enumerate(ids):
            left = ids[k - 1] if k > 0 else 4
            right = ids[k + 1] if k + 1 < len(ids) else 4
            rows.append(_signature(token, 0.0)
                        + 0.5 * _signature(left, 2.3)
                        + 0.5 * _signature(right, 3.1))
        return [torch.tensor(np.stack(rows))]


def pairing(family):
    """Position -> partner, from the stems the pipeline itself parses."""
    mapping = {}
    for stem in parse_stems(family["dot_bracket"], family["sequence"]):
        for i, j in stem:
            mapping[i] = j
            mapping[j] = i
    return mapping


def score(adapter_class, families):
    adapter = adapter_class()
    precisions, chances = [], []
    for family in families:
        adapter.structure = pairing(family)
        result = run_phase6(adapter, [family], compute_null=True)
        entry = result["per_rna"][family["name"]]
        if entry.get("skipped") or entry.get("h3_chance_fraction") is None:
            continue
        precisions.append(entry["h3_precision"]["fraction"])
        chances.append(entry["h3_chance_fraction"])
    return np.array(precisions), np.array(chances)


def main():
    families = [f for f in load_rfam_families() if f["name"] not in QUARANTINED]
    print(f"  {len(families)} non-quarantined families\n")
    print(f"  {'adapter':10s} {'n':>4s} {'precision':>10s} {'chance':>8s} {'excess':>8s}")
    out = {}
    for adapter_class in (Oracle, Local):
        precisions, chances = score(adapter_class, families)
        out[adapter_class.name] = {
            "n": len(precisions),
            "precision": float(precisions.mean()),
            "chance": float(chances.mean()),
        }
        print(f"  {adapter_class.name:10s} {len(precisions):4d} "
              f"{precisions.mean():10.3f} {chances.mean():8.3f} "
              f"{precisions.mean() - chances.mean():+8.3f}")

    path = REPO / "results" / "probe_positive_control.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"\n  wrote {path.relative_to(REPO)}")
    print("  A working probe reads near 1.0 for the oracle and near its own\n"
          "  chance rate for the local model.")


if __name__ == "__main__":
    main()
