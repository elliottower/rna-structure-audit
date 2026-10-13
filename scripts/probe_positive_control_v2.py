"""Can the Rung 3 probe detect partner specificity in a system known to have it?

Several models score at their own chance rate. That is consistent with two
readings a reader cannot separate: the models do not resolve base-pair partners,
or the probe cannot detect it. Nothing in the panel distinguishes them, because
every arm is a model under test.

Supersedes `probe_positive_control.py`, which ran the superseded Phase 6 script
and compared precision with a chance rate over every eligible pair. This runs
`run_phase6` of `phase6_compensatory_mutation_v3.py`, the script every model row
now comes from, and reads the chance rate on the pairs precision is read on. It
also records each construction's matched-decoy rate against all decoys.

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

Run inside the panel's own image by `scripts/modal_probe_control.py`.
"""

import json
from pathlib import Path

import numpy as np
import torch

from phase6_compensatory_mutation_v3 import load_rfam_families, run_phase6
from stem_pairs import parse_stems

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
    """Per-family precision, matched chance rate and matched-decoy rate."""
    adapter = adapter_class()
    per_family = {}
    for family in families:
        adapter.structure = pairing(family)
        result = run_phase6(adapter, [family], compute_null=True)
        entry = result["per_rna"][family["name"]]
        if entry.get("skipped") or entry.get("h3_chance_fraction_matched") is None:
            continue
        decoy = entry["decoys"]["family_all"]
        per_family[family["name"]] = {
            "precision": entry["h3_precision"]["fraction"],
            "chance_matched": entry["h3_chance_fraction_matched"],
            # One layer, so index 0; None where the family has no matched decoy.
            "decoy_rate_all": decoy["rate"][0],
            "decoy_pairs_all": decoy["pairs"][0],
        }
    return per_family


def run():
    """Both constructions on the non-quarantined panel, as one result structure."""
    families = [f for f in load_rfam_families() if f["name"] not in QUARANTINED]
    out = {"script": "scripts/probe_positive_control_v2.py",
           "phase6": "scripts/phase6_compensatory_mutation_v3.py",
           "n_families_loaded": len(families), "d_model": D_MODEL, "adapters": {}}
    for adapter_class in (Oracle, Local):
        per_family = score(adapter_class, families)
        precisions = np.array([v["precision"] for v in per_family.values()])
        chances = np.array([v["chance_matched"] for v in per_family.values()])
        decoys = np.array([v["decoy_rate_all"] for v in per_family.values()
                           if v["decoy_rate_all"] is not None and v["decoy_pairs_all"] >= 3])
        out["adapters"][adapter_class.name] = {
            "n": len(per_family),
            "precision": float(precisions.mean()),
            "chance_matched": float(chances.mean()),
            "decoy_families": int(len(decoys)),
            "decoy_rate_all": float(decoys.mean()) if len(decoys) else None,
            "per_family": per_family,
        }
    return out
