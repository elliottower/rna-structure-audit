"""The random seed a family's null is drawn from.

Both rungs build their null within a single family: Rungs 1-2 permute stem/loop
labels within nucleotide strata of that family's own sequence, and Rung 3
deranges partner assignments within that family's own stems. Neither null is
defined over the panel, so neither should change when the panel does.

Both violated that. `phases_1_to_5.py` seeded from `42 + rna_idx * 1000`, the
family's position in the loaded list, and `phase6_compensatory_mutation.py`
seeded one global generator that every family then drew from in panel order.
Withdrawing a record shifted the stream every later family received, so families
whose annotation never changed silently got different nulls.

Seeding from the family name fixes both: a family's null draws are reproducible
across runs and unaffected by which other families are present.
"""

import hashlib

import numpy as np


def family_rng(name: str) -> np.random.Generator:
    """A generator determined by the family name and nothing else.

    blake2b rather than hash(), which is salted per process and would not
    reproduce across runs.
    """
    digest = hashlib.blake2b(name.encode("utf-8"), digest_size=8).digest()
    return np.random.default_rng(int.from_bytes(digest, "big"))
