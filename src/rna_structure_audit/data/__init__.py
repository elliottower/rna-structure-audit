"""Bundled Rfam family data for the benchmark."""

import json
from pathlib import Path

DATA_DIR = Path(__file__).parent


def load_families(include_withdrawn: bool = False) -> list[dict]:
    """Load the bundled Rfam families the benchmark scores.

    Each family is a dict with keys: name, sequence, dot_bracket, rfam_id.

    Five of the 52 curated families were withdrawn on annotation review and
    carry an ``excluded`` key giving the reason. They ship so that a result
    computed before the review can still be reproduced, and they are skipped
    unless asked for: scoring a model on an annotation known to be wrong
    produces a number nobody should compare against the published panel.
    """
    families = []
    for p in sorted(DATA_DIR.glob("*.json")):
        with open(p) as f:
            family = json.load(f)
        if "excluded" in family and not include_withdrawn:
            continue
        families.append(family)
    return families


def withdrawn_families() -> dict[str, dict]:
    """The withdrawn families and why, so a caller can report what it skipped.

    Each value is the family's ``excluded`` record: ``reason``, the ``rule``
    applied, and where the decision is ``logged``.
    """
    out = {}
    for p in sorted(DATA_DIR.glob("*.json")):
        with open(p) as f:
            family = json.load(f)
        if "excluded" in family:
            out[family["name"]] = family["excluded"]
    return out
