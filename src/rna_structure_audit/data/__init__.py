"""Bundled Rfam family data for the benchmark."""

import json
from pathlib import Path

DATA_DIR = Path(__file__).parent


def load_families() -> list[dict]:
    """Load all bundled Rfam families.

    Each family is a dict with keys: name, sequence, dot_bracket, rfam_id.
    """
    families = []
    for p in sorted(DATA_DIR.glob("*.json")):
        with open(p) as f:
            families.append(json.load(f))
    return families
