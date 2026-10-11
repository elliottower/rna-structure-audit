"""The moved structure functions must equal the superseded script's on every family."""

import json
from pathlib import Path

import phase6_compensatory_mutation as superseded
import stem_pairs

FAMILIES = Path(__file__).resolve().parents[1] / "data" / "rfam_families"


def test_stems_and_eligible_pairs_match_the_superseded_script():
    paths = sorted(FAMILIES.glob("*.json"))
    assert paths
    for path in paths:
        family = json.loads(path.read_text())
        sequence, structure = family["sequence"], family["dot_bracket"]
        if len(sequence) != len(structure):
            continue
        stems = stem_pairs.parse_stems(structure, sequence)
        assert stems == superseded.parse_stems(structure, sequence), family["name"]
        assert (stem_pairs.get_eligible_pairs(sequence, stems)
                == superseded.get_eligible_pairs(sequence, stems)), family["name"]
