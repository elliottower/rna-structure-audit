"""The confirmatory set is decided by the frozen registration, not by a literal.

Ten scripts carry the quarantined families as a hardcoded set. The point of
reading them out of `PREREGISTRATION_PHASE6_V2.md` is that an edit to any one of
those copies -- or to the frozen document -- stops being invisible, so the tests
that matter are the ones that fail when a copy drifts and when the document's two
statements of the commitment stop agreeing.
"""

import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

_spec = importlib.util.spec_from_file_location(
    "registered_quarantine", REPO / "scripts" / "registered_quarantine.py")
rq = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rq)


def test_every_hardcoded_copy_matches_the_frozen_document():
    disagreeing = {path.name: sorted(value)
                   for path, value in rq.literals().items()
                   if value != rq.QUARANTINE}
    assert disagreeing == {}


def test_the_quarantine_is_a_subset_of_the_panel():
    assert rq.QUARANTINE <= rq.panel_names()


def _doctored(tmp_path, replacement):
    text = (REPO / "PREREGISTRATION_PHASE6_V2.md").read_text()
    path = tmp_path / "doctored.md"
    path.write_text(text.replace("tRNA_Ala_human", replacement, 1))
    return path


def test_the_two_statements_must_name_the_same_families(tmp_path):
    """Replacing the name in whichever sentence comes first leaves the two
    disagreeing, which is the failure a single literal cannot express."""
    with pytest.raises(ValueError, match="same"):
        rq.registered_quarantine(_doctored(tmp_path, "U1_snRNA"))


def test_a_name_the_panel_does_not_contain_is_not_read_as_a_family(tmp_path):
    with pytest.raises(ValueError, match="same"):
        rq.registered_quarantine(_doctored(tmp_path, "tRNA_Ala_martian"))


def test_a_document_quarantining_nothing_is_an_error(tmp_path):
    text = (REPO / "PREREGISTRATION_PHASE6_V2.md").read_text()
    for name in sorted(rq.QUARANTINE):
        text = text.replace(name, "some_family")
    path = tmp_path / "empty.md"
    path.write_text(text)
    with pytest.raises(ValueError, match="no quarantined family"):
        rq.registered_quarantine(path)
