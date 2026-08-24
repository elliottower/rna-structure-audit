"""Every anchor in the v13 patch still matches the document it edits.

The patch is a list of verbatim strings lifted from paper_v12.tex. An anchor
that stopped matching -- a rewrapped line, a changed number, an edit applied
twice -- is invisible until the run lands and the patch is finally allowed to
execute, which is the worst moment to discover it. `build` is separated from
`main` so these run without the generated tables or the results on disk.
"""

import importlib.util
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
V12 = (REPO / "paper" / "paper_v12.tex").read_text()

_spec = importlib.util.spec_from_file_location(
    "patch_v13", REPO / "paper" / "patches" / "patch_v13_repaired_panel.py")
patch = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(patch)


@pytest.fixture(scope="module")
def v13():
    return patch.build(V12)


def test_every_edit_anchor_matches_v12_exactly_once():
    # Reported one at a time. build() aborts on the first, which hides the rest
    # behind whichever happens to come earliest in the list.
    wrong = {label: V12.count(old) for label, old, _ in patch.EDITS
             if V12.count(old) != 1}
    assert not wrong


def test_the_five_results_tables_are_gone_from_the_body(v13):
    for label, generated in patch.GENERATED_TABLES:
        assert "\\label{" + label + "}" not in v13
        assert "\\input{generated/" + generated + "}" in v13


def test_the_tables_that_are_not_regenerated_survive(v13):
    # tab:models, tab:htt and tab:quarantine carry no repaired-panel number and
    # have no generator, so a swap that took them too would delete content
    # nothing replaces.
    for label in ("tab:models", "tab:htt", "tab:quarantine"):
        assert "\\label{" + label + "}" in v13


def test_a_swapped_table_takes_its_whole_environment(v13):
    assert v13.count("\\begin{table}") == v13.count("\\end{table}")
    assert v13.count("\\begin{tabular}") == v13.count("\\end{tabular}")


def test_the_panel_macros_are_defined_before_the_abstract_uses_them(v13):
    definition = v13.index("\\input{generated/panel_description}")
    assert definition < v13.index("\\begin{abstract}")
    assert definition < v13.index("\\panelN{}")


def test_the_results_macros_are_defined_before_the_body_uses_them(v13):
    assert v13.index("\\input{generated/results_macros}") < v13.index("\\provenanceRuns{}")


def test_no_generated_file_is_inputted_twice(v13):
    for name in patch.REQUIRED_GENERATED:
        assert v13.count("\\input{generated/" + name[: -len(".tex")] + "}") == 1


def test_the_hand_carried_class_census_is_gone(v13):
    for stale in ("tRNAs (7 families)", "rRNAs (4)", "riboswitches (8)",
                  "CRISPR repeats (3)", "snRNAs (5)", "miRNA precursors (5)"):
        assert stale not in v13


def test_the_paper_does_not_narrate_its_own_corrections(v13):
    for changelog in ("An earlier analysis", "reverses this interpretation",
                      "in a prior run", "Corrected from", "originally passed"):
        assert changelog not in v13


def test_both_gates_reach_the_reader(v13):
    assert "registered as 7 at $N = 32$" in v13
    assert "the gate is 8 at the repaired $N$" in v13


def test_a_rewrapped_anchor_aborts_the_build():
    label, old, _ = patch.EDITS[0]
    with pytest.raises(SystemExit, match="0 matches"):
        patch.build(V12.replace(old, old.replace("\n", " "), 1))


def test_replace_table_refuses_a_label_it_cannot_find_once():
    with pytest.raises(SystemExit, match="0 labels"):
        patch.replace_table(V12, "tab:does_not_exist", "whatever")


def test_the_patch_introduces_macros_and_names_them_consistently():
    used = patch.macros_used()
    assert used, "the patch replaces hand-carried numbers with macros; it uses none"
    assert used == {name for name in used if name.startswith(("panel", "provenance"))}, \
        "an unexpected macro family: the two checks below only cover these two"


def test_every_panel_macro_is_one_the_panel_generator_defines():
    # panel_description.tex is on disk already -- it is computed from
    # data/rfam_families/ and needs no GPU run.
    defined = set(re.findall(r"\\newcommand\{\\([a-zA-Z]+)\}",
                             (REPO / "paper/generated/panel_description.tex").read_text()))
    assert {n for n in patch.macros_used() if n.startswith("panel")} <= defined


def test_every_provenance_macro_is_one_the_results_generator_writes():
    # Not on disk until the runs land, so the check is against the literal
    # macro("name", ...) calls that write it.
    source = (REPO / "scripts/generate_results_tables.py").read_text()
    written = set(re.findall(r'macro\("([a-zA-Z]+)"', source))
    assert {n for n in patch.macros_used() if n.startswith("provenance")} <= written
