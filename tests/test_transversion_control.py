"""The transversion alphabet must break the pair it substitutes into.

The control exists to answer whether the Watson-Crick substitution is a real
structural disruption, so the alphabet it replaces COMPLEMENT with has one job:
whatever it puts at a stem position must not pair with the partner that position
already had. An alphabet that fails this silently turns the control into a
second copy of the experiment it is controlling.
"""

import phases_1_to_5 as phases
from modal_repaired_panel import TRANSVERSION_COMPLEMENT, result_dir

WATSON_CRICK = phases.COMPLEMENT


def test_no_substitution_leaves_a_watson_crick_pair():
    for nucleotide, partner in WATSON_CRICK.items():
        substituted = TRANSVERSION_COMPLEMENT[nucleotide]
        assert WATSON_CRICK[substituted] != partner, (
            f"{nucleotide}->{substituted} still pairs with {partner}")


def test_every_substitution_changes_the_nucleotide():
    # run_mutation_sensitivity skips a position whose substitute equals the
    # nucleotide it replaces, so a fixed point would drop that position from
    # both the stem and the loop mean without saying so.
    for nucleotide, substituted in TRANSVERSION_COMPLEMENT.items():
        assert substituted != nucleotide


def test_the_alphabet_covers_what_the_panel_contains():
    assert set(TRANSVERSION_COMPLEMENT) == set(WATSON_CRICK)


def test_the_control_differs_from_the_experiment_at_every_nucleotide():
    for nucleotide in WATSON_CRICK:
        assert TRANSVERSION_COMPLEMENT[nucleotide] != WATSON_CRICK[nucleotide]


def test_the_control_does_not_write_where_the_run_it_controls_wrote():
    # _run_model empties a directory whose stored stamp differs from the one it
    # is about to write, and the control's stamp differs on the commit. Sharing
    # the directory would delete the Watson-Crick result the control exists to
    # be compared against -- an A100 run, on models no cheap stack loads.
    for model in ("rnafm", "evo", "nt_untrained"):
        assert result_dir(model, transversion=True) != result_dir(model, transversion=False)


def test_an_ordinary_run_still_writes_under_its_own_name():
    # generate_results_tables.load() reads results/repaired_panel/<key>/ by name
    # rather than by globbing, so a renamed directory is not a missing file, it
    # is a model that silently drops out of every table.
    for model in ("rnafm", "evo", "nt_untrained"):
        assert result_dir(model, transversion=False) == model
