"""No paper-facing number may be read from the mis-sized exceedance field.

`exceeds_null_primary` compares a maximum over L layers against a null evaluated
at the single layer that maximum selected, so its size is 1 - 0.95^L rather than
0.05 -- measured at 0.49 to 0.60 on randomly initialized weights, which carry no
partner specificity (D15). `exceeds_null_conservative` maxes over layers on the
null side too and sits at nominal.

Both are stored and the registration requires both to be reported wherever they
disagree, which they do in every model. What must not happen is a published
count silently reading the primary field again, so this names the files allowed
to touch it and fails on any other.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
FIELDS = ("exceeds_null_primary", "families_exceeding_null_primary")

# Scripts whose job is to compute, compare or explain both variants. Everything
# else in `scripts/` either reports a number or generates a table, and must read
# the conservative field.
ALLOWED = {
    # Writes both fields.
    "phase6_compensatory_mutation.py",
    # Measure the primary variant's size, or compare runs; reading it is the point.
    "audit_null_selection_bias.py",
    "audit_ps_noise_floor.py",
    "audit_run_to_run_drift.py",
    "compare_phase6_repair.py",
    # Superseded pipelines, kept for provenance and not run.
    "phase6_compensatory_mutation_v1.py",
    "phase6_compensatory_mutation_v2_preflight.py",
    # Merges shards field-by-field; it carries the name to copy the value.
    "merge_phase6_families.py",
    # Writes its own result files with its own field of the same name.
    "multi_seq_ps.py",
    # An audit of the deposited files, which carry the primary field.
    "audit_duplicate_families.py",
}

# Runners that produced earlier deposits and are superseded by
# `modal_repaired_panel.py`. They print a count to their own log rather than to
# the manuscript, and editing them would change scripts whose output is already
# on record. They are listed one by one rather than matched by prefix, so a new
# runner is caught.
SUPERSEDED_RUNNERS = {
    "modal_all_phases.py",
    "modal_phase6.py",
    "modal_phase6_caduceus.py",
    "modal_phase6_compat.py",
    "modal_phase6_dnabert2.py",
    "modal_phase6_flashattn_models.py",
    "modal_phase6_synthetic_covariation.py",
    "modal_phase6_untrained_all.py",
    "modal_phase6_untrained_caduceus.py",
    "modal_phase6_untrained_dnabert2_evo.py",
    "modal_phase6_untrained_ernierna.py",
    "modal_phase6_untrained_evo.py",
    "modal_phase6_untrained_rinalmo.py",
}
ALLOWED |= SUPERSEDED_RUNNERS


def offenders():
    """Files naming a primary field without also naming its conservative twin.

    The rule is not that the primary field is untouchable -- the registration
    requires both counts wherever they disagree, so a reporter that prints both
    is doing the right thing. The rule is that the primary count never appears
    alone, because alone it reads as a result rather than as the sensitivity
    analysis it is.
    """
    found = []
    for path in sorted(REPO.glob("scripts/*.py")):
        if path.name in ALLOWED:
            continue
        source = path.read_text()
        for field in FIELDS:
            twin = field.replace("primary", "conservative")
            if field in source and twin not in source:
                lines = [n for n, line in enumerate(source.splitlines(), start=1)
                         if field in line]
                found.append((path.name, lines[0], field))
    return found


def test_no_script_reports_the_mis_sized_count_without_the_calibrated_one():
    found = offenders()
    assert not found, (
        "these name a field whose test size is 0.49 to 0.60 without naming the "
        "calibrated one beside it:\n"
        + "\n".join(f"  {name}:{number}  {field}" for name, number, field in found))


def test_the_guard_would_notice_a_new_reader():
    """The allowlist must not have grown to cover everything."""
    scanned = [p for p in REPO.glob("scripts/*.py") if p.name not in ALLOWED]
    assert len(scanned) > 20, (
        f"only {len(scanned)} scripts are being checked, so the allowlist has "
        "swallowed the guard")


def test_both_variants_are_stored_so_the_conservative_one_can_be_read():
    source = (REPO / "scripts/phase6_compensatory_mutation.py").read_text()
    for field in ("exceeds_null_conservative", "families_exceeding_null_conservative"):
        assert re.search(rf'["\']{field}["\']', source), (
            f"{field} is not written, so nothing can report it")
