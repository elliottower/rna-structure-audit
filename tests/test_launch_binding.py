"""A launch must be refused unless the files about to run are the files that were tested."""

import json

import pytest

import modal_decoy_panel as launcher

PASSED = {"tests": 14, "failures": 0, "errors": 0, "skipped": 0, "passed": 14}


@pytest.fixture
def repo(tmp_path, monkeypatch):
    for path in launcher.TESTED_FILES:
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(f"contents of {path}\n")
    record = tmp_path / "results" / "audit" / "model_side_tests.json"
    record.parent.mkdir(parents=True)
    monkeypatch.setattr(launcher, "REPO", tmp_path)
    monkeypatch.setattr(launcher, "TEST_RECORD", record)
    return tmp_path


def write_record(counts=PASSED, returncode=0):
    launcher.TEST_RECORD.write_text(json.dumps(full_record(counts, returncode)))


def full_record(counts=PASSED, returncode=0):
    return {"returncode": returncode, "counts": counts,
            "sha256": launcher.file_hashes(launcher.TESTED_FILES),
            "recorded": "2026-10-11T00:00:00+00:00", "head_commit": "0" * 40,
            "libraries": {"torch": "2.6.0"}}


def test_a_passing_record_over_unchanged_files_is_accepted(repo):
    write_record()
    assert launcher.require_tested_files()["counts"] == PASSED


def test_no_record_is_refused(repo):
    with pytest.raises(launcher.UntestedLaunch, match="no test record"):
        launcher.require_tested_files()


@pytest.mark.parametrize("path", launcher.TESTED_FILES)
def test_one_changed_byte_in_any_tested_file_is_refused(repo, path):
    write_record()
    (repo / path).write_text(f"contents of {path}\n ")
    with pytest.raises(launcher.UntestedLaunch, match=path):
        launcher.require_tested_files()


@pytest.mark.parametrize("field", ["failures", "errors", "skipped"])
def test_a_run_with_any_test_not_passing_is_refused(repo, field):
    write_record({**PASSED, field: 1})
    with pytest.raises(launcher.UntestedLaunch, match="did not pass in full"):
        launcher.require_tested_files()


def test_a_run_that_collected_nothing_is_refused(repo):
    write_record({"tests": 0, "failures": 0, "errors": 0, "skipped": 0, "passed": 0})
    with pytest.raises(launcher.UntestedLaunch, match="did not pass in full"):
        launcher.require_tested_files()


def test_a_nonzero_exit_is_refused_whatever_the_counts_say(repo):
    write_record(returncode=1)
    with pytest.raises(launcher.UntestedLaunch, match="did not pass in full"):
        launcher.require_tested_files()


def test_the_launcher_and_both_model_side_files_are_among_the_tested_files():
    for path in ("scripts/modal_decoy_panel.py",
                 "scripts/phase6_compensatory_mutation_v3.py",
                 "scripts/rnafm_adapter_v2.py"):
        assert path in launcher.TESTED_FILES


def test_every_tested_and_registered_file_exists_in_the_repository():
    root = launcher.Path(launcher.__file__).resolve().parent.parent
    for path in launcher.TESTED_FILES + launcher.REGISTERED_FILES:
        assert (root / path).is_file(), path


def test_the_panel_hash_moves_when_one_record_changes(repo):
    families = repo / "data" / "rfam_families"
    families.mkdir(parents=True)
    (families / "a.json").write_text("{}")
    before = launcher.panel_tree_hash()
    (families / "a.json").write_text("{ }")
    assert launcher.panel_tree_hash() != before


@pytest.mark.parametrize("missing", ["returncode", "counts", "sha256", "recorded",
                                     "head_commit", "libraries"])
def test_a_record_missing_any_field_is_refused_as_untested(repo, missing):
    record = full_record()
    del record[missing]
    launcher.TEST_RECORD.write_text(json.dumps(record))
    with pytest.raises(launcher.UntestedLaunch, match="not one --test-only wrote"):
        launcher.require_tested_files()


def test_a_record_that_is_not_json_is_refused_as_untested(repo):
    launcher.TEST_RECORD.write_text("passed")
    with pytest.raises(launcher.UntestedLaunch, match="not one --test-only wrote"):
        launcher.require_tested_files()


def test_flag_checks_come_before_the_manifest_is_written():
    source = launcher.Path(launcher.__file__).read_text()
    main = source[source.index("def main("):]
    assert main.index("--transversion runs the mutation stage") < main.index("manifest = {")
    assert main.index("--ablate-bias zeroes") < main.index("manifest = {")
    assert main.index("Unknown models") < main.index("manifest = {")


def test_the_manifest_binds_the_whole_test_record_by_hash():
    source = launcher.Path(launcher.__file__).read_text()
    assert '"test_record_sha256": hashlib.sha256(TEST_RECORD.read_bytes()).hexdigest()' in source
