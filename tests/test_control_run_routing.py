"""A randomized-weight control must reach the same code paths as its trained twin.

Every stage in `phases_1_to_5` dispatches on the key it is handed: which
tokenizer the model has, whether the architecture has attention at all, and how
NT's 6-mers and DNABERT-2's BPE tokens map back onto nucleotides. A key
carrying an `_untrained` suffix is in none of those sets and nothing raises, so
a control runs the character-tokenizer path, reports a transformer as having no
attention, and drops the offset -- and the file it writes looks like every
other file. The runner therefore keys the stages on the base model and suffixes
only the output, which is what `modal_phase6_untrained_all.py` does.
"""

import inspect
import re

import modal_repaired_panel as runner
import phases_1_to_5 as phases
from phase6_compensatory_mutation import ADAPTER_OFFSETS, NON_CHARACTER_TOKENIZERS

STAGES = ("run_mutation_sensitivity", "run_attention_contact", "run_structure_probing")

DISPATCH_SETS = {
    "RNA_TOKENIZER_MODELS": phases.RNA_TOKENIZER_MODELS,
    "ATTENTION_MODELS": phases.ATTENTION_MODELS,
    "NON_CHARACTER_TOKENIZERS": NON_CHARACTER_TOKENIZERS,
    "ADAPTER_OFFSETS": set(ADAPTER_OFFSETS),
}


def test_a_suffixed_key_silently_loses_a_transformers_attention():
    suffixed = phases.run_attention_contact(None, "nt_untrained", [])
    assert suffixed["skipped"] == "no attention (SSM architecture)", (
        "the hazard is that this returns quietly rather than raising")
    assert "skipped" not in phases.run_attention_contact(None, "nt", [])


def test_no_dispatch_set_admits_a_control_under_its_own_name():
    for name in runner.UNTRAINED_MODELS:
        base = runner.base_model(name)
        assert base in runner.TRAINED_MODELS, name
        for label, members in DISPATCH_SETS.items():
            assert name not in members, (
                f"{label} lists {name}; the runner keys on the base model, so a "
                "membership added under the suffix would never be consulted")


def test_the_runner_hands_every_stage_the_base_key():
    source = inspect.getsource(runner._run_model)
    for stage in STAGES:
        call = re.search(rf"{stage}\(adapter, (\w+),", source)
        assert call is not None, f"{stage} is no longer called positionally"
        assert call.group(1) == "key", (
            f"{stage} is handed {call.group(1)}; a control's stages must be keyed "
            "on base_model(model_name), not on the run name")
    assert re.search(r"^\s*key = base_model\(model_name\)$", source, re.M)


def test_every_control_that_carries_attention_shares_its_twins_answer():
    for name in runner.UNTRAINED_MODELS:
        base = runner.base_model(name)
        keyed = phases.run_attention_contact(None, base, [])
        assert ("skipped" in keyed) == (base not in phases.ATTENTION_MODELS)
