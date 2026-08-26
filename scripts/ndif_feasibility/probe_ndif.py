"""Probe what NDIF will actually run, for the expanded RNA model panel.

The question: can the Rung 1-3 protocol (per-position embedding extraction under
single-nucleotide perturbation) run on NDIF instead of rented GPUs? Three things
have to hold, and each is probed separately.

  1. remote execution works with this API key
  2. NDIF provisions models absent from its published catalogue
  3. encoder / masked-LM models -- which is every RNA foundation model -- survive
     both provisioning and intervention-graph execution

Every probe is declared as data in PROBES below, so the emitted JSON carries the
repo id, the input string, the hook expression and the automodel for each row
alongside its result. A row and its provenance travel together.

Two findings drive the design:

* `LanguageModel._remoteable_model_key` serializes only {repo_id, revision}, so
  `automodel` never crosses the wire and the server builds every model as
  AutoModelForCausalLM. `shim=True` carries it inside the key instead;
  `_remoteable_from_model_key` merges the key JSON into constructor kwargs and
  TransformersMixin resolves a string automodel via getattr(modeling_auto, ...),
  so provisioning needs no server change.

* t5 is the control that separates "the shim broke it" from "encoders are
  broken". BERT registers with AutoModelForCausalLM as BertLMHeadModel, so it
  reaches the server through the stock unshimmed path.

Usage:
    python scripts/ndif_feasibility/probe_ndif.py
    python scripts/ndif_feasibility/probe_ndif.py --only t5_bert_stock_path_control
    python scripts/ndif_feasibility/probe_ndif.py --list
"""

import argparse
import json
import os
import traceback
from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional

from dotenv import load_dotenv
from huggingface_hub import HfApi

import nnsight
from nnsight import CONFIG, LanguageModel
from nnsight.intervention.backends.remote import RemoteBackend

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "ndif_feasibility"

# NDIF access was set up in the lookback audit repo; the key lives there.
DEFAULT_ENV = Path.home() / "Documents" / "GitHub" / "lookback-validity-audit" / ".env"

TEXT = "The capital of France is"
MASKED_TEXT = "The capital of France is [MASK]."
PROTEIN_WT = "MKTVRQERLKSIVRILERSKEPVSGAQ"
PROTEIN_MUT = "MKTVRQERLKSIVAILERSKEPVSGAQ"  # R -> A at residue 14


# ── the automodel shim ───────────────────────────────────────────────────────

_ID_CACHE: dict[str, str] = {}


def _keyed_with_automodel(self) -> str:
    """LanguageModel._remoteable_model_key, but carrying the automodel.

    to_model_key() derives the server-side import path from
    __func__.__module__ and __func__.__qualname__.split('.')[0], so both are
    reassigned below to keep the key pointing at the real LanguageModel class
    rather than at this module.
    """
    if self.repo_id not in _ID_CACHE:
        _ID_CACHE[self.repo_id] = HfApi().model_info(self.repo_id).id
    return json.dumps(
        {
            "repo_id": _ID_CACHE[self.repo_id],
            "revision": self.revision,
            "automodel": self.automodel.__name__,
        }
    )


_keyed_with_automodel.__module__ = "nnsight.modeling.language"
_keyed_with_automodel.__qualname__ = "LanguageModel._remoteable_model_key"

_STOCK_MODEL_KEY = LanguageModel._remoteable_model_key


def set_shim(enabled: bool) -> None:
    LanguageModel._remoteable_model_key = (
        _keyed_with_automodel if enabled else _STOCK_MODEL_KEY
    )


# ── probe specification ──────────────────────────────────────────────────────


@dataclass
class Probe:
    """One remote-execution attempt, fully described by its own fields."""

    name: str
    repo_id: str
    inputs: list[str]
    hook: str  # expression evaluated with `m` bound to the model
    why: str
    automodel: Optional[str] = None  # None -> nnsight default (AutoModelForCausalLM)
    shim: bool = False
    trace_kwargs: dict[str, Any] = field(default_factory=dict)
    expect: str = "ok"  # "ok" or "error" -- what the finding predicts
    postprocess: Optional[Callable[[list], dict]] = None


def _perturbation_summary(saved: list) -> dict:
    """Wild-type vs single-residue mutant, per-position L2 shift at the last layer."""
    hs_wt, hs_mut = saved
    shift = (hs_wt[-1].float() - hs_mut[-1].float()).norm(dim=-1)[0]
    return {
        "n_layers": len(hs_wt),
        "layer_shape": list(hs_wt[0].shape),
        "per_position_l2_shift": [round(float(x), 4) for x in shift],
        "argmax_position": int(shift.argmax()),
    }


PROBES: list[Probe] = [
    Probe(
        name="t1_hot_catalogued",
        repo_id="meta-llama/Llama-3.2-1B-Instruct",
        inputs=[TEXT],
        hook="m.model.layers[8].output[0]",
        why="remote execution works at all on a HOT catalogued model",
    ),
    Probe(
        name="t2_cold_catalogued",
        repo_id="openai-community/gpt2-medium",
        inputs=[TEXT],
        hook="m.transformer.h[6].output[0]",
        why="COLD catalogued model cold-starts on demand",
    ),
    Probe(
        name="t3_uncatalogued_decoder",
        repo_id="EleutherAI/pythia-160m",
        inputs=[TEXT],
        hook="m.gpt_neox.layers[6].output[0]",
        why="model absent from GET /status is provisioned on request, no shim",
    ),
    Probe(
        name="t4_masked_lm_stock_path",
        repo_id="facebook/esm2_t12_35M_UR50D",
        automodel="AutoModelForMaskedLM",
        inputs=[PROTEIN_WT],
        hook="m.esm.encoder.layer[6].nns_output[0]",
        why="masked LM through the stock path; automodel is dropped on the wire",
        expect="error",
    ),
    Probe(
        name="t4a_masked_lm_shimmed",
        repo_id="facebook/esm2_t12_35M_UR50D",
        automodel="AutoModelForMaskedLM",
        shim=True,
        inputs=[PROTEIN_WT],
        hook="m.esm.encoder.layer[6].nns_output[0]",
        why="shim fixes provisioning; does graph execution follow?",
        expect="error",
    ),
    Probe(
        name="t4b_second_encoder_family",
        repo_id="google-bert/bert-base-uncased",
        automodel="AutoModelForMaskedLM",
        shim=True,
        inputs=[MASKED_TEXT],
        hook="m.bert.encoder.layer[6].nns_output[0]",
        why="unrelated encoder family, to show the failure is not ESM-specific",
        expect="error",
    ),
    Probe(
        name="t5_bert_stock_path_control",
        repo_id="google-bert/bert-base-uncased",
        inputs=[MASKED_TEXT],
        hook="m.bert.encoder.layer[6].nns_output[0]",
        why="CONTROL: BERT registers as BertLMHeadModel so it needs no shim. "
        "If this fails too, the shim is exonerated and encoders are the fault.",
        expect="error",
    ),
    Probe(
        name="t6_all_layer_extraction",
        repo_id="facebook/esm2_t12_35M_UR50D",
        automodel="AutoModelForMaskedLM",
        shim=True,
        inputs=[PROTEIN_WT, PROTEIN_MUT],
        hook="m.output.hidden_states",
        trace_kwargs={"output_hidden_states": True},
        why="the operation the RNA protocol needs: all-layer embeddings under a "
        "point mutation, mirroring ModelAdapter.get_all_layer_embeddings",
        expect="error",
        postprocess=_perturbation_summary,
    ),
]


# ── runner ───────────────────────────────────────────────────────────────────


def setup(env_path: Path, read_timeout: float) -> None:
    load_dotenv(env_path)
    if "NDIF_KEY" not in os.environ:
        raise RuntimeError(f"NDIF_KEY not found in {env_path}")
    CONFIG.APP.REMOTE_LOGGING = False
    CONFIG.set_default_api_key(os.environ["NDIF_KEY"])
    RemoteBackend.CONNECT_TIMEOUT = 30.0
    RemoteBackend.READ_TIMEOUT = read_timeout


def build_model(probe: Probe) -> LanguageModel:
    kwargs = {}
    if probe.automodel is not None:
        from transformers import models as _tm  # noqa: F401
        from transformers.models.auto import modeling_auto

        kwargs["automodel"] = getattr(modeling_auto, probe.automodel)
    return LanguageModel(probe.repo_id, **kwargs)


def run_probe(probe: Probe) -> dict:
    set_shim(probe.shim)
    model = build_model(probe)
    saved = []
    for text in probe.inputs:
        # The .save() must happen inside the trace, the append outside it:
        # nnsight compiles the with-body into a graph rather than executing it,
        # so mutating a local container in there silently does nothing.
        with model.trace(text, remote=True, **probe.trace_kwargs):
            value = eval(probe.hook, {"m": model}).save()  # noqa: S307
        saved.append(value)

    out: dict[str, Any] = {"model_key_sent": model.to_model_key()}
    if probe.postprocess is not None:
        out.update(probe.postprocess(saved))
    else:
        first = saved[0]
        out["shape"] = list(first.shape)
        out["dtype"] = str(first.dtype)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--env", type=Path, default=DEFAULT_ENV)
    ap.add_argument("--read-timeout", type=float, default=300.0)
    ap.add_argument("--only", nargs="*", default=None, help="probe names to run")
    ap.add_argument("--list", action="store_true", help="print probes and exit")
    args = ap.parse_args(argv)

    if args.list:
        for p in PROBES:
            print(f"{p.name:<30} {p.repo_id:<45} expect={p.expect}")
        return 0

    setup(args.env, args.read_timeout)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"probe_{stamp}.json"

    record: dict[str, Any] = {
        "timestamp_utc": stamp,
        "client": {"nnsight": nnsight.__version__, "ndif_host": CONFIG.API.HOST},
        "probes": {},
    }

    for probe in PROBES:
        if args.only and probe.name not in args.only:
            continue

        spec = asdict(probe)
        spec.pop("postprocess", None)
        entry: dict[str, Any] = {"spec": spec}

        try:
            entry["status"] = "ok"
            entry["result"] = run_probe(probe)
        except Exception as exc:  # noqa: BLE001 -- the error text is the datum
            entry["status"] = "error"
            entry["error_type"] = type(exc).__name__
            entry["error"] = str(exc)
            entry["traceback"] = traceback.format_exc()

        entry["matches_expectation"] = entry["status"] == probe.expect
        record["probes"][probe.name] = entry

        # write after every probe: a hang or a kill must not lose finished work
        out_path.write_text(json.dumps(record, indent=2))

        flag = "OK  " if entry["status"] == "ok" else "FAIL"
        detail = (
            json.dumps(entry["result"])[:80]
            if entry["status"] == "ok"
            else entry["error"].strip().splitlines()[-1][:80]
        )
        print(f"  {flag} {probe.name:<30} {detail}")

    set_shim(False)

    n = len(record["probes"])
    agree = sum(1 for v in record["probes"].values() if v["matches_expectation"])
    print(f"\n{agree}/{n} probes matched the documented expectation")
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
