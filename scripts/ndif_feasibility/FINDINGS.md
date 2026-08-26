# Can the expanded RNA panel run on NDIF?

Reproduce with `scripts/ndif_feasibility/probe_ndif.py`; machine-readable output in
`results/ndif_feasibility/probe_<timestamp>.json`. Client nnsight 0.7.0 against
`https://api.ndif.us`.

## Result

| probe | shim | outcome |
|---|---|---|
| t1 HOT catalogued model (Llama-3.2-1B) | no | ok — `(6, 2048)` bfloat16 |
| t2 COLD catalogued model (gpt2-medium) | no | ok — `(5, 1024)`, cold-started on demand |
| t3 uncatalogued decoder (pythia-160m) | no | ok — `(5, 768)` |
| t4 masked LM, stock path (ESM-2 35M) | no | provisioning refused: `Unrecognized configuration class EsmConfig for this kind of AutoModel: AutoModelForCausalLM` |
| t4a masked LM, `automodel` in key (ESM-2 35M) | yes | provisions, then `RemoteException: name 'hooked_output' is not defined` |
| t4b second encoder family (bert-base-uncased) | yes | same failure |
| t5 encoder through stock path (bert-base-uncased) | **no** | same failure |
| t6 all-layer extraction + point mutation | yes | same failure |

## What each result establishes

**The published catalogue is not the constraint.** t3 provisions
`EleutherAI/pythia-160m`, which appears nowhere in the 118 deployments returned by
`GET https://api.ndif.us/status`, with no shim and no prior request. Arbitrary
HuggingFace decoder models are served on demand today.

**The stock remote path hardcodes a causal-LM head.**
`LanguageModel._remoteable_model_key` serializes `{repo_id, revision}` only, so the
`automodel` argument never crosses the wire and the server builds every model as
`AutoModelForCausalLM`. Every RNA foundation model is a masked LM, so every one of
them fails at t4.

**That part is fixable client-side.** `_remoteable_from_model_key` merges the key
JSON into the constructor kwargs, and `TransformersMixin.__init__` resolves a string
`automodel` through `getattr(modeling_auto, automodel)`. Adding `automodel` to the
key is sufficient to provision a masked LM, with no server change — t4a gets past
provisioning. `to_model_key` derives the server-side import path from
`__func__.__module__` and `__func__.__qualname__.split('.')[0]`, so a patch must set
both to keep the key pointing at the real `LanguageModel`.

**Encoder models then fail during graph execution, and the shim is not the cause.**
t5 is the control: BERT registers with `AutoModelForCausalLM` as `BertLMHeadModel`,
so it reaches the server through the unmodified path — and fails identically. Two
unrelated encoder families (ESM-2, BERT) fail the same way, shimmed and unshimmed,
at every hook point tried (`embeddings.output`, `lm_head.output`,
`encoder.nns_output[0]`, `encoder.layer[i].nns_output[0]`, `output.hidden_states`).
Decoder models at the same hook depth succeed. `hooked_output` is an nnsight
internal, not a name in any code here.

## Consequence for the panel

Not usable for the 26-model expansion until encoder remote execution works. The
panel stays on rented GPUs, which is where `scripts/modal_*.py` already runs it, and
is the right tool regardless: the largest addition (AIDO.RNA 1.6B) fits on one
GPU, so NDIF's tensor-parallel capacity buys nothing here.

Worth reporting upstream. The failure is one line, reproduces on a 110M public model
in isolation, and blocks every encoder architecture — BERT, ESM, and every RNA and
protein language model — from remote execution.

## Second-order note

Even with remote execution fixed, models loaded through the `multimolecule` package
(RiNALMo, ERNIE-RNA, SpliceBERT, UTR-LM, RNA-FM, RNABERT, RNAErnie, RNA-MSM,
3UTRBERT — nine of the expanded panel, including both Rung 3 passers) declare
architectures like `RiNALMoForPreTraining` with no `auto_map`. They need
`multimolecule` installed in the serving image, not `trust_remote_code`. That is a
smaller ask than it looks: one pip dependency, no arbitrary repo code.

## Local reproduction, no cluster needed

`nnsight_repro/test_envoy_overloaded_mount.py` reproduces the encoder failure
offline on two tiny public models in ~14 s:

```
3 failed, 3 passed, 2 skipped
```

Passing: decoder layer envoy is the base `Envoy` class and pickles cleanly, and
`.nns_output` stays mounted on the encoder (so any fix must keep that).
Failing: the encoder's synthesized class is named `Envoy.Preserved` — a dot in
`__name__`, never bound in its defining module, and therefore unpicklable:

```
PicklingError: Can't pickle <class 'nnsight.intervention.envoy.Envoy.Preserved'>:
attribute lookup Envoy.Preserved on nnsight.intervention.envoy failed
```

The two remote tests are skipped without `NDIF_KEY` and fail with the
`hooked_output` error when it is set.

## Fix submitted upstream

[ndif-team/nnsight#700](https://github.com/ndif-team/nnsight/pull/700) carries both
fixes and the test file above.

`_handle_overloaded_mount` now memoizes the synthesized class on
`(base class, mount point)` and binds it into the module that defines the base.
The attributes it installs are fully determined by that pair, so one shared class
is equivalent to one per instance, and per-instance state stays in `__dict__`. A
first attempt that only renamed the class was not sufficient: one class is created
per colliding envoy instance, so a single module-level name collides across them
(`PicklingError: it's not the same object as ...`). Memoizing is what makes the
name stable.

The second commit adds `automodel` to the remote model key, omitted when it is the
default so keys for existing causal-LM deployments stay byte-identical.

Until that lands and NDIF deploys it, the panel runs on rented GPUs regardless.
Even with remote execution working, `multimolecule` still has to be present in the
serving image for nine of the models.

## The fix, verified on Linux against real RNA models

`modal_verify_patch.py` installs one checkout, then runs identical checks on
`main` and on a branch carrying both fixes, switching with `git checkout -f`
between them. `main` is the negative control: a suite that passes on the patch
proves nothing unless the same suite fails without it. Raw output in
`results/ndif_feasibility/modal_verify_branches.json`.

| | `main` (87f4dad) | both fixes (e9ae3ef) |
|---|---|---|
| nnsight `tests/test_tiny.py` | pass | pass |
| encoder repro suite | 3 failed, 3 passed | **6 passed** |
| RNABERT envoy class | `Envoy.Preserved` | `Envoy__nns_output` |
| SpliceBERT envoy class | `Envoy.Preserved` | `Envoy__nns_output` |
| picklable | **no** | **yes** |
| distinct classes across 6 layers | 6 | 1 |

Their own suite passes on both branches, so the change is not a regression.

The computation is untouched. Per-position L2 shift at the last layer under a
single-nucleotide substitution is identical to four decimal places on both
branches, for both models, and peaks at the mutated position:

```
RNABERT     positions 8-12   0.3140  0.0797  0.8089  0.2187  0.1411   argmax 10
SpliceBERT  positions 8-12   3.3550  2.4945  8.5443  2.3225  2.1158   argmax 10
```

Tracing through `.nns_output` returns `(20, 120)` and `(20, 512)` on both.

RNA-FM fails on both branches with `ValueError: vocab_size (28) must be 26 when
codon=False`, a multimolecule/transformers version incompatibility unrelated to
this change. The panel does not load RNA-FM through multimolecule.

What this cannot show: the NDIF server runs its own nnsight, so `hooked_output`
will persist there until the fix is deployed. Verified here is that the patch
removes the proximate cause, passes the maintainers' suite, and leaves numerical
output unchanged on the models this panel uses.
