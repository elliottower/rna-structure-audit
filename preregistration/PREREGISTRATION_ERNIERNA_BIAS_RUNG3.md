# Does ERNIE-RNA's partner specificity come from its pairwise bias buffer or from its learned weights?

Status: FROZEN. Commit SHA: `d095d091895eed118d0cdeeffad7cc088e6b0e50`

Frozen before either ablated cell was computed. The two intact cells it
compares against were computed before this registration was written, and
their values are stated in Foreknowledge.

## Description

ERNIE-RNA scores the highest per-pair precision in the panel. On the 35 families
that carry a derangement chance rate, which is the analysis set here, it scores
0.888 against a chance rate of 0.124, an excess of 0.764. Its randomly
initialized control scores 0.490 against a chance rate of 0.241, an excess of
0.249. A control whose weights are random should carry no learned partner
information.

`multimolecule.ErnieRnaModel` holds `pairwise_bias_map`, a vocabulary-square
buffer registered `persistent=False`. After one forward pass it contains exactly
six non-zero cells: A-U and U-A at 2.0, C-G and G-C at 3.0, G-U and U-G at 0.8.
That is a Watson-Crick and wobble table, hardcoded rather than learned, and it is
added to attention through `pairwise_bias_proj`. `_randomize` in
`scripts/modal_repaired_panel.py` iterates `named_parameters()`, and a buffer is
not a parameter, so the table survives randomization untouched. Only
`pairwise_bias_proj`, 96 parameters, is randomized.

This registration asks whether that table is what the untrained control scores on.

## Hypotheses

**H1.** Ablating the bias reduces the randomly initialized model's excess over
its own derangement chance rate.

**H2.** Ablation costs the trained model less excess than it costs the untrained
model, because the trained model has learned weights to fall back on.

**H3.** The trained model's precision after ablation still exceeds its own
post-ablation derangement chance rate, so ERNIE-RNA's partner specificity does
not rest entirely on the table.

H1 is the primary question: without it the untrained control's excess of 0.249
has no attributed cause and the trained-versus-untrained comparison for ERNIE-RNA
is uninterpretable in either direction. H2 and H3 do not depend on it -- see the
manipulation check below, which replaces gating them on H1.

Each hypothesis compares a model against **its own** derangement chance rate,
computed in the same run on the same weights. A chance rate is a property of the
model, not a constant: across this panel it ranges from 0.000 to 0.316, so an
ablated model's rate cannot be assumed equal to the intact model's.

## Inference criteria

Write excess as `E = precision - chance`, both measured per arm.

| hypothesis | holds when |
|---|---|
| H1 | 95% CI on `E(untrained intact) - E(untrained ablated)` excludes zero |
| H2 | 95% CI on the difference of those two differences excludes zero, favouring the untrained arm |
| H3 | 95% CI on `E(trained ablated)` excludes zero |

Each interval is 10,000 bootstrap resamples of families, drawn once per replicate
and applied to both arms, because the arms cover the same panel. The point
estimate is reported beside every interval.

**What the design resolves, stated as power rather than as a threshold.** The
same bootstrap over the two intact arms gives a half-width of 0.126 on the
difference of excesses (`scripts/estimate_precision_resolution.py`). That is an
estimate from the closest paired comparison available before the run -- trained
against untrained -- and not from the arms the hypotheses are defined over, which
do not exist yet. Those two arms differ from each other far more than intact will
differ from ablated, so their family-level correlation need not be the same and
0.126 is a guide rather than a measurement of this design. A true
effect below that is not distinguishable from zero here, so an interval covering
zero is reported as "smaller than this design resolves" and never as evidence of
no effect. No pass/fail margin is registered: a threshold set at the detection
limit would find a true effect of that size about half the time, and pairing it
with a CI condition would test nearly the same event twice.

**Manipulation check: ablated precision must fall below intact precision in both
arms.** If it does not, the ablation did not take effect and no hypothesis is
evaluated. This replaces gating H2 and H3 on H1: the buffer explaining part but
not all of the advantage would fail H1 while leaving H2 and H3 both answerable,
and that is a likely outcome rather than a degenerate one.

**All hypotheses are void if either ablated arm carries a chance rate for fewer
than 33 of the 35 families in the analysis set**, which would mean the ablation
broke the forward pass rather than removed a prior. Both intact arms score 36 of
36 families, 35 of which carry a derangement null.

**The analysis set is the 35 families for which the derangement null is defined.**
The registered null uses only stems with at least three eligible pairs
(`PREREGISTRATION_PHASE6_V2.md`, step 5), and a family where no stem reaches
three has no null and therefore no chance rate. Exactly one non-quarantined
scored family is excluded on that rule: `c_di_GMP_riboswitch`, which the pipeline
records as `null_available: false` with 8 eligible pairs spread across its stems.
The set is pinned by the sha256 of its sorted family names,
`10bc3b36b3fe0d08`, so it can be checked rather than reconstructed.

**If an ablated arm's set differs from the intact 35, the comparison is taken on
the intersection**, and the size of the intersection is reported. The bootstrap
draws one family index per replicate and applies it to both arms, so a comparison
across arms with different sets would not be paired and the interval would not
mean what it says.

## An outcome that strengthens the model under examination

If H1 fails and H3 holds, ERNIE-RNA's partner specificity is neither the table
nor an artifact: it is learned, the untrained control's excess of 0.249 comes
from something else, and ERNIE-RNA's result is stronger than the difference of
excesses, 0.515, suggests. That outcome is reported whichever way it lands.

## Sampling plan

**Existing data.** Two of four cells are computed and will not be re-run: trained
intact (precision 0.888, chance 0.124, excess 0.764) and untrained intact (0.490,
0.241, excess 0.249), both `results/repaired_panel_v3`, commit `1350db1`, over
the 35-family analysis set.

**Data collection procedure.** Two new cells on the multimolecule image, one
container each: trained with the bias ablated, and randomly initialized with the
bias ablated. Ablation zeros `pairwise_bias_map` and every parameter of
`pairwise_bias_proj`, **and sets `model._inited = True`**. That flag is
load-bearing: `ErnieRnaEmbeddings` rebuilds the buffer on the first forward pass
when `_inited` is false, as a documented workaround for transformers v5 leaving
non-persistent buffers on the meta device. An ablation that omits it is a no-op,
which is what `scripts/modal_ernierna_ablation.py` did in July -- its Rungs 1--2
numbers measure nothing and are not used here. Randomization runs before zeroing,
under `RANDOM_INIT_SEED = 42`, so the zeroed tensors stay zero.

Both new cells also re-run Rungs 1--2, so the ablated and intact Rungs 1--2
numbers exist at one commit under one arithmetic. The July numbers are replaced,
not compared against: they predate the float64 metric, the token-span repair and
this panel.

**Sample size.** 47 analyzed families, 38 passing the registered Rung 3 filters,
36 outside the quarantine, 35 carrying a derangement null and therefore a chance
rate. The 35 are the analysis set, since every hypothesis here reads a chance
rate. It resolves a difference of excesses of 0.126 and cannot resolve one of
0.05.

## Foreknowledge of data or evidence

Substantial, and stated because it constrains what this registration can claim.

Known before writing: over the analysis set, trained precision 0.888 and
untrained 0.490, with derangement chance rates 0.124 and 0.241, giving excesses
of 0.764 and 0.249 whose difference is 0.515. Over all 36 non-quarantined scored
families the precisions are 0.864 and 0.476. The buffer is non-persistent, is
rebuilt on first forward, holds the six canonical pairs, and is identical across
model instantiations once built. Before that first forward pass it holds
uninitialized memory, which is what made an earlier inspection read NaNs and
values near 1e18; that state is transient and never reaches a measurement.

Not known: any Rung 3 number for either ablated cell, and whether the July
ablation would have changed Rungs 1--2 had it taken effect.

## Analysis plan

Precision is the mean over non-quarantined scored families of
`h3_precision.fraction`. Chance is the mean of `h3_chance_fraction` from the same
run. Intervals are 10,000 bootstrap resamples over families, matching
`N_BOOTSTRAP` and `SEED` elsewhere in this project. No new statistic is defined.

## Other

Headings of the OSF schema not answered above are N/A for one reason: this is a
computational ablation on fixed inputs.

Blinding — N/A; the analysis is deterministic given the weights.
Randomization — N/A; no units are assigned to conditions.
Missing data — N/A; a family scores or is recorded as skipped with its reason.
Exclusion criteria — N/A; the registered Rung 3 filters apply unchanged.
