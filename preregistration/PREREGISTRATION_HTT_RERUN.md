# Does repairing the pipeline change what the HTT case study shows?

Status: FROZEN. Commit SHA: `08de258e24794bc7c35f41f521ccf783c2874219`

Frozen before the HTT fragments are recomputed. The numbers in the manuscript's
HTT table were produced on 2026-07-13, before the token-span repair (D11--D14),
the float64 metric, and the RNA-FM port fix (D17).

## Description

Five HTT exon 1 fragments with CAG repeat counts 17, 21, 36, 40 and 60, folded
with ViennaRNA, carry a case study whose claim is about probe design rather than
about these five sequences: four models reach at least 97% linear-probe accuracy
for stem-versus-loop classification while showing no composition-controlled
mutation sensitivity. CAG repeats contain no U, so the composition bias a probe
can exploit is extreme.

The section's numbers were computed by the pipeline that carried D11 through D17.
Whether its claim survives recomputation is not known.

## Hypotheses

**H1.** RNA-FM's mean stem--loop ratio changes by more than 0.1, because its
deposited value came from a port that added an unseeded random vector to every
embedding and dropped the model's final layer norm (D17).

**H2.** NT v2 reports a mean ratio and a families-above-null count where the
manuscript prints a tokenizer caveat, because the positions it scores are now
placed by the tokenizer rather than by a closed-form guess (D11, D13).

**H3.** The qualitative claim holds: at least three of the five models with
probing accuracy above 0.97 exceed the composition-preserving null in no more
than one of the five fragments.

H3 carries the section. H1 and H2 are predictions about which numbers move, and
neither changes what the section argues.

## Inference criteria

| hypothesis | holds when |
|---|---|
| H1 | \|repaired RNA-FM ratio − 1.794\| > 0.1 |
| H2 | NT v2 yields a finite ratio and a families-above-null count |
| H3 | at least three models with probing > 0.97 exceed the null in ≤ 1 of 5 fragments |

**Amendment, recorded before the run.** The disposition below replaces a single
withdrawal condition on H3, which was stricter than the claim requires. H3 is a
statement about probe design -- CAG repeats contain no U, so a linear probe has
an extreme composition bias to exploit -- and one model behaving differently
does not touch it.

| models with probing > 0.97 that stay clean | disposition |
|---|---|
| 3 or 4 of 4 | the section stands as written |
| 1 or 2 of 4 | the section is rewritten around the models that stay clean, and names the ones that do not |
| 0 of 4 | the section is withdrawn |

**Withdrawal requires that no high-accuracy model stays clean**, because then
high probing accuracy and composition-controlled signal coincide and the section
has nothing to illustrate. The main panel makes the same point across
\panelN{} families and ten models, so nothing is lost by removing it.

## An outcome that strengthens what is under examination

If H2 holds and NT v2's repaired numbers show composition-controlled signal on
HTT hairpins where the deposited run showed none, the section reports it. That
would weaken the argument that probing accuracy is uninformative here, and it is
registered so it cannot be quietly dropped.

## Sampling plan

Five fragments, nine models. Sequences are reconstructed as
`HTT_LEFT_FLANK + "CAG" * n + HTT_RIGHT_FLANK` from
`scripts/modal_htt_all_models.py`, which produced the deposited run; structures
are ViennaRNA MFE predictions, as before. The deposited artifacts record repeat
counts, sequence lengths and MFE values but not the sequences or dot-brackets,
so both are regenerated.

**Sample size.** Five sequences. No interval computed on five fragments
separates models, and none is registered. The section's claim is a mechanism
illustration and remains one whatever the numbers.

## Foreknowledge of data or evidence

The deposited numbers are known and are in the manuscript: RNA-FM ratio 1.794
with 2 of 5 fragments above null, RiNALMo 1.031 with 0 of 5, ERNIE-RNA 1.086
with 0 of 5, SpliceBERT 0.996 with 0 of 5, NT v2 reported as a tokenizer caveat.
Probing accuracies are 0.813, 0.985, 0.985 and 0.977.

Not known: any repaired value.

## Analysis plan

`run_mutation_sensitivity`, `run_attention_contact` and `run_structure_probing`
from `scripts/phases_1_to_5.py`, unchanged, on the five fragments. The same
functions the panel uses.

## Other

Headings of the OSF schema not answered above are N/A: this is a computational
re-run on regenerated inputs.

Blinding — N/A; the analysis is deterministic given the weights.
Randomization — N/A; no units are assigned to conditions.
Missing data — N/A; a fragment scores or is recorded as skipped with its reason.
Exclusion criteria — N/A; all five fragments are scored.
