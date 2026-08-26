# Open defects

Known problems in the deposited data, results and manuscripts that are not yet
fixed. Resolved items move to `docs/CHANGELOG.md`; registration consequences go
to `docs/DEVIATIONS.md`. Each entry names the script that detects it, so the
list can be re-derived rather than trusted.

---

## D1. Six families in Rungs 1-2 carry a structure that is not their sequence's

**Detector:** `scripts/audit_annotation_validity.py`

An annotation that belongs to its sequence pairs almost only Watson-Crick and GU.
Across the 41 records that name the Rfam seed member they came from, the worst
non-canonical fraction is 18%. Across the 11 records whose `source` is
`Phase 1 curated (PDB/Rfam)` and names no accession, nine run from 24% to 77%.
The separation is complete: every family with a defective annotation is a family
with no accession behind it.

Three of the nine are already excluded post hoc for unbalanced brackets
(`U2_snRNA_stem`, `hammerhead_ribozyme`, `RNaseP_specificity`) and two are
quarantined (`tRNA_Phe_yeast`, `tRNA_Ala_human`). Six are inside the reported
analysis:

| family | non-WC pairs | annotated pairs |
|---|---|---|
| `SRP_RNA_helix8` | 77% | 31 |
| `5S_rRNA_ecoli` | 70% | 30 |
| `HDV_ribozyme` | 58% | 26 |
| `TPP_riboswitch` | 57% | 28 |
| `IRES_HCV_domainII` | 55% | 29 |
| `SAM_riboswitch` | 52% | 31 |

All six are in the 49 families reported for Rungs 1-2, where the metric is
stem-versus-loop discrimination and the positions the annotation calls paired are
positions that cannot pair. None reaches Rung 3: the requirement of 15 canonical
WC pairs excludes them, so the confirmatory set at N = 32 never contained them.

This is the same defect as the four families corrected in `5fc8914`, which were
detectable only because three of them were byte-identical. Duplication was the
symptom; the missing accession is the cause.

**Consequence for the reported result.** Rungs 1-2 report that no model exceeds
the composition-controlled null. A null result computed partly over annotations
that do not describe their sequences is weaker evidence for that conclusion, not
stronger, and the direction of the effect is not knowable without re-annotating.

---

## D2. Eleven families have no accession behind them

**Detector:** `scripts/check_panel_composition.py`

Forty-one records name their Rfam seed member with an accession. Eleven carry
`source: "Phase 1 curated (PDB/Rfam)"` and no identifier, so the sequence cannot
be traced to a record. `IRES_HCV_domainII` additionally carries `rfam_id: null`.

`paper_v12.tex:229` states that every sequence is "drawn from the Rfam seed
alignment, selected as the seed sequence with highest bit score." For these
eleven that cannot be checked, and four of the five families named as microRNA
precursors were wrong until `5fc8914`. Several of the eleven are named as
structural fragments rather than families (`SRP_RNA_helix8`, `U2_snRNA_stem`,
`RNaseP_specificity`, `IRES_HCV_domainII`), which are plausibly PDB-derived and
would not have a seed accession — but the manuscript sentence does not allow for
that, and nothing in `data/` records which origin applies to which family.

---

## D3. The printed panel statistics match no version of the panel

**Detector:** `scripts/check_panel_composition.py`

`paper_v12.tex:219` prints lengths of 22-301 nt, median 89, IQR 56-132. The
deposited panel gives 30-387 nt, median 105, IQR 79.75-154.5, and the pre-fix
panel gave 30-387, median 113.5, IQR 85.75-154.5. The printed figures match
neither, so they predate the Rfam expansion.

This is the same failure as the stem/loop GC pair at `paper_v12.tex:77-81` and
`:405`, where 68.5% / 43.9% is the 12-family pilot figure from
`PREREGISTRATION_STRUCTURE_METRICS.md:79-80` and was never updated. Both are
pilot-era numbers surviving into a manuscript describing a 52-family panel.

**Fix.** `scripts/generate_panel_description.py` writes every one of these figures to `paper/generated/panel_description.tex` as macros, computed from the records `load_rfam_families` loads. The manuscript takes them by `\input` from v13 on, so a hand-carried number cannot survive the next correction to the panel. The current panel is 47 analyzed families of 52 curated: 30--387 nt, median 106, IQR 77.5--166; stem GC 59.1%, loop GC 43.0%, stem above loop in 40 of 47. The defect stays open until v13 is built.

---

## D4. Stored results were computed against superseded annotations

**Detector:** `scripts/audit_duplicate_families.py`, `scripts/scope_rerun.py`

Every file in `results/` predates `5fc8914`. `data/` and `results/` are
inconsistent until the four corrected families are re-evaluated: 4 families
against 43 stored model-configuration files. The other 48 families are unchanged
and their records stay valid, because every null in the pipeline is computed
within a family and no aggregate pools across families.

---

## D5. `r_b = -0.28, p = 0.265` does not reproduce

**Detector:** `scripts/check_reported_values.py`

H2's rank-biserial correlation is printed as -0.28 with p = 0.265 in
`paper_v10.tex:743,806`, `paper_v11.tex:750`, `paper_v12.tex:750` and the
submitted `rna-structure-audit_v14.tex:328`. `generate_table5_registered.py`
gives +0.04, p = 1.000 at both N = 32 and N = 34. Where -0.28 came from is not
recoverable from the stored results. The FAIL verdict is unchanged either way.

`scripts/audit_h2_disposition.py` reaches +0.040, p = 1.000 by an independent
path, reading per-model mean PS from the phase 6 result files rather than from
the table generator, so the discrepancy is in the printed value and not in one
script. The ranking it prints also shows why the effect is near zero: the
RNA-pretrained group holds the two highest positions and three of the four
lowest, so the groups are bimodal rather than shifted, and a rank test on five
against five cannot register that.

---

## D6. One DNABERT-2 run is unseeded

**Detector:** `scripts/audit_duplicate_families.py`

DNABERT-2 returned 0.9304, 0.9449 and 0.9034 for three byte-identical inputs, so
that run carries an unseeded stochastic component. No other model shows this.

---

## D7. Two attention correlations have no stored run

The eager-attention recomputation for NT v2 (0.322) and the attention run for
untrained RNA-FM (0.061) are printed with no deposited per-family file. The
Availability section names both rather than claiming full deposition.

---

## D8. `PREREGISTRATION_PHASE6_V2.md` refers to a section it does not contain

Line 29 points at "Families Pending". No such section exists. The eligibility
criteria are stated in full at lines 29-37 and 128-131, so the plan is complete
and only the cross-reference is dangling. The frozen document is not edited; see
`docs/DEVIATIONS.md`.

---

## D9. The printed class breakdown describes a different panel

**Detector:** `scripts/audit_panel_classes.py`

`paper_v12.tex:215-218` gives the composition as tRNAs (7), rRNAs (4), ribozymes
(6), riboswitches (8), snRNAs (5), cis-regulatory elements (9), miRNA precursors
(5), CRISPR repeats (3), other ncRNAs (5). Five of those counts are wrong, and
not as a matter of classification judgment — the filenames fix the membership:

| class | printed | deposited |
|---|---|---|
| tRNAs | 7 | 2 |
| rRNAs | 4 | 1 |
| riboswitches | 8 | 14 |
| CRISPR repeats | 3 | 1 |
| miRNA precursors | 5 | 4 |

The errors cancel and the list still sums to 52, which is how it survived
proofreading. Same failure as D3: a description of the pilot panel carried into a
manuscript about the expanded one.

**Fix.** Generated with D3, from the same script. The analyzed panel is riboswitches (14), cis-regulatory elements (8), other ncRNAs (8), ribozymes (6), miRNA precursors (4), snRNAs (4), tRNAs (2), and CRISPR repeats (1). Eight classes, not nine: `5S_rRNA_ecoli` was the only rRNA family and is withdrawn, so the panel carries no ribosomal RNA. `paper/generated/panel_table.tex` lists all 52 records with class, Rfam accession, seed member, length and disposition.

---

## D10. The Rung 1 null rejects more often than the 5% it is set at

**Detector:** `scripts/audit_untrained_false_positives.py`

Seven models were run with randomly initialized weights. A random network's
stem-versus-loop ratio is not knowledge of structure, so the rate at which those
runs exceed their own composition-preserving null estimates the null's
false-positive rate against a registered 5%. Pooled over 329 family-runs it is
30, or 9.1%; exact binomial, one-sided, p = 1.3e-3.

The excess is not spread evenly across the seven:

| stratum | exceedances | rate |
|---|---|---|
| the five character-level tokenizers | 14/235 | 0.060 |
| NT v2 | 10/47 | 0.213 |
| DNABERT-2 | 6/47 | 0.128 |

DNABERT-2 is D11 and its Rung 1 values mean nothing either way. NT v2 is not:
its expansion misplaces 1.2% of nucleotides and never compares the wrong pair.
What NT v2 has is a six-nucleotide token. The null permutes stem and loop labels
within nucleotide strata, which scatters them; the annotation does not, because
stems and loops are contiguous runs. An embedding that averages six consecutive
positions therefore separates the real labels from the permuted ones through run
structure alone, and the wider the token the larger the separation. Any
positional encoding or convolutional receptive field does the same thing more
weakly; the token width is what makes it visible at 21%.

The five character-level models sit at 0.060 with a one-sided 95% lower bound of
0.036. That is consistent with nominal and does not establish it: 235 trials do
not separate 0.05 from 0.07.

**Consequence for the reported result.** A null whose threshold is too low
inflates every exceedance count at Rungs 1 and 2, so the discrimination claims
weaken and the negative claims strengthen. H6b (RNA-FM exceeds the nucleotide
null in at least 8 families, observed 0/47) and H7 (at least one family survives
the dinucleotide null, observed 30 for ERNIE-RNA) are the registered decisions
that read these counts.

**Fix.** A null that preserves run structure: circular rotation of the label
vector, or a permutation restricted to preserve run lengths. Per-position cosine
distances are not stored -- the deposited files carry `best_ratio`, `best_layer`,
the two null thresholds and the stem and loop counts, and nothing per layer or
per position -- so recomputing against a different null needs the embeddings
again. It is one GPU pass over all 17 runs, and it is the same pass as D11 and
D12.

---

## D11. DNABERT-2's per-position embeddings are placed by a guess about its tokenizer

**Detector:** `scripts/audit_token_alignment.py`

Rungs 1 and 2 read the cosine distance at the mutated nucleotide. Eight of the
ten models emit one token per nucleotide, so that position is a row of the hidden
state. NT v2 and DNABERT-2 do not, and `scripts/phases_1_to_5.py:61-76` expands
their token embeddings back to nucleotides with two closed-form guesses: `i // 6`
for NT v2, which assumes the sequence tiles into non-overlapping 6-mers, and
`i * n_tokens / n_nucleotides` for DNABERT-2, which assumes every byte-pair token
is the same width. Recovering each token's true span by decoding it -- both
tokenizers emit literal nucleotide strings, and the reconstruction is checked
against the sequence -- gives, over the 47 analyzed families and their 6,364
positions:

| | NT v2 | DNABERT-2 |
|---|---|---|
| token widths | 6 nt (1042 tokens), 1 nt (112) | 1 to 8 nt |
| nucleotides given a token that does not contain them | 74 (1.2%) | 3,361 (52.8%) |
| mutations that change the token count | 0 | 1,509 (23.7%) |
| mutations comparing rows that describe different stretches | 0 | 2,923 (45.9%) |

The two models fail differently. NT v2 tiles the sequence into 6-mers and pads
the remainder with single-nucleotide tokens, so `i // 6` is exact inside the
tiling and clamps at the tail; a complement substitution never changes the token
count, so both arms carry the identical misplacement and the metric compares one
position of the molecule against itself at six-nucleotide resolution. The values
are attenuated.

DNABERT-2's byte-pair tokens run 1 to 8 nucleotides and a substitution
re-segments the sequence, so in nearly half its trials the metric subtracts a
mutant row from a wild-type row describing a different stretch of the molecule.
That is not error around a true value; there is no pair.

**Consequence for the reported result.** DNABERT-2's Rungs 1 and 2 are
uninterpretable rather than imprecise: mean ratio 1.003, 8 of 47 exceeding the
nucleotide null, 5 of those 8 surviving the dinucleotide null. Untrained
DNABERT-2's mean ratio of 8.416 is not explained by this defect and is not
explained by anything else either -- the expansion does not read the weights, so
it should corrupt both arms alike, and trained DNABERT-2 sits at 1.003. NT v2's
values carry a bounded attenuation and stand with the caveat.

Rung 3 does not call the expansion and already flags both models through
`NON_CHARACTER_TOKENIZERS` (`scripts/phase6_compensatory_mutation.py:51`). Rungs
1 and 2 have no such flag, so the caveat is applied where the defect is smallest
and omitted where it is largest.

**Fix.** Assign nucleotide *i* to the token whose character span contains it,
from the tokenizer's own offsets. Needs the embeddings, so it goes in the pass
with D10 and D12.

---

## D12. The attention-contact map for DNABERT-2 is built on the same 6-mer assumption

**Detector:** `scripts/audit_token_alignment.py`

The attention rung maps the other way, from the nucleotide contact map onto token
pairs, and `_aggregate_contacts_to_tokens` (`scripts/phases_1_to_5.py:79-89`)
gives token *t* the window `[6t, 6t+6)` for both non-character models. DNABERT-2
averages 4.7 nucleotides per token, so its token rows outrun the sequence:

| | NT v2 | DNABERT-2 |
|---|---|---|
| token rows | 1,154 | 1,361 |
| rows whose window starts past the end of the sequence | 74 (6.4%) | 281 (20.6%) |
| rows whose window shares no nucleotide with the token | 74 (6.4%) | 1,207 (88.7%) |
| mean overlap between window and token | 0.936 | 0.074 |

**Consequence for the reported result.** DNABERT-2's attention-contact
correlation of 0.199 is a correlation against a contact map that describes, for
89% of its rows, a part of the molecule the token does not cover. Two registered
hypotheses are decided on it and both currently read PASS: H21 (DNABERT-2 trained
rho < 0.20, observed 0.199, under the threshold by 0.001) and H20 (trained minus
untrained below 0.05, observed 0.009). The post hoc architectural bound in
`scripts/generate_results_tables.py` reads the same delta. NT v2's 0.242, which
decides H11, carries a mean window overlap of 0.936 and moves little.

The direction is worth stating plainly: correcting this can only put two PASS
verdicts at risk, and H21 passes by 0.001.

**Fix.** Build the token contact map from the token spans, in the pass with D10
and D11.

---

## D13. Rung 3 reads the hidden state at the nucleotide index, whatever the tokenizer emits

**Detector:** `scripts/audit_rung3_token_alignment.py`

`compute_delta_profiles` (`scripts/phase6_compensatory_mutation.py:153`) reads
`emb[k + offset]` for nucleotide *k*, with `offset = 0` for every model in the
panel. For the eight models that emit one token per nucleotide a row is a
position and the read is right. NT v2 puts six nucleotides in a token, so a
106-nucleotide family has about eighteen rows; the partner position *j* of a
stem pair sits in the second half of the molecule and is past the end of that
array. The bounds test then drops the pair, and a family that loses all of its
pairs is written out as `"no valid PS values"` -- the same `skipped` field a
registered filter writes, so the loss reads as a filter in the stored file.

Over the 38 of 47 analyzed families that clear the registered Rung 3 filters:

| | NT v2 | DNABERT-2 | RiNALMo |
|---|---|---|---|
| families keeping at least one pair | 4 | 8 | 38 |
| eligible pairs scored | 16 of 583 (2.7%) | 28 of 583 (4.8%) | 583 of 583 |
| median distance from a scored position to the nearer edge of the token read for it | 110 nt | 124 nt | 0 |
| largest such distance | 230 nt | 264 nt | 0 |

Recomputing which pairs clear the bounds test, from the two tokenizers and the
panel alone, reproduces the deposited scored/skipped partition exactly for both
models. The detector loads no weights.

**Consequence for the reported result.** NT v2's and DNABERT-2's Rung 3 values
are computed on 2.7% and 4.8% of their pairs, read from tokens describing a
different part of the molecule, over 4 and 8 families. Three confirmatory
hypotheses use them. H1 asks whether at least one model reaches mean PS > 0 with
enough families exceeding their own null; restoring coverage can only add
candidates. H2 compares RNA-pretrained against DNA-pretrained mean PS, and two of
the five DNA-pretrained models are these two, so the correction moves the
comparison the manuscript's claim rests on and can move it either way. H3's
per-pair precision is a mean over families that these two models barely populate.

The command-line entry point refuses both models unless `--allow-non-character`
is passed (`NON_CHARACTER_TOKENIZERS`, `scripts/phase6_compensatory_mutation.py:52`).
`scripts/modal_repaired_panel.py` calls `run_phase6` directly and does not reach
that guard, which is how the deposited values were produced.

**Fix.** Map each nucleotide to the row holding it, from the tokenizer's own
spans (`scripts/token_spans.py`), for adapters whose `token_resolution` is a
subword scheme. Written and under test; the values change only when the pass
that fixes D10 to D12 runs.

---

## D14. NT v2's last content token is discarded in every family

**Detector:** `scripts/audit_special_token_slices.py`

`NTAdapter.get_all_layer_embeddings` sliced `hs[0, 1:-1, :]` under a comment
reading "Strip CLS/EOS". NT v2's tokenizer emits `<cls>` and nothing after the
sequence, in all 47 families, so the trailing element of that slice is a content
token rather than a special one:

| | leading specials | trailing specials | `1:-1` |
|---|---|---|---|
| NT v2 | 1 in 47 of 47 | 0 in 47 of 47 | drops a content token |
| DNABERT-2 | 1 in 47 of 47 | 1 in 47 of 47 | correct |

Across the panel that removes 92 nucleotides' worth of embedding rows, in every
family. The four multimolecule RNA tokenizers cannot be loaded in either stack
the detector runs under; it names them as unmeasured rather than passing them.

**Consequence for the reported result.** Every NT v2 number in the paper reads a
hidden state one row short. Rungs 1 and 2 place positions into that array through
the D11 expansion, whose tail clamping already touches the same region. Rung 3
loses a further row on top of D13. The attention rung takes the same slice, and
NT v2's contact correlation of 0.242 -- the highest in the panel, and the value
H11 is decided on -- is computed on it.

**Fix.** Every adapter that slices now takes its bounds from its own
tokenization (`content_bounds`, `scripts/token_spans.py`) instead of a constant.
Correctness stops depending on knowing each wrapper's convention, and an
unexpected layout raises instead of silently misaligning. RNA-FM builds its ids
in this repository and keeps its literal slice, with the reason written at the
call site.

**A canary weakens.** The four multimolecule models' slices move from hardcoded
to derived in the same change. Their results are bit-identical across re-runs
only if their tokenizers bracket symmetrically, which could not be checked on the
machine available. NT v2's results will move; a move in the other four is the
measurement, not a regression.

---

## D15. The Rung 3 exceedance test has size 1 - 0.95^L, not 0.05

**Detector:** `scripts/audit_null_selection_bias.py`

`best_ps` is the maximum of mean PS over layers
(`phase6_compensatory_mutation.py:243`). `null_95th_primary` is the 95th
percentile of the derangement null evaluated at the one layer that maximum
selected (lines 341, 351). The observed statistic maximizes over L layers; the
threshold it is compared against is calibrated for one.

Were the layers independent that comparison would reject, under the null, with
probability 1 - 0.95^L. Adjacent transformer layers are strongly correlated, so
the formula is an upper bound and the effective number of independent tests is
below L. The gap grows with L, which is why RiNALMo at L = 34 is 23 points under
its bound while the L = 13 models sit within 2 points of theirs:

| control | L | 1 - 0.95^L | primary | conservative |
|---|---|---|---|---|
| SpliceBERT | 7 | 0.30 | 17/35 = 0.49 | 4/35 = 0.11 |
| UTR-LM | 7 | 0.30 | 10/35 = 0.29 | 7/35 = 0.20 |
| ERNIE-RNA | 13 | 0.49 | 18/35 = 0.51 | 2/35 = 0.06 |
| RNA-FM | 13 | 0.49 | 18/35 = 0.51 | 1/35 = 0.03 |
| RiNALMo | 34 | 0.83 | 21/35 = 0.60 | 0/35 = 0.00 |

Of 47 analyzed families, 38 pass the registered Rung 3 filters, 36 of those are
outside the quarantine, and 35 of those contain a stem with at least three
eligible pairs and so have a derangement null at all. The confirmatory N is 36
and the denominator above is 35.

Randomly initialized weights carry no partner specificity, so the primary column
measures the test's size rather than estimating it. The conservative variant
takes the maximum over layers on the null side too (line 348) and needs no
assumption about layer independence. It sits at or below nominal: RiNALMo
untrained at 0 of 35 is below 0.05 rather than at it, because a maximum over
correlated layers has a heavier upper tail than the statistic it is calibrating.

**The registration specifies the primary variant.**
`PREREGISTRATION_PHASE6_V2.md:79` holds the null at the selected layer "to avoid
inflating the null". Holding the null at one layer while the observed statistic
maximizes over L deflates the null relative to the statistic. The same paragraph
requires both counts to be reported when the two nulls disagree on which
families exceed threshold, and they disagree in every model.

**Consequence for the reported result.** H1(a) counts families exceeding the
primary null against a gate of `ceil(4 * 0.05 * 36) = 8`. Against a test whose
measured size is 0.49 to 0.60, a gate of 8 of 36 is 22% and below the noise
rate, so the criterion as registered separates nothing. RiNALMo reaches 31 of 35
under the conservative test and ERNIE-RNA 33 of 35, so H1 passes on a calibrated
test; it did not pass on a calibrated test before.

**Fix.** Report both counts, which the registration already requires, and
designate the conservative count as the one the hypothesis is decided on.
`exceeds_null_conservative` is stored for every family in every deposited run,
so no GPU pass is needed. The designation is an emphasis chosen after seeing
which way the two counts cut, and it is registered as such in `DEVIATIONS.md`
with its reason: the primary variant's size is measured at 0.49 to 0.60 on
models that cannot have learned anything, so it is not a test.

Same defect class as D10. Both are nulls rejecting above nominal, by different
mechanisms -- D10 through run structure the permutation destroys, D15 through a
maximum tested against a single-layer threshold.

---

## D16. The chance rate for H3's three-way comparison is not 1/3

**Detector:** `scripts/audit_null_selection_bias.py`

H3 counts pairs where `d_partner` exceeds `max(d_prev, d_next)` and tests the
fraction against 1/3, "partner is max among three adjacent positions by chance"
(`PREREGISTRATION_PHASE6_V2.md:109`). Randomly initialized weights have no
partner specificity, so their fraction measures that chance rate:

| control | per-pair precision |
|---|---|
| UTR-LM | 0.113 |
| ERNIE-RNA | 0.161 |
| RiNALMo | 0.175 |
| RNA-FM | 0.232 |
| SpliceBERT | 0.275 |

All five sit below 1/3, and they span 0.113 to 0.275. No single constant could
have been correct in the registration; the rate is a property of each model.

The three positions are not exchangeable. The partner `j` lies between `j_prev`
and `j_next`, so whichever neighbor is nearer the mutated position carries
whatever perturbation reaches that far along the sequence, and a maximum of two
neighbors is taken against one partner. That the deficit is systematic is
measured. That sequence distance produces it is not: line 252 computes `d_prev`
and `d_next` separately, but only their maximum enters `pair_details`, and
`pair_details` is dropped before the result is written, so neither reaches disk.

**Consequence for the reported result.** A null of 1/3 against a chance rate
near 0.15 makes H3 conservative. RiNALMo at 0.838 and ERNIE-RNA at 0.873 pass
with room. What changes is the middle of the table: RNA-FM at 0.238 and
SpliceBERT at 0.243 are at or below their own controls, at 0.232 and 0.275, and
SpliceBERT's trained model scores below its untrained one.

**Fix.** The derangement null already computes the quantity. For each deranged
partner it evaluates `d_partner - max(d_prev, d_next)` on a pairing that is
false by construction, so the fraction of deranged pairs where that is positive
is H3's chance rate -- per model, matched to the same weights, sequences and
layer, and available for all ten models rather than the seven with untrained
counterparts. `derangement_null` returns percentiles only and discards the 1000
per-derangement values, so this is not recoverable from the deposit. Storing
that fraction is one float per family; storing `d_prev` and `d_next` alongside
`d_adj` is one float per pair and settles the mechanism. Both go in the next
pass.

---

## D17. RNA-FM ran with a random vector added to every embedding, and without its final layer norm

**Detector:** `scripts/audit_run_to_run_drift.py`, and the loader's own report.

RNA-FM is the only adapter that builds a bare `BertModel` and fills it with
`load_state_dict(strict=False)`. Every other model comes from
`from_pretrained`, which refuses a checkpoint it cannot match. Three of 199
tensors were absent from the checkpoint and kept whatever `BertModel(config)`
drew for them, which is not seeded:

| tensor | read by the forward pass |
|---|---|
| `embeddings.token_type_embeddings.weight` | yes, at every position of every sequence |
| `pooler.dense.weight`, `pooler.dense.bias` | no |

RNA-FM has no token types and `get_all_layer_embeddings` passes no
`token_type_ids`, so `BertModel` defaulted them to zero and added row 0 of that
random matrix to every position. Initialized `normal_(0, 0.02)` at
`hidden_size = 640`, its norm is about 0.5. It was redrawn in every container.

Two further tensors, `encoder.layer_norm.weight` and `.bias`, were in the
checkpoint and matched nothing in the model. RNA-FM is ESM-architecture and
normalizes after the last encoder block. The remap rule requires
`"encoder.layer."` with the trailing dot, so `encoder.layer_norm` fell through
unrenamed, and `BertModel` has no slot for it. The normalization was dropped.

**How it surfaced.** Comparing stored per-layer PS profiles across two runs, every
model moved by exactly zero except the two whose code had changed -- and trained
RNA-FM, which moved 2.26e-04 and switched `best_layer` in 15 of 36 families. Its
untrained control moved by exactly zero. The control re-initializes every
parameter under `RANDOM_INIT_SEED`, which overwrites the unseeded draw; the
trained model does not. That asymmetry is what named the mechanism.

**Consequence for the reported result.** Every RNA-FM number this project has
deposited, in every rung, was computed with a random per-container offset on all
embeddings and without the model's final normalization. The numbers are not
reproducible and are not attributable to RNA-FM. H10 is a trained-versus-
untrained sign test on RNA-FM and reads both arms; RNA-FM also contributes a row
to H2's five RNA-pretrained models and to the H3 table.

**This defect runs against the reported conclusion.** The manuscript reports
RNA-FM as not resolving partners. A model carrying a random embedding offset and
missing a layer norm is a confound for that negative, not evidence for it. Every
other defect registered on 2026-08-24 ran the other way.

**Fix.** `token_type_embeddings` zeroed, which is what having no token types
means. `encoder.layer_norm` held separately and applied to the final hidden
state, where ESM applies it; the checkpoint carries its weights but not its
epsilon, so it takes the one the rest of this port uses. The `0.9` loaded
threshold was written to catch a remapping that failed wholesale and never fired
on three tensors out of 199; any missing key outside the pooler, and any
unexpected key at all, now raises.

---

## D18. DNABERT-2 disagrees with itself by more than its own signal

**Detector:** `scripts/report_resolution_and_chance.py`

`noop_floor` runs one unmutated sequence through a model twice and computes the
same statistic as `best_ps`, whose true value is then exactly zero. Whatever
comes back is the floor the measurement resolves.

| model | floor, median | floor, max | mean PS |
|---|---|---|---|
| DNABERT-2 | 1.39e-02 | 6.21e-02 | -0.0204 |
| DNABERT-2 untrained | 3.92e-06 | 3.70e-03 | |
| every other model | ~2e-17 | ~7e-17 | |

The floor is larger than the quantity being measured. Every other model in the
panel, trained and untrained, floors at the cosine's own float64 rounding, so
this is DNABERT-2 rather than the pipeline. Its Triton flash attention is the
candidate; `torch.use_deterministic_algorithms` is called with `warn_only=True`,
which warns rather than forces, so nothing stopped it.

`best_ps` is a maximum over layers. When the floor swamps the signal at every
layer, that selection returns the largest noise draw rather than the largest
effect, so the statistic is biased upward by an amount nobody can bound.

**Consequence for the reported result.** No DNABERT-2 Rung 3 number is a
measurement. H2 reads its mean PS as one of five DNA-pretrained values and H20 is
a DNABERT-2 attention contrast. Whether Rungs 1 and 2 are affected is not known:
the floor was measured on the Rung 3 statistic, and the same control has not been
run on the ratio.

**Fix.** Force `torch.use_deterministic_algorithms(True)` rather than warn. If it
raises on the flash-attention kernel, that is the finding -- the model cannot be
run deterministically in this configuration -- and it is reported rather than
routed around.

**Forcing it is not sufficient.** Determinism is now forced and
`CUBLAS_WORKSPACE_CONFIG=":4096:8"` is set in the image, and DNABERT-2 still
disagrees with itself. Two HTT runs at commits `0c9e668` and `e9a0ed4`, both
stamped `determinism: forced`, both on torch 2.4.0, over the same five
fragments:

| model | mean stem-loop ratio | probing accuracy |
|---|---|---|
| DNABERT-2 | 0.984470 -> 0.998167 (1.4%) | 0.617814 -> 0.611894 (0.96%) |
| Caduceus | 1.2314275527 -> 1.2314275321 (1.7e-08) | identical |
| the other seven | identical | identical |

Caduceus's difference is at float64 rounding. DNABERT-2's is a percent, on a
quantity whose whole range across the panel is a few percent from unity, and it
moves the probing accuracy as well as the ratio, so it is not confined to the
Rung 3 statistic where the floor was measured. Eight models reproducing exactly
under the same pins rules out the run-to-run environment as the cause.

Its bundled `flash_attn_triton.py` is the remaining candidate and has not been
tested in isolation. Until it is, DNABERT-2's numbers carry a reproducibility
caveat wherever they appear, including Rungs 1 and 2.

---

## D19. An intervention applied before the first forward pass is silently discarded

**Detector:** `scripts/modal_inspect_multimol_buffers.py`

`transformers` v5 leaves a non-persistent buffer on the meta device after
`from_pretrained`, so it materializes as whatever memory held. Every
multimolecule model in this panel repairs its own buffers on first use, and none
of them is left damaged:

| model | non-persistent buffers | garbage before a forward | still garbage after |
|---|---|---|---|
| RiNALMo | 34 | 33 | 0 |
| UTR-LM | 7 | 4 | 0 |
| SpliceBERT | 1 | 0 | 0 |
| ERNIE-RNA | 3 | 2 | 0 |

No deposited result is contaminated by this. What it does break is any
intervention written between the load and the first forward pass, because the
repair overwrites it.

`scripts/modal_ernierna_ablation.py` is that case. It zeros `pairwise_bias_map`
and then calls `run_mutation_sensitivity`, whose first forward pass reaches
`if not self._inited:` in `ErnieRnaEmbeddings` and rebuilds the buffer with the
six canonical pairing values. The deposited run at
`results/ernierna_ablation/ernierna_no_attn_bias_20260717_104015` therefore
ablated nothing, and its numbers measure the unmodified model.

**Consequence for the reported result.** Nothing in the manuscript reads that
run. It was cited during this audit as evidence that the pairing bias barely
affects Rungs 1 and 2, and that reading is withdrawn: the comparison was between
one model and itself.

There is an accidental measurement in it. The intact model at 47 families under
the current code gives a mean best ratio of 1.2758; the same model at 52 families
under July's code and float32 metric gives 1.2348. The 0.041 gap is what the
panel change, the token-span repair and the arithmetic move a Rung 1 number when
the model is held fixed.

**Fix.** An ablation on a multimolecule model sets `model._inited = True` after
mutating, or patches the builder. `scripts/modal_ernierna_ablation.py` is left as
deposited, since its output is on record; the registered replacement carries the
flag.

---

## D20. A randomly initialized control is not a null model

**Detector:** `scripts/audit_null_selection_bias.py`,
`scripts/modal_inspect_multimol_buffers.py`

`_randomize` destroys every learned weight and keeps the architecture, which is
what makes the control a control. It iterates `named_parameters()`. A buffer is
not a parameter, so anything a model holds in one survives untouched, and the
"untrained" model keeps whatever that buffer encodes.

Measured against each model's own derangement chance rate, the controls are not
at chance:

| randomly initialized control | chance | observed | excess |
|---|---|---|---|
| ERNIE-RNA | 0.241 | 0.476 | +0.236 |
| RiNALMo | 0.232 | 0.327 | +0.095 |
| RNA-FM | 0.230 | 0.314 | +0.084 |
| SpliceBERT | 0.194 | 0.204 | +0.010 |
| UTR-LM | 0.127 | 0.135 | +0.008 |
| NT v2 | 0.000 | 0.000 | 0.000 |
| DNABERT-2 | 0.000 | 0.000 | 0.000 |

ERNIE-RNA's cause is identified: `pairwise_bias_map` is a non-persistent buffer
holding a Watson-Crick and wobble table -- A-U 2.0, C-G 3.0, G-U 0.8 -- which
randomization cannot reach. A model told which bases pair, with random weights on
top, is not a model that knows nothing about pairing.

**RNA-FM's excess is not carried by a buffer.** It holds two non-persistent
buffers, `position_ids` and `token_type_ids`, both integer, so there is no
non-parameter float state for an effect to live in.

**RiNALMo's reseeding is inconclusive by design.** Redrawing its 33 rotary
`inv_freq` tables drops its excess from +0.105 to +0.029, but `inv_freq` is a
deterministic geometric ladder, `1/base^(2i/d)`, rather than a table of learned
values. Replacing it with noise corrupts a well-formed position code instead of
removing information, and a drop is equally consistent with the buffers carrying
the effect and with a malformed position code degrading everything downstream.
The same manipulation on ERNIE-RNA drops +0.249 to +0.047, reproducing the
zeroing ablation through a different intervention, because there the buffer is a
lookup table and randomizing it does remove information.

**A candidate mechanism for both, and it is one this register already carries.**
Per-pair precision is read at the layer that maximizes mean PS, so a model with
more layers has more draws to select from. Across randomly initialized
character-level controls, with ERNIE-RNA set aside because its buffer is
identified, the excess is monotone in layer count:

| control | layers | excess over its own chance rate |
|---|---|---|
| RiNALMo | 34 | +0.105 |
| RNA-FM | 13 | +0.093 |
| SpliceBERT | 7 | +0.016 |
| UTR-LM | 7 | +0.011 |

Neither declared depth nor parameter count predicts it (Spearman +0.32 and
+0.27, both far from significance at n = 7). Layer count does, and it is the
asymmetry D15 identified in the exceedance count appearing in the precision
statistic: a maximum over layers compared against something not selected the
same way.

**The cause, and the candidate was stated at the wrong grain.** The layer is
chosen per family: 35 independent argmax draws, one for each family in the Rung 3
analysis set, not one draw for the model. Choosing once for the whole panel
instead leaves the trained models where they were and removes almost all of both
residual excesses (`scripts/layer_selection_audit.py`, which raises if a stored
`h3_precision` and its recomputation from `precision_per_layer` disagree):

| model | per-family layer | one panel layer | held-out |
|---|---|---|---|
| ERNIE-RNA trained | 0.888 | 0.887 | 0.885 [0.814, 0.955] |
| RiNALMo trained | 0.872 | 0.881 | 0.882 [0.822, 0.941] |
| ERNIE-RNA untrained | 0.490 | 0.479 | 0.478 [0.399, 0.554] |
| RiNALMo untrained | 0.337 | 0.158 | 0.158 [0.106, 0.208] |
| RNA-FM untrained | 0.323 | 0.143 | 0.143 [0.091, 0.193] |
| SpliceBERT untrained | 0.210 | 0.172 | 0.171 [0.121, 0.222] |
| UTR-LM untrained | 0.139 | 0.132 | 0.131 [0.081, 0.180] |

The held-out column chooses the layer on half the families and reports precision
on the other half, over 2,000 splits, and is selection-free by construction. It
reproduces the panel-layer column to three decimals in every row, so the two
disagree with the per-family column for the same reason.

Against each model's own chance rate the held-out excesses are +0.761
(ERNIE-RNA), +0.780 (RiNALMo), +0.238 (ERNIE-RNA untrained), **-0.074** (RiNALMo
untrained), **-0.087** (RNA-FM untrained), -0.023 (SpliceBERT untrained) and
+0.004 (UTR-LM untrained). Four of the five controls whose cause was unidentified
sit at or below their own chance rates once the layer is not chosen per family.
The one that does not is ERNIE-RNA, whose buffer is identified above and whose
ablation removes it.

The monotonicity in layer count above is the same fact seen from outside: a
34-layer control draws the per-family maximum from 34 profiles and a 7-layer
control from 7.

**The whole panel, and it is the whole of every mid-table excess.** Six further
models were re-run phase-6-only to store `precision_per_layer`; each reproduced
its prior numbers, and UTR-LM's six families differing at relative 1e-6 to 1e-4
are the only drift, explained by its earlier run predating forced determinism.
Excess over each model's own chance rate, per-family against held-out:

| model | layers | per-family | held-out |
|---|---|---|---|
| RiNALMo | 34 | +0.770 | +0.780 |
| ERNIE-RNA | 13 | +0.764 | +0.761 |
| ERNIE-RNA untrained | 13 | +0.249 | +0.237 |
| Caduceus | 17 | +0.122 | **-0.216** |
| RNA-FM | 13 | +0.109 | **-0.135** |
| RiNALMo untrained | 34 | +0.105 | -0.075 |
| RNA-FM untrained | 13 | +0.093 | -0.086 |
| Evo | 32 | +0.040 | -0.000 |
| UTR-LM | 7 | +0.016 | -0.039 |
| SpliceBERT untrained | 7 | +0.016 | -0.022 |
| UTR-LM untrained | 7 | +0.011 | +0.004 |
| SpliceBERT | 7 | +0.009 | -0.034 |
| HyenaDNA | 4 | -0.010 | -0.027 |

Every model except RiNALMo, ERNIE-RNA and ERNIE-RNA's bias-carrying control
sits at or below its own chance rate once the layer is not chosen per family.
Caduceus goes from 0.340 to 0.002, which is what a 17-layer noise profile looks
like when the maximum is taken 35 times and when it is not. The four models at
precision 0.000 or 0.007 cannot move and were not re-run.

The manuscript abstract reported Caduceus and RNA-FM as exceeding chance. Both
are corrected in `paper/rna-structure-audit_v16.tex`, which reports per-family,
panel-layer and held-out precision side by side.

**Consequence for the reported result.** Every trained-versus-untrained
comparison reads the untrained arm as a floor. H10, H14, H17 and H19 are sign
tests defined that way, and H11 and H20 are attention contrasts. Where the
untrained arm carries an architectural prior, the comparison measures learning
plus prior against prior, not learning against nothing, and the gap understates
or overstates depending on which side the prior helps.

It also removes the baseline D15 used. That entry measured the primary
exceedance test's size as 0.49 to 0.60 by treating the untrained controls as
null. The selection asymmetry D15 identified -- a maximum over L layers tested
against a single-layer threshold -- is a fact about the code and stands
independently. The number attached to it does not.

**Fix.** Not settled. `scripts/modal_inspect_multimol_buffers.py` lists what each
model holds outside its parameters, which is the first step. An
architecture-neutral null would have to remove the prior as well as the weights,
which is what the ERNIE-RNA ablation registered in
`preregistration/PREREGISTRATION_ERNIERNA_BIAS_RUNG3.md` tests for one model. For
the rest, the derangement null is the only baseline in this project that is
matched to a model's own weights, sequences and layer, and it should be preferred
wherever a chance rate is needed.

---

## D21. The HTT embedding-distance column measures sequence length

**Detector:** `scripts/generate_arm_tables.py`, which reports the Spearman
correlation at the embedding layer beside the one at the final layer.

The HTT case study reports cosine distance of each fragment's mean-pooled
embedding from the CAG17 fragment, and reads it as whether a model separates
normal from pathogenic alleles. The manuscript's sentence is that RNA-FM
"produces embedding distances from wild-type below 1e-4 for all repeat counts,
failing to distinguish normal from pathogenic alleles", and that RiNALMo "shows
the strongest separation (0.068 at CAG60)".

Recomputed in float64 at every layer, distance rises monotonically with repeat
count for all eight models with more than one layer, at Spearman +1.00 -- and
it already does so at the **embedding layer**, before any transformer block has
run. The five fragments have lengths 235, 247, 292, 304 and 364 and differ in
composition accordingly, so a mean-pooled token embedding separates them by
construction. A monotone distance profile is evidence about the fragments, not
about the representation.

Two consequences. The column cannot rank models by whether they encode the
expansion, so "RiNALMo shows the strongest separation" is a statement about
where the measure was read. And selecting the layer by this correlation, which
an earlier version of `_htt_distances` did, chooses among ties: every layer
scores +1.00.

**What survives.** The final-layer reading is still worth reporting, because one
model departs from the pattern there: RNA-FM collapses to 2.2e-07 in float64,
which is a real degeneracy of its final representation rather than a float32
artifact. The deposited column, computed on 2026-07-13 by a script no longer in
the repository, is reproduced exactly by the final-layer float64 recomputation
for ERNIE-RNA, RiNALMo, SpliceBERT and UTR-LM, which recovers that script's
specification.

**Fix.** Report the final layer and the embedding layer side by side, and state
that the embedding layer already reaches +1.00. Applied in
`paper/generated/htt_table.tex`.

---

## D22. Three parameter counts in the models table are not the checkpoints that ran

**Detector:** the HuggingFace API's `safetensors.total` for each model id in
`multi_model_audit.py`, which reports the parameter count of the exact
checkpoint loaded.

The models table reports a parameter count per model, and the manuscript quotes
those counts in three places to argue that scale does not order the panel.
Checked against the checkpoints the adapters actually load:

| model | checkpoint | table | actual |
|---|---|---|---|
| RNA-FM | `multimolecule/rnafm` | 99M | 99,523,469 |
| RiNALMo | `multimolecule/rinalmo-giga` | 650M | 650,917,159 |
| ERNIE-RNA | `multimolecule/ernierna` | 86M | 86,261,884 |
| SpliceBERT | `multimolecule/splicebert` | 19M | 19,719,708 |
| UTR-LM | `multimolecule/utrlm-te_el` | 1.2M | 1,210,494 |
| NT v2 | `InstaDeepAI/nucleotide-transformer-v2-50m-multi-species` | 56M | 55,905,164 |
| **HyenaDNA** | `LongSafari/hyenadna-tiny-1k-seqlen-hf` | **5.4M** | **450,712** |
| **Caduceus** | `kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16` | **14M** | **7,725,312** |
| Evo | `togethercomputer/evo-1-8k-base` | 7B | 6,452,781,248 |
| DNABERT-2 | `zhihan1996/DNABERT-2-117M` | 117M | no safetensors metadata |

HyenaDNA is out by 12x: the checkpoint loaded is the smallest release, and 5.4M
matches none of the published HyenaDNA sizes. Caduceus is out by 1.8x. Neither
checkpoint carries a parameter count in its name, so neither number was an
author designation -- both were wrong.

Evo and DNABERT-2 are different: 7B and 117M are what their authors call them,
and the names are how readers recognize the models, so those are kept. Evo's
checkpoint holds 6.45B, which is the usual gap between a model's name and its
released weights. DNABERT-2 publishes only a `pytorch_model.bin` with no
metadata, so its count is the authors' designation and is unverified here.

**Consequence for the reported result.** The direction of the argument is
unchanged and its force increases: the sentence "Caduceus at 14M outscores three
of the five RNA-pretrained models" is stronger at 7.7M, and the panel's range is
0.45M to 7B rather than 1.2M to 7B.

**Fix.** Applied in `scripts/generate_results_tables.py`, the v16 models table
and the workshop paper. A generator that asserted its parameter counts against
the API at build time would have caught this; it does not, and the counts remain
hand-entered.

---

## D23. The token-span repair treats an IUPAC ambiguity code as a special token

**Detector:** the multi-sequence arm, which is the first stage to score raw
seed-alignment sequences rather than curated families.

`scripts/token_spans.py` finds which rows of a hidden state carry sequence by
asking which token strings are made of nucleotides:

    NUCLEOTIDES = set("ACGTUN")
    is_content = [bool(p) and set(p) <= NUCLEOTIDES for p in pieces]

A token is content when its characters are a subset of that set, and anything
else is a special token. Rfam seed alignments carry the rest of IUPAC -- W, Y
and H appear in this panel's replicates -- so a position holding one of them is
classified as a special token, and because it sits between content tokens the
contiguity guard raises:

    TokenizationMismatch: special tokens ['W'] sit between content tokens

Six of ten models died on it at the same sequence, `Corona_5UTR_seq0`, the 21st
of 233. The four that survived do not call `content_bounds`: RNA-FM slices
`hidden[0, 1:-1, :]` directly and HyenaDNA, Caduceus and Evo are character-level
with no special tokens. The models running the repaired code are exactly the
models that failed, which is the shape of a repair introducing its own defect.

Both spellings occur. multimolecule tokenizers emit the letter (`'W'`);
DNABERT-2 emits `[UNK]` and NT v2 `<unk>`.

**Why it is not repaired in place.** For a character tokenizer an interior
unknown covers one nucleotide and the row could be kept. For NT v2's 6-mer and
DNABERT-2's byte-pair vocabularies a single unknown token covers a whole k-mer,
and which nucleotides that row holds cannot be read back from the tokenization.
Widening `NUCLEOTIDES` would fix four models and silently invent a span mapping
for two.

**Fix.** Sequences carrying characters outside `ACGTUN` are excluded from the
multi-sequence arm for every model, and the excluded names travel in the result
file as `excluded_for_ambiguity`. That is 6 of 233 sequences across 3 of 47
families -- Corona_5UTR (3 of 5), preQ1_riboswitch (2 of 5), THF_riboswitch
(1 of 5) -- and no family loses all its replicates. Excluding uniformly matters
more than excluding few: the arm compares within-family variance across models,
and a set that differs by model is not a comparison.

**Scope.** No curated panel family contains a character outside ACGU, so Rungs
1--3, the transversion, synthetic and HTT arms are untouched. Any future use of
this pipeline on raw alignment sequences meets the same wall.

---

## Checked, clean

- Every numeric literal in the manuscript against a stored source. The existing
  checks (`verify_paper_rung12_figures.py`, `verify_paper_phase6_figures.py`,
  `audit_table12_sources.py`, `audit_table5_aggregation.py`,
  `audit_attention_rho.py`) cover the tables and the Phase 6 figures. Prose
  numbers outside those tables have no systematic check, and D3 was found by
  hand.
- Completeness of `results/`: which model x family x rung cells exist and which
  are absent.
- Whether the class composition at `paper_v12.tex:215-218` (7 tRNAs, 4 rRNAs,
  6 ribozymes, ...) still describes the deposited panel.
