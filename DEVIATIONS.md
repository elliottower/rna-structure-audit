# Analysis Deviations and Corrections

Preregistrations governing this work, with sha256 of each frozen document:

| document | sha256 (first 16) |
|---|---|
| `PREREGISTRATION.md` (Phase 1) | `8b3b63b8407a4b2f` |
| `docs/PREREGISTRATION_EXPANDED_RFAM.md` (Phase 2) | `91eee9994c73ec8b` |
| `PREREGISTRATION_PHASE6_COMPENSATORY_MUTATION.md` (Phase 6 v1) | `17aac5a98732ba38` |
| `PREREGISTRATION_PHASE6_V2.md` | `d89be39cddf6f047` |
| `PREREGISTRATION_V1.md` | `94ac3d4ea41acaaa` |
| `preregistration/PREREGISTRATION_PHASE6_UNTRAINED_RINALMO.md` | `3966fd1156117b55` |

---

---

## 2026-08-24 — RNA-FM's port added a random vector to every embedding and dropped its final layer norm

**Registrations:** `docs/PREREGISTRATION_EXPANDED_RFAM.md` H10, the RNA-FM
trained-versus-untrained sign test. RNA-FM also carries a row in H2 and H3
(`PREREGISTRATION_PHASE6_V2.md`).

**Written before RNA-FM was recomputed.**

Three of 199 tensors were absent from the checkpoint and kept an unseeded random
initialization. One of them, `token_type_embeddings`, is read at every position
of every sequence, because RNA-FM has no token types and the forward pass
supplies none, so `BertModel` defaulted them to zero and added row 0 -- a vector
of norm about 0.5 -- to every embedding. It was redrawn in every container.
Separately, the checkpoint's `encoder.layer_norm` matched nothing in the model
and RNA-FM's final normalization was dropped. Details and the detector are in
`docs/OPEN_DEFECTS.md`, D17.

**Which registered decisions move.** Every RNA-FM number in every rung. H10 reads
both arms of a trained-versus-untrained comparison and the trained arm was not
reproducible. RNA-FM's rows in the H2 ranking and the H3 table move with it.

**The direction of this correction.** It runs against what the work reports.
RNA-FM is reported as not resolving partners, and a model with a random
embedding offset and no final layer norm is a confound for that negative rather
than evidence for it. Every other defect registered today ran toward the
reported conclusions; this one does not, and it is the reason the asymmetry noted
in the entry above is a property of where the errors were rather than of how
they were looked for.

**Cost.** RNA-FM's three cells -- trained, untrained control, transversion
control -- are recomputed against the repaired port. They carry a later commit
than the rest of the panel and are stamped with it.

---

## 2026-08-24 — Two null calibrations: the Rung 3 exceedance test has size 0.5, and H3's chance rate is not 1/3

**Registrations:** `PREREGISTRATION_PHASE6_V2.md`, H1(a) at line 96, H3 at line
110, and the null construction at line 79.

**Written before the affected numbers were recomputed.**

**D15.** `best_ps` is a maximum of mean PS over L layers.
`null_95th_primary` is the 95th percentile of the derangement null evaluated at
the single layer that maximum selected. The registration specifies this at line
79, reasoning that holding the null at the selected layer avoids inflating it;
the reasoning is inverted, because deflating the null relative to the statistic
inflates the rejection rate. Measured on randomly initialized weights, which
carry no partner specificity, the test's size is 0.49 to 0.60 rather than 0.05.
The conservative variant, which takes the maximum over layers on the null side
as well, sits at or below nominal on the same controls: 0.06, 0.03, 0.11, 0.20
and 0.00.

**D16.** H3 tests the fraction of pairs where the partner's perturbation exceeds
the larger of its two stem neighbors' against 1/3. On the same controls that
fraction is 0.113 to 0.275, never 1/3, and the spread means no single constant
could have been the right null. The partner sits between the two neighbors, so a
maximum of two is taken against one and the comparison is not exchangeable.

**Which registered decisions move.** H1(a) is decided on the conservative count,
which the registration already requires to be reported wherever the two nulls
disagree; they disagree in every model. The gate stands at
`ceil(4 * 0.05 * 36) = 8`. H3 is reported against each model's own randomly
initialized control where one exists, alongside the registered test against 1/3.
Designating the conservative count as the deciding one is an emphasis chosen
after both counts were seen, and is recorded here for that reason: the primary
variant's size is measured at 0.49 to 0.60 on models that cannot have learned
anything, so it is not a test.

**Cost.** `exceeds_null_conservative` is stored for every family in every
deposited run, so H1(a) needs no recomputation. H3's per-model chance rate is
computable from the derangement null, which already evaluates the same statistic
on pairings false by construction, but `derangement_null` returns percentiles
and discards the per-derangement values. That fraction, and `d_prev` and
`d_next` separately from their maximum, are added to the next pass.

**The direction of these corrections.** Every defect registered today runs
toward the conclusions this work reports: D13 restored pairs that could only
raise the DNA-pretrained side of H2, D15 removes a criterion that admitted noise
at half the families, and D16 shows the H3 threshold was too demanding rather
than too lax. That asymmetry is what one expects when the errors are in the
tests rather than in the data — a mis-sized test is found by running it on
controls that cannot have the effect, and controls only ever fail in one
direction. The order of operations is the only guarantee offered against
motivated correction: every defect above was registered here before the affected
numbers were recomputed, and the commits carry the sequence.

---

## 2026-08-24 — Two further defects: Rung 3 scored 2.7% of NT v2's pairs, and NT v2 lost an embedding row in every family

**Registrations:** `PREREGISTRATION_PHASE6_V2.md` (H1, H2, H3),
`docs/PREREGISTRATION_EXPANDED_RFAM.md` (H11).

**Written before the re-run.** Both defects were found while preparing the pass
that fixes D10 to D12, and no model has been run against either fix. The
measurements, the registered decisions each defect touches and the direction each
correction can move them are fixed here in advance of any value computed under
them. `scripts/audit_rung3_token_alignment.py` and
`scripts/audit_special_token_slices.py` reproduce every number below and load no
model weights.

**D13. Rung 3 reads the hidden state at the nucleotide index for every model.**
`compute_delta_profiles` reads `emb[k + offset]` for nucleotide *k* with
`offset = 0` throughout. Eight models emit one token per nucleotide and the read
is right. NT v2 emits about one row per six nucleotides, so the partner position
of a stem pair is past the end of the array; the pair is dropped by the bounds
test and a family that loses all of its pairs is recorded as `"no valid PS
values"`, the same field a registered filter writes. Of the 38 families clearing
the registered Rung 3 filters, NT v2 keeps 4 and DNABERT-2 keeps 8. NT v2's
reported Rung 3 rests on 16 of 583 eligible pairs and DNABERT-2's on 28, each
read from a token whose nearer edge is a median of 110 and 124 nucleotides from
the position it stands for. Recomputing the bounds test from the two tokenizers and
the panel reproduces the deposited scored/skipped partition exactly. The
corrected read maps each nucleotide to its row through the tokenizer's spans.

**D14. NT v2's adapter stripped a content token as though it were a separator.**
NT v2's tokenizer emits `<cls>` and nothing after the sequence, in all 47
families, and `NTAdapter.get_all_layer_embeddings` sliced `hs[0, 1:-1, :]`. The
last content token of every family is discarded, 92 nucleotides across the panel,
in all rungs including attention. DNABERT-2 brackets symmetrically and is
unaffected. Every adapter now derives its slice from its own tokenization.

**Which registered decisions move.** H1 (at least one model reaching mean PS > 0
with enough families exceeding their own derangement null), H2 (RNA-pretrained
mean PS above DNA-pretrained, rank-biserial > 0.5) and H3 (per-pair precision
above 1/3) are all computed over models that include NT v2 and DNABERT-2, and
D13 restores those two models from 4 and 8 families to their full qualifying set.
H11 (NT v2 attention correlation above 0.15 trained and below 0.05 untrained;
observed 0.242 and 0.246) reads a hidden state D14 shortens. No other registered
decision reads either model's Rung 3.

**The direction of these corrections.** D13 runs against what this work reports.
Two of the five DNA-pretrained models in H2 are the two whose coverage was near
zero, and the registration already anticipated that multi-nucleotide tokens would
dilute their PS -- the defect is not dilution but near-total loss, and restoring
the pairs can only raise the DNA-pretrained side of a comparison the manuscript
argues goes the other way. H1 becomes easier to pass, since restoring coverage
adds candidates and removes none. D14's direction is not predictable: an extra
row changes NT v2's per-layer distances in both arms of every ratio.

**Cost.** Neither defect needs a pass of its own. Both live in the embedding
read, so both are fixed in the single pass over the 17 runs that D10 to D12
already require.

---

## 2026-08-24 — Three defects in the Rungs 1-2 pipeline; the affected numbers are held for one re-run

**Registrations:** `docs/PREREGISTRATION_EXPANDED_RFAM.md` (H6b, H7, H11),
`docs/PREREG_PHASE4_EXPANDED_MODELS.md` (H20, H21).

**Written before the re-run.** No model has been run against a corrected null or
a corrected token assignment. The defects, the fixes, which registered decisions
each one touches and which way each correction can move them are fixed here in
advance of any value computed under them. `scripts/audit_untrained_false_positives.py`
and `scripts/audit_token_alignment.py` reproduce every number below and load no
model weights.

**D10. The Rung 1 null rejects at 9.1% where it is set at 5%.** Seven models were
run with randomly initialized weights, so the rate at which those runs exceed
their own composition-preserving null is a false-positive rate. Pooled: 30 of 329
family-runs, 9.1%, exact binomial one-sided p = 1.3e-3. The null permutes stem
and loop labels within nucleotide strata, which scatters them, while the
annotation puts them in contiguous runs; a model whose embedding pools
neighboring positions separates the real assignment from the permuted one for
that reason alone. NT v2, whose token spans six nucleotides, sits at 10 of 47.
The five character-level models sit at 14 of 235, one-sided 95% lower bound
0.036, which is consistent with 5% and does not establish it. The corrected null
preserves run structure, by circular rotation or by a permutation restricted to
preserve run lengths.

**D11. DNABERT-2's per-position embeddings are placed by a guess.** `phases_1_to_5.py:61-76`
maps nucleotide *i* to token `i * n_tokens / n_nucleotides`, which assumes
byte-pair tokens are all the same width. They run 1 to 8 nucleotides, and a
complement substitution re-segments the sequence, so 45.9% of mutation trials
subtract a wild-type row from a mutant row describing a different stretch of the
molecule. The corrected assignment uses the tokenizer's own character offsets.
NT v2's `i // 6` is exact inside its tiling, misplaces 1.2% of nucleotides at the
tail, and never mispairs the two arms.

**D12. The attention contact map is built on the same assumption.**
`_aggregate_contacts_to_tokens` gives token *t* the window `[6t, 6t+6)` for both
non-character models. For DNABERT-2, 1,207 of 1,361 token rows share no
nucleotide with the token they stand for and mean window overlap is 0.074.

**Which registered decisions move.** H6b (RNA-FM exceeds the nucleotide null in
at least 8 families; observed 0 of 47) and H7 (at least one family survives the
dinucleotide null; observed 30 for ERNIE-RNA) read exceedance counts and are
decided against the D10 null. H20 (DNABERT-2 trained minus untrained attention
below 0.05; observed 0.009) and H21 (DNABERT-2 trained rho below 0.20; observed
0.199) are decided on the D12 contact map. H11 (NT v2 attention above 0.15
trained and below 0.05 untrained; observed 0.242 and 0.246) reads NT v2's, whose
mean window overlap is 0.936. No registered hypothesis is decided on DNABERT-2's
Rung 1 or Rung 2 values; those appear in the tables only.

**The direction of these corrections.** A null that rejects too often inflates
exceedance counts, so correcting D10 weakens the discrimination claims and
strengthens the negative ones, which is the direction that favors what this work
reports. D12 runs the other way: both hypotheses it touches currently read PASS,
and H21 passes by 0.001. Nothing about the order in which the three were found
distinguishes the honest case from the motivated one, and the entry is written
before the re-run for the same reason the entry above it was.

**The pattern across the day's corrections is itself worth recording.** Of the
defects found on 2026-08-24, most move a number in the direction the manuscript
argues for. A run of same-signed corrections is what motivated cleaning produces,
and it is also what a manuscript that overclaimed produces when it is audited.
The record cannot settle which, so it states the run rather than leaving it to be
noticed: every correction is committed with its detector, each detector reads
data rather than model output where that is possible, and each entry precedes the
computation it licenses.

**Cost, and why it is one pass.** Per-position cosine distances are not stored.
The deposited Rung 1 files carry `best_ratio`, `best_layer`, the two null
thresholds, and the stem and loop counts, with nothing per layer or per position,
so neither the null nor the token assignment can be corrected from stored values.
All three defects need the embeddings, so all three are fixed in one GPU pass
over the 17 runs. That pass stores per-position distances, which makes a future
change of null a local recomputation.

**Amendment to the post hoc attention bound.** `scripts/generate_results_tables.py`
carried a bound of 0.025 on the trained-minus-untrained attention correlation,
below which a correlation is read as architectural rather than learned. The bound
is not registered; it was raised to 0.05 on 2026-08-24, after the deltas were
seen, which is the configuration in which a threshold is chosen to fit the data.
Two of the seven deltas it covers are DNABERT-2's, which D12 makes
uninterpretable. The manuscript reports the seven deltas and the widest gap as a
factor of the baseline it moved from, and makes no bounded claim, until the
re-run supplies a contact map the deltas can be computed against.

---

## 2026-08-24 — Unattributed annotations repaired; confirmatory set moves to N = 36 and the H1(a) gate to 8

**Registration:** `PREREGISTRATION_PHASE6_V2.md` (2026-07-13, SHA `c19aa59`).

**Written before the re-run.** No model has been run against any repaired
annotation. This entry fixes the repair rule, the disposition of every affected
record, the resulting N, the resulting gate and the reporting plan in advance of
any result computed under them. `scripts/repair_annotations.py` reproduces every
disposition below and writes `docs/annotation_repair_manifest.json`.

**What was found.** A dot-bracket annotation that belongs to its sequence pairs
almost only Watson-Crick and GU; one copied onto a different sequence pairs
whatever sits at those positions. `scripts/audit_annotation_validity.py` computes
that fraction for all 52 families. The 41 records naming an Rfam seed member
reach at most 18% non-canonical pairs. Of the 11 naming none, nine run 24-77%.
The two groups do not overlap, and the criterion reads no model output, so the
detection is independent of any outcome.

**The rule**, applied in order to all 11 records naming no seed member:

1. Annotation at most 25% non-canonical — **no repair**; it describes its own sequence.
2. Sequence is an exact seed member — **adopt structure**; the sequence is kept byte for byte.
3. No member of the family projects to a structure under 5% non-canonical — **drop**.
4. The name asserts a subdomain — **drop**; every seed member is a whole molecule.
5. The name asserts an organism the seed cannot supply — **defer**; a substitute would leave the name false.
6. Closest clean member differs in length by more than 30% — **drop**.
7. Otherwise **adopt member**, the one closest in length among members at most 5% non-canonical.

Rule 6's threshold sits in an empty band: the candidates produce length
differences of 0, 1, 2, 18, 71, 210 and 269 percent, so every threshold between
18 and 71 partitions them identically. Rules 4 and 5 are the only ones reading
anything but the data, and the claim each name makes is written out in the script.

**Dispositions.**

| family | disposition | length | non-canonical | Rung 3 after |
|---|---|---|---|---|
| SAM_riboswitch | adopt structure | 119 → 119 | 52% → 3% | qualifies |
| HDV_ribozyme | adopt member | 85 → 87 | 58% → 0% | qualifies |
| TPP_riboswitch | adopt member | 104 → 103 | 57% → 4% | qualifies |
| hammerhead_ribozyme | adopt member | 38 → 45 | 44% → 0% | no (12 WC pairs) |
| 5S_rRNA_ecoli | defer | 120 | 70% | — |
| IRES_HCV_domainII | drop | 86 | 55% | — |
| SRP_RNA_helix8 | drop | 90 | 77% | — |
| RNaseP_specificity | drop | 94 | 58% | — |
| U2_snRNA_stem | drop | 83 | 54% | — |
| tRNA_Ala_human | no repair | 75 | 24% | quarantined |
| tRNA_Phe_yeast | no repair | 76 | 0% | quarantined |

`5S_rRNA_ecoli` is deferred rather than substituted because RF00001's seed carries
no *E. coli* member: the stored sequence is 120 nt and shares 92% of its 8-mers
with a seed member, so it is genuine 5S rRNA, and replacing it would leave the
family name asserting an organism the record no longer holds. Repairing it needs
an annotation for the sequence actually stored, which a seed alignment cannot
supply.

**The rule is applied to all eleven, not to the six inside the analysis.** The
three families excluded post hoc for unbalanced brackets were excluded on the
strength of the same defective annotations, so that exclusion was not independent
of the defect; `hammerhead_ribozyme` returns to the Rungs 1-2 panel under a valid
annotation. The two quarantined families are a separate matter: `PREREGISTRATION_PHASE6_V2.md`
lines 13 and 147 quarantine them because pilot values were observed under a weaker
metric before the registration was written, which is foreknowledge and not data
quality. `tRNA_Phe_yeast`'s annotation is 0% non-canonical. Both quarantines stand
on their original grounds.

**Consequence for the registered set.** Applying the frozen eligibility criteria
to the repaired annotations, `HDV_ribozyme`, `SAM_riboswitch` and
`TPP_riboswitch` become eligible — their fabricated structures pair so
non-canonically that they carried 6, 13 and 8 Watson-Crick pairs against a
threshold of 15, and their repaired structures carry 22, 32 and 25. With
`mir_21_precursor`, entering under the correction recorded below, the confirmatory
set moves from N = 32 to N = 36. The Rungs 1-2 panel moves from 49 families to 47:
`5S_rRNA_ecoli`, `IRES_HCV_domainII` and `SRP_RNA_helix8` leave, `hammerhead_ribozyme`
returns.

**No registered rule changed.** The eligibility criteria — at least 15 canonical
WC pairs, stems of at least 3 consecutive WC pairs, the two terminal pairs of each
stem excluded, at least 5 eligible interior pairs remaining — are as frozen, and
`scripts/repair_annotations.py` applies them through the analysis package's own
`_parse_stems` and `_get_eligible_pairs` rather than restating them.
`PREREGISTRATION_PHASE6_V2.md:29` states in advance that "the exact count of
eligible families depends on data quality," and the document nowhere enumerates
the confirmatory families by name.

**The registered decision threshold moves.** H1 condition (a) is the only
confirmatory criterion carrying N: exceedances must reach `ceil(4 × 0.05 × N)`.
That is `ceil(6.4) = 7` at N = 32 and `ceil(7.2) = 8` at N = 36. H1(b), H2 and H3
fix their thresholds independently of N. `scripts/check_registered_n_sensitivity.py`
recomputes the gate at all three counts.

**The direction of the correction.** Every repair moves a family into eligibility
and none out, and the gate rises. A rising gate makes H1 harder to satisfy, which
is the direction favoring the negative reading this work reports elsewhere. That
is the configuration in which motivated data cleaning would appear, and no
argument made afterwards distinguishes it from the honest case. What is offered
instead is the order of operations: the criterion that flagged the defect reads no
model output, the repair rule was written before it was applied and applied
uniformly, and this entry precedes the re-run in the commit history.

**Amendment, same day, before the run: the nulls are reseeded.** Each null is
built inside one family — Rungs 1-2 permute stem and loop labels within
nucleotide strata of that family's own sequence, Rung 3 deranges partner
assignments within that family's own stems — so the distribution a family's
threshold is drawn from does not depend on which other families are present. The
draw did. `phases_1_to_5.py` seeded from `42 + rna_idx * 1000`, the family's
position in the loaded list, and `phase6_compensatory_mutation.py` seeded one
global generator that every family then drew from in panel order. Withdrawing five
records shifts the stream every later family receives, so families whose
annotation never changed would have received different nulls. `scripts/family_seed.py`
now derives each family's generator from its name, which makes the implementation
match the definition: a family's null is reproducible across runs and unaffected
by panel membership. The cost is that every family's null is redrawn once at this
run, so no per-family null in the re-run is comparable digit for digit with the
deposited one.

**What the correction can and cannot move.** H1(a) does not fail under the repair,
and with the nulls redrawn the claim rests on margins rather than on invariance.
`scripts/check_h1a_margin.py` measures, for every scored family, the gap between
its PS and its primary null 95th percentile, against the distance between the two
registered null constructions — a far larger perturbation than reseeding one of
them. Two families for RiNALMo and none for ERNIE-RNA have a gap smaller than that
distance. Were every one of them to cross the wrong way, the counts would be 26
and 28 against the raised gate of 8. No repair removes a scored family, since none
of the eleven carries a scored Rung 3 record, and entering families can only add.
H1(a) holds with a margin of at least 18 families.

**H1(b) is open, and is reported whichever way it lands.** It is a one-sample
Wilcoxon across families, so the three entering families can move it in either
direction and nothing above bounds it. The manuscript reports the statistic, its
p-value and the verdict against the registered alpha of 0.0167 whether or not the
test passes.

**Reporting.** The manuscript reports the analysis as run on the repaired data,
N = 36, and states the H1(a) verdict at both the registered gate of 7 and the
raised gate of 8 in the main text. A supplementary panel table gives all 52
families with Rfam accession, seed member and length, and records the disposition
of each repaired record. This supersedes the reporting decision in the entry
below, which was taken when the correction moved N by one and left the gate
unchanged: a correction that leaves the decision threshold in place is a defect in
the inputs and the corrected number stands alone, while a correction that moves
the threshold raises a robustness question the reader is entitled to see answered.

**Residual exposure.** Families enter a confirmatory set after results for the
others were known. The eligibility rule was frozen, is mechanical, and is applied
by committed code; this entry predates the entering families being scored; the
gate change is shown not to affect H1(a) by an argument that does not depend on
their values. None of that makes the sequence invisible, and it is stated rather
than argued away. `5S_rRNA_ecoli` is left unrepaired and out of the panel, and it
was the panel's only rRNA family, so the analyzed panel now spans eight classes
rather than nine and carries no ribosomal RNA. Ribosomal RNA is the structural
class with the deepest experimental annotation, and no result here speaks to it.
The confirmatory set carries no transfer RNA either. The panel's two tRNA
families are the two quarantined for pilot foreknowledge, so the 36 families the
Rung 3 hypotheses are decided on contain neither of the two classes whose folds
are best characterized; Rungs 1 and 2, where the quarantine does not apply, do
include both. `scripts/registered_quarantine.py` reads the quarantine out of the
frozen registration rather than repeating it, checks that the document's two
statements of it agree, and reports what it leaves.
The panel composition is generated from the records by
`scripts/generate_panel_description.py` rather than described by hand.

---

## 2026-08-24 — Phase 6 confirmatory set moves from 32 families to 33

**Registration:** `PREREGISTRATION_PHASE6_V2.md` (2026-07-13, SHA `c19aa59`).

**What changed.** Commit `5fc8914` replaced the annotations for
`mir_122_precursor`, `mir_155_precursor`, `mir_21_precursor` and
`mir_let7_precursor`. Before it, the first three carried one byte-identical
sequence and dot-bracket — `K02350.1/1-119`, an RF00001 (5S ribosomal RNA) seed
member — under three microRNA names, and the fourth carried a hand-made 72-nt
hairpin also labeled RF00001. All four now carry their own Rfam seed members
(RF00684, RF00731, RF00658, RF00027) and all 52 families are distinct.

**Consequence for the registered set.** Applying the frozen eligibility criteria
to the corrected annotations, `mir_21_precursor` becomes eligible: 16 canonical
Watson-Crick pairs, 3 stems, 7 eligible interior pairs, against thresholds of 15
and 5. Eligible families move from 34 to 35 and the confirmatory set from N = 32
to N = 33. The other three corrected families were eligible before and remain
eligible. The quarantine is untouched: `tRNA_Phe_yeast` and `tRNA_Ala_human` were
not among the four.

**No registered rule changed.** The eligibility criteria — at least 15 canonical
WC pairs, stems of at least 3 consecutive WC pairs, the two terminal pairs of each
stem excluded, at least 5 eligible interior pairs remaining — are as frozen.
`scripts/scope_rerun.py` applies them through the analysis package's own
`_parse_stems` and `_get_eligible_pairs` rather than restating them, so the
eligibility determination is made by the same code that made it at N = 32.
`PREREGISTRATION_PHASE6_V2.md:29` states in advance that "the exact count of
eligible families depends on data quality," and the document nowhere enumerates
the 32 families by name.

**The registered decision threshold does not move.** H1 condition (a) is the only
confirmatory criterion carrying N: exceedances must reach `ceil(4 × 0.05 × N)`.
That is `ceil(6.4) = 7` at N = 32 and `ceil(6.6) = 7` at N = 33. H1(b), H2 and H3
fix their thresholds independently of N. `scripts/check_registered_n_sensitivity.py`
recomputes both.

**Timing.** This entry is written before any model has been run against the
corrected annotations. `mir_21_precursor` has no computed PS value under the
preregistered metric at the time of writing, so the registration's closing
statement — "PS is fully a priori across all families" — still holds for the
family entering the set.

**Reporting.** The manuscript reports the analysis as run on the deposited data:
N = 33, derived from the registered criteria. The superseded count is not printed.
The registration's both-counts convention at line 79 covers the primary and
independently-max'd nulls, which are two defensible readings of the same data; a
count computed from annotations that were wrong is not a second reading, and
printing it beside the correct one would ask a reader to adjudicate a defect in
the inputs. This record and `docs/CHANGELOG.md` carry the correction.

**Residual exposure.** A family enters a confirmatory set after results for the
other 32 were known. The eligibility rule was frozen, is mechanical, and is
applied by committed code, and this entry predates the family being scored; none
of that makes the sequence invisible, and it is stated rather than argued away.

**Defect in the frozen document.** `PREREGISTRATION_PHASE6_V2.md:29` refers to a
"Families Pending" section that the document does not contain. The eligibility
criteria are stated in full at lines 29–37 and in the kill criteria at lines
128–131, so nothing is missing from the plan; the cross-reference points at a
section that was never written. The frozen file is left as it stands.

**Status of the stored results.** Every file in `results/` was computed against
the superseded annotations and none has been re-run. Until the four families are
re-evaluated, `data/` and `results/` are inconsistent. See `docs/CHANGELOG.md`
for what moves and by how much, and `scripts/audit_duplicate_families.py` to
reproduce it.

---

## Data corrections (2026-08-23)

**Three family files were duplicates of a fourth.** `mir_122_precursor.json`,
`mir_155_precursor.json` and `mir_let7_precursor.json` in `data/rfam_families/`
contained byte-identical sequence and structure, all labeled `RF00001` (5S rRNA)
with source `K02350.1/1-119`. They were counted as three independent families in
every exceedance count, every bootstrap resample over families, and every binomial
test. The defect is visible in the stored results: six models return bit-identical
`best_ratio` across the three (Evo 0.9703185608691253; RNA-FM 1.26115305684212).

Corrected by rebuilding each from its own Rfam seed alignment: let-7 `RF00027`,
mir-122 `RF00684`, mir-155 `RF00731`. Protocol, hypotheses and decision thresholds
are unchanged; only the input data is corrected.

**mir-21 carried the wrong accession and an unbalanced structure.**
`mir_21_precursor.json` held a genuine mir-21 sequence but was labeled `RF00001`
and its dot-bracket had 33 opening against 29 closing brackets, which placed it
among the three families excluded post hoc for unbalanced annotations. Rebuilt from
`RF00658`.

**Consequence for reported counts.** With mir-21 valid, the annotation filter now
excludes two families (`U2_snRNA_stem`, `hammerhead_ribozyme`) rather than three.
All four rebuilt families clear the Rung 3 gate (>= 15 Watson-Crick pairs, >= 3 per
stem), so the qualifying panel changes and all model results require recomputation.
The previously reported values (RiNALMo mean PS 0.2150, ERNIE-RNA 0.1204, both on
31/34 gate-passing families) are superseded.

**Sequence-structure misalignment under investigation.** An audit of canonical
pairing across all 52 families finds eight where more than a quarter of annotated
pairs are non-canonical: `SRP_RNA_helix8` (0.77), `5S_rRNA_ecoli` (0.70),
`HDV_ribozyme` (0.58), `TPP_riboswitch` (0.57), `IRES_HCV_domainII` (0.55),
`U2_snRNA_stem` (0.54), `SAM_riboswitch` (0.52), `hammerhead_ribozyme` (0.44).
Ribozymes and riboswitches do contain genuine non-canonical pairs, so the lower
values are not by themselves evidence of error; rates near 0.7 are not explicable
that way. Stem and loop labels at Rungs 1 and 2 derive from these annotations.
Resolution pending.

## Registered analyses not reported in v10

The following are required by the registrations and are absent from v10. They are
reported in the revision.

- **Conservative (max-over-layers) null.** `PREREGISTRATION_PHASE6_V2.md` requires
  that where the primary and conservative nulls disagree on which families exceed
  threshold, both counts are reported. They disagree substantially (RNA-FM 4 vs 0;
  SpliceBERT 14 vs 3; UTR-LM 13 vs 3; HyenaDNA 10 vs 5). Only the primary is
  reported.
- **Quarantine.** The registration sets N = 32 confirmatory families after
  quarantining `tRNA_Phe_yeast` and `tRNA_Ala_human`, and excludes them from all
  confirmatory hypotheses. The manuscript uses 34 throughout and does not mention
  the quarantine.
- **Registered hypotheses absent from the summary table.** H6b (RNA-FM exceeds the
  null in >= 9 families) fails and is not listed; H8 (probing margin >= 0.02 for at
  least one model) passes and is not listed.
- **Ten further hypotheses from two registrations absent entirely.**
  `docs/PREREG_PHASE3_RNA_PRETRAINED.md` registers H12-H15 and
  `docs/PREREG_PHASE4_EXPANDED_MODELS.md` registers H16-H21. None appears in the
  hypothesis summary. Three fail against the paper's own values: H13 (UTR-LM
  attention rho > 0.15; observed 0.046-0.052), H20 (|DNABERT-2 trained minus
  untrained| < 0.05; observed 0.257), H21 (DNABERT-2 trained rho < 0.20; observed
  0.257). The Discussion features the H20 result as a finding without noting that
  it falsifies the prediction.
- **Wilcoxon component of H1_6.** The registration requires a one-sample Wilcoxon
  signed-rank test alongside the family count. Only counts are reported.
- **Registered exploratory analyses not run or not reported.** E2 (guanine-cytosine
  versus adenine-uracil pair type; the values are stored in every result file), E3
  (forward/reverse asymmetry), and E5, the loop-mutation control, which distinguishes
  partner specificity from non-specific perturbation propagation along the backbone.

## Analyses not covered by any registration

The following are reported in the manuscript and are post hoc. The manuscript
describes all evaluations as preregistered; that statement is incorrect and is
being changed.

- Transversion control
- Synthetic covariation control. `PREREGISTRATION_PHASE6_V2.md` states of this
  experiment: "Distinguishing these would require synthetic sequences with known
  structure but no evolutionary covariation. This is noted as a limitation, **not
  tested**."
- Mutual-information baseline against Rfam seed alignments
- Within-family coefficient-of-variation study. The manuscript cites a
  "preregistered threshold of 10%"; no registration contains such a threshold.
- ERNIE-RNA attention-bias zeroing ablation

## Reporting error: Table 5 does not use the specified computation

`run_phase6` computes `mean_best_ps` over families that are not skipped, not
quarantined, and that pass the positive-control gate -- the exclusion the Methods
state. Table 5 reports a mean over all non-skipped families, gate failures
included. The two computations are compared by `scripts/check_reported_values.py`:

| model | pipeline (as specified) | Table 5 | |
|---|---|---|---|
| RiNALMo | 0.2265 (N=29) | 0.2150 | |
| ERNIE-RNA | 0.1170 (N=30) | 0.1204 | |
| Evo | 0.001052 (N=10) | 0.0011 | |
| SpliceBERT | 0.000350 (N=21) | 0.0002 | |
| RNA-FM | 0.0000488 (N=8) | 0.0001 | |
| UTR-LM | 0.0000137 (N=16) | 0.00001 | |
| HyenaDNA | 0.00000097 (N=25) | 0.0003 | 307x |
| ERNIE-RNA untrained | 0.0000000255 (N=10) | 0.000000095 | |
| NT v2 | 0.0 (N=1) | 0.001 | pipeline is exactly zero |

The ordering below the top two models does not survive the correction. Table 5
and Figure 3B rank NT v2 above SpliceBERT, RNA-FM, UTR-LM and HyenaDNA; on the
pipeline's values NT v2 is last among trained models. Discussion of the middle of
the panel was written against the incorrect ordering.

Tables are regenerated from the pipeline's own output.

## Conservative null

`scripts/check_reported_values.py` also prints families exceeding the conservative
(max-over-layers) null, which `PREREGISTRATION_PHASE6_V2.md` requires be reported
alongside the primary wherever the two disagree:

| model | primary | conservative |
|---|---|---|
| ERNIE-RNA | 28 | 28 |
| RiNALMo | 28 | 26 |
| SpliceBERT | 14 | 3 |
| UTR-LM | 13 | 3 |
| HyenaDNA | 10 | 5 |
| RNA-FM | 4 | 0 |
| ERNIE-RNA untrained | 9 | 4 |

## Untrained controls: run, not reported

`data/gpu_results/phase6_untrained/` holds Phase 6 untrained controls for nine
models. None had been committed to any repository, and the manuscript reports
only ERNIE-RNA. Mean perturbation specificity, untrained:

| model | mean PS | N | exceeding null |
|---|---|---|---|
| RiNALMo | 8.0e-09 | 6 | 3 |
| RNA-FM | 5.7e-09 | 3 | 2 |
| Caduceus | -5.5e-09 | 30 | 15 |
| SpliceBERT | -7.2e-09 | 11 | 4 |
| Nucleotide transformer v2 | -9.9e-09 | 1 | 0 |
| HyenaDNA | -1.3e-08 | 18 | 7 |
| Evo | -1.4e-08 | 2 | 0 |
| UTR-LM | -3.2e-08 | 23 | 8 |
| DNABERT-2 | -7.1e-03 | 1 | 0 |

`PREREGISTRATION_PHASE6_UNTRAINED_RINALMO.md` registers **H_null: untrained
RiNALMo produces mean PS indistinguishable from zero (|PS| < 0.001)**. Observed
8.0e-09. H_null holds by five orders of magnitude, and the registration's
requirement to "report alongside the ERNIE-RNA untrained control in the paper"
is unmet in v10.

Every architecture, without learned weights, sits within 3.2e-08 of zero. The
exceedance counts in the right-hand column are the negative-threshold artifact
described under the derangement null: models scoring at zero clear a null whose
95th percentile is below zero.
