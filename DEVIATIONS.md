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

**What the correction can and cannot move.** H1(a) cannot flip. Every null here is
a within-family derangement, so a family's exceedance flag does not depend on
which other families are present, and no repair removes a scored family — none of
the eleven carries a scored Rung 3 record. The stored counts are 28 of 29 for
RiNALMo and 28 of 30 for ERNIE-RNA, both clearing the raised gate of 8 by 20, and
entering families can only raise them. H1(b) is a Wilcoxon across families and has
no such argument: the entering families can move it either way, and its outcome is
left open.

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
