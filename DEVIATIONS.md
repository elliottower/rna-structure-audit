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

## Data corrections (2026-08-23)

**Three family files were duplicates of a fourth.** `mir_122_precursor.json`,
`mir_155_precursor.json` and `mir_let7_precursor.json` in `data/rfam_families/`
contained byte-identical sequence and structure, all labelled `RF00001` (5S rRNA)
with source `K02350.1/1-119`. They were counted as three independent families in
every exceedance count, every bootstrap resample over families, and every binomial
test. The defect is visible in the stored results: six models return bit-identical
`best_ratio` across the three (Evo 0.9703185608691253; RNA-FM 1.26115305684212).

Corrected by rebuilding each from its own Rfam seed alignment: let-7 `RF00027`,
mir-122 `RF00684`, mir-155 `RF00731`. Protocol, hypotheses and decision thresholds
are unchanged; only the input data is corrected.

**mir-21 carried the wrong accession and an unbalanced structure.**
`mir_21_precursor.json` held a genuine mir-21 sequence but was labelled `RF00001`
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
