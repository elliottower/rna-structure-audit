# Changelog

Corrections to the data, the code, and the reported figures, newest first.
Each entry names what changed, which files settled it, and which reported
numbers move. `docs/provenance.md` carries the longer narrative.

## 2026-08-26 — 0.2.0: the package scores the corrected panel

**`load_families()` returns 47 families, not 52.** Five were withdrawn on
annotation review and the bundled copies carried no exclusion marker, so
`evaluate()` scored families whose names assert something the sequence is not —
a stem of U2 snRNA that is the whole snRNA, a specificity domain of RNase P that
is the whole RNA. They still ship, and `load_families(include_withdrawn=True)`
returns the full curated set so a result computed before the review reproduces.
`withdrawn_families()` reports what was skipped and why.

**Eight bundled annotations were the pre-correction ones.** `HDV_ribozyme`,
`SAM_riboswitch`, `TPP_riboswitch`, `hammerhead_ribozyme` and the four `mir_*`
precursors differed from `data/rfam_families/` in sequence, dot-bracket or both.
A user scoring through the package and an author scoring through the analysis
tree got different panels under the same names. All 13 differing files are
synced.

**Grades change for four models.** `_grade()` is unchanged; the panel it reads
is not.

| model | 0.1.0 | 0.2.0 | why |
|---|---|---|---|
| Nucleotide Transformer v2 | D | B | 8 families survive the dinucleotide null |
| HyenaDNA | D | B | 7 survive |
| Evo | D | B | 6 survive |
| Caduceus | D | C | 8 exceed the Rung 1 null, 4 survive Rung 2 |

ERNIE-RNA and RiNALMo stay A; SpliceBERT, UTR-LM, RNA-FM and DNABERT-2 stay D.
A B means families survive the dinucleotide null, which is weaker than partner
specificity and is not evidence that a model resolves pairing.

**Reported parameter counts corrected.** HyenaDNA is 450,712 parameters, not
5.4M — the bundled adapter loads the smallest release. Caduceus is 7,725,312,
not 14M. Evo and DNABERT-2 keep their authors' designations.

`docs/OPEN_DEFECTS.md` D24 carries the detectors.

## 2026-08-24 — four family annotations corrected; stored results not yet re-run

**Four of the 52 evaluation families carried wrong annotations, and three of
them were the same sequence.** Until commit `5fc8914` (2026-08-23),
`mir_122_precursor`, `mir_155_precursor` and `mir_let7_precursor` all held one
byte-identical sequence and dot-bracket — `K02350.1/1-119`, a seed member of
**RF00001, 5S ribosomal RNA** — under three microRNA names.
`mir_21_precursor` held a hand-made 72-nt hairpin, also labelled RF00001. The
evaluation set therefore contained 50 distinct sequences, not 52, and one 5S
rRNA sequence entered every aggregate three times.

All four now carry genuine Rfam seed members: RF00684 (mir-122), RF00731
(mir-155), RF00658 (mir-21), RF00027 (let-7). All 52 families are distinct.

**The results in `results/` predate the correction and have not been re-run.**
`data/rfam_families/` and `results/` are therefore inconsistent until the
Phase 1–6 runs are repeated. `scripts/audit_duplicate_families.py` reports what
the Rung 3 aggregates become when the duplicates are collapsed to one
representative and when they are dropped; `scripts/verify_composition_figures.py`
recomputes the stem/loop GC enrichment from the corrected annotations.

| quantity | as reported | duplicates collapsed | duplicates dropped |
|---|---|---|---|
| RiNALMo mean PS | 0.2098 (N = 32) | 0.1913 (N = 30) | 0.1811 (N = 29) |
| RiNALMo exceeding null / H3 precision | 28/29, 0.874 | 26/27, 0.865 | 25/26, 0.859 |
| ERNIE-RNA mean PS | 0.1136 (N = 32) | 0.1084 (N = 30) | 0.1056 (N = 29) |
| ERNIE-RNA exceeding null / H3 precision | 28/30, 0.870 | 26/28, 0.860 | 25/27, 0.855 |
| stem / loop GC content | 59.7% / 45.6%, +14.1 pp, 42 of 52 | — | 58.1% / 45.6%, +12.5 pp, 40 of 52 |

Both collapsing and dropping are sensitivities.
`PREREGISTRATION_PHASE6_V2.md` fixes N = 32, and changing N is a registration
matter rather than an arithmetic one.

**No verdict moves.** H1 (at least 7 gate-passing families exceeding the null)
and H3 (per-pair precision above 1/3) still hold for RiNALMo and ERNIE-RNA,
both remain roughly 30x above every other model, and the count of models
encoding pairing partners is unchanged. ERNIE-RNA's Rung 1 count of 31 of 52
families exceeding the composition null is 29 of 50 distinct sequences.

**Also observed.** DNABERT-2 returns three different Rung 1 ratios (0.9304,
0.9449, 0.9034) for the three byte-identical inputs, so that run carries an
unseeded stochastic component.

## 2026-08-24 — figures and analysis scripts recovered

`figures/` was gitignored in the pre-rewrite repository, so `figure2_composition`
and `figure3_partner_specificity` were never committed and the submitted
manuscript did not build from a clean checkout. Both figures are now tracked,
along with fourteen files that existed only in the old working tree:
`plot_figure2.py`, `plot_figure3.py`, `multi_seq_ps.py`, `parse_stockholm.py`,
`generate_supplementary_tables.py`, six Modal wrappers for the multi-sequence
and untrained-control runs, and `results/rinalmo_untrained_phase6_ps.json`.

## 2026-08-24 — H2₆ rank-biserial

`r_b = -0.28, p = 0.265` appears in `paper_v10.tex` through `paper_v12.tex` and
in the submitted manuscript. It does not reproduce from the stored Table 5
means under either sign convention at either N = 32 or N = 34;
`scripts/generate_table5_registered.py` gives `r_b = +0.04, p = 1.000` in all
four cases. Where the printed value came from is not recoverable from the
stored results. The FAIL verdict is unchanged.

## 2026-08-23 — bootstrap confidence intervals reconciled

Two scripts wrote confidence intervals for Rungs 1–2 and seeded one generator
consumed in different orders, so a few bounds diverged by about 0.001.
`scripts/compute_bootstrap_cis.py` is the source of record;
`paper/compute_cis.py` moved to `scripts/superseded/` because it assigned
p-values from a bucketed lookup on effect size before correcting them. Six
bounds moved; no point estimate or verdict did. The BH-corrected counts
(31 to 28, 18 to 15) were removed: nothing in the pipeline produces them, and
they are not computable from stored results, which record a binary exceedance
per family against its own 95th-percentile null rather than the permutation
distribution.

## 2026-08-23 — Table 5 aggregation

`paper_v10.tex` computed Table 5 over all 34 eligible families.
`PREREGISTRATION_PHASE6_V2.md` fixes N = 32 after quarantining
`tRNA_Phe_yeast` and `tRNA_Ala_human` as pilot data. The pipeline applied the
quarantine correctly and only the table was assembled wrong. No verdict moved.
