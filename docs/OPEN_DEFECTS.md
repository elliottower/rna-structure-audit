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
