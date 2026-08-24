# Registered-plan deviations and disclosures

Departures from a frozen preregistration, and data corrections that change a
registered quantity without changing a registered rule. Newest first. Frozen
preregistrations are never edited; everything that would amend one is recorded
here instead.

---

## 2026-08-24 — Phase 6 confirmatory set moves from 32 families to 33

**Registration:** `PREREGISTRATION_PHASE6_V2.md` (2026-07-13, SHA `c19aa59`).

**What changed.** Commit `5fc8914` replaced the annotations for
`mir_122_precursor`, `mir_155_precursor`, `mir_21_precursor` and
`mir_let7_precursor`. Before it, the first three carried one byte-identical
sequence and dot-bracket — `K02350.1/1-119`, an RF00001 (5S ribosomal RNA) seed
member — under three microRNA names, and the fourth carried a hand-made 72-nt
hairpin also labelled RF00001. All four now carry their own Rfam seed members
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
