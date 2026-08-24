# Provenance

Dated record of what was measured, what changed, and which file settled it.
Nothing here belongs in the manuscript.

## 2026-08-24 — Table 5 was computed over 34 families; the registration specifies 32

`PREREGISTRATION_PHASE6_V2.md` (frozen 2026-07-13, commit `891d6af` after the
history rewrite) states: "Eligible families: 34 of 52 families pass the 15-WC-pair
threshold. After quarantine (2 families), N = 32 confirmatory families," and "All
confirmatory hypotheses exclude quarantined families." The quarantined families are
`tRNA_Phe_yeast` and `tRNA_Ala_human`, both used with RNA-FM while the
perturbation-specificity metric was being developed.

Table 5 in `paper/paper_v10.tex` was built over all 34.
`scripts/audit_table5_aggregation.py` establishes this rather than inferring it: the
printed gate numerators reproduce `pass/34` for all nine models that report one, and
the printed H3 precisions reproduce the mean of per-family fractions over
gate-passing families out of 34 for all eight. Neither reproduces at 32. Matching is
done at each literal's own printed precision, since a value shown to one significant
figure cannot be matched with a percentage tolerance.

**The pipeline applied the quarantine correctly.** Every Phase 6 results file stores
`families_quarantined: ["tRNA_Ala_human", "tRNA_Phe_yeast"]`, carries a per-family
`quarantined` flag, and its stored `mean_best_ps` equals the non-quarantined
gate-passing mean. The defect is in how the manuscript table was assembled from
`per_rna`, not in the runs. No model was re-run and no GPU work was done:
`per_rna` already stores each family's `best_ps`, `positive_control`,
`exceeds_null_primary` and `h3_precision`.

`scripts/generate_table5_registered.py` recomputes every Phase 6 quantity over the
32. Mean PS is taken over all 32 rather than over the 32 that also clear the positive
control: the gate is a separate screening criterion with its own column, H1's count
criterion is already evaluated among gate-passing families, and the gated mean would
be a second methodological change the manuscript never made. The gated mean is
printed beside it as a sensitivity. `paper/paper_v11.tex` is generated from v10 by
`paper/patches/patch_v11_registered_phase6.py` (18 exact-match-once edits, with
assertions that no figure computed over 34 survives); v10 is untouched.
`scripts/verify_paper_phase6_figures.py` re-derives all 11 Table 5 rows and 15 prose
figures from `results/` and compares them to the .tex.

**No verdict moves.** H1₆ and H3₆ pass, H2₆ fails, and the untrained control still
collapses. Restoring both quarantined families would also change no verdict, which
the new Appendix A reports (RiNALMo 30/31 and H3 0.882; ERNIE-RNA 29/31 and H3 0.874;
H2₆ fails either way).

### Two figures that were already wrong in v10

- **"two or more orders of magnitude below the two leaders."** At the v10 values this
  was 0.2150/0.0034 = 63× and 0.1204/0.0034 = 35×; at the registered values it is 60×
  and 32×. v11 says "at least a factor of 30" in all three places it appears.
- **H2₆ rank-biserial.** v10 prints `rb = −0.28, p = 0.265`. It does not reproduce
  from the Table 5 means under either sign convention, at either N = 32 or N = 34;
  both give `rb = +0.04, p = 1.000` from `mannwhitneyu` over 5 RNA and 5 DNA models.
  The FAIL verdict is unchanged, and v11 prints the reproducible figure. Where v10's
  −0.28 came from is not recoverable from the stored results.

### Order of operations

`scripts/check_freeze_order.py` confirms no run in this repository is dated before
the preregistration governing it, across all five phases. The check reads timestamps
from directory names, not file names, which is where they live.
