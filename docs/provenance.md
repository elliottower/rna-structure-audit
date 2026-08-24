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

## 2026-08-24 — two scripts wrote `results/bootstrap_cis.json`, and the paper quoted both

`scripts/compute_bootstrap_cis.py` and `paper/compute_cis.py` bootstrapped the same
per-family ratios, drew from one seeded generator in different orders, and wrote the
same output filename in incompatible schemas. Running the second overwrote the
artifact the first produced with a flat list no other script could read. Both were in
the tree and neither said which was authoritative.

`scripts/audit_table12_sources.py` attributes each of the 42 interval bounds in
Tables 1 and 2 to its producer, comparing at each literal's own printed precision and
separately against a bound widened outward. The two bootstraps agree to about 0.001,
so 36 bounds are consistent with either. The six that discriminate settle it:

| bound | printed in v11 | `compute_bootstrap_cis` | `compute_cis` |
|---|---|---|---|
| RiNALMo, Rung 1 lower | 1.11 | 1.1045 | 1.1049 |
| RNA-FM, Rung 1 lower | 1.55 | 1.5566 | **1.5500** |
| RNA-FM, Rung 1 upper | 1.98 | 1.9714 | **1.9771** |
| Evo, Rung 1 upper | 1.59 | 1.5824 | **1.5860** |
| ERNIE-RNA untrained, Rung 1 lower | 1.63 | 1.6393 | **1.6336** |
| DNABERT-2, Rung 2 upper | 0.57 | 0.7143 | **0.5714** |

Five reproduce from `compute_cis.py` under round-to-nearest and from the committed
artifact only if the bound is widened outward, which the other 36 bounds are not. So
`paper/compute_cis.py` produced the printed tables. RiNALMo's lower bound reproduces
from neither: both scripts put it at 1.104, and 1.11 rounds a lower bound up, which
narrows the interval.

DNABERT-2's is the one large disagreement and it is a property of the quantity.
Bootstrapping 7 binary outcomes puts the 97.5th percentile on a seven-point lattice,
so two draw orders land a whole family apart. The point estimate, 2 of 7, is the same
either way.

**Resolution.** `scripts/compute_bootstrap_cis.py` is the source of record and
`results/bootstrap_cis.json` is its output. `paper/compute_cis.py` moved to
`scripts/superseded/compute_cis.py`, its output path changed so it can no longer
overwrite the artifact, its docstring stating why it is retired. `paper_v12.tex` is
generated from v11 by `paper/patches/patch_v12_reproducible_intervals.py` and moves
six bounds by at most 0.01, four of them by one unit in the last printed place. No
point estimate, family count, or verdict changes; RiNALMo's interval still excludes
1.0. `scripts/verify_paper_rung12_figures.py` re-derives all 22 rows from the
artifact and passes.

### The Benjamini-Hochberg counts came from nowhere

v11's Limitations section reported "ERNIE-RNA: 31 → 28; RiNALMo: 18 → 15" under
Benjamini-Hochberg at FDR 0.05. Neither script produces those numbers.
`compute_bootstrap_cis.py` computes no p-values at all. `compute_cis.py` gives
31 → 31 and 18 → 0, and its p-values are not measured but assigned from a lookup
table keyed on effect size: 0.001 above 2× the null threshold, 0.005 above 1.5×,
0.01 above 1.2×, 0.03 otherwise. Benjamini-Hochberg over those is arithmetic on
invented inputs.

The correction is not computable from what this repository stores. Each result file
records a binary exceedance against each family's own 95th-percentile null and not
the permutation distribution behind it, so no continuous per-family p-value exists to
correct. v12 says so and rests the global control on the binomial test already
reported in the same paragraph, which is computed from the exceedance counts and
needs no p-value per family.

### Path repair in the Rung 1–2 chain

`scripts/compute_bootstrap_cis.py` and `paper/generate_figures.py` read from absolute
paths inside `rna-structure-awareness` and `causal-rna`, the second of which no
longer exists. Five inputs the bootstrap needs were outside the tree:
`rnafm_expanded_20260716.json`, `nt_expanded_20260716.json`,
`hyenadna_expanded_20260716.json`, `caduceus_expanded_20260716.json` and
`evo_expanded_20260716.json`. Every file shared between the old and new trees was
compared byte for byte before any copy. The five are now in
`data/gpu_results/expanded_rfam/`, the paths are relative to the repository root, and
re-running the script reproduces the committed artifact with zero numeric
differences — `timestamp` and the recorded input paths are the only fields that move.

### Two attention correlations have no stored run

`scripts/audit_attention_rho.py` reads every file in the tree that parses as JSON,
computes the mean of `attention_trained.per_rna[*].best_corr` over scored families
the way `scripts/compute_bootstrap_cis.py` does, and attributes each value Table 1
prints. Candidates are found by parsing rather than by globbing `*.json`: three
result files under `data/gpu_results/expanded_rfam_rerun/` carry no extension, and
they are 52-family reruns.

Seven of the nine printed correlations resolve. DNABERT-2's 0.257 is
`dnabert2_expanded_20260713_200727.json` over 52 families, matching its footnote.
Two do not:

- **NT v2, 0.322.** The footnote reports it as corrected from 0.230 under the SDPA
  backend, recomputed with eager attention across 52 families. Both NT files here
  give 0.2301 over 49 scored families and are identical to each other. No
  aggregation of either reaches 0.322: mean over all entries, mean of absolute
  values, median, max, sum over 52 and mean of the top half give 0.2301, 0.2301,
  0.2033, 0.6547, 0.2168 and 0.2988.
- **RNA-FM untrained, 0.061.** No file in the tree carries an attention block for
  this model at all.

Neither is in `rna-structure-awareness`, the tree this repository was rewritten
from, which was searched the same way. The measurements are not recoverable, and
0.230 is not restorable in place of 0.322 because the footnote reports it as a
backend artifact. v12 narrows the Availability section to name both rather than
leaving "All ... per-family result files" standing. `audit_attention_rho.py` holds
the two in `KNOWN_UNRECORDED` and fails on any gap not listed there, or on either
of these reproducing later.
