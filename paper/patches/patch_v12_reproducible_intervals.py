"""Build paper_v12.tex from paper_v11.tex.

Run:  uv run --no-project --python 3.12 python \
          paper/patches/patch_v12_reproducible_intervals.py

Eight edits. Six reconcile the Rung 1 and Rung 2 intervals to
results/bootstrap_cis.json, the seventh removes a multiple-testing correction no
script in this repository produces, and the eighth narrows an availability claim
to what the deposit holds.

Two scripts bootstrap the same per-family ratios and differ only in the order
they draw from one seeded generator, so their intervals disagree in the third
decimal. scripts/compute_bootstrap_cis.py writes results/bootstrap_cis.json,
which the repository ships. paper/compute_cis.py writes the same filename with an
incompatible flat schema, and it produced the printed tables:
scripts/audit_table12_sources.py attributes every bound in Tables 1 and 2 to its
producer, and every one of the 22 rows matches compute_cis.py under
round-to-nearest at its printed precision, with one exception noted below.
Neither bootstrap is the more correct estimate. The paper is reconciled to the
committed artifact so that one file backs every interval it prints.

Edits 1--5 change six bounds by at most 0.01, four of them by one unit in the
last printed place:

  RiNALMo             lo  1.11 -> 1.10   (artifact 1.1045)
  RNA-FM              lo  1.55 -> 1.56   (artifact 1.5566)
  RNA-FM              hi  1.98 -> 1.97   (artifact 1.9714)
  Evo                 hi  1.59 -> 1.58   (artifact 1.5824)
  ERNIE-RNA untrained lo  1.63 -> 1.64   (artifact 1.6393)
  DNABERT-2 retention hi  0.57 -> 0.71   (artifact 0.7143)

RiNALMo's lower bound is the exception: 1.11 is reachable from neither script,
since both put the bound at 1.104, and rounding a lower bound up narrows the
interval. The bound still excludes 1.0, so the elevation claim it supports is
unchanged.

DNABERT-2's retention bound is the one large difference, and it is a property of
the quantity rather than of either script: bootstrapping seven binary outcomes
puts the 97.5th percentile on a seven-point lattice, so two draw orders land a
whole family apart. The point estimate, 2 of 7, is unaffected. NT v2's
[0.57, 1.00] is a different 0.5714 on which both scripts agree, and is left
alone.

Edit 7. The Limitations section quotes Benjamini-Hochberg counts of 31 -> 28 and
18 -> 15. No script produces them: compute_cis.py gives 31 -> 31 and 18 -> 0
from p-values it assigns by a lookup table on effect size, and
compute_bootstrap_cis.py computes no p-values at all. A correction across
families is not computable from what the result files store, which is a binary
exceedance per family and no permutation distribution. The paragraph now says so
and rests the global control on the binomial test reported two sentences
earlier.

Edit 8. The Availability section claims "All analysis scripts, per-family result
files, 52 Rfam structures, and preregistration documents" are deposited. Two
attention correlations in Table 1 have no per-family output here, in the tree this
repository was rewritten from, or under any other aggregation of the files that
are here: NT v2's 0.322 and untrained RNA-FM's 0.061, both recomputations the
table footnotes already describe. scripts/audit_attention_rho.py establishes this
by reading every file in the tree that parses as JSON. The two are named rather
than the claim being left standing, and neither figure changes: 0.230 is the
value the deposited NT v2 run gives, and the footnote already reports it as the
superseded SDPA measurement.

No verdict moves. H1, H2 and H3 at every rung are unchanged, as are all point
estimates and every family count.

v11 is not modified.
"""

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "paper" / "paper_v11.tex"
DST = REPO / "paper" / "paper_v12.tex"

BH_OLD = r"""Applying Benjamini-Hochberg correction at FDR $= 0.05$
within each model reduces the counts modestly for top models
(ERNIE-RNA: 31 $\to$ 28; RiNALMo: 18 $\to$ 15) and
zeros out floor-level models entirely. The main text reports
uncorrected per-family counts because the permutation null already
controls the per-test type~I error rate, but the BH-adjusted counts
confirm that the separation between the top two models and the rest
survives correction."""

BH_NEW = r"""The main text reports uncorrected per-family counts. A
false-discovery-rate correction across families requires a $p$-value
per family; the pipeline records each family's exceedance against its
own 95th-percentile null and nothing finer, and recovering a
continuous $p$-value would require the full permutation distribution,
which the result files do not store. The binomial test above serves as
the global control, treating each family as one $\alpha = 0.05$ test
and asking whether the observed exceedance count is compatible with
the false-positive floor."""

EDITS: list[tuple[str, str, str]] = [
    (
        "Table 1, RiNALMo: lower bound, printed 1.11 and reachable from neither script",
        r"RiNALMo (650M, RNA)     & 1.128 & [1.11, 1.15] & 18/49",
        r"RiNALMo (650M, RNA)     & 1.128 & [1.10, 1.15] & 18/49",
    ),
    (
        "Rung 1 prose: RiNALMo's interval, quoted again",
        r"RiNALMo's $[1.11, 1.15]$ exclude",
        r"RiNALMo's $[1.10, 1.15]$ exclude",
    ),
    (
        "Table 1, RNA-FM: both bounds",
        r"RNA-FM (99M, RNA)       & 1.749 & [1.55, 1.98] & 2/49",
        r"RNA-FM (99M, RNA)       & 1.749 & [1.56, 1.97] & 2/49",
    ),
    (
        "Table 1, Evo: upper bound",
        r"Evo (7B, DNA)           & 1.408 & [1.26, 1.59] & 10/49",
        r"Evo (7B, DNA)           & 1.408 & [1.26, 1.58] & 10/49",
    ),
    (
        "Table 1, ERNIE-RNA untrained: lower bound",
        r"ERNIE-RNA untrained     & 1.848 & [1.63, 2.11] & 3/49",
        r"ERNIE-RNA untrained     & 1.848 & [1.64, 2.11] & 3/49",
    ),
    (
        "Table 2, DNABERT-2: retention interval",
        r"DNABERT-2 (117M, DNA)   &  7 &  2 & 29\% & [0.00, 0.57] \\",
        r"DNABERT-2 (117M, DNA)   &  7 &  2 & 29\% & [0.00, 0.71] \\",
    ),
    (
        "Limitations: the unreproducible Benjamini-Hochberg counts",
        BH_OLD,
        BH_NEW,
    ),
    (
        "Availability: the two attention correlations the deposit does not contain",
        "under a CC-BY-4.0 license.",
        "under a CC-BY-4.0 license. Two attention correlations in\n"
        "Table~\\ref{tab:rung1} come from runs whose per-family output is not\n"
        "among the deposited files: the eager-attention recomputation for\n"
        "NT~v2 and the attention run for untrained RNA-FM.",
    ),
]


def main() -> int:
    text = SRC.read_text()
    failures = [(label, text.count(old)) for label, old, _ in EDITS
                if text.count(old) != 1]
    if failures:
        print("ABORT -- these targets did not match exactly once:\n")
        for label, count in failures:
            print(f"  [{count} matches] {label}")
        return 1
    for label, old, new in EDITS:
        text = text.replace(old, new, 1)
        print(f"  applied: {label}")

    for stale in ["[1.11, 1.15]", "[1.55, 1.98]", "[1.26, 1.59]", "[1.63, 2.11]",
                  "[0.00, 0.57]", r"31 $\to$ 28", r"18 $\to$ 15",
                  "Benjamini-Hochberg", "BH-adjusted", "FDR $= 0.05$"]:
        assert stale not in text, f"a superseded figure survived: {stale}"

    assert text.count("[1.10, 1.15]") == 2, "RiNALMo's interval is quoted twice"
    assert "[0.57, 1.00]" in text, "NT v2's own 0.5714 bound was not to be touched"

    # Every point estimate and family count the two tables report is unchanged.
    for unchanged in ["& 1.128 &", "& 1.749 &", "& 1.408 &", "& 1.848 &",
                      "18/49", "31/49", "2/49", "30 of 31 first-order",
                      "RiNALMo retains 16 of 18", "DNABERT-2 retains 2 of 7 (29\\%)",
                      "$p < 10^{-20}$", "$p < 10^{-10}$", "2.45 are expected"]:
        assert unchanged in text, f"a figure that should not have moved is gone: {unchanged}"

    DST.write_text(text)
    print(f"\n{len(EDITS)} edits applied.\nwrote {DST.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
