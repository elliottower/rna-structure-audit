"""Build paper_v11.tex from paper_v10.tex: Phase 6 as preregistered, new SHAs, new URL.

Run:  uv run --no-project --python 3.12 python \
          paper/patches/patch_v11_registered_phase6.py

PREREGISTRATION_PHASE6_V2.md quarantines tRNA_Phe_yeast and tRNA_Ala_human as
pilot data and fixes the confirmatory set at N = 32. v10 reports Phase 6 over
all 34 eligible families: audit_table5_aggregation.py shows the printed gate
numerators reproduce pass/34 for all nine models and the printed H3 precisions
reproduce mean-of-fractions over gate-passing families out of 34 for all eight
reporting one. Every Phase 6 figure here comes from
generate_table5_registered.py, which recomputes from stored per-family values.

Three further edits are independent of the quarantine:

- The repository moved to github.com/elliottower/rna-structure-audit.
- Filtering the weight-geometry work out of the history renumbered every
  commit, so the three preregistration SHAs are remapped by map_prereg_shas.py.
- "Two orders of magnitude" overstates the separation between the two leaders
  and the third model, and did so at v10's own numbers as well (0.2150/0.0034 =
  63x, 0.1204/0.0034 = 35x). The separation is stated as a factor throughout.

v10 is not modified.
"""

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "paper" / "paper_v10.tex"
DST = REPO / "paper" / "paper_v11.tex"

TABLE_ROWS_OLD = r"""\textbf{RiNALMo} (650M)    & RNA & \textbf{0.2150} & 31/34 & \textbf{30/31} & \textbf{0.882} \\
\textbf{ERNIE-RNA} (86M)  & RNA & \textbf{0.1204} & 31/34 & \textbf{29/31} & 0.874 \\
ERNIE-RNA untrained        & RNA & $9.5 \times 10^{-8}$ & 10/34 & 9/10 & --- \\
Caduceus (14M)             & DNA & 0.0034 & 16/34 & 13/16 & 0.416 \\
Evo (7B)                   & DNA & 0.0011 & 10/34 & 4/10 & 0.278 \\
SpliceBERT (19M)           & RNA & 0.0002 & 22/34 & 15/22 & 0.276 \\
HyenaDNA (5.4M)            & DNA & 0.0003 & 27/34 & 10/27 & 0.083 \\
RNA-FM (99M)               & RNA & 0.0001 & 8/34 & 4/8 & 0.325 \\
UTR-LM ($\sim$2M)          & RNA & 0.00001 & 16/34 & 13/16 & 0.350 \\
NT~v2 (56M)                & DNA & 0.001 & 1/4 & 0/1 & --- \\
DNABERT-2 (117M)           & DNA & $-0.016$ & 0/8 & --- & --- \\"""

TABLE_ROWS_NEW = r"""\textbf{RiNALMo} (650M)    & RNA & \textbf{0.2098} & 29/32 & \textbf{28/29} & \textbf{0.874} \\
\textbf{ERNIE-RNA} (86M)  & RNA & \textbf{0.1136} & 30/32 & \textbf{28/30} & 0.870 \\
ERNIE-RNA untrained        & RNA & $7.3 \times 10^{-8}$ & 10/32 & 9/10 & 0.351 \\
Caduceus (14M)             & DNA & 0.0035 & 16/32 & 13/16 & 0.416 \\
Evo (7B)                   & DNA & 0.0015 & 10/32 & 4/10 & 0.278 \\
SpliceBERT (19M)           & RNA & 0.0002 & 21/32 & 14/21 & 0.265 \\
HyenaDNA (5.4M)            & DNA & 0.0003 & 25/32 & 10/25 & 0.070 \\
RNA-FM (99M)               & RNA & $9.7 \times 10^{-5}$ & 8/32 & 4/8 & 0.325 \\
UTR-LM ($\sim$2M)          & RNA & $7.0 \times 10^{-6}$ & 16/32 & 13/16 & 0.350 \\
NT~v2 (56M)                & DNA & 0.0014 & 1/4 & 0/1 & 0.000 \\
DNABERT-2 (117M)           & DNA & $-0.0161$ & 0/8 & --- & --- \\"""

EDITS: list[tuple[str, str, str]] = [
    (
        "Abstract: PS values, H3 precision, and the separation as a factor",
        """which position pairs with which (PS~$= 0.120$ and $0.215$,
precision~$> 87\\%$). The remaining eight fall two or more orders of
magnitude below, with Caduceus the only one exceeding the 1/3
partner-is-max baseline.""",
        """which position pairs with which (PS~$= 0.114$ and $0.210$,
precision~$\\approx 87\\%$). The remaining eight fall at least a factor
of 30 below, with Caduceus the only one exceeding the 1/3
partner-is-max baseline.""",
    ),
    (
        "Abstract: repository URL",
        """Code, data, and preregistration documents at
\\url{https://github.com/elliottower/rna-structure-awareness}; Zenodo""",
        """Code, data, and preregistration documents at
\\url{https://github.com/elliottower/rna-structure-audit}; Zenodo""",
    ),
    (
        "Methods: preregistration SHAs remapped onto the filtered history",
        """Preregistration documents (Phase~1 SHA \\texttt{694b43b}, Phase~2
\\texttt{bd4b3fd}, Phase~6 \\texttt{c19aa59}) are deposited with the
code.""",
        """Preregistration documents (Phase~1 SHA \\texttt{a207535}, Phase~2
\\texttt{ae6712e}, Phase~6 \\texttt{891d6af}) are deposited with the
code.""",
    ),
    (
        "Methods: the quarantine and the confirmatory N",
        """the derangement null; 34 families qualify (18 excluded by pair count,
3 by the annotation filter above).""",
        """the derangement null; 34 families qualify (18 excluded by pair count,
3 by the annotation filter above). Two of the 34, tRNA\\_Phe\\_yeast and
tRNA\\_Ala\\_human, were used while the perturbation-specificity metric was
being developed, and the preregistration excludes both from confirmatory
analysis, leaving $N = 32$. Their values are reported in
Appendix~\\ref{app:quarantine}.""",
    ),
    (
        "Figure 2 caption: separation stated as a factor",
        """RiNALMo and ERNIE-RNA separate from the remaining models by two
orders of magnitude.""",
        """RiNALMo and ERNIE-RNA separate from the remaining models by at
least a factor of 30.""",
    ),
    (
        "Table 5 caption: separation, quarantine, and the confirmatory N",
        """\\caption{Perturbation specificity across models. RiNALMo and
ERNIE-RNA separate from the pack by two orders of magnitude.
34 of 52 families qualify ($\\geq 15$ WC pairs with $\\geq 3$
per stem); NT~v2 and DNABERT-2 evaluate fewer due to tokenization.}""",
        """\\caption{Perturbation specificity across models. RiNALMo and
ERNIE-RNA separate from the pack by at least a factor of 30.
34 of 52 families qualify ($\\geq 15$ WC pairs with $\\geq 3$
per stem); two are quarantined as pilot data, leaving $N = 32$
confirmatory families. Gate is the number passing the stem-versus-loop
positive control; $>$~null and H3 precision are computed among those.
NT~v2 and DNABERT-2 evaluate fewer families due to tokenization.}""",
    ),
    ("Table 5 rows", TABLE_ROWS_OLD, TABLE_ROWS_NEW),
    (
        "Rung 3 prose: leader figures and the separation from Caduceus",
        """RiNALMo achieves the highest perturbation specificity (PS $= 0.215$),
with 30 of 31 gate-passing families exceeding the derangement null
and a partner-is-max precision of 0.882.
ERNIE-RNA follows at PS $= 0.120$, with 29/31 families exceeding the
null and H3 precision of 0.874. Both are two orders
of magnitude above the third model and far exceed the
$1/3$ chance baseline ($p \\ll 0.001$, binomial test).""",
        """RiNALMo achieves the highest perturbation specificity (PS $= 0.210$),
with 28 of 29 gate-passing families exceeding the derangement null
and a partner-is-max precision of 0.874.
ERNIE-RNA follows at PS $= 0.114$, with 28/30 families exceeding the
null and H3 precision of 0.870. Both are more than an order
of magnitude above the third model and far exceed the
$1/3$ chance baseline ($p \\ll 0.001$, binomial test).""",
    ),
    (
        "Rung 3 prose: Caduceus, SpliceBERT, HyenaDNA",
        """Caduceus (14M, DNA, BiMamba SSM) ranks third at PS $= 0.003$, with
13/16 gate-passing families exceeding the null and H3 precision of
0.416---above chance but two orders of magnitude below the leaders.
The remaining models cluster near zero PS. SpliceBERT
exceeds the null in 15 of 22 gate-passing families, but its PS
magnitude (0.0002) is three orders of magnitude below the two
leaders. HyenaDNA's H3 fraction (0.083) falls \\emph{below} chance,""",
        """Caduceus (14M, DNA, BiMamba SSM) ranks third at PS $= 0.004$, with
13/16 gate-passing families exceeding the null and H3 precision of
0.416---above chance but 32-fold below ERNIE-RNA and 60-fold below
RiNALMo. The remaining models cluster near zero PS. SpliceBERT
exceeds the null in 14 of 21 gate-passing families, but its PS
magnitude (0.0002) is roughly three orders of magnitude below the two
leaders. HyenaDNA's H3 fraction (0.070) falls \\emph{below} chance,""",
    ),
    (
        "Untrained control: PS value and gate denominators",
        """$\\text{PS} = 9.5 \\times 10^{-8}$, six orders of magnitude below
the trained model. Only 10 of 34 families pass the positive control
gate (versus 31/34 trained).""",
        """$\\text{PS} = 7.3 \\times 10^{-8}$, six orders of magnitude below
the trained model. Only 10 of 32 families pass the positive control
gate (versus 30/32 trained).""",
    ),
    (
        "Hypothesis summary: H1_6 and H3_6 counts",
        """H1$_6$ & $\\geq 1$ model: PS $> 0$, $\\geq 7$ fam $>$ null & RiNALMo: 30/31; ERNIE-RNA: 29/31 & \\textbf{PASS} \\\\
H2$_6$ & RNA PS $>$ DNA PS (rb $> 0.5$) & rb $= -0.28$, $p = 0.265$ & FAIL \\\\
H3$_6$ & Partner-is-max $> 1/3$ & RiNALMo: 0.882; ERNIE-RNA: 0.874 & \\textbf{PASS} \\\\""",
        """H1$_6$ & $\\geq 1$ model: PS $> 0$, $\\geq 7$ fam $>$ null & RiNALMo: 28/29; ERNIE-RNA: 28/30 & \\textbf{PASS} \\\\
H2$_6$ & RNA PS $>$ DNA PS (rb $> 0.5$) & rb $= +0.04$, $p = 1.000$ & FAIL \\\\
H3$_6$ & Partner-is-max $> 1/3$ & RiNALMo: 0.874; ERNIE-RNA: 0.870 & \\textbf{PASS} \\\\""",
    ),
    (
        "Hypothesis summary footnote: the RNA/DNA groups are indistinguishable",
        """higher PS than three RNA models, reversing the effect direction.}""",
        """higher PS than three RNA models, leaving the two groups
indistinguishable.}""",
    ),
    (
        "Synthetic covariation control: natural-side values",
        """ERNIE-RNA's mean PS drops from 0.120 (natural) to 0.035 (synthetic),
a 71\\% reduction. RiNALMo drops from 0.215 to 0.061, a 71\\%
reduction.""",
        """ERNIE-RNA's mean PS drops from 0.114 (natural) to 0.034 (synthetic),
a 70\\% reduction. RiNALMo drops from 0.210 to 0.060, a 71\\%
reduction.""",
    ),
    (
        "Synthetic covariation control: residual stated as a factor",
        """PS---still orders of magnitude above the weak-signal
models---indicates""",
        """PS---still roughly an order of magnitude above the weak-signal
models---indicates""",
    ),
    (
        "Availability: repository URL",
        """preregistration documents (SHAs listed in Section~3.1) are available at
\\url{https://github.com/elliottower/rna-structure-awareness} and""",
        """preregistration documents (SHAs listed in Section~3.1) are available at
\\url{https://github.com/elliottower/rna-structure-audit} and""",
    ),
    (
        "Figure 1 ladder, Rung 3 box: the registered leader means",
        "{\\footnotesize 2 of 10 models pass: RiNALMo (PS\\,=\\,0.215),\n"
        "  ERNIE-RNA (PS\\,=\\,0.120)}",
        "{\\footnotesize 2 of 10 models pass: RiNALMo (PS\\,=\\,0.210),\n"
        "  ERNIE-RNA (PS\\,=\\,0.114)}",
    ),
    (
        "Conclusion: the registered leader means and H3 precisions",
        "pairs with which (RiNALMo, $\\text{PS} = 0.215$, 88.2\\% precision;\n"
        "ERNIE-RNA, $\\text{PS} = 0.120$, 87.4\\% precision)",
        "pairs with which (RiNALMo, $\\text{PS} = 0.210$, 87.4\\% precision;\n"
        "ERNIE-RNA, $\\text{PS} = 0.114$, 87.0\\% precision)",
    ),
]

# The quarantined families are reported, as the registration requires: "Their PS
# values under the preregistered metric will be computed and reported as
# exploratory pilot observations. If either family would have changed the
# outcome of a confirmatory test, this is disclosed."
APPENDIX = r"""
\section{Quarantined pilot families}
\label{app:quarantine}

tRNA\_Phe\_yeast and tRNA\_Ala\_human were used while the
perturbation-specificity metric was being developed. Coupling-ratio values
under a weaker variant of the metric---partner perturbation against
unpaired-position perturbation, rather than the within-stem comparison used
here---were observed for both families with RNA-FM before the Phase~6
preregistration was frozen. The preregistration excludes both from every
confirmatory test and requires their values to be reported separately
(Table~\ref{tab:quarantine}).

\begin{table}[htbp]
\centering
\caption{Perturbation specificity on the two quarantined families, and the
effect of restoring them. Restoring both families changes no confirmatory
verdict: H1$_6$ and H3$_6$ pass either way, and H2$_6$ fails either way.}
\label{tab:quarantine}
\small
\begin{tabular}{lrrrr}
\toprule
\textbf{Model} & \textbf{tRNA\_Phe} & \textbf{tRNA\_Ala} & \textbf{PS ($N$=32)} & \textbf{PS ($N$=34)} \\
\midrule
RiNALMo                    & 0.5445 & 0.0521 & 0.2098 & 0.2150 \\
ERNIE-RNA                  & 0.2805 & 0.1784 & 0.1136 & 0.1204 \\
Caduceus                   & $-1.2 \times 10^{-8}$ & 0.0044 & 0.0035 & 0.0034 \\
Evo                        & $-0.0034$ & $-0.0079$ & 0.0015 & 0.0011 \\
SpliceBERT                 & 0.0007 & $-5.0 \times 10^{-8}$ & 0.0002 & 0.0002 \\
HyenaDNA                   & $-7.2 \times 10^{-8}$ & $-4.0 \times 10^{-8}$ & 0.0003 & 0.0003 \\
RNA-FM                     & $1.8 \times 10^{-8}$ & 0.0001 & $9.7 \times 10^{-5}$ & $9.5 \times 10^{-5}$ \\
UTR-LM                     & $1.9 \times 10^{-5}$ & $-2.0 \times 10^{-8}$ & $7.0 \times 10^{-6}$ & $7.2 \times 10^{-6}$ \\
\bottomrule
\end{tabular}
\end{table}

Restoring both families raises the two leaders---RiNALMo from 0.2098 to 0.2150
and ERNIE-RNA from 0.1136 to 0.1204---and leaves every other model within
rounding of its registered value. tRNA\_Phe\_yeast carries the increase for both
leaders; tRNA\_Ala\_human sits below RiNALMo's mean and above ERNIE-RNA's. The
registered analysis is the lower of the two for both leaders.
"""


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

    anchor = "\\section*{Declarations}"
    assert text.count(anchor) == 1, "cannot place the appendix"
    text = text.replace(anchor, "\\appendix\n" + APPENDIX + "\n" + anchor, 1)
    print("  applied: appendix reporting the quarantined families")

    # No Phase 6 figure computed over 34 may survive. Bare counts like 30/31
    # are not checked: Rung 2 reports a legitimate 30/31 over its own 49-family
    # set, and the appendix quotes the N=34 means on purpose.
    for stale in ["31/34", "22/34", "27/34", "16/34", "10/34", "8/34",
                  r"\textbf{0.2150}", r"\textbf{0.1204}", "PS $= 0.215$",
                  "PS $= 0.120$", "rb $= -0.28$", "0.882; ERNIE-RNA: 0.874",
                  r"PS\,=\,0.215", r"PS\,=\,0.120",
                  r"$\text{PS} = 0.215$", r"88.2\% precision"]:
        assert stale not in text, f"a figure computed over 34 survived: {stale}"
    assert "rna-structure-awareness" not in text, "old repository URL survived"
    for sha in ["694b43b", "bd4b3fd", "c19aa59"]:
        assert sha not in text, f"pre-rewrite SHA survived: {sha}"
    # The two leaders and the untrained control must still read as before.
    for kept in ["0.2098", "0.1136", "28/29", "28/30", "891d6af",
                 "app:quarantine"]:
        assert kept in text, f"{kept} missing"

    DST.write_text(text)
    print(f"\n{len(EDITS) + 1} edits applied.\nwrote {DST.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
