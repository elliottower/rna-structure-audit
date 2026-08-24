r"""Write the manuscript's results tables from the stamped repaired-panel re-run.

Run:  PYTHONPATH=scripts uv run --no-project --with numpy --with scipy \
          --with tqdm --python 3.12 python scripts/generate_results_tables.py

Every number in Rungs 1-3, the attention trained-versus-untrained table and the
hypothesis summary comes from `results/repaired_panel/<model>/`, and each of
those directories carries a stamp naming the commit, the panel hash and the
library versions the numbers were produced under. The tables are written to
`paper/generated/` and reach the manuscript through \input, so the next
correction to the data corrects the manuscript instead of leaving a number typed
into the body for someone to find later.

The hypothesis summary is assembled from the registrations here for the first
time. paper_v12's version was carried by hand: it printed eight rows where the
registrations define twenty-one hypotheses, omitting H6b and H8
(docs/PREREGISTRATION_EXPANDED_RFAM.md), H12-H15
(docs/PREREG_PHASE3_RNA_PRETRAINED.md) and H16-H21
(docs/PREREG_PHASE4_EXPANDED_MODELS.md), and it printed RNA-FM's stem/loop ratio
against H6, whose registered quantity is the trained-to-untrained ratio. Each
row below names its registration and is decided by an evaluator that reads the
stored per-family values.

H6, H14 and H19 are registered as a "mean trained/untrained ratio". That wording
admits the mean of per-family ratios and the ratio of per-family means, and the
registration does not choose. The table prints the mean of per-family ratios,
which is what "mean ... ratio across N families" says; results_macros.tex
carries both, so a reader can see what the choice costs.

H1 is Phase 1's underpowered version of H6, decided on twelve families that no
longer exist as a set: five of the seven records still naming Phase 1
provenance are withdrawn, and the two that remain are the quarantined pilot
families. It is reported at its Phase 1 value against the panel it was decided
on rather than recomputed against a panel it was never registered for.
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

from compute_bootstrap_cis import N_BOOTSTRAP, SEED, bootstrap_ci
from generate_table5_registered import QUARANTINE, mean, render_ps, summarize_entries

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results" / "repaired_panel"
GENERATED = REPO / "paper" / "generated"
FAMILIES = REPO / "data" / "rfam_families"

# Asserted, not derived: a generator that reports whatever it happens to find
# cannot tell a repaired panel from a half-fetched one.
PANEL_N = 47
PANEL_CURATED = 52

# Phase 1, docs/PREREGISTRATION_EXPANDED_RFAM.md line 13: "H1 (trained/untrained
# ratio >= 2.0): FAIL at 1.74x". Twelve families, five of them now withdrawn.
PHASE1_H1_RATIO = 1.74
PHASE1_H1_N = 12

# (key, short name, parameter count, domain, macro-safe name).
MODELS = [
    ("ernierna", "ERNIE-RNA", "86M", "RNA", "ErnieRNA"),
    ("rinalmo", "RiNALMo", "650M", "RNA", "RiNALMo"),
    ("rnafm", "RNA-FM", "99M", "RNA", "RNAFM"),
    ("utrlm", "UTR-LM", r"$\sim$2M", "RNA", "UTRLM"),
    ("splicebert", "SpliceBERT", "19M", "RNA", "SpliceBERT"),
    ("nt", "NT~v2", "56M", "DNA", "NTvTwo"),
    ("dnabert2", "DNABERT-2", "117M", "DNA", "DNABERTTwo"),
    ("hyenadna", "HyenaDNA", "5.4M", "DNA", "HyenaDNA"),
    ("caduceus", "Caduceus", "14M", "DNA", "Caduceus"),
    ("evo", "Evo", "7B", "DNA", "Evo"),
]

# Every model carrying a registered trained-versus-untrained hypothesis, plus
# UTR-LM, whose untrained row the attention table needs and no hypothesis names.
# HyenaDNA, Caduceus and Evo have no attention and no such hypothesis, so an
# untrained run of them would serve nothing that is reported.
UNTRAINED_KEYS = ["ernierna", "rinalmo", "rnafm", "utrlm", "splicebert", "nt",
                  "dnabert2"]

RNA_KEYS = [k for k, _, _, domain, _ in MODELS if domain == "RNA"]
DNA_KEYS = [k for k, _, _, domain, _ in MODELS if domain == "DNA"]

BY_KEY = {key: (short, size, domain, macro) for key, short, size, domain, macro in MODELS}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load(key: str) -> dict:
    """The two result files and the stamp for one model directory."""
    directory = RESULTS / key
    phases = directory / f"{key}_phases_1_to_5.json"
    phase6 = directory / f"{key}_phase6_ps.json"
    stamp = directory / "stamp.json"
    missing = [p for p in (phases, phase6, stamp) if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"{key}: missing " + ", ".join(str(p.relative_to(REPO)) for p in missing))
    return {
        "phases": json.loads(phases.read_text()),
        "phase6": json.loads(phase6.read_text())["results"],
        "stamp": json.loads(stamp.read_text()),
    }


# The stamp fields that determine the numbers. A run disagreeing on any of them
# analyzed different data or drew different randomness, and a table mixing it
# with the others is a table no reader can attribute.
IDENTICAL_FIELDS = ("panel_sha256", "n_families", "withdrawn", "seeding")


def load_all() -> dict:
    """Every model that ran, with the stamps checked against each other.

    The commit is deliberately not in `IDENTICAL_FIELDS`. No stack loads all ten
    models, so the runs are launched in groups and a later group carries a later
    commit; requiring one commit would mean re-running every finished model to
    change a string. What it costs is that the reader has to be told which cell
    came from which commit, so `provenance_table.tex` prints the commit, the
    library stack and the device for all seventeen runs.
    """
    keys = [key for key, *_ in MODELS] + [f"{key}_untrained" for key in UNTRAINED_KEYS]
    runs = {key: load(key) for key in keys}
    stamps = {key: run["stamp"] for key, run in runs.items()}
    for field in IDENTICAL_FIELDS:
        values = {key: json.dumps(stamp[field], sort_keys=True)
                  for key, stamp in stamps.items()}
        if len(set(values.values())) != 1:
            lines = "\n".join(f"    {key:22s} {value}" for key, value in sorted(values.items()))
            raise ValueError(f"models disagree on stamp field {field!r}:\n{lines}")
    any_stamp = next(iter(stamps.values()))
    if any_stamp["n_families"] != PANEL_N:
        raise ValueError(f"the run loaded {any_stamp['n_families']} families, "
                         f"not the repaired panel's {PANEL_N}")
    if len(any_stamp["withdrawn"]) != PANEL_CURATED - PANEL_N:
        raise ValueError(f"{len(any_stamp['withdrawn'])} withdrawn records, "
                         f"not {PANEL_CURATED - PANEL_N}")
    commits = sorted({stamp["commit"] for stamp in stamps.values()})
    print(f"  {len(runs)} runs on panel {any_stamp['panel_sha256'][:12]}, "
          f"{any_stamp['n_families']} families, "
          f"{len(commits)} commit(s): {', '.join(c[:7] for c in commits)}")
    return runs


def scored(section: dict) -> dict:
    """Per-family entries that hold a result rather than a skip reason."""
    if section.get("skipped"):
        return {}
    return {name: body for name, body in section.get("per_rna", {}).items()
            if isinstance(body, dict) and "skipped" not in body}


# ---------------------------------------------------------------------------
# Rung 1 and Rung 2
# ---------------------------------------------------------------------------


def mutation_stats(run: dict, rng) -> dict:
    """Mean ratio with its bootstrap interval, and null exceedance counts."""
    entries = scored(run["phases"]["mutation_trained"])
    ratios = np.array([body["best_ratio"] for body in entries.values()])
    exceeds_nuc = [name for name, body in entries.items() if body.get("exceeds_nuc_null")]
    survives = [name for name in exceeds_nuc if entries[name].get("exceeds_dinuc_null")]
    interval = bootstrap_ci(ratios, np.mean, rng) if len(ratios) else None
    retention = None
    if exceeds_nuc:
        flags = np.array([float(name in survives) for name in exceeds_nuc])
        retention = bootstrap_ci(flags, np.mean, rng)
    return {
        "n_scored": len(ratios),
        "mean_ratio": float(np.mean(ratios)) if len(ratios) else None,
        "ci": interval,
        "exceeds_nuc": len(exceeds_nuc),
        "survives_dinuc": len(survives),
        "retention": retention,
        "per_family": {name: body["best_ratio"] for name, body in entries.items()},
    }


def attention_mean(run: dict) -> float | None:
    """Mean best-head-best-layer Spearman rho, or None for an SSM."""
    entries = scored(run["phases"]["attention_trained"])
    return mean([body["best_corr"] for body in entries.values()]) if entries else None


def probing_accuracy(run: dict) -> float | None:
    return run["phases"]["probing"].get("best_accuracy")


def fmt_ci(interval: dict | None, digits: int = 2) -> str:
    if interval is None:
        return "---"
    return f"[{interval['ci_lower']:.{digits}f}, {interval['ci_upper']:.{digits}f}]"


def rung1_table(runs: dict, rng) -> str:
    rows = []
    for key, short, size, domain, _macro in MODELS:
        stats_ = mutation_stats(runs[key], rng)
        rho = attention_mean(runs[key])
        accuracy = probing_accuracy(runs[key])
        rows.append(
            f"{short} ({size}, {domain})".ljust(26)
            + f" & {stats_['mean_ratio']:.3f} & {fmt_ci(stats_['ci'])}"
            + f" & {stats_['exceeds_nuc']}/{stats_['n_scored']}"
            + (f" & {rho:.3f}" if rho is not None else " & ---")
            + (f" & ${accuracy - 0.5:+.3f}$" if accuracy is not None else " & ---")
            + r" \\")
    untrained = []
    for key in UNTRAINED_KEYS:
        short = BY_KEY[key][0]
        run = runs[f"{key}_untrained"]
        stats_ = mutation_stats(run, rng)
        rho = attention_mean(run)
        accuracy = probing_accuracy(run)
        untrained.append(
            f"{short} untrained".ljust(26)
            + f" & {stats_['mean_ratio']:.3f} & {fmt_ci(stats_['ci'])}"
            + f" & {stats_['exceeds_nuc']}/{stats_['n_scored']}"
            + (f" & {rho:.3f}" if rho is not None else " & ---")
            + (f" & ${accuracy - 0.5:+.3f}$" if accuracy is not None else " & ---")
            + r" \\")
    n_scored = mutation_stats(runs["ernierna"], rng)["n_scored"]
    return r"""\begin{table}[htbp]
\centering
\caption{Mutation sensitivity across ten models on the repaired panel
($N = \panelN{}$ Rfam families). Bootstrap 95\% CIs on the mean ratio
resample per-family ratios (""" + f"$B = {N_BOOTSTRAP:,}".replace(",", "{,}") + r"""$).
At $\alpha = 0.05$ per family, """ + f"{0.05 * n_scored:.1f} of {n_scored}" + r"""
families are expected to exceed the null by chance. Attention $\rho$ is the mean
best-head-best-layer Spearman correlation between symmetrized attention and the
contact map, and is undefined for the three state-space models. Probing $\Delta$
is balanced accuracy minus 0.5. Untrained rows randomize every weight at a fixed
seed and load through the same adapter as the trained row above them.}
\label{tab:rung1}
\small
\begin{tabular}{lrlccc}
\toprule
\textbf{Model} & \textbf{Mean ratio} & \textbf{95\% CI} & \textbf{Fam $>$ null} & \textbf{Attn $\rho$} & \textbf{Probing $\Delta$} \\
\midrule
""" + "\n".join(rows) + "\n\\midrule\n" + "\n".join(untrained) + r"""
\bottomrule
\end{tabular}
\end{table}
"""


def rung2_table(runs: dict, rng) -> str:
    trained = sorted(
        ((key, mutation_stats(runs[key], rng)) for key, *_ in MODELS),
        key=lambda pair: -pair[1]["exceeds_nuc"])
    rows = []
    for key, stats_ in trained:
        short, size, domain, _macro = BY_KEY[key]
        retention = (f"{100 * stats_['survives_dinuc'] / stats_['exceeds_nuc']:.0f}\\%"
                     if stats_["exceeds_nuc"] else "---")
        rows.append(f"{short} ({size}, {domain})".ljust(26)
                    + f" & {stats_['exceeds_nuc']} & {stats_['survives_dinuc']}"
                    + f" & {retention} & {fmt_ci(stats_['retention'])} \\\\")
    untrained = []
    for key in UNTRAINED_KEYS:
        stats_ = mutation_stats(runs[f"{key}_untrained"], rng)
        retention = (f"{100 * stats_['survives_dinuc'] / stats_['exceeds_nuc']:.0f}\\%"
                     if stats_["exceeds_nuc"] else "---")
        untrained.append(f"{BY_KEY[key][0]} untrained".ljust(26)
                         + f" & {stats_['exceeds_nuc']} & {stats_['survives_dinuc']}"
                         + f" & {retention} & {fmt_ci(stats_['retention'])} \\\\")
    return r"""\begin{table}[htbp]
\centering
\caption{Dinucleotide-stratified null across all models on the repaired panel.
``Fam $>$ nuc'' counts families exceeding the nucleotide-stratified null;
``Survive dinuc'' counts how many of those also exceed the
dinucleotide-stratified null, which is computed only where the first-order null
is exceeded. Retention is the second count over the first, and its bootstrap
95\% CI resamples the families in the first.}
\label{tab:rung2}
\small
\begin{tabular}{lrrll}
\toprule
\textbf{Model} & \textbf{Fam $>$ nuc} & \textbf{Survive dinuc} & \textbf{Retention} & \textbf{95\% CI} \\
\midrule
""" + "\n".join(rows) + "\n\\midrule\n" + "\n".join(untrained) + r"""
\bottomrule
\end{tabular}
\end{table}
"""


def attention_table(runs: dict) -> str:
    """Trained against untrained rho, for the seven models that have attention."""
    rows = []
    for key in UNTRAINED_KEYS:
        trained = attention_mean(runs[key])
        untrained = attention_mean(runs[f"{key}_untrained"])
        if trained is None or untrained is None:
            continue
        delta = trained - untrained
        reading = ("Architectural" if abs(delta) < 0.025
                   else ("Learned" if delta > 0 else "Untrained exceeds trained"))
        rows.append((trained, f"{BY_KEY[key][0]:<12} & {trained:.3f} & {untrained:.3f}"
                     f" & {reading} ($\\Delta = {delta:+.3f}$) \\\\"))
    body = "\n".join(row for _rho, row in sorted(rows, key=lambda pair: -pair[0]))
    return r"""\begin{table}[htbp]
\centering
\caption{Attention-contact Spearman correlation, trained against randomly
initialized weights, for the seven models in the panel that expose attention.
$\Delta$ within $\pm 0.025$ is read as architectural: the correlation is present
before any pretraining. Both columns are means over the repaired panel, computed
with eager attention.}
\label{tab:attn_trained_untrained}
\small
\begin{tabular}{@{}lccl@{}}
\toprule
\textbf{Model} & \textbf{Trained} & \textbf{Untrained} & \textbf{Interpretation} \\
\midrule
""" + body + r"""
\bottomrule
\end{tabular}
\end{table}
"""


# ---------------------------------------------------------------------------
# Rung 3
# ---------------------------------------------------------------------------


def rung3_stats(run: dict) -> dict:
    """The registered Table 5 aggregation, plus what the hypothesis tests need."""
    entries = {name: body for name, body in run["phase6"]["per_rna"].items()
               if isinstance(body, dict) and "skipped" not in body
               and isinstance(body.get("best_ps"), (int, float))}
    summary = summarize_entries(entries)
    nonq = {name: body for name, body in entries.items() if name not in QUARANTINE}
    passing = {name: body for name, body in nonq.items()
               if body.get("positive_control", {}).get("pass")}
    summary["ps_values"] = [body["best_ps"] for body in nonq.values()]
    summary["exceed_conservative"] = sum(
        1 for body in passing.values() if body.get("exceeds_null_conservative") is True)
    summary["h3_counts"] = (
        sum(body["h3_precision"]["partner_max_count"] for body in passing.values()
            if isinstance(body.get("h3_precision"), dict)
            and body["h3_precision"].get("n")),
        sum(body["h3_precision"]["n"] for body in passing.values()
            if isinstance(body.get("h3_precision"), dict)
            and body["h3_precision"].get("n")),
    )
    return summary


def gate_threshold(n: int) -> int:
    """The registered count criterion, ceil(4 * 0.05 * N), at the panel's N."""
    return math.ceil(4 * 0.05 * n)


def h1_six(stats_: dict) -> tuple[bool, bool, float]:
    """The two registered conditions, at both the registered gate and the panel's.

    Returns (passes at the repaired gate, passes at the registered 7, Wilcoxon p).
    """
    values = np.array(stats_["ps_values"])
    if len(values) < 2 or not np.any(values):
        return False, False, 1.0
    p_value = float(stats.wilcoxon(values, alternative="greater").pvalue)
    positive = stats_["mean_ps"] is not None and stats_["mean_ps"] > 0 and p_value < 0.0167
    return (positive and stats_["exceed"] >= gate_threshold(stats_["eligible"]),
            positive and stats_["exceed"] >= 7,
            p_value)


def rung3_table(runs: dict, rung3: dict) -> str:
    order = sorted((key for key, *_ in MODELS),
                   key=lambda key: -(rung3[key]["mean_ps"] or 0.0))
    rows = []
    for key in order:
        stats_ = rung3[key]
        short, size, domain, _macro = BY_KEY[key]
        passes, _registered, _p = h1_six(stats_)
        label = f"{short} ({size})"
        ps_cell = "---" if stats_["mean_ps"] is None else render_ps(stats_["mean_ps"])
        null_cell = "---" if not stats_["gate"] else f"{stats_['exceed']}/{stats_['gate']}"
        conservative = ("---" if not stats_["gate"]
                        else f"{stats_['exceed_conservative']}/{stats_['gate']}")
        h3_cell = "---" if stats_["h3"] is None else f"{stats_['h3']:.3f}"
        if passes:
            label = f"\\textbf{{{short}}} ({size})"
            ps_cell, null_cell = f"\\textbf{{{ps_cell}}}", f"\\textbf{{{null_cell}}}"
        rows.append(f"{label:<28} & {domain} & {ps_cell} & {stats_['gate']}/"
                    f"{stats_['eligible']} & {null_cell} & {conservative} & {h3_cell} \\\\")
    untrained = []
    for key in UNTRAINED_KEYS:
        stats_ = rung3[f"{key}_untrained"]
        short, size, domain, _macro = BY_KEY[key]
        ps_cell = "---" if stats_["mean_ps"] is None else render_ps(stats_["mean_ps"])
        null_cell = "---" if not stats_["gate"] else f"{stats_['exceed']}/{stats_['gate']}"
        conservative = ("---" if not stats_["gate"]
                        else f"{stats_['exceed_conservative']}/{stats_['gate']}")
        h3_cell = "---" if stats_["h3"] is None else f"{stats_['h3']:.3f}"
        untrained.append(f"{short + ' untrained':<28} & {domain} & {ps_cell} & "
                         f"{stats_['gate']}/{stats_['eligible']} & {null_cell} & "
                         f"{conservative} & {h3_cell} \\\\")
    eligible = rung3["ernierna"]["eligible"]
    return r"""\begin{table}[htbp]
\centering
\caption{Perturbation specificity on the repaired panel. """ + (
        f"{eligible + len(QUARANTINE)} of \\panelN{{}} families qualify "
        r"($\geq 15$ WC pairs with $\geq 3$ per stem); the two pilot families "
        f"are quarantined, leaving $N = {eligible}$ confirmatory families. Mean "
        r"PS is taken over those $N$. Gate counts families passing the "
        r"stem-versus-loop positive control, and the two null columns are "
        r"computed among gate-passing families: $>$~null is the primary "
        r"within-stem derangement null and $>$~cons.\ the conservative variant, "
        r"reported together because the registration defines both. H3 precision "
        r"is the mean of per-family partner-is-max fractions. Bold rows satisfy "
        r"H1$_6$ at the repaired gate of "
        f"{gate_threshold(eligible)}. NT~v2 and DNABERT-2 qualify on fewer "
        r"families because their tokenizers resolve short families into too few "
        r"tokens to perturb.}") + r"""
\label{tab:rung3}
\small
\begin{tabular}{llrcccc}
\toprule
\textbf{Model} & \textbf{Domain} & \textbf{Mean PS} & \textbf{Gate} & \textbf{$>$ null} & \textbf{$>$ cons.} & \textbf{H3 prec.} \\
\midrule
""" + "\n".join(rows) + "\n\\midrule\n" + "\n".join(untrained) + r"""
\bottomrule
\end{tabular}
\end{table}
"""


# ---------------------------------------------------------------------------
# Registered hypotheses
# ---------------------------------------------------------------------------


def sign_test(runs: dict, key: str) -> tuple[int, int, float]:
    """Families where the trained ratio beats the untrained one, and a binomial p."""
    trained = mutation_stats(runs[key], np.random.default_rng(SEED))["per_family"]
    untrained = mutation_stats(runs[f"{key}_untrained"],
                               np.random.default_rng(SEED))["per_family"]
    shared = sorted(set(trained) & set(untrained))
    wins = sum(1 for name in shared if trained[name] > untrained[name])
    p_value = stats.binomtest(wins, len(shared), 0.5, alternative="greater").pvalue
    return wins, len(shared), float(p_value)


def ratio_readings(runs: dict, key: str) -> tuple[float, float]:
    """Mean of per-family trained/untrained ratios, and ratio of the two means."""
    trained = mutation_stats(runs[key], np.random.default_rng(SEED))["per_family"]
    untrained = mutation_stats(runs[f"{key}_untrained"],
                               np.random.default_rng(SEED))["per_family"]
    shared = sorted(name for name in set(trained) & set(untrained)
                    if untrained[name] > 0)
    per_family = float(np.mean([trained[name] / untrained[name] for name in shared]))
    of_means = float(np.mean([trained[name] for name in shared])
                     / np.mean([untrained[name] for name in shared]))
    return per_family, of_means


def verdict(passes: bool) -> str:
    return r"\textbf{PASS}" if passes else "FAIL"


def hypothesis_rows(runs: dict, rung3: dict) -> list[tuple[str, str, str, str, str]]:
    """(tag, criterion, panel, result, verdict) for every registered hypothesis."""
    rng = np.random.default_rng(SEED)
    mutation = {key: mutation_stats(runs[key], rng) for key in runs}
    rows = []

    rows.append(("H1", r"RNA-FM trained/untrained ratio $\geq 2.0$",
                 f"Phase 1, $N = {PHASE1_H1_N}$",
                 f"{PHASE1_H1_RATIO:.2f}$\\times$", verdict(False)))

    for tag, key in [("H6", "rnafm"), ("H14", "rinalmo"), ("H19", "splicebert")]:
        per_family, _of_means = ratio_readings(runs, key)
        # H19 is registered in the passing direction: SpliceBERT is expected to
        # fall below 2.0, and the registration says so explicitly.
        below = tag == "H19"
        passes = (per_family < 2.0) if below else (per_family >= 2.0)
        comparator = "$<$" if below else r"$\geq$"
        rows.append((tag, f"{BY_KEY[key][0]} trained/untrained ratio {comparator} 2.0",
                     r"repaired, $N = \panelN{}$",
                     f"{per_family:.2f}$\\times$", verdict(passes)))

    rnafm = mutation["rnafm"]
    threshold = math.ceil(0.17 * rnafm["n_scored"])
    rows.append(("H6b", f"RNA-FM exceeds nuc.\\ null in $\\geq \\lceil 0.17 N \\rceil$ "
                 f"$= {threshold}$ families", r"repaired, $N = \panelN{}$",
                 f"{rnafm['exceeds_nuc']}/{rnafm['n_scored']}",
                 verdict(rnafm["exceeds_nuc"] >= threshold)))

    best_dinuc = max(((key, mutation[key]["survives_dinuc"]) for key, *_ in MODELS),
                     key=lambda pair: pair[1])
    rows.append(("H7", r"$\geq 1$ family survives the dinuc.\ null",
                 r"repaired, $N = \panelN{}$",
                 f"{best_dinuc[1]} fam ({BY_KEY[best_dinuc[0]][0]})",
                 verdict(best_dinuc[1] >= 1)))

    probes = [(key, probing_accuracy(runs[key])) for key, *_ in MODELS]
    best_probe = max((pair for pair in probes if pair[1] is not None),
                     key=lambda pair: pair[1])
    rows.append(("H8", r"best probing accuracy $>$ 0.5 by $\geq 0.02$",
                 r"repaired, $N = \panelN{}$",
                 f"{best_probe[1]:.3f} ({BY_KEY[best_probe[0]][0]})",
                 verdict(best_probe[1] - 0.5 >= 0.02)))

    for tag, key in [("H10", "rnafm"), ("H15", "rinalmo"), ("H17", "ernierna")]:
        wins, total, p_value = sign_test(runs, key)
        passes = wins >= 0.75 * total and p_value < 0.01
        rows.append((tag, f"{BY_KEY[key][0]} trained $>$ untrained in "
                     r"$\geq 75\%$ of families", r"repaired, $N = \panelN{}$",
                     f"{wins}/{total} = {100 * wins / total:.0f}\\%", verdict(passes)))

    nt_trained = attention_mean(runs["nt"])
    nt_untrained = attention_mean(runs["nt_untrained"])
    rows.append(("H11", r"NT~v2 attn.\ $> 0.15$ trained and $< 0.05$ untrained",
                 r"repaired, $N = \panelN{}$",
                 f"{nt_trained:.3f} / {nt_untrained:.3f}",
                 verdict(nt_trained > 0.15 and nt_untrained < 0.05)))

    for tag, key, bound, passes_below in [("H12", "rinalmo", 0.10, True),
                                          ("H13", "utrlm", 0.15, False),
                                          ("H16", "ernierna", 0.10, True),
                                          ("H18", "splicebert", 0.10, True),
                                          ("H21", "dnabert2", 0.20, True)]:
        rho = attention_mean(runs[key])
        comparator = "$<$" if passes_below else "$>$"
        passes = (rho < bound) if passes_below else (rho > bound)
        rows.append((tag, f"{BY_KEY[key][0]} attn.\\ $\\rho$ {comparator} {bound:.2f}",
                     r"repaired, $N = \panelN{}$", f"{rho:.3f}", verdict(passes)))

    delta = abs(attention_mean(runs["dnabert2"]) - attention_mean(runs["dnabert2_untrained"]))
    rows.append(("H20", r"$|$DNABERT-2 trained $-$ untrained attn.$| < 0.05$",
                 r"repaired, $N = \panelN{}$", f"{delta:.3f}", verdict(delta < 0.05)))

    eligible = rung3["ernierna"]["eligible"]
    gate = gate_threshold(eligible)
    leaders = []
    for key, *_ in MODELS:
        passes, registered, _p = h1_six(rung3[key])
        if passes or registered:
            leaders.append((key, rung3[key], passes, registered))
    leaders.sort(key=lambda item: -item[1]["mean_ps"])
    if leaders:
        detail = "; ".join(f"{BY_KEY[key][0]}: {stats_['exceed']}/{stats_['gate']}"
                           for key, stats_, _p, _r in leaders[:2])
        any_repaired = any(passes for _k, _s, passes, _r in leaders)
    else:
        detail, any_repaired = "no model", False
    rows.append((r"H1$_6$", f"$\\geq 1$ model: PS $> 0$, $\\geq {gate}$ fam $>$ null",
                 f"confirmatory, $N = {eligible}$", detail, verdict(any_repaired)))

    rna = [rung3[key]["mean_ps"] or 0.0 for key in RNA_KEYS]
    dna = [rung3[key]["mean_ps"] or 0.0 for key in DNA_KEYS]
    u_stat, p_value = stats.mannwhitneyu(rna, dna, alternative="two-sided")
    rank_biserial = 2 * u_stat / (len(rna) * len(dna)) - 1
    rows.append((r"H2$_6$", r"RNA PS $>$ DNA PS (rb $> 0.5$)",
                 f"confirmatory, $N = {eligible}$",
                 f"rb $= {rank_biserial:+.2f}$, $p = {p_value:.3f}$",
                 verdict(rank_biserial > 0.5)))

    h3_leaders = sorted(((key, rung3[key]) for key, *_ in MODELS),
                        key=lambda pair: -(pair[1]["h3"] or 0.0))[:2]
    h3_pass = False
    for key, stats_ in h3_leaders:
        count, total = stats_["h3_counts"]
        if total and stats.binomtest(count, total, 1 / 3,
                                     alternative="greater").pvalue < 0.0167:
            h3_pass = True
    rows.append((r"H3$_6$", r"partner-is-max $> 1/3$",
                 f"confirmatory, $N = {eligible}$",
                 "; ".join(f"{BY_KEY[key][0]}: {stats_['h3']:.3f}"
                           for key, stats_ in h3_leaders),
                 verdict(h3_pass)))
    return rows


def hypotheses_table(runs: dict, rung3: dict) -> str:
    rows = hypothesis_rows(runs, rung3)
    body = "\n".join(f"{tag:<8} & {criterion} & {panel} & {result} & {mark} \\\\"
                     for tag, criterion, panel, result, mark in rows)
    eligible = rung3["ernierna"]["eligible"]
    return r"""\begin{table}[htbp]
\centering
\caption{Every registered hypothesis, with the panel it is decided on. Rungs~1
and~2 are decided on the repaired panel of \panelN{} families; Rung~3 on the
""" + str(eligible) + r""" confirmatory families that qualify after quarantine.
H1 is Phase~1's underpowered form of H6 and is reported at the twelve-family
value it was decided on, because five of the twelve records are withdrawn and
the panel it names no longer exists. H1$_6$ is a single decision requiring both
a positive mean PS by Wilcoxon signed-rank ($p < 0.0167$) and the count
criterion; the count criterion was registered as 7 at $N = 32$ and is
""" + str(gate_threshold(eligible)) + r""" at the repaired $N$, and both are
reported in the text. Sources: \texttt{PREREGISTRATION\_EXPANDED\_RFAM}
(H1, H6, H6b, H7, H8, H10, H11), \texttt{PREREG\_PHASE3\_RNA\_PRETRAINED}
(H12--H15), \texttt{PREREG\_PHASE4\_EXPANDED\_MODELS} (H16--H21),
\texttt{PREREGISTRATION\_PHASE6\_V2} (H1$_6$--H3$_6$).}
\label{tab:hypotheses}
\small
\begin{tabular}{lllll}
\toprule
\textbf{Hyp.} & \textbf{Criterion} & \textbf{Panel} & \textbf{Result} & \textbf{Verdict} \\
\midrule
""" + body + r"""
\bottomrule
\end{tabular}
\end{table}
"""


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def stack(stamp: dict) -> str:
    """The library versions that could move a number, as one cell."""
    libs = stamp["libraries"]
    parts = [f"torch {libs['torch']}", f"transformers {libs['transformers']}"]
    for name in ("multimolecule", "flash-attn", "mamba-ssm"):
        if libs.get(name):
            parts.append(f"{name} {libs[name]}")
    return ", ".join(parts)


def provenance_table(runs: dict) -> str:
    """One row per container, so every cell in the results tables is attributable."""
    rows = []
    for key, short, _size, _domain, _macro in MODELS:
        for name, weights in ((key, "pretrained"),
                              (f"{key}_untrained", "randomized")):
            if name not in runs:
                continue
            stamp = runs[name]["stamp"]
            rows.append(" & ".join([
                short if weights == "pretrained" else f"{short} (untrained)",
                weights,
                stamp["device"],
                r"\texttt{" + stamp["commit"][:7] + "}",
                r"\footnotesize " + stack(stamp),
            ]) + r" \\")

    any_stamp = next(iter(runs.values()))["stamp"]
    return "\n".join([
        r"% Written by scripts/generate_results_tables.py. Do not edit.",
        r"\begin{table}[htbp]",
        r"\centering",
        r"\small",
        r"\caption{Provenance for each run behind the results tables. Every run "
        r"analyzes the same panel (SHA-256 \texttt{" + any_stamp["panel_sha256"][:12] +
        r"}, " + str(any_stamp["n_families"]) + " families, " +
        str(len(any_stamp["withdrawn"])) + r" withdrawn) under the same "
        r"family-derived seeding. No single library stack loads all ten models, "
        r"so the stack is given per run, and so is the commit.}",
        r"\label{tab:provenance}",
        r"\begin{tabular}{lllll}",
        r"\toprule",
        r"Model & Weights & Device & Commit & Library stack \\",
        r"\midrule",
        *rows,
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ])


# ---------------------------------------------------------------------------
# Macros for the prose
# ---------------------------------------------------------------------------


def macros(runs: dict, rung3: dict) -> str:
    rng = np.random.default_rng(SEED)
    lines = [r"% Written by scripts/generate_results_tables.py. Do not edit."]

    def macro(name: str, value: str) -> None:
        lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")

    eligible = rung3["ernierna"]["eligible"]
    macro("confirmatoryN", str(eligible))
    macro("gateThreshold", str(gate_threshold(eligible)))
    macro("gateRegistered", "7")
    macro("bootstrapB", f"{N_BOOTSTRAP:,}".replace(",", "{,}"))

    any_stamp = next(iter(runs.values()))["stamp"]
    macro("panelHashShort", r"\texttt{" + any_stamp["panel_sha256"][:12] + "}")
    macro("provenanceRuns", str(len(runs)))
    macro("provenanceCommits",
          str(len({run["stamp"]["commit"] for run in runs.values()})))

    for key, _short, _size, _domain, name in MODELS:
        stats_ = mutation_stats(runs[key], rng)
        macro(f"ratio{name}", f"{stats_['mean_ratio']:.3f}")
        macro(f"exceed{name}", f"{stats_['exceeds_nuc']}")
        macro(f"dinuc{name}", f"{stats_['survives_dinuc']}")
        accuracy = probing_accuracy(runs[key])
        if accuracy is not None:
            macro(f"probe{name}", f"{accuracy:.3f}")
        rho = attention_mean(runs[key])
        if rho is not None:
            macro(f"attn{name}", f"{rho:.3f}")
        three = rung3[key]
        if three["mean_ps"] is not None:
            macro(f"ps{name}", render_ps(three["mean_ps"]))
        macro(f"gate{name}", f"{three['gate']}/{three['eligible']}")
        if three["gate"]:
            macro(f"null{name}", f"{three['exceed']}/{three['gate']}")
            macro(f"cons{name}", f"{three['exceed_conservative']}/{three['gate']}")
        if three["h3"] is not None:
            macro(f"hthree{name}", f"{three['h3']:.3f}")

    macro("familiesScored", str(mutation_stats(runs["ernierna"], rng)["n_scored"]))
    macro("expectedFalsePositives",
          f"{0.05 * mutation_stats(runs['ernierna'], rng)['n_scored']:.1f}")

    # Both readings of the ambiguous "mean trained/untrained ratio", so the
    # choice the table makes is visible rather than buried in this script.
    for key in ["rnafm", "rinalmo", "splicebert"]:
        per_family, of_means = ratio_readings(runs, key)
        macro(f"tuPerFamily{BY_KEY[key][3]}", f"{per_family:.2f}")
        macro(f"tuOfMeans{BY_KEY[key][3]}", f"{of_means:.2f}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------


def main() -> int:
    runs = load_all()
    rung3 = {key: rung3_stats(run) for key, run in runs.items()}
    rng = np.random.default_rng(SEED)
    GENERATED.mkdir(parents=True, exist_ok=True)
    written = {
        "rung1_table.tex": rung1_table(runs, rng),
        "rung2_table.tex": rung2_table(runs, rng),
        "rung3_table.tex": rung3_table(runs, rung3),
        "attention_table.tex": attention_table(runs),
        "provenance_table.tex": provenance_table(runs),
        "hypotheses_table.tex": hypotheses_table(runs, rung3),
        "results_macros.tex": macros(runs, rung3),
    }
    for name, text in written.items():
        (GENERATED / name).write_text(text)
        print(f"  wrote paper/generated/{name}  ({len(text.splitlines())} lines)")

    print("\nRegistered hypotheses:")
    for tag, criterion, panel, result, mark in hypothesis_rows(runs, rung3):
        plain = mark.replace(r"\textbf{", "").replace("}", "")
        print(f"  {tag:<8} {plain:<5} {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
