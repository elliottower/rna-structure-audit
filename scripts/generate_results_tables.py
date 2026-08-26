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
import zlib
from pathlib import Path

import numpy as np
from scipy import stats

from compute_bootstrap_cis import N_BOOTSTRAP, SEED, bootstrap_ci
from generate_table5_registered import QUARANTINE, mean, render_ps, summarize_entries

REPO = Path(__file__).resolve().parents[1]
# The panel every table reads. `repaired_panel` is the deposited run and predates
# the Rung 3 token-alignment repair, the float64 metric and the measured
# resolution floor; `repaired_panel_v3` is the pass that carries all three, at one
# commit per cell.
RESULTS = REPO / "results" / "repaired_panel_v3"
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
    ("utrlm", "UTR-LM", "1.2M", "RNA", "UTRLM"),
    ("splicebert", "SpliceBERT", "19M", "RNA", "SpliceBERT"),
    ("nt", "NT~v2", "56M", "DNA", "NTvTwo"),
    ("dnabert2", "DNABERT-2", "117M", "DNA", "DNABERTTwo"),
    ("hyenadna", "HyenaDNA", "0.45M", "DNA", "HyenaDNA"),
    ("caduceus", "Caduceus", "7.7M", "DNA", "Caduceus"),
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

# How far a trained attention-contact correlation may sit from its randomly
# initialized counterpart and still be read as architectural. The absolute delta
# is the wrong scale on its own -- 0.04 on NT~v2's 0.242 is a different claim
# from 0.04 on RiNALMo's 0.064 -- so the prose quotes the widest gap and the
# factor it represents alongside this bound rather than resting on it.
ATTENTION_ARCHITECTURAL_DELTA = 0.05

BY_KEY = {key: (short, size, domain, macro) for key, short, size, domain, macro in MODELS}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def shown(path: Path) -> str:
    """Relative to the repo where it is under it, absolute where it is not."""
    return str(path.relative_to(REPO)) if path.is_relative_to(REPO) else str(path)


def load(key: str) -> dict:
    """The two result files and the stamp for one model directory."""
    directory = RESULTS / key
    phases = directory / f"{key}_phases_1_to_5.json"
    phase6 = directory / f"{key}_phase6_ps.json"
    stamp = directory / "stamp.json"
    missing = [p for p in (phases, phase6, stamp) if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"{key}: missing " + ", ".join(shown(p) for p in missing))
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


def load_transversion() -> dict:
    """Every model's transversion control, from the directories beside the runs.

    The control writes to `<model>_transversion/` rather than into the run it
    controls: `_run_model` empties a directory whose stored stamp disagrees with
    the one it is about to write, and the control is launched from a later
    commit than the Watson-Crick runs, so sharing a directory would have the
    control delete the numbers it exists to be compared against.

    The substitution alphabet travels in the payload as well as in the directory
    name, and is checked here on two counts. A substitution that leaves a
    nucleotide unchanged is dropped by `run_mutation_sensitivity` from both the
    stem and the loop mean without saying so, and two files naming different
    alphabets would put two substitutions in one column.
    """
    controls = {}
    for key, *_ in MODELS:
        path = RESULTS / f"{key}_transversion" / f"{key}_transversion.json"
        if not path.exists():
            raise FileNotFoundError(f"{key}: missing {shown(path)}")
        payload = json.loads(path.read_text())
        alphabet = payload.get("complement")
        if not alphabet:
            raise ValueError(f"{key}: the transversion file does not name its alphabet")
        fixed = sorted(n for n, sub in alphabet.items() if n == sub)
        if fixed:
            raise ValueError(
                f"{key}: the transversion alphabet leaves {', '.join(fixed)} "
                "unchanged, so those positions leave both means silently")
        controls[key] = {"alphabet": alphabet, "mutation": payload["mutation_trained"]}
    alphabets = {json.dumps(c["alphabet"], sort_keys=True) for c in controls.values()}
    if len(alphabets) != 1:
        raise ValueError(f"{len(alphabets)} transversion alphabets across the panel; "
                         "the control column would mix two substitutions")
    return controls


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


def mutation_all(runs: dict) -> dict:
    """`mutation_stats` once per run, on a bootstrap stream fixed by the key.

    Threading one generator through the renderers ties every interval to the
    order the renderers happen to be called in, and the same model's CI then
    differs between the table that prints it and the macro that quotes it --
    v12 was assembled that way, with four independent streams over the same
    families. Seeding from the key makes an interval a property of the model.
    """
    return {key: mutation_stats(run, np.random.default_rng(
                [SEED, zlib.crc32(key.encode())]))
            for key, run in runs.items()}


def attention_mean(run: dict) -> float | None:
    """Mean best-head-best-layer Spearman rho, or None for an SSM."""
    entries = scored(run["phases"]["attention_trained"])
    return mean([body["best_corr"] for body in entries.values()]) if entries else None


def attention_per_family(run: dict) -> dict:
    """Best-head-best-layer rho for each family, for trained-untrained deltas."""
    return {name: body["best_corr"]
            for name, body in scored(run["phases"]["attention_trained"]).items()}


def probing_accuracy(run: dict) -> float | None:
    return run["phases"]["probing"].get("best_accuracy")


def probing_layer(run: dict) -> int | None:
    return run["phases"]["probing"].get("best_layer")


# Balanced accuracy on a two-class target, which is what the probe reports.
PROBE_CHANCE = 0.5


def fmt_ci(interval: dict | None, digits: int = 2) -> str:
    if interval is None:
        return "---"
    return f"[{interval['ci_lower']:.{digits}f}, {interval['ci_upper']:.{digits}f}]"


def rung1_table(runs: dict, mut: dict) -> str:
    rows = []
    for key, short, size, domain, _macro in MODELS:
        stats_ = mut[key]
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
        stats_ = mut[f"{key}_untrained"]
        rho = attention_mean(run)
        accuracy = probing_accuracy(run)
        untrained.append(
            f"{short} untrained".ljust(26)
            + f" & {stats_['mean_ratio']:.3f} & {fmt_ci(stats_['ci'])}"
            + f" & {stats_['exceeds_nuc']}/{stats_['n_scored']}"
            + (f" & {rho:.3f}" if rho is not None else " & ---")
            + (f" & ${accuracy - 0.5:+.3f}$" if accuracy is not None else " & ---")
            + r" \\")
    n_scored = mut["ernierna"]["n_scored"]
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


def rung2_table(mut: dict) -> str:
    trained = sorted(
        ((key, mut[key]) for key, *_ in MODELS),
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
        stats_ = mut[f"{key}_untrained"]
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
        reading = ("Architectural" if abs(delta) < ATTENTION_ARCHITECTURAL_DELTA
                   else ("Learned" if delta > 0 else "Untrained exceeds trained"))
        rows.append((trained, f"{BY_KEY[key][0]:<12} & {trained:.3f} & {untrained:.3f}"
                     f" & {reading} ($\\Delta = {delta:+.3f}$) \\\\"))
    body = "\n".join(row for _rho, row in sorted(rows, key=lambda pair: -pair[0]))
    return r"""\begin{table}[htbp]
\centering
\caption{Attention-contact Spearman correlation, trained against randomly
initialized weights, for the seven models in the panel that expose attention.
$\Delta$ within $\pm """ + f"{ATTENTION_ARCHITECTURAL_DELTA:g}" + r"""$ is read as architectural: the
correlation is present before any pretraining. Both columns are means over the
repaired panel, computed with eager attention.}
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


def transversion_table(mut: dict, controls: dict) -> str:
    """Watson-Crick against transversion mean ratio, one row per model."""
    rows = []
    for key, short, size, domain, _macro in MODELS:
        watson = mut[key]["mean_ratio"]
        control = transversion_ratio(controls[key])
        rows.append((abs(control - watson) / watson,
                     f"{short} ({size}, {domain})".ljust(26)
                     + f" & {watson:.3f} & {control:.3f}"
                     + f" & ${100 * (control - watson) / watson:+.1f}$\\%"
                     + r" \\"))
    body = "\n".join(row for _deviation, row in sorted(rows, key=lambda pair: -pair[0]))
    alphabet = next(iter(controls.values()))["alphabet"]
    written = ", ".join(f"{nuc}$\\to${sub}" for nuc, sub in sorted(alphabet.items())
                        if nuc != "T")
    return r"""\begin{table}[htbp]
\centering
\caption{Mutation sensitivity under the Watson-Crick substitution and under a
substitution that leaves no position paired with the partner it had (""" + written + r"""; T
follows U). Rows are ordered by the size of the change. A model whose ratio
depends on the substituted nucleotide being a complement rather than on the
pairing being broken would move between the two columns.}
\label{tab:transversion}
\small
\begin{tabular}{lrrr}
\toprule
\textbf{Model} & \textbf{Watson-Crick} & \textbf{Transversion} & \textbf{Change} \\
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
    # Where PS peaks, among the families that pass the gate: a peak at layer 0
    # or 1 is what a positional-encoding artifact looks like, and the count is
    # what says whether the peak is the model's or one family's.
    summary["peak_layers"] = [body["best_layer"] for body in passing.values()
                              if isinstance(body.get("best_layer"), int)]
    summary["n_ps_layers"] = max(
        (len(body["per_layer_ps"]) for body in entries.values()
         if isinstance(body.get("per_layer_ps"), list)), default=0)
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


def sign_test(mut: dict, key: str) -> tuple[int, int, float]:
    """Families where the trained ratio beats the untrained one, and a binomial p."""
    trained = mut[key]["per_family"]
    untrained = mut[f"{key}_untrained"]["per_family"]
    shared = sorted(set(trained) & set(untrained))
    wins = sum(1 for name in shared if trained[name] > untrained[name])
    p_value = stats.binomtest(wins, len(shared), 0.5, alternative="greater").pvalue
    return wins, len(shared), float(p_value)


def ratio_readings(mut: dict, key: str) -> tuple[float, float]:
    """Mean of per-family trained/untrained ratios, and ratio of the two means."""
    trained = mut[key]["per_family"]
    untrained = mut[f"{key}_untrained"]["per_family"]
    shared = sorted(name for name in set(trained) & set(untrained)
                    if untrained[name] > 0)
    per_family = float(np.mean([trained[name] / untrained[name] for name in shared]))
    of_means = float(np.mean([trained[name] for name in shared])
                     / np.mean([untrained[name] for name in shared]))
    return per_family, of_means


def domain_effect(rung3: dict) -> tuple[float, float]:
    """H2$_6$: rank-biserial and Mann-Whitney p for RNA against DNA mean PS.

    Computed here rather than twice, because the body and the hypothesis table
    both quote it and v12 printed two different values for the same test.
    """
    rna = [rung3[key]["mean_ps"] or 0.0 for key in RNA_KEYS]
    dna = [rung3[key]["mean_ps"] or 0.0 for key in DNA_KEYS]
    u_stat, p_value = stats.mannwhitneyu(rna, dna, alternative="two-sided")
    return 2 * u_stat / (len(rna) * len(dna)) - 1, float(p_value)


def verdict(passes: bool) -> str:
    return r"\textbf{PASS}" if passes else "FAIL"


def hypothesis_rows(runs: dict, rung3: dict,
                    mutation: dict) -> list[tuple[str, str, str, str, str]]:
    """(tag, criterion, panel, result, verdict) for every registered hypothesis."""
    rows = []

    rows.append(("H1", r"RNA-FM trained/untrained ratio $\geq 2.0$",
                 f"Phase 1, $N = {PHASE1_H1_N}$",
                 f"{PHASE1_H1_RATIO:.2f}$\\times$", verdict(False)))

    for tag, key in [("H6", "rnafm"), ("H14", "rinalmo"), ("H19", "splicebert")]:
        per_family, _of_means = ratio_readings(mutation, key)
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
        wins, total, p_value = sign_test(mutation, key)
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

    rank_biserial, p_value = domain_effect(rung3)
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


def hypotheses_table(runs: dict, rung3: dict, mut: dict) -> str:
    rows = hypothesis_rows(runs, rung3, mut)
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


def peak_layer(stats_: dict) -> tuple[int, int, int] | None:
    """The layer PS most often peaks at, how many families peak there, and of how many."""
    layers = stats_["peak_layers"]
    if not layers:
        return None
    modal = max(set(layers), key=layers.count)
    return modal, layers.count(modal), len(layers)


def transversion_ratio(control: dict) -> float:
    """Mean best ratio over the families the control scored."""
    return mean([body["best_ratio"] for body in scored(control["mutation"]).values()])


def claim_failures(runs: dict, rung3: dict, mut: dict, controls: dict) -> list[str]:
    """Every sentence in the body the repaired panel contradicts.

    A macro carries a number, not an ordering. The body says ERNIE-RNA exceeds
    the null in more families than any other model, that RNA-FM has the highest
    mean ratio, that RiNALMo and ERNIE-RNA separate from the third model on
    perturbation specificity; a re-run that reverses any of those leaves the
    sentence standing beside a table that contradicts it, with every number in
    the sentence correct. The claims are re-read here from the files the tables
    come from, and a flip aborts with what the ordering now is, so the sentence
    is rewritten against data rather than discovered by a reader.
    """
    ps = {key: (rung3[key]["mean_ps"] or 0.0) for key, *_ in MODELS}
    ranked_ps = sorted(ps, key=lambda key: -ps[key])
    attn = {key: attention_mean(runs[key]) for key in UNTRAINED_KEYS}
    rna_attn = {key: value for key, value in attn.items() if key in RNA_KEYS}

    failures = []

    def require(holds: bool, sentence: str, observed: str) -> None:
        if not holds:
            failures.append(f"  the body says: {sentence}\n     the panel says: {observed}")

    def short(key: str) -> str:
        return BY_KEY[key][0]

    top_exceed = max(mut, key=lambda key: mut[key]["exceeds_nuc"] if key in ps else -1)
    require(top_exceed == "ernierna",
            "ERNIE-RNA exceeds the null in more families than any other model",
            f"{short(top_exceed)} does, in {mut[top_exceed]['exceeds_nuc']}")

    top_ratio = max(ps, key=lambda key: mut[key]["mean_ratio"])
    require(top_ratio == "rnafm", "RNA-FM has the highest mean ratio",
            f"{short(top_ratio)} does, at {mut[top_ratio]['mean_ratio']:.3f}")

    top_probe = max((key for key in ps if probing_accuracy(runs[key]) is not None),
                    key=lambda key: probing_accuracy(runs[key]))
    require(top_probe == "ernierna", "ERNIE-RNA has the strongest probing signal",
            f"{short(top_probe)} does, at {probing_accuracy(runs[top_probe]):.3f}")

    for key in ("ernierna", "rinalmo"):
        require(mut[key]["ci"]["ci_lower"] > 1.0,
                f"{short(key)}'s ratio CI excludes 1.0",
                f"it is {fmt_ci(mut[key]['ci'])}")
    require(mut["dnabert2"]["ci"]["ci_lower"] < 1.0 < mut["dnabert2"]["ci"]["ci_upper"],
            "DNABERT-2's ratio CI straddles 1.0",
            f"it is {fmt_ci(mut['dnabert2']['ci'])}")

    top_attn = max(attn, key=lambda key: attn[key])
    require(top_attn == "nt", "NT~v2 has the highest attention-contact correlation",
            f"{short(top_attn)} does, at {attn[top_attn]:.3f}")
    require(attention_mean(runs["nt_untrained"]) >= attn["nt"],
            "untrained NT~v2 matches or exceeds trained NT~v2 on attention",
            f"trained {attn['nt']:.3f} against untrained "
            f"{attention_mean(runs['nt_untrained']):.3f}")
    require(min(rna_attn, key=lambda key: rna_attn[key]) == "rnafm",
            "RNA-FM ranks last on attention among the RNA models",
            f"{short(min(rna_attn, key=lambda key: rna_attn[key]))} does")
    for key in RNA_KEYS:
        delta = attn[key] - attention_mean(runs[f"{key}_untrained"])
        require(abs(delta) < ATTENTION_ARCHITECTURAL_DELTA,
                f"{short(key)}'s trained attention stays within "
                f"{ATTENTION_ARCHITECTURAL_DELTA:g} of untrained",
                f"the delta is {delta:+.3f}")

    require(ranked_ps[:3] == ["rinalmo", "ernierna", "caduceus"],
            "RiNALMo, ERNIE-RNA and Caduceus are the top three on perturbation "
            "specificity, in that order",
            "the order is " + ", ".join(short(key) for key in ranked_ps[:3]))
    third = ps[ranked_ps[2]]
    require(third > 0 and ps[ranked_ps[1]] / third >= 10,
            "the two leaders stand more than an order of magnitude above the third model",
            f"the factor is {ps[ranked_ps[1]] / third:.1f}" if third > 0
            else f"{short(ranked_ps[2])} has non-positive mean PS")
    require(ps["dnabert2"] < 0, "DNABERT-2's mean PS is negative",
            f"it is {ps['dnabert2']:.4f}")
    require(rung3["hyenadna"]["h3"] is not None and rung3["hyenadna"]["h3"] < 1 / 3,
            "HyenaDNA's H3 fraction falls below the 1/3 chance baseline",
            f"it is {rung3['hyenadna']['h3']}")
    require(ps["ernierna"] / max(rung3["ernierna_untrained"]["mean_ps"] or 0.0, 1e-30) >= 1e5,
            "untrained ERNIE-RNA sits orders of magnitude below trained on PS",
            f"trained {ps['ernierna']:.4g} against untrained "
            f"{rung3['ernierna_untrained']['mean_ps']:.4g}")
    for key in ("evo", "caduceus"):
        layer = probing_layer(runs[key])
        require(layer is not None and layer <= 1,
                f"{short(key)}'s probing accuracy peaks at layer 0 or 1",
                f"it peaks at layer {layer}")
    require(mut["evo"]["mean_ratio"] > mut["caduceus"]["mean_ratio"],
            "Evo reaches a higher mean ratio than Caduceus",
            f"Evo {mut['evo']['mean_ratio']:.3f}, "
            f"Caduceus {mut['caduceus']['mean_ratio']:.3f}")

    watson_crick = {key: mut[key]["mean_ratio"] for key, *_ in MODELS}
    for key, control in controls.items():
        require(scored(control["mutation"]).keys() == set(mut[key]["per_family"]),
                f"{short(key)}'s control covers the families its run covered",
                f"control {len(scored(control['mutation']))} families, "
                f"run {len(mut[key]['per_family'])}")
        require(transversion_ratio(control) > 0,
                f"{short(key)}'s transversion ratio is a number to report",
                f"it is {transversion_ratio(control)}")
    ranked_probe = sorted((key for key in ps if probing_accuracy(runs[key]) is not None),
                          key=lambda key: -probing_accuracy(runs[key]))
    require(ranked_probe[1] == "rinalmo", "RiNALMo has the second-strongest probing signal",
            f"{short(ranked_probe[1])} does")

    def retention(key: str) -> float:
        stats_ = mut[key]
        return (stats_["survives_dinuc"] / stats_["exceeds_nuc"]
                if stats_["exceeds_nuc"] else 0.0)

    require(mut["rnafm"]["exceeds_nuc"] == 0,
            "RNA-FM exceeds the nucleotide-stratified null in no family",
            f"it exceeds it in {mut['rnafm']['exceeds_nuc']}")
    for key in ("ernierna", "rinalmo"):
        require(retention(key) >= 0.85,
                f"{short(key)} keeps almost all of its first-order survivors",
                f"it keeps {mut[key]['survives_dinuc']} of {mut[key]['exceeds_nuc']}")
    thin = [key for key, *_ in MODELS if mut[key]["exceeds_nuc"] <= 8]
    require(len(thin) == 7,
            "seven of the ten models exceed the first-order null in eight "
            "families or fewer, so their retention fractions rest on "
            "single-digit denominators",
            f"{len(thin)} do: " + ", ".join(short(key) for key in thin))

    beaten = [key for key in RNA_KEYS if ps[key] < ps["caduceus"]]
    require(sorted(beaten) == sorted(["splicebert", "rnafm", "utrlm"]),
            "Caduceus exceeds SpliceBERT, RNA-FM and UTR-LM on perturbation specificity",
            "it exceeds " + (", ".join(short(key) for key in beaten) or "no RNA model"))

    for key in ("evo", "caduceus"):
        require(transversion_ratio(controls[key]) > 0,
                f"{short(key)}'s control ran", "it did not")

    return failures


def check_claims(runs: dict, rung3: dict, mut: dict, controls: dict) -> None:
    """Abort rather than print a superlative the table beside it contradicts."""
    failures = claim_failures(runs, rung3, mut, controls)
    if failures:
        raise ValueError("the repaired panel contradicts the body:\n"
                         + "\n".join(failures))


# ---------------------------------------------------------------------------
# Macros for the prose
# ---------------------------------------------------------------------------


def macros(runs: dict, rung3: dict, mut: dict, controls: dict) -> str:
    lines = [r"% Written by scripts/generate_results_tables.py. Do not edit."]

    def macro(name: str, value: str) -> None:
        lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")

    eligible = rung3["ernierna"]["eligible"]
    macro("confirmatoryN", str(eligible))

    # The panel narrows three times and every count has a reason, so each is
    # generated rather than typed: 52 curated records, 5 withdrawn for
    # annotations that did not match their sequences, 47 analyzed. Of those, the
    # families passing the registered Rung 3 filters, then the non-quarantined
    # ones, then the ones whose stems admit a derangement null and therefore a
    # chance rate. Prose that types any of these drifts from the tables.
    per_rna = json.loads(
        (RESULTS / "ernierna" / "ernierna_phase6_ps.json").read_text()
    )["results"]["per_rna"]
    scored = [name for name, e in per_rna.items() if e.get("best_ps") is not None]
    non_quarantined = [n for n in scored if n not in QUARANTINE]
    with_chance = [n for n in non_quarantined
                   if per_rna[n].get("h3_chance_fraction") is not None]
    macro("rungThreeQualifying", str(len(scored)))
    macro("rungThreeQuarantined", str(len(scored) - len(non_quarantined)))
    macro("chanceRateN", str(len(with_chance)))
    macro("gateThreshold", str(gate_threshold(eligible)))
    macro("gateRegistered", "7")
    macro("bootstrapB", f"{N_BOOTSTRAP:,}".replace(",", "{,}"))

    any_stamp = next(iter(runs.values()))["stamp"]
    macro("panelHashShort", r"\texttt{" + any_stamp["panel_sha256"][:12] + "}")
    macro("provenanceRuns", str(len(runs)))
    macro("provenanceCommits",
          str(len({run["stamp"]["commit"] for run in runs.values()})))
    # Counted from the versions themselves, not from the four image definitions
    # in modal_repaired_panel.py: two images sharing a pin set would make the
    # count of images wrong and the count of stacks right.
    macro("provenanceStacks",
          str(len({stack(run["stamp"]) for run in runs.values()})))

    for key, _short, _size, _domain, name in MODELS:
        stats_ = mut[key]
        macro(f"ratio{name}", f"{stats_['mean_ratio']:.3f}")
        macro(f"ci{name}", fmt_ci(stats_["ci"]))
        macro(f"exceed{name}", f"{stats_['exceeds_nuc']}")
        macro(f"dinuc{name}", f"{stats_['survives_dinuc']}")
        if stats_["exceeds_nuc"]:
            macro(f"retain{name}", f"{stats_['survives_dinuc']} of "
                  f"{stats_['exceeds_nuc']}")
            macro(f"retainPct{name}",
                  f"{100 * stats_['survives_dinuc'] / stats_['exceeds_nuc']:.0f}\\%")
        accuracy = probing_accuracy(runs[key])
        if accuracy is not None:
            macro(f"probe{name}", f"{accuracy:.3f}")
            macro(f"probeDelta{name}", f"{accuracy - PROBE_CHANCE:+.3f}")
            macro(f"probeLayer{name}", str(probing_layer(runs[key])))
        rho = attention_mean(runs[key])
        if rho is not None:
            macro(f"attn{name}", f"{rho:.3f}")
        three = rung3[key]
        if three["mean_ps"] is not None:
            macro(f"ps{name}", render_ps(three["mean_ps"]))
        macro(f"gate{name}", f"{three['gate']}/{three['eligible']}")
        macro(f"gateCount{name}", str(three["gate"]))
        macro(f"eligible{name}", str(three["eligible"]))
        if three["gate"]:
            macro(f"null{name}", f"{three['exceed']}/{three['gate']}")
            macro(f"cons{name}", f"{three['exceed_conservative']}/{three['gate']}")
        if three["h3"] is not None:
            macro(f"hthree{name}", f"{three['h3']:.3f}")
        peak = peak_layer(three)
        if peak is not None:
            modal, at_modal, total = peak
            macro(f"psPeak{name}", str(modal))
            macro(f"psLayers{name}", str(three["n_ps_layers"]))
            macro(f"psPeakCount{name}", f"{at_modal} of {total}")

    # The untrained rows the prose quotes directly: the false-positive check in
    # Rung 1 and the architecture-without-training paragraph in Rung 3.
    for key in UNTRAINED_KEYS:
        name = BY_KEY[key][3]
        stats_ = mut[f"{key}_untrained"]
        macro(f"exceedUntrained{name}", f"{stats_['exceeds_nuc']}")
        macro(f"ratioUntrained{name}", f"{stats_['mean_ratio']:.3f}")
        if stats_["exceeds_nuc"]:
            macro(f"retainUntrained{name}", f"{stats_['survives_dinuc']} of "
                  f"{stats_['exceeds_nuc']}")
        untrained_rho = attention_mean(runs[f"{key}_untrained"])
        if untrained_rho is not None and attention_mean(runs[key]) is not None:
            macro(f"attnUntrained{name}", f"{untrained_rho:.3f}")
            macro(f"attnDelta{name}",
                  f"{attention_mean(runs[key]) - untrained_rho:+.3f}")
            # Whether the mean delta is one family or all of them.
            trained_family = attention_per_family(runs[key])
            untrained_family = attention_per_family(runs[f"{key}_untrained"])
            shared = sorted(set(trained_family) & set(untrained_family))
            deltas = [trained_family[n] - untrained_family[n] for n in shared]
            higher = sum(1 for value in deltas if value > 0)
            macro(f"attnFamilyDelta{name}", f"{mean(deltas):+.3f}")
            macro(f"attnHigherTrained{name}", str(higher))
            macro(f"attnHigherUntrained{name}", str(len(deltas) - higher))
        three = rung3[f"{key}_untrained"]
        if three["mean_ps"] is not None:
            macro(f"psUntrained{name}", render_ps(three["mean_ps"]))
        macro(f"gateUntrained{name}", f"{three['gate']}/{three['eligible']}")
        if three["gate"]:
            macro(f"nullUntrained{name}", f"{three['exceed']}/{three['gate']}")
        wins, total, _p = sign_test(mut, key)
        macro(f"sign{name}", f"{wins}/{total}")
        macro(f"signPct{name}", f"{100 * wins / total:.0f}\\%")

    rank_biserial, p_value = domain_effect(rung3)
    macro("domainRankBiserial", f"{rank_biserial:+.2f}")
    macro("domainP", f"{p_value:.3f}")

    macro("familiesScored", str(mut["ernierna"]["n_scored"]))
    macro("expectedFalsePositives", f"{0.05 * mut['ernierna']['n_scored']:.1f}")

    # Both readings of the ambiguous "mean trained/untrained ratio", so the
    # choice the table makes is visible rather than buried in this script.
    for key in ["rnafm", "rinalmo", "splicebert"]:
        per_family, of_means = ratio_readings(mut, key)
        macro(f"tuPerFamily{BY_KEY[key][3]}", f"{per_family:.2f}")
        macro(f"tuOfMeans{BY_KEY[key][3]}", f"{of_means:.2f}")

    # The transversion control, against the Watson-Crick run family for family.
    alphabet = next(iter(controls.values()))["alphabet"]
    macro("transversionAlphabet",
          ", ".join(f"{nuc}$\\to${sub}" for nuc, sub in sorted(alphabet.items())
                    if nuc != "T"))
    deviations = {}
    for key, control in controls.items():
        name = BY_KEY[key][3]
        ratio = transversion_ratio(control)
        macro(f"trans{name}", f"{ratio:.3f}")
        deviation = 100 * (ratio - mut[key]["mean_ratio"]) / mut[key]["mean_ratio"]
        macro(f"transDelta{name}", f"{deviation:+.1f}\\%")
        deviations[key] = abs(deviation)
    widest = max(deviations, key=lambda key: deviations[key])
    macro("transversionWidest", BY_KEY[widest][0])
    macro("transversionWidestDeviation", f"{deviations[widest]:.1f}\\%")
    macro("transversionBound", f"{math.ceil(deviations[widest]):d}\\%")

    # The attention bound, and the widest gap it covers. An absolute delta says
    # nothing without the baseline it moved from, so the factor goes out too.
    macro("attnBound", f"{ATTENTION_ARCHITECTURAL_DELTA:g}")
    gaps = {}
    for key in UNTRAINED_KEYS:
        trained, untrained = attention_mean(runs[key]), attention_mean(runs[f"{key}_untrained"])
        if trained is not None and untrained is not None and untrained > 0:
            gaps[key] = (abs(trained - untrained), trained / untrained)
    attn_widest = max(gaps, key=lambda key: gaps[key][0])
    macro("attnWidest", BY_KEY[attn_widest][0])
    macro("attnWidestFactor", f"{gaps[attn_widest][1]:.1f}")

    # The separation the overview figure quotes, rather than a round number
    # chosen once and left to drift.
    ps_rank = sorted((rung3[key]["mean_ps"] or 0.0 for key, *_ in MODELS), reverse=True)
    macro("psSeparation", f"{ps_rank[1] / ps_rank[2]:.0f}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------


def main() -> int:
    runs = load_all()
    controls = load_transversion()
    rung3 = {key: rung3_stats(run) for key, run in runs.items()}
    mut = mutation_all(runs)
    check_claims(runs, rung3, mut, controls)
    GENERATED.mkdir(parents=True, exist_ok=True)
    written = {
        "rung1_table.tex": rung1_table(runs, mut),
        "rung2_table.tex": rung2_table(mut),
        "rung3_table.tex": rung3_table(runs, rung3),
        "attention_table.tex": attention_table(runs),
        "provenance_table.tex": provenance_table(runs),
        "hypotheses_table.tex": hypotheses_table(runs, rung3, mut),
        "transversion_table.tex": transversion_table(mut, controls),
        "results_macros.tex": macros(runs, rung3, mut, controls),
    }
    for name, text in written.items():
        (GENERATED / name).write_text(text)
        print(f"  wrote paper/generated/{name}  ({len(text.splitlines())} lines)")

    print("\nRegistered hypotheses:")
    for tag, criterion, panel, result, mark in hypothesis_rows(runs, rung3, mut):
        plain = mark.replace(r"\textbf{", "").replace("}", "")
        print(f"  {tag:<8} {plain:<5} {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
