"""Compute mutual information between paired positions in Rfam seed alignments.

Non-neural baseline for partner specificity (PS): if paired positions show
high MI in evolutionary alignments, the PS signal might come from sequence
covariation rather than learned structure.

For each Rfam family:
  1. Parse the Stockholm seed alignment to get aligned sequences and SS_cons.
  2. Identify paired columns from SS_cons.
  3. Compute MI(i, j) for each paired column pair across the alignment.
  4. Compute MI for an equal-sized sample of unpaired column pairs.
  5. Compare per-family mean MI with per-family PS from Phase 6 results.

Usage:
    cd /path/to/causal-rna
    uv run --with numpy --with scipy --with tqdm python scripts/compute_mi_baseline.py
"""

import json
import re
import string
from pathlib import Path

import numpy as np
from scipy import stats
from tqdm import tqdm

ROOT = Path(__file__).resolve().parent.parent
SEED_DIR = ROOT / "data" / "rfam_seeds"
FAMILY_DIR = ROOT / "data" / "rfam_families"
RESULTS_DIR = ROOT / "results"
OUTPUT_PATH = RESULTS_DIR / "mi_baseline.json"

# Phase 6 PS result files (most recent per model)
PS_RESULT_FILES = {
    "rinalmo": RESULTS_DIR / "rinalmo_phase6_ps.json",
    "ernierna": RESULTS_DIR / "audit" / "phase6_ernierna_20260714_060520" / "ernierna_phase6_ps.json",
}

NUCLEOTIDES = ["A", "C", "G", "U", "-"]
NUC_TO_IDX = {n: i for i, n in enumerate(NUCLEOTIDES)}
N_STATES = len(NUCLEOTIDES)
GAP_FRACTION_THRESHOLD = 0.5
MIN_SEQUENCES = 5

RNG = np.random.default_rng(42)

# Stockholm SS_cons bracket types: each opener maps to its closer
BRACKET_PAIRS = {"<": ">", "(": ")", "[": "]", "{": "}"}
# Pseudoknot levels use uppercase/lowercase letter pairs (A/a, B/b, ...)
PSEUDOKNOT_OPENERS = set(string.ascii_uppercase)
PSEUDOKNOT_CLOSERS = set(string.ascii_lowercase)
UNPAIRED_CHARS = set(":,._-~")


def load_family_rfam_map() -> dict[str, str]:
    """Map family name -> Rfam ID from family JSON files."""
    mapping = {}
    for f in sorted(FAMILY_DIR.iterdir()):
        if f.suffix == ".json":
            with open(f) as fh:
                data = json.load(fh)
            rfam_id = data.get("rfam_id")
            if rfam_id and rfam_id != "None":
                mapping[data["name"]] = rfam_id
    return mapping


def parse_stockholm_alignment(sto_path: Path) -> tuple[list[str], str]:
    """Parse Stockholm file, return (list of aligned sequences, SS_cons string).

    Sequences are uppercased with T->U. Gap characters (.-~) are normalized to '-'.
    """
    sequences: dict[str, str] = {}
    ss_cons = ""

    with open(sto_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#=GC SS_cons"):
                ss_cons += line.split(None, 2)[2]
            elif line.startswith("#") or line.startswith("//") or not line.strip():
                continue
            else:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    name, seq = parts
                    if name in sequences:
                        sequences[name] += seq
                    else:
                        sequences[name] = seq

    # Normalize sequences: uppercase, T->U, gaps to '-'
    normalized = []
    for seq in sequences.values():
        s = seq.upper().replace("T", "U")
        s = re.sub(r"[.~]", "-", s)
        normalized.append(s)

    return normalized, ss_cons


def parse_pairs_from_ss_cons(ss_cons: str) -> list[tuple[int, int]]:
    """Extract base pair positions (i, j) from SS_cons string.

    Handles standard brackets (<>, (), [], {}) and pseudoknot annotations
    (uppercase/lowercase letter pairs like A/a, B/b).
    Returns list of (i, j) tuples where i < j.
    """
    pairs = []

    # Standard bracket types
    for opener, closer in BRACKET_PAIRS.items():
        stack = []
        for col, ch in enumerate(ss_cons):
            if ch == opener:
                stack.append(col)
            elif ch == closer:
                if stack:
                    i = stack.pop()
                    pairs.append((i, col))

    # Pseudoknot letter pairs (A/a, B/b, ...)
    letter_stacks: dict[str, list[int]] = {}
    for col, ch in enumerate(ss_cons):
        if ch in PSEUDOKNOT_OPENERS:
            letter_stacks.setdefault(ch, []).append(col)
        elif ch in PSEUDOKNOT_CLOSERS:
            upper = ch.upper()
            if upper in letter_stacks and letter_stacks[upper]:
                i = letter_stacks[upper].pop()
                pairs.append((i, col))

    # Ensure i < j
    pairs = [(min(i, j), max(i, j)) for i, j in pairs]
    return sorted(pairs)


def get_unpaired_columns(ss_cons: str, paired_cols: set[int]) -> list[int]:
    """Return column indices that are unpaired (not in any base pair)."""
    unpaired = []
    for col, ch in enumerate(ss_cons):
        if col not in paired_cols:
            unpaired.append(col)
    return unpaired


def column_vectors(sequences: list[str], col: int) -> np.ndarray:
    """Extract nucleotide indices at a given column across all sequences.

    Returns array of shape (n_seqs,) with values in [0, N_STATES-1].
    Unknown characters map to gap index.
    """
    indices = np.empty(len(sequences), dtype=np.int32)
    for i, seq in enumerate(sequences):
        ch = seq[col] if col < len(seq) else "-"
        indices[i] = NUC_TO_IDX.get(ch, NUC_TO_IDX["-"])
    return indices


def compute_mi(col_i: np.ndarray, col_j: np.ndarray) -> float:
    """Compute mutual information between two column vectors.

    Uses base-4 log so MI is in bits of nucleotide information.
    Treats gaps as a 5th state.
    """
    n = len(col_i)
    if n == 0:
        return 0.0

    # Joint counts
    joint = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    for a, b in zip(col_i, col_j):
        joint[a, b] += 1.0
    joint /= n

    # Marginals
    p_i = joint.sum(axis=1)
    p_j = joint.sum(axis=0)

    mi = 0.0
    for a in range(N_STATES):
        for b in range(N_STATES):
            if joint[a, b] > 0 and p_i[a] > 0 and p_j[b] > 0:
                mi += joint[a, b] * np.log(joint[a, b] / (p_i[a] * p_j[b]))

    # Convert from nats to base-4 bits
    mi /= np.log(4.0)
    return float(mi)


def compute_normalized_mi(col_i: np.ndarray, col_j: np.ndarray) -> float:
    """Compute normalized mutual information NMI = MI / log(min(|X|, |Y|)).

    Normalized to [0, 1]. Uses base-4 log throughout.
    """
    mi = compute_mi(col_i, col_j)
    if mi <= 0.0:
        return 0.0

    # Entropy of each column
    n = len(col_i)
    counts_i = np.bincount(col_i, minlength=N_STATES).astype(np.float64)
    counts_j = np.bincount(col_j, minlength=N_STATES).astype(np.float64)
    p_i = counts_i / n
    p_j = counts_j / n

    h_i = -np.sum(p_i[p_i > 0] * np.log(p_i[p_i > 0])) / np.log(4.0)
    h_j = -np.sum(p_j[p_j > 0] * np.log(p_j[p_j > 0])) / np.log(4.0)

    denom = min(h_i, h_j)
    if denom <= 0.0:
        return 0.0

    return float(min(mi / denom, 1.0))


def gap_fraction(sequences: list[str], col: int) -> float:
    """Fraction of sequences with a gap at this column."""
    n_gaps = sum(1 for seq in sequences if col >= len(seq) or seq[col] == "-")
    return n_gaps / len(sequences)


def process_family(sto_path: Path) -> dict | None:
    """Compute MI statistics for one Rfam family.

    Returns dict with paired_mi, unpaired_mi, n_pairs, n_sequences, etc.
    or None if the family cannot be processed.
    """
    sequences, ss_cons = parse_stockholm_alignment(sto_path)

    if len(sequences) < MIN_SEQUENCES:
        return None

    alignment_len = len(ss_cons)
    # Verify all sequences match alignment length
    sequences = [s for s in sequences if len(s) == alignment_len]
    if len(sequences) < MIN_SEQUENCES:
        return None

    # Find valid columns (gap fraction below threshold)
    valid_cols = set()
    for col in range(alignment_len):
        if gap_fraction(sequences, col) <= GAP_FRACTION_THRESHOLD:
            valid_cols.add(col)

    # Parse paired positions
    all_pairs = parse_pairs_from_ss_cons(ss_cons)

    # Filter to pairs where both columns are valid
    paired_positions = [(i, j) for i, j in all_pairs if i in valid_cols and j in valid_cols]

    if len(paired_positions) < 3:
        return None

    # Compute MI for paired positions
    paired_mi_values = []
    paired_nmi_values = []
    for i, j in paired_positions:
        ci = column_vectors(sequences, i)
        cj = column_vectors(sequences, j)
        paired_mi_values.append(compute_mi(ci, cj))
        paired_nmi_values.append(compute_normalized_mi(ci, cj))

    # Get unpaired columns
    paired_cols = set()
    for i, j in all_pairs:
        paired_cols.add(i)
        paired_cols.add(j)
    unpaired_cols = [c for c in valid_cols if c not in paired_cols]

    # Sample unpaired pairs (same count as paired)
    n_sample = len(paired_positions)
    unpaired_mi_values = []
    unpaired_nmi_values = []

    if len(unpaired_cols) >= 2:
        n_possible = len(unpaired_cols) * (len(unpaired_cols) - 1) // 2
        n_sample = min(n_sample, n_possible)

        # Generate random unpaired pairs
        sampled = set()
        attempts = 0
        while len(sampled) < n_sample and attempts < n_sample * 10:
            idx = RNG.choice(len(unpaired_cols), size=2, replace=False)
            pair = (min(unpaired_cols[idx[0]], unpaired_cols[idx[1]]),
                    max(unpaired_cols[idx[0]], unpaired_cols[idx[1]]))
            sampled.add(pair)
            attempts += 1

        for i, j in sampled:
            ci = column_vectors(sequences, i)
            cj = column_vectors(sequences, j)
            unpaired_mi_values.append(compute_mi(ci, cj))
            unpaired_nmi_values.append(compute_normalized_mi(ci, cj))

    return {
        "n_sequences": len(sequences),
        "alignment_length": alignment_len,
        "n_valid_columns": len(valid_cols),
        "n_paired_positions": len(paired_positions),
        "n_unpaired_sampled": len(unpaired_mi_values),
        "paired_mi_mean": float(np.mean(paired_mi_values)),
        "paired_mi_std": float(np.std(paired_mi_values)),
        "paired_mi_median": float(np.median(paired_mi_values)),
        "paired_nmi_mean": float(np.mean(paired_nmi_values)),
        "paired_nmi_std": float(np.std(paired_nmi_values)),
        "unpaired_mi_mean": float(np.mean(unpaired_mi_values)) if unpaired_mi_values else None,
        "unpaired_mi_std": float(np.std(unpaired_mi_values)) if unpaired_mi_values else None,
        "unpaired_nmi_mean": float(np.mean(unpaired_nmi_values)) if unpaired_nmi_values else None,
        "unpaired_nmi_std": float(np.std(unpaired_nmi_values)) if unpaired_nmi_values else None,
        "mi_ratio": (float(np.mean(paired_mi_values)) / float(np.mean(unpaired_mi_values))
                     if unpaired_mi_values and np.mean(unpaired_mi_values) > 0 else None),
    }


def load_ps_results(path: Path) -> dict[str, float]:
    """Load per-family best_ps from a Phase 6 result file.

    Returns {family_name: best_ps} for non-skipped families.
    """
    with open(path) as f:
        data = json.load(f)

    ps = {}
    per_rna = data.get("results", {}).get("per_rna", {})
    for name, v in per_rna.items():
        if v.get("skipped"):
            continue
        if "best_ps" in v:
            ps[name] = v["best_ps"]
    return ps


def main():
    # Load family -> Rfam ID mapping
    family_rfam = load_family_rfam_map()
    rfam_to_family = {v: k for k, v in family_rfam.items()}

    print(f"Loaded {len(family_rfam)} families with Rfam IDs")

    # Process each family
    mi_results = {}
    for family_name in tqdm(sorted(family_rfam.keys()), desc="Computing MI"):
        rfam_id = family_rfam[family_name]
        sto_path = SEED_DIR / f"{rfam_id}.sto"

        if not sto_path.exists():
            print(f"  {family_name}: no seed alignment at {sto_path}")
            continue

        result = process_family(sto_path)
        if result is None:
            print(f"  {family_name} ({rfam_id}): skipped (too few sequences or pairs)")
            continue

        mi_results[family_name] = result

    print(f"\nProcessed {len(mi_results)} families")

    # Summary statistics across all families
    all_paired_mi = [r["paired_mi_mean"] for r in mi_results.values()]
    all_unpaired_mi = [r["unpaired_mi_mean"] for r in mi_results.values()
                       if r["unpaired_mi_mean"] is not None]
    all_paired_nmi = [r["paired_nmi_mean"] for r in mi_results.values()]
    all_unpaired_nmi = [r["unpaired_nmi_mean"] for r in mi_results.values()
                        if r["unpaired_nmi_mean"] is not None]

    print(f"\nMI summary (base-4 log):")
    print(f"  Paired positions:   mean={np.mean(all_paired_mi):.4f}  std={np.std(all_paired_mi):.4f}")
    print(f"  Unpaired positions: mean={np.mean(all_unpaired_mi):.4f}  std={np.std(all_unpaired_mi):.4f}")
    print(f"  Wilcoxon paired>unpaired: ", end="")
    # Use families that have both paired and unpaired
    matched_paired = []
    matched_unpaired = []
    for r in mi_results.values():
        if r["unpaired_mi_mean"] is not None:
            matched_paired.append(r["paired_mi_mean"])
            matched_unpaired.append(r["unpaired_mi_mean"])
    if len(matched_paired) >= 5:
        stat, pval = stats.wilcoxon(matched_paired, matched_unpaired, alternative="greater")
        print(f"W={stat:.1f}, p={pval:.2e}")
    else:
        stat, pval = None, None
        print("too few families")

    print(f"\nNMI summary:")
    print(f"  Paired positions:   mean={np.mean(all_paired_nmi):.4f}  std={np.std(all_paired_nmi):.4f}")
    print(f"  Unpaired positions: mean={np.mean(all_unpaired_nmi):.4f}  std={np.std(all_unpaired_nmi):.4f}")

    # Compare MI with PS from Phase 6 results
    correlations = {}
    for model_name, ps_path in PS_RESULT_FILES.items():
        if not ps_path.exists():
            print(f"\n{model_name}: PS result file not found at {ps_path}")
            continue

        ps_values = load_ps_results(ps_path)
        print(f"\n{model_name}: loaded PS for {len(ps_values)} families")

        # Match families present in both MI and PS
        common_families = sorted(set(mi_results.keys()) & set(ps_values.keys()))
        print(f"  Common families: {len(common_families)}")

        if len(common_families) < 5:
            print("  Too few common families for correlation")
            continue

        mi_vec = [mi_results[f]["paired_mi_mean"] for f in common_families]
        nmi_vec = [mi_results[f]["paired_nmi_mean"] for f in common_families]
        ps_vec = [ps_values[f] for f in common_families]

        rho_mi, p_mi = stats.spearmanr(mi_vec, ps_vec)
        rho_nmi, p_nmi = stats.spearmanr(nmi_vec, ps_vec)

        print(f"  Spearman(MI, PS):  rho={rho_mi:.3f}, p={p_mi:.3e}")
        print(f"  Spearman(NMI, PS): rho={rho_nmi:.3f}, p={p_nmi:.3e}")

        correlations[model_name] = {
            "n_families": len(common_families),
            "families": common_families,
            "spearman_mi_ps": {"rho": float(rho_mi), "p_value": float(p_mi)},
            "spearman_nmi_ps": {"rho": float(rho_nmi), "p_value": float(p_nmi)},
            "per_family": {
                f: {"mi": mi_vec[i], "nmi": nmi_vec[i], "ps": ps_vec[i]}
                for i, f in enumerate(common_families)
            },
        }

    # Build output
    output = {
        "description": "Mutual information baseline for partner specificity",
        "method": {
            "log_base": 4,
            "states": NUCLEOTIDES,
            "gap_fraction_threshold": GAP_FRACTION_THRESHOLD,
            "min_sequences": MIN_SEQUENCES,
            "normalization": "MI / min(H(X), H(Y))",
        },
        "summary": {
            "n_families_processed": len(mi_results),
            "paired_mi_mean": float(np.mean(all_paired_mi)),
            "paired_mi_std": float(np.std(all_paired_mi)),
            "unpaired_mi_mean": float(np.mean(all_unpaired_mi)),
            "unpaired_mi_std": float(np.std(all_unpaired_mi)),
            "paired_nmi_mean": float(np.mean(all_paired_nmi)),
            "unpaired_nmi_mean": float(np.mean(all_unpaired_nmi)),
            "wilcoxon_paired_gt_unpaired": {
                "statistic": float(stat) if stat is not None else None,
                "p_value": float(pval) if pval is not None else None,
            },
        },
        "per_family": mi_results,
        "correlations_with_ps": correlations,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
