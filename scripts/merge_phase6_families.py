"""Merge recomputed families into an existing Phase 6 result file.

Only four families changed when the corrupted miRNA precursor files were
replaced (let-7 RF00027, mir-122 RF00684, mir-155 RF00731, mir-21 RF00658).
Per-family Phase 6 computation is independent -- run_phase6 loops over families
with no cross-family state -- so those four can be recomputed alone and merged
into the stored results for the other 48.

The aggregate fields are recomputed here. They must match run_phase6 exactly,
so `--verify` recomputes them from every stored file's own per_rna and compares
against the stored values. Run --verify and get a clean pass BEFORE trusting a
merge.

Usage:
    uv run python scripts/merge_phase6_families.py --verify
    uv run python scripts/merge_phase6_families.py \
        --base results/rinalmo_phase6_ps.json \
        --patch results/rerun/rinalmo_phase6_ps.json \
        --out   results/merged/rinalmo_phase6_ps.json
"""

import argparse
import json
import pathlib
import statistics

RECOMPUTED_FAMILIES = [
    "mir_let7_precursor",
    "mir_122_precursor",
    "mir_155_precursor",
    "mir_21_precursor",
]

# Fields run_phase6 derives from per_rna. Mirrors phase6_compensatory_mutation.py.
AGGREGATE_FIELDS = [
    "mean_best_ps",
    "median_best_ps",
    "families_exceeding_null_primary",
    "families_exceeding_null_conservative",
    "families_no_null",
    "families_total",
    "families_quarantined",
    "families_skipped",
    "families_failed_gate",
]


def compute_aggregates(per_rna, compute_null=True):
    """Reproduce the aggregate block at the end of run_phase6.

    active = not skipped, not quarantined, and positive_control passed.
    """
    out = {}
    active = {
        k: v for k, v in per_rna.items()
        if not v.get("skipped")
        and not v.get("quarantined")
        and v.get("positive_control", {}).get("pass", False)
    }
    if active:
        ps_values = [v["best_ps"] for v in active.values()]
        out["mean_best_ps"] = float(statistics.fmean(ps_values))
        out["median_best_ps"] = float(statistics.median(ps_values))
        if compute_null:
            out["families_exceeding_null_primary"] = sum(
                1 for v in active.values() if v.get("exceeds_null_primary") is True)
            out["families_exceeding_null_conservative"] = sum(
                1 for v in active.values() if v.get("exceeds_null_conservative") is True)
            out["families_no_null"] = [
                k for k, v in active.items() if v.get("null_available") is False]
    out["families_total"] = len(active)
    out["families_quarantined"] = [k for k, v in per_rna.items() if v.get("quarantined")]
    out["families_skipped"] = [k for k, v in per_rna.items() if v.get("skipped")]
    out["families_failed_gate"] = [
        k for k, v in per_rna.items()
        if not v.get("skipped")
        and not v.get("quarantined")
        and not v.get("positive_control", {}).get("pass", False)
    ]
    return out


def close(a, b, tol=1e-9):
    if isinstance(a, float) and isinstance(b, float):
        return abs(a - b) <= tol * max(1.0, abs(a), abs(b))
    if isinstance(a, list) and isinstance(b, list):
        return sorted(a) == sorted(b)
    return a == b


def verify(root):
    """Recompute aggregates from each stored file's own per_rna and compare."""
    files = sorted(pathlib.Path(root).rglob("*phase6*.json"))
    checked = failed = 0
    for f in files:
        try:
            d = json.load(open(f))
        except Exception:
            continue
        res = d.get("results")
        if not isinstance(res, dict) or "per_rna" not in res:
            continue
        checked += 1
        got = compute_aggregates(res["per_rna"], compute_null="families_exceeding_null_primary" in res)
        diffs = []
        for k in AGGREGATE_FIELDS:
            if k not in res and k not in got:
                continue
            if not close(got.get(k), res.get(k)):
                diffs.append(f"{k}: stored={res.get(k)!r} recomputed={got.get(k)!r}")
        rel = f.relative_to(root)
        if diffs:
            failed += 1
            print(f"  MISMATCH {rel}")
            for x in diffs:
                print(f"      {x}")
        else:
            print(f"  ok       {rel}  (mean_best_ps={res.get('mean_best_ps')})")
    print(f"\n  {checked} files checked, {failed} mismatched")
    return failed == 0


def merge(base_path, patch_path, out_path):
    base = json.load(open(base_path))
    patch = json.load(open(patch_path))
    bp = base["results"]["per_rna"]
    pp = patch["results"]["per_rna"]

    replaced, missing = [], []
    for fam in RECOMPUTED_FAMILIES:
        if fam in pp:
            bp[fam] = pp[fam]
            replaced.append(fam)
        else:
            missing.append(fam)
    if missing:
        raise SystemExit(f"patch is missing recomputed families: {missing}")

    unexpected = set(pp) - set(RECOMPUTED_FAMILIES)
    if unexpected:
        raise SystemExit(f"patch contains families outside the recomputed set: {sorted(unexpected)}")

    agg = compute_aggregates(bp, compute_null="families_exceeding_null_primary" in base["results"])
    for k in AGGREGATE_FIELDS:
        base["results"].pop(k, None)
    base["results"].update(agg)
    base["merged_from"] = {
        "base": str(base_path), "patch": str(patch_path),
        "families_replaced": replaced,
    }
    pathlib.Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    json.dump(base, open(out_path, "w"), indent=2)
    print(f"  replaced {len(replaced)} families -> {out_path}")
    print(f"  mean_best_ps {json.load(open(base_path))['results'].get('mean_best_ps')} -> {agg.get('mean_best_ps')}")
    print(f"  families_total {json.load(open(base_path))['results'].get('families_total')} -> {agg.get('families_total')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true",
                    help="recompute aggregates from stored per_rna and compare to stored values")
    ap.add_argument("--root", default="results")
    ap.add_argument("--base")
    ap.add_argument("--patch")
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.verify:
        raise SystemExit(0 if verify(a.root) else 1)
    if not (a.base and a.patch and a.out):
        raise SystemExit("need --base, --patch and --out (or --verify)")
    merge(a.base, a.patch, a.out)


if __name__ == "__main__":
    main()
