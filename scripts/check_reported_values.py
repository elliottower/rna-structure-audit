"""Compare the values reported in the manuscript against the pipeline's output.

The manuscript's Table 5 reports a mean perturbation specificity that the
pipeline did not produce. run_phase6 computes `mean_best_ps` over families that
are not skipped, not quarantined, and that pass the positive-control gate. The
manuscript reports a mean over all non-skipped families, which the Methods
explicitly exclude.

This script reads the stored `mean_best_ps` from each result file and prints it
beside the manuscript value, so the discrepancy is visible without opening a
JSON by hand. Update REPORTED after regenerating the tables.

Usage:
    uv run python scripts/check_reported_values.py
    uv run python scripts/check_reported_values.py --tolerance 0.05 --strict
"""

import argparse
import glob
import json
import sys

# Values as printed in Table 5 of the manuscript.
REPORTED = {
    "rinalmo": 0.2150,
    "ernierna": 0.1204,
    "ernierna_untrained": 9.5e-8,
    "caduceus": 0.0034,
    "evo": 0.0011,
    "hyenadna": 0.0003,
    "splicebert": 0.0002,
    "rnafm": 0.0001,
    "utrlm": 0.00001,
    "nt": 0.001,
    "dnabert2": -0.016,
}


def collect(pattern):
    found = {}
    for path in sorted(glob.glob(pattern, recursive=True)):
        try:
            blob = json.load(open(path))
        except Exception:
            continue
        results = blob.get("results")
        if not isinstance(results, dict) or "per_rna" not in results:
            continue
        model = blob.get("model")
        if model and model not in found:
            found[model] = {
                "mean": results.get("mean_best_ps"),
                "n_active": results.get("families_total"),
                "exceeding": results.get("families_exceeding_null_conservative"),
                "conservative": results.get("families_exceeding_null_conservative"),
                "path": path,
            }
    return found


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", default="results/**/*phase6*.json")
    parser.add_argument("--tolerance", type=float, default=0.10,
                        help="relative difference tolerated before flagging")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    found = collect(args.pattern)
    header = f"{'model':<22}{'pipeline':>14}{'N':>4}{'>null':>7}{'>cons':>7}{'manuscript':>14}   status"
    print(header)
    print("-" * len(header))

    flagged = []
    for model, row in sorted(found.items(), key=lambda kv: -(kv[1]["mean"] or -1)):
        pipeline, reported = row["mean"], REPORTED.get(model)
        status = ""
        if pipeline is None or reported is None:
            status = "no comparison"
        elif pipeline == 0 and reported != 0:
            status = "pipeline is exactly zero"
            flagged.append(model)
        elif pipeline != 0:
            relative = abs(reported - pipeline) / abs(pipeline)
            if relative > args.tolerance:
                status = f"differs by {relative:.0f}x" if relative > 1 else f"differs {relative:.0%}"
                flagged.append(model)
        shown_pipeline = "None" if pipeline is None else f"{pipeline:.3e}"
        shown_reported = "-" if reported is None else f"{reported:.3e}"
        print(f"{model:<22}{shown_pipeline:>14}{str(row['n_active']):>4}"
              f"{str(row['exceeding']):>7}{str(row['conservative']):>7}"
              f"{shown_reported:>14}   {status}")

    print(f"\n  {len(found)} models compared, {len(flagged)} outside tolerance")
    if flagged:
        print(f"  flagged: {', '.join(sorted(flagged))}")
    if args.strict and flagged:
        sys.exit(1)


if __name__ == "__main__":
    main()
