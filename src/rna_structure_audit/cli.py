"""Command-line interface for rna-structure-audit."""

import argparse
import importlib.util
import json
import sys
from pathlib import Path


def _load_adapter_from_file(path):
    spec = importlib.util.spec_from_file_location("user_adapter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for attr_name in dir(module):
        attr = getattr(module, attr_name)
        if (isinstance(attr, type)
                and hasattr(attr, "load")
                and hasattr(attr, "tokenize")
                and hasattr(attr, "get_all_layer_embeddings")
                and attr_name != "ModelAdapter"):
            return attr()
    raise ValueError(f"No ModelAdapter subclass found in {path}")


def main():
    parser = argparse.ArgumentParser(
        prog="rna-structure-audit",
        description="Three-rung benchmark for RNA/DNA foundation model structure awareness.",
    )
    parser.add_argument(
        "--adapter", required=True,
        help="Path to a .py file containing a ModelAdapter subclass.",
    )
    parser.add_argument(
        "--device", default="cpu",
        help="Torch device (default: cpu).",
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="Path to save full JSON results.",
    )
    parser.add_argument(
        "--rungs", default="1,2,3",
        help="Comma-separated rungs to run (default: 1,2,3).",
    )
    parser.add_argument(
        "--n-permutations", type=int, default=1000,
        help="Permutations for Rung 1/2 null (default: 1000).",
    )
    parser.add_argument(
        "--n-derangements", type=int, default=1000,
        help="Derangements for Rung 3 null (default: 1000).",
    )

    args = parser.parse_args()

    adapter = _load_adapter_from_file(args.adapter)
    rungs = tuple(int(r) for r in args.rungs.split(","))

    from rna_structure_audit.evaluate import evaluate

    results = evaluate(
        adapter,
        device=args.device,
        n_permutations=args.n_permutations,
        n_derangements=args.n_derangements,
        output_path=args.output,
        rungs=rungs,
    )

    report = results["report"]
    print()
    print(f"  RNA Structure Audit: {report['model']}")
    print(f"  {'=' * 45}")

    if "rung1" in report:
        r1 = report["rung1"]
        status = "PASS" if r1["families_exceeding_null"] > 5 else "MARGINAL" if r1["families_exceeding_null"] > 2 else "FAIL"
        print(f"  Rung 1: ratio {r1['mean_ratio']:.3f}, "
              f"{r1['families_exceeding_null']}/{r1['families_scored']} fam > null  "
              f"[{status}]")

    if "rung2" in report:
        r2 = report["rung2"]
        status = "PASS" if r2["surviving_dinuc"] > 5 else "FAIL"
        print(f"  Rung 2: {r2['surviving_dinuc']}/{r2['rung1_passing']} survive dinuc null  "
              f"[{status}]")

    if "rung3" in report:
        r3 = report["rung3"]
        status = "PASS" if r3["mean_ps"] > 0.05 and r3["mean_h3_precision"] > 0.5 else "FAIL"
        print(f"  Rung 3: PS = {r3['mean_ps']:.4f}, "
              f"precision {r3['mean_h3_precision']:.1%}  "
              f"[{status}]")

    if "grade" in report:
        print(f"  {'─' * 45}")
        print(f"  Grade: {report['grade']} — {report['grade_description']}")

    print()

    if args.output:
        print(f"  Full results saved to {args.output}")


if __name__ == "__main__":
    main()
