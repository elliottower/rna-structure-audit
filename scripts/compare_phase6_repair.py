"""What the token-span repair changed in Rung 3, model by model.

D13 had `compute_delta_profiles` read `emb[k]` for every model, so on a
subword tokenizer a nucleotide index addressed a token array. A pair whose
index ran past the token array was dropped, and a family that lost every pair
was written out as `"no valid PS values"` -- the same field a registered
filter writes, which is why the loss was invisible in the deposited results.

This compares a repaired run against the deposited one on the quantities the
registered hypotheses read: how many families were scored at all, the mean PS
over non-quarantined families that H1 and H2 use, the gate count H1 tests, and
the per-pair precision H3 tests.

    uv run --no-project --python 3.12 python scripts/compare_phase6_repair.py \
        --old results/repaired_panel --new <fetched dir> [--models nt,dnabert2]

Both arguments are directories laid out as `<key>/<key>_phase6_ps.json`. A
single fetched file may be passed to `--new` instead when one model is being
checked before the rest of a run lands.
"""

import argparse
import json
from pathlib import Path
from statistics import mean

REPO = Path(__file__).resolve().parents[1]

# Registered in PREREGISTRATION_PHASE6_V2.md. Both are tRNA families whose
# published structures the panel's dot-bracket annotation does not reproduce,
# so they are excluded from every aggregate a hypothesis reads.
QUARANTINED = {"tRNA_Phe_yeast", "tRNA_Ala_human"}


def load(path):
    """The per-family block of a phase 6 result file."""
    return json.loads(Path(path).read_text())["results"]["per_rna"]


def summarize(per_rna):
    """The four quantities H1, H2 and H3 are defined on.

    `scored` counts families carrying a PS value at all, which is the count D13
    moved. The mean is over non-quarantined scored families rather than the
    stored `mean_best_ps`, which is taken over families that passed the
    positive-control gate and so answers a different question.
    """
    kept = {name: entry for name, entry in per_rna.items()
            if name not in QUARANTINED and entry.get("best_ps") is not None}
    precisions = [entry["h3_precision"]["fraction"] for entry in kept.values()
                  if entry.get("h3_precision")]
    return {
        "scored": len(kept),
        "mean_ps": mean(e["best_ps"] for e in kept.values()) if kept else None,
        "exceeding_null": sum(1 for e in kept.values()
                              if e.get("exceeds_null_primary")),
        "precision": mean(precisions) if precisions else None,
    }


def fmt(value):
    return "--" if value is None else f"{value:+.4f}" if isinstance(value, float) else str(value)


def row(model, old, new):
    print(f"\n  {model}")
    print(f"    {'':22s} {'deposited':>12s} {'repaired':>12s}")
    for label, key in [("families scored", "scored"),
                       ("mean PS", "mean_ps"),
                       ("families over null", "exceeding_null"),
                       ("per-pair precision", "precision")]:
        print(f"    {label:22s} {fmt(old[key]):>12s} {fmt(new[key]):>12s}")


def per_family_diff(old, new):
    """Families whose stored block differs, and which field carries it.

    A model the repair does not touch should differ nowhere. One that differs
    in a single field while its deltas match is the case worth reading: the
    embeddings are the same and something downstream of them is not.
    """
    for name in sorted(set(old) | set(new)):
        before, after = old.get(name), new.get(name)
        if before == after:
            continue
        if before is None or after is None:
            print(f"      {name}: present in only one run")
            continue
        for field in sorted(set(before) | set(after)):
            if before.get(field) != after.get(field):
                print(f"      {name}.{field}: "
                      f"{before.get(field)!r} -> {after.get(field)!r}")


def resolve(directory, model, explicit_file):
    if explicit_file is not None:
        return explicit_file
    return Path(directory) / model / f"{model}_phase6_ps.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--old", default=str(REPO / "results/repaired_panel"))
    parser.add_argument("--new", required=True)
    parser.add_argument("--models", default="nt,dnabert2")
    parser.add_argument("--diff", action="store_true",
                        help="print every per-family field that differs")
    args = parser.parse_args()

    new_path = Path(args.new)
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    single = new_path if new_path.is_file() else None
    if single is not None and len(models) != 1:
        raise SystemExit("--new names one file, so --models must name one model")

    print(f"Rung 3, deposited vs repaired. Quarantined and excluded: "
          f"{', '.join(sorted(QUARANTINED))}")
    for model in models:
        old_per = load(resolve(args.old, model, None))
        new_per = load(resolve(args.new, model, single))
        row(model, summarize(old_per), summarize(new_per))
        if args.diff:
            per_family_diff(old_per, new_per)
    print()


if __name__ == "__main__":
    main()
