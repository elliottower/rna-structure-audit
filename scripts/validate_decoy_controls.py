"""Validation gate for the decoy comparison, on the panel's own families.

Two delta profiles with known answers are pushed through the functions the
model runs use. The partner oracle moves only the true partner: its decoy rate
must be exactly 1 in every family. The complementarity profile moves every
position carrying the partner's nucleotide equally: its decoy rate must be
exactly 0.5 in every family, for every scope and decoy set. The gate fails on
any other value.

Also recorded, as a description of the panel: the Rung 3 per-pair precision the
complementarity profile reaches, which is the share of deep-interior pairs
whose partner differs in nucleotide from both of its stem neighbours.

    uv run --no-project python scripts/validate_decoy_controls.py
"""

import json
from pathlib import Path

from count_matched_decoys import confirmatory_families
from decoy_comparison import (
    DECOY_SETS,
    SCOPES,
    complementarity_profile,
    oracle_profile,
    partner_beats_decoy,
    precision_at_every_layer,
)
from stem_pairs import is_deep_interior

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "audit" / "decoy_control_validation.json"

EXPECTED = {"oracle": 1.0, "complementarity": 0.5}
PROFILES = {"oracle": oracle_profile, "complementarity": complementarity_profile}


class GateFailure(RuntimeError):
    """A control with a known answer returned another value."""


def main() -> int:
    failures = []
    per_family = {}
    hits = total = 0
    for name, (sequence, eligible) in confirmatory_families().items():
        row = {}
        for control, build in PROFILES.items():
            profiles = [{0: build(sequence, pair)} for pair in eligible]
            row[control] = {
                "rung3_precision": precision_at_every_layer(eligible, profiles, 1)[0]}
            for scope in SCOPES:
                for decoy_set in DECOY_SETS:
                    result = partner_beats_decoy(
                        sequence, eligible, profiles, 1, scope, decoy_set)
                    rate = result["rate"][0]
                    row[control][f"{scope}_{decoy_set}"] = rate
                    if rate is not None and rate != EXPECTED[control]:
                        failures.append((name, control, scope, decoy_set, rate))
        per_family[name] = row
        deep = sum(is_deep_interior(pair) for pair in eligible)
        precision = row["complementarity"]["rung3_precision"]
        if precision is not None:
            hits += round(precision * deep)
            total += deep

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "expected": EXPECTED,
        "passed": not failures,
        "failures": failures,
        "complementarity_rung3_precision_pooled": hits / total if total else None,
        "deep_interior_pairs": total,
        "families": per_family,
    }, indent=2))
    print(f"wrote {OUT.relative_to(REPO)}")
    if failures:
        raise GateFailure(f"{len(failures)} control values off target: {failures[:5]}")
    print("gate passed; complementarity profile Rung 3 precision, pooled:",
          hits / total if total else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
