"""Validation gate for the decoy comparison, on the panel's own families.

Supersedes `validate_decoy_controls.py`, which accepted a rate of None in any
family and so could not catch an implementation that returned None everywhere.

Two delta profiles with known answers are pushed through the functions the
model runs use. For every family, scope and decoy set, the number of
contributing pairs must equal the structure-only count of
`results/audit/matched_decoy_counts.json`. Where that count is above zero the
partner oracle must score exactly 1 and the complementarity control exactly
0.5; where it is zero the condition is recorded as unavailable and the rate
must be None. Any other outcome fails the gate.

Also recorded, as a description of the panel: the Rung 3 per-pair precision of
the complementarity control, counted pair by pair.

    uv run --no-project python scripts/validate_decoy_controls_v2.py
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
)
from stem_pairs import is_deep_interior

REPO = Path(__file__).resolve().parents[1]
COUNTS = REPO / "results" / "audit" / "matched_decoy_counts.json"
OUT = REPO / "results" / "audit" / "decoy_control_validation_v2.json"

EXPECTED = {"oracle": 1.0, "complementarity": 0.5}
PROFILES = {"oracle": oracle_profile, "complementarity": complementarity_profile}


class GateFailure(RuntimeError):
    """A control with a known answer returned another value."""


def rung3_hits(eligible, profiles) -> tuple[int, int]:
    """(pairs whose partner beats both stem neighbours, deep-interior pairs)."""
    hits = total = 0
    for pair, profile in zip(eligible, profiles):
        if not is_deep_interior(pair):
            continue
        deltas = profile[0]
        total += 1
        hits += deltas[pair["j"]] > max(deltas[pair["j_prev"]], deltas[pair["j_next"]])
    return hits, total


def main() -> int:
    counts = json.loads(COUNTS.read_text())["families"]
    families = confirmatory_families()
    if set(families) != set(counts):
        raise GateFailure("the count file and the panel name different families")

    failures, per_family = [], {}
    available = unavailable = 0
    control_hits = control_total = 0
    for name, (sequence, eligible) in families.items():
        row = {}
        for control, build in PROFILES.items():
            profiles = [{0: build(sequence, pair)} for pair in eligible]
            hits, total = rung3_hits(eligible, profiles)
            row[control] = {"rung3_hits": hits, "rung3_pairs": total}
            if control == "complementarity":
                control_hits += hits
                control_total += total
            for scope in SCOPES:
                for decoy_set in DECOY_SETS:
                    tag = f"{scope}_{decoy_set}"
                    expected_pairs = counts[name][tag]["pairs"]
                    result = partner_beats_decoy(
                        sequence, eligible, profiles, 1, scope, decoy_set)
                    rate, pairs = result["rate"][0], result["pairs"][0]
                    row[control][tag] = {"rate": rate, "pairs": pairs}
                    if pairs != expected_pairs:
                        failures.append((name, control, tag, "pairs", pairs, expected_pairs))
                    elif expected_pairs == 0:
                        unavailable += 1
                        if rate is not None:
                            failures.append((name, control, tag, "rate", rate, None))
                    else:
                        available += 1
                        if rate != EXPECTED[control]:
                            failures.append(
                                (name, control, tag, "rate", rate, EXPECTED[control]))
        per_family[name] = row

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "expected": EXPECTED,
        "counts": str(COUNTS.relative_to(REPO)),
        "passed": not failures,
        "failures": failures,
        "conditions_evaluated": available,
        "conditions_unavailable": unavailable,
        "complementarity_control_rung3": {
            "hits": control_hits, "pairs": control_total,
            "precision_pooled": control_hits / control_total},
        "families": per_family,
    }, indent=2))
    print(f"wrote {OUT.relative_to(REPO)}")
    if failures:
        raise GateFailure(f"{len(failures)} control values off target: {failures[:5]}")
    print(f"gate passed: {available} conditions evaluated, {unavailable} unavailable; "
          f"complementarity control Rung 3: {control_hits}/{control_total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
