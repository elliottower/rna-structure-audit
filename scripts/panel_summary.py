"""Print the repaired panel's numbers, one line per model.

`generate_results_tables.py` refuses to write a table when the panel reverses a
sentence in the body, which leaves the prose to be rewritten with no tables to
read it from. This prints the same quantities the generator computes, through
the same loaders, so a sentence gets rewritten against the run rather than
against a recollection of it.
"""

from __future__ import annotations

import generate_results_tables as g

HEADINGS = ("model", "ratio", "95% CI", ">null", "dinuc", "attn", "probe",
            "PS", "gate", ">null", "H3", "trans")
WIDTHS = (22, 8, 16, 7, 8, 8, 8, 12, 8, 7, 8, 8)


def row(cells: tuple) -> str:
    return "".join(f"{cell:>{width}}" for cell, width in zip(cells, WIDTHS))


def main() -> int:
    runs = g.load_all()
    controls = g.load_transversion()
    mut = g.mutation_all(runs)
    rung3 = {key: g.rung3_stats(run) for key, run in runs.items()}

    print(row(HEADINGS))
    for key, short, _size, domain, _macro in g.MODELS:
        stats_, three = mut[key], rung3[key]
        interval = stats_["ci"]
        attention = g.attention_mean(runs[key])
        probing = g.probing_accuracy(runs[key])
        print(row((
            f"{short} ({domain})",
            f"{stats_['mean_ratio']:.3f}",
            f"[{interval['ci_lower']:.2f}, {interval['ci_upper']:.2f}]",
            stats_["exceeds_nuc"],
            f"{stats_['survives_dinuc']}/{stats_['exceeds_nuc']}",
            "---" if attention is None else f"{attention:.3f}",
            "---" if probing is None else f"{probing:.3f}",
            f"{three['mean_ps']:.2e}" if three["mean_ps"] is not None else "---",
            f"{three['gate']}/{three['eligible']}",
            three["exceed"],
            f"{three['h3']:.3f}" if three["h3"] is not None else "---",
            f"{g.transversion_ratio(controls[key]):.3f}",
        )))

    print("\nuntrained")
    for key in g.UNTRAINED_KEYS:
        stats_, three = mut[f"{key}_untrained"], rung3[f"{key}_untrained"]
        trained, untrained = g.attention_mean(runs[key]), \
            g.attention_mean(runs[f"{key}_untrained"])
        delta = ("---" if trained is None or untrained is None
                 else f"{trained - untrained:+.3f}")
        wins, total, _p = g.sign_test(mut, key)
        print(f"  {key:<12} ratio {stats_['mean_ratio']:.3f}  "
              f">null {stats_['exceeds_nuc']:>3}  "
              f"dinuc {stats_['survives_dinuc']}/{stats_['exceeds_nuc']:<3}  "
              f"attn {'---' if untrained is None else f'{untrained:.3f}'}  "
              f"delta {delta}  "
              f"PS {three['mean_ps']:.2e}  "
              f"gate {three['gate']}/{three['eligible']}  "
              f">null {three['exceed']}  "
              f"sign {wins}/{total}")

    trained_keys = [key for key, *_ in g.MODELS]
    ranked = sorted(trained_keys, key=lambda key: -(rung3[key]["mean_ps"] or 0.0))
    print("\nPS ranking:", ", ".join(
        f"{key} {rung3[key]['mean_ps']:.4g}" for key in ranked))
    print(f"second-to-third factor: "
          f"{rung3[ranked[1]]['mean_ps'] / rung3[ranked[2]]['mean_ps']:.1f}")
    for key in ranked[:3]:
        print(f"  {key}: peak layer {g.peak_layer(rung3[key])} of "
              f"{rung3[key]['n_ps_layers']}, "
              f"conservative {rung3[key]['exceed_conservative']}")

    print("\nfailures the build reports:")
    for line in g.claim_failures(runs, rung3, mut, controls) or ["  none"]:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
