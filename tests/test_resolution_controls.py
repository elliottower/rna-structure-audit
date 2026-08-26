"""The metric resolves below float32, and a deterministic model floors at zero.

D15 and D16 turned on quantities that were computed and thrown away. These
check the ones that were added back: the no-op floor, which measures resolution
instead of assuming it, the chance rate the derangement null already knew, and
the neighbour asymmetry that makes D16's mechanism testable.
"""

import numpy as np
import pytest
import torch

from phase6_compensatory_mutation import (
    cosine_distance, load_rfam_families, noop_floor, run_phase6,
)
from test_checkpoint_resume import StubAdapter


FLOAT32_ULP_NEAR_005 = 2.0 ** -28


def a_pair_about_five_percent_apart():
    """Two float32 vectors whose cosine distance lands in [2^-5, 2^-4).

    That is the interval the deposited values came from: every stored PS was an
    integer multiple of 2^-28, which is float32's spacing there, because the
    distances were float32 and their difference inherited their granularity.
    """
    # Rotating `a` by a fixed angle puts the distance at 1 - cos(theta) by
    # construction, so this needs no rejection sampling and cannot hang.
    theta = 0.32
    a = torch.rand(64, dtype=torch.float64) + 0.5
    w = torch.rand(64, dtype=torch.float64) - 0.5
    unit = a / a.norm()
    w = w - (w @ unit) * unit
    b = a * torch.cos(torch.tensor(theta)) + w / w.norm() * a.norm() * torch.sin(
        torch.tensor(theta))
    return a.float(), b.float()


def test_the_metric_carries_more_than_a_float32_could_hold():
    """A float32 in [2^-5, 2^-4) is an integer multiple of 2^-28 by construction.

    So a metric that subtracts in float32 can only ever return multiples of
    2^-28 at this magnitude, and PS -- a difference of two such numbers --
    inherits that floor. The repaired metric must land off that lattice.
    """
    off_lattice = 0
    for _ in range(200):
        a, b = a_pair_about_five_percent_apart()
        distance = cosine_distance(a, b)

        in_float32 = (1.0 - torch.nn.functional.cosine_similarity(
            a.unsqueeze(0), b.unsqueeze(0))).item()
        assert abs(in_float32 / FLOAT32_ULP_NEAR_005
                   - round(in_float32 / FLOAT32_ULP_NEAR_005)) < 1e-6, (
            "the float32 route must sit on the 2^-28 lattice, or this test is "
            "not measuring what the deposited values suffered from")

        quotient = distance / FLOAT32_ULP_NEAR_005
        off_lattice += abs(quotient - round(quotient)) > 1e-6

    assert off_lattice > 190, (
        f"only {off_lattice}/200 distances carry information finer than "
        "float32's spacing, so the metric is still quantized")


def test_a_deterministic_model_measures_a_floor_far_below_the_signal():
    families = sorted(load_rfam_families(), key=lambda f: len(f["sequence"]))
    adapter = StubAdapter()
    for family in families:
        result = run_phase6(adapter, [family], compute_null=False)
        entry = result["per_rna"][family["name"]]
        if entry.get("skipped"):
            continue
        floor = entry["noop_floor"]
        # The two passes are bitwise equal, so what is left is the cosine's own
        # rounding on each position, not a difference between the passes. That
        # residue is float64-sized; the signal RiNALMo carries is 0.2.
        assert abs(floor["best_ps_noop"]) < 1e-12, (
            "two forward passes of one sequence through a deterministic model "
            "must agree to far below anything the panel reports")
        assert all(abs(v) < 1e-12 for v in floor["per_layer_ps_noop"])
        return
    pytest.skip("no family survives the registered filters")


def test_the_floor_control_reads_the_same_rows_the_metric_does():
    families = sorted(load_rfam_families(), key=lambda f: len(f["sequence"]))
    adapter = StubAdapter()
    for family in families:
        result = run_phase6(adapter, [family], compute_null=False)
        entry = result["per_rna"][family["name"]]
        if entry.get("skipped"):
            continue
        assert (len(entry["noop_floor"]["per_layer_ps_noop"])
                == len(entry["per_layer_ps"])), (
            "the control must cover every layer the statistic is maximized over")
        return
    pytest.skip("no family survives the registered filters")


def test_the_derangement_reports_a_chance_rate_and_the_neighbors_an_asymmetry():
    families = sorted(load_rfam_families(), key=lambda f: len(f["sequence"]))
    adapter = StubAdapter()
    for family in families:
        result = run_phase6(adapter, [family], compute_null=True)
        entry = result["per_rna"][family["name"]]
        if entry.get("skipped") or not entry.get("null_available"):
            continue
        chance = entry["h3_chance_fraction"]
        assert chance is not None and 0.0 <= chance <= 1.0
        nearer = entry["nearer_neighbor_larger"]
        assert nearer is not None and 0.0 <= nearer <= 1.0
        return
    pytest.skip("no family has a derangement null")
