# Runbook: matched-decoy rerun

The commands that execute the run matrix of `PREREG.md`, in order. Each is run from the
repository root. Model lists are explicit because several arms take a flag that applies to
every model named in the same command.

## 1. Model-side tests

    modal run scripts/modal_decoy_panel.py --test-only

Writes `results/audit/model_side_tests.json`. Every later command refuses to start unless
that record shows every test passed with none skipped and the tested files still have the
hashes it names.

## 2. Seal and commit

Seal `scripts/phase6_compensatory_mutation_v3.py` and `scripts/rnafm_adapter_v2.py` with
`results seal`, write their hashes to the registration's log with `prereg log`, add the
`DEVIATIONS.md` entry, and commit. The launcher refuses a tree with uncommitted code.

## 3. Smoke checks

    modal run scripts/modal_decoy_panel.py --smoke-only
    modal run scripts/modal_decoy_panel.py --smoke-only --models ernierna_untrained --ablate-bias
    modal run scripts/modal_decoy_panel.py --smoke-only --models ernierna_untrained,rinalmo_untrained,rnafm_untrained --reseed-buffers

No tested file is edited after this point. An edit sends the run back to step 1.

## 4. The matrix

| arm | command |
|---|---|
| Phase 6, trained | `modal run --detach scripts/modal_decoy_panel.py --phase6-only --models rinalmo,ernierna,utrlm,splicebert,caduceus,hyenadna,evo,nt,dnabert2` |
| Phase 6, randomly initialized | `modal run --detach scripts/modal_decoy_panel.py --phase6-only --models ernierna_untrained,rinalmo_untrained,splicebert_untrained,utrlm_untrained,nt_untrained,dnabert2_untrained` |
| Phase 6, randomly initialized, bias ablated | `modal run --detach scripts/modal_decoy_panel.py --phase6-only --ablate-bias --models ernierna_untrained` |
| RNA-FM, every stage | `modal run --detach scripts/modal_decoy_panel.py --models rnafm,rnafm_untrained` |
| RNA-FM, transversion | `modal run --detach scripts/modal_decoy_panel.py --transversion --models rnafm` |
| RNA-FM, synthetic sequences | `modal run --detach scripts/modal_decoy_panel.py --synthetic --phase6-only --models rnafm` |
| RNA-FM, multi-sequence | `modal run --detach scripts/modal_decoy_panel.py --multi-seq --models rnafm` |
| RNA-FM, HTT fragments | `modal run --detach scripts/modal_decoy_panel.py --htt --models rnafm` |
| RNA-FM, randomly initialized, buffers reseeded | `modal run --detach scripts/modal_decoy_panel.py --reseed-buffers --phase6-only --models rnafm_untrained` |

Each launch writes `results/audit/launch_manifests/launch_<UTC time>.json`. Progress is read
from app state (`modal app list`, `modal app logs rna-decoy-panel`), not from the number of
files on the volume.

## 5. Fetch, gate, audit

Fetch the volume `rna-decoy-panel-results` to a new folder under `results/`, record the
outputs with `results run`, then:

    uv run python scripts/validate_decoy_controls_v2.py
    uv run python scripts/decoy_audit_v2.py --panel results/<fetched folder>

The audit refuses to run without a passing gate record, and stops if any of its ten folders
lacks exactly one Phase 6 file.
