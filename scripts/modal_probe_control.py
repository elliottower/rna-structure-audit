"""Run the probe's two constructed adapters inside the panel's own image.

`probe_positive_control_v2.py` needs torch and the Phase 6 script. This runs it
on Modal, in the multimolecule image `modal_decoy_panel.py` defines, on a CPU,
and writes what it returns to `results/audit/probe_positive_control_v2.json`.

No checkpoint inside the unit: the whole computation is two constructed
adapters over the panel's families with no model loaded, a few minutes of CPU,
and it writes one result. A run that dies is rerun from the start.

    modal run scripts/modal_probe_control.py
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import modal

from modal_decoy_panel import _library_versions, multimol_image

app = modal.App("rna-probe-control")
REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "results" / "audit" / "probe_positive_control_v2.json"


@app.function(image=multimol_image, timeout=3600)
def probe_control():
    import os
    os.chdir("/root/project")
    from probe_positive_control_v2 import run

    out = run()
    out["libraries"] = _library_versions()
    return out


@app.local_entrypoint()
def main():
    if OUT.exists():
        raise SystemExit(f"{OUT} exists; a result file is not overwritten")
    out = probe_control.remote()
    out["recorded"] = datetime.now(timezone.utc).isoformat()
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    for name, row in out["adapters"].items():
        print(f"  {name:8s} n={row['n']:3d} precision {row['precision']:.3f} "
              f"chance {row['chance_matched']:.3f} decoy rate {row['decoy_rate_all']}")
    print(f"wrote {OUT}")
