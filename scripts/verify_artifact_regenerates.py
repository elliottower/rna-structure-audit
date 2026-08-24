"""Does results/bootstrap_cis.json still reproduce from the inputs in this repo?

Run:  uv run --no-project --with numpy --with tqdm --python 3.12 python \
          scripts/verify_artifact_regenerates.py

Re-runs scripts/compute_bootstrap_cis.py into a temporary file and compares it to
the committed artifact field by field. Two fields are expected to differ and are
ignored: `timestamp`, and the `file` recorded for each model, which was an
absolute path into a repository that no longer exists before the inputs were
brought into this tree.

The bootstrap is seeded, so any numeric difference at all means an input changed.
"""

import json
import tempfile
from pathlib import Path

from paper_versions import REPO

import compute_bootstrap_cis

ARTIFACT = REPO / "results" / "bootstrap_cis.json"
IGNORED = {"timestamp", "file"}


def flatten(node, path: str = "") -> dict[str, object]:
    """Every leaf in the tree, keyed by its dotted path."""
    if isinstance(node, dict):
        return {k: v for key, child in node.items() if key not in IGNORED
                for k, v in flatten(child, f"{path}.{key}" if path else key).items()}
    if isinstance(node, list):
        return {k: v for i, child in enumerate(node)
                for k, v in flatten(child, f"{path}[{i}]").items()}
    return {path: node}


def main() -> int:
    committed = flatten(json.loads(ARTIFACT.read_text()))

    with tempfile.TemporaryDirectory() as tmp:
        candidate_path = Path(tmp) / "bootstrap_cis.json"
        compute_bootstrap_cis.main(candidate_path)
        candidate = flatten(json.loads(candidate_path.read_text()))

    only_committed = sorted(set(committed) - set(candidate))
    only_candidate = sorted(set(candidate) - set(committed))
    differing = sorted(k for k in set(committed) & set(candidate)
                       if committed[k] != candidate[k])

    print(f"\n{len(committed)} leaf values compared\n")
    for label, keys in [("only in the committed artifact", only_committed),
                        ("only in the fresh run", only_candidate),
                        ("differing", differing)]:
        if keys:
            print(f"{len(keys)} {label}:")
            for key in keys[:20]:
                print(f"  {key}: {committed.get(key)!r} vs {candidate.get(key)!r}")
            if len(keys) > 20:
                print(f"  ... and {len(keys) - 20} more")

    if only_committed or only_candidate or differing:
        return 1
    print("The committed artifact reproduces exactly from the inputs in this repo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
