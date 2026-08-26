"""RNA Structure Audit: three-rung benchmark for RNA/DNA foundation models."""

# Read from the installed distribution rather than restated here. Hardcoding it
# gives the version two sources of truth, and 0.2.0 shipped reporting itself as
# 0.1.0 because only pyproject.toml was bumped.
import importlib.metadata as _metadata

__version__ = _metadata.version("rna-structure-audit")
