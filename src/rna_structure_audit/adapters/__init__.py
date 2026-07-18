"""Built-in model adapters for common RNA/DNA foundation models.

These adapters require the `all-models` optional dependency:
    pip install rna-structure-audit[all-models]
"""

from rna_structure_audit.adapters.caduceus import CaduceusAdapter
from rna_structure_audit.adapters.dnabert2 import DNABERT2Adapter
from rna_structure_audit.adapters.ernierna import ERNIERNAAdapter
from rna_structure_audit.adapters.evo import EvoAdapter
from rna_structure_audit.adapters.hyenadna import HyenaDNAAdapter
from rna_structure_audit.adapters.nt import NTAdapter
from rna_structure_audit.adapters.rinalmo import RiNALMoAdapter
from rna_structure_audit.adapters.rnafm import RNAFMAdapter
from rna_structure_audit.adapters.splicebert import SpliceBERTAdapter
from rna_structure_audit.adapters.utrlm import UTRLMAdapter

ALL_ADAPTERS = [
    RNAFMAdapter,
    NTAdapter,
    HyenaDNAAdapter,
    EvoAdapter,
    CaduceusAdapter,
    RiNALMoAdapter,
    UTRLMAdapter,
    ERNIERNAAdapter,
    SpliceBERTAdapter,
    DNABERT2Adapter,
]

__all__ = [
    "CaduceusAdapter",
    "DNABERT2Adapter",
    "ERNIERNAAdapter",
    "EvoAdapter",
    "HyenaDNAAdapter",
    "NTAdapter",
    "RiNALMoAdapter",
    "RNAFMAdapter",
    "SpliceBERTAdapter",
    "UTRLMAdapter",
    "ALL_ADAPTERS",
]
