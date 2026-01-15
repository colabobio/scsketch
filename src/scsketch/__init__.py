"""
scSketch: Interactive exploration of single-cell embeddings with directional analysis.

scSketch is an interactive tool for exploring single-cell data embeddings (UMAP, tSNE, etc.)
in Python notebooks. It provides directional analysis capabilities to identify genes varying
along user-specified directions and integrates with Reactome for pathway exploration.
"""

from .scsketch import ScSketch
from .widgets import (
    CorrelationTable,
    PathwayTable,
    InteractiveSVG,
    Label,
    Div,
    GenePathwayWidget,
)
from .utils import Selection, Selections, Lasso, create_selection, fetch_pathways
from .analysis import compute_directional_analysis

__version__ = "0.1.0"

__all__ = [
    "ScSketch",
    "CorrelationTable",
    "PathwayTable",
    "InteractiveSVG",
    "Label",
    "Div",
    "GenePathwayWidget",
    "Selection",
    "Selections",
    "Lasso",
    "create_selection",
    "fetch_pathways",
    "compute_directional_analysis",
]
