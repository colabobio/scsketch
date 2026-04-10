"""Interactive scatterplot of gene expression vs projection."""

from pathlib import Path

from anywidget import AnyWidget
from traitlets import Dict, List, Unicode

_STATIC = Path(__file__).parent.parent / "static" / "widgets"


class GeneProjectionPlot(AnyWidget):
    """Interactive scatterplot of gene expression vs projection."""

    _esm = _STATIC / "gene_projection_plot.js"
    _css = _STATIC / "gene_projection_plot.css"

    data = List(Dict()).tag(sync=True)  # [{"projection": float, "expression": float}]
    gene = Unicode("").tag(sync=True)
