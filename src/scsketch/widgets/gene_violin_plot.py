"""Violin-like plot for selected vs background expression distributions."""

from pathlib import Path

from anywidget import AnyWidget
from traitlets import Dict, Unicode

_STATIC = Path(__file__).parent.parent / "static" / "widgets"


class GeneViolinPlot(AnyWidget):
    """Compact violin-like plot for selected vs background distributions.

    This widget expects binned densities, not raw points, to keep payloads small.
    """

    _esm = _STATIC / "gene_violin_plot.js"
    _css = _STATIC / "gene_violin_plot.css"

    data = Dict(default_value={}).tag(sync=True)
    gene = Unicode("").tag(sync=True)
