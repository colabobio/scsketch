"""Interactive correlation table widget for directional analysis results."""

from pathlib import Path

import traitlets
from anywidget import AnyWidget
from traitlets import Dict, List, Unicode

_STATIC = Path(__file__).parent.parent / "static" / "widgets"


class CorrelationTable(AnyWidget):
    _esm = _STATIC / "correlation_table.js"
    _css = _STATIC / "correlation_table.css"

    data = List(Dict()).tag(sync=True)
    columns = List(Unicode(), default_value=["Gene", "R", "p", "Selection"]).tag(sync=True)
    selected_gene = traitlets.Unicode("").tag(sync=True)
