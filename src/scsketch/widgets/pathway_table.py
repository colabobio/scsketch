"""Pathway table widget — displays a list of Reactome pathways."""

from pathlib import Path

import traitlets
from anywidget import AnyWidget
from traitlets import List, Unicode

_STATIC = Path(__file__).parent.parent / "static" / "widgets"


class PathwayTable(AnyWidget):
    _esm = _STATIC / "pathway_table.js"
    _css = _STATIC / "pathway_table.css"

    data = List(default_value=[]).tag(sync=True)
    selected_pathway = Unicode("").tag(sync=True)
