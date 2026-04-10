"""Divider widget for visual separation.

Adapted from the dimbridge notebook in jupyter-scatter:
https://github.com/flekschas/jupyter-scatter/blob/main/notebooks/dimbridge.ipynb
"""

from pathlib import Path

from anywidget import AnyWidget
from traitlets import Dict

_STATIC = Path(__file__).parent.parent / "static" / "widgets"


class Div(AnyWidget):
    _esm = _STATIC / "divider.js"
    _css = _STATIC / "divider.css"

    style = Dict({}).tag(sync=True)
