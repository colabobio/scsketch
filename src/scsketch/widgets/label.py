"""Label widget for selection display.

Adapted from the dimbridge notebook in jupyter-scatter:
https://github.com/flekschas/jupyter-scatter/blob/main/notebooks/dimbridge.ipynb
"""

from pathlib import Path

from anywidget import AnyWidget
from traitlets import Bool, Dict, Unicode

_STATIC = Path(__file__).parent.parent / "static" / "widgets"


class Label(AnyWidget):
    _esm = _STATIC / "label.js"
    _css = _STATIC / "label.css"

    name = Unicode("").tag(sync=True)
    style = Dict({}).tag(sync=True)
    focus = Bool(False).tag(sync=True)
