"""Interactive SVG viewer widget with zoom and pan capabilities."""

from pathlib import Path

from anywidget import AnyWidget
from traitlets import Unicode

_STATIC = Path(__file__).parent.parent / "static" / "widgets"


class InteractiveSVG(AnyWidget):
    _esm = _STATIC / "interactive_svg.js"
    _css = _STATIC / "interactive_svg.css"

    svg_content = Unicode("").tag(sync=True)
