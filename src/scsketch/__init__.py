import importlib.metadata
import pathlib
import anywidget
import traitlets

try:
    __version__ = importlib.metadata.version("scsketch")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"

class Widget(anywidget.AnyWidget):
    _esm = pathlib.Path(__file__).parent / "static" / "widget.js"
    _css = pathlib.Path(__file__).parent / "static" / "widget.css"
    value = traitlets.Int(0).tag(sync=True)

# Public view comes from directional.py
from .directional import view

__all__ = ["__version__", "Widget", "view"]