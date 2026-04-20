"""ScScatter — scSketch's customizable scatter plot component.

Two-layer subclass architecture
--------------------------------
* **Layer 1 — Python orchestrator** (this file):
  ``ScScatter`` subclasses ``jscatter.Scatter``.  It enforces scSketch-specific
  defaults, to add scSketch helper methods in the future, and overrides ``widget`` 
  so that the underlying anywidget instance is our own ``ScJupyterScatter`` instead of
  jscatter's ``JupyterScatter``.

* **Layer 2 — Frontend widget** (also this file):
  ``ScJupyterScatter`` subclasses ``jscatter.widget.JupyterScatter``.  Right now
  it is a transparent pass-through — identical JS/CSS as upstream.  When we need
  to implement our own toolbar, interaction modes or UI design changes we will swap
  ``_esm``/``_css`` to point at files in ``src/scsketch/static/``.

Upgrade path
------------
To customise the frontend:
1. Copy ``jscatter``'s ``bundle.js`` into ``src/scsketch/static/scscatter.js``
   (or author your own ES-module that re-exports from jscatter's bundle).
2. Create ``src/scscatter.css`` with any CSS overrides.
3. Point ``ScJupyterScatter._esm`` and ``._css`` at those files.
4. No changes needed in the rest of scSketch — ``ScScatter.widget`` already
   instantiates ``ScJupyterScatter``.
"""

from __future__ import annotations

import pathlib
from typing import Optional, Union, List

from jscatter.jscatter import (
    Scatter,
    get_scale,
    get_domain,
    normalize_annotations,
    sanitize_tooltip_properties,
    order_map,
    DEFAULT_TILE_SIZE,
)
from jscatter.utils import to_hex, to_scale_type
from jscatter.widget import JupyterScatter

# ---------------------------------------------------------------------------
# Layer 2 — Frontend widget
# ---------------------------------------------------------------------------

class ScJupyterScatter(JupyterScatter):
    """Subclass of jscatter's anywidget.

    Currently uses jscatter's own JS bundle unchanged.  Override ``_esm``
    and/or ``_css`` here when scSketch ships a customised frontend.

    Example (future)::

        _esm = pathlib.Path(__file__).parent / "static" / "scscatter.js"
        _css = pathlib.Path(__file__).parent / "static" / "scscatter.css"
    """
    # Inherits _esm from JupyterScatter (jscatter's bundle.js) — no change yet.
    pass


# ---------------------------------------------------------------------------
# Layer 1 — Python orchestrator
# ---------------------------------------------------------------------------

# These are the scSketch-specific defaults applied whenever ScScatter is
# instantiated without an explicit override.
_SCSKETCH_DEFAULTS: dict = {
    "axes": False,
    "background_color": "#111111",
    "tooltip": True,
    "legend": False,
}


class ScScatter(Scatter):
    """scSketch's scatter plot component.

    A drop-in replacement for ``jscatter.Scatter`` that:

    * Enforces scSketch visual defaults (dark background, no axes, etc.).
    * Instantiates ``ScJupyterScatter`` instead of ``JupyterScatter`` so the
      frontend can be customised independently of jscatter releases.
    * Provides scSketch-specific helper methods (added here as the need arises).

    Usage inside ScSketch::

        self.scatter = ScScatter(
            data=self.df,
            x="x",
            y="y",
            height=self.height,
            color_by=self.color_by_default,
            color_map=self.color_map_default,
            tooltip_properties=[...],
        )

    All ``jscatter.Scatter`` keyword arguments are accepted; scSketch defaults
    are applied *before* forwarding to the parent ``__init__``, so any explicit
    caller value wins.
    """

    def __init__(self, x, y, data=None, **kwargs):
        # Apply scSketch defaults for any key not explicitly supplied by caller.
        for key, default_value in _SCSKETCH_DEFAULTS.items():
            kwargs.setdefault(key, default_value)
        super().__init__(x, y, data=data, **kwargs)

    # ------------------------------------------------------------------
    # Override the lazy ``widget`` property to return ScJupyterScatter
    # ------------------------------------------------------------------

    @property
    def widget(self) -> ScJupyterScatter:  # type: ignore[override]
        """Return (and lazily create) the underlying ``ScJupyterScatter`` instance."""
        if self._widget is not None:
            return self._widget  # type: ignore[return-value]

        # Replicate the parent's lazy-init block, substituting ScJupyterScatter.
        # We forward all the same arguments as jscatter.Scatter.widget does.
        # Keep this in sync with jscatter when upgrading the dependency.
        self._widget = ScJupyterScatter(
            data=self._data,
            label_placement=self._label_placement,
            annotations=normalize_annotations(
                self._annotations, self._x_scale, self._y_scale
            ),
            axes=self._axes,
            axes_color=self.get_axes_color(),
            axes_grid=self._axes_grid,
            axes_labels=self.get_axes_labels(),
            background_color=self._background_color,
            background_image=self._background_image,
            camera_distance=self._camera_distance,
            camera_rotation=self._camera_rotation,
            camera_target=self._camera_target,
            camera_view=self._camera_view,
            camera_is_fixed=self._camera_is_fixed,
            color=order_map(self._color_map, self._color_order)
            if self._color_map
            else self._color,
            color_by=self.js_color_by,
            color_domain=get_domain(self, "color"),
            color_histogram=self._color_histogram,
            color_histogram_range=self.get_histogram_range("color"),
            color_hover=self._color_hover,
            color_scale=get_scale(self, "color"),
            color_selected=self._color_selected,
            color_title=self._color_by,
            connect=bool(self._connect_by),
            connection_color=order_map(
                self._connection_color_map, self._connection_color_order
            )
            if self._connection_color_map
            else self._connection_color,
            connection_color_by=self.js_connection_color_by,
            connection_color_hover=self._connection_color_hover,
            connection_color_selected=self._connection_color_selected,
            connection_opacity=order_map(
                self._connection_opacity_map, self._connection_opacity_order
            )
            if self._connection_opacity_map
            else self._connection_opacity,
            connection_opacity_by=self.js_connection_opacity_by,
            connection_size=order_map(
                self._connection_size_map, self._connection_size_order
            )
            if self._connection_size_map
            else self._connection_size,
            connection_size_by=self.js_connection_size_by,
            filter=self._filtered_points_idxs,
            height=self._height,
            labels=self._labels,
            label_shadow_color=to_hex(self._background_color)
            if self._label_shadow_color == "auto"
            else self._label_shadow_color,
            label_align=self._label_align,
            label_offset=self._label_offset,
            label_scale_function=self._label_scale_function,
            label_tile_size=self._label_placement.tile_size
            if self._label_placement
            else DEFAULT_TILE_SIZE,
            lasso_color=self._lasso_color,
            lasso_initiator=self._lasso_initiator,
            lasso_min_delay=self._lasso_min_delay,
            lasso_min_dist=self._lasso_min_dist,
            lasso_on_long_press=self._lasso_on_long_press,
            lasso_type=self._lasso_type,
            lasso_brush_size=self._lasso_brush_size,
            legend=self._legend,
            legend_color=self.get_legend_color(),
            legend_encoding=self.get_legend_encoding(),
            legend_position=self._legend_position,
            legend_size=self._legend_size,
            mouse_mode=self._mouse_mode,
            opacity=order_map(self._opacity_map, self._opacity_order)
            if self._opacity_map
            else self._opacity,
            opacity_by=self.js_opacity_by,
            opacity_domain=get_domain(self, "opacity"),
            opacity_histogram=self._opacity_histogram,
            opacity_histogram_range=self.get_histogram_range("opacity"),
            opacity_scale=get_scale(self, "opacity"),
            opacity_title=self._opacity_by,
            opacity_unselected=self._opacity_unselected,
            points=self.get_point_list(),
            reticle=self._reticle,
            reticle_color=self.get_reticle_color(),
            selection=self._selected_points_idxs,
            size=order_map(self._size_map, self._size_order)
            if self._size_map
            else self._size,
            size_by=self.js_size_by,
            size_domain=get_domain(self, "size"),
            size_histogram=self._size_histogram,
            size_histogram_range=self.get_histogram_range("size"),
            size_scale=get_scale(self, "size"),
            size_title=self._size_by,
            size_scale_function=self._size_scale_function,
            tooltip_enable=self._tooltip,
            tooltip_color=self.get_tooltip_color(),
            tooltip_properties=sanitize_tooltip_properties(
                self._data,
                ["x", "y", "color", "opacity", "size"],  # visual_properties
                self._tooltip_properties,
            ),
            tooltip_properties_non_visual_info=self._tooltip_properties_non_visual_info,
            tooltip_histograms=self._tooltip_histograms,
            tooltip_histograms_ranges=self._tooltip_histograms_ranges,
            tooltip_histograms_size=self._tooltip_histograms_size,
            tooltip_preview=self._tooltip_preview,
            tooltip_preview_type=self._tooltip_preview_type,
            tooltip_preview_text_lines=self._tooltip_preview_text_lines,
            tooltip_preview_image_background_color=self._tooltip_preview_image_background_color,
            tooltip_preview_image_position=self._tooltip_preview_image_position,
            tooltip_preview_image_size=self._tooltip_preview_image_size,
            tooltip_preview_image_height=self._tooltip_preview_image_height,
            tooltip_preview_audio_length=self._tooltip_preview_audio_length,
            tooltip_preview_audio_loop=self._tooltip_preview_audio_loop,
            tooltip_preview_audio_controls=self._tooltip_preview_audio_controls,
            tooltip_size=self._tooltip_size,
            width=self._width,
            x_domain=self._x_domain,
            x_histogram=self._x_histogram,
            x_histogram_range=self.get_histogram_range("x"),
            x_scale=to_scale_type(self._x_scale),
            x_scale_domain=self._x_scale_domain,
            x_title=self._x_by,
            y_domain=self._y_domain,
            y_histogram=self._y_histogram,
            y_histogram_range=self.get_histogram_range("y"),
            y_scale=to_scale_type(self._y_scale),
            y_scale_domain=self._y_scale_domain,
            y_title=self._y_by,
            zoom_animation=self._zoom_animation,
            zoom_on_filter=self._zoom_on_filter,
            zoom_on_selection=self._zoom_on_selection,
            zoom_padding=self._zoom_padding,
            zoom_to=self._zoom_to,
            regl_scatterplot_options=self._regl_scatterplot_options,
            transition_points_duration=self._transition_points_duration,
        )

        return self._widget  # type: ignore[return-value]
