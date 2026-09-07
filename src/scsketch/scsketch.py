"""Main ScSketch widget for interactive single-cell embedding exploration.

Architecture overview
---------------------
``ScSketch`` is a thin *orchestrator*:

* Data preparation  -> :mod:`scsketch._data`        (``build_embedding_df``)
* External HTTP     -> :mod:`scsketch._api`          (``fetch_pathways``, ``fetch_pathway_svg``)
* Diff. expression  -> :mod:`scsketch._diffexpr`     (``DiffExprEngine``)
* Widget building   -> :mod:`scsketch._ui`           (``build_controls``, ``UIControls``)
* Results display   -> :mod:`scsketch._results`      (``show_directional_results``, ...)
* Pure statistics   -> :mod:`scsketch.analysis`      (``test_direction``, ``lord_test``)

The UI composition and selection management patterns in this module are adapted from
the dimbridge notebook in jupyter-scatter by Fritz Lekschas:
https://github.com/flekschas/jupyter-scatter/blob/main/notebooks/dimbridge.ipynb
"""

from __future__ import annotations

import json
from html import escape
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from ipywidgets import GridBox, Layout
from jscatter.widgets import Button
from matplotlib import colormaps
from matplotlib.colors import to_hex
from scipy.spatial import ConvexHull

from jscatter import Line, Scatter, okabe_ito

from ._action_log import append_action as _append_action
from ._action_log import build_action_log_table as _build_action_log_table
from ._action_log import normalize_action_log as _normalize_action_log
from ._analysis import lord_test, test_direction
from ._data import build_embedding_df
from ._diffexpr import DiffExprEngine
from ._logging import LogLevel, configure_logging
from ._results import clear_results, show_diffexpr_results, show_directional_results
from ._scores import discovery_scores
from ._scatter import ScScatter
from ._session import (
    apply_playback_step as _apply_playback_step,
)
from ._session import (
    build_session_player as _build_session_player,
)
from ._session import (
    export_session as _export_session,
)
from ._session import (
    load_session as _load_session,
)
from ._session import (
    read_session as _read_session,
)
from ._session import (
    write_session as _write_session,
)
from ._ui import (
    UIControls,
    build_controls,
    clear_analysis_progress,
    set_analysis_progress,
)
from ._utils import (
    Lasso,
    Selection,
    Selections,
    points_in_polygon,
    split_line_equidistant,
)


def _uploaded_files(value):
    """Return uploaded file records across ipywidgets 7/8 value formats."""
    if isinstance(value, dict):
        return list(value.values())
    if isinstance(value, (list, tuple)):
        return list(value)
    return []


def _uploaded_content(uploaded) -> bytes:
    content = (
        uploaded.get("content")
        if isinstance(uploaded, dict)
        else getattr(uploaded, "content", None)
    )
    return bytes(content)


def _uploaded_name(uploaded) -> str:
    if isinstance(uploaded, dict):
        return str(uploaded.get("name") or "")
    return str(getattr(uploaded, "name", "") or "")


class ScSketch:
    """
    ScSketch: Interactive exploration of single-cell embeddings with directional analysis.

    This widget provides an interactive interface for exploring single-cell data embeddings,
    performing directional analysis to identify genes varying along user-specified directions,
    and exploring Reactome pathways.
    """

    def __init__(
        self,
        adata: AnnData,
        metadata_cols: Optional[List[str]] = None,
        color_by_default: str = "seurat_clusters",
        height: int = 720,
        background_color: str = "#111111",
        max_genes: int = 0,
        fdr_alpha: float = 0.05,
        verbosity: LogLevel = "warning",
        diffexpr_disk_cache_dir: str | Path | None = None,
        extra_views: dict[str, Scatter] | None = None,
        gene_annotation_species: str | int = "human",
        reactome_species: str | int = "human",
    ):
        """
        Initialize ScSketch widget.

        Args:
            adata: AnnData object with 'X_umap' in obsm and gene expression data
            metadata_cols: List of metadata column names for color encoding
            color_by_default: Default column to color by
            height: Height of the scatter plot in pixels
            background_color: Background color of the scatter plot
            max_genes: Maximum number of genes to preload into the plot DataFrame (0 for none).
                Directional analysis still uses all genes from ``adata`` regardless of this setting.
            fdr_alpha: False discovery rate alpha threshold for directional analysis
            verbosity: Logging verbosity level
            diffexpr_disk_cache_dir: Optional path to persist DE global stats cache to disk.
            extra_views: Optional mapping of labels to additional
                :class:`jscatter.Scatter` instances shown in the right panel
                when multi-view mode is enabled. Extra views are matched to the
                main scSketch view by row index.
            gene_annotation_species: Species filter passed to MyGene.info gene
                annotation lookups. Defaults to human.
            reactome_species: Species filter passed to Reactome pathway lookups.
                Defaults to human.
        """
        self.logger = configure_logging(verbosity)
        self.adata = adata
        self.metadata_cols = metadata_cols
        self.height = height
        self.background_color = background_color
        self.max_genes = max_genes
        self.fdr_alpha = fdr_alpha
        self.verbosity = verbosity
        self.gene_annotation_species = gene_annotation_species
        self.reactome_species = reactome_species
        self.extra_views: dict[str, Scatter] = dict(extra_views or {})
        self._multi_view_syncing = False
        self._extra_view_size = max(240, min(int(height), 420))
        self._gene_display_names = self._build_gene_display_names(adata)

        # -- Build DataFrame ------------------------------------------------
        result = build_embedding_df(
            adata=adata,
            metadata_cols=metadata_cols,
            color_by_default=color_by_default,
            max_genes=max_genes,
        )
        self.df: pd.DataFrame = result["df"]
        self.available_metadata_cols: list[str] = result["available_metadata_cols"]
        self.all_gene_sorted: list[str] = result["all_gene_sorted"]
        self.meta_cols_present: list[str] = result["meta_cols_present"]
        self.categorical_cols: list[str] = result["categorical_cols"]
        self.categorical_color_maps: dict[str, dict] = result["categorical_color_maps"]
        self.color_map_default = result["color_map_default"]
        self.color_by_default: str | None = result["color_by_default"]

        # -- Selection state ------------------------------------------------
        self.lasso = Lasso()
        self.selections = Selections()
        self.active_selection: Selection | None = None
        self.analysis_mode: str = "directional"  # "directional" | "differential"
        self.action_log: list[dict] = []
        self._action_log_paused = False
        self._history_dropdown_paused = False
        self._selection_archive: dict[str, Selection] = {}

        # -- LORD++ online state --------------------------------------------
        self.batch_results = None
        self.online_results = None

        # -- Colour cycling -------------------------------------------------
        self.available_colors = list(okabe_ito.copy())
        self.continuous_color_maps = [
            ["#00dadb", "#da00db"],
            ["#00dadb", "#a994dc", "#da00db"],
            ["#00dadb", "#8faddc", "#bd77dc", "#da00db"],
            ["#00dadb", "#7eb9dc", "#a994dc", "#c567dc", "#da00db"],
            ["#00dadb", "#72c0db", "#9aa3dc", "#b583dc", "#ca5cdb", "#da00db"],
            ["#00dadb", "#69c4db", "#8faddc", "#a994dc", "#bd77dc", "#cd54db", "#da00db"],
            ["#00dadb", "#62c7db", "#86b4dc", "#9e9fdc", "#b288dc", "#c16edc", "#cf4ddb", "#da00db"],
            ["#00dadb", "#5ccadb", "#7eb9dc", "#96a7dc", "#a994dc", "#b87fdc", "#c567dc", "#d048db", "#da00db"],
            ["#00dadb", "#57ccdb", "#78bddc", "#8faddc", "#a19ddc", "#b08bdc", "#bd77dc", "#c861db", "#d144db", "#da00db"],
        ]

        # -- Pending DE result (before selection is named) ------------------
        self._pending_diffexpr: dict | None = None

        # -- Differential expression engine --------------------------------
        disk_cache = Path(diffexpr_disk_cache_dir) if diffexpr_disk_cache_dir is not None else None
        self._de_engine = DiffExprEngine(
            adata=adata,
            disk_cache_dir=disk_cache,
        )

        # -- Build scatter plot --------------------------------------------
        self.scatter = ScScatter(
            data=self.df,
            x="x",
            y="y",
            # scSketch defaults (axes=False, background_color, tooltip=True,
            # legend=False) are applied by ScScatter automatically.
            background_color=self.background_color,
            height=self.height,
            color_by=self.color_by_default,
            color_map=self.color_map_default,
            tooltip_properties=[c for c in self.df.columns if c in self.meta_cols_present],
        )
        self.scatter.widget.color_selected = "#00dadb"
        self._configure_extra_views()

        # -- Build UI ------------------------------------------------------
        self._ctrl: UIControls = build_controls(
            scatter=self.scatter,
            available_metadata_cols=self.available_metadata_cols,
            all_gene_sorted=self.all_gene_sorted,
            color_by_default=self.color_by_default,
            max_genes=self.max_genes,
            df_columns=list(self.df.columns),
            extra_views=self.extra_views,
        )
        if self.scatter.widget.lasso_type == "freeform":
            self.analysis_mode = "differential"

        self._setup_handlers()
        self._sync_extra_view_color(self.color_by_default)
        self._apply_multi_view_visibility()

    # -- Logging shortcut --------------------------------------------------

    def _log(self, *args):
        self.logger.debug(" ".join(str(a) for a in args))

    def _record_action(
        self,
        action_type: str,
        label: str,
        *,
        selection: str | None = None,
        payload: dict | None = None,
    ) -> dict | None:
        """Append a user-facing audit event unless logging is paused."""
        if getattr(self, "_action_log_paused", False):
            return None
        if not hasattr(self, "action_log"):
            self.action_log = []
        entry = _append_action(
            self.action_log,
            action_type,
            label,
            selection=selection,
            analysis_mode=self.analysis_mode,
            payload=payload,
        )
        self._refresh_history_options()
        return entry

    def _archive_selection(self, selection: Selection | None) -> None:
        """Keep a snapshot reference for playback of selections later removed."""
        if selection is not None:
            self._selection_archive[selection.name] = selection

    # -- Progress helpers --------------------------------------------------

    def _set_progress(self, mode: str, step: int, total: int, message: str):
        set_analysis_progress(self._ctrl, mode, step, total, message)

    def _clear_progress(self, mode: str):
        clear_analysis_progress(self._ctrl, mode)

    # -- Results display ---------------------------------------------------

    def _show_directional_results(self, directional_results, *, initial_gene=None):
        ctrl = self._ctrl
        ctrl.compute_predicates.style = ""
        ctrl.compute_predicates.description = "Clear Results"
        ctrl.compute_predicates.on_click(self._clear_predicates)

        show_directional_results(
            directional_results,
            selections_predicates=ctrl.selections_predicates,
            pathway_table_container=ctrl.pathway_table_container,
            reactome_diagram_container=ctrl.reactome_diagram_container,
            df=self.df,
            adata=self.adata,
            active_selection=self.active_selection,
            on_gene_selected=self._handle_gene_selected,
            on_results_cleared=self._clear_predicates,
            log=self._log,
            on_pathway_selected=self._handle_pathway_selected,
            initial_gene=initial_gene,
            show_gene_details=self._show_gene_details_in_right_panel,
            gene_display_name=self._display_gene_name,
            gene_annotation_species=self.gene_annotation_species,
            reactome_species=self.reactome_species,
        )
        self._apply_multi_view_visibility()

    def _show_diffexpr_results(
        self,
        diff_results,
        selection_label,
        *,
        selected_indices=None,
        initial_gene=None,
    ):
        ctrl = self._ctrl
        show_diffexpr_results(
            diff_results,
            selection_label,
            selected_indices=selected_indices,
            selections_predicates=ctrl.selections_predicates,
            pathway_table_container=ctrl.pathway_table_container,
            reactome_diagram_container=ctrl.reactome_diagram_container,
            de_source_fn=self._de_engine.de_source,
            active_selection=self.active_selection,
            scatter=self.scatter,
            on_gene_selected=self._handle_gene_selected,
            log=self._log,
            initial_gene=initial_gene,
            show_gene_details=self._show_gene_details_in_right_panel,
            gene_display_name=self._display_gene_name,
            gene_annotation_species=self.gene_annotation_species,
        )
        self._apply_multi_view_visibility()

    @staticmethod
    def _build_gene_display_names(adata: AnnData) -> dict[str, str]:
        columns = (
            "gene_short_name",
            "gene_symbols",
            "gene_symbol",
            "gene_name",
            "gene_names",
            "symbol",
        )
        mapping: dict[str, str] = {}

        def add_var_labels(var_names, var) -> None:
            for column in columns:
                if column not in var.columns:
                    continue
                labels = var[column].astype(str)
                if not labels.str.strip().replace({"nan": ""}).any():
                    continue
                for gene_id, label in zip(var_names, labels):
                    label = str(label).strip()
                    if label and label.lower() != "nan":
                        mapping.setdefault(str(gene_id), label)
                return

        add_var_labels(adata.var_names, adata.var)
        raw = getattr(adata, "raw", None)
        if raw is not None and getattr(raw, "var", None) is not None:
            add_var_labels(raw.var_names, raw.var)
        return mapping

    def _display_gene_name(self, gene: str) -> str:
        return self._gene_display_names.get(str(gene), str(gene))

    def _clear_results_display(self, message: str | None = None):
        ctrl = self._ctrl
        clear_results(
            ctrl.selections_predicates,
            ctrl.pathway_table_container,
            ctrl.reactome_diagram_container,
            message=message,
        )

    def _gene_expression_color_map(self, steps: int = 256) -> list[str]:
        """Return a perceptually uniform sequential map for expression coloring."""
        steps = max(2, int(steps))
        cmap = colormaps["viridis"]
        return [to_hex(cmap(t)) for t in np.linspace(0.0, 1.0, steps)]

    def _show_gene_expression_caption(self, gene: str, display_gene: str | None = None):
        """Show the gene-expression color caption beneath the scatterplot."""
        colors = self._gene_expression_color_map(7)
        denom = max(1, len(colors) - 1)
        stops = ", ".join(
            f"{color} {int(round(i * 100 / denom))}%"
            for i, color in enumerate(colors)
        )
        safe_gene = escape(display_gene or gene)
        self._ctrl.gene_expression_caption.value = f"""
<div style="display:flex;align-items:center;gap:8px;flex-wrap:wrap;
            margin:4px 0 2px 58px;font:12px sans-serif;color:#333;">
  <span><b>{safe_gene}</b> expression</span>
  <span>Low exp.</span>
  <span style="display:inline-block;width:150px;height:18px;
               border:1px solid #bbb;background:linear-gradient(to right,{stops});">
  </span>
  <span>Max exp.</span>
</div>
"""
        self._ctrl.gene_expression_caption.layout.display = "block"

    def _hide_gene_expression_caption(self):
        self._ctrl.gene_expression_caption.layout.display = "none"
        self._ctrl.gene_expression_caption.value = ""

    def _color_embedding_by_gene(self, gene: str):
        """Recolor the full embedding by expression of the selected gene."""
        if not gene:
            return
        display_gene = self._display_gene_name(gene)

        if gene in self.df.columns:
            expr = pd.to_numeric(self.df[gene], errors="coerce").to_numpy(dtype=float)
        else:
            sub = self.adata[:, [gene]].X
            if sp.issparse(sub):
                sub = sub.toarray()
            expr = np.asarray(sub, dtype=float).ravel()

        finite = np.isfinite(expr)
        if not finite.any():
            self._log(f"[color] {gene!r}: no finite expression values; skipping recolor")
            return
        if not finite.all():
            expr = expr.copy()
            expr[~finite] = float(np.median(expr[finite]))

        vmin = float(expr.min())
        vmax = float(expr.max())
        if vmax <= vmin:
            vmax = vmin + 1e-6

        self.scatter.color(
            by=expr,
            map=self._gene_expression_color_map(),
            norm=(vmin, vmax),
            labeling={"minValue": "Low", "maxValue": "High", "variable": display_gene},
        )
        self._color_extra_views_by_gene_expression(
            gene, display_gene, expr, vmin, vmax
        )
        self._show_gene_expression_caption(gene, display_gene)
        self._refresh_scatter_non_spatial_points()
        self._log(
            f"[color] recolored embedding for gene {gene!r} ({display_gene!r}) "
            f"range=[{vmin:.4g}, {vmax:.4g}]"
        )

    def _handle_gene_selected(self, gene: str):
        """Record and handle a result-table gene click."""
        selection = None if self.active_selection is None else self.active_selection.name
        self._record_action(
            "select_gene",
            f"Selected gene {gene}",
            selection=selection,
            payload={"gene": gene},
        )
        self._color_embedding_by_gene(gene)

    def _handle_pathway_selected(self, pathway_id: str):
        """Record a Reactome pathway click."""
        selection = None if self.active_selection is None else self.active_selection.name
        self._record_action(
            "select_pathway",
            f"Selected pathway {pathway_id}",
            selection=selection,
            payload={"pathway_id": pathway_id},
        )

    def _refresh_scatter_non_spatial_points(self):
        """Force jscatter to sync changed non-spatial encodings to the frontend."""
        if self.scatter is None or self.scatter.widget is None:
            return
        try:
            self.scatter.widget.prevent_filter_reset = True
            self.scatter.widget.non_spatial_points_update = True
            self.scatter.widget.points = self.scatter.get_point_list()
        except Exception:
            self.logger.exception("Failed to refresh scatter point encodings")

    # -- Scatter annotations -----------------------------------------------

    def _update_annotations(self):
        scatter = self.scatter
        if scatter is None:
            return
        try:
            lasso_polygon = [] if self.lasso.polygon is None else [self.lasso.polygon]
            overlays = (
                self.selections.all_lassos()
                + self.selections.all_direction_guides()
                + lasso_polygon
            )
            scatter.annotations(overlays)
        except Exception:
            self.logger.exception("Failed to update annotations")

    # -- Event handlers ----------------------------------------------------

    def _setup_handlers(self):
        """Wire all event handlers to scatter and UI controls."""
        scatter = self.scatter
        scatter.widget.observe(
            self._lasso_selection_polygon_change_handler,
            names=["lasso_selection_polygon"],
        )
        scatter.widget.observe(self._selection_handler, names=["selection"])
        scatter.widget.observe(self._lasso_type_change_handler, names=["lasso_type"])
        if hasattr(scatter.widget, "lasso_brush_size"):
            scatter.widget.observe(
                self._lasso_brush_size_change_handler,
                names=["lasso_brush_size"],
            )
        if self.extra_views:
            scatter.widget.observe(
                self._main_selection_multi_view_handler,
                names=["selection"],
            )
            for view in self.extra_views.values():
                view.widget.observe(
                    self._extra_view_selection_handler,
                    names=["selection"],
                )

        ctrl = self._ctrl
        ctrl.selection_add.on_click(self._selection_add_handler)
        ctrl.compute_predicates.on_click(self._compute_predicates_handler)
        ctrl.compute_diffexpr.on_click(self._compute_diffexpr_handler)
        ctrl.color_by.observe(self._color_by_change_handler, names=["value"])
        ctrl.compute_predicates_between_selections.observe(
            self._compare_between_change_handler,
            names=["value"],
        )
        ctrl.diff_t_threshold.observe(self._diff_threshold_change_handler, names=["value"])
        ctrl.diff_p_threshold.observe(self._diff_threshold_change_handler, names=["value"])
        ctrl.history_dropdown.observe(
            self._history_dropdown_change_handler,
            names=["value"],
        )
        ctrl.session_save.on_click(self._session_save_handler)
        ctrl.session_upload.observe(self._session_upload_handler, names=["value"])
        ctrl.multi_view_toggle.observe(self._multi_view_toggle_handler, names=["value"])
        self._refresh_history_options()

    def _show_gene_details_in_right_panel(self) -> bool:
        ctrl = getattr(self, "_ctrl", None)
        if ctrl is None or not self.extra_views:
            return True
        return not bool(ctrl.multi_view_toggle.value)

    def _apply_multi_view_visibility(self) -> None:
        ctrl = getattr(self, "_ctrl", None)
        if ctrl is None or not self.extra_views:
            return

        if bool(ctrl.multi_view_toggle.value):
            ctrl.multi_view_container.layout.display = "flex"
            ctrl.pathway_table_container.layout.display = "none"
            ctrl.reactome_diagram_container.layout.display = "none"
        else:
            ctrl.multi_view_container.layout.display = "none"
            if ctrl.pathway_table_container.children:
                ctrl.pathway_table_container.layout.display = (
                    "block" if self.analysis_mode == "differential" else "flex"
                )

    def _multi_view_toggle_handler(self, change):
        enabled = bool(change["new"])
        self._apply_multi_view_visibility()
        self._record_action(
            "toggle_multi_view",
            f"Multi-view set to {enabled}",
            payload={"multi_view": enabled},
        )

    def _main_selection_multi_view_handler(self, change):
        if self._multi_view_syncing:
            return
        self._multi_view_syncing = True
        try:
            selection = self._clean_multi_view_selection(change["new"])
            if not selection and self.active_selection is not None:
                selection = self._clean_multi_view_selection(
                    self.active_selection.points
                )
            for view in self.extra_views.values():
                view.selection(selection)
        finally:
            self._multi_view_syncing = False

    def _extra_view_selection_handler(self, change):
        if self._multi_view_syncing:
            return
        self._multi_view_syncing = True
        try:
            self.scatter.selection(self._clean_multi_view_selection(change["new"]))
        finally:
            self._multi_view_syncing = False

    def _current_multi_view_selection(self) -> list[int]:
        if self.scatter is not None:
            selection = self._clean_multi_view_selection(self.scatter.selection())
            if selection:
                return selection
        if self.active_selection is not None:
            return self._clean_multi_view_selection(self.active_selection.points)
        return []

    def _sync_extra_views_to_active_selection(self, *, force: bool = False) -> None:
        if not self.extra_views:
            return
        points = self._current_multi_view_selection()
        self._multi_view_syncing = True
        try:
            for view in self.extra_views.values():
                if force and points:
                    view.selection([])
                view.selection(points)
        finally:
            self._multi_view_syncing = False

    @staticmethod
    def _clean_multi_view_selection(value) -> list[int]:
        if value is None:
            return []
        return np.asarray(value, dtype=int).tolist()

    def _sync_extra_view_color(self, color_by: str | None) -> None:
        if not color_by:
            return
        color_map = self.categorical_color_maps.get(color_by)
        for view in self.extra_views.values():
            data = getattr(view, "_data", None)
            if data is None or color_by not in data.columns:
                continue
            if color_map is not None:
                view.color(by=color_by, map=color_map)
            else:
                view.color(by=color_by, map="magma")
            self._refresh_extra_view_non_spatial_points(view)

    def _color_extra_views_by_gene_expression(
        self,
        gene: str,
        display_gene: str,
        expr: np.ndarray,
        vmin: float,
        vmax: float,
    ) -> None:
        if not self.extra_views:
            return
        for label, view in self.extra_views.items():
            data = getattr(view, "_data", None)
            if data is not None and len(data) != len(expr):
                self._log(
                    f"[multi-view] skipped gene coloring for {label!r}: "
                    f"expected {len(expr)} rows, found {len(data)}"
                )
                continue
            view.color(
                by=expr,
                map=self._gene_expression_color_map(),
                norm=(vmin, vmax),
                labeling={
                    "minValue": "Low",
                    "maxValue": "High",
                    "variable": display_gene,
                },
            )
            self._refresh_extra_view_non_spatial_points(view)
        self._sync_extra_views_to_active_selection(force=True)

    def _configure_extra_views(self) -> None:
        for label, view in self.extra_views.items():
            try:
                view.width(self._extra_view_size)
                view.height(self._extra_view_size)
            except Exception:
                self.logger.exception("Failed to size extra view %r", label)

    def _refresh_extra_view_non_spatial_points(self, view: Scatter) -> None:
        try:
            view.widget.prevent_filter_reset = True
            view.widget.non_spatial_points_update = True
            view.widget.points = view.get_point_list()
        except Exception:
            self.logger.exception("Failed to refresh extra-view point encodings")

    def _lasso_selection_polygon_change_handler(self, change):
        scatter = self.scatter
        if change["new"] is None:
            self.lasso.polygon = None
        else:
            points = np.asarray(change["new"], dtype=float).tolist()
            points.append(points[0])
            self.lasso.polygon = Line(points, line_color=scatter.widget.color_selected)
        self._update_annotations()

    def _selection_handler(self, change):
        ctrl = self._ctrl
        if len(change["new"]) > 0:
            ctrl.selection_add.disabled = False
            ctrl.selection_name.disabled = False
            ctrl.selection_name.placeholder = "Name selection..."
            new_index = (
                self.selections.selections[-1].index + 1
                if self.selections.selections
                else 1
            )
            ctrl.selection_name.value = f"Selection {new_index}"
            if self.analysis_mode == "differential":
                ctrl.compute_diffexpr.disabled = False
        else:
            ctrl.selection_add.disabled = True
            ctrl.selection_name.disabled = True
            ctrl.selection_name.placeholder = "Select some points..."
            ctrl.selection_name.value = ""
            if self.analysis_mode == "differential":
                ctrl.compute_diffexpr.disabled = self.active_selection is None

    def _lasso_type_change_handler(self, change):
        ctrl = self._ctrl
        if change["new"] == "freeform":
            self.analysis_mode = "differential"
            ctrl.directional_controls_box.layout.display = "none"
            ctrl.diff_controls_box.layout.display = "flex"
            ctrl.complete_add.children = (ctrl.add_controls,)
            ctrl.compute_diffexpr.disabled = (
                self.scatter is None or len(self.scatter.selection()) == 0
            ) and (self.active_selection is None)
            self._clear_results_display(
                "<em>Differential mode. Make a freeform selection to compute DE.</em>"
            )
        else:
            self.analysis_mode = "directional"
            ctrl.directional_controls_box.layout.display = "flex"
            ctrl.diff_controls_box.layout.display = "none"
            ctrl.complete_add.children = (ctrl.add_controls,)
            ctrl.compute_predicates.style = "primary"
            ctrl.compute_predicates.description = "Compute Directional Search"
            ctrl.compute_predicates.on_click(self._compute_predicates_handler)
            self._clear_results_display(None)
        self._update_annotations()
        self._record_action(
            "change_lasso_type",
            f"Changed lasso type to {change['new']}",
            payload={"lasso_type": change["new"]},
        )

    def _lasso_brush_size_change_handler(self, change):
        self._record_action(
            "change_brush_size",
            f"Changed brush size to {change['new']}",
            payload={"brush_size": change["new"]},
        )

    def _compare_between_change_handler(self, change):
        self._record_action(
            "toggle_compare_between_selections",
            f"Compare between selections set to {bool(change['new'])}",
            payload={"compare_between_selections": bool(change["new"])},
        )

    def _diff_threshold_change_handler(self, change):
        self._record_action(
            "change_diffexpr_threshold",
            "Changed differential-expression threshold",
            payload={
                "field": change["owner"].description,
                "value": change["new"],
            },
        )

    def _color_by_change_handler(self, change):
        new = change["new"]
        if new in self.categorical_color_maps:
            self.scatter.color(by=new, map=self.categorical_color_maps[new])
        else:
            self.scatter.color(by=new, map="magma")
        self._sync_extra_view_color(new)
        self._hide_gene_expression_caption()
        self._refresh_scatter_non_spatial_points()
        self._record_action(
            "change_color_by",
            f"Changed color by to {new}",
            payload={"color_by": new},
        )

    # -- Session UI controls -------------------------------------------------

    def _refresh_history_options(self):
        ctrl = getattr(self, "_ctrl", None)
        if ctrl is None or not hasattr(ctrl, "history_dropdown"):
            return

        rows = _normalize_action_log(getattr(self, "action_log", []))
        options = [("No recorded actions", None)]
        if rows:
            options = [("Choose action...", None)] + [
                (f"{entry['index']}: {entry['label']}", entry["index"] - 1)
                for entry in rows
            ]

        current = ctrl.history_dropdown.value
        valid_values = {value for _, value in options}
        next_value = current if current in valid_values else None

        self._history_dropdown_paused = True
        try:
            ctrl.history_dropdown.options = options
            ctrl.history_dropdown.disabled = not rows
            ctrl.history_dropdown.value = next_value
        finally:
            self._history_dropdown_paused = False

    def _apply_history_action(self, step_index: int) -> str:
        document = _export_session(self)
        step = _apply_playback_step(self, document, int(step_index))
        self._refresh_history_options()
        return step.get("_status") or step.get("label") or step.get("type") or ""

    def _history_dropdown_change_handler(self, change):
        if getattr(self, "_history_dropdown_paused", False):
            return
        step_index = change["new"]
        if step_index is None:
            return
        try:
            message = self._apply_history_action(int(step_index))
            self._set_session_status(f"<em>{escape(message)}</em>")
        except Exception as exc:
            self.logger.exception("Failed to restore history action")
            self._set_session_status(
                f"<em>Could not restore action: {escape(str(exc))}</em>"
            )

    def _set_session_status(self, html: str) -> None:
        if not html:
            self._ctrl.session_status.value = ""
            self._ctrl.session_status.layout.display = "none"
            return
        self._ctrl.session_status.value = (
            '<div style="max-width:100%;overflow-wrap:anywhere;'
            f'white-space:normal;">{html}</div>'
        )
        self._ctrl.session_status.layout.display = "block"

    def _session_save_handler(self, event):
        filename = self._ctrl.session_filename.value.strip()
        path = Path(filename or "scsketch-session.scsketch.json")
        try:
            document = self.export_session(path)
            n_actions = len(document.get("action_log") or [])
            self._set_session_status(
                f"<em>Saved {n_actions} actions to "
                f"<code>{escape(str(path))}</code>.</em>"
            )
        except Exception as exc:
            self.logger.exception("Failed to save session")
            self._set_session_status(
                f"<em>Could not save session: {escape(str(exc))}</em>"
            )

    def _session_upload_handler(self, change):
        files = _uploaded_files(change["new"])
        if not files:
            return
        uploaded = files[0]
        try:
            content = _uploaded_content(uploaded)
            document = json.loads(content.decode("utf-8"))
            warnings = self.load_session(document)
            name = _uploaded_name(uploaded) or "uploaded session"
            n_actions = len(getattr(self, "action_log", []))
            warning_text = f" {len(warnings)} warnings." if warnings else ""
            self._set_session_status(
                f"<em>Loaded {n_actions} actions from "
                f"<code>{escape(name)}</code>.{warning_text}</em>"
            )
        except Exception as exc:
            self.logger.exception("Failed to load session")
            self._set_session_status(
                f"<em>Could not load session: {escape(str(exc))}</em>"
            )

    # -- Selection management ----------------------------------------------

    def _add_selection_element(self, selection: Selection):
        """Create and register a labelled row in the sidebar for *selection*."""
        from .widgets import Label
        scatter = self.scatter
        ctrl = self._ctrl
        hex_color = to_hex(selection.color)
        selection_label_widget = Label(name=selection.name, style={"background": hex_color})
        selection_remove = Button(
            description="",
            tooltip="Remove Selection",
            icon="trash",
            width=36,
            background=hex_color,
            rounded=["top-right", "bottom-right"],
        )
        element = GridBox(
            [selection_label_widget, selection_remove],
            layout=Layout(grid_template_columns="1fr 40px"),
        )

        def focus_handler(change):
            if change["new"]:
                scatter.zoom(to=selection.points, animation=500, padding=2)
                self.active_selection = selection
                self._sync_extra_views_to_active_selection()
                self._record_action(
                    "focus_selection",
                    f"Focused selection {selection.name}",
                    selection=selection.name,
                    payload={"n_cells": int(len(selection.points))},
                )
                if self.analysis_mode == "differential":
                    if selection.cached_diffexpr is not None:
                        self._show_diffexpr_results(
                            selection.cached_diffexpr,
                            selection.name,
                            selected_indices=np.asarray(selection.points, dtype=int),
                        )
                    else:
                        self._clear_results_display(
                            f"<em>No cached differential results for "
                            f"<b>{selection.name}</b> yet. "
                            f"Make a freeform selection to compute DE.</em>"
                        )
                else:
                    if selection.cached_results is not None:
                        self._show_directional_results([selection.cached_results])
                    else:
                        self._clear_results_display(
                            f"<em>No cached results for <b>{selection.name}</b> yet. "
                            f"Click <b>Compute Directional Search</b>.</em>"
                        )
            else:
                scatter.zoom(to=None, animation=500, padding=0)

        selection_label_widget.observe(focus_handler, names=["focus"])

        def remove_handler(change):
            self._archive_selection(selection)
            ctrl.selections_elements.children = [
                e for e in ctrl.selections_elements.children if e != element
            ]
            self.selections.selections = [
                s for s in self.selections.selections if s != selection
            ]
            if self.active_selection is selection:
                self.active_selection = (
                    self.selections.selections[-1] if self.selections.selections else None
                )
                self._sync_extra_views_to_active_selection()
                if self.analysis_mode == "differential":
                    if (
                        self.active_selection is None
                        or self.active_selection.cached_diffexpr is None
                    ):
                        self._clear_results_display(
                            None
                            if self.active_selection is None
                            else (
                                f"<em>No cached differential results for "
                                f"<b>{self.active_selection.name}</b> yet.</em>"
                            )
                        )
                    else:
                        self._show_diffexpr_results(
                            self.active_selection.cached_diffexpr,
                            self.active_selection.name,
                            selected_indices=np.asarray(
                                self.active_selection.points, dtype=int
                            ),
                        )
                else:
                    if (
                        self.active_selection is None
                        or self.active_selection.cached_results is None
                    ):
                        self._clear_results_display(
                            None
                            if self.active_selection is None
                            else (
                                f"<em>No cached results for "
                                f"<b>{self.active_selection.name}</b> yet.</em>"
                            )
                        )
                    else:
                        self._show_directional_results(
                            [self.active_selection.cached_results]
                        )
            self._update_annotations()
            ctrl.compute_predicates.disabled = len(self.selections.selections) == 0
            self._record_action(
                "remove_selection",
                f"Removed selection {selection.name}",
                selection=selection.name,
                payload={"n_cells": int(len(selection.points))},
            )

        selection_remove.on_click(remove_handler)
        ctrl.selections_elements.children = ctrl.selections_elements.children + (element,)

    def _add_subdivided_selections(self):
        scatter = self.scatter
        df = self.df
        ctrl = self._ctrl

        lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon, dtype=float)
        lasso_mid = lasso_polygon.shape[0] // 2
        lasso_part_one = lasso_polygon[:lasso_mid, :]
        lasso_part_two = lasso_polygon[lasso_mid:, :][::-1]

        n_split_points = ctrl.selection_num_subdivisions.value + 1
        sub_lassos_one = split_line_equidistant(lasso_part_one, n_split_points)
        sub_lassos_two = split_line_equidistant(lasso_part_two, n_split_points)

        base_name = ctrl.selection_name.value or f"Selection {len(self.selections.selections) + 1}"
        color_map = self.continuous_color_maps[ctrl.selection_num_subdivisions.value]

        for i, part_one in enumerate(sub_lassos_one):
            polygon = np.vstack((part_one, sub_lassos_two[i][::-1]))
            idxs = np.where(points_in_polygon(df[["x", "y"]].values, polygon))[0]
            pts_xy = df.iloc[idxs][["x", "y"]].values
            hull = ConvexHull(pts_xy)
            hull_pts = np.vstack((pts_xy[hull.vertices], pts_xy[hull.vertices[0]]))
            spine = (part_one + sub_lassos_two[i]) / 2

            poly_list = polygon.astype(float).tolist()
            poly_list.append(poly_list[0])

            sel = Selection(
                index=len(self.selections.selections) + 1,
                name=f"{base_name}.{i + 1}",
                points=idxs,
                color=color_map[i],
                lasso=Line(poly_list, line_color=color_map[i], line_width=2),
                hull=Line(hull_pts.astype(float).tolist(), line_color=color_map[i], line_width=2),
                path=spine,
            )
            self.selections.selections.append(sel)
            self._archive_selection(sel)
            self._add_selection_element(sel)

    def _add_selection(self):
        scatter = self.scatter
        df = self.df
        ctrl = self._ctrl

        idxs = scatter.selection()
        pts = df.iloc[idxs][["x", "y"]].values
        hull = ConvexHull(pts)
        hull_pts = np.vstack((pts[hull.vertices], pts[hull.vertices[0]]))

        color = self.available_colors.pop(0)

        spine = None
        if scatter.widget.lasso_type == "brush":
            lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon)
            if lasso_polygon.shape[0] >= 2:
                if lasso_polygon.shape[0] % 2 == 1:
                    lasso_polygon = lasso_polygon[:-1]
                mid = lasso_polygon.shape[0] // 2
                spine = (
                    lasso_polygon[:mid, :] + lasso_polygon[mid:, :][::-1]
                ) / 2

        name = ctrl.selection_name.value or f"Selection {len(self.selections.selections) + 1}"

        lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon, dtype=float).tolist()
        lasso_polygon.append(lasso_polygon[0])

        sel = Selection(
            index=len(self.selections.selections) + 1,
            name=name,
            points=idxs,
            color=color,
            lasso=Line(lasso_polygon, line_color=color, line_width=2),
            hull=Line(hull_pts.astype(float).tolist(), line_color=color, line_width=2),
            path=spine,
        )
        self.selections.selections.append(sel)
        self._archive_selection(sel)
        self._add_selection_element(sel)

    def _selection_add_handler(self, event):
        ctrl = self._ctrl
        try:
            self.lasso.polygon = None

            if self.scatter.widget.lasso_type == "brush" and ctrl.selection_subdivide.value:
                self._add_subdivided_selections()
            else:
                self._add_selection()

            self.active_selection = (
                self.selections.selections[-1] if self.selections.selections else None
            )
            if self.analysis_mode == "differential":
                if self.active_selection is not None and self._pending_diffexpr is not None:
                    pts = np.sort(np.asarray(self.active_selection.points, dtype=int))
                    pend_pts = np.asarray(
                        self._pending_diffexpr.get("points", []), dtype=int
                    )
                    if pend_pts.size == pts.size and np.array_equal(pend_pts, pts):
                        pending = self._pending_diffexpr.get("results", []) or []
                        for entry in pending:
                            if isinstance(entry, dict):
                                entry["direction"] = self.active_selection.name
                        self.active_selection.cached_diffexpr = pending
                if (
                    self.active_selection is not None
                    and self.active_selection.cached_diffexpr is not None
                ):
                    self._show_diffexpr_results(
                        self.active_selection.cached_diffexpr,
                        self.active_selection.name,
                        selected_indices=np.asarray(
                            self.active_selection.points, dtype=int
                        ),
                    )
                else:
                    self._clear_results_display(
                        None
                        if self.active_selection is None
                        else (
                            f"<em>No cached differential results for "
                            f"<b>{self.active_selection.name}</b> yet.</em>"
                        )
                    )
            else:
                self._clear_results_display(
                    None
                    if self.active_selection is None
                    else (
                        f"<em>No cached results for <b>{self.active_selection.name}</b> yet. "
                        f"Click <b>Compute Directional Search</b>.</em>"
                    )
                )

            ctrl.compute_predicates.disabled = False
            self.scatter.selection([])
            self._sync_extra_views_to_active_selection()
            self._update_annotations()

            if len(self.selections.selections) > 1:
                ctrl.directional_controls_box.children = (
                    ctrl.compute_predicates_between_selections,
                    ctrl.compute_predicates,
                    ctrl.directional_progress_box,
                )
            else:
                ctrl.directional_controls_box.children = (
                    ctrl.compute_predicates,
                    ctrl.directional_progress_box,
                )
            if self.active_selection is not None:
                self._archive_selection(self.active_selection)
                self._record_action(
                    "save_selection",
                    f"Saved selection {self.active_selection.name}",
                    selection=self.active_selection.name,
                    payload={
                        "n_cells": int(len(self.active_selection.points)),
                        "lasso_type": self.scatter.widget.lasso_type,
                    },
                )
        except Exception:
            self.logger.exception("Error in _selection_add_handler")

    def _clear_predicates(self, event):
        ctrl = self._ctrl
        ctrl.compute_predicates.style = "primary"
        ctrl.compute_predicates.description = "Compute Directional Search"
        ctrl.compute_predicates.on_click(self._compute_predicates_handler)
        self._clear_results_display(None)

        if len(self.selections.selections) > 1:
            ctrl.directional_controls_box.children = (
                ctrl.compute_predicates_between_selections,
                ctrl.compute_predicates,
                ctrl.directional_progress_box,
            )
        else:
            ctrl.directional_controls_box.children = (
                ctrl.compute_predicates,
                ctrl.directional_progress_box,
            )
        self._record_action("clear_results", "Cleared visible results")

    # -- Directional analysis -----------------------------------------------

    def _compute_directional_analysis(self, df: pd.DataFrame, selections: Selections):
        """Run per-selection directional (correlation) analysis."""
        if not selections.selections:
            return []

        results = []
        for selection in selections.selections:
            selected_indices = selection.points
            selected_embeddings = df.iloc[selected_indices][["x", "y"]].values
            if selected_embeddings.shape[0] < 2:
                results.append([])
                continue

            path = None if selection.path is None else np.asarray(selection.path, dtype=float)
            if path is not None and path.shape[0] >= 2:
                start_point = path[0]
                end_point = path[-1]
            else:
                start_point = selected_embeddings[0]
                end_point = selected_embeddings[-1]

            v = end_point - start_point
            nv = np.linalg.norm(v)
            if nv <= 1e-15:
                results.append([])
                continue
            v = v / nv
            projections = np.dot(selected_embeddings - start_point, v)

            X_sel = self.adata.X[selected_indices, :]
            batch_result = test_direction(X_sel, projections)
            rs = batch_result["correlation"].astype(float)
            ps = batch_result["p_value"].astype(float)
            genes = list(self.adata.var_names)
            n_new = len(ps)

            prev_len = (
                0 if self.batch_results is None else len(self.batch_results["p_value"])
            )
            p_values = (
                ps
                if self.batch_results is None
                else np.concatenate([self.batch_results["p_value"], ps])
            )

            online_results_new = lord_test(p_values, self.online_results, alpha=self.fdr_alpha)
            self.online_results = online_results_new
            self.batch_results = {"p_value": p_values}

            alpha_chunk = online_results_new["alpha_i"][prev_len: prev_len + n_new]
            R_chunk = online_results_new["R"][prev_len: prev_len + n_new]

            correlations = [
                {
                    "attribute": gene,
                    "interval": (float(rs[j]), float(ps[j])),
                    "quality": abs(float(rs[j])),
                    "alpha_i": float(alpha_chunk[j]),
                    "reject": True,
                    "direction": selection.name,
                }
                for j, gene in enumerate(genes)
                if bool(R_chunk[j])
            ]
            results.append(correlations)

        return results

    def _compute_predicates_handler(self, event):
        ctrl = self._ctrl
        self._clear_progress("Directional")
        try:
            if not self.selections.selections:
                return
            ctrl.compute_predicates.disabled = True
            ctrl.compute_predicates.description = "Computing Directional Analysis..."
            self._set_progress("Directional", 1, 4, "Preparing selection")

            if (
                ctrl.compute_predicates_between_selections is not None
                and ctrl.compute_predicates_between_selections.value
            ):
                sels_for_run = self.selections
            else:
                target = self.active_selection or self.selections.selections[-1]
                sels_for_run = Selections(selections=[target])

            self._set_progress("Directional", 2, 4, "Computing correlations / p-values")
            directional_results = self._compute_directional_analysis(self.df, sels_for_run)
            self._set_progress("Directional", 3, 4, "Caching results")
            for sel, res in zip(sels_for_run.selections, directional_results):
                sel.cached_results = res
                self._record_action(
                    "compute_directional",
                    f"Computed directional results for {sel.name}",
                    selection=sel.name,
                    payload={
                        "result_count": int(len(res)),
                        "compare_between_selections": bool(
                            ctrl.compute_predicates_between_selections is not None
                            and ctrl.compute_predicates_between_selections.value
                        ),
                        "fdr_alpha": float(self.fdr_alpha),
                    },
                )
            self._set_progress("Directional", 4, 4, "Rendering")
            self._show_directional_results(directional_results)
        except Exception:
            self.logger.exception("Error in _compute_predicates_handler")
        finally:
            ctrl.compute_predicates.disabled = False
            self._clear_progress("Directional")

    # -- Differential expression --------------------------------------------

    def _compute_diffexpr_handler(self, event):
        ctrl = self._ctrl
        self._clear_progress("DE")
        try:
            ctrl.compute_diffexpr.disabled = True
            ctrl.compute_diffexpr.description = "Computing DE..."
            self._set_progress("DE", 1, 4, "Preparing selection")

            self._de_engine.t_threshold = float(ctrl.diff_t_threshold.value)
            self._de_engine.p_threshold = float(ctrl.diff_p_threshold.value)

            if self.scatter is not None and len(self.scatter.selection()) > 0:
                sel = np.asarray(self.scatter.selection(), dtype=int)
                label = "Current selection"
                self._set_progress("DE", 2, 4, "Computing statistics")
                res = self._de_engine.compute(sel, label)
                self._pending_diffexpr = {"points": np.sort(np.unique(sel)), "results": res}
                self._set_progress("DE", 3, 4, "Formatting results")
                self._show_diffexpr_results(res, label, selected_indices=np.unique(sel))
                self._record_action(
                    "compute_diffexpr",
                    "Computed DE for current unsaved selection",
                    payload={
                        "n_cells": int(len(np.unique(sel))),
                        "result_count": int(len(res)),
                        "t_threshold": float(ctrl.diff_t_threshold.value),
                        "p_threshold": float(ctrl.diff_p_threshold.value),
                    },
                )
            elif self.active_selection is not None:
                label = self.active_selection.name
                self._set_progress("DE", 2, 4, "Computing statistics")
                res = self._de_engine.compute(self.active_selection.points, label)
                self.active_selection.cached_diffexpr = res
                self._set_progress("DE", 3, 4, "Formatting results")
                self._show_diffexpr_results(
                    res,
                    label,
                    selected_indices=np.asarray(self.active_selection.points, dtype=int),
                )
                self._record_action(
                    "compute_diffexpr",
                    f"Computed DE for {label}",
                    selection=label,
                    payload={
                        "n_cells": int(len(self.active_selection.points)),
                        "result_count": int(len(res)),
                        "t_threshold": float(ctrl.diff_t_threshold.value),
                        "p_threshold": float(ctrl.diff_p_threshold.value),
                    },
                )
            self._set_progress("DE", 4, 4, "Done")
        finally:
            ctrl.compute_diffexpr.disabled = False
            ctrl.compute_diffexpr.description = "Compute DE"
            self._clear_progress("DE")

    # -- Public API ---------------------------------------------------------

    def get_genes(self, sel_name: str = "Selection 1") -> pd.DataFrame:
        """Return a DataFrame of significant genes for the named selection.

        Parameters
        ----------
        sel_name:
            Name of the selection to retrieve results from.

        Returns
        -------
        DataFrame with columns ``gene``, ``correlation``, ``p-value``,
        and ``discovery_score``, sorted descending by correlation. Empty
        DataFrame if the selection is not found or has no cached results.
        """
        for sel in self.selections.selections:
            if sel.name == sel_name and sel.cached_results is not None:
                data = [
                    {
                        "gene": entry["attribute"],
                        "correlation": entry["interval"][0],
                        "p-value": entry["interval"][1],
                    }
                    for entry in sel.cached_results
                ]
                df = pd.DataFrame(data)
                if df.empty:
                    return df
                df["discovery_score"] = discovery_scores(df["p-value"])
                return df.sort_values(by="correlation", ascending=False)
        return pd.DataFrame()

    def get_diffexpr_genes(self, sel_name: str = "Selection 1") -> pd.DataFrame:
        """Return differential-expression genes for the named selection.

        Parameters
        ----------
        sel_name:
            Name of the selection to retrieve differential-expression results from.

        Returns
        -------
        DataFrame with columns ``gene``, ``t-statistic``, ``p-value``,
        ``discovery_score``, and ``selection``, sorted descending by absolute
        t-statistic. Empty DataFrame if the selection is not found or has no
        cached DE results.
        """
        for sel in self.selections.selections:
            if sel.name == sel_name and sel.cached_diffexpr is not None:
                data = [
                    {
                        "gene": entry["attribute"],
                        "t-statistic": entry["interval"][0],
                        "p-value": entry["interval"][1],
                        "selection": entry.get("direction", sel.name),
                    }
                    for entry in sel.cached_diffexpr
                ]
                df = pd.DataFrame(data)
                if df.empty:
                    return df
                df["discovery_score"] = discovery_scores(df["p-value"])
                df["_abs_t"] = df["t-statistic"].abs()
                return (
                    df.sort_values(by="_abs_t", ascending=False)
                    .drop(columns=["_abs_t"])
                    .reset_index(drop=True)
                )
        return pd.DataFrame()

    def get_de_genes(self, sel_name: str = "Selection 1") -> pd.DataFrame:
        """Alias for :meth:`get_diffexpr_genes`."""
        return self.get_diffexpr_genes(sel_name)

    def export_session(self, path: str | Path | None = None) -> dict:
        """Export selections, analysis parameters, and cached results.

        If ``path`` is provided, the session document is also written as JSON.
        """
        if path is None:
            return _export_session(self)
        self._record_action(
            "export_session",
            f"Exported session to {path}",
            payload={"path": str(path)},
        )
        return _write_session(self, path)

    def load_session(self, session: str | Path | dict) -> list[str]:
        """Load a session document or JSON file into this ScSketch instance.

        Returns
        -------
        list[str]
            Non-fatal compatibility warnings, such as dataset fingerprint
            mismatches.
        """
        document = (
            _read_session(session) if isinstance(session, (str, Path)) else session
        )
        warnings = _load_session(self, document)
        self._refresh_history_options()
        return warnings

    def get_action_log(self) -> pd.DataFrame:
        """Return the recorded user-action log as a DataFrame."""
        return pd.DataFrame(getattr(self, "action_log", []))

    def show_action_log(self):
        """Return a clickable table of recorded user actions."""
        def apply_entry(entry):
            step_index = max(0, int(entry.get("index") or 1) - 1)
            return self._apply_history_action(step_index)

        return _build_action_log_table(self, apply_entry=apply_entry)

    def show_session_player(self, session: str | Path | dict | None = None):
        """Return a playback UI for a saved session.

        If ``session`` is omitted, playback is built from the current sketch
        state. Otherwise, pass a session document or path to a session JSON file.
        """
        if session is None:
            document = _export_session(self)
        else:
            document = (
                _read_session(session) if isinstance(session, (str, Path)) else session
            )
        return _build_session_player(self, document)

    def show(self):
        """Display the ScSketch widget."""
        return self._ctrl.ui
