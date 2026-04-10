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

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List

import ipywidgets as ipyw
from ipywidgets import GridBox, Layout
from jscatter import Scatter, okabe_ito, Line
from jscatter.widgets import Button
from matplotlib.colors import to_hex
from scipy.spatial import ConvexHull

from anndata import AnnData

from ._logging import LogLevel, configure_logging
from ._data import build_embedding_df
from ._diffexpr import DiffExprEngine
from ._ui import UIControls, build_controls, set_analysis_progress, clear_analysis_progress
from ._results import show_directional_results, show_diffexpr_results, clear_results
from .utils import (
    Lasso,
    Selection,
    Selections,
    points_in_polygon,
    split_line_equidistant,
)
from .analysis import (
    test_direction,
    lord_test,
)


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
        """
        self.logger = configure_logging(verbosity)
        self.adata = adata
        self.metadata_cols = metadata_cols
        self.height = height
        self.background_color = background_color
        self.max_genes = max_genes
        self.fdr_alpha = fdr_alpha
        self.verbosity = verbosity

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
        self.scatter = Scatter(
            data=self.df,
            x="x",
            y="y",
            background_color=self.background_color,
            axes=False,
            height=self.height,
            color_by=self.color_by_default,
            color_map=self.color_map_default,
            tooltip=True,
            legend=False,
            tooltip_properties=[c for c in self.df.columns if c in self.meta_cols_present],
        )
        self.scatter.widget.color_selected = "#00dadb"

        # -- Build UI ------------------------------------------------------
        self._ctrl: UIControls = build_controls(
            scatter=self.scatter,
            available_metadata_cols=self.available_metadata_cols,
            all_gene_sorted=self.all_gene_sorted,
            color_by_default=self.color_by_default,
            max_genes=self.max_genes,
            df_columns=list(self.df.columns),
        )
        if self.scatter.widget.lasso_type == "freeform":
            self.analysis_mode = "differential"

        self._setup_handlers()

    # -- Logging shortcut --------------------------------------------------

    def _log(self, *args):
        self.logger.debug(" ".join(str(a) for a in args))

    # -- Progress helpers --------------------------------------------------

    def _set_progress(self, mode: str, step: int, total: int, message: str):
        set_analysis_progress(self._ctrl, mode, step, total, message)

    def _clear_progress(self, mode: str):
        clear_analysis_progress(self._ctrl, mode)

    # -- Results display ---------------------------------------------------

    def _show_directional_results(self, directional_results):
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
            on_results_cleared=self._clear_predicates,
            log=self._log,
        )

    def _show_diffexpr_results(self, diff_results, selection_label, *, selected_indices=None):
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
            log=self._log,
        )

    def _clear_results_display(self, message: str | None = None):
        ctrl = self._ctrl
        clear_results(
            ctrl.selections_predicates,
            ctrl.pathway_table_container,
            ctrl.reactome_diagram_container,
            message=message,
        )

    # -- Scatter annotations -----------------------------------------------

    def _update_annotations(self):
        scatter = self.scatter
        if scatter is None:
            return
        try:
            lasso_polygon = [] if self.lasso.polygon is None else [self.lasso.polygon]
            overlays = self.selections.all_hulls() + lasso_polygon
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

        ctrl = self._ctrl
        ctrl.selection_add.on_click(self._selection_add_handler)
        ctrl.compute_predicates.on_click(self._compute_predicates_handler)
        ctrl.compute_diffexpr.on_click(self._compute_diffexpr_handler)
        ctrl.color_by.observe(self._color_by_change_handler, names=["value"])

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

    def _color_by_change_handler(self, change):
        new = change["new"]
        if new in self.categorical_color_maps:
            self.scatter.color(by=new, map=self.categorical_color_maps[new])
        else:
            self.scatter.color(by=new, map="magma")

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

            poly_list = polygon.astype(float).tolist()
            poly_list.append(poly_list[0])

            sel = Selection(
                index=len(self.selections.selections) + 1,
                name=f"{base_name}.{i + 1}",
                points=idxs,
                color=color_map[i],
                lasso=Line(poly_list),
                hull=Line(hull_pts.astype(float).tolist(), line_color=color_map[i], line_width=2),
            )
            self.selections.selections.append(sel)
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
                spine = (lasso_polygon[:mid, :] + lasso_polygon[mid:, :]) / 2

        name = ctrl.selection_name.value or f"Selection {len(self.selections.selections) + 1}"

        lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon, dtype=float).tolist()
        lasso_polygon.append(lasso_polygon[0])

        sel = Selection(
            index=len(self.selections.selections) + 1,
            name=name,
            points=idxs,
            color=color,
            lasso=Line(lasso_polygon),
            hull=Line(hull_pts.astype(float).tolist(), line_color=color, line_width=2),
            path=spine,
        )
        self.selections.selections.append(sel)
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

            v = selected_embeddings[-1] - selected_embeddings[0]
            v = v / np.linalg.norm(v)
            start_point = selected_embeddings[0]
            projections = np.array(
                [np.dot(pt - start_point, v) for pt in selected_embeddings]
            )

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
        sorted descending by correlation.  Empty DataFrame if the selection
        is not found or has no cached results.
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
                return pd.DataFrame(data).sort_values(by="correlation", ascending=False)
        return pd.DataFrame()

    def show(self):
        """Display the ScSketch widget."""
        return self._ctrl.ui
