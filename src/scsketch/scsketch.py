"""Main ScSketch widget for interactive single-cell embedding exploration.

The UI composition and selection management patterns in this module are adapted from
the dimbridge notebook in jupyter-scatter by Fritz Lekschas:
https://github.com/flekschas/jupyter-scatter/blob/main/notebooks/dimbridge.ipynb
"""

import base64
import numpy as np
import pandas as pd
import requests
from itertools import cycle
from typing import Optional, List
import ipywidgets as ipyw

from IPython.display import display, HTML
from ipywidgets import Checkbox, Dropdown, GridBox, HBox, Layout, IntText, FloatText, Text, VBox
from jscatter import Scatter, glasbey_light, link, okabe_ito, Line
from jscatter.widgets import Button
from matplotlib.colors import to_hex
from scipy.spatial import ConvexHull
import scipy.sparse as sp
import scipy.stats as ss
from anndata import AnnData

from ._logging import LogLevel, configure_logging

from .utils import (
    Lasso,
    Selection,
    Selections,
    find_equidistant_vertices,
    points_in_polygon,
    split_line_at_points,
    split_line_equidistant,    
    create_selection, fetch_pathways
)
from .widgets import (
    GenePathwayWidget, CorrelationTable, PathwayTable,
    InteractiveSVG, Label, Div, GeneViolinPlot
)
from .analysis import (
    compute_directional_analysis, test_direction, lord_test
)

_PROGRESS_SPINNER_HTML = (
    '<svg width="14" height="14" viewBox="0 0 50 50" '
    'style="vertical-align:-0.125em;margin-right:6px;">'
    '<circle cx="25" cy="25" r="20" fill="none" stroke="#999" stroke-width="5" '
    'stroke-linecap="round" stroke-dasharray="60 70">'
    '<animateTransform attributeName="transform" type="rotate" from="0 25 25" to="360 25 25" '
    'dur="0.8s" repeatCount="indefinite"/>'
    "</circle></svg>"
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
                Directional analysis still uses all genes from `adata` regardless of this setting.
            fdr_alpha: False discovery rate alpha threshold for directional analysis
            verbosity: Logging verbosity level
        """
        self.logger = configure_logging(verbosity)
        self.adata = adata
        self.metadata_cols = metadata_cols
        self.color_by_default = color_by_default
        self.height = height
        self.background_color = background_color
        self.max_genes = max_genes
        self.fdr_alpha = fdr_alpha
        self.verbosity = verbosity

        self.df: pd.DataFrame | None = None
        self.available_metadata_cols: list[str] = []
        self.meta_cols_present: list[str] = []
        self.categorical_cols: list[str] = []
        self.categorical_color_maps: dict[str, dict] = {}
        self.color_map_default = None

        self.scatter: Scatter | None = None
        self.lasso = Lasso()
        self.selections = Selections()    

        self.available_colors = list(okabe_ito.copy())
        self.continuous_color_maps = [
            ["#00dadb", "#da00db"],
            ["#00dadb", "#a994dc", "#da00db"],
            ["#00dadb", "#8faddc", "#bd77dc", "#da00db"],
            ["#00dadb", "#7eb9dc", "#a994dc", "#c567dc", "#da00db"],
            ["#00dadb", "#72c0db", "#9aa3dc", "#b583dc", "#ca5cdb", "#da00db"],
            ["#00dadb", "#69c4db", "#8faddc", "#a994dc", "#bd77dc", "#cd54db", "#da00db"],
            [
                "#00dadb",
                "#62c7db",
                "#86b4dc",
                "#9e9fdc",
                "#b288dc",
                "#c16edc",
                "#cf4ddb",
                "#da00db",
            ],
            [
                "#00dadb",
                "#5ccadb",
                "#7eb9dc",
                "#96a7dc",
                "#a994dc",
                "#b87fdc",
                "#c567dc",
                "#d048db",
                "#da00db",
            ],
            [
                "#00dadb",
                "#57ccdb",
                "#78bddc",
                "#8faddc",
                "#a19ddc",
                "#b08bdc",
                "#bd77dc",
                "#c861db",
                "#d144db",
                "#da00db",
            ],
        ]

        self.batch_results = None
        self.online_results = None

        self.debug_out: ipyw.Output | None = None
        self.active_selection: Selection | None = None
        self.analysis_mode: str = "directional"  # "directional" | "differential"
        self._pending_diffexpr: dict | None = None
        self._diffexpr_global_stats: dict | None = None

        self.selection_name: Text | None = None
        self.selection_add: Button | None = None
        self.selection_subdivide: Checkbox | None = None
        self.selection_num_subdivisions: IntText | None = None
        self.selection_subdivide_wrapper: HBox | None = None
        self.selections_elements: VBox | None = None
        self.selections_predicates: VBox | None = None
        self.selections_predicates_wrapper: VBox | None = None
        self.compute_predicates: Button | None = None
        self.compute_predicates_between_selections: Checkbox | None = None
        self.compute_predicates_wrapper: VBox | None = None
        self.directional_controls_box: VBox | None = None
        self.directional_progress_box: VBox | None = None
        self.directional_progress_label: ipyw.HTML | None = None
        self.directional_progress_bar: ipyw.IntProgress | None = None
        self.diff_t_threshold: FloatText | None = None
        self.diff_p_threshold: FloatText | None = None
        self.compute_diffexpr: Button | None = None
        self.diff_controls_box: VBox | None = None
        self.diff_progress_box: VBox | None = None
        self.diff_progress_label: ipyw.HTML | None = None
        self.diff_progress_bar: ipyw.IntProgress | None = None

        self.add_controls: GridBox | None = None
        self.complete_add: VBox | None = None

        self.color_by: Dropdown | None = None
        self.plot_wrapper: VBox | None = None
        self.sidebar: GridBox | None = None

        self.pathway_table_container: VBox | None = None
        self.reactome_diagram_container: VBox | None = None

        self.ui: VBox | None = None

        self._build_df()
        self._build_scatter()
        self._build_ui()
        self._build_layout()
        self._setup_handlers()

    def _log(self, *args):
        self.logger.debug(" ".join(str(a) for a in args))

    def _set_analysis_progress(self, mode: str, step: int, total: int, message: str):
        if mode == "Directional":
            box = self.directional_progress_box
            label = self.directional_progress_label
            bar = self.directional_progress_bar
        elif mode == "DE":
            box = self.diff_progress_box
            label = self.diff_progress_label
            bar = self.diff_progress_bar
        else:
            return

        if box is None or label is None or bar is None:
            return

        total = int(max(1, total))
        bar.max = total
        bar.value = int(min(max(step, 0), total))
        label.value = f"{_PROGRESS_SPINNER_HTML}{mode} ({bar.value}/{total}): {message}"
        box.layout.display = "flex"

    def _clear_analysis_progress(self, mode: str):
        if mode == "Directional":
            box = self.directional_progress_box
            label = self.directional_progress_label
            bar = self.directional_progress_bar
        elif mode == "DE":
            box = self.diff_progress_box
            label = self.diff_progress_label
            bar = self.diff_progress_bar
        else:
            return

        if box is None or label is None or bar is None:
            return

        label.value = ""
        bar.value = 0
        box.layout.display = "none"

    def _brush_spine_from_polygon(self, lasso_polygon: np.ndarray) -> np.ndarray | None:
        lasso_polygon = np.asarray(lasso_polygon, dtype=float)
        if lasso_polygon.ndim != 2 or lasso_polygon.shape[1] != 2 or lasso_polygon.shape[0] < 4:
            return None
        if lasso_polygon.shape[0] % 2 == 1:
            lasso_polygon = lasso_polygon[:-1]
        mid = lasso_polygon.shape[0] // 2
        if mid < 2:
            return None
        return (lasso_polygon[:mid, :] + lasso_polygon[mid:, :]) / 2

    def _arrow_lines(
        self,
        start: np.ndarray,
        end: np.ndarray,
        *,
        color: str,
        width: int = 3,
    ) -> list[Line]:
        start = np.asarray(start, dtype=float).ravel()
        end = np.asarray(end, dtype=float).ravel()
        if start.size != 2 or end.size != 2:
            return []

        v = end - start
        norm = float(np.linalg.norm(v))
        if not np.isfinite(norm) or norm <= 0.0:
            return []

        u = v / norm
        p = np.array([-u[1], u[0]], dtype=float)

        # Arrowhead size in data units: scale by embedding extent but cap to segment length.
        if self.df is not None:
            xs = np.asarray(self.df["x"], dtype=float)
            ys = np.asarray(self.df["y"], dtype=float)
            diag = float(np.hypot(float(np.nanmax(xs) - np.nanmin(xs)), float(np.nanmax(ys) - np.nanmin(ys))))
            base = diag * 0.03 if np.isfinite(diag) and diag > 0 else norm * 0.15
        else:
            base = norm * 0.15

        head_len = float(min(norm * 0.25, max(base, 1e-6)))
        head_wid = float(head_len * 0.55)

        tip = end
        left = tip - head_len * u + head_wid * p
        right = tip - head_len * u - head_wid * p

        shaft = Line([start.tolist(), tip.tolist()], line_color=color, line_width=width)
        head1 = Line([tip.tolist(), left.tolist()], line_color=color, line_width=width)
        head2 = Line([tip.tolist(), right.tolist()], line_color=color, line_width=width)
        return [shaft, head1, head2]

    def _arrow_lines_from_spine(self, spine: np.ndarray, *, color: str, width: int = 3) -> list[Line]:
        spine = np.asarray(spine, dtype=float)
        if spine.ndim != 2 or spine.shape[1] != 2 or spine.shape[0] < 2:
            return []

        points = spine.astype(float).tolist()
        shaft = Line(points, line_color=color, line_width=width)

        # Orient the arrowhead using the last segment direction, but size it using the overall arrow scale.
        tip = spine[-1]
        seg = spine[-1] - spine[-2]
        seg_norm = float(np.linalg.norm(seg))
        if not np.isfinite(seg_norm) or seg_norm <= 0.0:
            return [shaft]

        u = seg / seg_norm
        p = np.array([-u[1], u[0]], dtype=float)

        overall = spine[-1] - spine[0]
        overall_norm = float(np.linalg.norm(overall))
        if not np.isfinite(overall_norm) or overall_norm <= 0.0:
            overall_norm = seg_norm

        if self.df is not None:
            xs = np.asarray(self.df["x"], dtype=float)
            ys = np.asarray(self.df["y"], dtype=float)
            diag = float(np.hypot(float(np.nanmax(xs) - np.nanmin(xs)), float(np.nanmax(ys) - np.nanmin(ys))))
            base = diag * 0.03 if np.isfinite(diag) and diag > 0 else overall_norm * 0.15
        else:
            base = overall_norm * 0.15

        head_len = float(min(overall_norm * 0.35, max(base, 1e-6)))
        head_wid = float(head_len * 0.55)

        left = tip - head_len * u + head_wid * p
        right = tip - head_len * u - head_wid * p

        head1 = Line([tip.tolist(), left.tolist()], line_color=color, line_width=width)
        head2 = Line([tip.tolist(), right.tolist()], line_color=color, line_width=width)
        return [shaft, head1, head2]

    def _build_df(self):
        umap_df = pd.DataFrame(
            self.adata.obsm["X_umap"],
            columns=["x", "y"],
            index=self.adata.obs_names,
        )

        if self.metadata_cols is not None:
            available_metadata_cols = [
                col for col in self.metadata_cols if col in self.adata.obs.columns
            ]
            if len(available_metadata_cols) > 0:
                metadata_df = self.adata.obs[available_metadata_cols].copy(deep=False)
                coerce_cols = [
                    col
                    for col in available_metadata_cols
                    if pd.api.types.is_object_dtype(metadata_df[col])
                    or pd.api.types.is_categorical_dtype(metadata_df[col])
                ]
                if len(coerce_cols) > 0:
                    metadata_df = metadata_df.assign(
                        **{col: metadata_df[col].astype(str) for col in coerce_cols}
                    )
            else:
                self.logger.info("No requested metadata columns found; continuing without metadata.")
                available_metadata_cols = []
                metadata_df = pd.DataFrame(index=self.adata.obs_names)
        else:
            available_metadata_cols = []
            metadata_df = pd.DataFrame(index=self.adata.obs_names)
            self.logger.info("No metadata passed; continuing with UMAP + gene expression only.")

        all_gene_sorted = sorted(list(self.adata.var_names), key=lambda s: s.lower())
        gene_subset = all_gene_sorted[: self.max_genes] if self.max_genes > 0 else []
        if len(gene_subset) > 0:
            subX = self.adata[:, gene_subset].X
            if sp.issparse(subX):
                gene_exp_df = pd.DataFrame.sparse.from_spmatrix(
                    subX, columns=gene_subset, index=self.adata.obs_names
                )
            else:
                gene_exp_df = pd.DataFrame(
                    np.asarray(subX), columns=gene_subset, index=self.adata.obs_names
                )
        else:
            gene_exp_df = pd.DataFrame(index=self.adata.obs_names)

        df = pd.concat([umap_df, metadata_df, gene_exp_df], axis=1)
        df = df.loc[:, ~df.columns.duplicated()]

        meta_cols_present = [c for c in (self.metadata_cols or []) if c in df.columns]
        categorical_cols = [
            c
            for c in meta_cols_present
            if df[c].dtype == "object" or df[c].nunique(dropna=False) <= 30
        ]
        for c in categorical_cols:
            df[c] = df[c].astype(str)

        def _sorted_cats(x: pd.Series):
            vals = pd.Index(x.unique().astype(str))
            if vals.str.fullmatch(r"\d+").all():
                return sorted(vals, key=lambda s: int(s))
            return sorted(vals, key=str)

        categorical_color_maps = {
            c: dict(zip(_sorted_cats(df[c]), cycle(glasbey_light[1:])))
            for c in categorical_cols
        }            

        if self.color_by_default in categorical_cols:
            color_map_default = categorical_color_maps[self.color_by_default]
        else:
            priority_cats = [
                c
                for c in [
                    "cell_type",
                    "celltype.l1",
                    "celltype.l2",
                    "cell_population",
                    "seurat_clusters",
                    "leiden",
                    "louvain",
                    "clusters",
                ]
                if c in categorical_cols
            ]

            if len(priority_cats) > 0:
                self.color_by_default = priority_cats[0]
                color_map_default = categorical_color_maps[self.color_by_default]
            else:
                fallback_cont = next(
                    (c for c in ["n_genes", "total_counts", "pct_counts_mt"] if c in df.columns),
                    None,
                )
                self.color_by_default = fallback_cont
                color_map_default = None

        self.df = df
        self.available_metadata_cols = available_metadata_cols
        self.all_gene_sorted = all_gene_sorted
        self.meta_cols_present = meta_cols_present
        self.categorical_cols = categorical_cols
        self.categorical_color_maps = categorical_color_maps
        self.color_map_default = color_map_default

    def _build_scatter(self):
        if self.df is None:
            raise RuntimeError("No data is available to build the scatter plot.")
    
        scatter = Scatter(
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
        scatter.widget.color_selected = "#00dadb"
        self.scatter = scatter

    def _build_ui(self):
        """Build the UI components."""

        self.selection_name = Text(value="", placeholder="Select some points…", disabled=True)
        self.selection_name.layout.width = "100%"

        self.selection_add = Button(
            description="",
            tooltip="Save Selection",
            disabled=True,
            icon="plus",
            width=36,
            rounded=["top-right", "bottom-right"],
        )

        self.selection_subdivide = Checkbox(value=False, description="Subdivide", indent=False)
        self.selection_num_subdivisions = IntText(value=5, min=2, max=10, step=1, description="Parts")
        self.selection_subdivide_wrapper = HBox([self.selection_subdivide, self.selection_num_subdivisions])

        self.selections_elements = VBox(layout=Layout(grid_gap="2px"))

        self.selections_predicates = VBox(layout=Layout(overflow_y="auto", height="100%", grid_gap="6px"))
        self.selections_predicates_wrapper = VBox([self.selections_predicates], layout=Layout(height="100%"))

        self.compute_predicates = Button(
            description="Compute Directional Search",
            style="primary",
            disabled=True,
            full_width=True,
        )
        self.compute_predicates_between_selections = Checkbox(
            value=False, description="Compare Between Selections", indent=False
        )
        self.directional_progress_label = ipyw.HTML("")
        self.directional_progress_bar = ipyw.IntProgress(
            value=0,
            min=0,
            max=5,
            description="",
            layout=Layout(width="100%"),
        )
        self.directional_progress_box = VBox(
            [self.directional_progress_label, self.directional_progress_bar],
            layout=Layout(display="none", width="100%", grid_gap="2px"),
        )
        self.directional_controls_box = VBox([self.compute_predicates, self.directional_progress_box])

        self.diff_t_threshold = FloatText(value=2.0, step=0.5, description="|T| ≥")
        self.diff_t_threshold.layout.width = "100%"
        self.diff_t_threshold.style = {"description_width": "50px"}

        self.diff_p_threshold = FloatText(value=0.05, step=0.01, description="p ≤")
        self.diff_p_threshold.layout.width = "100%"
        self.diff_p_threshold.style = {"description_width": "35px"}

        diff_thresholds = GridBox(
            [self.diff_t_threshold, self.diff_p_threshold],
            layout=Layout(grid_template_columns="1fr 1fr", width="100%", min_width="0px", grid_gap="6px"),
        )
        self.compute_diffexpr = Button(
            description="Compute DE",
            style="primary",
            disabled=True,
            full_width=True,
        )
        self.diff_progress_label = ipyw.HTML("")
        self.diff_progress_bar = ipyw.IntProgress(
            value=0,
            min=0,
            max=5,
            description="",
            layout=Layout(width="100%"),
        )
        self.diff_progress_box = VBox(
            [self.diff_progress_label, self.diff_progress_bar],
            layout=Layout(display="none", width="100%", grid_gap="2px"),
        )
        self.diff_controls_box = VBox(
            [
                diff_thresholds,
                self.compute_diffexpr,
                self.diff_progress_box,
            ],
            layout=Layout(display="none"),
        )

        self.compute_predicates_wrapper = VBox([self.directional_controls_box, self.diff_controls_box])

        self.add_controls = GridBox(
            [self.selection_name, self.selection_add],
            layout=Layout(grid_template_columns="1fr 40px"),
        )
        self.complete_add = VBox([self.add_controls], layout=Layout(grid_gap="4px"))

        scatter = self.scatter
        if scatter is not None and scatter.widget.lasso_type == "freeform":
            self.analysis_mode = "differential"
            self.directional_controls_box.layout.display = "none"
            self.diff_controls_box.layout.display = "flex"

    def _build_layout(self):
        """Build the entire UI layout."""

        if self.df is None or self.scatter is None:
            raise RuntimeError("DataFrame or scatter plot not initialized")
        
        df = self.df
        scatter = self.scatter

        _cat_priority = ["cell_population", "seurat_clusters", "leiden", "louvain","clusters"]
        meta_cols = list(self.available_metadata_cols or [])
        meta_in_df = [c for c in _cat_priority if c in meta_cols]
        meta_rest = sorted([c for c in meta_cols if c not in _cat_priority], key=lambda s: s.lower())
        meta_ordered = meta_in_df + meta_rest
        meta_opts = [(c.replace("_", " ").title(), c) for c in meta_ordered]

        _qc_order = ["n_genes", "total_counts", "pct_counts_mt"]
        obs_opts = [(c.replace("_", " ").title(), c) for c in _qc_order if c in df.columns]

        _gene_sorted = sorted(list(self.adata.var_names), key=lambda s: s.lower())
        gene_options = [(g, g) for g in _gene_sorted[: self.max_genes]] if self.max_genes > 0 else []
        dropdown_options = meta_opts + obs_opts + gene_options

        self.color_by = Dropdown(
            options=dropdown_options,
            value=self.color_by_default,
            description="Color By:",
        )
        self.plot_wrapper = VBox(
            [scatter.show(), self.color_by],
            layout=Layout(width="100%", min_width="0px"),
        )

        self.pathway_table_container = VBox(
            [],
            layout=Layout(
                overflow_y="auto",
                height="400px",
                min_width="0px",
                border="1px solid #ddd",
                padding="10px",
                display="none",
            ),
        )

        self.reactome_diagram_container = VBox(
            [],
            layout=Layout(
                overflow="auto",
                min_height="0px",
                width="100%",
                max_width="100%",
                min_width="0px",
                padding="10px",
                display="none",
                border="1px solid #ccc",
            ),
        )

        self.sidebar = GridBox(
            [
                self.complete_add,
                self.selections_elements,
                self.selections_predicates_wrapper,
                self.compute_predicates_wrapper,
            ],
            layout=Layout(
                grid_template_rows="min-content max-content 1fr min-content",
                overflow_y="hidden",
                height="100%",
                min_width="0px",
                grid_gap="4px",
            ),
        )

        self.pathway_table_container.layout = Layout(
            overflow_y="auto",
            max_height="400px",
            border="1px solid #ddd",
            padding="10px",
            display="none",
        )

        top_layout = GridBox(
            [self.plot_wrapper, self.sidebar, self.pathway_table_container],
            layout=Layout(grid_template_columns="2fr 1fr 1fr", grid_gap="10px", height="auto"),
        )

        display(
            HTML(
                """
<style>
.jp-OutputArea-output, .jp-Cell-outputArea, .jp-Notebook {
    overflow: auto !important;
    max-height: none !important;
}
</style>
"""
            )
        )

        combined_gene_pathway_panel = GridBox(
            [
                VBox([self.sidebar], layout=Layout(overflow_y="auto", height="100%", min_width="0px")),
                VBox(
                    [self.pathway_table_container],
                    layout=Layout(overflow_y="auto", height="100%", min_width="0px"),
                ),
            ],
            layout=Layout(
                grid_template_columns="2fr 3fr",
                grid_gap="5px",
                align_items="flex-start",
                height="760px",
                min_width="0px",
            ),
        )

        top_layout_updated = GridBox(
            [self.plot_wrapper, combined_gene_pathway_panel],
            layout=Layout(
                grid_template_columns="2.5fr 2.5fr",
                grid_gap="10px",
                height="auto",
                align_items="flex-start",
                width="100%",
            ),
        )

        self.ui = VBox(
            [top_layout_updated, VBox([self.reactome_diagram_container], layout=Layout(width="100%"))],
            layout=Layout(
                grid_gap="10px",
                width="100%",
                overflow_y="auto",
                align_items="flex-start",
            ),
        )

    def _setup_handlers(self):
        """Set up all event handlers."""

        scatter = self.scatter
        if scatter is None:
            raise RuntimeError("scatter not initialized")

        scatter.widget.observe(self._lasso_selection_polygon_change_handler, names=["lasso_selection_polygon"])
        scatter.widget.observe(self._selection_handler, names=["selection"])
        scatter.widget.observe(self._lasso_type_change_handler, names=["lasso_type"])

        if self.selection_add is None or self.compute_predicates is None or self.color_by is None:
            raise RuntimeError("UI not initialized")

        self.selection_add.on_click(self._selection_add_handler)
        self.compute_predicates.on_click(self._compute_predicates_handler)
        if self.compute_diffexpr is not None:
            self.compute_diffexpr.on_click(self._compute_diffexpr_handler)
        self.color_by.observe(self._color_by_change_handler, names=["value"])

    def _update_annotations(self):
        """Update scatter plot annotations."""
        scatter = self.scatter
        if scatter is None:
            return
        try:
            lasso_polygon = [] if self.lasso.polygon is None else [self.lasso.polygon]
            overlays = self.selections.all_hulls() + lasso_polygon

            arrow_overlays: list[Line] = []
            if self.analysis_mode == "directional" and self.active_selection is not None and self.active_selection.path is not None:
                spine = np.asarray(self.active_selection.path, dtype=float)
                if spine.ndim == 2 and spine.shape[0] >= 2 and spine.shape[1] == 2:
                    arrow_overlays.extend(self._arrow_lines_from_spine(spine, color=to_hex(self.active_selection.color), width=3))

            if scatter.widget.lasso_type == "brush" and scatter.widget.lasso_selection_polygon is not None:
                spine = self._brush_spine_from_polygon(np.asarray(scatter.widget.lasso_selection_polygon))
                if spine is not None and spine.shape[0] >= 2:
                    arrow_overlays.extend(self._arrow_lines_from_spine(spine, color=scatter.widget.color_selected, width=3))

            overlays = overlays + arrow_overlays
            scatter.annotations(overlays)
        except Exception:
            self.logger.exception("Failed to update annotations")

    def _lasso_selection_polygon_change_handler(self, change):
        scatter = self.scatter
        if scatter is None:
            return
        if change["new"] is None:
            self.lasso.polygon = None
        else:
            points = np.asarray(change["new"], dtype=float).tolist()
            points.append(points[0])
            self.lasso.polygon = Line(points, line_color=scatter.widget.color_selected)
        self._update_annotations()

    def _add_selection_element(self, selection: Selection):
        scatter = self.scatter
        if scatter is None or self.selections_elements is None or self.compute_predicates is None:
            return

        hex_color = to_hex(selection.color)
        selection_name = Label(name=selection.name, style={"background": hex_color})
        selection_remove = Button(
            description="",
            tooltip="Remove Selection",
            icon="trash",
            width=36,
            background=hex_color,
            rounded=["top-right", "bottom-right"],
        )
        element = GridBox([selection_name, selection_remove], layout=Layout(grid_template_columns="1fr 40px"))

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
                            f"<em>No cached differential results for <b>{selection.name}</b> yet. "
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

        selection_name.observe(focus_handler, names=["focus"])

        def remove_handler(change):
            self.selections_elements.children = [e for e in self.selections_elements.children if e != element]
            self.selections.selections = [s for s in self.selections.selections if s != selection]
            if self.active_selection is selection:
                self.active_selection = self.selections.selections[-1] if self.selections.selections else None
                if self.analysis_mode == "differential":
                    if self.active_selection is None or self.active_selection.cached_diffexpr is None:
                        self._clear_results_display(
                            None
                            if self.active_selection is None
                            else (
                                f"<em>No cached differential results for <b>{self.active_selection.name}</b> yet. "
                                f"Make a freeform selection to compute DE.</em>"
                            )
                        )
                    else:
                        self._show_diffexpr_results(
                            self.active_selection.cached_diffexpr,
                            self.active_selection.name,
                            selected_indices=np.asarray(self.active_selection.points, dtype=int),
                        )
                else:
                    if self.active_selection is None or self.active_selection.cached_results is None:
                        self._clear_results_display(
                            None
                            if self.active_selection is None
                            else (
                                f"<em>No cached results for <b>{self.active_selection.name}</b> yet. "
                                f"Click <b>Compute Directional Search</b>.</em>"
                            )
                        )
                    else:
                        self._show_directional_results([self.active_selection.cached_results])
            self._update_annotations()
            self.compute_predicates.disabled = len(self.selections.selections) == 0

        selection_remove.on_click(remove_handler)
        self.selections_elements.children = self.selections_elements.children + (element,)

    def _add_subdivided_selections(self):
        if self.scatter is None or self.df is None:
            return
        if self.selection_num_subdivisions is None or self.selection_name is None:
            return

        scatter = self.scatter
        df = self.df
        selection_num_subdivisions = self.selection_num_subdivisions
        selection_name = self.selection_name

        lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon, dtype=float)
        lasso_mid = int(lasso_polygon.shape[0] / 2)
        lasso_part_one = lasso_polygon[:lasso_mid, :]
        lasso_part_two = lasso_polygon[lasso_mid:, :][::-1]

        n_split_points = selection_num_subdivisions.value + 1
        sub_lassos_part_one = split_line_equidistant(lasso_part_one, n_split_points)
        sub_lassos_part_two = split_line_equidistant(lasso_part_two, n_split_points)

        base_name = selection_name.value
        if len(base_name) == 0:
            base_name = f"Selection {len(self.selections.selections) + 1}"

        color_map = self.continuous_color_maps[selection_num_subdivisions.value]

        for i, part_one in enumerate(sub_lassos_part_one):
            polygon = np.vstack((part_one, sub_lassos_part_two[i][::-1]))
            idxs = np.where(points_in_polygon(df[["x", "y"]].values, polygon))[0]
            points = df.iloc[idxs][["x", "y"]].values
            hull = ConvexHull(points)
            hull_points = np.vstack((points[hull.vertices], points[hull.vertices[0]]))
            color = color_map[i]
            name = f"{base_name}.{i + 1}"

            lasso_polygon_list = polygon.astype(float).tolist()
            lasso_polygon_list.append(lasso_polygon_list[0])
            hull_points_list = hull_points.astype(float).tolist()

            selection = Selection(
                index=len(self.selections.selections) + 1,
                name=name,
                points=idxs,
                color=color,
                lasso=Line(lasso_polygon_list),
                hull=Line(hull_points_list, line_color=color, line_width=2),
            )
            self.selections.selections.append(selection)
            self._add_selection_element(selection)

    def _add_selection(self):
        if self.scatter is None or self.df is None:
            return
        if self.selection_name is None:
            return

        scatter = self.scatter
        df = self.df

        idxs = scatter.selection()
        points = df.iloc[idxs][["x", "y"]].values
        hull = ConvexHull(points)
        hull_points = np.vstack((points[hull.vertices], points[hull.vertices[0]]))
        hull_points_py = hull_points.astype(float).tolist()
        color = self.available_colors.pop(0)

        spine = None
        if scatter.widget.lasso_type == "brush":
            lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon)
            if lasso_polygon.shape[0] >= 2:
                if lasso_polygon.shape[0] % 2 == 1:
                    lasso_polygon = lasso_polygon[:-1]
                mid = lasso_polygon.shape[0] // 2
                spine = (lasso_polygon[:mid, :] + lasso_polygon[mid:, :]) / 2

        name = self.selection_name.value
        if len(name) == 0:
            name = f"Selection {len(self.selections.selections) + 1}"

        lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon, dtype=float)
        lasso_polygon_list = lasso_polygon.astype(float).tolist()
        lasso_polygon_list.append(lasso_polygon_list[0])

        selection = Selection(
            index=len(self.selections.selections) + 1,
            name=name,
            points=idxs,
            color=color,
            lasso=Line(lasso_polygon_list),
            hull=Line(hull_points_py, line_color=color, line_width=2),
            path=spine,
        )
        self.selections.selections.append(selection)
        self._add_selection_element(selection)

    def _selection_add_handler(self, event):
        if self.scatter is None or self.compute_predicates is None or self.compute_predicates_wrapper is None:
            return
        if self.compute_predicates_between_selections is None:
            return
        if self.selection_subdivide is None:
            return
        if self.directional_controls_box is None:
            return

        scatter = self.scatter
        selection_subdivide = self.selection_subdivide

        try:
            self.lasso.polygon = None

            if scatter.widget.lasso_type == "brush" and selection_subdivide.value:
                self._add_subdivided_selections()
            else:
                self._add_selection()

            self.active_selection = self.selections.selections[-1] if self.selections.selections else None
            if self.analysis_mode == "differential":
                if self.active_selection is not None and self._pending_diffexpr is not None:
                    pts = np.sort(np.asarray(self.active_selection.points, dtype=int))
                    pend_pts = np.asarray(self._pending_diffexpr.get("points", []), dtype=int)
                    if pend_pts.size == pts.size and np.array_equal(pend_pts, pts):
                        pending = self._pending_diffexpr.get("results", []) or []
                        for entry in pending:
                            if isinstance(entry, dict):
                                entry["direction"] = self.active_selection.name
                        self.active_selection.cached_diffexpr = pending
                if self.active_selection is not None and self.active_selection.cached_diffexpr is not None:
                    self._show_diffexpr_results(
                        self.active_selection.cached_diffexpr,
                        self.active_selection.name,
                        selected_indices=np.asarray(self.active_selection.points, dtype=int),
                    )
                elif self.active_selection is not None:
                    self._clear_results_display(
                        f"<em>No cached differential results for <b>{self.active_selection.name}</b> yet. "
                        f"Make a freeform selection to compute DE.</em>"
                    )
                else:
                    self._clear_results_display(None)
            else:
                if self.active_selection is not None:
                    self._clear_results_display(
                        f"<em>No cached results for <b>{self.active_selection.name}</b> yet. "
                        f"Click <b>Compute Directional Search</b>.</em>"
                    )
                else:
                    self._clear_results_display(None)

            self.compute_predicates.disabled = False
            scatter.selection([])
            self._update_annotations()

            if len(self.selections.selections) > 1:
                self.directional_controls_box.children = (
                    self.compute_predicates_between_selections,
                    self.compute_predicates,
                    self.directional_progress_box,
                )
            else:
                self.directional_controls_box.children = (self.compute_predicates, self.directional_progress_box)
        except Exception:
            self.logger.exception("Error updating compute_predicates UI state")

    def _selection_handler(self, change):
        if self.selection_add is None or self.selection_name is None:
            return

        if len(change["new"]) > 0:
            self.selection_add.disabled = False
            self.selection_name.disabled = False
            self.selection_name.placeholder = "Name selection…"
            new_index = 1
            if len(self.selections.selections) > 0:
                new_index = self.selections.selections[-1].index + 1
            self.selection_name.value = f"Selection {new_index}"
            if self.analysis_mode == "differential":
                if self.compute_diffexpr is not None:
                    self.compute_diffexpr.disabled = False
        else:
            self.selection_add.disabled = True
            self.selection_name.disabled = True
            self.selection_name.placeholder = "Select some points…"
            self.selection_name.value = ""
            if self.analysis_mode == "differential" and self.compute_diffexpr is not None:
                self.compute_diffexpr.disabled = self.active_selection is None

    def _clear_predicates(self, event):
        if self.compute_predicates is None or self.compute_predicates_wrapper is None:
            return
        if self.directional_controls_box is None:
            return

        self.compute_predicates.style = "primary"
        self.compute_predicates.description = "Compute Directional Search"
        self.compute_predicates.on_click(self._compute_predicates_handler)

        self._clear_results_display(None)

        if self.compute_predicates_between_selections is not None and len(self.selections.selections) > 1:
            self.directional_controls_box.children = (
                self.compute_predicates_between_selections,
                self.compute_predicates,
                self.directional_progress_box,
            )
        else:
            self.directional_controls_box.children = (self.compute_predicates, self.directional_progress_box)

    def _fetch_pathways(self, gene):
        url = f"https://reactome.org/ContentService/data/mapping/UniProt/{gene}/pathways?species=9606"
        try:
            response = requests.get(url)
            response.raise_for_status()
            pathways = response.json()
            return [{"Pathway": entry["displayName"], "stId": entry["stId"]} for entry in pathways]
        except requests.exceptions.RequestException as e:
            self.logger.warning("Error fetching Reactome pathways for %s: %s", gene, e)
            return []

    def _show_directional_results(self, directional_results):
        if self.compute_predicates is None or self.selections_predicates is None:
            return
        if self.pathway_table_container is None or self.reactome_diagram_container is None:
            return

        compute_predicates = self.compute_predicates
        selections_predicates = self.selections_predicates
        pathway_table_container = self.pathway_table_container
        reactome_diagram_container = self.reactome_diagram_container
        log = self._log
        adata = self.adata

        compute_predicates.style = ""
        compute_predicates.description = "Clear Results"
        compute_predicates.on_click(self._clear_predicates)

        all_results = []
        for i, result in enumerate(directional_results):
            for entry in result:
                if "reject" in entry and not entry["reject"]:
                    continue
                all_results.append(
                    {
                        "Gene": entry["attribute"],
                        "R": float(np.round(entry["interval"][0], 4)),
                        "p": f"{entry['interval'][1]:.3e}",
                        "Selection": entry.get("direction", f"Selection {i+1}"),
                    }
                )

        results_df = pd.DataFrame(all_results)
        results_df = results_df.dropna(subset=["R", "p"])
        results_df = results_df.sort_values(by="R", ascending=False).reset_index(drop=True)

        gene_table_widget = CorrelationTable(data=results_df.to_dict(orient="records"))
        pathway_table_widget = PathwayTable(data=[])

        pathway_table_container.layout.display = "none"
        reactome_diagram_container.layout.display = "none"

        from .widgets import GeneProjectionPlot
        from ipywidgets import HTML as _HTML

        gene_proj_plot = GeneProjectionPlot()
        pathway_msg = _HTML("")

        plot_box = VBox(
            [gene_proj_plot],
            layout=Layout(
                flex="0 0 auto",
                height="300px",
                max_height="300px",
                min_height="300px",
                padding="0px",
                overflow="visible",
                align_self="stretch",
                display="none",
            ),
        )

        pathway_box = VBox(
            [pathway_msg, pathway_table_widget],
            layout=Layout(flex="1 1 auto", overflow="auto"),
        )

        pathway_table_container.children = [plot_box, pathway_box]
        pathway_table_container.layout = Layout(
            display="flex",
            flex_direction="column",
            height="420px",
            max_height="420px",
            overflow="visible",
        )

        def on_gene_click(change):
            import traceback

            gene = change["new"]
            log(f"[UI] gene clicked: {gene!r}")
            try:
                pathways = self._fetch_pathways(gene)
                pathway_table_widget.data = pathways
                pathway_table_container.layout.display = "block"
                pathway_msg.value = (
                    f"<em>No Reactome pathways found for <b>{gene}</b>.</em>"
                    if len(pathways) == 0
                    else f"<b>Reactome pathways for {gene}</b>"
                )
                reactome_diagram_container.layout.display = "none"

                df = self.df
                selections = self.selections

                if len(selections.selections) == 0:
                    log("No selections available for gene projection plot.")
                    return

                sel = self.active_selection if self.active_selection is not None else selections.selections[-1]
                selected_indices = sel.points
                log(f"[plot] selection name={sel.name!r} n={len(selected_indices)}")

                X = df[["x", "y"]].to_numpy(dtype=float)
                Ux = X[:, 0] - X[:, 0].min()
                Uy = X[:, 1] - X[:, 1].min()
                px = float(np.ptp(Ux)) or 1.0
                py = float(np.ptp(Uy)) or 1.0
                U = np.c_[Ux / px, Uy / py]
                U_sel = U[selected_indices]

                dirv = U_sel[-1] - U_sel[0]
                L2 = float(np.dot(dirv, dirv))
                if L2 <= 1e-15:
                    pts = df.iloc[selected_indices][["x", "y"]].to_numpy(dtype=float)
                    v = pts[-1] - pts[0]
                    nv = np.linalg.norm(v)
                    if nv <= 1e-15:
                        log("Degenerate selection; skipping plot.")
                        return
                    v = v / nv
                    start = pts[0]
                    proj = np.dot(pts - start, v)
                    span = float(proj.max() - proj.min()) or 1.0
                    proj = (proj - proj.min()) / span
                else:
                    celv = U_sel - U_sel[0]
                    proj = (celv @ dirv) / (L2 + 1e-12)
                    proj = np.clip(np.nan_to_num(proj, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)

                if gene in df.columns:
                    expr = pd.to_numeric(df.iloc[selected_indices][gene], errors="coerce").to_numpy(dtype=float)
                    log(f"[plot] gene {gene!r} found in df columns (n={len(expr)})")
                else:
                    sub = adata[selected_indices, [gene]].X
                    if hasattr(sub, "toarray"):
                        sub = sub.toarray()
                    expr = np.asarray(sub).ravel().astype(float)
                    log(
                        f"[plot] gene {gene!r} NOT in df slice; loaded from adata for selection (n={len(expr)})"
                    )
                if not np.isfinite(expr).all():
                    finite = np.isfinite(expr)
                    if finite.any():
                        med = float(np.median(expr[finite]))
                        expr = np.where(finite, expr, med)
                    else:
                        log(f"Expression for {gene} has no finite values in the selection.")
                        return

                dx = proj.max() - proj.min()
                dy = expr.max() - expr.min()
                proj = np.linspace(0, 1, len(proj)) if dx <= 1e-12 else (proj - proj.min()) / (dx + 1e-12)
                expr = np.zeros_like(expr) if dy <= 1e-12 else (expr - expr.min()) / (dy + 1e-12)

                plot_data = [{"projection": float(p), "expression": float(e)} for p, e in zip(proj, expr)]
                gene_proj_plot.data = plot_data
                plot_box.layout.display = "block"
                gene_proj_plot.gene = gene

                log(
                    f"[plot] {gene}: n={len(plot_data)} pts "
                    f"proj=[{proj.min():.3f},{proj.max():.3f}] expr=[{expr.min():.3f},{expr.max():.3f}]"
                )
                log(f"[plot] first5={plot_data[:5]}")
            except Exception:
                self.logger.exception("Error handling gene click")

        interactive_svg_widget = InteractiveSVG()

        def on_pathway_click(change):
            pathway_id = change["new"]
            if not pathway_id:
                return

            svg_url = f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.svg"
            try:
                response = requests.get(svg_url)
                response.raise_for_status()
                svg_text = response.text.strip()
                if len(svg_text) < 50:
                    log("Empty SVG returned from Reactome")
                    return

                import base64

                svg_base64 = base64.b64encode(svg_text.encode("utf-8")).decode("utf-8")
                interactive_svg_widget.svg_content = svg_base64

                reactome_diagram_container.children = [interactive_svg_widget]
                reactome_diagram_container.layout.display = "block"

                log(f"[SVG] Loaded Reactome diagram for {pathway_id}")
            except Exception as e:
                self.logger.warning("Error fetching SVG diagram: %s", e)

        gene_table_widget.observe(on_gene_click, names=["selected_gene"])
        self.logger.debug("[UI] gene click observer attached")
        pathway_table_widget.observe(on_pathway_click, names=["selected_pathway"])
        selections_predicates.children = [gene_table_widget]
        log("Showing directional results...")

    def _compute_directional_analysis(self, df, selections: Selections):
        if len(selections.selections) == 0:
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
            projections = np.array([np.dot(pt - start_point, v) for pt in selected_embeddings])

            X_sel = self.adata.X[selected_indices, :]

            batch_result_new = test_direction(X_sel, projections)
            rs = batch_result_new["correlation"].astype(float)
            ps = batch_result_new["p_value"].astype(float)
            genes = list(self.adata.var_names)
            n_new = len(ps)

            prev_len = 0 if (self.batch_results is None) else len(self.batch_results["p_value"])
            p_values = ps if (self.batch_results is None) else np.concatenate([self.batch_results["p_value"], ps])

            online_results_new = lord_test(p_values, self.online_results, alpha=self.fdr_alpha)
            self.online_results = online_results_new
            self.batch_results = {"p_value": p_values}

            alpha_chunk = online_results_new["alpha_i"][prev_len : prev_len + n_new]
            R_chunk = online_results_new["R"][prev_len : prev_len + n_new]

            correlations = []
            for j, gene in enumerate(genes):
                if not bool(R_chunk[j]):
                    continue

                r = float(rs[j])
                p = float(ps[j])
                a = float(alpha_chunk[j])
                correlations.append(
                    {
                        "attribute": gene,
                        "interval": (r, p),
                        "quality": abs(r),
                        "alpha_i": a,
                        "reject": True,
                        "direction": selection.name,
                    }
                )
            results.append(correlations)

        return results

    def _compute_predicates_handler(self, event):
        if self.compute_predicates is None:
            return
        self._clear_analysis_progress("Directional")
        try:
            if len(self.selections.selections) == 0:
                return

            self.compute_predicates.disabled = True
            self.compute_predicates.description = "Computing Directional Analysis…"

            self._set_analysis_progress("Directional", 1, 4, "Preparing selection")

            if self.compute_predicates_between_selections is not None and self.compute_predicates_between_selections.value:
                sels_for_run = self.selections
            else:
                target = self.active_selection if self.active_selection is not None else self.selections.selections[-1]
                last_only = Selections(selections=[target])
                sels_for_run = last_only

            self._set_analysis_progress("Directional", 2, 4, "Computing correlations / p-values")
            directional_results = self._compute_directional_analysis(self.df, sels_for_run)
            self._set_analysis_progress("Directional", 3, 4, "Caching results")
            for sel, res in zip(sels_for_run.selections, directional_results):
                sel.cached_results = res
            self._set_analysis_progress("Directional", 4, 4, "Rendering")
            self._show_directional_results(directional_results)
        except Exception:
            import traceback

            traceback.print_exc()
        finally:
            self.compute_predicates.disabled = False
            self._clear_analysis_progress("Directional")

    def _lasso_type_change_handler(self, change):
        if self.complete_add is None or self.add_controls is None or self.selection_subdivide_wrapper is None:
            return
        if change["new"] == "freeform":
            self.analysis_mode = "differential"
            if self.directional_controls_box is not None:
                self.directional_controls_box.layout.display = "none"
            if self.diff_controls_box is not None:
                self.diff_controls_box.layout.display = "flex"
            self.complete_add.children = (self.add_controls,)
            if self.compute_diffexpr is not None:
                self.compute_diffexpr.disabled = (self.scatter is None or len(self.scatter.selection()) == 0) and (self.active_selection is None)
            self._clear_results_display("<em>Differential mode. Make a freeform selection to compute DE.</em>")
        elif change["new"] == "brush":
            self.analysis_mode = "directional"
            if self.directional_controls_box is not None:
                self.directional_controls_box.layout.display = "flex"
            if self.diff_controls_box is not None:
                self.diff_controls_box.layout.display = "none"
            self.complete_add.children = (self.add_controls,)
            if self.compute_predicates is not None:
                self.compute_predicates.style = "primary"
                self.compute_predicates.description = "Compute Directional Search"
                self.compute_predicates.on_click(self._compute_predicates_handler)
            self._clear_results_display(None)
        else:
            self.analysis_mode = "directional"
            if self.directional_controls_box is not None:
                self.directional_controls_box.layout.display = "flex"
            if self.diff_controls_box is not None:
                self.diff_controls_box.layout.display = "none"
            self.complete_add.children = (self.add_controls,)
            if self.compute_predicates is not None:
                self.compute_predicates.style = "primary"
                self.compute_predicates.description = "Compute Directional Search"
                self.compute_predicates.on_click(self._compute_predicates_handler)
            self._clear_results_display(None)
        self._update_annotations()

    def _compute_diffexpr_handler(self, event):
        if self.scatter is None:
            return
        if self.compute_diffexpr is None:
            return
        self._clear_analysis_progress("DE")
        try:
            self.compute_diffexpr.disabled = True
            self.compute_diffexpr.description = "Computing DE…"
            self._set_analysis_progress("DE", 1, 4, "Preparing selection")

            if len(self.scatter.selection()) > 0:
                sel = np.asarray(self.scatter.selection(), dtype=int)
                label = "Current selection"
                self._set_analysis_progress("DE", 2, 4, "Computing statistics")
                res = self._compute_diffexpr(sel, label)
                self._pending_diffexpr = {"points": np.sort(np.unique(sel)), "results": res}
                self._set_analysis_progress("DE", 3, 4, "Formatting results")
                self._show_diffexpr_results(res, label, selected_indices=np.unique(sel))
            elif self.active_selection is not None:
                label = self.active_selection.name
                self._set_analysis_progress("DE", 2, 4, "Computing statistics")
                res = self._compute_diffexpr(self.active_selection.points, label)
                self.active_selection.cached_diffexpr = res
                self._set_analysis_progress("DE", 3, 4, "Formatting results")
                self._show_diffexpr_results(res, label, selected_indices=np.asarray(self.active_selection.points, dtype=int))
            self._set_analysis_progress("DE", 4, 4, "Done")
        finally:
            self.compute_diffexpr.disabled = False
            self.compute_diffexpr.description = "Compute DE"
            self._clear_analysis_progress("DE")

    def _color_by_change_handler(self, change):
        if self.scatter is None:
            return
        new = change["new"]
        if new in self.categorical_color_maps:
            self.scatter.color(by=new, map=self.categorical_color_maps[new])
        else:
            self.scatter.color(by=new, map="magma")

    def _clear_results_display(self, message: str | None = None):
        if self.selections_predicates is not None:
            if message:
                self.selections_predicates.children = (ipyw.HTML(message),)
            else:
                self.selections_predicates.children = ()

        if self.pathway_table_container is not None:
            self.pathway_table_container.children = ()
            self.pathway_table_container.layout.display = "none"

        if self.reactome_diagram_container is not None:
            self.reactome_diagram_container.children = ()
            self.reactome_diagram_container.layout.display = "none"

    def _de_source(self):
        adata = self.adata
        if getattr(adata, "raw", None) is not None and getattr(adata.raw, "X", None) is not None:
            return adata.raw.X, list(adata.raw.var_names), "raw"
        return adata.X, list(adata.var_names), "X"

    def _ensure_diffexpr_global_stats(self):
        X, var_names, source = self._de_source()
        if (
            self._diffexpr_global_stats is not None
            and self._diffexpr_global_stats.get("source") == source
            and self._diffexpr_global_stats.get("n_obs") == int(X.shape[0])
            and self._diffexpr_global_stats.get("n_vars") == int(X.shape[1])
        ):
            return

        if hasattr(X, "sum"):
            total_sum = np.asarray(X.sum(axis=0)).ravel()
        else:
            total_sum = np.asarray(np.sum(X, axis=0)).ravel()

        if hasattr(X, "power"):
            total_sqsum = np.asarray(X.power(2).sum(axis=0)).ravel()
        else:
            total_sqsum = np.asarray(np.einsum("ij,ij->j", X, X)).ravel()

        self._diffexpr_global_stats = {
            "source": source,
            "n_obs": int(X.shape[0]),
            "n_vars": int(X.shape[1]),
            "var_names": var_names,
            "sum": total_sum.astype(float, copy=False),
            "sqsum": total_sqsum.astype(float, copy=False),
        }           

    def _compute_diffexpr(self, selected_indices: np.ndarray, selection_label: str) -> list[dict]:
        self._ensure_diffexpr_global_stats()
        stats = self._diffexpr_global_stats
        if stats is None:
            return []

        X, var_names, source = self._de_source()
        total_sum = stats["sum"]
        total_sqsum = stats["sqsum"]

        selected_indices = np.asarray(selected_indices, dtype=int)
        selected_indices = selected_indices[(selected_indices >= 0) & (selected_indices < int(X.shape[0]))]
        selected_indices = np.unique(selected_indices)

        n1 = int(selected_indices.shape[0])
        n = int(X.shape[0])
        n2 = n - n1
        if n1 < 2 or n2 < 2:
            return []

        X1 = X[selected_indices, :]
        if hasattr(X1, "sum"):
            sum1 = np.asarray(X1.sum(axis=0)).ravel().astype(float, copy=False)
        else:
            sum1 = np.asarray(np.sum(X1, axis=0)).ravel().astype(float, copy=False)

        if hasattr(X1, "power"):
            sqsum1 = np.asarray(X1.power(2).sum(axis=0)).ravel().astype(float, copy=False)
        else:
            sqsum1 = np.asarray(np.einsum("ij,ij->j", X1, X1)).ravel().astype(float, copy=False)

        sum2 = (total_sum - sum1).astype(float, copy=False)
        sqsum2 = (total_sqsum - sqsum1).astype(float, copy=False)

        mean1 = sum1 / n1
        mean2 = sum2 / n2

        var1 = (sqsum1 - (sum1 * sum1) / n1) / (n1 - 1)
        var2 = (sqsum2 - (sum2 * sum2) / n2) / (n2 - 1)
        var1 = np.maximum(var1, 0.0)
        var2 = np.maximum(var2, 0.0)
        std1 = np.sqrt(var1)
        std2 = np.sqrt(var2)

        t_stat, p_val = ss.ttest_ind_from_stats(
            mean1=mean1,
            std1=std1,
            nobs1=n1,
            mean2=mean2,
            std2=std2,
            nobs2=n2,
            equal_var=False,
        )

        t_stat = np.asarray(t_stat, dtype=float)
        p_val = np.asarray(p_val, dtype=float)
        ok = np.isfinite(t_stat) & np.isfinite(p_val)

        t_thr = float(self.diff_t_threshold.value) if self.diff_t_threshold is not None else 2.0
        p_thr = float(self.diff_p_threshold.value) if self.diff_p_threshold is not None else 0.05

        keep = ok & (np.abs(t_stat) >= t_thr) & (p_val <= p_thr)
        idx = np.where(keep)[0]
        if idx.size == 0:
            return []

        order = np.argsort(np.abs(t_stat[idx]))[::-1]
        idx = idx[order]
        idx = idx[:200]

        out: list[dict] = []
        for j in idx:
            out.append(
                {
                    "attribute": var_names[int(j)],
                    "interval": (float(t_stat[int(j)]), float(p_val[int(j)])),
                    "quality": float(abs(t_stat[int(j)])),
                    "direction": selection_label,
                }
            )
        return out

    def _show_diffexpr_results(
        self,
        diff_results: list[dict],
        selection_label: str,
        selected_indices: np.ndarray | None = None,
    ):
        if self.selections_predicates is None:
            return
        if self.pathway_table_container is None or self.reactome_diagram_container is None:
            return

        if len(diff_results) == 0:
            self._clear_results_display(
                f"<em>No differential results passed thresholds for <b>{selection_label}</b>.</em>"
            )
            return

        rows = []
        for entry in diff_results:
            t, p = entry["interval"]
            rows.append(
                {
                    "Gene": entry["attribute"],
                    "T": float(np.round(float(t), 4)),
                    "p": f"{float(p):.3e}",
                    "Selection": selection_label,
                }
            )

        results_df = pd.DataFrame(rows)
        results_df = results_df.dropna(subset=["T", "p"])
        results_df["_absT"] = results_df["T"].abs()
        results_df = results_df.sort_values(by="_absT", ascending=False).drop(columns=["_absT"]).reset_index(drop=True)

        gene_table_widget = CorrelationTable(
            data=results_df.to_dict(orient="records"),
            columns=["Gene", "T", "p", "Selection"],
        )

        from ipywidgets import HTML as _HTML

        violin = GeneViolinPlot()
        title = _HTML("")
        plot_box = VBox(
            [title, violin],
            layout=Layout(
                flex="0 0 auto",
                height="270px",
                max_height="270px",
                min_height="270px",
                padding="0px",
                overflow="hidden",
                align_self="stretch",
                display="none",
            ),
        )
        self.pathway_table_container.children = [plot_box]
        self.pathway_table_container.layout.display = "block"
        self.reactome_diagram_container.layout.display = "none"

        def on_gene_click(change):
            gene = change["new"]
            if not gene:
                return
            try:
                X, var_names, _ = self._de_source()
                if gene not in var_names:
                    return
                gene_idx = int(var_names.index(gene))

                if selected_indices is not None:
                    sel = np.unique(np.asarray(selected_indices, dtype=int))
                elif self.active_selection is not None:
                    sel = np.unique(np.asarray(self.active_selection.points, dtype=int))
                elif self.scatter is not None:
                    sel = np.unique(np.asarray(self.scatter.selection(), dtype=int))
                else:
                    return
                n_obs = int(X.shape[0])
                if sel.size < 1 or sel.size >= n_obs:
                    return

                max_sel = 5000
                max_bg = 5000
                rng = np.random.default_rng(0)

                if sel.size > max_sel:
                    sel_samp = rng.choice(sel, size=max_sel, replace=False)
                else:
                    sel_samp = sel

                mask = np.ones(n_obs, dtype=bool)
                mask[sel] = False
                bg_idx = np.flatnonzero(mask)
                if bg_idx.size > max_bg:
                    bg_samp = rng.choice(bg_idx, size=max_bg, replace=False)
                else:
                    bg_samp = bg_idx

                col_sel = X[sel_samp, gene_idx]
                col_bg = X[bg_samp, gene_idx]
                sel_vals = np.asarray(col_sel.toarray() if hasattr(col_sel, "toarray") else col_sel).ravel().astype(float)
                bg_vals = np.asarray(col_bg.toarray() if hasattr(col_bg, "toarray") else col_bg).ravel().astype(float)

                vmin = float(np.min([np.min(sel_vals) if sel_vals.size else 0.0, np.min(bg_vals) if bg_vals.size else 0.0]))
                vmax = float(np.max([np.max(sel_vals) if sel_vals.size else 0.0, np.max(bg_vals) if bg_vals.size else 0.0]))
                if not np.isfinite(vmin) or not np.isfinite(vmax):
                    return
                if vmax <= vmin:
                    vmax = vmin + 1e-6

                bins = np.linspace(vmin, vmax, 51)
                sel_hist, _ = np.histogram(sel_vals, bins=bins, density=True)
                bg_hist, _ = np.histogram(bg_vals, bins=bins, density=True)

                violin.gene = gene
                violin.data = {
                    "bins": bins.tolist(),
                    "selected": sel_hist.tolist(),
                    "background": bg_hist.tolist(),
                }
                title.value = f"<b>{gene}</b> — selected n={sel.size}, background n={n_obs - sel.size} (sampled)"
                plot_box.layout.display = "block"
            except Exception:
                self.logger.exception("Error rendering differential plot")

        gene_table_widget.observe(on_gene_click, names=["selected_gene"])
        self.selections_predicates.children = [gene_table_widget]

    def get_genes(self, sel_name="Selection 1"):
        """Get the genes in a named selection."""
        selections = self.selections.selections
        for sel in selections:
            if sel.name == sel_name:
                if sel.cached_results is not None:
                    gene_names = [entry["attribute"] for entry in sel.cached_results]
                    gene_corr = [entry["interval"][0] for entry in sel.cached_results]
                    gene_pval = [entry["interval"][1] for entry in sel.cached_results]
                    results = pd.DataFrame({'gene': gene_names, 'correlation': gene_corr, 'p-value': gene_pval})
                    results.sort_values(by=['correlation'], ascending=False, inplace=True)
                    return results
        return pd.DataFrame()  # Empty DataFrame if selection not found

    def show(self):
        """Display the ScSketch widget."""
        return self.ui
