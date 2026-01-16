"""Main ScSketch widget for interactive single-cell embedding exploration.

The UI composition and selection management patterns in this module are adapted from
the dimbridge notebook in jupyter-scatter by Fritz Lekschas:
https://github.com/flekschas/jupyter-scatter/blob/main/notebooks/dimbridge.ipynb
"""

import base64
import pandas as pd
import requests
from itertools import cycle
from typing import Optional, List

from anndata import AnnData
from IPython.display import display, HTML
from ipywidgets import Checkbox, Dropdown, GridBox, HBox, Layout, IntText, Text, VBox
from jscatter import Scatter, glasbey_light, okabe_ito, Line
from jscatter.widgets import Button
from matplotlib.colors import to_hex

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
            max_genes: Maximum number of genes to include (0 for all)
            verbosity: Logging verbosity level
        """
        self.logger = configure_logging(verbosity)
        self.adata = adata
        self.metadata_cols = metadata_cols
        self.color_by_default = color_by_default
        self.height = height
        self.background_color = background_color
        self.max_genes = max_genes        
        self.verbosity = verbosity        

        # Initialize color management
        self.all_colors = okabe_ito.copy()
        self.available_colors = [color for color in self.all_colors]

        # Initialize state management
        self.lasso = Lasso()
        self.selections = Selections()

        self._build_df()
        self._build_scatter()

        # Build the UI
        self._build_ui()
        self._setup_handlers()

    def _log(self, *args):
        self.logger.debug(" ".join(str(a) for a in args))

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
                metadata_df = self.adata.obs[available_metadata_cols].copy()
                for col in available_metadata_cols:
                    if pd.api.types.is_object_dtype(metadata_df[col]) or pd.api.types.is_categorical_dtype(metadata_df[col]):
                        metadata_df[col] = metadata_df[col].astype(str)
            else:
                self.logger.info("No requested metadata columns found; continuing without metadata.")
                available_metadata_cols = []
                metadata_df = pd.DataFrame(index=self.adata.obs_names)
        else:
            available_metadata_cols = []
            metadata_df = pd.DataFrame(index=self.adata.obs_names)
            self.logger.info("No metadata passed; continuing with UMAP + gene expression only.")

        all_gene_sorted = sorted(list(self.adata.var_names), key=lambda s: s.lower())
        if self.max_genes > 0:
            gene_subset = all_gene_sorted[: self.max_genes]
        else: 
            gene_subset = all_gene_sorted
        if len(gene_subset) > 0:
            subX = self.adata[:, gene_subset].X
            if hasattr(subX, "toarray"):
                subX = subX.toarray()
            gene_exp_df = pd.DataFrame(subX, columns=gene_subset, index=self.adata.obs_names)
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
        """Build the complete UI layout."""
        # Selection controls
        self.selection_name = Text(
            value="", placeholder="Select some points…", disabled=True
        )
        self.selection_name.layout.width = "100%"

        self.selection_add = Button(
            description="",
            tooltip="Save Selection",
            disabled=True,
            icon="plus",
            width=36,
            rounded=["top-right", "bottom-right"],
        )

        self.selection_subdivide = Checkbox(
            value=False, description="Subdivide", indent=False
        )
        self.selection_num_subdivisions = IntText(
            value=5, min=2, max=10, step=1, description="Parts"
        )
        self.selection_subdivide_wrapper = HBox(
            [self.selection_subdivide, self.selection_num_subdivisions]
        )

        self.selections_elements = VBox(layout=Layout(grid_gap="2px"))

        # Directional analysis controls
        self.selections_predicates = VBox(
            layout=Layout(
                top="4px", left="0px", right="0px", bottom="4px", grid_gap="4px"
            )
        )
        self.selections_predicates.add_class(
            "jupyter-scatter-scsketch-selections-predicates"
        )

        self.selections_predicates_wrapper = VBox(
            [self.selections_predicates],
            layout=Layout(height="100%"),
        )
        self.selections_predicates_wrapper.add_class(
            "jupyter-scatter-scsketch-selections-predicates-wrapper"
        )

        self.compute_predicates = Button(
            description="Compute Directional Search",
            style="primary",
            disabled=True,
            full_width=True,
        )

        self.compute_predicates_between_selections = Checkbox(
            value=False, description="Compare Between Selections", indent=False
        )

        self.compute_predicates_wrapper = VBox([self.compute_predicates])

        # Color by dropdown
        _cat_priority = ["cell_population", "seurat_clusters", "leiden", "clusters"]
        cat_in_df = [c for c in _cat_priority if c in self.categorical_cols]
        cat_rest = sorted([c for c in self.categorical_cols if c not in _cat_priority], key=lambda s: s.lower())
        cat_ordered = cat_in_df + cat_rest
        cat_opts = [(c.replace("_", " ").title(), c) for c in cat_ordered]

        _qc_order = ["n_genes", "total_counts", "pct_counts_mt"]
        obs_opts = [(c.replace("_", " ").title(), c) for c in _qc_order if c in self.df.columns]

        if self.max_genes > 0:
            gene_options = [(g, g) for g in self.all_gene_sorted[: self.max_genes]]
        else: 
            gene_options = [(g, g) for g in self.all_gene_sorted]
        dropdown_options = cat_opts + obs_opts + gene_options

        self.color_by = Dropdown(
            options=dropdown_options,
            value=self.color_by_default,
            description="Color By:",
        )

        # Pathway exploration widgets
        self.pathway_table_container = VBox(
            [],
            layout=Layout(
                overflow_y="auto",
                height="800px",
                border="1px solid #ddd",
                padding="10px",
                display="none",
            ),
        )

        self.reactome_diagram_container = VBox(
            [],
            layout=Layout(
                overflow_y="auto", height="800px", padding="10px", display="none"
            ),
        )

        self.interactive_svg_widget = InteractiveSVG()

        # Add controls layout
        self.add_controls = GridBox(
            [self.selection_name, self.selection_add],
            layout=Layout(grid_template_columns="1fr 40px"),
        )

        self.complete_add = VBox([self.add_controls], layout=Layout(grid_gap="4px"))

        # Sidebar layout
        self.sidebar = GridBox(
            [
                self.complete_add,
                self.selections_elements,
                self.selections_predicates_wrapper,
                self.compute_predicates_wrapper,
            ],
            layout=Layout(
                grid_template_rows="min-content max-content 1fr min-content",
                overflow_y="auto",
                height="800px",
                grid_gap="4px",
            ),
        )

        # Combined gene/pathway panel
        self.combined_gene_pathway_panel = GridBox(
            [
                VBox([self.sidebar], layout=Layout(overflow_y="auto", height="800px")),
                VBox(
                    [self.pathway_table_container],
                    layout=Layout(overflow_y="auto", height="800px"),
                ),
            ],
            layout=Layout(
                grid_template_columns="3fr 2fr", grid_gap="5px", height="800px"
            ),
        )

        # Plot wrapper
        self.plot_wrapper = VBox([self.scatter.show(), self.color_by])

        # Top layout
        self.top_layout = GridBox(
            [self.plot_wrapper, self.combined_gene_pathway_panel],
            layout=Layout(
                grid_template_columns="3fr 2fr", grid_gap="10px", height="800px"
            ),
        )

        # Final layout
        self.layout = VBox(
            [self.top_layout, self.reactome_diagram_container],
            layout=Layout(grid_gap="10px", width="100%", height="auto"),
        )

        # Inject CSS
        display(
            HTML("""
        <style>
        .jupyter-scatter-scsketch-selections-predicates {
            position: absolute !important;
        }
        .jupyter-scatter-scsketch-selections-predicates-wrapper {
            position: relative;
        }
        .jp-OutputArea-output, .jp-Cell-outputArea, .jp-Notebook {
            overflow: auto !important;
            max-height: none !important;
        }
        </style>
        """)
        )

    def _setup_handlers(self):
        """Set up all event handlers."""
        # Lasso selection polygon change handler
        self.scatter.widget.observe(
            self._lasso_selection_polygon_change_handler,
            names=["lasso_selection_polygon"],
        )

        # Selection handler
        self.scatter.widget.observe(self._selection_handler, names=["selection"])

        # Lasso type change handler
        self.scatter.widget.observe(
            self._lasso_type_change_handler, names=["lasso_type"]
        )

        # Color by change handler
        self.color_by.observe(self._color_by_change_handler, names=["value"])

        # Add selection button handler
        self.selection_add.on_click(self._selection_add_handler)

        # Compute predicates button handler
        self.compute_predicates.on_click(self._compute_predicates_handler)

    def _update_annotations(self):
        """Update scatter plot annotations."""
        lasso_polygon = [] if self.lasso.polygon is None else [self.lasso.polygon]
        self.scatter.annotations(self.selections.all_hulls() + lasso_polygon)

    def _lasso_selection_polygon_change_handler(self, change):
        """Handle lasso polygon changes."""
        if change["new"] is None:
            self.lasso.polygon = None
        else:
            points = change["new"].tolist()
            points.append(points[0])
            self.lasso.polygon = Line(
                points, line_color=self.scatter.widget.color_selected
            )
        self._update_annotations()

    def _selection_handler(self, change):
        """Handle selection changes."""
        if len(change["new"]) > 0:
            self.selection_add.disabled = False
            self.selection_name.disabled = False
            self.selection_name.placeholder = "Name selection…"
            new_index = 1
            if len(self.selections.selections) > 0:
                new_index = self.selections.selections[-1].index + 1
            self.selection_name.value = f"Selection {new_index}"
        else:
            self.selection_add.disabled = True
            self.selection_name.disabled = True
            self.selection_name.placeholder = "Select some points…"
            self.selection_name.value = ""

    def _lasso_type_change_handler(self, change):
        """Handle lasso type changes."""
        if change["new"] == "brush":
            self.complete_add.children = (
                self.add_controls,
                self.selection_subdivide_wrapper,
            )
        else:
            self.complete_add.children = (self.add_controls,)

    def _color_by_change_handler(self, change):
        """Handle color by dropdown changes."""
        selected_col = change["new"]
        if selected_col in self.categorical_color_maps:
            self.scatter.color(
                by=selected_col, map=self.categorical_color_maps[selected_col]
            )
        else:
            self.scatter.color(by=selected_col, map="plasma")

    def _add_selection(self):
        """Add a single selection."""
        idxs = self.scatter.selection()
        points = self.df.iloc[idxs][["x", "y"]].values
        color = self.available_colors.pop(0) if self.available_colors else okabe_ito[0]

        name = self.selection_name.value
        if len(name) == 0:
            name = f"Selection {len(self.selections.selections) + 1}"

        lasso_polygon = self.scatter.widget.lasso_selection_polygon

        selection = create_selection(
            index=len(self.selections.selections) + 1,
            name=name,
            points_indices=idxs,
            points_coords=points,
            lasso_polygon=lasso_polygon,
            color=color,
        )
        self.selections.selections.append(selection)
        self._add_selection_element(selection)

    def _add_selection_element(self, selection: Selection):
        """Add a selection element to the UI."""
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

        element = GridBox(
            [selection_name, selection_remove],
            layout=Layout(grid_template_columns="1fr 40px"),
        )

        def focus_handler(change):
            if change["new"]:
                self.scatter.zoom(to=selection.points, animation=500, padding=2)
            else:
                self.scatter.zoom(to=None, animation=500, padding=0)

        selection_name.observe(focus_handler, names=["focus"])

        def remove_handler(change):
            self.selections_elements.children = [
                e for e in self.selections_elements.children if e != element
            ]
            self.selections.selections = [
                s for s in self.selections.selections if s != selection
            ]
            self._update_annotations()
            self.compute_predicates.disabled = len(self.selections.selections) == 0

        selection_remove.on_click(remove_handler)

        self.selections_elements.children = self.selections_elements.children + (
            element,
        )

    def _selection_add_handler(self, event):
        """Handle add selection button click."""
        self.lasso.polygon = None
        self._add_selection()
        self.compute_predicates.disabled = False
        self.scatter.selection([])
        self._update_annotations()

        if len(self.selections.selections) > 1:
            self.compute_predicates_wrapper.children = (
                self.compute_predicates_between_selections,
                self.compute_predicates,
            )
        else:
            self.compute_predicates_wrapper.children = (self.compute_predicates,)

    def _compute_predicates_handler(self, event):
        """Handle compute predicates button click."""
        if self.compute_predicates.description == "Clear Results":
            self._clear_predicates()
            return

        if len(self.selections.selections) == 0:
            return

        self.compute_predicates.disabled = True
        self.compute_predicates.description = "Computing Directional Analysis…"

        # Compute directional correlations
        directional_results = compute_directional_analysis(
            self.df,
            self.selections,
            metadata_columns=["x", "y"] + self.available_metadata_cols,
        )

        # Display results
        self._show_directional_results(directional_results)

        self.compute_predicates.disabled = False

    def _clear_predicates(self):
        """Clear directional analysis results."""
        self.compute_predicates.style = "primary"
        self.compute_predicates.description = "Compute Directional Search"
        self.selections_predicates.children = ()
        self.pathway_table_container.layout.display = "none"
        self.reactome_diagram_container.layout.display = "none"

        if len(self.selections.selections) > 1:
            self.compute_predicates_wrapper.children = (
                self.compute_predicates_between_selections,
                self.compute_predicates,
            )
        else:
            self.compute_predicates_wrapper.children = (self.compute_predicates,)

    def _show_directional_results(self, directional_results):
        """Display directional analysis results."""
        self.compute_predicates.style = ""
        self.compute_predicates.description = "Clear Results"

        all_results = []
        for i, result in enumerate(directional_results):
            for entry in result:
                all_results.append(
                    {
                        "Gene": entry["attribute"],
                        "R": round(entry["interval"][0], 4),
                        "p": round(entry["interval"][1], 6),
                    }
                )

        # Convert to DataFrame and sort by absolute R-value
        results_df = pd.DataFrame(all_results)
        results_df = results_df.dropna(subset=["R", "p"])
        results_df = results_df.sort_values(by="R", ascending=False).reset_index(
            drop=True
        )

        # Create interactive table
        gene_table_widget = CorrelationTable(data=results_df.to_dict(orient="records"))
        pathway_table_widget = PathwayTable(data=[])

        self.pathway_table_container.layout.display = "none"
        self.reactome_diagram_container.layout.display = "none"

        def on_gene_click(change):
            gene = change["new"]
            pathways = fetch_pathways(gene)
            pathway_table_widget.data = pathways
            self.pathway_table_container.layout.display = (
                "block" if pathways else "none"
            )
            self.reactome_diagram_container.layout.display = "none"

        def on_pathway_click(change):
            pathway_id = change["new"]
            svg_url = (
                f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.svg"
            )

            try:
                response = requests.get(svg_url)
                response.raise_for_status()
                svg_content = response.text
                svg_base64 = base64.b64encode(svg_content.encode("utf-8")).decode(
                    "utf-8"
                )

                self.interactive_svg_widget.svg_content = svg_base64
                self.reactome_diagram_container.layout.display = "block"
                self.reactome_diagram_container.children = [self.interactive_svg_widget]

            except requests.exceptions.RequestException as e:
                print(f"Error fetching SVG diagram: {e}")

        gene_table_widget.observe(on_gene_click, names=["selected_gene"])
        pathway_table_widget.observe(on_pathway_click, names=["selected_pathway"])

        self.selections_predicates.children = [gene_table_widget]
        self.pathway_table_container.children = [
            HTML("<b>Reactome Pathways</b>"),
            pathway_table_widget,
        ]

    def show(self):
        """Display the ScSketch widget."""
        return self.layout
