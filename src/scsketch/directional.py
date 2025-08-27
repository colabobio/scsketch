import numpy as np
import pandas as pd
import traitlets
from dataclasses import dataclass, field
from itertools import cycle
from IPython.display import display, HTML
from ipywidgets import Checkbox, Dropdown, GridBox, HBox, Layout, IntText, Text, VBox
from jscatter import Scatter, glasbey_light, link, okabe_ito, Line
from jscatter.widgets import Button
from matplotlib.colors import to_hex
from scipy.spatial import ConvexHull
import scipy.stats as ss
import requests

from .widgets import (
    GenePathwayWidget, CorrelationTable, PathwayTable,
    InteractiveSVG, Label, Div
)

#Widget Composition - Finally, we're going to instantiate the scatter plot and all the other widgets and link them using their traits. The output is the UI you've been waiting for :)
def view(adata, metadata_cols=None, max_gene_options=50):
    """
    Visualize an AnnData object in scSketch.

    Parameters
    ----------
    adata: AnnData
        The annotated data matrix (must contain a UMAP in adata.obsm["X_umap"]).
    metadata_cols: list of str, optional
        List of obs columns to include as metadata (e.g., ['dpi', 'strain', ...]).
        If None, metadata will be skipped.
    """
    
    print("I am in view function")
    from IPython.display import display, HTML
    import numpy as np
    
    from collections import OrderedDict
    from dataclasses import dataclass, field
    
    from IPython.display import HTML
    from ipywidgets import Checkbox, Dropdown, GridBox, HBox, Layout, IntText, Text, VBox
    from itertools import cycle
    from jscatter import Scatter, glasbey_light, link, okabe_ito, Line
    from jscatter.widgets import Button
    from numpy import histogram, isnan
    from matplotlib.colors import to_hex
    from scipy.spatial import ConvexHull
####################################################################################################################################
    import pandas as pd
    
    # UMAP coordinates
    umap_df = pd.DataFrame(
        adata.obsm["X_umap"], 
        columns=["x", "y"], 
        index=adata.obs_names,
    )

    # Metadata (optional)
    if metadata_cols is not None:
        available_metadata_cols = [col for col in metadata_cols if col in adata.obs.columns]
        if len(available_metadata_cols) > 0:
            metadata_df = adata.obs[available_metadata_cols].copy()
            # cast to str for categorical handling
            for col in available_metadata_cols:
                metadata_df[col] = metadata_df[col].astype(str)
            print(f"Using metadata columns: {available_metadata_cols}")
        else:
            print("No requested metadata columns found, continuing without metadata.")
            available_metadata_cols = []
            metadata_df = pd.DataFrame(index=adata.obs_names)
           
    else:
        available_metadata_cols = []
        metadata_df = pd.DataFrame(index=adata.obs_names)
        print("No metadata passed, continuing with UMAP + gene expression only.")

    # Gene Expression
    gene_exp_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        columns = adata.var_names,
        index = adata.obs_names,
    )

    # Combine
    df = pd.concat([umap_df, metadata_df, gene_exp_df], axis=1)
    df = df.loc[:, ~df.columns.duplicated()]

    # Colors for selections (used later)
    all_colors = okabe_ito.copy()
    available_colors = [color for color in all_colors]

    #Continuous color ramps for subdivided selections (used later)
    continuous_color_maps = [
    ["#00dadb", "#da00db"],
    ["#00dadb", "#a994dc", "#da00db"],
    ["#00dadb", "#8faddc", "#bd77dc", "#da00db"],
    ["#00dadb", "#7eb9dc", "#a994dc", "#c567dc", "#da00db"],
    ["#00dadb", "#72c0db", "#9aa3dc", "#b583dc", "#ca5cdb", "#da00db"],
    ["#00dadb", "#69c4db", "#8faddc", "#a994dc", "#bd77dc", "#cd54db", "#da00db"],
    [
        "#00dadb","#62c7db","#86b4dc","#9e9fdc","#b288dc","#c16edc","#cf4ddb","#da00db",
    ],
    [
        "#00dadb","#5ccadb","#7eb9dc","#96a7dc","#a994dc","#b87fdc","#c567dc","#d048db","#da00db",
    ],
    [
        "#00dadb","#57ccdb","#78bddc","#8faddc","#a19ddc","#b08bdc","#bd77dc","#c861db","#d144db","#da00db",
    ],
]

    
    # Categorical color maps
    categorical_cols = [col for col in available_metadata_cols if col in df.columns]
    categorical_color_maps = {
        col:dict(zip(df[col].unique(), cycle(glasbey_light)))
        for col in categorical_cols
    }

    # Pick default coloring 
    color_by_name = "seurat_clusters" if "seurat_clusters" in df.columns else None

    scatter = Scatter(
        data = df,
        x = "x",
        y = "y",
        background_color = "#111111",
        axes = False,
        height = 720, 
        color_by = color_by_name,
        color_map = categorical_color_maps.get(color_by_name, "plasma"),
        tooltip = True, 
        legend = True,
        tooltip_properties = list(categorical_cols), # only columsn that actually exist
    )

    # If no metadata-based color, use a bright solid color for all points
    if color_by_name is None:
        # Ensure we're not using a colormap
        scatter.color(by=None, map=None)
        # Set a high-contrast RGBA (0..1); e.g., bright gold
        scatter.widget.color = (0.98, 0.82, 0.20, 1.0)
        
    
    # display(scatter.show())
    
    @dataclass
    class Selection:
        """Class for keeping track of a selection."""
    
        index: int
        name: str
        points: np.ndarray
        color: str
        lasso: Line
        hull: Line
    
    
    @dataclass
    class Selections:
        """Class for keeping track of selections."""
    
        selections: list[Selection] = field(default_factory=list)
    
        def all_points(self) -> np.ndarray:
            return np.unique(
                np.concatenate(
                    list(map(lambda selection: selection.points, self.selections))
                )
            )
    
        def all_hulls(self) -> list[Line]:
            return [s.hull for s in self.selections]
    
    
    @dataclass
    class Lasso:
        """Class for keeping track of the lasso polygon."""
    
        polygon: Line | None = None
    
    
    lasso = Lasso()
    selections = Selections()
    
    
    def update_annotations():
        lasso_polygon = [] if lasso.polygon is None else [lasso.polygon]
        scatter.annotations(selections.all_hulls() + lasso_polygon)
    
    
    def lasso_selection_polygon_change_handler(change):
        if change["new"] is None:
            lasso.polygon = None
        else:
            points = change["new"].tolist()
            points.append(points[0])
            lasso.polygon = Line(points, line_color=scatter.widget.color_selected)
        update_annotations()
    
    
    scatter.widget.observe(
        lasso_selection_polygon_change_handler, names=["lasso_selection_polygon"]
    )
    
    selection_name = Text(value="", placeholder="Select some points…", disabled=True)
    selection_name.layout.width = "100%"
    
    selection_add = Button(
        description="",
        tooltip="Save Selection",
        disabled=True,
        icon="plus",
        width=36,
        rounded=["top-right", "bottom-right"],
    )
    
    selection_subdivide = Checkbox(value=False, description="Subdivide", indent=False)
    
    selection_num_subdivisions = IntText(
        value=5,
        min=2,
        max=10,
        step=1,
        description="Parts",
    )
    
    selection_subdivide_wrapper = HBox([selection_subdivide, selection_num_subdivisions])
    
    selections_elements = VBox(layout=Layout(grid_gap="2px"))
    
    selections_predicates_css = """
    <style>
    .jupyter-scatter-dimbridge-selections-predicates {
        position: absolute !important;
    }
    
    .jupyter-scatter-dimbridge-selections-predicates-wrapper {
        position: relative;
    }
    </style>
    """
    
    display(HTML(selections_predicates_css))
    
    selections_predicates = VBox(
        layout=Layout(
            top="4px",
            left="0px",
            right="0px",
            bottom="4px",
            grid_gap="4px",
        )
    )
    selections_predicates.add_class("jupyter-scatter-dimbridge-selections-predicates")
    
    selections_predicates_wrapper = VBox(
        [selections_predicates],
        layout=Layout(
            height="100%",
        ),
    )
    selections_predicates_wrapper.add_class(
        "jupyter-scatter-dimbridge-selections-predicates-wrapper"
    )
    
    compute_predicates = Button(
        description="Compute Directional Search",
        style="primary",
        disabled=True,
        full_width=True,
    )
    
    compute_predicates_between_selections = Checkbox(
        value=False, description="Compare Between Selections", indent=False
    )
    
    compute_predicates_wrapper = VBox([compute_predicates])
    
    
    def add_selection_element(selection: Selection):
        hex_color = to_hex(selection.color)
    
        selection_name = Label(
            name=selection.name,
            style={"background": hex_color},
        )
    
        selection_remove = Button(
            description="",
            tooltip="Remove Selection",
            icon="trash",
            width=36,
            background=hex_color,
            rounded=["top-right", "bottom-right"],
        )
    
        element = GridBox(
            [
                selection_name,
                selection_remove,
            ],
            layout=Layout(grid_template_columns="1fr 40px"),
        )
    
        def focus_handler(change):
            if change["new"]:
                scatter.zoom(to=selection.points, animation=500, padding=2)
            else:
                scatter.zoom(to=None, animation=500, padding=0)
    
        selection_name.observe(focus_handler, names=["focus"])
    
        def remove_handler(change):
            selections_elements.children = [
                e for e in selections_elements.children if e != element
            ]
            selections.selections = [s for s in selections.selections if s != selection]
            update_annotations()
            compute_predicates.disabled = len(selections.selections) == 0
    
        selection_remove.on_click(remove_handler)
    
        selections_elements.children = selections_elements.children + (element,)
    
    
    def add_subdivided_selections():
        lasso_polygon = scatter.widget.lasso_selection_polygon
        lasso_points = lasso_polygon.shape[0]
    
        lasso_mid = int(lasso_polygon.shape[0] / 2)
        lasso_spine = (lasso_polygon[:lasso_mid, :] + lasso_polygon[lasso_mid:, :]) / 2
    
        lasso_part_one = lasso_polygon[:lasso_mid, :]
        lasso_part_two = lasso_polygon[lasso_mid:, :][::-1]
    
        n_split_points = selection_num_subdivisions.value + 1
    
        sub_lassos_part_one = split_line_equidistant(lasso_part_one, n_split_points)
        sub_lassos_part_two = split_line_equidistant(lasso_part_two, n_split_points)
    
        base_name = selection_name.value
        if len(base_name) == 0:
            base_name = f"Selection {len(selections.selections) + 1}"
    
        color_map = continuous_color_maps[selection_num_subdivisions.value]
    
        for i, part_one in enumerate(sub_lassos_part_one):
            polygon = np.vstack((part_one, sub_lassos_part_two[i][::-1]))
            idxs = np.where(points_in_polygon(df[["x", "y"]].values, polygon))[0]
            points = df.iloc[idxs][["x", "y"]].values
            hull = ConvexHull(points)
            hull_points = np.vstack((points[hull.vertices], points[hull.vertices[0]]))
            color = color_map[i]
            name = f"{base_name}.{i + 1}"
    
            lasso_polygon = polygon.tolist()
            lasso_polygon.append(lasso_polygon[0])
    
            selection = Selection(
                index=len(selections.selections) + 1,
                name=name,
                points=idxs,
                color=color,
                lasso=Line(lasso_polygon),
                hull=Line(hull_points, line_color=color, line_width=2),
            )
            selections.selections.append(selection)
            add_selection_element(selection)
    
    
    def add_selection():
        idxs = scatter.selection()
        points = df.iloc[idxs][["x", "y"]].values
        hull = ConvexHull(points)
        hull_points = np.vstack((points[hull.vertices], points[hull.vertices[0]]))
        color = available_colors.pop(0)
    
        name = selection_name.value
        if len(name) == 0:
            name = f"Selection {len(selections.selections) + 1}"
    
        lasso_polygon = scatter.widget.lasso_selection_polygon.tolist()
        lasso_polygon.append(lasso_polygon[0])
    
        selection = Selection(
            index=len(selections.selections) + 1,
            name=name,
            points=idxs,
            color=color,
            lasso=Line(lasso_polygon),
            hull=Line(hull_points, line_color=color, line_width=2),
        )
        selections.selections.append(selection)
        add_selection_element(selection)
    
    
    def selection_add_handler(event):
        lasso.polygon = None
    
        if scatter.widget.lasso_type == "brush" and selection_subdivide.value:
            add_subdivided_selections()
        else:
            add_selection()
    
        compute_predicates.disabled = False
    
        scatter.selection([])
        update_annotations()
    
        if len(selections.selections) > 1:
            compute_predicates_wrapper.children = (
                compute_predicates_between_selections,
                compute_predicates,
            )
        else:
            compute_predicates_wrapper.children = (compute_predicates,)
    
    
    selection_add.on_click(selection_add_handler)
    
    
    def selection_handler(change):
        if len(change["new"]) > 0:
            selection_add.disabled = False
            selection_name.disabled = False
            selection_name.placeholder = "Name selection…"
            new_index = 1
            if len(selections.selections) > 0:
                new_index = selections.selections[-1].index + 1
            selection_name.value = f"Selection {new_index}"
        else:
            selection_add.disabled = True
            selection_name.disabled = True
            selection_name.placeholder = "Select some points…"
            selection_name.value = ""
    
    
    scatter.widget.observe(selection_handler, names=["selection"])
    
    
    def clear_predicates(event):
        compute_predicates.style = "primary"
        compute_predicates.description = "Compute Predicates"
        compute_predicates.on_click(compute_predicates_handler)
    
        selections_predicates.children = ()
    
        if len(selections.selections) > 1:
            compute_predicates_wrapper.children = (
                compute_predicates_between_selections,
                compute_predicates,
            )
        else:
            compute_predicates_wrapper.children = (compute_predicates,)
    
    
    import ipywidgets as widgets
    from IPython.display import display
    
    
    def fetch_pathways(gene):
        """Fetch Reactome pathways for a given gene symbol."""
        url = f"https://reactome.org/ContentService/data/mapping/UniProt/{gene}/pathways?species=9606"
        try:
            response = requests.get(url)
            response.raise_for_status()
            pathways = response.json()
            return [
                {"Pathway": entry["displayName"], "stId": entry["stId"]}
                for entry in pathways
            ]
        except requests.exceptions.RequestException as e:
            print(f"Error fetching Reactome pathways for {gene}: {e}")
            return []
    
    
    def fetch_pathway_image(pathway_id):
        """Fetch Reactome pathway diagram image."""
        url = f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.png"
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.content
        except requests.exceptions.RequestException as e:
            print(f"Error fetching pathway image: {e}")
            return None
    
    
    search_gene = widgets.Text()
    
    
    def show_directional_results(directional_results):
        # Display the results of the directional analysis as a table.
        # Args:directional_results (list): List of computed correlations from directional analysis.
    
        compute_predicates.style = ""
        compute_predicates.description = "Clear Results"
        compute_predicates.on_click(clear_predicates)  # Attach clear button
    
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
        # Filter out rows where 'R' or 'p' are NaN
        results_df = results_df.dropna(subset=["R", "p"])
        # Sort after removing NaNs
        results_df = results_df.sort_values(by="R", ascending=False).reset_index(drop=True)
    
        # create interactive table with click support
        # existing gene correlation table widget (already displayed):
        gene_table_widget = CorrelationTable(data=results_df.to_dict(orient="records"))
    
        # create new widgets explicitly for Reactome pathway table and diagram
        pathway_table_widget = PathwayTable(data=[])  # initially empty
    
        pathway_table_container.layout.display = "none"
        reactome_diagram_container.layout.display = "none"
    
        # link the selected gene in table to GenePathwayWidget
        # handlers for interactive selections
        def on_gene_click(change):
            gene = change["new"]
            print(f"[UI] gene clicked: {gene}")
            pathways = fetch_pathways(gene)
            pathway_table_widget.data = pathways
            pathway_table_container.layout.display = "block" if pathways else "none"
            reactome_diagram_container.layout.display = "none"
    
        import base64
        import requests
        from ipywidgets import HTML
    
        # instantiate the widget only once, outside the click handler
        interactive_svg_widget = InteractiveSVG()
    
        def on_pathway_click(change):
            pathway_id = change["new"]
            svg_url = (
                f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.svg"
            )
    
            try:
                response = requests.get(svg_url)
                response.raise_for_status()
                svg_content = response.text
                svg_base64 = base64.b64encode(svg_content.encode("utf-8")).decode("utf-8")
    
                interactive_svg_widget.svg_content = svg_base64
                reactome_diagram_container.layout.display = "block"
                reactome_diagram_container.children = [interactive_svg_widget]
    
            except requests.exceptions.RequestException as e:
                print(f"Error fetching SVG diagram: {e}")
    
        # connect handlers
        gene_table_widget.observe(on_gene_click, names=["selected_gene"])
        print("[UI] gene click observer attached")
        pathway_table_widget.observe(
            on_pathway_click, names=["selected_pathway"]
        )  # use selected_pathway traitlet
    
        # Show in the UI
        selections_predicates.children = [gene_table_widget]
    
        pathway_table_container.children = [
            widgets.HTML("<b>Reactome Pathways</b>"),
            pathway_table_widget,
        ]
    
        # reactome_diagram_container.children = [pathway_image_widget]
    
        print("Showing directional results...")
    
    
    #############################Part 1
    
    import numpy as np
    import scipy.stats as ss
    import pandas as pd
    
    
    def compute_directional_analysis(df, selections):
        # Computes the correlation of gene expression along a directional axis.
        # Args:
        # df (pd.DataFrame): Dataframe containing gene expression data and spatial coordinates.
        # selections (Selections): The selected points for directional analysis.
        # Returns:
        #     list: A list of dictionaries containing the computed correlations.
    
        if len(selections.selections) == 0:
            return []
    
        results = []
    
        for selection in selections.selections:
            selected_indices = selection.points
            selected_embeddings = df.iloc[selected_indices][["x", "y"]].values
    
            # Ensure we have at least two points for a valid direction vector
            if selected_embeddings.shape[0] < 2:
                continue
    
            # Compute direction vector
            v = selected_embeddings[-1] - selected_embeddings[0]
            v = v / np.linalg.norm(v)  # Normalize
    
            # Compute projections
            start_point = selected_embeddings[0]
            projections = np.array(
                [np.dot(pt - start_point, v) for pt in selected_embeddings]
            )

            base_drop = [
                "x", "y", 
                "dpi","strain","percent.cmv.log10","seurat_clusters","virus.presence",
                "Infection_localminima","UL123_define_infection","Infection_state","Infection_state_bkgd",
            ]
            extra_meta = available_metadata_cols if 'available_metadata_cols' in locals() else []
            # Get gene expression data
            columns_to_drop = [col for col in set(base_drop).union(extra_meta) if col in df.columns]
            
            selected_expression = df.iloc[selected_indices].drop(columns=columns_to_drop, errors="ignore")
    
            # Compute correlations
            correlations = []
            for gene in selected_expression.columns:
                r, p = ss.pearsonr(projections, selected_expression[gene])
                correlations.append(
                    {
                        "attribute": gene,
                        "interval": (r, p),
                        "quality": abs(
                            r
                        ),  # Use absolute correlation as a measure of quality
                    }
                )
    
            results.append(correlations)
    
        return results
    
    
    ######################Part 2
    
    
    def compute_predicates_handler(event):
        if len(selections.selections) == 0:
            return
    
        compute_predicates.disabled = True
        compute_predicates.description = "Computing Directional Analysis…"
    
        # Compute directional correlations
        directional_results = compute_directional_analysis(df, selections)
    
        # Display in a table instead of histogram
        show_directional_results(directional_results)
    
        compute_predicates.disabled = False
        from IPython.display import display
        import ipywidgets as widgets
    
        debug_output = widgets.Output()
        display(debug_output)
        with debug_output:
            print("Running directional analysis...")
    
    
    compute_predicates.on_click(compute_predicates_handler)
    
    
    compute_predicates.on_click(compute_predicates_handler)
    
    add = GridBox(
        [
            selection_name,
            selection_add,
        ],
        layout=Layout(grid_template_columns="1fr 40px"),
    )
    
    complete_add = VBox([add], layout=Layout(grid_gap="4px"))
    
    
    def lasso_type_change_handler(change):
        if change["new"] == "brush":
            complete_add.children = (add, selection_subdivide_wrapper)
        else:
            complete_add.children = (add,)
    
    
    scatter.widget.observe(lasso_type_change_handler, names=["lasso_type"])

    # Dynamic Color By menu from available metadata + genes (limit to keep UI responsive)
    metadata_cols_lower = [c.lower() for c in available_metadata_cols]
    gene_options = [(g, g) for g in list(adata.var_names)[:max_gene_options]]

    # Initial selection matches the initial scatter coloring
    initial_color_value = color_by_name if color_by_name is not None else (gene_options[0][1] if gene_options else None)
    from ipywidgets import Dropdown

    color_by = Dropdown(
        options=(
            [("Seurat Clusters", "seurat_clusters")] if "seurat_clusters" in df.columns else []
        )
        + [(c.capitalize(), c) for c in metadata_cols_lower if c not in ["x", "y", "seurat_clusters"]]
        + ([("— Genes —", None)] if gene_options else [])
        + gene_options,
        value=initial_color_value,
        description="Color By:",
    )

    def color_by_change_handler(change):
        selected_col = change["new"]
        if selected_col is None:
            return # ignore the divider
        if selected_col in categorical_color_maps:
            scatter.color(by = selected_col, map = categorical_color_maps[selected_col])
        else:
            scatter.color(by = selected_col, map = "plasma")
    color_by.observe(color_by_change_handler, names = ["value"])
 
    
    # Main scatterplot and color selection
    plot_wrapper = VBox([scatter.show(), color_by])
    
    pathway_table_container = VBox(
        [],
        layout=Layout(
            overflow_y="auto",
            height="400px",
            border="1px solid #ddd",
            padding="10px",
            display="none",
        ),
    )
    
    reactome_diagram_container = VBox(
        [], layout=Layout(overflow_y="auto", height="400px", padding="10px", display="none")
    )
    
    # Sidebar with selection controls
    sidebar = GridBox(
        [
            complete_add,
            selections_elements,
            selections_predicates_wrapper,
            compute_predicates_wrapper,
        ],
        layout=Layout(
            # grid_template_rows='min-content max-content 1fr min-content',
            grid_template_rows="min-content max-content 1fr min-content",
            overflow_y="auto",
            height="800px",
            grid_gap="4px",
            # height='100%',
        ),
    )
    
    # Pathway table (right panel)
    pathway_table_container.layout = Layout(
        overflow_y="auto",
        height="800px",
        border="1px solid #ddd",
        padding="10px",
        display="none",  # initially hidden until gene selection
    )
    
    # Pathway diagram (bottom panel)
    reactome_diagram_container.layout = Layout(
        overflow_y="auto",
        height="800px",
        border="1px solid #ddd",
        padding="10px",
        display="none",  # initially hidden until pathway selection
    )
    
    # Combine top three panels
    top_layout = GridBox(
        [
            plot_wrapper,
            sidebar,
            pathway_table_container,
        ],
        layout=Layout(
            grid_template_columns="2fr 1fr 1fr",
            grid_gap="10px",
            height="auto",
        ),
    )
    
    from IPython.display import display, HTML
    
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
    
    # Final combined layout
    combined_gene_pathway_panel = GridBox(
        [
            VBox([sidebar], layout=Layout(overflow_y="auto", height="800px")),
            VBox(
                [pathway_table_container], layout=Layout(overflow_y="auto", height="800px")
            ),
        ],
        layout=Layout(
            grid_template_columns="3fr 2fr",  # Gene table 60% and pathway table 40%
            grid_gap="5px",
            # overflow='hidden',
            height="800px",
        ),
    )
    
    # Update the top-level GridBox to include only two columns now
    top_layout_updated = GridBox(
        [
            plot_wrapper,
            combined_gene_pathway_panel,  # combined gene/pathway panel
        ],
        layout=Layout(
            grid_template_columns="3fr 2fr",  # Scatterplot 60%, combined panel 40%
            grid_gap="10px",
            # overflow='hidden',
            height="800px",
        ),
    )
    
    # Final updated layout with pathway diagram at the bottom
    final_layout_updated = VBox(
        [
            top_layout_updated,
            reactome_diagram_container,
        ],
        layout=Layout(grid_gap="10px", width="100%", height="auto"),
    )
    
    # Display the final layout
    display(final_layout_updated)