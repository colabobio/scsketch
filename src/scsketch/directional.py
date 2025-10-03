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

from typing import Optional 
import ipywidgets as ipyw

from matplotlib.path import Path

def find_equidistant_vertices(vertices: np.ndarray, n_points: int) -> np.ndarray:
    seg = np.diff(vertices, axis=0)
    seg_len = np.linalg.norm(seg, axis=1)
    cum = np.concatenate(([0.0], np.cumsum(seg_len)))
    total = cum[-1]
    targets = np.linspace(0.0, total, n_points)
    out = np.zeros((n_points, 2))
    for i, t in enumerate(targets):
        j = max(0, min(np.searchsorted(cum, t, side="right")-1, len(seg)-1))
        alpha = (t - cum[j]) / (seg_len[j] if seg_len[j] else 1.0)
        out[i] = vertices[j] + alpha * seg[j]
    return out

def split_line_at_points(vertices: np.ndarray, split_points: np.ndarray) -> list[np.ndarray]:
    if len(split_points) < 2:
        return []
    return [np.vstack([split_points[i], split_points[i+1]]) for i in range(len(split_points)-1)]

def split_line_equidistant(vertices: np.ndarray, n_points: int) -> list[np.ndarray]:
    return split_line_at_points(vertices, find_equidistant_vertices(vertices, n_points))

def points_in_polygon(points: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    return Path(polygon).contains_points(points)

# Inside view(), where Selection is defined
@dataclass
class Context:
    df: pd.DataFrame
    selections: object
    scatter: object
    ui: ipyw.Widget
    pathway_table_container: ipyw.Widget
    reactome_diagram_container: ipyw.Widget

_last_context: Optional[Context] = None

def get_context() -> Optional[Context]:
    """Return the most recent scSketch view context (or None if not created yet)."""
    return _last_context
###############################################################################################################################
def test_direction(X, projection):
    rs = np.corrcoef(projection, X, rowvar=False)[0, 1:]

    n = len(X)
    T = -np.abs(rs * np.sqrt(n - 2)) / np.sqrt(1 - (rs**2))
    return {"correlation": rs, "p_value": ss.t.cdf(T, df=n - 2) * 2}


def lord_test(pval, initial_results=None, gammai=None, alpha=0.05, w0=0.005):
    """"
    This is a translation of "version 1" under:

    https://github.com/bioc/onlineFDR/blob/devel/src/lord.cpp
    
    The only changes are that we don't recompute threhsolds for hypotheses that
    we have already seen. This only necessary because we may continue testing
    for many directions.
    """
    N = len(pval)

    if gammai is None:
        gammai = (
            0.07720838
            * np.log(np.maximum(np.arange(1, N + 2), 2))
            / (np.arange(1, N + 2) * np.exp(np.sqrt(np.log(np.arange(1, N + 2)))))
        )

    # setup variables, substituting previous results if needed
    alphai = np.zeros(N)
    R = np.zeros(N, dtype=bool)
    tau = []
    if initial_results is not None:
        N0 = len(initial_results["p_value"])
        alphai[range(N0)] = initial_results["alpha_i"]
        R[range(N0)] = initial_results["R"]
        tau = initial_results["tau"]
    else:
        N0 = 1
        alphai[0] = gammai[0] * w0
        R[0] = pval[0] <= alphai[0]
        if R[0]:
            tau.append(0)

    # compute lord thresholds iteratively
    K = int(np.sum(R))
    for i in range(N0, N):
        if K <= 1:
            if R[i - 1]:
                tau = [i - 1]
            Cjsum = sum(gammai[i - tau[j] - 1] for j in range(K))
            alphai[i] = w0 * gammai[i] + (alpha - w0) * Cjsum
        else:
            if R[i - 1]:
                tau.append(i - 1)
            tau2 = tau[1:]
            Cjsum = sum(gammai[i - tau2[j] - 1] for j in range(K - 1))
            alphai[i] = (
                w0 * gammai[i] + (alpha - w0) * gammai[i - tau[0] - 1] + alpha * Cjsum
            )

        if pval[i] <= alphai[i]:
            R[i] = True
            K += 1

    return {"p_value": pval, "alpha_i": alphai, "R": R, "tau": tau}

###############################################################################################################################

#Widget Composition - Finally, we're going to instantiate the scatter plot and all the other widgets and link them using their traits. The output is the UI you've been waiting for :)
def view(adata, metadata_cols=None, max_gene_options=50, fdr_alpha=0.05):
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
                if pd.api.types.is_object_dtype(metadata_df[col]) or pd.api.types.is_categorical_dtype(metadata_df[col]):
                    metadata_df[col] = metadata_df[col].astype(str)
                #else leave numerics as numerics
                
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

    # Determine categorical vs continuous from the metadata we actually have
    meta_cols_present = [c for c in (metadata_cols or []) if c in df.columns]

    categorical_cols = [
        c for c in meta_cols_present
        if df[c].dtype == "object" or df[c].nunique(dropna=False) <= 30
    ]
    # ensure string labels for categories
    for c in categorical_cols:
        df[c] = df[c].astype(str)

    def _sorted_cats(x: pd.Series):
        vals = pd.Index(x.unique().astype(str))
        # numeric sort if all labels are integers; else lexicographic
        if vals.str.fullmatch(r"\d+").all():
            return sorted(vals, key=lambda s: int(s))
        return sorted(vals, key=str)
    
    categorical_color_maps = {
        c: dict(zip(_sorted_cats(df[c]), cycle(glasbey_light[1:])))
        for c in categorical_cols
    }

    # Prefer a categorical default. Otherwise fall back to a continuous obs metric.
    priority_cats = [c for c in ["cell_population", "seurat_clusters", "leiden", "clusters"] if c in categorical_cols]
    if len(priority_cats) > 0:
        color_by = priority_cats[0]
        color_map = categorical_color_maps[color_by]
    else:
        fallback_cont = next((c for c in ["n_genes","total_counts","pct_counts_mt"] if c in df.columns), None)
        color_by = fallback_cont
        color_map = None
    # --- 3) Now instantiate Scatter with this default ---
    scatter = Scatter(
        data=df,
        x="x",
        y="y",
        background_color="#111111",
        axes=False,
        height=720,
        color_by=color_by,
        color_map=color_map,
        tooltip=True,
        legend=True,
        tooltip_properties=[c for c in df.columns if c in meta_cols_present],
    )

    # Colors for selections 
    all_colors = okabe_ito.copy()
    available_colors = [color for color in all_colors]

    #Continuous color ramps for subdivided selections 
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

    scatter.widget.color_selected = "#00dadb"
    
    @dataclass
    class Selection:
        """Class for keeping track of a selection."""
    
        index: int
        name: str
        points: np.ndarray
        color: str
        lasso: Line
        hull: Line
        path: np.ndarray | None = None 
    
    
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
        try:
            lasso_polygon = [] if lasso.polygon is None else [lasso.polygon]
            overlays = selections.all_hulls() + lasso_polygon
            scatter.annotations(overlays)
        except Exception as e:
            with debug_out:
                import traceback; traceback.print_exc()
    
    def lasso_selection_polygon_change_handler(change):
        if change["new"] is None:
            lasso.polygon = None
        else:
            points = change["new"].tolist()
            points.append(points[0]) #closes loop
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

     # --- DEBUG PANEL: visible print area under the plot ---
    import ipywidgets as widgets
    debug_out = widgets.Output(layout=widgets.Layout(border='1px solid #444', max_height='140px', overflow_y='auto'))
    display(debug_out)
    
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
                # path=lasso_spine,
            )
            selections.selections.append(selection)
            add_selection_element(selection)
    
    def add_selection():
        idxs = scatter.selection()
        points = df.iloc[idxs][["x", "y"]].values
        hull = ConvexHull(points)
        hull_points = np.vstack((points[hull.vertices], points[hull.vertices[0]]))
        color = available_colors.pop(0)
        
        # Build brush spine (midline of polygon) 
        spine = None
        if scatter.widget.lasso_type == "brush":
            lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon)
            if lasso_polygon.shape[0] >= 2:
                if lasso_polygon.shape[0] % 2 == 1:
                    lasso_polygon = lasso_polygon[:-1]
                mid = lasso_polygon.shape[0] // 2
                spine = (lasso_polygon[:mid, :] + lasso_polygon[mid:, :]) / 2
    
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
            hull=Line(hull_points, line_color=color, line_width=3),
            # path=lasso_spine,
        )
        selections.selections.append(selection)
        add_selection_element(selection)
    
    
    def selection_add_handler(event):
        try:        
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
        except Exception as e:
            import traceback; traceback.print_exc()
    
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
                #If Online FDR fields are present, hide non-significant genes
                if "reject" in entry and not entry["reject"]:
                    continue
                
                all_results.append(
                    {
                        "Gene": entry["attribute"],
                        "R": float(np.round(entry["interval"][0], 4)),
                        "p": f"{entry['interval'][1]:.3e}",
                        # Extra fields are appended to the row (widget will ignore extra columns)
                        # "alpha_i": float(entry.get("alpha_i", float("nan"))),
                        # "reject": bool(entry.get("reject", False)),ß
                        "Selection": entry.get("direction", f"Selection {i+1}"),
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
    
    batch_results = None
    online_results = None
    
    def compute_directional_analysis(df, selections):
        # Computes the correlation of gene expression along a directional axis.
        # Args:
        # df (pd.DataFrame): Dataframe containing gene expression data and spatial coordinates.
        # selections (Selections): The selected points for directional analysis.
        # Returns:
        #     list: A list of dictionaries containing the computed correlations.
        #Computes correlation along the selection's direction and appends Online FDR fields.
        #Only genes passing Online FDR are returned.

        nonlocal batch_results, online_results #keep LORD++ state across button clicks
        
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

            # Vectorized per-gene correlation and p-values for THIS selection only
            batch_result_new = test_direction(selected_expression.values, projections) # {'correlation': r, 'p_value':p}
            rs = batch_result_new["correlation"].astype(float)
            ps = batch_result_new["p_value"].astype(float)
            genes = list(selected_expression.columns)
            n_new = len(ps)

            #Append to the running Online FDR stream and compute new thresholds for the new tail
            prev_len = 0 if (batch_results is None) else len(batch_results["p_value"])
            p_values = ps if (batch_results is None) else np.concatenate([batch_results["p_value"], ps])
            
            online_results_new = lord_test(p_values, online_results, alpha=fdr_alpha)
            online_results = online_results_new
            batch_results = {"p_value": p_values} #keep the accumulated stream in the same variable name

            #Extract the chunk belonging to THIS selection
            alpha_chunk = online_results_new["alpha_i"][prev_len:prev_len + n_new]
            R_chunk = online_results_new["R"][prev_len:prev_len + n_new]

            #Build output rows: append alpha_i/reject and filter to keep only significant genes
            
            # Compute correlations
            correlations = []
            
            for j, gene in enumerate(genes):
                if not bool(R_chunk[j]):
                    continue # hide non-significant genes

                r = float(rs[j])
                p = float(ps[j])
                a = float(alpha_chunk[j])
                
                correlations.append(
                    {
                        "attribute": gene,
                        "interval": (r, p),     #(correlation, p-value)
                        "quality": abs(r),
                        "alpha_i": a,          #Online FDR threshold used for this gene
                        "reject": True,        #passed Online FDR 
                        "direction": selection.name, #keep which selection this came from
                    }
                )
            
            results.append(correlations)
    
        return results
    
    
    ######################Part 2
    
    
    def compute_predicates_handler(event):
        try:
            
            if len(selections.selections) == 0:
                return
        
            compute_predicates.disabled = True
            compute_predicates.description = "Computing Directional Analysis…"
    
            if compute_predicates_between_selections.value:
                #compare mode: use all saved selections 
                sels_for_run = selections
            else:
                #default: only the most recent selection
                #build a temporary Selections with the last one
                last_only = Selections(selections=[selections.selections[-1]])
                sels_for_run = last_only
    
            #Compute directional correlations for the chosen selection(s)
            directional_results = compute_directional_analysis(df, sels_for_run)
    
            #Show only what we just computed (i.e., last selection if not comparing)
            show_directional_results(directional_results)
            
        except Exception:
            import traceback; traceback.print_exc()
        finally:
            compute_predicates.disabled = False
            #compute_predicates.description = "Compute Directional Search" #optional restore
    
    
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

    # --- Build dropdown options (categoricals + QC first, then genes) ---

    # 1) Prioritize categorical labels in a fixed order, then any other categoricals (alpha)
    _cat_priority = ["cell_population", "seurat_clusters", "leiden", "clusters"]
    cat_in_df = [c for c in _cat_priority if c in categorical_cols]
    
    # any remaining categoricals (not in priority list), sorted by label
    cat_rest = sorted([c for c in categorical_cols if c not in _cat_priority],
                      key=lambda s: s.lower())
    
    cat_ordered = cat_in_df + cat_rest
    cat_opts = [(c.replace("_", " ").title(), c) for c in cat_ordered]
    
    # 2) QC metrics in a sensible fixed order
    _qc_order = [ "n_genes", "total_counts", "pct_counts_mt" ]
    obs_opts = [(c.replace("_", " ").title(), c) for c in _qc_order if c in df.columns]
    
    # 3) Genes: case-insensitive alphabetical, limited by max_gene_options
    _gene_sorted = sorted(list(adata.var_names), key=lambda s: s.lower())
    gene_options = [(g, g) for g in _gene_sorted[:max_gene_options]]
    
    # 4) Final: categoricals + QC first, then a separator, then genes
    dropdown_options = (
        cat_opts
        + obs_opts
        # + ([("— Genes —", None)] if gene_options else [])
        + gene_options
    )
    
    from ipywidgets import Dropdown
    color_by = Dropdown(
        options=dropdown_options,
        value=color_by,   # the default chosen earlier via priority_cats/fallback
        description="Color By:",
    )
    
    def color_by_change_handler(change):
        new = change["new"]
        if new in categorical_color_maps:
            scatter.color(by=new, map=categorical_color_maps[new])
        else:
            scatter.color(by=new, map="magma")
        #categorical (clusters) -> bright glasbey map built at init
        #continuous (genes/other numeric) -> magma 
        # cmap = color_map if (color_map is not None and new == "seurat_clusters") else "magma"
        # scatter.color(by=new, map=cmap)
        
    # Switch palettes: categorical → Glasbey, otherwise → magma    
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

    global _last_context
    ctx = Context(
        df=df,
        selections=selections,
        scatter=scatter,
        ui=final_layout_updated,
        pathway_table_container=pathway_table_container,
        reactome_diagram_container=reactome_diagram_container,
    )
    _last_context = ctx
    return ctx


