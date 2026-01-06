import numpy as np
import pandas as pd
import traitlets
from dataclasses import dataclass
from itertools import cycle
from IPython.display import display, HTML
from ipywidgets import Checkbox, Dropdown, GridBox, HBox, Layout, IntText, Text, VBox
from jscatter import Scatter, glasbey_light, link, okabe_ito, Line
from jscatter.widgets import Button
from matplotlib.colors import to_hex
from scipy.spatial import ConvexHull
import scipy.stats as ss
import requests

from .utils import (
    Lasso,
    Selection,
    Selections,
    find_equidistant_vertices,
    points_in_polygon,
    split_line_at_points,
    split_line_equidistant,
)
from .widgets import (
    GenePathwayWidget, CorrelationTable, PathwayTable,
    InteractiveSVG, Label, Div
)

from typing import Optional 
import ipywidgets as ipyw

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


#Widget Composition - Finally, we're going to instantiate the scatter plot and all the other widgets and link them using their traits. The output is the UI you've been waiting for :)
def _legacy_view(adata, metadata_cols=None, max_gene_options=50, fdr_alpha=0.05):
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
    
    # print("I am in view function")
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
                
            # print(f"Using metadata columns: {available_metadata_cols}")
        else:
            print("No requested metadata columns found, continuing without metadata.")
            available_metadata_cols = []
            metadata_df = pd.DataFrame(index=adata.obs_names)
           
    else:
        available_metadata_cols = []
        metadata_df = pd.DataFrame(index=adata.obs_names)
        print("No metadata passed, continuing with UMAP + gene expression only.")

    # # Gene Expression
    # gene_exp_df = pd.DataFrame(
    #     adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
    #     columns = adata.var_names,
    #     index = adata.obs_names,
    # )

    # # Combine
    # df = pd.concat([umap_df, metadata_df, gene_exp_df], axis=1)
    # df = df.loc[:, ~df.columns.duplicated()]

    # --- Gene Expression (thin slice for UI; avoid full densify) ---
    _gene_sorted = sorted(list(adata.var_names), key=lambda s: s.lower())
    gene_subset = _gene_sorted[:max_gene_options]
    
    if len(gene_subset) > 0:
        subX = adata[:, gene_subset].X
        if hasattr(subX, "toarray"):
            subX = subX.toarray()
        gene_exp_df = pd.DataFrame(subX, columns=gene_subset, index=adata.obs_names)
    else:
        gene_exp_df = pd.DataFrame(index=adata.obs_names)
    
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
    # priority_cats = [c for c in ["cell_population", "seurat_clusters", "leiden", "clusters"] if c in categorical_cols]
    priority_cats = [c for c in ["cell_type", "celltype.l1", "celltype.l2",
                             "cell_population", "seurat_clusters", "leiden", "clusters"]
                 if c in categorical_cols]
    if len(priority_cats) > 0:
        color_by = priority_cats[0]
        color_map = categorical_color_maps[color_by]
    else:
        fallback_cont = next((c for c in ["n_genes","total_counts","pct_counts_mt"] if c in df.columns), None)
        color_by = fallback_cont
        color_map = None

    scatter = Scatter(
        data=df,
        x="x", y="y",
        height=720,
        axes=False,
        background_color="#111111",
        color_by=color_by,
        color_map=color_map,
        tooltip=True,
        legend=False,               
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
            # points = change["new"].tolist()
            points = np.asarray(change["new"], dtype=float).tolist()
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

    selections_predicates = VBox(
    layout=Layout(overflow_y="auto", height="100%", grid_gap="6px")
    )
    
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

    def log(*args):
        with debug_out:
            print(*args)
    
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
       
        # lasso_polygon = scatter.widget.lasso_selection_polygon
        lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon, dtype = float)
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
           
            # lasso_polygon = polygon.tolist()
            lasso_polygon = polygon.astype(float).tolist()
            lasso_polygon.append(lasso_polygon[0])
            hull_points = hull_points.astype(float).tolist()

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
        hull_points_py = hull_points.astype(float).tolist()
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

        # lasso_polygon = scatter.widget.lasso_selection_polygon.tolist()
        lasso_polygon = np.asarray(scatter.widget.lasso_selection_polygon, dtype = float)
        lasso_polygon = lasso_polygon.tolist()
        lasso_polygon.append(lasso_polygon[0])
    
        selection = Selection(
            index=len(selections.selections) + 1,
            name=name,
            points=idxs,
            color=color,
            lasso=Line(lasso_polygon),
            hull=Line(hull_points_py, line_color=color, line_width=3),
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
    
        from .widgets import GeneProjectionPlot
        from ipywidgets import VBox, Layout, HTML as _HTML
    
        gene_proj_plot = GeneProjectionPlot()
        pathway_msg = _HTML("")  # default title - empty by default, nothing visible

        # --- put plot (top) + pathway table (bottom) into the right column ---
        from ipywidgets import VBox, Layout, HTML as _HTML
        
        plot_box = VBox(
            [gene_proj_plot],
            layout=Layout(
                flex="0 0 auto",
                height="300px",
                max_height="300px",
                min_height="300px",
                #border="1px solid #4aa3ff",
                padding="0px",
                overflow="visible",
                align_self="stretch",
                display="none",
            ),
        )
        
        pathway_box = VBox(
            [pathway_msg, pathway_table_widget],
            layout=Layout(
                flex="1 1 auto",
                # max_height="70%",
                overflow="auto",             # scroll if long
            ),
        )
        
        # Make the container a 50/50 vertical split
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
                # --- pathways UI ---
                pathways = fetch_pathways(gene)
                pathway_table_widget.data = pathways
                pathway_table_container.layout.display = "block"
                pathway_msg.value = (
                    f"<em>No Reactome pathways found for <b>{gene}</b>.</em>"
                    if len(pathways) == 0
                    else f"<b>Reactome pathways for {gene}</b>"
                )
                reactome_diagram_container.layout.display = "none"
        
                # --- context/selection guards ---
                ctx = get_context()
                if ctx is None:
                    log("No context found for projection plot")
                    return
                df = ctx.df
                selections = ctx.selections
        
                if len(selections.selections) == 0:
                    log("No selections available for gene projection plot.")
                    return
        
                sel = selections.selections[-1]
                selected_indices = sel.points
                log(f"[plot] selection name={sel.name!r} n={len(selected_indices)}")
        
                # --- build projection along selection (global normalized) ---
                X = df[["x", "y"]].to_numpy(dtype=float)
                Ux = (X[:, 0] - X[:, 0].min())
                Uy = (X[:, 1] - X[:, 1].min())
                
                px = float(np.ptp(Ux)) or 1.0
                py = float(np.ptp(Uy)) or 1.0
                
                U = np.c_[Ux/px, Uy/py]
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
        
                # --- expression vector (fallback to adata if gene not in df slice) ---
                if gene in df.columns:
                    expr_all = pd.to_numeric(df[gene], errors="coerce").to_numpy(dtype=float)
                    log(f"[plot] gene {gene!r} found in df columns")
                else:
                    sub = adata[:, [gene]].X
                    if hasattr(sub, "toarray"):
                        sub = sub.toarray()
                    expr_all = np.asarray(sub).ravel().astype(float)
                    log(f"[plot] gene {gene!r} NOT in df slice; loaded from adata (len={len(expr_all)})")
        
                expr = expr_all[selected_indices]
                # sanitize numerics
                if not np.isfinite(expr).all():
                    finite = np.isfinite(expr)
                    if finite.any():
                        med = float(np.median(expr[finite]))
                        expr = np.where(finite, expr, med)
                    else:
                        log(f"Expression for {gene} has no finite values in the selection.")
                        return
        
                # normalize safely
                dx = proj.max() - proj.min()
                dy = expr.max() - expr.min()
                proj = np.linspace(0, 1, len(proj)) if dx <= 1e-12 else (proj - proj.min()) / (dx + 1e-12)
                expr = np.zeros_like(expr)           if dy <= 1e-12 else (expr - expr.min()) / (dy + 1e-12)
        
                plot_data = [{"projection": float(p), "expression": float(e)} for p, e in zip(proj, expr)]
                gene_proj_plot.data = plot_data
                plot_box.layout.display = "block"  # show plot when gene is clicked
                gene_proj_plot.gene = gene
        
                log(f"[plot] {gene}: n={len(plot_data)} pts "
                    f"proj=[{proj.min():.3f},{proj.max():.3f}] expr=[{expr.min():.3f},{expr.max():.3f}]")
                log(f"[plot] first5={plot_data[:5]}")
        
            except Exception:
                # show full traceback in the debug panel so nothing is silent
                with debug_out:
                    traceback.print_exc()        

            # Normalize safely so max != min (prevents flat or NaN scaling in JS)
            dx = proj.max() - proj.min()
            dy = expr.max() - expr.min()
            
            if dx <= 1e-12:
                proj = np.linspace(0, 1, len(proj))
            else:
                proj = (proj - proj.min()) / (dx + 1e-12)
            
            if dy <= 1e-12:
                expr = np.zeros_like(expr)
            else:
                expr = (expr - expr.min()) / (dy + 1e-12)
            
            log(f"[plot] {gene}: n={proj.size}, proj=({proj.min():.4f},{proj.max():.4f}), expr=({expr.min():.4f},{expr.max():.4f})")
            
            # Push to widget
            plot_data = [{"projection": float(p), "expression": float(e)}
                         for p, e in zip(proj, expr)]
            gene_proj_plot.data = plot_data
            gene_proj_plot.gene = gene
        
            # Show under the table
            # selections_predicates.children = [gene_table_widget, gene_proj_plot]
        
            # (Optional) Debug line
            log(f"[plot] {gene}: n={len(plot_data)}")
            log(f"[plot-debug] sending {len(plot_data)} points to widget for {gene}")
            log(f"[plot-debug] first 5 points: {plot_data[:5]}")
            
        import base64
        import requests
        from ipywidgets import HTML
        from ipywidgets import Image

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
                with debug_out:
                    print("Error fetching SVG diagram:", e) 
                    
        # connect handlers
        gene_table_widget.observe(on_gene_click, names=["selected_gene"])
        print("[UI] gene click observer attached")
        pathway_table_widget.observe(
            on_pathway_click, names=["selected_pathway"]
        )  # use selected_pathway traitlet
    
        # Show in the UI
        selections_predicates.children = [gene_table_widget]
    
        log("Showing directional results...")
    
    
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
        [],
        layout=Layout(
            overflow="auto",
            min_height="0px",
            width="100%",
            max_width="100%",
            padding="10px",
            display="none",
            border="1px solid #ccc",
        ),
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
        max_height="400px",
        border="1px solid #ddd",
        padding="10px",
        display="none",  # initially hidden until gene selection
    )
    
    # Pathway diagram (bottom panel)
    # reactome_diagram_container.layout = Layout(
    #     overflow_y="auto",
    #     height="800px",
    #     border="1px solid #ddd",
    #     padding="10px",
    #     display="none",  # initially hidden until pathway selection
    # )
    
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
            grid_template_columns="2fr 3fr",  # Gene table 60% and pathway table 40%
            grid_gap="5px",
            # overflow='hidden',
            align_items="flex-start",
            height="auto",
            max_height="500px",
        ),
    )
    
    # Update the top-level GridBox to include only two columns now
    top_layout_updated = GridBox(
        [
            plot_wrapper,
            combined_gene_pathway_panel,  # combined gene/pathway panel
        ],
        layout=Layout(
            grid_template_columns="2.5fr 2.5fr",  # Scatterplot 60%, combined panel 40%
            grid_gap="10px",
            # overflow='hidden',
            height="auto",
            align_items="flex-start",
            # max_height="550px",
        ),
    )
    
    # Final updated layout with pathway diagram at the bottom
    final_layout_updated = VBox(
        [
            top_layout_updated,
            VBox([reactome_diagram_container], layout=Layout(width="100%")),
        ],
        layout=Layout(
            grid_gap="10px", 
            width="100%",
            # max_height="90vh", #cap to visible window height
            overflow_y="auto", #scroll inside, not the whole notebook
            align_items="flex-start",
        ),
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


class _ScSketchDirectionalView:
    """
    Object-oriented wrapper around the existing scSketch `view()` logic.

    This refactor is intended to keep widget behavior identical while organizing the code
    into: dataframe construction, UI construction, and handler wiring, with state stored on
    the instance.
    """

    def __init__(self, adata, metadata_cols=None, max_gene_options=50, fdr_alpha=0.05):
        self.adata = adata
        self.metadata_cols = metadata_cols
        self.max_gene_options = max_gene_options
        self.fdr_alpha = fdr_alpha

        self.df: pd.DataFrame | None = None
        self.available_metadata_cols: list[str] = []
        self.meta_cols_present: list[str] = []
        self.categorical_cols: list[str] = []
        self.categorical_color_maps: dict[str, dict] = {}
        self.color_by_default: str | None = None
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
        if self.debug_out is None:
            return
        with self.debug_out:
            print(*args)

    def _build_df(self):
        adata = self.adata

        umap_df = pd.DataFrame(
            adata.obsm["X_umap"],
            columns=["x", "y"],
            index=adata.obs_names,
        )

        if self.metadata_cols is not None:
            available_metadata_cols = [
                col for col in self.metadata_cols if col in adata.obs.columns
            ]
            if len(available_metadata_cols) > 0:
                metadata_df = adata.obs[available_metadata_cols].copy()
                for col in available_metadata_cols:
                    if pd.api.types.is_object_dtype(metadata_df[col]) or pd.api.types.is_categorical_dtype(metadata_df[col]):
                        metadata_df[col] = metadata_df[col].astype(str)
            else:
                print("No requested metadata columns found, continuing without metadata.")
                available_metadata_cols = []
                metadata_df = pd.DataFrame(index=adata.obs_names)
        else:
            available_metadata_cols = []
            metadata_df = pd.DataFrame(index=adata.obs_names)
            print("No metadata passed, continuing with UMAP + gene expression only.")

        _gene_sorted = sorted(list(adata.var_names), key=lambda s: s.lower())
        gene_subset = _gene_sorted[: self.max_gene_options]
        if len(gene_subset) > 0:
            subX = adata[:, gene_subset].X
            if hasattr(subX, "toarray"):
                subX = subX.toarray()
            gene_exp_df = pd.DataFrame(subX, columns=gene_subset, index=adata.obs_names)
        else:
            gene_exp_df = pd.DataFrame(index=adata.obs_names)

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
            color_by_default = priority_cats[0]
            color_map_default = categorical_color_maps[color_by_default]
        else:
            fallback_cont = next(
                (c for c in ["n_genes", "total_counts", "pct_counts_mt"] if c in df.columns),
                None,
            )
            color_by_default = fallback_cont
            color_map_default = None

        self.df = df
        self.available_metadata_cols = available_metadata_cols
        self.meta_cols_present = meta_cols_present
        self.categorical_cols = categorical_cols
        self.categorical_color_maps = categorical_color_maps
        self.color_by_default = color_by_default
        self.color_map_default = color_map_default

    def _build_scatter(self):
        df = self.df
        if df is None:
            raise RuntimeError("df was not built")

        scatter = Scatter(
            data=df,
            x="x",
            y="y",
            height=720,
            axes=False,
            background_color="#111111",
            color_by=self.color_by_default,
            color_map=self.color_map_default,
            tooltip=True,
            legend=False,
            tooltip_properties=[c for c in df.columns if c in self.meta_cols_present],
        )
        scatter.widget.color_selected = "#00dadb"
        self.scatter = scatter

    def _build_ui(self):
        self.debug_out = ipyw.Output(
            layout=ipyw.Layout(border="1px solid #444", max_height="140px", overflow_y="auto")
        )
        display(self.debug_out)

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
        self.compute_predicates_wrapper = VBox([self.compute_predicates])

        self.add_controls = GridBox(
            [self.selection_name, self.selection_add],
            layout=Layout(grid_template_columns="1fr 40px"),
        )
        self.complete_add = VBox([self.add_controls], layout=Layout(grid_gap="4px"))

    def _build_layout(self):
        if self.df is None or self.scatter is None:
            raise RuntimeError("View not initialized")

        df = self.df
        scatter = self.scatter

        _cat_priority = ["cell_population", "seurat_clusters", "leiden", "clusters"]
        cat_in_df = [c for c in _cat_priority if c in self.categorical_cols]
        cat_rest = sorted([c for c in self.categorical_cols if c not in _cat_priority], key=lambda s: s.lower())
        cat_ordered = cat_in_df + cat_rest
        cat_opts = [(c.replace("_", " ").title(), c) for c in cat_ordered]

        _qc_order = ["n_genes", "total_counts", "pct_counts_mt"]
        obs_opts = [(c.replace("_", " ").title(), c) for c in _qc_order if c in df.columns]

        _gene_sorted = sorted(list(self.adata.var_names), key=lambda s: s.lower())
        gene_options = [(g, g) for g in _gene_sorted[: self.max_gene_options]]
        dropdown_options = cat_opts + obs_opts + gene_options

        self.color_by = Dropdown(
            options=dropdown_options,
            value=self.color_by_default,
            description="Color By:",
        )
        self.plot_wrapper = VBox([scatter.show(), self.color_by])

        self.pathway_table_container = VBox(
            [],
            layout=Layout(
                overflow_y="auto",
                height="400px",
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
                overflow_y="auto",
                height="800px",
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
                VBox([self.sidebar], layout=Layout(overflow_y="auto", height="800px")),
                VBox([self.pathway_table_container], layout=Layout(overflow_y="auto", height="800px")),
            ],
            layout=Layout(
                grid_template_columns="2fr 3fr",
                grid_gap="5px",
                align_items="flex-start",
                height="auto",
                max_height="500px",
            ),
        )

        top_layout_updated = GridBox(
            [self.plot_wrapper, combined_gene_pathway_panel],
            layout=Layout(
                grid_template_columns="2.5fr 2.5fr",
                grid_gap="10px",
                height="auto",
                align_items="flex-start",
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
        self.color_by.observe(self._color_by_change_handler, names=["value"])

    def _update_annotations(self):
        scatter = self.scatter
        if scatter is None:
            return
        try:
            lasso_polygon = [] if self.lasso.polygon is None else [self.lasso.polygon]
            overlays = self.selections.all_hulls() + lasso_polygon
            scatter.annotations(overlays)
        except Exception:
            if self.debug_out is not None:
                with self.debug_out:
                    import traceback

                    traceback.print_exc()

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
            else:
                scatter.zoom(to=None, animation=500, padding=0)

        selection_name.observe(focus_handler, names=["focus"])

        def remove_handler(change):
            self.selections_elements.children = [e for e in self.selections_elements.children if e != element]
            self.selections.selections = [s for s in self.selections.selections if s != selection]
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

        scatter = self.scatter
        selection_subdivide = self.selection_subdivide

        try:
            self.lasso.polygon = None

            if scatter.widget.lasso_type == "brush" and selection_subdivide.value:
                self._add_subdivided_selections()
            else:
                self._add_selection()

            self.compute_predicates.disabled = False
            scatter.selection([])
            self._update_annotations()

            if len(self.selections.selections) > 1:
                self.compute_predicates_wrapper.children = (
                    self.compute_predicates_between_selections,
                    self.compute_predicates,
                )
            else:
                self.compute_predicates_wrapper.children = (self.compute_predicates,)
        except Exception:
            import traceback

            traceback.print_exc()

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
        else:
            self.selection_add.disabled = True
            self.selection_name.disabled = True
            self.selection_name.placeholder = "Select some points…"
            self.selection_name.value = ""

    def _clear_predicates(self, event):
        if self.compute_predicates is None or self.compute_predicates_wrapper is None:
            return

        self.compute_predicates.style = "primary"
        self.compute_predicates.description = "Compute Predicates"
        self.compute_predicates.on_click(self._compute_predicates_handler)

        if self.selections_predicates is not None:
            self.selections_predicates.children = ()

        if self.compute_predicates_between_selections is not None and len(self.selections.selections) > 1:
            self.compute_predicates_wrapper.children = (
                self.compute_predicates_between_selections,
                self.compute_predicates,
            )
        else:
            self.compute_predicates_wrapper.children = (self.compute_predicates,)

    def _fetch_pathways(self, gene):
        url = f"https://reactome.org/ContentService/data/mapping/UniProt/{gene}/pathways?species=9606"
        try:
            response = requests.get(url)
            response.raise_for_status()
            pathways = response.json()
            return [{"Pathway": entry["displayName"], "stId": entry["stId"]} for entry in pathways]
        except requests.exceptions.RequestException as e:
            print(f"Error fetching Reactome pathways for {gene}: {e}")
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
        debug_out = self.debug_out
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

                ctx = get_context()
                if ctx is None:
                    log("No context found for projection plot")
                    return
                df = ctx.df
                selections = ctx.selections

                if len(selections.selections) == 0:
                    log("No selections available for gene projection plot.")
                    return

                sel = selections.selections[-1]
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
                    expr_all = pd.to_numeric(df[gene], errors="coerce").to_numpy(dtype=float)
                    log(f"[plot] gene {gene!r} found in df columns")
                else:
                    sub = adata[:, [gene]].X
                    if hasattr(sub, "toarray"):
                        sub = sub.toarray()
                    expr_all = np.asarray(sub).ravel().astype(float)
                    log(f"[plot] gene {gene!r} NOT in df slice; loaded from adata (len={len(expr_all)})")

                expr = expr_all[selected_indices]
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
                if debug_out is not None:
                    with debug_out:
                        traceback.print_exc()

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
                if debug_out is not None:
                    with debug_out:
                        print("Error fetching SVG diagram:", e)

        gene_table_widget.observe(on_gene_click, names=["selected_gene"])
        print("[UI] gene click observer attached")
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
                continue

            v = selected_embeddings[-1] - selected_embeddings[0]
            v = v / np.linalg.norm(v)
            start_point = selected_embeddings[0]
            projections = np.array([np.dot(pt - start_point, v) for pt in selected_embeddings])

            base_drop = [
                "x",
                "y",
                "dpi",
                "strain",
                "percent.cmv.log10",
                "seurat_clusters",
                "virus.presence",
                "Infection_localminima",
                "UL123_define_infection",
                "Infection_state",
                "Infection_state_bkgd",
            ]
            columns_to_drop = [col for col in set(base_drop).union(self.available_metadata_cols) if col in df.columns]
            selected_expression = df.iloc[selected_indices].drop(columns=columns_to_drop, errors="ignore")

            batch_result_new = test_direction(selected_expression.values, projections)
            rs = batch_result_new["correlation"].astype(float)
            ps = batch_result_new["p_value"].astype(float)
            genes = list(selected_expression.columns)
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
        try:
            if len(self.selections.selections) == 0:
                return

            self.compute_predicates.disabled = True
            self.compute_predicates.description = "Computing Directional Analysis…"

            if self.compute_predicates_between_selections is not None and self.compute_predicates_between_selections.value:
                sels_for_run = self.selections
            else:
                last_only = Selections(selections=[self.selections.selections[-1]])
                sels_for_run = last_only

            directional_results = self._compute_directional_analysis(self.df, sels_for_run)
            self._show_directional_results(directional_results)
        except Exception:
            import traceback

            traceback.print_exc()
        finally:
            self.compute_predicates.disabled = False

    def _lasso_type_change_handler(self, change):
        if self.complete_add is None or self.add_controls is None or self.selection_subdivide_wrapper is None:
            return
        if change["new"] == "brush":
            self.complete_add.children = (self.add_controls, self.selection_subdivide_wrapper)
        else:
            self.complete_add.children = (self.add_controls,)

    def _color_by_change_handler(self, change):
        if self.scatter is None:
            return
        new = change["new"]
        if new in self.categorical_color_maps:
            self.scatter.color(by=new, map=self.categorical_color_maps[new])
        else:
            self.scatter.color(by=new, map="magma")

    def render(self) -> Context:
        if self.ui is None:
            raise RuntimeError("UI not built")
        display(self.ui)

        global _last_context
        ctx = Context(
            df=self.df,
            selections=self.selections,
            scatter=self.scatter,
            ui=self.ui,
            pathway_table_container=self.pathway_table_container,
            reactome_diagram_container=self.reactome_diagram_container,
        )
        _last_context = ctx
        return ctx


def view(adata, metadata_cols=None, max_gene_options=50, fdr_alpha=0.05):
    """Public API: construct and display the scSketch UI, returning the view Context."""
    return _ScSketchDirectionalView(
        adata,
        metadata_cols=metadata_cols,
        max_gene_options=max_gene_options,
        fdr_alpha=fdr_alpha,
    ).render()
