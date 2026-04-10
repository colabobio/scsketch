"""UI widget construction for ScSketch.

``build_controls()`` creates every ipywidget control used in the sidebar and
returns them bundled in a ``UIControls`` dataclass so callers can reference
them by name.  ``build_layout()`` assembles those controls plus the scatter
plot into the final ``VBox`` tree that ``ScSketch.show()`` returns.

No business logic lives here — only widget instantiation and layout wiring.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import ipywidgets as ipyw
from ipywidgets import (
    Checkbox,
    Dropdown,
    FloatText,
    GridBox,
    HBox,
    IntText,
    Layout,
    Text,
    VBox,
)
from IPython.display import display, HTML
from jscatter import Scatter
from jscatter.widgets import Button


# ── Progress spinner HTML (shared) ──────────────────────────────────────────

PROGRESS_SPINNER_HTML = (
    '<svg width="14" height="14" viewBox="0 0 50 50" '
    'style="vertical-align:-0.125em;margin-right:6px;">'
    '<circle cx="25" cy="25" r="20" fill="none" stroke="#999" stroke-width="5" '
    'stroke-linecap="round" stroke-dasharray="60 70">'
    '<animateTransform attributeName="transform" type="rotate" from="0 25 25" to="360 25 25" '
    'dur="0.8s" repeatCount="indefinite"/>'
    "</circle></svg>"
)


# ── Dataclass holding every widget reference ─────────────────────────────────

@dataclass
class UIControls:
    """All ipywidget instances created by :func:`build_controls`."""

    # ── Selection add row ────────────────────────────────────────────────────
    selection_name: Text
    selection_add: Button
    selection_subdivide: Checkbox
    selection_num_subdivisions: IntText
    selection_subdivide_wrapper: HBox
    selections_elements: VBox
    selections_predicates: VBox
    selections_predicates_wrapper: VBox
    add_controls: GridBox
    complete_add: VBox

    # ── Directional analysis controls ────────────────────────────────────────
    compute_predicates: Button
    compute_predicates_between_selections: Checkbox
    directional_progress_label: ipyw.HTML
    directional_progress_bar: ipyw.IntProgress
    directional_progress_box: VBox
    directional_controls_box: VBox

    # ── Differential expression controls ────────────────────────────────────
    diff_t_threshold: FloatText
    diff_p_threshold: FloatText
    compute_diffexpr: Button
    diff_progress_label: ipyw.HTML
    diff_progress_bar: ipyw.IntProgress
    diff_progress_box: VBox
    diff_controls_box: VBox

    # ── Combined wrapper shown in sidebar ────────────────────────────────────
    compute_predicates_wrapper: VBox

    # ── Plot controls ────────────────────────────────────────────────────────
    color_by: Dropdown

    # ── Result display containers ────────────────────────────────────────────
    pathway_table_container: VBox
    reactome_diagram_container: VBox

    # ── Top-level layout widget ───────────────────────────────────────────────
    plot_wrapper: VBox
    sidebar: GridBox
    ui: VBox


def build_controls(
    scatter: Scatter,
    available_metadata_cols: List[str],
    all_gene_sorted: List[str],
    color_by_default: Optional[str],
    max_genes: int,
    df_columns: List[str],
) -> UIControls:
    """Create and return all UI controls for ScSketch.

    Parameters
    ----------
    scatter:
        The :class:`jscatter.Scatter` instance.
    available_metadata_cols:
        Metadata column names that exist in the DataFrame.
    all_gene_sorted:
        All gene names, alphabetically sorted.
    color_by_default:
        Initial value for the *Color By* dropdown.
    max_genes:
        How many genes were added to the DataFrame (0 = none).
    df_columns:
        Column names of the working DataFrame (for QC options in dropdown).

    Returns
    -------
    :class:`UIControls` with every widget populated and wired into its local
    sub-layout.  Event handlers are **not** attached here — see
    :meth:`ScSketch._setup_handlers`.
    """

    # ── Selection add row ────────────────────────────────────────────────────
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
    selection_num_subdivisions = IntText(value=5, min=2, max=10, step=1, description="Parts")
    selection_subdivide_wrapper = HBox([selection_subdivide, selection_num_subdivisions])

    selections_elements = VBox(layout=Layout(grid_gap="2px"))
    selections_predicates = VBox(layout=Layout(overflow_y="auto", height="100%", grid_gap="6px"))
    selections_predicates_wrapper = VBox([selections_predicates], layout=Layout(height="100%"))

    add_controls = GridBox(
        [selection_name, selection_add],
        layout=Layout(grid_template_columns="1fr 40px"),
    )
    complete_add = VBox([add_controls], layout=Layout(grid_gap="4px"))

    # ── Directional analysis controls ────────────────────────────────────────
    compute_predicates = Button(
        description="Compute Directional Search",
        style="primary",
        disabled=True,
        full_width=True,
    )
    compute_predicates_between_selections = Checkbox(
        value=False, description="Compare Between Selections", indent=False
    )
    directional_progress_label = ipyw.HTML("")
    directional_progress_bar = ipyw.IntProgress(
        value=0, min=0, max=5, description="", layout=Layout(width="100%")
    )
    directional_progress_box = VBox(
        [directional_progress_label, directional_progress_bar],
        layout=Layout(display="none", width="100%", grid_gap="2px"),
    )
    directional_controls_box = VBox([compute_predicates, directional_progress_box])

    # ── Differential expression controls ────────────────────────────────────
    diff_t_threshold = FloatText(value=2.0, step=0.5, description="|T| ≥")
    diff_t_threshold.layout.width = "100%"
    diff_t_threshold.style = {"description_width": "50px"}

    diff_p_threshold = FloatText(value=0.05, step=0.01, description="p ≤")
    diff_p_threshold.layout.width = "100%"
    diff_p_threshold.style = {"description_width": "35px"}

    diff_thresholds = GridBox(
        [diff_t_threshold, diff_p_threshold],
        layout=Layout(grid_template_columns="1fr 1fr", width="100%", min_width="0px", grid_gap="6px"),
    )
    compute_diffexpr = Button(
        description="Compute DE",
        style="primary",
        disabled=True,
        full_width=True,
    )
    diff_progress_label = ipyw.HTML("")
    diff_progress_bar = ipyw.IntProgress(
        value=0, min=0, max=5, description="", layout=Layout(width="100%")
    )
    diff_progress_box = VBox(
        [diff_progress_label, diff_progress_bar],
        layout=Layout(display="none", width="100%", grid_gap="2px"),
    )
    diff_controls_box = VBox(
        [diff_thresholds, compute_diffexpr, diff_progress_box],
        layout=Layout(display="none"),
    )

    compute_predicates_wrapper = VBox([directional_controls_box, diff_controls_box])

    # ── Switch to DE mode if scatter already in freeform mode ───────────────
    if scatter.widget.lasso_type == "freeform":
        directional_controls_box.layout.display = "none"
        diff_controls_box.layout.display = "flex"

    # ── Color-by dropdown ────────────────────────────────────────────────────
    _cat_priority = ["cell_population", "seurat_clusters", "leiden", "louvain", "clusters"]
    meta_in_priority = [c for c in _cat_priority if c in available_metadata_cols]
    meta_rest = sorted(
        [c for c in available_metadata_cols if c not in _cat_priority], key=lambda s: s.lower()
    )
    meta_ordered = meta_in_priority + meta_rest
    meta_opts = [(c.replace("_", " ").title(), c) for c in meta_ordered]

    _qc_order = ["n_genes", "total_counts", "pct_counts_mt"]
    obs_opts = [(c.replace("_", " ").title(), c) for c in _qc_order if c in df_columns]

    gene_options = [(g, g) for g in all_gene_sorted[:max_genes]] if max_genes > 0 else []
    dropdown_options = meta_opts + obs_opts + gene_options

    color_by = Dropdown(
        options=dropdown_options,
        value=color_by_default,
        description="Color By:",
    )

    # ── Result display containers ────────────────────────────────────────────
    pathway_table_container = VBox(
        [],
        layout=Layout(
            overflow_y="auto",
            max_height="400px",
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
            min_width="0px",
            padding="10px",
            display="none",
            border="1px solid #ccc",
        ),
    )

    # ── Sidebar ──────────────────────────────────────────────────────────────
    sidebar = GridBox(
        [complete_add, selections_elements, selections_predicates_wrapper, compute_predicates_wrapper],
        layout=Layout(
            grid_template_rows="min-content max-content 1fr min-content",
            overflow_y="hidden",
            height="100%",
            min_width="0px",
            grid_gap="4px",
        ),
    )

    # ── Plot wrapper ─────────────────────────────────────────────────────────
    plot_wrapper = VBox(
        [scatter.show(), color_by],
        layout=Layout(width="100%", min_width="0px"),
    )

    # ── Suppress Jupyter output-area overflow clipping ───────────────────────
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

    # ── Top-level layout ─────────────────────────────────────────────────────
    combined_gene_pathway_panel = GridBox(
        [
            VBox([sidebar], layout=Layout(overflow_y="auto", height="100%", min_width="0px")),
            VBox(
                [pathway_table_container],
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

    top_layout = GridBox(
        [plot_wrapper, combined_gene_pathway_panel],
        layout=Layout(
            grid_template_columns="2.5fr 2.5fr",
            grid_gap="10px",
            height="auto",
            align_items="flex-start",
            width="100%",
        ),
    )

    ui = VBox(
        [top_layout, VBox([reactome_diagram_container], layout=Layout(width="100%"))],
        layout=Layout(
            grid_gap="10px",
            width="100%",
            overflow_y="auto",
            align_items="flex-start",
        ),
    )

    return UIControls(
        selection_name=selection_name,
        selection_add=selection_add,
        selection_subdivide=selection_subdivide,
        selection_num_subdivisions=selection_num_subdivisions,
        selection_subdivide_wrapper=selection_subdivide_wrapper,
        selections_elements=selections_elements,
        selections_predicates=selections_predicates,
        selections_predicates_wrapper=selections_predicates_wrapper,
        add_controls=add_controls,
        complete_add=complete_add,
        compute_predicates=compute_predicates,
        compute_predicates_between_selections=compute_predicates_between_selections,
        directional_progress_label=directional_progress_label,
        directional_progress_bar=directional_progress_bar,
        directional_progress_box=directional_progress_box,
        directional_controls_box=directional_controls_box,
        diff_t_threshold=diff_t_threshold,
        diff_p_threshold=diff_p_threshold,
        compute_diffexpr=compute_diffexpr,
        diff_progress_label=diff_progress_label,
        diff_progress_bar=diff_progress_bar,
        diff_progress_box=diff_progress_box,
        diff_controls_box=diff_controls_box,
        compute_predicates_wrapper=compute_predicates_wrapper,
        color_by=color_by,
        pathway_table_container=pathway_table_container,
        reactome_diagram_container=reactome_diagram_container,
        plot_wrapper=plot_wrapper,
        sidebar=sidebar,
        ui=ui,
    )


# ── Progress helpers ─────────────────────────────────────────────────────────

def set_analysis_progress(
    controls: UIControls, mode: str, step: int, total: int, message: str
) -> None:
    """Update the progress bar and label for ``mode`` (``"Directional"`` or ``"DE"``)."""
    if mode == "Directional":
        box, label, bar = controls.directional_progress_box, controls.directional_progress_label, controls.directional_progress_bar
    elif mode == "DE":
        box, label, bar = controls.diff_progress_box, controls.diff_progress_label, controls.diff_progress_bar
    else:
        return

    total = int(max(1, total))
    bar.max = total
    bar.value = int(min(max(step, 0), total))
    label.value = f"{PROGRESS_SPINNER_HTML}{mode} ({bar.value}/{total}): {message}"
    box.layout.display = "flex"


def clear_analysis_progress(controls: UIControls, mode: str) -> None:
    """Hide and reset the progress bar for ``mode``."""
    if mode == "Directional":
        box, label, bar = controls.directional_progress_box, controls.directional_progress_label, controls.directional_progress_bar
    elif mode == "DE":
        box, label, bar = controls.diff_progress_box, controls.diff_progress_label, controls.diff_progress_bar
    else:
        return

    label.value = ""
    bar.value = 0
    box.layout.display = "none"
