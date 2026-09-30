"""UI widget construction for ScSketch.

``build_controls()`` creates every ipywidget control used in the sidebar and
returns them bundled in a ``UIControls`` dataclass so callers can reference
them by name.  ``build_layout()`` assembles those controls plus the scatter
plot into the final ``VBox`` tree that ``ScSketch.show()`` returns.

No business logic lives here — only widget instantiation and layout wiring.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import ipywidgets as ipyw
from IPython.display import HTML, display
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
from jscatter.widgets import Button

from jscatter import Scatter

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


# ── Theme (light chrome, dark plot) ─────────────────────────────────────────

_THEME_CSS_PATH = Path(__file__).parent / "static" / "scsketch.css"

# Simple placeholder mark; swap for the final logo asset when available.
_LOGO_HTML = (
    '<div class="sc-logo">'
    '<svg viewBox="0 0 26 26" xmlns="http://www.w3.org/2000/svg" aria-hidden="true">'
    '<rect width="26" height="26" rx="7" fill="#111111"/>'
    '<path d="M5 17c3-8 6-8 8-3s5 3 8-5" fill="none" stroke="#7fb2ff" '
    'stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>'
    "</svg><span>scSketch</span></div>"
)


def _theme_css() -> str:
    """Return the scSketch theme stylesheet ('' if the file is missing)."""
    try:
        return _THEME_CSS_PATH.read_text(encoding="utf-8")
    except OSError:
        return ""


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
    history_dropdown: Dropdown
    gene_expression_caption: ipyw.HTML
    session_filename: Text
    session_save: ipyw.Button
    session_upload: ipyw.FileUpload
    session_status: ipyw.HTML

    # ── Result display containers ────────────────────────────────────────────
    multi_view_toggle: ipyw.ToggleButtons
    multi_view_container: VBox
    pathway_table_container: VBox
    reactome_diagram_container: VBox

    # ── Top-level layout widget ───────────────────────────────────────────────
    plot_wrapper: VBox
    sidebar: VBox
    ui: VBox


def build_controls(
    scatter: Scatter,
    available_metadata_cols: List[str],
    all_gene_sorted: List[str],
    color_by_default: Optional[str],
    max_genes: int,
    df_columns: List[str],
    extra_views: Optional[Dict[str, Scatter]] = None,
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
    extra_views:
        Optional mapping of display labels to additional :class:`jscatter.Scatter`
        instances shown in the right detail panel when multi-view mode is enabled.

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
    compute_predicates_wrapper.add_class("sc-compute")

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
        description="Color by",
    )
    color_by.layout.width = "auto"
    color_by.layout.min_width = "0px"
    color_by.style = {"description_width": "auto"}
    color_by.add_class("sc-colorby")

    history_dropdown = Dropdown(
        options=[("No recorded actions", None)],
        value=None,
        description="History",
        disabled=True,
    )
    history_dropdown.layout.min_width = "0px"
    history_dropdown.style = {"description_width": "0px"}
    history_dropdown.add_class("sc-history")

    gene_expression_caption = ipyw.HTML(
        "",
        layout=Layout(display="none", width="100%", min_width="0px"),
    )
    gene_expression_caption.add_class("sc-caption")
    # ── Session import/export controls ───────────────────────────────────────
    session_filename = Text(
        value="scsketch-session.scsketch.json",
        description="",
    )
    session_filename.layout.min_width = "0px"
    session_filename.add_class("sc-session-name")
    session_save = ipyw.Button(
        description="Save",
        icon="",
        tooltip="Save scSketch session JSON",
        layout=Layout(width="auto", min_width="0px"),
    )
    session_save.add_class("sc-btn")
    session_save.add_class("sc-btn-save")
    session_upload = ipyw.FileUpload(
        accept=".json,.scsketch.json",
        multiple=False,
        description="Upload Session",
        icon="",
        tooltip="Load scSketch session JSON",
        layout=Layout(width="auto", min_width="0px"),
    )
    session_upload.add_class("sc-btn")
    session_upload.add_class("sc-btn-upload")
    session_status = ipyw.HTML(
        "",
        layout=Layout(
            display="none",
            width="100%",
            min_width="0px",
            overflow="hidden",
        ),
    )

    # ── Result display containers ────────────────────────────────────────────
    extra_views = extra_views or {}
    has_extra_views = len(extra_views) > 0
    multi_view_toggle = ipyw.ToggleButtons(
        options=[("OFF", False), ("ON", True)],
        value=has_extra_views,
        description="Multi-view",
        disabled=not has_extra_views,
        layout=Layout(width="auto", display="flex" if has_extra_views else "none"),
        style={"description_width": "80px"},
    )
    extra_view_children: list[ipyw.Widget] = []
    for label, view in extra_views.items():
        extra_view_children.extend(
            [
                ipyw.HTML(f"<b>{label}</b>"),
                view.show(),
            ]
        )
    multi_view_container = VBox(
        extra_view_children,
        layout=Layout(
            display="flex" if has_extra_views else "none",
            flex_direction="column",
            height="100%",
            min_height="0px",
            width="100%",
            min_width="0px",
            overflow="auto",
            grid_gap="8px",
        ),
    )
    pathway_table_container = VBox(
        [],
        layout=Layout(
            overflow_y="auto",
            max_height="400px",
            border="1px solid #e3e3df",
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
            border="1px solid #e3e3df",
        ),
    )

    # ── Theme ────────────────────────────────────────────────────────────────
    # Injected once per build; scoped to ``.sc-root`` so nothing leaks out.
    display(
        HTML(
            "<style>\n"
            + _theme_css()
            + "\n.jp-OutputArea-output, .jp-Cell-outputArea, .jp-Notebook {\n"
            "    overflow: auto !important;\n    max-height: none !important;\n}\n"
            "</style>"
        )
    )

    def _section_title(text: str) -> ipyw.HTML:
        return ipyw.HTML(f'<div class="sc-section-title">{text}</div>')

    # ── Top bar: logo | session name + Save + Upload | History ──────────────
    logo = ipyw.HTML(_LOGO_HTML)
    session_label = ipyw.HTML('<span class="sc-session-label">Session name</span>')
    session_group = HBox(
        [session_label, session_filename, session_save, session_upload],
        layout=Layout(min_width="0px"),
    )
    session_group.add_class("sc-session")
    top_bar = HBox([logo, session_group, history_dropdown])
    top_bar.add_class("sc-topbar")

    # ── Plot (dark) ──────────────────────────────────────────────────────────
    plot_wrapper = VBox(
        [color_by, scatter.show(), gene_expression_caption, session_status],
        layout=Layout(width="100%", min_width="0px", overflow="hidden"),
    )
    plot_wrapper.add_class("sc-plot")

    # ── Sidebar: SELECTIONS | GENE | compute ─────────────────────────────────
    selections_section = VBox(
        [_section_title("Selections"), complete_add, selections_elements],
        layout=Layout(flex="0 0 auto"),
    )
    divider = VBox([], layout=Layout(width="100%"))
    divider.add_class("sc-divider")
    genes_section = VBox(
        [_section_title("Gene"), selections_predicates_wrapper],
        layout=Layout(flex="1 1 0px", min_height="0px", overflow="auto"),
    )
    genes_section.add_class("sc-sidebar-genes")
    sidebar = VBox(
        [selections_section, divider, genes_section, compute_predicates_wrapper],
        layout=Layout(height="100%", min_width="0px", overflow="hidden"),
    )
    sidebar.add_class("sc-sidebar")

    workspace = GridBox(
        [plot_wrapper, sidebar],
        layout=Layout(
            grid_template_columns="minmax(0, 1fr) minmax(300px, 410px)",
            width="100%",
            min_width="0px",
            align_items="stretch",
        ),
    )
    workspace.add_class("sc-workspace")

    # ── Results dock: gene details (plot + pathways) and Reactome diagram ────
    dock_details = VBox(
        [multi_view_toggle, multi_view_container, pathway_table_container],
        layout=Layout(flex="1 1 0px", min_width="0px", grid_gap="8px"),
    )
    dock_row = HBox([dock_details], layout=Layout(width="100%", min_width="0px"))
    dock_row.add_class("sc-dock-row")
    dock = VBox(
        [dock_row, reactome_diagram_container],
        layout=Layout(width="100%", min_width="0px"),
    )
    dock.add_class("sc-dock")

    # ── Top-level layout ─────────────────────────────────────────────────────
    ui = VBox(
        [top_bar, workspace, dock],
        layout=Layout(width="100%", min_width="0px"),
    )
    ui.add_class("sc-root")

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
        history_dropdown=history_dropdown,
        gene_expression_caption=gene_expression_caption,
        session_filename=session_filename,
        session_save=session_save,
        session_upload=session_upload,
        session_status=session_status,
        multi_view_toggle=multi_view_toggle,
        multi_view_container=multi_view_container,
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
