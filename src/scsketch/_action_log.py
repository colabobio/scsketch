"""Structured user-action logging for ScSketch sessions."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from html import escape
from typing import Any

import ipywidgets as ipyw
import numpy as np

ACTION_LOG_SCHEMA_VERSION = "1.0"


def append_action(
    action_log: list[dict[str, Any]],
    action_type: str,
    label: str,
    *,
    selection: str | None = None,
    analysis_mode: str | None = None,
    payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Append and return one JSON-serializable action-log entry."""
    entry = {
        "index": len(action_log) + 1,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "type": action_type,
        "label": label,
        "selection": selection,
        "analysis_mode": analysis_mode,
        "payload": _jsonable(payload or {}),
    }
    action_log.append(entry)
    return entry


def normalize_action_log(entries) -> list[dict[str, Any]]:
    """Return action-log entries with stable indexes and JSON-safe payloads."""
    normalized = []
    for i, raw in enumerate(entries or [], start=1):
        entry = dict(raw)
        entry["index"] = int(entry.get("index") or i)
        entry.setdefault("timestamp", "")
        entry.setdefault("type", "unknown")
        entry.setdefault("label", entry["type"])
        entry.setdefault("selection", None)
        entry.setdefault("analysis_mode", None)
        entry["payload"] = _jsonable(entry.get("payload") or {})
        normalized.append(entry)
    return normalized


def action_log_from_steps(steps) -> list[dict[str, Any]]:
    """Build action-log entries from older session playback steps."""
    action_log = []
    for step in steps or []:
        append_action(
            action_log,
            step.get("type") or "unknown",
            step.get("label") or step.get("type") or "Unknown action",
            selection=step.get("selection"),
            analysis_mode=_mode_for_action(step.get("type")),
            payload={"source": "playback_step"},
        )
    return action_log


def build_action_log_table(sketch, *, apply_entry=None) -> ipyw.VBox:
    """Return a clickable action-log table widget."""
    status = ipyw.HTML("")
    rows_box = ipyw.VBox(layout=ipyw.Layout(grid_gap="2px"))

    def render_rows() -> None:
        rows = normalize_action_log(getattr(sketch, "action_log", []))
        if not rows:
            rows_box.children = (ipyw.HTML("<em>No recorded actions yet.</em>"),)
            return

        header = _row(
            [
                "#",
                "Time",
                "Action",
                "Selection",
                "Details",
            ],
            header=True,
        )
        rendered_rows = [header]
        for entry in rows:
            rendered_rows.append(_action_row(sketch, entry, status, apply_entry))
        rows_box.children = tuple(rendered_rows)

    refresh = ipyw.Button(description="Refresh", icon="refresh", tooltip="Refresh log")
    refresh.on_click(lambda _event: render_rows())
    title = ipyw.HTML("<b>scSketch Action Log</b>")
    render_rows()
    return ipyw.VBox(
        [ipyw.HBox([title, refresh]), rows_box, status],
        layout=ipyw.Layout(
            border="1px solid #ddd",
            padding="8px",
            margin="8px 0",
            width="100%",
        ),
    )


def apply_action_log_entry(sketch, entry: dict[str, Any]) -> str:
    """Apply an action-log row as far as possible to the current sketch UI."""
    action_type = entry.get("type")
    global_message = _apply_global_action(sketch, entry)
    if global_message is not None:
        return global_message

    selection_name = entry.get("selection")
    selection = _find_selection(sketch, selection_name)

    if selection is None:
        return f"Action {entry.get('index')} has no saved selection to restore."

    sketch.active_selection = selection
    _restore_selection_visuals(sketch, selection)

    mode = entry.get("analysis_mode") or _mode_for_action(action_type)
    if mode:
        _set_mode(sketch, mode)

    can_render = getattr(sketch, "_ctrl", None) is not None
    if action_type in {"compute_directional", "show_directional_results"}:
        if can_render and selection.cached_results is not None and hasattr(
            sketch, "_show_directional_results"
        ):
            sketch._show_directional_results([selection.cached_results])
            return f"Showing directional results for {selection.name}."
    elif action_type in {"compute_diffexpr", "show_diffexpr_results"}:
        if can_render and selection.cached_diffexpr is not None and hasattr(
            sketch, "_show_diffexpr_results"
        ):
            sketch._show_diffexpr_results(
                selection.cached_diffexpr,
                selection.name,
                selected_indices=np.asarray(selection.points, dtype=int),
            )
            return f"Showing differential results for {selection.name}."
    elif action_type in {"save_selection", "focus_selection"}:
        _show_cached_results(sketch, selection, mode)
        return f"Restored selection {selection.name}."
    elif action_type == "select_gene":
        gene = entry.get("payload", {}).get("gene")
        _show_cached_results(sketch, selection, mode, initial_gene=gene)
        if gene and getattr(sketch, "_ctrl", None) is None and hasattr(
            sketch, "_color_embedding_by_gene"
        ):
            sketch._color_embedding_by_gene(gene)
            return f"Restored selection {selection.name} and gene {gene}."
        if gene:
            return f"Restored selection {selection.name} and gene {gene}."
        return f"Restored selection {selection.name}."
    elif action_type == "select_pathway":
        pathway_id = entry.get("payload", {}).get("pathway_id")
        if pathway_id:
            _show_cached_results(sketch, selection, mode)
            if _show_pathway_svg(sketch, pathway_id):
                return f"Showing pathway {pathway_id} for {selection.name}."
            return f"Restored selection {selection.name}; pathway was {pathway_id}."

    return f"Restored selection {selection.name}."


def _apply_global_action(sketch, entry: dict[str, Any]) -> str | None:
    action_type = entry.get("type")
    payload = entry.get("payload") or {}
    ctrl = getattr(sketch, "_ctrl", None)
    scatter = getattr(sketch, "scatter", None)

    if action_type == "change_lasso_type":
        lasso_type = payload.get("lasso_type")
        if scatter is not None and lasso_type:
            _with_paused_log(
                sketch,
                lambda: setattr(scatter.widget, "lasso_type", lasso_type),
            )
        if lasso_type == "freeform":
            _set_mode(sketch, "differential")
        elif lasso_type:
            _set_mode(sketch, "directional")
        return f"Restored lasso type {lasso_type}."

    if action_type == "change_brush_size":
        brush_size = payload.get("brush_size")
        if scatter is not None and brush_size is not None:
            _with_paused_log(
                sketch,
                lambda: setattr(scatter.widget, "lasso_brush_size", brush_size),
            )
        return f"Restored brush size {brush_size}."

    if action_type == "change_color_by":
        color_by = payload.get("color_by")
        if ctrl is not None and color_by is not None:
            _with_paused_log(sketch, lambda: setattr(ctrl.color_by, "value", color_by))
        return f"Restored color by {color_by}."

    if action_type == "toggle_compare_between_selections":
        value = bool(payload.get("compare_between_selections"))
        if ctrl is not None:
            _with_paused_log(
                sketch,
                lambda: setattr(
                    ctrl.compute_predicates_between_selections,
                    "value",
                    value,
                ),
            )
        return f"Restored compare between selections to {value}."

    if action_type == "toggle_multi_view":
        value = bool(payload.get("multi_view"))
        if ctrl is not None and hasattr(ctrl, "multi_view_toggle"):
            _with_paused_log(
                sketch,
                lambda: setattr(ctrl.multi_view_toggle, "value", value),
            )
        if hasattr(sketch, "_apply_multi_view_visibility"):
            sketch._apply_multi_view_visibility()
        return f"Restored multi-view to {value}."

    if action_type == "change_diffexpr_threshold":
        value = payload.get("value")
        field = str(payload.get("field") or "")
        if ctrl is not None and value is not None:
            target = ctrl.diff_t_threshold if "T" in field else ctrl.diff_p_threshold
            _with_paused_log(sketch, lambda: setattr(target, "value", value))
        return "Restored differential-expression threshold."

    if action_type == "clear_results":
        if hasattr(sketch, "_clear_results_display"):
            sketch._clear_results_display(None)
        return "Cleared visible results."

    if action_type in {"export_session", "load_session"}:
        return str(entry.get("label") or action_type)

    return None


def _action_row(
    sketch,
    entry: dict[str, Any],
    status: ipyw.HTML,
    apply_entry,
) -> ipyw.GridBox:
    action_button = ipyw.Button(
        description=str(entry.get("type") or "action"),
        tooltip=str(entry.get("label") or ""),
        layout=ipyw.Layout(width="100%"),
    )

    def on_click(_event) -> None:
        message = (
            apply_entry(entry)
            if apply_entry is not None
            else apply_action_log_entry(sketch, entry)
        )
        status.value = f"<em>{escape(message)}</em>"

    action_button.on_click(on_click)

    return _row(
        [
            str(entry.get("index") or ""),
            _short_time(entry.get("timestamp")),
            action_button,
            entry.get("selection") or "",
            entry.get("label") or "",
        ]
    )


def _row(cells: list[Any], *, header: bool = False) -> ipyw.GridBox:
    widgets = []
    for cell in cells:
        if isinstance(cell, ipyw.Widget):
            widgets.append(cell)
        else:
            safe = escape(str(cell))
            value = f"<b>{safe}</b>" if header else safe
            widgets.append(ipyw.HTML(value))
    return ipyw.GridBox(
        widgets,
        layout=ipyw.Layout(
            grid_template_columns="48px 130px 180px 160px minmax(240px, 1fr)",
            grid_gap="6px",
            align_items="center",
            width="100%",
        ),
    )


def _short_time(timestamp: str | None) -> str:
    if not timestamp:
        return ""
    return str(timestamp).replace("T", " ").split("+", maxsplit=1)[0]


def _find_selection(sketch, name: str | None):
    if not name:
        return None
    for selection in getattr(sketch.selections, "selections", []):
        if selection.name == name:
            return selection
    return None


def _restore_selection_visuals(sketch, selection) -> None:
    scatter = getattr(sketch, "scatter", None)
    if scatter is not None:
        try:
            _with_paused_log(
                sketch,
                lambda: scatter.selection(np.asarray(selection.points, dtype=int)),
            )
            scatter.zoom(to=selection.points, animation=500, padding=2)
        except Exception:
            if hasattr(sketch, "logger"):
                sketch.logger.exception("Failed to restore action-log selection")
    if hasattr(sketch, "_update_annotations"):
        sketch._update_annotations()


def _show_cached_results(
    sketch,
    selection,
    mode: str | None,
    *,
    initial_gene: str | None = None,
) -> None:
    if getattr(sketch, "_ctrl", None) is None:
        return

    def render() -> None:
        if mode == "differential" and selection.cached_diffexpr is not None:
            sketch._show_diffexpr_results(
                selection.cached_diffexpr,
                selection.name,
                selected_indices=np.asarray(selection.points, dtype=int),
                initial_gene=initial_gene,
            )
        elif selection.cached_results is not None:
            sketch._show_directional_results(
                [selection.cached_results],
                initial_gene=initial_gene,
            )

    _with_paused_log(sketch, render)


def _show_pathway_svg(sketch, pathway_id: str) -> bool:
    ctrl = getattr(sketch, "_ctrl", None)
    if ctrl is None:
        return False
    try:
        from ._api import fetch_pathway_svg
        from .widgets import InteractiveSVG

        svg_b64 = fetch_pathway_svg(pathway_id)
        if svg_b64 is None:
            return False
        widget = InteractiveSVG()
        widget.svg_content = svg_b64
        ctrl.reactome_diagram_container.children = [widget]
        ctrl.reactome_diagram_container.layout.display = "block"
        return True
    except Exception:
        if hasattr(sketch, "logger"):
            sketch.logger.exception("Failed to restore Reactome pathway diagram")
        return False


def _set_mode(sketch, mode: str) -> None:
    sketch.analysis_mode = mode
    ctrl = getattr(sketch, "_ctrl", None)
    if ctrl is None:
        return
    if mode == "differential":
        ctrl.directional_controls_box.layout.display = "none"
        ctrl.diff_controls_box.layout.display = "flex"
    else:
        ctrl.directional_controls_box.layout.display = "flex"
        ctrl.diff_controls_box.layout.display = "none"


def _with_paused_log(sketch, fn) -> None:
    previous = getattr(sketch, "_action_log_paused", False)
    sketch._action_log_paused = True
    try:
        fn()
    finally:
        sketch._action_log_paused = previous


def _mode_for_action(action_type: str | None) -> str | None:
    if action_type in {"compute_diffexpr", "show_diffexpr_results"}:
        return "differential"
    if action_type in {"compute_directional", "show_directional_results"}:
        return "directional"
    return None


def _jsonable(value):
    return json.loads(json.dumps(value))
