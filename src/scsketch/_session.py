"""Session-log import/export helpers for ScSketch."""

from __future__ import annotations

import hashlib
import json
from html import escape as _escape_html
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import ipywidgets as ipyw
import numpy as np

from jscatter import Line

from ._action_log import (
    action_log_from_steps,
    apply_action_log_entry,
    normalize_action_log,
)
from ._utils import Selection, Selections

SESSION_SCHEMA_VERSION = "1.0"


def export_session(sketch) -> dict[str, Any]:
    """Return a JSON-serializable session document for ``sketch``."""
    adata = sketch.adata
    return {
        "schema_version": SESSION_SCHEMA_VERSION,
        "scsketch_version": _scsketch_version(),
        "dataset": _dataset_fingerprint(adata),
        "initial_config": _initial_config(sketch),
        "state": {
            "analysis_mode": getattr(sketch, "analysis_mode", "directional"),
            "active_selection": (
                None
                if getattr(sketch, "active_selection", None) is None
                else sketch.active_selection.name
            ),
            "differential": _differential_state(sketch),
            "multi_view": _multi_view_state(sketch),
        },
        "selections": [
            _serialize_selection(selection, adata)
            for selection in sketch.selections.selections
        ],
        "selection_archive": [
            _serialize_selection(selection, adata)
            for selection in _archived_selections(sketch)
        ],
        "steps": _steps_for_selections(sketch.selections.selections),
        "action_log": normalize_action_log(getattr(sketch, "action_log", [])),
    }


def write_session(sketch, path: str | Path) -> dict[str, Any]:
    """Write ``sketch`` session state to ``path`` and return the document."""
    document = export_session(sketch)
    path = Path(path)
    path.write_text(json.dumps(document, indent=2), encoding="utf-8")
    return document


def read_session(path: str | Path) -> dict[str, Any]:
    """Read a session document from ``path``."""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_session(sketch, session: dict[str, Any]) -> list[str]:
    """Load a session document into ``sketch``.

    Returns a list of warnings for non-fatal compatibility issues.
    """
    warnings = validate_session(sketch, session)
    selections = []
    for raw in session.get("selections", []):
        selections.append(_deserialize_selection(raw, sketch.adata))

    sketch.selections = Selections(selections=selections)
    sketch._selection_archive = {
        selection.name: selection
        for selection in _deserialize_selection_archive(session, sketch.adata)
    }
    active_name = session.get("state", {}).get("active_selection")
    sketch.active_selection = _find_selection(selections, active_name)
    if sketch.active_selection is None and selections:
        sketch.active_selection = selections[-1]

    state = session.get("state", {})
    if "analysis_mode" in state:
        sketch.analysis_mode = state["analysis_mode"]
    sketch.action_log = normalize_action_log(
        session.get("action_log") or action_log_from_steps(session.get("steps") or [])
    )
    _restore_widget_state(sketch, state)

    _refresh_loaded_session_ui(sketch)
    return warnings


def build_playback_steps(session: dict[str, Any]) -> list[dict[str, Any]]:
    """Return ordered playback steps for a session document."""
    action_log = normalize_action_log(session.get("action_log") or [])
    if action_log:
        return action_log

    steps = session.get("steps") or []
    if steps:
        return list(steps)

    selections = [_selection_stub(raw) for raw in session.get("selections", [])]
    return _steps_for_selections(selections)


def apply_playback_step(
    sketch,
    session: dict[str, Any],
    step_index: int,
) -> dict[str, Any]:
    """Apply one playback step to ``sketch`` and return that step."""
    steps = build_playback_steps(session)
    if not steps:
        raise ValueError("Session log does not contain any playback steps.")
    if step_index < 0 or step_index >= len(steps):
        raise IndexError("Playback step index is out of range.")

    step = steps[step_index]
    if _is_action_log_step(step):
        selections = _selections_through_action(
            session,
            sketch.adata,
            step_index,
            step,
        )
        sketch.selections = Selections(
            selections=selections
        )
        _restore_widget_state(sketch, session.get("state", {}))
        message = apply_action_log_entry(sketch, step)
        _refresh_selection_sidebar(sketch)
        return {**step, "_status": message}

    selections = _selections_through_step(session, sketch.adata, step)
    sketch.selections = Selections(selections=selections)
    sketch.active_selection = _find_selection(selections, step.get("selection"))
    if sketch.active_selection is None and selections:
        sketch.active_selection = selections[-1]

    step_type = step.get("type")
    if step_type == "show_diffexpr_results":
        sketch.analysis_mode = "differential"
    else:
        sketch.analysis_mode = "directional"

    _restore_widget_state(sketch, session.get("state", {}))
    _refresh_loaded_session_ui(sketch)
    _show_playback_step(sketch, step)
    return step


def build_session_player(sketch, session: dict[str, Any]) -> ipyw.VBox:
    """Build a small ipywidgets playback panel for a session document."""
    warnings = validate_session(sketch, session)
    steps = build_playback_steps(session)

    title = ipyw.HTML("<b>scSketch Session Playback</b>")
    step_label = ipyw.HTML("")
    status = ipyw.HTML("")
    warning_box = ipyw.HTML(_format_warnings(warnings))
    prev_button = ipyw.Button(description="Previous", icon="chevron-left")
    next_button = ipyw.Button(description="Next", icon="chevron-right")

    state = {"index": 0}

    def render() -> None:
        if not steps:
            step_label.value = "No playback steps found."
            status.value = ""
            prev_button.disabled = True
            next_button.disabled = True
            return

        index = state["index"]
        step = apply_playback_step(sketch, session, index)
        step_label.value = (
            f"<b>Step {index + 1} of {len(steps)}</b>: "
            f"{_escape_html(step.get('label') or step.get('type') or '')}"
        )
        if step.get("_status"):
            status.value = f"<em>{_escape_html(step['_status'])}</em>"
        else:
            status.value = _step_status(step)
        prev_button.disabled = index == 0
        next_button.disabled = index == len(steps) - 1

    def previous(_event) -> None:
        state["index"] = max(0, state["index"] - 1)
        render()

    def next_step(_event) -> None:
        state["index"] = min(len(steps) - 1, state["index"] + 1)
        render()

    prev_button.on_click(previous)
    next_button.on_click(next_step)
    render()

    return ipyw.VBox(
        [
            title,
            ipyw.HBox([prev_button, next_button]),
            step_label,
            status,
            warning_box,
        ],
        layout=ipyw.Layout(
            border="1px solid #ddd",
            padding="8px",
            margin="8px 0",
            width="100%",
        ),
    )


def validate_session(sketch, session: dict[str, Any]) -> list[str]:
    """Return warnings describing session/dataset compatibility issues."""
    warnings = []
    schema_version = session.get("schema_version")
    if schema_version != SESSION_SCHEMA_VERSION:
        warnings.append(
            f"Session schema version {schema_version!r} may not be compatible "
            f"with {SESSION_SCHEMA_VERSION!r}."
        )

    current = _dataset_fingerprint(sketch.adata)
    recorded = session.get("dataset", {})
    for key in ("n_obs", "n_vars", "obs_names_hash", "var_names_hash"):
        if recorded.get(key) != current.get(key):
            warnings.append(f"Dataset {key} differs from the session log.")
    if (
        recorded.get("embedding_hash") is not None
        and recorded.get("embedding_hash") != current.get("embedding_hash")
    ):
        warnings.append("Dataset embedding coordinates differ from the session log.")
    return warnings


def _scsketch_version() -> str:
    try:
        return version("scsketch")
    except PackageNotFoundError:
        return "unknown"


def _initial_config(sketch) -> dict[str, Any]:
    return {
        "metadata_cols": list(getattr(sketch, "metadata_cols", []) or []),
        "color_by_default": getattr(sketch, "color_by_default", None),
        "height": getattr(sketch, "height", None),
        "background_color": getattr(sketch, "background_color", None),
        "max_genes": getattr(sketch, "max_genes", None),
        "fdr_alpha": getattr(sketch, "fdr_alpha", None),
    }


def _differential_state(sketch) -> dict[str, Any]:
    ctrl = getattr(sketch, "_ctrl", None)
    if ctrl is None:
        return {}
    return {
        "t_threshold": getattr(ctrl.diff_t_threshold, "value", None),
        "p_threshold": getattr(ctrl.diff_p_threshold, "value", None),
    }


def _multi_view_state(sketch) -> bool | None:
    ctrl = getattr(sketch, "_ctrl", None)
    if ctrl is None or not hasattr(ctrl, "multi_view_toggle"):
        return None
    return bool(ctrl.multi_view_toggle.value)


def _dataset_fingerprint(adata) -> dict[str, Any]:
    embedding = None
    if hasattr(adata, "obsm") and "X_umap" in adata.obsm:
        embedding = np.asarray(adata.obsm["X_umap"])
    return {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obs_names_hash": _hash_strings(adata.obs_names),
        "var_names_hash": _hash_strings(adata.var_names),
        "embedding_key": "X_umap" if embedding is not None else None,
        "embedding_hash": None if embedding is None else _hash_array(embedding),
    }


def _hash_strings(values) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def _hash_array(values: np.ndarray) -> str:
    arr = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(arr.dtype).encode("ascii"))
    digest.update(str(arr.shape).encode("ascii"))
    digest.update(arr.tobytes())
    return digest.hexdigest()


def _serialize_selection(selection: Selection, adata) -> dict[str, Any]:
    points = np.asarray(selection.points, dtype=int)
    return {
        "index": int(selection.index),
        "name": selection.name,
        "color": _jsonable(selection.color),
        "points_indices": points.tolist(),
        "points_obs_names": [str(adata.obs_names[i]) for i in points],
        "lasso_polygon": _line_points(selection.lasso),
        "hull": _line_points(selection.hull),
        "path": _array_or_none(selection.path),
        "cached_results": _jsonable(selection.cached_results),
        "cached_diffexpr": _jsonable(selection.cached_diffexpr),
    }


def _steps_for_selections(selections) -> list[dict[str, Any]]:
    steps = []
    for selection in selections:
        steps.append(
            {
                "type": "restore_selection",
                "selection": selection.name,
                "label": f"Restore selection: {selection.name}",
            }
        )
        if getattr(selection, "cached_results", None) is not None:
            steps.append(
                {
                    "type": "show_directional_results",
                    "selection": selection.name,
                    "label": f"Show directional results: {selection.name}",
                }
            )
        if getattr(selection, "cached_diffexpr", None) is not None:
            steps.append(
                {
                    "type": "show_diffexpr_results",
                    "selection": selection.name,
                    "label": f"Show differential results: {selection.name}",
                }
            )
    return steps


def _selection_stub(raw: dict[str, Any]):
    return type(
        "SessionSelectionStub",
        (),
        {
            "name": raw.get("name", ""),
            "cached_results": raw.get("cached_results"),
            "cached_diffexpr": raw.get("cached_diffexpr"),
        },
    )()


def _deserialize_selection(raw: dict[str, Any], adata) -> Selection:
    points = _restore_points(raw, adata)
    color = raw.get("color") or "#000000"
    lasso_polygon = raw.get("lasso_polygon") or []
    hull = raw.get("hull") or lasso_polygon
    return Selection(
        index=int(raw.get("index", 0)),
        name=str(raw.get("name", "")),
        points=points,
        color=color,
        lasso=Line(lasso_polygon, line_color=color, line_width=2),
        hull=Line(hull, line_color=color, line_width=2),
        path=None if raw.get("path") is None else np.asarray(raw["path"], dtype=float),
        cached_results=raw.get("cached_results"),
        cached_diffexpr=raw.get("cached_diffexpr"),
    )


def _restore_points(raw: dict[str, Any], adata) -> np.ndarray:
    obs_names = raw.get("points_obs_names") or []
    if obs_names:
        index_by_name = {str(name): i for i, name in enumerate(adata.obs_names)}
        mapped = [index_by_name.get(str(name)) for name in obs_names]
        if all(index is not None for index in mapped):
            return np.asarray(mapped, dtype=int)

    points = np.asarray(raw.get("points_indices", []), dtype=int)
    if points.size and (points.min() < 0 or points.max() >= int(adata.n_obs)):
        raise ValueError(
            f"Selection {raw.get('name')!r} has out-of-bounds cell indices."
        )
    return points


def _selections_through_step(
    session: dict[str, Any],
    adata,
    step: dict[str, Any],
) -> list[Selection]:
    target_name = step.get("selection")
    selections = []
    for raw in session.get("selections", []):
        selections.append(_deserialize_selection(raw, adata))
        if raw.get("name") == target_name:
            break
    return selections


def _selections_through_action(
    session: dict[str, Any],
    adata,
    step_index: int,
    step: dict[str, Any],
) -> list[Selection]:
    action_log = normalize_action_log(session.get("action_log") or [])
    saved_names = []
    removed_names = set()
    for entry in action_log[: step_index + 1]:
        entry_type = entry.get("type")
        selection = entry.get("selection")
        if entry_type == "save_selection" and selection:
            removed_names.discard(selection)
            saved_names.append(selection)
        elif entry_type == "remove_selection" and selection:
            removed_names.add(selection)

    if step.get("selection") and step["selection"] not in saved_names:
        saved_names.append(step["selection"])

    names = [name for name in saved_names if name not in removed_names]
    raw_selections = session.get("selection_archive") or session.get("selections", [])
    raw_by_name = {
        raw.get("name"): raw
        for raw in raw_selections
        if raw.get("name") is not None
    }
    return [
        _deserialize_selection(raw_by_name[name], adata)
        for name in names
        if name in raw_by_name
    ]


def _is_action_log_step(step: dict[str, Any]) -> bool:
    return "payload" in step or "timestamp" in step


def _show_playback_step(sketch, step: dict[str, Any]) -> None:
    active = getattr(sketch, "active_selection", None)
    if active is None or getattr(sketch, "_ctrl", None) is None:
        return

    step_type = step.get("type")
    if step_type == "restore_selection":
        if hasattr(sketch, "_clear_results_display"):
            sketch._clear_results_display(
                f"<em>Playback restored <b>{_escape_html(active.name)}</b>.</em>"
            )
    elif step_type == "show_directional_results":
        if active.cached_results is not None and hasattr(
            sketch, "_show_directional_results"
        ):
            sketch._show_directional_results([active.cached_results])
    elif step_type == "show_diffexpr_results":
        if active.cached_diffexpr is not None and hasattr(
            sketch, "_show_diffexpr_results"
        ):
            sketch._show_diffexpr_results(
                active.cached_diffexpr,
                active.name,
                selected_indices=np.asarray(active.points, dtype=int),
            )


def _format_warnings(warnings: list[str]) -> str:
    if not warnings:
        return ""
    items = "".join(f"<li>{_escape_html(warning)}</li>" for warning in warnings)
    return f"<b>Compatibility warnings</b><ul>{items}</ul>"


def _step_status(step: dict[str, Any]) -> str:
    selection = _escape_html(str(step.get("selection") or ""))
    step_type = step.get("type")
    if step_type == "restore_selection":
        return f"Selection active: <b>{selection}</b>"
    if step_type == "show_directional_results":
        return f"Showing cached directional results for <b>{selection}</b>."
    if step_type == "show_diffexpr_results":
        return f"Showing cached differential results for <b>{selection}</b>."
    return ""


def _line_points(line) -> list[list[float]]:
    if line is None:
        return []
    for attr in ("points", "vertices", "data"):
        value = getattr(line, attr, None)
        if value is not None:
            return np.asarray(value, dtype=float).tolist()
    for attr in ("_points", "_vertices", "_data"):
        value = getattr(line, attr, None)
        if value is not None:
            return np.asarray(value, dtype=float).tolist()
    if isinstance(line, (list, tuple, np.ndarray)):
        return np.asarray(line, dtype=float).tolist()
    raise TypeError(f"Cannot serialize line object of type {type(line).__name__}.")


def _array_or_none(values) -> list[list[float]] | None:
    if values is None:
        return None
    return np.asarray(values, dtype=float).tolist()


def _jsonable(value):
    if value is None:
        return None
    return json.loads(json.dumps(value))


def _find_selection(selections: list[Selection], name: str | None) -> Selection | None:
    if name is None:
        return None
    for selection in selections:
        if selection.name == name:
            return selection
    return None


def _restore_widget_state(sketch, state: dict[str, Any]) -> None:
    ctrl = getattr(sketch, "_ctrl", None)
    if ctrl is None:
        return

    differential = state.get("differential", {})
    if differential.get("t_threshold") is not None:
        ctrl.diff_t_threshold.value = differential["t_threshold"]
    if differential.get("p_threshold") is not None:
        ctrl.diff_p_threshold.value = differential["p_threshold"]

    if getattr(sketch, "analysis_mode", "directional") == "differential":
        ctrl.directional_controls_box.layout.display = "none"
        ctrl.diff_controls_box.layout.display = "flex"
    else:
        ctrl.directional_controls_box.layout.display = "flex"
        ctrl.diff_controls_box.layout.display = "none"

    if state.get("multi_view") is not None and hasattr(ctrl, "multi_view_toggle"):
        was_paused = getattr(sketch, "_action_log_paused", False)
        sketch._action_log_paused = True
        try:
            ctrl.multi_view_toggle.value = bool(state["multi_view"])
        finally:
            sketch._action_log_paused = was_paused
        if hasattr(sketch, "_apply_multi_view_visibility"):
            sketch._apply_multi_view_visibility()


def _refresh_selection_sidebar(sketch) -> None:
    ctrl = getattr(sketch, "_ctrl", None)
    if ctrl is None:
        return

    ctrl.selections_elements.children = ()
    for selection in sketch.selections.selections:
        sketch._add_selection_element(selection)

    n_selections = len(sketch.selections.selections)
    ctrl.compute_predicates.disabled = n_selections == 0
    if n_selections > 1:
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

    if hasattr(sketch, "_update_annotations"):
        sketch._update_annotations()


def _refresh_loaded_session_ui(sketch) -> None:
    ctrl = getattr(sketch, "_ctrl", None)
    if ctrl is None:
        return

    _refresh_selection_sidebar(sketch)

    active = getattr(sketch, "active_selection", None)
    if active is None:
        if hasattr(sketch, "_clear_results_display"):
            sketch._clear_results_display(None)
    elif getattr(sketch, "analysis_mode", "directional") == "differential":
        if active.cached_diffexpr is not None and hasattr(
            sketch, "_show_diffexpr_results"
        ):
            sketch._show_diffexpr_results(
                active.cached_diffexpr,
                active.name,
                selected_indices=np.asarray(active.points, dtype=int),
            )
        elif hasattr(sketch, "_clear_results_display"):
            sketch._clear_results_display(
                f"<em>No cached differential results for <b>{active.name}</b> yet.</em>"
            )
    elif active.cached_results is not None and hasattr(
        sketch, "_show_directional_results"
    ):
        sketch._show_directional_results([active.cached_results])
    elif hasattr(sketch, "_clear_results_display"):
        sketch._clear_results_display(
            f"<em>No cached results for <b>{active.name}</b> yet.</em>"
        )


def _archived_selections(sketch) -> list[Selection]:
    selections: list[Selection] = []
    seen = set()
    archive = getattr(sketch, "_selection_archive", {}) or {}
    for selection in archive.values():
        if selection is not None and selection.name not in seen:
            selections.append(selection)
            seen.add(selection.name)
    for selection in getattr(sketch.selections, "selections", []):
        if selection.name not in seen:
            selections.append(selection)
            seen.add(selection.name)
    return selections


def _deserialize_selection_archive(
    session: dict[str, Any],
    adata,
) -> list[Selection]:
    raw_selections = session.get("selection_archive") or session.get("selections", [])
    return [_deserialize_selection(raw, adata) for raw in raw_selections]
