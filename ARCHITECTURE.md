# scSketch architecture

This document describes how the interactive scSketch UI is structured and where state lives.

## High-level flow

- The user calls `scsketch.view(adata, ...)`.
- `view()` constructs `_ScSketchDirectionalView` and calls `.render()`.
- `.render()` displays a composed ipywidgets UI and returns a `Context` (also stored as the module-level “last context”).

## Core modules

- `src/scsketch/directional.py`
  - Owns the interactive UI (“view”), selection management, and directional analysis.
- `src/scsketch/utils.py`
  - Defines `Selection`, `Selections`, `Lasso`, and small geometry helpers.
- `src/scsketch/widgets/*`
  - AnyWidget/ipywidget components used by the UI (label rows, tables, plots).

## Key runtime state

All UI state lives on a `_ScSketchDirectionalView` instance:

- `self.selections` (`Selections`)
  - The saved selections (each a `Selection` with `points`, `name`, `color`, etc.).
- `self.active_selection` (`Selection | None`)
  - The currently “active” selection.
  - Default behavior: when a new selection is saved, it becomes active (i.e., “latest selection” remains the default).
  - Clicking a saved selection in the sidebar also makes it active.
- `Selection.cached_results` (`list[dict] | None`)
  - Cached per-selection directional results.
  - Stored so users can switch between selections and re-open the previously computed gene list without recomputing.

## Compute vs render

- Compute happens in `_ScSketchDirectionalView._compute_directional_analysis(df, selections)`:
  - Returns a list of per-selection result lists (one entry per selection).
- Rendering happens in `_ScSketchDirectionalView._show_directional_results(directional_results)`:
  - Builds the gene table widget from results and wires the gene-click handlers.

## UI behavior rules

- **Compute target**:
  - If “Compare Between Selections” is enabled, compute runs over all saved selections.
  - Otherwise, compute runs over the active selection (fallback to the latest selection if none is active).
- **Clear Results**:
  - Clears the visible results panel only.
  - Does not delete `Selection.cached_results`.
- **Selection click**:
  - Zooms to the selection and activates it.
  - If cached results exist, restores them immediately; otherwise shows a “No cached results yet” message.

## Layout notes

- The UI uses `GridBox` for the main layout and relies on `min_width="0px"` for grid children to prevent CSS grid “min-content” sizing from shrinking the scatter plot column.
- The sidebar uses `grid_template_rows="min-content max-content 1fr min-content"` with a fixed overall panel height so the middle section can scroll while keeping the compute button visible.
