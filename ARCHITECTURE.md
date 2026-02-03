# scSketch architecture

This document describes how the interactive scSketch UI is structured and where state lives.

## High-level flow

- The user calls `sketch = ScSketch(adata=adata,...` to construct an ScSketch oject.
- `sketch.show()` displays a composed ipywidgets UI ready to use.

## Core modules

- `src/scsketch/scsketch.py`
  - Owns the interactive UI (“view”), selection management, and directional analysis.
- `src/scsketch/utils.py`
  - Defines `Selection`, `Selections`, `Lasso`, and small geometry helpers.
- `src/scsketch/widgets/*`
  - AnyWidget/ipywidget components used by the UI (label rows, tables, plots).

## Key runtime state

All UI state lives on a `ScSketch` instance:

- `self.df` (`pd.DataFrame`)
  - The main table passed to `jscatter.Scatter`.
  - Always contains `x`, `y` (UMAP coords) and any requested metadata columns.
  - May contain a small, optional set of preloaded gene-expression columns when `max_genes > 0` (to enable gene coloring in the dropdown).
  - Does **not** contain all gene expression by default (to avoid huge memory use on large datasets).
  - When `max_genes == 0`, no gene-expression columns are preloaded into `self.df` (analysis still uses all genes from `adata`).
- `self.selections` (`Selections`)
  - The saved selections (each a `Selection` with `points`, `name`, `color`, etc.).
- `self.active_selection` (`Selection | None`)
  - The currently “active” selection.
  - Default behavior: when a new selection is saved, it becomes active (i.e., “latest selection” remains the default).
  - Clicking a saved selection in the sidebar also makes it active.
- `Selection.cached_results` (`list[dict] | None`)
  - Cached per-selection directional results.
  - Stored so users can switch between selections and re-open the previously computed gene list without recomputing.
- `Selection.cached_diffexpr` (`list[dict] | None`)
  - Cached per-selection differential-expression results (Welch t-test; selected vs background).
  - Stored so users can switch between saved selections and re-open the previous DE gene list without recomputing.
- `self.analysis_mode` (`"directional"` | `"differential"`)
  - Current analysis mode.
  - `lasso_type == "freeform"` enables differential mode and hides the directional compute UI.
  - Other lasso types use directional mode.
- `self._pending_diffexpr` (`dict | None`)
  - Temporary cache for differential results computed from an unsaved “current selection”.
  - When the user saves that selection, the cached DE results are attached to the created `Selection.cached_diffexpr`.

## Compute vs render

- Directional compute happens in `ScSketch._compute_directional_analysis(df, selections)`:
  - Uses `df[["x","y"]]` for geometry/projections.
  - Pulls gene-expression from `adata.X` for the selected cells, so it can analyze all genes without preloading them into `self.df`.
  - Uses a sparse-aware correlation implementation (`test_direction`) to avoid densifying `adata.X` when it is sparse.
  - Returns a list of per-selection result lists (one entry per selection).
- Rendering happens in `ScSketch._show_directional_results(directional_results)`:
  - Builds the gene table widget from results and wires the gene-click handlers.
  - Gene click renders a projection-vs-expression plot for the active selection; expression is loaded for the selected cells only (to avoid loading 1M+ values).
- Differential compute happens in `ScSketch._compute_diffexpr(selected_indices, selection_label)`:
  - Uses `adata.raw.X` if present, else `adata.X`.
  - Compares selected cells vs all non-selected cells using Welch t-test computed from summary stats.
- Differential rendering happens in `ScSketch._show_diffexpr_results(diff_results, ...)`:
  - Renders a table with `T` and `p` (Sciviewer-like).
  - Gene click renders a compact violin-like plot of selected vs background (sampled for plotting).

## Progress indicators

- `demo.ipynb` dataset download uses `urllib.request.urlretrieve(..., reporthook=...)` to drive an `ipywidgets.IntProgress` (0–100%) while downloading; when the file already exists, the notebook shows a one-line “Found, skipping download” status instead of a persistent full bar.
- scSketch compute feedback is step-based (not byte/gene-level):
  - Directional compute advances coarse steps inside `ScSketch._compute_predicates_handler`.
  - Differential compute advances coarse steps inside `ScSketch._compute_diffexpr_handler`.
  - Progress widgets live next to the compute buttons, include an animated spinner while active, and are hidden when idle.
- scSketch directional selections (brush mode) show a visual direction arrow overlay:
  - The arrow is derived from the brush “spine” (the midpoint line through the brush polygon).
  - This arrow is a visual cue only; directional compute is unchanged.

## UI behavior rules

- **Compute target**:
  - If “Compare Between Selections” is enabled, compute runs over all saved selections.
  - Otherwise, compute runs over the active selection (fallback to the latest selection if none is active).
- **Differential mode** (`lasso_type == "freeform"`):
  - Directional controls are hidden and DE controls are shown (thresholds, `Compute DE`).
- **Subdivide / Parts UI**:
  - The subdivide controls are currently hidden (selections are saved as a single selection).
- **Clear Results**:
  - Clears the visible results panel only.
  - Does not delete `Selection.cached_results`.
- **Selection click**:
  - Zooms to the selection and activates it.
  - In directional mode: restores `Selection.cached_results` if present.
  - In differential mode: restores `Selection.cached_diffexpr` if present.
  - Otherwise shows a “No cached results yet” message for the current mode.

## Layout notes

- The UI uses `GridBox` for the main layout and relies on `min_width="0px"` for grid children to prevent CSS grid “min-content” sizing from shrinking the scatter plot column.
- The sidebar uses `grid_template_rows="min-content max-content 1fr min-content"` with a fixed overall panel height so the middle section can scroll while keeping the compute button visible.
