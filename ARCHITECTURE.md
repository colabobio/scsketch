# scSketch architecture

This document describes how the interactive scSketch UI is structured and where state lives.

## High-level flow

- The user calls `sketch = ScSketch(adata=adata, ...)` to construct a ScSketch object.
- `sketch.show()` displays a composed ipywidgets UI ready to use.

## Module map

`ScSketch` (`scsketch.py`) is a thin **orchestrator** — it owns selection state and wires
event handlers, but delegates all other responsibilities to single-purpose modules:

| Module | Responsibility |
|---|---|
| `scsketch.py` | Orchestrator: selection management, event wiring, `show()` |
| `_data.py` | `build_embedding_df()` — builds `df` from AnnData; no widgets |
| `_api.py` | All external HTTP — Reactome and MyGene.info |
| `_diffexpr.py` | `DiffExprEngine` — stateful global-stats cache + Welch t-test |
| `_ui.py` | `build_controls()` / `UIControls` — widget construction and layout only |
| `_results.py` | `show_directional_results()`, `show_diffexpr_results()` — result panel wiring |
| `_analysis.py` | Pure statistics: `test_direction`, `lord_test`, numba-accelerated DE kernels |
| `_utils.py` | `Selection`, `Selections`, `Lasso`, geometry helpers, `create_selection` |
| `_logging.py` | `configure_logging()` — module-scoped logger, `LogLevel` type alias |
| `_cli.py` | `scsketch demo` CLI entry point (downloads + launches demo notebook via `uv`) |
| `widgets/` | Eight AnyWidget components (tables, plots, SVG viewer, labels) |

## Key runtime state

All UI state lives on a `ScSketch` instance:

- `self.df` (`pd.DataFrame`)
  - Built by `_data.build_embedding_df()` from the AnnData object.
  - Always contains `x`, `y` (UMAP coords) and any requested metadata columns.
  - May contain a small, optional set of preloaded gene-expression columns when `max_genes > 0` (to enable gene coloring in the dropdown).
  - Does **not** contain all gene expression by default (to avoid huge memory use on large datasets).
  - When `max_genes == 0`, no gene-expression columns are preloaded into `self.df` (analysis still uses all genes from `adata`).
- `self.controls` (`UIControls`)
  - Dataclass holding every named `ipywidgets` widget reference, returned by `_ui.build_controls()`.
  - Keeps the orchestrator decoupled from widget construction details.
- `self.de_engine` (`DiffExprEngine`)
  - Owns the per-dataset global `sum`/`sqsum` cache and all Welch t-test logic.
  - Extracted from `ScSketch` so differential-expression computation has no widget dependencies.
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
  - Uses a sparse-aware correlation implementation (`_analysis.test_direction`) to avoid densifying `adata.X` when it is sparse.
  - Returns a list of per-selection result lists (one entry per selection).
- Rendering happens in `_results.show_directional_results(...)`:
  - Stateless function — takes results + widget refs, builds the gene table, and wires gene-click handlers.
  - Gene click renders a `GeneProjectionPlot` widget for the active selection; expression is loaded for the selected cells only.
- Differential compute is owned by `DiffExprEngine` (`_diffexpr.py`):
  - Uses `adata.raw.X` if present, else `adata.X`.
  - Compares selected cells vs all non-selected cells using Welch t-test computed from summary stats.
  - Maintains a per-dataset cache of global summary statistics (`sum` / `sqsum`) to support fast repeated DE queries.
    - Optionally, these global stats can be persisted to disk via `ScSketch(diffexpr_disk_cache_dir=...)` to avoid
      recomputing them across notebook sessions on large datasets.
    - For SciPy CSR matrices, global `sum`/`sqsum` are computed in one pass over CSR storage (to avoid materializing
      `X.power(2)` for large sparse matrices).
  - For SciPy CSR matrices, per-selection `sum`/`sqsum` are computed directly from CSR storage using
    `_analysis.diffexpr_sum_sqsum_selected_csr` (Numba-accelerated when the `[fast]` extra is installed;
    pure-NumPy fallback otherwise).
- Differential rendering happens in `_results.show_diffexpr_results(...)`:
  - Stateless function — renders a table with `T` and `p`, wires gene-click to a `GeneViolinPlot` widget.

## Widgets

All eight widgets in `src/scsketch/widgets/` follow the same pattern:

```python
_STATIC = Path(__file__).parent.parent / "static" / "widgets"
_esm = _STATIC / "<name>.js"   # anywidget reads at render time
_css = _STATIC / "<name>.css"  # anywidget injects into shadow DOM
```

This means JS and CSS can be edited live during development (with `ANYWIDGET_HMR=1`) without restarting the kernel.

| Widget class | Purpose |
|---|---|
| `CorrelationTable` | Gene correlation results (directional mode) |
| `PathwayTable` | Reactome pathway list for a selected gene |
| `InteractiveSVG` | Clickable Reactome pathway diagram |
| `GeneProjectionPlot` | Scatter plot: projection vs expression (directional gene click) |
| `GeneViolinPlot` | Selected vs background expression distribution (DE gene click) |
| `GenePathwayWidget` | Combined gene search + pathway table |
| `Label` | Styled section header |
| `Div` | Horizontal divider |

## Progress indicators

- `demo.ipynb` dataset download uses `urllib.request.urlretrieve(..., reporthook=...)` to drive an `ipywidgets.IntProgress` (0–100%) while downloading; when the file already exists, the notebook shows a one-line "Found, skipping download" status instead of a persistent full bar.
- scSketch compute feedback is step-based (not byte/gene-level):
  - Directional compute advances coarse steps inside `ScSketch._compute_predicates_handler`.
  - Differential compute advances coarse steps inside `ScSketch._compute_diffexpr_handler`.
  - Progress widget references are held in `UIControls` (`directional_progress_box`, `diff_progress_box`).
  - `_ui.set_analysis_progress()` / `_ui.clear_analysis_progress()` update them; `ScSketch` calls these helpers so it never manipulates widget internals directly.
  - Progress boxes include an animated SVG spinner while active and are hidden when idle.

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
