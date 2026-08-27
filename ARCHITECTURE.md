# scSketch architecture

This document describes how the interactive scSketch UI is structured and where state lives.

## High-level flow

- The user calls `sketch = ScSketch(adata=adata, ...)` to construct a ScSketch object.
- `sketch.show()` displays a composed ipywidgets UI ready to use.

## Manuscript revision analyses

Revision-specific notebooks live in `Manuscript/revision_analyses/`. They are
analysis artifacts rather than package runtime code:

- `directional_differential_video_demo.ipynb` records the core interactive
  scSketch workflow for the supplementary video.
- `umap_pca_linked_view_suppfig.ipynb` builds the linked UMAP/PCA diagnostic
  figure used to show that a saved scSketch selection can be inspected in a
  second projection.
- `paga_monocle_tutorial_baseline.ipynb` runs a scSketch-independent Scanpy
  PAGA/DPT baseline on the Monocle tutorial AnnData counts and metadata, exports
  cluster/pseudotime diagnostics, produces a PAGA/DPT gene ranking for later
  comparison, and optionally visualizes top pseudotime-associated genes along an
  inferred PAGA path with `scanpy.pl.paga_path`. It also includes optional
  bottom cells for launching scSketch on the PAGA-initialized Scanpy UMAP and
  exporting that session for downstream comparison. After scSketch selections
  are exported, the notebook can map each saved selection to Leiden/PAGA
  clusters, infer a connected PAGA path, rank genes by DPT pseudotime within
  that path, and write per-selection overlap summaries. A final summary section
  condenses these outputs into a reviewer-facing CSV table and PNG/PDF figure.
- `paga_scsketch_comparison.ipynb` runs a Leiden-based PAGA analysis from the
  AnnData neighbor graph, maps a saved scSketch session selection onto the PAGA
  clusters, infers a connected PAGA path through the selected clusters, and
  compares a PAGA/DPT path gene ranking against cached scSketch directional
  genes.

These notebooks consume files in `Manuscript/revision_analyses/data/` and write
derived CSV/PNG outputs under `Manuscript/revision_analyses/outputs/`. They
should remain self-contained so reviewer-response analyses are reproducible
without changing the public `scsketch` API.

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
| `_action_log.py` | Structured user-action logging and clickable audit-log table |
| `_session.py` | Versioned session-log export/import, dataset fingerprints, selection rehydration |
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
  - `Selection.points` is the authoritative saved cell set, copied from `jscatter.selection()` when the user clicks save.
  - Saved selection overlays render the original lasso/brush outline (`Selection.lasso`), not the convex hull of the selected points.
  - Brush selections also render a simple direction guide from `Selection.path`: a spine plus small arrowhead pointing from the brush start toward the brush end.
  - `Selection.hull` is still stored as derived geometry for compatibility, but it is not used to decide which cells belong to a selection.
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
- `self.action_log` (`list[dict]`)
  - Structured audit trail of meaningful user actions recorded by `ScSketch` event handlers.
  - Entries include an index, timestamp, action type, human-readable label, optional selection name, analysis mode, and JSON-safe payload.
  - The action log is exported in session JSON and restored on load.
  - `ScSketch.show_action_log()` renders it as a clickable table; clicking a row jumps to that action in the same reconstructed timeline used by session playback.
  - The main UI also exposes the same jump behavior through the `History:` dropdown beside `Color By:`.
- `self._selection_archive` (`dict[str, Selection]`)
  - Keeps saved selection snapshots even if they are later removed from the visible selection list.
  - Exported as `selection_archive` so playback can reconstruct historical selection state while moving backward and forward.
- `self.extra_views` (`dict[str, jscatter.Scatter]`)
  - Optional mapping supplied through `ScSketch(extra_views={"PCA": pca_scatter})`.
  - Extra views are caller-built `jscatter.Scatter` instances, so users can provide PCA, tSNE, PHATE, diffusion-map, or other 2D views without scSketch computing those embeddings.
  - Extra views are assumed to have the same row order as the main `adata` / `self.df`.
  - When present, the right detail panel starts in multi-view mode and displays the extra views instead of the gene-detail plot.
  - Extra views are normalized to a compact square size before rendering in the right panel.
  - Main-scatter selection and extra-view selections are synchronized by integer point index.
  - After the user saves a selection, scSketch clears the transient main-scatter selection but keeps the active saved selection highlighted in extra views.
  - Metadata color changes propagate to extra views when their data contains the selected metadata column; categorical columns reuse scSketch's category-to-color map.
  - Gene-result clicks recolor extra views with the same expression vector and continuous color scale used in the main embedding.
- `self.gene_annotation_species` (`str | int`)
  - Species filter passed to MyGene.info for gene annotation lookups.
  - Defaults to `"human"` so existing notebooks keep the same behavior unless a caller opts into another species.
- `self.reactome_species` (`str | int`)
  - Species filter passed to Reactome pathway lookups.
  - Defaults to `"human"`; internally the Reactome helper maps `"human"` / `"Homo sapiens"` to taxon `9606` to preserve the previous hard-coded behavior.

## Compute vs render

- Directional compute happens in `ScSketch._compute_directional_analysis(df, selections)`:
  - Uses `df[["x","y"]]` for geometry/projections.
  - For brush selections, uses the stored brush centerline (`selection.path`) from stroke start to stroke end as the projection axis; otherwise falls back to the selected-cell endpoint order.
  - The brush path is reconstructed from the stored brush outline by splitting the outline polygon into two halves, reversing the second half, and averaging the matched left/right boundary points into a centerline.
    ```text
    outline order:   A1 -> A2 -> A3 -> A4 -> B4 -> B3 -> B2 -> B1
    pairing used:    (A1,B1), (A2,B2), (A3,B3), (A4,B4)
    centerline pts:    *       *       *       *
    ```
  - Pulls gene-expression from `adata.X` for the selected cells, so it can analyze all genes without preloading them into `self.df`.
  - Uses a sparse-aware correlation implementation (`_analysis.test_direction`) to avoid densifying `adata.X` when it is sparse.
  - Returns a list of per-selection result lists (one entry per selection).
- Rendering happens in `_results.show_directional_results(...)`:
  - Stateless function — takes results + widget refs, builds the gene table, and wires gene-click handlers.
  - The visible gene table reports a `Discovery Score`, an integer 0-10 ranking score derived from the nominal p-value. The score is intended for prioritizing genes/features during exploratory analysis and is not interpreted as a valid post-selection p-value.
  - Result tables display readable gene labels from symbol-like `adata.var` / `adata.raw.var` columns when available, but keep the original `adata.var_names` identifier in a hidden `_gene_id` field for click handling and expression lookup.
  - Gene click recolors the main embedding by full-dataset expression using a blue-to-green continuous gradient; captions and color legend labels use the readable display name when available while expression lookup still uses the original gene ID.
  - Gene click fetches a cached MyGene.info annotation (`symbol`, `name`, `summary`, IDs) via `_api.fetch_gene_description()` and renders it above the pathway list with a scrollable summary area. The lookup uses `ScSketch(gene_annotation_species=...)`, which defaults to `"human"`. Ensembl-like IDs use MyGene's direct gene endpoint, WormBase IDs (`WBGene...`) use a `wormbase:` query, and other IDs are treated as symbols. If the internal gene ID lookup fails and a readable `adata.var` display label exists, the results panel retries the annotation lookup with that display label.
  - The same gene click also queries Reactome through `_api.fetch_pathways()` using `ScSketch(reactome_species=...)`, which defaults to `"human"` / Reactome taxon `9606`.
  - The gene click also renders a `GeneProjectionPlot` widget for the active selection using the same path-based projection direction as the analysis; expression is loaded for the selected cells only.
- Differential compute is owned by `DiffExprEngine` (`_diffexpr.py`):
  - Uses `adata.raw.X` if present, else `adata.X`.
  - Compares selected cells vs all non-selected cells using Welch t-test computed from summary stats.
  - Returns all genes passing the active `|T|` and `p` thresholds, sorted by `|T|`; there is no hardcoded top-N cap.
  - Maintains a per-dataset cache of global summary statistics (`sum` / `sqsum`) to support fast repeated DE queries.
    - Optionally, these global stats can be persisted to disk via `ScSketch(diffexpr_disk_cache_dir=...)` to avoid
      recomputing them across notebook sessions on large datasets.
    - For SciPy CSR matrices, global `sum`/`sqsum` are computed in one pass over CSR storage (to avoid materializing
      `X.power(2)` for large sparse matrices).
  - For SciPy CSR matrices, per-selection `sum`/`sqsum` are computed directly from CSR storage using
    `_analysis.diffexpr_sum_sqsum_selected_csr` (Numba-accelerated when the `[fast]` extra is installed;
    pure-NumPy fallback otherwise).
- Differential rendering happens in `_results.show_diffexpr_results(...)`:
  - Stateless function — renders a table with `T` and `Discovery Score`, wires gene-click to embedding recoloring, cached MyGene.info annotation rendering using `ScSketch(gene_annotation_species=...)` with the same internal-ID-then-display-label fallback as directional results, and a `GeneViolinPlot` widget.

## Multi-view mode

- `ScSketch(extra_views=...)` enables an optional multi-view panel.
- The public API expects a dictionary mapping display labels to prebuilt `jscatter.Scatter` instances:
  ```python
  sketch = ScSketch(
      adata=adata,
      metadata_cols=["cell_type"],
      color_by_default="cell_type",
      extra_views={"PCA": pca_scatter},
  )
  ```
- `_ui.build_controls()` adds a right-panel `Multi-view` OFF/ON segmented toggle only when extra views are provided.
- When the toggle is ON:
  - The right panel shows the extra scatter view(s).
  - Saved active selections remain highlighted in the extra view(s), even after the main view clears its transient selected-point state.
  - Gene clicks still recolor the main embedding.
  - Gene clicks also recolor the extra scatter view(s) by the same gene expression values, which supports visual comparison of expression gradients across embeddings.
  - Gene projection, violin, pathway table, and Reactome detail panels are hidden so the extra view stays visible.
- When the toggle is OFF:
  - The right panel returns to the existing gene-detail behavior.
  - Directional gene clicks can show `GeneProjectionPlot`.
  - Differential gene clicks can show `GeneViolinPlot`.
- Multi-view synchronization and gene-expression coloring are intentionally index-based for the first implementation; callers should build extra views from the same cells in the same order as the `AnnData` passed to `ScSketch`.

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

## Public result export

- `ScSketch.get_genes(selection_name)` exports cached directional results from `Selection.cached_results` with `gene`, `correlation`, `p-value`, and `discovery_score` columns.
- `ScSketch.get_diffexpr_genes(selection_name)` exports cached differential-expression results from `Selection.cached_diffexpr` with `gene`, `t-statistic`, `p-value`, `discovery_score`, and `selection` columns.
- `ScSketch.get_de_genes(selection_name)` is a short alias for `get_diffexpr_genes(...)`.

## Session logs

- `ScSketch.export_session()` returns a JSON-serializable session document.
- `ScSketch.export_session(path)` writes that document to disk.
- `ScSketch.load_session(session_or_path)` loads a session document or JSON file and returns non-fatal compatibility warnings.
- `ScSketch.get_action_log()` returns the recorded action log as a DataFrame.
- `ScSketch.show_action_log()` returns a clickable audit-log table for the current sketch.
- `ScSketch.show_session_player(session_or_path)` returns a small ipywidgets playback panel for stepping through a saved session.
- The main UI includes a session filename field plus `Save` and `Load` controls. `Save` writes the session JSON to the notebook working directory or user-provided path; `Load` accepts an uploaded `.scsketch.json` file and restores it.
- Session logs are versioned with `schema_version == "1.0"` and are owned by `_session.py`.
- A session log stores:
  - Dataset fingerprints: `n_obs`, `n_vars`, `obs_names_hash`, `var_names_hash`, and an `X_umap` coordinate hash when present.
  - Initial ScSketch config: metadata columns, default color, height, background, `max_genes`, and directional `fdr_alpha`.
  - UI analysis state: active selection, analysis mode, multi-view toggle state, and DE thresholds when widgets are available.
  - Saved selections: name, index, color, selected cell indices, selected `obs_names`, lasso polygon, hull, path, and cached directional/DE results.
  - Selection archive: all selections needed to replay history, including selections that were later deleted from the current visible state.
  - Action log entries recorded while the user works: selection save/focus/remove, mode/color/brush-size changes, threshold changes, compare toggles, computes, result clears, gene/pathway clicks, and session exports.
  - Playback steps synthesized from saved selections and cached result availability.
- On load, cell identity is restored by `obs_names` first and falls back to saved integer indices when names are unavailable.
- Loading a session replaces the current saved selections, restores selection overlays, restores active selection, and re-renders cached results when a full widget instance is available.
- Playback is intentionally conservative:
  - `restore_selection` steps restore selections up to the named selection and make that selection active.
  - `show_directional_results` steps restore the relevant selection and show its cached directional results.
  - `show_diffexpr_results` steps restore the relevant selection and show its cached DE results.
  - `select_gene` steps restore the relevant selection, recolor the embedding, and replay the gene-specific detail panel when cached results are available.
  - Action-log playback reconstructs the visible selection list at each step from the archived selections, so rewinding before a `remove_selection` action can bring that selection back.
  - The clickable action-log table and built-in `History:` dropdown use the same reconstruction path as the session player, so row clicks, dropdown jumps, and step-by-step playback stay consistent.
  - If an older session log has no explicit `steps` field, `_session.py` synthesizes the same step sequence from saved selections.
- Compatibility warnings are informational. They do not block load, because a user may intentionally replay a session against a reordered or closely related AnnData object.
- Live event recording is intentionally limited to reproducibility-relevant UI actions. Low-level pointer movement while drawing a lasso is not logged separately; the saved selection stores the resulting points and shape.

## Layout notes

- The UI uses `GridBox` for the main layout and relies on `min_width="0px"` for grid children to prevent CSS grid “min-content” sizing from shrinking the scatter plot column.
- The sidebar uses `grid_template_rows="min-content max-content 1fr min-content"` with a fixed overall panel height so the middle section can scroll while keeping the compute button visible.
