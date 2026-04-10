"""Results display helpers for ScSketch.

Functions in this module take analysis results plus the widget references
they need and build/wire the interactive results panels.  They contain no
stateful class logic — all state is passed in explicitly.
"""

from __future__ import annotations

import logging
from typing import Callable, Optional

import numpy as np
import pandas as pd
import ipywidgets as ipyw
from ipywidgets import HTML as _HTML, Layout, VBox

from anndata import AnnData

from ._api import fetch_pathways, fetch_pathway_svg
from .widgets import (
    CorrelationTable,
    PathwayTable,
    InteractiveSVG,
    GeneProjectionPlot,
    GeneViolinPlot,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Directional analysis results
# ─────────────────────────────────────────────────────────────────────────────

def show_directional_results(
    directional_results: list,
    *,
    selections_predicates: ipyw.VBox,
    pathway_table_container: ipyw.VBox,
    reactome_diagram_container: ipyw.VBox,
    df: pd.DataFrame,
    adata: AnnData,
    active_selection,           # Selection | None
    on_results_cleared: Callable,
    log: Callable,
) -> None:
    """Populate the sidebar with gene correlation results and wire pathway interactions.

    Parameters
    ----------
    directional_results:
        List-of-lists returned by
        :meth:`ScSketch._compute_directional_analysis`.
    selections_predicates:
        The ``VBox`` that holds the gene table widget.
    pathway_table_container:
        The container ``VBox`` for the pathway table / projection plot.
    reactome_diagram_container:
        The container ``VBox`` for the Reactome SVG diagram.
    df:
        The working DataFrame (UMAP + metadata).
    adata:
        The AnnData object (for gene expression when not in ``df``).
    active_selection:
        Currently focused selection (may be ``None``).
    on_results_cleared:
        Callback invoked when the user clicks *Clear Results*.
    log:
        Debug logging callable (e.g. ``self._log``).
    """
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
                    "Selection": entry.get("direction", f"Selection {i + 1}"),
                }
            )

    results_df = pd.DataFrame(all_results).dropna(subset=["R", "p"])
    results_df = results_df.sort_values(by="R", ascending=False).reset_index(drop=True)

    gene_table_widget = CorrelationTable(data=results_df.to_dict(orient="records"))
    pathway_table_widget = PathwayTable(data=[])

    pathway_table_container.layout.display = "none"
    reactome_diagram_container.layout.display = "none"

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
    pathway_box = VBox([pathway_msg, pathway_table_widget], layout=Layout(flex="1 1 auto", overflow="auto"))
    pathway_table_container.children = [plot_box, pathway_box]
    pathway_table_container.layout = Layout(
        display="flex",
        flex_direction="column",
        height="420px",
        max_height="420px",
        overflow="visible",
    )

    def on_gene_click(change):
        gene = change["new"]
        log(f"[UI] gene clicked: {gene!r}")
        try:
            pathways = fetch_pathways(gene)
            pathway_table_widget.data = pathways
            pathway_table_container.layout.display = "block"
            pathway_msg.value = (
                f"<em>No Reactome pathways found for <b>{gene}</b>.</em>"
                if len(pathways) == 0
                else f"<b>Reactome pathways for {gene}</b>"
            )
            reactome_diagram_container.layout.display = "none"

            if not active_selection or len(active_selection.points) == 0:
                log("No selections available for gene projection plot.")
                return

            sel = active_selection
            selected_indices = sel.points
            log(f"[plot] selection name={sel.name!r} n={len(selected_indices)}")

            X_umap = df[["x", "y"]].to_numpy(dtype=float)
            Ux = X_umap[:, 0] - X_umap[:, 0].min()
            Uy = X_umap[:, 1] - X_umap[:, 1].min()
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
                proj = np.dot(pts - pts[0], v)
                span = float(proj.max() - proj.min()) or 1.0
                proj = (proj - proj.min()) / span
            else:
                celv = U_sel - U_sel[0]
                proj = (celv @ dirv) / (L2 + 1e-12)
                proj = np.clip(np.nan_to_num(proj, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)

            if gene in df.columns:
                expr = pd.to_numeric(df.iloc[selected_indices][gene], errors="coerce").to_numpy(dtype=float)
                log(f"[plot] gene {gene!r} found in df columns (n={len(expr)})")
            else:
                sub = adata[selected_indices, [gene]].X
                if hasattr(sub, "toarray"):
                    sub = sub.toarray()
                expr = np.asarray(sub).ravel().astype(float)
                log(f"[plot] gene {gene!r} NOT in df; loaded from adata (n={len(expr)})")

            if not np.isfinite(expr).all():
                finite = np.isfinite(expr)
                if finite.any():
                    expr = np.where(finite, expr, float(np.median(expr[finite])))
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
            log(f"[plot] {gene}: n={len(plot_data)} proj=[{proj.min():.3f},{proj.max():.3f}]")
        except Exception:
            logger.exception("Error handling gene click")

    interactive_svg_widget = InteractiveSVG()

    def on_pathway_click(change):
        pathway_id = change["new"]
        if not pathway_id:
            return
        svg_b64 = fetch_pathway_svg(pathway_id)
        if svg_b64 is None:
            return
        interactive_svg_widget.svg_content = svg_b64
        reactome_diagram_container.children = [interactive_svg_widget]
        reactome_diagram_container.layout.display = "block"
        log(f"[SVG] Loaded Reactome diagram for {pathway_id}")

    gene_table_widget.observe(on_gene_click, names=["selected_gene"])
    pathway_table_widget.observe(on_pathway_click, names=["selected_pathway"])
    selections_predicates.children = [gene_table_widget]
    log("Showing directional results...")


# ─────────────────────────────────────────────────────────────────────────────
# Differential expression results
# ─────────────────────────────────────────────────────────────────────────────

def show_diffexpr_results(
    diff_results: list[dict],
    selection_label: str,
    *,
    selected_indices: Optional[np.ndarray],
    selections_predicates: ipyw.VBox,
    pathway_table_container: ipyw.VBox,
    reactome_diagram_container: ipyw.VBox,
    de_source_fn: Callable,      # () -> (X, var_names, source_tag)
    active_selection,            # Selection | None
    scatter,                     # Scatter | None
    log: Callable,
) -> None:
    """Populate the sidebar with differential expression results and violin plots.

    Parameters
    ----------
    diff_results:
        List of dicts from :meth:`DiffExprEngine.compute`.
    selection_label:
        Human-readable name for the result set.
    selected_indices:
        Cell indices of the foreground group (may be ``None``).
    selections_predicates:
        VBox that will contain the gene table.
    pathway_table_container:
        Container for the violin/histogram plot.
    reactome_diagram_container:
        Container for the Reactome diagram (hidden for DE mode).
    de_source_fn:
        Callable returning ``(X, var_names, source_tag)`` — typically
        ``DiffExprEngine.de_source``.
    active_selection:
        Currently focused selection.
    scatter:
        The jscatter Scatter instance (for fallback selection lookup).
    log:
        Debug logging callable.
    """
    if not diff_results:
        _clear_results(
            selections_predicates,
            pathway_table_container,
            reactome_diagram_container,
            message=f"<em>No differential results passed thresholds for <b>{selection_label}</b>.</em>",
        )
        return

    rows = [
        {
            "Gene": e["attribute"],
            "T": float(np.round(float(e["interval"][0]), 4)),
            "p": f"{float(e['interval'][1]):.3e}",
            "Selection": selection_label,
        }
        for e in diff_results
    ]
    results_df = pd.DataFrame(rows).dropna(subset=["T", "p"])
    results_df["_absT"] = results_df["T"].abs()
    results_df = results_df.sort_values(by="_absT", ascending=False).drop(columns=["_absT"]).reset_index(drop=True)

    gene_table_widget = CorrelationTable(
        data=results_df.to_dict(orient="records"),
        columns=["Gene", "T", "p", "Selection"],
    )

    violin = GeneViolinPlot()
    title = _HTML("")
    plot_box = VBox(
        [title, violin],
        layout=Layout(
            flex="0 0 auto",
            height="270px",
            max_height="270px",
            min_height="270px",
            padding="0px",
            overflow="hidden",
            align_self="stretch",
            display="none",
        ),
    )
    pathway_table_container.children = [plot_box]
    pathway_table_container.layout.display = "block"
    reactome_diagram_container.layout.display = "none"

    def on_gene_click(change):
        gene = change["new"]
        if not gene:
            return
        try:
            X, var_names, _ = de_source_fn()
            if gene not in var_names:
                return
            gene_idx = int(var_names.index(gene))

            if selected_indices is not None:
                sel = np.unique(np.asarray(selected_indices, dtype=int))
            elif active_selection is not None:
                sel = np.unique(np.asarray(active_selection.points, dtype=int))
            elif scatter is not None:
                sel = np.unique(np.asarray(scatter.selection(), dtype=int))
            else:
                return

            n_obs = int(X.shape[0])
            if sel.size < 1 or sel.size >= n_obs:
                return

            max_samp = 5000
            rng = np.random.default_rng(0)
            sel_samp = rng.choice(sel, size=min(sel.size, max_samp), replace=False)
            mask = np.ones(n_obs, dtype=bool)
            mask[sel] = False
            bg_idx = np.flatnonzero(mask)
            bg_samp = rng.choice(bg_idx, size=min(bg_idx.size, max_samp), replace=False)

            def _extract(rows):
                col = X[rows, gene_idx]
                return np.asarray(col.toarray() if hasattr(col, "toarray") else col).ravel().astype(float)

            sel_vals = _extract(sel_samp)
            bg_vals = _extract(bg_samp)

            vmin = float(min(sel_vals.min() if sel_vals.size else 0.0, bg_vals.min() if bg_vals.size else 0.0))
            vmax = float(max(sel_vals.max() if sel_vals.size else 0.0, bg_vals.max() if bg_vals.size else 0.0))
            if not (np.isfinite(vmin) and np.isfinite(vmax)):
                return
            if vmax <= vmin:
                vmax = vmin + 1e-6

            bins = np.linspace(vmin, vmax, 51)
            sel_hist, _ = np.histogram(sel_vals, bins=bins, density=True)
            bg_hist, _ = np.histogram(bg_vals, bins=bins, density=True)

            violin.gene = gene
            violin.data = {
                "bins": bins.tolist(),
                "selected": sel_hist.tolist(),
                "background": bg_hist.tolist(),
            }
            title.value = (
                f"<b>{gene}</b> — selected n={sel.size}, "
                f"background n={n_obs - sel.size} (sampled)"
            )
            plot_box.layout.display = "block"
        except Exception:
            logger.exception("Error rendering differential plot")

    gene_table_widget.observe(on_gene_click, names=["selected_gene"])
    selections_predicates.children = [gene_table_widget]


# ─────────────────────────────────────────────────────────────────────────────
# Clear helper
# ─────────────────────────────────────────────────────────────────────────────

def clear_results(
    selections_predicates: ipyw.VBox,
    pathway_table_container: ipyw.VBox,
    reactome_diagram_container: ipyw.VBox,
    message: Optional[str] = None,
) -> None:
    """Clear all result panels, optionally showing a status ``message``."""
    _clear_results(selections_predicates, pathway_table_container, reactome_diagram_container, message)


def _clear_results(
    selections_predicates: ipyw.VBox,
    pathway_table_container: ipyw.VBox,
    reactome_diagram_container: ipyw.VBox,
    message: Optional[str] = None,
) -> None:
    if message:
        selections_predicates.children = (ipyw.HTML(message),)
    else:
        selections_predicates.children = ()
    pathway_table_container.children = ()
    pathway_table_container.layout.display = "none"
    reactome_diagram_container.children = ()
    reactome_diagram_container.layout.display = "none"
