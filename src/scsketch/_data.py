"""Utilities for building the embedding DataFrame from an AnnData object.

This module is purely concerned with data preparation — no widgets, no UI.
"""

from __future__ import annotations

import logging
from itertools import cycle
from typing import Optional, List

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from jscatter import glasbey_light

logger = logging.getLogger(__name__)


def build_embedding_df(
    adata: AnnData,
    metadata_cols: Optional[List[str]],
    color_by_default: str,
    max_genes: int,
) -> dict:
    """Build the working DataFrame and colour-map metadata from an AnnData object.

    Parameters
    ----------
    adata:
        AnnData object with ``obsm["X_umap"]`` and gene expression in ``X``.
    metadata_cols:
        Column names in ``adata.obs`` to include in the DataFrame.
    color_by_default:
        Preferred column to colour by initially.
    max_genes:
        Maximum number of genes to include as columns (0 = none).

    Returns
    -------
    dict with keys:

    * ``df`` – combined DataFrame (UMAP + metadata + gene expression)
    * ``available_metadata_cols`` – list of metadata cols actually present
    * ``all_gene_sorted`` – alphabetically sorted gene names
    * ``meta_cols_present`` – requested metadata cols found in df
    * ``categorical_cols`` – subset of meta cols treated as categorical
    * ``categorical_color_maps`` – ``{col: {value: hex_color}}``
    * ``color_map_default`` – colour map dict for the default column (or ``None``)
    * ``color_by_default`` – (potentially updated) column name to colour by
    """
    # ── UMAP coordinates ────────────────────────────────────────────────────
    umap_df = pd.DataFrame(
        adata.obsm["X_umap"],
        columns=["x", "y"],
        index=adata.obs_names,
    )

    # ── Metadata columns ────────────────────────────────────────────────────
    if metadata_cols is not None:
        available_metadata_cols = [
            col for col in metadata_cols if col in adata.obs.columns
        ]
        if len(available_metadata_cols) > 0:
            metadata_df = adata.obs[available_metadata_cols].copy(deep=False)
            coerce_cols = [
                col
                for col in available_metadata_cols
                if pd.api.types.is_object_dtype(metadata_df[col])
                or pd.api.types.is_categorical_dtype(metadata_df[col])
            ]
            if coerce_cols:
                metadata_df = metadata_df.assign(
                    **{col: metadata_df[col].astype(str) for col in coerce_cols}
                )
        else:
            logger.info("No requested metadata columns found; continuing without metadata.")
            available_metadata_cols = []
            metadata_df = pd.DataFrame(index=adata.obs_names)
    else:
        available_metadata_cols = []
        metadata_df = pd.DataFrame(index=adata.obs_names)
        logger.info("No metadata passed; continuing with UMAP + gene expression only.")

    # ── Gene expression columns ─────────────────────────────────────────────
    all_gene_sorted = sorted(list(adata.var_names), key=lambda s: s.lower())
    gene_subset = all_gene_sorted[:max_genes] if max_genes > 0 else []
    if gene_subset:
        subX = adata[:, gene_subset].X
        if sp.issparse(subX):
            gene_exp_df = pd.DataFrame.sparse.from_spmatrix(
                subX, columns=gene_subset, index=adata.obs_names
            )
        else:
            gene_exp_df = pd.DataFrame(
                np.asarray(subX), columns=gene_subset, index=adata.obs_names
            )
    else:
        gene_exp_df = pd.DataFrame(index=adata.obs_names)

    # ── Combine ──────────────────────────────────────────────────────────────
    df = pd.concat([umap_df, metadata_df, gene_exp_df], axis=1)
    df = df.loc[:, ~df.columns.duplicated()]

    meta_cols_present = [c for c in (metadata_cols or []) if c in df.columns]
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

    # ── Default colour column ────────────────────────────────────────────────
    color_map_default = None
    if color_by_default in categorical_cols:
        color_map_default = categorical_color_maps[color_by_default]
    else:
        priority_cats = [
            c
            for c in [
                "cell_type",
                "celltype.l1",
                "celltype.l2",
                "cell_population",
                "seurat_clusters",
                "leiden",
                "louvain",
                "clusters",
            ]
            if c in categorical_cols
        ]
        if priority_cats:
            color_by_default = priority_cats[0]
            color_map_default = categorical_color_maps[color_by_default]
        else:
            fallback_cont = next(
                (c for c in ["n_genes", "total_counts", "pct_counts_mt"] if c in df.columns),
                None,
            )
            color_by_default = fallback_cont
            color_map_default = None

    return {
        "df": df,
        "available_metadata_cols": available_metadata_cols,
        "all_gene_sorted": all_gene_sorted,
        "meta_cols_present": meta_cols_present,
        "categorical_cols": categorical_cols,
        "categorical_color_maps": categorical_color_maps,
        "color_map_default": color_map_default,
        "color_by_default": color_by_default,
    }
