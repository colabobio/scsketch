#!/usr/bin/env python3
"""
Performance retest harness for scSketch.

This script benchmarks:
  - Directional Search compute: correlation+p-values (test_direction)
  - Online FDR update: LORD (lord_test), timed separately
  - Differential Expression compute: Welch t-test from summary stats

It intentionally excludes:
  - Reactome/network calls
  - Widget/UI rendering
  - Gene-click plots
"""

from __future__ import annotations

import argparse
import csv
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import anndata as ad
import scipy.sparse as sp
import scipy.stats as ss

try:
    import numba as nb
except Exception:  # pragma: no cover
    nb = None

from scsketch.analysis import (
    diffexpr_sum_sqsum_csr_global,
    diffexpr_sum_sqsum_selected_csr,
    lord_test,
    test_direction,
)


@dataclass(frozen=True)
class DatasetSpec:
    label: str
    path: Path
    prefer_obsm: str | None = None


@dataclass(frozen=True)
class ExportedSelection:
    dataset: str
    name: str
    indices: np.ndarray


DATASETS: list[DatasetSpec] = [
    DatasetSpec(label="pbmc3k", path=Path("pbmc3k.h5ad"), prefer_obsm="X_umap"),
    DatasetSpec(label="pbmc161k", path=Path("nygc_pbmc_161k_lite.h5ad"), prefer_obsm="X_wnn.umap"),
    DatasetSpec(label="1.2m", path=Path("1M_20260121.h5ad"), prefer_obsm="X_umap"),
]

SELECTION_SIZES: dict[str, list[int]] = {
    "pbmc3k": [500, 1000, 2000],
    "pbmc161k": [20661, 25531, 50000],
    "1.2m": [25531, 143672, 376364],
}


def _npz_str(z: dict, key: str) -> str | None:
    if key not in z:
        return None
    v = z[key]
    if isinstance(v, np.ndarray) and v.shape == ():
        return str(v.item())
    return str(v)


def load_exported_selection(path: Path) -> ExportedSelection:
    with np.load(path, allow_pickle=False) as z:
        idx = np.asarray(z["indices"], dtype=int)
        dataset = _npz_str(z, "dataset") or ""
        name = _npz_str(z, "name") or path.stem

    if not dataset:
        raise ValueError(f"{path}: missing required field 'dataset' (e.g. 'pbmc161k').")
    if idx.ndim != 1:
        raise ValueError(f"{path}: expected 1D indices array, got shape={idx.shape}.")
    return ExportedSelection(dataset=dataset, name=name, indices=idx)


def take_rows(X, row_indices: np.ndarray):
    """
    Take rows from X in the given order.

    Note: h5py.Dataset fancy indexing requires increasing indices, so we reorder
    for IO and then permute back to the requested order.
    """
    idx = np.asarray(row_indices, dtype=int)
    if idx.size == 0:
        return X[idx, :]

    if sp.issparse(X):
        return X[idx, :]

    if type(X).__module__.startswith("h5py"):
        idx_sorted = np.sort(idx)
        X_sorted = X[idx_sorted, :]
        pos = np.searchsorted(idx_sorted, idx)
        return X_sorted[pos, :]

    return X[idx, :]


def _pick_embedding_key(adata: ad.AnnData, prefer: str | None) -> str:
    obsm_keys = set(getattr(adata, "obsm", {}).keys())
    default = adata.uns.get("default_embedding", None)

    if isinstance(default, str) and default in obsm_keys:
        return default
    if prefer and prefer in obsm_keys:
        return prefer
    for k in ("X_wnn.umap", "X_umap"):
        if k in obsm_keys:
            return k
    raise KeyError(f"No suitable 2D embedding found in .obsm. Available: {sorted(obsm_keys)}")


def _principal_axis(coords: np.ndarray) -> np.ndarray:
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError(f"Expected coords shape (n,2), got {coords.shape}")
    centered = coords - coords.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    v = vt[0]
    v = v / (np.linalg.norm(v) + 1e-12)
    return v


def select_cells_along_axis(
    coords: np.ndarray,
    n_select: int,
    *,
    axis_quantiles: tuple[float, float] = (0.05, 0.95),
) -> np.ndarray:
    """
    Deterministically pick `n_select` cells near a principal-axis "trajectory".

    The selection is:
      - restricted to the middle quantile range along the principal axis, and
      - chosen as the closest points to that axis (perpendicular distance),
      - ordered along the axis so the selection has a consistent direction.
    """
    coords = np.asarray(coords, dtype=float)
    n_obs = int(coords.shape[0])
    n_select = int(min(max(1, n_select), n_obs))

    v = _principal_axis(coords)
    mu = coords.mean(axis=0, keepdims=True)
    centered = coords - mu

    t = centered @ v
    q0, q1 = axis_quantiles
    t0 = float(np.quantile(t, q0))
    t1 = float(np.quantile(t, q1))
    if not np.isfinite(t0) or not np.isfinite(t1) or t1 <= t0:
        t0, t1 = float(np.min(t)), float(np.max(t))

    candidates = np.flatnonzero((t >= t0) & (t <= t1))
    if candidates.size == 0:
        candidates = np.arange(n_obs, dtype=int)

    tc = t[candidates]
    xc = centered[candidates]
    # Perpendicular distance to the axis: ||x - (x·v)v||
    proj = np.outer(tc, v)
    resid = xc - proj
    dist2 = np.einsum("ij,ij->i", resid, resid)

    if candidates.size <= n_select:
        chosen = candidates
    else:
        part = np.argpartition(dist2, n_select - 1)[:n_select]
        chosen = candidates[part]

    # Order from low -> high along the axis to induce a direction.
    order = np.argsort(t[chosen], kind="mergesort")
    return chosen[order].astype(int, copy=False)


def projections_for_selection(coords: np.ndarray, selected_indices: np.ndarray) -> np.ndarray:
    coords = np.asarray(coords, dtype=float)
    sel = np.asarray(selected_indices, dtype=int)
    pts = coords[sel]
    if pts.shape[0] < 2:
        return np.zeros((pts.shape[0],), dtype=float)

    v = pts[-1] - pts[0]
    nv = float(np.linalg.norm(v))
    if not np.isfinite(nv) or nv <= 1e-15:
        return np.zeros((pts.shape[0],), dtype=float)
    v = v / nv
    start = pts[0]
    return (pts - start) @ v


def de_global_stats(X) -> dict:
    n_obs, n_vars = map(int, X.shape)
    if type(X).__module__.startswith("anndata._core.sparse_dataset"):
        raise RuntimeError(
            "Backed sparse matrices (anndata sparse_dataset) are not supported by this benchmark "
            "in the current environment. Load the .h5ad in-memory (do not use backed='r')."
        )
    if sp.issparse(X):
        if sp.isspmatrix_csr(X):
            total_sum, total_sqsum = diffexpr_sum_sqsum_csr_global(X, backend="auto")
        else:
            total_sum = np.asarray(X.sum(axis=0)).ravel()
            total_sqsum = np.asarray(X.power(2).sum(axis=0)).ravel()
    elif type(X).__module__.startswith("h5py"):
        # Chunked reduction to avoid materializing large dense arrays and to
        # support h5py.Dataset, which does not behave like a NumPy array in all ops.
        total_sum = np.zeros((n_vars,), dtype=float)
        total_sqsum = np.zeros((n_vars,), dtype=float)
        chunk = 8192
        for start in range(0, n_obs, chunk):
            stop = min(n_obs, start + chunk)
            block = np.asarray(X[start:stop, :], dtype=float)
            total_sum += block.sum(axis=0)
            total_sqsum += np.square(block).sum(axis=0)
    else:
        arr = np.asarray(X, dtype=float)
        if arr.ndim != 2:
            raise ValueError(f"Expected X to be 2D, got shape={arr.shape}")
        total_sum = np.asarray(arr.sum(axis=0)).ravel()
        total_sqsum = np.asarray(np.square(arr).sum(axis=0)).ravel()

    return {
        "n_obs": n_obs,
        "n_vars": n_vars,
        "sum": total_sum.astype(float, copy=False),
        "sqsum": total_sqsum.astype(float, copy=False),
    }


def de_compute_from_stats(
    X,
    *,
    selected_indices: np.ndarray,
    total_sum: np.ndarray,
    total_sqsum: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    sel = np.unique(np.asarray(selected_indices, dtype=int))
    sel = sel[(sel >= 0) & (sel < int(X.shape[0]))]

    n1 = int(sel.size)
    n = int(X.shape[0])
    n2 = n - n1
    if n1 < 2 or n2 < 2:
        return np.zeros((int(X.shape[1]),), dtype=float), np.ones((int(X.shape[1]),), dtype=float)

    if sp.isspmatrix_csr(X):
        sum1, sqsum1 = diffexpr_sum_sqsum_selected_csr(X, sel, backend="auto")
        sum1 = np.asarray(sum1, dtype=float).ravel()
        sqsum1 = np.asarray(sqsum1, dtype=float).ravel()
    else:
        X1 = take_rows(X, sel)
        if hasattr(X1, "sum"):
            sum1 = np.asarray(X1.sum(axis=0)).ravel().astype(float, copy=False)
        else:
            sum1 = np.asarray(np.sum(X1, axis=0)).ravel().astype(float, copy=False)

        if hasattr(X1, "power"):
            sqsum1 = np.asarray(X1.power(2).sum(axis=0)).ravel().astype(float, copy=False)
        else:
            sqsum1 = np.asarray(np.square(np.asarray(X1, dtype=float)).sum(axis=0)).ravel().astype(float, copy=False)

    sum2 = (total_sum - sum1).astype(float, copy=False)
    sqsum2 = (total_sqsum - sqsum1).astype(float, copy=False)

    mean1 = sum1 / n1
    mean2 = sum2 / n2

    var1 = (sqsum1 - (sum1 * sum1) / n1) / (n1 - 1)
    var2 = (sqsum2 - (sum2 * sum2) / n2) / (n2 - 1)
    var1 = np.maximum(var1, 0.0)
    var2 = np.maximum(var2, 0.0)

    std1 = np.sqrt(var1)
    std2 = np.sqrt(var2)

    t_stat, p_val = ss.ttest_ind_from_stats(
        mean1=mean1,
        std1=std1,
        nobs1=n1,
        mean2=mean2,
        std2=std2,
        nobs2=n2,
        equal_var=False,
    )

    t_stat = np.nan_to_num(np.asarray(t_stat, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    p_val = np.nan_to_num(np.asarray(p_val, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
    return t_stat, p_val


def _nnz_of(X) -> int | None:
    try:
        if sp.issparse(X):
            return int(X.nnz)
        nnz = getattr(X, "nnz", None)
        if nnz is not None:
            return int(nnz)
        return None
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser(prog="perf_retest")
    parser.add_argument(
        "--out",
        default="benchmarks/perf_results.csv",
        help="Output CSV path (default: benchmarks/perf_results.csv)",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=[d.label for d in DATASETS],
        help=f"Datasets to run (default: {', '.join(d.label for d in DATASETS)})",
    )
    parser.add_argument(
        "--selection-npz",
        action="append",
        default=[],
        help=(
            "Optional: path to a .npz exported from a notebook containing "
            "`dataset` (label) and `indices` (1D int array). Can be passed multiple times."
        ),
    )
    parser.add_argument(
        "--no-synthetic",
        action="store_true",
        help="If set, do not run the built-in synthetic selections (only exported selections).",
    )
    args = parser.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    selected_specs = [d for d in DATASETS if d.label in set(args.datasets)]
    if not selected_specs:
        raise SystemExit("No datasets selected.")

    exported: list[ExportedSelection] = []
    for p in args.selection_npz:
        exported.append(load_exported_selection(Path(p)))

    rows: list[dict] = []
    for spec in selected_specs:
        if not spec.path.exists():
            raise FileNotFoundError(
                f"Missing dataset file: {spec.path}. Place it in the repo root before benchmarking."
            )

        adata = ad.read_h5ad(spec.path)
        embed_key = _pick_embedding_key(adata, spec.prefer_obsm)
        coords = np.asarray(adata.obsm[embed_key])
        if coords.shape[1] != 2:
            raise ValueError(f"{spec.label}: embedding {embed_key!r} is not 2D (shape={coords.shape}).")

        X = adata.X
        n_obs, n_vars = map(int, adata.shape)
        x_type = type(X).__name__

        # One-time DE global stats per dataset (mirrors scSketch caching behavior).
        t0 = time.perf_counter()
        de_stats = de_global_stats(X)
        de_global_s = time.perf_counter() - t0

        runs: list[tuple[str, str, np.ndarray]] = []
        if not bool(args.no_synthetic):
            for n_select in SELECTION_SIZES[spec.label]:
                sel = select_cells_along_axis(coords, n_select)
                runs.append(("synthetic", f"synthetic_n{int(sel.size)}", sel))

        for s in exported:
            if s.dataset == spec.label:
                runs.append(("exported", s.name, s.indices))

        if not runs:
            continue

        for source, name, sel in runs:
            sel = np.asarray(sel, dtype=int)
            proj = projections_for_selection(coords, sel)
            X_sel = take_rows(X, sel)

            # Directional: correlation + p-values.
            t0 = time.perf_counter()
            td = test_direction(X_sel, proj)
            dir_corr_s = time.perf_counter() - t0

            # LORD update (separate timing).
            t0 = time.perf_counter()
            _ = lord_test(np.asarray(td["p_value"], dtype=float), initial_results=None, alpha=0.05, w0=0.005)
            lord_s = time.perf_counter() - t0

            # Differential expression for this selection using precomputed global sums.
            t0 = time.perf_counter()
            t_stat, p_val = de_compute_from_stats(
                X,
                selected_indices=sel,
                total_sum=de_stats["sum"],
                total_sqsum=de_stats["sqsum"],
            )
            de_s = time.perf_counter() - t0

            rows.append(
                {
                    "dataset": spec.label,
                    "file": str(spec.path),
                    "embedding": embed_key,
                    "selection_source": source,
                    "selection_name": name,
                    "n_obs": n_obs,
                    "n_vars": n_vars,
                    "X_type": x_type,
                    "X_nnz_total": _nnz_of(X),
                    "n_selected": int(sel.size),
                    "X_sel_type": type(X_sel).__name__,
                    "X_sel_nnz": _nnz_of(X_sel),
                    "directional_corr_p_s": dir_corr_s,
                    "directional_lord_s": lord_s,
                    "de_global_stats_s": de_global_s,
                    "de_compute_s": de_s,
                }
            )

    fieldnames = list(rows[0].keys()) if rows else []
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote {len(rows)} rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
