"""Differential expression computation engine.

This module is UI-free.  ``DiffExprEngine`` wraps all the stateful
global-stats caching and Welch t-test logic that was previously
embedded in ``ScSketch``.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.stats as ss
from anndata import AnnData

from .analysis import (
    diffexpr_sum_sqsum_selected_csr,
    diffexpr_sum_sqsum_csr_global,
)

logger = logging.getLogger(__name__)


class DiffExprEngine:
    """Stateful helper for computing differential expression.

    Holds the per-dataset global sum / sum-of-squares cache so that
    repeated DE queries do not recompute the full-matrix pass.

    Parameters
    ----------
    adata:
        AnnData object whose expression matrix will be queried.
    t_threshold:
        |T| threshold for reporting a gene as differentially expressed.
    p_threshold:
        p-value threshold for reporting a gene.
    disk_cache_dir:
        Optional directory to persist the global stats to disk.
    """

    def __init__(
        self,
        adata: AnnData,
        t_threshold: float = 2.0,
        p_threshold: float = 0.05,
        disk_cache_dir: Optional[Path] = None,
    ) -> None:
        self.adata = adata
        self.t_threshold = t_threshold
        self.p_threshold = p_threshold
        self.disk_cache_dir = disk_cache_dir
        self._global_stats: dict | None = None

    # ── Data source ──────────────────────────────────────────────────────────

    def de_source(self) -> tuple:
        """Return ``(X, var_names, source_tag)`` preferring ``adata.raw``."""
        adata = self.adata
        if getattr(adata, "raw", None) is not None and getattr(adata.raw, "X", None) is not None:
            return adata.raw.X, list(adata.raw.var_names), "raw"
        return adata.X, list(adata.var_names), "X"

    # ── Disk cache helpers ───────────────────────────────────────────────────

    def _cache_path(self, *, source: str, n_obs: int, n_vars: int, var_names: list[str]) -> Path | None:
        if self.disk_cache_dir is None:
            return None
        h = hashlib.sha256()
        h.update(source.encode("utf-8"))
        h.update(b"\0")
        h.update(str(n_obs).encode("utf-8"))
        h.update(b"\0")
        h.update(str(n_vars).encode("utf-8"))
        h.update(b"\0")
        for name in var_names:
            h.update(str(name).encode("utf-8"))
            h.update(b"\0")
        key = h.hexdigest()[:16]
        return self.disk_cache_dir / f"diffexpr_global_stats_{key}.npz"

    # ── Global stats (lazy, cached) ──────────────────────────────────────────

    def ensure_global_stats(self) -> None:
        """Compute (or load from disk) per-gene sum and sum-of-squares across all cells."""
        X, var_names, source = self.de_source()

        # Already up-to-date?
        gs = self._global_stats
        if (
            gs is not None
            and gs.get("source") == source
            and gs.get("n_obs") == int(X.shape[0])
            and gs.get("n_vars") == int(X.shape[1])
        ):
            return

        cache_path = self._cache_path(
            source=source,
            n_obs=int(X.shape[0]),
            n_vars=int(X.shape[1]),
            var_names=var_names,
        )

        # Try loading from disk cache
        if cache_path is not None and cache_path.exists():
            try:
                with np.load(cache_path, allow_pickle=False) as z:
                    cached_sum = np.asarray(z["sum"], dtype=float).ravel()
                    cached_sqsum = np.asarray(z["sqsum"], dtype=float).ravel()
                    cached_n_obs = int(z["n_obs"])
                    cached_n_vars = int(z["n_vars"])
                    cached_source = str(z["source"])
                if (
                    cached_source == source
                    and cached_n_obs == int(X.shape[0])
                    and cached_n_vars == int(X.shape[1])
                    and cached_sum.shape[0] == cached_n_vars
                    and cached_sqsum.shape[0] == cached_n_vars
                ):
                    self._global_stats = {
                        "source": source,
                        "n_obs": int(X.shape[0]),
                        "n_vars": int(X.shape[1]),
                        "var_names": var_names,
                        "sum": cached_sum,
                        "sqsum": cached_sqsum,
                        "disk_cache": str(cache_path),
                    }
                    return
            except Exception:
                logger.exception("Failed to load DE disk cache; recomputing global stats.")

        # Compute from scratch
        if sp.isspmatrix_csr(X):
            total_sum, total_sqsum = diffexpr_sum_sqsum_csr_global(X, backend="auto")
        else:
            total_sum = np.asarray(X.sum(axis=0) if hasattr(X, "sum") else np.sum(X, axis=0)).ravel().astype(float, copy=False)
            total_sqsum = (
                np.asarray(X.power(2).sum(axis=0) if hasattr(X, "power") else np.einsum("ij,ij->j", X, X))
                .ravel()
                .astype(float, copy=False)
            )

        self._global_stats = {
            "source": source,
            "n_obs": int(X.shape[0]),
            "n_vars": int(X.shape[1]),
            "var_names": var_names,
            "sum": total_sum,
            "sqsum": total_sqsum,
        }

        # Persist to disk
        if cache_path is not None:
            try:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                tmp = cache_path.with_name(cache_path.name + ".tmp")
                with tmp.open("wb") as f:
                    np.savez(
                        f,
                        source=source,
                        n_obs=int(X.shape[0]),
                        n_vars=int(X.shape[1]),
                        sum=self._global_stats["sum"],
                        sqsum=self._global_stats["sqsum"],
                    )
                tmp.replace(cache_path)
                self._global_stats["disk_cache"] = str(cache_path)
            except Exception:
                logger.exception("Failed to write DE disk cache; continuing without it.")

    # ── Main computation ─────────────────────────────────────────────────────

    def compute(
        self,
        selected_indices: np.ndarray,
        selection_label: str,
        *,
        t_threshold: float | None = None,
        p_threshold: float | None = None,
    ) -> list[dict]:
        """Run Welch's t-test between selected and background cells.

        Parameters
        ----------
        selected_indices:
            Integer cell indices forming the foreground group.
        selection_label:
            Human-readable name stored in each result dict.
        t_threshold:
            Override for ``|T|`` threshold (uses ``self.t_threshold`` if ``None``).
        p_threshold:
            Override for p-value threshold (uses ``self.p_threshold`` if ``None``).

        Returns
        -------
        List of dicts with keys ``attribute``, ``interval``, ``quality``, ``direction``.
        Sorted descending by ``|T|``, capped at 200 genes.
        """
        self.ensure_global_stats()
        stats = self._global_stats
        if stats is None:
            return []

        X, var_names, _ = self.de_source()
        total_sum: np.ndarray = stats["sum"]
        total_sqsum: np.ndarray = stats["sqsum"]

        t_thr = t_threshold if t_threshold is not None else self.t_threshold
        p_thr = p_threshold if p_threshold is not None else self.p_threshold

        selected_indices = np.asarray(selected_indices, dtype=int)
        selected_indices = selected_indices[
            (selected_indices >= 0) & (selected_indices < int(X.shape[0]))
        ]
        selected_indices = np.unique(selected_indices)

        n1 = int(selected_indices.shape[0])
        n = int(X.shape[0])
        n2 = n - n1
        if n1 < 2 or n2 < 2:
            return []

        # Per-group sums
        if sp.isspmatrix_csr(X):
            sum1, sqsum1 = diffexpr_sum_sqsum_selected_csr(X, selected_indices, backend="auto")
            sum1 = np.asarray(sum1, dtype=float).ravel()
            sqsum1 = np.asarray(sqsum1, dtype=float).ravel()
        else:
            X1 = X[selected_indices, :]
            sum1 = np.asarray(X1.sum(axis=0) if hasattr(X1, "sum") else np.sum(X1, axis=0)).ravel().astype(float, copy=False)
            sqsum1 = (
                np.asarray(X1.power(2).sum(axis=0) if hasattr(X1, "power") else np.einsum("ij,ij->j", X1, X1))
                .ravel()
                .astype(float, copy=False)
            )

        sum2 = (total_sum - sum1).astype(float, copy=False)
        sqsum2 = (total_sqsum - sqsum1).astype(float, copy=False)

        mean1 = sum1 / n1
        mean2 = sum2 / n2

        var1 = np.maximum((sqsum1 - (sum1 * sum1) / n1) / (n1 - 1), 0.0)
        var2 = np.maximum((sqsum2 - (sum2 * sum2) / n2) / (n2 - 1), 0.0)

        t_stat, p_val = ss.ttest_ind_from_stats(
            mean1=mean1,
            std1=np.sqrt(var1),
            nobs1=n1,
            mean2=mean2,
            std2=np.sqrt(var2),
            nobs2=n2,
            equal_var=False,
        )
        t_stat = np.asarray(t_stat, dtype=float)
        p_val = np.asarray(p_val, dtype=float)

        ok = np.isfinite(t_stat) & np.isfinite(p_val)
        keep = ok & (np.abs(t_stat) >= t_thr) & (p_val <= p_thr)
        idx = np.where(keep)[0]
        if idx.size == 0:
            return []

        idx = idx[np.argsort(np.abs(t_stat[idx]))[::-1]][:200]

        return [
            {
                "attribute": var_names[int(j)],
                "interval": (float(t_stat[int(j)]), float(p_val[int(j)])),
                "quality": float(abs(t_stat[int(j)])),
                "direction": selection_label,
            }
            for j in idx
        ]
