"""Directional analysis functions for scSketch."""

import numpy as np
import scipy.sparse as sp
import scipy.stats as ss

try:
    import numba as nb
except Exception:  # pragma: no cover
    nb = None


def _default_gammai(N: int) -> np.ndarray:
    return (
        0.07720838
        * np.log(np.maximum(np.arange(1, N + 2), 2))
        / (np.arange(1, N + 2) * np.exp(np.sqrt(np.log(np.arange(1, N + 2)))))
    )


def compute_directional_analysis(df, selections, metadata_columns=None):
    """
    Computes the correlation of gene expression along a directional axis.

    Args:
        df (pd.DataFrame): Dataframe containing gene expression data and spatial coordinates.
        selections (Selections): The selected points for directional analysis.
        metadata_columns (list, optional): List of metadata columns to exclude from analysis.

    Returns:
        list: A list of dictionaries containing the computed correlations.
    """
    if metadata_columns is None:
        metadata_columns = ["x", "y"]

    if len(selections.selections) == 0:
        return []

    results = []

    for selection in selections.selections:
        selected_indices = selection.points
        selected_embeddings = df.iloc[selected_indices][["x", "y"]].values

        # Ensure we have at least two points for a valid direction vector
        if selected_embeddings.shape[0] < 2:
            continue

        # Compute direction vector
        v = selected_embeddings[-1] - selected_embeddings[0]
        v = v / np.linalg.norm(v)  # Normalize

        # Compute projections
        start_point = selected_embeddings[0]
        projections = np.array(
            [np.dot(pt - start_point, v) for pt in selected_embeddings]
        )

        # Get gene expression data
        columns_to_drop = [col for col in metadata_columns if col in df.columns]
        selected_expression = df.iloc[selected_indices].drop(columns=columns_to_drop)

        # Compute correlations
        correlations = []
        for gene in selected_expression.columns:
            r, p = ss.pearsonr(projections, selected_expression[gene])
            correlations.append(
                {
                    "attribute": gene,
                    "interval": (r, p),
                    "quality": abs(r),  # Use absolute correlation as a measure of quality
                }
            )

        results.append(correlations)

    return results

def test_direction(X, projection):
    """Compute per-gene correlation/p-values against a 1D projection.

    This is implemented to avoid densifying `X` when it is sparse (common for scRNA-seq).
    """
    projection = np.asarray(projection, dtype=float).ravel()
    n = int(projection.shape[0])

    n_genes = int(getattr(X, "shape", (n, 0))[1])
    if n < 3 or n_genes == 0:
        return {
            "correlation": np.zeros(n_genes, dtype=float),
            "p_value": np.ones(n_genes, dtype=float),
        }

    pc = projection - float(np.mean(projection))
    ss_pc = float(np.dot(pc, pc))
    if not np.isfinite(ss_pc) or ss_pc <= 0.0:
        return {
            "correlation": np.zeros(n_genes, dtype=float),
            "p_value": np.ones(n_genes, dtype=float),
        }

    std_pc = float(np.sqrt(ss_pc / (n - 1)))

    if sp.issparse(X):
        if int(X.shape[0]) != n:
            raise ValueError("X and projection must have the same number of rows.")

        sum_x = np.asarray(X.sum(axis=0)).ravel().astype(float, copy=False)
        sum_x2 = np.asarray(X.power(2).sum(axis=0)).ravel().astype(float, copy=False)
        var_x = (sum_x2 - (sum_x * sum_x) / n) / (n - 1)
        var_x = np.maximum(var_x, 0.0)
        std_x = np.sqrt(var_x)

        cov = (np.asarray(X.T @ pc).ravel().astype(float, copy=False)) / (n - 1)
    else:
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        if int(X.shape[0]) != n:
            raise ValueError("X and projection must have the same number of rows.")

        std_x = np.std(X, axis=0, ddof=1)
        cov = (pc @ X) / (n - 1)

    denom = std_pc * std_x
    with np.errstate(divide="ignore", invalid="ignore"):
        rs = cov / denom

    rs = np.nan_to_num(rs, nan=0.0, posinf=0.0, neginf=0.0)
    rs = np.clip(rs, -1.0, 1.0)

    df = n - 2
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        T = -np.abs(rs * np.sqrt(df)) / np.sqrt(np.maximum(1.0 - (rs**2), 1e-300))
        p = ss.t.cdf(T, df=df) * 2

    p = np.nan_to_num(p, nan=1.0, posinf=1.0, neginf=1.0)
    p = np.clip(p, 0.0, 1.0)
    return {"correlation": rs, "p_value": p}

def _lord_test_python(pval, initial_results=None, gammai=None, alpha=0.05, w0=0.005):
    """"
    This is a translation of "version 1" under:

    https://github.com/bioc/onlineFDR/blob/devel/src/lord.cpp
    
    The only changes are that we don't recompute threhsolds for hypotheses that
    we have already seen. This only necessary because we may continue testing
    for many directions.
    """
    N = len(pval)

    if gammai is None:
        gammai = _default_gammai(N)

    # setup variables, substituting previous results if needed
    alphai = np.zeros(N)
    R = np.zeros(N, dtype=bool)
    tau = []
    if initial_results is not None:
        N0 = len(initial_results["p_value"])
        alphai[range(N0)] = initial_results["alpha_i"]
        R[range(N0)] = initial_results["R"]
        tau = initial_results["tau"]
    else:
        N0 = 1
        alphai[0] = gammai[0] * w0
        R[0] = pval[0] <= alphai[0]
        if R[0]:
            tau.append(0)

    # compute lord thresholds iteratively
    K = int(np.sum(R))
    for i in range(N0, N):
        if K <= 1:
            if R[i - 1]:
                tau = [i - 1]
            Cjsum = sum(gammai[i - tau[j] - 1] for j in range(K))
            alphai[i] = w0 * gammai[i] + (alpha - w0) * Cjsum
        else:
            if R[i - 1]:
                tau.append(i - 1)
            tau2 = tau[1:]
            Cjsum = sum(gammai[i - tau2[j] - 1] for j in range(K - 1))
            alphai[i] = (
                w0 * gammai[i] + (alpha - w0) * gammai[i - tau[0] - 1] + alpha * Cjsum
            )

        if pval[i] <= alphai[i]:
            R[i] = True
            K += 1

    return {"p_value": pval, "alpha_i": alphai, "R": R, "tau": tau}


if nb is not None:

    @nb.njit(cache=True)
    def _lord_update_numba(pval, alphai, R, gammai, alpha, w0, tau, tau_len, start_i, K):
        N = pval.shape[0]
        for i in range(start_i, N):
            if K <= 1:
                if R[i - 1]:
                    tau[0] = i - 1
                    tau_len = 1
                Cjsum = 0.0
                if K == 1:
                    Cjsum = gammai[i - tau[0] - 1]
                alphai[i] = w0 * gammai[i] + (alpha - w0) * Cjsum
            else:
                if R[i - 1]:
                    # Tau can lag K by 1 if the last rejection was at i-1 (end-of-call behavior).
                    if tau_len < K:
                        tau[tau_len] = i - 1
                        tau_len += 1
                Cjsum = 0.0
                for j in range(1, K):
                    Cjsum += gammai[i - tau[j] - 1]
                alphai[i] = (
                    w0 * gammai[i]
                    + (alpha - w0) * gammai[i - tau[0] - 1]
                    + alpha * Cjsum
                )

            if pval[i] <= alphai[i]:
                R[i] = True
                K += 1
        return tau_len, K


def lord_test(pval, initial_results=None, gammai=None, alpha=0.05, w0=0.005, *, backend="auto"):
    """
    LORD++ online FDR control.

    Parameters
    ----------
    pval:
        1D array of p-values.
    initial_results:
        Previous output from `lord_test` (for incremental updates).
    gammai:
        LORD++ gamma sequence. If None, a default sequence is used.
    backend:
        "auto" (default), "python", or "numba".
    """
    if backend not in {"auto", "python", "numba"}:
        raise ValueError("backend must be one of: 'auto', 'python', 'numba'.")

    pval = np.asarray(pval, dtype=float).ravel()
    N = int(pval.shape[0])
    if gammai is None:
        gammai = _default_gammai(N)
    else:
        gammai = np.asarray(gammai, dtype=float).ravel()

    want_numba = backend == "numba" or (backend == "auto" and nb is not None and N >= 5000)
    if not want_numba or nb is None:
        return _lord_test_python(
            pval, initial_results=initial_results, gammai=gammai, alpha=alpha, w0=w0
        )

    alphai = np.zeros(N, dtype=float)
    R = np.zeros(N, dtype=np.bool_)
    tau_list: list[int] = []
    if initial_results is not None:
        N0 = int(len(initial_results["p_value"]))
        alphai[:N0] = np.asarray(initial_results["alpha_i"], dtype=float)
        R[:N0] = np.asarray(initial_results["R"], dtype=np.bool_)
        tau_list = list(initial_results.get("tau", []))
    else:
        N0 = 1
        alphai[0] = float(gammai[0]) * float(w0)
        R[0] = bool(pval[0] <= alphai[0])
        if R[0]:
            tau_list = [0]

    K = int(np.sum(R))
    tau = np.empty(N, dtype=np.int64)
    tau_len = int(len(tau_list))
    if tau_len:
        tau[:tau_len] = np.asarray(tau_list, dtype=np.int64)

    tau_len, _ = _lord_update_numba(
        pval,
        alphai,
        R,
        gammai,
        float(alpha),
        float(w0),
        tau,
        tau_len,
        int(N0),
        int(K),
    )
    return {"p_value": pval, "alpha_i": alphai, "R": R, "tau": tau[:tau_len].tolist()}
