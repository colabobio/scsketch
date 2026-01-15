"""Directional analysis functions for scSketch."""

import numpy as np
import scipy.stats as ss


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
    rs = np.corrcoef(projection, X, rowvar=False)[0, 1:]

    n = len(X)
    T = -np.abs(rs * np.sqrt(n - 2)) / np.sqrt(1 - (rs**2))
    return {"correlation": rs, "p_value": ss.t.cdf(T, df=n - 2) * 2}

def lord_test(pval, initial_results=None, gammai=None, alpha=0.05, w0=0.005):
    """"
    This is a translation of "version 1" under:

    https://github.com/bioc/onlineFDR/blob/devel/src/lord.cpp
    
    The only changes are that we don't recompute threhsolds for hypotheses that
    we have already seen. This only necessary because we may continue testing
    for many directions.
    """
    N = len(pval)

    if gammai is None:
        gammai = (
            0.07720838
            * np.log(np.maximum(np.arange(1, N + 2), 2))
            / (np.arange(1, N + 2) * np.exp(np.sqrt(np.log(np.arange(1, N + 2)))))
        )

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

