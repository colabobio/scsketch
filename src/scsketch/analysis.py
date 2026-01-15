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
