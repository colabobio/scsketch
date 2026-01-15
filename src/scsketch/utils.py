"""Utility functions and data structures for scSketch."""

import numpy as np
from dataclasses import dataclass, field
from jscatter import Line
from scipy.spatial import ConvexHull


@dataclass
class Selection:
    """Class for keeping track of a selection."""

    index: int
    name: str
    points: np.ndarray
    color: str
    lasso: Line
    hull: Line


@dataclass
class Selections:
    """Class for keeping track of selections."""

    selections: list[Selection] = field(default_factory=list)

    def all_points(self) -> np.ndarray:
        """Get all points from all selections."""
        return np.unique(
            np.concatenate(
                list(map(lambda selection: selection.points, self.selections))
            )
        )

    def all_hulls(self) -> list[Line]:
        """Get all convex hulls from all selections."""
        return [s.hull for s in self.selections]


@dataclass
class Lasso:
    """Class for keeping track of the lasso polygon."""

    polygon: Line | None = None


def create_selection(
    index: int,
    name: str,
    points_indices: np.ndarray,
    points_coords: np.ndarray,
    lasso_polygon: np.ndarray,
    color: str,
) -> Selection:
    """
    Create a Selection object from raw data.

    Args:
        index: Selection index
        name: Selection name
        points_indices: Indices of selected points in the dataframe
        points_coords: XY coordinates of selected points
        lasso_polygon: Lasso polygon coordinates
        color: Selection color

    Returns:
        Selection object
    """
    hull = ConvexHull(points_coords)
    hull_points = np.vstack((points_coords[hull.vertices], points_coords[hull.vertices[0]]))

    lasso_polygon_list = lasso_polygon.tolist()
    lasso_polygon_list.append(lasso_polygon_list[0])

    return Selection(
        index=index,
        name=name,
        points=points_indices,
        color=color,
        lasso=Line(lasso_polygon_list),
        hull=Line(hull_points, line_color=color, line_width=2),
    )


def fetch_pathways(gene: str) -> list[dict]:
    """
    Fetch Reactome pathways for a given gene symbol.

    Args:
        gene: Gene symbol

    Returns:
        List of pathway dictionaries with 'Pathway' and 'stId' keys
    """
    import requests

    url = f"https://reactome.org/ContentService/data/mapping/UniProt/{gene}/pathways?species=9606"
    try:
        response = requests.get(url)
        response.raise_for_status()
        pathways = response.json()
        return [
            {"Pathway": entry["displayName"], "stId": entry["stId"]}
            for entry in pathways
        ]
    except requests.exceptions.RequestException as e:
        print(f"Error fetching Reactome pathways for {gene}: {e}")
        return []
