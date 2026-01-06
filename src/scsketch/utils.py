from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from jscatter import Line
from matplotlib.path import Path


def find_equidistant_vertices(vertices: np.ndarray, n_points: int) -> np.ndarray:
    seg = np.diff(vertices, axis=0)
    seg_len = np.linalg.norm(seg, axis=1)
    cum = np.concatenate(([0.0], np.cumsum(seg_len)))
    total = cum[-1]
    targets = np.linspace(0.0, total, n_points)
    out = np.zeros((n_points, 2))
    for i, t in enumerate(targets):
        j = max(0, min(np.searchsorted(cum, t, side="right") - 1, len(seg) - 1))
        alpha = (t - cum[j]) / (seg_len[j] if seg_len[j] else 1.0)
        out[i] = vertices[j] + alpha * seg[j]
    return out


def split_line_at_points(vertices: np.ndarray, split_points: np.ndarray) -> list[np.ndarray]:
    if len(split_points) < 2:
        return []
    return [
        np.vstack([split_points[i], split_points[i + 1]])
        for i in range(len(split_points) - 1)
    ]


def split_line_equidistant(vertices: np.ndarray, n_points: int) -> list[np.ndarray]:
    return split_line_at_points(vertices, find_equidistant_vertices(vertices, n_points))


def points_in_polygon(points: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    return Path(polygon).contains_points(points)


@dataclass
class Selection:
    """Class for keeping track of a selection."""

    index: int
    name: str
    points: np.ndarray
    color: str
    lasso: Line
    hull: Line
    path: np.ndarray | None = None


@dataclass
class Selections:
    """Class for keeping track of selections."""

    selections: list[Selection] = field(default_factory=list)

    def all_points(self) -> np.ndarray:
        return np.unique(
            np.concatenate(list(map(lambda selection: selection.points, self.selections)))
        )

    def all_hulls(self) -> list[Line]:
        return [s.hull for s in self.selections]


@dataclass
class Lasso:
    """Class for keeping track of the lasso polygon."""

    polygon: Line | None = None

