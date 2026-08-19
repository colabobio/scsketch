"""Ranking score helpers for displaying exploratory results."""

from __future__ import annotations

import numpy as np


def discovery_scores(p_values) -> np.ndarray:
    """Return integer 0-10 Discovery Scores derived from nominal p-values."""
    p = np.asarray(p_values, dtype=float)
    if p.size == 0:
        return np.asarray([], dtype=float)

    p = np.clip(p, np.finfo(float).tiny, 1.0)
    signal = -np.log10(p)
    lo = float(np.min(signal))
    hi = float(np.max(signal))

    if hi <= lo:
        return np.full(signal.shape, 10, dtype=int)

    return np.rint(10.0 * (signal - lo) / (hi - lo)).astype(int)
