import numpy as np

from scsketch.analysis import lord_test


def _pvals_with_forced_rejections(n: int) -> np.ndarray:
    p = np.ones((n,), dtype=float)
    # Force rejections at stable, non-borderline positions.
    for idx in (0, 10, 50, 123, n // 2, n - 2, n - 1):
        if 0 <= idx < n:
            p[idx] = 0.0
    return p


def test_lord_numba_matches_python_full_run():
    p = _pvals_with_forced_rejections(2000)
    res_py = lord_test(p, backend="python")
    res_nb = lord_test(p, backend="numba")

    assert np.allclose(res_py["alpha_i"], res_nb["alpha_i"], rtol=0.0, atol=1e-12)
    assert np.array_equal(res_py["R"], res_nb["R"])
    assert res_py["tau"] == res_nb["tau"]


def test_lord_numba_incremental_matches_full_run():
    p = _pvals_with_forced_rejections(2000)
    p1, p2 = p[:1000], p[1000:]

    res1 = lord_test(p1, backend="numba")
    res_inc = lord_test(np.concatenate([p1, p2]), initial_results=res1, backend="numba")
    res_full = lord_test(np.concatenate([p1, p2]), backend="numba")

    assert np.allclose(res_inc["alpha_i"], res_full["alpha_i"], rtol=0.0, atol=1e-12)
    assert np.array_equal(res_inc["R"], res_full["R"])
    assert res_inc["tau"] == res_full["tau"]

