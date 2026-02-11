import numpy as np
import scipy.sparse as sp

from scsketch.analysis import diffexpr_sum_sqsum_selected_csr


def _baseline_sum_sqsum(X: sp.csr_matrix, sel: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    X1 = X[sel, :]
    sum1 = np.asarray(X1.sum(axis=0)).ravel().astype(np.float64, copy=False)
    sqsum1 = np.asarray(X1.power(2).sum(axis=0)).ravel().astype(np.float64, copy=False)
    return sum1, sqsum1


def test_diffexpr_sum_sqsum_selected_csr_matches_baseline_python():
    rng = np.random.default_rng(0)
    X = sp.random(
        2000,
        300,
        density=0.05,
        format="csr",
        data_rvs=lambda n: rng.normal(size=n).astype(np.float32),
        random_state=0,
    )

    sel = rng.choice(X.shape[0], size=500, replace=False)
    sel = np.concatenate([sel, sel[:10], np.array([-1, -2, X.shape[0] + 10])])

    sum_b, sq_b = _baseline_sum_sqsum(X, np.unique(sel[(sel >= 0) & (sel < X.shape[0])]))
    sum_p, sq_p = diffexpr_sum_sqsum_selected_csr(X, sel, backend="python")

    assert np.allclose(sum_b, sum_p, rtol=1e-6, atol=1e-6)
    assert np.allclose(sq_b, sq_p, rtol=1e-6, atol=1e-6)


def test_diffexpr_sum_sqsum_selected_csr_matches_baseline_numba_when_available():
    try:
        import numba  # noqa: F401
    except Exception:
        return

    rng = np.random.default_rng(1)
    X = sp.random(
        4000,
        500,
        density=0.03,
        format="csr",
        data_rvs=lambda n: rng.normal(size=n).astype(np.float32),
        random_state=1,
    )

    sel = rng.choice(X.shape[0], size=1200, replace=False)

    sum_b, sq_b = _baseline_sum_sqsum(X, sel)
    sum_n, sq_n = diffexpr_sum_sqsum_selected_csr(X, sel, backend="numba")

    assert np.allclose(sum_b, sum_n, rtol=1e-6, atol=1e-6)
    assert np.allclose(sq_b, sq_n, rtol=1e-6, atol=1e-6)

