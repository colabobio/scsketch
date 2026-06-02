import numpy as np
from anndata import AnnData

from scsketch._diffexpr import DiffExprEngine


def test_diffexpr_compute_returns_all_passing_genes_without_200_cap():
    rng = np.random.default_rng(0)
    n_cells = 40
    n_genes = 260
    selected = np.arange(20)

    X = rng.normal(loc=0.0, scale=0.4, size=(n_cells, n_genes))
    X[selected, :] += 2.0
    adata = AnnData(
        X=X,
        var={"gene": [f"GENE{i}" for i in range(n_genes)]},
    )
    adata.var_names = adata.var["gene"]

    engine = DiffExprEngine(adata, t_threshold=0.0, p_threshold=1.0)
    results = engine.compute(selected, "Selection 1")

    assert len(results) == n_genes
    assert results == sorted(
        results,
        key=lambda entry: abs(entry["interval"][0]),
        reverse=True,
    )
