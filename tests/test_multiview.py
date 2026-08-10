import numpy as np
import pandas as pd
from anndata import AnnData

from jscatter import Scatter
from scsketch._action_log import apply_action_log_entry
from scsketch.scsketch import ScSketch


def _adata():
    adata = AnnData(
        X=np.array(
            [
                [0.0, 4.0],
                [1.0, 3.0],
                [2.0, 2.0],
                [3.0, 1.0],
            ]
        ),
        obs={"cluster": ["A", "B", "A", "C"]},
    )
    adata.obs_names = [f"cell_{i}" for i in range(4)]
    adata.var_names = ["GeneA", "GeneB"]
    adata.obsm["X_umap"] = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ]
    )
    return adata


def _extra_view(adata):
    data = pd.DataFrame(
        {
            "PC1": [0.0, 1.0, 2.0, 3.0],
            "PC2": [3.0, 2.0, 1.0, 0.0],
            "cluster": adata.obs["cluster"].astype(str).to_numpy(),
        },
        index=adata.obs_names,
    )
    return Scatter(data=data, x="PC1", y="PC2", color_by="cluster")


def test_extra_views_enable_multi_view_panel_and_share_color_map():
    adata = _adata()
    pca = _extra_view(adata)

    sketch = ScSketch(
        adata=adata,
        metadata_cols=["cluster"],
        color_by_default="cluster",
        extra_views={"PCA": pca},
    )

    assert sketch._ctrl.multi_view_toggle.value is True
    assert sketch._ctrl.multi_view_container.layout.display == "flex"
    assert sketch._ctrl.pathway_table_container.layout.display == "none"
    assert pca.width() == 420
    assert pca.height() == 420
    assert pca._color_map_order == list(sketch.categorical_color_maps["cluster"])


def test_extra_view_selection_syncs_with_main_scatter():
    adata = _adata()
    pca = _extra_view(adata)
    sketch = ScSketch(
        adata=adata,
        metadata_cols=["cluster"],
        color_by_default="cluster",
        extra_views={"PCA": pca},
    )

    sketch._main_selection_multi_view_handler({"new": [0, 2]})
    assert np.asarray(pca.selection()).tolist() == [0, 2]

    sketch._extra_view_selection_handler({"new": [1, 3]})
    assert np.asarray(sketch.scatter.selection()).tolist() == [1, 3]


def test_gene_expression_coloring_syncs_to_extra_views():
    adata = _adata()
    pca = _extra_view(adata)
    sketch = ScSketch(
        adata=adata,
        metadata_cols=["cluster"],
        color_by_default="cluster",
        extra_views={"PCA": pca},
    )

    sketch._color_embedding_by_gene("GeneA")

    assert pca._color_by == "Custom Color Data"
    assert pca._color_labeling["variable"] == "GeneA"
    np.testing.assert_allclose(pca._color_data.to_numpy(), [0.0, 1.0, 2.0, 3.0])


def test_multi_view_toggle_action_restores_panel_state():
    adata = _adata()
    pca = _extra_view(adata)
    sketch = ScSketch(
        adata=adata,
        metadata_cols=["cluster"],
        color_by_default="cluster",
        extra_views={"PCA": pca},
    )

    message = apply_action_log_entry(
        sketch,
        {"type": "toggle_multi_view", "payload": {"multi_view": False}},
    )

    assert message == "Restored multi-view to False."
    assert sketch._ctrl.multi_view_toggle.value is False
    assert sketch._ctrl.multi_view_container.layout.display == "none"


def test_session_round_trip_restores_multi_view_toggle_state():
    adata = _adata()
    sketch = ScSketch(
        adata=adata,
        metadata_cols=["cluster"],
        color_by_default="cluster",
        extra_views={"PCA": _extra_view(adata)},
    )
    sketch._ctrl.multi_view_toggle.value = False

    session = sketch.export_session()
    target = ScSketch(
        adata=adata,
        metadata_cols=["cluster"],
        color_by_default="cluster",
        extra_views={"PCA": _extra_view(adata)},
    )
    target.load_session(session)

    assert session["state"]["multi_view"] is False
    assert target._ctrl.multi_view_toggle.value is False
    assert target._ctrl.multi_view_container.layout.display == "none"
