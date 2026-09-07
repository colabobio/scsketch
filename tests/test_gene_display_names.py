import ipywidgets as ipyw
import numpy as np
import pandas as pd
from anndata import AnnData

from jscatter import Scatter
from scsketch._results import show_directional_results
from scsketch.scsketch import ScSketch


def test_gene_display_names_prefer_gene_short_name():
    adata = AnnData(
        X=np.zeros((2, 2)),
        var=pd.DataFrame(
            {"gene_short_name": ["nduo-6", ""]},
            index=["WBGene00010957", "WBGene00000000"],
        ),
    )

    labels = ScSketch._build_gene_display_names(adata)

    assert labels["WBGene00010957"] == "nduo-6"
    assert "WBGene00000000" not in labels


def test_directional_result_table_displays_label_and_keeps_gene_id():
    adata = AnnData(X=np.zeros((2, 1)))
    adata.var_names = ["WBGene00010957"]
    df = pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]})
    results_box = ipyw.VBox()

    show_directional_results(
        [
            [
                {
                    "attribute": "WBGene00010957",
                    "interval": (0.8, 0.001),
                    "reject": True,
                    "direction": "Selection 1",
                }
            ]
        ],
        selections_predicates=results_box,
        pathway_table_container=ipyw.VBox(),
        reactome_diagram_container=ipyw.VBox(),
        df=df,
        adata=adata,
        active_selection=None,
        on_gene_selected=None,
        on_results_cleared=lambda *_: None,
        log=lambda *_: None,
        gene_display_name=lambda gene: "nduo-6",
    )

    table = results_box.children[0]
    assert table.data[0]["Gene"] == "nduo-6"
    assert table.data[0]["_gene_id"] == "WBGene00010957"


def test_directional_gene_click_uses_configured_species(monkeypatch):
    adata = AnnData(X=np.zeros((2, 1)))
    adata.var_names = ["unc-13"]
    df = pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]})
    results_box = ipyw.VBox()
    calls = []

    def fake_fetch_gene_description(gene, *, species="human"):
        calls.append(("annotation", gene, species))
        return {"symbol": gene, "name": "fake gene"}

    def fake_fetch_pathways(gene, *, species="human"):
        calls.append(("reactome", gene, species))
        return []

    monkeypatch.setattr(
        "scsketch._results.fetch_gene_description",
        fake_fetch_gene_description,
    )
    monkeypatch.setattr("scsketch._results.fetch_pathways", fake_fetch_pathways)

    show_directional_results(
        [
            [
                {
                    "attribute": "unc-13",
                    "interval": (0.8, 0.001),
                    "reject": True,
                    "direction": "Selection 1",
                }
            ]
        ],
        selections_predicates=results_box,
        pathway_table_container=ipyw.VBox(),
        reactome_diagram_container=ipyw.VBox(),
        df=df,
        adata=adata,
        active_selection=None,
        on_gene_selected=None,
        on_results_cleared=lambda *_: None,
        log=lambda *_: None,
        gene_annotation_species=6239,
        reactome_species=6239,
    )

    table = results_box.children[0]
    table.selected_gene = "unc-13"

    assert calls == [
        ("annotation", "unc-13", 6239),
        ("reactome", "unc-13", 6239),
    ]


def test_directional_gene_click_retries_description_with_display_name(monkeypatch):
    adata = AnnData(X=np.zeros((2, 1)))
    adata.var_names = ["WBGene00010957"]
    df = pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]})
    results_box = ipyw.VBox()
    details_box = ipyw.VBox()
    calls = []

    def fake_fetch_gene_description(gene, *, species="human"):
        calls.append(("annotation", gene, species))
        if gene == "nduo-6":
            return {"symbol": "nduo-6", "name": "NADH:ubiquinone oxidoreductase"}
        return None

    def fake_fetch_pathways(gene, *, species="human"):
        calls.append(("reactome", gene, species))
        return []

    monkeypatch.setattr(
        "scsketch._results.fetch_gene_description",
        fake_fetch_gene_description,
    )
    monkeypatch.setattr("scsketch._results.fetch_pathways", fake_fetch_pathways)

    show_directional_results(
        [
            [
                {
                    "attribute": "WBGene00010957",
                    "interval": (0.8, 0.001),
                    "reject": True,
                    "direction": "Selection 1",
                }
            ]
        ],
        selections_predicates=results_box,
        pathway_table_container=details_box,
        reactome_diagram_container=ipyw.VBox(),
        df=df,
        adata=adata,
        active_selection=None,
        on_gene_selected=None,
        on_results_cleared=lambda *_: None,
        log=lambda *_: None,
        gene_display_name=lambda gene: "nduo-6",
        gene_annotation_species=6239,
        reactome_species=6239,
    )

    table = results_box.children[0]
    table.selected_gene = "WBGene00010957"

    assert calls == [
        ("annotation", "WBGene00010957", 6239),
        ("annotation", "nduo-6", 6239),
        ("reactome", "WBGene00010957", 6239),
    ]
    assert "nduo-6" in details_box.children[1].children[0].value


def test_gene_expression_coloring_uses_display_name_in_caption_and_extra_view():
    adata = AnnData(
        X=np.array([[0.0], [1.0], [2.0]]),
        obs=pd.DataFrame({"cluster": ["A", "B", "A"]}),
        var=pd.DataFrame(
            {"gene_short_name": ["nduo-6"]},
            index=["WBGene00010957"],
        ),
    )
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    pca = Scatter(
        data=pd.DataFrame(
            {
                "PC1": [0.0, 1.0, 2.0],
                "PC2": [2.0, 1.0, 0.0],
                "cluster": ["A", "B", "A"],
            },
            index=adata.obs_names,
        ),
        x="PC1",
        y="PC2",
        color_by="cluster",
    )
    sketch = ScSketch(
        adata=adata,
        metadata_cols=["cluster"],
        color_by_default="cluster",
        extra_views={"PCA": pca},
    )

    sketch._color_embedding_by_gene("WBGene00010957")

    assert "<b>nduo-6</b> expression" in sketch._ctrl.gene_expression_caption.value
    assert "WBGene00010957 expression" not in sketch._ctrl.gene_expression_caption.value
    assert pca._color_labeling["variable"] == "nduo-6"
