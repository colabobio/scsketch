from types import SimpleNamespace

from scsketch.scsketch import ScSketch


def test_get_diffexpr_genes_exports_cached_de_results():
    sketch = ScSketch.__new__(ScSketch)
    sketch.selections = SimpleNamespace(
        selections=[
            SimpleNamespace(
                name="Cluster 6 3dpi portion",
                cached_diffexpr=[
                    {
                        "attribute": "GENE_LOW",
                        "interval": (2.1, 0.01),
                        "direction": "Cluster 6 3dpi portion",
                    },
                    {
                        "attribute": "GENE_HIGH",
                        "interval": (-5.5, 0.001),
                        "direction": "Cluster 6 3dpi portion",
                    },
                ],
            )
        ]
    )

    genes = sketch.get_diffexpr_genes("Cluster 6 3dpi portion")

    assert genes.to_dict(orient="records") == [
        {
            "gene": "GENE_HIGH",
            "t-statistic": -5.5,
            "p-value": 0.001,
            "selection": "Cluster 6 3dpi portion",
        },
        {
            "gene": "GENE_LOW",
            "t-statistic": 2.1,
            "p-value": 0.01,
            "selection": "Cluster 6 3dpi portion",
        },
    ]


def test_get_de_genes_aliases_get_diffexpr_genes():
    sketch = ScSketch.__new__(ScSketch)
    sketch.selections = SimpleNamespace(
        selections=[
            SimpleNamespace(
                name="Selection 1",
                cached_diffexpr=[
                    {
                        "attribute": "GENE1",
                        "interval": (3.0, 0.02),
                        "direction": "Selection 1",
                    }
                ],
            )
        ]
    )

    assert sketch.get_de_genes("Selection 1").equals(
        sketch.get_diffexpr_genes("Selection 1")
    )


def test_get_diffexpr_genes_returns_empty_dataframe_without_cached_de():
    sketch = ScSketch.__new__(ScSketch)
    sketch.selections = SimpleNamespace(
        selections=[
            SimpleNamespace(name="Selection 1", cached_diffexpr=None),
        ]
    )

    assert sketch.get_diffexpr_genes("Selection 1").empty
