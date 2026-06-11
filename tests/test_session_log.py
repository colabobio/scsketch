from types import SimpleNamespace

import numpy as np
from anndata import AnnData

from jscatter import Line
from scsketch._session import apply_playback_step, build_playback_steps
from scsketch._utils import Selection, Selections
from scsketch.scsketch import ScSketch


def _adata(obs_names=("cell_a", "cell_b", "cell_c")):
    adata = AnnData(
        X=np.zeros((len(obs_names), 2)),
        obs={"cluster": ["0"] * len(obs_names)},
    )
    adata.obs_names = list(obs_names)
    adata.var_names = ["GeneA", "GeneB"]
    adata.obsm["X_umap"] = np.array(
        [[float(i), float(i + 1)] for i in range(len(obs_names))]
    )
    return adata


def _sketch(adata):
    sketch = ScSketch.__new__(ScSketch)
    sketch.adata = adata
    sketch.metadata_cols = ["cluster"]
    sketch.color_by_default = "cluster"
    sketch.height = 720
    sketch.background_color = "#111111"
    sketch.max_genes = 0
    sketch.fdr_alpha = 0.05
    sketch.analysis_mode = "directional"
    sketch.active_selection = None
    sketch._ctrl = None
    sketch.selections = Selections(
        selections=[
            Selection(
                index=1,
                name="Selection 1",
                points=np.array([0, 2]),
                color="#00dadb",
                lasso=Line([[0, 0], [1, 0], [1, 1], [0, 0]], line_color="#00dadb"),
                hull=Line([[0, 0], [1, 1], [0, 0]], line_color="#00dadb"),
                path=np.array([[0.0, 0.0], [1.0, 1.0]]),
                cached_results=[
                    {
                        "attribute": "GeneA",
                        "interval": (0.8, 0.001),
                        "quality": 0.8,
                        "reject": True,
                        "direction": "Selection 1",
                    }
                ],
            )
        ]
    )
    sketch.active_selection = sketch.selections.selections[0]
    return sketch


def _add_de_selection(sketch):
    sketch.selections.selections.append(
        Selection(
            index=2,
            name="Selection 2",
            points=np.array([1]),
            color="#da00db",
            lasso=Line([[2, 2], [3, 2], [2, 2]], line_color="#da00db"),
            hull=Line([[2, 2], [3, 2], [2, 2]], line_color="#da00db"),
            cached_diffexpr=[
                {
                    "attribute": "GeneB",
                    "interval": (3.2, 0.02),
                    "quality": 3.2,
                    "direction": "Selection 2",
                }
            ],
        )
    )


def test_export_session_includes_replayable_selection_state():
    sketch = _sketch(_adata())

    session = sketch.export_session()

    assert session["schema_version"] == "1.0"
    assert session["dataset"]["n_obs"] == 3
    assert session["dataset"]["n_vars"] == 2
    assert session["initial_config"]["metadata_cols"] == ["cluster"]
    assert session["state"]["active_selection"] == "Selection 1"

    selection = session["selections"][0]
    assert selection["name"] == "Selection 1"
    assert selection["points_indices"] == [0, 2]
    assert selection["points_obs_names"] == ["cell_a", "cell_c"]
    assert selection["lasso_polygon"] == [
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 0.0],
    ]
    assert selection["path"] == [[0.0, 0.0], [1.0, 1.0]]
    assert selection["cached_results"][0]["interval"] == [0.8, 0.001]
    assert session["steps"] == [
        {
            "type": "restore_selection",
            "selection": "Selection 1",
            "label": "Restore selection: Selection 1",
        },
        {
            "type": "show_directional_results",
            "selection": "Selection 1",
            "label": "Show directional results: Selection 1",
        },
    ]


def test_load_session_restores_points_by_obs_names_when_order_changes():
    source = _sketch(_adata())
    session = source.export_session()
    target = _sketch(_adata(obs_names=("cell_c", "cell_b", "cell_a")))
    target.selections = SimpleNamespace(selections=[])
    target.active_selection = None

    warnings = target.load_session(session)

    restored = target.selections.selections[0]
    assert restored.name == "Selection 1"
    assert restored.points.tolist() == [2, 0]
    assert restored.cached_results[0]["attribute"] == "GeneA"
    assert any("obs_names_hash" in warning for warning in warnings)


def test_session_can_round_trip_through_json_file(tmp_path):
    source = _sketch(_adata())
    path = tmp_path / "analysis.scsketch.json"

    written = source.export_session(path)

    target = _sketch(_adata())
    target.selections = SimpleNamespace(selections=[])
    target.active_selection = None
    warnings = target.load_session(path)

    assert warnings == []
    assert written["selections"][0]["name"] == "Selection 1"
    assert target.active_selection.name == "Selection 1"
    assert target.selections.selections[0].points.tolist() == [0, 2]


def test_build_playback_steps_synthesizes_legacy_session_steps():
    sketch = _sketch(_adata())
    _add_de_selection(sketch)
    session = sketch.export_session()
    session.pop("steps")

    steps = build_playback_steps(session)

    assert [step["type"] for step in steps] == [
        "restore_selection",
        "show_directional_results",
        "restore_selection",
        "show_diffexpr_results",
    ]


def test_apply_playback_step_restores_state_through_step():
    sketch = _sketch(_adata())
    _add_de_selection(sketch)
    session = sketch.export_session()
    target = _sketch(_adata())
    target.selections = SimpleNamespace(selections=[])
    target.active_selection = None

    step = apply_playback_step(target, session, 3)

    assert step["type"] == "show_diffexpr_results"
    assert target.analysis_mode == "differential"
    assert target.active_selection.name == "Selection 2"
    assert [selection.name for selection in target.selections.selections] == [
        "Selection 1",
        "Selection 2",
    ]


def test_show_session_player_returns_widget_and_applies_first_step():
    sketch = _sketch(_adata())
    _add_de_selection(sketch)
    session = sketch.export_session()
    target = _sketch(_adata())
    target.selections = SimpleNamespace(selections=[])
    target.active_selection = None

    player = target.show_session_player(session)

    assert player.__class__.__name__ == "VBox"
    assert target.active_selection.name == "Selection 1"
