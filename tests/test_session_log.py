from types import SimpleNamespace

import ipywidgets as ipyw
import numpy as np
from anndata import AnnData

from jscatter import Line
from scsketch._action_log import apply_action_log_entry
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
    sketch.logger = SimpleNamespace(exception=lambda *args, **kwargs: None)
    sketch.analysis_mode = "directional"
    sketch.active_selection = None
    sketch.action_log = []
    sketch._action_log_paused = False
    sketch._history_dropdown_paused = False
    sketch._selection_archive = {}
    sketch.scatter = None
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
    sketch._selection_archive = {
        selection.name: selection for selection in sketch.selections.selections
    }
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
    assert session["selection_archive"][0]["name"] == "Selection 1"
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


def test_recorded_action_log_is_exported_and_loaded():
    source = _sketch(_adata())
    source._record_action(
        "save_selection",
        "Saved selection Selection 1",
        selection="Selection 1",
        payload={"n_cells": 2},
    )
    session = source.export_session()
    target = _sketch(_adata())
    target.action_log = []

    target.load_session(session)

    assert session["action_log"][0]["type"] == "save_selection"
    assert session["action_log"][0]["selection"] == "Selection 1"
    assert target.action_log[0]["payload"] == {"n_cells": 2}


def test_action_log_row_can_restore_selection_state():
    sketch = _sketch(_adata())
    entry = sketch._record_action(
        "compute_directional",
        "Computed directional results for Selection 1",
        selection="Selection 1",
    )
    sketch.active_selection = None

    message = apply_action_log_entry(sketch, entry)

    assert message == "Restored selection Selection 1."
    assert sketch.active_selection.name == "Selection 1"


def test_select_gene_action_replays_gene_result_details():
    sketch = _sketch(_adata())
    calls = []
    sketch._ctrl = SimpleNamespace(
        directional_controls_box=SimpleNamespace(layout=SimpleNamespace(display="")),
        diff_controls_box=SimpleNamespace(layout=SimpleNamespace(display="")),
    )
    def show_directional_results(results, *, initial_gene=None):
        calls.append(("directional", results, initial_gene))

    sketch._show_directional_results = show_directional_results
    sketch._color_embedding_by_gene = lambda gene: calls.append(("gene", gene))
    entry = sketch._record_action(
        "select_gene",
        "Selected gene GeneA",
        selection="Selection 1",
        payload={"gene": "GeneA"},
    )

    message = apply_action_log_entry(sketch, entry)

    assert message == "Restored selection Selection 1 and gene GeneA."
    assert calls == [("directional", [sketch.active_selection.cached_results], "GeneA")]


def test_show_action_log_returns_widget():
    sketch = _sketch(_adata())
    sketch._record_action("save_selection", "Saved selection Selection 1")

    table = sketch.show_action_log()

    assert table.__class__.__name__ == "VBox"


def test_show_action_log_click_reconstructs_deleted_selection_timeline():
    sketch = _sketch(_adata())
    _add_de_selection(sketch)
    for selection in sketch.selections.selections:
        sketch._archive_selection(selection)
    sketch._record_action(
        "save_selection",
        "Saved selection Selection 1",
        selection="Selection 1",
    )
    sketch._record_action(
        "save_selection",
        "Saved selection Selection 2",
        selection="Selection 2",
    )
    removed = sketch.selections.selections.pop()
    sketch._archive_selection(removed)
    sketch._record_action(
        "remove_selection",
        "Removed selection Selection 2",
        selection="Selection 2",
    )

    table = sketch.show_action_log()
    rows_box = table.children[1]
    save_selection_2_row = rows_box.children[2]
    save_selection_2_button = save_selection_2_row.children[2]
    save_selection_2_button.click()

    assert [selection.name for selection in sketch.selections.selections] == [
        "Selection 1",
        "Selection 2",
    ]


def test_history_options_include_recorded_actions():
    sketch = _sketch(_adata())
    sketch._ctrl = SimpleNamespace(
        history_dropdown=ipyw.Dropdown(options=[("No recorded actions", None)]),
    )

    sketch._record_action(
        "save_selection",
        "Saved selection Selection 1",
        selection="Selection 1",
    )

    assert not sketch._ctrl.history_dropdown.disabled
    assert sketch._ctrl.history_dropdown.options == (
        ("Choose action...", None),
        ("1: Saved selection Selection 1", 0),
    )


def test_history_action_jump_reconstructs_deleted_selection():
    sketch = _sketch(_adata())
    _add_de_selection(sketch)
    for selection in sketch.selections.selections:
        sketch._archive_selection(selection)
    sketch._record_action(
        "save_selection",
        "Saved selection Selection 1",
        selection="Selection 1",
    )
    sketch._record_action(
        "save_selection",
        "Saved selection Selection 2",
        selection="Selection 2",
    )
    removed = sketch.selections.selections.pop()
    sketch._archive_selection(removed)
    sketch._record_action(
        "remove_selection",
        "Removed selection Selection 2",
        selection="Selection 2",
    )

    message = sketch._apply_history_action(1)

    assert message == "Restored selection Selection 2."
    assert [selection.name for selection in sketch.selections.selections] == [
        "Selection 1",
        "Selection 2",
    ]


def test_session_save_button_writes_json_file(tmp_path):
    sketch = _sketch(_adata())
    path = tmp_path / "saved.scsketch.json"
    sketch._ctrl = SimpleNamespace(
        history_dropdown=ipyw.Dropdown(options=[("No recorded actions", None)]),
        diff_t_threshold=ipyw.FloatText(value=2.0),
        diff_p_threshold=ipyw.FloatText(value=0.05),
        session_filename=ipyw.Text(value=str(path)),
        session_status=ipyw.HTML(""),
    )

    sketch._session_save_handler(None)

    assert path.exists()
    assert "Saved 1 actions" in sketch._ctrl.session_status.value


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


def test_build_playback_steps_prefers_action_log_when_present():
    sketch = _sketch(_adata())
    sketch._record_action(
        "save_selection",
        "Saved selection Selection 1",
        selection="Selection 1",
    )
    session = sketch.export_session()

    steps = build_playback_steps(session)

    assert [step["type"] for step in steps] == ["save_selection"]


def test_apply_playback_step_can_apply_action_log_entry():
    sketch = _sketch(_adata())
    sketch._record_action(
        "save_selection",
        "Saved selection Selection 1",
        selection="Selection 1",
    )
    session = sketch.export_session()
    target = _sketch(_adata())
    target.selections = SimpleNamespace(selections=[])
    target.active_selection = None

    step = apply_playback_step(target, session, 0)

    assert step["type"] == "save_selection"
    assert target.active_selection.name == "Selection 1"
    assert target.selections.selections[0].name == "Selection 1"


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


def test_apply_playback_step_can_rewind_before_removed_selection():
    sketch = _sketch(_adata())
    _add_de_selection(sketch)
    for selection in sketch.selections.selections:
        sketch._archive_selection(selection)
    sketch._record_action(
        "save_selection",
        "Saved selection Selection 1",
        selection="Selection 1",
    )
    sketch._record_action(
        "save_selection",
        "Saved selection Selection 2",
        selection="Selection 2",
    )
    removed = sketch.selections.selections.pop()
    sketch._archive_selection(removed)
    sketch._record_action(
        "remove_selection",
        "Removed selection Selection 2",
        selection="Selection 2",
    )
    session = sketch.export_session()
    target = _sketch(_adata())
    target.selections = SimpleNamespace(selections=[])
    target.active_selection = None

    step_before_delete = apply_playback_step(target, session, 1)
    assert step_before_delete["type"] == "save_selection"
    assert [selection.name for selection in target.selections.selections] == [
        "Selection 1",
        "Selection 2",
    ]

    delete_step = apply_playback_step(target, session, 2)
    assert delete_step["type"] == "remove_selection"
    assert [selection.name for selection in target.selections.selections] == [
        "Selection 1",
    ]
    assert {raw["name"] for raw in session["selection_archive"]} == {
        "Selection 1",
        "Selection 2",
    }


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
