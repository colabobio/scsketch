import numpy as np

from scsketch._scores import discovery_scores


def test_discovery_scores_normalize_nominal_p_values_to_zero_ten():
    scores = discovery_scores([0.05, 0.005, 0.0005])

    np.testing.assert_array_equal(scores, [0, 5, 10])
    assert np.issubdtype(scores.dtype, np.integer)


def test_discovery_scores_handles_identical_p_values():
    scores = discovery_scores([0.01, 0.01])

    np.testing.assert_array_equal(scores, [10, 10])
