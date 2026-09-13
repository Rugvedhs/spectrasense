"""Conformal guarantees are easy to break silently, so they are tested directly
against the property they promise."""

import numpy as np
import pytest

from ptgbench.conformal import (
    MondrianConformal,
    NormalizedConformal,
    SplitConformal,
    conformal_quantile,
    similarity_bins,
)


def test_quantile_uses_the_finite_sample_correction():
    scores = np.arange(1.0, 20.0)          # n = 19
    # ceil((19+1) * 0.9) = 18 -> the 18th smallest score.
    assert conformal_quantile(scores, 0.1) == 18.0


def test_quantile_is_infinite_when_calibration_is_too_small():
    # n = 8: ceil(9 * 0.9) = 9 > 8, so no finite quantile is valid.
    assert np.isinf(conformal_quantile(np.arange(1.0, 9.0), 0.1))


def test_quantile_of_empty_calibration_is_infinite():
    assert np.isinf(conformal_quantile(np.array([]), 0.1))


@pytest.mark.parametrize("alpha", [0.05, 0.1, 0.2])
def test_marginal_coverage_on_exchangeable_data(alpha):
    rng = np.random.default_rng(0)
    y = rng.normal(0, 50, 8000)
    pred = y + rng.normal(0, 20, 8000)
    predictor = SplitConformal(alpha).calibrate(y[:4000], pred[:4000])
    intervals = predictor.predict(pred[4000:])
    coverage = intervals.covers(y[4000:]).mean()
    assert coverage == pytest.approx(1 - alpha, abs=0.02)


def test_mondrian_gives_each_category_its_own_width():
    rng = np.random.default_rng(1)
    n = 6000
    category = rng.integers(0, 2, n)
    noise = np.where(category == 0, 5.0, 60.0)
    y = rng.normal(0, 30, n)
    pred = y + rng.normal(0, noise)

    predictor = MondrianConformal(0.1).calibrate(
        y[:3000], pred[:3000], categories_calib=category[:3000]
    )
    intervals = predictor.predict(pred[3000:], categories_test=category[3000:])
    widths = intervals.width
    easy = widths[category[3000:] == 0].mean()
    hard = widths[category[3000:] == 1].mean()
    assert hard > 3 * easy

    covered = intervals.covers(y[3000:])
    for c in (0, 1):
        assert covered[category[3000:] == c].mean() == pytest.approx(0.9, abs=0.03)


def test_mondrian_falls_back_for_an_unseen_category():
    """The family-holdout failure mode, in miniature."""
    rng = np.random.default_rng(2)
    y = rng.normal(0, 20, 400)
    pred = y + rng.normal(0, 5, 400)
    calib_categories = np.array(["seen"] * 400)
    predictor = MondrianConformal(0.1).calibrate(
        y, pred, categories_calib=calib_categories
    )
    predictor.predict(pred[:50], categories_test=np.array(["unseen"] * 50))
    assert predictor.n_fallback_ == 50


def test_normalised_conformal_widens_where_the_model_struggles():
    rng = np.random.default_rng(3)
    n = 6000
    X = rng.normal(0, 1, (n, 4))
    noise = np.where(X[:, 0] > 0, 60.0, 5.0)
    y = rng.normal(0, 30, n)
    pred = y + rng.normal(0, noise)

    predictor = NormalizedConformal(0.1, random_state=0)
    predictor.fit_difficulty(X[:2000], y[:2000], pred[:2000])
    predictor.calibrate(y[2000:4000], pred[2000:4000], X_calib=X[2000:4000])
    intervals = predictor.predict(pred[4000:], X_test=X[4000:])

    hard = intervals.width[X[4000:, 0] > 0].mean()
    easy = intervals.width[X[4000:, 0] <= 0].mean()
    assert hard > 2 * easy


def test_similarity_bins_are_ordered_and_bounded():
    values = np.array([0.0, 0.45, 0.55, 0.65, 0.75, 0.95])
    assert list(similarity_bins(values)) == [0, 1, 2, 3, 4, 5]
