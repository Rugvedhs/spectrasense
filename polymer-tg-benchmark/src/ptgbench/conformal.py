"""Split-conformal predictors and the conditional variants tested against them.

Split conformal prediction turns any point predictor into interval predictions
with a finite-sample *marginal* coverage guarantee, assuming only that
calibration and test data are exchangeable.  That assumption is exactly what a
family holdout breaks, and marginal coverage is a weak promise in any case: a
predictor can hit 90% overall while covering 99% of easy repeat units and 60% of
the hard ones.  Since screening decisions are made on the hard ones, the
quantity that matters is *conditional* coverage.

Three predictors are compared:

:class:`SplitConformal`        one global quantile of absolute residuals.
:class:`NormalizedConformal`   residuals scaled by a fitted difficulty estimate,
                               which widens intervals where the model expects to
                               struggle instead of everywhere at once.
:class:`MondrianConformal`     a separate quantile per category, restoring exact
                               coverage within each category by construction —
                               provided the category has calibration data.

The last proviso is the crux of the family-holdout result.  A Mondrian taxonomy
built on polymer family cannot help for an unseen family: by construction it has
no calibration members.  A taxonomy built on *structural novelty* — how far the
nearest training analogue is — is defined for any repeat unit, seen family or
not, and is therefore the variant that can survive the shift.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor


def conformal_quantile(scores: np.ndarray, alpha: float) -> float:
    """The finite-sample-valid conformal quantile of calibration scores.

    Uses the ``ceil((n+1)(1-alpha))``-th order statistic, the correction that
    makes coverage at least ``1-alpha`` for finite ``n`` rather than only
    asymptotically.  When the calibration set is too small for the requested
    level, the honest answer is an infinite interval, not a quietly narrow one.
    """
    scores = np.asarray(scores, dtype=float)
    scores = scores[np.isfinite(scores)]
    n = scores.size
    if n == 0:
        return np.inf
    k = int(np.ceil((n + 1) * (1.0 - alpha)))
    if k > n:
        return np.inf
    return float(np.sort(scores)[k - 1])


@dataclass
class Intervals:
    lower: np.ndarray
    upper: np.ndarray
    point: np.ndarray

    @property
    def width(self) -> np.ndarray:
        return self.upper - self.lower

    def covers(self, y_true: np.ndarray) -> np.ndarray:
        return (y_true >= self.lower) & (y_true <= self.upper)


class SplitConformal:
    """Standard split conformal with absolute-residual nonconformity."""

    name = "SCP"

    def __init__(self, alpha: float = 0.1):
        self.alpha = alpha

    def calibrate(self, y_calib: np.ndarray, pred_calib: np.ndarray, **_) -> "SplitConformal":
        self.q_ = conformal_quantile(np.abs(y_calib - pred_calib), self.alpha)
        return self

    def predict(self, pred_test: np.ndarray, **_) -> Intervals:
        return Intervals(pred_test - self.q_, pred_test + self.q_, pred_test)


class NormalizedConformal:
    """Locally weighted conformal: residuals divided by a difficulty estimate.

    The difficulty model is a random forest fitted to predict ``log`` absolute
    residuals of the point predictor on held-out data.  Working in log space
    keeps the estimate positive and stops a handful of extreme residuals from
    dominating the fit.  ``beta`` floors the scale so that a confidently-wrong
    region cannot collapse the interval to zero width.
    """

    name = "Normalised SCP"

    def __init__(self, alpha: float = 0.1, beta: float = 1.0, random_state: int = 0):
        self.alpha = alpha
        self.beta = beta
        self.random_state = random_state

    def fit_difficulty(
        self, X_diff: np.ndarray, y_diff: np.ndarray, pred_diff: np.ndarray
    ) -> "NormalizedConformal":
        residual = np.abs(y_diff - pred_diff)
        self.difficulty_ = RandomForestRegressor(
            n_estimators=200,
            min_samples_leaf=5,
            n_jobs=-1,
            random_state=self.random_state,
        ).fit(X_diff, np.log1p(residual))
        return self

    def _scale(self, X: np.ndarray) -> np.ndarray:
        return np.expm1(self.difficulty_.predict(X)) + self.beta

    def calibrate(
        self, y_calib: np.ndarray, pred_calib: np.ndarray, X_calib: np.ndarray = None, **_
    ) -> "NormalizedConformal":
        scale = self._scale(X_calib)
        self.q_ = conformal_quantile(np.abs(y_calib - pred_calib) / scale, self.alpha)
        return self

    def predict(self, pred_test: np.ndarray, X_test: np.ndarray = None, **_) -> Intervals:
        half = self.q_ * self._scale(X_test)
        return Intervals(pred_test - half, pred_test + half, pred_test)


class MondrianConformal:
    """Category-conditional conformal: one quantile per category.

    Categories with fewer than ``min_per_bin`` calibration points fall back to
    the global quantile.  The fallback is reported rather than hidden: a
    taxonomy that has to fall back for the very cases it was meant to protect has
    not solved the problem, and that is precisely what happens to a
    family-conditioned taxonomy under family holdout.
    """

    def __init__(self, alpha: float = 0.1, min_per_bin: int = 20, name: str = "Mondrian SCP"):
        self.alpha = alpha
        self.min_per_bin = min_per_bin
        self.name = name

    def calibrate(
        self,
        y_calib: np.ndarray,
        pred_calib: np.ndarray,
        categories_calib: np.ndarray = None,
        **_,
    ) -> "MondrianConformal":
        scores = np.abs(y_calib - pred_calib)
        self.global_q_ = conformal_quantile(scores, self.alpha)
        self.q_by_category_ = {}
        self.n_by_category_ = {}
        for category in np.unique(categories_calib):
            mask = categories_calib == category
            self.n_by_category_[category] = int(mask.sum())
            if mask.sum() >= self.min_per_bin:
                self.q_by_category_[category] = conformal_quantile(
                    scores[mask], self.alpha
                )
        return self

    def predict(
        self, pred_test: np.ndarray, categories_test: np.ndarray = None, **_
    ) -> Intervals:
        half = np.full(pred_test.shape, self.global_q_, dtype=float)
        self.n_fallback_ = 0
        for i, category in enumerate(categories_test):
            if category in self.q_by_category_:
                half[i] = self.q_by_category_[category]
            else:
                self.n_fallback_ += 1
        return Intervals(pred_test - half, pred_test + half, pred_test)


def similarity_bins(
    similarity: np.ndarray, edges: tuple[float, ...] = (0.4, 0.5, 0.6, 0.7, 0.8)
) -> np.ndarray:
    """Discretise nearest-training Tanimoto similarity into ordered bins.

    The edges bracket the range where structure-activity reasoning is usually
    held to break down; using fixed edges rather than data-driven quantiles keeps
    the taxonomy comparable across split regimes, which is required if the bins
    are to be used as a Mondrian category under distribution shift.
    """
    return np.digitize(similarity, np.asarray(edges, dtype=float))


class FamilyOutCalibration:
    """Calibration residuals produced by models that had never seen the family.

    Split conformal draws its calibration set at random from the training pool, so
    every calibration residual describes *interpolation*: the model had relatives of
    that repeat unit in training.  Under a family holdout the test residuals
    describe extrapolation instead, the exchangeability assumption fails, and the
    quantile is calibrated against the wrong population.

    This calibrator repairs the mismatch at its source.  The training pool is
    partitioned by family; for each held-out training family a model is refitted on
    the remainder and used to predict it.  The pooled residuals then describe the
    same kind of prediction the test set will demand — one made without local
    chemistry — so the resulting quantile is calibrated for extrapolation.

    The cost is the usual cross-conformal caveat: calibration residuals come from
    models fitted on slightly less data than the deployed one, which makes the
    resulting intervals mildly conservative rather than anti-conservative. That is
    the safe direction for a screening decision.
    """

    def __init__(self, model_factory, n_families: int = 8, min_family_size: int = 30):
        self.model_factory = model_factory
        self.n_families = n_families
        self.min_family_size = min_family_size

    def residual_pool(
        self,
        X: np.ndarray,
        y: np.ndarray,
        pool_idx: np.ndarray,
        families: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return the calibrated indices and the out-of-family predictions for them."""
        pool_families = families[pool_idx]
        counts = pd.Series(pool_families).value_counts()
        eligible = counts[counts >= self.min_family_size].index.tolist()
        chosen = eligible[: self.n_families]

        indices, predictions = [], []
        for family in chosen:
            held = pool_idx[pool_families == family]
            rest = pool_idx[pool_families != family]
            if len(held) == 0 or len(rest) < 50:
                continue
            model = self.model_factory()
            model.fit(X[rest], y[rest])
            indices.append(held)
            predictions.append(model.predict(X[held]))

        if not indices:
            return np.array([], dtype=int), np.array([], dtype=float)
        return np.concatenate(indices), np.concatenate(predictions)
