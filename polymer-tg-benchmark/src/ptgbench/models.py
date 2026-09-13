"""The regressor zoo, wired so that every data-dependent step is fitted in-fold.

Each entry is a scikit-learn pipeline beginning with :class:`DescriptorFilter`.
That placement is deliberate: deciding which descriptor columns are constant or
degenerate is a decision informed by data, so doing it once over the full table —
as is common in this literature — leaks information about the test structures
into the feature definition.  Inside a pipeline the filter sees only the
training fold.
"""

from __future__ import annotations

import numpy as np
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import (
    ExtraTreesRegressor,
    HistGradientBoostingRegressor,
    RandomForestRegressor,
)
from sklearn.kernel_ridge import KernelRidge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

from .features import DescriptorFilter

RANDOM_STATE = 0


def _tree_pipeline(estimator) -> Pipeline:
    # Trees are scale-invariant, so no scaler: it would only cost time and blur
    # the descriptor units that make importances interpretable.
    return Pipeline([("filter", DescriptorFilter()), ("model", estimator)])


def _scaled_pipeline(estimator) -> Pipeline:
    return Pipeline(
        [
            ("filter", DescriptorFilter()),
            ("scale", StandardScaler()),
            ("model", estimator),
        ]
    )


def build_model(name: str, random_state: int = RANDOM_STATE) -> Pipeline:
    """Instantiate one model by name."""
    if name == "median":
        return _tree_pipeline(DummyRegressor(strategy="median"))
    if name == "random_forest":
        return _tree_pipeline(
            RandomForestRegressor(
                n_estimators=500, min_samples_leaf=1, n_jobs=-1,
                random_state=random_state,
            )
        )
    if name == "extra_trees":
        return _tree_pipeline(
            ExtraTreesRegressor(
                n_estimators=500, min_samples_leaf=1, n_jobs=-1,
                random_state=random_state,
            )
        )
    if name == "hist_gbr":
        return _tree_pipeline(
            HistGradientBoostingRegressor(
                max_iter=500, learning_rate=0.06, max_leaf_nodes=31,
                early_stopping=False, random_state=random_state,
            )
        )
    if name == "svr":
        return _scaled_pipeline(SVR(C=100.0, epsilon=5.0, gamma="scale"))
    if name == "krr":
        return _scaled_pipeline(KernelRidge(alpha=1.0, kernel="rbf", gamma=None))
    raise ValueError(f"unknown model {name!r}")


MODEL_NAMES = ("median", "random_forest", "extra_trees", "hist_gbr", "svr")


def ensemble_std(pipeline: Pipeline, X: np.ndarray) -> np.ndarray:
    """Spread of individual tree predictions, as a cheap epistemic signal.

    Defined only for forest-style ensembles; returns ``None`` otherwise so
    callers can fall back to a fitted difficulty model.
    """
    model = pipeline.named_steps.get("model")
    if not hasattr(model, "estimators_"):
        return None
    transformed = X
    for step_name, step in pipeline.steps[:-1]:
        transformed = step.transform(transformed)
    per_tree = np.stack([tree.predict(transformed) for tree in model.estimators_])
    return per_tree.std(axis=0)
