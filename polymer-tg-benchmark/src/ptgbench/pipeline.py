"""Running one split end to end: fit, calibrate, predict, score.

Every training pool is divided three ways rather than two:

``fit``         trains the point predictor.
``difficulty``  supplies *out-of-sample* residuals for the normalised predictor's
                difficulty model.  Fitting that model on the point predictor's
                own training residuals would teach it that the model is
                everywhere accurate, which is exactly the wrong lesson.
``calibration`` supplies the conformal quantile.

All conformal variants share one point predictor and one calibration set, so any
difference between them comes from how they turn residuals into intervals and
not from a difference in the underlying regression.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .conformal import (
    Intervals,
    MondrianConformal,
    NormalizedConformal,
    SplitConformal,
    similarity_bins,
)
from .features import max_similarity_to_reference
from .metrics import (
    conditional_coverage,
    coverage_width_criterion,
    interval_metrics,
    point_metrics,
    sparsification_error,
    worst_subgroup_coverage,
)
from .models import build_model, ensemble_std
from .splits import Split


@dataclass
class SplitResult:
    point: pd.DataFrame
    intervals: pd.DataFrame
    conditional: pd.DataFrame
    predictions: pd.DataFrame


def three_way_split(
    train_idx: np.ndarray,
    seed: int,
    difficulty_fraction: float = 0.2,
    calib_fraction: float = 0.2,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    perm = rng.permutation(train_idx)
    n = len(perm)
    n_diff = max(1, int(round(n * difficulty_fraction)))
    n_calib = max(1, int(round(n * calib_fraction)))
    return (
        np.sort(perm[n_diff + n_calib :]),
        np.sort(perm[n_calib : n_calib + n_diff]),
        np.sort(perm[:n_calib]),
    )


def run_split(
    split: Split,
    model_name: str,
    X: np.ndarray,
    y: np.ndarray,
    fingerprints: list,
    families: np.ndarray,
    alpha: float = 0.1,
    seed: int = 0,
    with_conformal: bool = True,
) -> SplitResult:
    fit_idx, diff_idx, calib_idx = three_way_split(split.train_idx, seed)
    test_idx = split.test_idx

    model = build_model(model_name, random_state=seed)
    model.fit(X[fit_idx], y[fit_idx])

    pred_test = model.predict(X[test_idx])
    point = point_metrics(y[test_idx], pred_test)
    point.update(
        {
            "split": split.name,
            "regime": split.regime,
            "arm": split.arm,
            "group": split.group,
            "model": model_name,
            "n_train": len(split.train_idx),
            "n_fit": len(fit_idx),
        }
    )

    # Structural novelty is measured against the data the point model actually
    # saw, not the whole training pool, so it matches what the model can know.
    fit_fps = [fingerprints[i] for i in fit_idx]
    sim_test = max_similarity_to_reference([fingerprints[i] for i in test_idx], fit_fps)
    point["mean_nn_similarity"] = float(np.mean(sim_test))

    predictions = pd.DataFrame(
        {
            "split": split.name,
            "regime": split.regime,
            "arm": split.arm,
            "group": split.group,
            "model": model_name,
            "index": test_idx,
            "y_true": y[test_idx],
            "y_pred": pred_test,
            "nn_similarity": sim_test,
            "family": families[test_idx],
        }
    )

    if not with_conformal:
        empty = pd.DataFrame()
        return SplitResult(pd.DataFrame([point]), empty, empty, predictions)

    pred_calib = model.predict(X[calib_idx])
    pred_diff = model.predict(X[diff_idx])
    sim_calib = max_similarity_to_reference(
        [fingerprints[i] for i in calib_idx], fit_fps
    )

    bins_calib = similarity_bins(sim_calib)
    bins_test = similarity_bins(sim_test)

    predictors = [
        ("SCP", SplitConformal(alpha), {}, {}),
        (
            "Mondrian-similarity",
            MondrianConformal(alpha, name="Mondrian-similarity"),
            {"categories_calib": bins_calib},
            {"categories_test": bins_test},
        ),
        (
            "Mondrian-family",
            MondrianConformal(alpha, name="Mondrian-family"),
            {"categories_calib": families[calib_idx]},
            {"categories_test": families[test_idx]},
        ),
    ]

    normalised = NormalizedConformal(alpha, random_state=seed)
    normalised.fit_difficulty(X[diff_idx], y[diff_idx], pred_diff)
    predictors.append(
        (
            "Normalised SCP",
            normalised,
            {"X_calib": X[calib_idx]},
            {"X_test": X[test_idx]},
        )
    )

    interval_rows, conditional_rows = [], []
    for label, predictor, calib_kwargs, test_kwargs in predictors:
        predictor.calibrate(y[calib_idx], pred_calib, **calib_kwargs)
        iv = predictor.predict(pred_test, **test_kwargs)

        row = interval_metrics(iv, y[test_idx], alpha)
        row.update(
            {
                "split": split.name,
                "regime": split.regime,
                "arm": split.arm,
                "group": split.group,
                "model": model_name,
                "conformal": label,
                "cwc": coverage_width_criterion(iv, y[test_idx], alpha),
            }
        )
        if isinstance(predictor, MondrianConformal):
            row["fallback_fraction"] = predictor.n_fallback_ / max(1, len(test_idx))
        else:
            row["fallback_fraction"] = 0.0

        by_similarity = conditional_coverage(iv, y[test_idx], bins_test, alpha)
        row.update(
            {
                f"sim_{k}": v
                for k, v in worst_subgroup_coverage(by_similarity).items()
            }
        )
        by_family = conditional_coverage(iv, y[test_idx], families[test_idx], alpha)
        row.update(
            {f"fam_{k}": v for k, v in worst_subgroup_coverage(by_family).items()}
        )
        interval_rows.append(row)

        for frame, taxonomy in ((by_similarity, "similarity"), (by_family, "family")):
            frame = frame.copy()
            frame["taxonomy"] = taxonomy
            frame["conformal"] = label
            frame["split"] = split.name
            frame["regime"] = split.regime
            frame["model"] = model_name
            conditional_rows.append(frame)

        predictions[f"lower_{label}"] = iv.lower
        predictions[f"upper_{label}"] = iv.upper

    spread = ensemble_std(model, X[test_idx])
    if spread is not None:
        point.update(sparsification_error(y[test_idx], pred_test, spread))

    return SplitResult(
        pd.DataFrame([point]),
        pd.DataFrame(interval_rows),
        pd.concat(conditional_rows, ignore_index=True),
        predictions,
    )
