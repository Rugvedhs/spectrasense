"""Point-accuracy, interval and subgroup-coverage metrics.

Point metrics answer "how close on average"; interval metrics answer "how often
is the truth inside the stated range, and how wide did that range have to be".
Both are needed: an interval predictor can buy coverage with uselessly wide
intervals, and a point predictor with a good MAE can still be confidently wrong
on precisely the unfamiliar chemistry a screening campaign cares about.

``R^2`` is reported but never used alone to rank difficulty.  It is normalised by
the variance of whichever test set it is computed on, so a family with a narrow
Tg spread can show a deeply negative R^2 at a smaller absolute error than a
family with a wide spread — an artefact of the denominator, not a worse model.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .conformal import Intervals


def point_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    residual = y_true - y_pred
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    spearman = (
        float(stats.spearmanr(y_true, y_pred).statistic) if y_true.size > 2 else np.nan
    )
    return {
        "n": int(y_true.size),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        "spearman": spearman,
        "bias": float(np.mean(residual)),
    }


def interval_metrics(
    intervals: Intervals, y_true: np.ndarray, alpha: float = 0.1
) -> dict[str, float]:
    covered = intervals.covers(np.asarray(y_true, dtype=float))
    width = intervals.width
    finite = np.isfinite(width)
    return {
        "coverage": float(np.mean(covered)),
        "target_coverage": 1.0 - alpha,
        "coverage_error": float(np.mean(covered) - (1.0 - alpha)),
        "mean_width": float(np.mean(width[finite])) if finite.any() else np.inf,
        "median_width": float(np.median(width[finite])) if finite.any() else np.inf,
        "frac_infinite_width": float(np.mean(~finite)),
    }


def conditional_coverage(
    intervals: Intervals,
    y_true: np.ndarray,
    groups: np.ndarray,
    alpha: float = 0.1,
    min_size: int = 10,
) -> pd.DataFrame:
    """Coverage and mean width within each subgroup.

    Subgroups smaller than ``min_size`` are reported but flagged, because a
    coverage estimate from a handful of points is mostly binomial noise and
    should not drive a conclusion about conditional validity.
    """
    covered = intervals.covers(np.asarray(y_true, dtype=float))
    width = intervals.width
    rows = []
    for group in pd.unique(groups):
        mask = groups == group
        n = int(mask.sum())
        finite = np.isfinite(width[mask])
        rows.append(
            {
                "group": group,
                "n": n,
                "coverage": float(np.mean(covered[mask])),
                "coverage_error": float(np.mean(covered[mask]) - (1.0 - alpha)),
                "mean_width": float(np.mean(width[mask][finite]))
                if finite.any()
                else np.inf,
                "reliable": n >= min_size,
            }
        )
    return pd.DataFrame(rows).sort_values("coverage").reset_index(drop=True)


def worst_subgroup_coverage(
    conditional: pd.DataFrame, reliable_only: bool = True
) -> dict[str, float]:
    """Summarise how badly conditional coverage departs from the target.

    ``max_coverage_gap`` is the headline number: the largest shortfall below the
    nominal level across subgroups.  Marginal coverage can be exactly on target
    while this is severe, which is the failure mode the paper is about.
    """
    frame = conditional[conditional["reliable"]] if reliable_only else conditional
    if frame.empty:
        return {"worst_coverage": np.nan, "max_coverage_gap": np.nan, "n_groups": 0}
    return {
        "worst_coverage": float(frame["coverage"].min()),
        "max_coverage_gap": float(-frame["coverage_error"].min()),
        "mean_abs_coverage_error": float(frame["coverage_error"].abs().mean()),
        "n_groups": int(len(frame)),
    }


def coverage_width_criterion(
    intervals: Intervals, y_true: np.ndarray, alpha: float = 0.1, eta: float = 30.0
) -> float:
    """Penalised width: narrow intervals are only rewarded if coverage holds.

    Follows the usual CWC construction — mean width inflated by an exponential
    penalty when empirical coverage falls below the nominal level. Lower is
    better. ``eta`` sets how harshly undercoverage is punished.
    """
    stats_ = interval_metrics(intervals, y_true, alpha)
    coverage, width = stats_["coverage"], stats_["mean_width"]
    target = 1.0 - alpha
    penalty = np.exp(-eta * (coverage - target)) if coverage < target else 1.0
    return float(width * penalty)


def sparsification_error(
    y_true: np.ndarray, y_pred: np.ndarray, uncertainty: np.ndarray, n_steps: int = 20
) -> dict[str, float]:
    """Area under the sparsification error curve (AUSE).

    Repeatedly discard the points the model claims are most uncertain and watch
    the error on what remains.  A useful uncertainty ranks errors, so its curve
    should track the oracle curve obtained by discarding the truly worst points.
    AUSE is the gap between them; zero is a perfect ranking.
    """
    residual = np.abs(np.asarray(y_true, float) - np.asarray(y_pred, float))
    n = residual.size
    by_uncertainty = np.argsort(-np.asarray(uncertainty, float))
    by_error = np.argsort(-residual)

    fractions = np.linspace(0, 0.8, n_steps)
    model_curve, oracle_curve = [], []
    for f in fractions:
        keep = n - int(round(f * n))
        model_curve.append(residual[by_uncertainty[n - keep :]].mean())
        oracle_curve.append(residual[by_error[n - keep :]].mean())

    model_curve = np.asarray(model_curve)
    oracle_curve = np.asarray(oracle_curve)
    scale = model_curve[0] if model_curve[0] > 0 else 1.0
    return {
        "ause": float(np.trapezoid(model_curve - oracle_curve, fractions) / scale),
        "sparsification_gain": float((model_curve[0] - model_curve[-1]) / scale),
    }
