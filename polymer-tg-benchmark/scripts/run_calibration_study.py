"""Does the conformal calibration set itself have to be built for extrapolation?

Split conformal draws calibration points at random from the training pool, so
every calibration residual describes interpolation.  Under a family holdout the
test residuals describe extrapolation, and the quantile is calibrated against the
wrong population.  This study compares the standard random calibration set
against one built by holding families out *within* the training pool, under both
a global and a novelty-conditioned quantile.
"""

from __future__ import annotations

import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd

import _bootstrap  # noqa: F401
from ptgbench.conformal import (
    FamilyOutCalibration,
    MondrianConformal,
    SplitConformal,
    similarity_bins,
)
from ptgbench.features import max_similarity_to_reference
from ptgbench.metrics import (
    conditional_coverage,
    interval_metrics,
    worst_subgroup_coverage,
)
from ptgbench.models import build_model
from ptgbench.pipeline import three_way_split
from ptgbench.splits import family_holdout_splits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-dir", default="data/processed")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--model", default="extra_trees")
    parser.add_argument("--alpha", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--n-calibration-families", type=int, default=8)
    args = parser.parse_args()

    processed = Path(args.processed_dir)
    table = pd.read_parquet(processed / "structures.parquet")
    X = pd.read_parquet(processed / "descriptors.parquet").to_numpy()
    y = table["Tg"].to_numpy()
    families = table["family"].to_numpy()
    with open(processed / "fingerprints.pkl", "rb") as handle:
        fps = pickle.load(handle)

    splits = family_holdout_splits(families)
    print(f"family-holdout splits: {len(splits)}", flush=True)

    rows, conditional_rows = [], []
    started = time.time()
    for n_done, split in enumerate(splits, start=1):
        fit_idx, _, calib_idx = three_way_split(split.train_idx, args.seed)
        test_idx = split.test_idx

        model = build_model(args.model, random_state=args.seed)
        model.fit(X[fit_idx], y[fit_idx])
        pred_test = model.predict(X[test_idx])

        fit_fps = [fps[i] for i in fit_idx]
        sim_test = max_similarity_to_reference([fps[i] for i in test_idx], fit_fps)
        bins_test = similarity_bins(sim_test)

        # --- calibration source 1: the usual random draw from the training pool
        pred_calib = model.predict(X[calib_idx])
        sim_calib = max_similarity_to_reference([fps[i] for i in calib_idx], fit_fps)
        sources = {
            "random": (y[calib_idx], pred_calib, similarity_bins(sim_calib))
        }

        # --- calibration source 2: residuals from models blind to each family
        calibrator = FamilyOutCalibration(
            lambda: build_model(args.model, random_state=args.seed),
            n_families=args.n_calibration_families,
        )
        fo_idx, fo_pred = calibrator.residual_pool(X, y, fit_idx, families)
        if fo_idx.size:
            fo_sim = np.zeros(fo_idx.size, dtype=float)
            for family in np.unique(families[fo_idx]):
                mask = families[fo_idx] == family
                reference = [fps[i] for i in fit_idx if families[i] != family]
                fo_sim[mask] = max_similarity_to_reference(
                    [fps[i] for i in fo_idx[mask]], reference
                )
            sources["family_out"] = (y[fo_idx], fo_pred, similarity_bins(fo_sim))

        for source_name, (y_cal, p_cal, bins_cal) in sources.items():
            for label, predictor, calib_kwargs, test_kwargs in (
                ("SCP", SplitConformal(args.alpha), {}, {}),
                ("Mondrian-similarity",
                 MondrianConformal(args.alpha, name="Mondrian-similarity"),
                 {"categories_calib": bins_cal}, {"categories_test": bins_test}),
            ):
                predictor.calibrate(y_cal, p_cal, **calib_kwargs)
                iv = predictor.predict(pred_test, **test_kwargs)

                row = interval_metrics(iv, y[test_idx], args.alpha)
                row.update({
                    "split": split.name, "family": split.group,
                    "model": args.model, "conformal": label,
                    "calibration": source_name, "n_calibration": len(y_cal),
                    "n_test": len(test_idx),
                    "mean_nn_similarity": float(np.mean(sim_test)),
                })
                by_similarity = conditional_coverage(iv, y[test_idx], bins_test, args.alpha)
                row.update({f"sim_{k}": v
                            for k, v in worst_subgroup_coverage(by_similarity).items()})
                if isinstance(predictor, MondrianConformal):
                    row["fallback_fraction"] = predictor.n_fallback_ / max(1, len(test_idx))
                else:
                    row["fallback_fraction"] = 0.0
                rows.append(row)

                frame = by_similarity.copy()
                frame["split"] = split.name
                frame["conformal"] = label
                frame["calibration"] = source_name
                conditional_rows.append(frame)

        elapsed = time.time() - started
        print(f"[{n_done}/{len(splits)}] {split.group:20s} "
              f"eta={elapsed / n_done * (len(splits) - n_done) / 60:5.1f}min", flush=True)

    results = Path(args.results_dir)
    results.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(results / "calibration_study.csv", index=False)
    pd.concat(conditional_rows, ignore_index=True).to_csv(
        results / "calibration_study_conditional.csv", index=False)

    summary = pd.DataFrame(rows).groupby(["calibration", "conformal"]).agg(
        coverage=("coverage", "mean"),
        worst_similarity_coverage=("sim_worst_coverage", "mean"),
        similarity_gap=("sim_max_coverage_gap", "mean"),
        mean_width=("mean_width", "mean"),
        fallback=("fallback_fraction", "mean"),
        n_calibration=("n_calibration", "mean"),
    ).reset_index()
    summary.to_csv(results / "table_calibration_study.csv", index=False)
    print(summary.to_string(index=False))
    print(f"done in {(time.time() - started) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
