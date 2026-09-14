"""Turn raw stage results into the paper's tables and figures."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import _bootstrap  # noqa: F401
from ptgbench import figures as F
from ptgbench.figures import save, use_paper_style

ALPHA = 0.1
CWC_ETA = 30.0
REGIME_ORDER = ["random", "scaffold", "cluster", "family"]
MODEL_LABELS = {
    "median": "Median baseline",
    "svr": "SVR",
    "random_forest": "Random forest",
    "hist_gbr": "Hist. gradient boosting",
    "extra_trees": "Extra trees",
}
CONFORMAL_ORDER = ["SCP", "Normalised SCP", "Mondrian-family", "Mondrian-similarity"]
BIN_LABELS = ["<0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", ">=0.8"]


def read(results: Path, kind: str, stages: list[str]) -> pd.DataFrame:
    frames = []
    for stage in stages:
        path = results / f"{kind}_{stage}.csv"
        if path.exists():
            frames.append(pd.read_csv(path))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def bootstrap_ci(values: np.ndarray, n_boot: int = 10000, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    draws = rng.choice(values, size=(n_boot, values.size), replace=True).mean(axis=1)
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


# ---------------------------------------------------------------- tables ----
def table_point_by_regime(point: pd.DataFrame, out: Path) -> pd.DataFrame:
    frame = point[point["regime"].isin(REGIME_ORDER)]
    grouped = frame.groupby(["regime", "model"]).agg(
        n_splits=("mae", "size"),
        mae_mean=("mae", "mean"), mae_sd=("mae", "std"),
        rmse_mean=("rmse", "mean"), r2_mean_of_splits=("r2", "mean"),
        similarity=("mean_nn_similarity", "mean"),
    ).reset_index()
    grouped["regime"] = pd.Categorical(grouped["regime"], REGIME_ORDER, ordered=True)
    grouped = grouped.sort_values(["regime", "mae_mean"]).reset_index(drop=True)
    grouped.to_csv(out / "table_point_by_regime.csv", index=False)
    return grouped


def table_family_holdout(point: pd.DataFrame, out: Path) -> pd.DataFrame:
    family = point[point["regime"] == "family"].copy()
    if family.empty:
        return family
    reference = (
        point[point["regime"] == "random"].groupby("model")["mae"].mean().rename("mae_random")
    )
    family = family.merge(reference, on="model", how="left")
    family["deterioration"] = family["mae"] / family["mae_random"]
    keep = ["group", "model", "n", "mae", "rmse", "r2", "mae_random",
            "deterioration", "mean_nn_similarity"]
    family = family[keep].sort_values(["model", "deterioration"], ascending=[True, False])
    family.to_csv(out / "table_family_holdout.csv", index=False)
    return family


def table_matched(point_matched: pd.DataFrame, out: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if point_matched.empty:
        return point_matched, point_matched

    wide = point_matched.pivot_table(
        index=["group", "model", "split"], columns="arm", values="mae"
    ).reset_index()
    wide["split_rep"] = wide["split"].str.rsplit("_", n=1).str[0]
    wide = wide.groupby(["group", "model", "split_rep"], as_index=False)[
        ["informed", "control", "naive"]
    ].mean()

    wide["family_effect"] = wide["naive"] - wide["control"]
    wide["size_effect"] = wide["control"] - wide["informed"]
    wide["total_effect"] = wide["naive"] - wide["informed"]
    wide.to_csv(out / "table_matched_raw.csv", index=False)

    rows = []
    for (group, model), block in wide.groupby(["group", "model"]):
        row = {"group": group, "model": model, "n_repeats": len(block)}
        for col in ("informed", "control", "naive", "family_effect",
                    "size_effect", "total_effect"):
            row[f"{col}_mean"] = float(block[col].mean())
        rows.append(row)
    summary = pd.DataFrame(rows)

    stat_rows = []
    for model, block in wide.groupby("model"):
        per_family = block.groupby("group")[
            ["family_effect", "size_effect", "total_effect"]
        ].mean()
        for effect in ("family_effect", "size_effect", "total_effect"):
            values = per_family[effect].to_numpy()
            low, high = bootstrap_ci(values)
            wilcoxon = stats.wilcoxon(values) if len(values) > 5 else None
            stat_rows.append({
                "model": model, "effect": effect, "n_families": len(values),
                "mean": float(values.mean()), "ci_low": low, "ci_high": high,
                "median": float(np.median(values)),
                "wilcoxon_p": float(wilcoxon.pvalue) if wilcoxon else np.nan,
                "n_families_positive": int((values > 0).sum()),
            })
    statistics = pd.DataFrame(stat_rows)
    summary.to_csv(out / "table_matched_summary.csv", index=False)
    statistics.to_csv(out / "table_matched_statistics.csv", index=False)
    return summary, statistics


def table_conformal(intervals: pd.DataFrame, out: Path) -> pd.DataFrame:
    frame = intervals[intervals["model"] == "extra_trees"].copy()
    grouped = frame.groupby(["regime", "conformal"]).agg(
        n_splits=("coverage", "size"),
        coverage=("coverage", "mean"),
        worst_similarity_coverage=("sim_worst_coverage", "mean"),
        similarity_gap=("sim_max_coverage_gap", "mean"),
        worst_family_coverage=("fam_worst_coverage", "mean"),
        family_gap=("fam_max_coverage_gap", "mean"),
        mean_width=("mean_width", "mean"),
        fallback=("fallback_fraction", "mean"),
    ).reset_index()

    # CWC penalises undercoverage exponentially, so averaging per-split values
    # lets one badly-covered split dominate the mean by orders of magnitude.
    # Recomputing it from the aggregated coverage and width keeps the number
    # interpretable and comparable across regimes.
    deficit = (1 - ALPHA) - grouped["coverage"]
    grouped["cwc"] = grouped["mean_width"] * np.where(
        deficit > 0, np.exp(CWC_ETA * deficit), 1.0
    )
    grouped["regime"] = pd.Categorical(grouped["regime"], REGIME_ORDER, ordered=True)
    grouped = grouped.sort_values(["regime", "conformal"]).reset_index(drop=True)
    grouped.to_csv(out / "table_conformal.csv", index=False)
    return grouped


def table_applicability_domain(conditional: pd.DataFrame, out: Path) -> pd.DataFrame:
    frame = conditional[
        (conditional["taxonomy"] == "similarity") & (conditional["model"] == "extra_trees")
    ].copy()
    grouped = frame.groupby(["regime", "conformal", "group"]).apply(
        lambda b: pd.Series({
            "n": b["n"].sum(),
            "coverage": np.average(b["coverage"], weights=b["n"]),
            "mean_width": np.average(b["mean_width"], weights=b["n"]),
        }),
        include_groups=False,
    ).reset_index().rename(columns={"group": "similarity_bin"})
    grouped["similarity_label"] = [
        BIN_LABELS[int(b)] if int(b) < len(BIN_LABELS) else str(b)
        for b in grouped["similarity_bin"]
    ]
    grouped.to_csv(out / "table_applicability_domain.csv", index=False)
    return grouped


def table_gap_drivers(point: pd.DataFrame, structures: pd.DataFrame, out: Path) -> pd.DataFrame:
    """What predicts how badly a family transfers?

    Three candidate explanations are tested against the family-holdout
    deterioration: how structurally isolated the family is, how far its Tg
    distribution sits from the training median, and how large it is.  These are
    associations across twenty families, not a causal decomposition — the matched
    design is what supplies the causal claim.
    """
    family = point[point["regime"] == "family"]
    if family.empty:
        return pd.DataFrame()

    reference = point[point["regime"] == "random"].groupby("model")["mae"].mean()
    rows = []
    overall_median = structures["Tg"].median()
    for (group, model), block in family.groupby(["group", "model"]):
        if model == "median":
            continue
        members = structures[structures["family"] == group]["Tg"]
        others = structures[structures["family"] != group]["Tg"]
        rows.append({
            "group": group, "model": model,
            "deterioration": block["mae"].mean() / reference[model],
            "mae": block["mae"].mean(),
            "mean_nn_similarity": block["mean_nn_similarity"].mean(),
            "tg_median_shift": abs(members.median() - others.median()),
            "tg_iqr_ratio": (members.quantile(0.75) - members.quantile(0.25))
                            / (others.quantile(0.75) - others.quantile(0.25)),
            "n_members": len(members),
        })
    drivers = pd.DataFrame(rows)
    drivers.to_csv(out / "table_gap_drivers.csv", index=False)

    stat_rows = []
    for model, block in drivers.groupby("model"):
        for column in ("mean_nn_similarity", "tg_median_shift", "tg_iqr_ratio",
                       "n_members"):
            pearson = stats.pearsonr(block[column], block["deterioration"])
            spearman = stats.spearmanr(block[column], block["deterioration"])
            stat_rows.append({
                "model": model, "driver": column, "n_families": len(block),
                "pearson_r": float(pearson.statistic),
                "pearson_p": float(pearson.pvalue),
                "spearman_rho": float(spearman.statistic),
                "spearman_p": float(spearman.pvalue),
            })
    correlations = pd.DataFrame(stat_rows)
    correlations.to_csv(out / "table_gap_driver_correlations.csv", index=False)
    return correlations


def table_representation_ablation(results: Path, out: Path) -> pd.DataFrame:
    """Does the generalisation gap depend on the molecular representation?

    If swapping physically interpretable descriptors for raw Morgan bits left the
    family-holdout penalty unchanged, the penalty is a property of the data rather
    than of the feature set — which is the claim the paper needs to defend against
    "a better representation would fix it".
    """
    frames = []
    for stage in ("random", "family"):
        for suffix, label in (("", "descriptors"), ("_morgan", "morgan")):
            path = results / f"point_{stage}{suffix}.csv"
            if path.exists():
                block = pd.read_csv(path)
                block["representation"] = block.get("representation", label)
                block["representation"] = block["representation"].fillna(label)
                frames.append(block[block["model"] == "extra_trees"])
    if len(frames) < 3:
        return pd.DataFrame()

    frame = pd.concat(frames, ignore_index=True)
    summary = frame.groupby(["representation", "regime"]).agg(
        n_splits=("mae", "size"), mae=("mae", "mean"), mae_sd=("mae", "std"),
        r2=("r2", "mean"),
    ).reset_index()

    rows = []
    for representation, block in summary.groupby("representation"):
        by_regime = block.set_index("regime")["mae"]
        if {"random", "family"}.issubset(by_regime.index):
            rows.append({
                "representation": representation,
                "mae_random": by_regime["random"],
                "mae_family": by_regime["family"],
                "penalty_K": by_regime["family"] - by_regime["random"],
                "penalty_ratio": by_regime["family"] / by_regime["random"],
            })
    ablation = pd.DataFrame(rows)
    summary.to_csv(out / "table_representation_summary.csv", index=False)
    ablation.to_csv(out / "table_representation_ablation.csv", index=False)
    return ablation


def table_baseline_comparison(results: Path, out: Path) -> pd.DataFrame:
    """Additive group contribution against the machine-learned ensembles.

    The question is not only which is more accurate on a random split — the
    ensembles are — but which degrades less when the split forces novel backbone
    chemistry. An additive law is compositional, so it has a prior claim to
    transferring: a group it has seen in one family carries the same contribution
    in another. Whether that claim survives contact with the data is the point.
    """
    frames = []
    for stage in REGIME_ORDER:
        for suffix in ("", "_gc"):
            path = results / f"point_{stage}{suffix}.csv"
            if path.exists():
                frames.append(pd.read_csv(path))
    if not frames:
        return pd.DataFrame()

    frame = pd.concat(frames, ignore_index=True)
    frame = frame[frame["model"].isin(["group_contribution", "hist_gbr", "median"])]
    summary = frame.groupby(["model", "regime"]).agg(
        n_splits=("mae", "size"), mae=("mae", "mean"), mae_sd=("mae", "std"),
        r2=("r2", "mean"),
    ).reset_index()
    if summary.empty:
        return summary

    wide = summary.pivot_table(index="model", columns="regime", values="mae")
    available = [r for r in REGIME_ORDER if r in wide.columns]
    wide = wide[available]
    if {"random", "family"}.issubset(wide.columns):
        wide["penalty_K"] = wide["family"] - wide["random"]
        wide["penalty_ratio"] = wide["family"] / wide["random"]
    wide = wide.reset_index()

    summary.to_csv(out / "table_baseline_by_regime.csv", index=False)
    wide.to_csv(out / "table_baseline_comparison.csv", index=False)
    return wide


def table_pooled_r2(results: Path, out: Path) -> pd.DataFrame:
    """Pooled R2 alongside the mean of per-split R2, because they differ.

    Averaging R2 over family holdouts normalises each family by its own Tg
    variance, so a narrow-spread family can contribute a large negative value at
    a modest absolute error. The mean of those is not comparable to a
    random-split R2 and must not be quoted as though it were. Pooling the
    residuals over all held-out structures gives the comparable quantity.
    """
    rows = []
    for stage in REGIME_ORDER:
        path = results / f"predictions_{stage}.parquet"
        point_path = results / f"point_{stage}.csv"
        if not (path.exists() and point_path.exists()):
            continue
        predictions = pd.read_parquet(path)
        point = pd.read_csv(point_path)
        for model in predictions["model"].unique():
            block = predictions[predictions["model"] == model]
            y = block["y_true"].to_numpy()
            residual = y - block["y_pred"].to_numpy()
            denominator = ((y - y.mean()) ** 2).sum()
            rows.append({
                "regime": stage, "model": model,
                "n_predictions": len(block),
                "r2_pooled": 1 - (residual ** 2).sum() / denominator
                if denominator > 0 else np.nan,
                "r2_mean_of_splits": point[point["model"] == model]["r2"].mean(),
                "mae_pooled": float(np.abs(residual).mean()),
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "table_pooled_r2.csv", index=False)
    return frame


def table_band_mae(results: Path, out: Path) -> pd.DataFrame:
    """Absolute error resolved by nearest-training similarity band.

    The applicability-domain argument needs this: an interval predictor that
    applies one width across bands whose error differs by a factor of X is
    mis-specified by exactly that factor, and X has to be measured rather than
    asserted.
    """
    rows = []
    for stage in REGIME_ORDER:
        path = results / f"predictions_{stage}.parquet"
        if not path.exists():
            continue
        predictions = pd.read_parquet(path)
        predictions = predictions.assign(
            abs_error=(predictions["y_true"] - predictions["y_pred"]).abs(),
            band=np.digitize(predictions["nn_similarity"], [0.4, 0.5, 0.6, 0.7, 0.8]),
        )
        for (model, band), block in predictions.groupby(["model", "band"]):
            rows.append({
                "regime": stage, "model": model, "band": band,
                "band_label": BIN_LABELS[int(band)] if int(band) < len(BIN_LABELS)
                else str(band),
                "n": len(block), "mae": float(block["abs_error"].mean()),
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "table_band_mae.csv", index=False)

    ratios = []
    for (regime, model), block in frame.groupby(["regime", "model"]):
        reliable = block[block["n"] >= 30]
        if len(reliable) > 1:
            ratios.append({
                "regime": regime, "model": model,
                "mae_worst_band": reliable["mae"].max(),
                "mae_best_band": reliable["mae"].min(),
                "ratio": reliable["mae"].max() / reliable["mae"].min(),
            })
    ratio_frame = pd.DataFrame(ratios)
    ratio_frame.to_csv(out / "table_band_mae_ratio.csv", index=False)
    return ratio_frame


def table_model_comparison(results: Path, out: Path, reference: str = "hist_gbr") -> pd.DataFrame:
    """Paired tests between models on the repeated random splits.

    Both a paired t-test and a Wilcoxon signed-rank test are reported, because
    for the closest pair they disagree about significance at 0.05 and the
    manuscript's claim rests on which is believed. Repeated splits of one dataset
    are not independent draws, so the t-test is anti-conservative here; the
    honest reading is that the two leading models are not separated by an amount
    this design can resolve.
    """
    path = results / "point_random.csv"
    if not path.exists():
        return pd.DataFrame()
    point = pd.read_csv(path)
    if reference not in set(point["model"]):
        return pd.DataFrame()
    base = point[point["model"] == reference].sort_values("split")["mae"].to_numpy()

    rows = []
    for model in point["model"].unique():
        if model == reference:
            continue
        other = point[point["model"] == model].sort_values("split")["mae"].to_numpy()
        if len(other) != len(base):
            continue
        rows.append({
            "model": model, "reference": reference, "n_splits": len(base),
            "mae": float(other.mean()), "reference_mae": float(base.mean()),
            "delta_mae": float(other.mean() - base.mean()),
            "paired_t_p": float(stats.ttest_rel(other, base).pvalue),
            "wilcoxon_p": float(stats.wilcoxon(other, base).pvalue),
        })
    frame = pd.DataFrame(rows).sort_values("delta_mae").reset_index(drop=True)
    frame.to_csv(out / "table_model_comparison.csv", index=False)
    return frame


def table_width_matched_control(results: Path, out: Path,
                                model: str = "extra_trees") -> pd.DataFrame:
    """Can a uniformly inflated split-conformal interval do the same job?

    The novelty-conditioned predictor buys its conditional coverage partly with
    width, so the fair control is split conformal scaled by a single constant to
    the same mean width. That control is an *oracle*: the constant is computed
    from the novelty-conditioned widths on the test set, and nothing available at
    prediction time supplies it. It is reported anyway, because the honest
    comparison is against the strongest strawman rather than the most convenient
    one.

    A McNemar test on the paired coverage indicators settles the band the
    argument rests on; the two predictors see identical structures, so the
    discordant pairs are the whole of the evidence.
    """
    rows = []
    for regime in ("family", "cluster"):
        path = results / f"predictions_{regime}.parquet"
        if not path.exists():
            continue
        frame = pd.read_parquet(path)
        frame = frame[frame["model"] == model]
        if frame.empty or "lower_Mondrian-similarity" not in frame.columns:
            continue

        y = frame["y_true"].to_numpy()
        scp_width = (frame["upper_SCP"] - frame["lower_SCP"]).to_numpy()
        centre = ((frame["upper_SCP"] + frame["lower_SCP"]) / 2).to_numpy()
        mondrian = (
            (y >= frame["lower_Mondrian-similarity"])
            & (y <= frame["upper_Mondrian-similarity"])
        ).to_numpy()
        mondrian_width = (
            frame["upper_Mondrian-similarity"] - frame["lower_Mondrian-similarity"]
        ).to_numpy()
        band = np.digitize(frame["nn_similarity"].to_numpy(), [0.4, 0.5, 0.6, 0.7, 0.8])
        lowest = band == 0

        # Sweep the inflation factor, including the oracle value that equalises
        # mean width and the factor that would reach nominal in the lowest band.
        oracle_k = mondrian_width.mean() / scp_width.mean()
        for k in sorted({round(oracle_k, 4), 1.0, 1.2, 1.4, 1.6, 1.8, 1.93, 2.0, 2.2}):
            half = scp_width / 2 * k
            covered = (y >= centre - half) & (y <= centre + half)
            rows.append({
                "regime": regime, "model": model, "k": k,
                "is_oracle_k": abs(k - round(oracle_k, 4)) < 1e-9,
                "coverage_pooled": float(covered.mean()),
                "coverage_lowest_band": float(covered[lowest].mean()),
                "mean_width": float((2 * half).mean()),
                "mondrian_coverage_pooled": float(mondrian.mean()),
                "mondrian_coverage_lowest_band": float(mondrian[lowest].mean()),
                "mondrian_mean_width": float(mondrian_width.mean()),
                "n_lowest_band": int(lowest.sum()),
            })

        # McNemar at the oracle factor, in the band the claim is about.
        half = scp_width / 2 * oracle_k
        matched = (y >= centre - half) & (y <= centre + half)
        only_mondrian = int((mondrian[lowest] & ~matched[lowest]).sum())
        only_matched = int((~mondrian[lowest] & matched[lowest]).sum())
        discordant = only_mondrian + only_matched
        # Attach to the oracle row, not to whichever k the sweep ended on.
        oracle_row = next(
            r for r in rows if r["regime"] == regime and r["is_oracle_k"]
        )
        oracle_row["mcnemar_only_mondrian"] = only_mondrian
        oracle_row["mcnemar_only_width_matched"] = only_matched
        oracle_row["mcnemar_p"] = (
            float(stats.binomtest(only_mondrian, discordant, 0.5).pvalue)
            if discordant else np.nan
        )

    frame = pd.DataFrame(rows)
    frame.to_csv(out / "table_width_matched.csv", index=False)
    return frame


def table_bias_shrinkage(results: Path, structures: pd.DataFrame, out: Path,
                         model: str = "hist_gbr") -> pd.DataFrame:
    """Is the family-resolved bias chemistry, or regression to the training mean?

    A model fitted mostly on other families pulls its predictions toward the bulk
    of the training distribution, so any family whose Tg sits far from that bulk
    is biased toward it for reasons that have nothing to do with its backbone.
    That alternative has to be removed before a signed error can be read as
    evidence of a chemical mechanism.

    Regressing bias on the family's Tg offset separates the two. The families
    that move *against* the shrinkage line are the ones whose bias the offset
    cannot explain, and they are where the chemical argument actually lives.
    """
    path = results / "point_family.csv"
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    frame = frame[(frame["model"] == model)].dropna(subset=["group"]).copy()
    if frame.empty:
        return frame

    offsets = []
    for family in frame["group"]:
        member = structures["family"] == family
        offsets.append(
            structures.loc[member, "Tg"].median()
            - structures.loc[~member, "Tg"].median()
        )
    frame["tg_offset"] = offsets

    slope, intercept = np.polyfit(frame["tg_offset"], frame["bias"], 1)
    frame["bias_expected_from_offset"] = slope * frame["tg_offset"] + intercept
    frame["bias_residual"] = frame["bias"] - frame["bias_expected_from_offset"]

    correlation = stats.pearsonr(frame["tg_offset"], frame["bias"])
    frame["shrinkage_r"] = correlation.statistic
    frame["shrinkage_p"] = correlation.pvalue

    frame = frame[[
        "group", "n", "mae", "bias", "tg_offset", "bias_expected_from_offset",
        "bias_residual", "shrinkage_r", "shrinkage_p",
    ]].sort_values("bias_residual", ascending=False)
    frame.to_csv(out / "table_bias_shrinkage.csv", index=False)
    return frame


def table_matched_by_analogue(results: Path, out: Path) -> pd.DataFrame:
    """Does the matched effect survive where near-duplicates were never available?

    The informed arm keeps the rest of the family, so for some test structures it
    holds a very close analogue that the naive arm cannot have. If the family
    effect were only the loss of near-duplicates it would vanish once those
    structures are excluded. Splitting the effect by how close the informed arm's
    nearest training neighbour was tests exactly that.
    """
    path = results / "predictions_matched.parquet"
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_parquet(path)
    if "arm" not in frame.columns:
        return pd.DataFrame()

    informed = frame[frame["arm"] == "informed"][
        ["split", "model", "index", "nn_similarity"]
    ].rename(columns={"nn_similarity": "informed_similarity"})
    informed["pair"] = informed["split"].str.rsplit("_", n=1).str[0]

    merged = frame.copy()
    merged["pair"] = merged["split"].str.rsplit("_", n=1).str[0]
    merged = merged.merge(
        informed[["pair", "model", "index", "informed_similarity"]],
        on=["pair", "model", "index"], how="left",
    )
    merged["abs_error"] = (merged["y_true"] - merged["y_pred"]).abs()
    merged["analogue_band"] = pd.cut(
        merged["informed_similarity"], [0, 0.6, 0.8, 0.9, 1.01],
        labels=["<0.6", "0.6-0.8", "0.8-0.9", ">=0.9"],
    )

    rows = []
    for (model, band), block in merged.groupby(["model", "analogue_band"],
                                               observed=True):
        wide = block.pivot_table(index=["pair", "index"], columns="arm",
                                 values="abs_error")
        if not {"naive", "control"}.issubset(wide.columns):
            continue
        wide = wide.dropna(subset=["naive", "control"])
        if wide.empty:
            continue
        rows.append({
            "model": model, "analogue_band": str(band), "n": len(wide),
            "mae_informed": float(wide["informed"].mean())
            if "informed" in wide else np.nan,
            "mae_control": float(wide["control"].mean()),
            "mae_naive": float(wide["naive"].mean()),
            "family_effect": float((wide["naive"] - wide["control"]).mean()),
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "table_matched_by_analogue.csv", index=False)
    return frame


# --------------------------------------------------------------- figures ----
def fig_dataset(structures: pd.DataFrame, out: Path) -> None:
    """Family Tg distributions, horizontal so family names need no rotation."""
    order = structures.groupby("family")["Tg"].median().sort_values().index.tolist()
    data = [structures.loc[structures["family"] == f, "Tg"].to_numpy() for f in order]

    fig, ax = plt.subplots(figsize=(F.ONE_HALF_COLUMN, 4.1))
    parts = ax.boxplot(
        data, orientation="horizontal", patch_artist=True, widths=0.6,
        showfliers=False, showcaps=False,
        medianprops=dict(color=F.INK, linewidth=1.0),
        whiskerprops=dict(color=F.HAIRLINE, linewidth=0.5),
        boxprops=dict(linewidth=0.0),
    )
    for patch in parts["boxes"]:
        patch.set_facecolor(F.FILL)
        patch.set_edgecolor("none")

    ax.set_yticks(range(1, len(order) + 1))
    ax.set_yticklabels(order)
    ax.set_xlabel("Glass transition temperature (°C)")
    F.reference_line(ax, 0, orientation="v")
    F.yardstick(ax, axis="x")
    ax.set_ylim(0.3, len(order) + 0.7)

    # Counts as a right-hand column rather than floating annotations.
    right = ax.get_xlim()[1]
    ax.text(right, len(order) + 0.9, "n", fontsize=6.5, color=F.MUTED,
            ha="right", va="bottom", fontstyle="italic")
    for i, family in enumerate(order, start=1):
        ax.text(right, i, f"{len(structures[structures.family == family])}",
                fontsize=6.3, color=F.MUTED, ha="right", va="center")
    F.save(fig, out / "fig1_dataset")


def fig_generalization_gap(point: pd.DataFrame, out: Path) -> None:
    """Error against split regime: a trend, so a line rather than grouped bars."""
    frame = point[point["regime"].isin(REGIME_ORDER) & (point["model"] != "median")]
    models = [m for m in MODEL_LABELS if m in set(frame["model"]) and m != "median"]
    x = np.arange(len(REGIME_ORDER))

    fig, ax = plt.subplots(figsize=(F.ONE_HALF_COLUMN, 2.6))
    # Dodge each series slightly: undodged error bars sit on top of one another
    # at the wide-spread regimes and read as hatching rather than as uncertainty.
    offsets = np.linspace(-0.055, 0.055, len(models))
    finals = []
    for model, offset in zip(models, offsets):
        means, errs = [], []
        for regime in REGIME_ORDER:
            block = frame[(frame["regime"] == regime) & (frame["model"] == model)]["mae"]
            means.append(block.mean())
            errs.append(block.std() if len(block) > 1 else 0.0)
        colour = F.MODEL_COLORS[model]
        ax.errorbar(x + offset, means, yerr=errs, marker="o", color=colour,
                    capsize=0, elinewidth=0.5, alpha=0.95,
                    markeredgecolor="white", markeredgewidth=0.4, zorder=3)
        finals.append(means[-1])

    for model, y_label, y_data, offset in zip(
        models, F.spread_labels(finals, minimum_gap=2.0), finals, offsets
    ):
        colour = F.MODEL_COLORS[model]
        if abs(y_label - y_data) > 0.05:
            ax.plot([x[-1] + offset, x[-1] + 0.20], [y_data, y_label],
                    color=colour, linewidth=0.5, zorder=2)
        F.end_label(ax, x[-1], y_label, MODEL_LABELS[model], colour, dx=0.22)

    ax.set_xticks(x)
    ax.set_xticklabels([r.capitalize() for r in REGIME_ORDER])
    ax.set_xlim(-0.3, len(REGIME_ORDER) - 1 + 1.55)
    ax.set_ylabel("Test MAE (K)")
    ax.set_xlabel("Partitioning regime, in order of enforced structural novelty")
    F.yardstick(ax)
    F.save(fig, out / "fig2_generalization_gap")


def fig_matched(matched_raw: pd.DataFrame, statistics: pd.DataFrame, out: Path) -> None:
    if matched_raw.empty:
        return
    model = ("extra_trees" if "extra_trees" in set(matched_raw["model"])
             else matched_raw["model"].iloc[0])
    frame = matched_raw[matched_raw["model"] == model]
    per_family = frame.groupby("group")[["informed", "control", "naive"]].mean()
    per_family["family_effect"] = per_family["naive"] - per_family["control"]
    per_family = per_family.sort_values("family_effect")
    y = np.arange(len(per_family))

    fig, axes = plt.subplots(1, 2, figsize=(F.DOUBLE_COLUMN, 3.5),
                             gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.06})

    ax = axes[0]
    for i, (_, row) in enumerate(per_family.iterrows()):
        ax.plot([row["control"], row["naive"]], [i, i], color=F.HAIRLINE,
                linewidth=0.6, zorder=2, solid_capstyle="butt")
    for arm in ("informed", "control", "naive"):
        ax.scatter(per_family[arm], y, s=14, color=F.ARM_COLORS[arm], label=arm,
                   edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(per_family.index)
    ax.set_xlabel("MAE on the same held-out structures (K)")
    ax.set_ylim(-0.8, len(per_family) - 0.2)
    F.yardstick(ax, axis="x")
    ax.legend(loc="lower right", ncol=1)
    F.panel(ax, "a")

    ax = axes[1]
    ax.barh(y, per_family["family_effect"], color=F.VERMILLION, height=0.55,
            linewidth=0)
    ax.axvline(0, color=F.HAIRLINE, linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("MAE penalty from removing the family (K)")
    ax.set_ylim(-0.8, len(per_family) - 0.2)
    F.yardstick(ax, axis="x")
    F.panel(ax, "b")

    # Label only the extremes: the rest are readable off the axis.
    for i in (0, len(per_family) - 1):
        value = per_family["family_effect"].iloc[i]
        ax.text(value + 0.9, i, f"{value:+.1f}", fontsize=6.3, color=F.MUTED,
                va="center")

    row = statistics[(statistics["model"] == model)
                     & (statistics["effect"] == "family_effect")]
    if not row.empty:
        r = row.iloc[0]
        ax.text(0.97, 0.06,
                f"mean {r['mean']:+.1f} K  (95% CI {r['ci_low']:+.1f} to "
                f"{r['ci_high']:+.1f})\n{int(r['n_families_positive'])} of "
                f"{int(r['n_families'])} families worse",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6.3,
                color=F.INK)
    F.save(fig, out / "fig3_matched_design")


# Line style is a second channel alongside hue: under family holdout two methods
# coincide exactly, and identical solid lines would hide that rather than show it.
_CONFORMAL_DASH = {
    "SCP": (0, ()),
    "Normalised SCP": (0, (4, 2)),
    "Mondrian-family": (0, (1.2, 1.6)),
    "Mondrian-similarity": (0, ()),
}
_CONFORMAL_WIDTH = {"Mondrian-similarity": 1.6}


def _conformal_line(ax, x, values, method, marker="o"):
    ax.plot(x, values, marker=marker, color=F.CONFORMAL_COLORS[method],
            linestyle=_CONFORMAL_DASH[method],
            linewidth=_CONFORMAL_WIDTH.get(method, 1.1),
            markeredgecolor="white", markeredgewidth=0.4,
            label=method, zorder=3)


def fig_coverage(conformal_table: pd.DataFrame, out: Path) -> None:
    regimes = [r for r in REGIME_ORDER if r in set(conformal_table["regime"])]
    methods = [m for m in CONFORMAL_ORDER if m in set(conformal_table["conformal"])]
    x = np.arange(len(regimes))

    fig, axes = plt.subplots(1, 2, figsize=(F.DOUBLE_COLUMN, 2.7), sharey=True,
                             gridspec_kw={"wspace": 0.08})
    for ax, column, letter in ((axes[0], "coverage", "a"),
                               (axes[1], "worst_similarity_coverage", "b")):
        for method in methods:
            values = [
                conformal_table.loc[(conformal_table["regime"] == r)
                                    & (conformal_table["conformal"] == method),
                                    column].mean()
                for r in regimes
            ]
            _conformal_line(ax, x, values, method)
        F.reference_line(ax, 0.9, "nominal 0.90")
        ax.set_xticks(x)
        ax.set_xticklabels([r.capitalize() for r in regimes])
        ax.set_xlim(-0.3, len(regimes) - 0.7)
        ax.set_ylim(0.45, 1.0)
        F.yardstick(ax)
        F.panel(ax, letter)

    axes[0].set_ylabel("Empirical coverage")
    axes[0].set_xlabel("Partitioning regime")
    axes[1].set_xlabel("Partitioning regime")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.10))
    F.save(fig, out / "fig4_coverage")


def fig_coverage_vs_similarity(ad: pd.DataFrame, out: Path) -> None:
    regimes = [r for r in REGIME_ORDER if r in set(ad["regime"])]
    methods = [m for m in CONFORMAL_ORDER if m in set(ad["conformal"])]

    fig, axes = plt.subplots(1, len(regimes), figsize=(F.DOUBLE_COLUMN, 2.3),
                             sharey=True, gridspec_kw={"wspace": 0.08})
    axes = np.atleast_1d(axes)
    for ax, regime, letter in zip(axes, regimes, "abcd"):
        for method in methods:
            block = ad[(ad["regime"] == regime)
                       & (ad["conformal"] == method)].sort_values("similarity_bin")
            if block.empty:
                continue
            _conformal_line(ax, block["similarity_bin"], block["coverage"], method)
        F.reference_line(ax, 0.9)
        ax.set_xticks(range(len(BIN_LABELS)))
        ax.set_xticklabels(BIN_LABELS, rotation=90, fontsize=6)
        ax.set_ylim(0.45, 1.0)
        F.yardstick(ax)
        F.panel(ax, letter)
        ax.text(0.5, 1.04, regime.capitalize(), transform=ax.transAxes,
                fontsize=7, color=F.INK, ha="center", va="bottom")

    axes[0].set_ylabel("Coverage")
    fig.text(0.5, -0.16, "Nearest-training Tanimoto similarity", ha="center",
             fontsize=7.5)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.30))
    F.save(fig, out / "fig5_coverage_vs_similarity")


def fig_width_tradeoff(conformal_table: pd.DataFrame, out: Path) -> None:
    methods = [m for m in CONFORMAL_ORDER if m in set(conformal_table["conformal"])]
    markers = {"random": "o", "scaffold": "s", "cluster": "^", "family": "D"}

    fig, ax = plt.subplots(figsize=(F.SINGLE_COLUMN, 2.7))
    for method in methods:
        for _, row in conformal_table[conformal_table["conformal"] == method].iterrows():
            ax.scatter(row["mean_width"], row["worst_similarity_coverage"], s=20,
                       color=F.CONFORMAL_COLORS[method],
                       marker=markers.get(row["regime"], "o"),
                       edgecolor="white", linewidth=0.4, zorder=3)
    F.reference_line(ax, 0.9, "nominal 0.90")
    ax.set_xlabel("Mean interval width (K)")
    ax.set_ylabel("Worst similarity-band coverage")
    F.yardstick(ax)

    method_handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=3.2,
                                 color=F.CONFORMAL_COLORS[m], label=m)
                      for m in methods]
    regime_handles = [plt.Line2D([], [], marker=markers[r], linestyle="",
                                 markersize=3.2, color=F.MUTED, label=r.capitalize())
                      for r in markers if r in set(conformal_table["regime"])]
    ax.legend(handles=method_handles + regime_handles, loc="lower right",
              fontsize=6, ncol=1)
    F.save(fig, out / "fig6_width_tradeoff")


def fig_error_vs_similarity(predictions_dir: Path, out: Path) -> None:
    frames = []
    for stage in ("random", "family"):
        path = predictions_dir / f"predictions_{stage}.parquet"
        if path.exists():
            frames.append(pd.read_parquet(path))
    if not frames:
        return
    frame = pd.concat(frames, ignore_index=True)
    frame = frame[frame["model"] == "extra_trees"].copy()
    frame["abs_error"] = (frame["y_true"] - frame["y_pred"]).abs()
    frame["bin"] = np.digitize(frame["nn_similarity"], [0.4, 0.5, 0.6, 0.7, 0.8])

    fig, ax = plt.subplots(figsize=(F.SINGLE_COLUMN, 2.5))
    for regime in ("random", "family"):
        block = frame[frame["regime"] == regime]
        if block.empty:
            continue
        grouped = block.groupby("bin")["abs_error"].agg(["mean", "sem"])
        colour = F.REGIME_COLORS[regime]
        ax.errorbar(grouped.index, grouped["mean"], yerr=grouped["sem"], marker="o",
                    color=colour, capsize=1.6, elinewidth=0.5, capthick=0.5,
                    markeredgecolor="white", markeredgewidth=0.4,
                    label=f"{regime.capitalize()} split", zorder=3)
    ax.set_xticks(range(len(BIN_LABELS)))
    ax.set_xticklabels(BIN_LABELS, rotation=90)
    ax.set_xlabel("Nearest-training Tanimoto similarity")
    ax.set_ylabel("Mean absolute error (K)")
    F.yardstick(ax)
    ax.legend(loc="upper right")
    F.save(fig, out / "fig7_error_vs_similarity")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--processed-dir", default="data/processed")
    parser.add_argument("--figures-dir", default="paper/figures")
    args = parser.parse_args()

    use_paper_style()
    results = Path(args.results_dir)
    figures = Path(args.figures_dir)
    figures.mkdir(parents=True, exist_ok=True)

    stages = REGIME_ORDER
    point = read(results, "point", stages)
    intervals = read(results, "intervals", stages)
    conditional = read(results, "conditional", stages)
    point_matched = read(results, "point", ["matched"])
    structures = pd.read_parquet(Path(args.processed_dir) / "structures.parquet")

    if not point.empty:
        summary = table_point_by_regime(point, results)
        print(summary.to_string(index=False))
        table_family_holdout(point, results)
        correlations = table_gap_drivers(point, structures, results)
        if not correlations.empty:
            print(correlations.to_string(index=False))
        fig_generalization_gap(point, figures)
    fig_dataset(structures, figures)

    if not intervals.empty:
        conformal_table = table_conformal(intervals, results)
        print(conformal_table.to_string(index=False))
        fig_coverage(conformal_table, figures)
        fig_width_tradeoff(conformal_table, figures)
    if not conditional.empty:
        ad = table_applicability_domain(conditional, results)
        fig_coverage_vs_similarity(ad, figures)
    if not point_matched.empty:
        _, statistics = table_matched(point_matched, results)
        raw = pd.read_csv(results / "table_matched_raw.csv")
        print(statistics.to_string(index=False))
        fig_matched(raw, statistics, figures)
    fig_error_vs_similarity(results, figures)

    width_matched = table_width_matched_control(results, results)
    if not width_matched.empty:
        print("\nwidth-matched oracle control:")
        print(width_matched[width_matched["is_oracle_k"]][
            ["regime", "k", "coverage_pooled", "coverage_lowest_band",
             "mondrian_coverage_lowest_band", "mean_width", "mcnemar_p"]
        ].to_string(index=False))

    shrinkage = table_bias_shrinkage(results, structures, results)
    if not shrinkage.empty:
        print(f"\nbias vs family Tg offset: r={shrinkage['shrinkage_r'].iloc[0]:+.3f} "
              f"p={shrinkage['shrinkage_p'].iloc[0]:.4f}")
        print(shrinkage.head(4)[["group", "tg_offset", "bias",
                                 "bias_residual"]].to_string(index=False))

    analogue = table_matched_by_analogue(results, results)
    if not analogue.empty:
        print("\nmatched family effect by informed-arm analogue closeness:")
        print(analogue.to_string(index=False))

    pooled = table_pooled_r2(results, results)
    if not pooled.empty:
        print("\npooled vs mean-of-split R2:")
        print(pooled[pooled.model == "hist_gbr"][
            ["regime", "r2_pooled", "r2_mean_of_splits"]].to_string(index=False))

    bands = table_band_mae(results, results)
    if not bands.empty:
        print("\nerror ratio across similarity bands:")
        print(bands[bands.model == "hist_gbr"].to_string(index=False))

    comparison = table_model_comparison(results, results)
    if not comparison.empty:
        print("\nmodel comparison on repeated random splits:")
        print(comparison.to_string(index=False))

    baseline = table_baseline_comparison(results, results)
    if not baseline.empty:
        print("\nadditive baseline vs ensemble (MAE by regime):")
        print(baseline.to_string(index=False))

    ablation = table_representation_ablation(results, results)
    if not ablation.empty:
        print("\nrepresentation ablation (extra trees):")
        print(ablation.to_string(index=False))
    print(f"\nfigures written to {figures}")


if __name__ == "__main__":
    main()
