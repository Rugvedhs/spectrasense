"""Turn raw stage results into the paper's tables and figures."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import _bootstrap  # noqa: F401
from ptgbench.figures import (
    ARM_COLORS,
    CONFORMAL_COLORS,
    REGIME_COLORS,
    SEQUENTIAL,
    SERIES,
    TEXT_SECONDARY,
    save,
    use_paper_style,
)

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
        rmse_mean=("rmse", "mean"), r2_mean=("r2", "mean"),
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


# --------------------------------------------------------------- figures ----
def fig_dataset(structures: pd.DataFrame, out: Path) -> None:
    order = structures.groupby("family")["Tg"].median().sort_values().index.tolist()
    data = [structures.loc[structures["family"] == f, "Tg"].to_numpy() for f in order]
    counts = [len(d) for d in data]

    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    parts = ax.boxplot(
        data, orientation="vertical", patch_artist=True, widths=0.62, showfliers=False,
        medianprops=dict(color="white", linewidth=1.2),
        whiskerprops=dict(color=TEXT_SECONDARY, linewidth=0.7),
        capprops=dict(color=TEXT_SECONDARY, linewidth=0.7),
    )
    for patch in parts["boxes"]:
        patch.set_facecolor(SERIES[0])
        patch.set_edgecolor("white")
        patch.set_linewidth(0.8)

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=38, ha="right")
    ax.set_ylabel("Glass transition temperature (°C)")
    ax.axhline(0, color=TEXT_SECONDARY, linewidth=0.6, linestyle=":")
    top = ax.get_ylim()[1]
    for i, n in enumerate(counts, start=1):
        ax.text(i, top, f"{n}", ha="center", va="bottom", fontsize=6, color=TEXT_SECONDARY)
    ax.set_ylim(top=top * 1.12)
    ax.set_title("Reported $T_g$ by repeat-unit family (structure level, n above each box)")
    save(fig, out / "fig1_dataset")


def fig_generalization_gap(point: pd.DataFrame, out: Path) -> None:
    frame = point[point["regime"].isin(REGIME_ORDER) & (point["model"] != "median")]
    models = [m for m in MODEL_LABELS if m in frame["model"].unique() and m != "median"]
    width = 0.8 / len(REGIME_ORDER)

    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    for j, regime in enumerate(REGIME_ORDER):
        means, errs = [], []
        for model in models:
            block = frame[(frame["regime"] == regime) & (frame["model"] == model)]["mae"]
            means.append(block.mean())
            errs.append(block.std() if len(block) > 1 else 0.0)
        positions = np.arange(len(models)) + j * width - 0.4 + width / 2
        bars = ax.bar(positions, means, width * 0.9, yerr=errs, capsize=2,
                      color=REGIME_COLORS[regime], label=regime.capitalize(),
                      edgecolor="white", linewidth=0.8,
                      error_kw=dict(ecolor=TEXT_SECONDARY, elinewidth=0.7))
        ax.bar_label(bars, fmt="%.0f", fontsize=5.5, padding=1, color=TEXT_SECONDARY)

    ax.set_xticks(np.arange(len(models)))
    ax.set_xticklabels([MODEL_LABELS[m] for m in models])
    ax.set_ylabel("Test MAE (°C)")
    ax.set_title("Prediction error rises as the split forces structural novelty")
    ax.legend(ncol=4, loc="upper left")
    save(fig, out / "fig2_generalization_gap")


def fig_matched(matched_raw: pd.DataFrame, statistics: pd.DataFrame, out: Path) -> None:
    if matched_raw.empty:
        return
    model = "extra_trees" if "extra_trees" in matched_raw["model"].unique() else matched_raw["model"].iloc[0]
    frame = matched_raw[matched_raw["model"] == model]
    per_family = frame.groupby("group")[["informed", "control", "naive"]].mean()
    per_family["family_effect"] = per_family["naive"] - per_family["control"]
    per_family = per_family.sort_values("family_effect")

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [1.15, 1]})

    ax = axes[0]
    y = np.arange(len(per_family))
    for arm in ("informed", "control", "naive"):
        ax.scatter(per_family[arm], y, s=26, color=ARM_COLORS[arm], label=arm.capitalize(),
                   edgecolor="white", linewidth=0.8, zorder=3)
    for i, (_, row) in enumerate(per_family.iterrows()):
        ax.plot([row["control"], row["naive"]], [i, i], color=TEXT_SECONDARY,
                linewidth=0.7, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(per_family.index, fontsize=6.5)
    ax.set_xlabel("MAE on the same held-out structures (°C)")
    ax.set_title("Matched arms, identical test set")
    ax.legend(ncol=3, loc="lower right")

    ax = axes[1]
    bars = ax.barh(y, per_family["family_effect"], color=SERIES[1], edgecolor="white",
                   linewidth=0.8, height=0.68)
    ax.bar_label(bars, fmt="%+.1f", fontsize=6, padding=2, color=TEXT_SECONDARY)
    ax.axvline(0, color=TEXT_SECONDARY, linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("MAE penalty from removing the family (°C)")
    ax.set_title("Family effect, size held constant")

    row = statistics[(statistics["model"] == model) & (statistics["effect"] == "family_effect")]
    if not row.empty:
        r = row.iloc[0]
        ax.text(0.98, 0.02,
                f"mean {r['mean']:+.1f} °C\n95% CI [{r['ci_low']:+.1f}, {r['ci_high']:+.1f}]\n"
                f"{int(r['n_families_positive'])}/{int(r['n_families'])} families worse",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6,
                color=TEXT_SECONDARY)
    save(fig, out / "fig3_matched_design")


def fig_coverage(conformal_table: pd.DataFrame, out: Path) -> None:
    regimes = [r for r in REGIME_ORDER if r in conformal_table["regime"].unique()]
    methods = [m for m in CONFORMAL_ORDER if m in conformal_table["conformal"].unique()]
    width = 0.8 / len(methods)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharex=True)
    for ax, column, title in (
        (axes[0], "coverage", "Marginal coverage: holds on random splits, fails under shift"),
        (axes[1], "worst_similarity_coverage",
         "Least-familiar band: fails even on random splits"),
    ):
        for j, method in enumerate(methods):
            values = [
                conformal_table.loc[
                    (conformal_table["regime"] == r) & (conformal_table["conformal"] == method),
                    column,
                ].mean()
                for r in regimes
            ]
            positions = np.arange(len(regimes)) + j * width - 0.4 + width / 2
            bars = ax.bar(positions, values, width * 0.9, color=CONFORMAL_COLORS[method],
                          label=method, edgecolor="white", linewidth=0.8)
            # Rotated: four near-equal bars per group collide with horizontal labels.
            ax.bar_label(bars, fmt="%.2f", fontsize=5.5, padding=2,
                         color=TEXT_SECONDARY, rotation=90)
        ax.axhline(0.9, color="#e34948", linewidth=1.0, linestyle="--")
        # Below the line and hard left, where no bar label can reach it.
        ax.text(-0.48, 0.893, "nominal 90%", fontsize=6, color="#e34948",
                ha="left", va="top")
        ax.set_xticks(np.arange(len(regimes)))
        ax.set_xticklabels([r.capitalize() for r in regimes])
        ax.set_ylim(0.3, 1.10)
        ax.set_title(title, fontsize=8)
    axes[0].set_ylabel("Empirical coverage")
    axes[0].legend(ncol=2, loc="lower left")
    save(fig, out / "fig4_coverage")


def fig_coverage_vs_similarity(ad: pd.DataFrame, out: Path) -> None:
    regimes = [r for r in REGIME_ORDER if r in ad["regime"].unique()]
    methods = [m for m in CONFORMAL_ORDER if m in ad["conformal"].unique()]

    fig, axes = plt.subplots(1, len(regimes), figsize=(7.4, 2.6), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, regime in zip(axes, regimes):
        for method in methods:
            block = ad[(ad["regime"] == regime) & (ad["conformal"] == method)]
            block = block.sort_values("similarity_bin")
            if block.empty:
                continue
            ax.plot(block["similarity_bin"], block["coverage"], marker="o",
                    color=CONFORMAL_COLORS[method], label=method,
                    markeredgecolor="white", markeredgewidth=0.7)
        ax.axhline(0.9, color="#e34948", linewidth=1.0, linestyle="--")
        ax.set_xticks(range(len(BIN_LABELS)))
        ax.set_xticklabels(BIN_LABELS, rotation=45, ha="right", fontsize=6)
        ax.set_title(regime.capitalize())
        ax.set_xlabel("Nearest-training Tanimoto")
    axes[0].set_ylabel("Coverage")
    axes[0].set_ylim(0.3, 1.02)
    axes[-1].legend(loc="lower right")
    fig.suptitle("Conditional coverage against structural novelty", y=1.04)
    save(fig, out / "fig5_coverage_vs_similarity")


def fig_width_tradeoff(conformal_table: pd.DataFrame, out: Path) -> None:
    methods = [m for m in CONFORMAL_ORDER if m in conformal_table["conformal"].unique()]
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    markers = {"random": "o", "scaffold": "s", "cluster": "^", "family": "D"}
    for method in methods:
        block = conformal_table[conformal_table["conformal"] == method]
        for _, row in block.iterrows():
            ax.scatter(row["mean_width"], row["worst_similarity_coverage"],
                       s=42, color=CONFORMAL_COLORS[method],
                       marker=markers.get(row["regime"], "o"),
                       edgecolor="white", linewidth=0.9, zorder=3)
    ax.axhline(0.9, color="#e34948", linewidth=1.0, linestyle="--")
    ax.set_xlabel("Mean interval width (°C)")
    ax.set_ylabel("Worst similarity-band coverage")
    ax.set_title("Wider is not automatically safer")
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=CONFORMAL_COLORS[m],
                          markeredgecolor="white", label=m) for m in methods]
    handles += [plt.Line2D([], [], marker=markers[r], linestyle="", color=TEXT_SECONDARY,
                           label=r.capitalize()) for r in markers if r in
                conformal_table["regime"].unique()]
    ax.legend(handles=handles, ncol=2, fontsize=6, loc="lower right")
    save(fig, out / "fig6_width_tradeoff")


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

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    for regime, colour in (("random", SERIES[0]), ("family", SERIES[3])):
        block = frame[frame["regime"] == regime]
        if block.empty:
            continue
        grouped = block.groupby("bin")["abs_error"].agg(["mean", "count", "sem"])
        ax.errorbar(grouped.index, grouped["mean"], yerr=grouped["sem"], marker="o",
                    color=colour, label=f"{regime.capitalize()} split", capsize=2,
                    markeredgecolor="white", markeredgewidth=0.7)
        for b, row in grouped.iterrows():
            ax.annotate(f"{int(row['count'])}", (b, row["mean"]), textcoords="offset points",
                        xytext=(0, 7), ha="center", fontsize=5.5, color=TEXT_SECONDARY)
    ax.set_xticks(range(len(BIN_LABELS)))
    ax.set_xticklabels(BIN_LABELS, rotation=45, ha="right")
    ax.set_xlabel("Nearest-training Tanimoto similarity")
    ax.set_ylabel("Mean absolute error (°C)")
    ax.set_title("Error grows as the nearest training analogue recedes")
    ax.legend()
    save(fig, out / "fig7_error_vs_similarity")


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

    ablation = table_representation_ablation(results, results)
    if not ablation.empty:
        print("\nrepresentation ablation (extra trees):")
        print(ablation.to_string(index=False))
    print(f"\nfigures written to {figures}")


if __name__ == "__main__":
    main()
