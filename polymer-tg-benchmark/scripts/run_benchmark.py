"""Run one benchmark stage over the cached features and write tidy result tables.

Stages are separate so a long run can be resumed, and so the expensive
matched-pair design can use a reduced model set without re-running the headline
comparison.
"""

from __future__ import annotations

import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd

import _bootstrap  # noqa: F401
from ptgbench import splits as S
from ptgbench.pipeline import run_split

MAIN_MODELS = ("median", "svr", "random_forest", "hist_gbr", "extra_trees")
MATCHED_MODELS = ("extra_trees", "hist_gbr")


def load_cache(processed_dir: str, representation: str = "descriptors"):
    """Load the cached table and the chosen model input.

    ``morgan`` swaps the RDKit descriptor block for the raw fingerprint bits.  The
    fingerprint matrix is always returned regardless, because clustering and the
    structural-novelty coordinate are defined in fingerprint space no matter what
    the regressor consumes.
    """
    processed = Path(processed_dir)
    table = pd.read_parquet(processed / "structures.parquet")
    matrix = np.load(processed / "fingerprint_matrix.npy")
    with open(processed / "fingerprints.pkl", "rb") as handle:
        fps = pickle.load(handle)

    if representation == "descriptors":
        X = pd.read_parquet(processed / "descriptors.parquet").to_numpy()
    elif representation == "morgan":
        X = matrix.astype(float)
    else:
        raise ValueError(f"unknown representation {representation!r}")
    return table, X, table["Tg"].to_numpy(), matrix, fps


def build_splits(stage: str, table, matrix, n_repeats: int, seed: int) -> list[S.Split]:
    smiles = table["canonical_psmiles"].tolist()
    families = table["family"].to_numpy()
    n = len(table)

    if stage == "random":
        return S.random_splits(n, n_repeats=n_repeats, seed=seed)
    if stage == "scaffold":
        return S.scaffold_splits(smiles, n_repeats=n_repeats, seed=seed)
    if stage == "cluster":
        clusters = S.fingerprint_clusters(matrix)
        return S.cluster_splits(clusters, n_repeats=n_repeats, seed=seed)
    if stage == "family":
        return S.family_holdout_splits(families)
    if stage == "matched":
        return S.matched_pair_splits(families, n_repeats=n_repeats, seed=seed)
    raise ValueError(f"unknown stage {stage!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True,
                        choices=["random", "scaffold", "cluster", "family", "matched"])
    parser.add_argument("--processed-dir", default="data/processed")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--n-repeats", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--representation", default="descriptors",
                        choices=["descriptors", "morgan"])
    parser.add_argument("--models", nargs="*", default=None,
                        help="override the model list for this stage")
    parser.add_argument("--suffix", default="",
                        help="appended to output filenames, to keep an ablation "
                             "run from overwriting the headline results")
    args = parser.parse_args()

    table, X, y, matrix, fps = load_cache(args.processed_dir, args.representation)
    families = table["family"].to_numpy()

    splits = build_splits(args.stage, table, matrix, args.n_repeats, args.seed)
    models = args.models or (
        MATCHED_MODELS if args.stage == "matched" else MAIN_MODELS
    )
    print(f"stage={args.stage}  splits={len(splits)}  models={len(models)}", flush=True)

    point, intervals, conditional, predictions = [], [], [], []
    started = time.time()
    total = len(splits) * len(models)
    done = 0

    for split in splits:
        for model_name in models:
            result = run_split(
                split, model_name, X, y, fps, families,
                alpha=args.alpha, seed=args.seed,
            )
            result.point["representation"] = args.representation
            point.append(result.point)
            if not result.intervals.empty:
                intervals.append(result.intervals)
                conditional.append(result.conditional)
            predictions.append(result.predictions)
            done += 1
            elapsed = time.time() - started
            print(
                f"[{done}/{total}] {split.name:42s} {model_name:14s} "
                f"mae={result.point['mae'].iloc[0]:7.2f} "
                f"eta={elapsed / done * (total - done) / 60:5.1f}min",
                flush=True,
            )

    results = Path(args.results_dir)
    results.mkdir(parents=True, exist_ok=True)
    pd.concat(point, ignore_index=True).to_csv(
        results / f"point_{args.stage}{args.suffix}.csv", index=False)
    if intervals:
        pd.concat(intervals, ignore_index=True).to_csv(
            results / f"intervals_{args.stage}{args.suffix}.csv", index=False)
        pd.concat(conditional, ignore_index=True).to_csv(
            results / f"conditional_{args.stage}{args.suffix}.csv", index=False)
    pd.concat(predictions, ignore_index=True).to_parquet(
        results / f"predictions_{args.stage}{args.suffix}.parquet")
    print(f"done in {(time.time() - started) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
