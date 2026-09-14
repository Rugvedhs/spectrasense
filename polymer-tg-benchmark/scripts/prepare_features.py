"""Compute and cache the descriptor block, fingerprints and structure table.

Descriptor computation dominates the cost of a cold run, so it is done once here
and every experiment stage reads the cache.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

import _bootstrap  # noqa: F401  (puts src/ on the path)
from ptgbench.data import family_summary, load_structure_level
from ptgbench.features import compute_descriptors, compute_fingerprints
from ptgbench.groupcontrib import compute_group_counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--out-dir", default="data/processed")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--fp-radius", type=int, default=2)
    parser.add_argument("--fp-bits", type=int, default=2048)
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results = Path(args.results_dir)
    results.mkdir(parents=True, exist_ok=True)

    table, report = load_structure_level(args.raw_dir)
    print(report.to_frame().to_string(index=False))
    print("attachment points:", report.attachment_point_counts)

    smiles = table["canonical_psmiles"].tolist()
    descriptors = compute_descriptors(smiles)
    matrix, fps = compute_fingerprints(smiles, args.fp_radius, args.fp_bits)
    groups, masses = compute_group_counts(smiles)
    groups["molar_mass"] = masses

    table.to_parquet(out / "structures.parquet")
    descriptors.to_parquet(out / "descriptors.parquet")
    groups.to_parquet(out / "group_counts.parquet")
    np.save(out / "fingerprint_matrix.npy", matrix)
    with open(out / "fingerprints.pkl", "wb") as handle:
        pickle.dump(fps, handle)

    report.to_frame().to_csv(results / "table_curation.csv", index=False)
    family_summary(table).to_csv(results / "table_families.csv", index=False)
    pd.DataFrame(
        sorted(report.attachment_point_counts.items()),
        columns=["n_attachment_points", "n_structures"],
    ).to_csv(results / "table_attachment_points.csv", index=False)

    print(f"\nstructures {table.shape}  descriptors {descriptors.shape}  fp {matrix.shape}")
    print(family_summary(table).to_string(index=False))


if __name__ == "__main__":
    main()
