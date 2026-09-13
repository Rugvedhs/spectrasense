"""Partitioning schemes, from optimistic to deliberately adversarial.

Four regimes are used, in increasing order of the structural novelty they force
the model to confront:

``random``   repeated random partitions; measures interpolation inside
             well-populated chemistry.
``scaffold`` Bemis-Murcko scaffold partitions; no scaffold is shared between
             train and test.
``cluster``  agglomerative clusters in Morgan space are held out whole, so test
             repeat units have no close training analogue.
``family``   an entire polymer family is withheld; the model must extrapolate to
             a linkage chemistry it has never seen.

The fifth entry point, :func:`matched_pair_splits`, is not another difficulty
level but a control design.  Comparing a family-holdout score against a random
split score confounds two things: the test *population* changes (families differ
in how hard they intrinsically are) and the training *set* changes.  The matched
design fixes the test set and varies only the training pool, and adds a
size-matched arm so that the effect of losing a family is separated from the
effect of simply training on less data.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.cluster import AgglomerativeClustering


@dataclass(frozen=True)
class Split:
    """One train/test partition, with the metadata needed to interpret it."""

    name: str
    regime: str
    train_idx: np.ndarray
    test_idx: np.ndarray
    # Populated only by the matched-pair design.
    arm: str | None = None
    group: str | None = None

    def __post_init__(self) -> None:
        overlap = np.intersect1d(self.train_idx, self.test_idx)
        if overlap.size:
            raise ValueError(
                f"split {self.name!r} leaks {overlap.size} structures between "
                "train and test"
            )


def murcko_scaffold(psmiles: str) -> str:
    """Generic Bemis-Murcko scaffold, with attachment points stripped.

    Dummy atoms are removed first: retaining them would make almost every repeat
    unit its own scaffold and reduce the split to a random one.
    """
    mol = Chem.MolFromSmiles(psmiles)
    if mol is None:
        return ""
    editable = Chem.RWMol(mol)
    for atom in sorted(
        (a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0), reverse=True
    ):
        editable.RemoveAtom(atom)
    stripped = editable.GetMol()
    try:
        Chem.SanitizeMol(stripped)
        scaffold = MurckoScaffold.GetScaffoldForMol(stripped)
        return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaffold))
    except Exception:
        return ""


def random_splits(
    n: int, n_repeats: int = 10, test_fraction: float = 0.2, seed: int = 0
) -> list[Split]:
    """Repeated random partitions.

    Repetition matters: a single random split gives no sense of how much of a
    model ranking is sampling noise, and much of the published spread between
    tree ensembles on this dataset is within that noise.
    """
    rng = np.random.default_rng(seed)
    n_test = int(round(n * test_fraction))
    splits = []
    for r in range(n_repeats):
        perm = rng.permutation(n)
        splits.append(
            Split(
                name=f"random_{r}",
                regime="random",
                train_idx=np.sort(perm[n_test:]),
                test_idx=np.sort(perm[:n_test]),
            )
        )
    return splits


def _grouped_holdout_splits(
    groups: np.ndarray, regime: str, n_repeats: int, test_fraction: float, seed: int
) -> list[Split]:
    """Assign whole groups to the test side until the target fraction is met."""
    rng = np.random.default_rng(seed)
    unique = np.unique(groups)
    n = len(groups)
    target = n * test_fraction

    splits = []
    for r in range(n_repeats):
        order = rng.permutation(unique)
        test_groups, size = [], 0
        for g in order:
            if size >= target:
                break
            test_groups.append(g)
            size += int((groups == g).sum())
        mask = np.isin(groups, test_groups)
        splits.append(
            Split(
                name=f"{regime}_{r}",
                regime=regime,
                train_idx=np.flatnonzero(~mask),
                test_idx=np.flatnonzero(mask),
            )
        )
    return splits


def scaffold_splits(
    psmiles: list[str], n_repeats: int = 10, test_fraction: float = 0.2, seed: int = 0
) -> list[Split]:
    scaffolds = np.array([murcko_scaffold(s) for s in psmiles])
    return _grouped_holdout_splits(scaffolds, "scaffold", n_repeats, test_fraction, seed)


def fingerprint_clusters(
    fingerprint_matrix: np.ndarray, n_clusters: int = 60
) -> np.ndarray:
    """Agglomerative clusters in Jaccard space over Morgan bits.

    Average linkage on Jaccard distance groups repeat units that share
    substructural environments, which is the property a cluster holdout needs:
    holding out a cluster removes a *neighbourhood*, not scattered points.
    """
    bits = fingerprint_matrix.astype(bool)
    n = bits.shape[0]
    intersection = bits.astype(np.float32) @ bits.astype(np.float32).T
    counts = bits.sum(axis=1).astype(np.float32)
    union = counts[:, None] + counts[None, :] - intersection
    distance = 1.0 - np.divide(
        intersection, union, out=np.zeros_like(intersection), where=union > 0
    )
    np.fill_diagonal(distance, 0.0)
    model = AgglomerativeClustering(
        n_clusters=min(n_clusters, n), metric="precomputed", linkage="average"
    )
    return model.fit_predict(distance)


def cluster_splits(
    clusters: np.ndarray,
    n_repeats: int = 10,
    test_fraction: float = 0.2,
    seed: int = 0,
) -> list[Split]:
    return _grouped_holdout_splits(
        clusters, "cluster", n_repeats, test_fraction, seed
    )


def family_holdout_splits(
    families: np.ndarray, min_size: int = 40
) -> list[Split]:
    """Leave-one-family-out, for every family large enough to score reliably."""
    counts = pd.Series(families).value_counts()
    eligible = counts[counts >= min_size].index.tolist()

    splits = []
    for family in sorted(eligible):
        mask = families == family
        splits.append(
            Split(
                name=f"family_{family}",
                regime="family",
                train_idx=np.flatnonzero(~mask),
                test_idx=np.flatnonzero(mask),
                group=family,
            )
        )
    return splits


def matched_pair_splits(
    families: np.ndarray,
    min_size: int = 40,
    test_fraction: float = 0.3,
    n_repeats: int = 5,
    seed: int = 0,
) -> list[Split]:
    """Fixed test set, three training pools differing only in what is removed.

    For each eligible family ``F`` a test set ``T`` is drawn from ``F``.  Three
    models then predict that *same* ``T``:

    ``informed``   trained on everything except ``T`` — so the rest of ``F`` is
                   available and the model has seen the linkage chemistry.
    ``naive``      trained on everything except ``F`` — the family is absent
                   entirely, which is the extrapolation case.
    ``control``    trained on ``informed`` minus a random block of other-family
                   structures of the same size as ``F \\ T``.

    ``naive`` versus ``control`` is the contrast that matters: both pools are the
    same size, so a difference between them is attributable to *which* chemistry
    was removed rather than to how much data was removed.
    """
    rng = np.random.default_rng(seed)
    counts = pd.Series(families).value_counts()
    eligible = sorted(counts[counts >= min_size].index.tolist())
    n = len(families)

    splits: list[Split] = []
    for family in eligible:
        member_idx = np.flatnonzero(families == family)
        n_test = max(1, int(round(len(member_idx) * test_fraction)))

        for r in range(n_repeats):
            test_idx = np.sort(rng.choice(member_idx, size=n_test, replace=False))
            rest_of_family = np.setdiff1d(member_idx, test_idx)

            informed = np.setdiff1d(np.arange(n), test_idx)
            naive = np.setdiff1d(np.arange(n), member_idx)

            outside = np.setdiff1d(informed, rest_of_family)
            n_drop = min(len(rest_of_family), len(outside) - 1)
            dropped = rng.choice(outside, size=n_drop, replace=False)
            control = np.setdiff1d(informed, dropped)

            for arm, train_idx in (
                ("informed", informed),
                ("naive", naive),
                ("control", control),
            ):
                splits.append(
                    Split(
                        name=f"matched_{family}_{r}_{arm}",
                        regime="matched",
                        train_idx=np.sort(train_idx),
                        test_idx=test_idx,
                        arm=arm,
                        group=family,
                    )
                )
    return splits


def calibration_split(
    train_idx: np.ndarray, calib_fraction: float = 0.25, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """Carve a conformal calibration set out of a training pool.

    Split conformal needs calibration residuals from data the point predictor
    never saw; reusing training residuals would understate them and silently
    break the coverage guarantee.
    """
    rng = np.random.default_rng(seed)
    perm = rng.permutation(train_idx)
    n_calib = max(1, int(round(len(train_idx) * calib_fraction)))
    return np.sort(perm[n_calib:]), np.sort(perm[:n_calib])
