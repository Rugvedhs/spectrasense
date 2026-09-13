"""Splits are where leakage hides."""

import numpy as np
import pytest

from ptgbench.splits import (
    Split,
    calibration_split,
    family_holdout_splits,
    matched_pair_splits,
    murcko_scaffold,
    random_splits,
    scaffold_splits,
)


def test_split_rejects_overlap():
    with pytest.raises(ValueError, match="leaks"):
        Split("bad", "random", np.array([1, 2, 3]), np.array([3, 4]))


def test_random_splits_partition_exactly():
    splits = random_splits(1000, n_repeats=4, test_fraction=0.2, seed=0)
    for split in splits:
        assert len(split.train_idx) + len(split.test_idx) == 1000
        assert not set(split.train_idx) & set(split.test_idx)
        assert len(split.test_idx) == 200


def test_scaffold_splits_share_no_scaffold():
    smiles = [
        "*CC(*)c1ccccc1", "*CC(*)c1ccc(C)cc1", "*OCC*", "*OCCC*",
        "*NC(=O)CCC(=O)N*", "*NC(=O)CCCC(=O)N*", "*CC(*)C(=O)OC",
        "*CC(*)C(=O)OCC", "*O[Si](C)(C)*", "*c1ccc(*)cc1",
    ]
    for split in scaffold_splits(smiles, n_repeats=3, seed=0):
        train = {murcko_scaffold(smiles[i]) for i in split.train_idx}
        test = {murcko_scaffold(smiles[i]) for i in split.test_idx}
        assert not train & test


def test_family_holdout_removes_the_whole_family():
    families = np.array(["A"] * 60 + ["B"] * 50 + ["C"] * 10)
    splits = family_holdout_splits(families, min_size=40)
    assert {s.group for s in splits} == {"A", "B"}      # C is too small
    for split in splits:
        assert set(families[split.train_idx]).isdisjoint({split.group})
        assert set(families[split.test_idx]) == {split.group}


def test_matched_arms_share_a_test_set_and_match_on_size():
    """The control that makes the family effect interpretable."""
    families = np.array(["A"] * 80 + ["B"] * 80 + ["C"] * 80)
    splits = matched_pair_splits(families, min_size=40, n_repeats=2, seed=0)

    by_key: dict[tuple, dict] = {}
    for split in splits:
        key = (split.group, split.name.rsplit("_", 1)[0])
        by_key.setdefault(key, {})[split.arm] = split

    assert by_key
    for arms in by_key.values():
        assert set(arms) == {"informed", "control", "naive"}
        test_sets = [tuple(a.test_idx) for a in arms.values()]
        assert len(set(test_sets)) == 1, "arms must predict identical structures"
        assert len(arms["control"].train_idx) == len(arms["naive"].train_idx)
        assert len(arms["informed"].train_idx) > len(arms["naive"].train_idx)
        # The naive arm must not have seen the family at all.
        assert split_free_of(families, arms["naive"], arms["naive"].group)
        # The informed arm must have seen it.
        assert not split_free_of(families, arms["informed"], arms["informed"].group)


def split_free_of(families, split, group) -> bool:
    return group not in set(families[split.train_idx])


def test_calibration_split_is_disjoint():
    train = np.arange(500)
    fit, calib = calibration_split(train, calib_fraction=0.25, seed=0)
    assert not set(fit) & set(calib)
    assert len(fit) + len(calib) == 500
    assert len(calib) == 125
