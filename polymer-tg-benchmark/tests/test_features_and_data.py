"""Feature handling must not let test structures shape the feature definition."""

import numpy as np
import pandas as pd
import pytest

from ptgbench.data import canonicalize
from ptgbench.features import (
    DescriptorFilter,
    compute_descriptors,
    compute_fingerprints,
    max_similarity_to_reference,
)


def test_canonicalisation_collapses_equivalent_writings():
    assert canonicalize("*CC(*)c1ccccc1") == canonicalize("c1ccccc1C(*)C*")


def test_canonicalisation_reports_failure():
    assert canonicalize("this is not a molecule") is None


def test_descriptor_filter_is_fitted_on_training_data_only():
    train = pd.DataFrame({"constant": [1.0] * 10, "varies": np.arange(10.0),
                          "dup": np.arange(10.0)})
    # 'constant' varies in the test block, but the filter never saw that.
    test = pd.DataFrame({"constant": np.arange(10.0), "varies": np.arange(10.0),
                         "dup": np.arange(10.0)})
    flt = DescriptorFilter().fit(train)
    assert "constant" not in flt.columns_
    assert len(flt.columns_) == 1          # 'dup' is an exact duplicate
    assert flt.transform(test).shape == (10, 1)


def test_descriptor_filter_imputes_rather_than_drops_rows():
    train = pd.DataFrame({"a": np.arange(10.0), "b": np.arange(10.0) * 2})
    flt = DescriptorFilter().fit(train)
    test = pd.DataFrame({"a": [np.nan, 1.0], "b": [np.inf, 2.0]})
    out = flt.transform(test)
    assert out.shape == (2, 2)
    assert np.isfinite(out).all()


def test_descriptors_are_finite_after_log_scaling():
    """Ipc overflows to inf for larger repeat units without the log transform."""
    big = "*" + "c1ccc(-c2ccc(" * 3 + "C" + ")cc2)cc1" * 3 + "*"
    frame = compute_descriptors(["*CC*", big])
    assert "Ipc" in frame.columns
    assert np.isfinite(frame["Ipc"]).all()


def test_similarity_is_one_for_an_identical_structure():
    _, fps = compute_fingerprints(["*CC(*)c1ccccc1", "*OCC*"])
    sims = max_similarity_to_reference([fps[0]], [fps[0], fps[1]])
    assert sims[0] == pytest.approx(1.0)


def test_similarity_falls_with_structural_distance():
    _, fps = compute_fingerprints(
        ["*CC(*)c1ccccc1", "*CC(*)c1ccc(C)cc1", "*O[Si](C)(C)*"]
    )
    near = max_similarity_to_reference([fps[0]], [fps[1]])[0]
    far = max_similarity_to_reference([fps[0]], [fps[2]])[0]
    assert near > far
