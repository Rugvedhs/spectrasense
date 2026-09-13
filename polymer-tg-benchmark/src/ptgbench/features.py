"""Molecular representations for repeat units.

Two representations are used throughout:

* the full RDKit descriptor block, as a physically interpretable feature set, and
* binary Morgan fingerprints, used both as an alternative model input and as the
  similarity space in which structural novelty is measured.

Attachment points (``*``) are retained.  They mark where the chain continues, so
they carry real information about the linkage chemistry, and Morgan environments
centred near them describe exactly the backbone motifs that set Tg.

Descriptor *filtering* deliberately does not live here.  Dropping constant,
degenerate or duplicated columns is a data-dependent decision, so it belongs
inside a pipeline fitted on the training partition only — see
:class:`DescriptorFilter`.  Computing the full block once and filtering later is
what keeps the evaluation free of preprocessing leakage.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from sklearn.base import BaseEstimator, TransformerMixin

RDLogger.DisableLog("rdApp.*")

# Ipc grows factorially with graph size and overflows to inf for larger repeat
# units; the log form is the standard remedy and is monotone in the original.
_LOG_SCALE_DESCRIPTORS = ("Ipc",)


def descriptor_names() -> list[str]:
    return [name for name, _ in Descriptors.descList]


def compute_descriptors(psmiles: list[str]) -> pd.DataFrame:
    """Full RDKit descriptor block, one row per repeat unit."""
    names = descriptor_names()
    calc = dict(Descriptors.descList)
    rows = []
    for smi in psmiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append([np.nan] * len(names))
            continue
        row = []
        for name in names:
            try:
                value = calc[name](mol)
            except Exception:
                value = np.nan
            row.append(value)
        rows.append(row)

    frame = pd.DataFrame(rows, columns=names, dtype=float)
    for name in _LOG_SCALE_DESCRIPTORS:
        if name in frame.columns:
            frame[name] = np.log10(np.clip(frame[name].to_numpy(), 1.0, None))
    return frame.replace([np.inf, -np.inf], np.nan)


def morgan_generator(radius: int = 2, n_bits: int = 2048):
    return GetMorganGenerator(radius=radius, fpSize=n_bits)


def compute_fingerprints(
    psmiles: list[str], radius: int = 2, n_bits: int = 2048
) -> tuple[np.ndarray, list]:
    """Return a dense bit matrix and the RDKit fingerprint objects.

    The fingerprint objects are kept because RDKit's bulk Tanimoto routines are
    far faster than recomputing similarity from the dense matrix.
    """
    gen = morgan_generator(radius, n_bits)
    fps, matrix = [], np.zeros((len(psmiles), n_bits), dtype=np.uint8)
    for i, smi in enumerate(psmiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fps.append(None)
            continue
        fp = gen.GetFingerprint(mol)
        fps.append(fp)
        arr = np.zeros((n_bits,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        matrix[i] = arr
    return matrix, fps


def max_similarity_to_reference(
    query_fps: list, reference_fps: list, batch: int = 512
) -> np.ndarray:
    """Nearest-neighbour Tanimoto similarity of each query to a reference set.

    This is the structural-novelty coordinate used for applicability-domain
    analysis and for the Mondrian conformal taxonomy: a repeat unit whose closest
    training analogue is remote is one the model has no local evidence for.
    """
    reference = [f for f in reference_fps if f is not None]
    if not reference:
        return np.zeros(len(query_fps), dtype=float)

    out = np.zeros(len(query_fps), dtype=float)
    for i, fp in enumerate(query_fps):
        if fp is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, reference)
        out[i] = max(sims) if sims else 0.0
    return out


class DescriptorFilter(BaseEstimator, TransformerMixin):
    """Drop non-finite, constant and duplicated descriptor columns.

    Fitted on the training partition only.  Columns that survive are recorded so
    the same selection is applied unchanged to validation and test data; any
    residual non-finite value in unseen data is replaced by the training median
    rather than discarded, so no test row is ever silently dropped.
    """

    def __init__(self, variance_threshold: float = 0.0):
        self.variance_threshold = variance_threshold

    def fit(self, X, y=None):
        frame = pd.DataFrame(X).astype(float).replace([np.inf, -np.inf], np.nan)
        finite = frame.columns[frame.notna().all(axis=0)]
        frame = frame[finite]

        variances = frame.var(axis=0, ddof=0)
        keep = variances[variances > self.variance_threshold].index
        frame = frame[keep]

        # Exact duplicates carry no extra information and destabilise
        # tree-importance attribution; keep the first occurrence of each.
        deduped = frame.T.drop_duplicates().T

        self.columns_ = list(deduped.columns)
        self.medians_ = deduped.median(axis=0).to_numpy(dtype=float)
        self.n_features_in_ = frame.shape[1] if hasattr(frame, "shape") else None
        return self

    def transform(self, X):
        frame = pd.DataFrame(X).astype(float).replace([np.inf, -np.inf], np.nan)
        frame = frame.reindex(columns=self.columns_)
        # copy=True: a reindexed single-dtype frame can hand back a read-only
        # view, which the imputation below writes into.
        values = np.array(frame.to_numpy(dtype=float), copy=True)
        missing = np.isnan(values)
        if missing.any():
            values[missing] = np.take(self.medians_, np.where(missing)[1])
        return values

    def get_feature_names_out(self, input_features=None):
        return np.asarray(self.columns_, dtype=object)
