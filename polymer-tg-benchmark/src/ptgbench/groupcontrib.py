"""A group-contribution baseline in the Van Krevelen additive form.

Referees in polymer journals reasonably ask how a machine-learned model compares
against the additive schemes it proposes to replace. This implements the classical
form

    Tg = Yg / M,    Yg = sum_i n_i Y_i

where ``M`` is the repeat-unit molar mass, ``n_i`` the count of structural group
``i``, and ``Y_i`` its molar glass-transition contribution.

The group contributions ``Y_i`` are **fitted to the training partition** rather
than taken from a published table. That choice is deliberate and matters for the
comparison: it gives the additive scheme the same data advantage the machine-learned
models get, so any difference between them reflects the functional form rather than
the calibration set. It also means the baseline is subject to the same split
regimes, which is the point — additivity is a compositional assumption, and whether
compositional models extrapolate to unseen backbone chemistry better or worse than
ensembles is exactly the question the split regimes are built to answer.

Working in ``Yg = Tg * M`` space rather than on Tg directly is what makes the model
additive: contributions sum over groups, and dividing by molar mass at the end
restores an intensive temperature. Absolute temperature is used throughout, since
an additive law expressed in Celsius would make the contributions depend on an
arbitrary zero.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.linear_model import Ridge

RDLogger.DisableLog("rdApp.*")

KELVIN_OFFSET = 273.15

# Element and topology counts, on top of RDKit's functional-group fragments.
# These carry the backbone skeleton that fragment counts alone would miss.
_ELEMENTS = ("C", "N", "O", "S", "Si", "P", "F", "Cl", "Br", "I")


def _fragment_functions():
    return [(name, fn) for name, fn in Descriptors.descList if name.startswith("fr_")]


def group_count_names() -> list[str]:
    names = [name for name, _ in _fragment_functions()]
    names += [f"n_{element}" for element in _ELEMENTS]
    names += [
        "n_aromatic_atoms", "n_rings", "n_rotatable", "n_hbond_donors",
        "n_hbond_acceptors", "n_heavy_atoms", "n_backbone_atoms", "intercept",
    ]
    return names


def compute_group_counts(psmiles: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    """Integer group counts per repeat unit, plus the molar mass of each.

    The trailing constant column plays the role of a chain-end / baseline term:
    without it every contribution is forced through the origin, which fits the
    additive law badly for small repeat units.
    """
    from .families import backbone_atoms

    fragments = _fragment_functions()
    rows, masses = [], []
    for smiles in psmiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            rows.append([0.0] * len(group_count_names()))
            masses.append(np.nan)
            continue

        row = []
        for _, fn in fragments:
            try:
                row.append(float(fn(mol)))
            except Exception:
                row.append(0.0)

        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        row += [float(symbols.count(e)) for e in _ELEMENTS]
        row += [
            float(sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())),
            float(mol.GetRingInfo().NumRings()),
            float(Descriptors.NumRotatableBonds(mol)),
            float(Descriptors.NumHDonors(mol)),
            float(Descriptors.NumHAcceptors(mol)),
            float(mol.GetNumHeavyAtoms()),
            float(len(backbone_atoms(mol))),
            1.0,
        ]
        rows.append(row)
        masses.append(Descriptors.MolWt(mol))

    frame = pd.DataFrame(rows, columns=group_count_names(), dtype=float)
    return frame, np.asarray(masses, dtype=float)


class GroupContributionRegressor(BaseEstimator, RegressorMixin):
    """Additive group-contribution regressor, fitted in molar (Yg) space.

    The design matrix must carry repeat-unit molar mass in its final column; the
    remaining columns are group counts. Ridge regularisation is used because
    fragment counts are strongly collinear — many groups co-occur in every member
    of a family — and an unregularised fit puts wild compensating values on rare
    groups, which is precisely the behaviour that would collapse under a family
    holdout for the wrong reason.
    """

    def __init__(self, alpha: float = 1.0):
        self.alpha = alpha

    @staticmethod
    def _split(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        X = np.asarray(X, dtype=float)
        return X[:, :-1], X[:, -1]

    def fit(self, X, y):
        counts, mass = self._split(X)
        mass = np.where(np.isfinite(mass) & (mass > 0), mass, 1.0)
        molar_target = (np.asarray(y, dtype=float) + KELVIN_OFFSET) * mass
        self.model_ = Ridge(alpha=self.alpha, fit_intercept=False)
        self.model_.fit(counts, molar_target)
        return self

    def predict(self, X):
        counts, mass = self._split(X)
        mass = np.where(np.isfinite(mass) & (mass > 0), mass, 1.0)
        kelvin = self.model_.predict(counts) / mass
        # An additive fit can produce physically impossible temperatures for
        # repeat units unlike anything it was fitted on. Clipping at absolute
        # zero keeps the error finite without hiding that the prediction is poor.
        return np.clip(kelvin, 1.0, None) - KELVIN_OFFSET

    @property
    def contributions_(self) -> np.ndarray:
        return self.model_.coef_
