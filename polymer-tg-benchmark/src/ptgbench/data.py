"""Loading, auditing and structure-level aggregation of the repeat-unit records.

The raw table is a record-level export: the same repeat unit can appear several
times with Tg values from different primary sources.  Modelling on record level
would let identical inputs straddle a train/test boundary and would silently
weight frequently reported polymers more heavily, so everything downstream works
on a *structure-level* table with one row per canonical repeat unit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from .families import assign_families

RDLogger.DisableLog("rdApp.*")

RAW_FILENAME = "raw_tg_7208.csv"


@dataclass
class CurationReport:
    """Counts recorded at each curation step, for the data-provenance table."""

    n_raw_records: int = 0
    n_unique_raw_smiles: int = 0
    n_parse_failures: int = 0
    n_canonical_structures: int = 0
    n_aggregated_groups: int = 0
    n_records_absorbed: int = 0
    n_disagreeing_groups: int = 0
    max_within_group_range: float = 0.0
    attachment_point_counts: dict[int, int] = field(default_factory=dict)

    def to_frame(self) -> pd.DataFrame:
        rows = [
            ("Raw records", self.n_raw_records),
            ("Unique raw SMILES", self.n_unique_raw_smiles),
            ("RDKit parse failures", self.n_parse_failures),
            ("Canonical structures retained", self.n_canonical_structures),
            ("Structures aggregated from >1 record", self.n_aggregated_groups),
            ("Records absorbed by aggregation", self.n_records_absorbed),
            ("Aggregated groups with disagreeing Tg", self.n_disagreeing_groups),
            ("Max within-structure Tg range (K)", self.max_within_group_range),
        ]
        return pd.DataFrame(rows, columns=["step", "value"])


def canonicalize(psmiles: str) -> str | None:
    """Canonical SMILES for a repeat unit, or ``None`` if RDKit rejects it."""
    mol = Chem.MolFromSmiles(psmiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def _count_attachment_points(psmiles: str) -> int:
    mol = Chem.MolFromSmiles(psmiles)
    if mol is None:
        return 0
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0)


def load_structure_level(
    raw_dir: str | Path, temperature_unit: str = "C"
) -> tuple[pd.DataFrame, CurationReport]:
    """Build the structure-level modelling table.

    Duplicate canonical structures are collapsed to their median Tg — the same
    rule the upstream curation used for multi-source records, and the robust
    choice given that disagreements between primary sources reach a hundred
    kelvin or more.  Aggregation provenance is kept in ``n_records`` and
    ``tg_range`` so that later analysis can condition on measurement spread.

    Returns the table and a :class:`CurationReport` documenting every step.
    """
    raw_path = Path(raw_dir) / RAW_FILENAME
    raw = pd.read_csv(raw_path)

    report = CurationReport(
        n_raw_records=len(raw),
        n_unique_raw_smiles=raw["psmiles"].nunique(),
    )

    raw = raw.copy()
    raw["canonical_psmiles"] = [canonicalize(s) for s in raw["psmiles"]]
    report.n_parse_failures = int(raw["canonical_psmiles"].isna().sum())
    raw = raw.dropna(subset=["canonical_psmiles"]).reset_index(drop=True)

    grouped = raw.groupby("canonical_psmiles")["Tg_C"]
    table = pd.DataFrame(
        {
            "Tg_C": grouped.median(),
            "n_records": grouped.size(),
            "tg_min": grouped.min(),
            "tg_max": grouped.max(),
        }
    ).reset_index()
    table["tg_range"] = table["tg_max"] - table["tg_min"]

    multi = table["n_records"] > 1
    report.n_canonical_structures = len(table)
    report.n_aggregated_groups = int(multi.sum())
    report.n_records_absorbed = int(table.loc[multi, "n_records"].sum() - multi.sum())
    report.n_disagreeing_groups = int((table["tg_range"] > 0).sum())
    report.max_within_group_range = float(table["tg_range"].max())

    counts = pd.Series(
        [_count_attachment_points(s) for s in table["canonical_psmiles"]]
    )
    report.attachment_point_counts = counts.value_counts().sort_index().to_dict()
    table["n_attachment_points"] = counts.to_numpy()

    table["family"] = assign_families(table["canonical_psmiles"].tolist())

    if temperature_unit.upper() == "K":
        table["Tg"] = table["Tg_C"] + 273.15
    else:
        table["Tg"] = table["Tg_C"]

    table = table.sort_values("canonical_psmiles", ignore_index=True)
    return table, report


def family_summary(table: pd.DataFrame) -> pd.DataFrame:
    """Per-family count and Tg distribution, ordered by descending size."""
    g = table.groupby("family")["Tg"]
    out = pd.DataFrame(
        {
            "n": g.size(),
            "tg_mean": g.mean(),
            "tg_std": g.std(),
            "tg_min": g.min(),
            "tg_median": g.median(),
            "tg_max": g.max(),
        }
    )
    return out.sort_values("n", ascending=False).reset_index()
