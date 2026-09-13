"""Emit supplementary tables that document the family classifier's chemistry.

Table S1 is the classifier's validation set: every case is a polymer whose family
is not in dispute, so a reader can check the taxonomy against their own chemical
judgement rather than taking it on trust.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import _bootstrap  # noqa: F401
from ptgbench.families import assign_family

VALIDATION_CASES = [
    ("Poly(ethylene)", "*CC*", "Polyolefins"),
    ("Poly(propylene)", "*CC(C)*", "Polyolefins"),
    ("Poly(styrene)", "*CC(*)c1ccccc1", "Polystyrenes"),
    ("Poly(methyl acrylate)", "*CC(*)C(=O)OC", "Polyacrylics"),
    ("Poly(methyl methacrylate)", "*CC(C)(*)C(=O)OC", "Polyacrylics"),
    ("Poly(butyl acrylate)", "*CC(*)C(=O)OCCCC", "Polyacrylics"),
    ("Poly(acrylonitrile)", "*CC(*)C#N", "Polyacrylics"),
    ("Poly(acrylamide)", "*CC(*)C(=O)N", "Polyacrylics"),
    ("Poly(vinyl chloride)", "*CC(*)Cl", "Polyhalo-olefins"),
    ("Poly(tetrafluoroethylene)", "*C(F)(F)C(F)(F)*", "Polyhalo-olefins"),
    ("Poly(vinyl alcohol)", "*CC(*)O", "Polyvinyls"),
    ("Poly(vinyl acetate)", "*CC(*)OC(C)=O", "Polyvinyls"),
    ("Poly(ethylene oxide)", "*OCC*", "Polyethers"),
    ("Poly(ethylene terephthalate)", "*OC(=O)c1ccc(C(=O)OCC*)cc1", "Polyesters"),
    ("Aliphatic polyester", "*OC(=O)CCCCC(=O)OCC*", "Polyesters"),
    ("Aliphatic polyamide", "*NC(=O)CCCCC(=O)N*", "Polyamides"),
    ("Polyurethane", "*NC(=O)OCCO*", "Polyurethanes"),
    ("Polyurea", "*NC(=O)N*", "Polyureas"),
    ("Bisphenol-A polycarbonate", "*OC(=O)Oc1ccc(C(C)(C)c2ccc(*)cc2)cc1", "Polycarbonates"),
    ("Poly(dimethylsiloxane)", "*O[Si](C)(C)*", "Polysiloxanes"),
    ("Aromatic polyimide", "*N(*)C(=O)c1ccc2c(c1)C(=O)N(*)C2=O", "Polyimides"),
    ("Poly(p-phenylene)", "*c1ccc(*)cc1", "Polyphenylenes"),
    ("Poly(butadiene)", "*CC=CC*", "Polydienes"),
    ("Poly(ethylene sulfide)", "*SCC*", "Polysulfides"),
    ("Polyanhydride", "*C(=O)OC(=O)CC*", "Polyanhydrides"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="results")
    args = parser.parse_args()

    rows = []
    for name, psmiles, expected in VALIDATION_CASES:
        assigned = assign_family(psmiles)
        rows.append({
            "polymer": name,
            "psmiles": psmiles,
            "expected_family": expected,
            "assigned_family": assigned,
            "agrees": assigned == expected,
        })
    frame = pd.DataFrame(rows)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "table_S1_family_validation.csv", index=False)

    n_agree = int(frame["agrees"].sum())
    print(frame.to_string(index=False))
    print(f"\n{n_agree}/{len(frame)} reference polymers assigned as expected")
    if n_agree != len(frame):
        raise SystemExit("family classifier disagrees with a reference assignment")


if __name__ == "__main__":
    main()
