export type StructureType = "rocksalt" | "zincblende" | "wurtzite" | "diamond";

export type Family =
  | "IV-VI rocksalt (Pb-chalcogenide)"
  | "IV-VI rocksalt (Sn-chalcogenide)"
  | "III-V zincblende"
  | "II-VI zincblende"
  | "II-VI wurtzite"
  | "Group IV (diamond)";

export interface Component {
  symbol: string;
  frac: number; // site occupancy fraction, sums to 1 within a sublattice
}

export interface Material {
  id: string;
  formula: string;
  displayName: string;
  cations: Component[];
  anions: Component[];
  structure: StructureType;
  family: Family;
  latticeConstantAng: number;
  actualEgEv: number;
  note: string;
}

// Curated set of known NIR-photodetector-adjacent chemistries, spanning the
// families the spec calls out: PbS/PbSe (IV-VI), InGaAs-adjacent (III-V),
// CdTe-family (II-VI), plus group-IV anchors (Si, Ge) for contrast.
export const MATERIALS: Material[] = [
  {
    id: "PbS", formula: "PbS", displayName: "Lead sulfide",
    cations: [{ symbol: "Pb", frac: 1 }], anions: [{ symbol: "S", frac: 1 }],
    structure: "rocksalt", family: "IV-VI rocksalt (Pb-chalcogenide)",
    latticeConstantAng: 5.936, actualEgEv: 0.41,
    note: "Classic colloidal-quantum-dot NIR/SWIR photodetector material.",
  },
  {
    id: "PbSe", formula: "PbSe", displayName: "Lead selenide",
    cations: [{ symbol: "Pb", frac: 1 }], anions: [{ symbol: "Se", frac: 1 }],
    structure: "rocksalt", family: "IV-VI rocksalt (Pb-chalcogenide)",
    latticeConstantAng: 6.117, actualEgEv: 0.27,
    note: "Narrower gap than PbS; tuned via quantum confinement for SWIR imaging.",
  },
  {
    id: "PbTe", formula: "PbTe", displayName: "Lead telluride",
    cations: [{ symbol: "Pb", frac: 1 }], anions: [{ symbol: "Te", frac: 1 }],
    structure: "rocksalt", family: "IV-VI rocksalt (Pb-chalcogenide)",
    latticeConstantAng: 6.462, actualEgEv: 0.32,
    note: "Gap set by relativistic band inversion, not simple ionicity trends.",
  },
  {
    id: "SnTe", formula: "SnTe", displayName: "Tin telluride",
    cations: [{ symbol: "Sn", frac: 1 }], anions: [{ symbol: "Te", frac: 1 }],
    structure: "rocksalt", family: "IV-VI rocksalt (Sn-chalcogenide)",
    latticeConstantAng: 6.327, actualEgEv: 0.18,
    note: "Topological-crystalline-insulator parent compound; very narrow gap.",
  },
  {
    id: "InAs", formula: "InAs", displayName: "Indium arsenide",
    cations: [{ symbol: "In", frac: 1 }], anions: [{ symbol: "As", frac: 1 }],
    structure: "zincblende", family: "III-V zincblende",
    latticeConstantAng: 6.058, actualEgEv: 0.36,
    note: "High-mobility narrow-gap III-V, used in SWIR avalanche photodiodes.",
  },
  {
    id: "InSb", formula: "InSb", displayName: "Indium antimonide",
    cations: [{ symbol: "In", frac: 1 }], anions: [{ symbol: "Sb", frac: 1 }],
    structure: "zincblende", family: "III-V zincblende",
    latticeConstantAng: 6.479, actualEgEv: 0.17,
    note: "Narrowest-gap III-V in this set; mainly a mid-IR material at the SWIR edge.",
  },
  {
    id: "GaSb", formula: "GaSb", displayName: "Gallium antimonide",
    cations: [{ symbol: "Ga", frac: 1 }], anions: [{ symbol: "Sb", frac: 1 }],
    structure: "zincblende", family: "III-V zincblende",
    latticeConstantAng: 6.096, actualEgEv: 0.73,
    note: "Substrate and barrier layer for InGaAs/GaSb SWIR photodiode stacks.",
  },
  {
    id: "InGaAs-lattice-matched", formula: "In₀.₅₃Ga₀.₄₇As", displayName: "InGaAs (telecom, lattice-matched to InP)",
    cations: [{ symbol: "In", frac: 0.53 }, { symbol: "Ga", frac: 0.47 }], anions: [{ symbol: "As", frac: 1 }],
    structure: "zincblende", family: "III-V zincblende",
    latticeConstantAng: 5.869, actualEgEv: 0.75,
    note: "The incumbent commercial NIR detector material, used across the 900-1700 nm telecom band.",
  },
  {
    id: "InP", formula: "InP", displayName: "Indium phosphide",
    cations: [{ symbol: "In", frac: 1 }], anions: [{ symbol: "P", frac: 1 }],
    structure: "zincblende", family: "III-V zincblende",
    latticeConstantAng: 5.869, actualEgEv: 1.35,
    note: "Substrate/window layer for InGaAs detectors; near the upper edge of the target window.",
  },
  {
    id: "GaAs", formula: "GaAs", displayName: "Gallium arsenide",
    cations: [{ symbol: "Ga", frac: 1 }], anions: [{ symbol: "As", frac: 1 }],
    structure: "zincblende", family: "III-V zincblende",
    latticeConstantAng: 5.653, actualEgEv: 1.42,
    note: "Just above the target window; included as an upper-bound reference point.",
  },
  {
    id: "CdTe", formula: "CdTe", displayName: "Cadmium telluride",
    cations: [{ symbol: "Cd", frac: 1 }], anions: [{ symbol: "Te", frac: 1 }],
    structure: "zincblende", family: "II-VI zincblende",
    latticeConstantAng: 6.482, actualEgEv: 1.50,
    note: "Anchors the CdTe-family named explicitly in the spec's material scope.",
  },
  {
    id: "CdSe", formula: "CdSe", displayName: "Cadmium selenide",
    cations: [{ symbol: "Cd", frac: 1 }], anions: [{ symbol: "Se", frac: 1 }],
    structure: "wurtzite", family: "II-VI wurtzite",
    latticeConstantAng: 4.30, actualEgEv: 1.74,
    note: "Wurtzite polymorph; classic quantum-dot host, just past the target window.",
  },
  {
    id: "Si", formula: "Si", displayName: "Silicon",
    cations: [{ symbol: "Si", frac: 1 }], anions: [{ symbol: "Si", frac: 1 }],
    structure: "diamond", family: "Group IV (diamond)",
    latticeConstantAng: 5.431, actualEgEv: 1.12,
    note: "Ubiquitous visible/near-NIR detector material; the group-IV control case.",
  },
  {
    id: "Ge", formula: "Ge", displayName: "Germanium",
    cations: [{ symbol: "Ge", frac: 1 }], anions: [{ symbol: "Ge", frac: 1 }],
    structure: "diamond", family: "Group IV (diamond)",
    latticeConstantAng: 5.658, actualEgEv: 0.67,
    note: "Standard NIR photodiode material below 1550 nm; a second group-IV control case.",
  },
];

export const MATERIAL_BY_ID: Record<string, Material> = Object.fromEntries(
  MATERIALS.map((m) => [m.id, m])
);
