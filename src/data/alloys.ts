import { MATERIAL_BY_ID, type Material } from "./materials";

// Pseudo-binary alloy systems for Sandbox mode's composition slider.
// Band gap follows the standard semiconductor-alloy interpolation with a
// quadratic bowing term: Eg(x) = (1-x)*Eg_A + x*Eg_B - b*x*(1-x).
export interface AlloySystem {
  id: string;
  label: string;
  endpointAId: string;
  endpointBId: string;
  bowingEv: number;
  varyingSublattice: "cation" | "anion";
  description: string;
}

export const ALLOY_SYSTEMS: AlloySystem[] = [
  {
    id: "InGaAs",
    label: "In₁₋ₓGaₓAs",
    endpointAId: "InAs",
    endpointBId: "GaAs",
    bowingEv: 0.477,
    varyingSublattice: "cation",
    description: "The commercial NIR-detector alloy family; x = 0.47 is the telecom composition lattice-matched to InP.",
  },
  {
    id: "PbSSe",
    label: "PbS₁₋ₓSeₓ",
    endpointAId: "PbS",
    endpointBId: "PbSe",
    bowingEv: 0.15,
    varyingSublattice: "anion",
    description: "Rocksalt IV-VI alloy; tunes SWIR cutoff wavelength across both lead-chalcogenide endpoints.",
  },
  {
    id: "CdSeTe",
    label: "CdSe₁₋ₓTeₓ",
    endpointAId: "CdSe",
    endpointBId: "CdTe",
    bowingEv: 0.3,
    varyingSublattice: "anion",
    description: "II-VI alloy spanning the wurtzite/zincblende boundary near the top of the target window.",
  },
];

export function alloyEgAt(system: AlloySystem, x: number): number {
  const a = MATERIAL_BY_ID[system.endpointAId].actualEgEv;
  const b = MATERIAL_BY_ID[system.endpointBId].actualEgEv;
  return (1 - x) * a + x * b - system.bowingEv * x * (1 - x);
}

// Builds a synthetic Material for an alloy composition so it can flow through
// the same featurization / prediction / rendering pipeline as the curated set.
export function alloyVirtualMaterial(system: AlloySystem, x: number): Material {
  const endA = MATERIAL_BY_ID[system.endpointAId];
  const endB = MATERIAL_BY_ID[system.endpointBId];
  const lerp = (p: number, q: number) => p * (1 - x) + q * x;

  const cations =
    system.varyingSublattice === "cation"
      ? [
          { symbol: endA.cations[0].symbol, frac: 1 - x },
          { symbol: endB.cations[0].symbol, frac: x },
        ]
      : endA.cations;

  const anions =
    system.varyingSublattice === "anion"
      ? [
          { symbol: endA.anions[0].symbol, frac: 1 - x },
          { symbol: endB.anions[0].symbol, frac: x },
        ]
      : endA.anions;

  return {
    id: `${system.id}-x${x.toFixed(2)}`,
    formula: system.label,
    displayName: `${system.label} (x = ${x.toFixed(2)})`,
    cations,
    anions,
    structure: endA.structure,
    family: endA.family,
    latticeConstantAng: lerp(endA.latticeConstantAng, endB.latticeConstantAng),
    actualEgEv: alloyEgAt(system, x),
    note: system.description,
  };
}
