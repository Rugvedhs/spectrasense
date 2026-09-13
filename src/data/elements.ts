// Simplified elemental reference data used for composition-only ("baseline") featurization.
// Electronegativity: Pauling scale. Radius: empirical covalent radius, picometers.
// Values are standard textbook approximations, adequate for feature engineering, not
// meant as a citable primary source.
export interface ElementRef {
  symbol: string;
  name: string;
  electronegativity: number;
  covalentRadiusPm: number;
  valenceElectrons: number;
  period: number;
  color: string; // render color in the structure viewer
}

export const ELEMENTS: Record<string, ElementRef> = {
  Pb: { symbol: "Pb", name: "Lead", electronegativity: 2.33, covalentRadiusPm: 146, valenceElectrons: 4, period: 6, color: "#5b5650" },
  Sn: { symbol: "Sn", name: "Tin", electronegativity: 1.96, covalentRadiusPm: 139, valenceElectrons: 4, period: 5, color: "#8f9ca3" },
  S:  { symbol: "S",  name: "Sulfur", electronegativity: 2.58, covalentRadiusPm: 105, valenceElectrons: 6, period: 3, color: "#d9b23c" },
  Se: { symbol: "Se", name: "Selenium", electronegativity: 2.55, covalentRadiusPm: 120, valenceElectrons: 6, period: 4, color: "#bf8f2e" },
  Te: { symbol: "Te", name: "Tellurium", electronegativity: 2.10, covalentRadiusPm: 140, valenceElectrons: 6, period: 5, color: "#8a7048" },
  In: { symbol: "In", name: "Indium", electronegativity: 1.78, covalentRadiusPm: 142, valenceElectrons: 3, period: 5, color: "#7d8fae" },
  Ga: { symbol: "Ga", name: "Gallium", electronegativity: 1.81, covalentRadiusPm: 122, valenceElectrons: 3, period: 4, color: "#6f9c92" },
  As: { symbol: "As", name: "Arsenic", electronegativity: 2.18, covalentRadiusPm: 119, valenceElectrons: 5, period: 4, color: "#7a689c" },
  Sb: { symbol: "Sb", name: "Antimony", electronegativity: 2.05, covalentRadiusPm: 139, valenceElectrons: 5, period: 5, color: "#8a6ba0" },
  P:  { symbol: "P",  name: "Phosphorus", electronegativity: 2.19, covalentRadiusPm: 107, valenceElectrons: 5, period: 3, color: "#c17a3e" },
  Cd: { symbol: "Cd", name: "Cadmium", electronegativity: 1.69, covalentRadiusPm: 144, valenceElectrons: 2, period: 5, color: "#a888b8" },
  Ge: { symbol: "Ge", name: "Germanium", electronegativity: 2.01, covalentRadiusPm: 122, valenceElectrons: 4, period: 4, color: "#5f8a7c" },
  Si: { symbol: "Si", name: "Silicon", electronegativity: 1.90, covalentRadiusPm: 111, valenceElectrons: 4, period: 3, color: "#5a7f8c" },
};

export type ElementSymbol = keyof typeof ELEMENTS;
