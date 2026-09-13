import { ELEMENTS } from "../data/elements";
import type { Material, StructureType } from "../data/materials";

export type Vec3 = [number, number, number];

export interface StructureAtom {
  position: Vec3;
  color: string;
  symbol: string;
  role: "cation" | "anion";
}

export interface StructureGeometry {
  atoms: StructureAtom[];
  bonds: [Vec3, Vec3][];
  cellEdges: [Vec3, Vec3][];
  scale: number;
}

const add = (a: Vec3, b: Vec3): Vec3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const scaleV = (a: Vec3, s: number): Vec3 => [a[0] * s, a[1] * s, a[2] * s];
const dist = (a: Vec3, b: Vec3) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);

// Fractional coordinates with a 0-component also appear at the 1-boundary, so
// conventional cells render as closed boxes (corner/edge/face-shared atoms),
// the standard textbook picture rather than a truncated interior slice.
function expandBoundaryImages(frac: Vec3): Vec3[] {
  const opts = frac.map((v) => (Math.abs(v) < 1e-6 ? [0, 1] : [v]));
  const out: Vec3[] = [];
  for (const x of opts[0]) for (const y of opts[1]) for (const z of opts[2]) out.push([x, y, z]);
  return out;
}

function fracToCartesian(frac: Vec3, basis: [Vec3, Vec3, Vec3]): Vec3 {
  return add(add(scaleV(basis[0], frac[0]), scaleV(basis[1], frac[1])), scaleV(basis[2], frac[2]));
}

// FCC sublattice base positions used by rocksalt / zincblende / diamond.
const FCC_BASE: Vec3[] = [
  [0, 0, 0],
  [0, 0.5, 0.5],
  [0.5, 0, 0.5],
  [0.5, 0.5, 0],
];

function assignSymbols(components: { symbol: string; frac: number }[], count: number): string[] {
  if (components.length === 1) return Array(count).fill(components[0].symbol);
  // Deterministically split `count` sites between two components by fraction.
  const secondCount = Math.round(count * components[1].frac);
  const symbols: string[] = [];
  for (let i = 0; i < count; i++) symbols.push(i < count - secondCount ? components[0].symbol : components[1].symbol);
  return symbols;
}

function buildSublattice(
  components: { symbol: string; frac: number }[],
  baseOffsets: Vec3[],
  basis: [Vec3, Vec3, Vec3],
  role: "cation" | "anion"
): StructureAtom[] {
  const symbols = assignSymbols(components, baseOffsets.length);
  const atoms: StructureAtom[] = [];
  baseOffsets.forEach((offset, i) => {
    const symbol = symbols[i];
    const color = ELEMENTS[symbol]?.color ?? "#999999";
    for (const frac of expandBoundaryImages(offset)) {
      atoms.push({ position: fracToCartesian(frac, basis), color, symbol, role });
    }
  });
  return atoms;
}

function cellWireframe(basis: [Vec3, Vec3, Vec3]): [Vec3, Vec3][] {
  const corners: Vec3[] = [];
  for (const x of [0, 1]) for (const y of [0, 1]) for (const z of [0, 1]) corners.push(fracToCartesian([x, y, z], basis));
  const idx = (x: number, y: number, z: number) => x * 4 + y * 2 + z;
  const edges: [Vec3, Vec3][] = [];
  for (const x of [0, 1]) for (const y of [0, 1]) for (const z of [0, 1]) {
    if (x === 0) edges.push([corners[idx(0, y, z)], corners[idx(1, y, z)]]);
    if (y === 0) edges.push([corners[idx(x, 0, z)], corners[idx(x, 1, z)]]);
    if (z === 0) edges.push([corners[idx(x, y, 0)], corners[idx(x, y, 1)]]);
  }
  return edges;
}

function nearestNeighborBonds(cations: StructureAtom[], anions: StructureAtom[]): [Vec3, Vec3][] {
  let minDist = Infinity;
  for (const c of cations) for (const a of anions) minDist = Math.min(minDist, dist(c.position, a.position));
  const cutoff = minDist * 1.06;
  const bonds: [Vec3, Vec3][] = [];
  for (const c of cations) {
    for (const a of anions) {
      if (dist(c.position, a.position) <= cutoff) bonds.push([c.position, a.position]);
    }
  }
  return bonds;
}

function cubicGeometry(m: Material, diamond: boolean): StructureGeometry {
  const scale = m.latticeConstantAng / 6; // normalizes typical 5.4-6.5 Ang cells to ~1 scene unit
  const basis: [Vec3, Vec3, Vec3] = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
  ];
  const cationOffsets = FCC_BASE;
  const anionOffsets: Vec3[] = diamond
    ? FCC_BASE.map((p) => add(p, [0.25, 0.25, 0.25]))
    : m.structure === "rocksalt"
    ? FCC_BASE.map((p) => add(p, [0.5, 0, 0]))
    : FCC_BASE.map((p) => add(p, [0.25, 0.25, 0.25]));

  const cationAtoms = buildSublattice(m.cations, cationOffsets, basis, "cation");
  const anionAtoms = diamond
    ? buildSublattice(m.cations, anionOffsets, basis, "anion") // same species, second FCC sublattice
    : buildSublattice(m.anions, anionOffsets, basis, "anion");

  return {
    atoms: [...cationAtoms, ...anionAtoms],
    bonds: nearestNeighborBonds(cationAtoms, anionAtoms),
    cellEdges: cellWireframe(basis),
    scale,
  };
}

function wurtziteGeometry(m: Material): StructureGeometry {
  const a = 1;
  const c = 1.62 * a;
  const basis: [Vec3, Vec3, Vec3] = [
    [a, 0, 0],
    [-a / 2, (a * Math.sqrt(3)) / 2, 0],
    [0, 0, c],
  ];
  const u = 0.375;
  const cationOffsets: Vec3[] = [
    [1 / 3, 2 / 3, 0],
    [2 / 3, 1 / 3, 0.5],
  ];
  const anionOffsets: Vec3[] = [
    [1 / 3, 2 / 3, u],
    [2 / 3, 1 / 3, 0.5 + u],
  ];
  const cationAtoms = buildSublattice(m.cations, cationOffsets, basis, "cation");
  const anionAtoms = buildSublattice(m.anions, anionOffsets, basis, "anion");
  const scale = m.latticeConstantAng / 4.3;

  return {
    atoms: [...cationAtoms, ...anionAtoms],
    bonds: nearestNeighborBonds(cationAtoms, anionAtoms),
    cellEdges: cellWireframe(basis),
    scale,
  };
}

export function buildStructureGeometry(m: Material): StructureGeometry {
  const kind: StructureType = m.structure;
  if (kind === "wurtzite") return wurtziteGeometry(m);
  if (kind === "diamond") return cubicGeometry(m, true);
  return cubicGeometry(m, false);
}
