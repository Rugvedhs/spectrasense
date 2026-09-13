import { ELEMENTS } from "../data/elements";
import { MATERIALS, type Material } from "../data/materials";

/**
 * IMPORTANT — read before wiring up real training data:
 *
 * Phase 1-3 of the build roadmap (Materials Project data pull, baseline
 * composition-only regressor, structure-aware corrected model) have not run
 * yet — this session was built without an MP API key. The two "models" below
 * are deterministic stand-ins that encode the same *qualitative* hypothesis
 * the real pipeline is meant to test: composition-only descriptors track
 * band gap reasonably for ordinary zincblende/diamond semiconductors but miss
 * the relativistic band-gap suppression in rocksalt Pb-chalcogenides, while a
 * structure-aware model does not. Replace `baselinePredict` / `correctedPredict`
 * with real model inference once Phase 2/3 land; nothing else in the app
 * depends on how the numbers are produced.
 */

export interface CompositionFeatures {
  meanElectronegativityDiff: number;
  meanValenceElectrons: number;
  meanPeriod: number;
  meanCovalentRadiusPm: number;
}

export function compositionFeatures(m: Material): CompositionFeatures {
  const cationEN = weightedMean(m.cations, (s) => ELEMENTS[s].electronegativity);
  const anionEN = weightedMean(m.anions, (s) => ELEMENTS[s].electronegativity);
  const allSites = [...m.cations, ...m.anions];
  return {
    meanElectronegativityDiff: Math.abs(anionEN - cationEN),
    meanValenceElectrons: weightedMean(allSites, (s) => ELEMENTS[s].valenceElectrons),
    meanPeriod: weightedMean(allSites, (s) => ELEMENTS[s].period),
    meanCovalentRadiusPm: weightedMean(allSites, (s) => ELEMENTS[s].covalentRadiusPm),
  };
}

function weightedMean(components: { symbol: string; frac: number }[], get: (s: string) => number): number {
  const totalFrac = components.reduce((a, c) => a + c.frac, 0) || 1;
  return components.reduce((a, c) => a + get(c.symbol) * c.frac, 0) / totalFrac;
}

// Deterministic pseudo-random in [-1, 1], seeded by a string — stands in for
// per-material model residual so the same material always yields the same
// "prediction" across renders without needing real inference.
function seededUnit(seed: string): number {
  let h = 2166136261;
  for (let i = 0; i < seed.length; i++) {
    h ^= seed.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  const t = (h >>> 0) / 4294967295;
  return t * 2 - 1;
}

const FAMILY_BASELINE_BIAS: Record<Material["family"], number> = {
  "III-V zincblende": 0.92,
  "II-VI zincblende": 0.90,
  "II-VI wurtzite": 0.90,
  "Group IV (diamond)": 0.88,
  "IV-VI rocksalt (Sn-chalcogenide)": 1.8,
  "IV-VI rocksalt (Pb-chalcogenide)": 2.6,
};

export function baselinePredict(m: Material): number {
  const bias = FAMILY_BASELINE_BIAS[m.family];
  const jitter = seededUnit(m.id + ":baseline") * 0.08 * m.actualEgEv;
  return Math.max(0.02, m.actualEgEv * bias + jitter);
}

export function correctedPredict(m: Material): number {
  const jitter = seededUnit(m.id + ":corrected") * 0.09 * m.actualEgEv;
  return Math.max(0.02, m.actualEgEv + jitter);
}

// A rough calibrated "uncertainty" band per model — wider where the model is
// systematically worse, used for the confidence indicator in the UI.
export function baselineUncertaintyEv(m: Material): number {
  const bias = FAMILY_BASELINE_BIAS[m.family];
  return Math.abs(bias - 1) * m.actualEgEv * 0.5 + 0.03;
}

export function correctedUncertaintyEv(m: Material): number {
  return 0.09 * m.actualEgEv + 0.02;
}

export function evToNm(eg: number): number {
  return 1239.84 / eg;
}

export interface Metrics {
  mae: number;
  rmse: number;
  // R² is statistically undefined when there isn't enough spread in the
  // reference values to define a variance to explain (n < 2, or a
  // single-material family) — null rather than a fabricated 0/1 in that case.
  r2: number | null;
  n: number;
}

export function computeMetrics(materials: Material[], predict: (m: Material) => number): Metrics {
  const n = materials.length;
  const actual = materials.map((m) => m.actualEgEv);
  const pred = materials.map(predict);
  const meanActual = actual.reduce((a, b) => a + b, 0) / n;

  let sumAbs = 0;
  let sumSq = 0;
  let ssRes = 0;
  let ssTot = 0;
  for (let i = 0; i < n; i++) {
    const err = pred[i] - actual[i];
    sumAbs += Math.abs(err);
    sumSq += err * err;
    ssRes += err * err;
    ssTot += (actual[i] - meanActual) ** 2;
  }
  return {
    mae: sumAbs / n,
    rmse: Math.sqrt(sumSq / n),
    r2: n < 2 || ssTot === 0 ? null : 1 - ssRes / ssTot,
    n,
  };
}

export function metricsByFamily(materials: Material[], predict: (m: Material) => number): Record<string, Metrics> {
  const families = Array.from(new Set(materials.map((m) => m.family)));
  const out: Record<string, Metrics> = {};
  for (const f of families) {
    out[f] = computeMetrics(materials.filter((m) => m.family === f), predict);
  }
  return out;
}

export const OVERALL_BASELINE_METRICS = computeMetrics(MATERIALS, baselinePredict);
export const OVERALL_CORRECTED_METRICS = computeMetrics(MATERIALS, correctedPredict);
