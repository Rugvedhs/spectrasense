# SpectraSense

An interactive materials-discovery sandbox for NIR-active sensing materials. Pick or build a
composition/structure and see predicted band gap from two models: a composition-only baseline and
a structure-aware corrected model. Both are compared against a reference value, with the gap
reported honestly as MAE / RMSE / R² rather than a single flattering number.

## Current build status

This is the **Phase 4 (front-end sandbox)** slice of the project, built ahead of Phases 1–3
because no Materials Project API key was available yet. Concretely:

- The 14-material curated set in [`src/data/materials.ts`](src/data/materials.ts) is hand-entered
  from standard reference band-gap values, not a live `mp-api` pull.
- `baselinePredict` and `correctedPredict` in [`src/lib/models.ts`](src/lib/models.ts) are
  **deterministic stand-ins**, not trained models. They encode the same qualitative hypothesis the
  real pipeline is meant to test: composition-only descriptors track band gap reasonably for
  ordinary zincblende/diamond semiconductors but miss the relativistic band-gap suppression in
  rocksalt Pb-chalcogenides. That keeps the app honest about *why* the comparison looks the way it
  does, without claiming a result it hasn't earned yet.
- Every metric shown in the app (R², MAE, RMSE, scatter plots) is computed live from those
  functions, never hardcoded, so swapping in real models changes the numbers everywhere at once.

Wiring in real data and real models (Phases 1-3 of the roadmap below) does not require touching
the UI layer. `Explore`, `Sandbox`, `Compare`, and `About` all consume `materials.ts` and
`models.ts` through the same interfaces a trained pipeline would fill in.

## Running it

```bash
npm install
npm run dev
```

## Architecture

```
src/
  data/
    elements.ts    element reference data (electronegativity, radius, valence, render color)
    materials.ts   curated 14-material reference set + Material/StructureType/Family types
    alloys.ts      pseudo-binary alloy systems (InGaAs, PbSSe, CdSeTe) for the Sandbox slider
  lib/
    models.ts      composition featurization, baseline/corrected predictors, MAE/RMSE/R² metrics
    crystal.ts      unit-cell geometry generator (rocksalt, zincblende, wurtzite, diamond cubic)
    format.ts       eV <-> nm conversion, target-window logic, display formatting
  components/
    TopNav.tsx           Explore / Sandbox / Compare / About mode switcher
    Explore.tsx          landing explainer + curated-material gallery
    Sandbox.tsx          material/alloy picker, structure viewer, live predictions
    Compare.tsx          baseline-vs-corrected scatter plots + per-family metrics table
    About.tsx            methodology write-up, validation numbers, honest limitations
    StructureViewer.tsx  three.js unit-cell renderer (spheres, bonds, cell wireframe)
    PredictionPanel.tsx  per-model prediction card with uncertainty band + NIR-sensing context
    SpectrumAxis.tsx      wavelength/eV axis figure used in the gallery and elsewhere
```

**Why plain `three.js` instead of `@react-three/fiber`:** `@react-three/fiber`'s current release
peer-depends on `react <19.3`, which conflicts with this project's React 19.3. Driving `three`
directly from a `useEffect` avoids the dependency conflict entirely and keeps the render loop
fully inspectable.

## Data layer (target, once an MP API key is available)

- Source: Materials Project API (`mp-api` Python client), ~150,000+ computed materials.
- Target property: band gap, as the tractable proxy for NIR absorption-edge behavior.
- Scope: narrow-gap semiconductors (0.5–1.5 eV), II-VI/III-V chalcogenides, known NIR
  photodetector families (PbS, PbSe, InGaAs-adjacent, CdTe-family).
- Pipeline: query → clean → deduplicate polymorphs → train/validation/test split, stratified by
  chemical family to prevent leakage.

## Modeling layer (target)

- **Baseline:** Magpie-style composition-only descriptors → gradient-boosted regressor (XGBoost).
- **Corrected:** adds coordination number, bond-length statistics, space group, or (stretch goal) a
  graph neural network over the crystal graph (CGCNN/MEGNet-style).
- Same train/test split for both, evaluated side by side. The comparison *is* the result.

## Build roadmap

| Phase | Scope | Status |
|---|---|---|
| 1. Data pipeline | Materials Project API integration, filtering, cleaning, splitting | not started |
| 2. Baseline model | Composition-feature regressor + evaluation | not started |
| 3. Corrected model | Structure-aware model + evaluation, side-by-side comparison | not started |
| 4. Front-end sandbox | React app, interactive controls, structure visualization | **this build** |
| 5. Validation writeup | Technical report with honest error analysis | not started |
| 6. Deployment + polish | Hosting, README, repo cleanup, demo recording | README done; hosting pending |

## Tech stack

- **Data/ML (target):** Python, `mp-api`, pandas, scikit-learn/XGBoost, PyTorch Geometric.
- **Front-end (this build):** React 19 + TypeScript + Vite, `three.js` for the crystal-structure
  viewer, no CSS framework. Design tokens are hand-written in `src/styles/global.css`.
- **Deployment (target):** static/serverless hosting for the front end, a small API layer for live
  predictions once real models exist.

## Honest limitations

- The R²/MAE gap shown in **Compare** is illustrative of the hypothesis, not yet a trained-model
  result (see "Current build status" above).
- A full GNN structure encoder is a stretch goal; the documented fallback is the hand-engineered
  structural feature set already used by the placeholder `correctedPredict`.
- Scope is deliberately narrow: one target property (band gap), one physically motivated material
  set, rather than broad and shallow coverage.
