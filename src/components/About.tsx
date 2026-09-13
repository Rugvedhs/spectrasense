import { OVERALL_BASELINE_METRICS, OVERALL_CORRECTED_METRICS } from "../lib/models";
import { MATERIALS } from "../data/materials";
import { fmtR2 } from "../lib/format";

export default function About() {
  return (
    <div className="view view--about">
      <section>
        <h1>Methodology</h1>
        <p>
          SpectraSense predicts semiconductor band gap — used here as the tractable proxy for NIR
          absorption-edge behavior — from two feature sets computed over a curated,{" "}
          {MATERIALS.length}-material set spanning IV-VI rocksalt, III-V/II-VI zincblende, II-VI
          wurtzite, and group-IV diamond-cubic chemistries relevant to NIR photodetectors.
        </p>
      </section>

      <section>
        <h2>Data layer</h2>
        <p>
          Production data source: the <strong>Materials Project API</strong> (<code>mp-api</code>),
          filtered to compositions with band gap in the 0.5–1.5 eV target window plus reference
          compounds just outside it, deduplicated by polymorph, split by chemical family to prevent
          leakage between train and test sets.
        </p>
        <p className="callout-note">
          <strong>Current build status:</strong> this preview ships without a Materials Project API
          key, so the {MATERIALS.length} materials above are hand-curated literature values rather
          than a live API pull, and the two models below are deterministic stand-ins for the trained
          regressors described in the spec (see <code>src/lib/models.ts</code>). Swapping in a real
          <code>mp-api</code> pull and trained models does not require touching the UI layer.
        </p>
      </section>

      <section>
        <h2>Modeling layer</h2>
        <h3>Baseline — composition only</h3>
        <p>
          Magpie-style elemental statistics (mean electronegativity difference, mean valence electron
          count, mean period, mean covalent radius) fed into a regressor with no knowledge of crystal
          structure.
        </p>
        <h3>Corrected — structure-aware</h3>
        <p>
          Adds coordination environment and bonding topology (rocksalt vs. zincblende vs. wurtzite vs.
          diamond, and the resulting nearest-neighbor geometry) — the information the baseline cannot
          see.
        </p>
      </section>

      <section>
        <h2>Validation — current numbers</h2>
        <div className="about-metrics">
          <div>
            <span className="about-metrics__label">Baseline R&sup2;</span>
            <span className="about-metrics__value">{fmtR2(OVERALL_BASELINE_METRICS.r2)}</span>
          </div>
          <div>
            <span className="about-metrics__label">Corrected R&sup2;</span>
            <span className="about-metrics__value">{fmtR2(OVERALL_CORRECTED_METRICS.r2)}</span>
          </div>
          <div>
            <span className="about-metrics__label">Baseline MAE</span>
            <span className="about-metrics__value">{OVERALL_BASELINE_METRICS.mae.toFixed(2)} eV</span>
          </div>
          <div>
            <span className="about-metrics__label">Corrected MAE</span>
            <span className="about-metrics__value">{OVERALL_CORRECTED_METRICS.mae.toFixed(2)} eV</span>
          </div>
        </div>
        <p>
          Full breakdown by material family — including where the baseline fails worst — is on the{" "}
          <strong>Compare</strong> page, computed live from the same functions each time the app renders.
        </p>
      </section>

      <section>
        <h2>Honest limitations</h2>
        <ul>
          <li>The gap between baseline and corrected R&sup2; shown here is illustrative of the hypothesis, not yet a trained-model result — that requires Phase 1-3 of the build roadmap against real Materials Project data.</li>
          <li>A full graph-neural-network structure encoder (CGCNN/MEGNet-style) is a stretch goal; the documented fallback is a hand-engineered structural feature set (coordination number, bond-length statistics) feeding the same gradient-boosted regressor.</li>
          <li>Scope is deliberately narrow — one target property, one physically motivated material scope — rather than broad and shallow.</li>
        </ul>
      </section>
    </div>
  );
}
