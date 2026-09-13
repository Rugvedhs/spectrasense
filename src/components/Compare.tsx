import { MATERIALS } from "../data/materials";
import {
  baselinePredict,
  correctedPredict,
  OVERALL_BASELINE_METRICS,
  OVERALL_CORRECTED_METRICS,
  metricsByFamily,
  type Metrics,
} from "../lib/models";
import { fmtR2 } from "../lib/format";

const AXIS_MAX = 1.9;
const SIZE = 300;
const PAD = 34;

function scale(v: number) {
  return PAD + (v / AXIS_MAX) * (SIZE - PAD * 1.5);
}

function ScatterChart({ title, predict, color }: { title: string; predict: (m: (typeof MATERIALS)[number]) => number; color: string }) {
  return (
    <div className="scatter">
      <div className="scatter__title">{title}</div>
      <svg viewBox={`0 0 ${SIZE} ${SIZE}`} role="img" aria-label={`${title}: predicted versus reference band gap scatter plot`}>
        <line x1={scale(0)} y1={SIZE - scale(0)} x2={scale(AXIS_MAX)} y2={SIZE - scale(AXIS_MAX)} stroke="var(--border-strong)" strokeDasharray="4 4" />
        <line x1={scale(0)} y1={SIZE - PAD} x2={scale(AXIS_MAX)} y2={SIZE - PAD} stroke="var(--border)" />
        <line x1={PAD} y1={SIZE - scale(0)} x2={PAD} y2={SIZE - scale(AXIS_MAX)} stroke="var(--border)" />
        {[0, 0.5, 1, 1.5].map((t) => (
          <g key={t} fontFamily="IBM Plex Mono, monospace" fontSize="9.5" fill="var(--text-faint)">
            <text x={scale(t)} y={SIZE - PAD + 14} textAnchor="middle">{t}</text>
            <text x={PAD - 8} y={SIZE - scale(t) + 3} textAnchor="end">{t}</text>
          </g>
        ))}
        <text x={SIZE / 2} y={SIZE - 4} textAnchor="middle" fontFamily="IBM Plex Mono, monospace" fontSize="10" fill="var(--text-muted)">reference (eV)</text>
        <text x={10} y={SIZE / 2} textAnchor="middle" fontFamily="IBM Plex Mono, monospace" fontSize="10" fill="var(--text-muted)" transform={`rotate(-90 10 ${SIZE / 2})`}>predicted (eV)</text>

        {MATERIALS.map((m) => (
          <circle key={m.id} cx={scale(m.actualEgEv)} cy={SIZE - scale(predict(m))} r="4.5" fill={color} fillOpacity="0.85" />
        ))}
      </svg>
    </div>
  );
}

function MetricsRow({ label, m }: { label: string; m: Metrics }) {
  const r2Display = fmtR2(m.r2);
  return (
    <tr>
      <td>{label}</td>
      <td>{m.mae.toFixed(3)}</td>
      <td>{m.rmse.toFixed(3)}</td>
      <td className={m.r2 !== null && m.r2 < 0 ? "neg" : ""}>{r2Display}</td>
      <td>{m.n}</td>
    </tr>
  );
}

export default function Compare() {
  const baselineByFamily = metricsByFamily(MATERIALS, baselinePredict);
  const correctedByFamily = metricsByFamily(MATERIALS, correctedPredict);
  const families = Object.keys(baselineByFamily);

  return (
    <div className="view view--compare">
      <section className="compare-headline">
        <h1>Baseline vs. corrected — the actual result</h1>
        <p>
          Composition-only features reach <strong>R&sup2; = {fmtR2(OVERALL_BASELINE_METRICS.r2)}</strong> against
          the {OVERALL_BASELINE_METRICS.n}-material reference set (MAE {OVERALL_BASELINE_METRICS.mae.toFixed(2)} eV).
          Adding structural information — coordination environment and bonding topology — brings that to{" "}
          <strong>R&sup2; = {fmtR2(OVERALL_CORRECTED_METRICS.r2)}</strong> (MAE {OVERALL_CORRECTED_METRICS.mae.toFixed(2)} eV).
          Figures below are computed live from the current model outputs, not hardcoded.
        </p>
      </section>

      <section className="scatter-row">
        <ScatterChart title="Baseline · composition-only" predict={baselinePredict} color="var(--warn)" />
        <ScatterChart title="Corrected · structure-aware" predict={correctedPredict} color="var(--accent)" />
      </section>

      <section className="metrics-table">
        <h2>Metrics by material family</h2>
        <div className="table-wrap">
          <table>
            <thead>
              <tr><th>Model / family</th><th>MAE (eV)</th><th>RMSE (eV)</th><th>R&sup2;</th><th>n</th></tr>
            </thead>
            <tbody>
              <MetricsRow label="Baseline — overall" m={OVERALL_BASELINE_METRICS} />
              <MetricsRow label="Corrected — overall" m={OVERALL_CORRECTED_METRICS} />
              {families.map((f) => (
                <MetricsRow key={`b-${f}`} label={`Baseline — ${f}`} m={baselineByFamily[f]} />
              ))}
              {families.map((f) => (
                <MetricsRow key={`c-${f}`} label={`Corrected — ${f}`} m={correctedByFamily[f]} />
              ))}
            </tbody>
          </table>
        </div>
        <p className="metrics-table__note">
          The rocksalt Pb-chalcogenide family carries the baseline's worst error — composition-only descriptors have no
          way to encode the relativistic band-inversion physics that suppresses their gap, which is exactly the
          structure-sensitivity this comparison was built to surface. R&sup2; reads <code>n/a</code> for families with
          only one member: with no spread in the reference values there is no variance for R&sup2; to explain, so MAE/RMSE
          are the honest metrics there, not a fabricated 0 or 1.
        </p>
      </section>
    </div>
  );
}
