import type { Material } from "../data/materials";
import {
  baselinePredict,
  correctedPredict,
  baselineUncertaintyEv,
  correctedUncertaintyEv,
  evToNm,
} from "../lib/models";
import { fmtEv, fmtNm, sensingContext, TARGET_MIN_EV, TARGET_MAX_EV } from "../lib/format";

interface Props {
  material: Material;
}

const TRACK_MIN = 0;
const TRACK_MAX = 2.2;
const pct = (v: number) => `${(Math.min(Math.max(v, TRACK_MIN), TRACK_MAX) / TRACK_MAX) * 100}%`;

function ModelRow({
  label,
  predicted,
  uncertainty,
  actual,
  tone,
}: {
  label: string;
  predicted: number;
  uncertainty: number;
  actual: number;
  tone: "baseline" | "corrected";
}) {
  const lo = predicted - uncertainty;
  const hi = predicted + uncertainty;
  const errAbs = Math.abs(predicted - actual);
  return (
    <div className={`model-row model-row--${tone}`}>
      <div className="model-row__head">
        <span className="model-row__label">{label}</span>
        <span className="model-row__value">
          {fmtEv(predicted)} <span className="model-row__nm">· {fmtNm(evToNm(predicted))}</span>
        </span>
      </div>
      <div className="model-row__track">
        <div className="model-row__window" style={{ left: pct(TARGET_MIN_EV), width: `calc(${pct(TARGET_MAX_EV)} - ${pct(TARGET_MIN_EV)})` }} />
        <div className="model-row__band" style={{ left: pct(lo), width: `calc(${pct(hi)} - ${pct(lo)})` }} />
        <div className="model-row__actual" style={{ left: pct(actual) }} title={`reference: ${fmtEv(actual)}`} />
        <div className="model-row__point" style={{ left: pct(predicted) }} />
      </div>
      <div className="model-row__foot">
        <span>&plusmn; {fmtEv(uncertainty)} uncertainty</span>
        <span>{fmtEv(errAbs)} from reference</span>
      </div>
    </div>
  );
}

export default function PredictionPanel({ material }: Props) {
  const baseline = baselinePredict(material);
  const corrected = correctedPredict(material);
  const baseUnc = baselineUncertaintyEv(material);
  const corrUnc = correctedUncertaintyEv(material);

  return (
    <div className="prediction-panel">
      <ModelRow label="Baseline · composition-only" predicted={baseline} uncertainty={baseUnc} actual={material.actualEgEv} tone="baseline" />
      <ModelRow label="Corrected · structure-aware" predicted={corrected} uncertainty={corrUnc} actual={material.actualEgEv} tone="corrected" />

      <div className="reference-row">
        <span>Reference value (Materials Project-style)</span>
        <span className="reference-row__value">{fmtEv(material.actualEgEv)} · {fmtNm(evToNm(material.actualEgEv))}</span>
      </div>

      <div className="explainer">
        <span className="explainer__label">Why this might work for NIR sensing</span>
        <p>
          At {fmtEv(material.actualEgEv)} ({fmtNm(evToNm(material.actualEgEv))}), {sensingContext(material.actualEgEv)}
        </p>
      </div>
    </div>
  );
}
