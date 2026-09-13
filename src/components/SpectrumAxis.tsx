import { evToNm } from "../lib/models";

export interface SpectrumMarker {
  label: string;
  eg: number;
  color: string;
}

interface Props {
  markers: SpectrumMarker[];
  compact?: boolean;
}

const LAMBDA_MIN = 380;
const LAMBDA_MAX = 2600;
const X0 = 20;
const X1 = 620;

function xForNm(nm: number): number {
  const clamped = Math.min(Math.max(nm, LAMBDA_MIN), LAMBDA_MAX);
  return X0 + ((clamped - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN)) * (X1 - X0);
}

const visibleEndX = xForNm(700);
const targetStartX = xForNm(evToNm(1.5));
const targetEndX = xForNm(evToNm(0.5));

export default function SpectrumAxis({ markers, compact }: Props) {
  const h = compact ? 62 : 96;
  const bandY = compact ? 14 : 20;
  const bandH = compact ? 12 : 16;

  return (
    <svg viewBox={`0 0 640 ${h}`} className="spectrum-axis" role="img" aria-label="Wavelength axis from 380 to 2600 nanometers with the 0.5 to 1.5 electronvolt target window highlighted">
      <rect x={X0} y={bandY} width={visibleEndX - X0} height={bandH} fill="var(--vis-swatch)" opacity="0.55" rx="2" />
      <rect x={visibleEndX} y={bandY} width={X1 - visibleEndX} height={bandH} fill="var(--nir-fill)" rx="2" />
      <rect x={targetStartX} y={bandY} width={targetEndX - targetStartX} height={bandH} fill="var(--nir-target-fill)" stroke="var(--accent)" strokeWidth="1" strokeDasharray="3 3" rx="2" />
      <line x1={X0} y1={bandY + bandH + 10} x2={X1} y2={bandY + bandH + 10} stroke="var(--border-strong)" strokeWidth="1" />

      {!compact && (
        <g fontFamily="IBM Plex Mono, monospace" fontSize="9.5" fill="var(--text-faint)">
          <text x={X0} y={bandY + bandH + 24}>380</text>
          <text x={targetStartX} y={bandY + bandH + 24} textAnchor="middle">827</text>
          <text x={targetEndX} y={bandY + bandH + 24} textAnchor="middle">2480</text>
          <text x={X1} y={bandY + bandH + 24} textAnchor="end">2600 nm</text>
        </g>
      )}

      {markers.map((m, i) => {
        const x = xForNm(evToNm(m.eg));
        const y = bandY - 6;
        return (
          <g key={i}>
            <line x1={x} y1={y} x2={x} y2={bandY + bandH} stroke={m.color} strokeWidth="2" />
            <circle cx={x} cy={y} r={compact ? 3 : 4} fill={m.color} />
            {!compact && (
              <text x={x} y={y - 8} textAnchor="middle" fontFamily="IBM Plex Mono, monospace" fontSize="10" fill={m.color} fontWeight={600}>
                {m.label}
              </text>
            )}
          </g>
        );
      })}
    </svg>
  );
}
