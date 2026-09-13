export function fmtR2(v: number | null, digits = 2): string {
  return v === null ? "n/a" : v.toFixed(digits);
}

export function fmtEv(v: number, digits = 2): string {
  return `${v.toFixed(digits)} eV`;
}

export function fmtNm(v: number, digits = 0): string {
  return `${v.toFixed(digits)} nm`;
}

export const TARGET_MIN_EV = 0.5;
export const TARGET_MAX_EV = 1.5;

export function inTargetWindow(eg: number): boolean {
  return eg >= TARGET_MIN_EV && eg <= TARGET_MAX_EV;
}

export function sensingContext(eg: number): string {
  if (eg < 0.2) return "deep short-wave IR territory — beyond typical NIR photodetector design, closer to thermal-imaging chemistries.";
  if (eg < 0.5) return "short-wave-IR (SWIR) range — matches extended-InGaAs and PbS/PbSe colloidal-QD detector designs.";
  if (eg <= 1.5) return "squarely inside the NIR sensing target window — the same range as InGaAs telecom photodiodes and CdTe-family absorbers.";
  if (eg <= 1.8) return "just past the target window, at the NIR/visible boundary — still relevant as an upper-bound reference point.";
  return "visible-range gap — outside the NIR sensing scope this sandbox is built around.";
}
