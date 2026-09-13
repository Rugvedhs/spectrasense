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
  if (eg < 0.2) return "this sits in deep short-wave IR territory, beyond typical NIR photodetector design and closer to thermal-imaging chemistries.";
  if (eg < 0.5) return "this falls in the short-wave-IR (SWIR) range, matching extended-InGaAs and PbS/PbSe colloidal-QD detector designs.";
  if (eg <= 1.5) return "this lands squarely inside the NIR sensing target window, the same range used by InGaAs telecom photodiodes and CdTe-family absorbers.";
  if (eg <= 1.8) return "this sits just past the target window at the NIR/visible boundary, still useful as an upper-bound reference point.";
  return "this has a visible-range gap, outside the scope this sandbox is built around.";
}
