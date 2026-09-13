import { MATERIALS } from "../data/materials";
import { baselinePredict, correctedPredict } from "../lib/models";
import { fmtEv } from "../lib/format";
import SpectrumAxis from "./SpectrumAxis";

interface Props {
  onSelect: (id: string) => void;
}

export default function Explore({ onSelect }: Props) {
  return (
    <div className="view view--explore">
      <section className="landing">
        <h1>Finding the next NIR sensing material, computationally</h1>
        <p className="landing__lede">
          Near-infrared sensing — allergen detection, produce ripeness, moisture and material
          identification — depends on a detector material whose band gap sits in the right window,
          roughly <strong>0.5–1.5 eV</strong>. Today that search is mostly a literature exercise: look up
          what's been made before, borrow its band gap, hope it transfers. SpectraSense asks whether a
          model trained on open computational materials data can do better — and whether it actually
          needs to know a material's <em>structure</em>, not just its <em>composition</em>, to get there.
        </p>
        <p className="landing__lede">
          The curated set below spans exactly the chemistries a physical NIR sensor build (like a
          PbS/PbSe-based allergen detector) would draw from: lead- and tin-chalcogenide rocksalts,
          III-V and II-VI zincblende semiconductors, and the group-IV control cases, Si and Ge.
        </p>
      </section>

      <section className="gallery">
        <div className="gallery__head">
          <h2>Explore mode</h2>
          <span className="gallery__count">{MATERIALS.length} curated materials</span>
        </div>
        <div className="gallery__grid">
          {MATERIALS.map((m) => {
            const baseline = baselinePredict(m);
            const corrected = correctedPredict(m);
            return (
              <button key={m.id} className="material-card" onClick={() => onSelect(m.id)}>
                <div className="material-card__top">
                  <span className="material-card__formula">{m.formula}</span>
                  <span className={`structure-badge structure-badge--${m.structure}`}>{m.structure}</span>
                </div>
                <div className="material-card__name">{m.displayName}</div>
                <div className="material-card__family">{m.family}</div>

                <SpectrumAxis
                  compact
                  markers={[
                    { label: "act", eg: m.actualEgEv, color: "var(--text-muted)" },
                    { label: "base", eg: baseline, color: "var(--warn)" },
                    { label: "corr", eg: corrected, color: "var(--accent)" },
                  ]}
                />

                <div className="material-card__stats">
                  <div><span>actual</span><strong>{fmtEv(m.actualEgEv)}</strong></div>
                  <div><span>baseline</span><strong>{fmtEv(baseline)}</strong></div>
                  <div><span>corrected</span><strong>{fmtEv(corrected)}</strong></div>
                </div>
              </button>
            );
          })}
        </div>
      </section>
    </div>
  );
}
