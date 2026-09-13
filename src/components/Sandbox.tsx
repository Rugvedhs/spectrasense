import { useEffect, useMemo, useState } from "react";
import { MATERIALS, MATERIAL_BY_ID } from "../data/materials";
import { ALLOY_SYSTEMS, alloyVirtualMaterial } from "../data/alloys";
import StructureViewer from "./StructureViewer";
import PredictionPanel from "./PredictionPanel";
import { ELEMENTS } from "../data/elements";

interface Props {
  initialMaterialId: string | null;
}

export default function Sandbox({ initialMaterialId }: Props) {
  const [pickMode, setPickMode] = useState<"curated" | "alloy">("curated");
  const [curatedId, setCuratedId] = useState<string>(initialMaterialId ?? MATERIALS[0].id);
  const [alloySystemId, setAlloySystemId] = useState<string>(ALLOY_SYSTEMS[0].id);
  const [x, setX] = useState(0.47);

  useEffect(() => {
    if (initialMaterialId) {
      setPickMode("curated");
      setCuratedId(initialMaterialId);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [initialMaterialId]);

  const alloySystem = ALLOY_SYSTEMS.find((s) => s.id === alloySystemId)!;
  const material = useMemo(() => {
    return pickMode === "curated" ? MATERIAL_BY_ID[curatedId] : alloyVirtualMaterial(alloySystem, x);
  }, [pickMode, curatedId, alloySystem, x]);

  const legendSymbols = useMemo(() => {
    const set = new Map<string, string>();
    for (const c of material.cations) set.set(c.symbol, ELEMENTS[c.symbol]?.color ?? "#999");
    for (const a of material.anions) set.set(a.symbol, ELEMENTS[a.symbol]?.color ?? "#999");
    return Array.from(set.entries());
  }, [material]);

  return (
    <div className="view view--sandbox">
      <div className="sandbox__controls">
        <div className="control-group">
          <span className="control-group__label">Pick a starting point</span>
          <div className="segmented">
            <button className={pickMode === "curated" ? "is-active" : ""} onClick={() => setPickMode("curated")}>Curated material</button>
            <button className={pickMode === "alloy" ? "is-active" : ""} onClick={() => setPickMode("alloy")}>Alloy composition</button>
          </div>
        </div>

        {pickMode === "curated" ? (
          <div className="control-group">
            <label className="control-group__label" htmlFor="material-select">Base composition / structure</label>
            <select id="material-select" value={curatedId} onChange={(e) => setCuratedId(e.target.value)}>
              {MATERIALS.map((m) => (
                <option key={m.id} value={m.id}>{m.displayName} — {m.formula}</option>
              ))}
            </select>
          </div>
        ) : (
          <>
            <div className="control-group">
              <label className="control-group__label" htmlFor="alloy-select">Alloy system</label>
              <select id="alloy-select" value={alloySystemId} onChange={(e) => setAlloySystemId(e.target.value)}>
                {ALLOY_SYSTEMS.map((s) => (
                  <option key={s.id} value={s.id}>{s.label}</option>
                ))}
              </select>
              <p className="control-group__hint">{alloySystem.description}</p>
            </div>
            <div className="control-group">
              <label className="control-group__label" htmlFor="x-slider">
                Composition x = {x.toFixed(2)}
              </label>
              <input id="x-slider" type="range" min={0} max={1} step={0.01} value={x} onChange={(e) => setX(parseFloat(e.target.value))} />
              <div className="control-group__endpoints">
                <span>x = 0 ({MATERIAL_BY_ID[alloySystem.endpointAId].formula})</span>
                <span>x = 1 ({MATERIAL_BY_ID[alloySystem.endpointBId].formula})</span>
              </div>
            </div>
          </>
        )}
      </div>

      <div className="sandbox__stage">
        <div className="structure-panel">
          <StructureViewer material={material} />
          <div className="structure-panel__legend">
            {legendSymbols.map(([symbol, color]) => (
              <span key={symbol} className="legend-chip">
                <span className="legend-chip__swatch" style={{ background: color }} />
                {symbol} — {ELEMENTS[symbol]?.name ?? symbol}
              </span>
            ))}
          </div>
          <div className="structure-panel__meta">
            <span>{material.structure}</span>
            <span>{material.family}</span>
            <span>a &asymp; {material.latticeConstantAng.toFixed(2)} &Aring;</span>
          </div>
        </div>

        <PredictionPanel material={material} />
      </div>
    </div>
  );
}
