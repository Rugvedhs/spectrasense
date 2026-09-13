import { useState } from "react";
import TopNav, { type Mode } from "./components/TopNav";
import Explore from "./components/Explore";
import Sandbox from "./components/Sandbox";
import Compare from "./components/Compare";
import About from "./components/About";

export default function App() {
  const [mode, setMode] = useState<Mode>("explore");
  const [selectedMaterialId, setSelectedMaterialId] = useState<string | null>(null);

  const goToSandbox = (id: string) => {
    setSelectedMaterialId(id);
    setMode("sandbox");
  };

  return (
    <div className="app">
      <TopNav mode={mode} onChange={setMode} />
      <main className="app__main">
        {mode === "explore" && <Explore onSelect={goToSandbox} />}
        {mode === "sandbox" && <Sandbox initialMaterialId={selectedMaterialId} />}
        {mode === "compare" && <Compare />}
        {mode === "about" && <About />}
      </main>
      <footer className="app__footer">
        <span>SpectraSense — interactive NIR materials-discovery sandbox</span>
        <span>Data: Materials Project–style reference values · Preview build, simulated models</span>
      </footer>
    </div>
  );
}
