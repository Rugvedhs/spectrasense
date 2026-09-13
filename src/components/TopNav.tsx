export type Mode = "explore" | "sandbox" | "compare" | "about";

const TABS: { id: Mode; label: string }[] = [
  { id: "explore", label: "Explore" },
  { id: "sandbox", label: "Sandbox" },
  { id: "compare", label: "Compare" },
  { id: "about", label: "About" },
];

interface Props {
  mode: Mode;
  onChange: (m: Mode) => void;
}

export default function TopNav({ mode, onChange }: Props) {
  return (
    <header className="top-nav">
      <div className="top-nav__brand">
        <span className="top-nav__mark" aria-hidden="true" />
        <span className="top-nav__name">SpectraSense</span>
        <span className="top-nav__badge">preview build · simulated models</span>
      </div>
      <nav className="top-nav__tabs" aria-label="View mode">
        {TABS.map((t) => (
          <button
            key={t.id}
            className={`top-nav__tab ${mode === t.id ? "is-active" : ""}`}
            onClick={() => onChange(t.id)}
            aria-current={mode === t.id ? "page" : undefined}
          >
            {t.label}
          </button>
        ))}
      </nav>
    </header>
  );
}
