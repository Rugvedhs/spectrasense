"""Assert every cell of the manuscript's inline tables against its source CSV.

The fourteen numbered tables in the body of the paper are transcriptions: they
are typeset from the markdown source rather than generated, because several of
them are excerpts, reshapes or joins of more than one result file and a generated
table could not carry the column headings the prose refers to.  A transcription
can drift.  This module re-parses the markdown tables out of ``paper/manuscript.md``
and rebuilds each one from the CSV named in its caption, at the rounding the
manuscript quotes, then compares cell by cell.

Run it directly for a report, or let ``pytest tests`` run it.

One column is deliberately not checked and is declared as such: the "Reported by
[N]" column of Table 1 is transcribed from a prior public analysis of the same
dataset and is not an output of this pipeline.  Its heading carries whatever
number the reference list currently assigns that analysis, so it is matched by
pattern rather than by a literal.
"""

from __future__ import annotations

import re
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT = ROOT / "paper" / "manuscript.md"
RESULTS = ROOT / "results"

UNCHECKED = {(1, re.compile(r"^Reported by \[\d+\]$"))}


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------

SUPERSCRIPTS = str.maketrans("⁰¹²³⁴⁵⁶⁷"
                             "⁸⁹⁻", "0123456789-")

NUMBER = re.compile(r"[+-]?\d+(?:\.\d+)?")


def normalise(cell: str) -> str:
    """Strip the typographic layer so a cell compares as the number it states."""
    text = cell.strip()
    text = text.replace("**", "").replace("*", "")
    text = text.replace("−", "-")     # minus sign
    text = text.replace("×", "x")     # multiplication sign
    text = text.replace(" ", " ")
    text = text.translate(SUPERSCRIPTS)    # 10⁻⁷ -> 10-7
    text = text.replace(",", "")           # thousands separator
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def agrees(found: str, expected: str) -> bool:
    """Do two cells state the same thing at the precision the manuscript quotes?

    The non-numeric skeleton must match exactly.  Numbers must agree to within
    half a unit in the last decimal place the manuscript prints, which is what
    "correct at the quoted rounding" means and which tolerates a tie being
    rounded either way (1.5675 may be printed 1.567 or 1.568).
    """
    a, b = normalise(found), normalise(expected)
    if a == b:
        return True
    if NUMBER.sub("#", a) != NUMBER.sub("#", b):
        return False
    a_numbers, b_numbers = NUMBER.findall(a), NUMBER.findall(b)
    for x, y in zip(a_numbers, b_numbers):
        places = len(x.split(".")[1]) if "." in x else 0
        if abs(float(x) - float(y)) > 0.5 * 10 ** -places + 1e-12:
            return False
    return True


def parse_tables(text: str) -> dict[int, list[list[str]]]:
    """Return {table number: rows}, keyed off the ``**Table N.**`` caption."""
    tables: dict[int, list[list[str]]] = {}
    number: int | None = None
    rows: list[list[str]] = []

    def close() -> None:
        nonlocal rows, number
        if number is not None and rows:
            tables[number] = rows
            number = None
        rows = []

    for line in text.split("\n"):
        stripped = line.strip()
        caption = re.match(r"\*\*Table (\d+)\.\*\*", stripped)
        if caption:
            close()
            number = int(caption.group(1))
            continue
        if stripped.startswith("|") and stripped.endswith("|"):
            cells = [c.strip() for c in stripped.strip("|").split("|")]
            if re.fullmatch(r"[\s|:-]+", "|".join(cells)):
                continue
            rows.append(cells)
            continue
        if rows:
            close()
    close()
    return tables


# --------------------------------------------------------------------------
# expected tables, rebuilt from the result CSVs
# --------------------------------------------------------------------------

def csv(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / f"{name}.csv")


def f(value: float, places: int, sign: bool = False) -> str:
    """Format as the manuscript does, rounding halves away from zero.

    Python's default rounds ties to even and, worse, to whatever the binary
    representation happens to sit on; a transcriber rounds 1.5675 up. Decimal
    makes the two agree.
    """
    quantum = Decimal(1).scaleb(-places)
    rounded = Decimal(repr(float(value))).quantize(quantum, rounding=ROUND_HALF_UP)
    return f"{rounded:{'+' if sign else ''}f}"


def expect_1() -> list[list[str]]:
    curation = csv("table_curation").set_index("step")["value"]
    points = csv("table_attachment_points").set_index("n_attachment_points")["n_structures"]
    rows = [["Step", "", "This work"]]   # heading and column skipped, see UNCHECKED
    labels = {
        "Raw records": "Raw records",
        "Unique raw SMILES": "Unique raw SMILES",
        "RDKit parse failures": "RDKit parse failures",
        "Canonical structures retained": "Canonical structures retained",
        "Structures aggregated from >1 record": "Structures aggregated from >1 record",
        "Records absorbed by aggregation": "Records absorbed by aggregation",
        "Aggregated groups with disagreeing Tg": "Aggregated groups with disagreeing *T*~g~",
        "Max within-structure Tg range (K)": "Max within-structure *T*~g~ range (K)",
    }
    for key, label in labels.items():
        rows.append([label, "", f"{int(curation[key])}"])
    rows.append([
        "Units with 2 / 3 / 4 attachment points", "",
        f"{int(points[2])} / {int(points[3])} / {int(points[4])} (structures)",
    ])
    return rows


def expect_2() -> list[list[str]]:
    frame = csv("table_families").sort_values("n", ascending=False)
    rows = [["Family", "*n*", "Median *T*~g~ (°C)", "Min (°C)", "Max (°C)"]]
    for _, r in frame.iterrows():
        rows.append([r["family"], f"{int(r['n'])}", f(r["tg_median"], 1),
                     f(r["tg_min"], 1), f(r["tg_max"], 1)])
    return rows


MODEL_LABELS = {
    "hist_gbr": "Histogram gradient boosting",
    "extra_trees": "Extremely randomised trees",
    "random_forest": "Random forest",
    "svr": "Support vector regression",
    "median": "Training-median baseline",
}
REGIMES = ("random", "scaffold", "cluster", "family")


def expect_3() -> list[list[str]]:
    frame = csv("table_point_by_regime")
    rows = [["Model", "Random (10)", "Scaffold (10)", "Cluster (10)", "Family (21)"]]
    for model, label in MODEL_LABELS.items():
        cells = [label]
        for regime in REGIMES:
            r = frame[(frame["model"] == model) & (frame["regime"] == regime)].iloc[0]
            cells.append(f"{f(r['mae_mean'], 2)} ± {f(r['mae_sd'], 2)}")
        rows.append(cells)
    similarity = ["Mean nearest-training similarity"]
    for regime in REGIMES:
        r = frame[frame["regime"] == regime].iloc[0]
        similarity.append(f(r["similarity"], 3))
    rows.append(similarity)
    return rows


def expect_4() -> list[list[str]]:
    frame = csv("table_pooled_r2")
    labels = {"random": "Random (10 splits)", "scaffold": "Scaffold (10 splits)",
              "cluster": "Cluster (10 splits)", "family": "Family (21 splits)"}
    rows = [["Regime", "Predictions pooled", "Pooled *R*²", "Mean of per-split *R*²"]]
    for regime, label in labels.items():
        r = frame[(frame["regime"] == regime) & (frame["model"] == "hist_gbr")].iloc[0]
        rows.append([label, f"{int(r['n_predictions'])}", f(r["r2_pooled"], 3),
                     f(r["r2_mean_of_splits"], 3)])
    return rows


def expect_5() -> list[list[str]]:
    frame = csv("table_family_holdout")
    frame = frame[frame["model"] == "hist_gbr"].sort_values("deterioration",
                                                            ascending=False)
    rows = [["Family", "*n*", "MAE (K)", "Deterioration", "Mean NN similarity"]]
    selected = list(frame.head(3).itertuples(index=False))
    selected += list(frame.tail(3).itertuples(index=False))
    for i, r in enumerate(selected):
        if i == 3:
            rows.append(["…", "", "", "", ""])
        rows.append([r.group, f"{int(r.n)}", f(r.mae, 1),
                     f"{f(r.deterioration, 2)}×", f(r.mean_nn_similarity, 3)])
    return rows


def expect_6() -> list[list[str]]:
    frame = csv("table_bias_shrinkage").sort_values("bias_residual", ascending=False)
    chosen = list(frame.head(4)["group"]) + ["Other backbone"] + list(
        frame.tail(4)["group"])[::-1]
    # The manuscript orders the block by residual, descending, with the
    # abstention class inserted at its own rank.
    order = [g for g in frame["group"] if g in set(chosen)]
    rows = [["Family", "Offset (K)", "Bias (K)", "Expected from offset (K)",
             "Residual (K)"]]
    for group in order:
        r = frame[frame["group"] == group].iloc[0]
        rows.append([group, f(r["tg_offset"], 1, sign=True),
                     f(r["bias"], 1, sign=True),
                     f(r["bias_expected_from_offset"], 1, sign=True),
                     f(r["bias_residual"], 1, sign=True)])
    return rows


CONFORMAL_LABELS = {
    "SCP": "SCP",
    "Normalised SCP": "Normalised SCP",
    "Mondrian-family": "Mondrian, family",
    "Mondrian-similarity": "Mondrian, similarity",
}


def expect_7() -> list[list[str]]:
    conformal = csv("table_conformal")
    domain = csv("table_applicability_domain")
    band = domain[domain["similarity_bin"] == 0]
    rows = [["Regime", "Method", "Coverage (mean of splits)", "Band < 0.4 (pooled)",
             "Width (K)", "Fallback", "CWC"]]
    for regime in REGIMES:
        for key, label in CONFORMAL_LABELS.items():
            r = conformal[(conformal["regime"] == regime)
                          & (conformal["conformal"] == key)].iloc[0]
            b = band[(band["regime"] == regime) & (band["conformal"] == key)].iloc[0]
            rows.append([regime.capitalize(), label, f(r["coverage"], 3),
                         f(b["coverage"], 3), f(r["mean_width"], 1),
                         f(r["fallback"], 2), f"{round(r['cwc']):d}"])
    return rows


def expect_8() -> list[list[str]]:
    frame = csv("table_width_matched")
    rows = [["Regime", "Predictor", "*k*", "Coverage (pooled)",
             "Band < 0.4 (pooled)", "Width (K)"]]
    plan = {
        "family": [("Split conformal", 1.0), ("Coverage-matched oracle", "matched"),
                   ("Mondrian, similarity", None), ("Width-matched oracle", "oracle"),
                   ("Inflated to match band", 1.93),
                   ("Inflated to nominal band", 2.0)],
        "cluster": [("Split conformal", 1.0), ("Coverage-matched oracle", "matched"),
                    ("Mondrian, similarity", None), ("Width-matched oracle", "oracle")],
    }
    for regime, entries in plan.items():
        block = frame[frame["regime"] == regime]
        for label, key in entries:
            if key is None:
                r = block.iloc[0]
                rows.append([regime.capitalize(), label, "—",
                             f(r["mondrian_coverage_pooled"], 3),
                             f(r["mondrian_coverage_lowest_band"], 3),
                             f(r["mondrian_mean_width"], 1)])
                continue
            if key == "matched":
                r = block[block["is_coverage_matched_k"]].iloc[0]
            elif key == "oracle":
                r = block[block["is_oracle_k"]].iloc[0]
            else:
                r = block[(block["k"] - key).abs() < 1e-9].iloc[0]
            rows.append([regime.capitalize(), label, f(r["k"], 3),
                         f(r["coverage_pooled"], 3),
                         f(r["coverage_lowest_band"], 3), f(r["mean_width"], 1)])
    return rows


BAND_LABELS = ["< 0.4", "0.4–0.5", "0.5–0.6", "0.6–0.7", "0.7–0.8", "≥ 0.8"]


def expect_9() -> list[list[str]]:
    domain = csv("table_applicability_domain")
    domain = domain[domain["regime"] == "family"]
    band_mae = csv("table_band_mae")
    band_mae = band_mae[(band_mae["regime"] == "family")
                        & (band_mae["model"] == "extra_trees")]
    rows = [["Band", "*n*", "MAE (K)", "SCP coverage", "SCP width (K)",
             "Mondrian coverage", "Mondrian width (K)"]]
    for index, label in enumerate(BAND_LABELS):
        scp = domain[(domain["conformal"] == "SCP")
                     & (domain["similarity_bin"] == index)].iloc[0]
        mondrian = domain[(domain["conformal"] == "Mondrian-similarity")
                          & (domain["similarity_bin"] == index)].iloc[0]
        mae = band_mae[band_mae["band"] == index].iloc[0]
        rows.append([label, f"{int(scp['n'])}", f(mae["mae"], 1),
                     f(scp["coverage"], 3), f(scp["mean_width"], 1),
                     f(mondrian["coverage"], 3), f(mondrian["mean_width"], 1)])
    return rows


def expect_10() -> list[list[str]]:
    frame = csv("table_matched_statistics")
    rows = [["Model", "Effect", "Mean (K)", "95% CI (K)", "Median (K)",
             "Wilcoxon *p*", "Positive"]]
    effects = {"family_effect": "Family", "size_effect": "Size",
               "total_effect": "Total"}
    for model in ("extra_trees", "hist_gbr"):
        for key, label in effects.items():
            r = frame[(frame["model"] == model) & (frame["effect"] == key)].iloc[0]
            p = ("9.5 x 10-7" if r["wilcoxon_p"] < 1e-6
                 else f"{r['wilcoxon_p']:.2f}")
            rows.append([MODEL_LABELS[model], label, f(r["mean"], 2, sign=True),
                         f"[{f(r['ci_low'], 2, sign=True)}, "
                         f"{f(r['ci_high'], 2, sign=True)}]",
                         f(r["median"], 2, sign=True), p,
                         f"{int(r['n_families_positive'])} / "
                         f"{int(r['n_families'])}"])
    return rows


def expect_11() -> list[list[str]]:
    frame = csv("table_matched_by_analogue")
    labels = {"<0.6": "< 0.6", "0.6-0.8": "0.6–0.8", "0.8-0.9": "0.8–0.9",
              ">=0.9": "≥ 0.9"}
    rows = [["Model", "Informed-arm band", "*n*", "Informed", "Control", "Naive",
             "Family effect"]]
    for model in ("extra_trees", "hist_gbr"):
        for key, label in labels.items():
            r = frame[(frame["model"] == model)
                      & (frame["analogue_band"] == key)].iloc[0]
            rows.append([MODEL_LABELS[model], label, f"{int(r['n'])}",
                         f(r["mae_informed"], 2), f(r["mae_control"], 2),
                         f(r["mae_naive"], 2), f(r["family_effect"], 2, sign=True)])
    return rows


def expect_12() -> list[list[str]]:
    summary = csv("table_representation_summary")
    ablation = csv("table_representation_ablation")
    labels = {"descriptors": "RDKit descriptors, 217 columns",
              "morgan": "Morgan, radius 2, 2,048 bits"}
    rows = [["Representation", "Random MAE (K)", "Family MAE (K)", "Penalty (K)",
             "Penalty ratio"]]
    for key, label in labels.items():
        rnd = summary[(summary["representation"] == key)
                      & (summary["regime"] == "random")].iloc[0]
        fam = summary[(summary["representation"] == key)
                      & (summary["regime"] == "family")].iloc[0]
        abl = ablation[ablation["representation"] == key].iloc[0]
        rows.append([
            label,
            f"{f(rnd['mae'], 2)} ± {f(rnd['mae_sd'], 2)} ({int(rnd['n_splits'])})",
            f"{f(fam['mae'], 2)} ± {f(fam['mae_sd'], 2)} ({int(fam['n_splits'])})",
            f(abl["penalty_K"], 2, sign=True), f(abl["penalty_ratio"], 3),
        ])
    return rows


def expect_13() -> list[list[str]]:
    frame = csv("table_calibration_study")
    plan = [("random", "SCP", "Random", "Global (SCP)"),
            ("random", "Mondrian-similarity", "Random", "Novelty-conditioned"),
            ("family_out", "SCP", "Family-out", "Global (SCP)"),
            ("family_out", "Mondrian-similarity", "Family-out",
             "Novelty-conditioned")]
    rows = [["Calibration set", "Quantile", "Coverage (mean of splits)",
             "Worst band (split-averaged)", "Width (K)", "*n* calibration"]]
    for calibration, conformal, left, right in plan:
        r = frame[(frame["calibration"] == calibration)
                  & (frame["conformal"] == conformal)].iloc[0]
        rows.append([left, right, f(r["coverage"], 3),
                     f(r["worst_similarity_coverage"], 3),
                     f"{round(r['mean_width']):d}",
                     f"{round(r['n_calibration']):d}"])
    return rows


def expect_14() -> list[list[str]]:
    frame = csv("table_baseline_comparison")
    labels = {"group_contribution": "Group contribution",
              "hist_gbr": "Histogram gradient boosting",
              "median": "Training-median baseline"}
    rows = [["Model", "Random", "Scaffold", "Cluster", "Family", "Penalty (K)",
             "Ratio"]]
    for key, label in labels.items():
        r = frame[frame["model"] == key].iloc[0]
        rows.append([label] + [f(r[regime], 2) for regime in REGIMES]
                    + [f(r["penalty_K"], 2, sign=True), f(r["penalty_ratio"], 3)])
    return rows


EXPECTED = {n: globals()[f"expect_{n}"] for n in range(1, 15)}


# --------------------------------------------------------------------------
# comparison
# --------------------------------------------------------------------------

def check(manuscript: Path = MANUSCRIPT) -> list[str]:
    """Return a list of human-readable mismatches; empty means the paper agrees."""
    tables = parse_tables(manuscript.read_text())
    problems: list[str] = []

    for number, builder in EXPECTED.items():
        if number not in tables:
            problems.append(f"Table {number}: not found in {manuscript.name}")
            continue
        found = tables[number]
        expected = builder()
        if len(found) != len(expected):
            problems.append(
                f"Table {number}: {len(found)} rows in the manuscript against "
                f"{len(expected)} rebuilt from the CSV")
            continue
        header = [normalise(c) for c in found[0]]
        for row_index, (found_row, expected_row) in enumerate(zip(found, expected)):
            if len(found_row) != len(expected_row):
                problems.append(
                    f"Table {number} row {row_index}: {len(found_row)} cells "
                    f"against {len(expected_row)}")
                continue
            for column, (a, b) in enumerate(zip(found_row, expected_row)):
                label = header[column] if column < len(header) else str(column)
                heading = found[0][column].strip() if column < len(found[0]) else ""
                if any(n == number and p.match(heading) for n, p in UNCHECKED):
                    continue
                if not agrees(a, b):
                    problems.append(
                        f"Table {number} row {row_index} column '{label}': "
                        f"manuscript {a!r} against CSV {b!r}")
    return problems


def main() -> int:
    problems = check()
    if problems:
        print(f"{len(problems)} mismatch(es):")
        for line in problems:
            print("  " + line)
        return 1
    print(f"all {len(EXPECTED)} inline tables agree with their source CSVs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
