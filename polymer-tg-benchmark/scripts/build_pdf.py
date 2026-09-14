"""Typeset the manuscript as a PDF, with figures and key tables placed inline.

The markdown source stays the single point of truth: this renders it rather than
holding a second copy of the text.  Figures are inserted after the paragraph that
first cites them, and the result tables are read from ``results/`` at build time,
so a rebuilt PDF always carries current numbers.

Fonts are DejaVu rather than the ReportLab built-ins.  The manuscript contains
minus signs, multiplication signs, superscripts and Greek letters, none of which
exist in the standard Type-1 fonts; they would render as black boxes.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

FONT_DIR = Path(matplotlib.__file__).parent / "mpl-data" / "fonts" / "ttf"

FIGURE_CAPTIONS = {
    1: "Reported glass transition temperature by repeat-unit family, ordered by "
       "median; boxes span the interquartile range, whiskers the 5th to 95th "
       "percentiles, and family sizes are listed at the right. The ordering is "
       "the one polymer chemistry predicts, which is evidence that the derived "
       "taxonomy is meaningful rather than merely self-consistent.",
    2: "Test MAE by model and split regime; error bars are the standard deviation "
       "across repeated splits. The spread between regimes dwarfs the spread "
       "between models.",
    3: "Matched-pair design. (a) The three training arms predicting an identical "
       "set of held-out structures; informed and control are almost coincident, "
       "while naive is displaced to the right in every family. (b) The family "
       "effect, naive minus control, with training-set size held constant.",
    4: "Coverage of nominally 90% conformal intervals. (a) Marginal coverage. "
       "(b) Coverage in the least-similar structural band. Line style encodes "
       "method alongside colour, so the exact coincidence of split conformal and "
       "the family-conditioned quantile under family holdout remains visible.",
    5: "Conditional coverage against nearest-training Tanimoto similarity, one "
       "panel per split regime. Split conformal degrades monotonically with "
       "structural novelty; the novelty-conditioned quantile stays near nominal.",
    6: "Worst similarity-band coverage against mean interval width. Marker colour "
       "is the conformal method, shape the split regime.",
    7: "Mean absolute error against nearest-training similarity, for random and "
       "family-holdout splits. Group sizes are annotated.",
}

TABLES = {
    "table_point_by_regime.csv": (
        "Point accuracy by split regime and model",
        ["regime", "model", "n_splits", "mae_mean", "mae_sd", "rmse_mean", "r2_mean"],
    ),
    "table_conformal.csv": (
        "Conformal coverage, width and subgroup validity",
        ["regime", "conformal", "coverage", "worst_similarity_coverage",
         "mean_width", "fallback"],
    ),
    "table_matched_statistics.csv": (
        "Matched-pair effects, with bootstrap confidence intervals",
        ["model", "effect", "n_families", "mean", "ci_low", "ci_high",
         "wilcoxon_p", "n_families_positive"],
    ),
    "table_family_holdout.csv": (
        "Family-holdout error and deterioration (histogram gradient boosting)",
        ["group", "n", "mae", "r2", "deterioration", "mean_nn_similarity"],
    ),
    "table_pooled_r2.csv": (
        "Pooled R2 against the mean of per-split R2, by regime and model",
        ["regime", "model", "n_predictions", "r2_pooled", "r2_mean_of_splits"],
    ),
    "table_baseline_comparison.csv": (
        "Group-contribution baseline against the ensemble, by regime",
        ["model", "random", "scaffold", "cluster", "family", "penalty_K",
         "penalty_ratio"],
    ),
}


def register_fonts() -> None:
    faces = {
        "Body": "DejaVuSerif.ttf",
        "Body-Bold": "DejaVuSerif-Bold.ttf",
        "Body-Italic": "DejaVuSerif-Italic.ttf",
        "Body-BoldItalic": "DejaVuSerif-BoldItalic.ttf",
        "Head": "DejaVuSans-Bold.ttf",
        "Mono": "DejaVuSansMono.ttf",
    }
    for name, filename in faces.items():
        pdfmetrics.registerFont(TTFont(name, str(FONT_DIR / filename)))
    pdfmetrics.registerFontFamily(
        "Body", normal="Body", bold="Body-Bold",
        italic="Body-Italic", boldItalic="Body-BoldItalic",
    )


def inline(text: str) -> str:
    """Convert the markdown subset used in the manuscript to ReportLab markup."""
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"~([A-Za-z0-9]+)~", r"<sub>\1</sub>", text)      # T~g~

    # Code spans are lifted out before emphasis is processed: their contents are
    # literal, and a span such as `*` would otherwise be eaten by the italic rule
    # and leave unbalanced markup behind.
    spans: list[str] = []

    def _stash(match: re.Match) -> str:
        spans.append(match.group(1))
        return f"\x00{len(spans) - 1}\x00"

    text = re.sub(r"`([^`]+)`", _stash, text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(
        r"\x00(\d+)\x00",
        lambda m: f'<font face="Mono" size="8">{spans[int(m.group(1))]}</font>',
        text,
    )
    # Superscripts that would otherwise rely on a single Unicode codepoint.
    text = text.replace("²", "<super>2</super>")
    text = re.sub(r"10⁻([0-9⁰¹²³⁴⁵⁶⁷⁸⁹]+)", lambda m: "10<super>-"
                  + m.group(1).translate(str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789"))
                  + "</super>", text)
    return text


def build_styles():
    styles = getSampleStyleSheet()
    base = dict(fontName="Body", fontSize=9.5, leading=13.5, alignment=TA_JUSTIFY)
    return {
        "title": ParagraphStyle("title", **{**base, "fontName": "Head", "fontSize": 16,
                                            "leading": 20, "alignment": TA_CENTER,
                                            "spaceAfter": 10}),
        "h1": ParagraphStyle("h1", **{**base, "fontName": "Head", "fontSize": 12.5,
                                      "leading": 16, "spaceBefore": 14, "spaceAfter": 6,
                                      "alignment": 0}),
        "h2": ParagraphStyle("h2", **{**base, "fontName": "Head", "fontSize": 10.5,
                                      "leading": 14, "spaceBefore": 10, "spaceAfter": 4,
                                      "alignment": 0}),
        "body": ParagraphStyle("body", **{**base, "spaceAfter": 6}),
        "abstract": ParagraphStyle("abstract", **{**base, "fontSize": 9,
                                                  "leading": 12.5, "leftIndent": 10,
                                                  "rightIndent": 10, "spaceAfter": 8}),
        "bullet": ParagraphStyle("bullet", **{**base, "leftIndent": 14,
                                              "bulletIndent": 4, "spaceAfter": 3}),
        "caption": ParagraphStyle("caption", **{**base, "fontSize": 8, "leading": 11,
                                                "spaceBefore": 3, "spaceAfter": 12,
                                                "textColor": colors.HexColor("#333333")}),
        "cell": ParagraphStyle("cell", fontName="Body", fontSize=6.6, leading=8.4),
        "cellhead": ParagraphStyle("cellhead", fontName="Body-Bold", fontSize=6.6,
                                   leading=8.4),
    }


def figure_flowable(number: int, figures_dir: Path, styles, width: float):
    path = figures_dir / f"fig{number}_*.png"
    matches = sorted(figures_dir.glob(f"fig{number}_*.png"))
    if not matches:
        return []
    from PIL import Image as PILImage

    with PILImage.open(matches[0]) as im:
        w, h = im.size
        dpi = im.info.get("dpi", (600, 600))[0] or 600

    # Figures are drawn at their intended printed width. Stretching them to the
    # text measure would enlarge the type inside them past the body text, which
    # is the usual reason a figure looks pasted in rather than typeset.
    natural = w / dpi * 72.0
    drawn = min(natural, width)
    image = Image(str(matches[0]), width=drawn, height=drawn * h / w)
    image.hAlign = "CENTER"
    caption = Paragraph(
        f"<b>Figure {number}.</b> {inline(FIGURE_CAPTIONS.get(number, ''))}",
        styles["caption"],
    )
    return [Spacer(1, 6), KeepTogether([image, caption])]


def table_flowable(path: Path, title: str, columns, styles, width: float):
    if not path.exists():
        return []
    frame = pd.read_csv(path)
    if path.name == "table_family_holdout.csv":
        frame = frame[frame["model"] == "hist_gbr"]
    columns = [c for c in columns if c in frame.columns]
    frame = frame[columns].copy()
    for column in frame.columns:
        if pd.api.types.is_float_dtype(frame[column]):
            frame[column] = frame[column].map(
                lambda v: "" if pd.isna(v) else (f"{v:.3g}" if abs(v) < 0.01 else f"{v:.3f}")
            )
    header = [Paragraph(c.replace("_", " "), styles["cellhead"]) for c in frame.columns]
    rows = [header] + [
        [Paragraph(str(v), styles["cell"]) for v in row]
        for row in frame.astype(str).itertuples(index=False)
    ]
    table = Table(rows, repeatRows=1, hAlign="LEFT", colWidths=width / len(frame.columns))
    table.setStyle(TableStyle([
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#bbbbbb")),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef2f7")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    return [
        Spacer(1, 8),
        Paragraph(f"<b>{inline(title)}</b>", styles["caption"]),
        table,
        Spacer(1, 10),
    ]


def parse_manuscript(text: str):
    """Yield ('kind', payload) blocks from the markdown source."""
    lines = text.split("\n")
    if lines and lines[0].strip() == "---":
        end = lines.index("---", 1)
        front = "\n".join(lines[1:end])
        title = re.search(r'title:\s*"(.*)"', front)
        yield "title", title.group(1) if title else "Manuscript"
        lines = lines[end + 1:]

    buffer: list[str] = []
    table: list[list[str]] = []

    def flush_table():
        """Drop the |---|---| separator row and hand back the cells."""
        rows = [r for r in table if not re.match(r"^[\s|:-]+$", "|".join(r))]
        return rows

    for raw in lines:
        line = raw.rstrip()

        if line.strip().startswith("|") and line.strip().endswith("|"):
            if buffer:
                yield "para", " ".join(buffer); buffer = []
            table.append([c.strip() for c in line.strip().strip("|").split("|")])
            continue
        if table:
            rows = flush_table()
            if rows:
                yield "table", rows
            table = []

        if not line.strip():
            if buffer:
                yield "para", " ".join(buffer)
                buffer = []
            continue
        if line.startswith("### "):
            if buffer:
                yield "para", " ".join(buffer); buffer = []
            yield "h2", line[4:]
        elif line.startswith("## "):
            if buffer:
                yield "para", " ".join(buffer); buffer = []
            yield "h1", line[3:]
        elif re.match(r"^\s*[*-]\s+", line) or re.match(r"^\s*\d+\.\s+", line):
            if buffer:
                yield "para", " ".join(buffer); buffer = []
            yield "bullet", re.sub(r"^\s*(?:[*-]|\d+\.)\s+", "", line)
        else:
            buffer.append(line.strip())
    if table:
        rows = flush_table()
        if rows:
            yield "table", rows
    if buffer:
        yield "para", " ".join(buffer)


def inline_table(rows, styles, width: float):
    """Render an inline markdown table with the same look as the appendix ones."""
    n_columns = max(len(r) for r in rows)
    rendered = []
    for index, row in enumerate(rows):
        row = row + [""] * (n_columns - len(row))
        style = styles["cellhead"] if index == 0 else styles["cell"]
        rendered.append([Paragraph(inline(c), style) for c in row])

    table = Table(rendered, repeatRows=1, hAlign="LEFT",
                  colWidths=width / n_columns)
    table.setStyle(TableStyle([
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#bbbbbb")),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef2f7")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return [Spacer(1, 6), table, Spacer(1, 8)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", default="paper/manuscript.md")
    parser.add_argument("--references", default="paper/references.md")
    parser.add_argument("--figures-dir", default="paper/figures")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--out", default="paper/polymer_tg_benchmark.pdf")
    args = parser.parse_args()

    register_fonts()
    styles = build_styles()

    doc = SimpleDocTemplate(
        args.out, pagesize=A4,
        leftMargin=20 * mm, rightMargin=20 * mm,
        topMargin=18 * mm, bottomMargin=18 * mm,
        title="Backbone chemistry and the reliability of machine-learned polymer Tg",
    )
    width = doc.width
    story: list = []
    placed: set[int] = set()
    in_abstract = False

    for kind, payload in parse_manuscript(Path(args.manuscript).read_text()):
        if kind == "title":
            story.append(Paragraph(inline(payload), styles["title"]))
            story.append(Spacer(1, 4))
        elif kind == "h1":
            in_abstract = payload.strip().lower() == "abstract"
            story.append(Paragraph(inline(payload), styles["h1"]))
        elif kind == "h2":
            in_abstract = False
            story.append(Paragraph(inline(payload), styles["h2"]))
        elif kind == "bullet":
            story.append(Paragraph(inline(payload), styles["bullet"], bulletText="•"))
        elif kind == "table":
            story.extend(inline_table(payload, styles, width))
        else:
            style = styles["abstract"] if in_abstract else styles["body"]
            story.append(Paragraph(inline(payload), style))
            for number in sorted(set(int(n) for n in re.findall(r"Figure (\d)", payload))):
                if number not in placed:
                    placed.add(number)
                    story.extend(figure_flowable(number, Path(args.figures_dir),
                                                 styles, width))

    # Any figure never cited still belongs in the document.
    for number in sorted(set(FIGURE_CAPTIONS) - placed):
        story.extend(figure_flowable(number, Path(args.figures_dir), styles, width))

    story.append(PageBreak())
    story.append(Paragraph("Appendix A. Result tables", styles["h1"]))
    for filename, (title, columns) in TABLES.items():
        story.extend(table_flowable(Path(args.results_dir) / filename, title,
                                    columns, styles, width))

    references = Path(args.references)
    if references.exists():
        story.append(PageBreak())
        for kind, payload in parse_manuscript(references.read_text()):
            if kind == "h1":
                story.append(Paragraph(inline(payload), styles["h1"]))
            elif kind == "h2":
                story.append(Paragraph(inline(payload), styles["h2"]))
            elif kind == "bullet":
                story.append(Paragraph(inline(payload), styles["bullet"],
                                       bulletText="•"))
            elif kind == "table":
                story.extend(inline_table(payload, styles, width))
            else:
                story.append(Paragraph(inline(payload), styles["body"]))

    doc.build(story)
    size = Path(args.out).stat().st_size / 1024
    print(f"wrote {args.out} ({size:.0f} KB)")


if __name__ == "__main__":
    main()
