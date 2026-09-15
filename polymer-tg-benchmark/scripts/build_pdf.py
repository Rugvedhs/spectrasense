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

# Figure captions and the appendix table caption live in the manuscript source,
# under "## Figure Captions" and "## Appendix A", so that a referee reading the
# submitted source can check them.  This script reads them from there; the only
# thing it keeps is the column selection for the appendix table, which is
# typesetting rather than prose.
CAPTION_SECTIONS = ("Figure Captions", "Appendix A")

APPENDIX_TABLES = {
    # Tables 1-14 are numbered, captioned and typeset inline from the manuscript
    # source.  The appendix carries only what an inline table deliberately does
    # not: the full twenty-one-family listing that Table 5 quotes six rows of.
    "table_family_holdout.csv": (
        "Table A1",
        ["group", "n", "mae", "rmse", "r2", "deterioration", "mean_nn_similarity"],
    ),
}


def extract_captions(text: str) -> tuple[dict[int, str], dict[str, str]]:
    """Pull the figure and appendix-table captions out of the manuscript source.

    A caption is a paragraph opening ``**Figure 3.**`` or ``**Table A1.**`` in
    one of the caption sections.  Returning them here keeps the markdown the
    single point of truth: there is no second copy to drift.
    """
    figures: dict[int, str] = {}
    tables: dict[str, str] = {}
    in_section = False
    for block in re.split(r"\n\s*\n", text):
        stripped = block.strip()
        if stripped.startswith("## "):
            heading = stripped[3:].strip()
            in_section = any(heading.startswith(s) for s in CAPTION_SECTIONS)
            continue
        if not in_section:
            continue
        match = re.match(r"\*\*Figure (\d+)\.\*\*\s*(.+)", stripped, re.S)
        if match:
            figures[int(match.group(1))] = " ".join(match.group(2).split())
            continue
        match = re.match(r"\*\*(Table A\d+)\.\*\*\s*(.+)", stripped, re.S)
        if match:
            tables[match.group(1)] = " ".join(match.group(2).split())
    return figures, tables


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


def figure_flowable(number: int, figures_dir: Path, styles, width: float,
                    captions: dict[int, str]):
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
        f"<b>Figure {number}.</b> {inline(captions.get(number, ''))}",
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

    manuscript = Path(args.manuscript).read_text()
    figure_captions, appendix_captions = extract_captions(manuscript)

    # The caption sections carry text that is typeset with the figure and with
    # the appendix table, so they are not emitted a second time as body prose.
    skipping = False

    for kind, payload in parse_manuscript(manuscript):
        if kind == "h1":
            skipping = any(payload.strip().startswith(s) for s in CAPTION_SECTIONS)
        if skipping and kind != "title":
            continue
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
                                                 styles, width, figure_captions))

    # Any figure never cited still belongs in the document.
    for number in sorted(set(figure_captions) - placed):
        story.extend(figure_flowable(number, Path(args.figures_dir), styles, width,
                                     figure_captions))

    story.append(PageBreak())
    story.append(Paragraph("Appendix A. Full per-family listing", styles["h1"]))
    for filename, (label, columns) in APPENDIX_TABLES.items():
        title = f"{label}. {appendix_captions.get(label, '')}"
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
