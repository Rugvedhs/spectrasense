"""Convert the manuscript from author-year citation keys to MDPI numeric style.

Polymers (MDPI) numbers references in order of first appearance in the text and
cites them as ``[1]``, ``[2,3]`` or ``[4-6]``.  The manuscript was drafted with
author-year keys such as ``[Tao2021]``.  This script performs the conversion once
and can then verify it, so the mapping is reproducible rather than hand-made:

    python scripts/renumber_citations.py            # convert in place
    python scripts/renumber_citations.py --check    # verify, change nothing

``--check`` asserts that every numeric citation in the manuscript resolves to an
entry in ``paper/references.md`` and that every entry is cited at least once.
It runs as part of the test suite.

The key-to-entry map below is the only hand-written part.  Each value is a string
that occurs in exactly one entry of ``paper/references.md``, so the map survives
the renumbering it performs and the script is idempotent: run a second time, it
finds no keys left and leaves both files alone.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT = ROOT / "paper" / "manuscript.md"
REFERENCES = ROOT / "paper" / "references.md"

# citation key -> a substring that identifies its entry in references.md uniquely
KEY_TO_ENTRY = {
    "Tao2021": "Benchmarking Machine Learning",
    "Casanola2024": "Casanola-Martin",
    "Kuenneth2023": "polyBERT",
    "Xu2023": "TransPolymer",
    "Jablonka2025": "PolyMetriX",
    "Xu2025": "POINT",
    "Babbar2023": "Explainability and Transferability",
    "Qiu2023": "Graph Neural Network",
    "Bicerano": "Bicerano, J.",
    "VanKrevelen2009": "Van Krevelen, D. W.",
    "Askadskii2003": "Askadskii",
    "Otsuka": "PoLyInfo: Polymer",
    "Akdogan2026": "Trust Beyond Accuracy",
    "Agrawal2026": "directed message-passing",
    "Norinder2014": "Introducing Conformal",
    "Papadopoulos2002": "Inductive Confidence",
    "Vovk2005": "Algorithmic Learning in a Random World",
    "Vovk2013": "Conditional validity of inductive",
    "Lei2018": "Distribution-Free Predictive Inference",
    "Papadopoulos2011": "Reliable prediction intervals",
    "Bostrom2020": "Mondrian conformal regressors",
    "Barber2021": "jackknife+",
    "Tibshirani2019": "Tibshirani, R. J.; Barber",
    "Bemis1996": "Bemis, G. W.",
    "Rogers2010": "Extended-Connectivity Fingerprints",
    "Sheridan2013": "Time-Split Cross-Validation",
    "Wu2018": "MoleculeNet",
    "Pedregosa2011": "Scikit-learn",
    "RDKit": "Open-source cheminformatics",
    "Teh": "Polymer Glass-Transition Temperature Prediction",
}

# Annotations naming a key no longer make sense once the keys are gone.  Written
# as patterns over flexible whitespace because the entries are wrapped.
ENTRY_CLEANUPS = [
    (r";\s+cited\s+as\s+\[Norinder2014\]\.", "."),
    (r";\s+cited\s+as\s+\[Vovk2013\]\.", "."),
    (r"cited\s+as\s+\[Tibshirani2019\]\s+and\s+discussed\s+in",
     "discussed in"),
    (r"\s+—\s+cited\s+as\s+\[Pedregosa2011\]\.", ""),
    (r"\s+—\s+cited\s+as\s+\[RDKit\]\.", ""),
    (r"dataset,\s+cited\s+as\s+\[Teh\]:", "dataset:"),
]

CITATION = re.compile(r"\[([A-Z][A-Za-z0-9]*(?:\s*,\s*[A-Z][A-Za-z0-9]*)*)\]")
NUMERIC_CITATION = re.compile(r"\[(\d+(?:\s*[,–-]\s*\d+)*)\]")
NOT_A_CITATION = {"AUTHOR"}


def split_entries(text: str) -> list[str]:
    """Return the reference entries, in file order, without their numbering."""
    entries: list[str] = []
    current: list[str] | None = None
    for line in text.split("\n"):
        if re.match(r"^\d+\.\s", line):
            if current is not None:
                entries.append("\n".join(current).rstrip())
            current = [re.sub(r"^\d+\.\s+", "", line)]
        elif line.startswith("#"):
            if current is not None:
                entries.append("\n".join(current).rstrip())
                current = None
        elif current is not None:
            if line.strip():
                current.append(line)
            else:
                entries.append("\n".join(current).rstrip())
                current = None
    if current is not None:
        entries.append("\n".join(current).rstrip())
    return entries


def format_group(numbers: list[int]) -> str:
    """MDPI grouping: [1], [1,2], [1-3] for a run of three or more."""
    numbers = sorted(numbers)
    parts: list[str] = []
    start = previous = numbers[0]
    for value in numbers[1:] + [None]:
        if value == previous + 1:
            previous = value
            continue
        if previous - start >= 2:
            parts.append(f"{start}–{previous}")
        elif previous != start:
            parts.extend([str(start), str(previous)])
        else:
            parts.append(str(start))
        if value is not None:
            start = previous = value
    return "[" + ",".join(parts) + "]"


def convert() -> int:
    manuscript = MANUSCRIPT.read_text()
    references = REFERENCES.read_text()

    entries = split_entries(references)
    if not entries:
        sys.exit("no numbered entries found in references.md")

    # Resolve every key to exactly one entry.
    key_to_index: dict[str, int] = {}
    for key, needle in KEY_TO_ENTRY.items():
        matches = [i for i, e in enumerate(entries) if needle in e]
        if len(matches) != 1:
            sys.exit(f"key {key!r} matches {len(matches)} entries via {needle!r}")
        key_to_index[key] = matches[0]

    # Number in order of first appearance.
    order: list[str] = []
    for match in CITATION.finditer(manuscript):
        for key in (k.strip() for k in match.group(1).split(",")):
            if key in NOT_A_CITATION:
                continue
            if key not in key_to_index:
                sys.exit(f"citation key {key!r} has no entry in references.md")
            if key not in order:
                order.append(key)

    if not order:
        print("manuscript already uses numeric citations; nothing to do")
        return 0

    uncited = set(KEY_TO_ENTRY) - set(order)
    if uncited:
        sys.exit(f"entries never cited: {sorted(uncited)}")

    number = {key: i + 1 for i, key in enumerate(order)}

    def replace(match: re.Match) -> str:
        keys = [k.strip() for k in match.group(1).split(",")]
        if any(k in NOT_A_CITATION for k in keys):
            return match.group(0)
        return format_group([number[k] for k in keys])

    MANUSCRIPT.write_text(CITATION.sub(replace, manuscript))

    rendered: list[str] = []
    for i, key in enumerate(order):
        label = f"{i + 1}. "
        lines = entries[key_to_index[key]].split("\n")
        indent = " " * len(label)
        rendered.append(
            label + lines[0]
            + "".join("\n" + indent + line.strip() for line in lines[1:])
        )
    body = "\n".join(rendered)
    for pattern, replacement in ENTRY_CLEANUPS:
        body = re.sub(pattern, replacement, body)
    header = (
        "# References\n\n"
        "Numbered in order of first appearance in the manuscript, as Polymers\n"
        "requires. Generated by `scripts/renumber_citations.py`.\n\n"
    )
    REFERENCES.write_text(header + body + "\n")

    for key in order:
        print(f"{number[key]:>3}  {key}")
    print(f"{len(order)} keys mapped to {len(entries)} entries")
    return 0


def check() -> int:
    manuscript = MANUSCRIPT.read_text()
    references = REFERENCES.read_text()
    entries = split_entries(references)
    problems: list[str] = []

    leftover = sorted({
        key.strip()
        for match in CITATION.finditer(manuscript)
        for key in match.group(1).split(",")
        if key.strip() not in NOT_A_CITATION
    })
    if leftover:
        problems.append(f"author-year keys still in the text: {leftover}")

    cited: set[int] = set()
    for match in NUMERIC_CITATION.finditer(manuscript):
        for part in match.group(1).split(","):
            part = part.strip()
            if "–" in part or re.fullmatch(r"\d+-\d+", part):
                low, high = re.split(r"[–-]", part)
                cited.update(range(int(low), int(high) + 1))
            else:
                cited.add(int(part))

    for value in sorted(cited):
        if not 1 <= value <= len(entries):
            problems.append(f"citation [{value}] has no reference entry")
    missing = sorted(set(range(1, len(entries) + 1)) - cited)
    if missing:
        problems.append(f"reference entries never cited: {missing}")

    if problems:
        for problem in problems:
            print(problem, file=sys.stderr)
        return 1
    print(f"all {len(entries)} references cited; every citation resolves")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="verify the numbering instead of rewriting it")
    args = parser.parse_args()
    sys.exit(check() if args.check else convert())


if __name__ == "__main__":
    main()
