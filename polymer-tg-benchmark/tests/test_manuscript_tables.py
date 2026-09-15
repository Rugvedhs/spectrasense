"""The inline tables in the paper must agree with the CSVs they cite.

The fourteen numbered tables are transcriptions, so Section 2.9's claim that the
paper's numbers come from the pipeline is only true if something checks them.
This is that something.
"""

import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import check_manuscript_tables as checker  # noqa: E402


@pytest.mark.skipif(not (checker.RESULTS / "table_conformal.csv").exists(),
                    reason="results/ not populated; run scripts/run_all.sh first")
def test_inline_tables_match_their_sources():
    problems = checker.check()
    assert not problems, "\n".join(problems)


def test_every_numbered_table_is_covered():
    """A table added to the paper without a checker entry should fail here."""
    found = set(checker.parse_tables(checker.MANUSCRIPT.read_text()))
    assert found == set(checker.EXPECTED), (
        f"tables in the manuscript: {sorted(found)}; "
        f"tables with a checker: {sorted(checker.EXPECTED)}"
    )
