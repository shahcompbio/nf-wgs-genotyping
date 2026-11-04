"""Shared pytest setup: put bin/ and tests/ on sys.path and build synthetic fixtures once."""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BIN = os.path.join(REPO, "bin")
for p in (BIN, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

SVTYPER_DATA = os.path.join(HERE, "data", "svtyper")


@pytest.fixture(scope="session")
def synth(tmp_path_factory):
    """Synthetic reference, BAM and caller VCFs with known truth (see make_fixtures.py)."""
    import make_fixtures

    out = tmp_path_factory.mktemp("synthetic")
    return make_fixtures.build_all(str(out))
