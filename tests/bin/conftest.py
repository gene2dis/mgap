import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin"))


@pytest.fixture
def data():
    return Path(__file__).resolve().parent / "data"
