import pathlib
import pytest


# Determine paths and make them available in tests using fixtures.
# Any function decorated with @pytest.fixture will be availabe as input
# variable in other test files.
@pytest.fixture
def rootdir():
    return pathlib.Path(__file__).parent.resolve()
