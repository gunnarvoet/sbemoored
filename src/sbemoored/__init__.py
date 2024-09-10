"""
Library for data processing of moored SBE instruments.

Data files for SBE56 and SBE37 sensors are quite different, therefore the
processing code has been separated into individual submodules.
"""

__all__ = ["sbe56", "sbe37"]
# __all__ = ["io", "sbe37"]

from . import sbe37, sbe56
