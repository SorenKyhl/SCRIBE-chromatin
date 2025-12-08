"""
Backward compatibility module.

This module has been renamed to `scribe.analysis`.
All imports from `scribe.epilib` will continue to work but are deprecated.

Use `from scribe.analysis import SimulationResult` instead of `from scribe.epilib import Sim`.
"""

import warnings

# Re-export everything from analysis for backward compatibility
from scribe.analysis import *  # noqa: F401, F403
from scribe.analysis import SimulationResult, Sim  # noqa: F401

warnings.warn(
    "scribe.epilib is deprecated. Use scribe.analysis instead. "
    "epilib.Sim is now analysis.SimulationResult.",
    DeprecationWarning,
    stacklevel=2,
)
