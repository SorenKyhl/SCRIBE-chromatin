"""
Path resolution for SCRIBE.

This module provides path utilities for locating project files and data directories.
It has no dependencies on other pylib modules to avoid circular imports.

Configuration priority for data directory:
    1. Explicit output_dir parameter (if provided)
    2. SCRIBE_DATA_DIR environment variable
    3. ~/.scribe/data/ (default user data directory)
"""

import os
from pathlib import Path


def get_project_root() -> Path:
    """Get the root directory of the SCRIBE project."""
    return Path(__file__).parent.parent


def get_data_dir(output_dir: str | None = None) -> Path:
    """
    Get the data directory for large files (Hi-C, ChIP-seq).

    Args:
        output_dir: Optional explicit path to use as data directory

    Returns:
        Path to the data directory (creates it if it doesn't exist)

    Environment Variables:
        SCRIBE_DATA_DIR: Override the default data directory location
    """
    if output_dir:
        data_dir = Path(output_dir)
    elif "SCRIBE_DATA_DIR" in os.environ:
        data_dir = Path(os.environ["SCRIBE_DATA_DIR"])
    else:
        data_dir = Path.home() / ".scribe" / "data"

    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


def get_defaults_dir() -> Path:
    """Get the directory containing default configuration files (inside pylib)."""
    return Path(__file__).parent / "defaults"
