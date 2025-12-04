"""
Default configuration and data paths for SCRIBE.

This module provides default simulation configurations and manages data paths
for Hi-C and ChIP-seq files. Large data files are stored in a configurable
location (not tracked by git).

Configuration priority for data directory:
    1. SCRIBE_DATA_DIR environment variable
    2. ~/.scribe/data/ (default user data directory)

To set up data files, either:
    1. Set SCRIBE_DATA_DIR to point to your existing data
    2. Run: python -m pylib.download_data --setup

Expected data directory structure:
    $SCRIBE_DATA_DIR/
    ├── hic/
    │   └── HCT116_auxin/
    │       └── *.hic
    └── chipseq/
        └── HCT116_hg19/
            ├── *.bigWig
            └── metadata.tsv
"""

import json
import logging
from pathlib import Path

from pylib.paths import get_data_dir, get_defaults_dir

logger = logging.getLogger(__name__)


# =============================================================================
# Default Configurations
# =============================================================================


def _load_json(path: Path) -> dict:
    """Load a JSON file."""
    with open(path) as f:
        return json.load(f)


def load_default_config(name: str = "config") -> dict:
    """
    Load a default simulation configuration.

    Args:
        name: Name of the config file (without .json extension).
              Default is "config" which loads pylib/defaults/config.json

    Returns:
        Configuration dictionary

    Raises:
        FileNotFoundError: If the config file is not found
    """
    config_path = get_defaults_dir() / f"{name}.json"
    if config_path.exists():
        return _load_json(config_path)

    raise FileNotFoundError(f"Config '{name}' not found in {get_defaults_dir()}")


def load_default_params() -> dict:
    """
    Load default maximum entropy optimization parameters.

    Returns:
        Parameters dictionary with keys: iterations, gamma, trust_region, etc.
    """
    params_path = get_defaults_dir() / "params.json"
    if params_path.exists():
        return _load_json(params_path)

    # Fallback to sensible defaults if file doesn't exist
    logger.warning("params.json not found, using built-in defaults")
    return {
        "iterations": 10,
        "gamma": 1.0,
        "trust_region": 1000,
        "equilib_sweeps": 10000,
        "production_sweeps": 50000,
        "parallel": 1,
        "method": "n",
    }


# =============================================================================
# Data File Access
# =============================================================================


def get_hic_path(cell_type: str) -> Path | None:
    """
    Get the path to a Hi-C file for a given cell type.

    Args:
        cell_type: Cell type identifier (e.g., "HCT116_auxin", "GM12878")

    Returns:
        Path to the .hic file, or None if not found

    Example:
        >>> hic_path = get_hic_path("HCT116_auxin")
        >>> if hic_path is None:
        ...     print("Data not found. Set SCRIBE_DATA_DIR or download data.")
    """
    data_dir = get_data_dir()

    # Check for .hic files in the cell type directory
    cell_dir = data_dir / "hic" / cell_type
    if cell_dir.exists():
        hic_files = list(cell_dir.glob("*.hic"))
        if hic_files:
            return hic_files[0]

    logger.warning(f"Hi-C file for '{cell_type}' not found. " f"Expected location: {cell_dir}/")
    return None


def get_chipseq_dir(cell_type: str, genome: str = "hg19") -> Path | None:
    """
    Get the directory containing ChIP-seq files for a given cell type.

    Args:
        cell_type: Cell type identifier (e.g., "HCT116")
        genome: Genome assembly (e.g., "hg19", "hg38")

    Returns:
        Path to the ChIP-seq directory, or None if not found
    """
    chipseq_dir = get_data_dir() / "chipseq" / f"{cell_type}_{genome}"

    if chipseq_dir.exists():
        return chipseq_dir

    logger.warning(
        f"ChIP-seq directory for '{cell_type}' ({genome}) not found. "
        f"Expected location: {chipseq_dir}"
    )
    return None


def get_cached_hic(cell_type: str, chrom: str, resolution: str = "20k") -> Path | None:
    """
    Get path to a cached/preprocessed Hi-C numpy array.

    Args:
        cell_type: Cell type identifier
        chrom: Chromosome number (as string)
        resolution: Resolution string (e.g., "20k")

    Returns:
        Path to the .npy file, or None if not found
    """
    cache_file = get_data_dir() / "cache" / f"{cell_type}_chr{chrom}_{resolution}.npy"
    if cache_file.exists():
        return cache_file
    return None


# =============================================================================
# Lazy-loaded default configurations
# =============================================================================

# Cache for loaded configs
_config_cache: dict | None = None
_params_cache: dict | None = None


def _get_config() -> dict:
    """Get default config (cached)."""
    global _config_cache
    if _config_cache is None:
        try:
            _config_cache = load_default_config("config")
        except FileNotFoundError:
            logger.warning("Default config.json not found, using minimal config")
            _config_cache = {
                "nbeads": 1024,
                "nSweeps": 50000,
                "seed": 12345,
                "bonded_on": True,
                "nonbonded_on": True,
            }
    return _config_cache.copy()


def _get_params() -> dict:
    """Get default params (cached)."""
    global _params_cache
    if _params_cache is None:
        _params_cache = load_default_params()
    return _params_cache.copy()


# =============================================================================
# Public API - Module-level defaults
# =============================================================================

# Load config and params at module import time
# These are copies to prevent accidental mutation
config = _get_config()
params = _get_params()


# For explicit function access
def get_config() -> dict:
    """Get a copy of the default simulation configuration."""
    return _get_config()


def get_params() -> dict:
    """Get a copy of the default optimization parameters."""
    return _get_params()


# Default genomic region parameters
res = 100000  # 100kb resolution
chrom = 2
start = 0
end = 120_000_000
size = 1024

# Data directory
data_dir = get_data_dir()


# =============================================================================
# Lazy-loaded pipeline objects (to avoid circular imports)
# =============================================================================

_chipseq_pipeline = None
_data_pipeline = None


def _get_chipseq_pipeline():
    """Get the default ChIP-seq processing pipeline (lazy loaded)."""
    global _chipseq_pipeline
    if _chipseq_pipeline is None:
        from pylib.chipseqPipeline import ChipseqPipeline, Normalize, Sigmoid, Smooth

        _chipseq_pipeline = ChipseqPipeline([Smooth(), Normalize(), Sigmoid()])
    return _chipseq_pipeline


def _get_data_pipeline():
    """Get the default data loader (lazy loaded)."""
    global _data_pipeline
    if _data_pipeline is None:
        from pylib.dataloader import DataLoader

        _data_pipeline = DataLoader(res, chrom, start, end, size)
    return _data_pipeline


class _LazyPipeline:
    """Lazy-loading wrapper for pipeline objects."""

    def __init__(self, getter):
        self._getter = getter
        self._obj = None

    def __getattr__(self, name):
        if self._obj is None:
            self._obj = self._getter()
        return getattr(self._obj, name)

    def __repr__(self):
        if self._obj is None:
            return f"<Lazy {self._getter.__name__}>"
        return repr(self._obj)


# These are lazy-loaded to avoid circular import issues
chipseq_pipeline = _LazyPipeline(_get_chipseq_pipeline)
data_pipeline = _LazyPipeline(_get_data_pipeline)


# =============================================================================
# Backwards compatibility - hic_paths dictionary
# =============================================================================


class _HicPaths(dict):
    """
    Dictionary-like object that resolves Hi-C paths on access.

    This provides backwards compatibility with code that uses:
        hic_path = default.hic_paths["HCT116_auxin"]
    """

    # Primary supported cell type
    _cell_types = ["HCT116_auxin"]

    def __getitem__(self, key: str) -> Path | None:
        path = get_hic_path(key)
        if path is None:
            raise KeyError(
                f"Hi-C data for '{key}' not found. "
                f"Set SCRIBE_DATA_DIR environment variable or download the data."
            )
        return path

    def get(self, key: str, default=None) -> Path | None:
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key: str) -> bool:
        return get_hic_path(key) is not None

    def keys(self):
        return self._cell_types


hic_paths = _HicPaths()


# =============================================================================
# Backwards compatibility - direct path variables (deprecated)
# =============================================================================

# These are kept for backwards compatibility but will return None if data not found
# Users should migrate to get_hic_path() and get_chipseq_dir()


def _lazy_hic_path(cell_type: str) -> Path | None:
    """Lazy loader for Hi-C paths."""
    return get_hic_path(cell_type)


class _LazyPath:
    """Descriptor that lazily resolves a path."""

    def __init__(self, getter):
        self.getter = getter
        self._value = None
        self._resolved = False

    def __get__(self, obj, objtype=None):
        if not self._resolved:
            self._value = self.getter()
            self._resolved = True
        return self._value


# Legacy path variables - these will be None if data not downloaded
HCT116_hic = get_hic_path("HCT116_auxin")
HCT116_chipseq = get_chipseq_dir("HCT116", "hg19")

# Optional: GM12878 (not required for basic usage)
GM12878_hic = get_hic_path("GM12878")

# Legacy cache paths
HCT116_hic_20k = {
    "2": get_cached_hic("HCT116", "2", "20k"),
    "10": get_cached_hic("HCT116", "10", "20k"),
}
HCT116_seqs_20k = get_data_dir() / "cache" / "seqs20k.npy"
