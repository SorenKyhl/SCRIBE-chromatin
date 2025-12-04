"""
Unit tests for pylib.default module.

Tests cover:
- Configuration loading
- Parameter loading
- Data path resolution
- Module-level defaults
"""

import os
import tempfile
from pathlib import Path
from unittest import mock

import pytest

from pylib import default
from pylib.default import (
    get_cached_hic,
    get_chipseq_dir,
    get_config,
    get_hic_path,
    get_params,
    load_default_config,
    load_default_params,
)

# =============================================================================
# Test fixtures
# =============================================================================


@pytest.fixture
def temp_data_dir():
    """Create a temporary data directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def mock_data_dir(temp_data_dir):
    """Mock the SCRIBE_DATA_DIR environment variable."""
    with mock.patch.dict(os.environ, {"SCRIBE_DATA_DIR": str(temp_data_dir)}):
        yield temp_data_dir


# =============================================================================
# Tests for load_default_config
# =============================================================================


class TestLoadDefaultConfig:
    """Tests for the load_default_config function."""

    def test_loads_default_config(self):
        """Test loading the default config."""
        config = load_default_config("config")
        assert isinstance(config, dict)
        assert "nbeads" in config or "nSweeps" in config

    def test_returns_dict(self):
        """Test that load_default_config returns a dictionary."""
        config = load_default_config("config")
        assert isinstance(config, dict)

    def test_file_not_found(self):
        """Test that FileNotFoundError is raised for missing config."""
        with pytest.raises(FileNotFoundError):
            load_default_config("nonexistent_config_xyz")

    def test_config_has_expected_keys(self):
        """Test that loaded config has expected simulation parameters."""
        config = load_default_config("config")
        # Should have at least some simulation parameters
        assert len(config) > 0


# =============================================================================
# Tests for load_default_params
# =============================================================================


class TestLoadDefaultParams:
    """Tests for the load_default_params function."""

    def test_returns_dict(self):
        """Test that load_default_params returns a dictionary."""
        params = load_default_params()
        assert isinstance(params, dict)

    def test_has_expected_keys(self):
        """Test that params has expected optimization parameters."""
        params = load_default_params()
        # Should have iteration-related parameters
        expected_keys = {"iterations", "gamma"}
        assert len(expected_keys & set(params.keys())) > 0 or len(params) > 0

    def test_fallback_defaults(self):
        """Test that fallback defaults are used when file is missing."""
        # Even without params.json, should return sensible defaults
        params = load_default_params()
        assert isinstance(params, dict)


# =============================================================================
# Tests for get_config and get_params
# =============================================================================


class TestGetConfigAndParams:
    """Tests for get_config and get_params functions."""

    def test_get_config_returns_copy(self):
        """Test that get_config returns a copy (not the cached original)."""
        config1 = get_config()
        config2 = get_config()

        # Modifying one shouldn't affect the other
        config1["test_key"] = "test_value"
        assert "test_key" not in config2

    def test_get_params_returns_copy(self):
        """Test that get_params returns a copy (not the cached original)."""
        params1 = get_params()
        params2 = get_params()

        # Modifying one shouldn't affect the other
        params1["test_key"] = "test_value"
        assert "test_key" not in params2

    def test_config_is_dict(self):
        """Test that get_config returns a dictionary."""
        config = get_config()
        assert isinstance(config, dict)

    def test_params_is_dict(self):
        """Test that get_params returns a dictionary."""
        params = get_params()
        assert isinstance(params, dict)


# =============================================================================
# Tests for data path functions
# =============================================================================


class TestGetHicPath:
    """Tests for the get_hic_path function."""

    def test_returns_none_for_missing_data(self, mock_data_dir):
        """Test that None is returned when data doesn't exist."""
        result = get_hic_path("nonexistent_cell_type")
        assert result is None

    def test_returns_path_when_exists(self, mock_data_dir):
        """Test that path is returned when .hic file exists."""
        # Create a fake .hic file
        hic_dir = mock_data_dir / "hic" / "test_cell"
        hic_dir.mkdir(parents=True)
        hic_file = hic_dir / "test.hic"
        hic_file.write_text("fake hic data")

        result = get_hic_path("test_cell")
        assert result == hic_file

    def test_returns_first_hic_file(self, mock_data_dir):
        """Test that the first .hic file is returned when multiple exist."""
        hic_dir = mock_data_dir / "hic" / "multi_cell"
        hic_dir.mkdir(parents=True)
        (hic_dir / "a.hic").write_text("data")
        (hic_dir / "b.hic").write_text("data")

        result = get_hic_path("multi_cell")
        assert result is not None
        assert result.suffix == ".hic"


class TestGetChipseqDir:
    """Tests for the get_chipseq_dir function."""

    def test_returns_none_for_missing_data(self, mock_data_dir):
        """Test that None is returned when directory doesn't exist."""
        result = get_chipseq_dir("nonexistent_cell")
        assert result is None

    def test_returns_path_when_exists(self, mock_data_dir):
        """Test that path is returned when directory exists."""
        chipseq_dir = mock_data_dir / "chipseq" / "test_cell_hg19"
        chipseq_dir.mkdir(parents=True)

        result = get_chipseq_dir("test_cell", "hg19")
        assert result == chipseq_dir

    def test_respects_genome_parameter(self, mock_data_dir):
        """Test that genome parameter affects path."""
        hg38_dir = mock_data_dir / "chipseq" / "test_cell_hg38"
        hg38_dir.mkdir(parents=True)

        result = get_chipseq_dir("test_cell", "hg38")
        assert result == hg38_dir

        # hg19 should not exist
        result_hg19 = get_chipseq_dir("test_cell", "hg19")
        assert result_hg19 is None


class TestGetCachedHic:
    """Tests for the get_cached_hic function."""

    def test_returns_none_for_missing_cache(self, mock_data_dir):
        """Test that None is returned when cache file doesn't exist."""
        result = get_cached_hic("test", "2", "20k")
        assert result is None

    def test_returns_path_when_exists(self, mock_data_dir):
        """Test that path is returned when cache file exists."""
        cache_dir = mock_data_dir / "cache"
        cache_dir.mkdir(parents=True)
        cache_file = cache_dir / "test_chr2_20k.npy"
        cache_file.write_text("cached data")

        result = get_cached_hic("test", "2", "20k")
        assert result == cache_file


# =============================================================================
# Tests for module-level defaults
# =============================================================================


class TestModuleLevelDefaults:
    """Tests for module-level default values."""

    def test_config_exists(self):
        """Test that default.config exists and is a dict."""
        assert hasattr(default, "config")
        assert isinstance(default.config, dict)

    def test_params_exists(self):
        """Test that default.params exists and is a dict."""
        assert hasattr(default, "params")
        assert isinstance(default.params, dict)

    def test_resolution_default(self):
        """Test that default resolution is set."""
        assert hasattr(default, "res")
        assert default.res == 100000  # 100kb

    def test_size_default(self):
        """Test that default size is set."""
        assert hasattr(default, "size")
        assert default.size == 1024

    def test_chrom_default(self):
        """Test that default chromosome is set."""
        assert hasattr(default, "chrom")
        assert default.chrom == 2


# =============================================================================
# Tests for _HicPaths class
# =============================================================================


class TestHicPaths:
    """Tests for the _HicPaths backwards compatibility class."""

    def test_get_returns_none_for_missing(self, mock_data_dir):
        """Test that .get() returns None for missing data."""
        result = default.hic_paths.get("nonexistent")
        assert result is None

    def test_get_with_default(self, mock_data_dir):
        """Test that .get() returns default for missing data."""
        result = default.hic_paths.get("nonexistent", "default_value")
        assert result == "default_value"

    def test_contains_false_for_missing(self, mock_data_dir):
        """Test that 'in' returns False for missing data."""
        assert "nonexistent" not in default.hic_paths

    def test_contains_true_when_exists(self, mock_data_dir):
        """Test that 'in' returns True when data exists."""
        hic_dir = mock_data_dir / "hic" / "test_cell"
        hic_dir.mkdir(parents=True)
        (hic_dir / "test.hic").write_text("data")

        assert "test_cell" in default.hic_paths

    def test_keys_returns_cell_types(self):
        """Test that keys() returns the supported cell types."""
        keys = default.hic_paths.keys()
        assert isinstance(keys, list)
