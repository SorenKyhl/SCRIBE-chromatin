"""
Unit tests for pylib.paths module.

Tests cover:
- Project root resolution
- Data directory resolution with env vars and defaults
- Defaults directory resolution
"""

import os
import tempfile
from pathlib import Path
from unittest import mock

import pytest

from pylib.paths import get_data_dir, get_defaults_dir, get_project_root


class TestGetProjectRoot:
    """Tests for the get_project_root function."""

    def test_returns_path(self):
        """Test that get_project_root returns a Path object."""
        result = get_project_root()
        assert isinstance(result, Path)

    def test_is_scribe_root(self):
        """Test that the returned path is the SCRIBE project root."""
        result = get_project_root()
        # Should contain key project files/directories
        assert (result / "pylib").exists()
        assert (result / "README.md").exists() or (result / "setup.py").exists()

    def test_consistent_results(self):
        """Test that repeated calls return the same path."""
        result1 = get_project_root()
        result2 = get_project_root()
        assert result1 == result2


class TestGetDataDir:
    """Tests for the get_data_dir function."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_explicit_output_dir(self, temp_dir):
        """Test that explicit output_dir takes precedence."""
        result = get_data_dir(output_dir=str(temp_dir))
        assert result == temp_dir

    def test_env_variable(self, temp_dir):
        """Test that SCRIBE_DATA_DIR environment variable is respected."""
        with mock.patch.dict(os.environ, {"SCRIBE_DATA_DIR": str(temp_dir)}):
            result = get_data_dir()
            assert result == temp_dir

    def test_explicit_overrides_env(self, temp_dir):
        """Test that explicit output_dir overrides env variable."""
        other_dir = temp_dir / "other"
        with mock.patch.dict(os.environ, {"SCRIBE_DATA_DIR": str(temp_dir)}):
            result = get_data_dir(output_dir=str(other_dir))
            assert result == other_dir

    def test_default_location(self):
        """Test default location when no env var or explicit dir."""
        env = os.environ.copy()
        env.pop("SCRIBE_DATA_DIR", None)

        with mock.patch.dict(os.environ, env, clear=True):
            result = get_data_dir()
            assert result == Path.home() / ".scribe" / "data"

    def test_creates_directory(self, temp_dir):
        """Test that the directory is created if it doesn't exist."""
        new_dir = temp_dir / "new_subdir" / "data"
        assert not new_dir.exists()

        result = get_data_dir(output_dir=str(new_dir))

        assert result.exists()
        assert result.is_dir()

    def test_returns_path_object(self, temp_dir):
        """Test that get_data_dir returns a Path object."""
        result = get_data_dir(output_dir=str(temp_dir))
        assert isinstance(result, Path)


class TestGetDefaultsDir:
    """Tests for the get_defaults_dir function."""

    def test_returns_path(self):
        """Test that get_defaults_dir returns a Path object."""
        result = get_defaults_dir()
        assert isinstance(result, Path)

    def test_points_to_defaults(self):
        """Test that the path points to the defaults directory inside pylib."""
        result = get_defaults_dir()
        assert result.name == "defaults"
        assert result.parent.name == "pylib"
        assert result.parent.parent == get_project_root()

    def test_defaults_dir_exists(self):
        """Test that the defaults directory exists."""
        result = get_defaults_dir()
        assert result.exists()
        assert result.is_dir()

    def test_contains_config_files(self):
        """Test that defaults directory contains expected config files."""
        result = get_defaults_dir()
        json_files = list(result.glob("*.json"))
        assert len(json_files) > 0, "Expected at least one JSON config file"
