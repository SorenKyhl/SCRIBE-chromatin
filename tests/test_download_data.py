"""
Unit tests for scribe.download_data module.

Tests cover:
- Data directory resolution
- File registry validation
- Download status checking
- URL validation
"""

import os
import tempfile
from pathlib import Path
from unittest import mock

import pytest

from scribe.download_data import (
    CHIPSEQ_FILES,
    HIC_FILES,
    check_data_status,
    download_file,
    setup_directory_structure,
)
from scribe.paths import get_data_dir

# =============================================================================
# Test fixtures
# =============================================================================


@pytest.fixture
def temp_data_dir():
    """Create a temporary directory for test data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def mock_env_data_dir(temp_data_dir):
    """Mock the SCRIBE_DATA_DIR environment variable."""
    with mock.patch.dict(os.environ, {"SCRIBE_DATA_DIR": str(temp_data_dir)}):
        yield temp_data_dir


# =============================================================================
# Tests for get_data_dir
# =============================================================================


class TestGetDataDir:
    """Tests for the get_data_dir function."""

    def test_explicit_output_dir(self, temp_data_dir):
        """Test that explicit output_dir takes precedence."""
        result = get_data_dir(output_dir=str(temp_data_dir))
        assert result == temp_data_dir
        assert result.exists()

    def test_env_variable(self, mock_env_data_dir):
        """Test that SCRIBE_DATA_DIR environment variable is respected."""
        result = get_data_dir()
        assert result == mock_env_data_dir

    def test_default_location(self, temp_data_dir):
        """Test default location when no env var or explicit dir."""
        # Clear the env var if set
        env = os.environ.copy()
        env.pop("SCRIBE_DATA_DIR", None)

        with mock.patch.dict(os.environ, env, clear=True):
            result = get_data_dir()
            assert result == Path.home() / ".scribe" / "data"

    def test_creates_directory(self, temp_data_dir):
        """Test that the directory is created if it doesn't exist."""
        new_dir = temp_data_dir / "new_subdir" / "data"
        result = get_data_dir(output_dir=str(new_dir))
        assert result.exists()
        assert result.is_dir()


# =============================================================================
# Tests for file registries
# =============================================================================


class TestFileRegistries:
    """Tests for HIC_FILES and CHIPSEQ_FILES registries."""

    def test_hic_files_structure(self):
        """Test that HIC_FILES has the expected structure."""
        assert isinstance(HIC_FILES, dict)
        assert len(HIC_FILES) > 0

        for name, info in HIC_FILES.items():
            assert isinstance(name, str)
            assert "description" in info
            assert "filename" in info
            assert "size_gb" in info
            assert isinstance(info["size_gb"], (int, float))

    def test_hic_files_have_urls(self):
        """Test that all HIC files have download URLs."""
        for name, info in HIC_FILES.items():
            assert "url" in info, f"Missing URL for {name}"
            assert info["url"].startswith("http"), f"Invalid URL for {name}"

    def test_chipseq_files_structure(self):
        """Test that CHIPSEQ_FILES has the expected structure."""
        assert isinstance(CHIPSEQ_FILES, dict)
        assert len(CHIPSEQ_FILES) > 0

        for name, info in CHIPSEQ_FILES.items():
            assert isinstance(name, str)
            assert "description" in info
            assert "files" in info
            assert isinstance(info["files"], dict)

    def test_chipseq_files_have_required_fields(self):
        """Test that all ChIP-seq files have required fields."""
        for cell_type, info in CHIPSEQ_FILES.items():
            for mark, file_info in info["files"].items():
                assert "accession" in file_info, f"Missing accession for {cell_type}/{mark}"
                assert "url" in file_info, f"Missing URL for {cell_type}/{mark}"
                assert "size_mb" in file_info, f"Missing size for {cell_type}/{mark}"
                assert file_info["url"].startswith("http"), f"Invalid URL for {cell_type}/{mark}"

    def test_chipseq_url_format(self):
        """Test that ChIP-seq URLs follow ENCODE format."""
        for _cell_type, info in CHIPSEQ_FILES.items():
            for mark, file_info in info["files"].items():
                url = file_info["url"]
                accession = file_info["accession"]
                # ENCODE URLs should contain the accession
                assert accession in url, f"URL doesn't contain accession for {mark}"
                # Should be bigWig format
                assert url.endswith(".bigWig"), f"URL doesn't end with .bigWig for {mark}"

    def test_hct116_has_expected_marks(self):
        """Test that HCT116 has the expected histone marks."""
        expected_marks = {"H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3", "H3K36me3", "H3K4me1"}

        assert "HCT116_hg19" in CHIPSEQ_FILES
        actual_marks = set(CHIPSEQ_FILES["HCT116_hg19"]["files"].keys())

        assert expected_marks == actual_marks, f"Missing marks: {expected_marks - actual_marks}"


# =============================================================================
# Tests for setup_directory_structure
# =============================================================================


class TestSetupDirectoryStructure:
    """Tests for the setup_directory_structure function."""

    def test_creates_expected_directories(self, mock_env_data_dir):
        """Test that all expected directories are created."""
        setup_directory_structure()

        expected_dirs = [
            mock_env_data_dir / "hic" / "HCT116_auxin",
            mock_env_data_dir / "chipseq" / "HCT116_hg19",
            mock_env_data_dir / "cache",
        ]

        for d in expected_dirs:
            assert d.exists(), f"Directory not created: {d}"
            assert d.is_dir()

    def test_creates_readme(self, mock_env_data_dir):
        """Test that README.md is created in data directory."""
        setup_directory_structure()

        readme_path = mock_env_data_dir / "README.md"
        assert readme_path.exists()

        content = readme_path.read_text()
        assert "SCRIBE" in content
        assert "hic" in content
        assert "chipseq" in content

    def test_idempotent(self, mock_env_data_dir):
        """Test that running setup twice doesn't cause errors."""
        setup_directory_structure()
        setup_directory_structure()  # Should not raise


# =============================================================================
# Tests for download_file
# =============================================================================


class TestDownloadFile:
    """Tests for the download_file function."""

    def test_creates_parent_directories(self, temp_data_dir):
        """Test that parent directories are created."""
        dest_path = temp_data_dir / "nested" / "path" / "file.txt"

        # Mock requests.get to avoid actual download
        mock_response = mock.Mock()
        mock_response.headers = {"content-length": "0"}
        mock_response.iter_content.return_value = []
        mock_response.raise_for_status = mock.Mock()

        with mock.patch("scribe.download_data.requests.get", return_value=mock_response):
            download_file("http://example.com/file.txt", dest_path)

        assert dest_path.parent.exists()

    def test_handles_download_error(self, temp_data_dir, capsys):
        """Test that download errors are handled gracefully."""
        dest_path = temp_data_dir / "file.txt"

        with mock.patch("scribe.download_data.requests.get", side_effect=Exception("Network error")):
            result = download_file("http://example.com/file.txt", dest_path)

        assert result is False
        captured = capsys.readouterr()
        assert "Error" in captured.out

    def test_successful_download(self, temp_data_dir):
        """Test successful file download."""
        dest_path = temp_data_dir / "file.txt"
        test_content = b"test content"

        mock_response = mock.Mock()
        mock_response.headers = {"content-length": str(len(test_content))}
        mock_response.iter_content.return_value = [test_content]
        mock_response.raise_for_status = mock.Mock()

        with mock.patch("scribe.download_data.requests.get", return_value=mock_response):
            result = download_file("http://example.com/file.txt", dest_path)

        assert result is True
        assert dest_path.exists()
        assert dest_path.read_bytes() == test_content


# =============================================================================
# Tests for check_data_status
# =============================================================================


class TestCheckDataStatus:
    """Tests for the check_data_status function."""

    def test_runs_without_error(self, mock_env_data_dir, capsys):
        """Test that check_data_status runs without errors."""
        check_data_status()

        captured = capsys.readouterr()
        assert "Hi-C Data Status" in captured.out
        assert "ChIP-seq Data Status" in captured.out

    def test_shows_missing_files(self, mock_env_data_dir, capsys):
        """Test that missing files are reported."""
        check_data_status()

        captured = capsys.readouterr()
        assert "Missing" in captured.out

    def test_shows_found_files(self, mock_env_data_dir, capsys):
        """Test that existing files are reported as found."""
        # Create a fake .hic file
        hic_dir = mock_env_data_dir / "hic" / "HCT116_auxin"
        hic_dir.mkdir(parents=True, exist_ok=True)
        (hic_dir / "test.hic").write_text("fake hic data")

        check_data_status()

        captured = capsys.readouterr()
        assert "Found" in captured.out


# =============================================================================
# Integration tests
# =============================================================================


class TestIntegration:
    """Integration tests for the download_data module."""

    def test_full_workflow(self, mock_env_data_dir, capsys):
        """Test the full workflow: setup -> status check."""
        # Setup directories
        setup_directory_structure()

        # Check status
        check_data_status()

        captured = capsys.readouterr()
        assert str(mock_env_data_dir) in captured.out

    def test_data_dir_consistency(self, mock_env_data_dir):
        """Test that get_data_dir returns consistent paths."""
        dir1 = get_data_dir()
        dir2 = get_data_dir()

        assert dir1 == dir2
        assert dir1 == mock_env_data_dir
