"""Integration tests for simulation pooling consistency.

These tests verify that the C++ simulation's contact map pooling matches
the Python pooling functions in scribe.hic.

These tests require the scribe_engine C++ extension to be built (run `make all`).
They are automatically skipped if the extension is not available.
"""

import shutil
import tempfile

import numpy as np
import pytest

# Check if simulation engine is available
try:
    from scribe import epilib, hic
    from scribe.ideal_chain import ideal_chain_simulation

    HAS_SIMULATION = True
except ImportError:
    HAS_SIMULATION = False

pytestmark = pytest.mark.skipif(
    not HAS_SIMULATION, reason="scribe_engine C++ extension not built (run `make all`)"
)


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for simulation output."""
    tmpdir = tempfile.mkdtemp(prefix="scribe_test_")
    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


class TestSimulationPoolingConsistency:
    """Test that simulation contact pooling matches Python pooling functions.

    This is the key integration test from test_suite/pooling/run.py.
    It verifies that:
    1. Running simulation with contact_resolution=2 and conservative pooling
       produces the same result as running at full resolution then applying hic.pool()
    2. Same for non-conservative pooling with hic.pool_sum()
    """

    @pytest.mark.slow
    def test_conservative_pooling_matches(self, temp_output_dir):
        """C++ conservative pooling should match Python hic.pool()."""
        import os

        os.chdir(temp_output_dir)

        nbeads = 256  # Smaller for faster test
        nsweeps = 5000

        # Run simulation without pooling (full resolution)
        sim_full = ideal_chain_simulation(nbeads)
        sim_full.config["nSweeps"] = nsweeps
        sim_full.config["conservative_contact_pooling"] = True
        sim_full.run("no_pooling")

        # Run simulation with C++ pooling (resolution=2)
        sim_pooled = ideal_chain_simulation(nbeads)
        sim_pooled.config["nSweeps"] = nsweeps
        sim_pooled.config["contact_resolution"] = 2
        sim_pooled.config["conservative_contact_pooling"] = True
        sim_pooled.run("pooling")

        # Load results
        result_full = epilib.Sim("no_pooling/")
        result_pooled = epilib.Sim("pooling/")

        # Apply Python pooling to full resolution result
        python_pooled = hic.pool(result_full.hic, 2)
        python_pooled_diag = epilib.get_diagonal(python_pooled)

        # Compare diagonals - they should match
        np.testing.assert_allclose(
            result_pooled.d,
            python_pooled_diag,
            rtol=1e-5,
            err_msg="C++ conservative pooling doesn't match Python hic.pool()",
        )

    @pytest.mark.slow
    def test_nonconservative_pooling_matches(self, temp_output_dir):
        """C++ non-conservative pooling should match Python hic.pool_sum()."""
        import os

        os.chdir(temp_output_dir)

        nbeads = 256  # Smaller for faster test
        nsweeps = 5000

        # Run simulation without pooling (full resolution)
        sim_full = ideal_chain_simulation(nbeads)
        sim_full.config["nSweeps"] = nsweeps
        sim_full.config["conservative_contact_pooling"] = False
        sim_full.run("no_pooling")

        # Run simulation with C++ non-conservative pooling
        sim_pooled = ideal_chain_simulation(nbeads)
        sim_pooled.config["nSweeps"] = nsweeps
        sim_pooled.config["contact_resolution"] = 2
        sim_pooled.config["conservative_contact_pooling"] = False
        sim_pooled.run("pooling_nonconservative")

        # Load results
        result_full = epilib.Sim("no_pooling/")
        result_pooled = epilib.Sim("pooling_nonconservative/")

        # Apply Python pooling to full resolution result
        python_pooled = hic.pool_sum(result_full.hic, 2)
        python_pooled_diag = epilib.get_diagonal(python_pooled)

        # Compare diagonals - they should match
        np.testing.assert_allclose(
            result_pooled.d,
            python_pooled_diag,
            rtol=1e-5,
            err_msg="C++ non-conservative pooling doesn't match Python hic.pool_sum()",
        )


class TestIdealChainSimulation:
    """Basic tests for ideal chain simulation."""

    @pytest.mark.slow
    def test_ideal_chain_runs(self, temp_output_dir):
        """Ideal chain simulation should run without errors."""
        import os

        os.chdir(temp_output_dir)

        sim = ideal_chain_simulation(128)
        sim.config["nSweeps"] = 1000
        sim.run("test_run")

        # Check output files exist
        assert (temp_output_dir / "test_run" / "hic.npy").exists() or os.path.exists(
            os.path.join(temp_output_dir, "test_run", "hic.npy")
        )

    @pytest.mark.slow
    def test_ideal_chain_produces_contact_map(self, temp_output_dir):
        """Ideal chain simulation should produce a valid contact map."""
        import os

        os.chdir(temp_output_dir)

        nbeads = 128
        sim = ideal_chain_simulation(nbeads)
        sim.config["nSweeps"] = 1000
        sim.run("test_run")

        result = epilib.Sim("test_run/")

        # Contact map should be symmetric
        np.testing.assert_array_almost_equal(result.hic, result.hic.T)

        # Contact map should be the right size
        assert result.hic.shape == (nbeads, nbeads)

        # Diagonal should have highest values (self-contacts)
        diag = np.diag(result.hic)
        off_diag = result.hic[~np.eye(nbeads, dtype=bool)]
        assert np.mean(diag) > np.mean(off_diag)
