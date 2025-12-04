"""Tests for scribe.hic module - Hi-C manipulation functions."""

import numpy as np
import pytest

from scribe import hic


class TestPool:
    """Tests for hic.pool() - conservative sum pooling."""

    def test_pool_reduces_size(self):
        """Pool should reduce matrix dimensions by the factor."""
        mat = np.ones((8, 8))
        result = hic.pool(mat, factor=2, normalize=False)
        assert result.shape == (4, 4)

    def test_pool_requires_even_division(self):
        """Pool should raise error if factor doesn't evenly divide size."""
        mat = np.ones((7, 7))
        with pytest.raises(AssertionError):
            hic.pool(mat, factor=2)

    def test_pool_symmetric_output(self):
        """Pool should produce symmetric output from symmetric input."""
        mat = np.random.rand(8, 8)
        mat = (mat + mat.T) / 2  # Make symmetric
        result = hic.pool(mat, factor=2, normalize=False)
        np.testing.assert_array_almost_equal(result, result.T)

    def test_pool_conserves_contacts(self):
        """Conservative pool should preserve total contact sum (approximately)."""
        # Create a symmetric matrix with known sum
        mat = np.random.rand(8, 8)
        mat = (mat + mat.T) / 2
        np.fill_diagonal(mat, np.random.rand(8))

        original_sum = np.sum(np.triu(mat))
        result = hic.pool(mat, factor=2, normalize=False)
        pooled_sum = np.sum(np.triu(result))

        # Should conserve upper triangle sum
        np.testing.assert_almost_equal(original_sum, pooled_sum, decimal=10)

    def test_pool_with_normalization(self):
        """Pool with normalize=True should produce normalized output."""
        mat = np.random.rand(8, 8) * 10
        mat = (mat + mat.T) / 2
        result = hic.pool(mat, factor=2, normalize=True)
        # With "ones" normalization, diagonal should all be 1
        np.testing.assert_array_almost_equal(np.diag(result), np.ones(4))


class TestPoolSum:
    """Tests for hic.pool_sum() - non-conservative sum pooling."""

    def test_pool_sum_reduces_size(self):
        """Pool sum should reduce matrix dimensions by the factor."""
        mat = np.ones((8, 8))
        result = hic.pool_sum(mat, factor=2, normalize=False)
        assert result.shape == (4, 4)

    def test_pool_sum_values(self):
        """Pool sum should sum all values in each block."""
        mat = np.ones((4, 4))
        result = hic.pool_sum(mat, factor=2, normalize=False)
        # Each 2x2 block of ones sums to 4
        expected = np.full((2, 2), 4.0)
        np.testing.assert_array_equal(result, expected)

    def test_pool_sum_with_normalization(self):
        """Pool sum with normalize=True should produce normalized output."""
        mat = np.random.rand(8, 8) * 10
        mat = (mat + mat.T) / 2
        result = hic.pool_sum(mat, factor=2, normalize=True)
        # With "ones" normalization, diagonal should all be 1
        np.testing.assert_array_almost_equal(np.diag(result), np.ones(4))


class TestPoolDiagonal:
    """Tests for hic.pool_diagonal() - diagonal pooling."""

    def test_pool_diagonal_reduces_size_by_half(self):
        """Pool diagonal should reduce matrix to half size."""
        mat = np.ones((8, 8))
        result = hic.pool_diagonal(mat, normalize=False)
        assert result.shape == (4, 4)

    def test_pool_diagonal_samples_diagonal_elements(self):
        """Pool diagonal should sample from diagonal positions."""
        mat = np.zeros((4, 4))
        # Set specific elements that should be sampled
        mat[0, 0] = 1
        mat[1, 1] = 2
        mat[2, 2] = 3
        mat[3, 3] = 4

        result = hic.pool_diagonal(mat, normalize=False)
        # Result[0,0] = mat[0,0] + mat[1,1] = 3
        # Result[1,1] = mat[2,2] + mat[3,3] = 7
        assert result[0, 0] == 3
        assert result[1, 1] == 7


class TestPoolSeqs:
    """Tests for hic.pool_seqs() - sequence pooling."""

    def test_pool_seqs_1d(self):
        """Pool seqs should work on 1D arrays."""
        seq = np.array([1.0, 2.0, 3.0, 4.0])
        result = hic.pool_seqs(seq, factor=2)
        expected = np.array([1.5, 3.5])  # Mean of pairs
        np.testing.assert_array_equal(result, expected)

    def test_pool_seqs_2d(self):
        """Pool seqs should work on 2D arrays (multiple sequences)."""
        seqs = np.array([[1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0]])
        result = hic.pool_seqs(seqs, factor=2)
        expected = np.array([[1.5, 3.5], [3.0, 7.0]])
        np.testing.assert_array_equal(result, expected)


class TestNormalizeHic:
    """Tests for hic.normalize_hic()."""

    def test_normalize_diagonal_all_ones(self):
        """Normalized Hi-C with default 'ones' method should have diagonal of 1s."""
        mat = np.random.rand(10, 10) * 100
        mat = (mat + mat.T) / 2
        result = hic.normalize_hic(mat.copy())
        np.testing.assert_array_almost_equal(np.diag(result), np.ones(10))

    def test_normalize_mean_method(self):
        """Normalized Hi-C with 'mean' method should have mean diagonal of 1."""
        mat = np.random.rand(10, 10) * 100
        mat = (mat + mat.T) / 2
        result = hic.normalize_hic(mat.copy(), method="mean")
        assert np.mean(np.diag(result)) == pytest.approx(1.0, rel=1e-5)

    def test_normalize_preserves_symmetry(self):
        """Normalization should preserve symmetry."""
        mat = np.random.rand(10, 10)
        mat = (mat + mat.T) / 2
        result = hic.normalize_hic(mat.copy())
        np.testing.assert_array_almost_equal(result, result.T)

    def test_normalize_invalid_method(self):
        """Normalization should raise error for invalid method."""
        mat = np.random.rand(5, 5)
        with pytest.raises(ValueError):
            hic.normalize_hic(mat, method="invalid")


class TestPoolingConsistency:
    """Tests that verify consistency between pooling methods.

    These tests verify the key property that:
    - Conservative pooling (pool) conserves total contacts
    - Non-conservative pooling (pool_sum) is a simple block sum
    """

    def test_conservative_vs_nonconservative_differ(self):
        """Conservative and non-conservative pooling should give different results."""
        mat = np.random.rand(8, 8)
        mat = (mat + mat.T) / 2

        conservative = hic.pool(mat, factor=2, normalize=False)
        nonconservative = hic.pool_sum(mat, factor=2, normalize=False)

        # They should differ (unless the matrix is very special)
        assert not np.allclose(conservative, nonconservative)

    def test_pooling_factor_4(self):
        """Pooling should work with factor > 2."""
        mat = np.random.rand(16, 16)
        mat = (mat + mat.T) / 2

        result = hic.pool(mat, factor=4, normalize=False)
        assert result.shape == (4, 4)

        result_sum = hic.pool_sum(mat, factor=4, normalize=False)
        assert result_sum.shape == (4, 4)
