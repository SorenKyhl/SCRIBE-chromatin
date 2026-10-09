"""Tests for the mixing diagnostics in scribe.mixing.

The estimators are checked on AR(1) processes, whose integrated
autocorrelation time is known exactly: tau = (1 + phi) / (1 - phi).
"""

import numpy as np
import pytest

from scribe import mixing


def ar1(phi, n, m=1, p=None, seed=0):
    """m chains of length n (optionally p independent components each)."""
    rng = np.random.default_rng(seed)
    shape = (m, n) if p is None else (m, n, p)
    eps = rng.standard_normal(shape)
    x = np.empty(shape)
    x[:, 0] = eps[:, 0] / np.sqrt(1 - phi**2)  # start in the stationary state
    for t in range(1, n):
        x[:, t] = phi * x[:, t - 1] + eps[:, t]
    return x


def tau_ar1(phi):
    return (1 + phi) / (1 - phi)


@pytest.mark.parametrize("phi", [0.0, 0.5, 0.9])
def test_integrated_time_ar1(phi):
    x = ar1(phi, 20000, m=8)
    assert mixing.integrated_time(x) == pytest.approx(tau_ar1(phi), rel=0.1)


def test_integrated_time_single_trace():
    x = ar1(0.8, 50000)[0]
    assert mixing.integrated_time(x) == pytest.approx(tau_ar1(0.8), rel=0.15)


def test_batch_means_time_ar1():
    x = ar1(0.9, 20000, m=8, seed=1)
    assert mixing.batch_means_time(x) == pytest.approx(tau_ar1(0.9), rel=0.2)


def test_replica_time_ar1():
    # many short replicas: the estimate has m - 1 degrees of freedom
    x = ar1(0.9, 2000, m=200, seed=2)
    assert mixing.replica_time(x) == pytest.approx(tau_ar1(0.9), rel=0.25)


def test_replica_time_detects_offset_chains():
    # chains offset from each other by more than their spread look slow
    x = ar1(0.0, 1000, m=10, seed=3) + np.arange(10)[:, None] * 5
    assert mixing.replica_time(x) > 100


def test_ess():
    x = ar1(0.9, 20000, m=4, seed=4)
    assert mixing.ess(x) == pytest.approx(x.size / tau_ar1(0.9), rel=0.1)


def test_multivariate_ess_iid_and_correlated():
    iid = ar1(0.0, 20000, m=2, p=5, seed=5)
    assert mixing.multivariate_ess(iid) == pytest.approx(40000, rel=0.15)
    slow = ar1(0.9, 20000, m=2, p=5, seed=6)
    assert mixing.multivariate_ess(slow) == pytest.approx(40000 / tau_ar1(0.9), rel=0.2)


def test_multivariate_ess_handles_redundant_columns():
    x = ar1(0.5, 20000, m=2, p=3, seed=7)
    # duplicated and constant columns must not make the covariance singular
    x = np.concatenate([x, x, np.zeros(x.shape[:2] + (1,))], axis=-1)
    assert mixing.multivariate_ess(x, n_components=3) == pytest.approx(
        40000 / tau_ar1(0.5), rel=0.2)


def test_zero_variance_raises():
    with pytest.raises(ValueError):
        mixing.integrated_time(np.ones(100))


def test_observable_times_and_summary():
    m, n = 4, 5000
    chains = []
    for k in range(m):
        chains.append({
            "sweep": np.arange(1, n + 1) * 10.0,
            "energy": np.stack([ar1(0.5, n, seed=10 + k)[0]] * 5, axis=1),
            "obs": np.stack([ar1(0.9, n, seed=20 + k)[0], np.ones(n)], axis=1),
            "structure": np.stack([ar1(0.95, n, seed=30 + k)[0]] * 3, axis=1),
        })
    times = mixing.observable_times(chains)
    # the constant observable column is skipped
    assert set(times["group"]) == {"energy", "obs", "structure"}
    assert len(times[times["group"] == "obs"]) == 1
    s = mixing.summarize(times, sample_interval=10, seconds_per_sweep=1e-3)
    assert s["tau_max_maxent"] == pytest.approx(10 * tau_ar1(0.9), rel=0.2)
    assert s["tau_rg2"] == pytest.approx(10 * tau_ar1(0.95), rel=0.25)
    assert s["tau_max_all"] == pytest.approx(s["tau_rg2"])
    assert s["seconds_rg2"] == pytest.approx(s["tau_rg2"] * 1e-3)


def test_write_last_frame(tmp_path):
    frames = []
    for f in range(3):
        frames += ["2", "atoms", f"0\t{f}.0\t0\t0", f"1\t{f}.5\t0\t0"]
    src = tmp_path / "out.xyz"
    src.write_text("\n".join(frames) + "\n")
    mixing.write_last_frame(src, tmp_path / "init.xyz")
    assert (tmp_path / "init.xyz").read_text().splitlines() == frames[-4:]
