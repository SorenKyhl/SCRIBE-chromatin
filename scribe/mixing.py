"""Mixing diagnostics for tuning the Monte Carlo move set.

Sweeps per second is not the quantity to optimize: what matters is how many
*independent* samples of the observables a run produces per CPU-second. For a
stationary trace with integrated autocorrelation time ``tau`` (in samples),

    Var(time average) = tau * Var(x) / n,      ESS = n / tau,

with the convention ``tau = 1 + 2 * sum_{t>=1} rho(t)`` (an uncorrelated trace
has tau = 1). Efficiency is then ``ESS / seconds``.

Estimators here:

- :func:`integrated_time`: Sokal's automatic windowing on the FFT
  autocorrelation function, averaged over chains (as in emcee).
- :func:`batch_means_time`: an independent cross-check from the variance of
  block averages.
- :func:`replica_time`: the variance of per-chain means across independent
  replicas; assumption-light but noisy (chi-square with m-1 dof).
- :func:`multivariate_ess`: Vats, Flegal & Jones (2019) multivariate ESS of a
  vector of observables, on its leading principal components.

The engine writes the traces these need when the config has
``"dump_mixing_stats": true`` (``structure.traj`` and ``move_stats.txt``, in
addition to the usual energy and observable trajectories); see
:func:`load_traces` and :func:`load_move_stats`.
"""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import pandas as pd

ENERGY_COLUMNS = ["bonded", "plaid", "diagonal", "boundary", "total"]
#: relative standard deviation below which a column counts as constant (the
#: engine prints 6 decimals, so parsed constants can have std ~1e-15)
CONSTANT_RTOL = 1e-9
STRUCTURE_COLUMNS = ["rg2", "ree2", "r2half"]


# ---------------------------------------------------------------- estimators


def autocorr_function(x: np.ndarray) -> np.ndarray:
    """Normalized autocorrelation function of a 1-D trace (FFT, biased)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if _is_constant(x):
        raise ValueError("trace has zero variance")
    f = np.fft.rfft(x - x.mean(), n=2 * n)
    acf = np.fft.irfft(f * np.conjugate(f))[:n]
    return acf / acf[0]


def _is_constant(x: np.ndarray, axis=None) -> np.ndarray:
    scale = np.maximum(np.abs(np.mean(x, axis=axis)), 1.0)
    return np.std(x, axis=axis) <= CONSTANT_RTOL * scale


def _as_chains(x) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[None, :]
    if x.ndim != 2:
        raise ValueError("expected a trace (n,) or chains (m, n)")
    return x


def integrated_time(x, c: float = 5.0, return_window: bool = False):
    """Integrated autocorrelation time, in samples.

    ``x`` is one trace (n,) or several chains of equal length (m, n); the
    chains' autocorrelation functions are averaged before windowing. The
    window M is the smallest M with M >= c * tau(M) (Sokal 1997). The estimate
    is reliable only when the trace is much longer than tau (n >~ 50 tau);
    check ``n / tau``. If no window satisfies the criterion, the window is the
    whole trace and the estimate is a lower bound.
    """
    chains = _as_chains(x)
    rho = np.mean([autocorr_function(ch) for ch in chains], axis=0)
    return integrated_time_from_acf(rho, c, return_window)


def integrated_time_from_acf(rho: np.ndarray, c: float = 5.0, return_window: bool = False):
    """Sokal-windowed integrated time from a normalized autocorrelation
    function (e.g. the mean of per-chain :func:`autocorr_function` results).

    Lets callers compute each chain's ACF once and reuse it, e.g. to
    bootstrap over chains by averaging resampled ACFs.
    """
    taus = 2.0 * np.cumsum(rho) - 1.0
    m = np.arange(len(taus)) < c * taus
    window = int(np.argmin(m)) if not m.all() else len(taus) - 1
    tau = float(taus[window])
    return (tau, window) if return_window else tau


def batch_means_time(x, n_batches: int | None = None) -> float:
    """Integrated time from the variance of block averages, in samples.

    tau ~= b * Var(block means) / Var(x) for block length b >> tau. Default
    block count: floor(sqrt(n)) per chain.
    """
    chains = _as_chains(x)
    m, n = chains.shape
    k = n_batches or max(int(np.sqrt(n)), 2)
    b = n // k
    if b < 1:
        raise ValueError("trace too short for batch means")
    blocks = chains[:, : k * b].reshape(m, k, b).mean(axis=2)
    # deviations from each chain's own mean, so replica offsets don't count
    blocks = blocks - blocks.mean(axis=1, keepdims=True)
    var_blocks = (blocks**2).sum() / (m * (k - 1))
    var_x = chains.var(axis=1, ddof=1).mean()
    return float(b * var_blocks / var_x)


def replica_time(x) -> float:
    """Integrated time from the spread of per-chain means, in samples.

    tau ~= n * Var(chain means) / Var(x). Needs independent, equilibrated
    chains; with m chains the estimate has m - 1 degrees of freedom.
    """
    chains = _as_chains(x)
    m, n = chains.shape
    if m < 2:
        raise ValueError("need at least two chains")
    var_means = chains.mean(axis=1).var(ddof=1)
    var_x = chains.var(axis=1, ddof=1).mean()
    return float(n * var_means / var_x)


def ess(x, c: float = 5.0) -> float:
    """Effective sample size of one trace or several chains (total)."""
    chains = _as_chains(x)
    return chains.size / integrated_time(chains, c)


def multivariate_ess(x, n_components: int | None = 10, n_batches: int | None = None) -> float:
    """Multivariate ESS (Vats, Flegal & Jones 2019) of a vector observable.

    ``x`` has shape (n, p) for one chain or (m, n, p) for several. The
    observables are standardized and projected on their ``n_components``
    leading principal components (None: all), because with many correlated
    columns the batch-means covariance is singular. Then

        mESS = N * (det(Lambda) / det(Sigma)) ** (1/p),

    with Lambda the sample covariance and Sigma the batch-means estimate of
    the asymptotic covariance. Columns with zero variance are dropped.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim == 2:
        x = x[None]
    m, n, p = x.shape
    flat = x.reshape(m * n, p)
    sd = flat.std(axis=0)
    keep = ~_is_constant(flat, axis=0)
    z = (x[..., keep] - flat[:, keep].mean(axis=0)) / sd[keep]
    if n_components is not None and n_components < z.shape[-1]:
        _, _, vt = np.linalg.svd(z.reshape(m * n, -1), full_matrices=False)
        z = z @ vt[:n_components].T
    q = z.shape[-1]
    k = n_batches or max(int(np.sqrt(n)), 2)
    b = n // k
    if m * (k - 1) <= q:
        raise ValueError(f"{m * k} batches are too few for {q} components")
    blocks = z[:, : k * b].reshape(m, k, b, q).mean(axis=2)
    blocks = blocks - blocks.mean(axis=1, keepdims=True)
    sigma = b * np.einsum("mki,mkj->ij", blocks, blocks) / (m * (k - 1))
    zc = z - z.mean(axis=1, keepdims=True)
    lam = np.einsum("mni,mnj->ij", zc, zc) / (m * (n - 1))
    _, logdet_lam = np.linalg.slogdet(lam)
    _, logdet_sig = np.linalg.slogdet(sigma)
    return float(m * n * np.exp((logdet_lam - logdet_sig) / q))


# ------------------------------------------------------------------- loading


def _read_traj(path: Path) -> np.ndarray:
    df = pd.read_csv(path, sep=r"\s+", header=None)
    return df.dropna(axis=1).to_numpy(dtype=float)


def load_traces(data_out) -> dict[str, np.ndarray]:
    """Per-sample traces written by one engine run (columns, without sweep).

    Returns a dict with ``"sweep"`` and whichever groups exist:
    ``"energy"`` (bonded, plaid, diagonal, boundary, total), ``"obs"``
    (plaid observables), ``"diag_obs"`` and ``"structure"`` (Rg^2, Ree^2,
    R^2(N/2), from ``dump_mixing_stats``).
    """
    d = Path(data_out)
    out: dict[str, np.ndarray] = {}
    for key, name in [("energy", "energy.traj"), ("obs", "observables.traj"),
                      ("diag_obs", "diag_observables.traj"),
                      ("structure", "structure.traj")]:
        f = d / name
        if f.exists() and f.stat().st_size > 0:
            a = _read_traj(f)
            out.setdefault("sweep", a[:, 0])
            out[key] = a[:, 1:]
    return out


def load_move_stats(data_out) -> pd.DataFrame:
    """Per-move-type attempts, acceptances, beads moved, summed squared
    displacement and time (``move_stats.txt``), plus derived columns."""
    df = pd.read_csv(Path(data_out) / "move_stats.txt", sep="\t", index_col=0)
    df["acceptance"] = df["accepts"] / df["attempts"].where(df["attempts"] > 0)
    df["mean_beads"] = df["beads_moved"] / df["accepts"].where(df["accepts"] > 0)
    # expected squared jump distance per attempted move and per second
    df["esjd_per_move"] = df["sq_disp_nm2"] / df["attempts"].where(df["attempts"] > 0)
    df["esjd_per_second"] = df["sq_disp_nm2"] / df["seconds"].where(df["seconds"] > 0)
    return df


# ------------------------------------------------------------------- summary


def observable_times(chains: list[dict[str, np.ndarray]], burn: int = 0,
                     c: float = 5.0) -> pd.DataFrame:
    """Integrated time (in samples) of every observable column.

    ``chains`` are :func:`load_traces` results of independent replicas of one
    condition; ``burn`` samples are dropped from each. Columns that are
    constant in any chain are skipped. One row per column, with the Sokal,
    batch-means and (for >= 2 chains) replica estimates and n/tau.
    ``converged`` is False when the Sokal window reached the end of the
    trace (tau is then a lower bound: run longer).
    """
    n = min(len(ch["sweep"]) for ch in chains) - burn
    if n <= 0:
        raise ValueError("burn-in longer than the traces")
    rows = []
    for group in ("energy", "obs", "diag_obs", "structure"):
        if not all(group in ch for ch in chains):
            continue
        ncol = chains[0][group].shape[1]
        for j in range(ncol):
            x = np.array([ch[group][burn:burn + n, j] for ch in chains])
            if np.any(_is_constant(x, axis=1)):
                continue
            name = {"energy": ENERGY_COLUMNS, "structure": STRUCTURE_COLUMNS}.get(group)
            col = name[j] if name and j < len(name) else str(j)
            tau, window = integrated_time(x, c, return_window=True)
            rows.append({
                "group": group, "column": col, "tau": tau, "window": window,
                "converged": bool(window < n - 1),
                "tau_batch": batch_means_time(x),
                "tau_replica": replica_time(x) if len(chains) > 1 else np.nan,
                "n_over_tau": n / tau,
            })
    return pd.DataFrame(rows)


def summarize(times: pd.DataFrame, sample_interval: int,
              seconds_per_sweep: float | None = None) -> dict:
    """Condense :func:`observable_times` into the tuning objectives.

    Times are converted from samples to sweeps (``sample_interval`` =
    ``dump_stats_frequency``) and, with ``seconds_per_sweep``, to seconds.
    Reports the worst (max) time per group and over the maxent observables
    (``obs`` + ``diag_obs``), and the worst over everything.
    """
    out: dict = {}
    t = times.assign(tau_sweeps=times["tau"] * sample_interval)
    for group, g in t.groupby("group"):
        out[f"tau_max_{group}"] = float(g["tau_sweeps"].max())
        out[f"tau_median_{group}"] = float(g["tau_sweeps"].median())
    for col in STRUCTURE_COLUMNS + ["total"]:
        r = t[t["column"] == col]
        if len(r):
            out[f"tau_{col}"] = float(r["tau_sweeps"].iloc[0])
    maxent = t[t["group"].isin(["obs", "diag_obs"])]
    if len(maxent):
        out["tau_max_maxent"] = float(maxent["tau_sweeps"].max())
    out["tau_max_all"] = float(t["tau_sweeps"].max())
    out["min_n_over_tau"] = float(t["n_over_tau"].min())
    out["all_converged"] = bool(t["converged"].all())
    if seconds_per_sweep:
        out["seconds_per_sweep"] = seconds_per_sweep
        for k in [k for k in out if k.startswith("tau_")]:
            out[k.replace("tau_", "seconds_", 1)] = out[k] * seconds_per_sweep
    return out


# ------------------------------------------------------------------- running


def run_chain(root, config: dict, seqs, init_xyz, seed: int, sweeps: int,
              name: str = "data_out", **overrides) -> Path:
    """Run one chain from a stored configuration; returns its output folder.

    ``init_xyz`` is an .xyz configuration (e.g. a previous run's final
    ``output.xyz`` frame) that the run starts from, so mixing is measured
    without burn-in from a random coil. ``overrides`` are config keys, e.g.
    move-set parameters (``n_pivot=20, exp_decay_pivot=50.0``). Mixing stats
    are switched on. ``root`` must not exist.
    """
    from scribe.scribe_sim import ScribeSim

    root = Path(root)
    cfg = copy.deepcopy(config)
    cfg.update({"load_configuration": True,
                "load_configuration_filename": "init.xyz",
                "nSweeps": int(sweeps), "seed": int(seed),
                "dump_mixing_stats": True})
    cfg.update(overrides)
    sim = ScribeSim(root=str(root), config=cfg, seqs=seqs, randomize_seed=False)
    sim.setup()
    write_last_frame(init_xyz, root / "init.xyz")
    sim.run(name)
    return root / name


def write_last_frame(xyz, dest) -> None:
    """Write the last frame of a (possibly multi-frame) .xyz file to ``dest``."""
    if Path(dest).exists() and Path(dest).resolve() == Path(xyz).resolve():
        raise ValueError("source and destination are the same file")
    lines = Path(xyz).read_text().splitlines()
    frame = int(lines[0]) + 2
    if len(lines) % frame:
        raise ValueError(f"{xyz}: not a sequence of {frame - 2}-bead frames")
    Path(dest).write_text("\n".join(lines[-frame:]) + "\n")
