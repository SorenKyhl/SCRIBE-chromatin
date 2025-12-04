"""
Energy calculation utilities for SCRIBE simulations.

This module provides functions to calculate the energy matrices and terms
used in the TICG (Theoretically Informed Coarse Grain) simulation.

The energy consists of:
- Plaid (off-diagonal) interactions: E_IJ based on χ_IJ parameters
- Diagonal interactions: distance-dependent χ(d) parameters

Key matrices:
- S: Raw energy matrix from sequences and chi parameters
- E: Symmetric energy matrix (S + S.T - diag)
- D: Diagonal energy contribution from distance-dependent chi
"""

import math

import numpy as np


def calculate_all_energy(config, x, chi):
    """Calculate all energy components from sequences and parameters.

    Args:
        config: Simulation configuration dict
        x: Sequence matrix of shape (nbeads, nspecies)
        chi: Interaction parameter matrix of shape (nspecies, nspecies)

    Returns:
        tuple: (S, D, E, ED) where:
            - S: Raw plaid energy matrix
            - D: Diagonal energy matrix
            - E: Symmetric plaid energy
            - ED: Net energy (plaid + diagonal)
    """
    diag_chis_continuous = calculate_diag_chi_step(config)
    D = calculate_D(diag_chis_continuous)
    E, S = calculate_E_S(x, chi)
    ED = calculate_net_energy(S, D)
    return S, D, E, ED


def calculate_net_energy(S, D):
    """Calculate net energy by combining plaid and diagonal contributions.

    Args:
        S: Plaid energy matrix
        D: Diagonal energy matrix

    Returns:
        ED: Combined symmetric energy matrix
    """
    S_sym = (S + S.T) / 2
    SD = S_sym + D + np.diag(np.diagonal(D.copy()))
    ED = s_to_E(SD)
    return ED


def calculate_E_S(x, chi):
    """Calculate both E and S matrices from sequences and chi.

    Args:
        x: Sequence matrix of shape (nbeads, nspecies)
        chi: Interaction parameters of shape (nspecies, nspecies)

    Returns:
        tuple: (E, S) energy matrices, or (None, None) if inputs are None
    """
    if x is None or chi is None:
        return None, None
    S = calculate_S(x, chi)
    E = s_to_E(S)
    return E, S


def calculate_E(x, chi):
    """Calculate symmetric energy matrix from sequences and chi.

    Args:
        x: Sequence matrix of shape (nbeads, nspecies)
        chi: Interaction parameters of shape (nspecies, nspecies)

    Returns:
        E: Symmetric energy matrix
    """
    S = calculate_S(x, chi)
    E = s_to_E(S)
    return E


def s_to_E(S):
    """Convert raw energy matrix S to symmetric energy E.

    Args:
        S: Raw (upper triangular) energy matrix

    Returns:
        E: Symmetric energy matrix (S + S.T - diag(S))
    """
    if S is None:
        return None

    return S + S.T - np.diag(np.diagonal(S).copy())


def calculate_S(x, chi):
    """Calculate raw energy matrix S = x @ chi @ x.T.

    The chi matrix should be upper triangular.

    Args:
        x: Sequence matrix of shape (nbeads, nspecies)
        chi: Upper triangular interaction parameters

    Returns:
        S: Raw energy matrix of shape (nbeads, nbeads)
    """
    assert len(chi.shape) == 2, f"chi has shape {chi.shape}"
    m, k = x.shape
    assert m > k, f"x has shape {x.shape}, try x.T"
    # zero lower triangle (double check)
    chi = np.triu(chi)

    try:
        S = x @ chi @ x.T
    except ValueError:
        print("x", x, x.shape)
        print("chi", chi, chi.shape)
        raise
    return S


def calculate_diag_chi_step(config, diag_chi=None):
    """Expand binned diagonal chi values to per-bead step function.

    The diagonal chi parameters are stored in bins for efficiency.
    This function expands them to a value for each genomic distance.

    Args:
        config: Simulation configuration
        diag_chi: Binned diagonal chi values (uses config if None)

    Returns:
        Array of length nbeads with chi value for each distance
    """
    m = config["nbeads"]
    if diag_chi is None:
        diag_chi = config["diag_chis"]
    diag_bins = len(diag_chi)

    if diag_bins == m:
        return diag_chi

    if "diag_start" in config.keys():
        diag_start = config["diag_start"]
    else:
        diag_start = 0

    if "diag_cutoff" in config.keys():
        diag_cutoff = config["diag_cutoff"]
    else:
        diag_cutoff = m

    if "dense_diagonal_on" in config.keys():
        dense = config["dense_diagonal_on"]
    else:
        dense = False

    if dense:
        if "n_small_bins" in config.keys():
            n_small_bins = config["n_small_bins"]
            small_binsize = config["small_binsize"]
            big_binsize = config["big_binsize"]
        else:
            # soren compatibility
            n_small_bins = int(config["dense_diagonal_loading"] * diag_bins)
            n_big_bins = diag_bins - n_small_bins
            m_eff = diag_cutoff - diag_start  # number of beads with nonzero interaction
            dividing_line = m_eff * config["dense_diagonal_cutoff"]
            small_binsize = int(dividing_line / (n_small_bins))
            big_binsize = int((m_eff - dividing_line) / n_big_bins)
    else:
        raise NotImplementedError("non-dense diagonal chis not implemented")

    diag_chi_step = np.zeros(m)
    for d in range(diag_cutoff):
        if d < diag_start:
            continue
        d_eff = d - diag_start
        if dense:
            dividing_line = n_small_bins * small_binsize

            if d_eff > dividing_line:
                bin = n_small_bins + math.floor((d_eff - dividing_line) / big_binsize)
            else:
                bin = math.floor(d_eff / small_binsize)
        else:
            binsize = m / diag_bins
            bin = int(d_eff / binsize)
        diag_chi_step[d] = diag_chi[bin]

    return diag_chi_step


def calculate_D(diag_chi_continuous):
    """Build diagonal energy matrix from per-distance chi values.

    Creates a symmetric matrix where D[i,j] = chi(|i-j|).

    Args:
        diag_chi_continuous: Array of chi values for each genomic distance

    Returns:
        D: Symmetric diagonal energy matrix of shape (m, m)
    """
    m = len(diag_chi_continuous)
    D = np.zeros((m, m))
    for d in range(m):
        rng = np.arange(m - d)
        D[rng, rng + d] = diag_chi_continuous[d]
        D[rng + d, rng] = diag_chi_continuous[d]

    return D
