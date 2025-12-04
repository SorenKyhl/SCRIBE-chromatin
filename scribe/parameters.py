import copy

import numpy as np

from scribe import default


def get_config(
    nbeads=None,
    config=default.config,
    grid_bond_ratio=0.95,
    base="gaussian-5k",
    scale="onethird",
    scale_diag_bins=False,
):
    """
    calculates physical parameters for a simulation with nbeads beads,
    if nbeads is not specified, take from config

    nbeads (int): number of beads in simulation
    config (optional):  default config file.
    grid_bond_ratio (float): grid size is defined by this ratio: bond_length / grid_size
    base = ["gaussian", "persistent", "gaussian-5k", "persistent-5k"]
        base parameters, from which to scale
        if gaussian, chain is lp = 16.5, b = 16.5
        if persistent, chain is lp = 50, b = 16.5
        if gaussian-5k, the gaussian chain is gaussian renormalized up to 5kbp/bead resolution
        if persistent-5k, the persistent chain is gaussian renormalized up to 5kbp/bead resolution
    scale (str): scaling method.
    scale_diag_bins: if on, keep the same diagonal bins, just scale them up to match the new nbeads.
                        if off, generate new bins using np.geomspace. doing this at different nbeads
                                will not necessarily generate the same nubmer of diagonal bins.
    """
    config = copy.deepcopy(config)
    nbeads_orig = config["nbeads"]

    if nbeads is None:
        nbeads = config["nbeads"]
    else:
        # scale config.chis to new size specified by nbeads
        config["chis"] = (config["nbeads"] / nbeads * np.array(config["chis"])).tolist()
        config["diag_chis"] = (config["nbeads"] / nbeads * np.array(config["diag_chis"])).tolist()

    # if config["diagonal_binning"]:
    # bin_boundaries = np.geomspace(1,config["nbeads"], 32)
    # config["diagonal_bin_boundaries"] = sorted(list(set(bin_boundaries)))
    # config["diag_chis"] = list(np.zeros(len(config["diagonal_bin_boundaries"])))
    # x = np.array(config["diagonal_bin_boundaries"])
    # config["diagonal_bin_boundaries"] = np.ceil(x * int(nbeads/config["nbeads"])).tolist()

    # if nbeads > 10241:
    #    config['contact_resolution'] = 5

    assert base in ["gaussian", "gaussian-5k", "persistent", "persistent-5k", "load"]
    if base == "gaussian":
        baseb = 16.5
        baseg = grid_bond_ratio * baseb
        baseN = 512000
        basev = 520
    elif base == "persistent":
        baseb = 100
        baseg = grid_bond_ratio * baseb
        baseN = 84480
        basev = 3151.52
    elif base == "gaussian-5k":
        baseb = 82.5
        baseg = grid_bond_ratio * baseb
        baseN = 20480
        basev = 13000
    elif base == "persistent-5k":
        baseb = 203.10
        baseg = grid_bond_ratio * baseb
        baseN = 20480
        basev = 13000
    elif base == "load":
        baseb = config["bond_length"]
        baseg = config["grid_size"]
        baseN = config["nbeads"]
        basev = config["beadvol"]
    else:
        raise ValueError(
            "base must be: ['gaussian' | 'gaussian-5k' | 'persistent' | 'persistent-5k]"
        )

    if scale == "onethird":
        factor = (nbeads / baseN) ** (-1 / 3)
        config["nbeads"] = nbeads
        config["bond_length"] = baseb * factor
        config["grid_size"] = baseg * factor
        config["beadvol"] = basev * baseN / nbeads
        config["diag_cutoff"] = nbeads
    elif scale == "gaussian":
        factor = (nbeads / baseN) ** (-1 / 2)
        config["nbeads"] = nbeads
        config["bond_length"] = baseb * factor
        config["beadvol"] = basev * baseN / nbeads
        config["diag_cutoff"] = nbeads
    else:
        raise ValueError("scale must be onethird or gaussian")

    if config["diagonal_binning"]:
        if scale_diag_bins:
            x = np.array(config["diagonal_bin_boundaries"])
            config["diagonal_bin_boundaries"] = np.ceil(x * int(nbeads / nbeads_orig)).tolist()
        else:
            bin_boundaries = np.geomspace(1, config["nbeads"], 32, dtype=int).tolist()
            config["diagonal_bin_boundaries"] = sorted(list(set(bin_boundaries)))
            config["diag_chis"] = np.zeros(len(config["diagonal_bin_boundaries"])).tolist()

    return config
