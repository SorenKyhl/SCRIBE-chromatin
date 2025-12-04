from os import PathLike

import numpy as np

from pylib import utils


class Config:
    def __init__(self, config):
        self.set_config(config)

    def __getitem__(self, key):
        "operates just like a normal config file"
        return self.config[key]

    def set_config(self, config):
        """load config from path
        config: [dict, filepath]
        """
        if isinstance(config, dict):
            self.config = config
        elif isinstance(config, PathLike):
            self.config = utils.load_json(config)
        else:
            raise ValueError("usage: Config( dict | str )")

    def flatten_chis(self):
        """returns 1-D list: [plaid_chis, diag_chis]"""
        indices = np.triu_indices(self.config["nspecies"])
        plaid_chis = np.array(self.config["chis"])[indices]
        diag_chis = np.array(self.config["diag_chis"])
        flat_chis = np.hstack((plaid_chis, diag_chis))
        return flat_chis

    def chis_to_matrix(self, flat_chis):
        """converts flattened plaid chis to matrix"""
        nspecies = self.config["nspecies"]
        X = np.zeros((nspecies, nspecies))
        X[np.triu_indices(nspecies)] = flat_chis
        X = X + X.T - np.diag(np.diag(X))
        return X

    def split_chis(self, allchis):
        """splits all chis into: [plaid_chis, diag_chis]"""
        nplaidchis = len(allchis) - len(self.config["diag_chis"])
        plaid_chis_flat, diag_chis = np.split(allchis, [nplaidchis])
        return plaid_chis_flat, diag_chis

    def set_chis(self, allchis):
        """takes 1d vector of all chis and updates config chi parameters"""
        plaid_chis_flat, diag_chis = self.split_chis(allchis)
        self.config["chis"] = self.chis_to_matrix(plaid_chis_flat).tolist()
        self.config["diag_chis"] = diag_chis.tolist()

    def scale_chis(self, scale_factor):
        return NotImplementedError
