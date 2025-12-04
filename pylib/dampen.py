import numpy as np

from pylib import epilib as ep
from pylib import utils
from pylib.maxent import Maxent


def dampen(gamma, me_path="me-1024", it=0):
    it0 = ep.Sim(me_path + f"/iteration{it}/production_out")
    params = utils.load_json(me_path + "/resources/params.json")
    config = utils.load_json(me_path + f"/iteration{it}/config.json")
    gthic = np.load(me_path + "/resources/experimental_hic.npy")
    params["gamma"] = gamma

    root = f"dampen{gamma}"
    me = Maxent(root, params, config, it0.seqs, gthic, dampen_first_step=True)
    me.fit()


def main():
    # dampen(0.125, "me-1024", 0)
    dampen(1, "dampen0.125", 4)


if __name__ == "__main__":
    main()
