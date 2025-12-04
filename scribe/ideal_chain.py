import sys

from scribe import parameters
from scribe.scribe_sim import ScribeSim

""" ideal_chain
module for running ideal chain simulations
"""


def ideal_chain_simulation(nbeads: int, grid_bond_ratio: float | None = None):
    """return simulation object with only bonded interactions

    Args:
        nbeads: number of beads in simulation
        grid_bond_ratio: ratio between grid size and polymer bond length

    Returns:
        ScribeSim object for simulating ideal chain
    """
    if grid_bond_ratio is None:
        config = parameters.get_config(nbeads)
    else:
        config = parameters.get_config(nbeads, grid_bond_ratio=grid_bond_ratio)

    config["nonbonded_on"] = False
    config["load_bead_types"] = False
    sim = ScribeSim(root=f"ideal-chain-{nbeads}", config=config, seqs=None)
    return sim


if __name__ == "__main__":
    if len(sys.argv) > 1:
        nbeads = int(sys.argv[1])
    else:
        raise ValueError("usage: ideal_chain.py nbeads")

    sim = ideal_chain_simulation(nbeads)
    sim.run()
