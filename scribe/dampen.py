"""
Dampening utilities for maximum entropy optimization.

This module provides functions to restart maximum entropy optimization
with a different learning rate (gamma) to fine-tune convergence.

When optimization oscillates or overshoots, reducing gamma can help
the parameters converge more smoothly to the optimal values.
"""

import numpy as np

from scribe import analysis
from scribe import utils
from scribe.maxent import Maxent


def dampen(gamma, me_path="me-1024", it=0):
    """Restart maximum entropy optimization with a new learning rate.

    This function loads a previous maxent run and continues optimization
    with a different gamma (learning rate) value. Useful for fine-tuning
    when the original gamma causes oscillation.

    Args:
        gamma: New learning rate for optimization
        me_path: Path to the previous maxent run directory
        it: Iteration number to restart from

    Example:
        >>> # If optimization at gamma=1.0 oscillates, try dampening:
        >>> dampen(0.5, "maxent_output", it=5)
    """
    it0 = analysis.SimulationResult(me_path + f"/iteration{it}/production_out")
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
