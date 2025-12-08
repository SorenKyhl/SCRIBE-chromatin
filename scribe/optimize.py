import copy
import functools
import os.path as osp

import numpy as np
from bayes_opt import BayesianOptimization
from scipy import optimize
from sklearn.metrics import mean_squared_error

from scribe import analysis, default, utils
from scribe.scribe_sim import ScribeSim

"""
module for optimizing simulation parameters

Context:
    The diagonal interactions determined by maximum entropy optimization are sensitive
    to the size of the grid used in TICG simulations. This is especially true for diagonal
    interaction parameters governing beads separated by short contour lengths (small s).

    Often it is useful to modify the grid size in such a way that the prbability of
    contact between nearest neighbors along the chain (i.e. p(s=1)) for a chain with no
    nonbonded interactions is identical to the probability observed in experiment. In
    this case, the maximum entropy method will need to make only minor
    adjustments to the diagonal interactions at short contour lengths.

    In the converse situation, large attractive or repulsive interactions are needed to recapitulate
    the experimental p(s) curve for small s, and sometimes the experimental probabilities are
    not achievable, because nearest neighbor beads cannot be brought into more (or less)
    frequent contact simply becuase the grid spacing is too large or too small relative
    to the bond length. In this case, the maximum entropy method does not converge.

    To ameliroate these issues, the grid size can be optimized prior to maximum entropy,
    to ensure stable behavior of the maximum entropy optimization.
"""


def nearest_neighbor_contact_error(grid_bond_ratio, sim_engine, gthic):
    """calculate the difference between simulated and experimental nearest-neighbor contact
    probability.

    Args:
        grid_bond_ratio [float]: grid size expressed as a ratio of the bond length.
        sim_engine [ScribeSim]: object to call simulation runs
        gthic [ndarray]: ground truth hic

    Returns:
        error [float]: difference between simulated and experimental contact frequency
                        for nearest-neighbor beads: p_simulated(1) - p_experiment(1)
    """
    try:
        nearest_neighbor_contact_error.counter += 1
    except AttributeError:
        nearest_neighbor_contact_error.counter = 1

    output = f"iteration{nearest_neighbor_contact_error.counter}"
    print(f"optimizing grid size, {output}")
    sim_engine.config["grid_size"] = sim_engine.config["bond_length"] * grid_bond_ratio

    sim_engine.run(output)

    output_dir = osp.join(sim_engine.root, output)
    sim_analysis = analysis.SimulationResult(output_dir, maxent_analysis=False)
    p1_sim = sim_analysis.d[1]
    p1_exp = analysis.get_diagonal(gthic)[1]
    error = p1_sim - p1_exp
    print(f"error, {error}")
    return error


def nearest_neighbor_contact_error_bond(grid_bond_ratio, sim_engine, gthic):
    """calculate the difference between simulated and experimental nearest-neighbor contact
    probability.

    Args:
        grid_bond_ratio [float]: grid size expressed as a ratio of the bond length.
        sim_engine [ScribeSim]: object to call simulation runs
        gthic [ndarray]: ground truth hic

    Returns:
        error [float]: difference between simulated and experimental contact frequency
                        for nearest-neighbor beads: p_simulated(1) - p_experiment(1)
    """
    try:
        nearest_neighbor_contact_error.counter += 1
    except AttributeError:
        nearest_neighbor_contact_error.counter = 1

    output = f"iteration{nearest_neighbor_contact_error.counter}"
    print(f"optimizing grid size, {output}")
    sim_engine.config["bond_length"] = sim_engine.config["grid_size"] / grid_bond_ratio

    sim_engine.run(output)

    output_dir = osp.join(sim_engine.root, output)
    sim_analysis = analysis.SimulationResult(output_dir, maxent_analysis=False)
    p1_sim = sim_analysis.d[1]
    p1_exp = analysis.get_diagonal(gthic)[1]
    error = p1_sim - p1_exp
    print(f"error, {error}")
    return error


def optimize_grid_size(config, gthic, low_bound=0.5, high_bound=1.5, root="optimize-grid-size"):
    """tune grid size until simulated nearest neighbor contact probability
    is equal to the same probability derived from the ground truth hic matrix.

    Args:
        config [json]: simulation config file. grid_size will be overwritten,
                            so it's initial value doesn't matter
        gthic [ndarray]: ground truth hic map. used as the target for optimization
    Returns:
        optimal_grid_size [float]: size of grid cell that matches first neighbor contact probability
    """
    config = copy.deepcopy(config)

    config["load_configuration"] = False
    config["nonbonded_on"] = False
    config["load_bead_types"] = False

    sim_engine = ScribeSim(root, config, seqs=None)

    try:
        result = optimize.brentq(
            nearest_neighbor_contact_error,
            low_bound,
            high_bound,
            args=(sim_engine, gthic),
            xtol=1e-3,
            maxiter=10,
        )
    except ValueError:
        low_bound /= 2
        high_bound *= 2
        result = optimize.brentq(
            nearest_neighbor_contact_error,
            low_bound,
            high_bound,
            args=(sim_engine, gthic),
            xtol=1e-3,
            maxiter=10,
        )

    # optimizer returns the grid_to_bond ratio... have to convert to real units.
    optimal_grid_size = result * config["bond_length"]
    return optimal_grid_size


def optimize_bond_length(config, gthic, low_bound=0.5, high_bound=1.5, root="optimize-bond-length"):
    """tune grid size until simulated nearest neighbor contact probability
    is equal to the same probability derived from the ground truth hic matrix.

    Args:
        config [json]: simulation config file. grid_size will be overwritten,
                            so it's initial value doesn't matter
        gthic [ndarray]: ground truth hic map. used as the target for optimization
    Returns:
        optimal_grid_size [float]: size of grid cell that matches first neighbor contact probability
    """
    config = copy.deepcopy(config)

    config["load_configuration"] = False
    config["nonbonded_on"] = False
    config["load_bead_types"] = False

    sim_engine = ScribeSim(root, config, seqs=None)

    try:
        result = optimize.brentq(
            nearest_neighbor_contact_error_bond,
            low_bound,
            high_bound,
            args=(sim_engine, gthic),
            xtol=1e-3,
            maxiter=10,
        )
    except ValueError:
        low_bound /= 2
        high_bound *= 2
        result = optimize.brentq(
            nearest_neighbor_contact_error_bond,
            low_bound,
            high_bound,
            args=(sim_engine, gthic),
            xtol=1e-3,
            maxiter=10,
        )

    # optimizer returns the grid_to_bond ratio... have to convert to real units.
    optimal_bond_length = config["grid_size"] / result
    return optimal_bond_length


def stiffness_root_error(hic, gthic):
    """calculate error for the purposes of optimizing stiffness.

    for lack of a better option, just try to match the p(s)
    at the  4th bead for N=1024
    """
    nbeads = len(hic)
    index = int(4 * nbeads / 1024) - 1
    # index = int(6*nbeads/1024) - 1
    # index = int(100*nbeads/1024) - 1
    return analysis.get_diagonal(hic)[index] - analysis.get_diagonal(gthic)[index]

    gthic_d = analysis.get_diagonal(gthic)
    hic_d = analysis.get_diagonal(hic)

    id_start = int(200 * nbeads / 1024)
    id_end = int(800 * nbeads / 1024)
    return np.mean(gthic_d[id_start:id_end]) - np.mean(hic_d[id_start:id_end])
    # return np.mean(analysis.get_diagonal(hic) - analysis.get_diagonal(gthic))


def stiffness_bayes_error(hic, gthic):
    """error used for bayesian optimization"""
    sim = analysis.get_diagonal(hic)
    exp = analysis.get_diagonal(gthic)
    return mean_squared_error(sim, exp, squared=False)


def simulate_stiffness_error(k_angle, sim_engine, gthic, method):
    """simulate and calculate the difference between simulated and experimental p(s)
    diagonal probability as a function of the chain stiffness defined by the angle potential.

    Args:
        k_angle [float]: grid size expressed as a ratio of the bond length.
        sim_engine [ScribeSim]: object to call simulation runs
        gthic [ndarray]: ground truth hic

    Returns:
        error [float]: difference between simulated and experimental p(s) diagonal probability
    """
    try:
        simulate_stiffness_error.counter += 1
    except AttributeError:
        simulate_stiffness_error.counter = 1

    output = f"iteration{simulate_stiffness_error.counter}"
    print(f"optimizing chain stiffness, {output}")
    sim_engine.config["k_angle"] = k_angle

    sim_engine.run(output)

    output_dir = osp.join(sim_engine.root, output)
    sim_analysis = analysis.SimulationResult(output_dir, maxent_analysis=False)

    if method == "bayes":
        # bayes method maximizes, so return negative of error
        error = -stiffness_bayes_error(sim_analysis.hic, gthic)
        print(f"bayes error, {error}")
    else:
        error = stiffness_root_error(sim_analysis.hic, gthic)
        print(f"root error, {error}")
    return error


def optimize_stiffness(
    config, gthic, low_bound=0, high_bound=1, method="notbayes", root="optimize-stiffness"
):
    """tune angle stiffness until simulated p(s) diagonal probabity
    is equal to the same probability derived from the ground truth hic matrix.

    Args:
        config [json]: simulation config file. grid_size will be overwritten,
                            so it's initial value doesn't matter
        gthic [ndarray]: ground truth hic map. used as the target for optimization
    Returns:
        optimal_k_angle[float]: angle stiffness
    """
    config = copy.deepcopy(config)
    config["load_configuration"] = False
    config["nonbonded_on"] = False
    config["load_bead_types"] = False

    sim_engine = ScribeSim(root, config, seqs=None, gthic=gthic)

    if method == "bayes":
        pbounds = {"k_angle": (low_bound, high_bound)}

        black_box = functools.partial(
            simulate_stiffness_error, sim_engine=sim_engine, gthic=gthic, method="bayes"
        )

        optimizer = BayesianOptimization(
            f=black_box,
            pbounds=pbounds,
            verbose=2,
            random_state=1,
        )

        optimizer.maximize(
            init_points=2,
            n_iter=20,
        )

        return optimizer.max["params"]["k_angle"]  # pyright: ignore

    else:
        # use rootfinding
        method = "notbayes"
        result = optimize.brentq(
            simulate_stiffness_error,
            low_bound,
            high_bound,
            args=(sim_engine, gthic, method),
            xtol=1e-2,
            maxiter=10,
        )
        return result


if __name__ == "__main__":
    """directory set up:
    config.json
    experimental_hic.npy [optional - if doesn't exist, will load from default pipeline.
    """
    config = utils.load_json("config.json")
    config["nSweeps"] = 10000

    gthic_path = "experimental_hic.npy"
    if osp.exists(gthic_path):
        gthic = np.load(gthic_path)
    else:
        pipe = default.data_pipeline.resize(config["nbeads"])
        gthic = pipe.load_hic(default.HCT116_hic)

    optimal_grid_size = optimize_grid_size(config, gthic)
    print("optimal grid size is:", optimal_grid_size)
