import logging

import matplotlib.pyplot as plt
import numpy as np

from scribe import default, epilib, hic, parameters, utils
from scribe.config import Config
from scribe.ideal_chain import ideal_chain_simulation
from scribe.maxent import Maxent
from scribe.optimize import optimize_bond_length, optimize_grid_size, optimize_stiffness
from scribe.pysim import Pysim


def plot_stiffness_error(ideal_small, ideal_large, gthic_big, pool_fn):
    """
    ideal_small: small sim with stiffness
    ideal_large: big sim without stiffness
    gthic_big: ground truth hic at big scale
    """
    factor = int(len(gthic_big) / len(ideal_small.hic))
    gthic_small = pool_fn(gthic_big, factor)

    factor = int(len(ideal_large.hic) / len(ideal_small.hic))
    id_pooled = pool_fn(ideal_large.hic, factor)

    ratio_pooled = epilib.get_diagonal(id_pooled) / epilib.get_diagonal(gthic_small)
    ratio_sim = ideal_small.d / epilib.get_diagonal(gthic_small)

    ratio2 = ratio_sim / ratio_pooled

    error = np.mean(ratio2[100:700])

    # optional:
    # epilib.plot_diagonal(gthic2k, 'k', scale='loglog')
    # epilib.plot_diagonal(id2k.d, label="ideal 2k")

    # epilib.plot_diagonal(gthic1k, 'k')
    epilib.plot_diagonal(ideal_small.d, label="stiff optimal", scale="loglog")
    epilib.plot_diagonal(id_pooled, label="pooled ideal 2k")
    plt.legend()
    plt.title("p(s) stiff vs pooled")
    plt.savefig("stiff_p.png")
    plt.close()

    plt.figure()
    epilib.plot_diagonal(ratio_sim, label="stiff", scale="loglog")
    epilib.plot_diagonal(ratio_pooled, label="pooled")
    plt.legend()
    plt.title("p(s) ratio stiff vs pooled")
    plt.savefig("stiff_ratios.png")
    plt.close()

    plt.figure()
    epilib.plot_diagonal(ratio2, scale="loglog")
    epilib.plot_diagonal(np.ones(1024), "k")
    plt.title(f"ratio of ratios, mean: {error}")
    plt.savefig("stiff_ratio_ratios.png")
    plt.close()


def tune_stiffness(
    nbeads_large,
    nbeads_small,
    pool_fn,
    grid_bond_ratio,
    method,
    large_contact_pooling_factor,
    match_ideal_large_grid=False,
    optimize_bond=False,
):
    """optimize chain stiffness of small system to match p(s) curve of large system

    simualte ideal chain at large scale and pool large ideal hic down to small scale,
    tune stiffness of ideal chain at small scale so that
    diagonal probability of small simulation matches that of the pooled large simulation.

    match_ideal_large_grid: tune small simulation to match first point of pooled ideal large chain
        (this is not the same as the first point of pooled large gthic)
    """
    data_out = "data_out"

    # run large ideal simulation
    try:
        ideal_chain_large = ideal_chain_simulation(nbeads_large, grid_bond_ratio)
        if nbeads_large >= 10240:
            ideal_chain_large.config["nSweeps"] = 25000
        ideal_chain_large.config["contact_resolution"] = large_contact_pooling_factor
        ideal_chain_large.run(data_out)
        large_hic = np.loadtxt(ideal_chain_large.root / data_out / "contacts.txt")
    except FileExistsError:
        large_out = "ideal-chain-" + str(nbeads_large) + "/" + data_out
        large_hic = np.loadtxt(large_out + "/contacts.txt")

    # get config for small ideal chain simulation
    ideal_chain_small = ideal_chain_simulation(nbeads_small, grid_bond_ratio)
    small_config = ideal_chain_small.config

    factor = int(len(large_hic) / nbeads_small)
    # don't use nbeads_large, because the large_hic might be pooled during simulation
    large_ideal_hic_pooled = pool_fn(large_hic, factor)

    if match_ideal_large_grid:
        # tune small simulation to match first point of pooled ideal large chain
        # (this is not the same as the first point of pooled large gthic)

        # tune grid size
        root = "optimize-grid-match-ideal-large"
        try:
            optimal_grid_size = optimize_grid_size(small_config, large_ideal_hic_pooled, root=root)
        except FileExistsError:
            optimal_grid_size = utils.load_json(f"{root}/config.json")["grid_size"]

        small_config["grid_size"] = optimal_grid_size

    if optimize_bond:
        try:
            new_bond_length = optimize_bond_length(
                small_config, large_ideal_hic_pooled, low_bound=0.5, high_bound=1.5
            )
        except FileExistsError:
            config_bond_opt = utils.load_json("optimize-bond-length/config.json")
            new_bond_length = config_bond_opt["bond_length"]

        small_config["bond_length"] = new_bond_length

    k_angle_opt = optimize_stiffness(
        small_config, large_ideal_hic_pooled, low_bound=0, high_bound=5, method=method
    )

    if match_ideal_large_grid:
        return k_angle_opt, optimal_grid_size
    if optimize_bond:
        return k_angle_opt, new_bond_length
    else:
        return k_angle_opt


def scaleup(
    nbeads_large,
    nbeads_small,
    pool_fn,
    method="notbayes",
    pool_large=True,
    zerodiag=False,
    match_ideal_large_grid=False,
    optimize_bond=False,
    cell="HCT116_auxin",
    drop_diag_inds=None,
):
    """optimize chis on small system, and scale up parameters to large system

    in order for the chi parameters to be transferrable from the coarse system to the fine system,
    they have to have the same ideal chain behavior. In other words, an ideal chain simulation
    of the fine system - when pooled to the coarse resolution - should have the same p(s) curve as
    the coarse ideal chain system.

    First, the grid size is tuned so that p(s=1/N) matches for both the coarse and fine grain systems.


    requires tuning the grid size and stiffness at small scale,
    in order for the chi parameters to be transferrable
    """

    if pool_large:
        large_contact_pooling_factor = int(nbeads_large / nbeads_small)
    else:
        large_contact_pooling_factor = 1

    config_small = parameters.get_config(nbeads_small)
    gthic_small = hic.load_hic(nbeads_small, pool_fn, cell=cell)

    if drop_diag_inds is not None:
        plaid_chis = config_small["nspecies"]
        start = int(plaid_chis * (plaid_chis + 1) / 2)
        diag_chis = len(config_small["diag_chis"])
        all_inds = np.arange(start + diag_chis, dtype=int)
        skip = 5
        opt_inds = np.hstack((all_inds[:start], all_inds[start + skip :]))
    else:
        opt_inds = None

    if pool_large:
        if not config_small["conservative_contact_pooling"]:
            raise ValueError("conservative contact pooling must be turned on if pooling large hic")
        gthic_large = gthic_small
    else:
        gthic_large = hic.load_hic(nbeads_large, pool_fn, cell=cell)

    seqs_large = hic.load_seqs(nbeads_large, 10, cell=cell)
    seqs_small = hic.load_seqs(nbeads_small, 10, cell=cell)

    # tune grid size
    try:
        optimal_grid_size = optimize_grid_size(config_small, gthic_small)
    except FileExistsError:
        optimal_grid_size = utils.load_json("optimize-grid-size/config.json")["grid_size"]

    config_small["grid_size"] = optimal_grid_size
    grid_bond_ratio = (
        optimal_grid_size / config_small["bond_length"]
    )  # for later, when getting large sim config

    # tune stiffness
    try:
        if match_ideal_large_grid:
            k_angle_opt, small_optimal_grid_size = tune_stiffness(
                nbeads_large,
                nbeads_small,
                pool_fn,
                grid_bond_ratio,
                method,
                large_contact_pooling_factor,
                match_ideal_large_grid,
            )
        if optimize_bond:
            k_angle_opt, new_bond_length = tune_stiffness(
                nbeads_large,
                nbeads_small,
                pool_fn,
                grid_bond_ratio,
                method,
                large_contact_pooling_factor,
                match_ideal_large_grid,
                optimize_bond,
            )
        else:
            k_angle_opt = tune_stiffness(
                nbeads_large,
                nbeads_small,
                pool_fn,
                grid_bond_ratio,
                method,
                large_contact_pooling_factor,
                match_ideal_large_grid,
            )
    except FileExistsError:
        stiff_opt_config = utils.load_json("optimize-stiffness/config.json")
        k_angle_opt = stiff_opt_config["k_angle"]

        if optimize_bond:
            new_bond_length = stiff_opt_config["bond_length"]  # if manually tuning

        if match_ideal_large_grid:
            small_optimal_grid_size = stiff_opt_config["grid_size"]

    config_small["k_angle"] = k_angle_opt

    if optimize_bond:
        config_small["bond_length"] = new_bond_length

    if match_ideal_large_grid:
        logging.info(f"optimal small grid size is : {small_optimal_grid_size}")
        config_small["grid_size"] = small_optimal_grid_size

    # plot results
    final_it_stiff = utils.get_last_iteration("optimize-stiffness")
    ideal_small = epilib.Sim(final_it_stiff)
    ideal_large = epilib.Sim(f"ideal-chain-{str(nbeads_large)}/data_out")
    plot_stiffness_error(ideal_small, ideal_large, gthic_large, pool_fn)

    # maxent at small size
    goals = epilib.get_goals(gthic_small, seqs_small, config_small)
    params = default.params
    params["goals"] = goals

    me_root = "me-" + str(nbeads_small)
    try:
        me = Maxent(
            me_root, params, config_small, seqs_small, gthic_small, optimize_indices=opt_inds
        )
        me.fit()
    except FileExistsError:
        pass

    # production at large simulation size
    final_it = utils.get_last_iteration(me_root)
    config_opt = Config(final_it / "config.json")

    config_large = parameters.get_config(
        nbeads_large, config_opt.config, grid_bond_ratio=grid_bond_ratio, scale_diag_bins=True
    )

    def scale_chis(config, scaling_ratio):
        config["chis"] = (scaling_ratio * np.array(config["chis"])).tolist()
        config["diag_chis"] = (scaling_ratio * np.array(config["diag_chis"])).tolist()
        return config

    # correct chis for difference in grid size due to matching large ideal pooled
    if match_ideal_large_grid:
        scaling_ratio = (optimal_grid_size / small_optimal_grid_size) ** 3
        config_large = scale_chis(config_large, scaling_ratio)

    config_large["contact_resolution"] = large_contact_pooling_factor
    config_large["k_angle"] = 0
    config_large["angles_on"] = False

    if zerodiag:
        config_large["diag_chis"][1] = 0

    sim_large_root = f"final-{nbeads_large}"
    sim_large = Pysim(sim_large_root, config_large, seqs_large, gthic=gthic_large)
    sim_large.run_eq(10000, 50000, 7)


if __name__ == "__main__":
    scaleup(2048, 1024, hic.pool_sum)
