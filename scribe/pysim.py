import time
from collections.abc import Sequence
from multiprocessing import Process
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

from scribe import utils
from scribe.pyticg import Sim
from scribe.utils import cat, cd, copy_last_snapshot

PathLike = Union[Path, str]
ArrayLike = Union[list, np.ndarray]

"""
pysim
"""


class Pysim:
    """wrapper around pyticg"""

    def __init__(
        self,
        root: PathLike,
        config: dict,
        seqs: ArrayLike | None,
        gthic=None,
        randomize_seed: bool = True,
        mkdir: bool = True,
        setup_needed: bool = True,
    ):
        self.set_root(root, mkdir)
        self.set_config(config)
        self.seqs = seqs
        self.setup_needed = setup_needed  # should be true unless instantiated using from_directory

        if gthic is not None:
            np.save(self.root / "experimental_hic.npy", gthic)

        if randomize_seed:
            self.randomize_seed()

    @classmethod
    def from_directory(cls, sim_dir: PathLike, new_root: PathLike | None = None):
        """construct a simulation object from a directory that's already set up
        (i.e. contains config and sequence files)

        dir: simulation directory from which to initialize
        root: new simulation root directory
        """
        sim_dir = Path(sim_dir).absolute()
        config = utils.load_json(sim_dir / "config.json")
        if config["load_bead_types"]:
            with utils.cd(sim_dir):
                seqs = utils.load_sequences(config)
        else:
            seqs = None

        if new_root is None:
            return cls(sim_dir, config, seqs, mkdir=False, setup_needed=False)
        else:
            return cls(new_root, config, seqs, mkdir=True, setup_needed=True)

    def set_root(self, root: PathLike, mkdir: bool):
        """set the root of the simulation. and (optionally) create the directory

        Args:
            root: simulation root, where all simulation outputs will be stored
            mkdir: if true, will make a new directory at the root path
        """

        self.root = Path(root)
        if mkdir:
            self.root.mkdir(exist_ok=False)

    def set_config(self, config: dict | str):
        """set the config, either by assignment or by loading from a json file"""
        if isinstance(config, dict):
            self.config = config
        else:
            self.config = utils.load_json(config)

    def randomize_seed(self):
        """sets seed for random number generator used in Monte Carlo simulation"""
        self.config["seed"] = np.random.randint(10000)

    def setup(self):
        """write simulation inputs in simulation root directory, if setup_needed"""
        # needs to happen regardless of setup_needed setting, because seed is randomized
        utils.write_json(self.config, Path(self.root, "config.json"))

        # further setup needed - write seuences to directory
        if self.setup_needed:
            if self.seqs is None:
                return

            self.seqs = np.array(self.seqs)
            if self.seqs.ndim > 1:
                for i, seq in enumerate(self.seqs):
                    self.write_sequence(seq, Path(self.root, f"pcf{i+1}.txt"))
            else:
                self.write_sequence(self.seqs, Path(self.root, "pcf1.txt"))

    def flatten_chis(self):
        """restructure chi parameters into a 1-dimensional list
        Returns: [1 x n] list = [plaid_chis, diag_chis]
        """
        if self.config["plaid_on"]:
            indices = np.triu_indices(self.config["nspecies"])
            plaid_chis = np.array(self.config["chis"])[indices]
            diag_chis = np.array(self.config["diag_chis"])
            flat_chis = np.hstack((plaid_chis, diag_chis))
            return flat_chis
        else:
            # no plaid chis, just return diag
            return self.config["diag_chis"]

    def chis_to_matrix(self, flat_chis):
        """restructure 1-D list of plaid chis into matrix
        Returns: [n x n] chi matrix
        """
        nspecies = self.config["nspecies"]
        X = np.zeros((nspecies, nspecies))
        X[np.triu_indices(nspecies)] = flat_chis
        X = X + X.T - np.diag(np.diag(X))
        return X

    def split_chis(self, allchis):
        """splits 1-D list of chis into: plaid_chis, diag_chis"""
        nplaidchis = len(allchis) - len(self.config["diag_chis"])
        plaid_chis_flat, diag_chis = np.split(allchis, [nplaidchis])
        return plaid_chis_flat, diag_chis

    def set_chis(self, allchis):
        """takes 1d vector of all chis and updates config chi parameters"""
        if self.config["plaid_on"]:
            plaid_chis_flat, diag_chis = self.split_chis(allchis)
            self.config["chis"] = self.chis_to_matrix(plaid_chis_flat).tolist()
            self.config["diag_chis"] = diag_chis.tolist()
        else:
            self.config["diag_chis"] = allchis

    def load_observables(self, jacobian=False):
        """load observable trajectories from simulation output.
        return mean of observables throughout the simulation,
        and (optionally) the jacobian of the observable matrix
        """
        if self.data_out is None:
            raise ValueError("data out member variable has not been initialized")

        obs_files = []
        if self.config["plaid_on"]:
            obs_files.append(self.root / self.data_out / "observables.traj")
        obs_files.append(self.root / self.data_out / "diag_observables.traj")

        df_total = pd.DataFrame()
        for file in obs_files:
            df = pd.read_csv(file, delimiter="\t", header=None)
            df = df.dropna(axis=1)
            df = df.drop(df.columns[0], axis=1)
            df_total = pd.concat((df_total, df), axis=1)

        self.obs = df_total.mean().values

        if jacobian:
            self.jac = df_total.cov().values  # pyright: ignore
            return self.obs, self.jac
        else:
            return self.obs

    def run(self, name=None):
        """run simulation"""
        self.data_out = name
        self.setup()
        with cd(self.root):
            if name:
                engine = Sim(self.data_out)  # output to dir: name
            else:
                engine = Sim()  # output to dir: data_out
                self.data_out = "data_out"  # don't set earlier...
                # in this case, use default constructor so engine.run() pipes to stdout
            engine.run()

    def run_eq(
        self,
        equilibrium_sweeps: int,
        production_sweeps: int,
        parallel_simulations: int = 1,
    ):
        """run equilibration followed by production simulation
        the production run can be executed in parallel.

        Args:
            equilibrium_sweeps: number of equilibrium simulation sweeps
            production_sweeps: number of production simulation sweeps (per core)
            parallel_simulations: number of parallel production simulations to execute
        """
        equilibration_dir = "equilibration"
        production_dir = "production_out"

        # equilibration
        self.config["nSweeps"] = equilibrium_sweeps
        self.config["load_configuration"] = False
        self.run(equilibration_dir)

        # production. copy the last structure to initialize production simulations
        equilibrated_structure = "equilibrated.xyz"
        copy_last_snapshot(
            xyz_in=(self.root / equilibration_dir / "output.xyz"),
            xyz_out=(self.root / equilibrated_structure),
            nbeads=self.config["nbeads"],
        )

        self.config["nSweeps"] = production_sweeps
        self.config["load_configuration"] = True
        self.config["load_configuration_filename"] = equilibrated_structure

        if parallel_simulations == 1:
            self.run(production_dir)
        else:
            self.run_parallel(production_dir, parallel_simulations)

    def run_parallel(self, name: str, cores: int):
        """run production, using parallel simulations with different initial seeds

        each parallel core dumps output to their own directory. after completion,
        the individual simulation data are aggregated into a final production directory

        Args:
            name: name of output file for aggregated simulation data
            cores: number of parallel simulations to execute
        """
        self.data_out = name
        processes = []
        # create processes
        for i in range(cores):
            target = self.run
            args = (f"core{i}",)
            processes.append(Process(target=target, args=args))

        # run simulations
        for p in processes:
            self.randomize_seed()
            p.start()  # run simulation for process p
            sleeptime = int(
                self.config["nbeads"] / 1000
            )  # larger simulations take longer to move over
            time.sleep(sleeptime)  # engine needs time to load config.json before next iteration

        # wait to join
        for p in processes:
            p.join()

        # getting segfaults, possibly due to all processes retaining their memory. does this fix it?
        for p in processes:
            p.terminate()

        self.aggregate_production_files()

    def aggregate_production_files(self):
        """aggregate simulation data from each core into final production folder.

        each parallel core dumps its output to a separate folder
        after all simulations are finished...
        - concatenate all observables and energies into one trajectory file
        - sum all contact maps to produce one aggregate contact map
        """
        Path(self.root / "production_out").mkdir()

        # aggregate files are concatenated together to form one trajectory
        aggregate_files = [
            "observables.traj",
            "diag_observables.traj",
            "energy.traj",
            "output.xyz",
        ]
        for file in aggregate_files:
            cat(
                self.root / "production_out" / file,
                self.root.glob("core*/" + str(file)),
            )

        # contact maps are summed together to form one contact map
        contact_files = list(self.root.glob("core*/contacts.txt"))
        self.combine_contactmaps(
            contact_files, output_file=self.root / "production_out/contacts.txt"
        )

        production = Pysim(self.root / "production_out", self.config, self.seqs, mkdir=False)
        production.setup()  # save config, seqs to production_out

        # TODO: add delete old files option

    def combine_contactmaps(
        self, contact_files: Sequence[PathLike], output_file: PathLike | None = None
    ):
        """combines (sums) multiple contact maps from separate parallel simulations

        Args:
            contact_files: list of files containing contact maps to be combined
            output_file: name for output combined contact map.
        Returns:
            contact map, or saves contact map to file if specified
        """
        combined = np.array(
            pd.read_csv(contact_files[0], header=None, sep=" ", dtype=int)
        )  # much faster than np.loadtxt
        for file in contact_files[1:]:
            tt = np.array(
                pd.read_csv(file, header=None, sep=" ", dtype=int)
            )  # much faster than np.loadtxt
            combined += tt
            # combined += np.loadtxt(file)

        if output_file:
            np.savetxt(output_file, combined, fmt="%d", delimiter=" ")
        else:
            return combined

    def write_sequence(self, sequence: ArrayLike, path: PathLike):
        """write sequence of polymer bead types to disk at specified path"""
        np.savetxt(path, sequence, fmt="%.8f", delimiter=" ")


def reduce_sequences(sim, inds):
    """
    returns a Pysim object with a subset of the original sequences
    """
    inds_complement = list(set(np.arange(sim.config["nspecies"])) - set(inds))  # indices to delete

    sim.config["nspecies"] = len(inds)
    sim.config["bead_type_files"] = np.delete(
        sim.config["bead_type_files"], inds_complement
    ).tolist()
    sim.config["chis"] = utils.delete_2d(sim.config["chis"], inds_complement).tolist()
    return sim
