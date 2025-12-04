import functools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams.update({"font.size": 18})


import scribe.maxent as me
from scribe import epilib as ep


class MaxentPipeline:
    """
    Pipeline for running maximum entropy optimizations.
    
    Args:
        name [str]: name of output directory
        gthic [ndarray]: ground truth hi-c contact map: optimization target
        config [dir]: simulation configuration
        params [dir]: maxent optimization parameters
        load_first [bool]: if true, load sequences and goals rather than compute them
        seqs_method [callable]: operates on gthic to return bead type sequences
        goals_method [callable]: operates on gthic, seqs, and config to return maxent goals
    """

    def __init__(
        self,
        name,
        gthic,
        config,
        params,
        load_first=False,
        seqs_method=functools.partial(ep.get_sequences, k=10),  # default is pca with k=10
        goals_method=ep.get_goals,
        analysis_on=True,
        initial_chis=None,
    ):
        self.name = name
        self.root = Path.cwd() / self.name
        self.gthic = gthic
        self.config = config
        self.params = params
        self.seqs_method = seqs_method
        self.goals_method = goals_method
        self.load_first = load_first
        self.analysis_on = analysis_on
        self.initial_chis = initial_chis

    def get_seqs(self):
        if self.load_first:
            self.seqs = self.load_sequences()
        else:
            self.seqs = self.seqs_method(self.gthic)

    def get_goals(self):
        if self.load_first:
            self.goals = self.load_goals()
        else:
            self.goals = self.goals_method(self.gthic, self.seqs, self.config)

    def fit(self):
        self.get_seqs()
        self.get_goals()
        self.params["goals"] = self.goals

        optimizer = me.Maxent(
            root=self.root,
            params=self.params,
            config=self.config,
            seqs=self.seqs,
            gthic=self.gthic,
            analysis_on=self.analysis_on,
            initial_chis=self.initial_chis,
        )

        optimizer.fit()

    def load_sequences(self):
        seqs = []
        for file in self.config["bead_type_files"]:
            print("loading", file)
            seqs.append(np.loadtxt(file))
        seqs = np.array(seqs)
        return seqs

    def load_goals(self):
        goals_plaid = np.loadtxt("obj_goal.txt")
        goals_diag = np.loadtxt("obj_goal_diag.txt")
        return np.hstack((goals_plaid, goals_diag))
