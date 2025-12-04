import numpy as np

from scribe import default, epilib
from scribe.plot_contactmap import plot_contactmap
from scribe.pysim import Pysim

# load simulation configuration file
config = default.config
# config = utils.load_json("config.json")

config["nSweeps"] = 10000  # modify desired config settings

# calculate bead labels "sequences" from ground truth hic
gthic = np.load("experimental_hic.npy")
seqs = epilib.get_sequences(gthic, k=10)

# ... or load them if they're already on disk
"""
dir_containing_seqs = "."
with utils.cd(dir_containing_seqs):
    utils.load_sequences(config)
"""

# construct simulation object
sim = Pysim(root=".", config=config, seqs=seqs, mkdir=False)

# run single short production simulation
sim.run("output_dir")  # if output_dir is none, defaults to "data_out"

# run equilibration, followed by production simulations
"""
sim.run_eq(eq_sweeps=10000, prod_sweeps=50000, parallel=1)
# if parallel = n (>1), it will run n production simulations in "core{n}"
# and aggregate the results in "production_out"
"""

# optional analysis-plot contactmap
plot_contactmap("output_dir")
