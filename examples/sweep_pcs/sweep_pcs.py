import functools

import numpy as np

from scribe import default
from scribe import epilib as ep
from scribe.maxent_pipeline import MaxentPipeline

"""
optimizes chromosome structure for HCT116-chr2-0Mbp-102.4Mbp using 1024 beads.
uses principal component analysis to assign bead types.
sweeps through the number of PCs, from 2 to 10
"""

# config = json.load(open("config.json"))
# config['bond_length'] = 30
config = default.config

gthic = np.load("experimental_hic.npy")


params = default.params
# params = json.load(open("params.json"))
# params["equilib_sweeps"] = 10000
# params["production_sweeps"] = 50000
# params["parallel"] = 7
# params["trust_region"] = 1000
# params["iterations"] = 12

for i in range(2, 3):
    seqs_method = functools.partial(ep.get_sequences, k=i)
    name = "quick-" + str(i)
    pipe = MaxentPipeline(name, gthic, config, params, seqs_method=seqs_method, load_first=False)
    pipe.fit()
