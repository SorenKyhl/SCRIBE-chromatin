"""
Minimal example: run a maxent optimization using synthetic data.
No downloads required.
"""

import numpy as np
from scribe import default
from scribe.maxent import Maxent
from scribe.data_pipeline import SyntheticDataPipeline
from scribe import epilib as ep

# Use synthetic data pipeline (no downloads needed)
pipeline = SyntheticDataPipeline(nbeads=1024)

# Start with the full default config and override only what's needed
config = default.config.copy()
config["nSweeps"] = 10

# Get default params for maxent optimization
params = default.params.copy()
params["parallel"] = 1  # Single process for quick testing
params["iterations"] = 2  # Just 2 iterations for demo
params["equilib_sweeps"] = 100
params["production_sweeps"] = 100

# Load synthetic data
gthic = pipeline.load_hic()
seqs = pipeline.load_chipseq()

# Compute goals from sequences and config
goals = ep.get_goals(gthic, seqs, config)
params["goals"] = goals

maxent = Maxent(
    root="synthetic_maxent",
    params=params,
    config=config,
    seqs=seqs,
    gthic=gthic,
)
maxent.fit()
print("Synthetic maxent optimization completed.")
