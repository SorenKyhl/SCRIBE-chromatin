"""
Minimal example: run a single forward simulation using synthetic data.
No downloads required.
"""

from scribe import default
from scribe.scribe_sim import ScribeSim
from scribe.data_pipeline import SyntheticDataPipeline

# Use synthetic data pipeline (no downloads needed)
pipeline = SyntheticDataPipeline(nbeads=1024)

# Start with the full default config and override only what's needed
config = default.config.copy()
config["nSweeps"] = 10

sim = ScribeSim(
    root="synthetic_simulation",
    config=config,
    seqs=pipeline.load_chipseq(),
)

sim.run()
print("Synthetic forward simulation completed.")
