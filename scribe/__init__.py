"""
SCRIBE: Simulation of Chromatin with Reinforced Inference and Bayesian Estimation

A Python package for chromatin simulation and analysis.

Key classes:
    - ScribeSim: Dispatch and run polymer simulations
    - SimulationResult: Analyze completed simulation results
    - Maxent: Maximum entropy optimization of interaction parameters
    - DataPipeline: Load Hi-C and ChIP-seq data by cell type

Example:
    >>> from scribe import ScribeSim, SimulationResult
    >>> sim = ScribeSim(root="output", config=config, seqs=sequences)
    >>> sim.run_eq(equilibrium_sweeps=10000, production_sweeps=50000)
    >>> result = SimulationResult("output/production_out")
    >>> result.plot_contactmap()
"""

# Core simulation classes
from scribe.scribe_sim import ScribeSim
from scribe.analysis import SimulationResult

# High-level data loading
from scribe.data_pipeline import DataPipeline, SyntheticDataPipeline

# Optimization
from scribe.maxent import Maxent

__all__ = [
    "ScribeSim",
    "SimulationResult", 
    "DataPipeline",
    "SyntheticDataPipeline",
    "Maxent",
]
