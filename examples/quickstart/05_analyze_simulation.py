#!/usr/bin/env python3
"""
SCRIBE Quickstart Example 5: Simulation Analysis

This script demonstrates how to:
1. Load a completed simulation
2. Analyze energy convergence
3. Compare predicted vs experimental contact maps
4. Calculate quantitative metrics (SCC, Pearson)

Prerequisites:
    - Run 03_run_simulation.py first

Output:
    - Analysis plots in the simulation directory
"""

import numpy as np
from scipy.stats import pearsonr

from scribe.analysis import SimulationResult, get_SCC
from scribe.analysis_pipeline import compare_analysis, sim_analysis


def main():
    # Load a completed simulation's outputs (production_out holds the
    # aggregated contacts.txt, energy.traj, etc. that SimulationResult reads)
    try:
        sim = SimulationResult("output/production_out", maxent_analysis=False)
    except Exception as e:
        print(f"Could not load simulation from 'output/production_out': {e}")
        print("Run 03_run_simulation.py first.")
        return

    print("Loaded simulation from: output/production_out")
    print(f"Contact map shape: {sim.hic.shape}")

    # Basic analysis: energy convergence, contact map visualization
    print("\nRunning basic analysis...")
    sim_analysis(sim)

    # Compare to experimental Hi-C (ground truth)
    try:
        experimental_hic = np.load("experimental_hic.npy")
        sim.gthic = experimental_hic

        print("\nComparing to experimental Hi-C...")
        compare_analysis(sim)  # Generates comparison plots

        # Quantitative metrics
        scc = get_SCC(sim.hic, experimental_hic)
        pearson_r, _ = pearsonr(sim.hic.flatten(), experimental_hic.flatten())

        print("\nQuantitative comparison:")
        print(f"  - SCC (Stratum-adjusted Correlation): {scc:.3f}")
        print(f"  - Pearson r: {pearson_r:.3f}")

    except FileNotFoundError:
        print("\nexperimental_hic.npy not found - skipping comparison analysis.")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
