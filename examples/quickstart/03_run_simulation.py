#!/usr/bin/env python3
"""
SCRIBE Quickstart Example 3: Run a Single Simulation

This script demonstrates how to:
1. Load polymer sequences and configuration
2. Run a forward simulation
3. Visualize the resulting contact map

Prerequisites:
    - Run 02_process_data.py first (or have chipseq_sequences.npy)

Output:
    - output/: Simulation output directory
    - output/hic.png: Contact map visualization
"""

import numpy as np

from pylib import default
from pylib.plot_contactmap import plot_contactmap
from pylib.pysim import Pysim


def main():
    # Load default configuration (contains interaction parameters χ)
    config = default.config.copy()

    # Load polymer sequences (epigenetic mark occupancies from ChIP-seq)
    try:
        sequences = np.load("chipseq_sequences.npy")
    except FileNotFoundError:
        print("chipseq_sequences.npy not found.")
        print("Run 01_process_data.py first, or download example data.")
        return

    print(f"Loaded sequences: {sequences.shape}")
    print(f"Config nbeads: {config.get('nbeads', 'not set')}")

    # Create simulation: sequences define bead identities, config defines χ parameters
    sim = Pysim(root="output", config=config, seqs=sequences)

    # Run equilibration + production to generate ensemble of 3D structures
    print("Running simulation...")
    print("  - Equilibration: 10,000 sweeps")
    print("  - Production: 50,000 sweeps")
    sim.run_eq(eq_sweeps=10000, prod_sweeps=50000)

    # Visualize the resulting contact map (averaged over ensemble)
    print("Generating contact map visualization...")
    plot_contactmap("output")

    print("\nSimulation complete!")
    print("Output saved to: output/")


if __name__ == "__main__":
    main()
