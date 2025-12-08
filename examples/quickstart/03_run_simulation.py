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

import shutil
from pathlib import Path

import numpy as np

from scribe import default
from scribe.plot_contactmap import plot_contactmap
from scribe.scribe_sim import ScribeSim


def main():
    # Load default configuration (contains interaction parameters χ)
    config = default.config.copy()

    # Load polymer sequences (epigenetic mark occupancies from ChIP-seq)
    try:
        sequences = np.load("chipseq_sequences.npy")
    except FileNotFoundError:
        print("chipseq_sequences.npy not found.")
        print("Run 02_process_data.py first, or download example data.")
        return

    print(f"Loaded sequences: {sequences.shape}")
    nspecies, nbeads = sequences.shape
    print(f"Sequences: {nspecies} species, {nbeads} beads")

    # Update config to match our sequences
    config["nspecies"] = nspecies
    config["nbeads"] = nbeads
    
    # Resize chi matrix to match nspecies (default is 10x10, we have 6 species)
    # Initialize with zeros - these would be optimized via maxent
    config["chis"] = [[0.0] * nspecies for _ in range(nspecies)]
    
    # Update bead type files list to match nspecies
    config["bead_type_files"] = [f"pcf{i+1}.txt" for i in range(nspecies)]
    
    print(f"Config: nspecies={config['nspecies']}, nbeads={config['nbeads']}")

    # Clean up previous output if it exists
    output_dir = Path("output")
    if output_dir.exists():
        print(f"Removing existing output directory: {output_dir}")
        shutil.rmtree(output_dir)

    # Create simulation: sequences define bead identities, config defines χ parameters
    sim = ScribeSim(root="output", config=config, seqs=sequences)

    # Run equilibration + production to generate ensemble of 3D structures
    print("Running simulation...")
    print("  - Equilibration: 10,000 sweeps")
    print("  - Production: 50,000 sweeps")
    sim.run_eq(equilibrium_sweeps=10000, production_sweeps=50000)

    # Visualize the resulting contact map (averaged over ensemble)
    print("Generating contact map visualization...")
    plot_contactmap("output/production_out")

    print("\nSimulation complete!")
    print("Output saved to: output/")


if __name__ == "__main__":
    main()
