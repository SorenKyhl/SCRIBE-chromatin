#!/usr/bin/env python3
"""
SCRIBE Full Quickstart: Data → Simulation

This script runs the complete workflow from data download to simulation:
1. Download Hi-C and ChIP-seq data (if not present)
2. Process data into polymer sequences
3. Run a forward simulation
4. Analyze and visualize results

This is a self-contained example that demonstrates the full SCRIBE pipeline
for generating 3D chromatin structures from epigenetic data.

Usage:
    python full_simulation.py

Output:
    - chipseq_sequences.npy: Processed polymer sequences
    - experimental_hic.npy: Hi-C contact map
    - simulation_output/: Simulation results and visualizations
"""

import subprocess
import sys

import numpy as np

from pylib import default
from pylib.datapipeline import DataPipeline
from pylib.plot_contactmap import plot_contactmap
from pylib.pysim import Pysim


def check_and_download_data():
    """Check if data exists, download if needed."""
    pipeline = DataPipeline(cell="HCT116_auxin", chrom=2, nbeads=1024)
    status = pipeline.status()

    if not status["hic"]["available"] or not status["chipseq"]["available"]:
        print("=" * 60)
        print("Data not found. Downloading...")
        print("=" * 60)
        print("\nThis will download ~36 GB of data.")
        print("Press Ctrl+C to cancel, or wait to continue...\n")

        try:
            import time

            time.sleep(3)
        except KeyboardInterrupt:
            print("\nDownload cancelled.")
            sys.exit(0)

        # Download data
        subprocess.run([sys.executable, "-m", "pylib.download_data", "--all"], check=True)

    return pipeline


def process_data(pipeline):
    """Process Hi-C and ChIP-seq data into simulation inputs."""
    print("\n" + "=" * 60)
    print("Step 1: Processing Data")
    print("=" * 60)

    # Load Hi-C
    print("\nLoading Hi-C contact map...")
    hic = pipeline.load_hic()
    print(f"  Shape: {hic.shape}")

    # Load ChIP-seq
    print("\nLoading ChIP-seq data...")
    sequences = pipeline.load_chipseq()
    seq_array = pipeline.load_chipseq_array()
    print(f"  Tracks: {list(sequences.keys())}")
    print(f"  Shape: {seq_array.shape}")

    # Save
    np.save("chipseq_sequences.npy", seq_array)
    np.save("experimental_hic.npy", hic)
    print("\nSaved: chipseq_sequences.npy, experimental_hic.npy")

    return hic, seq_array


def run_simulation(sequences):
    """Run forward simulation."""
    print("\n" + "=" * 60)
    print("Step 2: Running Simulation")
    print("=" * 60)

    config = default.config.copy()

    print("\nConfiguration:")
    print(f"  - Beads: {config.get('nbeads', sequences.shape[0])}")
    print(f"  - Sequence types: {sequences.shape[1]}")

    sim = Pysim(root="simulation_output", config=config, seqs=sequences)

    print("\nRunning simulation...")
    print("  - Equilibration: 10,000 sweeps")
    print("  - Production: 50,000 sweeps")
    sim.run_eq(eq_sweeps=10000, prod_sweeps=50000)

    return sim


def analyze_results(sim, experimental_hic):
    """Analyze simulation results."""
    print("\n" + "=" * 60)
    print("Step 3: Analysis")
    print("=" * 60)

    # Generate contact map visualization
    print("\nGenerating contact map...")
    plot_contactmap("simulation_output")

    # Calculate metrics
    from scipy.stats import pearsonr

    from pylib.epilib import SCC

    scc = SCC(sim.hic, experimental_hic)
    pearson_r, _ = pearsonr(sim.hic.flatten(), experimental_hic.flatten())

    print("\nComparison to experimental Hi-C:")
    print(f"  - SCC: {scc:.3f}")
    print(f"  - Pearson r: {pearson_r:.3f}")


def main():
    print("=" * 60)
    print("SCRIBE Full Quickstart: Data → Simulation")
    print("=" * 60)

    # Check/download data
    pipeline = check_and_download_data()
    print("\nData status:")
    print(f"  {pipeline.status()}")

    # Process data
    hic, sequences = process_data(pipeline)

    # Run simulation
    sim = run_simulation(sequences)

    # Analyze
    analyze_results(sim, hic)

    print("\n" + "=" * 60)
    print("Complete!")
    print("=" * 60)
    print("\nOutputs:")
    print("  - simulation_output/: Simulation results")
    print("  - simulation_output/hic.png: Contact map")
    print("  - chipseq_sequences.npy: Polymer sequences")
    print("  - experimental_hic.npy: Experimental Hi-C")


if __name__ == "__main__":
    main()
