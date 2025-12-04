#!/usr/bin/env python3
"""
SCRIBE Full Quickstart: Data → Maximum Entropy Optimization

This script runs the complete workflow from data download to parameter learning:
1. Download Hi-C and ChIP-seq data (if not present)
2. Process data into polymer sequences
3. Run maximum entropy optimization to learn χ parameters
4. Analyze convergence and learned parameters

This is a self-contained example that demonstrates how SCRIBE learns
the epigenetic interaction parameters that reproduce experimental Hi-C.

Usage:
    python full_maxent.py

Output:
    - chipseq_sequences.npy: Processed polymer sequences
    - experimental_hic.npy: Hi-C contact map
    - maxent_output/: Optimization results
    - maxent_output/chis.npy: Learned χ parameters
"""

import subprocess
import sys
from pathlib import Path

import numpy as np

from scribe import default
from scribe.data_pipeline import DataPipeline
from scribe.maxent import Maxent


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
        subprocess.run([sys.executable, "-m", "scribe.download_data", "--all"], check=True)

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


def run_maxent(sequences, experimental_hic):
    """Run maximum entropy optimization."""
    print("\n" + "=" * 60)
    print("Step 2: Maximum Entropy Optimization")
    print("=" * 60)

    config = default.config.copy()
    params = default.params.copy()

    print("\nConfiguration:")
    print(f"  - Beads: {config.get('nbeads', sequences.shape[0])}")
    print(f"  - Sequence types: {sequences.shape[1]}")
    print("\nOptimization parameters:")
    print(f"  - Iterations: {params.get('iterations', 'default')}")
    print(f"  - Gamma (learning rate): {params.get('gamma', 'default')}")

    me = Maxent(
        root="maxent_output",
        params=params,
        config=config,
        seqs=sequences,
        gthic=experimental_hic,
    )

    print("\nStarting optimization...")
    print("This will iteratively:")
    print("  1. Run simulations with current χ")
    print("  2. Compare predicted vs experimental Hi-C")
    print("  3. Update χ to reduce difference")
    print()

    me.fit()

    return me


def analyze_results():
    """Analyze optimization results."""
    print("\n" + "=" * 60)
    print("Step 3: Analysis")
    print("=" * 60)

    output_dir = Path("maxent_output")

    # Load results
    if (output_dir / "chis.npy").exists():
        chis = np.load(output_dir / "chis.npy")
        print("\nLearned χ parameters:")
        print(f"  - Shape: {chis.shape}")
        print(f"  - Final χ range: [{chis[-1].min():.3f}, {chis[-1].max():.3f}]")

    if (output_dir / "SCC.txt").exists():
        scc = np.loadtxt(output_dir / "SCC.txt")
        print("\nTraining progress (SCC):")
        print(f"  - Initial: {scc[0]:.3f}")
        print(f"  - Final: {scc[-1]:.3f}")
        print(f"  - Improvement: {scc[-1] - scc[0]:.3f}")


def main():
    print("=" * 60)
    print("SCRIBE Full Quickstart: Data → Maximum Entropy")
    print("=" * 60)

    # Check/download data
    pipeline = check_and_download_data()
    print("\nData status:")
    print(f"  {pipeline.status()}")

    # Process data
    hic, sequences = process_data(pipeline)

    # Run maximum entropy
    run_maxent(sequences, hic)

    # Analyze
    analyze_results()

    print("\n" + "=" * 60)
    print("Complete!")
    print("=" * 60)
    print("\nOutputs:")
    print("  - maxent_output/: Optimization results")
    print("  - maxent_output/chis.npy: Learned χ parameters")
    print("  - maxent_output/SCC.txt: Correlation trajectory")
    print("  - chipseq_sequences.npy: Polymer sequences")
    print("  - experimental_hic.npy: Experimental Hi-C")


if __name__ == "__main__":
    main()
