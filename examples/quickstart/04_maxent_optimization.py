#!/usr/bin/env python3
"""
SCRIBE Quickstart Example 4: Maximum Entropy Optimization

This script demonstrates how to:
1. Load experimental Hi-C and polymer sequences
2. Set up maximum entropy optimization
3. Learn χ parameters that reproduce the Hi-C contact map

Prerequisites:
    - Run 02_process_data.py first

Output:
    - maxent_output/: Optimization output directory
    - maxent_output/chis.npy: Learned χ parameters
    - maxent_output/SCC.txt: Stratum-adjusted correlation over iterations
"""

import numpy as np

from pylib import default
from pylib.maxent import Maxent


def main():
    # Load experimental Hi-C contact map (training target)
    try:
        hic_experimental = np.load("experimental_hic.npy")
    except FileNotFoundError:
        print("experimental_hic.npy not found.")
        print("Run 02_process_data.py first.")
        return

    # Load polymer sequences from ChIP-seq (defines bead identities)
    try:
        sequences = np.load("chipseq_sequences.npy")
    except FileNotFoundError:
        print("chipseq_sequences.npy not found.")
        print("Run 02_process_data.py first.")
        return

    print(f"Hi-C shape: {hic_experimental.shape}")
    print(f"Sequences shape: {sequences.shape}")

    # Set up maximum entropy optimization
    config = default.config.copy()
    params = default.params.copy()

    print("\nMaximum entropy parameters:")
    print(f"  - Iterations: {params.get('iterations', 'default')}")
    print(f"  - Gamma (learning rate): {params.get('gamma', 'default')}")

    me = Maxent(
        root="maxent_output",
        params=params,
        config=config,
        seqs=sequences,  # Input: epigenetic sequences
        gthic=hic_experimental,  # Target: experimental Hi-C to match
    )

    # Run optimization: learns χ parameters that reproduce Hi-C
    print("\nStarting maximum entropy optimization...")
    print("This will iteratively:")
    print("  1. Run simulations with current χ")
    print("  2. Compare predicted vs experimental Hi-C")
    print("  3. Update χ to reduce difference")
    print()

    me.fit()

    print("\nOptimization complete!")
    print("Output saved to: maxent_output/")
    print("  - chis.npy: Learned interaction parameters")
    print("  - SCC.txt: Correlation trajectory")


if __name__ == "__main__":
    main()
