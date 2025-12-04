#!/usr/bin/env python3
"""
SCRIBE Quickstart Example 6: High-Level Pipeline

This script demonstrates how to:
1. Use the Pipeline class for systematic experiments
2. Sweep over different sequence representations
3. Compare results across parameter settings

Prerequisites:
    - Run 02_process_data.py first

Output:
    - pc_1/, pc_2/, ..., pc_10/: Output directories for each run
"""

import functools

import numpy as np

from scribe import default
from scribe import epilib as ep
from scribe.maxent_pipeline import MaxentPipeline


def main():
    # Load data and config
    try:
        experimental_hic = np.load("experimental_hic.npy")
    except FileNotFoundError:
        print("experimental_hic.npy not found.")
        print("Run 02_process_data.py first.")
        return

    config = default.config.copy()
    params = default.params.copy()

    print(f"Hi-C shape: {experimental_hic.shape}")
    print("\nSweeping over number of principal components (k=1 to 10)")
    print("Each run will derive 2k sequences from Hi-C via PCA\n")

    # Sweep over different numbers of principal components
    for k in range(1, 11):
        print(f"--- Running with k={k} (2k={2*k} sequences) ---")

        seqs_method = functools.partial(ep.get_sequences, k=k)
        pipe = MaxentPipeline(
            name=f"pc_{k}",
            gthic=experimental_hic,
            config=config,
            params=params,
            seqs_method=seqs_method,  # Derives 2k sequences from Hi-C PCA
        )
        pipe.fit()  # Runs full maximum entropy optimization

    print("\nSweep complete!")
    print("Results saved to: pc_1/, pc_2/, ..., pc_10/")


if __name__ == "__main__":
    main()
