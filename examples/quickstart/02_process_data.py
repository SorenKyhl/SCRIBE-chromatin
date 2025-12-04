#!/usr/bin/env python3
"""
SCRIBE Quickstart Example 2: Data Processing

This script demonstrates how to:
1. Load Hi-C contact maps using the high-level DataPipeline
2. Load ChIP-seq data with automatic processing
3. Save data for simulation

Prerequisites:
    - Download data first: python -m scribe.download_data --all

Output:
    - chipseq_sequences.npy: Polymer sequences from ChIP-seq
    - experimental_hic.npy: Hi-C contact map
"""

import numpy as np

from scribe.datapipeline import DataPipeline


def main():
    # Create a pipeline for HCT116 cell line data
    # This automatically finds data in ~/.scribe/data/
    pipeline = DataPipeline(
        cell="HCT116_auxin",  # Cell line/condition
        chrom=2,  # Chromosome number
        nbeads=1024,  # Number of polymer beads
    )

    # Check data availability
    status = pipeline.status()
    print("Data status:", status)

    if not status["hic"]["available"]:
        print("\nHi-C data not found. Run: python -m scribe.download_data --hic")
        return

    if not status["chipseq"]["available"]:
        print("\nChIP-seq data not found. Run: python -m scribe.download_data --chipseq")
        return

    # Load Hi-C contact map (with automatic pooling and caching)
    print("\nLoading Hi-C contact map...")
    hic = pipeline.load_hic()
    print(f"Hi-C shape: {hic.shape}")

    # Load ChIP-seq data (automatically finds all bigWig files)
    print("\nLoading ChIP-seq data...")
    sequences = pipeline.load_chipseq()
    print(f"Loaded {len(sequences)} ChIP-seq tracks: {list(sequences.keys())}")

    # Convert to array and save for simulation
    seq_array = pipeline.load_chipseq_array()
    print(f"Sequence array shape: {seq_array.shape}")

    np.save("chipseq_sequences.npy", seq_array)
    np.save("experimental_hic.npy", hic)

    print("\nSaved:")
    print("  - chipseq_sequences.npy")
    print("  - experimental_hic.npy")


if __name__ == "__main__":
    main()
