#!/usr/bin/env python3
"""
SCRIBE Quickstart Example 2: Data Processing

This script demonstrates how to:
1. Load Hi-C contact maps using the high-level DataPipeline
2. Load ChIP-seq data with automatic processing
3. Save data for simulation
4. Optionally visualize the processed data

Prerequisites:
    - Download data first: python -m scribe.download_data --all

Usage:
    python 02_process_data.py           # Process data only
    python 02_process_data.py --plot    # Process and visualize

Output:
    - chipseq_sequences.npy: Polymer sequences from ChIP-seq
    - chipseq_marks.npy: Names of the ChIP-seq marks
    - experimental_hic.npy: Hi-C contact map
    - (with --plot) hic_contactmap.png, chipseq_heatmap.png, chipseq_tracks.png
"""

import argparse

import numpy as np

from scribe.data_pipeline import DataPipeline


def process_data():
    """Load and save Hi-C and ChIP-seq data."""
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
        return None, None, None

    if not status["chipseq"]["available"]:
        print("\nChIP-seq data not found. Run: python -m scribe.download_data --chipseq")
        return None, None, None

    # Load Hi-C contact map (with automatic pooling and caching)
    print("\nLoading Hi-C contact map...")
    hic = pipeline.load_hic()
    print(f"Hi-C shape: {hic.shape}")

    # Load ChIP-seq data (automatically finds all bigWig files)
    print("\nLoading ChIP-seq data...")
    sequences = pipeline.load_chipseq()
    mark_names = list(sequences.keys())
    print(f"Loaded {len(sequences)} ChIP-seq tracks: {mark_names}")

    # Convert to array and save for simulation
    seq_array = pipeline.load_chipseq_array()
    print(f"Sequence array shape: {seq_array.shape}")

    np.save("chipseq_sequences.npy", seq_array)
    np.save("chipseq_marks.npy", np.array(mark_names))
    np.save("experimental_hic.npy", hic)

    print("\nSaved:")
    print("  - chipseq_sequences.npy")
    print("  - chipseq_marks.npy")
    print("  - experimental_hic.npy")

    return hic, seq_array, mark_names


def visualize_data(hic, seqs, marks):
    """Create visualization plots of the processed data."""
    import matplotlib.pyplot as plt
    from scribe import epilib

    print("\nGenerating plots...")

    # Plot 1: Hi-C contact map using epilib
    epilib.plot_contactmap(hic)
    plt.title("Experimental Hi-C Contact Map (chr2)")
    plt.xlabel("Genomic position (100kb bins)")
    plt.ylabel("Genomic position (100kb bins)")
    plt.savefig("hic_contactmap.png", dpi=150, bbox_inches="tight")
    print("  Saved: hic_contactmap.png")
    plt.close()

    # Plot 2: ChIP-seq tracks as heatmap
    fig, ax = plt.subplots(figsize=(14, 4))

    # Normalize sequences for visualization
    seqs_norm = seqs.copy()
    for i in range(len(seqs_norm)):
        s = seqs_norm[i]
        seqs_norm[i] = (s - s.min()) / (s.max() - s.min() + 1e-10)

    im = ax.imshow(seqs_norm, aspect="auto", cmap="viridis")
    ax.set_xlabel("Genomic position (100kb bins)")
    ax.set_ylabel("ChIP-seq track")
    ax.set_yticks(range(len(marks)))
    ax.set_yticklabels(marks)
    ax.set_title("ChIP-seq Signal Tracks (normalized)")
    plt.colorbar(im, ax=ax, label="Normalized signal")
    plt.savefig("chipseq_heatmap.png", dpi=150, bbox_inches="tight")
    print("  Saved: chipseq_heatmap.png")
    plt.close()

    # Plot 3: ChIP-seq tracks as line plots
    fig, axes = plt.subplots(len(marks), 1, figsize=(14, 2 * len(marks)), sharex=True)
    if len(marks) == 1:
        axes = [axes]

    colors = plt.cm.tab10(np.linspace(0, 1, len(marks)))

    for i, ax in enumerate(axes):
        ax.fill_between(range(seqs.shape[1]), seqs[i], alpha=0.7, color=colors[i])
        ax.plot(seqs[i], linewidth=0.5, color=colors[i])
        ax.set_ylabel(marks[i], fontsize=10)
        ax.set_xlim(0, seqs.shape[1])
        ax.set_ylim(0, seqs[i].max() * 1.1)

    axes[-1].set_xlabel("Genomic position (100kb bins)")
    fig.suptitle("ChIP-seq Signal Tracks", fontsize=14)
    plt.tight_layout()
    plt.savefig("chipseq_tracks.png", dpi=150, bbox_inches="tight")
    print("  Saved: chipseq_tracks.png")
    plt.close()

    # Print statistics
    print("\n--- Data Statistics ---")
    print(f"Hi-C:")
    print(f"  Shape: {hic.shape}")
    print(f"  Min: {hic.min():.4f}, Max: {hic.max():.4f}, Mean: {hic.mean():.4f}")
    print(f"  Non-zero fraction: {(hic > 0).sum() / hic.size:.2%}")

    print(f"\nChIP-seq:")
    print(f"  Shape: {seqs.shape}")
    for i, mark in enumerate(marks):
        print(
            f"  {mark}: min={seqs[i].min():.3f}, max={seqs[i].max():.3f}, mean={seqs[i].mean():.3f}"
        )


def main():
    parser = argparse.ArgumentParser(description="Process Hi-C and ChIP-seq data for SCRIBE")
    parser.add_argument("--plot", action="store_true", help="Generate visualization plots")
    args = parser.parse_args()

    # Process data
    hic, seqs, marks = process_data()

    if hic is None:
        return

    # Optionally visualize
    if args.plot:
        visualize_data(hic, seqs, marks)


if __name__ == "__main__":
    main()
