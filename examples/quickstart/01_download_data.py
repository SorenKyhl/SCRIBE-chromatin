#!/usr/bin/env python3
"""
SCRIBE Quickstart Example 1: Download Data

This script demonstrates how to:
1. Check what data is available/missing
2. Download Hi-C and ChIP-seq data from ENCODE
3. Verify the download

Data is stored in ~/.scribe/data/ by default.
Set SCRIBE_DATA_DIR environment variable to use a custom location.

Output:
    ~/.scribe/data/hic/HCT116_auxin/*.hic  (~29 GB)
    ~/.scribe/data/chipseq/HCT116_hg19/*.bigWig  (~6.8 GB)
"""

import subprocess
import sys

from scribe import default
from scribe.paths import get_data_dir


def check_status():
    """Check what data is available."""
    print("=" * 60)
    print("Data Status")
    print("=" * 60)
    print(f"\nData directory: {get_data_dir()}")

    # Check Hi-C
    hic_path = default.get_hic_path("HCT116_auxin")
    if hic_path:
        print(f"\n✓ Hi-C data found: {hic_path}")
    else:
        print("\n✗ Hi-C data not found")

    # Check ChIP-seq
    chipseq_dir = default.get_chipseq_dir("HCT116")
    if chipseq_dir:
        bigwig_files = list(chipseq_dir.glob("*.bigWig"))
        print(f"✓ ChIP-seq data found: {chipseq_dir}")
        print(f"  {len(bigwig_files)} bigWig files")
    else:
        print("✗ ChIP-seq data not found")

    return hic_path is not None and chipseq_dir is not None


def download_all():
    """Download all data using the download_data module."""
    print("\n" + "=" * 60)
    print("Downloading Data")
    print("=" * 60)
    print("\nThis will download ~36 GB of data:")
    print("  - Hi-C: ~29 GB")
    print("  - ChIP-seq: ~6.8 GB")
    print("\nPress Ctrl+C to cancel...\n")

    try:
        import time

        time.sleep(3)
    except KeyboardInterrupt:
        print("\nDownload cancelled.")
        return False

    # Run the download
    subprocess.run([sys.executable, "-m", "scribe.download_data", "--all"], check=True)
    return True


def main():
    print("=" * 60)
    print("SCRIBE Quickstart: Download Data")
    print("=" * 60)

    # Check current status
    all_present = check_status()

    if all_present:
        print("\n" + "=" * 60)
        print("All data is already downloaded!")
        print("=" * 60)
        print("\nYou can proceed to the next example:")
        print("  python 02_process_data.py")
        return

    # Offer to download
    print("\n" + "-" * 60)
    response = input("Download missing data? [y/N]: ").strip().lower()

    if response == "y":
        success = download_all()
        if success:
            print("\n" + "=" * 60)
            print("Download complete!")
            print("=" * 60)
            check_status()
    else:
        print("\nTo download later, run:")
        print("  python -m scribe.download_data --all")
        print("\nOr download separately:")
        print("  python -m scribe.download_data --hic       # Hi-C only")
        print("  python -m scribe.download_data --chipseq   # ChIP-seq only")


if __name__ == "__main__":
    main()
