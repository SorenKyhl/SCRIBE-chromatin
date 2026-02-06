#!/usr/bin/env python3
"""
Download sample data for SCRIBE.

This script downloads Hi-C and ChIP-seq data files required for running
SCRIBE simulations and examples.

Usage:
    python -m scribe.download_data [--hic] [--chipseq] [--all] [--status] [--output-dir DIR]

Examples:
    # Check which data files are available
    python -m scribe.download_data --status

    # Download all data (Hi-C + ChIP-seq, ~36 GB total)
    python -m scribe.download_data --all

    # Download Hi-C data only (~29 GB)
    python -m scribe.download_data --hic

    # Download ENCODE ChIP-seq data only (6 marks, ~6.8 GB)
    python -m scribe.download_data --chipseq-encode

    # Download all histone ChIP-seq marks (12 marks, ~12.9 GB)
    python -m scribe.download_data --chipseq-histone

Note:
    Data files are large (several GB). Ensure you have sufficient disk space.

Data Sources:
    - Hi-C data: GEO (GSE104333)
    - ChIP-seq data: ENCODE
"""

import argparse
import os
import sys
from pathlib import Path

import requests
from tqdm import tqdm

from scribe.paths import get_data_dir

# =============================================================================
# Data file registry
# =============================================================================

# Hi-C files from GEO/4DN Data Portal
# GEO supplementary files follow pattern: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSExxxxx/suppl/
HIC_FILES = {
    "HCT116_auxin": {
        "description": "HCT116 cohesin-depleted (Rao et al. 2017)",
        "geo_accession": "GSE104333",
        "filename": "GSE104333_Rao-2017-treated_6hr_combined_30.hic",
        "size_gb": 29.2,
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104333/suppl/GSE104333_Rao-2017-treated_6hr_combined_30.hic",
    },
}

# ChIP-seq files from ENCODE
# HCT116 histone modifications aligned to hg19
# URLs from ENCODE portal: https://www.encodeproject.org/
# All files are signal p-value bigWig from pooled replicates

def _encode_mark(accession: str, experiment: str, size_mb: int) -> dict:
    """Build a mark entry from its ENCODE accession."""
    return {
        "accession": accession,
        "experiment": experiment,
        "url": f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.bigWig",
        "size_mb": size_mb,
    }

# Single registry of all HCT116 histone ChIP-seq marks (hg19, signal p-value)
_HCT116_MARKS = {
    "H3K27ac":  _encode_mark("ENCFF445BLD", "ENCSR661KMA", 1020),
    "H3K27me3": _encode_mark("ENCFF977LFP", "ENCSR810BDB", 1200),
    "H3K36me3": _encode_mark("ENCFF121SMM", "ENCSR091QXP", 1100),
    "H3K4me1":  _encode_mark("ENCFF205ZML", "ENCSR161MXP", 1150),
    "H3K4me2":  _encode_mark("ENCFF736AWH", "ENCSR794ULT", 1000),
    "H3K4me3":  _encode_mark("ENCFF176NSX", "ENCSR333OPW", 1034),
    "H3K79me2": _encode_mark("ENCFF954SUA", "ENCSR494CCN", 1000),
    "H3K9ac":   _encode_mark("ENCFF960LSV", "ENCSR093SHE", 1000),
    "H3K9me2":  _encode_mark("ENCFF679MXB", "ENCSR555LYM", 1000),
    "H3K9me3":  _encode_mark("ENCFF014WPW", "ENCSR179BUC", 1360),
    "H4K20me1": _encode_mark("ENCFF757LRS", "ENCSR474DOV", 1000),
    "H2AFZ":    _encode_mark("ENCFF570WCD", "ENCSR227XNT", 1000),
}

# Which marks belong to which set
_ENCODE_MARK_NAMES = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"]
_ALL_MARK_NAMES = list(_HCT116_MARKS.keys())

CHIPSEQ_FILES = {
    "HCT116_hg19": {
        "description": f"HCT116 ENCODE histone modifications ({len(_ENCODE_MARK_NAMES)} marks)",
        "files": {m: _HCT116_MARKS[m] for m in _ENCODE_MARK_NAMES},
    },
    "HCT116_hg19_all": {
        "description": f"HCT116 all histone modifications — ENCODE ({len(_ALL_MARK_NAMES)} marks)",
        "files": {m: _HCT116_MARKS[m] for m in _ALL_MARK_NAMES},
    },
}


# =============================================================================
# Download utilities
# =============================================================================


def download_file(url: str, dest_path: Path, desc: str = "Downloading") -> bool:
    """
    Download a file with progress bar.

    Args:
        url: URL to download from
        dest_path: Destination file path
        desc: Description for progress bar

    Returns:
        True if successful, False otherwise
    """
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()

        total_size = int(response.headers.get("content-length", 0))

        dest_path.parent.mkdir(parents=True, exist_ok=True)

        with open(dest_path, "wb") as f:
            with tqdm(total=total_size, unit="B", unit_scale=True, desc=desc) as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
                    pbar.update(len(chunk))

        return True

    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return False


def print_manual_instructions():
    """Print instructions for manually downloading data."""
    print("\n" + "=" * 70)
    print("SCRIBE Data Download Instructions")
    print("=" * 70)

    data_dir = get_data_dir()
    print(f"\nData directory: {data_dir}")
    print("\nDue to licensing and size, some files must be downloaded manually.")

    print("\n" + "-" * 70)
    print("Hi-C Data")
    print("-" * 70)
    for name, info in HIC_FILES.items():
        print(f"\n{name}: {info['description']}")
        print(f"  Size: ~{info['size_gb']} GB")
        print(f"  Destination: {data_dir / 'hic' / name / info['filename']}")
        if info.get("url"):
            print("  Auto-download: python -m scribe.download_data --hic")

    print("\n" + "-" * 70)
    print("ChIP-seq Data")
    print("-" * 70)
    for name, info in CHIPSEQ_FILES.items():
        print(f"\n{name}: {info['description']}")
        marks = list(info.get("files", {}).keys())
        print(f"  Marks: {', '.join(marks)}")
        print(f"  Destination: {data_dir / 'chipseq' / name}/")
        print("  Auto-download: python -m scribe.download_data --chipseq")

    print("\n" + "-" * 70)
    print("Environment Variable")
    print("-" * 70)
    print("""
To use a custom data location, set SCRIBE_DATA_DIR:

    export SCRIBE_DATA_DIR=/path/to/your/data

Add this to your ~/.bashrc or ~/.zshrc to make it permanent.
""")


def check_data_status():
    """Check which data files are available."""
    data_dir = get_data_dir()

    print(f"\nData directory: {data_dir}\n")
    print("Hi-C Data Status:")
    print("-" * 40)

    for name, info in HIC_FILES.items():
        hic_dir = data_dir / "hic" / name
        hic_files = list(hic_dir.glob("*.hic")) if hic_dir.exists() else []
        if hic_files:
            status = "✓ Found"
        elif info.get("url"):
            status = "✗ Missing (can auto-download)"
        else:
            status = "✗ Missing (manual download)"
        print(f"  {name}: {status}")
        if hic_files:
            for f in hic_files:
                print(f"    → {f.name}")

    print("\nChIP-seq Data Status:")
    print("-" * 40)

    for name, info in CHIPSEQ_FILES.items():
        chipseq_dir = data_dir / "chipseq" / name
        bw_files = list(chipseq_dir.glob("*.bigWig")) if chipseq_dir.exists() else []
        expected = len(info.get("files", {}))
        if len(bw_files) >= expected and expected > 0:
            status = f"✓ Found ({len(bw_files)} files)"
        elif bw_files:
            status = f"⚠ Partial ({len(bw_files)}/{expected} files)"
        elif expected > 0:
            status = "✗ Missing (can auto-download)"
        else:
            status = "✗ Missing"
        print(f"  {name}: {status}")


def download_hic_data(cell_types: list = None, force: bool = False):
    """
    Download Hi-C data files.

    Args:
        cell_types: List of cell types to download, or None for all
        force: If True, re-download even if file exists

    Returns:
        True if all downloads successful, False otherwise.
    """
    data_dir = get_data_dir()

    if cell_types is None:
        cell_types = list(HIC_FILES.keys())

    print("\n" + "=" * 70)
    print("Downloading Hi-C Data")
    print("=" * 70)

    all_success = True
    downloaded_count = 0
    skipped_count = 0
    failed_count = 0

    for name in cell_types:
        if name not in HIC_FILES:
            print(f"Unknown cell type: {name}")
            continue

        info = HIC_FILES[name]
        if not info.get("url"):
            print(f"\n{name}: No download URL available (manual download required)")
            continue

        dest_dir = data_dir / "hic" / name
        dest_path = dest_dir / info["filename"]

        if dest_path.exists() and not force:
            print(f"\n{name}: Already exists at {dest_path}")
            print("  Use --force to re-download")
            skipped_count += 1
            continue

        print(f"\nDownloading {name}: {info['description']}")
        print(f"  Size: ~{info['size_gb']} GB")
        print(f"  URL: {info['url']}")
        print(f"  Destination: {dest_path}")

        dest_dir.mkdir(parents=True, exist_ok=True)

        if download_file(info["url"], dest_path, desc=info["filename"]):
            print("  ✓ Downloaded successfully")
            downloaded_count += 1
        else:
            print("  ✗ Download failed")
            failed_count += 1
            all_success = False
            # Remove partial file if it exists
            if dest_path.exists():
                dest_path.unlink()

    # Print summary
    print("\n" + "-" * 70)
    print("Hi-C Download Summary:")
    print(f"  Downloaded: {downloaded_count}")
    print(f"  Skipped (already exist): {skipped_count}")
    print(f"  Failed: {failed_count}")
    print("-" * 70)

    return all_success


def download_chipseq_data(cell_types: list | None = None, force: bool = False) -> bool:
    """
    Download ChIP-seq data from ENCODE.

    Args:
        cell_types: List of cell types to download. If None, downloads all.
        force: If True, re-download even if files exist.

    Returns:
        True if all downloads successful, False otherwise.
    """
    data_dir = get_data_dir()

    if cell_types is None:
        cell_types = list(CHIPSEQ_FILES.keys())

    print("\n" + "=" * 70)
    print("Downloading ChIP-seq Data from ENCODE")
    print("=" * 70)

    all_success = True
    downloaded_count = 0
    skipped_count = 0
    failed_count = 0
    failed_marks = []

    for cell_type in cell_types:
        if cell_type not in CHIPSEQ_FILES:
            print(f"Unknown cell type: {cell_type}")
            continue

        info = CHIPSEQ_FILES[cell_type]
        dest_dir = data_dir / "chipseq" / cell_type
        dest_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{cell_type}: {info['description']}")

        total_size = sum(f["size_mb"] for f in info["files"].values())
        print(f"  Total size: ~{total_size} MB ({len(info['files'])} files)")

        for mark, file_info in info["files"].items():
            filename = f"{file_info['accession']}.bigWig"
            dest_path = dest_dir / filename

            if dest_path.exists() and not force:
                print(f"  {mark}: Already exists ({filename})")
                skipped_count += 1
                continue

            print(f"  Downloading {mark} ({filename})...")

            if download_file(file_info["url"], dest_path, desc=f"{mark}"):
                print("    ✓ Downloaded successfully")
                downloaded_count += 1
            else:
                print("    ✗ Download failed")
                failed_count += 1
                failed_marks.append(mark)
                all_success = False
                # Remove partial file if it exists
                if dest_path.exists():
                    dest_path.unlink()

        # Create a simple metadata file
        metadata_path = dest_dir / "metadata.tsv"
        if not metadata_path.exists() or force:
            with open(metadata_path, "w") as f:
                f.write("File accession\tAssay\tBiosample\tTarget\n")
                for mark, file_info in info["files"].items():
                    f.write(
                        f"{file_info['accession']}\tChIP-seq\t{cell_type.split('_')[0]}\t{mark}\n"
                    )
            print("  Created metadata.tsv")

    # Print summary
    print("\n" + "-" * 70)
    print("ChIP-seq Download Summary:")
    print(f"  Downloaded: {downloaded_count}")
    print(f"  Skipped (already exist): {skipped_count}")
    print(f"  Failed: {failed_count}")
    if failed_marks:
        print(f"  Failed marks: {', '.join(failed_marks)}")
        all_success = False
    print("-" * 70)

    return all_success


def setup_directory_structure():
    """Create the expected directory structure."""
    data_dir = get_data_dir()

    # Create directories
    dirs = [
        data_dir / "hic" / "HCT116_auxin",
        data_dir / "chipseq" / "HCT116_hg19",
        data_dir / "chipseq" / "HCT116_hg19_all",
        data_dir / "cache",
    ]

    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        print(f"Created: {d}")

    # Create a README in the data directory
    readme_content = """# SCRIBE Data Directory

This directory contains large data files for SCRIBE simulations.

## Directory Structure

```
data/
├── hic/                    # Hi-C contact maps (.hic format)
│   └── HCT116_auxin/       # HCT116 cohesin-depleted cells
├── chipseq/                # ChIP-seq tracks (.bigWig format)
│   └── HCT116_hg19/        # HCT116 histone modifications
└── cache/                  # Preprocessed/cached files (.npy)
```

## Downloading Data

### Hi-C Data
- HCT116: Download from GEO (GSE104333)

### ChIP-seq Data
- Download from ENCODE (https://www.encodeproject.org/)
- Search for cell type + "histone ChIP-seq"
- Download bigWig files and metadata.tsv

## See Also
- SCRIBE documentation: docs/source/installation.rst
- Run `python -m scribe.download_data --help` for more options
"""

    readme_path = data_dir / "README.md"
    with open(readme_path, "w") as f:
        f.write(readme_content)
    print(f"Created: {readme_path}")


# =============================================================================
# Main
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Download and manage SCRIBE data files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--force", action="store_true", help="Force re-download even if files exist"
    )
    parser.add_argument(
        "--all", action="store_true", help="Download all available data (Hi-C + ChIP-seq)"
    )
    parser.add_argument("--hic", action="store_true", help="Download Hi-C data only")
    parser.add_argument(
        "--chipseq-encode", action="store_true",
        help="Download ENCODE histone ChIP-seq marks (6 marks: HCT116_hg19)"
    )
    parser.add_argument(
        "--chipseq-histone", action="store_true",
        help="Download all histone ChIP-seq marks (12 marks: HCT116_hg19_all)"
    )
    parser.add_argument(
        "--status", action="store_true", help="Check which data files are available"
    )
    parser.add_argument("--setup", action="store_true", help="Create the data directory structure")
    parser.add_argument("--output-dir", type=str, help="Custom output directory for data files")

    args = parser.parse_args()

    if args.output_dir:
        os.environ["SCRIBE_DATA_DIR"] = args.output_dir

    if args.status:
        check_data_status()
    elif args.setup:
        setup_directory_structure()
    elif args.all or args.hic or getattr(args, 'chipseq_encode', False) or getattr(args, 'chipseq_histone', False):
        # Download requested data
        if args.all or args.hic:
            download_hic_data(force=args.force)
        if args.all or getattr(args, 'chipseq_encode', False):
            download_chipseq_data(cell_types=["HCT116_hg19"], force=args.force)
        if getattr(args, 'chipseq_histone', False):
            download_chipseq_data(cell_types=["HCT116_hg19_all"], force=args.force)
    elif len(sys.argv) == 1:
        # No args - show status and help
        check_data_status()
        print("\n" + "=" * 70)
        print("To download all data automatically:")
        print("  python -m scribe.download_data --all")
        print("\nFor more options:")
        print("  python -m scribe.download_data --help")
        print("=" * 70)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
