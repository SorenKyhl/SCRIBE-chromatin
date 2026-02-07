"""
Low-level data loading utilities for SCRIBE.

This module provides the DataLoader class for loading raw Hi-C and ChIP-seq
data from .hic and .bigWig files into numpy arrays.

For high-level data access by cell type, use the DataPipeline class from
scribe.datapipeline instead.
"""

from pathlib import Path

import hicstraw
import numpy as np
import pandas as pd
import pyBigWig

from scribe import analysis
from scribe import hic as hiclib


def get_experiment_marks(directory):
    """
    Make mapping between experimental codes and epigenetic marks in specified directory.

    Args:
        directory: Directory containing .bigWig files and metadata.tsv

    Returns:
        dict: Maps file accession codes to epigenetic mark names
    """
    directory = Path(directory)
    metadata = pd.read_csv(directory / "metadata.tsv", sep="\t")

    # Support both "Experiment target" (ENCODE format) and "Target" (simplified format)
    if "Experiment target" in metadata.columns:
        target_col = "Experiment target"
    elif "Target" in metadata.columns:
        target_col = "Target"
    else:
        raise KeyError(
            f"metadata.tsv must have 'Experiment target' or 'Target' column. Found: {list(metadata.columns)}"
        )

    marks = metadata[target_col].apply(lambda s: s.split("-")[0])
    lookup_table = dict(zip(metadata["File accession"], marks))
    return lookup_table


class DataLoader:
    """Low-level loader for Hi-C and ChIP-seq files.

    This class handles the conversion of raw data files (.hic, .bigWig) into
    numpy arrays. It manages genomic coordinates, resolution, and data cleaning.

    For high-level access by cell type, use DataPipeline instead.

    Args:
        res: Resolution of contact map in base pairs per bin
        chrom: Chromosome number or name
        start: Start position in base pairs
        end: End position in base pairs
        size: Desired number of bins along one dimension of contact map

    Example:
        >>> loader = DataLoader(res=100000, chrom=2, start=0, end=102_400_000, size=1024)
        >>> hic = loader.load_hic("/path/to/file.hic")
        >>> chipseq = loader.load_bigWig("/path/to/H3K27ac.bigWig")
    """

    def __init__(self, res, chrom, start, end, size):
        self.res = int(res)
        self.start = start
        self.end = end
        self.size = size
        self.set_chrom(chrom)

        self.bigsize = self.size
        self.dropped_inds = []

    def set_chrom(self, chrom):
        """Set the chromosome for data loading."""
        self.chrom = str(chrom)
        self.chromstr = "chr" + str(self.chrom)

    def resize(self, newsize):
        """Resize the loader to a new number of bins."""
        factor = int(newsize / self.size)
        self.res = int(self.res / factor)
        self.size = newsize
        return self

    def load_hic(self, filename, KR=True, clean=True, rescale_method="ones"):
        """Load Hi-C directly from .hic file at the loader's resolution.

        This loads the contact map at the resolution specified by self.res.
        For more accurate coarse-graining, use load_hic_pooled() instead,
        which loads at high resolution and pools down.

        Args:
            filename: Path to .hic file
            KR: Use Knight-Ruiz normalization
            clean: Remove rows/columns with zero diagonal
            rescale_method: Method to rescale contact map ("mean" or "ones")

        Returns:
            Contact map of shape (self.size, self.size)
        """
        filename = str(filename)  # hicstraw doesn't accept pathlib objects
        hic = hicstraw.HiCFile(filename)
        contactmap = analysis.load_contactmap_hicstraw(
            hic, self.res, self.chrom, self.start, self.end, KR=KR
        )

        self.bigsize, _ = np.shape(
            contactmap
        )  # used for loading the correct number of chipseq bins

        if clean:
            contactmap, self.dropped_inds = analysis.clean_contactmap(contactmap)

        # set main diag to one (on average)
        if rescale_method:
            contactmap = hiclib.normalize_hic(contactmap, rescale_method)

        return contactmap[0 : self.size, 0 : self.size]

    def load_hic_pooled(
        self,
        filename,
        nbeads: int = None,
        highres_beads: int = 20480,
        pool_fn=None,
        cache: bool = False,
        cache_dir: Path = None,
        KR=True,
        clean=True,
        rescale_method="ones",
    ):
        """Load Hi-C from .hic file at high resolution and pool down.

        This method loads the contact map at high resolution (default 20480 beads)
        and then pools down to the requested number of beads. This preserves
        contact counts more accurately than loading directly at low resolution.

        Args:
            filename: Path to .hic file
            nbeads: Target number of beads. If None, uses self.size.
            highres_beads: Number of beads to load at high resolution (default 20480).
                Must be divisible by nbeads.
            pool_fn: Pooling function to use. If None, uses hic.pool() which
                handles diagonal elements correctly to avoid double-counting.
            cache: If True, cache the high-res Hi-C to disk for faster subsequent loads.
            cache_dir: Directory for cache files. If None, uses ~/.scribe/data/cache/
            KR: Use Knight-Ruiz normalization
            clean: Remove rows/columns with zero diagonal
            rescale_method: Method to rescale contact map ("mean" or "ones")

        Returns:
            Pooled contact map of shape (nbeads, nbeads)
        """
        import copy

        if pool_fn is None:
            pool_fn = hiclib.pool

        if nbeads is None:
            nbeads = self.size

        # Set up cache directory
        if cache:
            if cache_dir is None:
                from scribe.paths import get_data_dir

                cache_dir = get_data_dir() / "cache"
            cache_dir = Path(cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)

            # Create cache key from filename and parameters
            filename_path = Path(filename)
            cache_key = (
                f"{filename_path.stem}_chr{self.chrom}_{self.start}-{self.end}_{highres_beads}"
            )
            cache_file = cache_dir / f"{cache_key}.npy"

        # Check cache
        if cache and cache_file.exists():
            gthic = np.load(cache_file)
            # Restore metadata for ChIP-seq loading
            # Note: dropped_inds are not cached, so clean=True won't work correctly
            # when loading from cache. This is a limitation.
            self.bigsize = len(gthic)
            self.dropped_inds = []
        else:
            # Calculate what resolution we need for highres_beads
            region_size = self.end - self.start
            highres_resolution = region_size // highres_beads

            # Create a copy of loader configured for high-res loading
            highres_loader = copy.deepcopy(self)
            highres_loader.res = highres_resolution
            highres_loader.size = highres_beads

            # Load at high resolution
            gthic = highres_loader.load_hic(
                filename, KR=KR, clean=clean, rescale_method=rescale_method
            )

            # Store dropped indices for ChIP-seq loading compatibility
            self.dropped_inds = highres_loader.dropped_inds
            self.bigsize = highres_loader.bigsize

            # Save to cache
            if cache:
                np.save(cache_file, gthic)

        # Pool down to target resolution
        factor = len(gthic) // nbeads
        if factor < 1:
            raise ValueError(f"highres_beads ({highres_beads}) must be >= nbeads ({nbeads})")
        if len(gthic) % nbeads != 0:
            raise ValueError(f"highres_beads ({len(gthic)}) must be divisible by nbeads ({nbeads})")

        pooled = pool_fn(gthic, factor, normalize=True)
        return pooled

    def load_bigWig(self, filename, method="mean"):
        """Load ChIP-seq signal from bigWig file.

        Can load from local or remote files.

        Args:
            filename: Path to bigWig file
            method: Aggregation method ("mean" or "max")

        Returns:
            Signal array of shape (self.size,)
        """
        assert method in ["mean", "max"]
        filename = str(filename)  # pyBigWig doesn't accept pathlib objects
        bw = pyBigWig.open(filename)
        signal = bw.stats(self.chromstr, self.start, self.end, type=method, nBins=self.bigsize)
        signal = np.delete(signal, self.dropped_inds)
        return np.array(signal[0 : self.size])

    def load_wig(self, filename, method):
        """Load ChIP-seq from .wig file format.

        Can only load from local files.

        Args:
            filename: Path to .wig file
            method: Aggregation method ("mean" or "max")

        Returns:
            Signal array
        """
        assert method in ["mean", "max"]
        df = pd.read_csv(filename, sep="\t", names=["start", "end", "value"], skiprows=1)
        chip = analysis.bin_chipseq(df, self.res, method=method)
        chip = np.nan_to_num(chip)
        return chip

    def load_chipseq_from_files(self, filenames, method):
        """Load ChIP-seq from a list of files.

        Args:
            filenames: List of paths to ChIP-seq files (.bigWig or .wig)
            method: Aggregation method ("mean" or "max")

        Returns:
            dict: Maps mark names to signal arrays
        """
        seqs = {}
        for file in filenames:
            extension = file.suffix

            if extension == ".bigWig":
                name = file.name.split("_")[0]
                seqs[name] = self.load_bigWig(file, method)
            elif extension == ".wig":
                name = str(file).split("_")[-2]
                seqs[name] = self.load_wig(file, method)
            else:
                raise ValueError(f"Unsupported file extension: {extension}")

        return seqs

    def load_chipseq_from_directory(self, directory, method):
        """Load all ChIP-seq files from a directory.

        Expects a metadata.tsv file in the directory mapping file accessions
        to experiment targets.

        Args:
            directory: Directory containing .bigWig files and metadata.tsv
            method: Aggregation method ("mean" or "max")

        Returns:
            dict: Maps mark names to signal arrays
        """
        seqs = {}
        directory = Path(directory)
        filenames = sorted(directory.glob("*.bigWig"))
        lookup_table = get_experiment_marks(directory)

        for file in filenames:
            extension = file.suffix

            if extension == ".bigWig":
                name = lookup_table[file.stem]
                seqs[name] = self.load_bigWig(file, method)
            elif extension == ".wig":
                name = lookup_table[file.stem]
                seqs[name] = self.load_wig(file, method)
            else:
                raise ValueError(f"Unsupported file extension: {extension}")
            print(f"loaded {name}")

        return seqs
