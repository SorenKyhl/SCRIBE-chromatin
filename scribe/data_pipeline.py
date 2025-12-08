"""
High-level data pipeline for SCRIBE.

This module provides the DataPipeline class for loading Hi-C and ChIP-seq
data by cell type name, automatically resolving file paths from the
~/.scribe/data/ database.

For low-level file loading, use DataLoader from scribe.data_loader.
"""

from pathlib import Path

import numpy as np

from scribe import hic as hiclib
from scribe.data_loader import DataLoader
from scribe.paths import get_data_dir


class DataPipeline:
    """
    High-level data pipeline for loading Hi-C and ChIP-seq by cell type.

    This class provides a semantic interface to the SCRIBE data directory.
    Instead of specifying file paths, you specify cell types and the pipeline
    automatically resolves paths from ``~/.scribe/data/``.

    Args:
        cell (str): Cell type name (e.g., ``HCT116_auxin``, ``GM12878``)
        chrom (int or str): Chromosome number or name
        nbeads (int): Number of polymer beads (bins) for the simulation
        start (int, optional): Start position in base pairs (default: 0)
        end (int, optional): End position in base pairs (default: inferred from nbeads and resolution)
        assembly (str, optional): Genome assembly for ChIP-seq (default: ``hg19``)
        highres_beads (int, optional): Number of beads for high-resolution loading (default: 20480)
        cache (bool, optional): Whether to cache processed data (default: True)

    Example:
        >>> pipeline = DataPipeline(cell="HCT116_auxin", chrom=2, nbeads=1024)
        >>> hic = pipeline.load_hic()
        >>> chipseq = pipeline.load_chipseq()

    Data Directory Structure::

        ~/.scribe/data/
        ├── hic/
        │   └── {cell}/           # e.g., HCT116_auxin/
        │       └── *.hic
        ├── chipseq/
        │   └── {cell}_{assembly}/ # e.g., HCT116_hg19/
        │       ├── *.bigWig
        │       └── metadata.tsv
        └── cache/                 # Cached processed data
    """

    # Default genomic parameters
    DEFAULT_RESOLUTION = 100_000  # 100kb target resolution
    DEFAULT_START = 0
    # High-resolution for loading - must be a resolution available in .hic files
    # Common resolutions: 5000, 10000, 25000, 50000, 100000
    DEFAULT_HIGHRES_RESOLUTION = 5000  # 5kb - highest common resolution
    # Fraction of extra region to load to account for bins lost during cleaning
    # Typical cleaning removes ~4% of bins due to zeros on diagonal
    CLEANING_BUFFER = 0.10  # Load 10% extra to be safe

    def __init__(
        self,
        cell: str,
        chrom: int = 2,
        nbeads: int = 1024,
        start: int = None,
        end: int = None,
        assembly: str = "hg19",
        highres_resolution: int = None,
        cache: bool = True,
    ):
        self.cell = cell
        self.chrom = chrom
        self.nbeads = nbeads
        self.assembly = assembly
        self.cache = cache

        # High-resolution loading - use a fixed resolution that's available in .hic files
        self.highres_resolution = highres_resolution if highres_resolution is not None else self.DEFAULT_HIGHRES_RESOLUTION
        
        # Calculate genomic coordinates
        # We load extra region to account for bins dropped during cleaning
        self.start = start if start is not None else self.DEFAULT_START
        
        if end is not None:
            self.end = end
        else:
            # Load extra region to account for cleaning losses
            # Target: nbeads at 100kb = nbeads * 100kb base pairs
            # But we load (1 + buffer) * that to have enough after cleaning
            target_bp = self.nbeads * self.DEFAULT_RESOLUTION
            raw_end = self.start + int(np.ceil(target_bp * (1 + self.CLEANING_BUFFER)))
            # Align end to bufsize boundary (bufsize = 1000 * highres_resolution in epilib)
            bufsize = 1000 * self.highres_resolution
            self.end = ((raw_end // bufsize) + 1) * bufsize  # Round up to next buffer boundary
        
        # Calculate how many high-res beads we'll load
        self.highres_beads = (self.end - self.start) // self.highres_resolution

        # Data directory
        self.data_dir = get_data_dir()
        self.cache_dir = self.data_dir / "cache"

        # Internal loader (created lazily)
        self._loader: DataLoader | None = None

        # Cached data
        self._hic = None
        self._chipseq = None
        self._dropped_inds = []
        self._cleaned_size = self.highres_beads  # Size after cleaning

    @property
    def loader(self) -> DataLoader:
        """Get or create the internal DataLoader."""
        if self._loader is None:
            self._loader = DataLoader(
                res=self.highres_resolution,
                chrom=self.chrom,
                start=self.start,
                end=self.end,
                size=self.highres_beads,
            )
        return self._loader

    def _get_hic_path(self) -> Path | None:
        """Get the path to the Hi-C file for this cell type."""
        hic_dir = self.data_dir / "hic" / self.cell
        if not hic_dir.exists():
            return None
        hic_files = list(hic_dir.glob("*.hic"))
        if not hic_files:
            return None
        return hic_files[0]

    def _get_chipseq_dir(self) -> Path | None:
        """Get the ChIP-seq directory for this cell type."""
        # Try with assembly suffix first
        cell_base = self.cell.split("_")[0]  # HCT116_auxin -> HCT116
        chipseq_dir = self.data_dir / "chipseq" / f"{cell_base}_{self.assembly}"
        if chipseq_dir.exists():
            return chipseq_dir

        # Try exact cell name
        chipseq_dir = self.data_dir / "chipseq" / self.cell
        if chipseq_dir.exists():
            return chipseq_dir

        return None

    def _get_cache_path(self, data_type: str) -> Path:
        """Get the cache file path for a given data type."""
        cache_key = f"{self.cell}_chr{self.chrom}_{self.start}-{self.end}_{self.highres_beads}"
        return self.cache_dir / f"{cache_key}_{data_type}.npy"

    def load_hic(self, pool_fn=None, force_reload: bool = False) -> np.ndarray:
        """Load Hi-C contact map for this cell type.

        Loads at high resolution and pools down to nbeads. Results are cached
        to disk for faster subsequent loads.

        Args:
            pool_fn: Pooling function (default: hic.pool with diagonal handling)
            force_reload: If True, reload from source even if cached

        Returns:
            Contact map of shape (nbeads, nbeads)

        Raises:
            FileNotFoundError: If Hi-C data not found for this cell type
        """
        if self._hic is not None and not force_reload:
            return self._hic

        if pool_fn is None:
            pool_fn = hiclib.pool

        # Check cache first
        cache_path = self._get_cache_path("hic_highres")
        if self.cache and cache_path.exists() and not force_reload:
            gthic = np.load(cache_path)
            self._cleaned_size = len(gthic)
        else:
            # Find Hi-C file
            hic_path = self._get_hic_path()
            if hic_path is None:
                raise FileNotFoundError(
                    f"Hi-C data not found for cell type '{self.cell}'. "
                    f"Expected location: {self.data_dir / 'hic' / self.cell}/"
                )

            # Load at high resolution (cleaning happens inside loader)
            gthic = self.loader.load_hic(hic_path)
            self._dropped_inds = self.loader.dropped_inds
            self._cleaned_size = len(gthic)

            # Cache the high-res cleaned data
            if self.cache:
                self.cache_dir.mkdir(parents=True, exist_ok=True)
                np.save(cache_path, gthic)

        # Pool down to target nbeads
        # After cleaning, we need at least nbeads worth of high-res bins
        # The pooling factor should give us exactly nbeads
        if self._cleaned_size < self.nbeads:
            raise ValueError(
                f"Hi-C data has {self._cleaned_size} bins after cleaning, but nbeads={self.nbeads}. "
                f"Try reducing nbeads or specifying a larger genomic region with 'end' parameter."
            )
        
        # Find the largest pooling factor that gives at least nbeads
        factor = self._cleaned_size // self.nbeads
        usable_size = factor * self.nbeads
        
        # Trim to usable size (symmetric matrix)
        gthic = gthic[:usable_size, :usable_size]

        self._hic = pool_fn(gthic, factor, normalize=True)
        return self._hic

    def load_chipseq(
        self,
        method: str = "mean",
        marks: list = None,
        apply_processing: bool = True,
        force_reload: bool = False,
    ) -> dict:
        """Load ChIP-seq data for this cell type.

        Args:
            method: Aggregation method for bigWig files ("mean" or "max")
            marks: List of marks to load (default: all available)
            apply_processing: If True, apply Smooth -> Normalize -> Sigmoid pipeline
            force_reload: If True, reload from source even if cached

        Returns:
            dict: Maps mark names to signal arrays of shape (nbeads,)

        Raises:
            FileNotFoundError: If ChIP-seq data not found for this cell type
        """
        if self._chipseq is not None and not force_reload:
            return self._chipseq

        # Find ChIP-seq directory
        chipseq_dir = self._get_chipseq_dir()
        if chipseq_dir is None:
            cell_base = self.cell.split("_")[0]
            raise FileNotFoundError(
                f"ChIP-seq data not found for cell type '{self.cell}'. "
                f"Expected location: {self.data_dir / 'chipseq' / f'{cell_base}_{self.assembly}'}/"
            )

        # Need to load Hi-C first to get dropped indices
        if self._hic is None:
            try:
                self.load_hic()
            except FileNotFoundError:
                # If no Hi-C, load without dropped indices
                self._dropped_inds = []
                self._cleaned_size = self.highres_beads

        # Update loader with dropped indices
        self.loader.dropped_inds = self._dropped_inds
        self.loader.bigsize = self.highres_beads

        # Load ChIP-seq
        seqs_raw = self.loader.load_chipseq_from_directory(chipseq_dir, method)

        # Filter marks if specified
        if marks is not None:
            seqs_raw = {k: v for k, v in seqs_raw.items() if k in marks}

        # Pool sequences to target nbeads (same factor as Hi-C)
        factor = self._cleaned_size // self.nbeads if self._cleaned_size > self.nbeads else 1
        usable_size = factor * self.nbeads
        
        seqs_pooled = {}
        for name, seq in seqs_raw.items():
            # Handle None values
            seq = np.nan_to_num(seq, nan=0.0)
            # Trim to usable size and pool
            seq = seq[:usable_size]
            if factor > 1:
                seq = hiclib.pool_seqs(seq, factor)
            seqs_pooled[name] = seq

        # Apply processing pipeline
        if apply_processing:
            from scribe.chipseq_pipeline import ChipseqPipeline, Normalize, Sigmoid, Smooth

            pipeline = ChipseqPipeline([Smooth(), Normalize(), Sigmoid()])
            seqs_processed = {}
            for name, seq in seqs_pooled.items():
                seqs_processed[name] = pipeline.fit(seq)
            self._chipseq = seqs_processed
        else:
            self._chipseq = seqs_pooled

        return self._chipseq

    def load_chipseq_array(self, **kwargs) -> np.ndarray:
        """Load ChIP-seq as a numpy array instead of dict.

        Args:
            **kwargs: Arguments passed to load_chipseq()

        Returns:
            Array of shape (n_marks, nbeads)
        """
        seqs_dict = self.load_chipseq(**kwargs)
        return np.array(list(seqs_dict.values()))

    def load_seqs_from_hic(self, k: int = 6, smooth: bool = True) -> np.ndarray:
        """Derive polymer sequences from Hi-C via PCA.

        This extracts principal components from the Hi-C contact map to use
        as polymer sequences. This is useful for benchmarking or when
        ChIP-seq data is not available.

        Args:
            k: Number of principal components (will return 2k sequences)
            smooth: Whether to Gaussian smooth the Hi-C before PCA

        Returns:
            Array of shape (2k, nbeads)
        """
        from scribe import analysis

        # Load Hi-C (already at target resolution and trimmed to nbeads)
        hic = self.load_hic()

        if smooth:
            hic = hiclib.smooth_hic(hic)

        # Get sequences via PCA
        seqs = analysis.get_sequences(hic, k, randomized=True, correct_PCA=True)

        return seqs

    def status(self) -> dict:
        """Check data availability for this cell type.

        Returns:
            dict with 'hic' and 'chipseq' availability status
        """
        hic_path = self._get_hic_path()
        chipseq_dir = self._get_chipseq_dir()

        return {
            "cell": self.cell,
            "chrom": self.chrom,
            "nbeads": self.nbeads,
            "hic": {
                "available": hic_path is not None,
                "path": str(hic_path) if hic_path else None,
            },
            "chipseq": {
                "available": chipseq_dir is not None,
                "path": str(chipseq_dir) if chipseq_dir else None,
                "n_files": len(list(chipseq_dir.glob("*.bigWig"))) if chipseq_dir else 0,
            },
        }

    def __repr__(self) -> str:
        return (
            f"DataPipeline(cell='{self.cell}', chrom={self.chrom}, "
            f"nbeads={self.nbeads}, highres_beads={self.highres_beads}, "
            f"region={self.start}-{self.end})"
        )


# SyntheticDataPipeline for testing and examples without real data
class SyntheticDataPipeline(DataPipeline):
    """
    Data pipeline that returns synthetic (all-zero) Hi-C and ChIP-seq data.
    Useful for testing installation, configuration, and workflows without requiring downloads.
    """
    def __init__(self, nbeads=128, nspecies=10, **kwargs):
        self._nspecies = nspecies
        super().__init__(cell="SYNTHETIC", chrom=1, nbeads=nbeads, **kwargs)

    def load_hic(self):
        import numpy as np
        return np.zeros((self.nbeads, self.nbeads))

    def load_chipseq(self):
        """Return synthetic ChIP-seq data as 2D array (nspecies, nbeads)."""
        import numpy as np
        return np.zeros((self._nspecies, self.nbeads))


# Convenience function for quick data loading
def load_data(
    cell: str,
    chrom: int = 2,
    nbeads: int = 1024,
    **kwargs,
) -> tuple:
    """Convenience function to load Hi-C and ChIP-seq data.

    Args:
        cell: Cell type name
        chrom: Chromosome number
        nbeads: Number of polymer beads
        **kwargs: Additional arguments passed to DataPipeline

    Returns:
        tuple: (hic_array, chipseq_dict)
    """
    pipeline = DataPipeline(cell=cell, chrom=chrom, nbeads=nbeads, **kwargs)
    hic = pipeline.load_hic()
    chipseq = pipeline.load_chipseq()
    return hic, chipseq
