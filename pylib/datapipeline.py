"""
High-level data pipeline for SCRIBE.

This module provides the DataPipeline class for loading Hi-C and ChIP-seq
data by cell type name, automatically resolving file paths from the
~/.scribe/data/ database.

For low-level file loading, use DataLoader from pylib.dataloader.
"""

from pathlib import Path

import numpy as np

from pylib import hic as hiclib
from pylib.dataloader import DataLoader
from pylib.paths import get_data_dir


class DataPipeline:
    """High-level data pipeline for loading Hi-C and ChIP-seq by cell type.

    This class provides a semantic interface to the SCRIBE data directory.
    Instead of specifying file paths, you specify cell types and the pipeline
    automatically resolves paths from ~/.scribe/data/.

    Args:
        cell: Cell type name (e.g., "HCT116_auxin", "GM12878")
        chrom: Chromosome number or name
        nbeads: Number of polymer beads (bins) for the simulation
        start: Start position in base pairs (default: 0)
        end: End position in base pairs (default: inferred from nbeads and resolution)
        assembly: Genome assembly for ChIP-seq (default: "hg19")
        highres_beads: Number of beads for high-resolution loading (default: 20480)
        cache: Whether to cache processed data (default: True)

    Example:
        >>> pipeline = DataPipeline(cell="HCT116_auxin", chrom=2, nbeads=1024)
        >>> hic = pipeline.load_hic()
        >>> chipseq = pipeline.load_chipseq()

    Data Directory Structure:
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
    DEFAULT_RESOLUTION = 100_000  # 100kb
    DEFAULT_START = 0
    DEFAULT_HIGHRES_BEADS = 20480

    def __init__(
        self,
        cell: str,
        chrom: int = 2,
        nbeads: int = 1024,
        start: int = None,
        end: int = None,
        assembly: str = "hg19",
        highres_beads: int = None,
        cache: bool = True,
    ):
        self.cell = cell
        self.chrom = chrom
        self.nbeads = nbeads
        self.assembly = assembly
        self.cache = cache
        self.highres_beads = highres_beads or self.DEFAULT_HIGHRES_BEADS

        # Calculate genomic coordinates
        self.start = start if start is not None else self.DEFAULT_START
        if end is not None:
            self.end = end
        else:
            # Infer end from nbeads and default resolution
            self.end = self.start + (self.highres_beads * self.DEFAULT_RESOLUTION)

        # Resolution for high-res loading
        self.highres_resolution = (self.end - self.start) // self.highres_beads

        # Data directory
        self.data_dir = get_data_dir()
        self.cache_dir = self.data_dir / "cache"

        # Internal loader (created lazily)
        self._loader: DataLoader | None = None

        # Cached data
        self._hic = None
        self._chipseq = None
        self._dropped_inds = []
        self._bigsize = self.highres_beads

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
            self._bigsize = len(gthic)
        else:
            # Find Hi-C file
            hic_path = self._get_hic_path()
            if hic_path is None:
                raise FileNotFoundError(
                    f"Hi-C data not found for cell type '{self.cell}'. "
                    f"Expected location: {self.data_dir / 'hic' / self.cell}/"
                )

            # Load at high resolution
            gthic = self.loader.load_hic(hic_path)
            self._dropped_inds = self.loader.dropped_inds
            self._bigsize = self.loader.bigsize

            # Cache the high-res data
            if self.cache:
                self.cache_dir.mkdir(parents=True, exist_ok=True)
                np.save(cache_path, gthic)

        # Pool down to target nbeads
        factor = len(gthic) // self.nbeads
        if factor < 1:
            raise ValueError(
                f"highres_beads ({self.highres_beads}) must be >= nbeads ({self.nbeads})"
            )
        if len(gthic) % self.nbeads != 0:
            raise ValueError(
                f"highres_beads ({len(gthic)}) must be divisible by nbeads ({self.nbeads})"
            )

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
                self._bigsize = self.highres_beads

        # Update loader with dropped indices
        self.loader.dropped_inds = self._dropped_inds
        self.loader.bigsize = self._bigsize

        # Load ChIP-seq
        seqs_raw = self.loader.load_chipseq_from_directory(chipseq_dir, method)

        # Filter marks if specified
        if marks is not None:
            seqs_raw = {k: v for k, v in seqs_raw.items() if k in marks}

        # Pool sequences to target nbeads
        factor = self._bigsize // self.nbeads if self._bigsize > self.nbeads else 1
        seqs_pooled = {}
        for name, seq in seqs_raw.items():
            if factor > 1:
                seq = hiclib.pool_seqs(seq, factor)
            # Handle None values
            seq = np.nan_to_num(seq, nan=0.0)
            seqs_pooled[name] = seq

        # Apply processing pipeline
        if apply_processing:
            from pylib.chipseqPipeline import ChipseqPipeline, Normalize, Sigmoid, Smooth

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
        from pylib import epilib

        # Load Hi-C at high resolution
        hic_path = self._get_hic_path()
        if hic_path is None:
            raise FileNotFoundError(f"Hi-C data not found for cell type '{self.cell}'")

        gthic = self.loader.load_hic(hic_path)

        if smooth:
            gthic = hiclib.smooth_hic(gthic)

        # Get sequences via PCA
        seqs = epilib.get_sequences(gthic, k, randomized=True, correct_PCA=True)

        # Pool to target nbeads
        factor = seqs.shape[1] // self.nbeads
        if factor > 1:
            seqs = hiclib.pool_seqs(seqs, factor)

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
            f"nbeads={self.nbeads}, assembly='{self.assembly}')"
        )


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
