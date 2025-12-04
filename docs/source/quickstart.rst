Quick Start
===========

This guide walks you through the main workflows in SCRIBE.

Downloading Data
----------------

Before running examples, download the required Hi-C and ChIP-seq data:

.. code-block:: bash

   # Check what data is available/missing
   python -m scribe.download_data --status

   # Download all data (~36 GB: Hi-C + ChIP-seq)
   python -m scribe.download_data --all

   # Or download separately:
   python -m scribe.download_data --hic        # Hi-C only (~29 GB)
   python -m scribe.download_data --chipseq    # ChIP-seq only (~6.8 GB)

Data is stored in ``~/.scribe/data/`` by default. Set ``SCRIBE_DATA_DIR`` to use a custom location.

Loading Data by Cell Type
-------------------------

The high-level ``DataPipeline`` loads data by cell type, automatically finding files in ``~/.scribe/data/``:

.. code-block:: python

   from scribe.data_pipeline import DataPipeline
   import numpy as np

   # Create pipeline for HCT116 cell line data
   pipeline = DataPipeline(
       cell="HCT116_auxin",  # Cell line/condition
       chrom=2,               # Chromosome number
       nbeads=1024,           # Number of polymer beads
   )

   # Check data availability
   print(pipeline.status())

   # Load Hi-C contact map (with automatic pooling and caching)
   hic = pipeline.load_hic()

   # Load all ChIP-seq tracks as a dictionary
   sequences = pipeline.load_chipseq()
   print(f"Loaded tracks: {list(sequences.keys())}")

   # Or get as array directly
   seq_array = pipeline.load_chipseq_array()

   # Save for simulation
   np.save("chipseq_sequences.npy", seq_array)
   np.save("experimental_hic.npy", hic)


Low-Level Data Loading
----------------------

For custom file paths or fine-grained control, use ``DataLoader``:

.. code-block:: python

   from scribe.data_loader import DataLoader
   import numpy as np

   # Define genomic region explicitly
   loader = DataLoader(
       res=100000,          # 100 kbp resolution
       chrom=2,             # chromosome number
       start=0,             # start position (bp)
       end=102_400_000,     # end position (bp)
       size=1024            # number of polymer beads
   )

   # Load from specific file paths
   hic = loader.load_hic("/path/to/file.hic")
   sequences = loader.load_chipseq_from_directory("/path/to/chipseq/", method="mean")

   seq_array = np.stack(list(sequences.values()), axis=1)


Running a Simulation
--------------------

Run a forward simulation using epigenetic sequences and interaction parameters (χ) to generate an ensemble of 3D genome structures:

.. code-block:: python

   from scribe.scribe_sim import ScribeSim
   from scribe import default
   from scribe.plot_contactmap import plot_contactmap
   import numpy as np

   # Load default configuration (contains interaction parameters χ)
   config = default.config.copy()

   # Load polymer sequences (epigenetic mark occupancies from ChIP-seq)
   sequences = np.load("chipseq_sequences.npy")

   # Create simulation: sequences define bead identities, config defines χ parameters
   sim = ScribeSim(root="output", config=config, seqs=sequences)

   # Run equilibration + production to generate ensemble of 3D structures
   sim.run_eq(eq_sweeps=10000, prod_sweeps=50000)

   # Visualize the resulting contact map (averaged over ensemble)
   plot_contactmap("output")


Maximum Entropy Optimization
----------------------------

Optimize the Flory-Huggins χ interaction parameters to match experimental Hi-C contact maps. The maximum entropy framework iteratively runs simulations and updates χ until the predicted contact frequencies match the experimental data:

.. code-block:: python

   from scribe.maxent import Maxent
   from scribe import default

   # Load experimental Hi-C contact map (training target)
   hic_experimental = np.load("experimental_hic.npy")

   # Load polymer sequences from ChIP-seq (defines bead identities)
   sequences = np.load("chipseq_sequences.npy")

   # Set up maximum entropy optimization
   config = default.config.copy()
   params = default.params.copy()

   me = Maxent(
       root="maxent_output",
       params=params,
       config=config,
       seqs=sequences,       # Input: epigenetic sequences
       gthic=hic_experimental  # Target: experimental Hi-C to match
   )

   # Run optimization: learns χ parameters that reproduce Hi-C
   me.fit()


High-Level MaxentPipeline
-------------------------

The ``MaxentPipeline`` class is a high-level wrapper for spawning multiple maximum entropy training runs. Use it to systematically compare different sequence representations derived from Hi-C data (e.g., varying the number of principal components):

.. code-block:: python

   from scribe.maxent_pipeline import MaxentPipeline
   from scribe import epilib as ep
   from scribe import default
   import functools
   import numpy as np

   # Load data and config
   experimental_hic = np.load("experimental_hic.npy")
   config = default.config.copy()
   params = default.params.copy()

   # Sweep over different numbers of principal components
   for k in range(1, 11):
       seqs_method = functools.partial(ep.get_sequences, k=k)
       pipe = MaxentPipeline(
           name=f"pc_{k}",
           gthic=experimental_hic,
           config=config,
           params=params,
           seqs_method=seqs_method  # Derives 2k sequences from Hi-C PCA
       )
       pipe.fit()  # Runs full maximum entropy optimization
