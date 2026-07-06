Quick Start
===========

Run the Example Scripts
------------------------

The fastest way to get started is to run the self-contained scripts in
``examples/quickstart/``. They download data (if needed), process it, and run
a full simulation or maximum-entropy optimization end to end.

.. code-block:: bash

   cd examples/quickstart

   # Data -> forward simulation -> contact map comparison
   python full_simulation.py

   # Data -> maximum entropy optimization of chi parameters
   python full_maxent.py

Both scripts are safe to run from any working directory (they write their
outputs relative to the current directory), so you can copy them elsewhere,
or run them against a separate results directory, e.g.:

.. code-block:: bash

   mkdir -p /path/to/results/run1 && cd /path/to/results/run1
   python /path/to/SCRIBE-chromatin/examples/quickstart/full_simulation.py

The same directory also has a numbered set of scripts (``01_download_data.py``
through ``06_pipeline_sweep.py``) that break the workflow into individual
steps -- useful if you want to inspect or modify one stage without rerunning
everything. See ``examples/quickstart/README.md`` for the full list.

The rest of this page walks through the same steps as library calls, for
when you want to script something custom instead of running the example
files directly.

Downloading Data
----------------

SCRIBE uses experimental Hi-C and ChIP-seq data for training and validation. The data files are large and hosted externally. We provide helper functions to download an example dataset for HCT116 cells from ENCODE.

.. code-block:: bash

   # Check what data is available/missing
   python -m scribe.download_data --status

   # Download all data (~36 GB: Hi-C + ChIP-seq)
   python -m scribe.download_data --all

   # Or download separately:
   python -m scribe.download_data --hic              # Hi-C only (~29 GB)
   python -m scribe.download_data --chipseq-encode    # ENCODE ChIP-seq (6 marks, ~6.8 GB)
   python -m scribe.download_data --chipseq-histone   # All histone ChIP-seq (12 marks, ~12.9 GB)

Data is stored in ``~/.scribe/data/`` by default. Set ``SCRIBE_DATA_DIR`` to use a custom location.

Loading Data by Cell Type
-------------------------

The high-level ``DataPipeline`` loads raw Hi-C and ChIP-seq files, processes them (pooling, normalization), and prepares model-ready inputs:

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
   from scribe.analysis import SimulationResult, get_SCC
   from scribe.plot_contactmap import plot_contactmap
   from scribe import default
   import numpy as np

   # Load default configuration (contains interaction parameters χ)
   config = default.config.copy()

   # Load polymer sequences (epigenetic mark occupancies from ChIP-seq)
   sequences = np.load("chipseq_sequences.npy")

   # Create simulation: sequences define bead identities, config defines χ parameters
   sim = ScribeSim(root="output", config=config, seqs=sequences)

   # Run equilibration + production to generate ensemble of 3D structures
   sim.run_eq(equilibrium_sweeps=10000, production_sweeps=50000)

   # ScribeSim is a dispatcher only -- it does not hold the resulting contact
   # map. Load the finished run (in "output/production_out") as a
   # SimulationResult to inspect and visualize it.
   result = SimulationResult("output/production_out", maxent_analysis=False)

   # Visualize the resulting contact map (averaged over ensemble)
   plot_contactmap("output/production_out")

   # Compare to experimental data
   experimental_hic = np.load("experimental_hic.npy")
   scc = get_SCC(result.hic, experimental_hic)
   print(f"SCC: {scc:.3f}")


Maximum Entropy Optimization
----------------------------

Optimize the Flory-Huggins χ interaction parameters to match experimental Hi-C contact maps. The maximum entropy framework iteratively runs simulations and updates χ until the predicted contact frequencies match the experimental data:

.. code-block:: python

   from scribe.maxent import Maxent
   from scribe.analysis import get_goals
   from scribe import default
   import numpy as np

   # Load experimental Hi-C contact map (training target)
   hic_experimental = np.load("experimental_hic.npy")

   # Load polymer sequences from ChIP-seq (defines bead identities)
   sequences = np.load("chipseq_sequences.npy")

   # Set up maximum entropy optimization
   config = default.config.copy()
   params = default.params.copy()

   # Maxent requires target observables ("goals") derived from the
   # experimental Hi-C, sequences, and config
   params["goals"] = get_goals(hic_experimental, sequences, config)

   me = Maxent(
       root="maxent_output",
       params=params,
       config=config,
       seqs=sequences,       # Input: epigenetic sequences
       gthic=hic_experimental  # Target: experimental Hi-C to match
   )

   # Run optimization: learns χ parameters that reproduce Hi-C
   me.fit()

   # To continue an interrupted run later (resumes and finishes remaining
   # iterations automatically):
   # me = Maxent.resume("maxent_output")


Running a Simulation from a Converged Result
---------------------------------------------

Running maximum entropy optimization just to get realistic chi parameters
is expensive. The package ships a converged config + sequences from a
completed ``full_maxent.py`` run, so you can spawn a forward simulation
with maxent-optimized chis directly:

.. code-block:: python

   from scribe import default
   from scribe.scribe_sim import ScribeSim

   config, seqs = default.load_converged("hct116_auxin_maxent")
   sim = ScribeSim(root="output", config=config, seqs=seqs)
   sim.run_eq(equilibrium_sweeps=10000, production_sweeps=50000)

See ``examples/quickstart/07_run_from_converged.py`` for the full script,
and ``scribe/defaults/hct116_auxin_maxent/README.md`` for provenance of
this particular result.


High-Level MaxentPipeline
-------------------------

The ``MaxentPipeline`` class is a high-level wrapper for spawning multiple maximum entropy training runs. Unlike ``Maxent``, it computes ``goals`` for you internally. Use it to systematically compare different sequence representations derived from Hi-C data (e.g., varying the number of principal components):

.. code-block:: python

   from scribe.maxent_pipeline import MaxentPipeline
   from scribe import analysis
   from scribe import default
   import functools
   import numpy as np

   # Load data and config
   experimental_hic = np.load("experimental_hic.npy")
   config = default.config.copy()
   params = default.params.copy()

   # Sweep over different numbers of principal components
   for k in range(1, 11):
       seqs_method = functools.partial(analysis.get_sequences, k=k)
       pipe = MaxentPipeline(
           name=f"pc_{k}",
           gthic=experimental_hic,
           config=config,
           params=params,
           seqs_method=seqs_method  # Derives 2k sequences from Hi-C PCA
       )
       pipe.fit()  # Runs full maximum entropy optimization
