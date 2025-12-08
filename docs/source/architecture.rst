Software Architecture
=====================

SCRIBE is organized into a C++ simulation core with Python interfaces for ease of use.
The Python API follows a layered architecture separating **simulation dispatch**, **result inspection**, and **pipeline orchestration**.

Design Philosophy
-----------------

The API is structured around three core concepts:

- **ScribeSim** (Dispatcher): Handles the simulation lifecycle—configuration, execution, and checkpointing. Use this to *run* simulations.
- **SimulationResult** (Inspector): Loads completed simulation outputs for analysis and metric computation. Use this to *analyze* results.
- **Pipelines** (Orchestrators): Combine multiple steps into reusable workflows. ``analysis_pipeline`` orchestrates analysis; ``maxent_pipeline`` orchestrates optimization runs.

This separation allows you to:

1. Run simulations and analyze results independently
2. Swap analysis methods without changing simulation code
3. Build custom workflows by combining low-level components

Directory Structure
-------------------

.. code-block:: text

   SCRIBE-chromatin/
   ├── src/                      # C++ simulation engine (TICG core)
   ├── scribe/                   # Python interface and analysis tools
   │   ├── scribe_sim.py         # Simulation dispatcher (ScribeSim)
   │   ├── analysis.py           # Result inspection (SimulationResult, metrics)
   │   ├── analysis_pipeline.py  # Analysis workflows (sim_analysis, compare_analysis)
   │   ├── maxent.py             # Maximum entropy optimizer (Maxent)
   │   ├── maxent_pipeline.py    # MaxEnt workflow automation (MaxentPipeline)
   │   ├── data_pipeline.py      # High-level data loading by cell type
   │   ├── data_loader.py        # Low-level file loading (.hic, .bigWig)
   │   └── default.py            # Default configurations
   ├── examples/                 # Tutorial notebooks and scripts
   ├── defaults/                 # Default configuration files
   └── tests/                    # Unit and integration tests


Module Hierarchy
----------------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Module
     - Role
     - Description
   * - ``scribe_engine``
     - Engine
     - C++ extension (pybind11 wrapper) providing core TICG simulation routines
   * - ``scribe_sim.ScribeSim``
     - Dispatcher
     - Set up, configure, and run simulations. Main interface for launching jobs.
   * - ``analysis.SimulationResult``
     - Inspector
     - Load and analyze completed simulation outputs. Computes metrics and generates plots.
   * - ``analysis_pipeline``
     - Orchestrator
     - High-level analysis workflows: ``sim_analysis()``, ``compare_analysis()``
   * - ``maxent.Maxent``
     - Optimizer
     - Core maximum entropy optimization. Iteratively updates χ parameters.
   * - ``maxent_pipeline.MaxentPipeline``
     - Orchestrator
     - Maximum entropy workflow. Spawns and manages multiple training runs.
   * - ``data_pipeline.DataPipeline``
     - Loader
     - High-level data loading by cell type from ``~/.scribe/data/``.
   * - ``data_loader.DataLoader``
     - Loader
     - Low-level file loading from .hic and .bigWig files.


Core Components
---------------

C++ Simulation Engine (``src/``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simulation engine implements the Theoretically Informed Coarse Grain (TICG) method:

- **O(N) complexity** via field-theoretic energy functionals
- **Monte Carlo sampling** of polymer configurations
- **Contact map computation** from ensemble of structures

The engine is compiled into a Python extension module (``scribe_engine``) using pybind11.

Python API (``scribe/``)
------------------------

**ScribeSim** - Simulation Dispatcher

- Manages simulation setup and configuration
- Handles file I/O for coordinates, contact maps, and parameters
- Runs equilibration and production phases
- Use this to *launch* simulations

**SimulationResult** - Result Inspector

- Loads completed simulation outputs (contacts, energy, observables)
- Computes metrics: SCC, RMSE, diagonal statistics
- Generates analysis plots: contact maps, energy trajectories, χ matrices
- Use this to *analyze* simulation results after they complete

**analysis_pipeline** - Analysis Orchestrator

- High-level functions combining multiple analysis steps
- ``sim_analysis(sim)``: Full analysis of a single simulation
- ``compare_analysis(sim)``: Compare simulation to experimental Hi-C
- Saves publication-quality figures

**Maxent** - Maximum Entropy Optimizer

- Implements the iterative optimization loop
- Updates χ parameters based on gradient of contact frequency differences
- Tracks convergence metrics (loss, SCC, parameter updates)

**MaxentPipeline** - MaxEnt Orchestrator

- Coordinates data loading, simulation, and optimization
- Enables parameter sweeps (e.g., over number of principal components)
- Manages output directory structure

**DataPipeline** - High-level Data Loader

- Loads data by cell type name (e.g., "HCT116_auxin")
- Automatically finds files in ``~/.scribe/data/``
- Provides caching for processed Hi-C data
- Methods: ``load_hic()``, ``load_chipseq()``, ``status()``

**DataLoader** - Low-level File Loader

- Reads ENCODE bigWig files for ChIP-seq data
- Reads .hic files for Hi-C contact maps
- Handles genomic coordinate systems and resolution matching
- Use when you need custom file paths or fine-grained control


Configuration
-------------

SCRIBE uses JSON configuration files. See the :doc:`configuration` page for a complete reference of all parameters.

Default configurations are provided in ``defaults/`` and can be loaded via:

.. code-block:: python

   from scribe import default
   config = default.config.copy()
   params = default.params.copy()
