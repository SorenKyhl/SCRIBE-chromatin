Software Architecture
=====================

SCRIBE is organized into a C++ simulation core with Python interfaces for ease of use.

Directory Structure
-------------------

.. code-block:: text

   SCRIBE-chromatin/
   ├── src/                 # C++ simulation engine (TICG core)
   ├── scribe/               # Python interface and analysis tools
   │   ├── scribe_sim.py         # High-level simulation interface
   │   ├── maxent.py        # Maximum entropy optimizer
   │   ├── pipeline.py      # End-to-end workflow automation
   │   ├── datapipeline.py  # High-level data loading by cell type
   │   ├── dataloader.py    # Low-level file loading (.hic, .bigWig)
   │   ├── analysis.py      # Analysis and visualization
   │   └── default.py       # Default configurations
   ├── examples/            # Tutorial notebooks and scripts
   ├── defaults/            # Default configuration files
   └── scripts/             # Analysis and batch processing scripts


Module Hierarchy
----------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Module
     - Level
     - Description
   * - ``scribe_engine``
     - Low
     - C++ extension (pybind11 wrapper) providing core simulation routines
   * - ``scribe_sim.ScribeSim``
     - High
     - Simulation setup, execution, and I/O. Main interface for running simulations.
   * - ``maxent.Maxent``
     - Low
     - Core maximum entropy optimization. Iteratively updates χ parameters.
   * - ``pipeline.Pipeline``
     - High
     - ChIP-seq processing + optimization workflow. Spawns multiple training runs.
   * - ``datapipeline.DataPipeline``
     - High
     - High-level data loading by cell type from ``~/.scribe/data/``.
   * - ``dataloader.DataLoader``
     - Low
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

**scribe_sim** - High-level simulation interface

- Manages simulation setup and configuration
- Handles file I/O for coordinates, contact maps, and parameters
- Provides plotting and analysis methods

**Maxent** - Maximum entropy optimizer

- Implements the iterative optimization loop
- Updates χ parameters based on gradient of contact frequency differences
- Tracks convergence metrics (loss, SCC, parameter updates)

**Pipeline** - Workflow automation

- Coordinates data loading, simulation, and optimization
- Enables parameter sweeps (e.g., over number of principal components)
- Manages output directory structure

**DataPipeline** - High-level data loading

- Loads data by cell type name (e.g., "HCT116_auxin")
- Automatically finds files in ``~/.scribe/data/``
- Provides caching for processed Hi-C data
- Methods: ``load_hic()``, ``load_chipseq()``, ``status()``

**DataLoader** - Low-level file loading

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
