Installation
============

Requirements
------------

- Python 3.8+
- C++17 compiler (gcc 7+ or clang 5+)
- Conda (recommended) or pip

Quick Install (Recommended)
---------------------------

The easiest way to install SCRIBE is using conda:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/SorenKyhl/SCRIBE-chromatin.git
   cd SCRIBE-chromatin

   # Create and activate conda environment with all dependencies
   conda env create -f environment.yml
   conda activate scribe

   # Build C++ simulation engine and install Python API
   make all

This command does two things:

1. **Compiles the C++ simulation engine** (``src/``) into a Python extension module (``pyticg``) using pybind11
2. **Installs the Python API** (``scribe/``) which provides high-level interfaces for running simulations and maximum entropy optimization

Partial Builds
--------------

To reinstall just the Python API (after the C++ engine has been built once):

.. code-block:: bash

   make install

To rebuild only the C++ engine:

.. code-block:: bash

   make build

Key Dependencies
----------------

The conda environment (``environment.yml``) includes:

- ``numpy``, ``scipy`` - Numerical computing
- ``matplotlib`` - Visualization
- ``scikit-learn`` - Dimensionality reduction for sequences
- ``pyBigWig`` - Reading ChIP-seq bigWig files
- ``hic-straw`` - Reading Hi-C .hic files
- ``pybind11`` - C++/Python bindings
- ``pytest`` - Unit testing

Data Setup
----------

SCRIBE uses experimental Hi-C and ChIP-seq data for training and validation. 
These data files are large and hosted externally.

**Configuring the data directory:**

- Set the ``SCRIBE_DATA_DIR`` environment variable to specify a custom location
- Default location: ``~/.scribe/data/`` (no configuration needed)

.. code-block:: bash

   # Optional: Set custom data directory
   export SCRIBE_DATA_DIR=/path/to/your/data

**Quick Setup (Recommended):**

.. code-block:: bash

   # Check what data is available/missing
   python -m scribe.download_data --status

   # Download all data (~36 GB: Hi-C + ChIP-seq)
   python -m scribe.download_data --all

   # Or download separately:
   python -m scribe.download_data --download   # Hi-C only (~29 GB)
   python -m scribe.download_data --chipseq    # ChIP-seq only (~6.8 GB)

**Manual Download:**

If you prefer to download manually or need additional datasets:

- **Hi-C data**: Download from `GEO GSE104333 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104333>`_
  
  - Place ``.hic`` files in ``$SCRIBE_DATA_DIR/hic/``
  
- **ChIP-seq data**: Download from `ENCODE <https://www.encodeproject.org/>`_
  
  - Place ``.bigWig`` files in ``$SCRIBE_DATA_DIR/chipseq/``

The expected directory structure:

.. code-block:: text

   $SCRIBE_DATA_DIR/
   ├── hic/
   │   ├── HCT116_DMSO.hic
   │   ├── HCT116_auxin.hic
   │   └── ...
   └── chipseq/
       ├── metadata.tsv
       ├── ENCFF014WPW.bigWig
       └── ...
