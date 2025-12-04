SCRIBE Documentation
====================

.. image:: ../images/logo_scribe.png
   :alt: SCRIBE Logo
   :align: center
   :width: 500px

**SCRIBE** (Sequence-to-Chromatin structure via epigenetic Representation, Investigation, Benchmarking, and Editing) is a computational framework for **reading** the epigenetic code by predicting 3D chromatin structure from histone modifications, and **writing** novel sequences for in silico epigenetic engineering.

By combining coarse-grained polymer physics with maximum entropy optimization, SCRIBE enables:

- **Direct structure prediction** from ChIP-seq histone modification data
- **Quantitative investigation** of the histone code hypothesis
- **In silico epigenetic engineering** through computational knock-out/knock-in experiments
- **Benchmarking** of sequence-to-structure prediction models

.. image:: ../images/overview.png
   :alt: SCRIBE Overview
   :align: center
   :width: 800px

Overview
--------

The 3D organization of chromatin plays a critical role in gene regulation. SCRIBE models chromatin as a block copolymer where each monomer's identity is explicitly defined by patterns of histone post-translational modifications (PTMs) from experimental ChIP-seq data. This explicit representation of epigenetic marks is what enables both **reading** and **writing**:

- **Reading**: By learning the effective interaction parameters (Flory-Huggins χ) between epigenetic marks that best reproduce experimental Hi-C contact maps, SCRIBE reveals how the epigenetic code shapes 3D genome organization.
- **Writing**: Because epigenetic sequences are explicit inputs, users can design novel or modified sequences to predict how changes in histone modifications would alter chromatin structure, enabling in silico epigenetic engineering.

**Key features:**

- Nucleosome-resolution simulations of entire chromosomes
- O(N) complexity via Theoretically Informed Coarse Grain (TICG) field-theoretic methods
- No chromatin state-calling required: works directly with continuous ChIP-seq signals

.. note::
   This project is under active development.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   configuration
   analysis
   architecture

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   modules

Examples
--------

See the ``examples/`` directory for detailed tutorials:

- **single-simulation/** - Basic simulation setup and execution
- **chipseq_maxent/** - Training on ChIP-seq data
- **sweep_pcs/** - Parameter sweep over principal components
- **snippets/** - Useful code snippets for common tasks

Citation
--------

If you use SCRIBE in your research, please cite:

.. code-block:: bibtex

   @phdthesis{kyhl2025investigating,
     title={Investigating the epigenetic code through data-driven chromosome structure modeling},
     author={Kyhl, Soren},
     year={2025},
     school={University of Chicago},
     url={https://knowledge.uchicago.edu/record/16570}
   }

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
