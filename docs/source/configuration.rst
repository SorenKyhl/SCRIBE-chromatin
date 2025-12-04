Configuration Reference
=======================

SCRIBE simulations are configured via a JSON file (``config.json``) that specifies all simulation parameters. This page provides a comprehensive reference for all available configuration options.

.. contents:: Table of Contents
   :local:
   :depth: 2

Example Configuration
---------------------

Here is a minimal configuration file to get started:

.. code-block:: json

   {
     "nbeads": 1024,
     "nSweeps": 50000,
     "seed": 82353,
     
     "bonded_on": true,
     "nonbonded_on": true,
     "plaid_on": true,
     "diagonal_on": true,
     
     "bond_type": "gaussian",
     "bond_length": 30,
     "beadvol": 520,
     "grid_size": 28.7,
     
     "phi_chromatin": 0.06,
     "boundary_type": "sphere"
   }

System Parameters
-----------------

Basic System Setup
^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``nbeads``
     - int
     - Number of beads (monomers) in the polymer chain. Determines the genomic resolution of the simulation.
   * - ``seed``
     - int
     - Random number generator seed for reproducibility.
   * - ``nSweeps``
     - int
     - Total number of Monte Carlo sweeps to perform.
   * - ``nEquilibSweeps``
     - int
     - Number of equilibration sweeps before production (deprecated, use separate equilibration run).

Physical Parameters
^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``beadvol``
     - float
     - Volume of each bead in nm³. Default: 520 nm³.
   * - ``bond_length``
     - float
     - Equilibrium bond length between adjacent beads in nm. Default: 30 nm.
   * - ``grid_size``
     - float
     - Size of each grid cell in nm. Controls the resolution of the field-theoretic interactions.
   * - ``phi_chromatin``
     - float
     - Target chromatin volume fraction (0-1). Used to calculate simulation volume.
   * - ``target_volume``
     - float
     - Alternative to ``phi_chromatin``: directly specify simulation volume in μm³.

Boundary Conditions
^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``boundary_type``
     - string
     - Simulation boundary geometry. Options: ``"sphere"``, ``"spherical"``, ``"cube"``, ``"cubic"``, ``"spheroid"``.
   * - ``aspect_ratio``
     - float
     - For spheroid boundary: ratio of polar to equatorial radius. Required when ``boundary_type`` is ``"spheroid"``.

Energy Terms
------------

Main Energy Toggles
^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``bonded_on``
     - bool
     - Enable bonded (chain connectivity) energy terms.
   * - ``nonbonded_on``
     - bool
     - Enable non-bonded (field-theoretic) energy terms.
   * - ``plaid_on``
     - bool
     - Enable epigenetic (plaid/checkerboard) interaction energy based on histone modifications.
   * - ``diagonal_on``
     - bool
     - Enable genomic-distance-dependent (diagonal) interaction energy.
   * - ``constant_chi_on``
     - bool
     - Enable uniform attractive interaction between all beads.
   * - ``boundary_attract_on``
     - bool
     - Enable attractive interactions with the simulation boundary (LAD-like).

Bond Parameters
^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``bond_type``
     - string
     - Type of bond potential. Options: ``"gaussian"`` (recommended), ``"FENE"``, ``"DSS"``, ``"SC"``.
   * - ``k_bond``
     - float
     - Bond spring constant (only for SC bond type).

Angle Parameters
^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``angles_on``
     - bool
     - Enable angular (bending) energy terms for chain stiffness.
   * - ``k_angle``
     - float
     - Angular spring constant.
   * - ``theta_0``
     - float
     - Equilibrium bond angle in degrees. Default: 180 (straight chain).

Epigenetic (Plaid) Interactions
-------------------------------

SCRIBE's key feature is modeling chromatin interactions based on epigenetic marks (histone modifications).

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``nspecies``
     - int
     - Number of epigenetic species/marks (e.g., 5 for H3K27ac, H3K4me3, H3K27me3, H3K9me3, H3K36me3).
   * - ``load_bead_types``
     - bool
     - Whether to load epigenetic mark profiles from files.
   * - ``bead_type_files``
     - list[str]
     - List of files containing epigenetic mark profiles (one file per species). Each file has ``nbeads`` lines with mark intensity values.
   * - ``chis``
     - 2D array
     - Interaction matrix (nspecies × nspecies) specifying the strength of interactions between epigenetic marks. Symmetric matrix where ``chis[i][j]`` is the interaction strength between species i and j.
   * - ``load_chipseq``
     - bool
     - Alias for ``load_bead_types``.

Example epigenetic configuration:

.. code-block:: json

   {
     "plaid_on": true,
     "nspecies": 5,
     "load_bead_types": true,
     "bead_type_files": ["pcf1.txt", "pcf2.txt", "pcf3.txt", "pcf4.txt", "pcf5.txt"],
     "chis": [
       [0.5, 0.0, -0.2, 0.0, 0.0],
       [0.0, 0.3, 0.0, 0.0, 0.0],
       [-0.2, 0.0, 0.8, 0.1, 0.0],
       [0.0, 0.0, 0.1, 0.6, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.2]
     ]
   }

Diagonal (Genomic Distance) Interactions
----------------------------------------

Distance-dependent interactions capture the decay of contact probability with genomic separation.

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``diag_chis``
     - list[float]
     - Interaction strengths for each genomic distance bin. Length determines number of bins.
   * - ``diag_start``
     - int
     - Starting genomic distance (in beads) for diagonal interactions.
   * - ``diagonal_linear``
     - bool
     - Use linear binning for diagonal interactions. If false, uses dense binning scheme.
   * - ``dense_diagonal_on``
     - bool
     - Enable variable-width binning with finer resolution near the diagonal.
   * - ``dense_diagonal_cutoff``
     - float
     - Fraction of chain length for high-resolution diagonal region (e.g., 0.0625).
   * - ``dense_diagonal_loading``
     - float
     - Fraction of bins allocated to the high-resolution region (e.g., 0.5).

Matrix-Based Interactions
^^^^^^^^^^^^^^^^^^^^^^^^^

For advanced users, interaction strengths can be specified as full matrices:

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``lmatrix_on``
     - bool
     - Enable L-matrix (plaid) from file.
   * - ``lmatrix_filename``
     - string
     - Path to L-matrix file (nbeads × nbeads).
   * - ``dmatrix_on``
     - bool
     - Enable D-matrix (diagonal) from file.
   * - ``dmatrix_filename``
     - string
     - Path to D-matrix file (nbeads × nbeads).
   * - ``umatrix_on``
     - bool
     - Enable U-matrix (combined L+D) from file.
   * - ``umatrix_filename``
     - string
     - Path to U-matrix file.

Monte Carlo Move Parameters
---------------------------

Move Types
^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``displacement_on``
     - bool
     - Enable single-bead displacement moves.
   * - ``translation_on``
     - bool
     - Enable segment translation moves (recommended: true).
   * - ``crankshaft_on``
     - bool
     - Enable crankshaft rotation moves (recommended: true).
   * - ``pivot_on``
     - bool
     - Enable pivot rotation moves (recommended: true).
   * - ``rotate_on``
     - bool
     - Enable bead orientation rotation moves (only for DSS bonds).
   * - ``gridmove_on``
     - bool
     - Enable grid origin translation (recommended: true).

Move Parameters
^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``decay_length``
     - int
     - Characteristic decay length for selecting move partners. Controls the typical size of collective moves. Default: nbeads/2.

Density Control
---------------

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``density_cap_on``
     - bool
     - Enable maximum local density constraint to prevent unphysical overlap.
   * - ``phi_solvent_max``
     - float
     - Maximum allowed solvent volume fraction per grid cell (e.g., 0.5).
   * - ``compressibility_on``
     - bool
     - Enable compressibility energy term for more realistic density fluctuations.
   * - ``kappa``
     - float
     - Compressibility coefficient (only used if ``compressibility_on`` is true).

Output Parameters
-----------------

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``dump_frequency``
     - int
     - Save configuration and contact map every N sweeps.
   * - ``dump_stats_frequency``
     - int
     - Save energy and observables every N sweeps.
   * - ``dump_density``
     - bool
     - Save density profiles to file.
   * - ``dump_observables``
     - bool
     - Save observable trajectories. Default: true.
   * - ``print_acceptance_rates``
     - bool
     - Print MC move acceptance rates during simulation.
   * - ``profiling_on``
     - bool
     - Enable detailed timing profiling (for debugging/optimization).
   * - ``print_trans``
     - bool
     - Print detailed translation move information (for debugging).

Contact Map Parameters
----------------------

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``contact_resolution``
     - int
     - Number of beads per contact map bin. Default: 1 (bead-resolution contact map).
   * - ``track_contactmap``
     - bool
     - Save contact maps at each dump (creates separate files per timepoint).
   * - ``conservative_contact_pooling``
     - bool
     - Use conservative contact counting (avoids double-counting in pooled bins).
   * - ``double_count_main_diagonal``
     - bool
     - Whether to count self-contacts (main diagonal) in contact map.
   * - ``contact_bead_skipping``
     - bool
     - Skip intermediate beads when computing contacts (for coarse-graining).
   * - ``visit_tracking``
     - bool
     - Only count one contact per pixel per configuration (reduces sampling bias).
   * - ``update_contacts_distance``
     - bool
     - Use distance-based contact criterion instead of grid-based.
   * - ``distance_cutoff``
     - float
     - Distance threshold for contacts (only if ``update_contacts_distance`` is true).

Initial Configuration
---------------------

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``load_configuration``
     - bool
     - Load initial bead positions from file instead of generating random coil.
   * - ``load_configuration_filename``
     - string
     - Path to XYZ file containing initial configuration.

Parallelization
---------------

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``parallel``
     - bool
     - Enable parallel computation of energy terms.
   * - ``set_num_threads``
     - bool
     - Manually set number of threads (only if ``parallel`` is true).
   * - ``num_threads``
     - int
     - Number of threads to use (only if ``set_num_threads`` is true).
   * - ``cell_volumes``
     - bool
     - Use cell-based volume calculation (for parallelization).

Advanced Parameters
-------------------

Pseudobeads (Diagonal Interpolation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``diag_pseudobeads_on``
     - bool
     - Enable pseudobead interpolation for smoother diagonal energy calculation.

Constant Interaction
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``constant_chi``
     - float
     - Uniform interaction strength between all bead pairs (if ``constant_chi_on`` is true).

Boundary Attraction
^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``boundary_chi``
     - float
     - Strength of boundary attraction (if ``boundary_attract_on`` is true). Useful for modeling lamina-associated domains (LADs).

Complete Example
----------------

Here is a complete configuration file for a typical SCRIBE simulation:

.. code-block:: json

   {
     "nbeads": 1024,
     "seed": 82353,
     "nSweeps": 50000,
     "dump_frequency": 50000,
     "dump_stats_frequency": 100,
     
     "bonded_on": true,
     "nonbonded_on": true,
     "plaid_on": true,
     "diagonal_on": true,
     
     "bond_type": "gaussian",
     "bond_length": 30,
     "beadvol": 520,
     "grid_size": 28.7,
     
     "phi_chromatin": 0.06,
     "phi_solvent_max": 0.5,
     "density_cap_on": true,
     "boundary_type": "sphere",
     
     "displacement_on": false,
     "translation_on": true,
     "crankshaft_on": true,
     "pivot_on": true,
     "rotate_on": false,
     "gridmove_on": true,
     
     "nspecies": 5,
     "load_bead_types": true,
     "bead_type_files": ["pcf1.txt", "pcf2.txt", "pcf3.txt", "pcf4.txt", "pcf5.txt"],
     "chis": [
       [0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0]
     ],
     
     "diag_start": 0,
     "dense_diagonal_on": true,
     "diag_pseudobeads_on": true,
     "dense_diagonal_cutoff": 0.0625,
     "dense_diagonal_loading": 0.5,
     "diagonal_linear": true,
     "diag_chis": [0, 0, 0, 0, 0, 0, 0, 0],
     
     "contact_resolution": 1,
     "track_contactmap": false,
     "conservative_contact_pooling": true,
     
     "lmatrix_on": false,
     "umatrix_on": false,
     "dmatrix_on": false,
     
     "angles_on": false,
     "compressibility_on": false,
     "parallel": false,
     "profiling_on": false,
     "print_acceptance_rates": true,
     "load_configuration": false
   }
