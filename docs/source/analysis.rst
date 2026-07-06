Analysis
========

SCRIBE provides tools for analyzing simulation results and comparing predictions to experimental data.

The analysis API is organized into two modules:

- **analysis.py**: Low-level result inspection via ``SimulationResult`` class, plus metrics (``SCC``, ``get_diagonal``, etc.)
- **analysis_pipeline.py**: High-level workflows (``sim_analysis``, ``compare_analysis``) that orchestrate multiple analysis steps

Inspecting Simulation Results
-----------------------------

Use ``SimulationResult`` to load and inspect completed simulation outputs:

.. code-block:: python

   from scribe.analysis import SimulationResult, get_SCC
   from scipy.stats import pearsonr
   import numpy as np

   # Load a completed simulation's outputs. "production_out" is the
   # directory with the aggregated contacts.txt/energy.traj/etc -- pass
   # maxent_analysis=False if this run has no associated maxent goals.
   result = SimulationResult("output/production_out", maxent_analysis=False)

   # Access simulation data
   print(f"Contact map shape: {result.hic.shape}")
   print(f"Energy trajectory: {len(result.energy)} steps")
   print(f"Number of sequences: {result.k}")

   # Compute metrics against experimental data
   experimental_hic = np.load("experimental_hic.npy")
   scc = get_SCC(result.hic, experimental_hic)
   pearson_r, _ = pearsonr(result.hic.flatten(), experimental_hic.flatten())
   print(f"SCC: {scc:.3f}, Pearson r: {pearson_r:.3f}")

   # Generate plots
   result.plot_contactmap()
   result.plot_energy()

Analysis Pipelines
------------------

Use the high-level pipeline functions for comprehensive analysis workflows:

.. code-block:: python

   from scribe.analysis import SimulationResult
   from scribe.analysis_pipeline import sim_analysis, compare_analysis
   import numpy as np

   # sim_analysis/compare_analysis operate on a SimulationResult, not the
   # ScribeSim dispatcher used to launch the run -- ScribeSim doesn't
   # retain the finished contact map.
   sim = SimulationResult("output/production_out", maxent_analysis=False)

   # Set ground truth before calling sim_analysis: it plots the simulated
   # diagonal against gthic's even when not doing a full comparison.
   experimental_hic = np.load("experimental_hic.npy")
   sim.gthic = experimental_hic

   # Full analysis: energy, contact map, observables, χ parameters
   sim_analysis(sim)

   # Compare to experimental Hi-C (ground truth)
   compare_analysis(sim)  # Generates comparison plots (scatter, triangle, difference)

Output Files
^^^^^^^^^^^^

``sim_analysis`` generates:

- ``y.npy`` - Predicted contact map
- ``energy.png`` - Energy convergence plot
- ``chis.png`` - χ parameter matrix visualization
- Contact map plots

``compare_analysis`` generates:

- ``tri.png``, ``tri_log.png`` - Triangle comparison plots
- ``scatter.png`` - Scatter plot of predicted vs experimental contacts
- ``diff.png`` - Difference map


Maximum Entropy Analysis
------------------------

Analyze convergence and learned parameters from a completed maximum entropy optimization:

.. code-block:: python

   from scribe.maxent import Maxent
   import numpy as np
   import matplotlib.pyplot as plt

   # Load a completed (or partially completed) maxent run from disk
   me = Maxent.from_directory("maxent_output")

   # Or, to resume a partially completed run and finish remaining iterations:
   # me = Maxent.resume("maxent_output")

   # Learned χ parameters over iterations are tracked directly on the
   # Maxent object, and also saved to disk each iteration
   print(f"Plaid chis shape: {me.track_plaid_chis.shape}")
   print(f"Diagonal chis shape: {me.track_diag_chis.shape}")

   plaid_chis = np.load("maxent_output/plaid_chis.npy")   # (iterations, n_plaid_params)
   diag_chis = np.load("maxent_output/daig_chis.npy")     # (iterations, n_diag_params)
   final_plaid_chis = plaid_chis[-1]
   final_diag_chis = diag_chis[-1]

   # SCC trajectory across iterations
   scc_trajectory = np.loadtxt("maxent_output/SCC.txt")

   # Plot training progress
   plt.figure()
   plt.plot(scc_trajectory)
   plt.xlabel("Iteration")
   plt.ylabel("SCC (Stratum-adjusted Correlation)")
   plt.title("Maximum Entropy Training Progress")
   plt.savefig("training_progress.png")

Convergence Metrics
^^^^^^^^^^^^^^^^^^^

The maximum entropy optimization tracks several metrics:

- **Loss** - The objective function being minimized (difference between predicted and experimental contact frequencies)
- **SCC** - Stratum-adjusted Correlation Coefficient, a Hi-C-specific similarity metric
- **RMSE** - Root Mean Square Error between contact maps
- **Parameter convergence** - L2 norm of χ parameter updates between iterations
