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

   from scribe.analysis import SimulationResult, SCC
   from scipy.stats import pearsonr
   import numpy as np

   # Load a completed simulation's outputs
   result = SimulationResult("output/production_out")

   # Access simulation data
   print(f"Contact map shape: {result.hic.shape}")
   print(f"Energy trajectory: {len(result.energy)} steps")
   print(f"Number of sequences: {result.k}")

   # Compute metrics against experimental data
   experimental_hic = np.load("experimental_hic.npy")
   scc = SCC(result.hic, experimental_hic)
   pearson_r, _ = pearsonr(result.hic.flatten(), experimental_hic.flatten())
   print(f"SCC: {scc:.3f}, Pearson r: {pearson_r:.3f}")

   # Generate plots
   result.plot_contactmap()
   result.plot_energy()

Analysis Pipelines
------------------

Use the high-level pipeline functions for comprehensive analysis workflows:

.. code-block:: python

   from scribe.scribe_sim import ScribeSim
   from scribe.analysis_pipeline import sim_analysis, compare_analysis
   from scribe.analysis import SCC
   import numpy as np

   # Load simulation via dispatcher (knows output structure)
   sim = ScribeSim.from_directory("output")

   # Full analysis: energy, contact map, observables, χ parameters
   sim_analysis(sim)

   # Compare to experimental Hi-C (ground truth)
   experimental_hic = np.load("experimental_hic.npy")
   sim.gthic = experimental_hic
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

   # Load a completed maxent run
   me = Maxent(root="maxent_output", load=True)

   # Plot optimization convergence (loss and parameter updates)
   me.plot_convergence()  # Saves loss.png and param_convergence.png

   # Visualize learned χ parameters over iterations
   me.plot_plaid_chis(legend=True)  # Track χ_IJ evolution
   me.plot_diag_chis()               # Track diagonal parameters

   # Load final χ matrix and SCC trajectory
   final_chis = np.load("maxent_output/chis.npy")[-1]  # Final χ parameters
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
