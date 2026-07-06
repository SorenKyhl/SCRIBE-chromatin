#!/usr/bin/env python3
"""
SCRIBE Quickstart Example 7: Run a Simulation from a Converged Maxent Result

Maximum entropy optimization (see 04_maxent_optimization.py / full_maxent.py)
is expensive to rerun just to get a simulation with realistic chi parameters.
This script shows how to spawn a fresh forward simulation directly from a
converged result shipped with the package, skipping optimization entirely.

This script demonstrates how to:
1. Load a converged config + sequences via scribe.default.load_converged()
2. Run a forward simulation with those maxent-optimized chi parameters
3. Visualize the resulting contact map

Output:
    - converged_output/: Simulation output directory
    - contactmap.png: Contact map visualization
"""

import shutil
from pathlib import Path

from scribe import default
from scribe.plot_contactmap import plot_contactmap
from scribe.scribe_sim import ScribeSim


def main():
    # config: maxent-optimized chi parameters (chis/diag_chis)
    # seqs: the (n_marks, nbeads) ChIP-seq sequence array that config was
    #       trained against -- see scribe/defaults/hct116_auxin_maxent/README.md
    config, seqs = default.load_converged("hct116_auxin_maxent")
    print(f"Loaded converged config + sequences: {seqs.shape[0]} marks, {seqs.shape[1]} beads")
    print(f"Marks: {config.get('mark_names')}")

    # Clean up previous output if it exists
    output_dir = Path("converged_output")
    if output_dir.exists():
        print(f"Removing existing output directory: {output_dir}")
        shutil.rmtree(output_dir)

    # Create simulation: converged chis, no maxent optimization needed
    sim = ScribeSim(root=str(output_dir), config=config, seqs=seqs)

    print("\nRunning simulation...")
    print("  - Equilibration: 10,000 sweeps")
    print("  - Production: 50,000 sweeps")
    sim.run_eq(equilibrium_sweeps=10000, production_sweeps=50000)

    # Visualize the resulting contact map (averaged over ensemble)
    print("\nGenerating contact map visualization...")
    plot_contactmap(str(output_dir / "production_out"))

    print("\nSimulation complete!")
    print(f"Output saved to: {output_dir}/")


if __name__ == "__main__":
    main()
