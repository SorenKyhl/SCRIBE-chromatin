# SCRIBE Quickstart Examples

This directory contains step-by-step example scripts demonstrating the SCRIBE workflow.

## Prerequisites

1. Install SCRIBE:
   ```bash
   conda env create -f environment.yml
   conda activate scribe
   make all
   ```

## Individual Scripts

Run these in order to understand each step:

| Script | Description |
|--------|-------------|
| `01_download_data.py` | Download Hi-C and ChIP-seq data from ENCODE |
| `02_process_data.py` | Load Hi-C and ChIP-seq data, save as numpy arrays |
| `03_run_simulation.py` | Run a forward simulation with default χ parameters |
| `04_maxent_optimization.py` | Learn χ parameters via maximum entropy |
| `05_analyze_simulation.py` | Analyze and compare simulation results |
| `06_pipeline_sweep.py` | Sweep over different sequence representations |
| `07_run_from_converged.py` | Run a simulation from a shipped, pre-converged maxent result (no optimization needed) |

## Full Workflows

Self-contained scripts that run the complete workflow:

| Script | Description |
|--------|-------------|
| `full_simulation.py` | Data → Simulation (downloads data if needed) |
| `full_maxent.py` | Data → Maximum Entropy Optimization |

## Quick Start

For the fastest start, run one of the full workflow scripts:

```bash
# Run a forward simulation
python full_simulation.py

# Or learn parameters with maximum entropy
python full_maxent.py
```

Maximum entropy optimization is expensive to rerun just to get realistic chi
parameters. If you just want a simulation with maxent-optimized chis without
running the optimization yourself, use the pre-converged result shipped with
the package:

```bash
python 07_run_from_converged.py
```

## Output Files

After running the examples:

- `chipseq_sequences.npy` - Polymer sequences from ChIP-seq
- `experimental_hic.npy` - Hi-C contact map
- `output/` or `simulation_output/` - Simulation results
- `maxent_output/` - Maximum entropy optimization results
  - `plaid_chis.npy`, `daig_chis.npy` - Learned interaction parameters (plaid and diagonal)
  - `SCC.txt` - Correlation trajectory
