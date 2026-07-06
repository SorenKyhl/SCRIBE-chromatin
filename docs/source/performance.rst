Engine Performance
==================

This page is a living record of performance work on the C++ Monte Carlo engine
(``src/``): how to measure it, what has been optimized and why, and what
remains. Update it whenever you profile or optimize the engine.

Benchmarking and profiling
---------------------------

Use ``benchmarks/benchmark_engine.py``. It drives the engine through the shipped
converged HCT116 config (:func:`scribe.default.load_converged`) with a **fixed
seed** and no seed randomization, so the Monte Carlo trajectory is
bitwise-identical across builds. This makes it both a benchmark and a
correctness check: any behavior-preserving change must reproduce the same final
acceptance rate and bead-move counts -- only sweeps/sec should move.

.. code-block:: bash

   # Throughput (sweeps/sec), 3 repeats
   python benchmarks/benchmark_engine.py --sweeps 5000 --repeat 3

   # Per-category time breakdown (enables the engine's scope timers)
   python benchmarks/benchmark_engine.py --sweeps 400 --profile

The ``--profile`` mode reports true per-move-category times (translation,
crankshaft, pivot, grid move, ...). For a function-level view, attach a
sampling profiler to a running benchmark, e.g. on macOS::

   python benchmarks/benchmark_engine.py --sweeps 40000 &
   sample $! 6

.. note::
   The engine reads ``config.json`` from the current working directory, so the
   benchmark writes into a temporary directory and runs there. The reported
   throughput excludes Python/setup overhead only loosely; compare relative
   numbers from the same machine.

Key finding: the hot path is allocation-bound
----------------------------------------------

The single most important observation: for the canonical converged config, the
Monte Carlo hot path was dominated by heap allocation, not floating-point
compute. Rebuilding with ``-O3 -march=native`` over the default ``-O2`` gave no
measurable speedup, while a sampling profile showed ``malloc``/``free``/
``memcpy`` dominating the innermost energy routine. Consequently the wins came
from removing per-move and per-cell allocations, and the build flags were left
unchanged (``-march=native`` also hurts binary portability for no benefit here).

Optimization log
----------------

All changes below are behavior-preserving: on the fixed-seed benchmark the
trajectory stayed bitwise-identical (overall acceptance 82.2284%). Throughput
figures are 5000-sweep runs at the default ``-O2`` on the same machine, so they
are only meaningful relative to each other.

.. list-table::
   :header-rows: 1
   :widths: 26 44 14 16

   * - Change
     - What / why
     - Result
     - Commit
   * - Reuse ``getDiagEnergy`` scratch + pass ``diag_chis`` by ref
     - The diagonal energy (hottest inner loop, run per cell per energy eval)
       heap-allocated a ``std::vector<int>`` every call and took its 28-element
       ``diag_chis`` **by value**. Reuse a ``thread_local`` buffer; pass by
       const-ref. Also removed dead ``std::chrono`` timing in
       ``getNonBondedEnergy`` and reused the crank/pivot ``old_positions``
       buffers.
     - ~842 -> ~1127
     - ``0fd6d17``
   * - Flat buffers + generation stamp for flagged cells
     - Each move built a ``std::unordered_set<Cell*>`` (flagged cells) and
       ``std::unordered_map`` (bead swaps), which node-allocate per element.
       Replace with reusable ``std::vector`` buffers on ``Sim``; deduplicate
       flagged cells with a per-``Cell`` ``uint64`` generation stamp
       (``flagCell``/``beginFlagging``) instead of a hash set. Energy functions
       now take ``const std::vector<Cell*>&``; ``Grid`` gained
       ``active_cells_vec``.
     - ~1127 -> ~1250
     - ``dbf6275``
   * - Grid move checks only the density cap
     - ``MCmove_grid`` is not a Metropolis move: it re-grids and accepts unless
       a cell overflows the density cap. It computed the full nonbonded energy
       (density cap + plaid + diagonal) over every active cell every sweep, but
       only the density-cap term can reach the rejection threshold. Compute
       just ``densityCapEnergy``.
     - ~1250 -> ~1500
     - ``0960bb7``
   * - Skip inert orientation work; drop exceptions in crank/pivot
     - Bead orientation (``Bead::u``) only affects the energy for DSS bonds;
       for gaussian/harmonic bonds it is inert, yet crankshaft/pivot saved and
       quaternion-rotated it per bead. Gate on an ``orientation_active`` flag.
       Also replaced exception-based rejection (``throw "rejected"``) with a
       plain flag.
     - ~1500 -> ~1650
     - ``b1d7e8a``

Cumulatively this is roughly a **1.9x** speedup (~842 -> ~1650 sweeps/sec),
with an identical simulation trajectory.

The scope-timer profiler was also fixed as part of this work: the per-category
``Timer`` objects in ``Sim::MC()`` shared one block scope (explicit destructor
calls were commented out), so each reported cumulative time-to-end-of-sweep.
Each category now has its own block scope, plus a ``gridmove`` timer, so
``--profile`` reports true per-category splits.

Current breakdown
-----------------

After the changes above, a typical ``--profile`` run of the converged config
(translation and crankshaft have equal move counts, ``n_trans == n_crank ==
decay_length``; ``n_pivot`` is 10x fewer) looks like:

=============  ===========
Category       % of moves
=============  ===========
translating    ~60%
cranking       ~28%
gridmove       ~9%
pivoting       ~3%
=============  ===========

Translation dominates because a translation displaces every bead in its segment
by the full step, so more beads cross grid-cell boundaries (more flagged cells,
more energy evaluation) than a crankshaft rotation of the same-length segment,
whose beads sit close to the rotation axis and barely move. Both are now
limited by the shared per-cell energy computation
(``Cell::getEnergy`` + ``Cell::getDiagEnergy``).

Remaining opportunities
-----------------------

Not yet done; roughly in increasing order of risk/effort:

- **``meshBeads`` resets every cell.** The grid move (and its revert) reset all
  ``(L+1)^3`` cells and re-insert every bead each sweep. Resetting only occupied
  or active cells would cut most of the remaining grid-move cost.
- **Grid-move frequency.** The grid move runs every sweep to suppress
  discretization artifacts. Reducing its frequency is a modeling decision, not a
  pure optimization -- it changes results -- so it is left to the user.
- **``Cell::getEnergy`` division hoisting.** The plaid inner loop
  (``ntypes`` x ``ntypes``) does a division per iteration; hoisting it changes
  floating-point rounding (not bitwise-identical), so it needs statistical
  validation.
- **``getDiagEnergy`` early exit.** Sorting a cell's bead indices would let the
  pairwise loop break once the genomic separation exceeds ``diag_cutoff``,
  pruning distant pairs. Marginal for sparse cells; helps dense ones.
- **Incremental / delta energy (large refactor).** Each move recomputes full
  old and new energies over all flagged cells from scratch, though only a few
  beads changed. Computing just the energy delta of the moved beads would be a
  substantial win but is a significant rework of the cell-based energy model.
