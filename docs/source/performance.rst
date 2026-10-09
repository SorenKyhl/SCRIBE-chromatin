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

   # The engine sources are #included into src/pybind_Sim.cpp, so setuptools
   # does not see edits to them. Force a rebuild before every measurement.
   python setup.py build_ext --inplace --force

   # Throughput (sweeps/sec), 3 repeats
   python benchmarks/benchmark_engine.py --sweeps 5000 --repeat 3

   # Per-category time breakdown (enables the engine's scope timers)
   python benchmarks/benchmark_engine.py --sweeps 400 --profile

Matching acceptance rates is a weak check. To show that a change is
behavior-preserving, keep the outputs of a baseline build and of the changed
build and byte-compare them. ``log.log`` contains wall-clock times, so exclude
it:

.. code-block:: bash

   python benchmarks/benchmark_engine.py --sweeps 3000 --repeat 1 --keep /tmp/before
   # ... apply the change, rebuild with --force ...
   python benchmarks/benchmark_engine.py --sweeps 3000 --repeat 1 --keep /tmp/after
   diff -r -x log.log /tmp/before /tmp/after   # no output == identical

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
trajectory stayed bitwise-identical (overall acceptance 82.2284%). The last two
rows change floating-point rounding, so identity is not guaranteed for them;
see *Rounding-changing optimizations* below for how they were validated. Throughput
figures are 5000-sweep runs at the default ``-O2`` on the same machine, so they
are only meaningful relative to each other.

.. note::
   **Re-verified 2026-10-08 (Apple M3, -O2).** The first four rows were
   recorded in July, and some of their figures did not reproduce. The
   pre-optimization engine (``d4421f7``) runs at **~610** sweeps/sec, not
   ~842. The engine after ``b1d7e8a`` runs at **~1570**, not ~1650. That makes
   the July work a ~2.6x speedup, not ~1.9x. The bitwise claim does hold. With
   a byte comparison of every output file (energy, observables, diagonal
   observables, contacts, final xyz) at 3000 sweeps, ``d4421f7`` and
   ``2abf125`` are identical. The no-gain-from-``-O3 -march=native`` finding
   also still holds, even now that the hot path is compute-bound. The rows
   from ``Cell::contains`` onward were measured in the same session (3000-sweep
   runs), each with the byte-comparison check.

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
   * - ``Cell::contains`` as a flat vector
     - Each cell's ``std::unordered_set<Bead*>`` node-allocated on every insert
       and freed on every erase or clear: in ``meshBeads`` every sweep, and in
       every local move's ``moveIn``/``moveOut``. Use a ``std::vector`` with
       swap-and-pop erase (cells hold ~3 beads). This was previously thought
       to need statistical validation because it changes iteration order. It
       doesn't: every consumer of ``contains`` on the default path is
       order-independent (integer pair counts per diagonal bin, ``typenums``
       for plaid, integer contact counts), and the outputs are byte-identical.
     - ~1570 -> ~1870
     - ``d2db0ac``
   * - Per-cell energy cache
     - Each move computes the flagged cells' nonbonded energy before and after.
       A cell's plaid and diagonal energies depend only on its contents and
       volume, so ``Cell`` caches them and ``moveIn``/``moveOut``/``reset``/
       volume updates invalidate the cache. After an accepted move, a cell's
       cached "new" energy is the next move's "old" energy. The cached value is
       exactly what recomputation would give, so results are byte-identical.
       This assumes ``chis``/``diag_chis`` are fixed for the life of a ``Sim``.
     - ~1870 -> ~2170
     - ``e31dd74``
   * - Fixed-size ``unit_vec``
     - Took and returned a dynamic ``Eigen::MatrixXd``, so it heap-allocated
       on every translation, rotation, pivot and grid move.
     - ~+1%
     - ``c97efa4``
   * - Skip angles when ``k_angle == 0``
     - The converged config has ``angles_on`` with zero stiffness, so each move
       evaluated four angle energies that were exactly zero.
     - ~2190 -> ~2280
     - ``36f7c82``
   * - One pass for plaid observables
     - ``saveObservables`` made 78 passes (one per species pair) over ~1300
       active cells every 10 sweeps. Now it makes one pass, with the same
       per-pair summation order.
     - ~2280 -> ~2470
     - ``c488f8c``
   * - Plaid energy from an incrementally maintained ``S n``
     - ``sum_{i<=j} chi_ij n_i n_j`` equals ``n.(S n)/2`` with ``S`` the
       symmetric matrix from the upper triangle of ``chis`` (diagonal
       doubled). Each bead carries ``S d`` (``Bead::chi_d``) and each cell keeps
       ``S n`` (``Cell::chi_n``) up to date in ``moveIn``/``moveOut``, so
       ``getEnergy`` is a 12-term dot product with one division instead of 78
       terms with a division each. Rounding-changing.
     - ~2500 -> ~3200
     - ``48cc9b7``
   * - Diagonal energy from a per-separation weight table
     - The energy is linear in the per-bin pair counts, so ``getDiagEnergy``
       sums ``diag_chis[bin(d)] * nbonds(d)`` per pair from a table indexed by
       genomic separation instead of binning into ``diag_phis`` and dotting.
       ``diag_phis`` are filled only for the observables
       (``Cell::updateDiagPhis``, same integer counts as before).
       Rounding-changing.
     - ~3200 -> ~4080
     - ``fae31e6``

Cumulatively the engine is now about **6.7x** faster than before this work
(~610 -> ~4080 sweeps/sec on an M3). The byte-identical changes from
``Cell::contains`` through the observables pass give ~1.57x over ``2abf125``
(~1570 -> ~2470, interleaved runs); the two rounding-changing energy changes
give another ~1.63x (``bb32a69`` measured at ~2500 on 2026-10-09, -> ~4080).

Rounding-changing optimizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The plaid ``S n`` and diagonal weight-table changes compute the same energies
in a different floating-point order, so they could in principle flip an
accept/reject decision whose random number lands within rounding error
(~1e-13 relative) of the acceptance probability. A rough estimate puts that
at about once per 1e12 moves, and the checks agree. Against
``bb32a69``, all of the following were **byte-identical** in every output
file (energy, observables, diagonal observables, extra, contacts, xyz):

- the benchmark config, seed 12345, 5000 sweeps, and seeds 1-4 at 20000
  sweeps each;
- the benchmark config with ``gridmove_on: false`` (20000 sweeps), which
  never re-meshes and so lets ``Cell::chi_n`` accumulate incremental rounding
  for the whole run, with ``diagonal_on: false``, and with ``smatrix_on``;
- the ideal chain (no plaid/diagonal), three contact-pooling modes;
- the research repo's regression tier 0 (stored-configuration energies, zero
  difference at printed precision vs. reference) and tier 1 (2000-sweep
  trajectory, identical to the baseline engine, acceptance 82.4308%).

Per the research tier definitions, tier 2 (ensemble statistics) is only needed
when tier 1 diverges, so it was not run. If a future config does diverge,
that is the expected consequence of rounding, not by itself a regression.


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
translating    ~52%
cranking       ~33%
gridmove       ~12%
pivoting       ~3%
=============  ===========

(Before the ``Cell::contains`` change, gridmove was ~14%, not the ~9% this
table used to show. It fell to ~5% after that change and rose back to ~12%
with the ``S n`` change, because the per-sweep re-mesh's ``moveIn`` now also
updates ``chi_n``, while the moves themselves got cheaper.) The scope timers cover only the moves. The
``dump_stats_frequency`` output (``saveEnergy``, ``updateContacts``,
``saveObservables``) is not timed and is roughly another 10-15% of wall time.

Translation dominates because a translation displaces every bead in its segment
by the full step, so more beads cross grid-cell boundaries (more flagged cells,
more energy evaluation) than a crankshaft rotation of the same-length segment,
whose beads sit close to the rotation axis and barely move. A sampling
profile shows no measurable ``malloc``/``free`` left in the hot path. After
the energy changes, the per-cell energy computation is no longer dominant:
the move routines' own bookkeeping (``MCmove_translate`` /
``MCmove_crankshaft`` bodies: cell lookups, flagging, bead updates) is about
40% of samples, ``Grid::diagEnergy`` (diagonal pair loop, inlined) about 15%,
``getNonBondedEnergy`` (plaid dot product and density cap, inlined) about 11%,
``meshBeads`` ~6%, and ``exp``/``log`` ~8%.

Remaining opportunities
-----------------------

The bitwise-preserving opportunities are largely used up, and so are the
cheap rounding-changing ones in the energy. Rounding-changing work can't be
proven correct by byte comparison, but in practice it stays byte-identical
at printed precision (see *Rounding-changing optimizations*), which is the
first thing to check. The engine is now limited by move bookkeeping rather
than energy arithmetic.

Bitwise-preserving, small:

- **Keep output files open.** Each stats dump ``fopen``/``fclose``-es four or
  five files. ``open``/``close``/``write`` syscalls are ~3% of samples. Keeping
  the ``FILE*`` open with an ``fflush`` per dump would remove most of that.

Tried and rejected (measured, not worth it):

- **Keep the energy cache across the grid move's re-mesh.** ``meshBeads``
  invalidates every cell each sweep, which is why only about half of "old"
  energy lookups hit the cache. Re-validating cells whose membership is
  unchanged and whose rebuilt ``typenums`` are bit-identical made the engine
  ~4% *slower*. Most occupied cells are touched within a sweep, so their
  incrementally updated ``typenums`` rarely bit-match the rebuilt ones. An
  incremental re-mesh (move only beads whose cell index changed) would avoid
  this, but it changes ``typenums`` rounding, so it needs statistical
  validation.
- **Skip** ``exp`` **for downhill moves** (``dU >= 0`` always accepts, since
  ``uniform()`` is in [0, 1)). Exact, but no measurable change.
- **Resetting only occupied cells in** ``meshBeads``: ~1% at most (see the git
  history of this page).
- **Incremental diagonal pair sum.** Keeping each cell's diagonal pair sum up
  to date in ``moveIn``/``moveOut`` (O(k) per bead) instead of recomputing it
  per energy evaluation (O(k^2)) made the engine ~5% *slower*: with ~3 beads
  per cell the pair loop is short, and the incremental update also runs on
  rejected moves and on every bead of the per-sweep re-mesh.
- **Per-separation bin/count tables alone** (bitwise-identical version of the
  diagonal change, keeping the ``diag_phis`` dot product): ~+0.5%. The cost
  was the per-call zeroing and dotting of ``diag_phis``, not the binning.

Larger, unexplored:

- **Move bookkeeping.** The translate/crankshaft bodies (cell lookup per
  bead, flagging, swap records, ``moveIn``/``moveOut`` and undo on rejection)
  are now the largest cost. Profile them at line level before changing
  anything.
- **Incremental re-mesh.** ``meshBeads`` rebuilds every cell each sweep;
  moving only beads whose cell index changed would cut it (and keep more
  energy-cache entries valid), at the cost of more incremental rounding in
  ``typenums``/``chi_n``.
- **getDiagEnergy early exit** (sorting a cell's bead ids to stop at
  ``diag_cutoff``) does nothing for the converged config, whose
  ``diag_cutoff`` equals ``nbeads``. Only worth it for configs with a short
  cutoff and dense cells.

Modeling decisions, not optimizations:

- **Grid-move frequency.** The grid move runs every sweep to suppress
  discretization artifacts. Reducing its frequency changes results.

Not a performance issue, but noticed while verifying the grid-move change:
``MCmove_grid`` rejects only when the summed density-cap energy reaches
``9999999999``, while one over-full cell contributes ``99999999 * phi``
(~5e7). A grid move is therefore rejected only if more than ~100 cells
overflow at once. For the converged config this never matters (no cell came
near the cap in 3000 grid moves), but the two constants look inconsistent.
