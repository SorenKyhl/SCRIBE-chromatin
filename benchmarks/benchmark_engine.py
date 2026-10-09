#!/usr/bin/env python3
"""Reproducible throughput / profiling harness for the SCRIBE C++ engine.

This drives the engine through the shipped converged HCT116 config with a
*fixed* seed and no seed randomization, so the Monte Carlo trajectory is
bitwise-identical across builds. That makes it a correctness check as well as
a benchmark: any behavior-preserving optimization must reproduce the same
final acceptance rate and bead-move counts -- only sweeps/sec should change.

Usage:
    python benchmark_engine.py [--sweeps N] [--repeat K] [--profile]

    --sweeps N    number of MC sweeps per run (default 5000)
    --repeat K    run the benchmark K times (default 3), reporting each
    --profile     enable the engine's per-category scope timers and print an
                  aggregated breakdown instead of a throughput number
    --keep DIR    copy the first run's output directory (energy/observable
                  traces, contacts, final xyz) to DIR, so two builds can be
                  byte-compared with `diff -r` -- a much stronger check than
                  the acceptance rate alone

The engine sources are #included into src/pybind_Sim.cpp, so setuptools does
not notice edits to them: rebuild with `python setup.py build_ext --inplace
--force` before benchmarking a change.

See docs/source/performance.rst for the optimization log this harness backs.
"""
import argparse
import re
import shutil
import tempfile
import time
from collections import defaultdict
from pathlib import Path

from scribe import default
from scribe.scribe_sim import ScribeSim

SEED = 12345


def _make_sim(sweeps: int, profile: bool, root: Path) -> ScribeSim:
    config, seqs = default.load_converged("hct116_auxin_maxent")
    config["seed"] = SEED
    config["nSweeps"] = sweeps
    config["load_configuration"] = False  # start from random coil; self-contained
    if profile:
        config["profiling_on"] = True
    if root.exists():
        shutil.rmtree(root)
    # randomize_seed=False keeps the fixed seed above, for reproducibility
    return ScribeSim(root=str(root), config=config, seqs=seqs, randomize_seed=False)


def run_throughput(sweeps: int, repeat: int, workdir: Path,
                   keep: Path | None = None) -> None:
    for i in range(repeat):
        root = workdir / f"run{i}"
        sim = _make_sim(sweeps, profile=False, root=root)
        t0 = time.perf_counter()
        sim.run("out")
        elapsed = time.perf_counter() - t0
        accept = _grep(root / "out" / "log.log", r"overall acceptance rate: (\S+)")
        if keep is not None and i == 0:
            if keep.exists():
                shutil.rmtree(keep)
            shutil.copytree(root / "out", keep)
        print(f"[{i + 1}/{repeat}] {sweeps} sweeps in {elapsed:6.3f} s"
              f"  =>  {sweeps / elapsed:7.1f} sweeps/sec"
              f"   (acceptance {accept})")


def run_profile(sweeps: int, workdir: Path) -> None:
    root = workdir / "profile"
    sim = _make_sim(sweeps, profile=True, root=root)
    sim.run("out")
    log = (root / "out" / "log.log").read_text()

    totals: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    for m in re.finditer(r"^(.*?) took (\d+) microseconds", log, re.M):
        totals[m.group(1)] += int(m.group(2))
        counts[m.group(1)] += 1

    move_total = sum(v for k, v in totals.items() if k != "Initializing")
    print(f"\nPer-category scope-timer totals over {sweeps} sweeps:")
    print(f"{'category':16s} {'seconds':>10s} {'calls':>8s} {'% of moves':>11s}")
    for k in sorted(totals, key=lambda k: -totals[k]):
        pct = 100 * totals[k] / move_total if k != "Initializing" and move_total else 0.0
        print(f"{k:16s} {totals[k] / 1e6:10.3f} {counts[k]:8d} {pct:10.1f}%")


def _grep(path: Path, pattern: str) -> str:
    if not path.exists():
        return "?"
    m = re.search(pattern, path.read_text())
    return m.group(1) if m else "?"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sweeps", type=int, default=5000)
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--profile", action="store_true")
    p.add_argument("--keep", type=Path, default=None)
    args = p.parse_args()

    with tempfile.TemporaryDirectory(prefix="scribe_bench_") as tmp:
        workdir = Path(tmp)
        if args.profile:
            run_profile(args.sweeps, workdir)
        else:
            run_throughput(args.sweeps, args.repeat, workdir, args.keep)


if __name__ == "__main__":
    main()
