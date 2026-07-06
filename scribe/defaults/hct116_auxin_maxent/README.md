# hct116_auxin_maxent

Converged maximum-entropy result for HCT116 (auxin-treated), chromosome 2,
1024 beads, 12 ChIP-seq histone marks (see `mark_names` in `config.json`).

- `config.json` -- final iteration's simulation config, i.e. `default.config`
  with `chis`/`diag_chis` replaced by the maxent-optimized values. Produced
  by `examples/quickstart/full_maxent.py` (10 iterations, final SCC ~0.72
  vs experimental Hi-C). `load_configuration` is reset to `False` so this
  config starts a fresh simulation rather than resuming from a checkpoint
  that isn't shipped with the package.
- `seqs.npy` -- the (12, 1024) ChIP-seq sequence array the optimization was
  run against. Sequence identities (not chi values) never change during
  maxent, so this is the same array as `chipseq_sequences.npy` from that run.

Load both together with `scribe.default.load_converged("hct116_auxin_maxent")`.
