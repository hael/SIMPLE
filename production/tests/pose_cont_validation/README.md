# Continuous 3-D pose-refinement Oracle validation

`run_oracle_validation.sh` is the single Phase 10 runtime entry point. Run it
only after the authorized Oracle compilation and environment setup. It runs all
retained focused cases, the pose and neutral mother suites, the PCG mother
suite, and the `pcg_recon` and `pcg_priors` regressions.

The runner changes to `~/Projects` before every SIMPLE executable invocation.
It creates one `continuous_3D_pose_validation_<timestamp>` package directly
below that directory. The package contains independent command and exit-status
records, executable hashes, raw logs, machine-readable test output, source
snapshots and hashes, generated matrix copies, environment data, aggregate
analysis, and a checksum manifest. An individual failure does not stop later
independent tests.

The PCG mother suite has one frozen known failure. Phase 10 accepts status 1
only when `halfset_fsc` reproduces the exact Phase 2 disposition: best lambda
10 for both halves, PCG raw L2 `0.6612755/0.6608928`, conventional-gridding raw
L2 `0.6518738/0.6388823`, and the same scale-sensitive assertion. Any different
result makes the package fail.

After the workflow exports `SIMPLE_PATH` and updates `PATH`, run:

```bash
cd "$HOME/Projects"
"$HOME/Projects/hael_SIMPLE-rsync-test/production/tests/pose_cont_validation/run_oracle_validation.sh"
```

The entry point does not compile, normalize source, synchronize files, invoke
`pcg_pose_polish`, or call a removed pose-polisher strategy.
