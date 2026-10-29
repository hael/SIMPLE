# SIMPLE resource planner

Run in a shell with installed SIMPLE available through `PATH` and `SIMPLE_PATH`.
Python 3.9+ and NumPy are required; the launcher discovers installed interpreters.
Use `--fixture-python` if NumPy is in a separate interpreter. No build is required.

```bash
./scripts/hardware_resource_planner.py                      # Run quick calibration
./scripts/hardware_resource_planner.py calibrate            # Preview only
./scripts/hardware_resource_planner.py calibrate --full --execute
./scripts/hardware_resource_planner.py calibrate --help
```

Both presets cover motion correction, CTF estimation, picking, solve2D and solve3D.

| Preset | Workload | Measured repetitions | Thread candidates |
| --- | --- | --- | --- |
| Quick | 4K movies, 100 particles, 128-pixel boxes | 1 | 4 and 8 where available |
| Full | 4K/8K movies, 2,000 particles, 128-pixel boxes | 5 | Powers of two plus usable CPU limit |

Each configuration also has a warmup. Quick solve3D starts from a supplied volume;
full solve3D runs de novo and uses a longer schedule. Full runs can take hours.
Single-node full motion/CTF runs also test 1/2/4 workers where CPU and input counts
allow. Fixed-seed synthetic inputs support scientific checks; parallel results and
timings are not bitwise deterministic. Recommendations remain synthetic estimates.

## Run options and output

```bash
./scripts/hardware_resource_planner.py calibrate --execute \
  --workflows motion_correct ctf_estimate --threads 4 8 \
  --workdir /usr/local/cache/simple-calibration/measurements \
  --output /usr/local/cache/simple-calibration/profile.json
```

Defaults write to `output/resource_calibration/standard/` relative to the current
directory. Choose fresh output paths: existing profiles and work directories are
refused. The work filesystem must allow execution of generated worker scripts.

The console and JSON report settings, timings, quality status and failures.
`measurements/` retains commands, logs and validation evidence. Preparation and
validation are excluded from commander timings. Repeated image outputs are removed
by default; use `--no-discard-images` to keep them. Generated results should not be
committed. See `--help` for workload, timeout and disk-budget overrides.

SIMPLE does **not** read or apply the JSON automatically (`auto_apply=false`).
Only runs passing scientific checks produce recommendations. Distributed memory
estimates conservatively sum process peaks, rather than measuring concurrent RSS.

## Slurm

Log into the cluster, load SIMPLE, and specify its partition:

```bash
./scripts/hardware_resource_planner.py calibrate --full --execute \
  --slurm csbdevel --slurm-cpus 16 --slurm-memory 65536 \
  --workdir ./cluster-calibration --output ./cluster-profile.json
```

The tool submits from the current host; it does not choose an SSH login host.
CPU and memory requests apply per node (memory in MiB). Optional flags include
`--slurm-account` and `--slurm-time`; default wall limits are one hour for quick
runs and two days for full runs.

Full Slurm runs compare one and two nodes. `--node-counts 1` limits the sweep to
one node; `--node-counts 1 2 4` also tests four. Shared executable storage is
required. Native workers run inside the allocation and the master merges results.
Each trial uses one worker per node and a fixed partition count equal to the
largest requested node count. Picking stays single-node.

Recommendations require successful worker placement and scientific output checks.
Placement uses Slurm node identities, which may differ from OS hostnames. JSON
includes per-node results, speedup, efficiency and node-hours. Failed node trials
return a nonzero exit status even if other trials succeed. Bounded one/two-node
runs have been validated; full-workload and four-node scaling remain unqualified.

## Optional storage comparison

```bash
./scripts/hardware_resource_planner.py calibrate --execute \
  --io-dirs /usr/local/data /usr/local/cache
```

Defaults: 256 MiB sequential data, 128 files of 64 KiB, three repetitions per
directory. Override with `--io-size-mib`, `--io-small-files` and `--io-repeats`.
Directories must exist and have sufficient free space. Probes remove only their
own temporary files; forced termination may leave a `.simple-io-*` directory.

The table and JSON include seconds, MiB/s and variability. Writes include file and
directory `fsync`; reads include checksum verification. First reads request cache
eviction, but cold reads are not guaranteed. Payload blocks repeat after 16 MiB,
so deduplication can affect results. These probes do not measure SIMPLE workflow
speedups or automatically choose scratch storage. On Slurm they run only on the
coordinator node and do not establish whether paths are shared. Failed directories
are reported while remaining probes continue; failures produce a nonzero exit.

## Development

Keep `resource_planner/`, `memory/`, and the memory-estimator scripts beside the
entry point when copying the tool. Workflow adapters live in `workflows/`, input
generators in `fixtures/`, and Python tests in `tests/`.

```bash
PYTHONPATH=scripts python3 -m unittest discover -s scripts/resource_planner/tests
```

Set `SIMPLE_FIXTURE_PYTHON` to a NumPy-capable interpreter for fixture replay tests.
See [remaining work](../../doc/implementation_notes/planned/hardware_resource_planning_unfinished.md)
for qualification limits and future native integration.
