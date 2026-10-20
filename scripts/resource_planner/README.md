# SIMPLE resource planner

The public entry point is `scripts/hardware_resource_planner.py`. Run it in
the shell where SIMPLE is installed and available through `PATH` and
`SIMPLE_PATH`; no build directory is required. The launcher requires Python
3.9 or newer and can select an installed newer interpreter if `python3` is
older. Fixture generation requires NumPy; use `--fixture-python` when it is
installed in a separate interpreter.

From the repository root:

```bash
./scripts/hardware_resource_planner.py discover
./scripts/hardware_resource_planner.py calibrate            # Preview only
./scripts/hardware_resource_planner.py                      # Execute standard calibration
./scripts/hardware_resource_planner.py calibrate --help
```

The standard workflows are motion correction, CTF estimation, picking,
solve2D and solve3D. Particle workflows import Python-generated clean molecular
projection stacks directly (`ctf=no`); solve3D additionally prepares a serial
solve2D result. They do not run motion correction, CTF estimation or picking
as prerequisites. Inputs use fixed seeds and bundled
molecular coordinates. Parallel solver output and timings need not be bitwise
identical. Optional older analytic fixtures use `--fixture-source analytic`.

## Selecting a run

```bash
./scripts/hardware_resource_planner.py calibrate \
  --workflows motion_correct ctf_estimate --threads 4 8 \
  --workdir output/resource_calibration/desktop/measurements \
  --output output/resource_calibration/desktop/profile.json --execute
```

Default output is `output/resource_calibration/standard/`, relative to the
launch directory. `profile.json` contains workflow status and qualified
settings; `measurements/` contains plans, individual run logs, validations
and timing evidence. Existing work directories and profile files are refused
to preserve earlier runs. Choose a new run directory for another calibration.
These generated files are not source code and should not be committed.

The solve3D candidate limit defaults to 1800 seconds. Use
`--solve3d-timeout SECONDS` to override it. An explicit `--timeout` applies
to solve3D too unless its dedicated override is supplied. One full desktop
solve3D run took about 17 minutes in the observed fixture; warm-ups and
repetitions can make a full campaign take hours. Preparation and validation
are excluded from the reported commander timings.

## Slurm

Log into the desired cluster, load SIMPLE there, and run the same entry point:

```bash
./scripts/hardware_resource_planner.py calibrate \
  --workflows motion_correct ctf_estimate \
  --slurm csbdevel --slurm-cpus 16 --slurm-memory 65536 \
  --slurm-time 01:00:00 --execute
```

The tool submits from the current machine; it does not select an SSH host.
The administrator supplies the partition and, when needed, `--slurm-account`.
Discovery occurs inside a single-node allocation. Request a larger Slurm
wall-time limit for full solve3D sweeps. Multi-node scaling is not qualified.

## Package layout

- `cli.py`: parsing and command dispatch.
- `hardware.py`: host inventory and usable CPU selection.
- `calibration.py`: preparation, candidate sweeps and evidence collection.
- `execution.py`, `supervise.py`: process execution, limits and cleanup.
- `slurm.py`: allocation metadata and batch submission.
- `profiles.py`, `reporting.py`: JSON profiles and console output.
- `campaign.py`, `adapters.py`: explicit custom campaigns and command contracts.
- `workflows/`: per-workflow adapters, shared mechanics and molecular preparation.
- `fixtures/`: independently executable NumPy helpers, validation and bundled data.
- `examples/campaign.json`: custom campaign template; supply installed executable
  paths and your own fixture directory and preparation/validation scripts.
- `tests/`: Python-only checks with fake native commands.

Keep the package beside the entry point when copying the tool. The analytic
motion fixture also imports the existing `scripts/memory/` helpers; include
that directory when copying the scripts. The retired desktop wrappers are no
longer public commands. Use `calibrate --workflows ...` instead.

## Validation and current limits

```bash
PYTHONPATH=scripts python3 -m unittest discover \
  -s scripts/resource_planner/tests -p 'test_*.py'
# Optional separate NumPy interpreter for the fixture replay check:
SIMPLE_FIXTURE_PYTHON=/usr/bin/python3.6 PYTHONPATH=scripts python3 -m unittest \
  discover -s scripts/resource_planner/tests -p 'test_*.py'
```

The JSON is evidence for future SIMPLE integration: `auto_apply=false`.
SIMPLE does not currently read it automatically. Synthetic quality gates
must pass before recommendations are exported. Current observed limitations
include CTF profile export losing a temporary input, picking accuracy, and
solve2D/solve3D scientific qualification. Reorganization does not fix these.
See the living design record in
[`hardware_resource_planning_completed.md`](../../doc/implementation_notes/completed/hardware_resource_planning_completed.md).
