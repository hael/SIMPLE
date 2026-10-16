# Simulated Validation Workflows

`simple_test_exec` provides built-in molecular workflows that generate their
own simulation data, run production processing stages, and check results
against known truth. Use these workflows to validate a processing chain.

For simulation commands that operate on your own inputs, see
[reproject](reproject.md), [simulate_movie](simulate_movie.md), and
[simulate_particles](simulate_particles.md).

## Choose a validation workflow

| Goal | Command | Starting point and result |
| --- | --- | --- |
| Process simulated movies through initial 3D reconstruction | `simple_test_exec test=simulated_workflow suite=6vxx` | Embedded molecular model; movies, a SIMPLE project, 2D classes, and a validated 3D map |
| Process synthetic micrographs through 2D classification | `simple_test_exec test=mini_stream suite=6vxx` | Embedded molecular model; aligned average micrographs, CTF estimates, extracted particles, and validated 2D classes |

## Before running

Use a SIMPLE installation with the environment configured as described in
[Installation](../installation.md). Keep `simple_test_exec`, `simple_exec`,
and `simple_private_exec` on `PATH`, because workflow stages can launch
worker processes.

Run the examples in a dedicated working directory with the local queue:

```shell
mkdir simulation_runs
cd simulation_runs
export SIMPLE_QSYS=local
```

The molecular coordinates are embedded in SIMPLE; the built-in examples need
no coordinate download or experimental dataset. They use fixed seeds and
fixture settings.

## Built-in molecular workflows

### Movies through initial 3D reconstruction

List the accepted suites, then run either molecular system:

```shell
simple_test_exec test=simulated_workflow suite=list
simple_test_exec test=simulated_workflow suite=6vxx
simple_test_exec test=simulated_workflow suite=1jxy
```

| Suite | Embedded molecule | Symmetry in this workflow |
| --- | --- | --- |
| `6vxx` | SARS-CoV-2 spike, 6VXX | `c3` |
| `1jxy` | Beta-galactosidase, 1JYX | `d2` |

The beta-galactosidase suite is named `1jxy` in the command line and output
names, while the embedded coordinate model uses the PDB identifier **1JYX**.

Each suite performs these stages:

1. Generate a centered 192 × 192 × 192 truth volume at 1.3 Å/pixel.
2. Generate 100 spiral projections and distribute them among ten movies,
   with ten projections per movie. Each movie has sixteen 1536 × 1536 frames.
3. Create a SIMPLE project, import the movies, correct motion, and estimate CTF.
4. Pick particles and extract them into 192-pixel boxes.
5. Run `solve2D`, then `solve3D`, and compare the reconstructed map to the truth.

The default picker is `picker=segdiam`. To exercise reference picking, run:

```shell
simple_test_exec test=simulated_workflow suite=6vxx picker=new
```

`picker=new` uses the generated reprojections as picking references and also
checks the generated reference dimensions, sampling, and physical size.
Both picker choices are accepted for either molecular suite.

The `1jxy` suite additionally runs `solve3D` from the same `solve2D` project
with `nparts=2`, using two local worker processes. This second run has its own
execution directory and final-map validation. The suite fixes `nthr=1`
internally for its processing stages.

### Outputs and validation

The run retains its data in `test_simulated_workflow_6vxx/` or
`test_simulated_workflow_1jxy/` beneath the directory where you launched it.
The truth volume is under `0_pdb2mrc/`, and projections are under
`1_reproject/`. The movie-generation directory contains numbered files such
as `simulate_movie_001.mrc`, `optimal_movie_average_001.mrc`, and
`simulate_movie_params_001.txt`. The subsequent project and processing
directories contain the extracted particles, class averages, and final maps.

**Rerunning the same suite removes and recreates its
`test_simulated_workflow_<suite>/` directory.** Copy results you want to retain
before rerunning, including when switching picker choices.

A successful run reports `PASS: simulated_workflow ...` with docking
correlations, whole-volume correlation, and the truth-FSC resolution. The
final-map checks require:

- A 192-pixel cubic map with 1.3 Å/pixel sampling, within 0.01 Å/pixel.
- Whole-volume Pearson correlation of at least 0.80 after docking against
  the truth in both possible hands and selecting the better fit.
- Masked FSC = 0.143 resolution against the truth of 40 Å or better.

The distributed `1jxy` result must pass the same checks. These are regression
thresholds for the synthetic fixture. The reported truth FSC compares the
reconstruction with the generating model; it is distinct from a half-map FSC.

### Synthetic micrographs through mini-stream

For a workflow that starts from micrographs and ends with 2D classes:

```shell
simple_test_exec test=mini_stream suite=list
simple_test_exec test=mini_stream suite=6vxx
simple_test_exec test=mini_stream suite=1jxy
```

Omitting `suite=` runs both molecular suites. Each suite generates six
micrographs from eight projections per micrograph, using eight movie frames
at 1.3 Å/pixel. The `6vxx` images are 1024 × 1024 pixels; the `1jxy` images
are 1280 × 1280 pixels. Projection generation uses `c3` for `6vxx` and `c1`
for this mini-stream `1jxy` fixture.

The simulator's aligned averages become `synthetic_micrograph_001.mrc`, etc.
The workflow saves matching `synthetic_micrograph_truth_001.txt` files,
removes the temporary movie stacks, and writes `mini_stream_micrographs.txt`
for input to `mini_stream`.

Results remain in a `test_mini_stream_<suite>_<process-id>/` directory,
including the `mini_stream.simple` project and its class-average outputs.
The checks compare recovered defocus with the simulated values (maximum
error 0.10 µm), require a picked-particle count between 40% and 150% of the
48 placed particles, and verify extraction counts, class assignments,
populations, shape ranks, and finite, nonconstant populated class averages.
The particle-count check does not measure coordinate-matched picking recall.
Successful validation reports `PASS: mini_stream ...` in the log.

## Related validation suites

`simple_test_exec test=simulate_particles` validates reprojection and
particle simulation against an embedded 6VXX fixture.
`simple_test_exec test=stream_preproc` exercises simulated movies through
streaming preprocessing and local worker jobs.

The molecular workflows also have CTest entries named
`simulated_workflow_6vxx`, `simulated_workflow_1jxy`, `mini_stream_6vxx`, and
`mini_stream_1jxy`. They carry the `highlevel` label and run only when
selected explicitly. For the broader suite catalogue and instructions for
running registered tests, see the
[Test environment policy](../policies/test_environment_policy.md).
