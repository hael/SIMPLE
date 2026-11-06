# Glossary

The names SIMPLE's code uses today for its recurring concepts, and what they mean. One concept often
has several names; new code uses the first name given. This page records the code as it is. The domain
language of the [domain-driven design proposal](refactoring_notes/planned/next_release/staged_domain_driven_design_refactor.md)
(its section 4) is a separate, future vocabulary. The folder guides (`.github/skills/simple-main-*`) state
the rules around each term.

## Workflows

| Term | Meaning |
| --- | --- |
| `refine3D` | 3D refinement: align particle images to reference volumes and reconstruct new ones, iteration by iteration |
| `solve3D` | De novo map determination: ab initio 3D reconstruction coupled with initial 3D refinement, run as stages of `refine3D` (`src/main/solve/`). Named `abinitio3D` before 2026-10-03 |
| `refine2D` | 2D classification: align particles to class averages and average them again; the 2D counterpart of `refine3D`. Named `cluster2D` before 2026-10-03 |
| `solve2D` | The 2D counterpart of `solve3D`: stages of `refine2D` from no references |
| stream | The streaming pipeline (`simple_stream`): preprocessing, picking and 2D classification as data arrive |
| sieve | The stream's particle sieve: chunks of particles classified in 2D, bad classes rejected |
| pass | One `refine3D_exec` or `refine2D_exec`: one iteration's search over the sampled particles, in one process or one part |
| master, worker, part | A distributed run splits the particles into `nparts` parts; each worker (`part=` given) processes one part and writes part files; the master merges them and alone writes the project |

## Data and records

| Term | Meaning |
| --- | --- |
| project, `projfile`, `sp_project` | The `.simple` file and its in-memory type: segments of records |
| segment, `oritype` | A table of the project: `mic` (micrographs), `stk` (particle stacks), `ptcl2D` and `ptcl3D` (one row per particle, for 2D and 3D work), `cls2D` and `cls3D` (one row per class average), `out` (registered outputs), `optics`, `projinfo` |
| `ori`, `oris` | One record (a row) and a table of records; particle records have fixed named slots (`I_*`), others a typed hash |
| `state` | The 3D label of a particle (`ptcl3D`) or its selection (`ptcl2D`); state 0 excludes the row (rejected) |
| `class`, `proj` | The 2D class of a particle; a projection direction (an index into the direction grid). 0 means unassigned |
| even/odd half: `eo`, `ieo`, `iseven` | The two independent halves of the data behind the gold-standard FSC. The key `'eo'` is 0 even, 1 odd, -1 none; loops use `ieo` (1 even, 2 odd) and `iseven` (logical). `get_half` returns 0 or 1 and stops on -1 |
| `sampled`, `updatecnt` | Fractional update bookkeeping on particle rows: the round in which a row was last drawn, and how many times it has been drawn |
| `update_frac` | The fraction of particles a pass updates (fractional update); with `trail_rec=yes` the previous reconstruction is carried as a weighted average (trailing reconstruction) |
| cavgs | Class averages (2D) |
| `pinds` | Particle indices: row numbers of the particle segment |

## Images and Fourier space

| Term | Meaning |
| --- | --- |
| `box`, `smpd` | Image side length in pixels, and the sampling distance (pixel size) in Ångström per pixel |
| `box_crop`, `smpd_crop` | The same for the Fourier-cropped (downscaled) images and volumes a search or reconstruction works on; `box*smpd == box_crop*smpd_crop` |
| `box_croppd`, `OSMPL_PAD_FAC` | The padded box of gridding (oversampling factor `OSMPL_PAD_FAC`) |
| `ldim` | Image dimensions `[nx, ny, nz]` (`nz` = 1 for 2D) |
| `lp`, `hp` | Low-pass and high-pass limits in Ångström |
| `kfromto` | The band of Fourier shells (indices) a search uses: `kfromto(1)` (from `hp`) to `kfromto(2)` (from `lp`) |
| `mskdiam` | Mask diameter in Ångström |
| FSC | Fourier shell correlation between the two half maps, per shell; resolution is read at 0.143 (or 0.5) |
| CTF | Contrast transfer function of the microscope: defocus `dfx`, `dfy`, astigmatism angle `angast`, phase shift |
| `pgrp` | Point-group symmetry (`c1`, `cn`, `dn`, `t`, `o`, `i`) |
| polar transforms, `pftc` | Particle and reference images resampled on polar Fourier grids for fast rotational search (`polarft_calc`) |
| `cftc` | The Cartesian Fourier calculator of continuous pose refinement (`refine=cont`) |
| `nspace` | The number of projection directions of the search grid (`nspace_sub` of its coarse subspace) |
| `trs` | The shift search range: maximum half-width in pixels |

## Search and scoring

| Term | Meaning |
| --- | --- |
| `refine` | The search mode of a pass: 3D `shc`, `neigh`, `greedy`, `prob`, `prob_neigh`, `cont`; 2D `snhc`, `snhc_smpl`, `prob`, ... |
| `objfun` | The score: `euclid` (noise-weighted Euclidean distance, needs `sigma2`) or `cc` (correlation) |
| `sigma2` | The noise power spectrum per even/odd half that weights the Euclidean score (`euclid_sigma2`) |
| probability table | The candidates and scores of `refine=prob*` from which a particle's pose or class is drawn (`probabilistic/`) |
| keyed draw | A random draw whose value depends only on the run seed, the iteration, the particle and the purpose (`rnd_keyed`), so it does not depend on threads, batches or parts; `SIMPLE_SEED` sets the seed |
| `nthr`, `nthr_glob`, team | The thread count the parameters give, its process-wide copy, and the OpenMP team a parallel region runs |

## Code layers

| Term | Meaning |
| --- | --- |
| `params`, `parameters` | The typed parameters of a run, built once from the command line (`params%new(cline)`) |
| `cline`, `cmdline` | The command line: keys and values given to a program |
| `build`, `builder` | The run's toolbox: the project, the orientation field (`spproj_field`), images, calculators |
| `ui` | Program and parameter descriptions (help, GUI) |
| exec router | Maps a program name to its commander (`src/main/exec`) |
| commander | A program's front door: normalises the command line, then runs a strategy (`execute(cline)`) |
| strategy | How a program runs (shared memory, distributed master, worker) and the master's steps; `run_<x>(cline)` is a strategy routine's entry |
| fast gate | `scripts/run_fast_gate.sh`: static checks and the fifteen unit suites, about 30 seconds |

## Name conventions

| Pattern | Meaning |
| --- | --- |
| `_8` | Double-precision variant of a routine (`gen_corr_for_rot_8`) |
| `l_` | A logical flag (`l_distr_worker_glob`, `params%l_update_frac`) |
| `_glob` | A process-wide module variable (`nthr_glob`, `CWD_GLOB`) |
| `new` / `kill` | Constructor and destructor of a type; every `new` has its `kill` |
| `_exec` | A step on the caller's `(params, build)` (`refine3D_exec`, `volassemble_exec`) |
| generic and specific | A type-bound generic (`set`, `add`) dispatches by argument types to specifics named by what distinguishes them (`set_real`, `add_img`) |
