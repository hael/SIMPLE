# Policies

A policy states the contract the code keeps today and what a change must
preserve. It describes the software as it is; the history of a change belongs in
`doc/refactoring_notes/` and `doc/implementation_notes/`. Skills in
`.github/skills/` summarise a subsystem for a developer at work; where a skill
and a policy differ, the policy and the code rule, and the skill is the one to
fix.

A policy opens with a scope statement and, when it is a draft or carries an open
validation, says so in its first lines.

## Before you add or change any step

| Policy | Read before |
|---|---|
| [layering_and_project_ownership_policy.md](layering_and_project_ownership_policy.md) | adding or changing a commander, a strategy routine, a worker step, or any code that reads or writes an `sp_project` |
| [compile_time_policy.md](compile_time_policy.md) | adding a type, a module-level `use`, a CMake target or a compile option |
| [test_environment_policy.md](test_environment_policy.md) | adding, moving or deleting a test, or changing what the fast gate runs |
| [ui_layer_policy.md](ui_layer_policy.md) | adding or changing a program, a parameter or its metadata (`src/main/ui`) |

## Workflows

| Area | Policy | Read before changing |
|---|---|---|
| 2D | [2D/solve2D_policy.md](2D/solve2D_policy.md) | `solve2D`, `refine2D`, the 2D matcher and strategy, class-average assembly |
| 2D | [2D/class_average_bootstrap_policy.md](2D/class_average_bootstrap_policy.md) | stochastic class-average expansion before `solve3D_cavgs` |
| 2D | [2D/particle_cache_policy.md](2D/particle_cache_policy.md) | the downscaled particle disk cache (`cache=yes`) |
| 3D | [3D/refine3D_policy.md](3D/refine3D_policy.md) | `refine3D`: the strategy, the matcher, the master's project, assembly, finalization |
| 3D | [3D/refine3D_auto_policy.md](3D/refine3D_auto_policy.md) | the automated wrapper around `refine3D` |
| 3D | [3D/solve3D_policy.md](3D/solve3D_policy.md) | the staged particle-based de novo map determination |
| 3D | [3D/solve3D_cavgs_policy.md](3D/solve3D_cavgs_policy.md) | the class-average route to a first map |
| 3D | [3D/solve3D_addon_policy.md](3D/solve3D_addon_policy.md) | extending a solved map with particles it did not contain |
| 3D | [3D/solve3D_chunks_policy.md](3D/solve3D_chunks_policy.md) | solving a map from all particles in balanced chunks |
| 3D | [3D/solve3D_pyramid_policy.md](3D/solve3D_pyramid_policy.md) | multistate 3D by class-average assignment, offline |
| 3D | [3D/reconstruct3D_pcg_policy.md](3D/reconstruct3D_pcg_policy.md) | the PCG reconstruction backend |
| 3D | [3D/automasking_policy.md](3D/automasking_policy.md) | automasking in refinement, `solve3D` and `volassemble` |
| 3D | [3D/classify3D_refs_policy.md](3D/classify3D_refs_policy.md) | reference-guided 3D classification |
| 3D | [3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md](3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md) | multi-state alignment and reconstruction peak memory |
| 3D | [NU/nonuniform_filtering_policy.md](NU/nonuniform_filtering_policy.md) | nonuniform filtering and its handoffs |
| Heterogeneity | [heterogeneity/refine3D_states_policy.md](heterogeneity/refine3D_states_policy.md) | particle-level conformational state refinement |
| Heterogeneity | [heterogeneity/refine3D_states_cavgs_policy.md](heterogeneity/refine3D_states_cavgs_policy.md) | sorting sub-class averages into states |
| Sampling | [importance_sampling_fractional_update_policy.md](importance_sampling_fractional_update_policy.md) | sampled updates, probabilistic candidate sampling, fractional class-average restoration, trailing reconstruction |
| Streaming | [stream/README.md](stream/README.md) | any stream stage, the pool, the IPC, restart |
| Sieving | [sieving_and_rejection/index.md](sieving_and_rejection/index.md) | particle sieve, class-average rejection, size compatibility |

## Data, I/O and infrastructure

| Policy | Read before changing |
|---|---|
| [project_stack_indexing_policy.md](project_stack_indexing_policy.md) | the indexing between project rows and physical particle stacks; `merge_projects`, `prune_project`, `fix_projfile` |
| [particle_read_path_policy.md](particle_read_path_policy.md) | how particle images are read from stack files |
| [coarray_execution_policy.md](coarray_execution_policy.md) | the coarray queue-system backend |
| [persistent_worker_policy.md](persistent_worker_policy.md) | the persistent-worker queue-system overlay |

## Preprocessing, optics and numerics

| Policy | Read before changing |
|---|---|
| [motion_correction_policy.md](motion_correction_policy.md) | reference-free motion correction |
| [motion_gain_analysis_policy.md](motion_gain_analysis_policy.md) | gain-orientation analysis and gain generation (used by the stream's preprocessing stage) |
| [phase_shift_ctf_policy.md](phase_shift_ctf_policy.md) | phase-shift CTF modelling |
| [microscope_parameters_policy.md](microscope_parameters_policy.md) | `kv`, `cs`, `fraca` and the optics model |
| [segdiampick_mics_policy.md](segdiampick_mics_policy.md) | `segdiampick_mics`, the diameter estimate of the `mini_stream` validation program (the stream picks with `segdiam_bin_picker`, see [stream/reference_generation_policy.md](stream/reference_generation_policy.md)) |
| [KB_Interpolation_Policy.md](KB_Interpolation_Policy.md) | Kaiser-Bessel interpolation on oversampled Fourier grids |
| [fsc_area_score_policy.md](fsc_area_score_policy.md) | the `fsc_area_score` diagnostic |

## Gaps

What no policy covers yet is tracked in
`doc/refactoring_notes/planned/policy_set_audit_2026-10-21.md`.
