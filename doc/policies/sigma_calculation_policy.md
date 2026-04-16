# Sigma Calculation Policy

## Scope

This document defines the current policy for Euclidean sigma calculation and use in SIMPLE.

It applies to workflows that run with `objfun='euclid'`, especially:

- `cluster2D`
- `abinitio2D`
- `refine3D`
- reconstruction paths that consume Euclidean sigma spectra

## Public policy

Sigma calculation is a Euclidean-objective feature.

- If `objfun='euclid'`, SIMPLE computes and persists sigma spectra.
- If `objfun!='euclid'`, the Euclidean sigma path is inactive.
- Sigma application during restoration/reconstruction is gated separately by `ml_reg='yes'`.

The effective reconstruction gate is `params%l_ml_reg`, which is only true when:

- `ml_reg='yes'`
- `cc_objfun == OBJFUN_EUCLID`

This is enforced in [src/main/simple_parameters.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_parameters.f90:1761).

## Ownership and persistence

The runtime owner of sigma state is `builder%esig`, the `euclid_sigma2` object held by the builder. That object owns:

- `sigma2_noise`: group-expanded per-particle spectra used by alignment/reconstruction
- `sigma2_part`: per-particle spectra written during the current search/update pass
- `sigma2_groups`: group spectra loaded from the STAR representation

See [src/main/simple_builder.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_builder.f90:30) and [src/main/simple_euclid_sigma2.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_euclid_sigma2.f90:15).

Persisted sigma data has two forms:

- partition-local binary files: `sigma2_part_<part>.dat`
- grouped STAR files: `sigma2_group_iter_<iter>.star`

The binary format is defined in [src/fileio/simple_sigma2_binfile.f90](/Users/elmlundho/src/SIMPLE/src/fileio/simple_sigma2_binfile.f90:1).

## Grouping policy

`sigma_est` controls grouping:

- `group`: group by particle `stkind`
- `global`: use one shared group per even/odd half-set

This is normalized into `params%l_sigma_glob` in [src/main/simple_parameters.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_parameters.f90:1726).

When group STAR files are loaded, the selected group spectrum is expanded back onto each active particle according to:

- particle `eo`
- particle `stkind` unless `l_sigma_glob`

See [src/main/simple_euclid_sigma2.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_euclid_sigma2.f90:171).

## Calculation lifecycle

### 1. Bootstrap spectra

If the expected grouped sigma STAR file for the current starting iteration does not exist, SIMPLE bootstraps sigma with `calc_pspec`.

In shared-memory `refine3D`, this check is performed against `sigma2_group_iter_<which_iter>.star` before entering the main loop; if missing, `calc_pspec` is executed after enforcing even/odd partitioning and, for `startit==1`, unit particle weights. See [src/main/strategies/parallelization/simple_refine3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_refine3D_strategy.f90:186).

`calc_pspec` computes per-particle noise power spectra from masked particle images and writes `init_pspec_part*.dat`. In the bootstrap assembly step:

- spectra are averaged per even/odd and group
- the average image power spectrum is subtracted
- negative values are repaired by neighbor propagation
- grouped STAR output is written
- partition-local `sigma2_part_<part>.dat` files are rewritten with the grouped values

See [src/main/strategies/parallelization/simple_calc_pspec_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_calc_pspec_strategy.f90:130) and [src/main/commanders/simple/simple_commanders_euclid_distr.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_euclid_distr.f90:15).

### 2. Per-particle update during search

During Euclidean 2D/3D matching, sigma is recalculated after the particle orientation/shift update.

- 2D uses `refkind='class'`
- 3D uses `refkind='proj'`

The residual is formed in polar Fourier space after shift, rotation, and CTF application, and the per-shell squared residual energy is normalized by `2 * pftsz`. See [src/main/pftc/simple_polarft_corr.f90](/Users/elmlundho/src/SIMPLE/src/main/pftc/simple_polarft_corr.f90:559).

Call sites:

- [src/main/strategies/search/simple_strategy2D_matcher.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/search/simple_strategy2D_matcher.f90:157)
- [src/main/strategies/search/simple_strategy3D_matcher.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/search/simple_strategy3D_matcher.f90:145)

At the end of the search pass, the updated per-particle spectra are written back to the partition binary file by `build%esig%write_sigma2`.

### 3. Iteration consolidation

The policy is that grouped STAR files are the durable handoff artifact between iterations and runs.

Current implementation details differ slightly by workflow:

- `cluster2D` consolidates after the iteration and writes `which_iter + 1` for the next iteration. See [src/main/strategies/parallelization/simple_cluster2D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_cluster2D_strategy.f90:296).
- shared-memory `refine3D` consolidates `which_iter` at iteration entry, then writes `which_iter + 1` again at finalization so a subsequent run can restart from grouped sigma state. See [src/main/strategies/parallelization/simple_refine3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_refine3D_strategy.f90:251) and [src/main/strategies/parallelization/simple_refine3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_refine3D_strategy.f90:362).
- distributed `refine3D` also writes `which_iter + 1` at finalization. See [src/main/strategies/parallelization/simple_refine3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_refine3D_strategy.f90:873).

The authoritative consolidation implementation is `exec_calc_group_sigmas` in [src/main/commanders/simple/simple_commanders_euclid.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_euclid.f90:41).

That implementation:

- reads all partition `sigma2_part_<part>.dat` files
- groups by `eo` and either `stkind` or one global bucket
- computes weighted group averages using particle weight `w`
- writes `sigma2_group_iter_<which_iter>.star`

This is intentionally different from bootstrap `calc_pspec_assemble`, which uses unweighted counts plus average-spectrum subtraction and negative repair.

## Consumption policy

Grouped sigma spectra are loaded into `builder%esig`, expanded per particle, and then consumed only on ML-regularized restoration/reconstruction paths.

For 3D reconstruction, `params%l_ml_reg` causes the sigma groups to be loaded before `calc_3Drec`. See [src/main/commanders/simple/simple_commanders_rec.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_rec.f90:67) and [src/main/strategies/parallelization/simple_rec3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_rec3D_strategy.f90:149).

For 2D class restoration, the same grouped sigma state is loaded before Fourier-plane generation. See [src/main/class/simple_classaverager_restore.f90](/Users/elmlundho/src/SIMPLE/src/main/class/simple_classaverager_restore.f90:79).

The Fourier-plane policy is:

- upsample sigma from the cropped shell grid to the padded reconstruction grid
- divide both the complex Fourier sample and its CTF power term by sigma

See [src/main/image/simple_image_ctf.f90](/Users/elmlundho/src/SIMPLE/src/main/image/simple_image_ctf.f90:203) and [src/main/image/simple_image_ctf.f90](/Users/elmlundho/src/SIMPLE/src/main/image/simple_image_ctf.f90:263).

## Current architectural boundary

The current boundary is acceptable and should be preserved:

- `simple_euclid_sigma2`: sigma object model and STAR/binary translation
- matcher/search code: per-particle sigma updates
- Euclid commanders/strategies: bootstrap and iteration-level consolidation
- reconstruction/restoration code: sigma consumption when `l_ml_reg`

## Recommended evolution

Near-term improvements should preserve numerical behavior:

- make iteration-number semantics explicit in one helper so `cluster2D` and `refine3D` use the same documented convention
- separate bootstrap-only policy from steady-state consolidation policy, because they intentionally use different averaging rules
- add regression coverage for restart/resume behavior using preexisting `sigma2_group_iter_<iter>.star`

The main rule to preserve is simple: grouped STAR sigma files are the cross-iteration contract, while partition binary sigma files are working state for the current pass.
