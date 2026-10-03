# Trailing reconstruction without the finished-halfmap blend

Date: 2026-10-03

Status: proposed; implementation not started. Release 4 cleanup item C6.

Validation level: static source inspection only. No source code was changed,
compiled or executed while preparing this note. Line references are to the
working tree of 2026-10-03 and will drift.

This is the single living design record for this refactor. Update it as the
implementation and validation land rather than creating companion plans.

## 1. Ruling

The developer's ruling of 2026-10-03:

> Trailing needs to be done consistently across using the newly implemented
> approach. No trailing should ever happen on finished halfmaps. If the
> solution to correct weighting and consistent application is one
> reconstruction, I am ok with it.

"The newly implemented approach" is the accumulator-domain trailing chain:
blended, unregularized even/odd Fourier sums and sampling densities at
full-dataset mass, restored once after the blend. "Finished halfmaps" are
restored maps after sampling-density correction, regularization and
deapodization. The ruling removes every path that blends such maps.

## 2. Current behaviour

### 2.1 The accumulator chain (kept)

For `trail_rec=yes`, `volassemble` keeps one chain per state
(`trailrec_stateNN_{even,odd}` plus rho files plus the `trailrec_stateNN.txt`
manifest). With the realized fraction `f` of the current partials and the
applied update weight `u` (`ufrac_trec` for a single state, otherwise `f`), the
population rule (`population_blend_weights`) scales the current partials by
`u/f` and the chain by `(1-u)N/M`. One restoration of the blended sums yields
halves whose current-map coefficient is exactly `u`, and the FSC is estimated
after the blend. When `u >= 0.99`, nothing is blended: the chain is rewritten
from the current partials at full mass. Code:
`simple_commanders_rec_distr.f90::restore_state_from_parts`,
`blend_trailing_accumulators`. The PCG backend keeps an equivalent raw
accumulator chain pair (`refine3D_pcg_trail_accum_fname`) with the same
population rule (`simple_rec3D_pcg_strategy.f90::set_chain_blend_weights`).

### 2.2 The bootstrap blend (to be removed)

When no valid chain exists and `u < 0.99`, both backends fall back to a
volume-domain bootstrap:

- **Gridding** (`simple_commanders_rec_distr.f90`):
  - `restore_eos_and_write_fsc` reads the previous half maps from
    `vol1..N` on the command line (`read_previous_halfmaps`). Their FSC drives
    the ML regularization of the new halves.
  - `trail_restored_halves_if_needed` then writes
    `u*new + (1-u)*previous` for the restored halves, the NU base and auxiliary
    inputs, and, under `lp` set, the merged volume.
  - `blend_trailing_accumulators` meanwhile seeds the chain from the current
    partials scaled by `1/f`.
- **PCG** (`simple_rec3D_pcg_strategy.f90::execute_rec3D_pcg_distributed_master`):
  - with `l_bootstrap`, the FSC pair is the previous shipped pair
    (`load_previous_state_halves`);
  - `blend_bootstrap_half` blends the solved halves, the ML pair and the
    solvent pair with the previous pair;
  - the support provenance of the blend is combined from both contributions;
  - the chain is seeded from the current raw accumulators (around line 1084).
  - `trail_bootstrap_states` reports the bootstrap per state to
    `filter_pcg_nonuniform_maps`, which only checks the array's size.

### 2.3 When the bootstrap runs

The bootstrap needs both an absent or invalid chain and `u < 0.99`. In a fresh
run, the first trailing iteration has `f = 1`: the updated pool is the current
sample. The chain is then seeded with nothing blended. The bootstrap therefore
runs only when the updated pool is larger than the current sample and no valid
chain exists:

- a refine3D run with `trail_rec=yes` that starts on a project that already has
  update history but has no chain in its directory (for example probabilistic
  modes, which keep `updatecnt`; or `startit > 1`);
- refine3D_auto and other commanders that switch `trail_rec` on mid-run
  (`simple_commanders_refine3D.f90`, the per-iteration `l_trail_rec` updates);
- a chain discarded on validation (particle population, state layout, larger
  grid, physical extent, or corrupt or mixed-generation files);
- a state whose realized fraction is about zero (`realized_update_frac < 0.001`
  routes that state to the bootstrap path).

solve3D normally avoids the bootstrap: the full reconstruction at each stage
boundary seeds the chain with `trail_seed=yes` before a trailing stage
(`simple_solve3D_utils.f90::calc_rec`). So does the distributed refine3D start
with PCG and `objfun=cc` when no volume is given
(`simple_refine3D_strategy.f90::distr_initialize`).

### 2.4 Why it is inconsistent

- **Two different blends in one trailing run.** The bootstrap iteration blends
  finished maps, after regularization and deapodization, which are not linear
  operations. Every later iteration blends raw sums before restoration. The
  results are not the same estimator.
- **A different FSC source.** On the bootstrap iteration, the regularization
  FSC comes from the previous half maps, not from the data being restored.
- **The chain loses the previous model.** The bootstrap output contains `1-u`
  of the previous maps, but the chain seeded at the same moment holds only the
  current sample, scaled by `1/f` to full mass. From the next iteration on,
  the previous model's information is gone from the chain, and the chain claims
  the mass of `N` particles with the noise of `n = f*N`.
- **A hidden command-line dependency.** The bootstrap requires the previous
  volumes on the command line; volassemble stops if `vol<state>` is missing.

## 3. Proposed design

### 3.1 Principle

A blend (`u < 0.99`) always uses a valid chain. When the chain that a blend
needs is missing, it is created by one full reconstruction, never by blending
maps.

### 3.2 Seeding at assembly time

1. **Decide after matching, in the strategy.** After matching (the sample, and
   thus `f`, is known only then), the refine3D strategy decides before
   assembly. Both the shared-memory and the distributed strategy do this, at
   the point where they now call `volassemble` or `assemble_refine3D_pcg`.
   For each populated state:
   - compute the updated pool `N`, the sample `n`, `f` and `u`, exactly as
     volassemble does (`get_group_update_counts`, `get_state_update_fracs`,
     `ufrac_trec`);
   - check the backend's chain for validity.

   Seeding is needed when any populated state has `u < 0.99` and no valid
   chain.
2. **Seed with one full reconstruction.** When seeding is needed, the strategy
   runs reconstruct3D in-process instead of assembling the partials. Its
   command line is the iteration command line with:
   - `prg=reconstruct3D` and `mkdir=no`;
   - `trail_rec`, `update_frac`, `ufrac_trec` and `fillin` deleted;
   - `trail_seed=yes`;
   - `which_iter` set.

   reconstruct3D already supports this on both backends: its partials come from
   `sample4rec`, which takes every active row with `updatecnt > 0`, at its
   current pose, including this iteration's updates. That is exactly the
   population the chain represents, so the chain is written at full mass
   (`TC_NSEED` from `get_state_rec_pops`). The distributed strategy's command
   line carries `nparts`, so the seed is distributed as well.
3. **Outputs.** The seed reconstruction writes the iteration's output volumes,
   half maps and FSC files under the names assembly would have used. The rest
   of the iteration (volume naming, postprocess, NU low-pass handoff,
   convergence, sigma2 commit) proceeds unchanged. On this iteration the map is
   the full reconstruction, so there is no lag and no weighting decision; from
   the next iteration on, every blend uses the chain.
4. **Cost.** One extra full particle pass in the seeding iteration only. The
   matcher's partials of that iteration go unused (they are written before `f`
   is known). This is a deliberate, approved exception to the single-read I/O
   contract.
5. **Sigma2 consistency.** The canonical sigma2 commit happens after assembly,
   so the seed reconstruction reads the same committed state as the matcher's
   partials did. The seed uses the run's objective (`objfun`), not the forced
   `cc` of the start-up reconstruction.

### 3.3 Backend changes

- **Gridding (`simple_commanders_rec_distr.f90`):**
  - remove `read_previous_halfmaps`, `trail_restored_halves_if_needed`, the
    bootstrap branch of `restore_eos_and_write_fsc` (the FSC is always
    estimated post-blend), the `vol_prev_even`, `vol_prev_odd` and `vol_merged`
    arguments, and the bootstrap seeding comment;
  - in `blend_trailing_accumulators`, keep the `u >= 0.99` branch that rewrites
    the chain from the current partials;
  - make "trail_rec, `u < 0.99`, no valid chain" a hard error naming the
    strategy contract;
  - for a state with `realized_update_frac < 0.001` and a valid chain, carry
    the chain unchanged and restore from it, as the solve3D_addon cohort path
    already does.
- **Chain validity for both callers.** Move the chain validation
  (`validate_trail_chain` and the helpers it uses) out of the contained scope
  into a public module procedure that has no side effects. Its result is the
  validity plus the generation and stored mass, and volassemble discards the
  set on failure as now. The same procedure serves the strategy's decision.
- **PCG (`simple_rec3D_pcg_strategy.f90`):**
  - remove the `l_bootstrap` branches: the previous FSC pair,
    `load_previous_state_halves`, `blend_bootstrap_half`, the support
    combination for blended pairs, and the bootstrap chain seeding;
  - a missing chain pair where a blend is needed becomes a hard error;
  - expose a side-effect-free chain-validity check, built from
    `discard_stale_trail_chain_pair` and the pair-existence test;
  - remove the `trail_bootstrap_states` argument, and the `l_trail_bootstrap`
    argument of `filter_pcg_nonuniform_maps` together with its callers in
    `simple_refine3D_strategy.f90` and `simple_rec3D_strategy.f90`.
- **Strategies (`simple_refine3D_strategy.f90`):** a shared private helper
  takes the seeding decision (3.2.1) and builds and runs the seed command line
  (3.2.2). The in-memory strategy writes the particle segment before seeding,
  so that reconstruct3D reads the current poses from the project file.

### 3.4 What does not change

- The population rule, the `u/f` scaling, the chain manifest and the
  integrity checks.
- The carry-over of complete chains on `continue=yes`.
- solve3D stage-boundary seeding.
- The solve3D_addon frozen paths, which already require a seeded chain.
- The `u >= 0.99` path.

## 4. Alternative considered: seeding at the start of the iteration

Seeding before matching would let the first trailing iteration blend with the
proper weight `u` against a full reconstruction at the previous poses. It is
not recommended: reconstruct3D writes the same state volumes, half maps and FSC
files that the iteration is about to use as matching references, so the seed
would silently replace them. The decision would also have to be made before `f`
is known, so seeds would be taken in iterations that need no blend.

## 5. Documents and skills to update

- `.github/skills/simple-frac-update-trailing/SKILL.md` and
  `references/frac-update-contract.md`: replace the bootstrap paragraph with
  the seeding contract.
- `.github/skills/simple-refine3d/SKILL.md`: one line on the seeding exception
  to the single-read contract.
- `doc/policies/importance_sampling_fractional_update_policy.md` and
  `doc/policies/3D/reconstruct3D_pcg_policy.md`: remove the bootstrap blend
  and the lag-one FSC pair; describe the seed.
- `doc/policies/3D/refine3D_policy.md`: the strategy's assembly step.
- The release 4 inventory: set C6 to done when this lands.

## 6. Tests

- `src/main/image/simple_accum_blend_tester.f90`
  (`run_bootstrap_then_update`, "post-bootstrap effective update weight equals
  realized fraction f"): replace with a seed test. A full-reconstruction seed
  of mass `N`, followed by an update with `f` and `u`, restores with
  current-map coefficient `u`. Keep the chain-mode recurrence cases.
- A strategy-level unit test of the seeding decision on synthetic oris:
  - `f = 1` and no chain: no seed;
  - `f < 1` and a valid chain: no seed;
  - `f < 1` and no chain: seed;
  - `f` about 0 and a valid chain: no seed, chain carried.
- Run `scripts/check_test_registry.py` after the test changes.

## 7. Validation plan (developer runs)

On the Dell, with ground truth per the validation memory (orientation errors
against truth for small noisy sets, FSC for larger ones):

1. **refine3D, `trail_rec=yes`, `update_frac=0.2`, gridding**, started from a
   project with update history and no chain. Expect one seeding iteration (log
   line), then accumulator blends. Compare orientation error and FSC per
   iteration with the same run before the change.
2. **The same run with `rec_backend=pcg`.**
3. **refine3D_auto switching `trail_rec` on mid-run.** Expect exactly one seed.
4. **A multi-state run where one state receives no sample in some iteration.**
   Expect that state's chain to be carried unchanged.
5. **solve3D end to end.** Expect no seeding iteration, because the stage
   boundaries seed, and maps unchanged within noise.

## 8. Open questions

- **Re-seeding scope.** reconstruct3D seeds all states together, so a seeding
  iteration also replaces the valid chains of the other states with full
  reconstructions. That is acceptable because seeding is rare, but it could be
  narrowed per state later.
- **States with active particles but none ever updated.** Such a state gets no
  rows from `sample4rec` and therefore no chain. Its next sampled iteration has
  `f = 1` and seeds itself through the `u >= 0.99` path, so no loop arises, but
  this should be confirmed in the multi-state run.
- **Seed log and accounting.** Whether the seeding iteration should be marked
  in the convergence output (`TRAIL_REC_UPDATE_FRACTION`), so the jump to a
  full update is visible in run reports.
