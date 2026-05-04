# Importance Sampling and Fractional Update Policy

This document records durable workflow contracts for sampled particle updates,
probabilistic candidate sampling, fractional class-average restoration, and
trailing reconstruction in `abinitio2D`, `cluster2D`, `abinitio3D`, and
`refine3D`. It is policy, not a line-by-line implementation map.

## 1. Core Model

SIMPLE has two sampling layers that must remain separate:

1. outer fractional-update sampling chooses which particles participate in the
   current iteration
2. inner importance sampling chooses which reference, orientation, or in-plane
   candidates are explored for those participating particles

The outer subset is recorded in the project through `sampled` and `updatecnt`.
Downstream restoration and reconstruction consume that recorded state. They must
not infer participation from the nominal command-line `update_frac` alone.

Probabilistic pre-alignment is a sample-once-and-reuse path: the pre-alignment
commander chooses the outer subset, probability-table workers reuse it, and the
matcher reuses it again for the hard particle update.

## 2. Ownership

`simple_commanders_abinitio2D.f90` owns `abinitio2D` orchestration: defaults,
stage execution, final fill-in, and final class-average generation.

`simple_abinitio2D_controller.f90` owns the 2D stage policy: `NPTCLS2SAMPLE_2D`,
`nsample` override handling, stage-local `update_frac`, search-mode transitions,
and the rule that stage 1 may sample particles without fractionally restoring
previous class averages.

`simple_commanders_abinitio.f90` and `simple_abinitio_controller.f90` own 3D
stage scheduling: dynamic `update_frac`, `fillin`, `frac_best`, `balance`,
`trail_rec`, and transitions between `shc_smpl`, `prob`, and `prob_neigh`.

`simple_matcher_smpl_and_lplims.f90` owns the shared outer subset-selection
helpers for 2D and 3D. This is where full update, random sampling,
update-count-biased sampling, class-balanced sampling, fill-in sampling, and
subset reproduction are dispatched.

`simple_oris_sampling.f90` and `simple_oris_getters.f90` own the bookkeeping:
`sampled`, `updatecnt`, exact subset reproduction, global realized update
fraction, and class-local realized update fractions.

`simple_commanders_prob.f90` owns probabilistic pre-alignment orchestration:
sampling the outer subset once, writing it to the project, running table
generation, aggregating probability-table outputs, and writing the assignment
artifact.

`simple_eul_prob_tab*.f90` owns inner candidate importance sampling. These
modules may sample references, orientations, neighbors, or in-plane candidates
inside the active particle subset, but they must not choose a new particle
subset.

`simple_strategy2D_matcher.f90` and `simple_strategy3D_matcher.f90` own
particle-domain search on the active subset, assignment consumption, pose or
class updates, sigma updates during search, and writing partition-local
reconstruction or class-average inputs.

The classaverager modules own 2D class-average restoration and assembly.
`commander_volassemble` owns 3D volume assembly and trailing reconstruction.
These layers consume sampled-update state; they do not own particle selection.

## 3. Bookkeeping Contracts

`sampled` marks the current sampling round. All particles with the latest
`sampled` value belong to the current active subset.

`updatecnt` tracks cumulative update history. Count-biased and fill-in paths use
it to prefer under-updated or never-updated active particles.

`sample4update_reprod` is the only correct way to reuse a previously selected
probabilistic subset. A probability-table worker or downstream matcher must not
silently resample when a probabilistic pre-step has already sampled the subset.

`get_update_frac` returns the global realized update fraction for 3D trailing
from the current `sampled` round, active particles, and particles with
`updatecnt > 0`.

`get_class_update_fracs` returns per-class realized update fractions for 2D
class-average carry-over. It uses active particles, current class assignments,
the latest `sampled` round, and `updatecnt > 0`.

The nominal `update_frac` is a target used by sampling. The realized fraction in
`simple_oris` is the downstream restoration and trailing contract.

## 4. Abinitio2D and Cluster2D

`abinitio2D` uses a fixed run-local target sample size:

- default: `NPTCLS2SAMPLE_2D = 200000`
- override: `nsample=<integer>`

The stage controller converts that target into:

```text
update_frac_2D = min(1.0, real(min(nptcls_eff, nsample_target_2D)) / real(nptcls_eff))
```

where `nptcls_eff` is the number of active particles with `state > 0`. If the
target covers almost all active particles, the stage command omits
`update_frac` and naturally becomes a full update.

Current stage policy:

- stage 1 uses the sampled-update machinery when needed, but fractional
  class-average carry-over is disabled
- while `startit == 1`, `sample_ptcls4update2D` keeps the initial subset sticky
  by reproducing it after the first random draw
- later non-probabilistic iterations use `sample4update_cnt`, which is
  stochastic but biased toward particles with lower `updatecnt`
- probabilistic stages use `prob_align2D` to sample once, then `prob_tab2D` and
  `cluster2D_exec` reproduce the same subset
- final fill-in is assignment-only for active particles with `updatecnt == 0`;
  it uses existing class averages and does not restore new class-average sums

Fractional 2D restoration is class-local. `cavger_init_online` reads or centers
previous partial sums when fractional update is active, obtains per-class
realized fractions through `get_class_update_fracs`, and weights previous
even/odd class sums and CTF-squared sums independently for each class. This is
the 2D analogue of respecting independently updated objects in 3D.

Distributed cleanup must preserve class-average partial sums while fractional
restoration still needs them as carry-over input. Assignment and distance
artifacts are per-iteration handoffs and may be removed before the next
iteration writes replacements.

## 5. Abinitio3D and Refine3D

The 3D controller derives the outer update policy from explicit `nsample`,
explicit `update_frac`, dynamic `nsample_start`/`nsample_stop`, or default sample
count bounds. The resulting update fraction is capped by `UPDATE_FRAC_MAX`.

Current high-level ab initio stage policy:

- early stages use `shc_smpl`
- later stages use `prob`
- final neighborhood stages use `prob_neigh`
- final active stages may switch to `fillin`, except where the multi-state
  policy disables it

`sample_ptcls4update3D` applies the normal 3D subset policy:

- if fractional update is off, select all active particles
- if `balance=yes`, use class-balanced sampling
- otherwise use update-count-biased sampling

`sample_ptcls4fillin` is a separate late-stage coverage policy. Its purpose is
to update particles with insufficient history, not to preserve the normal
balanced or count-biased exploration distribution.

The 3D matcher writes partial reconstructions from the active subset. Volume
assembly then restores volumes, calculates FSCs, postprocesses references, and
applies trailing reconstruction when requested. Trailing uses explicit
`ufrac_trec` if provided; otherwise it consumes the realized fraction from
`get_update_frac`.

## 6. Probabilistic Pre-Alignment

Probabilistic pre-alignment is not a second outer sampler.

The workflow is:

1. choose the outer subset through the normal 2D or 3D sampling helper
2. write the sampled project state
3. run probability-table generation only for that subset
4. aggregate table outputs into one assignment artifact
5. reproduce the same subset in the matcher
6. perform the hard particle update

`simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`, and
`simple_eul_prob_tab2D.f90` perform candidate-level importance sampling inside
that subset. They may use score-derived candidate distributions,
`angle_sampling`, `greedy_sampling`, or neighborhood sampling, but the selected
particle set is already fixed before they run.

## 7. Restoration and Assembly

2D class-average restoration consumes class-local realized update fractions:

- previous class contribution: `1 - rho(class)`
- current class contribution: the new partial sums for that class

Classes with no active updated particles keep a zero realized fraction. Classes
with full sampled participation replace previous sums.

3D volume assembly consumes the global realized update fraction for trailing:

- previous volume contribution: `1 - update_frac_trail_rec`
- current volume contribution: `update_frac_trail_rec`

Neither class-average restoration nor volume assembly should make new particle
sampling decisions. If a restoration or assembly change requires a different
subset policy, that policy belongs in the commander/controller/sampling-helper
layer and must be reflected in `sampled` and `updatecnt`.

## 8. Invariants

- Outer particle sampling happens before probabilistic table generation.
- Probabilistic table workers and downstream matchers reproduce the same subset.
- Candidate importance sampling never changes the particle subset.
- `sampled` remains the current-round marker.
- `updatecnt` remains cumulative update history.
- Downstream restoration uses realized update state, not only nominal
  `update_frac`.
- Stage 1 of `abinitio2D` may be sampled but must not fractionally carry over
  previous class-average sums.
- 2D fractional class-average restoration remains class-local.
- Final `abinitio2D` fill-in remains assignment-only unless the policy is
  explicitly changed.
- `volassemble` and the classaverager remain consumers of sampled-update state,
  not producers of particle-selection policy.

## 9. Review Checklist

For sampling, probabilistic alignment, class-average restoration, or volume
assembly changes, check:

- Does the outer subset get selected exactly once for a probabilistic
  pre-alignment iteration?
- Do table workers and matchers reuse the recorded subset through
  `sample4update_reprod`?
- Is candidate-level importance sampling kept separate from particle-level
  subset selection?
- Are `sampled` and `updatecnt` updated consistently before downstream
  restoration or trailing consumes them?
- Does 2D restoration use class-local realized fractions?
- Does 3D trailing consume the realized or explicit trailing fraction?
- Are shared-memory and distributed paths preserving the same scientific
  workflow and artifact contracts?
