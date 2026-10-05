# Class-Average State, Trailing-Chain Mass and Reconstruct3D Partials Refactoring

Date: 2026-09-05

Revised: 2026-09-11 after the canonical sigma2 cutover; 2026-09-28 after a
source review of 2D carry-over, and again on 2026-09-28 after the same review
of the 3D trailing chains.

Status: Changes 1 to 6 implemented on 2026-09-28, rebased onto `5fa147fc7` and
validated on 2026-09-29 (uncommitted work tree, for review). Decisions are in
Section 9; implementation findings and validation results are in Section 12.

Purpose: single living design and validation record for:

- removing an invalid reconstruct3D continuation check;
- retiring a dead streaming reconciliation path;
- correcting the class update fraction in distributed 2D;
- moving 2D class-average carry-over out of partition-shaped files into one
  blend at the assembly owner;
- making the 2D and 3D fractional blends keep the sampling mass of the
  population they represent.

Related contracts:

- [Canonical Sigma2 State Refactoring](../completed/canonical_sigma2_state_refactoring.md)
- [Solve2D Policy](../../policies/2D/solve2D_policy.md)
- [Class-Average Bootstrap Policy](../../policies/2D/class_average_bootstrap_policy.md)
- [Importance Sampling and Fractional Update Policy](../../policies/importance_sampling_fractional_update_policy.md)
- [Sampling and Fractional Updates](../../algorithms/sampling_and_fractional_updates.md)
- [Refine2D Class Averaging](../../algorithms/refine2d_class_averaging.md)
- [Refine3D Policy](../../policies/3D/refine3D_policy.md)
- [Reconstruct3D PCG Policy](../../policies/3D/reconstruct3D_pcg_policy.md)

## 1. Decision Summary

Six independently reviewable changes:

1. **Remove the gridding continuation file-count check** in refine3D. It
   counts transient files and rejects every valid multi-state continuation
   (Section 3.7).
2. **Retire the dead streaming chunk reconciliation.** `stream_chunk%read`
   and its `average_into` have no caller (Section 3.6).
3. **Correct the distributed class update fraction now.** Workers compute
   the fraction with a global denominator and a partition-local numerator,
   so previous class sums decay by about `1 - f/nparts` per iteration
   instead of `1 - f` (Section 3.2). Interim fix: each worker counts over its
   own particle range, pinned by a unit test. Change 4 replaces it.
4. **Blend 2D carry-over once, at the assembly owner, under the population
   rule.** Workers accumulate only the current iteration from zero. The
   owner sums the parts and reads one partless previous set that records
   the population it represents. It computes the counts from the merged
   project, then shifts, scales and adds once, writes the new set, and
   restores from it. This removes the dependence on `nparts`, gives
   distributed runs the shared-memory result, and keeps sampling mass equal
   to the represented population (Sections 3.3 and 4.1).
5. **Apply the same population rule to the 3D trailing chains**, gridding
   and PCG. Today the chains lose mass whenever the represented population
   changes (Section 3.3). Each chain records the population it represents.
6. **Make `sample4rec`'s "nothing updated yet" decision global.** Today each
   worker decides from its own range (Section 3.13).

Withdrawn from the 2026-09-11 plan, because they give sigma2-grade guarantees
to state that is an approximation by construction (Section 4.2):

- per-particle layout and checkpoint digests;
- generation-scoped worker payloads with coverage validation;
- project registration of a `cavg_state` path;
- directory sync and the crash-order trace between project and state writes;
- a full rebuild whenever particles are activated, deactivated or reordered.

A stale or partly stale aggregate is the error the recurrence already accepts
when particles change class or state, and it decays with every blend. The
recovery path for anything the owner cannot read or match is one full update
in 2D (Section 5.3), and the existing re-seed in 3D.

No public CLI selector will be introduced. Each workflow keeps one production
carry-over path, as in the sigma2 rollout.

## 2. Lessons Carried Forward from Canonical Sigma2

| Sigma2 lesson | Application here |
| --- | --- |
| Separate authoritative state from derived views. | Persist unregularized Fourier numerator and CTF-squared (2D) or density (3D) sums. Restored averages, maps and FSC/FRC stay derived outputs. |
| Durable identity must be scientific, not scheduler-shaped. | The 2D previous set has no part number and is independent of `nparts`, `numlen` and execution mode, as the 3D chains already are. |
| Workers never touch committed state. | Workers write current-only contributions. Only the owner reads, blends and writes carried state. |
| Publish a scientific unit together. | The four 2D arrays are published together (Section 5.2), through a temporary name and a rename. The 3D chain keeps its manifest-last scheme. |
| Audit every handoff. | Class centering, the pool's crop-box upsample, chain seeding, cleanup and continuation are part of the change. |
| Cut over vertically and remove the legacy path. | Shared memory, distributed 2D, the pool and the 3D chains are migrated in turn, and legacy part-file carry-over reads are removed. |

Not copied: sigma2 stores exact per-particle records, so it can validate
identity row by row and rebuild exactly. A class or state aggregate has no
smaller exact source. It is the state of a stochastic recurrence, so row-level
identity checks would guard a precision the stored data do not have.

## 3. Current-System Review Findings

Verified against `master` on 2026-09-28, by reading the code.

### 3.1 2D carry-over is partition-shaped (latent)

In `cavger_init_online` (`src/main/class/simple_classaverager_restore.f90`)
each worker reads the four `cavgs_*_partN` and `ctfsqsums_*_partN` files for
its own `part`, scales them by `1 - f(c)`, accumulates its current
particles, and rewrites the same files. `cavger_assemble_sums_from_parts`
then sums the part files. Shared memory writes the same files with part 1.

Shrinking `nparts` loses the removed parts' previous mass. Growing it, or
changing `numlen` (for example from 9 to 10 parts), asks for files that do
not exist. Cleanup (`cleanup_distributed_iteration_artifacts` in
`simple_cluster2D_strategy.f90`, `qsys_cleanup`) enumerates only the current
`nparts`, so files from a wider layout survive.

This defect is latent. `nparts` is constant within every current workflow:
abinitio2D stages, the streaming pool (`nparts_pool = nparts`) and chunks,
which start fresh. It triggers on `continue=yes` with a different `nparts`,
or on a switch between shared and distributed mode.

### 3.2 Distributed 2D workers use the wrong class update fraction (live)

When the sample is drawn inside the worker, `sample_ptcls4update2D` samples
only `[fromp, top]` and stamps the chosen rows with the global maximum
`sampled` marker plus one. `cavger_init_online` then calls
`get_class_update_fracs` (`src/main/ori/simple_oris_getters.f90`) on the
whole project, which every worker reads in full:

- the numerator counts rows of class `c` carrying the new marker, and only
  this worker's rows carry it;
- the denominator counts updated active rows of class `c` in all
  partitions.

Each worker's fraction is therefore its own share of the class's sample,
about `f(c)/nparts`. Summed over workers, previous mass decays by about
`1 - f(c)/nparts` per iteration instead of `1 - f(c)`. With
`update_frac = 0.1` and four parts, that is 0.975 instead of 0.9: old
contributions take about 27 iterations to halve instead of about 7, and
total sampling mass grows above that of the dataset.

Affected: every distributed iteration in which workers draw a new sample
(`startit > 1`, `sample4update_cnt`). That covers:

- the streaming pool (`refine=snhc_smpl` set in
  `simple_stream_p06_pool2D_new.f90`, `startit = pool_iter`);
- abinitio2D stages below `PROBREFINE_STAGE` (`snhc_smpl`, chosen in
  `set_cluster2D_stage_search_policy`) when `update_frac` is active;
- distributed `cluster2D` with `refine=snhc_smpl` and `update_frac`.

Not affected:

- shared memory, where one process sees the whole sample;
- `prob` and `prob_snhc`, where `prob_align2D` samples the whole project on
  the master and workers call `sample4update_reprod`;
- the `startit = 1` sticky path, which reproduces the sample and does not
  read carry-over;
- 3D, whose fractions are computed on the master (Section 3.10).

Consequence: shared memory is the reference behaviour. For affected runs the
current distributed output is not a valid baseline, and Changes 3 and 4 alter
distributed results on purpose.

### 3.3 The blends lose mass when the represented population changes (2D and 3D)

Both blends use the same realized fraction per group (2D class, 3D state),
computed after sampling from active rows with `updatecnt > 0`:

```text
N = active rows of the group with updatecnt > 0
n = those carrying the current sampled marker
f = n / N
```

- **2D:** `new = current + (1 - f) * previous`.
- **3D:** `new = (u/f) * current + (1 - u) * chain`, with `u = f` unless
  `ufrac_trec` overrides it (single-state only).

With the default `u = f`, the 3D current scale is 1 and the chain weight is
`1 - f`, so the 3D blend is exactly the 2D recurrence. The mass argument in
`blend_trailing_accumulators` (`simple_commanders_rec_distr.f90`),
`(u/f)*(f*D) + (1-u)*D = D`, assumes that the chain already represents the
same `N` particles the fraction is measured against. That fails whenever
the represented population changes between blends:

- **First-time particles.** Sampling increments `updatecnt` before the
  fraction is computed, so a particle sampled for the first time counts in
  `n` and `N` although the stored sums contain nothing from it.
- **Deactivation** removes rows from `N` while their contributions stay in
  the stored sums.
- **Re-activation** returns rows to `N` whose contributions were scaled
  away.

Example: a group whose stored sums represent 100 particles. This round
resamples 10 of them and samples 30 first-time particles, so `N = 130` and
`n = 40`.

| Rule | Previous weight | Mass after blend (population 130) |
| --- | --- | --- |
| Current: `1 - f` | 1 - 40/130 = 0.69 | 40 + 69 = 109 |
| Population rule (Section 4.1): `(N - n)/M` | 90/100 = 0.9 | 40 + 90 = 130 |

[Refine2D Class Averaging](../../algorithms/refine2d_class_averaging.md)
states that the 2D blend keeps the full dataset's sampling mass, and
[Sampling and Fractional Updates](../../algorithms/sampling_and_fractional_updates.md)
states it for the 3D chain. Both statements hold only when the represented
population is unchanged.

Where the population changes:

- **2D:** every streaming-pool iteration (appends, and class selections in
  `update_match_class_states_in_pool` switch rows on and off), and the early
  abinitio2D iterations (`sample4update_cnt` favours low `updatecnt`).
- **3D:** trailing iterations before sampling has covered every active
  particle, and any cleanup that deactivates rows. Appended rows change the
  row count, which already discards and re-seeds the chain (Section 3.12).

The size of the 3D effect depends on how many first-time particles trailing
stages see. Measuring it is part of the validation (Section 10).

### 3.4 Class centering shifts the 2D carry-over

`simple_matcher_pftc_prep.f90` calls `cavger_shift_partial_eosum` for each
class whose centering offset exceeds `CENTHRESH`. Each worker shifts its own
previous sums by an offset computed from the shared merged reference. Once
the owner holds the only previous set, it must apply that shift once, before
blending.

### 3.5 The class update fraction is aggregate policy

`apply_weights2cavgs` applies one scalar weight to all four arrays of a
class. This is an intentionally class-local stochastic update. It does not
identify which old Fourier contributions belong to which particles, and it
cannot move old mass when particles change class or leave the active set.
This refactor keeps that structure; only the scalar changes (Section 4.1).

### 3.6 Streaming `average_into` is dead

`stream_chunk%read` (`src/main/stream/simple_stream_chunk.f90`) sums each
chunk's part files, divides by `nparts_chunk`, and writes partless files.
Nothing calls the type-bound procedure, and nothing reads its four outputs;
the `chunks(...)%read_*` calls in `simple_projfile_utils.f90` are
`sp_project` methods. The division affects no live recurrence.

### 3.7 The refine3D continuation count check is invalid

On `continue=yes` with `update_frac` and the gridding backend,
`simple_refine3D_strategy.f90` requires `nparts * 4` files matching
`refine3D_partial_rec_glob` (`*recvol_state*part*`) in the previous
directory. The check is wrong in three ways:

- the pattern matches every state, so `nstates > 1` gives
  `4 * nparts * nstates` files and a valid continuation fails;
- a state with no particles in a part writes no partial (Section 3.8), so the
  count can be short;
- `remove_partial_rec_files` deletes these transient files before every
  round, and nothing on the continuation path reads them.

Remove the check; do not make reconstruct3D partials durable to satisfy it.

### 3.8 Missing gridding partials are zero by design

`simple_matcher_3Drec` writes no partial for a state with no particles in a
worker part, and writes both halves for a populated state even when one half
is empty. An absent state/part pair in `read_gridding_pair_accumulators`
therefore means a zero contribution. No reconstruct3D manifest or
explicit-empty marker is needed.

### 3.9 The PCG raw payload needs no change

The PCG raw `(B, D)` payload is versioned, published through a temporary
file, and reduced in ascending order. Only the PCG trailing chain's blend
weights change (Change 5).

### 3.10 3D fractions are computed on the master, from labels after the search

`get_state_update_fracs` is called in `exec_volassemble`
(`determine_trailing_update_fraction`, `simple_commanders_rec_distr.f90`) and
in the distributed PCG master (`simple_rec3D_pcg_strategy.f90`), after the
merge. Section 3.2 does not apply. The master sees the state labels after
the search, which the population rule handles without the labels from
before the search (Section 4.1). The PCG master also checks that its raw
particle counts match the gridding sampling bookkeeping
(`count_state_sampling`); that check must keep holding.

### 3.11 3D fractional update without trailing carries nothing

With `update_frac` and `trail_rec=no`, `exec_volassemble` reconstructs from
the current partials only. The `l_update_frac` branch in
`read_gridding_pair_accumulators` only zero-pads a smaller previous grid when
the trailing chain is read. There is no carry-over and nothing to change.

### 3.12 3D chain seeding and identity

- Stage-boundary seeds (`trail_seed`) come from a full reconstruction whose
  particles `sample4rec` selects: active rows with `updatecnt > 0`, which is
  `N` at the seed.
- The bootstrap seed writes the current partials scaled by `1/f`, which is
  also mass `N`.
- The gridding chain manifest (`write_trail_chain_set`,
  `validate_trail_chain`) records box, sampling, total row count, state
  layout, generation and component sizes, but not the population the chain
  represents. A row-count change discards and re-seeds the chain.
- The PCG chain carries its own provenance (`pcg_chain_provenance`).

Change 5 adds the represented population per state to both.

### 3.13 `sample4rec` decides coverage per worker range (latent)

`sample4rec` (`simple_oris_sampling.f90`) reconstructs active rows with
`updatecnt > 0` if its range contains any. Otherwise it reconstructs every
active row in the range. Workers call it with their own `[fromp, top]`
(`simple_rec3D_strategy.f90`, `simple_commanders_rec.f90`). A worker whose
range holds no updated row, while other ranges do, would therefore add
never-aligned particles to the reconstruction. Random sampling makes this
rare, but it is the same class of error as Section 3.2: a global decision
taken from a local view. The decision must use the whole project.

## 4. Scientific State Contract

### 4.1 Population rule

For each group `g` (2D class or 3D state), the stored sums record `M(g)`,
the population they represent. The owner reads `N(g)` and `n(g)` from the
merged project (Section 3.3) and applies

```text
f = n / N                         realized fraction (unchanged definition)
u = f, or ufrac_trec (3D, single state)
s = u / f                         current scale; 1 by default
w = (1 - u) * N / M               previous weight
new = s * current + w * shift(previous)
M  <- N                           recorded with the new sums
```

With the default `u = f`, `s = 1` and `w = (N - n) / M`. Here `N - n` is the
number of active, previously updated particles not resampled this round.
Unsampled particles keep their class or state, so this count is the same
under labels taken before or after the search. The mass after the blend is
`n + (N - n) = N`, the represented population, whatever joined or left.
`shift` is the 2D class centering shift (Section 3.4); 3D has none.

Properties:

- First-time particles are in `n` and `N` but not in `M`, so they add mass
  without displacing old mass.
- Deactivated rows leave `N`, so their share of old mass leaves too.
- Resampled particles and particles that moved group are in `n`, so their
  old contributions decay by `w`.
- With no change in population, `w = 1 - f`: the rule reproduces the current
  recurrence exactly.
- `n = 0` gives `w = (N - n)/M`, which is 1 unless rows left.
- `u = 1` replaces the stored sums.
- `w` can exceed 1 when deactivated rows return; see D2.

For the 3D `ufrac_trec` override, the current-map coefficient stays `u` and
the mass stays `N`.

**2D stored meaning.** Four unregularized Fourier accumulators per class:
even and odd numerators, even and odd CTF-squared sums. They are captured
before CTF-density correction, ML-prior attachment, inverse FFT,
low-resolution even/odd insertion, gridding correction or any other
restoration step. Workers and the shared-memory matcher accumulate the
current sample into zeroed arrays. The owner needs nothing from the workers
beyond their current sums and centering offsets; it computes the counts from
the merged project, which cluster2D merges before class-average assembly.

**3D stored meaning.** Unchanged: the even/odd Fourier sums and densities of
the gridding chain, and the raw `(B, D)` accumulators of the PCG chain.

### 4.2 Approximation boundary

The recurrence is not an exact sum over the current class or state
assignments after particles move or leave the active set. The population
rule keeps the mass right, not the membership: old contributions are removed
in proportion, not particle by particle. Exact transport would need
per-particle records with storage and I/O comparable to four particle-image
banks, and is out of scope. The 2D and 3D policy documents and the API
comments must say so. Terms such as "exact class sum" or "mass transport"
must not be used for the aggregate.

## 5. Ownership, Files and Lifecycle

### 5.1 Ownership

- The class-average domain (`simple_classaverager` and its restore
  submodule) owns reading, shifting, blending and writing the 2D previous
  set. In distributed runs the owner is the `cavgassemble` step
  (`cavger_assemble_sums_from_parts`); in shared memory it is the matcher's
  restoration finalisation.
- `exec_volassemble` owns the gridding chain blend, and the distributed PCG
  master owns the PCG chain blend, as today.
- The counts `N(g)` and `n(g)` come from one shared helper in `ori`, used by
  both 2D and 3D, which keeps `get_class_update_fracs` and
  `get_state_update_fracs` consistent.
- Workers never read carried state.

### 5.2 Files

- **2D previous set:** fixed names without a part number, in the run
  directory, as the current part files are today. It holds the four arrays
  and, per class, `M(c)`. The four arrays are published together: either
  one file with a small header and four sections, or four files carrying
  the same generation number, which the reader checks. Every write goes to a
  temporary name and is renamed into place.
- **2D current contributions:** one per worker per iteration, named
  distinctly from the legacy `cavgs_*_partN` and `ctfsqsums_*_partN` names
  so that legacy carry-over files are never mistaken for current
  contributions. They carry the worker's centering offsets. The master
  deletes them before launching workers and after the blend, like assignment
  files.
- **3D chains:** the gridding manifest and the PCG chain provenance gain
  `M(s)` and a format version. A chain without it (written by an older
  build) is discarded and re-seeded by the existing path.

### 5.3 2D header check and full update

The 2D previous set records `box_crop`, `smpd_crop`, class count, kinds and a
format version. If the set is missing, unreadable, internally inconsistent,
or disagrees with the run on any of these, the owner runs that iteration as a
full update: all particles, no carry-over, the path `startit = 1` already
takes. Neither pool class count nor abinitio2D `box_crop` changes within a
run, so this check adds no full updates to current workflows.

### 5.4 2D lifecycle

| Event | Behaviour |
| --- | --- |
| `nparts`, `numlen` or shared/distributed change | Accept; nothing depends on the layout. |
| Particles appended (pool) | Accept; first-time particles join `N` without displacing old mass. |
| Activation or deactivation (pool class selection) | Accept; the population rule removes or restores mass in proportion. |
| Class moves inside an update | Accept, as today. |
| Pool crop-box upsample | The owner pads the one previous set once (`cavger_pad_partial_sums` on one set, not per part). |
| Class count or any other geometry change | Full update. |
| Missing, unreadable or inconsistent previous set | Full update. |
| Crash between the state write and the project write | Accept what is on disk. Either side is at most one iteration out of step, and the population rule re-anchors the mass on the next blend. |

## 6. 2D Iteration Protocol

Distributed:

1. The master deletes stale current-contribution files and launches workers.
2. Each worker samples, searches, accumulates its current sums from zero,
   and writes its current-contribution file with its centering offsets.
3. After the normal worker barrier and the project merge, the owner checks
   that every part's file exists, sums the parts in ascending order, checks
   that the workers' centering offsets agree, and computes `N(c)` and `n(c)`
   from the merged project.
4. The owner reads the previous set (or selects a full update), shifts,
   blends with the population rule, writes the new set with `M(c) = N(c)`,
   restores class averages and FRCs, and writes the project.
5. The master deletes the current-contribution files.

Shared memory runs steps 2 to 4 in one process, without contribution files.

## 7. Compatibility and Scope Boundaries

- Legacy `cavgs_*_partN.mrc` and `ctfsqsums_*_partN.mrc` files are never read
  as carry-over. A 2D run continued from an old directory starts with one
  full update. No converter is planned.
- 3D chains from an older build are re-seeded once (Section 5.2).
- No public persistence selector, reconstruct3D manifest, explicit-empty
  markers or PCG raw-payload change is part of this work.
- Iteration class-average stacks, FRCs, JPEGs, restored maps and their names
  and retention are unchanged.
- Sampling code still owns the selected subset. Class-average and volume
  assembly code consume `sampled` and `updatecnt`; they do not choose or
  change the sample.

## 8. Rollout Plan

1. **Independent cleanups (Changes 1, 2 and 6).** Remove the continuation
   count check. Recheck that nothing calls `stream_chunk%read` or reads its
   outputs, then remove it. Make `sample4rec` decide coverage over the whole
   project, with a unit test.
2. **Interim 2D fraction fix (Change 3).** Give `get_class_update_fracs` an
   optional particle range. Distributed workers pass `[fromp, top]`, so each
   worker scales its partition-local sums by its partition-local fraction;
   shared memory is unchanged. Add a `simple_oris_tester` unit test. Save
   this state as a separately shippable patch.
3. **Shared counting helper and population rule.** One `ori` helper returns
   `N(g)` and `n(g)` for class or state labels. A pure blend-weight function
   returns `s` and `w`. Unit tests cover the cases in Section 10.
4. **Owner-side 2D carry-over (Change 4).** Shared memory and distributed
   2D: current-only workers, centering offsets, the previous set with
   `M(c)`, the population rule, cleanup. Remove the interim range argument
   if nothing else uses it.
5. **Streaming pool.** One previous set in the pool directory; the crop-box
   upsample pads it once.
6. **3D chains (Change 5).** Gridding and PCG chain blends use the population
   rule; the manifest and provenance record `M(s)`; seeds record `M(s) = N`.
7. **Cleanup and documents.** Remove the part-file carry-over reads, the
   per-part padding loop and the cleanup branches that preserve part files.
   Update the documents listed under Cutover audit in Section 10.

## 9. Decisions

- **D1, accepted 2026-09-28:** keep the class- and state-local aggregate
  recurrence; exact per-particle contributions are not required.
- **D2, accepted 2026-09-28:** adopt the population rule of Section 4.1 for
  both 2D and 3D. The realized fraction `f = n/N` keeps its documented
  definition, and each stored set records the population it represents. The
  rule changes results wherever the represented population changes
  (Section 3.3) and nowhere else. It supersedes the earlier proposal to
  count only particles updated before this round, which fixed first-time
  particles but not deactivation, and needed the labels from before the
  search. `w > 1` is allowed when deactivated rows return, since that keeps
  the mass equal to the represented population.
- **D3, accepted 2026-09-28:** shared memory is the 2D reference behaviour.
  Distributed results in the paths listed in Section 3.2 change on purpose.
- **D4, decided 2026-09-28 (implementation):** one file. The previous set is
  `cavg_state.bin` and each worker contribution `cavg_contrib_part<N>.bin`
  (`simple_cavg_sums`): a header (magic, format version, kind, class count,
  `box_crop`, `smpd_crop`, array shape), the per-class metadata (`M(c)`, or the
  centering offsets and accumulated populations), the four arrays and a
  closing byte count, written to `<name>.tmp` and renamed into place. Reason:
  one rename publishes the four arrays together, so a reader can never see a
  mixed set and no cross-file generation check is needed; a truncated or
  foreign file fails the size and magic checks and selects a full update.
- **Baselines:** before Changes 4 and 5, record pre-change runs (restored
  outputs, populations, FRC/FSC, and the carried mass per represented
  particle) as comparison targets.

## 10. Validation Matrix

### Cleanups

- `continue=yes` with fractional gridding refinement succeeds for
  `nstates = 1` and `nstates > 1`.
- Source and generated-index searches find no call to `stream_chunk%read`
  and no reader of its partless outputs; chunk memoization and pool import
  still use the converged project handoff.
- `sample4rec`: a unit test where only one range holds updated rows returns
  no never-updated row for any range.

### Interim fraction fix

- Unit test: in a two-partition fixture where only one range carries the new
  marker, the range-restricted fraction equals that partition's own
  fraction, and the whole-project call reproduces the old diluted value.

### Population rule (unit level, 2D and 3D)

- Mass conservation over a scripted sequence of blends with unit
  per-particle mass: first-time particles, resampling, group moves,
  deactivation, re-activation and appends. The carried mass equals `N(g)`
  after every blend, and the old rule's drift is reproduced as a control.
- Without a population change, the rule equals the current recurrence to
  single-precision tolerance.
- `f = 0`, `f = 1`, `n = 0`, `M = 0`, empty groups, and `ufrac_trec`, where
  the current-map coefficient stays `u`.
- 2D owner blend: for given sums, offsets and counts, splitting the current
  sums into 1, 2 or 4 parts gives the same result up to summation order, and
  the centering shift is applied once.
- 2D previous-set I/O: round trip, header refusals, a mixed set and a
  missing set each select a full update, and a write is never observed
  half-done.
- 3D chain manifest and PCG provenance: `M(s)` round trip; a chain without
  it is re-seeded.
- The existing trailing-reconstruction blend sub-suite
  (`test=unit_image suite=trailing_reconstruction_blend`) is extended with
  population changes.

### Workflow level

- Fast gate, the affected `unit_*` areas, `lib_reconstruction`,
  `lib_stream`, `pcg_recon`, `simulated_workflow_1jxy`,
  `simulated_workflow_6vxx`, `mini_stream_1jxy`, `mini_stream_6vxx` and the
  `abinitio3D_addon` gate pass, apart from failures shown to be pre-existing
  on the baseline build.
- Real data, before against after:
  - 2D: shared memory against distributed `nparts=4`, including an
    `update_frac` run with `refine=snhc_smpl`;
  - 3D: a single-state trailing run and a three-state run;
  - the add-on cohort chain.
- Real-data metrics:
  - carried mass per represented particle, per iteration: stable after the
    change; before it, distributed 2D drifts and populations that change
    lose mass;
  - first-time share of each 3D trailing sample;
  - populations and FRC/FSC resolutions;
  - map and class-average agreement.
- Changing `nparts` or execution mode across a 2D `continue=yes` neither
  crashes nor loses or multiplies previous mass.
- The pool upsample pads one set and agrees with the old per-part padding
  summed over parts.
- A run killed between the state write and the project write restarts and
  completes.

### Cutover audit

- No worker reads carried state, and no carry-over path reads `*_partN`
  class-average files.
- These documents describe the new contract:
  - `abinitio2D_policy.md`: the artifact list, the carry-over paragraph and
    the checklist bullet on preserving partial sums;
  - `importance_sampling_fractional_update_policy.md`: the
    `get_class_update_fracs` text, the distributed-cleanup paragraph, and
    the 3D trailing contribution and chain weights;
  - `cluster2d_class_averaging.md` and `sampling_and_fractional_updates.md`:
    the blend formulas and the mass statements;
  - the comments in `blend_trailing_accumulators` and the PCG chain blend;
  - the generated code maps.
- The `simple-frac-update-trailing` skill states the chain weights. Its
  update is proposed for approval, not applied.

## 11. Acceptance Criteria

- The `nparts * 4` reconstruct3D continuation check is gone, with no broader
  reconstruct3D change.
- The dead `stream_chunk%read` path is removed and not replaced by another
  chunk-part reconciliation file.
- `sample4rec` decides coverage from the whole project.
- Distributed 2D results follow the shared-memory recurrence.
- One partless 2D previous set is the only durable source of class-average
  carry-over, and its four arrays are published together.
- Carried 2D and 3D sums keep the mass of the population they represent
  through appends, first-time samples, deactivation and re-activation.
- Changing partition count or execution mode neither invalidates nor scales
  carry-over.
- The approximation boundary and the lifecycle rules are documented in the
  owning policies.
- Legacy part-file carry-over reads, per-part padding and the cleanup
  branches that preserve part files are gone.

## 12. Validation Record

- 2026-09-05: initial source review confirmed partition-shaped class-average
  carry-over and the invalid gridding continuation check. Missing gridding
  state/part files were confirmed to be zero contributions by design, so the
  proposed reconstruct3D transaction layer was dropped. No implementation,
  compilation or runtime validation.
- 2026-09-11: planning review after the canonical sigma2 cutover. The plan
  moved to one container, integrity and generation checks, project and
  checkpoint identity, and explicit lifecycle rules. It prohibited a public
  rollout selector and made rebuild and streaming ownership
  pre-implementation gates. Source inspection found no caller of
  `stream_chunk%read`, so item 2 became dead-code removal.
- 2026-09-28: source review of 2D against `master`, by reading code only.
  Findings:
  - distributed workers dilute the class update fraction to about `f/nparts`
    whenever they draw their own sample (Section 3.2);
  - first-time particles inflate the fraction against the documented mass
    statement (Section 3.3);
  - the carry-over is shifted by class centering (Section 3.4);
  - the continuation count check also fails every multi-state continuation
    (Section 3.7);
  - the `nparts` persistence defect is latent in current workflows;
  - pool class count and abinitio2D `box_crop` are fixed within a run.

  The plan was reduced to an owner-side blend. Per-particle digests,
  generation-scoped payloads, project registration, directory sync and
  rebuild-on-activation-change were withdrawn, and an interim fraction fix
  was added.
- 2026-09-28: source review of the 3D trailing chains, by reading code only.
  Findings:
  - with `u = f` the 3D chain blend is the 2D recurrence, so both lose mass
    when the represented population changes (Section 3.3);
  - 3D fractions are computed on the master (Section 3.10);
  - trailing-free fractional 3D carries nothing (Section 3.11);
  - chain seeds represent `N`, but chains do not record it (Section 3.12);
  - `sample4rec` decides coverage per worker range (Section 3.13).

  The first-time fix proposed earlier the same day was replaced by the
  population rule (Section 4.1), which also covers deactivation and needs no
  counts from workers. Changes 5 and 6 were added, and D2 was restated.
- 2026-09-28: Hans accepted D2, with `w > 1` allowed.
- 2026-09-28: implementation of Changes 1 to 6 (work tree on top of
  `75ffa8fef`). Code:
  - Change 1: the `nparts * 4` count check is gone from
    `simple_refine3D_strategy.f90`; `refine3D_partial_rec_glob` lost its only
    caller and was deleted.
  - Change 2: source and generated-index searches found no caller of
    `stream_chunk%read` and no reader of its partless outputs; the routine,
    its binding and the import it alone needed are deleted.
  - Change 6: `sample4rec` decides "any active row updated" over the whole
    project; `get_state_rec_pops` gives the per-state population that rule
    reconstructs (the stage-boundary seed's `M`).
  - Change 3 was implemented as specified (`get_class_update_fracs` with a
    worker range, unit test) and shipped as a separate patch, then removed
    with Change 4, which left `get_class_update_fracs` without a caller.
  - Section 5.1 helper: `oris%get_group_update_counts(label, ngroups, N, n)`;
    the pure weight function is `population_blend_weights`, first in a module
    of its own and, since 2026-09-29, part of `simple_oris` (see below).
  - Change 4: workers and the shared-memory matcher accumulate from zero;
    workers write contributions; the owner (`cavger_commit_carryover` in
    shared memory, `cavger_assemble_sums_from_parts` in `cavgassemble`) sums,
    checks offsets and carry-over mode, shifts the previous set once, blends
    and publishes. Legacy part-file reads, `apply_weights2cavgs`, the
    per-worker shift and the cleanup branch that kept part files are gone.
  - Change 5 (pool): `cavger_pad_carried_sums` pads the one carried set.
  - Change 5 (3D): gridding manifest version 2 records `M(s)`
    (`simple_trail_chain_manifest`); the PCG chain identity is
    `pcgtrail-v3`, whose header particle count per half is its represented
    population. Both blends use state-level population-rule weights; seeds
    record `M = N` (bootstrap and replacement) or the `sample4rec` population
    (stage-boundary seed).
  - Logging: one `>>> CAVG CARRY-OVER` line per 2D iteration (N, n, M before
    and after, w average/min/max, s, carried mass per represented particle);
    one `TRAILING ACCUMULATOR BLEND` (gridding) or `PCG TRAILING BLEND` line
    per state and iteration (w, s, the former weight 1 - u, N, n, M before
    and after, first-time rows of the sample).

  Findings and choices where the plan was silent or inexact:
  - The represented population recorded after a blend is `s*n + w*M`. It is
    `N` whenever `M > 0` or `n = N`; with no previous mass (`M = 0`) it is the
    mass actually stored, `n`, not `N`. An iteration without carry-over
    records the accumulated population of each class.
  - A consequence of the population rule, not a defect: when a carried set
    holds fewer particles than the rows it must stand for, `w` exceeds 1. The
    first carry-over iteration of `abinitio2D` is the common case: stage 1
    restores from its sticky subset only (`M = n`), and stage 2 then scales
    that set to stand for all previously updated, unsampled rows
    (`w = (N - n)/M`), where the former rule used `1 - f`.
  - Section 5.3 says the owner runs a full update when the previous set is
    unusable, but by the time the owner runs, the workers have already
    sampled. The decision is therefore taken by the master before the
    iteration (`carryover_needs_full_update` in the cluster2D strategies:
    `update_frac = 1` for that iteration, in the job description and the
    `prob_align2D` command line); the owner only warns if the set vanished in
    between and restores from the current sample.
  - The restored averages use the represented even/odd populations (active,
    updated rows of each class) when the owner blends, and the accumulated
    populations otherwise. Before, shared memory used the accumulated sample
    and distributed runs all active rows, so a class with carried sums but no
    sampled particle was zeroed in shared memory only.
  - The class-centering offset is computed from the merged reference and the
    project before the search, so workers agree; the owner warns and applies
    part 1's offsets if they differ by more than 1e-3 pixel.
  - The note's Section 3.3 example rounds: the former rule gives
    40 + (90/130)*100 = 109.23, not 109.
  - Two private routines of the restore submodule were already without a
    caller (`calc_class_center_shift`, `shift_stack_slice2D`) and were deleted
    with the code they duplicated.
  - The tracked `doc/code_overview/fortran-indexes` carry absolute paths of
    the owner's checkout; regenerating them on another machine rewrites every
    path, so they were not regenerated here (they list the deleted
    `refine3D_partial_rec_glob` and `average_into`).
  - `scripts/check_descr.py` reports `simple_refine_motion_model_strategy.f90`
    on the base tree too (pre-existing).
  - Not done: the killed-between-writes restart of Section 10 was not run as a
    workflow; the one-file publication through a rename and the refusal of a
    truncated or temporary set are unit-tested instead.
- 2026-09-29: rebased onto `5fa147fc7` (15 upstream commits) and retested on
  the Dell against a `5fa147fc7` baseline build. The diff applied cleanly: only
  two documents overlapped, and no upstream code calls a removed routine.
  - Unit tests: fast gate 13/13, with the new sub-suites run and passing
    (update blend 69, class-average carry-over 38, trailing chain identity 17
    and trailing-reconstruction blend 39 checks).
  - Library and workflow gates: `lib_reconstruction`, `lib_stream`,
    `pcg_recon`, `simulated_workflow_1jxy`, `simulated_workflow_6vxx`,
    `mini_stream_1jxy` and `abinitio3D_addon` pass. `mini_stream_6vxx` fails
    its picking check ("fewer than 40 percent of placed particles were
    picked") on the baseline too: pre-existing, and it fails before any changed
    code runs.
  - 2D, HolJunk base set (33016 particles, 100 classes, sampled fraction 0.24
    to 0.25). Carried CTF-squared mass over ten cluster2D continuation
    iterations:
    - baseline, shared memory: flat;
    - baseline, distributed: 1.368e13 to 2.747e13, +137 % from the abinitio2D
      end, as Section 3.2 predicts;
    - new code: flat (ratio 1.0000) in shared memory, distributed, after a
      `continue=yes` from 4 to 3 parts, and after switching from shared memory
      to distributed.

    Final class statistics are within run-to-run spread: population-weighted
    mean FRC0.143 21.0 to 22.2 A, medians equal or one shell apart, and the
    direction of the small differences reversed from the 2026-09-28 runs.
  - 3D, bgal single state (5513 particles, fraction 0.25): under the old rule
    the carried mass per represented particle fell to 0.80 when first-time
    particles joined; the new code keeps it at 1.0000. Final FSC 0.5/0.143 is
    4.41/3.84 A on both builds. The add-on reports are identical on both builds:
    3.75 A, +6 shells, IMPROVED, correlation 0.982 with the base map.
  - 3D, HolJunk three states:
    - abinitio3D, which does not trail (see the finding below): FSC0.143 21 to
      24 A on both builds;
    - trailing test: a ten-iteration refine3D continuation (`trail_rec=yes
      update_frac=0.3 refine=prob_state lp=10`) run by both builds from the
      baseline's abinitio3D result, so both start from the same input;
    - 39 % of each sample were first-time particles;
    - old rule: chain weights 0.50 to 0.71, and carried mass per particle down
      to 0.78;
    - new rule: `w` 0.70 to 1.00, and carried mass per particle 1.0000;
    - gridding, final FSC 0.5/0.143: baseline 27.6/10.0, 27.6/10.0, 10.2/9.7 A;
      new code 27.6/10.0, 10.2/10.0, 16.3/10.0 A. FSC0.143 is capped by
      `lp=10`, and there is no consistent difference;
    - PCG, final FSC0.143: baseline 25.6, 22.4, 21.1 A; new code 17.1, 10.0,
      10.0 A. The gain holds in all three states but was measured once per
      condition.
  - Finding for the owner: multi-state abinitio3D never trails.
    - With `nstates > 1`, abinitio3D runs `multivol_mode=independent` with
      `NSTAGES_INDEPENDENT = 5` stages, while `set_refine3D_trailrec_policy`
      enables trailing from `TRAILREC_STAGE_MULTI = NSTAGES = 8`.
    - The multi-state add-on runs of 2026-09-27 logged no chain blend either.
    - For several states, the population rule therefore acts only in refine3D
      with `trail_rec=yes` and in `multivol_mode=docked`.
  - Still not run: the killed-between-writes restart of Section 10.
- 2026-09-29, after review: `population_blend_weights` moved from its own
  module (`simple_update_blend`) into the `oris` class module. It is declared
  public in `simple_oris` beside the type and implemented in the
  `simple_oris_sampling` submodule, next to the `sampled`/`updatecnt`
  bookkeeping whose counts it consumes. It stays a free elemental procedure,
  because the gridding blend knows `M` only inside `restore_state_from_parts`,
  after the counts were taken in `exec_volassemble`. Its tests joined the
  `orientation collection` sub-suite of `unit_ori`, and the separate
  `update blend` sub-suite is gone. The empty submodule `simple_oris_weights`
  (retired particle-weight hooks) was deleted. A bounds-checked build exposed
  an out-of-range read in the control case of `test_trail_rec_population`
  (the loop index was used after the loop); the helper now takes the counts as
  arguments, and the control is pinned to its exact value.
