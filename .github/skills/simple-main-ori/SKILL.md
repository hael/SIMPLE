---
name: simple-main-ori
description: Use when working in SIMPLE's src/main/ori subsystem, including orientation records (ori) and collections (oris), the per-particle named slots and how keys are typed and saved, the Euler-angle, shift, label and even/odd conventions, point-group symmetry and the asymmetric unit, sampled/updatecnt bookkeeping and the sample draws, transforms, distances and their units, text and binary record I/O, statistics, neighbourhoods and the population weights of trailing reconstruction.
---

# `src/main/ori`: orientation records, collections and point groups

Owns `ori` (a record), `oris` (a collection: each project segment, `build%eulspace`), `sym` (point
groups, asymmetric unit, reference spiral), the Euler utilities, the `key=value` parser
(`simple_sauron`) and the mechanics of the sample draws. Particle slots: `src/defs/simple_defs_ori.f90`.
It does not decide which particles a pass samples (`strategies/search/simple_matcher_smpl_and_lplims.f90`),
the even/odd halves (`project/simple_eo_partition.f90`) or project I/O (`fileio/simple_binoris.f90`).

## Rules

Silent ones first: nothing reports these.

- **Particle records save only their slots.** A particle-typed record (`is_ptcl=.true.`: `os_ptcl2D`,
  `os_ptcl3D`) keeps named keys in real slots and others in hash tables, but `binoris` writes particle
  segments as the slots alone: a non-slot key (`z`, `mirr`, a new name) would be lost at the next
  project write, so a Debug build stops there (`check_particle_slots` in `simple_binoris`); producers strip
  what they do not keep (`strip_nonslot_keys`; `nonslot_key` names the first such key) and a key read back
  after a write is a slot (`frame`, slot 56). Other segments are hash-typed, saved as text. `=` into an existing record keeps the target's
  typing. A new particle field is a slot (`references/particle_slots.md`); never renumber or reuse one.
- **Absent reads as zero; presence depends on the typing.** `get`/`get_int` of a missing key give 0
  (permanent, by ruling: no stop, no `get_or`), `get_str` ''. On a particle record the slots listed
  in `oriparam_isthere` (class, proj, eo, corr, dfx, dfy, stkind, pind, ...) report absent at 0; the
  others (state, e1-e3, x, y, sampled, updatecnt, ...) are always present. So eo=0 (even) reads as
  absent: on a particle field `get_neven` counts the even rows only once a row is odd, and
  `get_nevenodd() == 0` means "no halves yet"; test one row with `get_eo`, and parse into a
  hash-typed record where presence matters (`read_ctfparams_state_eo`). `oris%isthere` is true if
  any record has the key.
- **Values are reals typed by their text.** Integers are stored as reals (exact to 2^24). Text that
  reads as a number parses as a real: `set(key, '12')` comes back a real after a write and read.
- **Euler angles** are ZYZ in degrees (`e1, e2, e3` = RELION rot, tilt, psi), `R = Rz(e3) Ry(e2)
  Rz(e1)` with the passive `rotmat` (`euler2m`). Compose with `compose`/`euler_compose` (`R2 R1`), never
  by adding angles. `set_euler` goes through the matrix: e1, e3 in [0,360), e2 in [0,180], and within
  0.026 degrees of a pole e3 is folded into e1 and set to 0. A 2D record's in-plane angle is `e3` with
  e1 = e2 = 0: set it with `e3set`, never `set_euler`.
- **Distance units differ.** `.euldist.`, `.inplrotdist.`, `geodesic_dist_trace`, `oris%euldist`:
  radians. `min_euldist`, `find_angres`, the angle of `nearest_proj_neighbors`: degrees.
- **Shifts** `x, y` are pixels of the particle's stack in the particle frame; an image is centred by
  shifting it by minus `[x, y]`. Cropping and resampling: `references/conventions.md`.
- **Labels.** `state` 0 excludes a row (`ptcl2D`'s selects, `ptcl3D`'s is the 3D label); `class` and
  `proj` (a direction of `eulspace`) are unassigned at 0; `eo` is 0 even, 1 odd, -1 none (`reject`); a half sum
  takes it through `get_half`, which stops on -1 (class averages, 3D reconstruction). `get_pop(s)` and
  `get_pinds` skip state 0; `get_n(label)` (largest label, at least 1) and `get_noris` count every row.
- **Halves** of particle fields have one maker, `ensure_field_halves`, with one label in `os_ptcl2D`
  and `os_ptcl3D`; `partition_eo` (row parity) is for class-average fields and time series.
- **Sampling.** Each drawing `sample4update_*` (all but `_reprod`) adds 1 to `updatecnt` of its rows
  and stamps `sampled` with the round: the largest `sampled` among active rows, +1 when
  `incr_sampled`. `_reprod` changes nothing and returns the latest round's rows. The round and
  `sample4rec`'s coverage read the whole field, not `fromto`: give them the full field. The project
  owner draws and writes; a worker only reproduces.
- **Reallocating** makes row indices stale: the calls in `references/conventions.md`.
- **Symmetry.** `sym%apply(o, isym, osym)` gives `R_o R_isym`, other keys kept (op 1 is the identity).
  The asymmetric unit (`get_eullims`) includes the mirror half. Icosahedral operators have one setting,
  SPIDER's (a 5-fold axis on z), everywhere: the builder's `pgrpsyms`, the symmetry search and analysis of
  `volume/`, `symmetrize_map`; the asymmetric unit assumes it. An unsupported group
  calls `stop`: check with `is_valid_pointgroup`.
- **OpenMP**: one row per thread; read strings with `get_static` (`get_str` allocates).

## Where to start

`simple_ori.f90`, `simple_oris.f90` (interface; bodies in its submodules), `simple_sym.f90`,
`simple_ori_utils.f90`. Stack indexing: `doc/policies/project_stack_indexing_policy.md`. Conventions
against RELION: the table in `doc/implementation_notes/completed/ghostbuster_operators.md` 5.2.

## Tests

- `build/bin/simple_test_exec test=unit_ori` (about 1.5 s; `suite=` orientation,
  orientation_collection, symmetry, euler_shift).
- Slots: `test=unit_project suite=binoris`; halves: `suite=eo_partition`; pose frame, hand and
  point-group conventions: `test=unit_numerics suite=truth_metrics`; sampling: `test=pcg_frac_update`.
- Pose or label behaviour: `test=simulated_workflow` (map correlation is a robustness profile: one run
  shows breakage, not improvement) and `test=multistate_workflows suite=states`.

## See also

`doc/policies/importance_sampling_fractional_update_policy.md`; skills `simple-main-project`,
`simple-solve3d-importance-sampling`, `simple-frac-update-trailing`, `simple-main-star`.
