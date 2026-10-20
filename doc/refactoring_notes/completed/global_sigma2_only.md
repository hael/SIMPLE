# Global sigma2 only: removing `sigma_est=group`

Status: completed 2026-10-13 by an unattended run on the Oracle Linux test
machine (section 12), uncommitted for review; the record of every phase is
in section 13 and the review report is `global_sigma2_only_report.md`
beside this file. Written as a refactoring plan on 2026-10-12. Rulings of
the same day:

- `sigma2_convert` is retired.
- The renames of section 5 are accepted.
- The incremental store reduction is decided by measurement (section 6).
- The stream is part of this refactoring.

The outline was revised the same day after a review; every claim below was
checked against the code. Rulings on the revision, also of the same day:

- The `solve3D` manifest key `sigma_est` leaves within manifest version 2,
  with no second version bump (section 4).
- Existing global stores with the STAR-seed provenance are refused and
  rebuilt from particle power spectra, like stores written under `group`
  (section 4).
- The science validation's full re-read of the store moves to the tests
  (section 6).
- The streaming commit and the block-copying merges are part of this
  refactoring. The sparse ranges wait for a measurement on the ribosome
  project, which the test machine does not hold (section 6).

## 1. Decision

- `sigma_est=group` is no longer supported. Under it, each particle stack
  had its own noise spectrum per half.
- Every workflow estimates one global noise spectrum per half: the average
  of the per-particle residual spectra. A half is the even or the odd
  half-set.
- The `sigma_est` key goes, with every place that sets, passes or checks it.
- The per-particle residual spectra (records) stay in the canonical sigma2
  store (`sigma2_state.bin`). A fractional update replaces only some
  records, and the global average must reflect every particle's latest
  residual (`doc/refactoring_notes/completed/canonical_sigma2_state_refactoring.md`,
  section 2).
- Release 4 keeps no backwards compatibility except for `.simple` project
  files.

## 2. What group estimation costs today

The choice between two policies is threaded through the code, and several
expensive structures exist only because the spectrum could differ by stack.

- **A noise table per particle.** At load, `read_groups`
  (`simple_euclid_sigma2.f90`) copies each particle's group spectrum into
  `sigma2_noise(shells, fromp:top)`.
  - Under the global policy every column is one of two spectra.
  - At box 336 (168 shells) that is 67 MB per 100,000 particles in every
    process that loads it, and 504 MB at 750,000.
  - Shared-memory runs hold the whole project. FLEX and `cls_expansion`
    widen the range to the whole project even in distributed runs.
  - Reconstruction consumers (`load_sigma2_groups`) then copy the table
    again into `sigma2_part` (`allocate_ptcls`), and never use the copy.
  - Matchers hold the partition's records (`read_part`) next to the table.
    They update those records and write them back, so the records stay.
- **Per-particle indexing in the hot loops.**
  - The polar correlation kernels (`simple_polarft_corr.f90`, 19
    references) and the particle self-terms (`simple_polarft_memo.f90`, 4)
    compute `real(k)/sigma2_noise(k,iptcl)` per shell, per evaluation.
  - The access is contiguous: shells vary fastest for a fixed particle.
    The cost is a division per shell where a precomputed weight would do,
    and the size of the table.
  - `wsqsum_ptcl_now` (`simple_polarft_memo.f90`) recomputes a particle's
    weighted self-term over all shells on every call, because the memoised
    value may be stale relative to the table.
  - The GPU path builds and uploads a weight table `wks(shells, nptcls)`
    per batch.
- **Per-particle spectra in reconstruction and the other consumers.**
  - `gen_fplane4rec` (`simple_image_ctf.f90`) receives the particle's
    column. For every particle it allocates and resamples it to the plane
    grid, rebuilds its shell lookup table, and divides by the spectrum or
    its square root per pixel.
  - The PCG path (`simple_rec3D_pcg_strategy.f90`) resamples every
    particle's column into `sig2(0:R, nptcls)`. In the state-weighted path
    it then divides each column by the particle's weight.
  - Class-average restoration (`simple_classaverager_restore.f90`) passes
    each particle's column to `gen_fplane4rec`.
  - The Cartesian slots copy the column per batch
    (`prep_cart_batch_of_build`, `simple_matcher_ptcl_batch.f90`), and each
    slot of `simple_cartft_calc.f90` allocates and validates its own copy.
  - `cls_expansion` fills a weight table `w(ncoeff, nptcls)` per class.
- **Class-average restoration reads records it never uses.**
  `cavger_read_euclid_sigma2` calls `read_part`, which reads the partition's
  per-particle records. Restoration only needs the spectra.
- **The store carries a grouping policy.**
  - The header records the policy and the group count. The grouped section
    holds 2 × the number of groups spectra.
  - Every reduction, validation and commit takes an array of group
    identifiers built from `stkind` over all particles.
- **The policy check exists in three copies.** Each validates the store,
  takes the layout digest and, under `group`, recounts stacks from
  `stkind`:
  - `canonical_sigma2_consumable` (`simple_sigma2_files.f90`);
  - the reconstruction start (`simple_rec3D_strategy.f90`);
  - the `solve2D` checkpoint resume (`simple_commanders_solve2D.f90`).

  The chunk merge (`concatenate_canonical_states` in
  `merge_chunk_projfiles`) and the `merge_projects` merge
  (`merge_selected_canonical_field`), both in `simple_projfile_utils.f90`,
  build the group arrays the same way. The chunk merge also serves project
  and class-average commanders and the sieve, not only the stream.
- **Group-only failure modes and workarounds.**
  - A selection that empties one half of a stack makes the reduction refuse
    ("empty even/odd half"). Under the global policy only a selection that
    empties one half of the whole project does.
  - FLEX carries a fallback that detects this and switches to global
    (`ensure_canonical_sigma_state`).
  - `cls_expansion` probes both policies and flips `params%sigma_est` to
    whichever store it finds (`prepare_flex_sigma2`).
- **Command-line plumbing.**
  - About 35 lines set `sigma_est=global` on a child command line, delete
    it, require it or warn about it.
  - Seven programs show it in their UI.
  - The `solve3D` manifest lists it as a recorded and replayed key.

## 3. Target design

- **In memory** (`euclid_sigma2`):
  - two global spectra, `sigma2(shells, 2)`, one per half;
  - the derived arrays each consumer needs, per half and computed once:
    `k/sigma2` for the polar kernels, `1/sqrt(sigma2)` for the Cartesian
    slots, `1/sigma2` for `cls_expansion`. Reconstruction keeps its own
    copies on the plane grid (section 6). Each is a few hundred floats, so
    they stay separate instead of being derived from one another in the
    kernels;
  - the per-particle records of the pass's partition, only for a pass that
    writes residuals (the matchers).
- **Consumers** address the spectrum or a derived array by the particle's
  half.
  - The polar calculator already knows each particle's half (`iseven`).
    The other consumers read it from the orientation (`eo`).
  - The reduction already requires every active particle's half to be 0 or
    1, and the lookup refuses any other value. Today `read_groups` would
    index out of bounds for an active particle without a half, while
    class-average restoration files such a particle under even.
- **On disk.**
  - The store keeps its version-two layout. The grouping word must be
    global and the group count 1, and the grouped section holds the two
    global spectra.
  - Existing global stores, the default everywhere, stay valid, except
    those whose provenance is the STAR seed (section 4).
  - A store written under `group` or with the STAR-seed provenance fails
    validation and is rebuilt from particle power spectra by the one
    rebuild entry below. Today that rebuild does not run everywhere.
- **One rebuild entry and one identity check.**
  - `canonical_sigma2_consumable(project, os, box, smpd, message)` stays
    the identity check.
  - The rebuild moves into one routine in `simple_sigma2_bootstrap.f90`,
    the module that already owns the bootstrap.
    `ensure_sigma2_for_iteration` is the starting point; today it reads
    only the `ptcl3D` field.
  - Today each workflow rolls its own rebuild: refine3D's initialisation,
    the reconstruction start, the final reconstruction, the `solve2D`
    checkpoint resume, and `ensure_sigma2_for_iteration` for `solve3D` and
    FLEX. Three paths read or update the store with no rebuild in front:
    - `refine2D` checks nothing. Its first `prepare_canonical_sigma_update`
      validates the store and stops on failure.
    - `make_cavgs` opens the store through `cavger_read_euclid_sigma2`
      without a check.
    - `load_sigma2_groups` checks and stops on failure, without
      rebuilding. It stays a hard check behind an entry that has rebuilt.
  - Every workflow entry that consumes or updates the store calls the
    central routine. Readers inside a pass (the matcher load, class-average
    restoration in a matcher) stay unchecked, because their master has run
    it.
  - This lands before the store refuses the stack policy (section 7).
  - The final reconstruction is the exception (decided in Phase 1). Its
    rebuild is `bootstrap_rec3D`, which seeds from particle power spectra
    and then runs a residual pass at the final sampling. A seed alone
    would skip that pass, so it keeps its own decision between
    `reconstruct3D` and `bootstrap_rec3D` on the one identity check.
- **The power-spectrum bootstrap.**
  - The bootstrap (`calc_pspec`) already samples globally: at most 25,000
    particles, split between the halves and drawn at random within each.
  - Its workers return per-half sums instead of a per-particle file. The
    normalisation stays exactly as it is (section 6).

## 4. What goes

### Parameters, UI and command lines

- **Parameters.** `sigma_est` and `l_sigma_glob` in `simple_parameters.f90`,
  their registration in `simple_parameters_parse.f90`, and their derivation
  in `simple_parameters_phases.f90`.
- **The UI parameter** `sigma_est` in `simple_ui_params_common.f90`, and its
  use by these programs:
  - `refine3D` and `refine3D_auto` (`simple_ui_refine3D.f90`);
  - `refine3D_nano` (`single_ui_nano3D.f90`);
  - `flex_pca` and `refine3D_states` (`simple_ui_heterogeneity.f90`);
  - `sigma2_convert` (`simple_ui_other.f90`).
- **`solve2D`'s own `sigma_est` input** (`new_solve2D_descriptor` in
  `simple_ui_refine2D.f90`). `solve2D_chunks` has no input, but its
  commander forces the value.
- **Every set, delete, require or warn of `sigma_est`:**
  - `simple_commanders_refine3D.f90`, eight sites, including the
    `bootstrap_rec3D` warning;
  - `simple_commanders_solve2D.f90`, three sites;
  - `simple_commanders_solve3D.f90`, two sites;
  - `simple_commanders_validate.f90`, two sites;
  - `single_commanders_nano3D.f90`;
  - `simple_external_reference_pose_initialization.f90`;
  - `simple_sigma2_bootstrap.f90`;
  - `simple_final_rec.f90`;
  - `simple_refine3D_strategy.f90`;
  - `simple_refine3D_states_controller.f90`;
  - `simple_flex_pca_project_gateway.f90`;
  - `simple_cls_expansion_strategy.f90`.
- **The stream's three sites:**
  - `simple_stream_pool2D.f90` passes it on;
  - `simple_stream_stage_initial_analysis.f90` and
    `simple_commanders_stream_p06_pool2D.f90` set it.
- **The stream's dead `imgkind=sigma2` loop.** The initial analysis copies an
  `imgkind=sigma2` output entry under a balanced name
  (`simple_stream_stage_initial_analysis.f90`, around line 1714). Nothing
  writes such entries any more, so the loop goes.
- **The `solve3D` manifest key.** `sigma_est` comes out of
  `MANIFEST_INPUT_KEYS` and `MANIFEST_REPLAY_KEYS`
  (`simple_solve3D_manifest.f90`), within manifest version 2 (ruled).
  - The manifest records the keys given at entry, before `solve3D` forces
    `sigma_est=global`.
  - `solve3D` does not offer `sigma_est` in its UI. But the parser accepts
    every vocabulary key for every program, so a hand-typed key would be
    recorded.
  - The reader refuses unknown input keys and any other schema version.
    Version 2 (the sigma2 artifact identified by its generation and valid
    file structure) is not yet in use, so the key leaves without a second bump.
    A manifest that recorded the key is refused as carrying an unknown key,
    and that solution can no longer be extended.
- **The memory benchmark script** (`scripts/memory/benchmark_solve2d.py`)
  passes and records `sigma_est`.

### The canonical store

- **`simple_sigma2_state_file.f90`.**
  - `SIGMA2_GROUP_STACK` goes.
  - `validate_header_values` accepts only the global policy with one group,
    and only the power-spectrum and residual provenances.
  - `sigma2_state_init_header` loses its `ngroups` and `grouping` arguments
    and writes the fixed values.
  - `SIGMA2_RANGE_DIRECT` is never written: the header's range-mode word is
    always the local mode. The constant goes, and validation accepts only
    the local mode.
- **`simple_sigma2_state.f90`.**
  - `sigma2_state_reduce_groups`, `sigma2_state_validate_science` and
    `sigma2_state_commit` lose the group identifier array and the policy
    switch.
  - `sigma2_state_validate_identity` loses `expected_grouping` and
    `expected_ngroups`.
- **The rebuild code of each workflow** (refine3D initialisation,
  reconstruction start, final reconstruction, `solve2D` resume) becomes a
  call of the central rebuild entry (section 3).
  `canonical_sigma2_consumable` loses its `l_sigma_glob` argument, and
  `sigma2_estimate_available` and `ensure_sigma2_for_iteration`
  (`simple_sigma2_bootstrap.f90`) lose `l_sigma_glob` and `sigma_est`.
- **`simple_solve3D_utils.f90`.** The unused imports of the two policy
  constants go, and so do the two `l_sigma_glob` passes.
- **The group-count and group-array code** in:
  - `exec_calc_group_sigmas` (`simple_commanders_euclid.f90`);
  - `exec_calc_pspec_assemble` (`simple_commanders_euclid_distr.f90`);
  - both merges in `simple_projfile_utils.f90`.

### In memory and in the consumers

- **`simple_euclid_sigma2.f90`.**
  - The `stkind` branch of `read_groups` goes, and so does the loop that
    copies the global spectrum into every particle column.
  - `sigma2_noise(shells, fromp:top)` and `sigma2_groups(2, ngroups, shells)`
    become the two global spectra and their derived arrays.
  - `allocate_ptcls` goes: reconstruction consumers need no records.
    `read_part` stays, for the matchers only.
  - The dead fields `pinds` and `micinds`, and the uncalled `write_info`,
    go.
- **The polar calculator.**
  - `assign_sigma2_noise` and the pointer in `simple_polarft_calc.f90` and
    `simple_polarft_core.f90` take the two spectra.
  - The kernels in `simple_polarft_corr.f90` and `simple_polarft_memo.f90`
    index by half.
  - The `sigma_contrib`-absent branch of `gen_sigma_contrib` goes. It writes
    a particle's residual into the noise table, and no caller takes it.
- **Reconstruction.**
  - `prep_imgs4rec` and `gen_rec_plane` (`simple_matcher_3Drec.f90`),
    `simple_rec3D_pcg_strategy.f90`, `simple_rec3D_service.f90` and
    `simple_commanders_rec.f90` pass the half's spectrum.
  - `gen_fplane4rec` (implementation in `simple_image_ctf.f90`, interface
    in `simple_image.f90`) accepts it already resampled to the plane grid
    (section 6).
  - In the PCG path the per-particle matrix `sig2(0:R, nptcls)` goes
    entirely. Every particle of one accumulation call has the same half,
    and the state weighting only scales the spectrum:
    `sigma2_i(k) = sigma2_half(k)/w_i` (`simple_rec3D_pcg_strategy.f90`).
  - So `prep_particles` (`simple_reconstructor_pcg.f90`) takes one
    resampled inverse spectrum for the half and one weight per particle.
    Where it divides by `sigma2_i` today it multiplies by the inverse
    spectrum and `w_i`. It forms `sqrt(w_i)` times the inverse square root
    only where it builds the transfer amplitude and the whitened
    observation.
- **Class averages.** `simple_classaverager_restore.f90` drops `read_part`
  and passes the half's spectrum.
- **Cartesian slots.**
  - `cartft_calc` holds the two spectra and their inverse square roots,
    validated once. A slot records its particle's half instead of a copy
    of the spectrum.
  - `prep_cart_batch_of_build` (`simple_matcher_ptcl_batch.f90`) no longer
    copies a column per particle.
  - The slots reuse their observation and transfer arrays across batches
    instead of reallocating them, and look the shell up in a table instead
    of taking a square root per pixel.
  - The residual's diagnostics are already computed only on request.
- **FLEX.**
  - The group fallback in `ensure_canonical_sigma_state` goes.
  - The per-particle checks and the unit-spectrum fallback in
    `load_and_validate_sigma` work on the two spectra. The widening of
    `fromp:top` to the whole project for the load goes: there is no
    per-particle table to cover.
  - The stale message pointing at `sigma2_group*.star` files goes.
  - The `allocated(sigma2_noise)` check in
    `simple_flex_reconstructor_latent_ops.f90` follows the new storage.
- **`cls_expansion`.**
  - The two-policy probe in `prepare_flex_sigma2` goes, and so does its
    widening of `fromp:top` to the whole project.
  - The weight table becomes `w(ncoeff, 2)` with a half label per particle,
    and the fit routines of `simple_flex_cls_expansion.f90` take it in that
    form. Without a usable store, the class-residual weights are the same
    for both halves.

### The STAR converter, retired

`sigma2_convert` imports a grouped STAR file as a lossy seed, expanded over
the project's stacks, and exports the stored spectra as a grouped STAR.
Under release 4 it is an old-format reader and writer, and it goes entirely:

- `commander_sigma2_convert` and `exec_sigma2_convert`
  (`simple_commanders_euclid.f90`);
- its exec entry (`simple_exec_other.f90`) and its UI program
  (`simple_ui_other.f90`);
- the `sigma_action` key: its parameter, its registration, and its check in
  `simple_parameters_phases.f90`;
- the `SIGMA2_PROV_STAR_SEED` provenance (`simple_sigma2_state_file.f90`),
  as ruled below;
- the STAR helpers in `simple_euclid_sigma2.f90`, which exist for it alone:
  `write_groups_starfile`, `read_sigma2_groups_file`, `read_sigma2_groups`,
  `init_from_group_header`, the `LENSTR` constant and the STAR wrapper
  import.

**Existing global stores with the STAR-seed provenance are refused and
rebuilt (ruled).**

- Only a store imported by `sigma2_convert` and never updated since carries
  it, because every matcher commit sets the provenance to residual.
- Validation refuses the provenance value, and the central rebuild entry
  rebuilds such a store from particle power spectra. The imported seed is
  lost. This follows the release-4 rule.

## 5. Renames

With one policy, the "group" names describe something that no longer
exists. The renames are accepted, and in the repository's practice they come
without aliases.

| Today | Proposed |
| --- | --- |
| `calc_group_sigmas` (private program and commander) | `commit_sigma2` |
| `load_sigma2_groups` | `load_sigma2` |
| `read_groups` (`euclid_sigma2`) | `read_global` |
| `sigma2_state_reduce_groups` | `sigma2_state_reduce` |
| `sigma2_state_read_groups`, `sigma2_state_write_groups` | `sigma2_state_read_global`, `sigma2_state_write_global` |
| header field `grouped_offset` | `global_offset` (same word on disk; checksum word 18 stays reserved) |
| `sigma2_group_iter` | `sigma2_commit_iter` |

## 6. Optimisations

Each item states what changes and how to measure it. A gain is claimed only
from SIMPLE's benchmark files on the same run, before and after.

- **Memory.**
  - Matchers keep their partition's records, because they update and write
    them. They lose the per-particle noise table: 168 floats per particle
    at box 336.
  - Reconstruction consumers lose both the table and the record copy
    (`allocate_ptcls`). Class-average restoration loses both the table and
    the records it reads. FLEX and `cls_expansion` lose both over the whole
    project.
  - Measure with the peak RSS that the refine2D and refine3D matcher
    benchmark files already record.
- **Weights computed once.**
  - The kernels read `k/sigma2` for the particle's half from a small array
    instead of dividing per shell per evaluation.
  - With the spectra fixed for the pass and each particle's half known at
    memoisation, `memoize_sqsum_ptcl` forms the weighted self-term once.
    `wsqsum_ptcl_now` becomes a lookup instead of a loop over the shells.
    This holds only if memoisation runs after the spectra are assigned; the
    order is checked when the change is made.
  - Measure with the "alignment search" timing of
    `REFINE2D_BENCH_ITERnnn.txt` and the refine3D matcher files.
- **Folding the weight into memoised transforms, to benchmark.**
  - Folding it into the reference transforms (`ft_ref`, `ft_ref2`) misses
    the common case. When the reference is shifted, the kernel transforms
    the shifted reference afresh and never reads `ft_ref`.
  - Folding it into the particle transforms (`ft_ptcl_ctf`, `ft_ctf2`)
    covers shifted and unshifted references alike. They are memoised per
    particle, whose half is known, and the weight multiply leaves the inner
    loop.
  - The same particle arrays feed the correlation kernels (`gen_corrs`, the
    normalised-correlation gradient). Folded arrays therefore serve only
    passes that score with the Euclidean objective alone. Anything else
    would need a second copy of the particle memo, which is large.
  - The other Euclidean kernels that read the weight directly (the
    gradients, `gen_many_euclids`, the GPU kernel) are checked one by one.
- **The GPU path.** Two weight vectors replace the `wks(shells, nptcls)`
  table built and uploaded per batch.
- **Reconstruction planes, cached for the grid's lifetime.**
  - The two spectra are resampled to the plane grid once, together with the
    shell lookup table and the inverse and inverse-square-root spectra.
  - The cache is keyed by the store generation, the plane grid (padded box
    and pixel size) and the reconstruction mode. It is rebuilt only when one
    of them changes.
  - `gen_fplane4rec` then does no allocation, resample or division per
    particle. It multiplies by the cached arrays.
  - Class-average restoration uses the same cache, and the PCG path takes
    its per-half inverse spectrum from it (section 4).
- **Class-average restoration buckets particles by class.**
  - `cavger_update_sums` loops over every class and, inside, over every
    particle of the batch, skipping the other classes' particles. Bucketing
    the batch by class once makes it one pass.
  - This is not a sigma2 change, but it sits in code this refactoring
    edits. Measure with the "class averaging" timing of the refine2D files.
- **No unused reads.**
  - Class-average restoration no longer reads the partition's records.
  - Reconstruction consumers no longer allocate a record copy.
- **Fewer copies in the residual path.**
  - `calc_sigma2` allocates a contribution vector per particle, and
    `store_contribution` copies it into the record column. The residual
    routines write into the column directly.
  - `read_part` reads the records into a temporary and copies them, and
    `write_sigma2` copies them into another temporary. Reading straight
    into the record array and writing straight from it removes both copies.
- **The power-spectrum bootstrap.**
  - Today each partition writes a file with every one of its rows: zeros
    except the sampled particles, which are scaled.
  - The master reads each file twice. The first read sums by half. The
    second only learns the partition's row range, before writing the seed
    into the candidate.
  - At 750,000 particles that is 504 MB written and about 1 GB read, for
    25,000 sampled spectra. The master then writes the seeded candidate and
    reads it twice more, for the reduction and the validation.
  - Instead, the workers write per-half sums of the rows they would have
    written, and the master takes the row ranges from the project.
  - The normalisation stays exactly as it is:
    - a sampled spectrum is scaled by `nactive/(2 nsample)` of its half,
      then sanitised: non-finite values become 1.0, and a spectrum with no
      value above `DTINY` becomes 1.0 throughout;
    - the master divides each half's sum by the half's number of active
      particles, sampled or not, which gives the sampled mean over two for
      every spectrum that was not floored;
    - it subtracts the average image's power over two, then repairs
      non-positive values (`remove_negative_sigmas`);
    - the average image is the sum of the sampled images whose spectra were
      valid, divided by the number of active particles in both halves (at the
      merge, by the number of images summed instead: section 15).
  - The selection stays the master's shared selection file, drawn per half.
    Drawing the sample inside the workers would change which particles are
    sampled, and with them the seed.
  - The master writes the seeded records in blocks. Every active row of a
    half holds the same spectrum, so the global spectra are known exactly
    and the reduction and validation re-reads go; the seed's validity is
    checked once instead of per row. When the committed prefix of an
    appended project is retained, its records are summed while they are
    copied.
  - Summation order changes, so the result is compared within tolerance,
    not bit for bit.
  - `sigma2_binfile` then leaves this path. Nothing else uses it, so the
    type and its re-export from `simple_pftc_srch_api.f90` go too.
- **Store I/O per commit.**
  - Today one matcher iteration reads the store five times and writes it
    three times:
    - each worker reads its partition's records (`read_part`) and writes
      them all back as a range file, updated or not;
    - the master copies the committed store into the candidate
      (`sigma2_state_create_candidate`), merges the range files into it
      (`sigma2_state_merge_local_ranges`), reduces it, and validates it
      (`sigma2_state_validate_science`), the last two each a full read;
    - each range merge opens the candidate, writes its rows and syncs the
      whole file (`sigma2_state_write_particles`). The candidate creation,
      the reduction and the publication sync again.
  - The copy is wasted. The merge requires the ranges to cover every row
    exactly, so every copied byte is overwritten.
  - **A streaming commit, version two kept (ruled: part of this
    refactoring).**
    - The master creates a fresh candidate at full size and fills it in row
      blocks. Each row comes from a range file if one covers it, otherwise
      from the committed store.
    - Every active, valid row is added to two `real64` sums and counts as
      it passes. The master then writes the global spectra and syncs once,
      through one file handle.
    - With today's ranges, which cover every row, nothing is read from the
      committed store. The iteration then reads the store twice and writes
      it twice (the workers' record read and range write, the master's
      merge), with no format change.
    - The science validation's full re-read moves to the tests (ruled).
      The commit still checks coverage, counts and finite global spectra.
    - Range files stay until publication, so a failed commit can be rerun.
      Deleting each as it is consumed would save disk space and lose that.
  - **Sparse ranges for fractional updates (ruled: after a measurement on
    the ribosome project, not in this refactoring).**
    - Workers write only the rows they updated, with their row indices, and
      stop reading records they do not change. The range file gains a row
      list. It is a transient file, so the store's format does not change.
    - The streaming commit takes the other rows from the committed store.
      An iteration that updates a fraction f of the rows then reads the
      store about once and writes it about 1 + f times, instead of twice
      each.
    - No switch between dense and sparse is needed: a worker that updated
      every row writes them all.
    - The polish overlay (`read_pending_range`) then reads the discrete
      pass's sparse range.
  - **Incremental sums (a format version three) are not planned.**
    - The global section would keep `real64` sums and counts per half. A
      commit would subtract the replaced rows and add the new ones instead
      of summing every row.
    - With the streaming commit the sums cost nothing. The candidate must
      be a complete file before its atomic publication, so every row passes
      through the commit anyway. Version three would pay only if the copy
      bypassed the process (a file-system copy), which is not portable.
    - Should measurements ever call for it: the sums are valid only for the
      active set and half labels they were formed under. The store would
      record an FNV digest of the two (SIMPLE has no project generation
      counter). Any change of the digest forces a full reduction, and so
      does a fixed number of commits, so rounding in the running sums
      cannot accumulate. Every existing store would be rebuilt once.
  - **Measure first.** Whether the sparse ranges are worth doing depends
    on how much of an iteration the store transaction takes after the
    streaming commit. The transaction is timed in its three places, in
    SIMPLE's benchmark format:
    - the candidate preparation (the copy), in the refine2D and refine3D
      strategies;
    - the commit phases inside `commit_sigma2`: the range merge, the
      reduction, the science validation and the publication;
    - the workers' record read and range write, in the matcher benchmark
      files.
  - **Today's coverage.** The refine3D strategy file has a "group-sigma
    consolidation" line, but it times only the shared-memory commit. The
    distributed commit and the candidate preparation are untimed. The
    refine2D strategy writes no strategy file, and gains one in the same
    format.
  - **Context.** Each file records the store's particle count, shell count
    and size in bytes, and the number of rows the update replaced.
  - **Test cases.** A large project on a spinning disk (the ribosome, with
    751,729 particle rows) and one on fast storage. The test machine holds
    only beta-galactosidase, whose store is a few megabytes, so the run
    adds the timings and records them there; the ribosome measurement is
    the maintainer's.
  - **Decision.** The streaming commit is made in this refactoring (ruled).
    The sparse ranges follow if the ribosome measurement after it still
    shows the transaction as a visible share of an iteration.
- **Merges.** The chunk merge (`concatenate_canonical_states`) and
  `merge_selected_canonical_field` read each source store whole into memory
  and write it through its own open and sync, then reduce and validate by
  re-reading. Instead they copy in fixed blocks through one handle, add to
  the global sums while copying, and sync once.

## 7. Order of work

Each step leaves the tree compiling. The maintainer compiles and runs the
gates.

1. **Timing first.** Add the transaction timings of section 6. Then record
   a baseline on the test cases before anything else changes: the
   transaction timings, plus the peak memory, alignment-search,
   class-averaging and reconstruction timings the benchmark files already
   hold.
2. **One rebuild entry** (section 3), while both policies still exist.
   - The central routine in `simple_sigma2_bootstrap.f90`, called by every
     workflow entry that consumes or updates the store, `refine2D` and
     `make_cavgs` included.
   - Calls of it replace the rebuild code of refine3D, the reconstruction
     start and the `solve2D` resume; the final reconstruction keeps its
     choice of `bootstrap_rec3D` on the identity check (section 3).
3. **Parameters, UI and command lines,** the stream's included.
   - Remove `sigma_est`, `l_sigma_glob` and every set, delete, require and
     warn.
   - Every caller of the identity check passes no policy.
   - The `solve3D` manifest key as ruled.
4. **The store.**
   - The fixed policy in the header, the group arrays out of reduction,
     validation and commit, `SIGMA2_RANGE_DIRECT` out.
   - The merges, `calc_pspec` assembly and `calc_group_sigmas`.
   - The converter's retirement, and the STAR-seed provenance as ruled.
5. **The in-memory model.**
   - Two spectra and their derived arrays in `euclid_sigma2`, then the
     consumers: the polar calculator, reconstruction, PCG with one spectrum
     per half and scalar weights, class averages, the Cartesian slots by
     half, FLEX, and `cls_expansion` with `w(ncoeff, 2)`.
   - The renames of section 5 land here.
6. **Kernel and reconstruction optimisations**, one at a time with
   before-and-after benchmark files: the precomputed weights and the
   memoised self-term, the GPU weights, the reconstruction-plane cache, the
   class bucketing of restoration, the residual-path copies. Folding the
   weight into the particle transforms is benchmarked last, for
   Euclidean-only passes.
7. **The bootstrap's per-half sums.**
8. **The streaming commit and the block-copying merges** (ruled), with the
   science validation's full re-read moved to the tests.
9. **Sparse ranges,** later and outside this refactoring, if the ribosome
   measurement after step 8 still shows the transaction as a visible share
   of an iteration.

## 8. Tests

- **`simple_sigma2_state_tester.f90`.**
  - The per-stack policy case
    (`test_policy('group', SIGMA2_GROUP_STACK, ...)`) goes. A new case
    checks that a header with the stack policy is refused.
  - The streaming commit: the committed global spectra equal a full
    reduction of the records; an uncovered or doubly covered row is
    refused; a failed publication leaves the committed store and the range
    files in place.
- **`simple_sigma2_bootstrap_tester.f90`.**
  - The `sigma_est` argument and its assertions go.
  - The central rebuild: a stack-policy header, a wrong layout and a
    missing file each trigger it; a valid global store does not.
  - The per-half sums reproduce today's seed spectra within tolerance,
    including a floored spectrum and an appended project's retained prefix.
- **`simple_cartft_calc_tester.f90`.** Its store fixture loses the group
  arrays and `l_sigma_glob`. Slots by half give the same objectives and
  residuals as slots with copies, within tolerance.
- **`simple_pcg_halfset_tester.f90`.** Its cases call `prep_particles`
  without a spectrum today. A new case checks that one spectrum per half
  with scalar weights gives the same accumulations as the per-particle
  matrix, within tolerance.
- **`cls_expansion`** has no tester. The fit with `w(ncoeff, 2)` is checked
  once against the fit with the full table when the change is made.
- **`simple_solve3D_manifest_tester.f90`.** A manifest that records
  `sigma_est` is refused as carrying an unknown key, within version 2.
- **`simple_flex_pca_application_tester.f90`.** Its `sigma_est` line goes.
- **Testers that use the per-particle table, `allocate_ptcls` or the
  `sigma2_state_init_header` arguments.** They follow the new storage and
  signatures:
  - `simple_strategy3D_cont_tester.f90`;
  - `simple_strategy3D_cont_1jyx_tester.f90`;
  - `simple_pftc_inplane_tester.f90`;
  - `simple_solve3D_manifest_tester.f90`;
  - `simple_cartft_calc_tester.f90`;
  - `simple_commanders_test_highlevel_frac.f90`;
  - the `sigma2_estimate_available` calls in
    `simple_sigma2_state_tester.f90`.
- **Renamed procedures** in testers follow the renames.
- **Lightweight gates.** `scripts/check_test_registry.py` runs whenever a
  suite table or test UI changes.
- **Scientific check.** A `solve2D` and a `refine3D` run on a project whose
  store is already global must give the same results before and after
  steps 2 to 5, within the existing numerical tolerance. Steps 6 to 8
  change the order of floating-point operations, so they are compared
  within tolerance as well.

## 9. Documents

- **Living documents.**
  - The forced-setting lines `sigma_est=global` in `solve3D_policy.md`,
    `refine3D_auto_policy.md`, `solve3D_cavgs_policy.md` and
    `refine3D_states_cavgs_policy.md` go.
  - "grouping policy" and "grouped sigma" wording goes from
    `solve3D_cavgs_policy.md`, `refine3D_policy.md`, `solve3D_policy.md` and
    `reconstruct3D_pcg_policy.md`.
  - The example command line in `doc/how2s/memory_estimator.md` goes.
  - The refusal of `sigma_est` in `refine3D_states_cavgs_policy.md` goes.
  - The `sigma_est=global` entry in
    `doc/policies/stream/reference_generation_policy.md` goes.
  - So do the mentions in the planned notes `single_pcg_integration_completed.md` and
    `refine3D_states_cavgs.md`.
  - Stale comments about grouping go from `simple_flex_pca_strategy.f90`,
    `simple_final_rec.f90`, `simple_sigma2_bootstrap.f90` and
    `simple_sigma2_files.f90`.
- **The skill file** `.github/skills/simple-refine3d/references/bayesian-model.md`
  lists `sigma_est`. It is edited only with your approval.
- **History documents** stay as written, including
  `canonical_sigma2_state_refactoring.md`. A completion note for this
  refactoring records the change.
- **Generated indexes** under `doc/code_overview/` are regenerated, not
  edited.

## 10. Consequences

- **Anyone who ran `sigma_est=group`** gets the global spectrum. That is the
  decision.
- **A project whose store was written under `group`** is rebuilt from
  particle power spectra at the entry of the next workflow that uses it,
  once the rebuild is central (step 2). Its residual records are lost, and
  the next run starts from image power.
- **A store with the STAR-seed provenance** is rebuilt the same way, and
  its imported seed is lost.
- **The stream** changes in its three `sigma_est` lines and loses the dead
  `imgkind=sigma2` loop.
- **`sigma2_convert` users,** if any, lose STAR import and export.
- **FLEX and `cls_expansion`** lose their policy fallbacks. The per-stack
  empty-half failure they worked around cannot occur with one group. A
  selection whose active particles all fall in one half of the project
  still makes the reduction refuse, as it does under the global policy
  today.
- **The bootstrap** gives the same seed up to summation order.

## 11. Open after the rulings

Waiting for measurements:

- **The sparse ranges** wait for the maintainer's ribosome measurement
  after the streaming commit (section 6). The incremental store reduction
  of the earlier ruling became the format version three of section 6,
  which is not planned because the streaming commit makes its sums free.

Needing approval:

- **The skill file** `.github/skills/simple-refine3d/references/bayesian-model.md`
  (section 9). The run proposes its new text in the report and does not
  edit it.

Decided at the merge (section 15):

- **The bootstrap's average image** divided the sum of at most 25,000
  sampled images by the number of all active particles. Up to 25,000 active
  particles that is the sample mean. Above, it shrank the subtracted power
  by the square of the sampled fraction: about 1/900 at 750,000 particles.
  The run preserved it; the merge divides by the number of images summed.

## 12. Execution

This plan is carried out as an unattended run on the Oracle Linux test
machine (arun run `global_sigma2_only`), in five phases. Nothing is
committed; the maintainer reviews the complete diff when all phases are
done. Every phase records itself in section 13.

Common to every phase:

- **Data.** Beta-galactosidase (5,513 particles, box 256 at 1.275 Å per
  pixel, `mskdiam=180 pgrp=d2`, 50 selected class averages), copied into a
  fresh directory under the run's scratch directory for each run; the
  originals are never written.
- **Execution settings.** `nparts=4 nthr=6` for distributed runs and
  `nthr=24` for shared-memory runs, identical in every phase. Timed runs
  use an out-of-tree Release build with tests built; nothing else heavy
  runs while a run is timed.
- **Validation.** Beta-galactosidase is large enough for the Fourier shell
  correlation (FSC). Every `solve3D` run uses `nsample=2500` and starts
  from the one common `solve3D_cavgs` result with `cavg_ini_ext=yes`.
  The simulated workflow tests are judged by their ground-truth
  orientation and shift errors.
- **Acceptance of a repeated case.** The final FSC 0.143 resolution of
  every state lies within the larger of the Phase 0 spread and 3 % of its
  Phase 0 value. Ground-truth errors lie within their Phase 0 values plus
  the same margin.
- **Acceptance of a repeated `solve2D` case** (cases E; the maintainer's
  ruling of 2026-10-12 replaces the class-count and median rule, which the
  unchanged code itself missed). Every run completes. Two measures of its
  last iteration are judged against the four `solve2D` runs on unchanged
  code (Phase 0's E1 and E2, and the two reruns of the base build made in
  Phase 2):
  - the mean score (the last `SCORE [0,1]` average of the run's log) is no
    lower than 0.97 times the lowest of the four, 0.389, so at least 0.377;
  - the mean class resolution over the populated classes (the last
    `RESOLUTION @ FSC=0.143` average of the run's log, taken over the
    classes of the active particles, so empty classes do not count) is no
    worse than 1.03 times the worst of the four, 11.151 Å, so at most
    11.49 Å.
  The number of populated classes and the median class resolution (from the
  final ranked class list) are reported, not judged. The maintainer's second
  ruling of 2026-10-12 confirms the particle-weighted reading; the unweighted
  mean over the final class list is reported, not judged.
- **Gates.** The fast gate (run by `./compile_debug.sh`), the CTest unit
  suites of the touched areas, and `python3 scripts/check_test_registry.py .`
  whenever tests or their registration change.
- **Not in this run.** The generated indexes under `doc/code_overview/`
  (the maintainer regenerates them), the skill file of section 9 (its new
  text goes into the report), and the sparse ranges of section 6.

### 12.1 Phase 0: timings and baseline

1. Find how SIMPLE is built on this machine and build the tree as it is
   (`./compile_debug.sh`, which also runs the fast gate), plus the Release
   build for timed runs.
2. Add the transaction timings of section 6, the only source change of this
   phase, in SIMPLE's benchmark conventions (`L_BENCH_GLOB`, `tic`/`toc`,
   no new test arrays):
   - the refine3D strategy file: the candidate preparation and the
     distributed commit, next to the existing "group-sigma consolidation"
     line of the shared-memory commit;
   - a refine2D strategy benchmark file in the same format, with the
     candidate preparation and the commit;
   - the commit's phases inside `calc_group_sigmas` (range merge,
     reduction, science validation, publication), in a benchmark file of
     their own (`SIGMA2_COMMIT_BENCH.txt`);
   - the workers' record read and range write, in the refine2D and refine3D
     matcher benchmark files;
   - the context of each: the store's particle count, shell count and size
     in bytes, and the number of rows the update replaced.
3. Run the baseline cases on the Release build with the timings:
   - **S.** `solve3D_cavgs`, single state; its final project and maps are
     kept as the common start.
   - **A.** `solve3D nsample=2500 cavg_ini_ext=yes` from S, distributed,
     twice (the spread).
   - **B.** `refine3D` from A's final project with `objfun=euclid` and
     `update_frac=0.2`, five iterations, distributed, twice.
   - **C.** B's command line in shared memory, once.
   - **D.** `reconstruct3D rec_backend=pcg` from B's final project, once.
   - **E.** `solve2D` on the particles, distributed, twice; then
     `make_cavgs` with `ml_reg=yes` on its final project, once.
   - **G.** A copy of the particle project with a store written under
     `sigma_est=group` (one `calc_pspec` run), kept for the rebuild checks
     of Phases 1 and 2.
   For every case record the exact command line, the wall time, the final
   FSC 0.5 and 0.143 resolutions per state (for E the number of classes and
   the median class resolution), the peak memory and the alignment-search,
   class-averaging and reconstruction timings of the benchmark files, and
   the transaction timings.
4. Run the simulated workflow tests that use the Euclidean objective and
   record their ground-truth errors.

Exit criteria:

- the build and the fast gate pass;
- the transaction timings appear in the benchmark files of B, C and E;
- every case is recorded with its command line in section 13, with the
  spreads of A, B and E;
- S and G are kept for the later phases.

### 12.2 Phase 1: one rebuild entry, parameters, UI and command lines

Steps 2 and 3 of section 7: the central rebuild routine of section 3 and
its callers, while both policies still exist; then `sigma_est` and
`l_sigma_glob` out of parameters, UI, command lines and the stream, and the
manifest key out within version 2.

Exit criteria:

- the build, the fast gate and the registry check pass;
- the rebuild test of section 8 passes;
- outside history documents, the tree holds no `sigma_est` and no
  `l_sigma_glob`;
- `refine2D` (two iterations) and `make_cavgs ml_reg=yes` on copies of G
  rebuild the store at entry and complete.

### 12.3 Phase 2: the store and the in-memory model

Steps 4 and 5 of section 7: the fixed policy in the store, the converter's
retirement and the STAR-seed provenance as ruled, `SIGMA2_RANGE_DIRECT`
out, the merges and the bootstrap assembly without group arrays; then the
two spectra and their derived arrays, every consumer of section 4 and the
renames of section 5.

Exit criteria:

- the Release and the debug builds, the fast gate and the registry check
  pass, with the tests of section 8 for the store, the Cartesian slots and
  the PCG backend;
- outside history documents, the tree holds none of `SIGMA2_GROUP_STACK`,
  `SIGMA2_PROV_STAR_SEED`, `SIGMA2_RANGE_DIRECT`, `sigma2_convert`,
  `sigma_action`, `allocate_ptcls`, `read_sigma2_groups` or the old names of
  section 5;
- `refine2D`, `refine3D`, `reconstruct3D` and `make_cavgs` on copies of G
  rebuild the store at entry and complete;
- cases A to E repeated: the acceptance holds, and the peak memory of each
  is recorded against Phase 0;
- the simulated workflow tests meet the acceptance.

### 12.4 Phase 3: kernel and reconstruction optimisations, bootstrap

Steps 6 and 7 of section 7, one optimisation at a time:

- implement it, build, and measure it on the case that exercises it (B for
  the alignment search, E for class averaging, A and D for reconstruction),
  against the build without it;
- keep it if it is not slower beyond the run-to-run spread; folding the
  weight into the particle transforms is kept only if it is measurably
  faster on a pass that scores with the Euclidean objective alone;
- record each before and after in section 13, kept or dropped.

Then the bootstrap's per-half sums, with the test of section 8.

Exit criteria:

- the build, the fast gate and the registry check pass;
- every optimisation has its before-and-after entry;
- cases A to E repeated on the build with every kept optimisation: the
  acceptance holds.

### 12.5 Phase 4: the streaming commit, documents and close

1. Step 8 of section 7: the streaming commit and the block-copying merges,
   the science validation's full re-read moved to the tests, with the tests
   of section 8.
2. The documents of section 9, except the skill file and the generated
   indexes.
3. Cases A to E repeated: the acceptance holds, and the transaction timings
   are compared with Phase 0. On beta-galactosidase the store is a few
   megabytes, so also report how many full passes over the store an
   iteration makes, before and after.
   The maintainer's ruling of 2026-10-13 trims this step to what the
   streaming commit can change:
   - cases A, B, C and E run once each with the Phase 0 command lines, one
     at a time, and the acceptance of section 12 applies to the single runs;
   - D and F are not repeated (their Phase 3 results stand); the
     block-copying merges are covered by their tests;
   - after the timed cases, `simulated_workflow_1jxy`,
     `simulated_workflow_6vxx`, `cont_refine3D_1jxy` and `solve3D_addon`
     run, each judged by its own checks; `flex_pca_blobs` is not rerun;
   - the transaction timings and the full passes over the store per
     iteration are compared with Phase 0 on A, B, C and E.
   A second ruling of the same day accepts `simulated_workflow_6vxx` as
   flaky for now: it runs once in this phase and is recorded, not judged,
   and the report lists it as an open follow-up.
4. Move this plan to `doc/refactoring_notes/completed/` with `mv`, and write
   the review report beside it. The report reads on its own: what changed
   and why in plain words; the evidence against Phase 0 in tables; every
   decision taken during the run and where it is recorded; the proposed
   text for the skill file; what is left open (the sparse ranges, the
   bootstrap's average image). SIMPLE terms are explained and acronyms
   spelled out on first use; no run-internal labels.

Exit criteria:

- the build, the fast gate and the registry check pass;
- the acceptance of step 3 holds;
- steps 2 and 4 are done.

### File table

| File | Change | Phases |
| --- | --- | --- |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90`, `simple_refine2D_strategy.f90` | Transaction timings (Phase 0); rebuild entry and policy removal; renamed commit commander; streaming commit's candidate preparation. | 0 1 2 4 |
| `src/main/strategies/search/simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90`, `simple_matcher_ptcl_batch.f90`, `simple_strategy3D_cont.f90` | Worker timings (Phase 0); the half's spectrum; Cartesian slots by half; residual path. | 0 2 3 |
| `src/main/commanders/simple/simple_commanders_euclid.f90` | Commit timings (Phase 0); policy and converter removal; rename to `commit_sigma2`; streaming commit. | 0 1 2 4 |
| `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` | `sigma_est`, `l_sigma_glob`; `sigma_action` with the converter. | 1 2 |
| `src/main/ui/simple_ui_params_common.f90`, `src/main/ui/simple/simple_ui_refine2D.f90`, `simple_ui_refine3D.f90`, `simple_ui_heterogeneity.f90`, `simple_ui_other.f90`, `src/main/ui/single/single_ui_nano3D.f90` | The UI parameter and its uses; the converter's program. | 1 2 |
| `src/main/ui/simple_private_prgs.f90` | The renamed private program. | 2 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90`, `simple_commanders_solve2D.f90`, `simple_commanders_solve3D.f90`, `simple_commanders_validate.f90` | Every set, delete, require or warn; rebuild entry; renames. | 1 2 |
| `src/main/commanders/single/single_commanders_nano3D.f90` | Policy plumbing. | 1 2 |
| `src/main/commanders/simple/simple_commanders_refine2D.f90` | The unused imports go. | 1 2 |
| `src/main/commanders/simple/simple_commanders_denoise.f90` | The `cls_expansion` commander's own rebuild (`calc_pspec` when the store file is missing) goes: the `refine2D` it runs next rebuilds through the one entry. | 1 |
| `src/main/commanders/simple/simple_commanders_rec.f90` | The half's spectrum; renamed loader. | 2 |
| `src/main/commanders/simple/simple_commanders_euclid_distr.f90` | Bootstrap assembly without group arrays; per-half sums. | 1 2 3 |
| `src/main/simple_final_rec.f90`, `simple_external_reference_pose_initialization.f90`, `simple_refine3D_states_controller.f90` | Policy plumbing; rebuild entry; stale comments. | 1 2 4 |
| `src/main/solve/simple_solve3D_utils.f90` | Policy plumbing and imports. | 1 2 |
| `src/main/solve/simple_solve3D_manifest.f90` | The manifest key, within version 2; stale comments on the grouped section. | 1 4 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | The `rec3D_backends` description still says sigma2 group files are symlinked in; the test reads the project's registered store. | 4 |
| `src/main/solve/simple_solve3D_manifest_tester.f90` | Its store fixture follows the fixed-policy header and the renamed global spectra. | 2 |
| `src/utils/qsys/simple_qsys_ctrl.f90`, `simple_qsys_ctrl_tester.f90` | Ruling: a distributed master stops with an error when a part fails or vanishes, instead of waiting. | 2 |
| `src/main/simple_cavgs_work_project.f90`, `simple_cavgs_work_project_tester.f90`, `doc/policies/heterogeneity/refine3D_states_cavgs_policy.md` | Ruling: a sub-class half without members is written with state 0 and becomes an inactive row. | 2 |
| `src/main/commanders/test/simple_commanders_test_highlevel_flex.f90` | Ruling: the gate checks that every blank sub-class half of cls_expansion carries state 0. | 2 |
| `src/main/image/simple_fplane_grid.f90` (new), `src/main/flex/fit/simple_flex_reconstructor_latent_ops.f90` | The reconstruction-plane cache: the memoised plane grid and its per-half spectra, and FLEX's direct caller of `gen_rec_plane`. | 3 |
| `src/main/sigma2/simple_sigma2_pspec_sums.f90` (new) | The power-spectrum bootstrap's per-half partition sums and seed spectra, shared by the workers and the master. | 3 |
| `src/main/sigma2/simple_sigma2_bootstrap.f90` | The central rebuild routine; policy arguments; stale comments. | 1 2 3 4 |
| `src/main/sigma2/simple_euclid_sigma2.f90` | Two spectra and derived arrays; STAR helpers, dead fields; residual path; range writes. | 1 2 3 4 |
| `src/main/sigma2/simple_sigma2_state.f90` | Group arrays out; renames; streaming commit. Its tester follows the policy removal in Phase 1. | 1 2 4 |
| `src/fileio/simple_sigma2_state_file.f90` | Fixed policy, provenances, range mode; renamed fields; a particle reader that fills the caller's array (residual path); streaming commit primitives. | 2 3 4 |
| `src/fileio/simple_sigma2_files.f90` | One identity check; renamed loader; stale comments; the polar weights refreshed after a load. | 1 2 3 4 |
| `src/fileio/simple_projfile_utils.f90` | Merges without group arrays; block-copying merges. | 2 4 |
| `src/fileio/simple_sigma2_binfile.f90` | Deleted with the per-partition bootstrap files. | 3 |
| `src/main/strategies/parallelization/simple_rec3D_strategy.f90`, `simple_cls_expansion_strategy.f90`, `simple_make_cavgs_strategy.f90` | Rebuild entry; policy probe out; the half's spectrum. | 1 2 |
| `src/main/strategies/parallelization/simple_rec3D_service.f90` | The half's spectrum; renamed loader. | 2 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | One spectrum per half and scalar weights; plane cache. | 2 3 |
| `src/main/strategies/parallelization/simple_calc_pspec_strategy.f90` | Per-half sums. | 3 |
| `src/main/strategies/parallelization/simple_flex_pca_strategy.f90` | Stale comment. | 4 |
| `src/main/strategies/search/simple_matcher_3Drec.f90` | The half's spectrum; plane cache. | 2 3 |
| `src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_core.f90`, `simple_polarft_corr.f90`, `simple_polarft_memo.f90` | Spectra by half; precomputed weights; memoised self-term; GPU weights; folding. | 2 3 |
| `src/main/cftc/simple_cartft_pose_opt_tester.f90` | Its Euclidean slots follow the Cartesian slots by half (the calculator holds the spectra). | 2 |
| `src/main/cftc/simple_cartft_calc.f90` | Slots by half. Its tester's store fixture follows the policy removal in Phase 1. | 1 2 3 |
| `src/main/cftc/simple_cartft_calc_tester.f90` | Its store fixture and transaction cases follow the streaming commit, which replaces the copied candidate and the separate merge, reduction and commit. | 4 |
| `src/main/image/simple_image_ctf.f90`, `simple_image.f90` | `gen_fplane4rec` with the cached plane-grid spectra. | 2 3 |
| `src/main/class/simple_classaverager_restore.f90`, `simple_classaverager.f90`, `simple_classaverager_core.f90` | No record read; the half's spectrum; class bucketing. | 2 3 |
| `src/main/volume/simple_reconstructor_pcg.f90`, `simple_reconstructor.f90` | `prep_particles` with one spectrum and weights; comments. | 2 3 |
| `src/main/flex/run/simple_flex_pca_project_gateway.f90`, `simple_flex_pca_application.f90` | Group fallback, range widening and stale message out; the polar weights refreshed after the unit spectra. | 1 2 3 |
| `src/main/flex/fit/simple_flex_reconstructor_latent_ops.f90`, `simple_flex_pca_plane_cache.f90` | The new storage; comments. | 2 |
| `src/main/flex/simple_flex_cls_expansion.f90` | `w(ncoeff, 2)` with a half label. | 2 |
| `src/main/stream/pool2D/simple_stream_pool2D.f90`, `src/main/stream/stages/simple_stream_stage_initial_analysis.f90`, `src/main/commanders/stream/simple_commanders_stream_p06_pool2D.f90` | The three `sigma_est` lines and the dead `imgkind=sigma2` loop. | 1 |
| `src/main/exec/simple_exec_other.f90`, `src/main/apis/simple_private_exec_api.f90`, `production/simple_private_exec_driver.f90` | The converter's exec entry; the renamed private program. | 2 |
| `src/main/apis/simple_pftc_srch_api.f90` | The `sigma2_binfile` re-export. | 3 |
| `src/main/simple_builder.f90` | Comments on the sigma2 member and `gen_fplane4rec`. | 2 3 |
| `src/main/strategies/search/simple_strategy3D_cont_1jyx_tester.f90`, `src/main/pftc/simple_pftc_inplane_tester.f90`, `src/main/flex/run/simple_flex_pca_application_tester.f90`, `src/main/flex/simple_flex_pca_tester.f90` | Testers that follow the new storage and signatures. | 1 2 3 |
| `src/main/volume/simple_pcg_halfset_tester.f90`, `simple_frozen_accum_tester.f90`, `simple_trail_chain_manifest_tester.f90`, `src/main/ctf/simple_ctf_tester.f90` | The PCG weights test; signatures of `prep_particles` and `gen_fplane4rec`. | 2 3 |
| `src/main/commanders/test/simple_commanders_test_highlevel_frac.f90`, `simple_commanders_test_highlevel_addon.f90`, `simple_commanders_test_highlevel_pcg.f90` | High-level tests that follow the new storage and signatures. | 1 2 3 4 |
| `scripts/memory/benchmark_solve2d.py` | `sigma_est` out. | 1 |
| `doc/policies/3D/solve3D_policy.md`, `refine3D_auto_policy.md`, `solve3D_cavgs_policy.md`, `refine3D_policy.md`, `reconstruct3D_pcg_policy.md` | Section 9; their `sigma_est` lines already in Phase 1, and the old names of section 5 in Phase 2, whose exits require the tree to hold neither outside history documents. | 1 2 4 |
| `doc/policies/heterogeneity/refine3D_states_cavgs_policy.md`, `doc/policies/stream/reference_generation_policy.md`, `doc/how2s/memory_estimator.md` | Section 9; their `sigma_est` lines already in Phase 1, whose exit requires the tree to hold no `sigma_est` outside history documents. | 1 4 |
| `doc/implementation_notes/completed/single_pcg_integration_completed.md`, `refine3D_states_cavgs.md` | Section 9; their `sigma_est` lines already in Phase 1, whose exit requires the tree to hold no `sigma_est` outside history documents. | 1 4 |
| `doc/refactoring_notes/completed/global_sigma2_only.md`, `global_sigma2_only_report.md` | The plan moved on completion, and the report. | 4 |

## 13. Progress

### Phase 0: timings and baseline (2026-10-12, done)

**Source change.** Only the transaction timings of section 12.1, step 2.
They write in SIMPLE's benchmark format (a context section, then a timings
section), only when the global benchmark flag `L_BENCH_GLOB` is set (it is
on), with `tic`/`toc`:

- `simple_commanders_euclid.f90`: the sigma2 commit (`calc_group_sigmas`)
  times its four phases (range merge, reduction, science validation,
  publication) and writes `SIGMA2_COMMIT_BENCH_ITERnnn.txt`. The file name
  carries the iteration like every other per-iteration benchmark file, so a
  run keeps one per iteration; refine2D passes its commit the iteration the
  committed store feeds (the matcher iteration plus one), refine3D the
  matcher iteration. To time validation and publication apart, the commit
  calls `sigma2_state_validate_science` and `sigma2_state_publish` itself,
  which is all `sigma2_state_commit` does.
- `simple_refine3D_strategy.f90`: `REFINE3D_STRATEGY_BENCH_ITERnnn.txt`
  gains the candidate preparation (both execution modes), and its existing
  "group-sigma consolidation" line now also times the distributed commit.
- `simple_refine2D_strategy.f90`: a new `REFINE2D_STRATEGY_BENCH_ITERnnn.txt`
  with the candidate preparation, the commit and the iteration's total.
- `simple_matcher_ptcl_batch.f90`, `simple_strategy2D_matcher.f90`,
  `simple_strategy3D_matcher.f90`: the workers' record read (`read_part`),
  the load of the per-particle noise table (`read_groups`) and the range
  write (`write_sigma2`), in `REFINE2D_BENCH_ITERnnn.txt` (written by
  partition 1 only, as before) and `REFINE3D_BENCH_ITERnnn_PARTppp.txt`.
- Context in every file: the store's particle rows, shells and size in
  bytes; the rows the update replaced (the active particles carrying the
  current sampling marker in the strategy files, the particles the
  partition updated in the matcher files, the rows merged in the commit
  file). Timings of the sigma2 lines are in milliseconds (three decimals),
  because the store of the test data is a few megabytes.

**Build.** On this machine SIMPLE builds with GCC 15.2 (gcc-toolset-15),
CMake 3.26.5 and OpenMP, `OPENBLAS_CORETYPE=Haswell`; the fast gate's
checkers need Python 3.12 as `python3`. `./compile_debug.sh` on the base
passed the fast gate (15 of 15 suites); with the timing code it builds
without compiler warnings and passes the fast gate again. Timed runs use an
out-of-tree Release build with tests (`BUILD_TESTS=ON`, `USE_ARCHOPT=ON`,
installed into its own tree), the configuration of the maintainer's
`~/src/SIMPLE/build`; its only warnings are the five pre-existing linker
warnings "requires executable stack".

Evidence is under `/home/elmlundho/agent_runs/global_sigma2_only/scratch`:
`run_case.sh` (the exact command lines and how each run directory is
prepared), `queue.sh` (the order), `case_summary.sh`, `cls2D_summary.sh`,
`bench_summary.py`, and in `logs/` per case `case_<name>.log` (full output,
each line prefixed by its time in seconds, `/usr/bin/time -v` at the end),
`summary_<name>.txt`, `queue_phase0.log` and `bench_tables_phase0.txt` (all
benchmark files condensed). Every benchmark file of every case is kept in
`scratch/keep/bench/phase0/<case>/`. The cases ran one at a time on the
otherwise idle 24-core machine (load below 1 at each start, apart from S,
which started as the builds ended).

**Data and kept projects.** Each run starts from a fresh copy, named
`bgal.simple`, of `/data2/bgal_addon/control/bgal_addon.simple`
(beta-galactosidase after 2D classification: 5,513 particles in one stack,
box 256 at 1.275 Å per pixel, 50 selected class averages) or of a kept
project. A project refers to its maps by absolute path inside the run
directory that wrote them, and to its sigma2 store by a path relative to
the project's directory, so `run_case.sh` recreates the directories a case
reads from `scratch/keep/` before it runs. Kept for later phases:
`keep/S/1_solve3D_cavgs` (the common `solve3D` start), `keep/A/1_solve3D`
(input of B and C, final map `recvol_state01.mrc`), `keep/B/1_refine3D`
(input of D), `keep/E/1_solve2D` (input of the `make_cavgs` run) and
`keep/G/G` (the store written under `sigma_est=group`; copy it to
`scratch/runs/G` before use, its project records that directory).

**Command lines.** Every run sets `SIMPLE_PATH` to the Release build and
puts its `bin` and `scripts` first on `PATH`. Later phases repeat these
exactly (`run_case.sh <name> <kind>`):

| Case | Starts from | Command line |
| --- | --- | --- |
| S | the control project | `simple_exec prg=solve3D_cavgs projfile=bgal.simple mskdiam=180 pgrp=d2 nparts=4 nthr=6` |
| A1, A2 | `keep/S/1_solve3D_cavgs/bgal.simple` | `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6` |
| B1, B2 | `keep/A/1_solve3D/bgal.simple` | `simple_exec prg=refine3D projfile=bgal.simple vol1=<scratch>/runs/A1/1_solve3D/recvol_state01.mrc mskdiam=180 pgrp=d2 objfun=euclid update_frac=0.2 maxits=5 nparts=4 nthr=6` |
| C | as B | `simple_exec prg=refine3D projfile=bgal.simple vol1=<same> mskdiam=180 pgrp=d2 objfun=euclid update_frac=0.2 maxits=5 nthr=24 mkdir=no` |
| D | `keep/B/1_refine3D/bgal.simple` | `simple_exec prg=reconstruct3D projfile=bgal.simple mskdiam=180 pgrp=d2 rec_backend=pcg nparts=4 nthr=6` |
| E1, E2 | the control project | `simple_exec prg=solve2D projfile=bgal.simple mskdiam=180 ncls=100 nparts=4 nthr=6` |
| F | inside a copy of `keep/E/1_solve2D` | `simple_exec prg=make_cavgs projfile=bgal.simple ml_reg=yes mkdir=no nparts=4 nthr=6` |
| G | the control project | `simple_private_exec prg=calc_pspec projfile=bgal.simple sigma_est=group nparts=4 nthr=6` |

Decisions taken for the command lines:

- `solve2D` requires `ncls`; 100 is the user guide's example (about 55
  particles per class).
- `calc_pspec` is a private program; it runs distributed like the others.
- **C runs with `mkdir=no`.** With the default `mkdir=yes`, shared-memory
  `refine3D` failed in its first iteration: it rebuilds its parameters each
  iteration from a command line that still holds `mkdir=yes`, creates a
  nested `1_refine3D/1_refine3D` and no longer finds its sigma2 store
  (`logs/case_C_failed_mkdir.log`). The `refine3D_auto` family resets
  `mkdir=no` after its own parameter parsing; plain shared-memory `refine3D`
  does not. This defect of the base is outside this refactoring and is
  reported, not fixed; `mkdir=no` runs the same computation in the case
  directory.
- **F runs with `mkdir=no` in a copy of E1's final directory.** In a new
  directory `make_cavgs` finds no store, because the store path is relative
  to the project's directory and `make_cavgs` reads the store without the
  rebuild the other workflows run (section 3). Every worker stopped and the
  distributed master then waited forever, a second defect of the base
  reported here (`logs/case_F_failed_nostore.log`). The central rebuild
  entry of Phase 1 removes the first; from Phase 1 on, F could also run in a
  new directory, but it repeats this command line so the cases compare.
- Every `refine3D`, `reconstruct3D` and `solve2D` run starts in a new
  directory, so it rebuilds its store from particle power spectra at entry
  ("rebuilding canonical state: sigma2 state file does not exist").

**Results.** Wall time and peak memory are from `/usr/bin/time -v` (peak
memory is the largest resident set size of any one process of the run).
The resolutions are the Fourier shell correlation (FSC) at 0.5 and 0.143 of
the final reconstruction: `RESOLUTION_FINAL_STATE01` for S and A,
`RESOLUTION_STATE01_ITER005` (the last iteration) for B and C, and
`RESOLUTION_STATE01` for D. With `update_frac=0.2` each refine3D iteration
reconstructs from its 20 % sample only, hence the low FSC 0.5 of B and C.
For E, the number of classes holding particles and their median resolution
come from the final project's `cls2D` field. Every run exited with status
0.

| Case | Wall time | Peak memory | FSC 0.5 / 0.143 (Å), or classes and median class resolution |
| --- | --- | --- | --- |
| S | 11:02 | 3.07 GB | 7.77 / 6.94 |
| A1 | 15:00 | 3.19 GB | 4.35 / 3.84 |
| A2 | 15:09 | 3.21 GB | 4.41 / 3.84 |
| B1 | 3:10 | 3.03 GB | 8.82 / 4.41 |
| B2 | 3:10 | 3.03 GB | 8.82 / 4.35 |
| C | 2:11 | 3.70 GB | 8.82 / 4.60 |
| D | 1:43 | 6.59 GB | 4.35 / 3.80 |
| E1 | 3:13 | 0.41 GB | 87 classes, 12.09 Å |
| E2 | 3:15 | 0.41 GB | 88 classes, 11.87 Å |
| F | 0:03.4 | 0.31 GB | (class averages only) |
| G | 0:00.8 | 0.12 GB | (store only) |

Spreads (the difference between the two runs): A, 0.06 Å at FSC 0.5 and
0.00 Å at 0.143; B, 0.00 Å and 0.06 Å; E, one class and 0.22 Å. The
acceptance rule of section 12 takes the larger of these and 3 % of the
Phase 0 value.

G's store has the version-two layout with one group (the project has one
stack) but the per-stack grouping word, and the power-spectrum provenance:
what the rebuild checks of Phases 1 and 2 must refuse.

**Benchmark timings** (seconds, summed over the iterations; for the
matcher lines, the slowest partition of each iteration). The store of
these cases has 5,513 rows of 128 shells, 2.82 MB; a refine3D iteration
updates 1,103 rows (`update_frac=0.2`), a `solve2D` iteration all of them.

| Line | A1 | A2 | B1 | B2 | C | E1 | E2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| iterations | 43 | 43 | 5 | 5 | 5 | 27 | 27 |
| matcher peak memory (GB) | 2.72 | 2.72 | 2.75 | 2.75 | 3.53 | 0.35 | 0.35 |
| alignment search (3D: orientation search) | 12.83 | 12.80 | 129.73 | 128.03 | 88.09 | 11.35 | 11.70 |
| class averaging | | | | | | 20.09 | 21.19 |
| matcher partial reconstruction | 42.52 | 42.11 | 8.90 | 8.88 | 9.61 | | |
| volume assembly (`VOLASSEMBLE` total) | 111.11 | 113.74 | 17.37 | 17.58 | 16.66 | | |
| workers' record read | 0.037 | 0.036 | 0.005 | 0.005 | 0.006 | 0.016 | 0.015 |
| workers' noise-table load | 0.107 | 0.042 | 0.013 | 0.013 | 0.003 | 0.001 | 0.001 |
| workers' range write | 1.225 | 1.296 | 0.146 | 0.098 | 0.084 | 0.224 | 0.293 |
| candidate preparation | 0.641 | 0.692 | 0.031 | 0.031 | 0.024 | 0.145 | 0.143 |
| commit, as the strategy sees it | 1.512 | 1.475 | 0.245 | 0.241 | 0.167 | 0.964 | 0.935 |
| commit: range merge | 0.635 | 0.628 | 0.079 | 0.078 | 0.026 | 0.434 | 0.423 |
| commit: reduction | 0.178 | 0.169 | 0.020 | 0.021 | 0.020 | 0.121 | 0.117 |
| commit: science validation | 0.033 | 0.020 | 0.005 | 0.004 | 0.003 | 0.017 | 0.017 |
| commit: publication | 0.254 | 0.246 | 0.027 | 0.027 | 0.031 | 0.174 | 0.168 |
| commit: total inside the commander | 1.093 | 1.075 | 0.130 | 0.129 | 0.082 | 0.744 | 0.728 |

The commit as the strategy sees it includes starting the commander (reading
the project, building its toolbox). On these few megabytes the whole sigma2
transaction of a B iteration (candidate, workers' read and write, commit)
takes under 0.1 s of about 38 s, and of an E iteration about 0.05 s of 7 s;
the ribosome measurement of section 6 is the maintainer's. D writes no
benchmark file (distributed PCG `reconstruct3D`); its wall time and peak
memory above are its reconstruction measure. F's class-average assembly
(`CAVGASSEMBLE_BENCH_ITER001.txt`) took 0.36 s.

**Simulated workflow tests.** The high-level CTest entries that score with
the Euclidean objective, run one at a time on the Release build
(`cd scratch/build_release && ctest -R '^<name>$' -V`, with `SIMPLE_PATH`
and `PATH` set to it; log `logs/highlevel_phase0.log`). Two of them judge
the final map against the simulated truth, one the per-particle poses:

| Test | Result | Time | Measures |
| --- | --- | --- | --- |
| `simulated_workflow_1jxy` (`solve2D` then `solve3D`) | passed | 276 s | whole-volume correlation 0.899 (floor 0.80), truth FSC 0.143 at 4.46 Å (limit 40 Å) |
| `simulated_workflow_6vxx` | passed | 808 s | whole-volume correlation 0.877, truth FSC 0.143 at 5.43 Å |
| `solve3D_addon` | **failed on the base** | 12 s | see below |
| `cont_refine3D_1jxy` | passed | 508 s | ground-truth pose errors, below |

`solve3D_addon` stops in its first `solve3D`: "balance=class needs a 2D
solution (classes in ptcl2D); the project has none: run with balance=none"
(`simple_commanders_solve3D.f90`). The test's simulated project has no 2D
classes and its command line does not set `balance=none`. This failure is
in the base and unrelated to sigma2; the test gives no baseline until the
test or the default is fixed, which is the maintainer's decision and not
part of this run.

`cont_refine3D_1jxy` (1,000 simulated 1JYX particles, starting poses the
truth perturbed by 15° and 2 pixels; median errors against the truth, the
gate's floor in brackets). Stage a uses the Euclidean objective, b the
correlation objective, c repeats a distributed in two parts, e
`refine3D_auto`, f the polish, g and h the external-reference pose
initialisation (pose errors up to the D2 symmetry mates):

| Stage | Rotation error (°) | Shift error (px) | Fraction within the 8 Å basin |
| --- | --- | --- | --- |
| a (Euclidean) | 0.908 [1.5] | 0.067 [0.2] | 0.998 [0.95] |
| b (correlation) | 0.226 [1.5] | 0.057 [0.2] | 0.997 [0.95] |
| c (a, distributed) | 0.908 [1.5] | 0.067 [0.2] | 0.998 [0.95] |
| e (`refine3D_auto`) | 0.915 [7.5] | 0.117 [1.0] | 0.999 |

Symmetry-aware pose errors: f without polish 0.935°, with polish 0.629°
[1.035]; h polar 1.486° [11.54], h continued 0.398° [1.586] (shift
0.331 px); g 0.186° [1.586].

**Exit criteria.**

- The build and the fast gate pass: `./compile_debug.sh` with the final
  timing code builds without compiler warnings, passes the fast gate (15 of
  15) and installs (`logs/build_debug_final.log`; on the base,
  `logs/build_debug_base.log`); the Release build has no new warnings
  (`logs/build_release.log`, `logs/build_release_incr1.log`).
- The transaction timings appear in the benchmark files of B, C and E: the
  table above.
- Every case is recorded with its command line, with the spreads of A, B
  and E: above.
- S and G are kept for the later phases, with A, B and E: `scratch/keep/`.

### Phase 1: one rebuild entry, parameters, UI and command lines (2026-10-12)

**Inventory before editing** (the code at the start of Phase 1; "store" is
the canonical sigma2 store `sigma2_state.bin`, "identity check" is
`canonical_sigma2_consumable`, which validates the store and compares its
grid, particle layout and grouping with the project).

The identity check and its two wrappers:

- `canonical_sigma2_consumable(project, os, box, smpd, l_sigma_glob,
  message)` (`simple_sigma2_files.f90`); under `group` it recounts stacks
  from `stkind`. Called by `load_sigma2_groups` (same file, a hard check
  without rebuild), the refine3D strategy (`canonical_sigma2_needs_bootstrap`),
  the final reconstruction (`canonical_final_rec_needs_bootstrap`), the
  bootstrap module (`sigma2_estimate_available`), `cls_expansion`'s
  two-policy probe (`prepare_flex_sigma2`) and the `solve3D_addon`
  high-level test.
- `sigma2_estimate_available(projfile, box, smpd, l_sigma_glob)` reads the
  `ptcl3D` field only; called by `ensure_sigma2_for_iteration` and by the
  `solve3D_addon` frozen-project check (`simple_solve3D_utils.f90`).
- `ensure_sigma2_for_iteration(..., l_sigma_glob, label, l_bootstrapped,
  sigma_est)` rebuilds through `calc_pspec` when the store is not
  consumable; called by the `solve3D` stage start (`simple_solve3D_utils.f90`)
  and FLEX (`ensure_canonical_sigma_state`, after its group fallback).

Workflow-local rebuilds (a check, then `calc_pspec`):

- refine3D strategy, shared memory (`inmem_initialize`) and distributed
  (`distr_initialize`, two branches: continued and fresh runs), through
  `canonical_sigma2_needs_bootstrap`, which also requires `ptcl3D`.
- reconstruction start (`ensure_canonical_sigma_state` in
  `simple_rec3D_strategy.f90`), its own copy of the identity check with
  both policies; run by the master and, in distributed runs, by every
  worker.
- `solve2D` checkpoint resume (`ensure_resume_sigma_state`), its own copy of
  the identity check with both policies.
- final reconstruction (`simple_final_rec.f90`): the identity check decides
  between a plain `reconstruct3D` and `bootstrap_rec3D` (seed, bootstrap map,
  residual pass); it does not call `calc_pspec` itself.
- `cls_expansion` commander (`simple_commanders_denoise.f90`): runs
  `calc_pspec` when the store file is missing, then `refine2D`.

Unconditional seeds, not rebuilds (they always reseed and stay): the
`refine3D_auto` startup, `bootstrap_rec3D`, the external-reference pose
initialisation's sigma bootstrap, and `solve2D`'s first stage.

No check at all: `refine2D` (its first candidate preparation fails on a bad
store) and `make_cavgs` (`cavger_read_euclid_sigma2` reads the store when
`ml_reg=yes`).

`sigma_est` and `l_sigma_glob` outside the history documents: the
parameter, its registration and derivation (`simple_parameters.f90`,
`simple_parameters_parse.f90`, `simple_parameters_phases.f90`); the UI
parameter (`simple_ui_params_common.f90`) and its seven uses (`refine3D`,
`refine3D_auto`, `refine3D_nano`, `flex_pca`, `refine3D_states`,
`sigma2_convert`, and `solve2D`'s own input); the commanders `refine3D`
(eight sites, including the `bootstrap_rec3D` warning), `solve2D` (three),
`solve3D` (two), `validate` (two), `nano3D` (one), `euclid` (the converter's
policy checks), `euclid_distr` (the bootstrap's group count and header
policy); `simple_external_reference_pose_initialization.f90` (three),
`simple_sigma2_bootstrap.f90`, `simple_sigma2_files.f90`,
`simple_euclid_sigma2.f90` (the global/stack branch of `read_groups`),
`simple_final_rec.f90`, `simple_refine3D_strategy.f90` (two),
`simple_refine3D_states_controller.f90`, `simple_rec3D_strategy.f90`,
`simple_solve3D_utils.f90` (two), `simple_flex_pca_project_gateway.f90`
(the group fallback), `simple_cls_expansion_strategy.f90` (the probe), the
`solve3D` manifest keys, the stream's three sites, the testers
`simple_sigma2_bootstrap_tester.f90`, `simple_cartft_calc_tester.f90` and
`simple_flex_pca_application_tester.f90`, the memory benchmark script, and
the living documents of section 9.

The code matches the plan's description, with three additions the plan
does not name: the reconstruction workers run the rebuild check too; the
`cls_expansion` commander has a file-exists rebuild before its `refine2D`
pass; and `simple_euclid_sigma2.f90` and the two `euclid` commanders branch
on `l_sigma_glob`, so with the parameter gone they keep only their global
branch in this phase (the store layer itself changes in Phase 2).

**What changed.**

- **One rebuild entry.** `ensure_sigma2_for_iteration(template, projfile,
  oritype, iter, box, smpd, label, l_bootstrapped)` in
  `simple_sigma2_bootstrap.f90` checks the store with
  `sigma2_estimate_available(projfile, oritype, box, smpd)`. That function
  reads the project's `ptcl2D` or `ptcl3D` field (the store's particle layout
  is taken over the field of the consumer) and calls the one identity check.
  When the store is not consumable, the entry runs `calc_pspec` with a fresh,
  sparse command line (`prepare_pspec_cline`, which now sets the field
  explicitly). The identity check `canonical_sigma2_consumable(project, os,
  box, smpd, message)` takes no policy and requires the global grouping with
  one group.
- **Its callers.**
  - the refine3D strategy in both execution modes;
  - the reconstruction start (the master or the shared-memory run; the parts
    of a distributed run no longer run the check);
  - the `solve2D` checkpoint resume;
  - the `solve3D` stage start and FLEX, as before;
  - new: `refine2D` at entry (masters and shared memory, not the workers);
  - new: `make_cavgs` with `ml_reg=yes` (the master, or the shared-memory run
    after it has written the even/odd halves).
  The copies of the identity check in the reconstruction strategy and the
  `solve2D` resume, and refine3D's `canonical_sigma2_needs_bootstrap`, are
  gone. `refine2D` now writes the even/odd split before it builds its
  references, so a `make_cavgs` run inside the reference initialisation seeds
  on the halves `refine2D` uses. The `cls_expansion` commander's own
  file-exists rebuild is gone; the `refine2D` round it runs rebuilds.
- **The final reconstruction** keeps its own decision (recorded in sections 3
  and 7): a store that is not consumable sends it to `bootstrap_rec3D`, which
  seeds and then runs a residual pass at the final sampling. A seed alone
  would skip that pass. It uses the identity check without policy.
- **`sigma_est` and `l_sigma_glob` are gone** from the parameters, the UI
  parameter and its seven program uses, every commander, the stream (with
  the dead `imgkind=sigma2` loop in the initial analysis), the `solve3D`
  manifest's recorded and replayed keys (within version 2), the memory
  benchmark script and the living documents. `cls_expansion` loses its
  two-policy probe and FLEX its group fallback. Where code branched on
  `l_sigma_glob` (`read_groups`, the converter, the bootstrap assembly), only
  the global branch remains. The store layer itself keeps its group
  machinery until Phase 2.
- **Tests.**
  - `simple_sigma2_bootstrap_tester.f90`: the `sigma_est` argument and its
    assertions go. A new `test_rebuild_entry` writes a small project and
    stores. A valid global store is consumable and the rebuild entry leaves
    it alone. A store with the stack policy, of another particle layout or of
    another grid, and a missing store, are each not consumable, which is
    exactly what makes the entry rebuild.
  - `simple_solve3D_manifest_tester.f90`: a manifest of version 2 that
    records `sigma_est` is refused as carrying an unknown input key.
  - Other testers follow the new signatures (`simple_sigma2_state_tester.f90`,
    the `solve3D_addon` high-level test) or lose a `sigma_est` or
    `l_sigma_glob` line (`simple_cartft_calc_tester.f90`,
    `simple_flex_pca_application_tester.f90`).
- **File table.** Two changes:
  - the living documents' `sigma_est` lines moved into this phase (their
    rows now read "1 4"), because this phase's exit requires the tree to hold
    no `sigma_est` outside history documents;
  - `simple_commanders_denoise.f90` was added for the `cls_expansion` rebuild.
- **Files changed in this phase** (51): the plan; the eight living
  documents; the memory benchmark script; and the Fortran sources and
  testers listed by `git status` under `src/` (the Phase 0 files among them
  only for the policy and rebuild changes above). `scripts/check_test_registry.py`
  was not affected.

**Evidence** (under `/home/elmlundho/agent_runs/global_sigma2_only/scratch/logs`):

- **Build and gates.** `./compile_debug.sh` builds without compiler warnings,
  passes the fast gate (15 of 15) and installs (`p1_build_debug_final.log`).
  `python3 scripts/check_test_registry.py .` passes (`p1_registry.log`). The
  Release build has no new warnings (`p1_build_release2.log`).
- **Unit suites.** The sigma2 bootstrap suite passes with 39 checks,
  including the rebuild test; the `solve3D` manifest suite passes with 125.
- **Rebuild on copies of G** (the store written under `sigma_est=group`).
  Each run recreates `scratch/runs/G` from `scratch/keep/G/G`, because G's
  project records that directory, and runs in it with `mkdir=no`:
  - `simple_private_exec prg=refine2D_distr projfile=bgal.simple mskdiam=180
    ncls=100 objfun=euclid maxits=2 mkdir=no nparts=4 nthr=6` (`refine2D` is
    the private program `refine2D_distr`; `prg=refine2D` is its worker). At
    entry: "canonical sigma2 grouping policy does not match", then a rebuild
    from particle power spectra. Two iterations, exit 0, 12 s
    (`case_G_r2D.log`).
  - `simple_exec prg=make_cavgs projfile=bgal.simple ml_reg=yes mkdir=no
    nparts=4 nthr=6`: the same refusal and rebuild, exit 0, 4 s
    (`case_G_mkcavgs.log`).
  - The rebuilt store has the global grouping word, one group, generation 2
    and the power-spectrum provenance.
- **High-level tests on the Release build** (not required by this phase;
  run because the rebuild entries moved):
  - `simulated_workflow_1jxy` passed, whole-volume correlation 0.899 and
    truth FSC 0.143 at 4.54 Å (Phase 0: 0.899, 4.46 Å);
  - `cont_refine3D_1jxy` passed;
  - `mini_stream_1jxy` passed;
  - `flex_pca_blobs` stops with a segmentation fault in its fixture setup
    (`phantom_layout`), before any sigma2 code. A Release build of the
    unchanged base, exported with `git archive HEAD`, fails identically
    (`flex_pca_blobs_base_segfault.log`). The failure is in the base and is
    reported, not fixed.
  Log: `highlevel_phase1.log`.

**Exit criteria.**

- The build, the fast gate and the registry check pass: above.
- The rebuild test of section 8 passes: `test_rebuild_entry`.
- Outside history documents the tree holds no `sigma_est` and no
  `l_sigma_glob`. Two exceptions remain:
  - the skill file `.github/skills/simple-refine3d/references/bayesian-model.md`,
    which section 9 leaves to the maintainer's approval;
  - the manifest tester's refusal case, which section 8 requires.
- `refine2D` (two iterations) and `make_cavgs ml_reg=yes` on copies of G
  rebuild the store at entry and complete: above.

**For the maintainer.**

- The rebuild check now runs at every `refine2D` entry, including each
  `solve2D` stage and each stream pool round. That check is a deep validation
  of the store: one full read.
- A store that cannot be consumed is now rebuilt, where `refine2D` used to
  stop. A defect in the stream's store merges would therefore show up as
  repeated reseeding, not as a stop.
- The stale comment about group structure in `simple_calc_pspec_strategy.f90`
  goes in Phase 3, which edits that file.

### Phase 2: the store and the in-memory model (2026-10-12 and 13; stopped on the 2D acceptance, resumed under the ruling, stopped on a second FLEX defect, resumed under the second ruling, stopped on blank sub-class averages, resumed under the third ruling, done)

**Consumers confirmed against the code.** Every consumer of section 4 was read before editing; the
plan's list matched, with these call-site notes:

- The polar calculator holds the particle's half per local slot (`iseven(i)`, `i = pinds(iptcl)`),
  so its column is `merge(1, 2, iseven(i))`. The two particle-batch builders memoized the particle's
  self-term before setting its half; the half is now set first.
- `gen_sigma_contrib` had an unused branch that wrote into the noise table; its output argument is
  now mandatory and the branch is gone.
- The class averager read the particle records and then overwrote them with the spectra; it now
  reads only the two spectra.
- The PCG (preconditioned conjugate gradient) reconstruction had a dead `whiten_observation`; it is
  deleted.
- `gen_fplane4rec` keeps its interface in this phase: callers pass their half's spectrum. The
  pre-resampled plane-grid input lands with the plane cache in Phase 3, which edits the same routine.

**What changed.**

- **The store.** The header validation accepts only one global spectrum per half (grouping word
  global, one group), the power-spectrum or residual provenance, and the local range mode. A stack
  grouping, the STAR-seed provenance and the direct range mode are refused with a message, and the
  rebuild entry then rebuilds. `sigma2_state_init_header` writes the fixed policy, with no grouping
  argument. The converter (`sigma2_convert`, its `sigma_action` parameter and its UI entry) is
  retired. The merges of `simple_projfile_utils.f90`, the bootstrap assembly, the reduction and the
  commit carry no group arrays.
- **The renames of section 5**, without aliases: `commit_sigma2` (private program, commander and
  stop message `SIMPLE_COMMIT_SIGMA2`), `load_sigma2`, `read_global`, `sigma2_state_reduce`,
  `sigma2_state_read_global` and `sigma2_state_write_global`, the header fields `global_offset` and
  `global_checksum` (same words on disk), and `sigma2_commit_iter`. The refine3D bench line reads
  "sigma2 commit".
- **The in-memory model.** `euclid_sigma2` holds `sigma2_noise(kfrom:kto, 2)`, the two spectra, with
  the public lookup `sigma2_half(eo)` (even is 1, odd is 2; any other value is refused). The
  per-particle noise table, `allocate_ptcls`, the STAR helpers and the dead fields are gone.
  `read_part` only fills the matchers' own records.
- **Consumers.**
  - Polar kernels, both memo routines and the offload path index by half.
  - The Cartesian calculator holds the two spectra, their inverse square roots and their validity
    once (`set_sigma2`), and each particle slot holds its half.
  - Reconstruction planes and class averages take their half's spectrum.
  - The PCG backend takes one spectrum per half plus scalar particle weights, using
    `sigma2_i = sigma2_half / w_i` (arrays `inv_sig2`, `rsqrt_sig2`, `wts`, `sqrt_wts`).
  - `cls_expansion` takes `w(ncoeff, 2)` and a half index per member.
  - The FLEX gateway checks the two spectra and no longer widens the particle range.
- **A Phase 1 defect fixed.** The rebuild entry read the project file from disk. In a new execution
  directory (`mkdir=yes`) the copied project still records the previous run's directory, so the
  entry found the previous run's store. Every consumer's builder resolves the store against its
  working directory and found none, and refine3D stopped at its first iteration. Phase 1's rebuild
  checks all ran with `mkdir=no`. `sigma2_estimate_available` now resolves the store against the
  working directory, as the builders do, so a new execution directory rebuilds its own store, as the
  base did (Phase 0's B1 and D logs show the bootstrap).
- **Cleanups from review.** The dead `swap_ptclsevenodd` is deleted. The testers delete their
  temporary stores. A few "group" comments are reworded.
- **Tests.**
  - The sigma2 state suite: a global transaction, plus a test that patches each retired header word
    (with a valid checksum) and expects its refusal. New public helpers `write_test_sigma2_store` and
    `patch_sigma2_header_word`.
  - The bootstrap suite: the stack-grouping case through a patched header word, plus
    `test_new_exec_dir`, which shows that a copied project in a new directory does not consume the
    store of the directory it records.
  - The Cartesian calculator suite: slots by half against slots with their own copy of the
    spectrum, and the round trip on the global spectra.
  - The PCG half-set suite: `test_half_spectrum_weights`. One accumulation with the half's spectrum
    and particle weights equals the sum of single-particle accumulations with `sigma2/w_i`, to a
    relative tolerance of 1e-4.
  - The manifest, pose optimizer, continuous in-plane and refine3D continuous testers follow the new
    interfaces.
- **File table.** Rows were added for `simple_cartft_pose_opt_tester.f90` and
  `simple_solve3D_manifest_tester.f90`. The row containing `reconstruct3D_pcg_policy.md` now reads
  "1 2 4", for its `load_sigma2_groups` paragraph.
- **Files changed in this phase** (52 at the end, all assigned by the file table; the four defect fixes added the queue controller and its tester, the FLEX application tester and the addon high-level test, below): the plan;
  `doc/policies/3D/reconstruct3D_pcg_policy.md`; the private-program driver, API and UI; the exec
  and UI files of the converter; the parameter files; and the Fortran sources and testers under
  `src/` listed by `git status`.

**Evidence** (under `/home/elmlundho/agent_runs/global_sigma2_only/scratch/logs` unless noted).

- **Builds and gates.**
  - `./compile_debug.sh` builds without compiler warnings (only the linker's executable-stack
    notes, as before), passes the fast gate (15 of 15) and installs (`p2_build_debug_final.log`).
  - The Release build is clean (`p2_build_release3.log`).
  - `python3 scripts/check_test_registry.py .` passes (`p2_registry.log`).
  - `git diff --summary` shows no mode changes.
- **Unit suites** (`scratch/p2test/`): sigma2 state 130 checks; sigma2 bootstrap 43; `solve3D`
  manifest 125; Cartesian calculator 167; pose optimizer 71; continuous in-plane 42; PCG half-set 43,
  including the new test; addon report docking 2.
  - On the final debug build, the unit tests `unit_cart_align3D`, `unit_pftc_align2D3D`,
    `unit_project` and `unit_reconstruction` pass, leaving no temporary stores
    (`scratch/p2test/final_*.log`).
  - The library reconstruction run's third suite, and the full debug run, stop with a segmentation
    fault in the FLEX phantom fixture (`phantom_layout`). This is the base failure recorded in
    Phase 1, before any sigma2 code.
- **Names.** Outside history documents (`doc/implementation_notes/completed/`) the tree holds none
  of the names this phase removes or renames.
- **Rebuild on copies of G** (Release build, `scratch/runs/G` recreated from `scratch/keep/G/G`,
  `mkdir=no`). Each logs "sigma2 state is not one global spectrum per half (a retired grouping
  policy)", rebuilds from particle power spectra and exits 0:
  - `refine2D` (`simple_private_exec prg=refine2D_distr ... maxits=2`), 12 s;
  - `refine3D` (`simple_exec prg=refine3D projfile=bgal.simple vol1=<A1 map> mskdiam=180 pgrp=d2
    objfun=euclid lp=8 maxits=1 mkdir=no nparts=4 nthr=6`), 33 s;
  - `reconstruct3D` (`mskdiam=180 pgrp=d2 mkdir=no nparts=4 nthr=6`), 16 s;
  - `make_cavgs ml_reg=yes`, 4 s.
  The first `refine3D` attempt, without `lp`, stopped because G has no resolution history
  (`case_p2G_r3D_failed_nolp.log`).
- **Cases A to F** repeated with the Phase 0 command lines on the new Release build, from
  `scratch/keep/`. Logs: `case_p2*.log` and `queue_phase2.log`; benchmark files in
  `scratch/keep/bench/phase2/`; benchmark comparison in
  `scratch/keep/bench/phase2/comparison_phase0_phase2.txt`.

| Case | Wall time | Peak memory (Phase 0) | Result (Phase 0) |
| --- | --- | --- | --- |
| A1 | 14:47 | 3.23 GB (3.19) | 4.41 / 3.84 Å (4.35 / 3.84) |
| A2 | 15:06 | 3.24 GB (3.21) | 4.35 / 3.80 Å (4.41 / 3.84) |
| B1 | 3:03 | 3.03 GB (3.03) | 8.82 / 4.35 Å (8.82 / 4.41) |
| B2 | 3:05 | 3.03 GB (3.03) | 8.82 / 4.35 Å (8.82 / 4.35) |
| C | 2:13 | 3.70 GB (3.70) | 8.16 / 4.66 Å (8.82 / 4.60) |
| D | 1:41 | 6.43 GB (6.59) | 4.35 / 3.80 Å (4.35 / 3.80) |
| E1 | 3:12 | 0.42 GB (0.41) | 85 classes, 14.84 Å (87, 12.09) |
| E2 | 3:15 | 0.43 GB (0.41) | 91 classes, 14.19 Å (88, 11.87) |
| F | 0:03.4 | | completes (0:03.4) |

  - The FSC (Fourier shell correlation) 0.143 values of A to D meet the acceptance (C is 1.3 %
    above Phase 0, under the 3 % margin).
  - The matcher's peak memory is unchanged.
  - The noise-table load falls from up to 0.107 s to under 1 ms, because only the two spectra are
    read; the particle record read falls by about a third. The rest of the sigma2 transaction is
    unchanged within noise.
  - The B orientation search took 121.5 and 123.8 s (Phase 0: 129.7 and 128.0).
  - The first B1 attempt is the defect described above (`case_p2B1_failed_stalecwd.log`).
- **The E acceptance of the plan failed, and the base failed it too** (the first stop). The global spectra at the end of E1 agree
  with Phase 0's to within 1 %, and every run's last iteration has the same mean score
  (0.389 ± 0.007) and a mean class resolution of 10.6 to 11.4 Å. Further runs of the same command
  line (`queue_spread.log`; base and Phase 1 trees built from `git archive HEAD` and the Phase 2
  start snapshot; populated classes and median class resolution):

| Build | Runs |
| --- | --- |
| Base, Phase 0 | 87, 12.09 Å; 88, 11.87 Å |
| Base, now | 90, 12.32 Å; 89, 10.20 Å |
| Phase 1 | 83, 9.60 Å; 87, 10.20 Å |
| Phase 2 | 85, 14.84 Å; 91, 14.19 Å; 89, 9.60 Å; 90, 12.63 Å |

  - The unchanged base misses the rule against its own Phase 0 runs (10.20 Å against 12.09 Å; 90
    classes against 87).
  - The rule's allowance (the Phase 0 spread of one class and 0.22 Å, or 3 %) is far below the run-
    to-run spread of these two numbers.
  - Two short `refine2D` runs of the same build agree on under 1 % of class labels, so the 2D
    assignments cannot be compared run by run.
- **Simulated workflow tests** (Release, `highlevel_phase2.log`):
  - `simulated_workflow_1jxy` passed: correlation 0.903, FSC 0.143 at 4.30 Å (Phase 0: 0.899,
    4.46 Å);
  - `simulated_workflow_6vxx` passed: 0.882 and 5.20 Å (Phase 0: 0.877, 5.43 Å);
  - `cont_refine3D_1jxy` passed, with its ground-truth errors within Phase 0's plus 3 %;
  - `solve3D_addon` fails at its start (`balance=class` needs a 2D solution), as on the base.

**Resumed under the maintainer's ruling of 2026-10-12.**

- **The `solve2D` acceptance applied** (rule recorded in section 12). The measures come from each run's
  log: the last `SCORE [0,1]` average, and the last `RESOLUTION @ FSC=0.143` average. That average is
  taken over the active particles, so it averages the class resolutions over the populated classes,
  weighted by population. The run directories were deleted after the first stop; each log's final
  ranked class list reproduces the `cls2D` counts and medians exactly.
  - Reference (base): mean scores 0.389 in all four runs, so the floor is 0.377. Mean class
    resolutions are 10.845 and 11.002 Å (Phase 0), and 11.151 and 10.825 Å (the reruns), so the
    limit is 1.03 × 11.151 = 11.49 Å.
  - Phase 2: E1 0.389, 10.783 Å; E2 0.389, 11.396 Å; E3 0.389, 10.763 Å; E4 0.389, 11.051 Å.
    **All four pass.** The reported, unjudged numbers are in the table above.
  - Reading of the ruling, for the maintainer to confirm: "mean class resolution over the populated
    classes" is taken as that log average, the measure quoted in the stop report the ruling answers.
    The unweighted mean over the final ranked class list behaves differently, because a few classes
    of two to four particles sit at the 326 Å floor and dominate it. On that measure the base gives
    18.28, 18.37, 23.03 and 18.36 Å (limit 23.72 Å). Phase 2's E2 (31.96 Å) would fail, and so would a
    Phase 1 run (29.25 Å), whose scoring is the base's.
- **Defect: shared-memory `refine3D` with `mkdir=yes` nested its directory.** Cause: every iteration
  rebuilds its parameters from the command line, which still held `mkdir=yes`, so each iteration
  created and entered a new execution directory and lost the sigma2 store. Fix: `inmem_initialize`
  (`simple_refine3D_strategy.f90`) sets `mkdir=no` on the command line once its parameters are
  built, as the distributed initialize does. Evidence: case C's command line with `mkdir=yes`
  (`run_case.sh` kind `refine3D_shmem_mkdir`) completes in 2:26 at 3.69 GB with one execution
  directory, five commits and 8.59 / 4.47 Å (`case_p2Cmk.log`, `queue_phase2fix.log`). Case C keeps
  its `mkdir=no` command line for the comparisons.
- **Defect: a distributed master waited forever when its parts died.** Cause: a partition script
  wrote no exit status, and the non-streaming scheduler only polled for each part's finished file,
  which a crashed or killed part never writes.
  - Fix (`simple_qsys_ctrl.f90`, file-table row added):
    - Partition scripts and array-job elements record themselves and write their exit status, with
      the streaming path's job record (`simple_qsys_job_record.f90`).
    - A part's earlier status and record are deleted before it is submitted.
    - `update_queue` marks a submitted part failed on a non-zero status. A part without a status gets
      the liveness check, which writes the lost status after two misses.
    - `schedule_jobs` and `schedule_array_jobs` then cancel the parts still running and stop with an
      error naming the failed parts, their exit status, their scripts and the shared output file.
    - The streaming path is unchanged.
  - Tests (`simple_qsys_ctrl_tester.f90`):
    - `test_failed_parts` checks the scripts' record and status lines, and that a non-zero or lost
      status fails a part while a zero status does not.
    - `test_failed_local_parts` submits real local part scripts whose executable fails and sees both
      fail within seconds.
    - The suite passes 60 checks.
  - Evidence: a distributed `refine3D` on case B's inputs, with part 2's executable killed with
    signal 9, stops 13 s into the run with "REFINE3D: 1 PART(S) FAILED; part 2 (exit status 137,
    script distr_simple_script_2); the output of every part is in .../SIMPLE_SUBPROC_OUTPUT". No part
    is left running (`case_killpart.log`, `killpart_killer.log`).
- **Defect: `flex_pca_blobs` and the fractional-reconstruction suite stopped with a segmentation
  fault in `phantom_layout`.** Cause, found with the debug build and the disassembly of the faulting
  address: the routine's first statement, `ntot = NPTCLS`, meant the module constant. Its optional
  dummy argument `nptcls` hides that constant (Fortran names are case-insensitive), so the statement
  read an absent argument. That has been so since commit `f83f12687`. Fix (tester only): the default
  is read as `FLEX_PHANTOM_NPTCLS`, the constant's public name. What FLEX computes is unchanged.
  - Behind the fault, the fractional-reconstruction suite then stopped in this phase's
    `sigma2_half` with a garbage half (96) from `gen_rec_plane`, although the particle's half was 0
    or 1 a line before. The half was computed by a function call inside the array-section subscript
    of an actual argument bound to `gen_fplane4rec`'s explicit-shape optional dummy, inside an
    OpenMP loop. It is now computed into a local variable first, there and at the three other call
    sites of that form (class averaging, the PCG strategy, the fractional high-level test).
    `sigma2_half`'s refusal now prints the value. The mathematics is unchanged, and no Release case
    run had met a refusal.
  - Evidence: the full debug library reconstruction run completes, with the PCG half-set (43
    checks), addon report docking (2) and fractional reconstruction (17) suites passing
    (`scratch/p2test/r_lib_reconstruction_dbg2.log`).
  - **`flex_pca_blobs` does not pass yet** (`highlevel_phase2_fixes.log`). Past the fixture, all its
    FLEX checks pass:
    - label accuracy 0.9995 (shared) and 0.9975 (distributed);
    - weighted state maps against the truth 0.990 and 0.991, above floors of 0.977.
    It then stops in its `refine3D_states_cavgs` stage, "zero right-hand side; nothing to solve with
    PCG", when a state has lost all its particles. The base, with only the fixture line fixed, fails
    the same way at iteration 3 for state 1, where Phase 2 fails at iteration 2 for state 2
    (`base_tree/flex_pca_blobs_base.log`, recorded in the journal; the tree is deleted). This is a
    second defect of the existing code, described under "For the maintainer".
- **Defect: the `solve3D_addon` test refused its own command line.** Cause: its simulated 2D
  stand-in had class labels but no class FRCs, so `solve3D` found no 2D solution for
  `balance=class`. Fix (`simple_commanders_test_highlevel_addon.f90`): the fixture writes class FRCs
  for its ten stand-in classes (a smooth fall-off through 0.5 at the gate's final low-pass limit;
  the runs force their low-pass range, so the values plan nothing) and registers them as the
  project's `frc2D`. Evidence: `solve3D_addon` passes on the Release build in 650 s.
  - Baseline for Phases 3 and 4: union map against truth 0.949 (floor 0.90); base map 0.951; union
    minus base −0.002 (floor −0.03); FSC 0.143 against truth 4.83 Å for both maps; report union
    against base 0.9989 (floor 0.90); cohort coverage 1.00; pose error, cohort against frozen,
    4.25° (floor 15); excess 0.50° (floor 5).
- **After the fixes**, on the final Release build (`highlevel_phase2_final.log`):
  - `simulated_workflow_1jxy`: 0.903 and 4.46 Å;
  - `cont_refine3D_1jxy`: passes;
  - `simulated_workflow_6vxx`: 0.884 and 5.09 Å;
  - case B1 rerun: 3:01, 3.03 GB, 8.82 / 4.35 Å, identical to before the fixes.
  Final checks:
  - `./compile_debug.sh`: no compiler warnings and the fast gate 15 of 15
    (`p2r_build_debug_final.log`);
  - the registry check passes (`p2r_registry.log`);
  - no mode changes;
  - 52 files changed in this phase, all assigned by the file table.

**Resumed under the maintainer's second ruling (2026-10-12).**

- The particle-weighted reading of the `solve2D` acceptance is confirmed (section 12), so E stands
  as applied above.
- **Defect: the distributed PCG reconstruction stopped on an empty state** (the fifth item).
  - Cause: the PCG master (`simple_rec3D_pcg_strategy.f90`) carried a state forward only when its
    trailing sample was below 0.001 of the state. A state left without particles reached the solve
    with a zero right-hand side and stopped the run. The gridding path drops such a state
    (`determine_dropped_states`), and the `refine3D_states_cavgs` policy says an empty state never
    fails the run.
  - Fix, the gridding rule:
    - Before the state loop, `find_empty_states` sums each state's particles over the raw
      accumulator headers of every part and both halves.
    - When another state has particles, a state with none is reported, dropped for the iteration,
      and its previous maps are carried forward, marked as carried as in the existing "no sample
      and no trailing chain" branch.
    - Without previous maps the solve stops with a message naming the state. When no state has
      particles, the existing refusal stays.
    - As on the gridding path, an add-on assembly and a `state=` run never drop a state, and a
      state that still has a trailing chain keeps the chain path.
  - Test (`simple_rec3D_service_tester.f90`, `test_empty_state_carried`): PCG `reconstruct3D` on
    the two-conformer phantom, then every particle moved to state 1 and run again in the same
    directory. The run completes and logs "STATE 2 HAS NO PARTICLES; DROPPING AND CARRYING
    PREVIOUS VOLUME FORWARD". State 2's map and half maps are unchanged byte for byte, and state 1
    is reconstructed again. The fractional-reconstruction suite passes 19 checks
    (`scratch/p2test/r_frac_empty.log`).
  - **`flex_pca_blobs` still fails, for a reason the rule does not cover**
    (`highlevel_phase2_empty.log`). In its `refine3D_states_cavgs` stage, iteration 2, state 1 is
    not empty: the raw headers count 3 rows in each half. Its right-hand side is zero all the same,
    because the rows it keeps are blank. `cls_expansion` writes 9 sub-class averages per half
    whose images are entirely zero: even 18, 48, 71, 77, 80, 92, 96, 101 and 110; odd 17, 72, 78,
    81, 90, 93, 95, 102 and 111. These score NaN and stay in the stage as active rows. After
    iteration 1, nearly every row with signal has moved to state 2 (174 and 175 rows), and states
    1 and 3 keep only a few rows, among them blank ones.
  - The base, with only the fixture line fixed, writes the same blank averages at the same
    indices, and fails at the same point with the same counts (3/3, 174/173, 3/4)
    (`flex_pca_blobs_base_ruling2.log`). So this is a third defect of the existing code, not of
    this run.
- After the change: `./compile_debug.sh` has no compiler warnings and passes the fast gate 15 of
  15 (`p2e_build_debug_final.log`). The registry check passes (`p2e_registry.log`). There are no
  mode changes, and 53 files changed in this phase, all assigned by the file table.

**Resumed under the maintainer's third ruling (2026-10-13).**

- **Defect: `cls_expansion` wrote a sub-class half without members as an active, blank average** (the
  sixth item).
  - Cause: the even and odd averages of a sub-class are restored from the members at even and odd
    positions of the class (`flex_cls_half_reproducibility`). A half with no weight in a sub-class
    came out as an all-zero image. The records gave every sub-class `state 1`, and the work project of
    `refine3D_states_cavgs` turned both halves into active rows, which scored NaN. On the phantom
    fixture that is 9 halves per stack.
  - Fix (`simple_cls_expansion_strategy.f90`):
    - `record_class` sums each sub-class's weight in its even half and in its odd half.
    - The sums travel to the master as columns seven and eight of the part class maps, as separation
      and reproducibility already do, and are written into the merged map.
    - `apply_split_project_updates` writes `state_even` and `state_odd` in `cls2D`: 0 for a half without
      members, 1 otherwise.
  - The work project (`simple_cavgs_work_project.f90`) makes a row with a half flagged 0 inactive, in
    `ptcl3D` and `ptcl2D`. That link is needed, although the ruling expected no downstream change: the
    rows took only the sub-class's `cls2D` state. A project without the flags reads as before, and the
    better-row rule is unchanged. The `refine3D_states_cavgs` policy describes it (section 2, the work
    project).
  - Tests:
    - `simple_cavgs_work_project_tester.f90`, `test_memberless_halves`: a flagged even or odd half is
      an inactive row and its other half stays active. The suite passes 117 checks
      (`scratch/p2test/r_cavgs_wp.log`).
    - `flex_pca_blobs` gains the gate check that every blank half carries state 0, and only those. On
      the Release build: 18 blank halves, 18 flagged, check passed.
- **`flex_pca_blobs`, out of this run's acceptance** (the ruling), run once on the Release build
  (`highlevel_phase2_blank.log`). **It still fails, on one check: a known failure of the existing
  code.**
  - 81 checks pass, among them:
    - truth label accuracy 0.9995 (shared memory) and 0.9975 (distributed);
    - leading eigenvolume against the truth 0.956 and 0.920;
    - refined state maps against the truth (weighted) 0.990 and 0.991, above floors of 0.977 and 0.976;
    - the `states_cavgs` checks: projects written and unchanged where they must be, and the blank-half
      check.
  - The empty-state rule now carries states 1 and 3 through the later iterations, in which they hold
    no rows.
  - The failing check is "two or three states are delivered". `refine3D_states_cavgs` ends with all
    180 sub-class averages in one state ("ONLY 1 STATE HOLDS SUB-CLASS AVERAGES"). The collapse is in
    the existing code: the base, with only the fixture line fixed, already had 174 and 173 rows in
    state 2 and only 3 or 4 in each other state at iteration 2 (`flex_pca_blobs_base_ruling2.log`).
    Following the ruling, it is not chased further.
- **Final checks:**
  - `./compile_debug.sh` has no compiler warnings and passes the fast gate 15 of 15
    (`p2b_build_debug_final.log`);
  - the registry check passes (`p2b_registry.log`);
  - there are no mode changes;
  - 57 files changed in this phase, all assigned by the file table;
  - no job is left running.

**Exit criteria.**

- Builds, fast gate, registry check and the section 8 tests for the store, the Cartesian slots and
  PCG: met.
- No retired or old names outside history documents: met.
- `refine2D`, `refine3D`, `reconstruct3D` and `make_cavgs` on copies of G rebuild at entry and
  complete: met.
- Cases A to E repeated, with peak memory recorded: A to D meet the FSC acceptance; E meets the
  ruling's `solve2D` acceptance (under the reading above).
- The simulated workflow tests meet the acceptance: met.
- The rulings' defects:
  - the directory nesting, the waiting master and the phantom fault are fixed, with their evidence;
  - `solve3D_addon` passes and is the baseline for Phases 3 and 4;
  - the library reconstruction run completes;
  - the PCG empty-state rule is in place and tested;
  - memberless sub-class halves are written with state 0.
- `flex_pca_blobs`, out of the acceptance by the third ruling, fails one check: the collapse into one
  state in `refine3D_states_cavgs`, which is in the existing code. It is recorded as a known failure.

**For the maintainer.**

- Known failure of the existing code: `flex_pca_blobs` fails because `refine3D_states_cavgs`
  collapses the phantom's sub-class averages into one state.
- A particle with `eo` outside 0 and 1 is now refused (`sigma2_half`), also inside parallel loops.
  The polar path alone maps `eo` of -1 to the even half, as before.
- Rebuilding a refused store restarts at generation 1. Before, the rebuild continued the old
  generation count.
- `cartft_calc%set_sigma2` replaces the spectra under slots that are already filled. Production
  sets them once per batch before filling; the tester helper sets them per slot.
- The per-half derived arrays of `euclid_sigma2`, the reuse of the slot arrays and the shell lookup
  table are Phase 3 work (section 6).

### Phase 3: kernel and reconstruction optimisations, bootstrap (2026-10-13)

**How each optimisation was measured.** Every optimisation was built twice on the Release
configuration from snapshot copies of the source: the tree without it (`scratch/src_ref` built in
`scratch/build_ref`) and the tree with it (`scratch/src_new` built in `scratch/build_new`). The two
ran the unchanged Phase 0 command lines alternately, in both orders, with nothing else heavy running
(`scratch/ab.sh`). The case that exercises the change was used: B (refine3D) for the alignment
search, E (solve2D) for class averaging, A and D for reconstruction. Single runs of the same build
differ by up to 5 % on this machine in a session (Phase 0 saw 1.3 %), so a change counts as faster
or slower only beyond that spread. The benchmark files are in
`scratch/keep/bench/phase3/<optimisation>/`; the queue logs in `scratch/logs/queue_phase3/`.

| Optimisation | Measured on | Without | With | Decision |
| --- | --- | --- | --- | --- |
| Shell weights computed once | B orientation search, four pairs | 117.2 / 123.5 / 127.9 / 127.5 s (mean 124.0) | 125.3 / 126.9 / 119.4 / 124.9 s (mean 124.1) | kept: no change beyond the spread |
| | E alignment search, two pairs | 11.88 / 11.39 s | 12.02 / 11.63 s | |
| GPU weight vectors | no GPU on this machine | | | coded, not measurable |
| Reconstruction-plane cache | E class averaging, two pairs | 20.16 / 23.12 s | 17.99 / 21.24 s | kept: faster on average, within the spread |
| | B partial reconstruction, two pairs | 8.93 / 8.89 s | 9.33 / 9.10 s | |
| | A partial reconstruction, one pair | 41.98 s | 42.96 s | the untouched volume assembly (111.1 / 114.6 s) and search (12.4 / 12.9 s) moved by the same 3 %: machine noise |
| | D (PCG reconstruct3D) wall time, two pairs | 1:49.2 / 1:47.9 | 1:48.6 / 1:47.5 | maps identical (FSC 0.143 3.795 Å) |
| Class bucketing in restoration | E class averaging, three pairs | 18.45 / 18.87 / 19.41 s (mean 18.91) | 17.69 / 19.18 / 17.89 s (mean 18.25) | kept: not slower |
| Residual-path copies | workers' record read, B and E, two pairs | B 0.004 / 0.001, E 0.006 / 0.002 s | 0.000 s in all four | kept; range write and peak memory unchanged |
| Cartesian shell table and slot reuse | library 1JYX pose-recovery suite, five pairs | mean 15.59 s (14.5 to 17.3) | mean 15.10 s (14.4 to 15.8) | kept; maps and scores identical |
| Folding the weight into the particle transforms | B orientation search, two pairs; E alignment search, three pairs | B 125.2 / 128.9, E 11.45 / 11.96 / 11.28 s | B 124.9 / 127.5, E 11.78 / 11.16 / 11.33 s | dropped: not measurably faster |
| Bootstrap per-half sums (step 7) | `calc_pspec` on the control project, two pairs | 0.75 / 0.74 s | 0.71 / 0.73 s | in (required); the stored spectra are bit-identical |

**What changed.**

- **Shell weights computed once (kept).** The polar calculator keeps the Euclidean weights
  k/sigma2(k) of the two spectra, in the single and double precision each kernel already used
  (`sigma2_wk`, `sigma2_wk8` in `simple_polarft_calc.f90`). `memoize_sigma2_weights`
  (`simple_polarft_memo.f90`) forms them, and the weighted self-term of every memoised particle,
  whenever the spectra are filled: at `assign_sigma2_noise` when the values are already there,
  after `read_global` in the matcher's `prep_sigmas_objfun`, in `load_sigma2` and after FLEX's unit
  spectra. `wsqsum_ptcl_now`, which recomputed the self-term on every evaluation, is now a lookup.
  The kernels use exactly the expressions they computed before, so scores are unchanged bit for
  bit. The routine does nothing on a calculator that was never constructed (`make_cavgs`,
  `reconstruct3D` attach spectra to one) or while the spectra are still the zeros
  `euclid_sigma2` allocates. The dead per-thread buffers `w_weights` and `sumsq_cache` are gone, and
  `reallocate_ptcls` now reallocates the per-shell square sums the refresh reads.
- **The GPU path (coded, not measurable here).** The OpenMP-offload benchmark routines upload the
  two weight vectors and each particle's half instead of a weight table per shell and particle, and
  read the memoised self-terms. The machine has no GPU and the offload build is off; a syntax check
  with the offload definitions enabled shows only the cuFFT interfaces that build provides.
- **Reconstruction-plane cache (kept).** A new type `fplane_grid` (`simple_fplane_grid.f90`) holds
  the shell lookup table of a padded Fourier plane and, per half, the spectrum resampled to it as an
  inverse and an inverse square root. `euclid_sigma2` keeps one (`fplane`, `prep_fplane_grid`),
  rebuilt only when the padded box or the lowest shell changes and dropped whenever the spectra are
  read again or killed; that invalidation stands in for the plan's store-generation key and also
  covers FLEX's unit spectra, which have no generation. `gen_fplane4rec` now takes the grid and the
  particle's half (its `kfromto` and raw-spectrum arguments are gone) and multiplies by the cached
  arrays: no allocation, resample, square-root table or division per particle. Gridding
  (`prep_imgs4rec`, `gen_rec_plane`), class-average restoration and FLEX's cached path build or
  fetch the grid once per call; the PCG worker reads its per-half inverse spectrum from the same
  grid at the padded shells 2n, which hold exactly sigma2(n), so `prep_particles` takes the inverse
  (`inv_sig2`) and the PCG maps are unchanged bit for bit. Gridding and restoration change in the last
  bits (a multiply by a stored reciprocal instead of a division).
- **Class bucketing in restoration (kept).** `cavger_update_sums` sorts the batch by class once
  (a stable counting sort) instead of scanning the whole batch for every class; each class's sums
  accumulate in the same order as before.
- **Fewer copies in the residual path (kept).** `read_part` reads the records straight into the
  record array through a new reader that fills the caller's array
  (`sigma2_state_read_particles_into`); the polar residual (`gen_sigma_contrib`) is written straight
  into the particle's record column; `write_sigma2` writes the records without a converted copy.
  The Cartesian residual keeps its allocating interface: it is formed from the calculator's own
  double-precision sums with the slot's shell bounds, twelve tester call sites rely on that, and the
  copy it saves (about a hundred values per particle) is negligible next to its per-pixel loop.
- **Cartesian slots (kept; section 4, deferred from Phase 2).** The calculator builds a shell table
  once per box and looks each pixel's shell up instead of taking a square root per pixel per
  evaluation (`in_shells` is gone); the residual takes the square root of the spectrum once per
  shell; a slot keeps its observation and transfer arrays across fills (every reader visits only
  the pixels of the slot's shells, which each fill writes). A new check in the batch-preparation
  test refills each slot under the other objective, half and band first and asserts the result
  equals a slot filled once.
- **Folding the weight into the particle transforms (dropped).** Benchmarked last, for passes that
  score with the Euclidean objective alone: the matchers folded each particle's shell weight into
  `ft_ptcl_ctf` and `ft_ctf2`, the Euclidean FFT kernels skipped the multiply and the correlation
  kernels refused folded transforms. It was not measurably faster, so it was not applied; the patch
  is kept as `scratch/keep/bench/phase3/O7_folding/folding_not_kept.patch`.
- **The bootstrap in per-half sums (step 7).** The `calc_pspec` workers write, per partition, the
  per-half sums of the scaled and sanitised spectra of their sampled particles (a new module,
  `simple_sigma2_pspec_sums.f90`) instead of a file with every row. The master adds the partitions'
  sums, divides each half by its number of active particles taken from the project, subtracts the
  average image's power over two and repairs non-positive values, all unchanged
  (`bootstrap_seed_spectra`; the average image is still divided by the number of active particles in
  both halves). It writes the candidate's records in blocks (the retained committed prefix of an
  appended project, then each active row with its half's seed), sums the valid active records into
  the global spectra while writing them, exactly as the reduction would, writes those and publishes.
  The reduction and science-validation re-reads are gone; publication keeps its file validation.
  The part files stay until publication. `sigma2_binfile` and its re-export from
  `simple_pftc_srch_api.f90` are deleted.
- **Tests.**
  - `simple_sigma2_bootstrap_tester.f90`, `test_pspec_sums`: the per-half sums through two partition
    files reproduce the row-wise seed spectra (a floored spectrum and a repaired shell included),
    within 1e-6; the global spectra summed while writing (a retained prefix with an invalid record
    included) equal a full reduction of the same records.
  - `simple_cartft_calc_tester.f90`: the refill check above.
  - Testers follow the new `gen_fplane4rec` and `prep_particles` signatures (CTF, frozen
    accumulator, PCG half-set, reconstructor, 1JYX continuous, high-level PCG and fractional).
- **Review.** A review of the diff found that forming the weights from the zero-filled spectra at
  assignment traps Debug builds (they trap on division by zero); fixed as above, and a Debug
  Euclidean `refine2D` now completes. Its smaller findings (refresh after `load_sigma2` and FLEX,
  the grid at `box_croppd`, part files kept until publication, a worker check of the half) are
  applied.
- **File table.** Rows were added for the plane-cache module with FLEX's caller, the bootstrap-sums
  module, and Phase 3 on `simple_sigma2_state_file.f90`, `simple_sigma2_files.f90` and the FLEX
  gateway.

**Evidence** (under `/home/elmlundho/agent_runs/global_sigma2_only/scratch`; logs in `logs/`).

- **Builds and gates.** `./compile_debug.sh` builds without compiler warnings (only the linker's
  executable-stack notes), passes the fast gate 15 of 15 and installs (`p3_build_debug_final.log`).
  `python3 scripts/check_test_registry.py .` passes (`p3_registry.log`). The Release build is clean
  (`p3_build_release_final.log`). `git diff --summary` shows only the deletion of
  `simple_sigma2_binfile.f90`, no mode change. 34 files changed in this phase, all assigned by the
  file table.
- **Unit suites** (Debug, `p3test/F_*.log`): `unit_pftc_align2D3D` (sigma2 state 130 checks, sigma2
  bootstrap 56 with `test_pspec_sums`), `unit_cart_align3D` (with the refill check),
  `unit_reconstruction`, `unit_project`, `unit_core` and the CTF suite (72) pass; the library
  reconstruction run (PCG half-set 43, fractional reconstruction 19) passes with the plane cache
  (`p3test/O3_lib_reconstruction.log`). A Debug Euclidean `refine2D` on a copy of G, with
  floating-point traps on, completes (`case_p3dbgG_r2D.log`).
- **Cases A to E on the build with every kept optimisation** (`queue_phase3/final.log`; benchmark
  files in `keep/bench/phase3/final/`, comparison with Phase 0 in
  `keep/bench/phase3/final/comparison_phase0_phase3.txt`):

| Case | Wall time | Peak memory (Phase 0) | Result (Phase 0) |
| --- | --- | --- | --- |
| A1 | 15:00 | 3.19 GB (3.19) | 4.35 / 3.84 Å (4.35 / 3.84) |
| A2 | 14:41 | 3.21 GB (3.21) | 4.41 / 3.84 Å (4.41 / 3.84) |
| B1 | 3:16 | 3.03 GB (3.03) | 8.82 / 4.35 Å (8.82 / 4.41) |
| B2 | 3:08 | 3.03 GB (3.03) | 8.82 / 4.41 Å (8.82 / 4.35) |
| C | 2:16 | 3.69 GB (3.70) | 8.82 / 4.60 Å (8.82 / 4.60) |
| D | 1:48 | 6.11 GB (6.59) | FSC 0.143 3.80 Å (3.80) |
| E1 | 3:07 | 0.41 GB (0.41) | score 0.389, mean class resolution 10.96 Å; 86 classes, median 9.89 Å |
| E2 | 3:07 | 0.42 GB (0.41) | score 0.389, 10.74 Å; 88 classes, median 11.48 Å |

  - The FSC 0.143 resolutions of A to D equal Phase 0's within its spread; E meets the `solve2D`
    acceptance of section 12 (score at least 0.377, mean class resolution at most 11.49 Å).
  - The matcher's peak memory is unchanged; the workers' record read and the noise-table load are
    below a millisecond; the store transaction is unchanged; D's peak memory fell by 0.48 GB.
  - The orientation search of B and C in these single runs is 1 to 6 % above Phase 0's, within the
    spread this machine showed in this session (117 to 135 s for B across the comparison runs); every
    optimisation was compared against its own reference and none was slower beyond the spread.
- **Simulated workflow tests on the final Release build** (`highlevel_phase3.log`):
  - `simulated_workflow_1jxy`: correlation 0.901, FSC 0.143 at 4.54 Å (Phase 0: 0.899, 4.46 Å);
  - `simulated_workflow_6vxx`: passes in this run at 0.868 and 6.09 Å (Phase 0: 0.877, 5.43 Å). Its
    runs are not reproducible (multithreaded reduction order changes the path from the first
    iteration of the 2D step), and a rerun converged to a wrong map (whole-volume correlation 0.54,
    truth FSC 49.9 Å). Repeated runs (`vxx_ab.log`, `scratch/vxx_ab.sh`) on the final build and on the
    tree as it was before this phase (`scratch/build_p2`, rebuilt from `state/phase_3.start`):
    the final build passed 6 of 8 runs (whole-volume correlation 0.868 to 0.885, truth FSC 4.62 to
    6.09 Å) and failed 2 (0.54 and 0.55, FSC 49.9 Å); the earlier tree passed 7 of 8 (0.867 to 0.888,
    4.71 to 5.31 Å) and failed 1 the same way (0.51, 83.2 Å). The failure, a de novo run converging to a
    wrong map, belongs to the existing code; Phase 3 changes no result except the last bits of gridding
    and restoration (a multiply by a stored reciprocal), and its passing runs are as good as before;
  - `cont_refine3D_1jxy`: passes, ground-truth errors as in Phase 0 (median rotation errors 0.91,
    0.23 and 0.91 degrees, 99.7 to 99.8 % in basin);
  - `solve3D_addon`: passes (union map against truth 0.938, FSC 5.03 Å, cohort pose error 9.7°;
    Phase 2 baseline 0.949, 4.83 Å, 4.25°; floors 0.90, 6.36 Å and 15°);
  - `flex_pca_blobs` (out of the acceptance by the third ruling): its FLEX checks and the blank-half
    check (18 blank halves, 18 flagged) pass; it fails in `refine3D_states_cavgs`, iteration 2, where
    the states collapse again (state 1 holds 341 rows, state 2 one even row and no odd row) and the
    distributed PCG refuses a state with only one half ("distributed PCG requires both halfsets").
    The empty-state rule of the second ruling covers a state without particles, not a one-sided one.
    This is the known failure of the existing code recorded in Phase 2, reached by another path.

**Exit criteria.**

- The build, the fast gate and the registry check pass: above.
- Every optimisation has its before-and-after entry: the table above (the GPU weights cannot run
  on this machine; folding was measured and dropped).
- Cases A to E repeated on the build with every kept optimisation: the acceptance holds.

**For the maintainer.**

- The plane cache is keyed on the padded box and the lowest shell and dropped whenever the spectra
  are reloaded, instead of the store generation (which `euclid_sigma2` does not hold, and which FLEX's
  unit spectra do not have).
- The Cartesian residual keeps its allocating interface (reason above).
- The sparse ranges and the streaming commit are Phase 4 and beyond; the bootstrap's store I/O gain
  only shows at scale (the ribosome measurement is the maintainer's).
- `flex_pca_blobs` stays a known failure: `refine3D_states_cavgs` collapses the phantom's states,
  and a state left with particles in one half only stops the distributed PCG solve.
- `simulated_workflow_6vxx` fails now and then on the existing code too (one run in eight before
  this phase, two in eight after), converging to a wrong map; it is not reproducible run to run.

### Phase 4: the streaming commit, documents and close (2026-10-13)

Validation follows the maintainer's two rulings of 2026-10-13, recorded in section 12.5 step 3.

**What changed.**

- **The streaming commit** (section 6).
  - The strategies' candidate preparation (`prepare_canonical_sigma_update` in the refine2D and refine3D
    strategies) now calls `sigma2_state_prepare_update(committed, nparts, numlen)`. It validates the
    committed store and removes the candidate and range files of the next generation that a crashed
    transaction may have left. It copies nothing.
  - The workers' range write (`write_sigma2`) takes the next generation and the particle layout from the
    committed store's header, so no candidate has to exist while the workers run.
  - The commit (`commit_sigma2`, `consolidate_canonical`) calls `sigma2_state_assemble_candidate` and then
    `sigma2_state_publish`, and deletes the range files after publication.
  - `sigma2_state_assemble_candidate` (`simple_sigma2_state.f90`) checks every range's header (generation,
    layout, shells, bounds), that the ranges lie in the schedule, cover it and overlap nowhere, before it
    writes anything. It then creates a fresh candidate of the next generation and fills it in blocks of
    4,096 rows through one handle, each row from its range or, outside the schedule, from the committed
    store. The valid rows of active particles are added to `real64` sums per half as they pass; the
    global spectra are their means, refused when a half is empty or a mean is not finite. The candidate
    is synced once.
  - New one-handle primitives in `simple_sigma2_state_file.f90`: `sigma2_state_open_candidate`,
    `sigma2_state_write_particles_unit`, `sigma2_state_close_candidate` (global spectra, checksum, one
    sync), `sigma2_state_open_store`, `sigma2_state_read_particles_unit`, `sigma2_state_open_local_range`
    and `sigma2_state_read_local_range_rows`. The path-based readers and writers are built on them.
  - The commit's benchmark file has the lines "streaming assembly", "publication" and "total time".
- **The block-copying merges.** `concatenate_canonical_states` and `merge_selected_canonical_field`
  (`simple_projfile_utils.f90`) call `sigma2_state_concatenate`, which copies the source stores in the same
  blocks into a fresh candidate, sums while copying and syncs once, then publish.
- **The science validation's full re-read moved to the tests** (ruled). `sigma2_state_reduce` stays: it
  is the full re-reading reduction the tests hold the streamed spectra against, and test fixtures use it
  to make a committed store. `sigma2_state_validate_science`, `sigma2_state_merge_local_ranges`,
  `sigma2_state_commit`, the copying branch of `sigma2_state_create_candidate` and
  `validate_copy_compatibility` had no other user and are deleted.
- **Tests** (section 8; `simple_sigma2_state_tester.f90`):
  - `test_global_transaction` on 4 and 9,000 rows (three blocks): the streamed candidate holds the ranges'
    rows, is the next generation with the residual provenance, its global spectra equal a full reduction
    (`require_full_reduction`, the former validation rule and tolerance), publication replaces the store
    and the ranges stay for the caller;
  - `test_unscheduled_rows_from_committed`: rows outside the schedule come from the committed store;
  - `test_recovery_guards`: an uncovered row, a doubly covered row and a range of another generation are
    refused, with the committed file unchanged byte for byte;
  - `test_failed_publication`: a candidate damaged after assembly is not published, the committed file is
    unchanged and the range file is still there; the commit reruns from it and publishes;
  - `test_update_preparation`: the stale candidate and range of the next generation are removed, and an
    update of a store that is not committed is refused;
  - `test_concatenation`: two stores (5,000 and 3 rows) concatenated in order, global spectra equal to a
    full reduction; sources with more rows than the candidate are refused;
  - `test_invalid_record_skip` and `test_estimate_available_on_disk` follow the streaming commit.
  - `simple_cartft_calc_tester.f90` (file-table row added): its fixture and the polish and round-trip
    cases follow the new commit.
- **Documents** (section 9). The forced `sigma_est` lines and the `memory_estimator`, reference-generation
  and planned-note mentions had already gone in Phase 1. Reworded now: "compatible grouped sigma
  estimates" and the manifest's "grouped-section identity" (`solve3D_policy.md`), "group-sigma
  consolidation" and the grouped sigma files (`refine3D_policy.md`), the "sigma-group loader" and
  "grouping policy" (`solve3D_cavgs_policy.md`), and "per-particle sigma2" (`reconstruct3D_pcg_policy.md`).
  Stale grouping comments are gone from `simple_flex_pca_strategy.f90`, `simple_final_rec.f90`,
  `simple_sigma2_bootstrap.f90` and `simple_sigma2_files.f90`, and two found on the way, with file-table
  rows: the `solve3D` manifest's comments and the `rec3D_backends` test description, which still spoke of
  symlinked sigma2 group files. The skill file and the generated indexes are untouched.
- **Files changed in this phase** (20 before the move, all assigned by the file table): the plan; four
  policies; `simple_sigma2_state.f90` and its tester, `simple_sigma2_state_file.f90`,
  `simple_euclid_sigma2.f90`, `simple_commanders_euclid.f90`, the refine2D and refine3D strategies,
  `simple_projfile_utils.f90`, `simple_cartft_calc_tester.f90`, and the five files with stale comments
  and the test UI. Then the plan was moved here and the report written beside it.

**Evidence** (under `/home/elmlundho/agent_runs/global_sigma2_only/scratch`; logs in `logs/`).

- **Builds and gates.** `./compile_debug.sh` builds without compiler warnings, passes the fast gate 15 of
  15 and installs (`p4_debug_build_final.log`; the first build after the commit change,
  `p4_debug_build1.log`). `python3 scripts/check_test_registry.py .` passes (`p4_registry.log`). The
  Release build is clean (`p4_release_build.log`). `git diff --summary` shows only Phase 3's deletion of
  `simple_sigma2_binfile.f90`, no mode change.
- **Unit suites** (Debug, `p4_unit_final.log`): `unit_pftc_align2D3D` (the sigma2 state suite with the
  cases above, and the sigma2 bootstrap suite; 490 checks), `unit_cart_align3D` (744, with the polish and
  round-trip cases on the streaming commit), `unit_project` and `unit_reconstruction` (228) pass.
- **Cases A, B, C and E, once each** (ruled), with the Phase 0 command lines on the Release build, one at
  a time (`queue_phase4.log`, `case_p4*.log`; benchmark files in `keep/bench/phase4/`, condensed with
  Phase 0's in `keep/bench/phase4/comparison_phase0_phase4.txt`):

| Case | Wall time (Phase 0) | Peak memory (Phase 0) | Result (Phase 0) |
| --- | --- | --- | --- |
| A1 | 15:11 (15:00) | 3.21 GB (3.19) | 4.35 / 3.89 Å (4.35 / 3.84) |
| B1 | 2:58 (3:10) | 3.03 GB (3.03) | 8.82 / 4.41 Å (8.82 / 4.41) |
| C | 2:04 (2:11) | 3.69 GB (3.70) | 8.59 / 4.53 Å (8.82 / 4.60) |
| E1 | 2:49 (3:13) | 0.43 GB (0.41) | score 0.389, mean class resolution 11.19 Å; 89 classes, median 11.66 Å (0.389, 10.85 Å; 87, 12.09 Å) |

  - The FSC 0.143 resolutions meet the acceptance: A1 is 0.05 Å (1.3 %) above Phase 0, C 0.07 Å below,
    under the 3 % margin; B1 equals Phase 0's B1. E1 meets the `solve2D` rule (score at least 0.377,
    mean class resolution at most 11.49 Å).
  - By construction the streaming commit sums the same rows in the same order as the reduction did, so it
    writes the global spectra the old commit wrote (not measured bit for bit); the runs differ from
    Phase 0 by the run-to-run variation of multithreaded alignment.
- **The transaction against Phase 0** (seconds summed over the run's iterations: A1 43, B1 and C 5, E1 27;
  matcher lines the slowest partition):

| Line | A1 P0 / P4 | B1 P0 / P4 | C P0 / P4 | E1 P0 / P4 |
| --- | --- | --- | --- | --- |
| candidate preparation | 0.641 / 0.000 | 0.031 / 0.000 | 0.024 / 0.000 | 0.145 / 0.000 |
| range merge + reduction + validation (P0); streaming assembly (P4) | 0.846 / 0.251 | 0.104 / 0.029 | 0.049 / 0.028 | 0.572 / 0.155 |
| publication | 0.254 / 0.265 | 0.027 / 0.030 | 0.031 / 0.027 | 0.174 / 0.177 |
| commit total in the commander | 1.093 / 0.525 | 0.130 / 0.061 | 0.082 / 0.057 | 0.744 / 0.339 |
| commit as the strategy sees it | 1.512 / 0.933 | 0.245 / 0.180 | 0.167 / 0.135 | 0.964 / 0.546 |
| workers' record read | 0.037 / 0.000 | 0.005 / 0.000 | 0.006 / 0.005 | 0.016 / 0.000 |
| workers' range write | 1.225 / 1.270 | 0.146 / 0.087 | 0.084 / 0.079 | 0.224 / 0.312 |

  The workers' range write is the same code in both; its variation is the machine's. Publication is
  dominated by its two syncs (file and directory) and is unchanged.
- **Full passes over the store per matcher iteration** (counted from the code each iteration runs;
  header and global-section reads are a few hundred bytes):
  - before: five reads (the workers' record read, the candidate copy, the range files in the merge, the
    reduction's re-read, the validation's re-read) and three writes (the workers' range files, the
    candidate copy, the merged rows), with seven file syncs at four parts (the copy, one per range merge,
    the global spectra, publication);
  - after: two reads (the workers' record read, the range files in the assembly) and two writes (the
    workers' range files, the candidate), with two syncs (the candidate, publication);
  - the merges of chunk and project stores: one read and one write with one sync, where they read and
    wrote each source with its own sync and then re-read the result twice.
- **Simulated workflow tests** on the Release build, side by side after the timed cases (ruled;
  `highlevel_phase4.log`, one log per test `highlevel_phase4_<test>.log`):
  - `simulated_workflow_1jxy` passes: whole-volume correlation 0.902, truth FSC 0.143 at 4.46 Å (Phase 0:
    0.899, 4.46 Å);
  - `simulated_workflow_6vxx` (recorded, not judged) passes: 0.880 and 5.31 Å (Phase 0: 0.877, 5.43 Å);
  - `cont_refine3D_1jxy` passes, its ground-truth errors as in Phase 0 (median rotation errors 0.91, 0.23
    and 0.91° for stages a, b and c, 99.7 to 99.8 % in the basin; polish 0.629° against 1.035);
  - `solve3D_addon` passes: union map against the truth 0.944, FSC 0.143 at 4.93 Å (floors 0.90 and
    6.03 Å), cohort pose error 10.05° (floor 15°), and an excess of the cohort's pose
    error over the frozen particles' of 4.76° against a limit of 5°. Because that is close to the limit,
    the test ran a second time, alone: excess 3.78°, cohort 9.66°, union map 0.941 and 4.93 Å, passing
    (`highlevel_phase4_solve3D_addon_rerun.log`). Two runs of the same build differ by 1°, so the excess
    varies from run to run; Phase 2 measured 0.50° and Phase 3 2.83°.
  - `flex_pca_blobs` is not rerun (ruled).
- **Disk.** The run directories of the four cases and the recreated input copies were deleted once their
  numbers were recorded; the benchmark files stay in `keep/bench/phase4/`.

**Exit criteria.**

- The build, the fast gate and the registry check pass: above.
- The acceptance of step 3 holds: A1, B1, C and E1 meet the section 12 acceptance; the four high-level
  tests pass; the transaction timings and the passes over the store are compared with Phase 0 above.
- Step 2: the section 9 documents are updated, except the skill file and the generated indexes.
- Step 4: this plan is moved to `doc/refactoring_notes/completed/`, and the report
  `global_sigma2_only_report.md` is written beside it.

**For the maintainer.**

- The candidate now exists only during the commit. Anything that expected a prepared candidate while the
  workers run (none is left in the tree) would have to read the committed header instead, as the workers
  do.
- `sigma2_state_reduce` is used only by tests and test fixtures now. It is kept as the independent full
  reduction the streamed spectra are checked against.
- The bootstrap master (`simple_commanders_euclid_distr.f90`, not assigned to this phase) still writes its
  records through `sigma2_state_write_particles`, one open and sync per block; the one-handle primitives
  would make that a single sync.
- The `solve3D_addon` pose-error excess came within 0.24° of its 5° limit in one run (above).

## 15. Merge into master (2026-10-13)

The run started from `7dbacf1de`. Master had moved on by fifteen commits when
its diff was brought into `~/src/SIMPLE` at `c6e344c6b`, two of them in the
same code. The merge keeps both sides as follows; section 13 records the run
as it ran and is not rewritten.

- **The store carries no checksums** (`b129120f4`, "Removed group and header
  checksums usage in sigma2 states files"). The header and global-spectra
  checksums, the deep validation and the range-file header checksum are gone;
  words 18 and 32 of the state header and word 16 of a range file stay in the
  version-two layout, are written as zero and are ignored on read. The
  manifest identifies the committed sigma2 state by its generation after a
  structural validation. The streaming commit therefore writes the global
  spectra without a header rewrite. In the tests, the retired-word checks of
  that commit are carried over (zero when written, a historical non-zero word
  accepted, a corrupted version refused), and the failed-publication test
  damages the candidate's header instead of its global spectra. The copied-
  project check of `test_estimate_available_on_disk` is kept beside the run's
  `test_new_exec_dir`; both cover the stale recorded directory, which master
  fixed with `update_projinfo` and the run by resolving the store against the
  process's working directory. Master's fix is kept: consumers resolve the
  store against the project file's directory (`update_projinfo`), which equals
  the working directory only when the project file sits there, and the copied-
  project check (a copy in another directory) fails under the run's version.
  `update_projinfo` also names the project after its file, as reading a
  project does, and the particle layout digest includes that name: the
  fixtures of the state and bootstrap testers are therefore named after their
  files (they used free names, which no project read from disk can carry).
- **`cls_expansion` keeps master's conditional sub-class averages**
  (`c8662fb2a`). Only the sigma2 interface is carried over: one consumability
  check and `load_sigma2` without widening the particle range, and the noise
  weight of a member taken from its half's spectrum. The run's
  `w(ncoeff, 2)` table is not carried over: the rewritten module keeps a
  weight per member. The state-0 flags for memberless halves are re-done in
  the rewritten module (the maintainer's ruling stands): a sub-class half
  without members, whose conditional average is only the class's half mean,
  is written with `state_even`/`state_odd` 0 in `cls2D`, counted from the
  hard labels and the half split; the counts travel as two columns of the
  class maps through the distributed merge. The `flex_pca_blobs` gate checks
  the flags against the membership. The commander's own `calc_pspec` rebuild
  was already gone from master.
- **The flex phantom fixture fix** of master (`d03750cd3`, `DEFAULT_NPTCLS`)
  is kept in place of the run's equivalent one-line fix.
- **The bootstrap's average image is the mean of the images summed**
  (the maintainer's decision at the merge). Each `calc_pspec` worker counts
  the images that enter its partition's Fourier sum (sampled, with a spectrum
  that was valid before sanitising) and writes the count into its per-half
  sums file, after the shell count (`SIMPLE_PSPECSUM2`; a count outside the
  partition's range is refused). The master adds the counts and divides the
  summed image by their total, not by the number of active particles. Up to
  25,000 active particles the seed changes only when a sampled spectrum was
  floored; above, the subtracted power is no longer shrunk by the square of
  the sampled fraction. The bootstrap tester checks the count through the
  files and its refusal.

## 14. Code map

| Area | Files |
| --- | --- |
| Parameters | `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` |
| UI | `src/main/ui/simple_ui_params_common.f90`, `simple/simple_ui_refine2D.f90`, `simple/simple_ui_refine3D.f90`, `simple/simple_ui_heterogeneity.f90`, `simple/simple_ui_other.f90`, `single/single_ui_nano3D.f90`, `simple_private_prgs.f90` |
| Sigma2 domain | `src/main/sigma2/simple_euclid_sigma2.f90`, `simple_sigma2_state.f90`, `simple_sigma2_bootstrap.f90`, `simple_sigma2_pspec_sums.f90` (Phase 3) |
| Store file layer | `src/fileio/simple_sigma2_state_file.f90`, `simple_sigma2_files.f90`, `simple_projfile_utils.f90` (`simple_sigma2_binfile.f90` deleted in Phase 3) |
| Commanders | `simple_commanders_euclid.f90`, `simple_commanders_euclid_distr.f90`, `simple_commanders_refine2D.f90` (imports `commander_calc_pspec` and the renamed commander but uses neither: the import goes), `simple_commanders_refine3D.f90`, `simple_commanders_solve2D.f90`, `simple_commanders_solve3D.f90`, `simple_commanders_validate.f90`, `simple_commanders_rec.f90`, `single_commanders_nano3D.f90` |
| Strategies | `simple_rec3D_strategy.f90`, `simple_rec3D_service.f90`, `simple_rec3D_pcg_strategy.f90`, `simple_refine2D_strategy.f90`, `simple_refine3D_strategy.f90`, `simple_calc_pspec_strategy.f90`, `simple_cls_expansion_strategy.f90`, `simple_make_cavgs_strategy.f90`, `simple_matcher_3Drec.f90`, `simple_matcher_ptcl_batch.f90` |
| Polar and Cartesian calculators | `src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_core.f90`, `simple_polarft_corr.f90`, `simple_polarft_memo.f90`, `src/main/cftc/simple_cartft_calc.f90` |
| Image and class averages | `src/main/image/simple_image_ctf.f90`, `simple_image.f90` (`gen_fplane4rec`), `simple_fplane_grid.f90` (Phase 3), `src/main/class/simple_classaverager_restore.f90` |
| Reconstruction backends | `src/main/volume/simple_reconstructor_pcg.f90` (`prep_particles`) |
| Workflows | `simple_final_rec.f90`, `simple_external_reference_pose_initialization.f90`, `simple_refine3D_states_controller.f90`, `src/main/solve/simple_solve3D_manifest.f90`, `simple_solve3D_utils.f90` |
| FLEX | `src/main/flex/run/simple_flex_pca_project_gateway.f90`, `src/main/flex/fit/simple_flex_reconstructor_latent_ops.f90`, `src/main/flex/simple_flex_cls_expansion.f90` |
| Stream | `simple_stream_pool2D.f90`, `simple_stream_stage_initial_analysis.f90`, `simple_commanders_stream_p06_pool2D.f90`, `doc/policies/stream/reference_generation_policy.md` |
| Exec and APIs | `src/main/exec/simple_exec_other.f90`, `src/main/apis/simple_private_exec_api.f90`, `src/main/apis/simple_pftc_srch_api.f90`, `production/simple_private_exec_driver.f90` |
| Tests | `simple_sigma2_state_tester.f90`, `simple_sigma2_bootstrap_tester.f90`, `simple_cartft_calc_tester.f90`, `simple_pcg_halfset_tester.f90`, `simple_flex_pca_application_tester.f90`, `simple_strategy3D_cont_tester.f90`, `simple_strategy3D_cont_1jyx_tester.f90`, `simple_pftc_inplane_tester.f90`, `simple_solve3D_manifest_tester.f90`, `simple_commanders_test_highlevel_frac.f90` |
| Scripts | `scripts/memory/benchmark_solve2d.py` |
