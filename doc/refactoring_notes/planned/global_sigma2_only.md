# Global sigma2 only: removing `sigma_est=group`

Status: refactoring plan, 2026-10-12. No code yet. It is carried out as an
unattended run on the Oracle Linux test machine (section 12). Rulings of
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
    Version 2 (the sigma2 artifact identified by its header and grouped
    section) is not yet in use, so the key leaves without a second bump.
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
| header fields `grouped_offset`, `group_checksum` | `global_offset`, `global_checksum` (same words on disk) |
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
      valid, divided by the number of active particles in both halves.
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
     start, the final reconstruction and the `solve2D` resume.
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
  - So do the mentions in the planned notes `single_pcg_integration.md` and
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

Not part of this refactoring, for a separate scientific decision:

- **The bootstrap's average image** divides the sum of at most 25,000
  sampled images by the number of all active particles. Up to 25,000 active
  particles that is the sample mean. Above, it shrinks the subtracted power
  by the square of the sampled fraction: about 1/900 at 750,000 particles.
  This refactoring preserves it.

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
  Phase 0 value. For `solve2D`, the number of classes and the median class
  resolution meet the same rule. Ground-truth errors lie within their
  Phase 0 values plus the same margin.
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
| `src/main/commanders/simple/simple_commanders_rec.f90` | The half's spectrum; renamed loader. | 2 |
| `src/main/commanders/simple/simple_commanders_euclid_distr.f90` | Bootstrap assembly without group arrays; per-half sums. | 1 2 3 |
| `src/main/simple_final_rec.f90`, `simple_external_reference_pose_initialization.f90`, `simple_refine3D_states_controller.f90` | Policy plumbing; rebuild entry; stale comments. | 1 2 4 |
| `src/main/solve/simple_solve3D_utils.f90` | Policy plumbing and imports. | 1 2 |
| `src/main/solve/simple_solve3D_manifest.f90` | The manifest key, within version 2. | 1 |
| `src/main/sigma2/simple_sigma2_bootstrap.f90` | The central rebuild routine; policy arguments; stale comments. | 1 2 3 4 |
| `src/main/sigma2/simple_euclid_sigma2.f90` | Two spectra and derived arrays; STAR helpers, dead fields; residual path; range writes. | 1 2 3 4 |
| `src/main/sigma2/simple_sigma2_state.f90` | Group arrays out; renames; streaming commit. Its tester follows the policy removal in Phase 1. | 1 2 4 |
| `src/fileio/simple_sigma2_state_file.f90` | Fixed policy, provenances, range mode; renamed fields; streaming commit primitives. | 2 4 |
| `src/fileio/simple_sigma2_files.f90` | One identity check; renamed loader; stale comments. | 1 2 4 |
| `src/fileio/simple_projfile_utils.f90` | Merges without group arrays; block-copying merges. | 2 4 |
| `src/fileio/simple_sigma2_binfile.f90` | Deleted with the per-partition bootstrap files. | 3 |
| `src/main/strategies/parallelization/simple_rec3D_strategy.f90`, `simple_cls_expansion_strategy.f90`, `simple_make_cavgs_strategy.f90` | Rebuild entry; policy probe out; the half's spectrum. | 1 2 |
| `src/main/strategies/parallelization/simple_rec3D_service.f90` | The half's spectrum; renamed loader. | 2 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | One spectrum per half and scalar weights; plane cache. | 2 3 |
| `src/main/strategies/parallelization/simple_calc_pspec_strategy.f90` | Per-half sums. | 3 |
| `src/main/strategies/parallelization/simple_flex_pca_strategy.f90` | Stale comment. | 4 |
| `src/main/strategies/search/simple_matcher_3Drec.f90` | The half's spectrum; plane cache. | 2 3 |
| `src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_core.f90`, `simple_polarft_corr.f90`, `simple_polarft_memo.f90` | Spectra by half; precomputed weights; memoised self-term; GPU weights; folding. | 2 3 |
| `src/main/cftc/simple_cartft_calc.f90` | Slots by half. Its tester's store fixture follows the policy removal in Phase 1. | 1 2 3 |
| `src/main/image/simple_image_ctf.f90`, `simple_image.f90` | `gen_fplane4rec` with the cached plane-grid spectra. | 2 3 |
| `src/main/class/simple_classaverager_restore.f90`, `simple_classaverager.f90`, `simple_classaverager_core.f90` | No record read; the half's spectrum; class bucketing. | 2 3 |
| `src/main/volume/simple_reconstructor_pcg.f90`, `simple_reconstructor.f90` | `prep_particles` with one spectrum and weights; comments. | 2 3 |
| `src/main/flex/run/simple_flex_pca_project_gateway.f90`, `simple_flex_pca_application.f90` | Group fallback, range widening and stale message out. | 1 2 |
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
| `doc/policies/3D/solve3D_policy.md`, `refine3D_auto_policy.md`, `solve3D_cavgs_policy.md`, `refine3D_policy.md`, `reconstruct3D_pcg_policy.md` | Section 9. | 4 |
| `doc/policies/heterogeneity/refine3D_states_cavgs_policy.md`, `doc/policies/stream/reference_generation_policy.md`, `doc/how2s/memory_estimator.md` | Section 9. | 4 |
| `doc/implementation_notes/planned/single_pcg_integration.md`, `refine3D_states_cavgs.md` | Section 9. | 4 |
| `doc/refactoring_notes/completed/global_sigma2_only.md`, `global_sigma2_only_report.md` | The plan moved on completion, and the report. | 4 |

## 13. Progress

No phase has run yet.

## 14. Code map

| Area | Files |
| --- | --- |
| Parameters | `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` |
| UI | `src/main/ui/simple_ui_params_common.f90`, `simple/simple_ui_refine2D.f90`, `simple/simple_ui_refine3D.f90`, `simple/simple_ui_heterogeneity.f90`, `simple/simple_ui_other.f90`, `single/single_ui_nano3D.f90`, `simple_private_prgs.f90` |
| Sigma2 domain | `src/main/sigma2/simple_euclid_sigma2.f90`, `simple_sigma2_state.f90`, `simple_sigma2_bootstrap.f90` |
| Store file layer | `src/fileio/simple_sigma2_state_file.f90`, `simple_sigma2_files.f90`, `simple_sigma2_binfile.f90`, `simple_projfile_utils.f90` |
| Commanders | `simple_commanders_euclid.f90`, `simple_commanders_euclid_distr.f90`, `simple_commanders_refine2D.f90` (imports `commander_calc_pspec` and the renamed commander but uses neither: the import goes), `simple_commanders_refine3D.f90`, `simple_commanders_solve2D.f90`, `simple_commanders_solve3D.f90`, `simple_commanders_validate.f90`, `simple_commanders_rec.f90`, `single_commanders_nano3D.f90` |
| Strategies | `simple_rec3D_strategy.f90`, `simple_rec3D_service.f90`, `simple_rec3D_pcg_strategy.f90`, `simple_refine2D_strategy.f90`, `simple_refine3D_strategy.f90`, `simple_calc_pspec_strategy.f90`, `simple_cls_expansion_strategy.f90`, `simple_make_cavgs_strategy.f90`, `simple_matcher_3Drec.f90`, `simple_matcher_ptcl_batch.f90` |
| Polar and Cartesian calculators | `src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_core.f90`, `simple_polarft_corr.f90`, `simple_polarft_memo.f90`, `src/main/cftc/simple_cartft_calc.f90` |
| Image and class averages | `src/main/image/simple_image_ctf.f90`, `simple_image.f90` (`gen_fplane4rec`), `src/main/class/simple_classaverager_restore.f90` |
| Reconstruction backends | `src/main/volume/simple_reconstructor_pcg.f90` (`prep_particles`) |
| Workflows | `simple_final_rec.f90`, `simple_external_reference_pose_initialization.f90`, `simple_refine3D_states_controller.f90`, `src/main/solve/simple_solve3D_manifest.f90`, `simple_solve3D_utils.f90` |
| FLEX | `src/main/flex/run/simple_flex_pca_project_gateway.f90`, `src/main/flex/fit/simple_flex_reconstructor_latent_ops.f90`, `src/main/flex/simple_flex_cls_expansion.f90` |
| Stream | `simple_stream_pool2D.f90`, `simple_stream_stage_initial_analysis.f90`, `simple_commanders_stream_p06_pool2D.f90`, `doc/policies/stream/reference_generation_policy.md` |
| Exec and APIs | `src/main/exec/simple_exec_other.f90`, `src/main/apis/simple_private_exec_api.f90`, `src/main/apis/simple_pftc_srch_api.f90`, `production/simple_private_exec_driver.f90` |
| Tests | `simple_sigma2_state_tester.f90`, `simple_sigma2_bootstrap_tester.f90`, `simple_cartft_calc_tester.f90`, `simple_pcg_halfset_tester.f90`, `simple_flex_pca_application_tester.f90`, `simple_strategy3D_cont_tester.f90`, `simple_strategy3D_cont_1jyx_tester.f90`, `simple_pftc_inplane_tester.f90`, `simple_solve3D_manifest_tester.f90`, `simple_commanders_test_highlevel_frac.f90` |
| Scripts | `scripts/memory/benchmark_solve2d.py` |
