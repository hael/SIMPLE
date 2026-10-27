# refine3D_local workflow proposal

Source status: implementation complete with the post-review stale-mask correction; administrator acceptance and integration remain pending.
Validation status: Linux compilation was reported successful after the latest correction. The local_startup suite passed 369/369 assertions, UI passed 927/927, and the fast gate passed 15/15 entries. The focused NU-handoff rerun after ownership extraction passed in 949.0 seconds before the stale-mask correction. Broader PCG/distributed and real-data validation remain as listed below.

The empty NU evidence envelope and final cFAR observation are recorded and deferred; they are not unfinished implementation. See [the latest validation and deferred-observation record](#deferred-nu-observations-and-administrator-review-handoff-2026-10-27) for evidence and limits.
Initial source review: `master` at `af13bc8286d5810d7e97523dfd4e013aa0c10951`,
19 October 2026. This is the living contract, review, implementation plan and validation record.

Introduce `prg=refine3D_local` as the single-state continuation workflow for continuous Cartesian
pose refinement. Extract the `pose_cont=only` route from `refine3D_auto`, retain the existing
`refine=cont` engine, and start directly from validated saved references and committed sigmas.
The principal saving is removal of unconditional startup particle passes. Keep the final
all-particle reconstruction until its equivalence to an existing map can actually be established.

## Review findings

1. **The continuous-only branch still pays for a fresh startup.**
   [`exec_refine3D_auto`](../../../src/main/commanders/simple/simple_commanders_refine3D.f90#L49)
   recognizes `pose_cont=only`, checks the stored poses and disables registration. However,
   lines 254–282 unconditionally call `calc_pspec` and `reconstruct3D`, then replace `vol1`
   with the new reconstruction. The discovered project/input volume is only the lag-one reference
   of that reconstruction. “USING PROJECT VOLUME” therefore does not mean the first local
   iteration directly reuses it. The unconditional power-spectrum pass also replaces an otherwise
   usable residual sigma estimate.

2. **The base engine already provides the required local search and sigma fallback.**
   [`strategy3D_cont`](../../../src/main/strategies/search/simple_strategy3D_cont.f90) refines
   rotation and shift from each stored pose, preserves state and half membership, and emits
   Cartesian residual sigmas. Both iteration strategies call
   [`ensure_sigma2_for_iteration`](../../../src/main/sigma2/simple_sigma2_bootstrap.f90), which
   seeds image-power sigmas only if the canonical state is not consumable. A separate local
   matcher, convergence rule or unconditional residual-preparation pass is unnecessary.

3. **A new directory needs a deliberate sigma handoff.**
   [`set_sigma2_state_path` / `get_sigma2_state_path`](../../../src/main/project/simple_sp_project_core.f90#L150)
   normally register the store relative to its owning project directory. A copied project alone
   does not carry that store. Simply deleting auto's `calc_pspec` call still causes the base
   strategy to rebuild it in a new directory. Resolve the source before the directory transition,
   validate it, copy the committed store into the destination and register that local copy.
   Never update the source store through an absolute path or a shared writable link.

4. **A merged volume is not a safe half-map handoff.**
   [`parameters` volume setup](../../../src/main/params/simple_parameters_phases.f90#L337)
   copies `vol1` into missing `_even` and `_odd` files for programs whose name contains `refine`.
   [`read_mask_filter_refvols`](../../../src/main/strategies/search/simple_matcher_refvol_utils.f90#L283)
   can also use the merged map for both halves. The new workflow must validate the real input pair
   before those fallbacks can run, including before its first `params%new`. Otherwise a missing
   half can silently become a nominally gold-standard input.

5. **Derived references and metadata need to travel together.**
   NU references are discovered by filenames and existence; the fresh-start NU LP lookup checks
   the even NU file and then reads the project's `lp`. FSCs are read under canonical names in
   the current directory; envelopes also use current-directory names. A `vol1` path alone does
   not transfer this state. File existence does not prove that an NU pair, FSC, mask and LP
   belong to the same reference generation. See
   [`try_nu_project_working_lp`](../../../src/main/strategies/search/simple_matcher_smpl_and_lplims.f90#L139).

6. **Reference preparation is already separate from reconstruction.**
   [`materialize_reprojection_model_from_volumes`](../../../src/main/strategies/search/simple_matcher_refvol_utils.f90#L165)
   prepares `cart_refvols_even.bin` and `cart_refvols_odd.bin` once per iteration. With `refine=cont`
   it generates no polar reprojection model. Reading, cropping, filtering and preparing Cartesian
   references should remain; those operations are not particle reconstruction. Rebuild these
   derived files after each new assembled map or matching-band change.

7. **Saved half-maps cannot initialize trailing sufficient statistics.**
   Gridding trails Fourier sums and sampling densities; PCG trails raw `(B,D)` statistics.
   Without a valid chain, assembly seeds from the current sample, and that iteration's map
   represents only that sample. Auto's startup reconstruction does not request `trail_seed`.
   Thus retaining that startup reconstruction is not itself a solution to fractional continuity.
   See the [trailing contract](../../../.github/skills/simple-frac-update-trailing/references/frac-update-contract.md).

8. **The final reconstruction is not generally a duplicate of the last iteration.**
   [`calc_final_rec`](../../../src/main/simple_final_rec.f90#L36) uses native sampling and all
   reconstruction-eligible particles. Iteration maps can be cropped, fractional or trailed.
   They also use the pre-iteration sigma generation; the new residual generation is committed
   after assembly. PCG final maps have a distinct cold-solve budget. Even a full native-grid
   iteration therefore does not automatically satisfy the final-map contract. If the registration
   box changed, the current final helper explicitly requests a sigma bootstrap, which includes
   two reconstructions and a polar `refine=sigma` pass. That pass searches no poses, but it is
   additional work worth avoiding on the normal reuse route.

## Proposed public contract

```text
simple_exec prg=refine3D_local projfile=<aligned-project.simple> mskdiam=<A> ...
```

The project must contain valid existing `ptcl3D` poses and half assignments. Use
`check_cont_seeds` before any particle pass, random-orientation fallback or half reassignment.
No global registration, CC pose acquisition, probabilistic neighborhood search or state search
runs in this workflow. All pose updates use `refine=cont`; `pose_cont=yes` is not a child setting.

- Keep base `refine3D refine=cont` for low-level callers, including existing solve/stream users.
- Keep `refine3D_auto pose_cont=no|yes`: its discrete workflow and optional continuous polish.
  Retire `pose_cont=only` with a clear message naming `refine3D_local`, without a compatibility alias.
- Initially retain auto's single-state scope, including `state=X` and frozen `m_estimator=flex`
  work-project selection through a shared helper. Do not silently collapse a multi-state project.
  Selection changes the row layout and may invalidate sigma/accumulator reuse.
- Keep Euclidean scoring, ML regularization, half-specific references, `center=no`,
  `combine_eo=no`, and the current NU/envelope defaults. Reject overrides that violate the local
  workflow, rather than silently replacing them. Keep `refine` internal and reject discrete-only
  options such as `regpass=yes`, `ref_pose_init=cc`, `nspace` and `nspace_sub`.
- Expose `athres_cont`, `trs`, `cont_conv_rot`, `cont_conv_shift`, `cont_conv_frac`, `nsample`,
  `minits`, `maxits`, filtering, masks and execution/backend controls. Honor an explicit `trs`;
  derive its default only when absent. The rotation/shift bounds apply around the pose at each
  iteration, not around the original run-entry pose; cumulative movement may exceed one bound.
- Retain the present sampling budget (25,000 by default), minimum three iterations and automatic
  four-updates-per-particle planning for the initial extraction. Preserve the existing explicit
  `maxits` behavior of also setting `minits`. Tune those separately after validation.
- Recommended new default: `autoscale=no`. Local high-resolution finishing uses native sampling
  unless the user requests a crop. This avoids automatically introducing the crop-to-native
  final sigma-bootstrap cycle. Expose `autoscale=yes` deliberately and report its finalization
  cost. This is a policy choice for the new command, not a change to auto's default.
- Prefer independent half-map modes (`nonuniform`, or the existing half-specific non-NU path).
  For this initial contract, reject `nonuniform_lpset`, explicit `lp` that activates merged-reference
  matching, and `combine_eo=yes`; offer `lpstop` to cap the band. Those existing merged-reference
  semantics remain available through base `refine3D` and need a separate decision if wanted here.

Starting a local workflow from a previous project's result is a fresh invocation, supported in
both execution modes. Do not implement it by forcing `continue=yes`: shared-memory base
refinement rejects that flag, while distributed continuation has additional chain requirements.
True interrupted-run resume retains its existing separate restrictions and must not repeat
startup seeding or reset convergence/sampling state as if it were a fresh invocation.

## Reference and sigma handoff

Choose the effective grid before planning work. Native sampling remains project-owned; pass only
an intentional `box_crop` to children. A compatible finer map can be Fourier-cropped without
reading particles. Validate physical extent, sampling, dimensions, coordinate frame and half
identity and amplitude convention before doing so; `read_and_crop` itself is not a provenance check. Never treat padding
a coarse map as recovered high-frequency information: cap its initial matching band accordingly.

Use explicit `vol1` when supplied, otherwise the selected project's registered unsharpened `vol`
and genuine half pair. An explicit input is a trusted aligned reference set, not a request for
pose acquisition; require its relationship to the project's half assignments. Do not run the
existing external-map normalizer on ordinary SIMPLE continuation products. Sharpened and
low-pass display products are not the primary source.

| Available material | Startup action | Particle work before iteration 1 |
| --- | --- | --- |
| Valid base pair, FSC, compatible NU pair/envelope/LP and sigmas | Stage/reference the bundle; prepare Cartesian references | None |
| Valid base pair and sigmas, missing or unverified derived products | Regenerate needed FSC/mask/NU products through volume-domain owners | None |
| Valid reference pair, missing or incompatible sigma store | Use `ensure_sigma2_for_iteration`; first continuous pass estimates residuals | One image-power pass only |
| Missing half, unreadable base pair or incompatible physical frame | Stop with the exact reason and preparation command to run | None silently |
| Deliberately changed reconstruction estimator/support | Reuse compatible seed references where valid; require reconstruction for products claiming the new estimator | Explicitly justified |

Missing NU products alone must never trigger a particle reconstruction. Factor a narrow
volume-domain preparation entry from existing half-map diagnostics and NU/mask owners if needed.
It must not attempt to recover unregularized halves from regularized ones, reapply an ML prior,
or pretend that multiplying a mask reproduces a PCG constrained solve. Where `_unfil` inputs
needed by the requested NU recipe are unavailable, use a clearly logged compatible half-map
fallback or fail with the missing requirement; do not silently change the scientific recipe.

The handoff consists of the merged map, base and (when required) `_unfil` half pairs, FSC with
grid/estimator identity, required mask/support provenance, and compatible NU products with their
matching LP. The canonical sigma store is validated separately by `canonical_sigma2_consumable`:
committed format, native grid, shell coverage and ordered particle/stack layout. Its current
`ptcl_layout_digest` covers lineage, stack references and row indices, not half assignments or
active membership. The handoff must separately verify preserved halves and compatible population
and noise-estimation settings; the existing validator alone does not establish those properties.
Source and destination must both pass validation after a store copy. A different `nparts` alone is not a
reason to discard a canonical store. Do not copy pending transactions or worker range files.

Record a small reference-set manifest at production/publication: generation, source project/state,
geometry, constituent file identities, half/pose/population provenance and relevant filter/support
settings. This is specific handoff metadata, not a general cache framework. Publish the record
last and reject mixed or incomplete sets. Reference and sigma generations are distinct: a valid
lagged reference need not have been reconstructed with the sigma generation that the next search
will consume.

For old projects without that record, accept their registered, geometrically valid base pair as
a legacy starting reference with that limitation logged. Do not infer derived-product freshness
or final-map equivalence from filenames or timestamps. Rebuild unverified derived products from
the chosen base inputs. Capture all source paths before `params%new` copies/chdirs; copy writable
artifacts into the run and preserve the source. Stage under unambiguous stems so stale sidecars
cannot be picked up accidentally, including with `mkdir=no`.

## Iteration and final reconstruction policy

Run the existing shared-memory/distributed strategy and `strategy3D_cont`. Preserve one particle
read per active matcher batch, its inline partial reconstruction, publication of updated poses and
sampling before file-backed assembly, assembly-owned filtering, and the post-assembly sigma commit.
Reference reuse removes startup work; ordinary self-consistent refinement still needs updated maps.
Do not turn the workflow into fixed-reference fitting by setting `volrec=no`.

For fractional refinement, carry a trailing chain only when its backend, grid, particle layout,
half assignment, state/weight identity and producer generation are compatible. Existing chain
validation is necessary but should not be assumed to establish every new cross-workflow invariant.
Do not manufacture a chain from finished half-maps or combine gridding and PCG accumulators.

Recommended behavior without a reusable chain: make the **first local refinement iteration a full
particle update**, assemble and seed its full-weight chain, then use the requested fractional budget.
This uses the already available reference to improve poses while collecting the necessary statistics,
instead of adding a separate reconstruction of unchanged poses. It is an explicit sampling-policy
change from the old route; log it, count it as iteration 1, and test its effect. Reset run-local
sampling/convergence bookkeeping consistently in both modes; do not let an old `sampled` generation
count toward new continuous convergence. Without this change, document that the first fractional
map discards the unsampled population's contribution.

Initially call the existing `calc_final_rec` once. On the native-grid reuse route, retain its direct
reconstruction using committed sigmas; no extra bootstrap is needed when its existing checks pass.
Use its current fallback when they fail. “Exclusively continuous” describes pose refinement: the
fallback's fixed-pose `refine=sigma` is not discrete pose search, but does currently build a polar
model. Do not claim the entire fallback is Cartesian-only. A future fixed-pose Cartesian residual
path can remove that dependency without abusing `refine=cont`, which would move poses again.

Defer skipping the final reconstruction until an explicit equivalence check covers all eligible
particles, final poses/weights, half partition, native grid, sigma generation, backend/solver budget,
regularization and support. FSC, unfiltered halves and final postprocessing products must belong to
that same result. Full sampling and matching filenames are insufficient. Any optimization belongs
in `simple_final_rec`, where output registration and postprocessing can be shared, not copied into
the new commander. A partial or trailed reference must never be relabeled a final all-particle map.

## Implementation sequence and ownership

1. **Extract the workflow boundary.** Add `commander_refine3D_local`, UI registration and execution
   routing beside auto. Share only genuinely common sampling and state-work-project preparation;
   do not duplicate the entire auto commander. Update parameter validation, help, UI tests and
   migrated `pose_cont=only` callers together. Selecting `prg=refine3D_local` is the explicit opt-in;
   ordinary auto/base calls retain their behavior apart from the intentional removal of `only`.
2. **Implement validated startup reuse.** Add the reference preflight before parameter-side half
   synthesis, capture source paths before directory creation, stage the bundle and committed
   sigma copy, and delegate missing-sigma recovery to its existing owner. Add the narrow
   volume-only derived-product preparation path and publication metadata. Preserve ownership:
   commander selects actions; project/sigma helpers manage their artifacts; volume modules do
   diagnostics/filtering; strategies run iterations.
3. **Close fractional continuity and execution parity.** Validate/copy compatible complete chains;
   otherwise execute the counted full first iteration and seed its chain on either backend.
   Preserve one sample across distributed workers and align fresh sampling resets in both modes.
4. **Use the shared ending and validate.** Retain `calc_final_rec`; exercise native and cropped
   routes. Update the base/auto policies and introduce the local policy when code implements it.
   Keep this note as the implementation and validation record. Final-map reuse is a subsequent
   optimization gated by the equivalence criteria above.

Use lean `only:` imports, allocatable heavy components and one parameter construction per commander.
No whole-`parameters` snapshot, new subset CMake target or separate test executable is needed.

## Validation and measurable acceptance criteria

Extend the existing unit suites and
[`cont_refine3D_1jxy`](../../../src/main/commanders/test/simple_commanders_test_highlevel_cont.f90).
Its auto-only stage and discrete-to-continuous stage are direct migration targets. Keep tests of
auto's `pose_cont=yes` polish and base `refine=cont`. Use fixed simulations and truth-based
pose/map assertions through `simple_test_utils`/the existing gate, not successful completion alone.

- **Fast-path accounting:** with a valid native reference bundle and sigma store in another run
  directory, zero startup `calc_pspec`, `bootstrap_rec3D` or `reconstruct3D` calls; zero discrete
  search/probability-table calls; one Cartesian reference preparation and one map assembly per
  refinement iteration; one justified final reconstruction. Assert calls, not only log text.
- **Integrity:** missing half fails before parsing can synthesize it; wrong extent, reordered rows,
  changed halves/weights, mixed NU generations and interrupted manifests take the stated failure
  or rebuild path. Source projects, maps and sigma stores remain unchanged after the new run.
- **Derived-only recovery:** missing NU/FSC/mask products invoke volume operations, with no stack
  reads. Reject incompatible support provenance; never count duplicated merged maps as two halves.
- **Sigma:** valid cross-directory and changed-partition reuse; corrupt/missing/layout-mismatched
  stores seed exactly once. First search consumes the expected generation, assembly uses the
  lagged model, and commit makes the next generation available.
- **Sampling:** a valid chain allows immediate fractional updates; absent/invalid chains cause the
  counted full first iteration and subsequent requested fractions on gridding and PCG. Verify
  realized fractions, persistent coverage, convergence streaks and no stale generation reuse.
- **Scientific behavior:** perturbed poses improve against simulation truth; discrete-to-local
  continuation preserves prior registration quality; state/half identities stay fixed; rejected
  transactions cannot count as stable. Exercise explicit and automatic iteration limits.
- **Finalization:** native, cropped, fractional, weighted-state and backend cases produce the
  existing native final products with correct FSC/support/sigma provenance. Test cropped sigma
  fallback explicitly; no claim of final-map reuse until the full equivalence gate exists.
- **Performance:** compare end-to-end wall time and particle-pass counts with the old continuous-only
  route on the same project/settings, then evaluate the proposed native default separately.
  Use existing `BENCH` reports for preparation, matching, assembly and sigma costs. Report measured
  savings and peak memory; do not infer a speedup percentage from removed calls.

Initial review performed by source inspection only. At that review, no Fortran source was changed, and no compilation,
linking or test binaries were run. Compilation and runtime validation remain with the maintainer, as required by the
repository instructions; the checks above are implementation acceptance criteria, not observed results.

## Implementation record: first preparatory change (2026-10-23)

User authorized the minimum unambiguous implementation pieces, with unresolved
policy choices reserved for review. Work is on feature/refine3d-local.

Implemented, not compiled: the existing continuous-only auto route now keeps an
explicit parsed trs instead of overwriting it with the automatic shift bound.
An absent trs still uses the existing derived default. Other auto routes retain
their previous behavior. This prepares the bound contract for extraction; it
does not introduce refine3D_local or retire pose_cont=only yet.

Validation: source inspection and diff whitespace check only. This worktree has
no build directory; compilation and runtime regression remain outstanding.
Required runtime cases: continuous-only explicit positive trs, explicit trs=0,
absent trs, and unchanged discrete auto controls.

Remaining implementation: workflow extraction and UI, reference preflight,
committed sigma handoff, derived-volume preparation, fractional continuity and
shared finalization. Review before proceeding: the proposed native-grid default,
reference-set manifest format, and supported restart/accumulator handoffs. No
hybrid search, scientific runs, commits or pushes were performed.

## Implementation record: reference preflight (2026-10-23)

Step 2 now has an early reference preflight on the existing `pose_cont=only`
route. The `refine3D_local` command has not yet been extracted. The preflight
runs before the first `params%new`, so parameter initialization cannot hide
missing input halves by copying the merged map.

The preflight reads the source project and merged/even/odd MRC headers. It
requires existing files, matching cubic dimensions and sampling, and physical
extent equal to the particle box times native sampling. A compatible coarser
grid is accepted. It reuses the continuous seed validator for positive state
and projection, finite angles and shifts, and explicit binary half labels.
Both half-sets must contain active particles. With `state=X`, these particle
checks cover only the selected state. No source project or reference is written.

Existing files establish that the pair predates startup. Headers cannot prove
independent half-map reconstruction or a common coordinate frame. Those remain
input provenance requirements. No map-content independence claim is made.

Regression assertions were added to the existing Cartesian strategy unit suite
(`test=unit_cart_align3D`): compatible native/coarse geometry, mismatched sizes,
sampling and extent, invalid/non-finite sampling, missing pose/half fields,
fractional half labels, and valid zero angles/shifts.

Observed validation: test registry, FLEX dependency DAG, source-header checks,
and diff whitespace check passed. Debug configuration with tests enabled passed.
Compilation stopped in unchanged `src/fileio/simple_posix.c` because Windows
pthread types/functions are unavailable. No test binary, fast gate, runtime
startup rejection case, or scientific dataset validation has passed for this
change. Logs are local `local-preflight-configure-short.log` and
`local-preflight-build-short.log`.

Build-only workaround: the MSYS generated make rules treated quoted
`SIMPLE_PATH=...` as a command. Only the two ignored generated `build.make` files
were changed to invoke `env`; repository CMake files were not changed. The short
Windows path was used for this build attempt. Do not treat that workaround as
validation of the standard unmodified Windows build.

Still required after the build blocker is resolved: run the Cartesian unit
suite and fast gate; exercise valid startup, each missing half, each geometry
mismatch, invalid stored seeds/half labels, selected-state validation, and
unchanged discrete auto startup. For each rejected input, verify nonzero exit
and no synthesized half-map or source-project write. Workflow extraction,
sigma reuse, reference preparation and finalization remain pending.

## Implementation record: local workflow boundary (2026-10-23)

Step 1 is implemented in source; Linux compilation and runtime acceptance are
pending. The user requires all compilation and executable tests on Linux.
Do not attempt Windows compilation again. The prior Windows failures above
are historical evidence, not a prerequisite for further source work.

Added `commander_refine3D_local`, `simple_exec` routing (including the existing
restart wrapper), and the `refine3D_local` UI/help registration. The two thin
commanders call one private `exec_refine3D_workflow` body with an explicit
continuous-only flag. Shared startup, sampling, state preparation, iteration
and finalization were retained without copying the auto workflow. The UI also
uses one input builder; local exposes bounds/convergence controls and does not
offer discrete registration, pose acquisition or a mode switch.

Auto retains `pose_cont=no|yes`. The commander and parameter validation reject
`pose_cont=only` with a message directing the caller to `refine3D_local`.
Local derives `refine=cont`, rejects projection-grid arguments and conflicting
workflow overrides, and requires `state=X` for multistate input. Its preflight
still precedes the single parameter construction. The preflight and work
project use one `select_state_rows` helper, so hard-label and frozen-FLEX
selections cover the same source rows. The continuous seed checker accepts an
optional row mask; its other callers retain their all-active-row behavior.
State work projects use the owning workflow's name.

Migrated both continuous-only stages of `cont_refine3D_1jxy`; retained their
existing metric names and directories for comparison. The perturbed-seed
stage now supplies its existing synthetic truth reference pair explicitly,
since that seed project has no registered output volume. Identical synthetic
reference halves are a fixture, not evidence of independent experimental halves.
The auto discrete/polish stages and base continuous stages remain in the gate.
Updated UI assertions and added a selected-row seed assertion to the existing
Cartesian suite. No CTest entry, executable or target was added.

Updated the auto/base policies and added the current local command policy.
The UI-default audit traces the delegated workflow, retaining program-specific
default generation. Startup `calc_pspec` and reconstruction passes still run;
saved-sigma reuse, derived-only reference preparation and fractional-chain
reuse remain step 2/3 work. The proposed native-grid default was not adopted.

Observed source-only checks: test registry, FLEX DAG, Fortran file headers,
diff whitespace, and UI-default generation passed. UI-default generation
reports 630 baselines and 110 program override sets. No compilation or test
executable was run for this extraction. Preserve the user's spacing.

### Linux handoff for this package

Sync this worktree, including the new local policy and all modified source,
test and UI files. Use the normal Linux environment and existing build.
Run the build first; run executable checks only after the build succeeds:

```sh
make -C build -j8
scripts/run_fast_gate.sh build
build/bin/simple_test_exec test=unit_ui
build/bin/simple_test_exec test=unit_cart_align3D
build/bin/simple_test_exec test=cont_refine3D_1jxy
build/bin/simple_test_exec test=flex_pca_blobs
```

The last gate also exercises the shared auto state/FLEX preparation. Run the
long gates in the background with separate logs. Retain nonzero exit status;
passing source checks is not evidence of a successful Linux build or run.
After these pass, run focused CLI cases for local routing/help, the retired
switch's migration error, conflicting local options, and the preflight failure
cases listed above. Check explicit/zero/absent `trs` and matched auto no/yes
controls. State/FLEX local continuation needs a dedicated runtime case as well.
Then test real data before any scientific or performance claim.

## Linux fast-gate review and seed corrections (2026-10-23)

Observed the supplied `run_fast_gate.log`: 13 of 15 suites passed in 48.01 s.
`unit_ui` passed. `unit_cart_align3D` stopped in `test_selected_seeds` while
copying an uninitialized `ori` (unallocated character hash). `unit_stream`
failed its assertion that assigned rows are valid continuous seeds.
These are failed validation gates; the local workflow is not runtime accepted.

Source tracing identified the stream failure: packed particle records always
have an `eo` slot, but `ori%isthere('eo')` uses nonzero presence semantics and
reports false for the valid even-half value 0. The new presence check must
apply only to hash-backed records. Raw finite/binary validation remains in
place for both representations. No stream code or generic orientation API
was changed.

Corrected the validator accordingly and initialized the selected-row test's
seed with `seed%new(is_ptcl=.true.)` before copying it. Added packed-record
assertions accepting both 0 and 1 and rejecting fractional/unassigned labels.
The existing hash-backed missing-label assertions remain. No source spacing
was reformatted. Source checks were rerun; compilation and executable tests
remain Linux-only and the fixes await that rerun.

Next Linux gate: sync these corrections, rebuild incrementally, and rerun
`scripts/run_fast_gate.sh build`. `unit_cart_align3D` and `unit_stream` must
both pass before proceeding to the high-level checks in the handoff above.

## Passed Linux results reviewed on 2026-10-26

The workflow-boundary phase is accepted against the supplied Linux logs.
These are observed test results, not Windows execution or real-data evidence.

| Evidence | Observed result |
| --- | --- |
| `run_fast_gate.log` | All 15 fast suites passed; 36.55 s wall time |
| `SIMPLE_TEST_unit_ui_20261023/simple_test_unit_ui_report.txt` | 894/894 assertions passed across five groups |
| `SIMPLE_TEST_unit_cart_align3D_20261023/simple_test_unit_cart_align3D_report.txt` | 808/808 assertions passed across five groups |
| `cont_refine3D_local_latest.log` | Final truth-validation PASS and normal stop; 2104.8 s |
| `flex_pca_blobs.log` | Final shared/distributed workflow PASS and normal stop; 2081.9 s |

The final two logs identify Git commit `c97106c4`; this stamp alone does not
identify uncommitted changes. The continuous log separately confirms the new
`REFINE3D_LOCAL` route for the perturbed and discrete-to-local stages. Its
perturbed local stage reached median rotation error 0.921536 degrees and shift
error 0.0888714 native pixels and passed iteration/sampling checks. The
continuation passed its symmetry-aware bound. FLEX passed state 1/3 weighted
and hard-label work-project checks and confirmed byte-identical parent project
content. These observations support boundary extraction and shared-helper
regression acceptance; they do not establish real-data quality or speedup.

The earlier incomplete continuous log used the old auto label and is not
acceptance evidence for the extracted local route. The final logs supersede
that snapshot and the earlier failed fast gate. The successful test fixture
directories were removed by normal test cleanup; reviewed logs were downloaded
to `C:/Users/hossainm7/Downloads/`. No map-level reanalysis is claimed.

Next package begins with focused startup-rejection assertions. Saved-reference
and committed-sigma reuse, first-full fractional continuity and real-data
validation remain pending. Preserve existing user whitespace; do not compile
or execute SIMPLE on Windows. Scientific runs must use the Linux environment.

## Next-package contract (2026-10-26)

User approved `autoscale=no` as the local default; auto retains `yes`. Explicit
caller settings take precedence. The focused rejection suite runs the real
local commander in forked test-dispatcher children, using synthetic headers and
saved poses but no particle stack. Each negative case must fail for its intended
reason before execution-directory creation, leave a missing half absent, and
preserve every source-project/map byte. Register it within the existing
Cartesian area suite; do not add a test executable or CTest entry.

Before removing startup work, preserve the reference/sigma handoff and trailing
contracts above. Do not reuse unverified NU sidecars or copy uncommitted sigma
transactions. Real-data project and scheduler details have been requested; no
real-data run or scientific result is implied by source implementation.

## Next-package implementation and Linux handoff (2026-10-26)

Recorded accepted boundary-phase results before these changes. The new code
is not covered by those earlier passed logs.

Implemented the focused `local startup` sub-suite of `unit_cart_align3D`, beside
the owning commander. It uses the existing Linux test dispatcher to invoke the
real local commander in isolated children. Twelve cases cover each missing
half, dimensions, sampling, extent, invalid projection, nonfinite angles and
shifts, fractional/unassigned halves, a single populated half and no active
particles. Assertions require the intended diagnostic and nonzero failure,
no newly created execution directory, unchanged input existence, and bytewise
unchanged source project/maps. Successful fixtures are cleaned; assertion
failures retain the fixture and child logs. No particle stack is created, so
these tests must not enter a particle pass. No new CTest entry or executable
was added. Run this sub-suite before the area's numerical suites.

Applied the user-approved local `autoscale=no` default in its thin commander;
auto keeps `yes`. Both retain explicit caller values. UI assertions cover the
generated defaults. The shared workflow no longer sets a competing default.

Prepared `scripts/submit_betagal_local.sh` from the existing data-testing Slurm
resource pattern: one node, local partition workers, 10 parts x 8 threads,
128G, with overrides. This is a short two-iteration continuation from an
existing aligned beta-gal project/reference pair, not movie processing or a
full fresh reconstruction campaign. It accepts exact input paths, refuses an
existing output directory, records command/executable/source-input hashes,
checks source inputs on success or failure, and retains all job/arm logs and
outputs. Optional `BASELINE_BUILD` runs the same local command against a saved
build from after boundary extraction but before startup-reuse changes. Without
that control it is a smoke run, not a speed comparison. It does not interpret
FSC or claim scientific acceptance. No job was submitted and no SIMPLE binary
was run while preparing it.

Observed checks: registry, FLEX DAG, leading Fortran descriptors, diff
whitespace and source-only UI-default generation passed. The generated lookup
confirms auto=yes and local=no. `bash -n` passed for the Slurm launcher; this
checks syntax only. New Fortran has not been compiled or executed on Windows
or Linux. Linux needs a CMake regeneration for the added tester source:

```sh
cmake -S . -B build
make -C build -j8
build/bin/simple_test_exec test=unit_cart_align3D suite=local_startup
build/bin/simple_test_exec test=unit_ui
scripts/run_fast_gate.sh build
```

Review the focused failure diagnostics before repeating the long scientific
gates. The generated `suite=local_startup` help and the existing area-suite
registration are updated together. Do not rerun the broad FLEX campaign just
for these additions.

Startup-reuse source audit: `canonical_sigma2_consumable` already validates
committed format/native grid/shell coverage/ordered layout. It does not cover
half/population provenance, so copying a registered path alone is insufficient.
The current local route still calls `prepare_external_init_vol`, which
normalizes an explicit pair, then unconditionally seeds sigmas and reconstructs.
That route must be replaced together with volume-only FSC/NU preparation and
source-path capture; merely skipping `calc_pspec`/`reconstruct3D` would leave
an inconsistent reference handoff. No such skip has been implemented here.
Shared-memory refinement already resets sampled/updatecnt at fresh start;
distributed parity and the full-first-iteration/trailing-seed transition still
need implementation and Linux checks. Saved-reference/sigma reuse and
fractional continuity remain the next source package after focused-test
validation, not accepted features. Real-data validation is prepared but not run.

Once that reuse package is ready and exact Linux inputs are selected, the
launcher entry is:

```sh
SIMPLE_PATH=/absolute/candidate/build bash scripts/submit_betagal_local.sh \
  /absolute/aligned/betagal.simple /absolute/aligned/recvol_state01.mrc \
  /absolute/new/betagal_local_check
```

Load the cluster's normal runtime modules first. The launcher defaults to
`csbdevel`, following the independent dataset launcher; override partition,
account, time, memory and worker counts through its documented environment
variables. Its source-map headers are validated by the real command; filename
checks alone do not establish aligned or independent experimental halves.

## Startup reuse and fractional continuity implementation (2026-10-26)

Source implementation is complete for the requested startup and first-update
paths; current changes have not been compiled or run. Auto retains its existing
startup sequence.

The local preflight now captures absolute map, FSC and committed-sigma paths
before `params%new` can create the run directory. After parameter setup, local
copies the merged/even/odd maps and support-provenance sidecar into its run
files and updates only the run project. A compatible registered FSC is copied;
otherwise it calculates FSC directly from the saved pair. It does not normalize
or write the source references. The matcher builds its Cartesian references
from the staged pair. Missing NU reference files use the existing FSC-filtered
half-map fallback for the first iteration; volume assembly produces the next
FSC, envelope-mask and NU outputs.

The source sigma file is staged only when its committed format/grid/shell/layout
checks pass both for the source project and the selected work project. The run
project points to its staged copy. Existing strategy initialization uses
`ensure_sigma2_for_iteration` to invoke `calc_pspec` only when that staged state
is absent or unusable; no startup reconstruction is run for local.

For fractional trailing, local checks the saved manifest, row and state layout,
weight identity, grid extent and component sizes before copying a complete
chain, with the manifest last. If that check fails, the existing `trail_seed`
handshake makes iteration 1 full in both shared and distributed execution; the
requested fraction is restored for later iterations. The existing terminal
`calc_final_rec` remains shared and unchanged. `cont_refine3D_1jxy` now asserts
that the no-chain first full update reaches every particle.

Validation observed in this worktree: `git diff --check` passed. No Windows
compile or executable test was run. A single Linux package remains required:
regenerate CMake for the new focused startup tester, build, run the focused
`local_startup` suite, `unit_ui`, the fast gate, and `cont_refine3D_1jxy`. After
that package passes, run the prepared short beta-gal Slurm continuation. The
current source has no observed Linux result for sigma reuse, chain reuse,
first-full sampling, or real data.

## FSC self-copy failure and temporary focused gate (2026-10-26)

Observed evidence: saved fixture `test_cont_refine3D_1jxy_4192080` retains
`h_continuation/fsc_state01.bin` at zero bytes. Its project registers the same
absolute FSC path that local startup uses as its destination with `mkdir=no`.
`simple_copy_file` opens its destination with `status="replace"` before reading
the input. Copying that FSC onto itself truncates it and causes the byte-buffer
read failure. This is a startup regression, not a scientific acceptance failure.

Required outcome: preserve an existing FSC when source and destination resolve
to the same path; regenerate missing or empty FSC from the validated half-maps.
Do not change the common copy routine or source fixture files. Temporarily run
only the local workflow in `cont_refine3D_1jxy`, including a second call in the
same directory without an explicit reference. Keep the full stages in source
for restoration after local validation passes. A focused pass is not a pass of
the complete gate.

Implementation approach: guard the FSC copy using the same absolute-path
comparison as map staging; use the existing half-map FSC fallback for an empty
registered file. Add a test-only `LOCAL_WORKFLOW_ONLY` flag, retain stage e and
its existing assertions, then continue from its saved project in place and
assert particle count, finite pose errors, iteration completion and nonempty
FSC. Leave other stages behind the flag. No Windows build or SIMPLE test run.
Implemented in the feature worktree: the FSC copy now skips identical absolute
paths and regenerates absent/empty registered FSC through the existing half-map
operation. The test's temporary `LOCAL_WORKFLOW_ONLY=.true.` runs stage e
(four iterations, nsample=500), then a two-iteration in-place continuation using
only the saved project and its registered references. Stages a-c, f-h and i-j
remain in source but are disabled. Set that flag to `.false.` after the local
checks pass to restore the full gate. The completion message explicitly refers
to selected stages, not acceptance of the full scientific gate.

Evidence review confirmed that the saved project's registered FSC resolves to
its own `h_continuation/fsc_state01.bin`; that file is zero bytes after the
failed self-copy. The downloaded fixture was not changed. Static inspection
found no trailing whitespace or continued fatal-macro arguments in the two
edited Fortran files. `git diff --check` could not complete reliably because
Windows returned an access error while launching/reading Git. No compilation or
SIMPLE test was attempted on Windows. Linux compilation, the focused test and
later restoration/re-execution of the complete gate remain pending.
In-place trailing-chain staging is also guarded now: each component and its
manifest is copied only when the absolute source and destination paths differ.
Cross-directory staging still copies all components before publishing the
manifest. This prevents the copy routine from truncating a reusable chain in
an in-place continuation.

Explicit crop handling now preserves a supplied `box_crop` even when
`autoscale=no`; child commands receive that crop and derive `smpd_crop` from
the native project grid. Automatic autoscaling remains unchanged. This source
change has not been compiled or run.

Reference-set provenance and the stronger sigma/chain reuse checks remain
pending. Their producer records must be published by gridding assembly, PCG
assembly and sigma commit, then checked by local staging; accepting a manifest
only in the local commander would not establish producer provenance.

## Observed Linux focused-local result (2026-10-26)

Reviewed `cont_refine3D_1jxy.log`: the selected local stages passed and ended
with `SIMPLE_TEST_CONT_REFINE3D_1JXY NORMAL STOP`, reported Git commit
`c97106c4`, and execution time 1639.8 s (27 min 19.8 s). The working directory
was removed by test cleanup; the retained log is the evidence. The reported
commit does not identify uncommitted source changes in the tested build.

Stage e bootstrapped missing sigmas, performed a full first update, seeded the
1000-particle chain, then used 500-particle updates with 0.5 trailing weights.
Its assertions passed: first-full coverage, iteration budget, sample size, no
polar model/probability table or discrete projection assignment, and recovery
against simulation truth. Median rotation error was 0.437455 degrees, median
shift error was 0.0584891 px, and the reported in-basin fraction was 1.0.

The subsequent in-place continuation reused committed sigmas and passed the
particle-count, finite-pose-error, iteration-budget and nonempty-FSC assertions.
It used full updates (1000 particles against the default target 25000), so it
did not exercise in-place fractional-chain staging. The log does not establish
validation of the latest explicit-crop edit. The complete gate's other stages
remain disabled and require restoration and a later full run.

The run logged empty NU evidence-envelope warnings and cFAR=0.0. These did not
fail its assertions; the focused pass does not resolve those diagnostics or
establish real-data quality. Reference-set provenance, stronger reuse checks,
fractional in-place chain-copy validation and real-data validation remain open.

## Provenance implementation contract (2026-10-26)

Publish a reference-set record after gridding/PCG has completed the state outputs.
Record the producer project/state, iteration, grid, ordered row identity, halves,
active population, poses, observation/CTF identity and filter/support settings.
Bind the record to constituent content digests. Mark publication in progress
before replacing files; an interrupted or corrupt recorded set is rejected.
Legacy base pairs remain admissible with a logged limitation, but their FSC
and derived products are not trusted by local startup.

Sigma commit and both accumulator-chain producers publish matching handoff
records, bound to the committed sigma file or the chain manifest and components.
Local sigma reuse checks rows, halves, population, observations and noise settings;
local chain reuse also checks saved poses and reconstruction settings. Reference
and sigma generations remain independent. Invalid/missing reuse provenance
causes sigma bootstrap or a full first update, never source-file deletion.
An accepted chain preserves represented-particle update history at fresh startup.

Implementation reuses the existing layout digest, content-digest primitives,
chain manifest, canonical sigma validator and fallback owners. Keep these checks
in a small shared handoff utility, with publication calls in the existing owners
and consumption in local startup. Add focused assertions to `local_startup`;
compile and run on Linux only. The user's latest Linux compilation passed before
this package; the 1639.8-second log belongs to the earlier build.

## Provenance source implementation and validation handoff (2026-10-26)

The preceding contract is implemented in `simple_refine3D_provenance` and the
existing gridding, PCG, sigma-commit and local-startup owners. Each anchor has
a versioned `.provenance` record. Producer project/state/iteration and ordered
constituent content identities define the reference generation; sigma and chain
records also bind their existing generation-bearing headers/manifests. The
record retains settings, grid, matching LP and particle identity digests.

Reference producers mark `.provenance.pending` before replacing outputs and
publish after NU/support outputs. Local startup rejects incomplete/mixed
recorded sets. Legacy bases rebuild FSC. Staged imported bases clear stale
stem-specific NU products and carry a `.local_import` marker: the matcher skips
unverified global lagged envelopes until assembly supplies new references.
Imported source maps and sidecars are not removed by staging.

Local sigma staging validates both copies and clears inherited registration
before deciding reuse, so rejection cannot bypass the existing bootstrap.
Local chain staging validates both copies, includes the handoff record before
publishing the manifest last, and sets the existing `trail_seed=no` handshake
on acceptance. Shared-memory initialization then preserves update history.
The common continuation copy also carries provenance and the PCG manifest.
Reference imports can have a different producer project; sigma and chain reuse
require the selected work-project identity. Missing old handoff records take
the documented fallback rather than claiming compatible provenance.

`local_startup` now adds tiny-file assertions for halves, active/represented
population, ordered rows, CTF/noise settings, lagged sigma poses, changed chain
poses/generation, same-size stale FSC, interrupted publication and legacy input.
The focused high-level continuation now requests `nsample=500` and asserts
fractional sampling and preserved update history. Its recovery thresholds and
six-iteration budget are unchanged. The earlier full-update continuation log
does not validate this path. Its six assembly phases totalled 761.4 s; simulation,
matching, bootstrap and two final reconstructions account for the remaining
runtime, without a separately measured breakdown.

Observed static checks: test registry, Fortran descriptions, FLEX DAG and
`git diff --check` passed. The new utility/tester also passed direct whitespace
and fatal-macro scans; a module-use graph check found no cycle through the new
utility. These are source checks only. No compilation or SIMPLE execution was
attempted on Windows. The user's successful latest Linux compilation predates
this provenance package.

Linux validation starts with CMake regeneration (the new module is globbed at
configuration), incremental build and the existing focused suite:

```sh
cmake -S . -B build
cmake --build build -j <jobs>
build/bin/simple_test_exec test=unit_cart_align3D suite=local_startup 2>&1 | tee build/local_startup.log
```

This generated-data check must report PASS and no failed assertions before the
long integration rerun. Then run the fast gate and the focused high-level gate.
The focused high-level log must show validated reference/sigma reuse, compatible
chain reuse at in-place startup, 50 percent continuation samples and the new
history assertions. PCG/distributed integration, full-gate restoration and real
beta-gal acceptance remain separate pending validation gates. No runtime speed
or scientific-quality claim is made for the provenance package.

## Review corrections: required outcome and implementation (2026-10-26)

The user authorized the three findings in the saved worktree review. This
correction implements the existing public and handoff contracts; it adds no
search algorithm or public option. Local must reject explicit `lp` and
`nonuniform_lpset` before parameter initialization. Auto keeps its choices.

A recorded continuation reference must agree with the source project's ordered
rows, half labels, active population and selected producer state. For explicitly
trusted external inputs with different row lineage, retain the existing external
reference contract and log that distinction. Same-lineage half changes are never
an external-input exception. Reference poses may lag; sigma generations remain
independent. Legacy geometry-only admission stays logged.

For a verified compatible NU bundle, stage its NU maps, local-resolution map
and recorded masks with same-file guards, and restore its recorded matching LP
to the run project. Preserve sources. Compare fixed preparation settings, not
the producer's evolving matching LP. Incompatible or unavailable derivatives
keep the logged FSC-filtered base fallback and imported-support marker.

Implementation: extend the existing provenance utility and its commander
consumer, restrict the shared UI builder only for local, and extend the existing
startup/UI suites. Test rejected options and changed reference halves through
the real public commander; test successful derivative copies, settings fallback,
source preservation, saved LP and in-place copies with tiny deterministic files.
Keep the temporary focused-gate flag unchanged. No Windows build or runtime run.

Linux evidence received before these corrections: `local_startup` passed 217
assertions in 3.64 s; the fast gate passed 15/15 suites in 50.29 s. The updated
focused integration is running. These results do not validate these corrections.

Source corrections are implemented. The public commander now rejects the two
merged-reference controls before initialization, and the local UI offers only
`none|nonuniform` plus the existing `lpstop` band cap. Reference validation
compares ordered layout, half and active-state digests and the selected producer
state. A zero-row call to the existing layout digest distinguishes project
lineage from row ordering, so a row reorder cannot use the external exception.

Verified compatible NU derivatives and recorded masks are copied by the volume
provenance utility with same-path guards. Startup restores the recorded LP to
the run particle field without setting explicit `lp` or enabling merged matching.
NU maps already contain their assembly-owned support. The base fallback still
suppresses unverified lagged envelopes. Preparation compatibility excludes the
producer's changing matching LP. The record format is version 2 (separate LP
and lineage digest); version 1 development records require regeneration. No
old-format reader was added. Legacy pairs without records remain admissible.

The existing startup suite now has 16 fatal-rejection cases, including both
merged-reference controls, changed recorded halves and a wrong producer state.
Additional deterministic assertions cover row ordering, active membership,
trusted external lineage, permitted pose lag, NU constituent copies and source
preservation, in-place copies, matching LP, changed settings, unavailable LP
and unrecorded derivatives. The local UI suite checks its filtering choices
and safe band control. The temporary high-level gate flag is unchanged.

Observed source checks after the corrections: test registry, Fortran descriptions,
FLEX DAG, UI-default generation to a temporary file, and `git diff --check` all
passed. No Fortran compilation or executable test ran on Windows. Linux must
compile these changes and run `unit_cart_align3D suite=local_startup`, `unit_ui`
and the fast gate before integration acceptance. The already-running focused
integration belongs to its earlier source snapshot and cannot validate these
corrections. Preserve and assess that result separately.

## Full integration gate restoration (2026-10-26)

Approved outcome: restore all `cont_refine3D_1jxy` stages and retain the new
in-place continuation assertions. Remove the temporary stage switch. Preserve
the original stages' order, thresholds and iteration budgets; append the
continuation regression after them, reading stage e's own saved project before
capturing its update-count baseline.

Observed downloaded Linux results before restoration: `local_startup` passed
331/331 assertions in 4.70 s; `unit_ui` passed 913/913 assertions; the fast
gate passed 15/15 suites in 44.46 s. Its startup failures are expected child
rejections with passing parent assertions. These results include the three
review corrections and the corrected `simple_copy_file` import.

The preceding focused integration passed in 1537.0 s. It verified reference,
sigma and fractional chain reuse, persistent update history and all six
in-place assertions. That run predates the three review corrections and record
version 2; it is not acceptance of the restored complete gate.

Restoration is implemented: the temporary switch and both stage guards are
removed. The original stages run in their previous order; the six in-place
assertions run afterward against stage e's reloaded project. Thresholds and
iteration budgets are unchanged. Test registry, Fortran descriptions, FLEX DAG
and diff whitespace checks passed. The restored complete gate has not been
compiled or run on Windows; Linux incremental build and full integration
execution remain pending.

## Positive NU handoff integration: contract and plan (2026-10-27)

Observed full-gate evidence: the downloaded `cont_refine3D_1jxy.log` passed
all restored stages and six in-place assertions in 5766.5 s. Both continuation
inputs came from completed workflows. With `automsk=yes`, the shared ending
intentionally reconstructs `filt_mode=none`, whereas local defaults to NU. The
reported preparation mismatch and base-half fallback are correct for that
handoff. This pass therefore does not cover positive compatible-NU reuse.

Approved outcome: add one focused production-context integration case to the
existing test entry. Use the existing seeded 1JYX fixture and one base
`refine3D refine=cont` iteration with NU assembly enabled. That base command
publishes the current NU bundle without calling the workflow final ending.
Start one `refine3D_local` iteration from this source project in a separate
execution directory. Keep production defaults, compatibility checks, final
reconstruction behavior, and the existing recovery thresholds unchanged.

Acceptance: the producer record validates; its matching LP and committed sigma
are present; all four staged NU/locres files match the producer contents; the
import-fallback marker is absent; the first local iteration's `LP_MATCHING`
agrees with the producer's `LP_NU_HANDOFF`; all particles receive continuous
work; producer project, reference constituents, producer record and sigmas
remain unchanged. Verify output half labels remain fixed. These assertions
check the handoff, not an additional scientific recovery threshold.

Implementation: reuse the current fixture, stage helper, production commanders,
file digest and statistics fields. Add `suite=nu_handoff` to run this case alone
without the polar/full stages. The default full gate also appends it after the
existing stages and in-place checks, preserving their random draws and budgets.
Share cleanup, distinguish focused/full success messages, and register help in
the existing test UI. Add no test executable, CTest entry or production hook.
Compile and run on Linux only.

Implemented the focused selector and appended case. It reuses the existing
base-continuous stage helper with NU/ML/envelope settings matching local
defaults. The source is captured before calling the public local commander
with `mkdir=yes`, exercising a real cross-directory project/reference/sigma
handoff. The first-iteration check compares `LP_MATCHING` with the producer's
`LP_NU_HANDOFF` at the actual Fourier grid. Constituent digests, the absent
import-fallback marker, positive continuous work, fixed halves, no discrete
artifacts and source preservation are assertions, not activation-log checks.

The selector runs after shared simulation/seed construction and before any
original full-gate stage. The default full gate appends the same case after
all existing assertions. Cleanup is shared; focused success is labelled as
focused. The test UI and policy describe `suite=nu_handoff`. No production
Fortran, CMake target or registry entry was changed.

Observed static validation: test registry, Fortran descriptions, FLEX DAG and
`git diff --check` passed. No Windows compilation or test execution was done.
Linux incremental compilation and the focused integration remain pending.
The earlier 5766.5-second full pass predates this new integration case.

After syncing and incrementally building on Linux, run only the new case:

```sh
simple_test_exec test=cont_refine3D_1jxy suite=nu_handoff 2>&1 | tee build/nu_handoff.log
```

Expected completion: `PASS: focused NU-bundle handoff integration` and normal
stop with every `nu_handoff_*` assertion passing. Keep the log for review.
This does not validate PCG/distributed local execution or real beta-gal data.

## Focused-suite parser correction (2026-10-27)

The downloaded `nu_handoff.log` contains stages a, b, c, e and f despite the
reported `suite=nu_handoff` command. Source tracing found an implementation
error: the selector reads `params%suite`, but `bind_input_carg` had no binding
for that existing field. Its default empty value therefore selects the full
gate. This is not evidence of an outdated binary or a failed positive NU
handoff; that case was not selected. Stop the unintended long run.

Correction contract: register the existing suite field through the normal
parameter parser, retain typed-parameter consumption, print the selected mode
before simulation, and add a fast constructor-level regression for default,
`nu_handoff` and `all` values. Use the existing startup suite and isolated
fixture; add no standalone test or alternative parser. Rebuild and pass the
fast regression before rerunning the focused integration.

The binding correction is implemented in the normal character registry. The
startup suite now has three constructor assertions for suite selection, and
the integration entry prints its selected suite before simulation. Registry,
Fortran-description, FLEX DAG and whitespace checks passed. No Windows build
or executable test was run. Linux must rebuild and pass the fast startup
regression before another focused run. The downloaded partial log is an
unintended full-gate execution, not positive NU-handoff evidence.

## Upstream integration without a feature commit (2026-10-27)

Saved and verified a complete recovery snapshot at `C:\Users\hossainm7\OneDrive - National Institutes of Health\Research\SIMPLE_worktrees\feature-refine3d-local-recovery-20261027-105615`. It includes the complete 2538-file working tree (untracked and ignored files included), original index, binary/full-index patches, original review export, per-file SHA-256 manifest and a verified Git bundle of the original HEAD history. No commit or stash was created.

Fetched canonical `origin/master` from `https://github.com/hael/SIMPLE.git`. The feature branch had no unique commits and was 45 commits behind. Fast-forwarded from `0b723362ce6354ca07038e4d35e4c3a9acfb0593` to `7484bf052de448597494fe28fe7d322b80a85d20`, then reapplied the saved local changes. The main master checkout was not modified.

Reviewed upstream changes in the 13 overlapping local files. Two textual conflicts were import-list overlaps: the exec router keeps both the upstream wave commander and local commander; the PCG strategy keeps both upstream image-array helpers and local provenance. All other local changes reapplied cleanly. The upstream PCG observation/accumulation refactor, posterior/annealing changes and added tests remain intact. Feature code, UI, suite parser binding and provenance remain present.

All 26 local paths remain uncommitted and unstaged; no unresolved index entries or conflict markers remain. The three untracked files are byte-identical to the snapshot. Before adding this record, non-overlap tracked contents were identical apart from Git LF/CRLF normalization in five files; the original byte representations remain in the archive. A module-use graph check found no cycle through the provenance utility. Test registry, description checks, FLEX DAG, generated UI defaults and diff whitespace checks passed.

No Windows compilation, SIMPLE execution or Linux run was performed. Earlier Linux passes predate this upstream integration; the suite-parser correction and focused positive NU handoff still require Linux validation. Regenerate/configure the Linux build for newly fetched sources, then build and run fast checks before `test=cont_refine3D_1jxy suite=nu_handoff`. Keep scientific thresholds unchanged when evaluating the new upstream behavior.

The adjacent Markdown review export now excludes this plan diff and links here directly. The complete original plan remains preserved in the recovery snapshot.

## Startup ownership correction (2026-10-27)

User approved correction of the two commander review findings. Preserve the
public contract above: no CLI, artifact, compatibility, numerical or final
reconstruction change. Keep early preflight before parameter initialization.

Implementation approach: move saved-reference staging and trailing-chain
validation/copying to the existing refine3D strategy. Each routine takes the
caller's loaded project, performs no project-file read/write, and returns its
reuse result. The commander remains the project owner and writes changed
segments before the child refinement starts. Move the existing raw half-pair
FSC calculation unchanged to the half-map diagnostics domain module. Remove
the preflight's redundant second reads of loaded stack and particle segments.
No new module, CMake target, parameter or persistent file is required.

Acceptance: preserve same-file copy guards, strict compatibility/provenance
checks, source preservation, NU matching-band transfer and full-first-update
fallback. Existing local_startup, fast and NU-handoff assertions stay unchanged.
Add a fast analytic check of the extracted FSC operation and a direct staging
check on the caller's project. Linux compilation/runtime remain required; the
pre-correction NU-handoff pass is the baseline, not validation of this move.

Observed before this correction: local_startup passed 334 assertions, unit_ui
passed 927 assertions and the fast gate passed all 15 entries. The focused
NU-handoff passed in 817.4 seconds: four staged constituent digests, first
matching LP (3.22759 A), 1000 continuous attempts, fixed half labels and 12
unchanged source files. Empty NU evidence support and final cFAR=0.0000 remain
separate scientific observations; this correction does not address them.

Source-checked correction: the two strategy routines now operate on a caller-owned
project, return reuse status, and do not read/write the project file. The
commander loads the required segments and publishes the updated particle/out
segments before the child run. Preflight retains its initial full read and
uses that object for sigma inspection instead of rereading stack/particle data.
The saved-half FSC sequence (read, FFT, FSC, binary write) is in
`simple_halfmap_diagnostics`; its masking/filtering behavior is unchanged.

Added direct checks to the existing local_startup tester: equal impulses give
FSC=+1 and opposite impulses give FSC=-1 on every shell (1e-6 tolerance);
reference staging updates the supplied project without a project file and
preserves source maps; both gridding and PCG chains accept in-place compatible
files without changing bytes and reject changed producer poses. Existing
constructor, rejection, provenance and integration assertions are retained.

Static registry, description, FLEX dependency, generated UI-default and
whitespace checks passed. Source review verified that only project loading,
publication and command-line selection remain with the commander for the moved
steps. No Windows build or SIMPLE test ran. Linux incremental build and
local_startup/fast checks must validate this correction; rerun the focused
NU-handoff after those pass to check production integration. The 817.4-second
pass above belongs to the source before this ownership correction.

### Chain fixture grid correction

The first Linux run after ownership extraction compiled successfully. UI
passed 927 assertions. local_startup passed 360 of 362; its new gridding and
PCG positive chain cases failed, and the fast gate passed 14 of 15 entries
in 50.93 seconds. The producer fixture left params%box_crop=0 while its
manifest used BOX=16. Publication therefore wrote an invalid grid, which
the unchanged provenance reader correctly rejected before compatibility
checks. The negative pose cases consequently did not test their intended
rejection reason.

Corrected the fixture's box_crop/smpd_crop to BOX/SMPD and added a production
provenance-consumer assertion before staging. A failed precondition records
an assertion failure and skips the dependent cases for that backend. No
production reader, validator, gate script or threshold changed. Static checks
are required; Linux local_startup and the fast gate remain pending for this
fixture correction. The large fast-gate log was CTest output-on-failure for
the failed area suite, not execution of the long integration stages.

## Deferred NU observations and administrator-review handoff (2026-10-27)

The user accepted deferring investigation of the empty NU evidence envelope
and final cFAR=0.0000. The source/log trace found no change to the NU filtering
algorithm, evidence-envelope algorithm, shared final-reconstruction policy or
cFAR scoring in this feature. Local changes concern validation, staging,
reference consumption and startup ownership. Do not treat these observations
as a reason to change NU thresholds or expand this refactoring's acceptance
criteria. Their scientific cause remains unresolved; defer them for separate
owner review unless evidence links a regression to this feature.

Retained evidence: `C:/Users/hossainm7/Downloads/nu_handoff.log`, focused
817.4-second run before the ownership extraction. Producer/local/final cFAR
values are 0.8331/0.9173/0.0000; corresponding global FSC=0.143 resolutions
are approximately 3.228/3.120/3.120 A. Both NU-producing iterations reported
an empty evidence field. With automsk=yes that field is diagnostic; the
density envelope is the intended active reference mask. The shared final
ending rebuilds at native sampling, uses the newly committed sigmas and
disables NU filtering under this mode. A good global FSC does not determine
the directional area ratio. The log lacks the directional curves and evidence
statistics needed to explain the printed zero or empty evidence. No defect or
scientific-quality conclusion is established by these observations alone.

The latest supplied Linux terminal output reports local_startup PASS,
364/364 assertions in 4.30 seconds. Both negative chain cases now identify
changed producer poses, after the positive fixture-provenance checks. The
downloaded fast-gate log reports 15/15 PASS in 49.41 seconds, including UI
and Cartesian alignment. The downloaded local_startup.log still contains the
earlier 362-check failure; use the later pasted terminal output for this
364-check result, not that stale file. Compilation was reported successful by
the user. No agent-side SIMPLE execution or Windows build occurred.

The focused `test=cont_refine3D_1jxy suite=nu_handoff` run after ownership
extraction passed in 949.0 seconds. It retained the verified NU handoff,
matching band, continuous attempts, fixed half labels and source-file checks.
It predates the stale-mask correction below. Broader PCG/distributed and
real-data gates remain as already stated; source review is not full
workflow/scientific acceptance.

### Post-review stale destination mask correction

Source review found that a valid NU reference record may omit an envelope
while a same-name mask from an earlier run remains in the destination. The
first volume assembly can then accept that stale file by geometry alone.
Reference staging now removes a destination automask or NU envelope if the
verified source record omits it. Recorded source masks and in-place copies
retain their existing file-identity protection. A local_startup assertion
starts with a stale destination envelope, republishes the source without that
constituent, and checks that staging removes only the unrecorded envelope.

The description, registry, FLEX dependency and whitespace static checks pass
after this correction. The user then reported a successful Linux compilation
and supplied the 2026-10-27 logs: local_startup passed 369/369 assertions in
7.87 seconds, five more than the preceding 364-check run; UI passed 927/927;
the fast gate passed all 15 entries in 74.48 seconds. Its 60-second budget
warning is a timing warning, not a failed test. These results validate the
focused correction and fast regression gate. The prior 949.0-second NU-handoff
pass predates this correction; broader PCG/distributed and real-data validation
remain separate acceptance work.
