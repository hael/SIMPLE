# refine3D_local workflow proposal

Status: proposed, not implemented. Source review: `master` at `af13bc8286d5810d7e97523dfd4e013aa0c10951`,
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

Review performed by source inspection only. No Fortran source was changed, and no compilation,
linking or test binaries were run. Compilation and runtime validation remain with the maintainer, as required by the
repository instructions; the checks above are implementation acceptance criteria, not observed results.
