# `flex_analysis` implementation plan

This is the forward-looking companion to
`doc/policies/flex_analysis_policy.md`.  It deliberately does not redefine
current behavior.  The policy is authoritative for behavior already wired into
the executable.

## Current boundary

The code has a registered Cartesian feature stack, a sparse gated diffusion
graph, embedding/medoid selection, and shared-memory local-linear pre-image
reconstruction.  With `nparts>1`, only feature registration is actually
scheduled on workers.  The master assembles features and performs graph,
embedding, and state reconstruction.

There are dormant worker handlers for graph rows and a legacy modal statistics
path.  Their existence must not be read as an implemented distributed
workflow.  Any work below must preserve the registered-image contract and use
normal `qsys_env` assignments; no custom worker-state format is permitted.

## 1. Distributed graph construction

Wire stage-2 graph workers into the master lifecycle only after the following
contract is implemented and tested:

1. Each worker receives its contiguous source-row assignment through the
   normal text assignment file.
2. Workers read the same complete registered residual-feature table and write
   only assigned neighbour rows.
3. The master verifies exact, non-overlapping row coverage before CSR assembly,
   then uses the shared graph normalization/eigensolver.
4. Shared and distributed graphs must agree in candidate counts, retained
   neighbour identities/distances, bandwidth, CSR structure, and embedding up
   to numerical eigenspace sign conventions.

The principal design question is memory: broadcasting/reading the full feature
table at every worker may erase the benefit.  Benchmark a read-only shared
table or an explicitly partitioned candidate exchange before selecting a
transport.  Do not re-register particles or reproject the mean in graph
workers.

## 2. Distributed local-linear pre-image reconstruction

The active estimator is per-state, per-voxel coupled local regression.  It
cannot use the existing legacy modal reduction unchanged: that protocol stores
a single shared density for independently fitted global residual modes,
whereas local linear needs the packed `(d+1) x (d+1)` normal-matrix moments for
each state.

The implementation should:

1. Assign registered-project row ranges and the required state kernels/target
   coordinates through ordinary assignment files.
2. Let each worker accumulate, for every assigned state, the RHS vector and
   packed normal-matrix volume using the existing coupled accumulation kernel.
3. Write the partial statistics using a versioned, validated transport owned by
   the projected-model/reconstructor layer, then sum corresponding statistics
   on the master and run the same coupled solve used in shared memory.
4. Compare shared and distributed state Fourier grids before output filtering;
   output maps alone are insufficient because filtering can mask a reduction
   error.

The reduction should reuse established reconstructor/statistics abstractions
where they fit.  It must not raw-dump expanded matrices with ad-hoc byte-offset
arithmetic.  Remove or isolate the unused legacy modal stage once the
local-linear reduction replaces it.

## 3. Local-linear validation still required

`simple_test_flex_projected_latent_model` currently exercises synthetic
local-constant reduction (T1), one-sided affine/curved response improvement
(T2), a rank-deficient finite-solve fallback (T5), and independent per-state
bandwidth selection (T6).  The remaining work is:

- T3: explicitly compare endpoint and interior behaviour;
- T4: apply local-linear `data_scales=w*g` and
  `density_scales=w*g*g^T` to the single-record/batched adjoint accumulation
  comparison, rather than testing the generic kernel alone;
- T7: specify and validate any frequency-dependent bandwidth before enabling
  one; none exists today;
- T8: add a synthetic statistical harness with half-set and permutation-null
  checks.

The current relative ridge is a numerical fallback, not an evidence-based
per-voxel uncertainty model.  A change to ridge strength or degeneracy policy
requires synthetic and real-data regression evidence.

## 4. Native-image reconstruction is a separate project

Current states reconstruct from the registered project.  A native-image path
requires an explicit, tested mapping from each registered row to its source
row and composition of the inverse removed `ptcl3D` in-plane transform with
the fixed projection orientation.  It must not use `ptcl2D` transforms.  Begin
only after the registered identity test and a native-vs-registered
reconstruction equivalence test agree at the selected support.

## 5. Filtering and uncertainty

The current code optionally transfers the input project's state-1 FSC filter
to every representative.  This is not a state-specific confidence estimate.
Before changing the filter policy, add state half-map reconstruction and define
which particle weights/assignments enter each half.  Posterior variance maps,
per-state FSC, and resolution-adaptive local bandwidths remain future work;
they must not be inferred from diffusion eigenvalues alone.

## 6. Performance and code cleanup

Benchmark the current shared and distributed-registration paths separately:

- registration/feature preparation time and stack I/O;
- graph memory, candidate counts, and distance-evaluation time;
- local-linear accumulation and coupled-solve time per state;
- peak memory as a function of `preimage_ndim`, box size, and state count.

Keep reuse boundaries intact: `transform_ptcls`/project-aware batched I/O for
registration, the shared diffusion graph engine for graph math, and the
reconstructor/coupled gridding path for state reconstruction.  Do not replace
these with flex-specific image readers, reprojection kernels, or stack formats.

After distributed local-linear reconstruction exists, delete or clearly
quarantine the unused global modal worker path so source structure no longer
implies a production capability that is absent.

Validate `preimage_mode` explicitly as `constant|linear|discrete` during flex
input validation. The discrete diagnostic must terminate after hard k-medoids
assignment and must not fall through to either pre-image estimator.
