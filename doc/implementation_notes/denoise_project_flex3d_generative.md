# Flex-generated 3D particle representatives for `denoise_project`

## Status and intent

This is a forward-looking implementation note. It does not describe behavior
that is currently available through `denoise_project`.

The goal is to add a 3D denoising mode to `denoise_project` that:

1. starts from a fixed `ptcl3D` registration and mean volume;
2. uses the `flex_analysis` registered-residual diffusion map to estimate a
   particle manifold;
3. fits a projection-aware generative 3D model over that manifold;
4. evaluates the forward projection of that model at every particle's latent
   coordinate and registered orientation;
5. writes those predictions as the one-to-one `stk_den` particle
   representatives expected by `objfun_den`; and
6. preserves the raw registered stack as the source of the Bayesian/Euclidean
   refinement term.

This gives `denoise_project` two scientific modes under one durable
dual-representation project contract:

- **2D:** the current class-registered, within-class diffusion-map pre-image
  path;
- **3D:** the proposed globally embedded, locally generative, tomographically
  consistent flex path.

The 3D mode is not a second refinement engine. It produces a denoised
particle-representation project. `refine3D` and `abinitio3D` continue to own
orientation search, Bayesian scoring, pose updates, and reconstruction.

## 1. Scientific model

Let:

- \(y_i\) be the observed particle;
- \(R_i\) be its fixed input `ptcl3D` orientation;
- \(t_i\) be its fixed input in-plane shift;
- \(z_i\) be its retained diffusion coordinate;
- \(\widehat V(z)\) be the flex generative 3D volume model; and
- \(P_R\) be the CTF-free forward projection operator.

The generated denoised representative is

```text
d_i = P_{R_i} V_hat(z_i).
```

The particle is first moved into the existing flex registered frame, so the
durable output project has zero stored in-plane angle and shift and retains
the registered projection direction. The generated `d_i` is written directly
in that same frame.

The existing hybrid objective then uses:

```text
score_i(R,t)
    = (1 - objfun_den_w) * score_raw_i(R,t)
    + objfun_den_w       * clamp(CC(reference(R,t), d_i), 0, 1).
```

The raw term remains the observation-driven Bayesian/Euclidean term. The flex
representative supplies a lower-noise, conformation-aware orientation cue.
The denoised term remains a cross-correlation term; this proposal does not
change its normalization or introduce a second noise model.

### 1.1 Global modal form

For a global modal model,

```text
V_hat(z_i) = V_mean + sum_q z_iq B_q,
```

and linearity of projection gives:

```text
d_i = P_Ri V_mean + sum_q z_iq P_Ri B_q.
```

This is already the algebra used by the legacy modal pre-image path when it
synthesizes a small number of target volumes. It is useful as an
implementation diagnostic and baseline, but it should not silently replace
the current production local-linear estimator.

### 1.2 Local-linear chart form

The intended production generator follows the current
`preimage_mode=linear` model. For chart \(c\), centered at medoid coordinate
\(z_c\),

```text
V_hat_c(z_i)
    = B_c,0 + sum_q=1:d (z_iq - z_cq) B_c,q.
```

Here:

- `B_c,0` is the local-linear intercept volume;
- `B_c,q` is the fitted tangent/slope volume for diffusion coordinate `q`;
- `d = min(preimage_ndim, nmodes)`; and
- particle `i` is evaluated only within an assigned or smoothly gated local
  chart.

The present local-linear reconstruction already solves for `d+1` components,
but only finalizes and writes component 1. Components `2:d+1` are destroyed.
The central model refactoring is to retain or immediately forward-project all
solved components before cleanup.

### 1.3 No particle-specific 3D volumes

The implementation must not materialize one 3D volume per particle. Projection
and latent synthesis commute:

```text
P_Ri V_hat_c(z_i)
    = P_Ri B_c,0
    + sum_q (z_iq - z_cq) P_Ri B_c,q.
```

For registered particles indexed on the `nspace` grid, the expensive work is
therefore bounded by:

```text
ncharts * (d + 1) * nspace
```

component projections. Particle representatives are inexpensive 2D linear
combinations of cached component projections.

## 2. Public `denoise_project` contract

### 2.1 Mode selection

Add an explicit public selector rather than overloading `graph`:

```text
denoise_mode=2d|3d
```

with `2d` as the compatibility default.

`graph=ori` is currently associated with the experimental SO(3)-mixture graph
path and must not become an implicit synonym for the flex generator. Graph
construction and output dimensionality are separate choices.

Recommended interpretation:

| Input | `denoise_mode=2d` | `denoise_mode=3d` |
| --- | --- | --- |
| registration source | `ptcl2D` classes | fixed `ptcl3D` |
| mean volume | not required | `vol1` required |
| `nspace` | optional/current use | required |
| embedding | per 2D class | global gated flex graph |
| pre-image | 2D Nyström residual | local generative 3D model |
| output | registered `stk` + `stk_den` | registered `stk` + generated `stk_den` |

### 2.2 Flex controls exposed by 3D mode

The 3D mode should reuse the established flex controls and defaults:

- `neigs`;
- `icm`;
- `k_nn`;
- `nang_nbrs`;
- `bandwidth_mode`;
- `bandwidth_tune`;
- `view_balance`;
- `npreimages`, interpreted as the local-chart count;
- `preimage_mode`, initially restricted to `linear` for the production
  generator;
- `preimage_ndim`;
- `lp`, the registered-residual graph-feature low-pass limit;
- `nspace`;
- `mskdiam`;
- `nparts`; and
- `nthr`.

The UI definitions should be shared with `flex_analysis` through a small
registration helper so descriptions, defaults, and validation cannot drift.
The two programs may override whether an input is required, but should not
maintain independent descriptions of the same scientific control.

### 2.3 Required 3D-mode validation

Before registration or graph work:

1. `vol1` must exist and be dimensionally compatible.
2. `ptcl2D` and `ptcl3D` must contain matching particle records.
3. At least three active `ptcl3D` particles must be selected.
4. Exactly one input state is supported initially.
5. Every active `ptcl3D%proj` must lie in `1:nspace`.
6. Selected stacks must have a uniform supported CTF convention.
7. `preimage_mode` must be exactly `linear` for the first production 3D path.
8. `npreimages >= 2`, `preimage_ndim >= 1`, and `neigs >= 1`.
9. The requested local dimension must be clamped to the retained embedding
   dimension and logged.

The input project must not be modified.

## 3. Durable frame, CTF, and project contracts

### 3.1 Registered frame

Reuse the flex registered-image contract:

1. phase-flip `ctf=yes` input once;
2. do not flip `ctf=flip` input again;
3. apply the inverse stored `ptcl3D` in-plane shift and angle exactly once;
4. rotate astigmatism orientation into the registered frame;
5. set registered `ptcl3D` in-plane angle and shift to zero;
6. retain the projection direction and particle identity mapping; and
7. never transform the generated image a second time before it is written.

The raw and generated stacks must have identical:

- particle counts;
- row ordering;
- box sizes;
- sampling;
- stack partitioning; and
- active/inactive particle semantics.

### 3.2 Generated-image transfer convention

The flex reconstruction is CTF-aware and estimates a deconvolved 3D density.
The generated representative supplied to `objfun_den` should be:

```text
d_i = P_Ri V_hat(z_i)
```

without applying particle-specific `abs(CTF)` amplitude modulation. This
matches the current denoised cross-correlation path, which correlates the
denoised particle representation with an unmodulated reference projection.
The raw Euclidean term continues to use its established CTF handling.

The existing cross-correlation path is also the authoritative solution to
scale/noise-normalization differences between empirical and generated
representatives. This project does not introduce a replacement normalization
scheme.

### 3.3 Dual-stack project

The final 3D-mode project must use the same durable stack keys as the current
2D mode:

```text
stk     = registered phase-flipped raw particles
stk_den = flex-generated CTF-free particle representatives
```

The output stack segment should continue to record `ctf=flip` for the
registered project. The generated source is selected explicitly through
`stk_den`; the project must not create a second orientation field or overload
the CTF flag to describe the generated images.

The final artifact must be directly consumable by:

```text
refine3D objfun=euclid objfun_den=yes ptcl_src=raw
```

without a new `refine3D` program mode.

### 3.4 Mapping back to the native project

The registered-to-native identity mapping must be durable. After refinement
of the denoised project, `map_params_from_den` should remain the supported
route for composing the refined registered-frame pose with the original
particle frame.

Add a 3D regression test for this composition rather than introducing a
flex-specific mapping command.

## 4. Ownership and refactoring boundaries

### 4.1 Public workflow ownership

`commander_denoise_project` remains a thin command entrypoint. It selects a
mode-specific strategy and owns no flex mathematics.

Recommended high-level structure:

```text
commander_denoise_project
  -> denoise_project strategy factory
     -> 2D class denoise strategy
     -> 3D flex-generative denoise strategy
```

The existing lifecycle remains:

```text
initialize -> execute -> finalize_run -> cleanup
```

Do not implement 3D mode by constructing a `flex_analysis` command line and
invoking the `flex_analysis` commander. Program-to-program calls would couple
artifact names, working-directory behavior, queue orchestration, and cleanup.

### 4.2 Reusable flex services

Extract the scientific stages currently coordinated privately by
`simple_flex_analysis_strategy` into reusable `src/main/flex` services:

1. selected-particle validation;
2. registered feature/project preparation;
3. gated graph construction;
4. embedding and ICM rank selection;
5. medoid/chart selection;
6. local-linear chart fitting;
7. chart-component forward projection; and
8. particle representative synthesis.

`flex_analysis` should call the same services and continue to write diagnostic
state volumes. `denoise_project` should call them and request particle
representatives. Neither workflow should parse the other workflow's text
outputs.

### 4.3 Proposed result types

Introduce explicit in-memory results rather than parallel allocatable
argument lists. Names are illustrative:

```text
type flex_embedding_result
    integer, allocatable :: pinds(:)
    integer, allocatable :: proj_ids(:)
    integer, allocatable :: medoids(:)
    integer, allocatable :: labels(:)
    real,    allocatable :: coords(:,:)
    real,    allocatable :: raw_coords(:,:)
    real,    allocatable :: spectral_coords(:,:)
    real,    allocatable :: state_kernel_weights(:,:)
    real,    allocatable :: state_bandwidths(:)
    real,    allocatable :: state_neff(:)
end type
```

and:

```text
type flex_local_chart_model
    real, allocatable :: center(:)
    real              :: bandwidth
    real              :: neff
    ! Intercept and tangent component ownership, or a projected cache handle.
end type

type flex_projected_model_cache
    integer :: nspace
    integer :: ncharts
    integer :: ndim
    ! Deterministic mapping:
    ! (chart, component, projection) -> cached 2D image.
end type
```

The existing `flex_embedding_result` used by other workflows should be
extended only if its current compact contract remains appropriate. Otherwise
introduce a richer internal result and keep the compact public result as a
derived view.

Every result type must own `kill`/cleanup behavior. Avoid returning live
builder pointers or relying on module-global allocatables.

### 4.4 Artifact naming

Reusable flex services must accept an artifact prefix or artifact-writer
policy. Hard-coded `flex_*` names are unsuitable when called from
`denoise_project`.

Suggested 3D-mode diagnostic prefix:

```text
denoise3d_*
```

The durable particle stacks may retain the established names:

```text
denoise_raw_particles[_partNNN].mrcs
denoise_particles[_partNNN].mrcs
```

Diagnostic coordinates, spectrum, graph summary, chart table, and projection
cache should use the `denoise3d_` prefix to avoid collisions with a separate
`flex_analysis` run in the same directory.

## 5. Local chart fitting refactor

### 5.1 Preserve the current estimator

For chart `c`, current local-linear accumulation forms:

```text
u_i       = z_i - z_c
g_i       = [1, u_i]
A         = sum_i w_ic g_i g_i^T CTF_i^2
b         = sum_i w_ic g_i D_i
components = solve(A, b).
```

The existing observability floor, packed normal-matrix representation,
relative ridge, CTF-aware plane preparation, and gridding correction remain
authoritative. The generator refactoring must not change the fitted intercept
volume when given identical inputs.

### 5.2 Finalize every solved component

After `solve_coupled_basis_exp`, component handling becomes:

1. finalize component 1 as the intercept;
2. finalize components `2:d+1` as tangent volumes using the same Fourier-grid,
   scaling, and gridding-correction conventions;
3. apply any linear output filter consistently to intercept and tangents, or
   apply it after 2D synthesis;
4. forward-project the components on the required orientation grid; and
5. release the 3D component volumes when their projected cache is complete.

This keeps peak memory close to one chart's `d+1` reconstructors plus the 2D
cache, rather than retaining all chart volumes simultaneously.

### 5.3 Chart assignment and evaluation

The initial implementation should use the existing k-medoid hard label to
select chart `c` for particle `i`. This is deterministic and avoids confusing
the column-normalized state reconstruction weights with particle membership
probabilities.

For each particle:

```text
delta_z = raw_coords(i,1:d) - chart(c)%center(1:d)
d_i     = projection(c,0,proj_i)
        + sum_q delta_z(q) * projection(c,q,proj_i).
```

Add diagnostics for:

- `norm(delta_z)` relative to chart bandwidth;
- fraction of particles outside the chart's well-supported radius;
- generated-image finite values and energy;
- intercept/tangent contribution ratio; and
- chart-boundary discontinuity.

A later smooth-chart version may blend predictions from multiple charts using
new row-normalized gating coefficients. It must not reuse
`state_kernel_weights` directly because those columns are normalized for
state reconstruction, not across states for particle interpolation.

### 5.4 Extrapolation control

The generated model should not extrapolate arbitrarily far along a poorly
supported tangent. The first implementation should provide a deterministic
guard, for example:

```text
delta_z_used = delta_z * min(1, radius_c / max(norm(delta_z), tiny)).
```

The default radius should be derived from the chart's weighted coordinate
support and logged. It must not be inferred from diffusion eigenvalue
magnitude alone.

The unclipped coordinate and applied scale should be written to the particle
map for validation.

## 6. Forward projection and particle synthesis

### 6.1 Projection orientation

The representative must be generated at the same registered `ptcl3D`
orientation that defined its flex observation. Prefer the exact stored
registered orientation.

If the production cache uses integer `proj` indices, add an explicit
equivalence test between:

- exact stored-orientation projection; and
- cached `nspace`-grid projection at `ptcl3D%proj`.

The grid cache may be used only when that difference is below a declared
tolerance at the active low-pass support. Otherwise cache by unique exact
orientation or project particle batches directly.

### 6.2 Projection cache

The cache index must be deterministic:

```text
index(chart, component, proj)
    = ((chart - 1) * (d + 1) + component) * nspace + proj
```

where `component=0:d`.

The cache manifest records:

- format version;
- box and sampling;
- `nspace`;
- number of charts;
- local dimension;
- coordinate convention;
- whether projections are CTF-free;
- applied low-pass/filter identity; and
- source model/project identity.

Use established image-stack I/O. Do not create an unversioned raw binary
projection format.

### 6.3 Generated stack writing

For each active particle row:

1. resolve raw particle identity and registered row;
2. resolve chart and latent displacement;
3. resolve projection image(s);
4. synthesize the real-space representative;
5. apply the same spatial support/mask convention expected by
   `objfun_den`;
6. verify finite, nonzero supported energy; and
7. write exactly one `stk_den` row.

Inactive particle rows follow the current dense-output convention and receive
blank raw and generated images while remaining `state=0`.

The finalizer verifies exact row coverage, no duplicates, and identical
raw/generated stack dimensions before writing the project.

## 7. Shared-memory and distributed execution

### 7.1 Shared-memory reference path

Implement and validate shared-memory 3D mode first:

```text
validate
  -> register particles and prepare residual features
  -> build graph and embedding
  -> select charts
  -> fit one chart at a time
  -> forward-project solved chart components
  -> synthesize all particle representatives
  -> write dual-stack project
```

This path defines the numerical reference for distributed execution.

### 7.2 Distributed stages

The current 2D distributed scheduler partitions 2D classes. That assignment
unit is invalid for global 3D flex mode. Use explicit 3D stages:

1. **registration/features:** particle-row assignments, reusing the current
   flex distributed registration contract;
2. **graph/embedding/chart fit:** initially master-owned, matching current
   flex production behavior;
3. **representative generation:** particle-row assignments reading a
   master-written projected-component cache; and
4. **finalization:** master validation and dual-project assembly.

Worker assignments must use normal text assignment files and the existing
`qsys_env` lifecycle. Stage identity must be explicit; do not infer it from
the presence of arbitrary temporary files.

The shared and distributed paths must agree on:

- selected particles and ordering;
- registered images;
- graph and retained embedding up to eigenvector sign conventions;
- medoids and hard chart labels;
- chart component Fourier grids;
- generated particle images; and
- final project metadata.

## 8. Cross-fitting and model-bias control

Because `d_i` is generated using an input orientation, the denoised
cross-correlation behaves partly as a learned pose prior centered on the
registration used to fit the model. This is intended, but must be controlled.

### 8.1 Half-crossed reconstruction

The first cross-fitted production option should:

1. keep a common registered embedding and chart definition;
2. fit chart component models separately from even and odd particles;
3. generate even-particle representatives from the odd model; and
4. generate odd-particle representatives from the even model.

This removes direct image contribution to the generated representative while
avoiding an immediate requirement for out-of-sample graph embedding.

Document this accurately as **cross-reconstructed**, not strictly
cross-embedded.

### 8.2 Strict cross-embedding

Strict held-out generation would additionally require:

- fitting each diffusion graph/embedding on a training half; and
- assigning held-out coordinates through a validated Nyström
  out-of-sample extension.

This is a later statistical extension and should not be implied by the
half-crossed reconstruction option.

### 8.3 Refinement schedule

Recommended initial workflow:

1. obtain a coarse raw-particle registration and mean map;
2. run `denoise_project denoise_mode=3d`;
3. run hybrid refinement with a moderate `objfun_den_w`;
4. optionally regenerate representatives from a later frozen registration;
5. reduce or disable `objfun_den` in late refinement; and
6. verify the final solution with raw-only refinement.

The denoised project is a frozen model prediction. It must not be regenerated
inside a particle search iteration from poses being updated by that same
search.

## 9. Diagnostics and durable metadata

Write sufficient metadata to reproduce every generated row:

```text
denoise3d_particle_map.txt
```

with at least:

- generated row;
- raw particle identity;
- registered particle row;
- projection index;
- half-set;
- chart;
- chart-center particle;
- raw diffusion coordinates;
- latent displacement used;
- extrapolation scale;
- chart bandwidth and effective sample size; and
- generated-image energy.

Retain flex graph diagnostics under the `denoise3d_` prefix:

- coordinate table;
- spectrum and selected rank;
- graph summary;
- chart/medoid table;
- chart support metrics; and
- projected-model cache manifest.

Diagnostic representative volumes may still be written at chart centers, but
they are not the particle-generation intermediate and should be optional.

## 10. Validation plan

### 10.1 Algebraic forward-model tests

1. **Global linearity**

   Compare projection of a synthesized global modal volume with the linear
   combination of projected mean/basis components.

2. **Local linearity**

   Compare projection of
   `B_c,0 + sum_q delta_z(q) B_c,q` with the corresponding combination of
   projected chart components.

3. **Chart-center identity**

   At `delta_z=0`, the generated projection must equal the projection of the
   current local-linear pre-image state volume.

4. **Zero-tangent identity**

   With tangent components forced to zero, every member of a chart must
   receive its intercept projection.

5. **Projection-grid equivalence**

   Compare cached `proj`-grid generation with exact stored-orientation
   projection.

### 10.2 Project and frame tests

1. Raw and generated stacks have identical row counts and dimensions.
2. Registered images have zero stored in-plane transforms.
3. Applying registration once agrees with the durable raw stack.
4. Generated representatives are in the same registered frame.
5. `map_params_from_den` composes a known synthetic 3D pose update back into
   the native project correctly.
6. A dual project is accepted by `refine3D objfun_den=yes` without metadata
   repair.

### 10.3 Reconstruction-regression tests

The local-linear refactoring must reproduce current chart-center volumes:

- before output filtering, compare Fourier grids;
- after filtering, compare real-space maps;
- verify intercept equality within numerical tolerance; and
- verify shared/distributed agreement.

This test is required before interpreting tangent volumes scientifically.

### 10.4 Statistical tests

On simulated heterogeneous data with known poses and latent states, compare:

- raw-only refinement;
- current 2D `denoise_project`;
- 3D intercept-only generation;
- 3D local-linear generation;
- same-half generation; and
- half-crossed generation.

Measure:

- orientation and shift error;
- latent-coordinate recovery;
- representative-to-clean-image correlation;
- reconstruction FSC against ground truth;
- half-map FSC;
- chart-boundary error;
- sensitivity to initial pose error; and
- convergence after the denoised term is disabled.

On real data, require:

- reproducible half-set improvement;
- no loss after raw-only cleanup;
- stable results under moderate changes to chart count and local dimension;
- no dependence on a single diffusion mode sign convention; and
- no high-resolution gain unsupported by raw half-map FSC.

### 10.5 Failure tests

Hard-fail on:

- missing mean volume;
- inconsistent particle fields;
- missing or invalid projection indices;
- unsupported CTF mixtures;
- empty graph or embedding;
- singular chart solve without an established finite fallback;
- missing projected cache rows;
- duplicate or incomplete particle output;
- nonfinite generated images; and
- project stack-count mismatch.

## 11. Implementation sequence

### Phase A: reusable model boundary

1. Introduce or extend an owned flex embedding result type.
2. Extract registration, graph, embedding, and chart-selection services from
   private strategy orchestration.
3. Add artifact-prefix control.
4. Keep `flex_analysis` outputs numerically unchanged.

Acceptance criterion: current `flex_analysis` tests and chart-center outputs
remain unchanged.

### Phase B: local generative chart components

1. Finalize local-linear tangent components after the coupled solve.
2. Add algebraic projection tests.
3. Add chart support and extrapolation diagnostics.
4. Preserve intercept identity with the current estimator.

Acceptance criterion: component synthesis and direct volume projection agree
within numerical tolerance.

### Phase C: shared-memory `denoise_mode=3d`

1. Add UI/parameter validation.
2. Reuse flex registration and fit the embedding/model.
3. Generate a particle-indexed `stk_den`.
4. Reuse the current dual-stack project writer after generalizing names and
   frame metadata.
5. Run an end-to-end `objfun_den` smoke test.

Acceptance criterion: the generated project runs through one shared-memory
`refine3D` iteration without conversion.

### Phase D: cross-reconstruction

1. Fit even/odd chart component models.
2. generate representatives from the opposite half;
3. add leakage and half-set diagnostics; and
4. compare against same-half generation.

Acceptance criterion: every active particle is generated without entering
its own reconstruction model.

### Phase E: distributed execution

1. Reuse distributed flex registration.
2. serialize the versioned projected-component cache.
3. distribute particle-row synthesis.
4. validate exact coverage and reduce/assemble output stacks.

Acceptance criterion: distributed and shared outputs agree within the declared
image and model tolerances.

### Phase F: workflow integration

1. Exercise the 3D mode from `abinitio3D` with the existing
   `objfun_den` controls.
2. define when the model is frozen and optionally regenerated;
3. define the denoised-weight annealing schedule; and
4. require a raw-only terminal refinement.

Acceptance criterion: demonstrated robustness gain without a raw-data FSC or
pose-accuracy regression.

## 12. Non-goals for the first implementation

The first 3D mode does not:

- update flex coordinates during a `refine3D` particle-search iteration;
- replace the raw Bayesian objective;
- treat the denoised cross-correlation as an independent calibrated
  likelihood;
- reconstruct one volume per particle;
- infer uncertainty from diffusion eigenvalues alone;
- provide strict held-out diffusion embeddings;
- implement frequency-dependent latent bandwidths;
- introduce a new native-frame image generator; or
- change the `objfun_den` normalization/scoring path.

## 13. Principal code locations

Expected primary edits:

- `src/main/ui/simple/simple_ui_denoise.f90`
  - add `denoise_mode=2d|3d`;
  - share flex-control registration.
- `src/main/params/simple_parameters*.f90`
  - add and validate the mode selector if no suitable existing parameter is
    reused.
- `src/main/strategies/parallelization/simple_denoise_project_strategy.f90`
  - split 2D and 3D scientific execution;
  - preserve the dual-project finalizer.
- `src/main/flex/simple_flex_analysis_strategy.f90`
  - move reusable scientific work out of private workflow-only routines.
- `src/main/flex/simple_flex_diffmap_rec3D.f90`
  - expose local chart fitting;
  - retain/finalize tangent components.
- `src/main/flex/simple_flex_projected_latent_model.f90`
  - own the fitted component and forward-projection kernels.
- a new focused flex module, for example
  `src/main/flex/simple_flex_generative_projector.f90`
  - cache component projections;
  - synthesize particle representatives;
  - write generation diagnostics.
- `src/main/commanders/simple/simple_commanders_denoise.f90`
  - keep the commander thin;
  - extend mapping tests rather than adding mapping logic here.
- `production/tests/`
  - add algebraic forward-model, project-contract, cross-reconstruction, and
    end-to-end hybrid-refinement tests.

No `refine3D` scientific change is required for the first implementation
beyond tests confirming that its existing `stk_den` and cross-correlation
contract accepts the generated representatives.

