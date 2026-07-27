# Dominant flex-density basin as a robust Bayesian 3D reconstruction estimator

## Status and intent

This is a forward-looking implementation note. It does not describe a mode
that is currently available in `refine3D`.

The goal is to couple the registered-residual diffusion manifold from
`flex_analysis` to Bayesian 3D refinement as a robust reconstruction
estimator. The first implementation deliberately does **not** sample or
reconstruct the full conformational landscape. It identifies the
view-balanced density basin with the largest integrated particle mass,
reconstructs the density peak of that basin, and uses the resulting even/odd
volumes as the next refinement references.

The intended outer iteration is:

```text
Bayesian pose update
        |
        v
registered-residual flex embedding
        |
        v
largest view-balanced density basin
        |
        v
local-linear even/odd pre-image at basin peak
        |
        v
normal volume assembly, filtering, and next Bayesian pose update
```

This is best viewed as a generalized robust EM procedure:

- the Bayesian matcher estimates pose variables from the current references;
- the flex graph estimates the dominant conformational population conditional
  on those poses; and
- reconstruction estimates the volume at that population's density peak.

The flex result is initially a **replacement reconstruction/reference
estimator**, not an additional same-iteration prior term. This avoids counting
the same particle observations twice while preserving the existing Bayesian
orientation objective.

This proposal is separate from
`denoise_project_flex3d_generative.md`. That note generates a one-to-one
denoised particle stack for the `objfun_den` correlation term. This note
changes which 3D density is reconstructed between refinement iterations.
Both mechanisms may eventually be combined, but they must first be validated
independently because otherwise the same flex model influences both alignment
and reconstruction.

## 1. Scientific definition

Let:

- \(y_i\) be particle \(i\);
- \(R_i,t_i\) be its current registered orientation and shift;
- \(z_i\in\mathbb R^d\) be its retained raw diffusion coordinate;
- \(p_i\) be its projection-direction bin;
- \(G=(V,E,W)\) be the symmetric registered-residual graph; and
- \(W_{ij}\) be the graph's raw Gaussian affinity.

The estimator has four distinct stages:

1. estimate particle poses by the current Bayesian refinement;
2. embed registered residuals using the current flex graph;
3. partition the graph into density-ascent basins and select the largest
   view-balanced basin;
4. reconstruct the conformation at that basin's density maximum.

### 1.1 Why the largest density basin

Three superficially similar targets have different meanings:

- A **global medoid** minimizes aggregate distance, but can lie between
  conformational modes.
- A **density maximum** selects the highest local peak, but a narrow,
  low-population state can have a higher peak than the prevalent state.
- The **largest density basin** selects the mode whose attraction basin has
  the greatest integrated population.

The third quantity is the required estimator. Within the selected basin, its
density maximum is the reconstruction target. The basin medoid remains useful
as a diagnostic but does not define the output state.

### 1.2 View-balanced particle mass

Preferred orientation must not make an oversampled view appear to be a
prevalent conformation. Reuse `projection_occupancy_weights`:

```text
omega_i = N / (N_occupied * n_{p_i}),
```

where \(n_{p_i}\) is the population of the occupied projection bin containing
particle \(i\). The weights have unit mean, and every occupied projection bin
has equal total mass.

These weights are used only to estimate conformational prevalence. They are
not reconstruction weights and must not be applied to the tomographic
denominator.

### 1.3 Graph density

For the first implementation, define the particle density directly from the
raw symmetric graph affinities:

```text
rho_i = omega_i + sum_{j in N(i)} W_ij * omega_j.
```

The first term is the \(K(0)=1\) self contribution. If
`view_balance=no`, set every `omega_i=1`.

Use `diffmap_graph%w`, not `diffmap_graph%wnorm`. Markov normalization is
appropriate for the diffusion operator but destroys the degree information
needed for a density estimate. The graph must therefore remain alive through
basin selection.

The initial implementation should not add a second density bandwidth. The
same raw kernel and bandwidth used to build the graph define `rho`. Record the
graph bandwidth in the diagnostics; if it is not currently persisted in
`diffmap_graph`, add it as explicit graph metadata.

### 1.4 Deterministic density-ascent basins

Partition the symmetric graph with a discrete watershed:

1. For every particle \(i\), inspect \(i\) and all graph neighbors.
2. Set `parent(i)` to the neighbor with the largest density strictly above
   `rho_i` by a numerical tolerance.
3. If no such neighbor exists, set `parent(i)=i`; this particle is a local
   density maximum.
4. Follow parent links to a maximum and path-compress the result.
5. Particles terminating at the same maximum belong to the same basin.

Determinism is required:

- compare densities using a documented relative tolerance;
- when densities are equal within tolerance, prefer the lower original
  particle index;
- assign basin identifiers by increasing peak particle index; and
- never let OpenMP traversal order decide a parent or basin label.

For basin \(b\), report:

```text
M_raw(b)  = count(i in b)
M_view(b) = sum_{i in b} omega_i
f_raw(b)  = M_raw(b)  / N
f_view(b) = M_view(b) / sum_i omega_i
```

The selected basin is:

```text
b_star = argmax_b M_view(b).
```

Ties are broken by larger raw mass, then higher peak density, then lower peak
particle index.

Both raw and view-balanced fractions are scientific outputs. A large
disagreement is not merely a warning about implementation quality; it can
indicate conformation-dependent preferred orientation and must remain visible
to the user.

### 1.5 Basin support diagnostics

For each basin, also compute:

- number and fraction of occupied projection bins;
- maximum fraction of its raw particles in any one projection bin;
- effective view-balanced population
  \((\sum_i\omega_i)^2/\sum_i\omega_i^2\);
- connected-component identifier and mass;
- peak density and peak particle;
- medoid particle and peak-to-medoid latent distance; and
- within-basin latent radius quantiles.

The largest basin remains the fixed selection rule. Quality gates determine
whether it is safe to use, not whether another basin should be selected.

If the selected basin fails minimum population or angular-support gates, the
iteration falls back to standard reconstruction. It must not silently choose
the second-largest basin.

## 2. Reconstructing the dominant conformation

Let \(b_\star\) be the selected basin and \(z_\star\) the raw diffusion
coordinate of its density peak. Use the existing local-linear pre-image model
centered at \(z_\star\):

```text
u_i = z_i - z_star
g_i = [1, u_i(1), ..., u_i(d)]
```

For each Fourier voxel, solve:

```text
A = sum_i w_i g_i g_i^T CTF_i^2
b = sum_i w_i g_i CTF_i^* P_Ri^* y_i
beta = (A + ridge)^-1 b
V_star = beta_0.
```

This is the same local-linear tomography already implemented by
`reconstruct_flex_diffmap_local_linear_states`. Only the target selection and
weight construction change.

### 2.1 Dedicated dominant-basin weights

Do not call `build_flex_preimage_kernel_weights` with `nstates=1`. Its current
multi-state density proxy is:

```text
q_i = sum_s exp(-||z_i-z_s||^2 / epsilon_s),
```

followed by division by `q_i`. With one state, that division cancels the
target kernel and produces nearly uniform weights.

Add a dedicated routine:

```fortran
build_flex_dominant_basin_weights( &
    raw_coords, basin_labels, selected_basin, peak_row, &
    regression_weights, responsibilities, bandwidth, neff )
```

For the first implementation:

```text
d_i^2 = ||z_i - z_star||^2
k_i   = exp(-d_i^2 / epsilon_star),  i in b_star
k_i   = 0,                            i outside b_star

w_i = k_i / sum_j k_j
r_i = k_i / max_j k_j.
```

Here:

- `w_i` is the normalized weight for the local-linear solve;
- `r_i` is an absolute \([0,1]\) dominant-state responsibility for a later
  weighted Bayesian reconstruction path.

The local bandwidth starts at the median positive \(d_i^2\) within the
selected basin and is inflated by the existing per-target factor of 1.5 until
the kernel effective sample size reaches a floor or the iteration cap is
reached. The density bandwidth used to find the basin and the reconstruction
bandwidth used to fit \(V_\star\) are different quantities and must be logged
separately.

The hard basin mask prevents a nearby secondary mode from bleeding into the
dominant-state regression across a density saddle. A genuinely unimodal
continuum has one basin and is unaffected.

Do not multiply `w_i` or `r_i` by `omega_i`. Projection occupancy weighting
defines which state is prevalent; tomographic sampling density and CTF
weighting already belong to reconstruction.

### 2.2 Three weights that must remain separate

The implementation will contain three conceptually different weights:

| Weight | Meaning | Normalization | Used by |
| --- | --- | --- | --- |
| `omega_i` | inverse projection-bin occupancy | mean one | density and basin mass |
| `w_i(z_star)` | local regression kernel | sums to one | flex local-linear pre-image |
| `r_i` | dominant-state responsibility | maximum one | future Bayesian weighted reconstruction |

They must have distinct names, types, output columns, and API arguments.
Reusing a generic `weights` array across these stages is too error-prone.

The existing `state_weights` columns sum to one because they are
reconstruction coefficients. They cannot be passed unchanged into a Bayesian
M-step: absolute scale affects its effective particle count and
regularization.

## 3. Coupling to the refinement loop

### 3.1 First production contract: replacement reference estimator

Add an explicit mode:

```text
robust_rec=none|flex_dominant
```

with `none` as the compatibility default.

When `robust_rec=flex_dominant` and the iteration is scheduled for volume
reconstruction:

1. complete the normal probabilistic/Bayesian pose update;
2. freeze the updated `ptcl3D` field;
3. run flex registered-residual feature extraction and graph embedding;
4. select the largest view-balanced density basin;
5. reconstruct its unfiltered even and odd local-linear pre-images;
6. pass those half-map products through the normal volume assembly,
   FSC, automasking, and filtering policy; and
7. use the assembled dominant-state volumes as the next iteration's
   references.

The basin estimator belongs after pose estimation and before volume assembly.
It must not be inserted inside `simple_strategy3D_matcher`; the matcher should
remain responsible for scoring and pose updates.

The final iteration of a production run should be standard raw-particle
reconstruction by default. The dominant-state estimator is used to stabilize
the pose/reference trajectory, while the terminal raw reconstruction provides
the least model-coupled reported half maps. A future explicit option may
retain the flex-dominant final map, but it must be labeled as a model-selected
output.

### 3.2 Schedule controls

The public controls should be small:

```text
robust_rec=none|flex_dominant
robust_rec_start=<iteration>
robust_rec_period=<positive integer>
robust_rec_terminal_raw=yes|no
robust_rec_strict=yes|no
```

Recommended defaults:

- start only after the current alignment resolution and angular sampling are
  adequate for meaningful registered residuals;
- update the basin every few refinement iterations rather than every
  iteration;
- reuse the last accepted basin estimator between scheduled flex updates only
  as the current reference, never by reapplying stale particle weights to new
  poses;
- use a standard terminal reconstruction; and
- fall back to standard reconstruction on failed quality gates unless
  `robust_rec_strict=yes`.

The usual flex controls (`k_nn`, `nang_nbrs`, `neigs`, `icm`,
`bandwidth_mode`, `bandwidth_tune`, `view_balance`, `preimage_ndim`, `lp`,
`nspace`, and `mskdiam`) should be registered once and reused. The dominant
mode has no `npreimages` control because basin selection is fixed.

### 3.3 Relationship to `objfun_den`

`objfun_den` already solves the scale/noise-normalization mismatch by using
cross-correlation to a denoised representative. No normalization changes are
needed here.

For initial production and validation:

- either use `robust_rec=flex_dominant`, or use flex-generated `stk_den` with
  `objfun_den=yes`;
- do not enable both by default.

The combined mode must be preceded by an ablation showing that flex influence
in both the alignment objective and reconstruction does not lock the
refinement to an early manifold error.

## 4. Gold-standard half-map contract

### 4.1 Required production behavior

Gold-standard independence is preserved only if each half selects and
reconstructs its dominant basin without using the opposite half's particle
images or map:

```text
even particles -> even flex graph -> even largest basin -> even V_star
odd particles  -> odd flex graph  -> odd largest basin  -> odd V_star
```

The even target and weights are used only for the even reconstruction; the
odd target and weights are used only for the odd reconstruction. Neither
half's map, high-resolution FSC, coordinates, peak particle, nor basin labels
may initialize or select the other half's result.

A common all-particle embedding followed by separate even/odd reconstruction
is useful as a diagnostic and development milestone, but is not the strict
production contract because the common graph transfers image-derived
neighborhood information between halves.

### 4.2 Half-basin compatibility

Independent halves can legitimately select different basins. Before accepting
the dominant estimator, compare them using only low-resolution,
non-circular diagnostics:

- raw and view-balanced basin mass;
- projection-bin support;
- low-dimensional basin radius and spectrum;
- low-pass map correlation below the alignment resolution; and
- orientation-distribution agreement.

Do not choose or match basins by maximizing high-frequency even/odd FSC.

If half compatibility fails, use standard reconstruction for both halves in
that iteration and record `fallback=half_basin_mismatch`. Never mix a
flex-dominant half with a standard half.

### 4.3 Output and FSC

Accepted dominant reconstructions must produce the same unfiltered even/odd
artifacts expected by `commander_volassemble`. FSC calculation, filtering,
automasking, nonuniform filtering, project `out` updates, and standard volume
names remain owned by the established volume-assembly path.

Do not compute a special flex FSC and silently substitute it for the
refinement FSC. Diagnostic flex-map correlations should have distinct names
and must not update project resolution metadata.

## 5. Architectural refactoring

The implementation should reuse flex algorithms, not call the
`flex_analysis` commander from `refine3D`.

### 5.1 Reusable flex engine

`run_flex_analysis` and related orchestration currently live in
`simple_flex_analysis_strategy.f90`. Extract a reusable service that accepts a
prepared project/registration context and returns data rather than writing
policy-specific final products:

```fortran
type :: flex_embedding_result
    integer, allocatable :: pinds(:)
    integer, allocatable :: proj_ids(:)
    real,    allocatable :: raw_coords(:,:)
    real,    allocatable :: coords(:,:)
    real,    allocatable :: eigvals(:)
    real,    allocatable :: view_weights(:)
    type(diffmap_graph)   :: graph
end type flex_embedding_result
```

Suggested module boundary:

```text
src/main/flex/simple_flex_diffmap_engine.f90
```

It should own:

- registered-residual feature preparation;
- projection-gated graph construction;
- optional occupancy weights;
- diffusion embedding; and
- graph/embedding diagnostics.

`flex_analysis`, `denoise_project`, and `refine3D` should be thin policy
clients of this engine.

### 5.2 Basin service

Add:

```text
src/main/flex/simple_flex_density_basin.f90
```

with a result object such as:

```fortran
type :: flex_density_basin_result
    real,    allocatable :: density(:)
    integer, allocatable :: parent(:)
    integer, allocatable :: labels(:)
    integer, allocatable :: peaks(:)
    integer              :: selected_basin = 0
    integer              :: selected_peak_row = 0
    integer              :: selected_peak_particle = 0
    real,    allocatable :: raw_mass(:)
    real,    allocatable :: view_mass(:)
    real,    allocatable :: view_fraction(:)
end type flex_density_basin_result
```

Public operations:

```fortran
select_largest_flex_density_basin(graph, pinds, proj_ids, view_weights, result)
build_flex_dominant_basin_weights(raw_coords, result, regression_weights, &
                                  responsibilities, bandwidth, neff)
```

The basin module must not depend on commanders, builders, image I/O, or global
project state. That makes its density, tie-breaking, and watershed behavior
unit-testable with synthetic graphs.

### 5.3 Reconstruction API

Refactor the current local-linear routine so the scientific core accepts an
explicit target coordinate and one weight vector:

```fortran
reconstruct_flex_diffmap_local_linear_target( &
    params, build, pinds, raw_coords, target_coord, regression_weights, &
    ndim, half_id, output_contract )
```

The existing multi-state routine can become a wrapper that loops over targets.
For the first dominant-basin implementation, `target_coord` is exactly the
density-peak particle coordinate. A later version may fit an off-particle
continuous mode without changing the reconstruction core.

The reusable solve should expose unfiltered Fourier/half-map products rather
than always finalizing a standalone `flex_diffmap_stateNN.mrc`. Refinement
policy, not the numerical pre-image solver, decides the durable filenames and
whether volume assembly follows.

### 5.4 Refinement integration point

In shared-memory refinement, add a robust-reconstruction policy step in
`inmem_execute_iteration` after:

```fortran
call refine3D_exec(...)
```

and before:

```fortran
call xvolassemble%execute(...)
```

The current `refine3D_exec` can emit reconstruction partials as part of the
matcher path. The refactoring must make pose estimation and reconstruction
dispatch explicit enough that a scheduled flex-dominant iteration can:

1. retain the updated orientations;
2. suppress or discard ordinary partial reconstruction products before
   assembly; and
3. emit compatible dominant-state even/odd products instead.

Do not let both partial sets reach `volassemble`.

The policy layer should decide among:

```text
standard partials
flex-dominant partials
standard fallback partials
```

while `commander_volassemble` remains unchanged wherever possible.

### 5.5 Shared and distributed execution

Implement shared-memory execution first to validate the estimator.

For distributed production, do not gather all registered particle images on
the master. Reuse the current distributed feature/neighbor pattern for
embedding, then distribute reconstruction accumulation. Each worker can emit
additive local-linear sufficient statistics for each half:

```text
B_q(k)     = sum_i w_i g_i(q) C_i^* P_Ri^* y_i
A_qr(k)    = sum_i w_i g_i(q) g_i(r) |C_i|^2.
```

The master reduces these arrays, solves the small coupled system per Fourier
voxel, and materializes the two half maps. This mirrors the existing partial
reconstruction architecture without transferring particle images.

The artifact manifest must identify:

- iteration;
- half;
- particle range/partition;
- embedding fingerprint;
- selected-basin fingerprint;
- latent dimension;
- target coordinate hash;
- reconstruction bandwidth;
- box and sampling;
- CTF convention; and
- completeness.

Assembly must reject mixed or incomplete fingerprints.

## 6. Longer-term Bayesian responsibility path

After the replacement estimator is validated, the dominant basin can be
expressed more directly as a robust Bayesian M-step using `r_i`:

```text
N(k) = sum_i r_i C_i^* P_Ri^* y_i
D(k) = sum_i r_i |C_i|^2
V(k) = N(k) / (D(k) + regularization).
```

This interprets off-basin particles as contamination rather than forcing
every particle to support one consensus volume. It generalizes naturally to
multiple responsibilities \(r_{is}\) later.

This path requires reconstruction accumulators that apply `r_i` to both the
data numerator and sampling-density denominator. Scaling only inserted image
values is incorrect.

The responsibility path is not required for the first implementation. The
local-linear dominant pre-image is scientifically preferable as the initial
robust estimator because it evaluates the conformation at the basin peak
rather than averaging the entire basin.

## 7. Optional empirical-Bayes prior, not first implementation

If a literal volume prior is later desired, it must use a stale,
half-specific flex estimate from a previous frozen outer iteration:

```text
p(V | V_flex_prev)
    proportional to
    exp[-1/2 sum_k tau(k) |V(k) - V_flex_prev(k)|^2]

V_MAP(k)
    = [N(k) + tau(k) V_flex_prev(k)]
      / [D(k) + tau(k)].
```

Required safeguards:

- `V_flex_prev` comes from the previous completed outer iteration;
- even uses only the previous even flex map and odd only the previous odd map;
- `tau(k)` is zero beyond a trusted low-resolution limit;
- prior strength decays as refinement becomes data-supported;
- the prior is disabled for the terminal raw reconstruction; and
- provenance records the source iteration and half.

Never use the opposite half as a prior. Never construct a flex map from the
current observations and add it as an independent same-iteration likelihood
term; that double-counts the data.

## 8. Artifacts and diagnostics

Every scheduled flex-dominant iteration should write:

```text
flex_dominant_iterNNN_summary.txt
flex_dominant_iterNNN_particles.txt
flex_dominant_iterNNN_basins.txt
flex_dominant_iterNNN_even_unfil.mrc
flex_dominant_iterNNN_odd_unfil.mrc
```

The particle table should contain at least:

```text
particle half proj_bin z1 ... zd density basin selected
view_weight regression_weight responsibility parent peak_particle
```

The basin table should contain:

```text
basin peak_particle raw_mass raw_fraction view_mass view_fraction
occupied_views max_view_fraction neff_view peak_density medoid_particle
```

The summary should record:

- all flex and robust-reconstruction parameters;
- graph size, connectivity, bandwidth, candidate-count diagnostics, and
  retained eigenvalues;
- selected basin and every tie-break decision;
- raw versus view-balanced occupancy;
- local regression bandwidth and effective sample size;
- half-basin compatibility metrics;
- accepted/fallback status and reason;
- input project and volume fingerprints; and
- exact output filenames used by the next refinement iteration.

Intermediate flex state maps must not be registered as the standard
refinement result unless the acceptance gates pass.

## 9. Failure and fallback policy

With `robust_rec_strict=no`, perform a complete standard reconstruction for
both halves if any of the following occurs:

- registered-residual feature or graph construction fails;
- the graph has an unacceptable disconnected-component structure;
- too few nontrivial diffusion modes remain;
- the largest basin is below the configured population or angular-support
  floor;
- the reconstruction kernel effective sample size is too small;
- the local-linear normal matrix is ill-conditioned beyond the ridge policy;
- even/odd dominant basins are incompatible;
- a distributed manifest is incomplete or inconsistent; or
- a generated half map is nonfinite or dimensionally invalid.

Fallback must be decided before volume assembly and logged with a stable
machine-readable reason. The iteration must never assemble one standard half
and one flex half.

With `robust_rec_strict=yes`, the same conditions are hard errors.

## 10. Validation plan

### 10.1 Unit tests

Add deterministic tests for:

1. one unimodal weighted graph producing one basin;
2. two modes where the broader, larger-mass basin wins over the taller,
   smaller peak;
3. preferred-orientation oversampling where raw and view-balanced selection
   differ;
4. equal-mass/equal-density tie-breaking by particle index;
5. disconnected graphs and singleton components;
6. invariance to CSR edge ordering and OpenMP thread count;
7. dominant-basin weights being zero outside the basin;
8. regression weights summing to one;
9. responsibilities having maximum one without column normalization; and
10. bandwidth inflation reaching the effective-sample floor independently.

### 10.2 Synthetic cryo-EM tests

Use simulated mixtures with known poses and CTFs:

- a dominant continuous state plus a minority alternative;
- a narrow minority mode with a taller density peak;
- conformation-dependent preferred orientation;
- increasing pose error;
- increasing noise; and
- a unimodal rigid control.

Measure:

- selected-basin precision/recall against latent truth;
- angular coverage;
- correlation/FSC to the true dominant conformation;
- pose error in the following Bayesian iteration;
- convergence rate;
- half-map independence; and
- fallback frequency.

The critical scientific regression is that the largest integrated basin, not
the tallest local peak or global medoid, is selected.

### 10.3 Real-data ablations

Compare:

```text
A. standard Bayesian refinement
B. standard refinement + objfun_den flex representatives
C. flex-dominant robust reconstruction
D. flex-dominant reconstruction + objfun_den
```

Run D only after B and C are independently stable. Compare orientation
stability, map interpretability, basin reproducibility across halves, FSC
without using it for selection, and terminal raw-reconstruction quality.

### 10.4 Restart and reproducibility

A restarted run must reproduce:

- particle ordering;
- projection occupancy weights;
- graph and eigenvector sign convention;
- density values;
- basin labels;
- selected peak;
- reconstruction weights; and
- output fingerprints.

Changing `nthr` or partition count must not change the selected basin except
within explicitly documented floating-point tolerances.

## 11. Phased implementation

### Phase 1: reusable selection and diagnostics

- Extract the reusable flex embedding service.
- Persist raw graph bandwidth metadata.
- Implement and unit-test deterministic density basins.
- Write basin/particle diagnostics.
- Run on frozen registrations without changing refinement output.

### Phase 2: shared-memory dominant half maps

- Add dedicated single-target basin weights.
- Refactor local-linear reconstruction to accept an explicit target.
- Independently embed and select even/odd basins.
- Produce and compare unfiltered dominant half maps.
- Keep the result diagnostic-only.

### Phase 3: refinement replacement estimator

- Add `robust_rec=flex_dominant` and schedule controls.
- Insert robust reconstruction between pose update and volume assembly.
- Add acceptance gates and all-or-nothing standard fallback.
- Preserve the normal FSC/postprocessing/project-output contract.
- Default to a terminal raw reconstruction.

### Phase 4: distributed execution

- Distribute flex feature and graph work using the existing flex partition
  pattern.
- Emit and reduce local-linear sufficient statistics by half.
- Add manifests, fingerprint validation, restart, and cleanup behavior.

### Phase 5: optional weighted Bayesian M-step

- Add absolute dominant responsibilities to numerator and denominator
  accumulation.
- Compare responsibility reconstruction with the local-linear peak estimator.
- Retain the local-linear estimator as the default unless validation supports
  replacement.

### Phase 6: optional combined and multi-state modes

- Calibrate combined `objfun_den` and robust reconstruction.
- Add stale half-specific empirical-Bayes priors only if they improve
  controlled tests.
- Generalize basin responsibilities to several conformations without changing
  the single-dominant default.

## 12. Initial scope and non-goals

The first production version is:

- single structural state in the refinement project;
- fixed largest view-balanced density-basin selection;
- density-peak local-linear reconstruction;
- independent even/odd analysis;
- shared-memory first;
- periodic robust reference estimation;
- standard terminal reconstruction; and
- deterministic standard fallback.

It does not initially:

- sample multiple conformations;
- switch among medoid, peak, or basin selectors;
- use flex coordinates directly in the Bayesian pose likelihood;
- add a same-iteration volume prior;
- combine automatically with `objfun_den`;
- select basins using FSC; or
- replace `commander_volassemble` and its postprocessing policy.

The key separation of concerns is:

```text
Bayesian refinement owns poses.
Flex owns registered-residual geometry and dominant-basin selection.
Local-linear tomography estimates the density at the basin peak.
Volume assembly owns half-map restoration, FSC, filtering, and durable output.
```

That separation marries manifold embedding with Bayesian 3D refinement while
keeping the statistical roles, half-map safeguards, and software ownership
clear.
