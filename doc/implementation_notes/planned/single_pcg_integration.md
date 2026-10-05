# SINGLE on the PCG backend: the atomic density prior

Implementation note, 2026-09-15. Design for review; nothing here is
implemented. Written against master `04ea41302`. Companion pages:
`pcg_backend_overview.md` (what the solver is today),
`reconstruct3D_pcg_policy.md`, `pcg_decision_log.md`.

## 0. Constraints that shape the design

SINGLE data sets are small (hundreds to a few thousand particles), the
objective is `cc` with `sigma_est=global`, `ml_reg=no` is forced
(`exec_refine3D_nano`), and the even/odd FSC carries no usable
information. Consequently, nothing below uses a half-map FSC, a
half-map noise estimate, or the ML/SSNR prior. The even/odd *particle
partition* is used once, for image-level cross-validation (§4.3), which
is a different thing from a half-map comparison.

The only regularization that has worked in SINGLE is ICM
(`image%ICM3D_eo`, `simple_image_filt.f90`). What it is, precisely: the
map is quantized to 256 grey levels, and three sweeps of coordinate
descent minimize, per voxel, `(x − y)²/σ² + λ Σ_{j∈N6(i)} (x_i − x_j)²`
with `λ = 0.1` and `σ² ≈ 10` (the eo local noise variance normalized to
mean 5, std 2, so only its spatial shape enters). That is the
conditional of the global Gaussian-MRF energy
`E(x) = ‖x − y‖²/σ² + λ xᵀ L x`, `L` the 6-neighbour graph Laplacian on
the voxel grid. Two properties matter: (i) it is a quadratic prior, so it
is exactly a precision term in the normal equations; (ii) the
quantization makes `λ` dimensionless relative to the map's dynamic
range, which is why one value transfers across data sets. Both are kept.

The simulated and reconstructed maps are not on a common scale. No
quantity below is a distance between them. The atomic model enters only
through a projector (scale-invariant by construction) and through
correlations; the one radial parameter that is fitted from the data is a
B factor, not an amplitude.

Current loop (`exec_autorefine3D_nano`): `refine3D_nano`
(`maxits_between` = 10 iterations, `trail_rec=yes`, `ufrac_trec=0.5`,
references = ICM-filtered eo pair) → `detect_atoms`
(`nanoparticle%identify_atomic_pos` → `_ATMS.pdb`, `_SIM.mrc` from
`atoms%convolve`, `_BIN.mrc`, `_CC.mrc`) → next `refine3D_nano` with
`vol1 = _SIM.mrc`. The simulated map enters the refinement only as the
starting reference of the next stage; the trailing reconstruction then
blends it away at a fixed fraction. That is the "fractional update" this
note replaces.

PCG state relevant here (`simple_reconstructor_pcg.f90`): the unknown
`x` is a real-space `box³` volume; `apply_normal(p)` = `mask · (H p) +
λ₀ p` with `H` the kernelized `AᴴWA`; the preconditioner is the Fourier
diagonal `1/(ρ + λ₀)`; `set_lambda_relative` expresses coefficients
relative to `data_scale` (mean of `ρ` over the six lowest native
shells, `update_lambda_from_density`); trailing is available at the
accumulator level (`add_raw_accum_weighted`). Because the iterate is in
real space, a real-space prior costs no extra FFTs.

## 1. The prior

The reconstruction is modeled as `x = G c + e`. `G` is the `N_vox ×
N_atoms` matrix whose column `g_i` is the density kernel of atom `i` at
its model position; `c` are free per-atom amplitudes; `e` is what the
model does not explain. The prior penalizes `e` and, independently,
roughness:

    Q(x) = λ_A · xᵀ (I − P_G) x  +  λ_L · xᵀ L x,       P_G = G (GᵀG)⁻¹ Gᵀ

`P_G` is the orthogonal projector onto the span of the atom kernels;
`I − P_G` is symmetric and idempotent, so `(I−P_G)ᵀ(I−P_G) = I − P_G`
and the term is already in normal-equation form. The regularized solve is

    (H + λ₀ I + λ_A (I − P_G) + λ_L L) x = b

symmetric positive semi-definite, solved by the existing CG. Properties:

- Scale invariance. `P_G` is invariant to rescaling any column of `G`
  and to rescaling `x`. The model contributes positions, element
  identity and relative kernel shape; nothing about amplitude, per atom
  or overall.
- Atomicity. Along `span(G)` there is no penalty: the data set every
  atom's amplitude. Off it, density is shrunk toward zero with weight
  `λ_A`. Where no kernel reaches, `P_G x = 0` and the term is plain
  `λ_A ‖x‖²`: vacuum between and around atoms is enforced by the same
  operator. As `λ_A → ∞` the solution is exactly `G ĉ`, a sum of atoms
  with least-squares amplitudes; `λ_A = 0` is the current base solve.
  One parameter spans the ladder from "no model" to "atoms are the
  unknowns"; no separate atomic-basis solver is needed.
- ICM inside the solver. `λ_L L` is the ICM energy as a precision term,
  applied against the particle data through `H` rather than against the
  gridded map, and minimized exactly rather than by three quantized
  sweeps. `L` is a Fourier multiplier
  `L̂(k) = 2 Σ_a (1 − cos 2π k_a / box)`, so it goes into the diagonal
  preconditioner exactly.

Version 1 uses one `λ_A` for the whole box. A real-space weight `W(r)`
(per-atom trust) would enter as `(I−P_G) W (I−P_G)`; it is deliberately
left out until the hard atom-inclusion rule of the M-step (§5) proves
insufficient.

## 2. The kernels and the one fitted radial parameter

`g_i` is the real-space profile `atoms%convolve` already uses: the
five-Gaussian electron scattering factor of element `Z_i`, each Gaussian
broadened by a B factor, evaluated within `cutoff = 8·smpd` of the atom
centre. `convolve` uses `B_min = (4·lp)²` with `lp = 2·smpd`. The prior
uses `B_i = B_min + B_data`, with one `B_data` fitted from the data per
stage. `B_data` is the only radial freedom, and it is what makes the
model's Fourier extent data-driven: the kernel's high-frequency content
is tied to its low-frequency content by the fitted shape, which is the
atomicity extrapolation the current SIM reference performs implicitly,
now with the width chosen by the data instead of by `2·smpd`.

Fit of `B_data` (`pcg_atomic_prior%fit_bfac_data(x_base)`): with `m`
the model map simulated at `B_min`, compute per shell `s`

    c_s = Re⟨m̂_s, x̂_s⟩ / ‖m̂_s‖²        (regression of map on model)

which is unbiased under map noise (noise is uncorrelated with `m`), and
the shell correlation `r_s = Re⟨m̂_s, x̂_s⟩ / (‖m̂_s‖ ‖x̂_s‖)`. Fit
`log c_s = log K − B_data k_s²/4` by weighted least squares over shells
with `r_s ≥ 0.5`, weights = shell voxel counts; `K` is discarded (the
projector does not see it). Clamp `B_data ≥ 0`. Using `r_s ≥ 0.5`
without noise correction is conservative in the right direction: map
noise depresses `r_s`, so the fit uses only shells where the model is
plainly supported. Per-atom B from the pdb (`betas`) is a later
refinement; version 1 uses one `B_data`.

Kernel tables: one radial table per `(Z, B)` pair on the crop grid;
footprints are precomputed as packed voxel-index lists per atom, clipped
to the solve support (columns of `G` are masked, so `P_G` commutes with
`mask_mul` and lives on the solve domain). Cost of a gather or scatter
is `N_atoms × footprint` ≈ `5·10³ × 2·10³ = 10⁷` operations, negligible
next to the FFTs of a CG iteration.

Gram matrix `GᵀG`: sparse, entries only for atoms whose footprints
overlap (centre distance `< 2·cutoff`). Solve `(GᵀG) c = Gᵀx` by CG
(inner solve, `N_atoms` unknowns, sparse matvec, ~20 iterations,
tolerance 1e-6) or by a banded Cholesky after ordering atoms along a
space-filling curve; the Gram factorization is built once per model
stage. `GᵀG` is positive definite as long as no two atoms coincide;
`detect_atoms` guarantees a minimum separation, and the builder rejects
a model that violates it.

## 3. Operator, preconditioner, convergence

`apply_normal(p)` gains one term, computed on the masked input and
masked on output:

    hp = mask · H (mask · p) + λ₀ p + λ_A (p − P_G p) + λ_L L p

`P_G p`: gather `Gᵀp`, inner solve, scatter. `L p`: 6-point stencil.
Neither touches Fourier space. Both are symmetric on the solve domain,
so CG's assumptions hold; the existing indefinite-stop guard remains the
safety net.

Preconditioner: `1/(ρ + λ₀ + λ_A + λ_L L̂(k))` on the padded lattice,
with `L̂` evaluated at native frequency `k/padf`. This is the exact
inverse of the operator on the orthogonal complement of `span(G)`; the
preconditioned operator differs from the identity by a term of rank at
most `N_atoms`. CG converges in a number of iterations governed by the
eigenvalue clustering of that low-rank part, not by `N_vox`. Expected
budget: `maxits_pcg` 20–30 at SINGLE box sizes, seconds per solve.
Report `ITS`, `RESID`, `MRES` as for the base solve; add the fraction of
`b`'s energy in `span(G)` as a one-line diagnostic.

Data term for SINGLE: `sig2 = 1` (unweighted; `objfun=cc` means the
alignment does not whiten either), `W = CTF` only. `data_scale` is
derived as now. `λ_A` and `λ_L` are stored as relative coefficients and
multiplied by `data_scale` in `update_lambda_from_density`, like
`lambda_rel`.

## 4. Strengths

Two numbers, both dimensionless relative to `data_scale`, neither from
half maps.

### 4.1 `λ_L` (ICM parity)

Calibrated once by matching the reference ICM produces. On two or three
reference sets, run `refine3D_nano` one iteration on gridding with
`icm=yes` and keep the filtered reference `r_ICM`; solve PCG with
`λ_A = 0` and `λ_L` on a log grid; take the `λ_L` maximizing
`corr(x(λ_L), r_ICM)` inside the mask. Expect one value to transfer,
for the same reason `λ = 0.1` transfers (both are dimensionless). Once
`λ_A` is in, `λ_L` is expected to shrink toward zero: `I − P_G` already
suppresses non-atomic structure and vacuum noise. Whether it is still
needed is decided by §4.3, not assumed.

### 4.2 `λ_A` default estimator (single map, no half maps)

The reconstruction noise is stationary in real space (no mask is
applied before the estimate), so its variance is measured in the empty
part of the box. On the base map `x_b` (the solve with `λ_A = 0`):

    N   = var(x_b) over the solvent shell: inside the spherical mask,
          outside the union of atom footprints dilated by cutoff
    T   = mean(x_b²) over the support V (union of footprints)
    S   = max(T − N, ε·N)                       signal power per voxel
    E   = mean((P_G x_b)²) over V               power explained by the model
    E'  = max(E − N·N_atoms/|V|, 0)             minus the noise that lands in span(G)
    C²  = min(E'/S, 0.98)                       model-explained fraction of signal
    λ_A,rel = (N/S) / (1 − C²)

Reading: the unexplained component has prior variance `(1 − C²) S`; the
data precision per voxel is `1/N`; the Wiener weight of the prior
relative to the data is their ratio. This is the sigma-A weight of
crystallographic phase combination, computed from map moments instead of
amplitude statistics. It has no free constant, and every input is on the
reconstruction's own scale, so the simulated map's scale never enters.
It is an approximation (shell dependence is collapsed into
`data_scale`); §4.3 is the check.

### 4.3 Validation of `λ_A` by image-level cross-validation

The quantity SINGLE refinement maximizes is the `cc` between particle
images and reprojections. Held-out `cc` is therefore the right criterion
for a prior, and it needs no FSC: solve on the even particles'
accumulator with `λ_A` on the grid `{0, ¼, ½, 1, 2, 4} × λ_A,rel`
(`λ_L` fixed), score

    CV(λ_A) = Σ_{i ∈ odd} cc( y_i , CTF_i · P(R_i) x_even(λ_A) ) + (even ↔ odd)

using `forward_plane` / `fourier_dot` of the reconstructor over the
`hp`–`lp` band of the run, and take the argmax. The even/odd partition
supplies independent images; no half-map comparison is made. Cost: six
solves per half per stage, seconds each. In development the scan runs
every model-building stage and the curve is logged; the default
estimator is accepted when the argmax sits at `1 × λ_A,rel` within a
factor of two across the benchmark sets. In production the scan is off
and the estimator is used, or `lambda_atm` is set explicitly.

## 5. Where it sits in the loop

`refine3D_nano` with `rec_backend=pcg atm_prior=yes pdbfile=<_ATMS.pdb>`:

1. Particle pass into raw accumulators (as now). Trailing stays at the
   accumulator level (`add_raw_accum_weighted`, `ufrac_trec`): it sums
   normal equations across iterations, which is a data-term operation and
   orthogonal to the prior. `vol1 = _SIM.mrc` is still accepted as the
   first alignment reference of a stage; it never enters the solve.
2. Base solve, `λ_A = 0`, `λ_L` as configured → `x_b` per half and
   merged. Fit `B_data` (§2) on the merged base map, build `G`,
   estimate `λ_A` (§4.2).
3. Regularized solve → `x_r`. References for alignment are `x_r`
   (masked as now). `icm` is forced to `no` when `atm_prior=yes` and
   logged; the ICM filtering block in `simple_matcher_refvol_utils.f90`
   is bypassed.
4. Products per iteration: `recvol_state01.mrc` = `x_r` (shipped, and
   the input of `detect_atoms`, as the ICM/trailed map is today);
   `recvol_state01_base.mrc` = `x_b`; `recvol_state01_res.mrc` =
   `x_b − G ĉ` with `ĉ = (GᵀG)⁻¹Gᵀx_b` (the residual map, computed on the
   base map so it is not shaped by the prior). A `_atmprior.txt` sidecar
   with `B_data`, `N`, `S`, `C²`, `λ_A,rel`, `λ_L`, `N_atoms`, CG stats.
5. `detect_atoms` runs on `x_r` unchanged in version 1. The residual map
   is diagnostic only until §7 step 4: peaks in `x_b − G ĉ` above
   `3·sqrt(N)` are missing-atom candidates, atoms with
   `ĉ_i < 0.2 · median(ĉ)` are removal candidates. Lattice fitting, atom
   validation and per-atom B stay where they are.
6. `autorefine3D_nano` sets `pdbfile` on `cline_refine3D_nano` after
   every `detect_atoms` and passes `atm_prior`, `lambda_atm`,
   `lambda_lap`, `bfac_data` through. Stage 1 (no model yet) runs with
   `λ_L` only.

The loop is then an EM iteration: the quadratic E-step is the
regularized solve; the M-step is `detect_atoms`, where the non-Gaussian
part of the prior (point process, lattice, minimum separation) lives.
The solver never sees anything but a quadratic form.

## 6. Code map

New: `src/main/volume/simple_pcg_atomic_prior.f90`, type
`pcg_atomic_prior`: fields — atom positions (crop voxel coordinates),
`Z`, `B`, per-`(Z,B)` radial tables, packed footprints (indices, values),
CSR Gram matrix and its factorization or inner-CG workspace, `λ_A`,
`λ_L`, support mask pointer. Procedures — `new/kill`,
`build(atoms, box, smpd, bfac_data, cutoff, mask)`, `gather(x) → Gᵀx`,
`scatter(c) → Gc`, `solve_gram(v) → (GᵀG)⁻¹v`, `project(x) → P_G x`,
`apply(x) → λ_A (x − P_G x) + λ_L L x`, `residual(x, c)`,
`fit_bfac_data(x_base)`, `estimate_lambda_atm(x_base, mask)`,
`laplacian_multiplier(k)`, `stats`.

`simple_reconstructor_pcg.f90`: `set_atomic_prior(prior)` (pointer,
`l_atomic_prior`); `apply_normal` adds `prior%apply(pm)` before the
final `mask_mul`; `build_precond` adds `λ_A + λ_L L̂(k)` to the
denominator; `update_lambda_from_density` scales the two relative
coefficients; `get_atomic_prior_stats`. The ML prior and
`shrink_by_ml_prior` are untouched and unused on this path.

`simple_parameters.f90` / phases: `atm_prior` (yes/no, default no),
`lambda_atm` (relative; `0` = use the §4.2 estimator), `lambda_lap`
(relative; default from the §4.1 calibration), `bfac_data` (override of
the fit), `atm_cv` (yes/no, the §4.3 scan). Validation:
`atm_prior=yes` requires `rec_backend=pcg` and `pdbfile`; forces
`icm=no`, `ml_reg=no`.

Reconstruction strategy (`simple_commanders_rec_distr.f90` PCG path):
base solve → prior build → regularized solve per half; the CV scan as a
subroutine called between them when `atm_cv=yes`; the three products
and the sidecar.

`single_commanders_nano3D.f90`: `pdbfile` handoff, parameter pass-through,
copying of the new products into `iteration_NN/` and `final_results/`.

`atoms%convolve`: factor the per-element `(a, b)` tables and `epot` into
a reusable kernel-table builder so the prior and the simulator share one
definition of the atom profile.

## 7. Order of work and the tests that gate each step

1. PCG in `refine3D_nano` with `λ_L L` only, `atm_prior=no`. Gate: ICM
   parity (§4.1) — a `λ_L` value that reproduces the ICM reference on
   the calibration sets; alignment convergence and the atom statistics
   from `detect_atoms` no worse than the gridding+ICM baseline over a
   full `autorefine3D_nano` run.
2. `pcg_atomic_prior` with unit tests on a synthetic box: `P_G` is a
   projector (`P_G² = P_G`, symmetric, `P_G G = G`); `apply` is symmetric
   (`⟨u, Q v⟩ = ⟨Q u, v⟩` to 1e-6); invariance of `P_G` to column
   scaling; Gram solve accuracy; `L̂` matches the stencil by FFT.
3. Full prior in the loop, `λ_A` from the estimator. Synthetic gate via
   `simulate_nanoparticle`: particles from a known model at realistic
   SNR and particle count; the prior fed a damaged model (10% atoms
   removed, 10% spurious added, all positions jittered 0.2 Å, `B` off by
   a factor 2, simulated map scaled by 10×). Pass: the missing atoms are
   the top peaks of the residual map; the spurious atoms get `ĉ_i ≈ 0`
   and no density in `x_r`; the fitted `B_data` is independent of the
   injected scale; `detect_atoms` on `x_r` recovers the true model with
   fewer errors than the baseline loop. Real gate: benchmark sets,
   baseline vs prior, on atom count, validity statistics, RMSD to the
   consensus model, and the held-out `cc` of §4.3 as the scalar summary.
4. §4.3 scan across the benchmark sets to confirm the estimator's
   calibration; decide whether `λ_L` stays.
5. `detect_atoms` consumes the residual map for add/drop; per-atom B in
   the kernels; per-atom trust `W(r)` only if step 3 shows the hard
   inclusion rule failing on a class of atoms (surface, low-Z).

## 8. Deliberately absent

No half-map FSC or half-map noise map anywhere (the eo partition is
used only as two sets of images in §4.3). No ML/SSNR prior. No Fourier-
shell precision inside the solve: the resolution dependence is carried by
the fitted kernel width, and the prior stays a real-space operator so the
solver's real-space iterate needs no extra FFTs. No distance between the
simulated and the reconstructed map: `P_G` is scale-free, `B_data` comes
from a regression slope, and `λ_A` from the reconstruction's own moments.
No per-atom soft weights in version 1.
