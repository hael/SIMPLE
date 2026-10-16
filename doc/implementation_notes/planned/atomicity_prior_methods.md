# The atomicity prior: methods, theory and outlook

Methods note, 2026-10-16. It describes the prior as implemented in `simple_pcg_atomic_prior.f90` and the PCG
reconstructor (commit `c59611ed1`, "atom prior for PCG"), relates it to Sayre's equation and to the GENFIRE
family of AET reconstructions, gives the FePt evidence, and looks ahead to biological specimens. The design and
its history are in `single_pcg_integration.md`; the solver is described in `pcg_backend_overview.md`.

## Summary

The reconstruction is the weighted least-squares solution of the projection data on a real-space support. The
prior adds one term. It treats the map as a sum of atom kernels at known positions, with free amplitudes, plus
a residual, and it penalizes only the residual. The penalty is the squared distance of the map from the span of
the atom kernels. It is quadratic, convex and scale-free, and it leaves the solver linear. One strength
parameter moves the solution from the prior-free map (zero) to the least-squares atomic model (infinity).

The prior is what Sayre's equation becomes when the atom positions are known and the amplitudes are free. It
is also a linear replacement for GENFIRE's positivity constraint, and its strong limit is the atomic-model
refinement that ends the AET pipelines. On the FePt tilt series it raised level-1 atom recall from 0.416 to 0.457
at unchanged precision and doubled the Fe atoms found by species discovery. At the higher strengths it made the
atoms rounder than in any map we have, the published GENFIRE map included. Up to the estimated strength it did
this while fitting the projections better than the prior-free map. On tilts
held out of the reconstruction it predicted as well as the best prior-free map, not better.

## 1. The reconstruction problem

### 1.1 Data and forward model

Image $i$ (a particle image in cryo-EM, a projection of a tilt series in AET) is modelled in Fourier space as a
central section of the volume transform:

$$
\hat b_i(\mathbf k) = C_i(\mathbf k)\, e^{-2\pi i\, \mathbf k\cdot \mathbf t_i}\, \hat x\!\left(R_i^{\mathsf T}\mathbf k\right) + \hat n_i(\mathbf k),
$$

with rotation $R_i$, shift $\mathbf t_i$, contrast transfer function $C_i$ (identically 1 for the annular
dark-field STEM tilt series) and noise $\hat n_i$ of variance $\sigma_i^2(|\mathbf k|)$. In SIMPLE the central
section is gathered from a 2× oversampled lattice with a Kaiser-Bessel kernel. Write $K_i$ for this linear map
from the volume to image $i$ and $W_i = \mathrm{diag}(1/\sigma_i^2)$ for the noise weights.

### 1.2 Weighted least squares on a support

The prior-free estimate minimizes

$$
\Phi_0(u) = \sum_i \big\| W_i^{1/2}\,(K_i\, P u - \hat b_i) \big\|^2 + \lambda_0 \|u\|^2 ,
$$

where $P$ is the indicator of the real-space solve support $\Omega$: a sphere, an automatic density envelope,
or a given mask, `pcg_mskfile`. $\lambda_0$ is a small absolute ridge. The normal equations are

$$
(P H P + \lambda_0 I)\, u = P r, \qquad H = \sum_i K_i^{\mathsf T} W_i K_i, \qquad r = \sum_i K_i^{\mathsf T} W_i \hat b_i ,
$$

and the shipped map is $x = w \odot u$, with $w$ the support window.

- $H$ is applied in kernel (Toeplitz) form. The particle pass accumulates, once, the kernel of $H$ on the padded
  lattice: the CTF² and noise weights deposited by the gridding kernel. After that, $Hp$ is a convolution done
  with a few FFTs, whatever the number of images. The kernel's interpolation envelope is divided out
  (deapodization), so $H$ acts on the map itself.
- The Fourier diagonal of $H$ is the sampling density $\rho(\mathbf k)$. Its mean over the native Nyquist sphere,
  zeros included, is $\bar d$ (`realspace_diag`). It stands for the diagonal of $H$ in real space, the data
  precision per voxel, up to a constant factor (the cube corners outside the sphere are left out).

### 1.3 Solver

Preconditioned conjugate gradients, from zero, on the system above.

- **Preconditioner.** The preconditioner is Fourier-diagonal:
  $m(\mathbf k) = 1/(\rho(\mathbf k) + \phi_s)$ on measured voxels and $m(\mathbf k) = 0$ elsewhere.
  - A voxel counts as measured when $\rho \ge 0.1\,\bar\rho_s$, with $\bar\rho_s$ the mean density over the
    sampled voxels of its shell.
  - $\phi_s$ is a shell-relative floor.
- **Support coupling.** The preconditioned residual is $z = P\,\mathcal F^{-1} m\, \mathcal F r$. The final
  multiplication by $P$ is a convolution in Fourier space. It is how unmeasured Fourier voxels (the missing
  wedge, the gaps between tilts) acquire values without ever being given a preconditioner gain.
- **Iterations.** Production runs few iterations: two in refinement, five for final maps. Early stopping is part
  of the prior-free regularization. On FePt the prior-free map predicts held-out tilts best at 10 iterations and
  gets worse after that (section 6).
- **Half maps.** Each half set of images is solved separately. The FSC of the two prior-free halves is the
  resolution authority.

## 2. The prior

### 2.1 Atom model

Given $N$ atom positions $\mathbf r_j$ (from atom detection on the prior-free map, section 3.1), let
$g_j(\mathbf r) = g(\mathbf r - \mathbf r_j)$ be a kernel centred on atom $j$, sampled on the voxel grid and
clipped to $\Omega$. The columns of $G \in \mathbb R^{V\times N}$ are the $g_j$. The model of the map is

$$
x = G c + e ,
$$

with amplitudes $c \in \mathbb R^N$ and a residual $e$. In version 1 every kernel is one isotropic Gaussian,
$g(\mathbf r) = \exp(-|\mathbf r|^2 / 2\sigma^2)$, with $B = 8\pi^2\sigma^2$; its transform falls as
$\exp(-B s^2/4)$ with $s$ the spatial frequency. The width is fitted on the prior-free map (section 3.2).

### 2.2 The penalty

Only the residual is penalized:

$$
Q(x) = \lambda_A \min_{c} \|x - Gc\|^2 = \lambda_A\, x^{\mathsf T} (I - P_G)\, x, \qquad P_G = G\,(G^{\mathsf T} G)^{-1} G^{\mathsf T}.
$$

$P_G$ is the orthogonal projector onto $\mathrm{span}(G)$. The identity follows from
$\|x - Gc\|^2 = \|(I-P_G)x\|^2 + \|P_G x - Gc\|^2$, whose second term vanishes at the minimizing $c$, and from
$I - P_G$ being symmetric and idempotent. Its properties:

- **Positive semi-definite.** $Q$ is convex and adds a symmetric PSD term to the normal operator, so CG applies
  unchanged.
- **Null space $\mathrm{span}(G)$.** Any combination of atom kernels costs nothing. The prior never pulls an
  atom's amplitude toward any value: species, partial occupancy and vacancies are left to the data.
- **Scale-free.** For any positive diagonal $D$, $(GD)\,(DG^{\mathsf T}GD)^{-1}(DG^{\mathsf T}) = P_G$. The model
  contributes positions and kernel shape, never a density scale, so no model map is ever compared with the
  reconstruction in amplitude.
- **Vacuum.** Where no kernel reaches, $P_G x = 0$ and $Q = \lambda_A \|x\|^2$ there: density between atoms is
  pulled to zero by the same term that leaves the atoms free. This is the part of positivity that matters for
  atom detection, kept linear.
- **Overlap handled exactly.** $(G^{\mathsf T}G)^{-1}$ separates overlapping kernels. The off-diagonal entries of
  $G^{\mathsf T}G$ are the overlaps $\langle g_a, g_b\rangle$ of neighbouring footprints. Nothing assumes atoms are
  resolved, but the conditioning degrades as they overlap (section 7).

### 2.3 Bayesian reading

Take the residual as Gaussian white noise of variance $\sigma_e^2$ and the amplitudes as unknown with a flat
prior. Integrating the amplitudes out gives

$$
p(x) \propto \int \exp\!\Big(-\tfrac{\|x - Gc\|^2}{2\sigma_e^2}\Big)\, dc \;\propto\; \exp\!\Big(-\tfrac{x^{\mathsf T}(I - P_G)\,x}{2\sigma_e^2}\Big),
$$

an improper Gaussian prior whose precision is $(I-P_G)/\sigma_e^2$. With the data likelihood of section 1.2,
whose precision is $H$ in the units of the normal equations, the maximum a posteriori map solves

$$
\big(P H P + \lambda_0 I + \lambda_A (I - P_G)\big)\, u = P r, \qquad \lambda_A = \frac{\sigma_n^2}{\sigma_e^2},
$$

where $\sigma_n^2$ is the noise variance per unit data weight. It is 1 when $W$ holds the inverse noise variances.
With unit weights, as in the FePt run, it is the unknown noise level, which the estimate of section 2.5 supplies.
The strength is a ratio of variances, as in any Wiener filter: data noise over the variance of the density the
atoms do not explain.

### 2.4 The two limits and what lies between

- **$\lambda_A = 0$.** The prior-free solve.
- **$\lambda_A \to \infty$.** The solution is forced into $\mathrm{span}(G)$, and the data term picks the
  amplitudes:

  $$
  x_\infty = G\hat c, \qquad \hat c = \arg\min_c \big\| W^{1/2}(K G c - \hat b) \big\|^2 \;\Longleftrightarrow\; G^{\mathsf T} H G\, \hat c = G^{\mathsf T} r .
  $$

  This is an atomic model whose amplitudes are fitted to the images, at fixed positions. Here and below $H$
  includes the small ridge $\lambda_0$, as the code's operator does.
- **In between, with white data.** Suppose $H = \bar d\, I$ (complete, uniform sampling), $P = I$, and neglect
  $\lambda_0$. Then $P_G$ commutes with the operator and the solution separates exactly:

  $$
  x_\lambda = P_G\, x_0 + \frac{1}{1 + \lambda_{\mathrm{rel}}}\,(I - P_G)\, x_0, \qquad \lambda_{\mathrm{rel}} = \lambda_A / \bar d ,
  $$

  with $x_0$ the prior-free map. The part of the map the atoms explain passes through untouched. The rest is
  shrunk by a Wiener factor. This is a real-space analogue of map combination in crystallography from a partial
  model: the model-explained part kept, the remainder weighted by its reliability.
- **In between, with incomplete data.** With a missing wedge, $H$ has a near null space: the unmeasured Fourier
  region. The prior removes it. Take $P = I$ and an exact null space. In the unmeasured directions the data say
  nothing, and stationarity of the objective requires the null-space component of $(I - P_G)x$ to vanish. So in
  the unmeasured region the transform of the solution equals that of $P_G x$, the atom model of the solution
  itself, whose amplitudes are set jointly with the measured part. With a support the statement holds
  approximately. In Fourier terms an element of $\mathrm{span}(G)$ is

  $$
  \hat x(\mathbf k) = \hat g(\mathbf k) \sum_{j=1}^{N} c_j\, e^{-2\pi i\, \mathbf k\cdot \mathbf r_j},
  $$

  defined for every $\mathbf k$ by $N$ numbers. Once the measured region constrains the $c_j$, the formula
  extrapolates into the wedge. The count is what makes this possible. FePt has about $2.3\times10^4$ atoms, each with one
  amplitude (four parameters with positions), against $68 \times 256^2 \approx 4.5\times10^6$ measured pixels.

### 2.5 Strength

The strength is set relative to the data precision per voxel:

$$
\lambda_A = \lambda_{\mathrm{rel}}\, \bar d .
$$

The earlier priors of the PCG backend anchored their strength to the mean data precision of the lowest shells.
Their effective strength then depended on how steeply the data precision falls across the band, which is a
property of each data set, and they did not carry across data sets (`pcg_priors_history.md`). A real-space term
competes with the real-space diagonal of $H$, for which $\bar d$ stands.

The estimate of $\lambda_{\mathrm{rel}}$ follows section 2.3. With $H \approx \bar d I$, the noise of a prior-free
map per voxel is $\mathcal N = \sigma_n^2/\bar d$, so $\lambda_{\mathrm{rel}} = \mathcal N/\sigma_e^2$. The residual
variance is the signal not explained by the atoms, $\sigma_e^2 = S(1 - C^2)$. All moments are means over the
support $\Omega$, computed on the prior-free pair $x_e$, $x_o$ and their average $x_b$:

$$
\mathcal N = \mathrm{var}\Big(\frac{x_e - x_o}{2}\Big), \quad
S = \max\big(\overline{x_b^2} - \mathcal N,\; 10^{-3}\mathcal N\big), \quad
E' = \max\Big(\overline{(P_G x_b)^2} - \mathcal N \frac{N}{|\Omega|},\; 0\Big), \quad
C^2 = \min(E'/S,\; 0.98),
$$

$$
\lambda_{\mathrm{rel}} = \frac{2\mathcal N / S}{1 - C^2}.
$$

- $\mathcal N$ is the noise of the merged map.
- $\mathcal N N/|\Omega|$ removes the noise that falls into the $N$-dimensional span of the kernels.
- The factor 2 is there because each half is solved separately, with twice the merged map's noise. The merged
  ratio $(\mathcal N/S)/(1-C^2)$ is logged beside it.

This is a sigma-A-like estimate, built from the map's own moments; no constant is tuned. Read's (1986) sigma-A
comes from a per-shell likelihood fit of observed against model structure factors; ours is one global
real-space ratio. The note's
acceptance rule is judged on held-out data: the estimate's held-out error must lie within 1% of the best on a
grid of strengths, on every test set. The rule is stated on the error, not the strength, because the optimum
is flat.

## 3. Ingredients as implemented

### 3.1 Positions

- **Source.** The positions are those of `detect_atoms` with species discovery (its `_species.pdb`), run on the
  prior-free merged map. PDB coordinates map to voxels as Å = 0-based voxel × sampling distance.
- **Cleaning.** Atoms whose centre voxel lies outside the support, and atoms closer than 0.5 Å to an earlier
  atom, are dropped.
- **Errs on inclusion.** A spurious kernel costs little: the data give it an amplitude near zero and the prior
  then does not penalize density there. A missing kernel pulls a real atom toward zero. So the model errs on
  inclusion, and the recovered atoms of discovery are part of it.

### 3.2 Kernel width

On the prior-free merged map $x_b$, a reference model $m = \sum_j g_{\sigma_f}(\mathbf r - \mathbf r_j)$ is built
with unit amplitudes and $\sigma_f = 1.5$ voxels, on a $5\sigma_f$ footprint. Per Fourier shell $s$:

$$
c_s = \frac{\operatorname{Re}\langle \hat m_s, \hat x_s\rangle}{\|\hat m_s\|^2}, \qquad
\log c_s = \log K - \frac{\Delta B}{4}\, s^2 .
$$

- **Why this slope.** The per-shell regression slope is unbiased under map noise. A map that is the model
  blurred by an extra $\Delta B$ gives exactly this line, for exact positions and equal amplitudes. With two
  species against a unit-amplitude model, $K$ varies somewhat with $s$.
- **The fit.** It runs over the shells where model and map correlate at 0.5 or more (at least three shells), and
  $K$ is discarded.
- **The kernel.** $B = B_f + \Delta B$, where $B_f$ is the reference model's B. It is reported as
  $B_{\mathrm{data}} = \max(0, B - B_{\min})$, with $B_{\min}$ the B of a one-voxel Gaussian.
- **Why a 1.5-voxel reference.** It is wide enough that the voxel lattice samples it without aliasing. Its
  untruncated footprint keeps a cut-off plateau out of its spectrum. The first version, which used kernels at
  $B_{\min}$ as the reference, recovered $B_{\mathrm{data}}$ 6% low in the unit test. Part of that came from the
  test's own truth, whose footprints were cut at 3σ and are now untruncated. The fit is now 0.4% off.
- **Footprints.** Kernel footprints are cubes of half-width $\mathrm{clamp}(\lceil 3\sigma\rceil, 4, 8)$ voxels.

On FePt: $B_{\mathrm{data}} = 17.2$ Å², $\sigma = 1.60$ voxels (0.60 Å), fitted over 24 shells; the total kernel
B is 28 Å². The position refinement of the FePt paper used B factors of 20 to 22 Å², but on a Kirkland form
factor, so the two numbers are not on the same baseline.

### 3.3 Applying $P_G$

- **The Gram matrix.** $G^{\mathsf T}G$ is sparse: only atoms whose footprint cubes overlap are coupled. It is
  built once per model from a cell grid of the footprint size and stored in compressed rows.
- **Applying it.** $P_G p$ is a gather ($G^{\mathsf T}p$), a solve with $G^{\mathsf T}G + \varepsilon D$
  ($D$ its diagonal, $\varepsilon = 10^{-7}$) by Jacobi-preconditioned CG in double precision to a relative
  residual of $10^{-10}$, and a scatter. The scatter runs serially, so the sum order is fixed.
- **Cost.** Each application of $\lambda_A(p - P_G p)$ is one gather and one scatter (about $N \times$ footprint
  operations each), an inner Gram CG of some tens of iterations over the non-zeros of $G^{\mathsf T}G$, and one
  update of the whole box. It needs no FFT.
- **Where it sits.** It enters the normal operator in `add_ridge`, after the support mask, where both operator
  modes (kernel and matrix-free) pass. Strength 0 attaches nothing and reproduces the prior-free solve bit for
  bit. On the command line `lambda_atm=0` means "use the estimate"; `atm_prior=no` is the prior-free solve.

### 3.4 Amplitudes

$\hat c$ of section 2.4 is computed per half by CG in atom space on $G^{\mathsf T}HG\,c = G^{\mathsf T}r$. It uses
the half's own normal operator ($\lambda_0$ included) and a Jacobi preconditioner on $\mathrm{diag}(G^{\mathsf T}G)$. Each iteration
is one application of $H$. The prior does not use $\hat c$, since $P_G$ is amplitude-free. $\hat c$ feeds the
residual map $x_b - G\hat c$ and is the intended intensity measure for species. In the numpy prototype on FePt,
amplitudes fitted this way separated Fe from Pt without supervision. Amplitudes measured on the map,
$(G^{\mathsf T}G)^{-1}G^{\mathsf T}x$, did not.

### 3.5 Solver with the prior

After the prior-free base pair, once per state:

1. Build the model, fit the width on the merged base map, and estimate $\lambda_{\mathrm{rel}}$ from the pair.
2. Per half, fit the amplitudes through the half's operator, before the prior is attached.
3. Attach the prior to both half operators.
4. Solve each half again from zero.

- **Preconditioner.** $\lambda_A/\mathrm{padsc}^2$ enters the denominator of the measured voxels only:
  $m = 1/(\rho + \phi_s + \lambda_A/\mathrm{padsc}^2)$.
  - Unmeasured voxels keep $m = 0$.
  - The alternative, the gain $1/\lambda_A$ on unmeasured voxels, stalled CG at small strengths and gave worse
    maps at every strength in the numpy prototype.
  - The wedge fills through the support coupling of section 1.3 instead.
- **Budget.** A fixed 30 iterations (`maxits_atm`), with a relative-residual stop available (`rtol_atm`). On FePt
  the CG residual is $1.5\times10^{-3}$ at 30 iterations and $1.2\times10^{-3}$ at 60.
  - Between 30 and 60 iterations, level-1 and discovery recall move by 1 to 2%.
  - The shape, Fe and species scores move by 4 to 11%, in both directions.
  - The map still changes by about 1% per iteration at 60: the wedge fills slowly through the support.
- **Iterations do not degrade it.** The prior-free solve gets clearly worse from 5 to 30 iterations. With the
  prior, 30 iterations are the best tested and 60 change the scores both ways.
- **Scope.** The prior requires the deapodized operator and excludes the ML prior, the solvent ridge and NU
  filtering in version 1.

### 3.6 Products and validation rules

- **Products.**
  - The regularized pair ships.
  - The prior-free pair is written under `_base` and carries the FSC.
  - The residual $x_b - G\hat c$ is written under `_res`.
  - A sidecar `_atmprior.txt` records the model, width, moments, strength and solver statistics.
- **Validation rules.** Both regularized halves are pulled toward the same atoms, so:
  - their FSC is inflated and is never used as evidence;
  - their difference is no noise reference: species discovery run with the regularized halves overcalls
    (6,298 atoms in its first stage on FePt) and must use the `_base` halves;
  - evidence comes from ground truth, from held-out images, or from a published model.

## 4. Relation to Sayre's equation

### 4.1 Sayre's equation

Sayre (1952) assumed a crystal of $N$ equal atoms that are resolved, so that their densities do not overlap:
$\rho(\mathbf r) = \sum_j g(\mathbf r - \mathbf r_j)$. Squaring changes the shape of each atom and nothing else:
$\rho^2(\mathbf r) = \sum_j g^2(\mathbf r - \mathbf r_j)$.

Let $f$ and $f^{(2)}$ be the transforms of $g$ and $g^2$, and
$Z(\mathbf h) = \sum_j e^{2\pi i\, \mathbf h\cdot\mathbf r_j}$ the structure factor of the positions. This section
uses the crystallographic sign convention, the opposite of section 2.4. Then
$F_{\rho}(\mathbf h) = f(\mathbf h) Z(\mathbf h)$ and $F_{\rho^2}(\mathbf h) = f^{(2)}(\mathbf h) Z(\mathbf h)$, and by
the convolution theorem $F_{\rho^2}(\mathbf h) = V^{-1}\sum_{\mathbf k} F(\mathbf k) F(\mathbf h - \mathbf k)$. Hence

$$
F(\mathbf h) = \frac{\theta(\mathbf h)}{V} \sum_{\mathbf k} F(\mathbf k)\, F(\mathbf h - \mathbf k), \qquad \theta = f / f^{(2)} .
$$

Two features made it one of the roots of direct methods. It needs no positions: it holds for any arrangement of equal
resolved atoms. And it couples every structure factor to products of others, so known ones constrain unknown
ones.

### 4.2 Both state atomicity, under opposite assumptions

In real space Sayre's equation reads $\rho = \check\theta * \rho^2$, with $\check\theta = \mathcal F^{-1}[\theta]$. For
Gaussian atoms $\theta = f/f^{(2)}$ decays with $s$, so $\check\theta$ is a smoothing kernel; it turns $g^2$ back into
$g$. Take a map in the span of ideal kernels (not clipped to a support) with non-overlapping footprints,
$\rho = \sum_j c_j g_j$. Then $\rho^2 = \sum_j c_j^2 g_j^2$ and $\check\theta * \rho^2 = \sum_j c_j^2 g_j$. Sayre's
equation therefore requires $\sum_j (c_j^2 - c_j)\, g_j = 0$. Because the kernels are linearly independent (the
condition under which $P_G$ exists), it holds if and only if $c_j \in \{0, 1\}$ for every $j$. If $\theta$ is scaled
to atoms of amplitude $a$, the condition is $c_j \in \{0, a\}$.

Non-overlap is needed only for the squaring step, and for Gaussians it holds only approximately. So:

- **Sayre's equation is atomicity with free positions and fixed, equal amplitudes.** The density is a sum of
  identical atoms, each present or absent, at positions nobody has to know. The constraint is quadratic in the
  density.
- **Our prior is atomicity with fixed positions and free amplitudes.** The density is a sum of kernels at known
  positions, with any amplitudes. The constraint is linear in the density.

The two are complementary relaxations of the statement "the density is a sum of atoms". Sayre trades knowledge of
the positions for a nonlinear constraint, which is why it can phase a structure ab initio and why it needs
equal, resolved atoms. Extensions to unequal atoms are approximate. We trade the ab initio power for linearity.
The positions come from the data, through atom detection on the prior-free map. In exchange the formulation:

- handles two or more species, partial occupancy and vacancies, since the amplitudes are free;
- separates overlapping atoms exactly through $(G^{\mathsf T}G)^{-1}$;
- keeps a residual for what the model lacks;
- has a convex optimum computed by the same linear solver as the prior-free map.

### 4.3 Filling unmeasured data

In crystallography Sayre's equation extends phases: measured amplitudes and some known phases determine
others. Density modification at atomic resolution goes further, extrapolating structure factors beyond the
measured resolution limit, the "free lunch" (Caliandro et al., 2005).

In tomography the phases are measured. What is missing is a region of Fourier space: the wedge outside the
tilt range, the gaps between tilts, and the space beyond Nyquist of the projections. The convolution in
Sayre's equation is what fills such regions. If $\mathbf h$ lies in the wedge, there are pairs
$\mathbf k$, $\mathbf h - \mathbf k$ both outside it, and their product predicts $F(\mathbf h)$. Our prior does the
same through the explicit structure factor
$\hat x(\mathbf k) = \hat g(\mathbf k)\sum_j c_j e^{-2\pi i\,\mathbf k\cdot\mathbf r_j}$ of section 2.4.

### 4.4 Partial-model phasing

Our prior relates to Sayre's equation as phasing from a partial model relates to direct methods. Once part of a
structure is known, the best map combines model and data with weights from the agreement between them (sigma-A,
Read 1986), and the free R value (Brünger, 1992) tells whether the model is helping or biasing. Our strength
estimate is that weight in real space (section 2.5). The held-out image check is the free R.

The loop of the SINGLE workflow is the same loop as automated model building (ARP/wARP; Perrakis et al.,
1999):

1. detect atoms on the prior-free map;
2. solve with the prior;
3. detect again.

## 5. Relation to GENFIRE and its successors

### 5.1 What GENFIRE does

GENFIRE (Pryor et al., 2017) produced the published FePt map. It works in two stages; we read both from the FePt
code.

- **Gridding.** Each projection is zero-padded four times. For every grid point within 0.05 grid units of a
  projection's central plane, the exact discrete transform of that projection is evaluated at the foot of the
  perpendicular, and values from several projections are averaged with inverse-distance weights. This gives
  measured values $m$ on a set $M$ of grid points; every other point is unmeasured.
- **Iteration.** With $\mathcal C$ the projection onto non-negative maps inside the support, 500 iterations of

  $$
  x_{n+1} = z_n - \mathcal F^{-1} M (\mathcal F z_n - m), \qquad z_n = \mathcal C(x_n),
  $$

  keeping the iterate with the smallest Fourier error at the measured points.

The update is a unit gradient step on the gridded misfit $\tfrac12\|M(\mathcal F z - m)\|^2$, taken from the
point that satisfies the constraints. GENFIRE is projected gradient descent. It is the error-reduction algorithm
of phase retrieval (Gerchberg and Saxton's, as generalized and read as steepest descent by Fienup, 1982), with a
tomographic data constraint.

Its successor RESIRE (Pham et al., 2023) minimizes the misfit to the projections themselves by fixed-step
gradient descent, with positivity after each step. Neither carries noise weights. The FePt pipeline then
traced atoms on the GENFIRE map and refined their positions against the projections. That refinement used one
form factor with a per-species amplitude and B factor, plus an embedded-atom potential as a restraint (Yang et
al., 2017).

### 5.2 The same problem, three ways

With the gridding approximation, PCG's operator becomes the Fourier-diagonal density $\rho$ and its right-hand
side the gridded backprojection $\tilde r$. Take $P = I$, with the support acting as a separate projection like
GENFIRE's $\mathcal C$. One preconditioned gradient step of unit size is then

$$
u_{n+1} = u_n - \mathcal F^{-1}\, \frac{\rho}{\rho + \phi}\, \big(\mathcal F u_n - \tilde r/\rho\big),
$$

which is GENFIRE's update with $M$ replaced by the smooth $\rho/(\rho+\phi)$ and without positivity. The prior-free
PCG solve improves on that iteration in four ways:

- conjugate directions with optimal steps;
- the exact coupled operator instead of the diagonal approximation;
- noise and CTF weights;
- a statistical regularizer.

What GENFIRE has and PCG lacked is a constraint that fills unmeasured regions: positivity.

The atomicity prior is that constraint in linear form. Both act in the missing wedge, and they differ in kind:

- **Positivity is a cone of dimension $V$.** For an object that fills its support, as a nanoparticle does, the
  set of positive maps consistent with the data is large. Positivity says nothing about density between atoms
  as long as it is positive. On FePt a numpy GENFIRE-style iteration on our data made the atoms rounder
  (elongation 2.02 to 1.90, 200 iterations) but lifted no more atoms above the background, and level-1 recall
  fell from 0.30 to 0.24 (`aet_codes_and_simple.md`, section 3.5). That test used the nominal tilts and a
  spherical support, which is why its recall is below the 0.416 of section 6.
- **Atomicity is a subspace of dimension $N \ll V$ plus a penalized residual.** It pulls interstitial density to
  zero, extrapolates into the wedge through the structure factor of known positions, and lets the data set every
  amplitude.

The strong limit of the prior is the AET pipelines' last step: an atomic model fitted to the projections. Ours
fits amplitudes at fixed positions; theirs fits positions with per-species amplitudes. The limit is reached in
the same linear solve that makes the map. So one quadratic problem spans their reconstruction (at $\lambda_A = 0$,
minus positivity) and their model refinement (at $\lambda_A \to \infty$, minus position updates), with a map at
every strength in between that is atomic where the data support atoms and free where they do not.

### 5.3 What we have not taken from them

- **Positions refined against the projections, and an interatomic potential as a restraint.** Position refinement
  is an open item (section 8). In the numpy prototype, Gauss-Newton position steps against the
  projections cut the false positives of the amplitude fit, which are real atoms whose detected position is off
  by 1 to 2 Å. A potential would write the physics into the answer: FePt chemical order is the quantity the 2017
  paper measured.
- **Positivity on top of the prior.** It would need a proximal step instead of plain CG. On FePt the prior at
  ½ × the estimate and above already gives rounder atoms than positivity did, though in a different set-up.

## 6. FePt evidence

### 6.1 Set-up

- **Data.** The FePt particle of Yang et al. (2017): 68 denoised projections from −65° to +64°, 0.3725 Å per
  pixel, padded to a 320 box.
- **Poses.** The refined poses of `refine_tomo` (one stage, 2 Å, two iterations), held fixed.
- **Support.** The density mask `tomo_mask.mrc`.
- **Base map.** Masked PCG, five iterations, which is what `refine_tomo` ships.
- **Model.** `detect_atoms` with discovery on the base map, 13,103 atoms.
- **Strength.** The estimate is $\lambda_{\mathrm{rel}} = 0.415$; the grid scales it by ⅛ to 2. All prior
  solves run 30 iterations from zero.
- **Scoring.** Against the published model of 23,196 atoms, with the same detector for every map:
  - level-1 recall and precision of `detect_atoms`;
  - surface atoms, those with fewer than 12 neighbours within 3.3 Å;
  - the discovery columns, which use the prior-free halves as noise reference;
  - Pt and Fe above the interstitial background, the share of atoms whose peak exceeds its 99th percentile;
  - elongation, the long-to-short axis ratio of the mean Pt atom;
  - the data error $|a\,Px - y|/|y|$ over the 68 projections, with one global scale.

The run's report is `~/agent_runs/runs/atm_prior_mac/report.md`.

### 6.2 Results

| Map | Level-1 recall (precision) | Surface | Discovery recall | Discovery Fe | Species agreement | Pt / Fe above background | Elongation | Data error |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| prior-free, 5 its | 0.416 (0.983) | 0.144 | 0.554 | 0.081 | 0.29 | 0.91 / 0.21 | 2.07 | 0.0690 |
| prior-free, 30 its | 0.372 (0.972) | 0.162 | 0.493 | 0.065 | 0.14 | 0.91 / 0.21 | 2.08 | 0.0691 |
| prior, ⅛ × estimate | 0.445 (0.983) | 0.191 | 0.603 | 0.150 | 0.44 | 0.93 / 0.27 | 1.99 | 0.0663 |
| prior, ¼ × | **0.457** (0.985) | **0.196** | **0.618** | 0.184 | 0.53 | 0.93 / 0.31 | 1.93 | **0.0658** |
| prior, ½ × | 0.450 (0.987) | 0.185 | 0.602 | 0.190 | 0.60 | 0.93 / 0.33 | 1.86 | 0.0662 |
| prior, 1 × | 0.442 (0.986) | 0.182 | 0.595 | 0.181 | 0.60 | 0.92 / 0.31 | 1.80 | 0.0674 |
| prior, 2 × | 0.434 (0.986) | 0.176 | 0.581 | 0.154 | 0.58 | 0.90 / 0.27 | **1.74** | 0.0696 |

For reference, their GENFIRE map gives level-1 recall 0.41 with the same detector and elongation 1.89.

Held out: every fourth tilt left out, maps from the other 51, model and estimate rebuilt from them (estimate
0.388). All held-out maps share the full-data `tomo_mask.mrc` as support, a small leak that is the same for
every map:

| Map | Held-out error | In-sample error |
| --- | --- | --- |
| prior-free, 3 / 5 / 10 / 20 / 30 its | 0.0917 / 0.0862 / **0.0852** / 0.0867 / 0.0890 | 0.0783 / 0.0711 / 0.0701 / 0.0702 / 0.0714 |
| prior, ⅛ / ¼ / ½ / 1 / 2 × estimate | 0.0867 / 0.0858 / 0.0854 / 0.0862 / 0.0881 | 0.0686 / 0.0678 / 0.0679 / 0.0690 / 0.0707 |

### 6.3 What the numbers say

- **More atoms, at no cost in fit.** At every strength from ⅛ to 1 × the estimate, the prior:
  - raises level-1 recall at unchanged or better precision;
  - raises discovery recall, with precision unchanged or better except at ⅛ × (0.979 against 0.981);
  - raises surface recall by a quarter to a third;
  - nearly doubles the Fe atoms discovery finds at ⅛ ×, and more than doubles them from ¼ × on;
  - fits the projections better than the prior-free map.

  The gain is not from iterations: the prior-free solve at 30 iterations is worse than at 5.
- **The strength trades detection against shape.** Detection peaks at ¼ × the estimate. Roundness improves
  steadily with strength. At 2 × the estimate the atoms are the roundest of the 30-iteration maps, rounder than
  in GENFIRE's map. But recall, surface recall and Fe calls have come back down, and the data error is above
  the prior-free map's.
- **Held-out tilts.** The best regularized map matches the best prior-free map on tilts it never saw, 0.0854
  against 0.0852, and the prior-free map reaches that only at its best stopping point. The 2 × map is 3.4%
  further from the unseen data than the best prior-free map, and 3.1% further than the best regularized one.
  The estimate lies 0.8% above the grid's best, within the 1% acceptance rule on this one test set.
- **Species are not solved.** The agreement column counts the matched atoms whose species call is right. The
  matched atoms are 91 to 96% Pt, so calling every one Pt would score 0.91 to 0.96. Discovery scores 0.29
  prior-free and 0.44 to 0.60 with the prior at 30 iterations (0.63 at 60): it still splits the Pt atoms into
  two classes instead of finding Fe. The data-fitted amplitudes of section 3.4, which the prototype showed to
  split Fe from Pt, are not yet written per atom or used by discovery.
- **The residual points at missing atoms.** Of the 10,343 reference atoms the model lacks, the same number of
  top peaks of $x_b - G\hat c$ falls on a missing atom in 61% of cases.

### 6.4 The 2 × map

The 2 × map is the most atomic-looking map, because it is the closest to the model. Of the 30-iteration maps it
has the roundest atoms (1.74), the highest Pt peak over interstitial noise (6.49, against 4.45 prior-free) and
the fewest negative voxels (6%). Yet by the share of Pt atoms above the 99th percentile of the interstitial
density it is the lowest (0.90, against 0.91 prior-free). That look is partly the prior's own: at that strength
more of what the data say beyond the model is suppressed. The held-out error and the falling recall agree on
this.

For analysis (detection, species, coordinates) the operating point is ¼ to ½ × the estimate, where held-out
prediction, detection and data fit are best together. The 2 × map is a fair figure of the structure, labelled
with its strength, in the way a sharpened map is a fair figure of a cryo-EM structure.

### 6.5 Independent checks of the same claims

A numpy prototype on the same data (nominal tilts, slab-separable 2D solves,
`~/agent_runs/runs/atm_prior_mac/numpy_proto/RESULTS.md`) found the same direction:

- the data fit was unchanged;
- with its own peak finder, recall at 90% precision rose from 0.72 to 0.89 and Pt above background from 0.29 to
  0.81 (partly circular, since its model came from the same data);
- elongation fell from 2.00 to about 1.45;
- the regularized maps predicted held-out tilts slightly better than the best prior-free map (0.117 against
  0.121).

On a synthetic tilt series with exact atoms, noise matched to the real map, an atom shape unlike the kernel, and
analytic projections, recall at 98% precision rose from 0.16 to 0.75. Amplitudes fitted to the projections, with
the positions refined too, called 93% of Fe atoms correctly at 97% precision.

Two differences from the Fortran run may explain its smaller held-out gain: the prototype used nominal tilts,
whose prior-free map over-fits worse, and it built its model from amplitudes and positions fitted to the data.
Neither has been tested in the Fortran run yet.

## 7. Biological specimens in cryo-EM

### 7.1 What carries over unchanged

- **The data term and the solver.** The operator already carries CTFs and noise weights. The prior is a
  real-space term and does not see them.
- **The rules.** The resolution comes from the prior-free pair, and the regularized halves are not compared with
  each other. The held-out check becomes easier: single-particle data sets have many more images than a tilt
  series. A free set of a few percent of the particles, held out of the reconstruction and predicted from it,
  is the free R of cryo-EM, and it costs nothing.
- **The strength.** The same anchor and estimate, which were built to carry across data sets.

### 7.2 What changes

- **Kernels.** The electron scattering factor of each element (the five-Gaussian fits of Peng et al., 1996),
  broadened by a B, replaces the single Gaussian. Since the prior is scale-free, only the kernel shape matters.
  A smooth B field, rather than one B or free per-atom B, follows the LocScale experience.
- **Solvent and the zero level.** Vitreous ice has density, and phase contrast does not measure the mean of the
  map. Pulling everything outside the atoms to zero would bias the molecular surface. A few unpenalized
  smooth basis vectors (a constant over the solvent region, or the solvent ridge's weight field) join $G$, as
  StatSTEM includes a constant background with its atomic columns.
- **Where the model comes from.**
  - A model built automatically from the prior-free map (for example with ModelAngelo; Jamali et al., 2024).
  - A deposited or predicted model fitted to it.
  - Partial models: only built residues enter $G$, and the residual keeps the rest.
- **Resolution.** At 1.2 to 1.5 Å (apoferritin; Nakane et al., 2020; Yip et al., 2020) atoms are separated and
  the formulation applies as on FePt. Between about 1.5 and 2.5 Å, kernels of B 30 to 80 Å² overlap.
  $G^{\mathsf T}G$ stays solvable, with the existing ridge and footprints that grow with B, but the span of $G$
  becomes smooth and the prior acts mainly on the highest frequencies. Beyond about 3 Å, atoms are the wrong
  basis. One Gaussian per residue, or a Gaussian mixture fitted to the map (Chen & Ludtke, 2021), keeps the form
  of the prior. Sayre's equation is the reminder of where the strong version holds: direct methods need data to
  about 1.2 Å.

### 7.3 Where it should help most

- **Missing cones.** The clearest change the prior made on FePt was to atom shape along the beam, the direction
  of the missing wedge (no measurement split its effect by Fourier region). Preferred orientation in
  single-particle cryo-EM produces the same defect: a cone of directions sampled thinly or not at all. The
  influenza hemagglutinin trimer is the standard case (Tan et al., 2017). An atomic model of the well-sampled
  directions predicts the poorly sampled ones through the same structure factor.
- **Electron cryo-tomography and subtomogram averaging.** Every tomogram has a missing wedge. Subtomogram averages
  of molecules with known structure could use the known model as $G$ to fill the wedge of each average.
- **Small and low-contrast particles.** Learned denoisers inside refinement (Blush; Kimanius et al., 2024) show
  how much a strong prior can lift alignment of small particles. An explicit model prior is the transparent
  counterpart wherever a model exists: the residual map shows what the model does not explain, and the held-out
  particles measure whether it helps.
- **Waters, ions and ligands at atomic resolution.** These are the weakest atoms. The prior suppresses the noise
  around them and the residual shows them. They are also where model bias is most dangerous (section 7.4).

### 7.4 Risks and how to keep them in check

- **Model bias.** Wherever $\lambda_A$ is large compared with the data precision, the map follows the model.
  The safeguards are those of crystallography:
  - build the model from the prior-free map;
  - judge every strength and every model change on held-out particles, never on the FSC of regularized halves;
  - show the residual map beside the regularized one;
  - add atoms to the model only when the held-out error falls.
- **Local strength.** Flexible regions and poorly built parts invite a per-region strength. The NU-evidence prior
  of the PCG backend did exactly that, worked on its test phantom, and was removed because its operating point
  did not carry across data sets. A local strength would have to be set per region by held-out prediction, not
  by a map statistic.
- **Heterogeneity.** A mixture of conformations is not a sum of atoms at one set of positions. The prior belongs
  after classification, per class, with the B field absorbing the remaining spread.

### 7.5 Relation to cryo-EM practice

- **Density modification for cryo-EM** (Terwilliger et al., 2020). It brought the crystallographic toolkit to
  half maps: solvent flattening and histogram matching, cross-validated between halves, without a model. It runs
  after reconstruction.
- **LocScale.** It scales local amplitudes to a model or pseudo-model reference after reconstruction.
- **Blush.** Its prior is learned and implicit, inside refinement.
- **This prior.** It is explicit and physical, inside the reconstruction. It is linear, so its optimum is unique
  and computed by the existing solver, and it shows its own residual.

### 7.6 First experiments in SIMPLE

1. An atomic-resolution apoferritin data set. The model is built from the prior-free map, a 5% free set of
   particles is held out, and the strength is scanned against held-out error. The test: waters and ions in the
   residual, held-out error at the estimate within 1% of the best.
2. A preferred-orientation data set (hemagglutinin). The test: the anisotropy of the map with the prior against
   tilted-data collections of the same specimen, which measure the missing cone.
3. The same code on SINGLE Pt particle data: the setting between AET and cryo-EM, with TEM kernels and many
   images.

## 8. Open items

- Model from amplitudes and positions fitted to the data, the prototype's version, against `detect_atoms`
  positions. This is the likely source of the remaining gap to the prototype.
- Per-atom amplitudes written out and used by species discovery.
- Position refinement against the images.
- Integration into `refine_tomo` and `refine3D_nano` / `autorefine3D_nano`.
- Peter Ercius's multislice simulations as exact truth for AET.
- The TEM kernel and the bulk-solvent basis for cryo-EM.

## References

- Brünger, A. T. (1992). Free R value: a novel statistical quantity for assessing the accuracy of crystal
  structures. *Nature* 355, 472–475.
- Caliandro, R., Carrozzini, B., Cascarano, G. L., De Caro, L., Giacovazzo, C. & Siliqi, D. (2005). Phasing at
  resolution higher than the experimental resolution. *Acta Cryst.* D61, 556–565.
- Chen, M. & Ludtke, S. J. (2021). Deep learning-based mixed-dimensional Gaussian mixture model for
  characterizing variability in cryo-EM. *Nat. Methods* 18, 930–936.
- De Backer, A., van den Bos, K. H. W., Van den Broek, W., Sijbers, J. & Van Aert, S. (2016). StatSTEM: an
  efficient approach for accurate and precise model-based quantification of atomic resolution electron
  microscopy images. *Ultramicroscopy* 171, 104–116.
- Fienup, J. R. (1982). Phase retrieval algorithms: a comparison. *Appl. Opt.* 21, 2758–2769.
- Jamali, K. et al. (2024). Automated model building and protein identification in cryo-EM maps. *Nature* 628,
  450–457.
- Kimanius, D. et al. (2024). Data-driven regularization lowers the size barrier of cryo-EM structure
  determination. *Nat. Methods* 21, 1216–1221.
- Nakane, T. et al. (2020). Single-particle cryo-EM at atomic resolution. *Nature* 587, 152–156.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive
  electron atomic scattering factors. *Acta Cryst.* A52, 257–276.
- Perrakis, A., Morris, R. & Lamzin, V. S. (1999). Automated protein model building combined with iterative
  structure refinement. *Nat. Struct. Biol.* 6, 458–463.
- Pham, M., Yuan, Y., Rana, A., Osher, S. & Miao, J. (2023). Accurate real space iterative reconstruction
  (RESIRE) algorithm for tomography. *Sci. Rep.* 13, 5624.
- Pryor, A. et al. (2017). GENFIRE: a generalized Fourier iterative reconstruction algorithm for high-resolution
  3D imaging. *Sci. Rep.* 7, 10409.
- Read, R. J. (1986). Improved Fourier coefficients for maps using phases from partial structures with errors.
  *Acta Cryst.* A42, 140–149.
- Sayre, D. (1952). The squaring method: a new method for phase determination. *Acta Cryst.* 5, 60–65.
- Tan, Y. Z. et al. (2017). Addressing preferred specimen orientation in single-particle cryo-EM through tilting.
  *Nat. Methods* 14, 793–796.
- Terwilliger, T. C., Ludtke, S. J., Read, R. J., Adams, P. D. & Afonine, P. V. (2020). Improvement of cryo-EM
  maps by density modification. *Nat. Methods* 17, 923–927.
- Yang, Y. et al. (2017). Deciphering chemical order/disorder and material properties at the single-atom level.
  *Nature* 542, 75–79.
- Yip, K. M., Fischer, N., Paknia, E., Chari, A. & Stark, H. (2020). Atomic-resolution protein structure
  determination by cryo-EM. *Nature* 587, 157–161.
