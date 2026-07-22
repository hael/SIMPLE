# Diffusion maps: lessons from a ManifoldEM review

## Status

**Implementation note — no code changes are implied by this document.**

This note records a systematic comparison between SIMPLE's diffusion-map
subsystem and the ManifoldEM Python package
(`/home/elmlundho/src/ManifoldEM/ManifoldEM`; original MATLAB by A. Dashti /
UWM 2016, Python port by H. Liao and E. Seitz / Columbia 2018–2020).  Its
purpose is to document algorithmic differences that could improve the SIMPLE
implementation and to provide exact pointers to both codebases so that any
future developer can trace the rationale without repeating the review.

---

## 1. SIMPLE files reviewed

| File | Role |
|---|---|
| `src/main/pca/simple_diff_map_graphs.f90` | CSR graph construction, Gaussian kernel, symmetric normalization |
| `src/main/pca/simple_diffusion_maps.f90` | Embedding driver, steerable DM, Nyström, eigengap selection |
| `src/main/pca/simple_diff_map_denoise.f90` | Graph-based denoising built on the embedding |

## 2. ManifoldEM files reviewed

| File | Role |
|---|---|
| `DMembeddingII.py` | Core DM pipeline: kernel construction, `slaplacian`, `sembedding`, `fergusonE`, `op` |
| `core.py` | `fergusonE` bandwidth estimator; `L2_distance`; `svdRF` |
| `calc_distance.py` | CTF-corrected pairwise distance computation per projection direction |
| `psi_analysis.py` | NLSA (delay-embedded DM) and movie generation |
| `manifoldTrimmingAuto.py` | Iterative radius-ball manifold trimming |
| `embedd.py` | Re-embedding after interactive trimming |
| `find_conformational_coords.py` | Global conformational coordinate propagation across projection directions |

---

## 3. Feature comparison

### 3.1 SIMPLE capabilities not present in ManifoldEM

- **Projection-gated kNN graph** — avoids O(N²) all-pairs comparisons by binning
  particles by projection direction and searching only neighbouring angular bins
  (`simple_diff_map_graphs.f90`, `build_gated_euclidean_knn_graph`,
  `find_gated_euclidean_neighbors_rows`).
- **Automatic eigengap dimensionality selection** — `auto_ndiff_from_eigengap`
  in `simple_diffusion_maps.f90`.  ManifoldEM always takes a fixed
  `nEigs = params.num_eigs`.
- **Nyström out-of-sample extension** — the optional `nystrom_coords` branch in
  `embed_graph` in `simple_diffusion_maps.f90`.  ManifoldEM does not
  implement out-of-sample extension.

> **Note**: SO(2)/SE(2) steerable diffusion maps (`embed_so2_graph`,
> `embed_se2_graph`, `connection_matvec`, `coalesce_coo_to_csr`,
> `pack_knn_to_csr`, and associated graph payload fields `theta`, `shift_x`,
> `shift_y`) were removed as unused code.

---

### 3.2 ManifoldEM capabilities not present in SIMPLE

The following sections describe each capability with exact file and line
references to ManifoldEM and to the corresponding SIMPLE location where a change
would be applied.

---

#### 3.2.1 Ferguson / tanh data-adaptive bandwidth (HIGH PRIORITY)

**ManifoldEM location**: `core.py`, function `fergusonE` (~line 182);
called from `DMembeddingII.py`, function `op` (~line 576):

```python
popt, logSumWij, resnorm, R_squared = fergusonE(np.sqrt(yVal), logEps)
sigma = tune * np.sqrt(2 * np.exp(-popt[1] / popt[0]))
```

`fergusonE` scans a range of log-scale bandwidths
(`logEps = np.arange(-150, 150.2, 0.2)`) and for each ε computes:

$$\log \sum_{ij} W_{ij}(\varepsilon) = \log \sum_{ij} \exp\!\Bigl(-\frac{d_{ij}^2}{2\varepsilon}\Bigr)$$

It then fits the tanh model:

$$f(x) = d + c \tanh(a x + b)$$

and extracts the inflection point as:

$$\sigma^* = \mathtt{tune} \cdot \sqrt{2 \exp(-b/a)}$$

This is the Coifman–Lafon optimal bandwidth selector (see Coifman & Lafon,
*Applied and Computational Harmonic Analysis* 21, 2006, §3.3).

**SIMPLE current approach**: `simple_diff_map_graphs.f90`,
`pack_scalar_knn_to_csr` (~line 450) and `pack_knn_to_csr` (~line 381):

```fortran
eps = median_positive(kth_d2)
if( eps < DTINY ) eps = max(sum(kth_d2) / real(max(size(kth_d2),1)), 1.e-6)
...
w = exp(-max(d2s(m,i), 0.) / eps)
```

The bandwidth is the median of the k-th-nearest-neighbour squared distances.
This is a reasonable heuristic but has no analytical guarantee of being in
the correct scale for the Laplace-Beltrami operator.

**Suggested change**: add an optional `bandwidth_mode` flag to
`build_euclidean_knn_graph` / `build_orientation_knn_graph`.  When set to
`'ferguson'`, collect all non-zero edge distances into a flat array after the
kNN search (before the CSR pack) and run the tanh scan to estimate `eps`.
Keep the median default so no existing caller is affected.  The scan itself
can be implemented compactly: the summation over distances at each trial ε is
an inner-product loop and the tanh fit is a 4-parameter nonlinear least
squares (e.g. Levenberg-Marquardt via `simple_linalg` or a dedicated new
utility).

---

#### 3.2.2 Coifman–Lafon α normalization (HIGH PRIORITY)

**ManifoldEM location**: `DMembeddingII.py`, function `slaplacian` (~line 97):

```python
d = np.array(l.sum(axis=0)).T
if options.alpha != 1:   # apply non-isotropic normalization
    d = d**options.alpha
yVal = yVal / (d[yRow].flatten('C') * d[yCol].flatten('C'))
```

ManifoldEM uses `alpha=1` by default throughout `op` (~line 594):

```python
alpha = 1   # Laplace-Beltrami operator
```

The three interpretations exposed by `slaplacian`:
- `alpha = 0` — plain graph Laplacian (sampling density is retained in the
  operator; what SIMPLE currently implements via its symmetric degree
  normalization in `normalize_diffmap_graph`).
- `alpha = 0.5` — Fokker-Planck diffusion.
- `alpha = 1` — Laplace-Beltrami operator (sampling density is divided out;
  eigenvalues reflect intrinsic geometry only).

**SIMPLE current approach**: `simple_diff_map_graphs.f90`,
`normalize_diffmap_graph` (~line 591):

```fortran
self%wnorm(p) = self%w(p) / sqrt(max(deg(i), DTINY) * max(deg(j), DTINY))
```

This is a pure symmetric degree normalization — it corresponds to `alpha = 0`
(graph Laplacian) applied to already-Gaussianized weights.  The
sampling-density correction that `alpha = 1` would provide is absent.

**Implication for cryo-EM**: because orientation coverage in cryo-EM datasets is
never uniform (preferred orientations, angular gaps), the plain graph Laplacian
embeds sampling density artefacts into the leading eigenvectors.  With
`alpha = 1` those artefacts are cancelled, and the eigenvectors reflect the
intrinsic molecular geometry of the particle set.

**Suggested change**: add an optional `alpha` argument (default `0`, preserving
current behaviour) to `normalize_diffmap_graph`.  When `alpha > 0`, compute a
pre-normalization step before the symmetric degree division:

```
1. Compute raw degree d_i = sum_j w_ij
2. For alpha > 0: rescale w_ij <- w_ij / (d_i^alpha * d_j^alpha)
3. Recompute rescaled degree and apply symmetric 1/sqrt(d_i d_j) normalisation
```

This is a two-pass change confined entirely to `normalize_diffmap_graph` and
does not affect graph construction or the eigensolver.

---

#### 3.2.3 Stationary measure / Riemannian density (MEDIUM PRIORITY)

**ManifoldEM location**: `DMembeddingII.py`, function `op` (~line 611):

```python
mu = v[:, 0]
mu = mu * mu   # Riemannian measure; sum(mu)=1
```

After `sembedding`, the zeroth eigenvector (trivial eigenfunction,
eigenvalue ≈ 1) is squared to give the stationary probability measure on the
manifold.  ManifoldEM stores `mu` and uses it in NLSA to weight images:

```python
# psi_analysis.py, _NLSA, ~line 57:
mu_psi = mu.reshape((-1,1)) * psi
A[ind4:ind5,:] = np.matmul(tmp, mu_psi)
```

**SIMPLE current approach**: `embed_graph` in `simple_diffusion_maps.f90`
(~line 139) skips eigenvector 0 entirely:

```fortran
do k = 1,ndiff_used
    j = nev - k          ! nev is ndiff_scan+1; k=1 takes the 2nd-largest eigenvalue
    coords(k,i) = evals(j) * evecs(i,j)
end do
```

The trivial eigenfunction is discarded without ever being exposed to callers.

**Suggested change**: add an optional output `stationary_measure` of shape
`(n)` to `embed_graph`.  When present, fill it with `evecs(:, nev)**2` before
the coordinate loop.  No existing caller is changed.  The measure can then be
used for density-weighted trajectory generation or averaging, if needed.

---

#### 3.2.4 Iterative manifold trimming (MEDIUM PRIORITY)

**ManifoldEM location**: `manifoldTrimmingAuto.py`, function `op` (~line 64):

```python
posPath1 = get_psiPath(psi, rad, 0)
cc = 0
while len(posPath1) < nS:
    cc += 1
    nS = len(posPath1)
    D1 = D[posPath1][:, posPath1]
    lamb, psi, sigma, mu, ... = DMembeddingII.op(D1, k, tune, 600000)
    posPath1 = posPath1[posPathInt]
```

The loop re-embeds only the subset of particles whose embedding coordinates
fall inside a radius ball centred on the manifold origin, iterating until the
included set stabilises.  The radius `rad` is a user parameter.

Also called from `embedd.py` after interactive point trimming:
```python
lamb, psi, sigma, mu, logEps, logSumWij, popt, R_squared = \
    DMembeddingII.op(D1, k, params.nlsa_tune, 60000)
```

**SIMPLE current approach**: no manifold trimming is implemented.  Outlier
particles remain in the graph and appear at the fringes of the embedding.

**Suggested change**: add a standalone utility function
`trim_embedding_by_radius` in `simple_diffusion_maps.f90` that accepts an
existing graph plus an embedding radius parameter and returns a logical mask of
kept particles.  Callers that want trimming run the masking loop, rebuild the
graph on the kept subset, and re-embed.  Keeping this as a utility (rather
than auto-running inside `embed_graph`) preserves the existing interface and
lets callers decide whether the cost of re-embedding is justified.

---

#### 3.2.5 Nonlinear Laplacian Spectral Analysis / delay embedding (LOW PRIORITY, HIGH IMPACT)

**ManifoldEM location**: `psi_analysis.py`, internal function `_NLSA`
(~line 37):

```python
ConD = np.zeros((num - ConOrder, num - ConOrder))
for i in range(ConOrder):
    Ind = range(i, num - ConOrder + i)
    ConD += DD[Ind][:, Ind]
lamb, psi, sigma, mu, ... = DMembeddingII.op(ConD, k, tune, 600000)
```

`ConD` is the delay-embedded distance matrix: a sum of `ConOrder` shifted
sub-matrices of the pairwise distance matrix `DD`.  This is a discrete
analogue of Takens delay embedding.  The DM eigenvalues and eigenfunctions of
`ConD` are smoother along the conformational trajectory than those of the
single-shot matrix `DD`, because noise that is uncorrelated across frames
averages out.

The NLSA eigenfunctions `psiC` are used to reconstruct CTF-deconvolved image
averages along the conformational trajectory:

```python
A[ind4:ind5,:] = np.matmul(tmp, mu_psi)
U, S, V = svdRF(A)
VX = np.matmul(V.T, psiC.T)
```

(`svdRF` is in `core.py` ~line 75; it is an SVD via eigendecomposition of
`A^T A` or `A A^T`.)

**SIMPLE current approach**: no delay-embedding or NLSA is implemented.
SIMPLE's conformational analysis operates on a single-shot embedding.

**Suggested change**: this is a significant addition.  The prerequisite is a
sorted distance matrix for each particle subset (i.e., particles must carry a
1-D ordering along the conformational coordinate before the delay matrix can
be constructed).  For a first prototype, the ordering could come from SIMPLE's
existing 1-D open-manifold fit (`fit_1D_open_manifold_3D` is an analogue in
ManifoldEM at `fit_1D_open_manifold_3D.py`).  The delay summation itself is
cheap once the ordering is known.

---

#### 3.2.6 Where Noise Weighting Lives In SIMPLE vs ManifoldEM (RELEVANT TO FLEX WORKFLOWS)

**ManifoldEM location**: `calc_distance.py`, `get_distance_CTF_local`
(~line 218):

```python
img_f = fft2(img)
CTF_i = CTF1[ind3, :, :]
img_f_wiener = img_f * (CTF_i / wiener_dom[i, :, :])
img = ifft2(img_f_wiener).real
```

Each particle image is Wiener-deconvolved with its per-particle CTF before
computing squared-Euclidean distances between images in the same projection
bin. The Wiener denominator `wiener_dom` accumulates `sum_i CTF_i^2 + 1/SNR`.

**Correction and SIMPLE-specific interpretation**:

The primary SIMPLE noise weighting is not ad hoc per-particle Wiener filtering
at distance-computation time. In the ML-regularized reconstruction path
(`ml_reg=yes`), weighting is applied in Fourier observation/operator assembly and
in volumetric regularization, using sigma2 and SSNR-derived priors.

In particular, the Fourier-plane generator supports sigma2-weighted assembly:

- `src/main/image/simple_image_ctf.f90` (`gen_fplane4rec`):
  `CTF*y/sigma2` and `CTF^2/sigma2` for reconstruction planes.
- The same routine also supports observation-model assembly with
  `y/sqrt(sigma2)` and `CTF/sqrt(sigma2)` transfer planes.

So the scientific question for SIMPLE is not "should we add ManifoldEM-style
particle-level SNR filtering before distance calculations?" The more relevant
question is:

> How do we preserve the same ML-style volumetric weighting/prior behavior in
> the flex_analysis latent residual-model reconstruction path?

**Current flex_analysis behavior**:

- `src/main/flex/simple_flex_analysis_strategy.f90` sets `ml_reg='no'` in
  `apply_defaults`.
- `src/main/flex/simple_flex_projected_latent_model.f90`
  (`prep_imgs4projected_model`) calls `gen_fplane4rec(..., observation_model=.true.)`
  without sigma2 input, so no sigma2 division is applied there.
- `src/main/flex/simple_flex_diffmap_rec3D.f90` uses
  `prepare_unfiltered_model_params` (currently only forces `lp=0.`), i.e. the
  pre-image reconstruction is intentionally unfiltered/unregularized relative to
  the `ml_reg=yes` refine3D path.

**Implication**:

The present flex_analysis path is internally coherent, but it does not mirror
the full `ml_reg=yes` volumetric weighting contract. If parity with ML
regularization is desired, the extension point is the latent-model plane
preparation and coupled-basis solve path, not distance-time per-particle
Wiener filtering.

**Implementation draft (optional ML-reg parity for flex_analysis)**:

1. Add an explicit flex toggle and keep current behavior as default.
   - Add a parameter (for example `flex_ml_reg`) with default `no`.
   - Keep `apply_defaults` in `src/main/flex/simple_flex_analysis_strategy.f90`
     backward-compatible (`ml_reg='no'` unless `flex_ml_reg=yes`).
   - When `flex_ml_reg=yes`, propagate `ml_reg='yes'` and set `l_ml_reg=.true.`
     in the flex reconstruction parameter copies.

2. Thread sigma2 into projected-model plane preparation.
   - Update `prep_imgs4projected_model` in
     `src/main/flex/simple_flex_projected_latent_model.f90` to optionally pass
     `sig2arr` into `gen_fplane4rec(..., observation_model=.true.)`.
   - Preserve the current code path when sigma2 is unavailable.
   - Keep the same Fourier contract (observation plane and transfer plane)
     while switching from unweighted to sigma2-weighted assembly.

3. Define a strict data-availability contract.
   - If `flex_ml_reg=yes` but required sigma2/noise-model inputs are missing,
     fail early with a clear error.
   - If `flex_ml_reg=no`, never require sigma2 (current behavior).

4. Keep regularization ownership in the reconstruction stage.
   - Do not introduce distance-time particle Wiener filtering in the DM metric.
   - Apply weighting where SIMPLE already expects it: Fourier-plane assembly,
     coupled basis accumulation, and volumetric solve/finalization.

5. Keep low-pass policy explicit.
   - `prepare_unfiltered_model_params` currently enforces `lp=0.` for model
     reconstruction in `src/main/flex/simple_flex_diffmap_rec3D.f90`.
   - If ML parity mode is enabled, decide explicitly whether this remains
     unfiltered for latent-basis estimation or follows refine3D-style policy.

6. Add targeted regression checks before broad rollout.
   - Extend `test_projected_model_plane_preparation` to run both modes:
     `flex_ml_reg=no` (legacy) and `flex_ml_reg=yes` (sigma2-weighted).
   - Keep algebraic identity checks in
     `test_fake_preimage_against_reconstruct3D` for both modes, with separate
     tolerances if needed.
   - Add one end-to-end flex fixture demonstrating that enabling
     `flex_ml_reg=yes` changes weighting behavior but preserves numerical
     stability and convergence.

7. Stage deployment to minimize risk.
   - Phase 1: plumb flags and sigma2 inputs, no default behavior change.
   - Phase 2: enable in controlled tests/benchmarks.
   - Phase 3: consider default changes only after quality and runtime review.

This path preserves SIMPLE's existing theory (ML-style weighting at
reconstruction time) while making flex_analysis capable of matching that
contract when explicitly requested.

---

## 4. Implementation priority summary

| Item | SIMPLE change location | Effort | Expected impact |
|---|---|---|---|
| α=1 Coifman-Lafon normalization | `normalize_diffmap_graph` in `simple_diff_map_graphs.f90` L591 | Low (one two-pass normalisation) | High — density-invariant geometry |
| Ferguson tanh bandwidth | `pack_scalar_knn_to_csr` in `simple_diff_map_graphs.f90` ~L380 | Medium (tanh scan + nonlinear fit) | Medium — principled kernel scale |
| Stationary measure output | `embed_graph` in `simple_diffusion_maps.f90` L139 | Low (optional output argument) | Low–medium |
| Iterative manifold trimming | New utility in `simple_diffusion_maps.f90` | Medium | Medium |
| NLSA / delay embedding | New section in `simple_diffusion_maps.f90` | High | High (if per-PrD movies wanted) |

The α normalization and Ferguson bandwidth are the lowest-risk and most
self-contained changes.  Both are confined to `simple_diff_map_graphs.f90` and
do not touch the eigensolver, the steerable path, or any caller.

---

## 5. References

- Coifman & Lafon, *Applied and Computational Harmonic Analysis* **21**, 5–30 (2006).
  Original DM paper; establishes the α normalisation and the tanh bandwidth selector.
- Ferguson, Kollias & Bharat, *Physical Review Letters* **104**, 118103 (2010).
  Source of the tanh/`fergusonE` approach used in ManifoldEM.
- Giannakis & Majda, *PNAS* **109**, 2222–2227 (2012).
  Introduces NLSA (the delay-embedding extension).
- Dashti et al., *Structure* **28**, 1120–1130 (2020).
  ManifoldEM application to cryo-EM conformational analysis.
