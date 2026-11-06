# Dose-series accumulation under a temporal attenuation model

Initial design, 3 November 2026, against SIMPLE `69c6c9d78` plus the working-tree
`refine_motion_model` changes. Producer integration is implemented and committed (section 10). The consumer,
`dose_filter3D`, is implemented in the working tree: its implementation and corrections are in section 11. The re-extraction of particles under the measured curve
(`reextract dose_weighting=measured`) is in section 12, its uncertainty and comparison corrections in 12.1, and
the validation to redo after these changes in section 13. This
note supersedes the `dose_filter3D_implementation_note.md` draft and incorporates the contract decisions of
section 9.

Working name: `dose_filter3D`. The voltage, box and doses come from the manifest, not the command line.

```text
simple_exec prg=dose_filter3D filetab=dose_series.txt smpd=<Å> mskdiam=<Å> \
    [hp=<Å>] [lp=<Å>] [bfac_model=line|free] [zerodose=yes|no] [outvol=dose_series.mrc] nthr=16
```

## 1 Scientific purpose and aim

A cryo-EM movie's later frames image a progressively more damaged specimen. The standard remedy is the
Grant & Grigorieff critical-dose curve — in SIMPLE `calc_dose_weights`
(`src/utils/filter/simple_estimate_ssnr.f90`), applied frame by frame in motion correction — which weights a
frame's amplitudes by $\exp\!\big(-D / 2N_e(k)\big)$ with

$$
N_e(k) = 0.245\,k^{-1.665} + 2.81 \quad \mathrm{e}/\text{Å}^2 .
\tag{1}
$$

One curve, with three constants, for every specimen ever imaged. It is a prior, not a measurement: the specimen,
its support and its ice decide how fast the structure actually fades.

**Aim.** Measure that decay on the data itself, from a series of per-dose-window reconstructions, and use the
measured curve in place of the generic constants — first to accumulate the windows optimally into one map, and
later to weight frames in motion correction.

Three products, in increasing order of what they commit to:

1. **A measurement with no model** — the relative attenuation of every window at every resolution shell, written
   before any premise is applied (section 4.3). A Guinier plot per window.
2. **A model** — one B-factor per window, those B-factors on a smooth curve of dose (sections 4.4-4.5), and the
   accumulation that curve implies (section 4.6).
3. **A comparison** — the model against the Grant & Grigorieff curve in the same estimator, as a cross-half
   prediction error band by band (section 5). It is descriptive: the attenuations are fitted on both halves, so
   it shows where the premise fails on these data but is not an independent test of it.

### 1.1 Relation to `uniform_filter3D`

The closest kin in SIMPLE is `uniform_filter3D` (`exec_uniform_filter3D` in `simple_commanders_resolest.f90`,
`estimate_lplim` in `src/utils/filter/simple_opt_filter.f90`): it filters the odd half with a Butterworth
low-pass at every candidate cutoff and keeps the cutoff whose filtered odd half lies closest to the even half
inside the mask, the squared distance being the cross-validation cost. `estimate_lplim3D` runs the same search
on refine3D's reference halves to set the low-pass limit; `nu_filt3D` runs it per voxel over a bank of cutoffs
with Potts smoothing of the labels.

The dose series is that search with $W$ half-map pairs instead of one, and the attenuation model in place of the
Butterworth cutoff as the filter family: the prediction of a window's half from the other halves stands where
the filtered odd half stood, the even/odd squared distance is the cost, and a window predicted from its own even
half alone recovers `uniform_filter3D`'s cost with the all-pass filter. Per shell the cost is a quadratic form in
the Gram matrix of the $2W$ halves (section 4.2), so a candidate costs a few numbers where `uniform_filter3D`
passes over the volume once per candidate.

One difference of purpose: `uniform_filter3D` delivers the *filtered* half (per shell its cost is minimised by
the FSC itself, the Wiener filter of one half); the dose series delivers the *unfiltered* accumulation and leaves
the filtering to `postprocess` or `enhance`.

## 2 Notation

Every symbol used in this note. Where a symbol collides with a conventional one, the clash is named.

### 2.1 Indices and geometry

| Symbol | Meaning |
|---|---|
| $w, v$ | Window index, $1 \le w,v \le W$; windows are ordered by increasing dose |
| $W$ | Number of dose windows in the series |
| $h$ | Half index, $h \in \{\mathrm{e}, \mathrm{o}\}$ for even and odd |
| $t$ | The *target* window in a cross-validation fold (section 5) |
| $\mathbf q$ | Integer Fourier coordinate on the half-volume lattice: $q_1 \ge 0$, and on the $q_1 = 0$ plane only $q_2 > 0$, or $q_2 = 0$ and $q_3 > 0$; the mates are conjugates |
| $N$ | Box width in pixels (**not** a noise symbol; noise is $\varepsilon$) |
| $d$ | Pixel size in Å |
| $s$ | Radial shell index, $s(\mathbf q) = \mathrm{nint}(\lVert\mathbf q\rVert)$, running $1 \le s \le N/2$ |
| $k_s$ | Spatial frequency of shell $s$, $k_s = s/(Nd)$ in Å$^{-1}$ |
| $n_s$ | Number of listed coefficients in shell $s$ |
| $n_{\mathrm{sh}}$ | Number of shells in the fitted band |
| $M(\mathbf x)$ | Real-space soft mask at the `mskdiam` radius |
| $f$ | Fraction of a shell's coefficients the mask leaves independent (section 4.2) |
| $n_s^{\mathrm{eff}}$ | Independent coefficient count of shell $s$, $n_s^{\mathrm{eff}} = f\,n_s$ |

### 2.2 Data and model

| Symbol | Meaning |
|---|---|
| $E_w(\mathbf q),\ O_w(\mathbf q)$ | Fourier coefficients of window $w$'s even and odd half maps |
| $Y_w^h(\mathbf q)$ | Either half of window $w$, when the distinction does not matter |
| $S(\mathbf q)$ | The structure, as window 1 sees it |
| $\varepsilon_w^h(\mathbf q)$ | Noise in half $h$ of window $w$ |
| $\sigma_w^2(s)$ | Variance of one half's noise coefficient in shell $s$ |
| $p_w(s)$ | Precision, $p_w(s) = 1/\sigma_w^2(s)$ |
| $P(s)$ | Signal power of window 1 in shell $s$, $P(s) = n_s^{-1}\sum_{\mathbf q\in s}\lvert S(\mathbf q)\rvert^2$ |
| $a_w(s)$ | Attenuation of window $w$ at shell $s$, relative to window 1 |
| $c_w$ | Scale of window $w$ (dimensionless); $c_1 \equiv 1$ |
| $B_w$ | B-factor of window $w$ in Å$^2$, relative to window 1; $B_1 \equiv 0$ |
| $D_w$ | Accumulated dose at window $w$'s centre, e/Å$^2$, in the scaled coordinate of section 3.3 |
| $B(D)$ | The dose curve, $B(D) = B_0 + \beta D$ |
| $B_0$ | Intercept of the dose curve, Å$^2$ |
| $\beta$ | Slope of the dose curve, Å$^2$ per e/Å$^2$ — **the headline measurement** |
| $\delta_w$ | Window $w$'s excess over the curve, $B_w = B(D_w) + \delta_w$ |
| $\tau$ | Scatter of $\delta_w$ across the windows, Å$^2$ |
| $B_w^\star$ | The B-factor actually used for window $w$, after shrinkage |
| $a_w^\star(s)$ | Attenuation from the shrunk parameters, $a_w^\star = c_w \exp(-B_w^\star k_s^2/4)$ |
| $N_e(k)$ | Critical exposure of equation (1), e/Å$^2$ (**not** the box width $N$) |

### 2.3 Statistics and output

| Symbol | Meaning |
|---|---|
| $\mathbf y(\mathbf q)$ | The $2W$-vector $(E_1,\dots,E_W,\,O_1,\dots,O_W)$ at coefficient $\mathbf q$ |
| $G(s)$ | $2W \times 2W$ real symmetric Gram matrix of shell $s$ |
| $G^{\mathrm{ee}}, G^{\mathrm{oo}}$ | Its within-half $W \times W$ diagonal blocks |
| $X_{wv}(s)$ | Its cross-half block, $X_{wv} = n_s^{-1}\sum_{\mathbf q\in s}\mathrm{Re}\big[E_w\overline{O_v}\big]$ |
| $C(s)$ | The symmetrised cross-half block, $C = \tfrac12(X + X^{\mathsf T})$ |
| $V_{wv}(s)$ | Variance of the estimate $C_{wv}(s)$ |
| $\mathrm{se}_w(s)$ | Standard error of $a_w(s)$ from section 4.3 |
| $\mathrm{se}_w$ | Standard error of $B_w$ from section 4.4 |
| $z_w$ | Residual of window $w$ against the dose curve, in units of its total error |
| $\omega_w(s)$ | Accumulation weight of window $w$ at shell $s$; $\sum_w \omega_w(s) = 1$ |
| $A_{\mathrm e}, A_{\mathrm o}, A$ | Accumulated even half, odd half, and their mean |
| $\bar a(s)$ | Attenuation of the accumulated map relative to window 1 |
| $n_{\mathrm{eff}}(s)$ | Kish's effective number of windows contributing to shell $s$ |
| $S_0(\mathbf q)$ | The zero-dose extrapolation |
| $\hat S_{\mathrm e}(\mathbf q)$ | Signal estimate built from the even halves alone (section 5) |
| $e_t(s)$ | Whitened cross-half prediction error for target window $t$ |
| $\overline{(\cdot)}$ | Complex conjugate |
| $(\cdot)^{\mathsf H}$ | Conjugate transpose |

A bar over a vector symbol is never used; $\bar a$ is a scalar per shell. $\sigma^2$ means sigma squared, not
sigma raised to a second power again.

## 3 The premise

### 3.1 The model

Between windows the structure does not change, it fades. Window $w$ sees the structure through one scale and one
B-factor:

$$
Y_w^h(\mathbf q) = a_w(s)\,S(\mathbf q) + \varepsilon_w^h(\mathbf q),
\qquad
a_w(s) = c_w \exp\!\left(-\frac{B_w\,k_s^2}{4}\right),
\qquad h \in \{\mathrm e, \mathrm o\}
\tag{2}
$$

and those B-factors lie on a smooth curve of the accumulated dose:

$$
B_w = B(D_w) + \delta_w,
\qquad
B(D) = B_0 + \beta D,
\qquad
\delta_w \sim \big(0,\ \tau^2\big).
\tag{3}
$$

The B-factor convention is the usual one on amplitudes, $\exp(-Bk^2/4)$ at $k = 1/\text{resolution}$. The
structure $S$ absorbs an overall scale and B-factor, so the gauge is

$$
c_1 = 1, \qquad B_1 = 0 :
\tag{4}
$$

$S$ is the structure *as window 1 sees it*, every $a_w$, $c_w$ and $B_w$ is relative to window 1, and the dose
curve extrapolates below it.

**Noise.** $\mathbb E[\varepsilon] = 0$ and $\mathrm{Var}(\varepsilon_w^h(\mathbf q)) = \sigma_w^2(s)$. Detector
noise is independent across windows, halves and coefficients: different windows are different movie frames, so
their shot noise is independent. What *is* shared between windows, within a half, is the structure and the
alignment error — the same particles at the same poses. That asymmetry is the whole basis of section 4.2.

### 3.2 Why this model

- **Why a B-factor.** A Gaussian in $k$ is the one-parameter description of a resolution-dependent amplitude
  loss, and the quantity `postprocess` and `enhance` already speak in.
- **Why a line in dose.** If damage accumulates smoothly, $B$ should be a smooth function of $D$. A line is the
  simplest thing the data can refute, and the residual scatter $\tau$ measures the refutation.
- **Why bother.** The model has $2W + 2$ parameters fitted from of order $10^8$ Fourier coefficients. Its weights
  are smooth in resolution by construction, where a free per-shell weight set would fit noise.
- **Closed form.** Under the model the optimal accumulation is the matched filter (section 4.6): no target
  window, no weight to choose per shell.
- **It extrapolates.** The structure at zero dose is a per-shell rescaling of the accumulated map.
- **It transfers.** The dose curve gives per-frame weights for motion correction, replacing the generic
  constants — a reconstruction from all frames at full time resolution with measured weights.

### 3.3 Orientation against Grant & Grigorieff

Their amplitude factor $\exp\!\big(-D/2N_e(k)\big)$ is equation (2) with $B = \beta D$ when

$$
\frac{D}{2N_e(k)} = \frac{\beta D\,k^2}{4}
\qquad\Longleftrightarrow\qquad
\beta = \frac{2}{N_e(k)\,k^2} .
\tag{5}
$$

Evaluating (5) at the constants of (1) gives $\beta = 2.2$ Å$^2$ per e/Å$^2$ at 2 Å, $4.2$ at 3 Å and $14$ at
10 Å. **Their curve is therefore not one B-factor**: it fades the low shells faster than a B-factor extrapolated
from the high shells would. Which description fits a given specimen is what sections 4.4 and 5 measure.

**Dose coordinate.** `calc_dose_weights` divides the accumulated dose by $1$ at 300 kV, $0.8$ at 200 kV and
$0.64$ at 100 kV — it *divides*, raising the effective dose at lower voltage. $D_w$ in this note is the dose in
that same scaled coordinate, so the fitted $\beta$ drops into `calc_dose_weights` unchanged and is directly
comparable to the equivalents of (5). The price is that $\beta$ is then not a physical specimen constant and is
not comparable across voltages; that is recorded as a limitation (section 8). At 300 kV the two coordinates
coincide. The input manifest records *physical* dose; the program applies the scaling.

## 4 Method

### 4.1 Inputs, window selection, noise and band

**Inputs.** $W$ raw half pairs from `refine_motion_model`'s interval reconstructions —
`recvol_state01_{even,odd}_unfil.mrc` in each `rec_frames_aa-bb/`, unfiltered and unmasked, state 1 only — at one
box $N$, one sampling $d$, one pose convention, with the same even/odd particle partition in every window. A
manifest written by the producer carries, per window: the two half-map paths, the inclusive fraction range
$[a,b]$, the physical dose bounds derived from the saved motion model's `dose_per_frame`, and a fingerprint of
the `eo` partition.

#### 4.1.1 The manifest

`filetab` names the manifest. It is not a bare `read_filetable` list: that helper reads every line
verbatim with no comment handling, so `dose_filter3D` parses the file itself. Comment lines start with
`#`; the first line is the version and is required. Half-map paths are relative to the manifest's own
directory, so a run directory can be moved. Doses are **physical** e/Å$^2$; the consumer applies the
voltage scaling of section 3.3.

```text
# SIMPLE dose series 1
# program            refine_motion_model
# state              1
# box                288
# smpd               1.06
# kv                 300.0
# nptcls             12043
# neven              6021
# nodd               6022
# eo_digest          6C4DD423214CC38D
# dose_per_fraction  1.2500000000E+00
# columns  from to dose_lo dose_hi even odd
 1  3  0.0000000000E+00  3.7500000000E+00  rec_frames_01-03/recvol_state01_even_unfil.mrc  rec_frames_01-03/recvol_state01_odd_unfil.mrc
 4  6  3.7500000000E+00  7.5000000000E+00  rec_frames_04-06/recvol_state01_even_unfil.mrc  rec_frames_04-06/recvol_state01_odd_unfil.mrc
```

`dose_lo` is the accumulated dose at the **start** of fraction `from`, `dose_hi` at the **end** of
fraction `to`; the producer takes both from the saved motion model's `dose_per_frame` and the original
fraction IDs. The window centre of section 3.3 is then
$D_w = \tfrac12(\texttt{dose\_lo} + \texttt{dose\_hi})$ after voltage scaling, and
$\Delta D_w = \texttt{dose\_hi} - \texttt{dose\_lo}$ is what section 8's deferred window-integrated
attenuation would need. `eo_digest` is a hash of the selected particles' `eo` vector, constant across
windows by construction; it is recorded once and only lets the consumer confirm that every window came
from one partition.

The consumer rejects: a missing or unknown version line; fewer than two windows after the pruning
below; non-increasing `from`; a `dose_hi` not greater than its `dose_lo`; non-increasing dose across
windows; a `kv` outside {300, 200, 100}; a missing half map; and any volume whose box, sampling or
cubic shape disagrees with the header.

**Window selection.** `refine_motion_model`'s last interval keeps its full `stepf` width, so when
$(\texttt{tof} - \texttt{fromf} + 1)$ is not a multiple of `stepf` it overlaps its predecessor and the two share
particle images — which would break the independence of (2). Drop the *predecessor*, keeping the full-width
highest-dose window. The survivors are disjoint but need not tile the exposure. Fewer than two surviving windows
is a hard error: below $W = 2$ the model is not identifiable.

**Lattice.** For the statistics each half is read, soft-masked by $M$ at the `mskdiam` radius about its own
solvent level, Fourier transformed, and the listed coefficients gathered into an $(n_{\mathrm{coeff}}, W)$
complex single-precision array per half. Removing the global mean instead would leave an offset under the soft
edge, which rings into every shell. For the accumulation of section 4.6 the halves are read again *as given*,
with only their mean removed; the two sets are never held at once.

**Noise.** Per window and shell, from the half difference:

$$
\hat\sigma_w^2(s) = \frac{1}{2\,n_s}\sum_{\mathbf q\in s}\big\lvert E_w(\mathbf q) - O_w(\mathbf q)\big\rvert^2,
\qquad
p_w(s) = \frac{1}{\hat\sigma_w^2(s)} .
\tag{6}
$$

> **Why (6) measures noise and nothing else.** Both halves see the same structure under the same attenuation, so
> the signal cancels exactly:
> $$E_w - O_w = \big(a_w S + \varepsilon_w^{\mathrm e}\big) - \big(a_w S + \varepsilon_w^{\mathrm o}\big) = \varepsilon_w^{\mathrm e} - \varepsilon_w^{\mathrm o}.$$
> With independent, equal-variance half noise $\mathrm{Var}(\varepsilon^{\mathrm e} - \varepsilon^{\mathrm o}) = 2\sigma_w^2$, hence the factor 2 in the denominator. This assumes the two halves have equal noise, which holds to the extent the even/odd particle split is balanced.
>
> **One spectrum suffices, not two.** The statistics run on masked halves and the accumulation on the halves as
> given, and those have different noise levels. But a soft mask multiplies in real space, so it convolves in
> Fourier, and where the noise power is locally smooth across the kernel width it scales *every* window's noise
> power by the same scalar $\langle M^2\rangle$. The weights of section 4.6 are normalised to sum one, so a
> factor common to all $p_w$ cancels: $\omega$ is invariant under $p_w \mapsto \gamma(s)\,p_w$. The masked
> spectrum therefore serves both roles — absolutely for the whitening of section 5, which is computed on the
> same masked quantities, and up to the cancelling factor for the weights. The approximation degrades only at
> the lowest few shells, where the noise power varies fast across the kernel.

**Band.** From the high-pass shell $s(\texttt{hp})$, `hp` defaulting to 20 Å, to $s(\texttt{lp})$ when `lp` is
given; otherwise up to the last shell where *any* window's own FSC,
$C_{ww}/\sqrt{G^{\mathrm{ee}}_{ww}G^{\mathrm{oo}}_{ww}}$, stays above $0.1$ — gaps allowed, at least six shells.
B-factors carry no information at the lowest shells, where a handful of coefficients is dominated by the mask
edge and the solvent: on a synthetic series a band starting at shell 1 let the 96 Å shell, 13 coefficients,
drive every $B_w$ negative.

> **Why "any window" and not the sum.** The early windows reach further in resolution than the full-dose sum
> does, because the sum's late windows add only noise there. Taking the band from the sum would cut it short
> exactly where this method has the most to gain.

### 4.2 The Gram matrices: sufficient statistics

One pass over the coefficients. With $\mathbf y(\mathbf q) = (E_1,\dots,E_W, O_1,\dots,O_W)$,

$$
G(s) = \frac{1}{n_s}\sum_{\mathbf q\in s}\mathrm{Re}\big[\mathbf y(\mathbf q)\,\mathbf y(\mathbf q)^{\mathsf H}\big]
\qquad (2W \times 2W,\ \text{real symmetric}),
\tag{7}
$$

with within-half blocks $G^{\mathrm{ee}}$, $G^{\mathrm{oo}}$ and cross-half block

$$
X_{wv}(s) = \frac{1}{n_s}\sum_{\mathbf q\in s}\mathrm{Re}\big[E_w(\mathbf q)\,\overline{O_v(\mathbf q)}\big],
\qquad
C = \tfrac12\big(X + X^{\mathsf T}\big),
\qquad
\mathbb E\big[C_{wv}(s)\big] = a_w(s)\,a_v(s)\,P(s).
\tag{8}
$$

> **Why the model is fitted on the cross-half block alone.** $E_w$ and $O_v$ are built from disjoint particle
> halves, so their noise and their alignment errors are independent and $\mathbb E[E_w\overline{O_v}]$ is
> unbiased for the signal product $a_w a_v P$. By contrast $\mathbb E[E_w\overline{E_v}]$ carries the half's
> shared error on top and is biased upward. Note that this also makes the *diagonal* $C_{ww}$ usable — it is the
> window's own FSC numerator — unlike the diagonal of an ordinary Gram matrix, which is pure power.

The variance $V_{wv}(s)$ is the noise of $C_{wv}(s)$ alone, given the structure. Per coefficient,
$\mathrm{Var}\,\mathrm{Re}[E_w\overline{O_v}] = \tfrac12\big(\lvert a_wS\rvert^2\sigma_v^2 + \lvert a_vS\rvert^2\sigma_w^2
+ \sigma_w^2\sigma_v^2\big)$, and with $a_w^2P = G^{\mathrm{ee}}_{ww} - \sigma_w^2$ every term is observable:

$$
\mathrm{Var}(X_{wv}) = \frac{G^{\mathrm{ee}}_{ww}\,\sigma_v^2 + G^{\mathrm{oo}}_{vv}\,\sigma_w^2 - \sigma_w^2\sigma_v^2}{2\,n_s^{\mathrm{eff}}},
\qquad
V_{wv} = \tfrac14\big(\mathrm{Var}(X_{wv}) + \mathrm{Var}(X_{vw})\big),\quad V_{ww} = \mathrm{Var}(X_{ww}).
\tag{8a}
$$

> **Why not the spread of the products.** The sample variance of $\mathrm{Re}[E_w\overline{O_v}]$ over a shell also
> holds $a_wa_v(\lvert S_{\mathbf q}\rvert^2 - P)$, the structure's own power fluctuating from coefficient to
> coefficient. That is no error: $P(s)$ is fitted per shell as the shell's realized power, and the fluctuation is
> common to every pair, so it cancels in the attenuations. Counting it inflates every error by about the amplitude
> SNR — 35-fold at high SNR on synthetic series. The off-diagonal $V_{wv}$ is a quarter of the summed variances
> because $C$ averages $X_{wv}$ and $X_{vw}$, two estimates with independent noise.

> **Why this independent count.** The mask correlates neighbouring Fourier coefficients, so a shell holds fewer
> independent samples than it holds coefficients. Masked noise is correlated by a kernel $\rho$, the transform of
> $M^2$, which is three-dimensional; a shell, one coefficient thick, cuts only a plane of it. By the projection-slice
> theorem that plane's $\sum\lvert\rho\rvert^2$ is $N^2\sum_{xy}\big(\sum_z M^2\big)^2 / \big(\sum M^2\big)^2$, so
> $$f = \frac{3\,\big(\sum_{\mathbf x} M^2\big)^2}{N^2\sum_{\text{axes}}\sum_{\text{plane}}\big(\sum_{\perp} M^2\big)^2},$$
> averaged over the three axes; for a ball of volume fraction $f_V$ it is $1.075\,f_V^{2/3}$, not $f_V$. The
> volume fraction is the count for a volume, not a shell: with SIMPLE's soft mask at `mskdiam` 66 Å in a 96 Å box,
> $\langle M^2\rangle = 0.142$ where a simulation of the masked noise products gives 0.314 and the formula 0.311.

> **Why one pass is enough.** Everything downstream is a quadratic form in $\mathbf y$. The cross-half prediction error of any
> weight set $\boldsymbol\omega$ predicting window $t$'s odd half from the even halves is
> $$\frac{1}{n_s}\sum_{\mathbf q\in s}\Big\lvert O_t(\mathbf q) - \sum_v \omega_v E_v(\mathbf q)\Big\rvert^2 = G^{\mathrm{oo}}_{tt} - 2\sum_v \omega_v X_{vt} + \boldsymbol\omega^{\mathsf T} G^{\mathrm{ee}} \boldsymbol\omega ,$$
> and the FSC of any accumulation is
> $$\mathrm{FSC}(s) = \frac{\boldsymbol\omega^{\mathsf T} C\, \boldsymbol\omega}{\sqrt{\big(\boldsymbol\omega^{\mathsf T} G^{\mathrm{ee}}\boldsymbol\omega\big)\big(\boldsymbol\omega^{\mathsf T} G^{\mathrm{oo}}\boldsymbol\omega\big)}} .$$
> So a candidate costs $O(W^2)$ arithmetic on $n_{\mathrm{sh}}$ small matrices, where `uniform_filter3D` passes
> over the volume once per candidate. Only section 4.6 touches the coefficients again.

### 4.3 Relative attenuation per shell, with no model

Per shell, the rank-one weighted least-squares fit of (8) with the gauge $a_1(s) = 1$:

$$
\min_{\,a_2 \dots a_W,\ P}\ \sum_{w \le v} \frac{\big[C_{wv}(s) - P(s)\,a_w(s)\,a_v(s)\big]^2}{V_{wv}(s)} .
\tag{9}
$$

Start from the leading eigenvector of the symmetric $C(s)$ scaled to $a_1 = 1$, then three alternating passes
(each $a_w$ given the others and $P$, then $P$ given all $a_w$). Output $a_w(s) \pm \mathrm{se}_w(s)$, the
standard errors from the sandwich covariance $H^{-1}\big(J^{\mathsf T}D\,\Sigma\,DJ\big)H^{-1}$ of (9), with $J$ the
Jacobian of $Pa_wa_v$, $D = \mathrm{diag}(1/V)$, $H = J^{\mathsf T}DJ$ and, to first order,

$$
C_{wv} - P a_w a_v = a_w\gamma_v + a_v\gamma_w + \varrho_{wv},
\qquad
\mathrm{Var}(\gamma_k) = \frac{P\sigma_k^2}{4\,n_s^{\mathrm{eff}}},
\qquad
\mathrm{Var}(\varrho_{wv}) = \frac{\sigma_w^2\sigma_v^2}{4\,n_s^{\mathrm{eff}}}\ \Big(\frac{\sigma_w^4}{2\,n_s^{\mathrm{eff}}}\ \text{for}\ w=v\Big).
\tag{9a}
$$

> **Why a sandwich.** $\gamma_k$, the projection of window $k$'s half noise onto the structure, enters every pair
> that holds $k$, so it does not average down over those pairs as the weights alone would assume: their curvature
> understates the error by about $\sqrt{W+1}$ at high SNR. $\Sigma$ is fixed by the quantities already estimated,
> and its diagonal is (8a).

> **What this is.** The Guinier diagnostic, written before any premise is applied. It asks only whether the
> windows differ by *one* multiplicative factor per shell — whether $C(s)$ is rank one. If they differ in more
> than one direction, which is structural change rather than fading, the rank-one fit will not describe $C$ and
> the residual shows it here, rather than being absorbed silently into a B-factor downstream.

### 4.4 One B-factor per window

Per window, over the band:

$$
\min_{c_w,\,B_w}\ \sum_s \frac{\Big[a_w(s) - c_w\exp\!\big(-B_w k_s^2/4\big)\Big]^2}{\mathrm{se}_w(s)^2} .
\tag{10}
$$

Start from the weighted regression of $\log a_w(s)$ on $k_s^2$ over the shells where $a_w(s)$ exceeds twice its
error, then refine (10) **in the linear domain** by Levenberg-Marquardt, keeping a step only when it lowers (10),
with standard errors $\mathrm{se}_w$ from the curvature, inflated by the mask's correlation across shells
(section 11.8). A window whose fit fails — no curvature, a non-finite objective, or a scale $c_w$ within five of
its own errors of zero — is **rejected** with a warning: it is left out
of the dose line and takes no weight in the model's accumulation (13), while the control and the flat mean keep
it. A window with no signal must take no weight, which the line's prediction at its dose would give it.

> **Why step acceptance.** Undamped Gauss-Newton steps diverged on a window carrying no signal: the objective rose
> from 10.67 to $3.5\times10^{45}$ with $B_w = -1927$, and under `bfac_model=free` that window took 99% of the
> weight at the band edge while the run reported normal completion. A window that is genuinely sharper than
> window 1 fits well and keeps its weight; the gate is on the fit, not on how much a window amplifies.

> **Why refine in $a$ and not in $\log a$.** Shells where the window has faded to noise have $a_w(s)$ near zero
> or negative. The logarithm is undefined or wildly weighted there. Fitting the residual in $a$ lets those
> shells contribute with the error they actually have.

Readout per window: the slope $-B_w/4$, the scale $c_w$, and the rms deviation of the shells from the line in
units of their errors. **A straight line in $k^2$ is the premise holding for that window.** A curvature
consistent across windows is the exponent question of section 8.

### 4.5 The dose curve

The weighted line $B(D) = B_0 + \beta D$ through the points $\big(D_w,\,B_w,\,\mathrm{se}_w\big)$ for
$w = 2 \dots W$, with the excess scatter $\tau^2$ by the usual random-effects (DerSimonian-Laird) estimate, its trace that of a weighted line (section 11.3) — zero when the windows
sit on the line within their errors — and per window the residual

$$
z_w = \frac{B_w - B(D_w)}{\sqrt{\mathrm{se}_w^2 + \tau^2}} .
\tag{11}
$$

> **Why window 1 is excluded from the fit.** $B_1 = 0$ is exact by the gauge (4), not measured: its standard
> error is identically zero, and both the weighted line and the shrinkage below divide by $\mathrm{se}^2$. It
> carries no information beyond fixing the origin. Fit on the $W-1$ real measurements and report the *predicted*
> $B(D_1)$, which should come out near zero, as the diagnostic in place of $z_1$.

With `bfac_model=line` the B-factor used for window $w$ is the shrinkage estimate

$$
B_w^\star = \frac{B_w/\mathrm{se}_w^2 + B(D_w)/\tau^2}{1/\mathrm{se}_w^2 + 1/\tau^2},
\qquad\text{with the limit}\quad \tau \to 0 \ \Longrightarrow\ B_w^\star = B(D_w).
\tag{12}
$$

> **Why shrink.** (12) is the posterior mean under $B_w \sim \mathcal N\big(B(D_w), \tau^2\big)$. It gives the
> line when the windows lie on it, and the window's own value when they do not: a motion-blurred first window
> keeps its own blur and shows up as a residual $z_1$, while the line is not pulled by it so long as $\tau$ is
> small. The $\tau\to0$ limit must be written explicitly — (12) is $0/0$ there, and that is precisely the case
> the premise predicts.

With `bfac_model=free`, $B_w^\star = B_w$. The scale is used as fitted, $c_w^\star = c_w$, since it is set by the
low shells where the signal is strong; its trend with dose is reported and should be flat.

### 4.6 The accumulation

With $a_w^\star(s) = c_w^\star \exp\!\big(-B_w^\star k_s^2/4\big)$, the weights and the accumulated halves are

$$
\omega_w(s) = \frac{a_w^\star(s)\,p_w(s)}{\sum_v a_v^\star(s)\,p_v(s)},
\qquad
A_{\mathrm e}(\mathbf q) = \sum_w \omega_w(s)\,E_w(\mathbf q),
\qquad
A = \tfrac12\big(A_{\mathrm e} + A_{\mathrm o}\big),
\tag{13}
$$

with $A_{\mathrm o}$ formed from the odd halves in the same way.

> **Why these weights.** Minimising $\sum_w \lvert Y_w - a_w S\rvert^2/\sigma_w^2$ over $S$ gives the matched
> filter under the model,
> $$\hat S = \frac{\sum_w a_w\,Y_w/\sigma_w^2}{\sum_w a_w^2/\sigma_w^2} ,$$
> and writing the delivered map as $A = \bar a\,\hat S$ with $\boldsymbol\omega$ normalised to sum one recovers
> (13) exactly. Note the weight goes as $a/\sigma^2$, **not** $a^2/\sigma^2$: weighting each window by its own
> SSNR — which is what a per-window FSC weighting would do — over-weights the early windows and is not the
> matched filter. Normalising to sum one rather than to the minimum-variance estimator keeps the output on the
> scale of a weighted mean of its inputs, with no per-shell amplification, and leaves the residual attenuation
> explicit and reported.

Per-shell diagnostics: Kish's effective number of windows (`kish_ess` in `src/utils/math/simple_stat.f90`), and
the attenuation of the accumulated map relative to window 1,

$$
n_{\mathrm{eff}}(s) = \frac{\big(\sum_w \omega_w(s)\big)^2}{\sum_w \omega_w(s)^2} = \frac{1}{\sum_w \omega_w(s)^2},
\qquad
\bar a(s) = \sum_w \omega_w(s)\,a_w^\star(s) .
\tag{14}
$$

Beyond the band, the weights of the last in-band shell continue. On request, the zero-dose extrapolation, applied
up to the band:

$$
S_0(\mathbf q) = \frac{A(\mathbf q)}{\bar a(s)\,\exp\!\big(B_0 k_s^2/4\big)} ,
\tag{15}
$$

where $B_0 < 0$ whenever the line reaches $B_1 = 0$ at $D_1 > 0$.

> **Why the halves stay separate.** The even output is a mean of even inputs only and the odd of odd inputs only,
> so the two share no noise and $\mathrm{FSC}(A_{\mathrm e}, A_{\mathrm o})$ is honest: it goes to `postprocess`
> or `enhance` like any half-map pair.
>
> **Why nothing is sharpened here.** A coefficient of the output is a weighted mean of the same coefficient over
> the windows; $\bar a(s)$ is reported so the residual attenuation is known; sharpening stays `postprocess`'s
> job. The zero-dose map (15) is a per-shell rescaling of $A$, so its FSC is unchanged by it — an interpretation
> product like `enhance`'s outputs, not a resolution claim, and `postprocess` would sharpen it a second time.

## 5 The control, and what would refute the premise

The standard answer on the same windows, in the same estimator, so the comparison is like for like: the
attenuation model of (2) with its parameters fixed by Grant & Grigorieff instead of measured,

$$
a_w^{\mathrm{GG}}(s) = \exp\!\left(-\frac{D_w - D_1}{2\,N_e(k_s)}\right),
\qquad c_w = 1,
\tag{16}
$$

with $N_e$ from (1) and the doses in the scaled coordinate of section 3.3. The difference $D_w - D_1$ appears
because the model's $a$ is relative to window 1 under the gauge (4) and the control must use the same gauge. Its
accumulation weights are (13) with $a^{\mathrm{GG}}$ in place of $a^\star$; `calc_dose_weights`' own
normalisation to unit sum of squares differs from (13) by a per-shell scalar, which changes no FSC.

For every target window $t$ and shell, the cross-half prediction error of predicting its odd half from the even
halves and vice versa, whitened by the target's own noise:

$$
e_t(s) = \frac{p_t(s)}{2\,n_s}\sum_{\mathbf q\in s}
\Big( \big\lvert O_t(\mathbf q) - \hat O_t(\mathbf q)\big\rvert^2
    + \big\lvert E_t(\mathbf q) - \hat E_t(\mathbf q)\big\rvert^2 \Big),
\tag{17}
$$

$$
\hat O_t(\mathbf q) = a_t(s)\,\hat S_{\mathrm e}(\mathbf q),
\qquad
\hat S_{\mathrm e}(\mathbf q) = \frac{\sum_v a_v(s)\,p_v(s)\,E_v(\mathbf q)}{\sum_v a_v(s)^2\,p_v(s)} ,
\tag{18}
$$

and $\hat E_t$ symmetrically from the odd halves. The attenuation $a$ in (17)-(18) is $a^\star$ for the model,
$a^{\mathrm{GG}}$ for the control, and $1$ for the flat mean — which is the full-dose reconstruction, the model
at $\beta = 0$. All three are quadratic forms in the Gram matrices of (7).

> **Why the target's own even half may be used to predict its odd half.** The two halves are built from disjoint
> particles, so there is no leakage; excluding it would only discard information.
>
> **Why the whitening.** A noise-free prediction leaves only the target's own noise, so
> $\mathbb E[e_t] = 1$; predicting a window from itself alone gives $\hat O_t = E_t$ and
> $\mathbb E[e_t] = 2$. The number is therefore comparable across shells, windows and specimens.
>
> **On honesty.** The model's error is mildly optimistic, since its $2W+2$ parameters were fitted on the same
> halves — a handful of numbers from of order $10^8$ coefficients. The error is therefore not held out, and
> section 6 reports its model-minus-control difference without a significance test. The control's is not optimistic at all: it has
> no fitted parameter.

Per resolution band, report the mean paired model-minus-control prediction-error difference. Negative values
indicate lower error for the measured model, positive values lower error for the control. This is a descriptive
comparison, not an independent validation or a calibrated significance test: attenuation parameters use both
halves, and masking correlates coefficients across shells. No z-score, significance threshold or categorical
winner is reported. Independent calibration data and reconstruction comparisons are needed to establish a gain.

The model's weights are delivered either way; the control's accumulation and FSC are written beside them.

**Four independent ways the premise can fail, each with its own readout.** Guinier lines curving consistently
(the exponent is not 2); the control winning the low shells (the loss has their shape there); $\tau$ large and
unstructured (the windows do not lie on a smooth dose curve); rank above one in (9) (the windows differ in more
than attenuation).

## 6 Products and readouts

`<outvol>` is `outvol` without its extension; `outvol` defaults to `dose_series.mrc`.

All products land in the execution directory unless `outvol` carries a path.

**Maps**

| File | Contents |
|---|---|
| `<outvol>_even.mrc`, `<outvol>_odd.mrc` | The accumulated halves $A_{\mathrm e}$, $A_{\mathrm o}$ of (13) under the measured model's weights; even from even inputs only, so they are an honest half pair for `postprocess` or `enhance` |
| `<outvol>.mrc` | Their mean $A$, the delivered map: the dose-weighted average of the windows, each shell weighted by (13) with weights summing to one; unfiltered, unsharpened, unmasked as the inputs were |
| `<outvol>_fsc.bin`, `<outvol>_fsc_masked.bin` | The FSC of the model's half pair, unmasked and under the soft mask of `mskdiam`, in SIMPLE's binary format |
| `<outvol>_control_even.mrc`, `_control_odd.mrc`, `_control.mrc`, `_control_fsc.bin`, `_control_fsc_masked.bin` | The same under the Grant & Grigorieff weights of section 5, for a like-for-like comparison |
| `<outvol>_fsc.txt` | Per shell: resolution and the FSC of the model's and the control's half pairs, each unmasked and masked |
| `<outvol>_zerodose.mrc` | With `zerodose=yes`: the extrapolation (15) of $A$ to zero dose; a per-shell rescaling with the same FSC, an interpretation product rather than a resolution claim |

**Tables**

| File | Contents |
|---|---|
| `<outvol>_model.txt` | Against the sharpest window (section 11.9), per window: dose $D_w$, scale $c_w$, $B_w \pm \mathrm{se}_w$, the line's $B(D_w)$, the $B_w^\star$ used, residual $z_w$, Guinier rms, and `used` (0 for a rejected window); last two lines the dose line $B_0$, $\beta$, $\mathrm{se}_\beta$, $\tau$ and the contrast line $\ln c_0$, $\gamma$, $\mathrm{se}_\gamma$, $\tau_c$ of section 12 |
| `<outvol>_curve.txt` | The dose curve `reextract` reads with `dose_weighting=measured` (section 12): voltage, physical dose per fraction, band, both lines, and per window its fractions, scaled centre dose, $c_w$, $B_w^\star$ and `used` |
| `<outvol>_attenuation.txt` | Per shell: every window's measured attenuation $a_w(s) \pm \mathrm{se}_w(s)$ relative to window 1 — the model-free Guinier diagnostic of (9) |
| `<outvol>_weights.txt` | Per shell: the model's weight $\omega_w(s)$ of every window, the effective number of windows $n_{\mathrm{eff}}(s)$, and the accumulation's residual attenuation $\bar a(s)$ |
| `<outvol>_weights_res.txt` | Per window: its dose and, at the shells nearest 20, 10, 8, 6 and 4 Å, the model's and the Grant & Grigorieff weights; NaN past Nyquist |
| `<outvol>_cv.txt` | Per shell: the whitened cross-half prediction error $e_t(s)$, averaged over $t$, for the model, the control and the flat mean, and the paired mean model-minus-control error difference (descriptive only) |
| `<outvol>_decay.txt` | Per band shell: the model-free critical exposure $N_e(s) \pm$ error, the lines' $1/(2\gamma + \beta k_s^2/2)$ (section 12), Grant & Grigorieff's $N_e(k_s)$, and the per-window model's through the same estimator (`ne_model_est`); then per window the measured attenuation, its error, the model's and theirs |
| `<outvol>_decay_steps.txt` | Per pair of neighbouring windows: its mean dose and, at the band shells nearest 20, 10, 8, 6 and 4 Å, the critical exposure between the two windows with its error and Grant & Grigorieff's; fading exponential in dose would make each column flat |

**Figures**

| File | Contents |
|---|---|
| `<outvol>_decay.png` | Left: attenuation against spatial frequency for up to five windows — points measured, solid the B-factor model, dashed Grant & Grigorieff. Middle: critical exposure per shell, measured without a model, against the per-window model through the same estimator (solid), the lines' (dotted) and Grant & Grigorieff's (dashed). Right: the critical exposure between neighbouring windows against their mean dose at 20, 10, 8, 6 and 4 Å, dashed Grant & Grigorieff's, which is flat |
| `<outvol>_weights_res.png` | The accumulation weight of each window against its dose at 20, 10, 8, 6 and 4 Å — solid the model, dashed Grant & Grigorieff, dotted equal weights |
| `<outvol>_plots.gp`, `<outvol>_plots_gnuplot.log` | The gnuplot commands that reproduce both figures from the tables, and gnuplot's log; a rendering failure warns and keeps the tables |

The per-shell critical exposure regresses $\ln a_w(s)$ on $D_w - D_1$ through the origin over the windows that
still resolve the shell, $N_e = -1/(2\,\mathrm{slope})$ as $\exp(-D/2N_e)$ has it, in the scaled dose
coordinate of section 3.3; a shell without decay beyond twice its error is left empty. It assumes no B-factor,
so where the points leave the model's curve the premise fails at that resolution.

Before the fit the log gives the mask's two factors: `INDEPENDENT FRACTION OF A SHELL UNDER THE MASK`, the $f$
of section 4.2, and `VARIANCE INFLATION ACROSS SHELLS UNDER THE MASK`, the factor of section 11.8 applied to the
Guinier errors.

Log lines, prefixed `Dose series:`. Per window: $c_w$, $B_w \pm \mathrm{se}$, $B(D_w)$, $B_w^\star$, $z_w$, the
Guinier rms. The line: $B_0$, $\beta \pm \mathrm{se}$, $\tau$, and $\beta$ against the equivalents of (5) at 2, 3
and 10 Å. The contrast line: $\ln c_0$, $\gamma \pm \mathrm{se}$, $\tau_c$. Per resolution band, in eight bands: $n_{\mathrm{eff}}$, $\bar a$, the three cross-half errors, and the
paired mean model-minus-control prediction-error difference, without a significance-based verdict. Shared
fitting data and cross-shell correlations preclude interpreting this difference as independent validation. The FSC 0.5
and 0.143 of the model, the control, the flat mean and window 1 alone, unmasked and under the soft sphere of
`mskdiam` with zero background, the convention reconstruct3D applies before its FSC.

## 7 Cost

The Gram pass of (7) is $2W^2 n_{\mathrm{coeff}}$ complex multiply-adds: at $W = 20$ and box 256 that is
$7\times10^9$, seconds with OpenMP over the coefficients. The fits of (9), (10) and (12), the control and every
FSC of section 5 are arithmetic on $n_{\mathrm{sh}}$ matrices of $2W \times 2W$. The accumulation (13) is
$W\,n_{\mathrm{coeff}}$.

Memory: the halves as complex single precision, $2W n_{\mathrm{coeff}} \times 8$ B $=$ 2.7 GB at $W = 20$,
box 256. The masked halves are held for the Gram pass only and the halves as given for the accumulation only, so
the peak is one set; Fourier-crop the halves to the box the band needs when memory is tight.

## 8 Standing assumptions and limitations

- **Shared poses.** The windows are reconstructed with the poses of the full-dose refinement, so the halves of
  different windows share the alignment errors of their half. Sections 4.2-4.5 and 5 use cross-half products
  only and are unaffected; the accumulated map is bounded by that refinement's alignment, like any reconstruction
  from it.
- **$\beta$ is in the scaled dose coordinate** of section 3.3: a drop-in for `calc_dose_weights`, but not a
  physical specimen constant and not comparable across voltages.
- **Window-centre dose.** $D_w$ is the dose at the window's centre, for both the model and the control. A window
  spans a dose range and its reconstruction averages frames whose attenuation differs; under the line the window's
  true attenuation is the mean of $\exp(-\beta D k^2/4)$ over its range,
  $$\exp\!\left(-\frac{\beta D_{\mathrm{lo}} k^2}{4}\right)\frac{1 - \exp\!\big(-\beta \Delta D k^2/4\big)}{\beta \Delta D k^2/4},$$
  not the value at the centre. This matters once $\beta\,\Delta D\,k^2/4$ passes about $0.1$ — a 3 e/Å$^2$ window
  at 3 Å with $\beta = 4$ gives $0.33$ — and curves a window's Guinier line slightly. Section 4.5 can fit the line
  with the window-integrated attenuation at no extra cost; section 4.4's per-window B-factor stays the model-free
  diagnostic. Deferred, and recorded here as a known bias in $\beta$.
- **float16 particle stacks.** Stacks are now written unconditionally in MRC mode 12, so every window is quantized
  identically and the attenuation ratios survive, but the high-resolution noise floor in (6) rises for all of them.
  To be checked once against a float32 control series, band by band.
- **One state** (state 1), **disjoint windows**, **one box and sampling**, **one even/odd particle partition**.
  The consumer sees only volumes, so it cannot verify the partition itself; the guarantee is the in-run check in
  `reconstruct_motion_fractions`, which throws if the selection or the `eo` labels change between intervals, plus
  the partition fingerprint in the manifest.
- **Everything is relative to window 1** under the gauge (4). Reported B-factors and attenuations are not absolute.
- **Deterministic.**

## 9 Open points

Settled with the maintainer, 3 November 2026: raw `_unfil` halves as the input (with `refine_motion_model` adjusted
to make that pair genuinely unfiltered and unmasked, leaving `_even`/`_odd` and therefore the exported FSC matrix
untouched); a producer-written manifest for the dose; dropping the window before an overlapping last one; state 1
only; the band rule of section 4.1; detector noise independent across windows with structure error shared; $\beta$
in the scaled dose coordinate; window-centre dose; window 1 excluded from the line fit.

Still open:

- **The exponent.** $a_w = c_w\exp\!\big(-B_w k^\gamma/4\big)$ with $\gamma$ a third global parameter
  ($\gamma = 2$ is the B-factor, $1.665$ the Grant & Grigorieff shape). Fit it only if the Guinier lines of the
  real series curve consistently, or the control wins the low shells. One parameter, the same (10).
- **The scale.** With equal windows and one reconstruction normalisation, $c_w$ should be 1 throughout. If the real
  series confirms it, fix it and drop a parameter per window.
- **`automsk`.** Whether the first release supports an automask, and if so which map it is built from.
- **The independent count.** Resolved as neither $\langle M\rangle$ nor $\langle M^2\rangle$: the plane formula of
  section 4.2, with a residual of about 10% recorded in section 11.1.
- **`bfac_model`.** Whether it is a CLI key at all, and whether the name reads too much like the existing `bfac`
  sharpening parameter.
- **Feeding $\beta$ back.** Section 1's third aim needs $\mathrm{twoNe} = 4/(\beta k^2)$ in place of
  $2(Ak^B + C)$ in `calc_dose_weights`, behind one key. In the first release, or named as the follow-up?
- **Per-window doses** for unequal windows and EER fractions, beyond what the manifest already carries.
- **Pre-exposure.** Dose received before the first frame shifts every $D_w$ by a constant: the slope $\beta$ is
  unaffected, but the zero-dose extrapolation (15) then refers to the start of the movie, not to the undamaged
  specimen.
- **Rank above one.** Where the Guinier lines curve and the exponent does not straighten them, the loss is more
  than attenuation; a rank-two extension of (9) is the next model, and the readout decides whether it is needed.
- **Local B-factors.** Damage is site-specific and one B-factor per window averages over the particle. `nu_filt3D`
  is the vehicle for the local version: a bank of dose curves per voxel in place of its bank of cutoffs, the same
  even/odd cost, the same Potts smoothing of the labels.


## 10 Producer integration

Status: implemented and validated in the working tree. Scope: make
`refine_motion_model` emit genuinely unmasked `_unfil` half maps and a complete, portable dose manifest,
and verify that `dose_filter3D` reads those products correctly. The existing attenuation fitting and
weight equations were preserved. Existing consumer changes, including its high-pass fitting bound,
were retained.

### 10.1 Unmasked half-map export

1. Add the internal yes/no parameter `unfil_unmasked`, default `no`, through the existing declaration,
   parser and validation modules. It has no application UI entry. Enabling it requires the gridding,
   nonregularized reconstruction route; reject incompatible use rather than silently ignoring it.
2. Set `unfil_unmasked=yes` on every child reconstruction in `reconstruct_motion_fractions`. Both
   reconstruction dispatch paths already forward the command line to `assemble_gridding` and
   `run_volassemble`; test that forwarding instead of adding another dispatch interface.
3. In the nonregularized branch of `restore_gridding_pair`, write each `_unfil` half after
   `restore_base` and `finalize_from_base`, before `mask3D_soft`. Write standard real-space float32 MRC
   volumes (mode 2). The corrected image can then be masked in place: no extra full-volume copy is needed.
4. Continue writing the standard even/odd maps and computing their FSC with the existing spherical
   support. Retain the existing merged reconstruction. Thus `_even_unfil`/`_odd_unfil` are unmasked,
   while standard maps, FSC binaries and the raw FSC surface retain their current scientific meaning.
   The export applies to populated states; the first dose manifest describes state 1 only.
5. Record `solve_support=none` and `solve_kind=gridding` for the unmasked halves. Extend the existing
   support-provenance writer/reader to represent `none` explicitly. Update the postprocessing check
   that currently treats the mere presence of a provenance sidecar as proof of support. Replacing an
   output must replace or remove its old provenance too. Existing sphere/density/NU cases remain covered.

### 10.2 Physical dose and selected-particle metadata

Collect metadata while the first copied interval project is already loaded in
`reconstruct_motion_fractions`, after `sample4rec` returns the actual reconstruction rows. Use the
state-1 subset of those rows, not every active particle in the source project.

- Use particle `stkind` and the existing project stack-to-micrograph mapping to identify distinct
  contributing micrographs. Read each contributing saved binary motion model once.
- The authoritative values are `motion_model%dose_per_frame` and `%voltage`. Require finite positive
  dose, a common supported voltage, sufficient saved frames, and a common dose per fraction within a
  named numerical tolerance. Do not average incompatible exposures or reject unused micrographs for
  metadata that cannot affect the state-1 maps.
- `dose_per_frame` already describes one decoded fraction, including EER grouping. Do not multiply it
  by `eer_fraction` or use `target_dose_per_frame` in its place.
- For original inclusive fraction indices `[a,b]`, compute in double precision:

  ```text
  dose_lo = (a - 1) * saved_dose_per_frame
  dose_hi = b       * saved_dose_per_frame
  ```

  These are physical e/Angstrom^2, without voltage scaling. If `fromf > 1`, retain the dose preceding
  the first exported interval. The saved model does not encode additional external pre-exposure.
- Verify that the state-1 selected particle rows and half labels are identical across intervals.
  Existing global row/half checks do not by themselves detect particles changing state membership.
- Record state-1 `nptcls`, `neven`, `nodd` and `eo_digest`. Reuse `string%to_fnv1a_hash64()` on an
  unambiguous serialization of ordered `(particle row, half label)` pairs. The digest records a
  producer-verified invariant; a consumer reading only volumes cannot independently prove it.

### 10.3 Shared manifest contract and reader validation

Introduce a small encapsulated `dose_series_manifest` type in
`src/main/volume/simple_dose_series_manifest.f90`, with private metadata/window storage, lifecycle,
read/write, validation and only the accessors the two callers need. Move the existing manifest/path
parsing helpers out of the dose-filter commander into that owner. Keep numerical fitting in
`simple_dose_series`; do not duplicate parsing rules or introduce another numerical fitting type.

The output is `dose_series.txt` at the reconstruction run root, using the version-1 format in section
4.1.1. It carries program, state, box, sampling, voltage, particle/half counts, partition digest and
physical dose per fraction. Each row has exactly the six columns `from to dose_lo dose_hi even odd`.
Paths are relative to the manifest directory; the reader also accepts absolute paths. Version 1 uses
space/tab separators and does not support whitespace inside path tokens. Detect unsupported or
truncated input explicitly. Required metadata and data rows must round-trip through the shared reader.

Validate the exact version marker, numeric parse status, finite values, ordered positive fraction
indices, positive dose spans, dose consistency, state 1, positive half populations and count consistency.
Keep all actual interval rows, including the final overlapping interval. Validate that any overlap is
only the supported final/predecessor case. The consumer drops that predecessor, then checks that all
survivors are disjoint and at least two remain. Record the actual selected intervals in its log.

Preflight every distinct even/odd file using `MrcImgHead` before calling `image%read`: dimensions,
real-space mode, all three cell/grid-derived pixel sizes, isotropy, common sampling, and consistent
axis/origin conventions. Producer outputs are mode 2; real-space float16 mode 12 can remain readable,
but complex Fourier modes are not valid inputs to this path. Check that each pair names distinct files.
Use named tolerances for header/dose comparisons, not scattered literals.

This preflight replaces an ineffective check: the current `prepare_half` compares the dimensions of
an image preallocated from the manifest, while image I/O retains that geometry and substitutes its
supplied pixel size in the in-memory header. No change to generic image I/O is needed. This work also
makes the newly recorded provenance available for checking supported producer outputs.

### 10.4 Publication and failure behavior

The master owns manifest publication. Accumulate rows in memory, verify the produced unmasked pairs
as intervals finish, and publish the final file only after every referenced reconstruction and the
required metadata have succeeded. Write a temporary file in the destination directory, close it and
rename it into place. A failed reconstruction must never leave a seemingly complete partial manifest;
handle any existing manifest explicitly so it cannot describe replaced or incomplete products.

The current application policy forbids dose-compatibility gates rejecting raw maps/FSC. Therefore
missing or heterogeneous dose/voltage, an unavailable state-1 half pair, an unsupported state-1 selection
change, or fewer than two disjoint usable windows prevents this new manifest from being published and
prints an explicit reason. Keep valid reconstructed maps, raw FSC outputs and retained stacks. Do not
invent dose, silently remove contributing particles or publish a subset as though the run completed.
Existing fatal reconstruction/partition corruption checks retain their current behavior. Manifest I/O
failures must be reported explicitly and must not leave a consumable final file.

### 10.5 File impact and implementation order

Implementation sequence:

1. Shared manifest type and contract tests, then connect the consumer to its reader/header checks.
2. Internal output parameter, child-command setting, pre-mask half-map writes and provenance support.
3. Producer metadata collection, interval verification and final manifest publication.
4. End-to-end producer/consumer validation and updates to the application algorithm/policy.

Expected production files:

- `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`,
  `simple_parameters_phases.f90`: internal output setting and validation.
- `src/main/strategies/parallelization/simple_refine_motion_model_strategy.f90`: enable export,
  collect authoritative metadata, verify interval identity and publish the manifest.
- `src/main/strategies/parallelization/simple_volassemble_strategy.f90`: write corrected halves before
  support, with the existing standard half-map/FSC path preserved.
- `src/main/volume/simple_halfmap_diagnostics.f90`: explicit unmasked support provenance.
- New `src/main/volume/simple_dose_series_manifest.f90`: shared file contract and header validation.
- `src/main/commanders/simple/simple_commanders_resolest.f90`: use the shared reader and validated inputs.
- `src/main/commanders/simple/simple_commanders_volops.f90`: distinguish `support=none` from applied
  support when interpreting provenance.

Tests extend the existing reconstruction/provenance suites and add the manifest tester to the existing
`unit_reconstruction` area, with producer metadata in `unit_motion` and export integration in
`lib_reconstruction`; no new executable or CMake target is needed. Register its suite in the normal test
registry/help locations. Update this living note, `doc/algorithms/refine_motion_model.md` and
`doc/policies/refine_motion_model_policy.md` as the implementation is completed.

Acceptance checks for the producer are listed in section 13.

### 10.6 Implementation notes

The reader inspects actual MRC headers before allocating image buffers; metadata and file errors return to the
producer as a manifest-publication diagnostic. The producer refuses an existing final manifest before writing new
reconstructions, validates headers at final publication and publishes through a temporary file, sync and atomic
rename, so no incomplete manifest is left. Canonical paths and symlink aliases cannot reuse a half-map file;
distinct hard links to one inode are not detected.

## 11 Consumer implementation

`simple_dose_series` (`src/utils/filter/`) is one encapsulated class holding the statistics, the fits, the three
accumulation schemes (model, control, flat), the paired cross-half comparison and the zero-dose rescaling. The
commander in `simple_commanders_resolest.f90` reads the manifest through the shared reader of section 10.3, gathers
and scatters the half-volume lattice and writes the products of section 6. Keys: `bfac_model` (`line|free`),
`zerodose` (`yes|no`), `hp` and `lp`; the box comes from the manifest. Departures from sections 4-6: the band
starts at `hp`, default 20 Å; the statistics copy is masked about its solvent level; every statistic is computed
after rescaling the coefficients to unit mean power, since an absolute guard on SIMPLE's FFT scale zeroed every FSC.

### 11.1 Error model

The errors of $a_w(s)$ use the noise-only pair variance (8a), not the sample variance of the cross products, which
holds the structure's own power fluctuation; the sandwich covariance (9a), since the projection of each window's
noise onto the structure is shared by every pair holding it; and the independent count of section 4.2. Masked
errors remain about 10% small, uniformly over shells and noise levels; about half is the coherence of the noise's
projection onto the structure, the rest is unexplained.

### 11.2 Guinier fits

Levenberg-Marquardt with step acceptance replaces undamped Gauss-Newton, which diverged on a window carrying no
signal. A fit that fails, or whose scale does not clear the detection gate of section 11.8, is rejected: its window
is left out of the line and of the model's weights, the model table marks it `used` 0 and the log says so.

### 11.3 Excess scatter about the dose line

$\tau^2$ is the DerSimonian-Laird estimate with the trace of a weighted line,
$\mathbb E[Q] = (k-2) + \tau^2\,\mathrm{tr}(P)$, $\mathrm{tr}(P) = \sum w -
\mathrm{tr}\big[(X^{\mathsf T}WX)^{-1}X^{\mathsf T}W^2X\big]$, $X = [\mathbf 1, D]$ (`line_tau2`). An intercept-only
trace left $\tau^2$ low by about $(k-2)/(k-1)$. The untruncated estimate is unbiased under the variance model;
clamping it at zero and taking the root bias $\tau$ with few windows.

### 11.4 Coefficient count

The count behind the unit-mean-power rescaling is formed in double precision: as a default integer it overflowed
past about box 590 with 20 windows.

### 11.5 Input guards

`dose_filter3D` stops when the band's noise is below $10^{-6}$ of the half-map power (one map given as both halves)
and when no window after the first can be fitted. `scripts/make_dose_series.py` refuses the same file, identical
data, and halves whose difference is below $10^{-6}$ of their power.

### 11.6 Figures

The decay figure and the weights figure of section 6 are rendered by gnuplot under Fortran control, as in
`refine_motion_model`; `posix_quote` in `simple_string_utils` quotes the paths. A rendering failure warns and keeps
the tables.

### 11.7 Test tools

`scripts/make_dose_series.py` builds a dose series from two experimental half maps with a known attenuation:
truths `model`, `integrated`, `fractions` (the mean over each window's fraction centres) and `gg`, with `--beta`,
`--gamma` (contrast fall) and `--b2` (quadratic B); `truth.txt` records the fraction-level $B$ and $c$ at each window
centre. `scripts/dose_filter3D_matrix.sh ROOT` runs the condition matrix (seeds, window widths, curvature with
same-seed unaveraged controls, contrast, voltage, Grant and Grigorieff truth, noise, layouts, options, guards, and
with `real` the re-extraction of real particles with its refusals); `scripts/dose_filter3D_matrix_summary.py ROOT`
tabulates truth against fit with a verdict per condition and the paired window-averaging bias.

### 11.8 Window errors across shells, and the detection gate

The mask correlates coefficients across neighbouring shells as well as within one, so a Guinier fit smooth over
many shells has a variance inflated by the sum of the correlations between shells: by Parseval and the
projection-slice theorem, per axis $a$,

$$
F = \frac{1}{3}\sum_a \frac{n_a \sum_{\mathbf x} M^4}{\sum_{\perp a}\big(\sum_a M^2\big)^2},
$$

1 without a mask, $(N/t + 2)/3$ for a slab of thickness $t$, about $4N/(3D)$ for a ball of diameter $D$.
`shell_correlation` computes it beside `shell_independence`, and both variances of each window's Guinier fit are
multiplied by it, which carries into the shrinkage, the residuals $z_w$ and the rejection gate. A window's scale
must clear five of its errors to count as measured (`C_DETECT`); the shells entering the log-domain start of the
fit keep twice their error.

### 11.9 Reports against the sharpest window

Window 1 is the gauge of the fit ($a_1 = 1$, $c_1 = 1$, $B_1 = 0$), but early beam-induced motion can leave it
blurrier than the next window, whose $B$ then reads negative against it. The model table and the log report
against the fitted window of least $B^\star$: $c_w/c_{\mathrm{ref}}$, $B_w - B^\star_{\mathrm{ref}}$,
$B(D_w) - B^\star_{\mathrm{ref}}$, $B^\star_w - B^\star_{\mathrm{ref}}$, $B_0 - B^\star_{\mathrm{ref}}$ and
$\ln c_0 - \ln c_{\mathrm{ref}}$, with the errors of section 12.1; the header names the reference. The shift leaves
every ratio between windows, $\beta$, $\gamma$, $\tau$, the weights, the maps and the dose curve as they were; the
fit, the attenuation and decay tables and the curve file keep window 1 as their gauge. The line's prediction
$B(D)$ can still be negative for early windows when the fading accelerates.

## 12 Re-extraction under the measured curve

`reextract` with `extractfrommov=yes dose_weighting=measured dose_curve=<outvol>_curve.txt` weights the movie
fractions of every particle by the measured curve instead of Grant and Grigorieff's.

**The contrast line.** The windows' scales fall with dose as their B-factors rise. $\ln c_w = \ln c_0 - \gamma
D_w$ is fitted over the measured windows $w \ge 2$ as the line of section 4.5 is: two passes, weights
$1/(\mathrm{se}_{\ln c,w}^2 + \tau_c^2)$ with $\mathrm{se}_{\ln c,w} = \mathrm{se}_{c,w}/c_w$ from the Guinier fit
and $\tau_c^2$ by the slope-aware trace of section 11.3. It serves only past the last window. A frequency-flat
fall is expected, not an artefact: a fading law in which $\ln a$ is not proportional to $k^2$ across the band
shows partly as an intercept; Grant and Grigorieff's curve, fitted the same way, shows one.

**The model's critical exposure.** Under both lines the model fades as $\exp(-(\gamma + \beta k^2/4)D)$, so its
critical exposure is $N_e = 1/(2\gamma + \beta k^2/2)$. Once the fading is not exponential in dose this and the
model-free curve measure different things: the lines' slope is the mean rate of windows 2 to $W$, the model-free
estimate a secant from window 1 weighted by $(a_w/\mathrm{se}_w)^2 (D_w - D_1)^2$, which leans on the early windows
at high frequency. The decay table and figure therefore also carry the per-window model ($c_w$, $B_w^\star$) read
by the same estimator, and `<outvol>_decay_steps.txt` the critical exposure between neighbouring windows, which is
flat when the fading is exponential in dose. Fading that accelerates mostly in $B$, with $B$ rising roughly as the
square of the dose, is what an uncorrected displacement growing linearly with dose gives
($B = 4\pi^2 v^2 \Delta D^2$); re-measuring after per-particle polishing would separate it from damage. The weights
of section 4.6 and of the re-extraction use each window's own $c_w$ and $B_w^\star$, not the lines, so a curve
measured before polishing does not hold after it.

**The file.** For example:

```
# SIMPLE dose curve 1
# program dose_filter3D
# units: dose in e/A^2 scaled as calc_dose_weights; B in A^2; band in A
# kv  2.0000000000000000E+002
# dose_per_fraction  1.2770832777023315E+000
# band  2.0423076702998234E+001  3.4038461171663723E+000
# line -6.6003819451131463E+001  7.6771854953970111E+000
# contrast  6.9855347773784701E-002  1.5744537108473067E-002
# nknots 8
# columns window from to dose c b used
1 1 3  2.3945311456918716E+000  1.0000000000000000E+000  0.0000000000000000E+000 1
2 4 6  7.1835934370756149E+000  9.4479002397246359E-001  6.0220362354875947E+000 1
...
8 22 24  3.5917967185378075E+001  5.9879111003567587E-001  2.2606929978079734E+002 1
```

`dose_per_fraction` is the manifest's physical dose; the band is in Å; `dose` is the scaled window centre;
`c` and `b` are $c_w$ and the $B_w^\star$ the model map used. A rejected window is written with `used` 0 and is
no knot. The reader (`simple_dose_curve`) refuses an unknown version, a missing or repeated header, a malformed
row, a window 1 off the gauge, doses that do not increase, fewer than two fitted windows, and a voltage other
than 100, 200 or 300 kV.

**Per fraction.** Fraction $f$ of a movie with physical dose $d$ per fraction sits at the scaled centre dose
$x_f = (f-\tfrac12)\,d/s_{\mathrm{kV}}$; a window's centre is the mean of its fractions'. Between fitted windows
$\ln c$ and $B$ are linear in dose, before the first window the first segment continues, and past the last the
two lines' slopes continue from it, clamped at zero so that an unmeasured fraction never gains on the last
measured one. Under the lines, equal-width windows need no correction for averaging over their fractions: each
carries the factor $\sinh(n\lambda\delta/2)/(n\sinh(\lambda\delta/2))$, $\lambda = \gamma + \beta k^2/4$, which
cancels against window 1, so the slopes measured at window centres are the per-fraction ones. With curvature,
interpolating between centres errs by $\tfrac12 B''\,\mathrm{Var}(x)$.

**Weights.** $\omega_f(k) = c_f e^{-B_f k^2/4} / \sqrt{\sum_g c_g^2 e^{-B_g k^2/2}}$, of unit squared sum at
every frequency like the Grant and Grigorieff mode, so the summed noise spectrum stays flat
(`image%apply_frame_bfactor_weights`, beside `apply_dose_weighing`). The curve was measured on plain sums of
fractions, so it also replaces the motion-correction frame weights: the weighted frames are summed equally.

**Checks.** The curve holds for one voltage and one dose rate, the producer's rule: a micrograph whose motion
model differs in either stops the run. A micrograph without a motion model stops it too, instead of falling back
to the integrated micrograph, which carries Grant and Grigorieff's weights. `dose_curve` is made absolute before
the run directory is entered, and the distributed master reads the file before it schedules a job.

### 12.1 Uncertainty and comparison corrections

The joint sandwich covariance of the shell attenuations is retained. For the through-origin critical-exposure
regression, with fixed fitted weights, the slope error is propagated as

$$
\operatorname{Var}(\hat m)=g^T\Sigma_a g,\qquad
g_w=\frac{w_w(D_w-D_1)}{a_w\sum_v w_v(D_v-D_1)^2}.
$$

This includes the common noisy reference. The neighboring-window log-ratio error includes the corresponding
negative twice-covariance term. First-order errors still assume the retained fitting observations and their
weights are fixed; the two-error decay gate now uses the propagated error.

Guinier B-factor covariance follows from each fit's sensitivity to its shell attenuations, with the existing
cross-shell variance inflation. The reporting shift uses the sensitivity of the reference's shrunk B-factor to
the measured B-factors, including the dose-line fit. The reported variance is that of
$B_w-B_{\mathrm{ref}}^\star$, with line weights and excess scatter held fixed. Only a genuine identical-estimator
difference has exactly zero error; the reference's measured B residual generally does not.

The comparison table's final column is now `error_model_minus_control`, the paired mean prediction-error
difference in whitened units. The log reports the same quantity per band and makes no significance-based
model/control/tie verdict. These comparisons remain descriptive until independent validation and a calibrated
covariance model are available.

Both stored curve voltage and movie voltage must be finite; movie dose per fraction must be finite and
positive before compatibility comparisons.

## 13 Validation

To be redone after the changes above. Unit suites: `unit_numerics suite=dose_series` and `suite=dose_curve`,
`unit_image suite=image`, `unit_motion suite=particle_extractor`, `unit_motion suite=fraction_dose_metadata`,
`unit_reconstruction suite=dose_series_manifest` and `suite=refine3d_provenance`; integration
`lib_reconstruction suite=dose_filter3d`, `suite=fraction_dose_export` and `suite=fractional_reconstruction`; the
fast gate. Then:

- Producer: standard maps and FSC unchanged by the exports; the `_unfil` halves keep density outside the sphere and
  reproduce the standard halves when masked; the manifest survives relocation; the final-overlap rule; malformed
  inputs refused; no manifest after a failed or inconsistent run; a real producer output accepted by
  `dose_filter3D`.
- `dose_filter3D` on synthetic series (section 11.7): $\beta$, $\gamma$ and $B_w$ recovered within their errors;
  error calibration over seeds; window averaging unbiased under exponential fading; guards stop.
- `dose_filter3D` on real series: the model and decay tables, the figures, the masked FSC against reconstruct3D's.
- `reextract`: Grant and Grigorieff against the measured curve, each reconstructed with fixed poses and compared
  by FSC; the refusals of section 12.
