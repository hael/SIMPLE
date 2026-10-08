# SINGLE without a species input: atom detection, Gaussian density model, species discovery and map simulation

Implementation note, 2026-10-07, revised 2026-10-08 after review. Written
against master `108c8c9d`. Companion page: `single_pcg_integration.md`
(section 11 says how the two fit together).

Status. The design is for review and is not approved as a replacement of the
present detection. What is approved is a small, opt-in prototype: the
foundations and the species-free level 1 of section 10 (phases 0 to 2, which
change nothing for an existing run) and the residual recovery and species
call behind an experimental key, producing diagnostics only (phase 3). The
prototype is developed alongside the existing detection machinery, which
keeps working unchanged throughout; section 10.0 states the rules. Promotion
of anything beyond diagnostics waits for the validation of section 10.5.
Section 15 records the review and what it changed.

Nothing in SIMPLE was compiled or run. The numbers in section 8 come from a
numpy emulation of the proposed pipeline on synthetic particles
(`emul_species.py`, which needs `emul_detect.py`; both are kept out of the
tree, like the other reference scripts), checked in section 8.8 by a second,
independent implementation (`codex_checks.py`). They show how the algorithm
behaves on idealised input and are not a measurement on data: the generator
produces the Gaussian atoms the detector assumes, so agreement between the two
says the design is coherent, not that it works on a reconstruction.

## 0. What is asked, and what constrains the design

The request: no `element` on the command line. The program finds the atoms,
finds how many kinds of atom there are and which atom is which, and simulates
a map from what it found. The first data are HAADF-STEM images (high-angle
annular dark-field scanning transmission electron microscopy, whose contrast
rises with atomic number) of Pt/Ni particles.

Two further points from the discussion shape the design.

- The B factor rises from the core to the surface, so anything that matches
  the density against a fixed shape has a different error rate at the surface
  than in the core. The density should be modelled by a generic Gaussian fit
  whose parameters can be examined as a function of radius. This is adopted:
  every atom gets its own fitted width, the species call does not use
  matching, and the radial profiles of the fit are a product (sections 1, 3.7,
  3.8, 7).
- An alternative is to specify the species, take the maximum over per-species
  correlation maps for the segmentation, and assign species with a Gaussian
  mixture model that includes the B factors. The mixture model with embedded
  B factors is what sections 3.6 and 3.7 are. The maximum over per-species
  maps does not change the segmentation, and section 5 says why and what
  specifying the species does buy.

"Species" in this note means an intensity class. HAADF contrast orders the
atoms by atomic number, so the data can say that there are two kinds of atom,
which is which, and what their intensity ratio is. They cannot name the
elements without a calibration, and nothing in the density model needs the
names (sections 5 and 12).

"Pt:Ni is 3-6:1" can be read as the atom ratio or as the intensity ratio.
The design does not depend on the reading: the emulation covers atom ratios of
3:1 and 6:1 crossed with intensity ratios of 3:1 and 6:1. Which one is meant
matters for what to expect from the data (Q1).

Constraints taken from the code:

- `refine3D_nano` has no species dependence
  (`single_commanders_nano3D.f90:269-295`). The species enters
  `autorefine3D_nano` only through `detect_atoms` and the `_SIM.mrc` it
  returns as `vol1`, and once upstream, in the `startvol.mrc` of
  `analysis2D_nano` (`single_commanders_nano2D.f90:196-203`).
- The simulated map is an alignment reference only. Trailing reconstruction
  blends accumulators, never volumes
  (`.github/skills/simple-frac-update-trailing/references/frac-update-contract.md`),
  so `_SIM.mrc` puts no density into the reconstruction. The intensities
  measured in the next map come from the particle images. (As read in
  `derive_sampling_settings`, `simple_parameters_phases.f90:504-512`,
  `l_trail_rec` is cleared when `update_frac` is above 0.99, and
  `autorefine3D_nano` sets `trail_rec=yes ufrac_trec=0.5` but not
  `update_frac`. If that reading is right, trailing is off in the default
  SINGLE loop. Not traced beyond that routine.)
- The volumes `refine3D_nano` ships are zero outside a soft sphere of radius
  `msk_crop` (`simple_commanders_rec_distr.f90:461-466`, `:638-682`), and
  `nanoparticle%new` masks again (`simple_nanoparticle.f90:271-278`). Any
  noise reference has to come from inside that sphere and away from the
  atoms.
- `image%phase_corr` divides by one global norm that depends on the map
  (`simple_image_calc.f90:1482-1485`). Correlation maps of different residuals
  are therefore not on one scale until they are standardised.
- `atoms%convolve` reads neither `beta` nor `occupancy` and silently skips
  any `Z` outside its table (`simple_atoms.f90:1235-1236`); `set_name` and
  `set_element` stop on an unknown symbol (`:727-743`, `:890-904`).
- A single-species data set must come out as it does today: one class, no
  atoms added.
- `discard_atoms` deletes an atom by clearing its component in `_BIN`,
  relabelling `_CC` and calling `find_centers`, which deallocates and
  reallocates `atominfo` (`simple_nanoparticle.f90:767-768`, `:1116`,
  `:1151`); `identify_atomic_pos` already has to recompute `valid_corr`
  after it for that reason (`:601-603`). Nothing fitted per atom survives
  pruning unless pruning comes first or works on the atom table. Its size
  test compares the longest extent of the component with half the atom
  radius (`:1132-1139`), so a component that is generated from the model
  can never fail it.
- The regression test of the present path is `single_atoms_stats`
  (`production/CMakeLists.txt:307-308`): it simulates a Pt particle in
  process (`simulate_nanoparticle`, box 160, `moldiam=20`), runs
  `detect_atoms element=Pt` and `atoms_stats`, and requires recall and
  precision of 0.90 or more within 1 A, a root-mean-square position error
  of at most 0.75 A and a correlation of at least 0.80 between `_SIM.mrc`
  and the input. The last recorded result is recall 1.000, precision 1.000,
  0.015 A, 1.000. The sub-suite `nanoparticle atoms` of `lib_single` that
  the first version of this note named no longer exists (moved 2026-09-25).

## 1. The model

    rho(r) = b(r) + sum_i A_i exp(-4 pi^2 |r - r_i|^2 / B_i)

`b` is a smooth background. Each atom is one isotropic Gaussian with its own
position, amplitude and width. Its integrated intensity is

    I_i = A_i (B_i / 4 pi)^(3/2)

The atoms form a mixture of `K` classes. Class `k` has an intensity `I_k`,
and an atom of class `k` has `I_i` near `I_k`. The width `B_i` belongs to the
atom and not to the class: it is free per atom, and it is expected to vary
with radius.

The reason for this split is physical. In an incoherent image the integral of
an atom is set by its scattering power, and displacement, the probe and the
reconstruction only spread it. So `I` carries the species and `B` carries the
environment. A statistic that mixes the two, such as the peak height, a
correlation with a fixed shape, or the amplitude at a fixed width, cannot
tell a light atom from a heavy atom that is smeared. Section 8.4 shows that
this is not a small effect: with widths that grow towards the surface, the
amplitude at a fixed width invents a species in a single-species particle.

Unknowns found from the map: the number of atoms, `r_i`, `B_i`, the number of
classes `K`, the class of each atom with its probability, and the `I_k`.

## 2. Length scales from the data

The element currently supplies an atom radius `r` and, through the lattice
table, a neighbour cutoff. Both follow from one measured number, the
nearest-neighbour distance `d_NN`: the median over atoms of the distance to
the nearest other centre, taken from the centres of the first binarisation.

- `r = 0.4 d_NN`. For Pt this gives 1.106 A against the tabulated 1.10 A, so
  `split_atoms`, `validate_atoms` and `discard_atoms` keep their behaviour.
- Neighbour cutoff `1.267 d_NN`, which is what `find_rMax` returns for FCC Pt
  (3.504 A). This assumes close packing. The general form is the midpoint
  between the first two peaks of the pair-distance histogram plus `0.15 r`;
  not emulated (Q5).
- The contact-score ceiling of `discard_atoms` (12 or 4, from the crystal
  system) is 12 in the prototype, the close-packed value that goes with the
  neighbour cutoff above (Q5), so that the pruning logic is unchanged.
  Deriving it from the data (the median contact score of the inner 30% of
  the atoms) belongs to 10.5.
- Fit and aperture radius `d_NN / 2`.
- Radial coordinate of an atom: its distance from the centre of the atom
  positions (`cendist`, already in `atom_stats`). For a particle that is far
  from spherical the depth below the surface is the better coordinate (Q3).

The distance tests of `split_atom` multiply squared voxel distances by `smpd`
where `smpd^2` is meant (`simple_nanoparticle.f90:914-969`), so the exclusion
radius written as `0.9 r` acts as 0.59 A for Pt at 0.358 A. The species-free
path should state it in Angstroms as a fraction of `d_NN`: 0.21 reproduces
today's effective value, 0.36 the written one. To be settled on the Pt test.

## 3. Pipeline

Nine parts inside `detect_atoms`; 3.9 gives the order they execute in (level
1, its pruning, then the rest). 3.1 and 3.9 are the present flow, with the
Gaussian kernel and the `d_NN` length scales substituted only when `element`
is absent. The rest is new and runs only under `discover_species=yes`.

### 3.1 Level 1: the present flow with a pseudo-atom

`phasecorr_one_atom`, the Otsu threshold ladder, the bisection,
`discard_small_ccs` and `split_atoms` run as now, with one change: the
template and the equal-atom simulation of `t2c` use the Gaussian of section 1
at a reference width `B_ref`, not an element. At 0.358 A the Pt template
correlates at 0.997 with a single Gaussian of sigma 0.42 A (`B` = 13.9 A^2),
which is the starting `B_ref`; with an element, `B_ref` is the Gaussian
fitted to that element's template. After the first fit (3.4) the median `B`
of the strong atoms is compared with `B_ref` and the ratio is reported
(`_species.txt`); re-running the detection at the fitted width is left to
10.5, because in the prototype level 1 and its products must not depend on
the opt-in key. Nothing is carried between calls of `detect_atoms`.

Level 1 finds the strongest class and nothing else. That is the behaviour to
preserve for single-species data, and the reason more levels are needed: the
Otsu gate sits at about a quarter of the strong peak, and in every
multi-species case of section 8 level 1 found the strongest class (all of it
in all cases but two, 395 and 394 of 396 there) and none of the weakest.

### 3.2 Noise reference

The noise region is the set of voxels inside the shipped sphere, at least the
soft-edge width from its rim, and farther than `1.5 d_NN` from every accepted
atom. It is rebuilt after each level. Three numbers come from it:

- the mean and standard deviation of each correlation map, which turn that
  map into a z map and remove the global norm of `phase_corr`;
- `sigma_n`, the standard deviation of the map itself;
- `s_A`, the standard deviation of the matched amplitude at `B_ref`, and
  `s_I`, that of the aperture intensity (3.5), both measured at 400 or more
  phantom sites placed at random in the region.

Measuring the noise of each statistic by applying that statistic to empty
positions makes no assumption about the colour of the noise. The emulation
used white noise, noise band-limited at 1 A, and strongly correlated noise,
with the same thresholds.

If the region is too small, because `mskdiam` hugs the particle, the program
stops and asks for a larger mask. The noise in the region is a lower bound on
the noise inside the particle, where alignment error adds to it.

**Half maps as the noise reference.** With `vol_even` and `vol_odd` given,
`sigma_n` and the mean and standard deviation that standardise each
correlation map come from the half-map difference (its noise is sqrt(2)
times that of a half and twice that of the average, so the full-map value is
the difference's divided by two), measured over the same region. The
difference is free of any density, so atoms that have not been found yet
cannot enter it. That matters: the region excludes 1.5 `d_NN` around every
accepted atom, which removes an undetected surface shell (its atoms sit
within one `d_NN` of the strong atoms) but not an undetected domain made only
of the weak species and thicker than that. In the check of section 8.8 such a
domain inside a mask of radius 14 to 16 A raised the region's standard
deviation by 36 to 40% (5% at 21 A); a robust spread removed only a quarter
to two fifths of that excess, because the matched filter spreads every atom
over a hundred voxels. The
direction of the error is conservative (a higher effective threshold, fewer
detections), and with a single map and no half maps it is accepted and
reported as the ratio of the region's robust to plain spread.

**The detection threshold is calibrated, not fixed.** A threshold on a z map
of local maxima is a multiple-testing problem: the number of noise maxima
above `k` grows with the volume searched and depends on the colour of the
noise, and the per-voxel Gaussian tail does not give it (it overpredicts the
count of maxima six times in the check of section 8.8). The region does. The
number of maxima of a smooth Gaussian field above `u` has the form
`C (u^2 - 1) exp(-u^2 / 2)` (the Euler characteristic density in three
dimensions), with `C` proportional to the volume searched. `C` is fitted to
the counts of region maxima above 2.5, 3.0 and 3.5, scaled by the ratio of
the search volume to the region volume, and `k` is the `u` at which the
expected number of false atoms in the search volume is a target of 0.2 per
particle. Scaling the region's own count by volume predicted the search
volume's count to within a factor 0.8 to 1.0 in the check; the calibrated `k`
came out at 5.0 for white and 1 A band-limited noise and 4.7 for strongly
correlated noise, so the fixed `k` = 5 of the first version of this note was
right for the emulated geometry and would be wrong for a larger search volume
or a different noise colour. The expected false count at the chosen `k` is a
reported diagnostic (section 7).

### 3.3 Residual levels

After level 1, repeat:

1. Subtract the fitted density of all accepted atoms from the map as read.
2. Matched-filter the residual with the same template; standardise to a z
   map with 3.2.
3. Candidates are local maxima of the z map above a threshold `k`
   (non-maximum suppression within `0.35 d_NN`) with at least `NVOX_THRESH` of
   the 27 surrounding voxels above `k`. The position is the centroid of the
   positive residual within `0.35 d_NN`.
4. Accept, fit each new atom on the residual, go to 1.

Two stages. Stage A uses the calibrated threshold `k_A` of 3.2 (5.0 on the
emulated geometry), searches the whole shipped sphere, and accepts any
candidate farther than `0.7 d_NN` from every accepted atom; it repeats until
it adds nothing. Stage B uses `k_B` = `k_A` - 1 and accepts a candidate only
if it passes a neighbour gate: no accepted atom closer than `0.85 d_NN`, and
at least three within `1.15 d_NN`. Its search volume is therefore the
accepted atoms dilated by `1.15 d_NN`. It also repeats until it adds nothing.
At most eight levels per stage. The gate is a switch (`nn_gate`, default
`yes`), and every atom records the stage and level that accepted it and its
detection z, so that a run with the gate off can be compared atom by atom.

Three choices here were forced by the emulation.

- Local maxima, not connected components. At a fixed low threshold the
  components of neighbouring atoms merge when the signal is strong, and the
  merged centroid lands between two atoms.
- No Otsu on the residual. It has no noise reference and can fall to the
  noise floor.
- The gate. Without it, `k_B` = 4 (its value on the emulated geometry,
  used for every number below) over the shipped sphere admits 2 to 4 noise
  maxima per particle (4.2 for white and band-limited noise, 2.3 for
  strongly correlated; 1.5 within the particle envelope; section 8.8,
  current variant of the levels; the first version of this note gave 6 to
  13 from an earlier variant), and in one single-species run of that earlier
  variant 13 of them were returned as a second species holding 2.4% of the
  atoms. That is how a species gets invented. With the gate, `k_B` = 4 gave
  0.2 to 0.4 false atoms per particle and raised the recall of weak, broad
  atoms from the 25% of stage A alone to 61% on a weak surface shell and 58%
  for weak adatoms, and from 59% to 78% on a weak-only domain.

**What the gate assumes and what it costs.** It is a local lattice prior: it
needs three found neighbours in a shell of 0.85 to 1.15 `d_NN`, so it
assumes a locally close-packed neighbourhood whose spacing is within 15% of
the particle's median. It uses only `d_NN`, so it does not assume one
lattice across the particle, but an amorphous or heavily strained surface
would fail it, and nothing in the emulation tests that. It cannot find a
weak atom with fewer than three accepted neighbours and a z below `k_A`; an
isolated adatom is found by stage A or not at all. A region made only of the
weak species is reached from its boundary inwards, one shell per level, and
in the check of section 8.8 the interior of such a domain was found to 99%.
Measured at the same `k`, the gate costs recall against an ungated search:
about six points on the surface shell and nine in the coordination 5 to 7
bin of the weak-only domain. Measured at the same false-positive rate, which
is the comparison that counts, the gated stage B at `k_B` = 4 beat an
ungated search at `k` = 4.5 by 7 to 15 points in every case and every
unsaturated coordination bin (the interior bin is 99% for both; 61 against
46% on the surface shell, 78
against 71% on the weak-only domain, 58 against 43% for weak adatoms of
coordination 2 to 4). The gate is kept, as a switch, on that evidence.

This is the one step that matches against a shape, and section 7 treats what
that means at the surface.

### 3.4 Density fit, stage 1: free amplitude

Three sweeps over all atoms. For atom `i` the data are the map minus the
background minus the fitted density of every other atom, inside a sphere of
radius `d_NN / 2`. The amplitude is linear; `B_i` is found by a
one-dimensional search that minimises

    SSE(B) / sigma_n^2 + ((ln B - m_i) / tau_i)^2

`SSE` is the sum of squared residuals over the sphere. `m_i` and `tau_i` are the median and the robust spread (floor 0.1) of `ln B`
over the strong atoms (amplitude above `10 s_A`) of the same radial shell,
within `d_NN` in radius; the global values are used when fewer than five
qualify. After each sweep the background is the residual smoothed with a
Gaussian of standard deviation `d_NN`.

A strong atom keeps its own width, since its data term dominates. A weak atom
cannot: with a free amplitude its width is not determined by its own density
(in a comparison of single-atom estimators made while drafting, not
reported here, the fitted sigma of an atom at 0.175 of
the strong amplitude ranged from 0.26 to 0.67 A around a true 0.35), so at
this stage it borrows the width of the strong atoms at its radius. Stage 2
removes the need for that.

Subtracting the neighbours is what allows the sphere to reach `d_NN / 2`
without a per-atom background term. An earlier variant, with the sphere of
`calc_isotropic_disp` (three quarters of that radius) and a free background
per atom, lost part of the integral of broad atoms to the background.

Positions were held fixed in the emulation. In the implementation each sweep
should also move the centre to the centroid of the atom's own residual.

### 3.5 Aperture intensity

    I_i = smpd^3 sum over |r - r_i| < d_NN/2 of [ map - b - sum_{j /= i} rho_j ] / F(d_NN / 2 sigma_i)

`F(t) = erf(t / sqrt 2) - sqrt(2 / pi) t exp(-t^2 / 2)` is the fraction of a
Gaussian inside `t` standard deviations, and `sigma_i^2 = B_i / 8 pi^2`. The
correction is 0.2% for the narrowest atoms of the emulation and 6% at twice
their variance. The factor `smpd^3` puts the voxel sum in the units of the
continuous integral `A_i (B_i / 4 pi)^(3/2)` of section 1, which the
stage-2 penalty of 3.7 compares it with; the first version of this note
omitted it (the mixture of 3.6 is scale-free, so the emulation did not see
the error, but the tie of 3.7 would have been off by `smpd^-3`, about 22 at
0.358 A). Every intensity in this note, in the products and in the tests is
in those units: the voxel sum of a pseudo-atom of intensity `q` is `q /
smpd^3` (section 10, phase 1).

This is the classification statistic. It is a plain sum, so it involves no
matching; it measures the quantity that carries the species; and it depends
on `B_i` only through the small correction. Its cost is noise: about 2.3
times that of the matched amplitude for white noise, and more when the noise
is concentrated at low frequency (case 8b of section 8.1).

### 3.6 Species discovery

Input: the `I_i` and `s_I`. A one-dimensional Gaussian mixture is fitted by
expectation-maximisation for `K` = 1, 2, 3, with every component variance
floored at `s_I^2`, from two deterministic starts (sorted quantiles, as in
`sortmeans`, and the `K-1` largest gaps of the sorted sample). Each fit gets
the Bayesian information criterion `BIC = -2 ln L + (3K - 1) ln N`, with `L`
the likelihood and `N` the number of atoms.

A `K` is admissible when every class holds at least `max(8, 0.02 N)` atoms and
adjacent classes are separated by `D >= 3`, with
`D = |mu_a - mu_b| / sqrt((var_a + var_b) / 2)`. `K` = 1 is always admissible.
The admissible `K` with the lowest BIC is taken. Each atom gets the posterior
probability of every class, which is stored, and is labelled by the largest.
Classes are numbered by decreasing intensity.

The labels use no spatial information. Whether the species are mixed,
ordered or segregated is the result, so no neighbour may influence a label.

`nspecies` fixes `K` when the user knows it, and a species list does the same
(section 5).

### 3.7 Density fit, stage 2: the B factor inside the mixture

With the classes known, every atom is refitted with its integral tied to its
class. For atom `i` of class `k` the one-dimensional search over `B` now
minimises

    SSE(A, B) / sigma_n^2 + ((A (B / 4 pi)^(3/2) - I_k) / tau_k)^2

with the amplitude given in closed form for each `B`. `tau_k` is the
intrinsic spread of the class, `sqrt(var_k - s_I^2)` with a floor of
`0.05 I_k`. Two rounds: fit, recompute the aperture intensities and the
classes, fit again.

The tie makes the width of a weak atom measurable. Its peak height is now
informative about its width, because a broader atom of the same species must
be lower by the cube of the width ratio. No width is borrowed from the strong
class any more, so a light species that is broader than the heavy one at the
same radius comes out broader (case R6 in section 8.2).

The cost is an assumption: all atoms of a class have the same integral to
within `tau_k`. An atom that is present part of the time will be fitted as a
broader atom of its class. The stage-1 widths are kept in the per-atom table
so that the two can be compared.

### 3.8 Radial profiles

The per-atom table (section 6) holds every fit parameter with the atom's
radius and coordination number, so any dependence can be plotted directly.
`detect_atoms` also writes the profiles: per class, in equal-count radial
shells of at least 20 atoms, the number of atoms, the mean radius, the width
(root mean square sigma and its standard error, and the same as `B`), the
mean intensity with its standard error, the mean coordination, and the
predicted detection signal-to-noise of section 7.

### 3.9 Pruning and validation

The order is: prune level 1, then fit it provisionally, then recover, then
prune the recovered atoms on the atom table, then fit and classify. No
per-atom fit survives a `discard_atoms` pass, so fitting starts after
level-1 pruning; the recovered atoms are pruned without `discard_atoms`,
which keeps their fits.

Level 1 is pruned as today: `discard_small_ccs`, `split_atoms`,
`validate_atoms` and `discard_atoms` run on the level-1 set, unchanged in
logic, with the length scales of section 2 in place of the element's. This
is the production path and its products (`_ATMS.pdb`, `_BIN`, `_CC`, `_MSK`,
`_SIM`) are written from it exactly as now. The first version of this note
placed pruning after the fits; the constraint of section 0 rules that out,
since `discard_atoms` rebuilds `atominfo` on every pass.

Residual-level atoms have no component at the level-1 threshold, so the
present pruning cannot apply to them and must not be imitated with
model-generated components (the size test would then be circular). They
are pruned on the atom table alone: the contact-score rule of
`discard_atoms` (fewer contacts than the threshold among the 15% of atoms
farthest from the centre), computed on the merged set, plus the exclusion
distances of 3.3, in a routine that touches neither `_BIN` nor `_CC`.

The stage-1 and stage-2 fits (3.4, 3.7), the aperture intensities (3.5) and
the species call (3.6) then run on the merged, pruned set. `valid_corr` of
the level-1 atoms is the one measured today against the equal-atom
simulation; for the recovered atoms the same correlation against a
simulation in which every atom has its class intensity and its own width is
written to the species table as a diagnostic. Once recovered atoms are
promoted into `_ATMS.pdb` (section 10.5), `valid_corr` moves to that
simulation for every atom, and two penalties of the present equal-atom,
equal-width simulation disappear: the one on a weak atom next to strong
neighbours, and the one on a broad atom at the surface.

## 4. Simulation

`atoms%convolve` gets a pseudo-atom branch: for a sentinel `Z` the kernel is
the Gaussian of section 1 with `B` from `beta` and amplitude
`q (4 pi / B)^(3/2)`, `q` from `occupancy`. An atom with `q` = 1 has unit
integrated intensity whatever its width.

At promotion (10.5; in the prototype `_SIM.mrc` is the present equal-atom
simulation) `_SIM.mrc` is rendered with the class intensity for `q` and the
fitted `B_i` of each atom. Against the noise-free generating density of the multi-species
cases this correlated at 0.99 or better, the same as rendering every atom
with its own fitted amplitude, so the class constraint costs nothing. One
width per class gave 0.94 to 0.99, and the equal-atom simulation used today
0.86 to 0.96 (section 8.1).

Because `occupancy` and `beta` are PDB columns, the PDB defines the
simulation. `simulate_nanoparticle pdbfile=` then is the species-free map
simulator, with per-atom amplitude and variance, and also the fixture
generator for the tests of section 10.

## 5. When the species are given

**The maximum over per-species correlation maps.** `phase_corr` removes the
scale of the template, and at this sampling the shapes of the element
templates are nearly the same: the Pt and Ni templates correlate at 0.998. In
the emulation the correlation map made with the Pt template and the one made
with the Ni template correlated at 0.9955, and the maximum of the two found
the same atoms as either: all of the strong class and, at an intensity ratio
of 6:1, none of the weak class (section 8.5). The weak atoms are missed
because the correlation map is linear in atom amplitude and the first
threshold sits above them. A template of the weak species has the same
shape, so it does not raise them. Residual levels do (3.3).

The same holds for widths: a template of the wrong width loses little. For an
atom with twice the variance of the template the matched signal falls by 8%.

**What a species list does buy.** Three things, none of them in the
segmentation: `K` is fixed, so model selection cannot go wrong; the classes
get element names, assigned by decreasing intensity to decreasing `Z`; and
`_ATMS.pdb` carries real element symbols. The class intensities and widths
still come from the data. They cannot come from the scattering-factor table,
which describes the electrostatic potential (Pt over Ni is 2.0 at the peak of
the templates and 1.7 in the integral) and not the HAADF signal.

**Proposal.** One implementation with two entries. Without `element` the
classes are discovered and named `X1`, `X2`, ... With a species list (the
existing `element` key already accepts two-element strings such as `CdSe`)
`K` is the length of the list and the classes take its names. The species
list belongs to promotion (10.5). In the prototype the discovery itself is
opt-in (`discover_species=yes`, section 9), `K` is fixed only through
`nspecies`, and leaving out `element` alone gives the present flow with the
Gaussian template and the `d_NN` length scales, and nothing more. Today a
two-element string only selects its first element for the template and for
every simulated atom (`simple_nanoparticle.f90:262-264`), so giving all
species has no effect on detection at present.

**The mixture model with B factors** is sections 3.6 and 3.7: class
intensities, per-atom class probabilities, and a width per atom estimated
inside the class model. A fuller form would compute the class probabilities
from the density likelihood of each atom under each class with `B` profiled
out, replacing the aperture sum. It uses the same model and the same code
path as stage 2; it was not emulated, and the aperture sum has the advantage
that it depends on the atom's own fit only through the small correction of
3.5.

## 6. Products

In the prototype the present products are written by the present path from
the level-1 atoms and are not changed: `_ATMS.pdb` (B column `valid_corr`,
as today), `_BIN.mrc`, `_CC.mrc`, `_MSK.mrc` and `_SIM.mrc` (equal atoms of
the reference width). Everything the prototype finds goes into three new
files, so that a run with `discover_species=yes` differs from one without it
only by their presence:

- `_species.csv`, one row per atom at full precision, level-1 and recovered
  atoms alike: position, radius, coordination number, stage and level that
  accepted it, detection z, amplitude, `B` of stage 1 and of stage 2,
  aperture intensity and its signal-to-noise, class, class probabilities,
  `valid_corr`, and whether the gate was needed to accept it.
- `_species_radial.csv`: the profiles of 3.8.
- `_species.txt`: `K`; class intensities, ratios and fractions; `D`; the BIC
  table with admissibility; the expected misclassification rate; `d_NN`,
  `B_ref`, the noise numbers and their source (half maps or region); the
  calibrated thresholds and the expected false count at each; atoms added
  per stage and level; the diagnostics of section 7.

Promotion, after the validation of section 10.5, adds the recovered atoms
and the class information to the present products:

- `_ATMS.pdb`. Element column: the class symbol (`X1`, `X2`, `X3`, registered
  in `get_element_Z_and_radius` with sentinel `Z`, precedent `CDSE`, 999; or
  the element names when a species list is given). Chain: `A`, `B`, `C` by
  class, so a class is one selection in a viewer. Occupancy: `I_k / I_1`.
  B column: `B_i`. Today the B column holds `valid_corr`
  (`simple_nanoparticle.f90:1821`); within SIMPLE only the unit tester reads
  it back, and `valid_corr_in_bfac_field.pdb` keeps it for `atoms_stats`.
- `_SIM.mrc` as in section 4.
- `_CC.mrc`, `_BIN.mrc`, `_MSK.mrc`. Atoms found at residual levels have no
  connected component of the level-1 threshold. Their components are derived
  from the data, never from the model: the connected region of the residual
  z map above `k_B / 2` around the accepted maximum, split between
  neighbouring maxima by a watershed on the z map (the merging of components
  at a low threshold is what moved detection to local maxima, so the split
  is required, not optional). The first version of this note generated the
  component from the fitted Gaussian; that would make the size test of
  `discard_atoms` circular and change what `atoms_stats` measures, and is
  withdrawn. `atoms_stats` reads `_CC.mrc` with `_ATMS.pdb`, so promoted
  atoms need a component. Not emulated.

## 7. Diagnostics and limits

- **Core and surface.** For an atom of integrated intensity `I` and width
  sigma the best possible detection signal-to-noise scales as
  `I sigma^(-3/2)`. An atom with twice the variance is detected at 0.59 of the
  signal-to-noise of the same atom in the core. That loss belongs to the data
  and no method avoids it; fitting the atom needs it to be found first. Using
  a core-width template on it costs a further 8%. A bank of two template
  widths was tried and dropped: it changed the number of atoms found by a few
  in either direction and admitted more noise peaks. So detection is less
  complete at the surface, for the light species first.
  The design does three things about it. The species call of an atom that
  was found does not depend on detection sensitivity, and neither does its
  fitted width. The predicted detection signal-to-noise is reported per class
  and radial shell, computed from the class intensity and the fitted width of
  that shell. And the shells where it falls below about 6.5 are flagged:
  there the emulation found 56-81% of the light atoms against 88-100% above
  (section 8.2), and the atoms that were found were the narrower ones, so the
  shell average of the light class is biased low in width.
- **Detection limit overall.** In the emulation a class whose found atoms had
  a median fitted amplitude of 6.6 or more times `s_A` was found to 88-100%;
  at 6.4 it was 84%, at 5.0 it was 65% and at 4.1 it was 45%. The same number
  can be read off an existing map before any code is written: the peak of the
  correlation map at the weakest visible atoms over its standard deviation
  away from the particle.
- **Missing species.** When a class is below the limit the program reports
  one class and cannot know better from intensities. The coordination of
  interior atoms can: with a quarter of the sites undetected the inner 30% of
  the atoms had a mean coordination of 8.9 where 12 is expected (case 8b).
  Reported as the interior coordination deficit.
- **What a class is.** An intensity class. A site that is occupied part of
  the time, or an atom that moved during acquisition, integrates to less than
  its species and can form a class of its own. The per-class coordination and
  radial profile are reported so that a class confined to the surface is
  visible. It is a flag and not a verdict, since surface segregation of the
  light species gives the same picture.
- **Marginal atoms.** Stage B occasionally admits a noise peak at a plausible
  surface site (section 8.6 gives the count). Each atom carries its detection
  z, its stage and level, and whether the gate was needed to accept it.
- **Expected false atoms.** The calibration of 3.2 gives, for each stage,
  the expected number of noise maxima above its threshold in its search
  volume. Both numbers are reported with the thresholds and the noise source
  (half-map difference or region), together with the ratio of the region's
  robust to plain spread, which flags a contaminated region.
- **Half-map agreement.** With `vol_even` and `vol_odd` given, the aperture
  intensities are measured at the final positions in each half map, and the
  label agreement and the intensity correlation are reported. The halves
  share their alignment, so this measures the effect of noise on the labels
  and not model bias. Not emulated.
- **Model bias.** Labels are found from scratch at every call; nothing is
  carried between iterations. Since `_SIM.mrc` only aligns (section 0), a
  wrong label can reach the next map only through the orientations. A direct
  test is cheap and worth having once `_SIM.mrc` carries the labels (after
  promotion, 10.5): render a random 10% of the atoms at the pooled intensity
  and check that they separate as well as the rest in the next iteration.
  Not emulated.

## 8. Emulation

Particle: FCC, a = 3.84 A, diameter 24 A, 528 atoms, box 128 at 0.358 A.
Gaussian atoms, sigma 0.35 A in the core with 5% scatter, integrated
intensity 1 for the strongest class. Two models of broadening, both at
conserved integral: "surface x2" doubles the variance of the 195 atoms with
coordination 9 or less; "radial x2" (or x3) lets the variance grow with the
square of the radius to twice (three times) its central value at the surface.
Noise is given as its standard deviation per voxel over the peak of a core
atom of the strongest class. One particle per case unless stated.

### 8.1 Detection, discovery and simulation

"Level 1" is the present flow. Counts are atoms found per class, strongest
first. "SNR" is the median fitted amplitude over `s_A` per class. The last
column is the correlation with the noise-free generating density of the
equal-atom simulation used today and of the proposed one.

| # | Case | Level 1 | All levels | False | K | Class intensities | SNR | Simulation: today, proposed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3:1 atoms, intensity 3:1, noise 0.03 | 396, 0 | 396/396, 132/132 | 1 | 2 | 1 : 0.330 | 79, 27 | 0.937, 0.998 |
| 2 | 3:1 atoms, intensity 6:1, noise 0.03 | 396, 0 | 396/396, 132/132 | 0 | 2 | 1 : 0.167 | 76, 12 | 0.903, 0.999 |
| 3 | 6:1 atoms, intensity 3:1, noise 0.05 | 453, 0 | 453/453, 75/75 | 0 | 2 | 1 : 0.336 | 46, 15 | 0.960, 0.997 |
| 4 | 6:1 atoms, intensity 6:1, noise 0.05 | 453, 0 | 453/453, 73/75 | 0 | 2 | 1 : 0.161 | 46, 7.8 | 0.943, 0.998 |
| 5 | as 2, noise 0.05, surface x2 | 396, 0 | 396/396, 116/132 | 0 | 2 | 1 : 0.169 | 42, 7.0 | 0.871, 0.997 |
| 6 | as 5, light species at the surface | 396, 0 | 396/396, 60/132 | 0 | 2 | 1 : 0.155 | 47, 4.1 | 0.918, 0.997 |
| 7 | as 5, light species in the core | 396, 0 | 396/396, 130/132 | 0 | 2 | 1 : 0.159 | 37, 8.6 | 0.862, 0.997 |
| 8 | as 2, noise 0.05 band-limited at 1 A | 396, 0 | 396/396, 132/132 | 0 | 2 | 1 : 0.169 | 49, 7.9 | 0.901, 0.998 |
| 8b | as 2, strongly correlated noise 0.05 | 395, 0 | 396/396, 1/132 | 0 | 1 | - | 12, 3.1 | 0.981, 0.984 |
| 9 | single species, noise 0.05 | 528 | 528/528 | 0 | 1 | - | 41 | 0.994, 0.998 |
| 10 | single species, noise 0.05, surface x2 | 528 | 528/528 | 0 | 1 | - | 45 | 0.952, 0.998 |
| 11 | as 10, noise band-limited at 1 A | 528 | 528/528 | 0 | 1 | - | 41 | 0.948, 0.998 |
| 12 | three species 2:1:1, intensities 1 : 0.5 : 0.167, noise 0.03 | 264, 132, 0 | all 528 | 0 | 3 | 1 : 0.495 : 0.166 | 79, 40, 13 | 0.877, 0.999 |
| R1 | as 2, noise 0.05, radial x2 | 396, 0 | 396/396, 111/132 | 0 | 2 | 1 : 0.173 | 33, 6.4 | 0.896, 0.996 |
| R2 | as R1, light species at the surface | 396, 0 | 396/396, 86/132 | 0 | 2 | 1 : 0.166 | 34, 5.0 | 0.909, 0.996 |
| R3 | single species, noise 0.05, radial x2 | 528 | 528/528 | 0 | 1 | - | 32 | 0.976, 0.997 |
| R4 | 3:1 atoms, intensity 3:1, noise 0.05, radial x3 | 394, 0 | 396/396, 128/132 | 1 | 2 | 1 : 0.314 | 25, 8.5 | 0.895, 0.995 |
| R5 | 6:1 atoms, intensity 6:1, noise 0.03, radial x2 | 453, 0 | 453/453, 74/75 | 0 | 2 | 1 : 0.170 | 52, 8.4 | 0.928, 0.999 |
| R6 | as R1 at noise 0.03, light species 1.3x the variance of the heavy one at every radius | 396, 0 | 396/396, 129/132 | 0 | 2 | 1 : 0.163 | 52, 6.6 | 0.881, 0.998 |
| R7 | single species, noise 0.05, radial x3 | 528 | 528/528 | 0 | 1 | - | 26 | 0.954, 0.996 |

Every atom that was found and matched to the generating model received the
right label, in every case. `D` was 10.1 or more wherever more than one class
was found. Class fractions equalled the fractions among the found atoms;
where weak atoms were missed (cases 5, 6, R1, R2) the light fraction is
underestimated accordingly (0.227, 0.132, 0.219 and 0.178 against 0.25).

In the multi-species cases the simulation with one width per class
correlated at 0.941-0.994 and the one with per-atom amplitudes at
0.995-0.999, against 0.995-0.999 for the proposed one.

### 8.2 Radial profiles

Case R6, where the light species is broader than the heavy one at every
radius. Widths are root-mean-square sigma per shell in Angstroms: fitted with
its standard error, then the generating value. "Stage 1" is the width before
the class tie (3.4).

| Class | Atoms | Mean radius (A) | Stage 1 | Stage 2 | Generating | Intensity / I_1 | Predicted SNR | Found in shell |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 80 | 5.3 | 0.383 | 0.384 +- 0.002 | 0.384 | 0.996 | 66 | 100% |
| 1 | 79 | 8.0 | 0.417 | 0.417 +- 0.003 | 0.417 | 1.002 | 59 | 100% |
| 1 | 79 | 9.6 | 0.448 | 0.448 +- 0.002 | 0.448 | 1.003 | 53 | 100% |
| 1 | 79 | 10.8 | 0.473 | 0.475 +- 0.003 | 0.477 | 1.001 | 48 | 100% |
| 1 | 79 | 11.7 | 0.489 | 0.492 +- 0.003 | 0.492 | 0.998 | 46 | 100% |
| 2 | 26 | 5.7 | 0.416 | 0.434 +- 0.007 | 0.438 | 0.174 | 9.0 | 100% |
| 2 | 26 | 8.4 | 0.447 | 0.475 +- 0.010 | 0.477 | 0.162 | 7.8 | 100% |
| 2 | 26 | 9.7 | 0.466 | 0.498 +- 0.007 | 0.516 | 0.161 | 7.3 | 100% |
| 2 | 26 | 10.7 | 0.466 | 0.504 +- 0.008 | 0.532 | 0.157 | 7.1 | 96% |
| 2 | 25 | 11.7 | 0.479 | 0.528 +- 0.008 | 0.548 | 0.161 | 6.6 | 92% |

What the profiles of all radial cases show:

- **Strongest class.** The fitted shell width was within 0.006 A (1.1%) of
  the generating value in every shell of every case, including the single
  species cases R3 and R7 and the factor-three gradient of R4 and R7.
- **Light class.** With stage 2 its own profile is recovered: in R6 it is
  broader than the heavy class in every shell, by 6% to 14% in sigma at
  similar radii where 14% was generated. Stage 1 alone pulls it to the heavy
  profile (first width column). Over all radial cases the fitted width of the
  light class was within 3% of the generating value in the inner shells, up
  to 6% low in outer shells where detection was complete or nearly so, and up
  to 10% low where it was not (R2, 56-70% found): there the atoms that are
  found are the narrower ones.
- **Intensity.** The class intensity is flat in radius for the strongest
  class (0.987 to 1.008). For the light class it is flat within its errors in
  the inner shells and falls by up to 10% in the outermost ones (0.157-0.161
  against 0.174 here).
- **Detection.** Over the light-class shells of R1, R2, R4, R5 and R6, a
  predicted signal-to-noise of 6.6 or more went with 88-100% of the atoms
  found, and 6.1 or less with 56-81%.

### 8.3 What the cases show

- Cases 1-4: both readings of "3-6:1", and both together, are within reach at
  these noise levels.
- Cases 5-7 and R1-R6: widths that grow towards the surface do not confuse
  the labels, wherever the light species sits, and the widths are recovered
  per class. They do limit detection: a light atom that is also broad is the
  hardest object (cases 6 and R2).
- Case 8: the thresholds hold when the noise is not white.
- Case 8b: when the light species is below the detection limit the answer is
  one class, with a coordination deficit. The aperture statistic is also at
  its worst here, because noise at low frequency adds coherently over the
  sphere.
- Cases 9-11, R3 and R7: no species is invented from broadening or from
  noise.
- Case 12: `K` is not tied to two.

### 8.4 Choice of the classification statistic

Three statistics were run through the same class discovery.

| Cases | Amplitude at fixed width | Integral from a single-atom fit with a width prior | Aperture intensity (3.5) |
| --- | --- | --- | --- |
| 1-4, 8, 12, R2 | right K, all labels right | right K, all labels right | right K, all labels right |
| R1, R5, R6 | right K, 98.9-99.2% | right K, 99.8-100% | right K, 100% |
| 5 | K = 3 | K = 2, 99.8% | K = 2, 100% |
| 6 | K = 3 | K = 2, 85.7% | K = 2, 100% |
| 7 | K = 3 | K = 1 | K = 2, 100% |
| R4 | K = 1: the light species is not seen | K = 1 | K = 2, 100% |
| 10, single species | K = 2: 37% of the atoms called a second species | K = 2 | K = 1 |
| 11, single species | K = 2 | K = 1 | K = 1 |
| 9, R3, R7, single species | K = 1 | K = 1 | K = 1 |
| 8b | K = 1: the light species was not detected | K = 1 | K = 1 |

The amplitude at a fixed width has the lowest noise and is the wrong
statistic: it reads a broad heavy atom as a lighter species. An earlier
review of the first draft (2026-10-07, before the one recorded in section
15) recommended it on the evidence of particles with uniform widths; this
table withdraws that recommendation.

### 8.5 Per-species correlation maps

The same synthetic maps, correlated with the Pt template, with the Ni
template, and the voxel-wise maximum of the two. "Gate" is the Otsu threshold
of the present flow.

| Case | Correlation of the two maps | Weak peak / strong peak | Gate / strong peak | Weak atoms above the gate: Pt map, Ni map, maximum |
| --- | --- | --- | --- | --- |
| 1 (intensity 3:1) | 0.9955 | 0.34 | 0.25 | 132, 132, 132 of 132 |
| 2 (intensity 6:1) | 0.9957 | 0.17 | 0.25 | 0, 0, 0 of 132 |
| 4 (intensity 6:1) | 0.9955 | 0.17 | 0.25 | 0, 0, 0 of 75 |

The three maps behave as one. At 3:1 the weak atoms clear the gate in all of
them, and the threshold criterion of the present flow then drops them: level 1
of case 1 returned none.

### 8.6 Repeats

Five seeds each for cases 2, 4, 6, 9, 10, R1, R2 and R3 (40 runs): the right
`K` in 40 of 40 and every label right. In the 15 single-species runs `K` was 1
throughout and one noise peak was accepted in one run. Recall of the light
species ranged over 73-75 of 75 (case 4), 60-74 of 132 (case 6), 94-107 of 132
(R1) and 75-89 of 132 (R2), and its intensity estimate over 0.162-0.187 where
0.167 was generated.

False atoms over all 60 runs of sections 8.1 and 8.6: five, one each in five
runs. The one that was examined (case 1) was accepted at stage B, at a surface
site with four accepted neighbours.

### 8.7 Not covered

`split_atoms`, `discard_atoms`, `validate_atoms`, the soft mask and the
iterated-conditional-modes denoising that a map may have been through
(`icm_imgfile`) were not emulated. Atoms are exact Gaussians, the background is zero,
positions are not refined, the particle is a sphere, and the noise is
stationary. The noise levels are guesses; what transfers is the
signal-to-noise scale of section 7. `D` = 3, the class-size floor and the 5%
floor on `tau_k` are starting values; the detection thresholds are
calibrated per map (3.2), and the target of 0.2 false atoms per particle is
a starting value.

### 8.8 Checks of 2026-10-08

A second implementation of the generator of section 8 and of the detector of
3.3 (`codex_checks.py`, numpy only, written without reference to
`emul_species.py`) was run to test the review's points about the thresholds
and the gate. Same box, sampling, template, lattice and diameter; its own cut
of the sphere gave 531 atoms, 210 of them with coordination 9 or less (528
and 195 in the emulation of 8.1); the level-1 set was taken as the strong
class at its true positions jittered by 0.05 A and fitted on the map, as the
emulation found it in every case; stage-1 fits with positions fixed and no
background; six noise maps per colour, five seeds per case.

Noise maxima above `k` in a pure-noise z map, standardised with the region's
mean and standard deviation, per map (search over the shipped sphere of
radius 21 A, 845,000 voxels; over the particle envelope of radius 14.7 A,
291,000 voxels; and inside the noise region, 233,000 voxels):

| Noise | k = 3.5: sphere, envelope, region | k = 4 | k = 4.5 | k = 5 | Gaussian tail x sphere voxels at k = 4 |
| --- | --- | --- | --- | --- | --- |
| white | 29.0, 10.3, 7.2 | 4.2, 1.5, 1.2 | 0.33, 0, 0 | 0, 0, 0 | 26.8 |
| band-limited at 1 A | 27.3, 10.2, 6.7 | 4.2, 1.5, 1.3 | 0.67, 0.17, 0 | 0, 0, 0 | 26.8 |
| strongly correlated (2.5 A) | 11.3, 4.2, 2.3 | 2.3, 0.67, 0.83 | 0.33, 0, 0.33 | 0, 0, 0 | 26.8 |

The region's count scaled by volume predicted the sphere's count at `k` = 4
to within a factor 0.98, 0.86 and 0.77. The calibration of 3.2 gave `k` =
5.01, 5.00 and 4.73 for a target of 0.2 false atoms per particle, with 0,
0 and 0.17 observed at those thresholds. For the tester of phase 1: the
white-noise region counts above 2.5, 3.0 and 3.5 were 154.5, 38.0 and 7.2
per map, the search-to-region volume ratio 3.628, the fitted `C` 2105
(the Euler-characteristic value for white noise smoothed by the template,
`V / FWHM^3 (4 ln 2)^(3/2) / (2 pi)^2`, is 1656), and the resulting `k`
4.99. `detect_peak_thres_fdr` of
`simple_segmentation.f90` (the false-discovery-rate routine: median and
robust spread of the peak values as the null, Benjamini-Hochberg selection
bounding the expected fraction of false discoveries at q = 0.05), applied to
all local maxima of the sphere, selected 0.2 to 0.3 noise peaks per map: the distribution of maxima
(mean z 1.84, spread 0.73) has a heavier right tail than its median and
spread describe.

Recall of the weak atoms and false atoms per particle at noise 0.05 with
surface broadening, three cases: the surface shell of case 6 (132 weak atoms
at intensity 1/6, coordination 5 to 9); a weak-only domain (118 atoms with
x > 4 A at intensity 1/6, 57 of them interior); 40 weak adatoms at
intensity 1/6 on sites of coordination 2 to 4 outside a single-species
particle. "Calibrated" is the threshold of 3.2 (5.0 within 0.02 on these
maps).

| Variant | Shell: recall, false | Domain: recall, false | Adatoms: recall, false | Domain interior (coordination 10-12) | Domain surface (5-7) |
| --- | --- | --- | --- | --- | --- |
| stage A only, k = 5 | 25.5%, 0.2 | 58.6%, 0 | 24.5%, 0 | 93% | 26% |
| A at k = 5, B gated at k = 4 | 60.6%, 0.4 | 77.6%, 0.4 | 57.5%, 0.2 | 99% | 55% |
| A at k = 5, B ungated at k = 4 | 66.2%, 3.8 | 82.4%, 3.2 | 68.0%, 3.8 | 100% | 64% |
| A at k = 5, B ungated at k = 4.5 | 45.9%, 0.6 | 71.0%, 0.6 | 43.0%, 0.4 | 99% | 43% |
| calibrated A, gated B at k = 4 | 60.6%, 0.4 | 77.6%, 0.4 | 57.0%, 0.2 | 99% | 55% |
| `detect_peak_thres_fdr`, q = 0.05, ungated | 55.5%, 1.6 | 76.8%, 1.4 | 52.5%, 0.4 | 99% | 54% |
| A at k = 5, B gated at k = 3.5 | 72.3%, 0.8 | 84.9%, 0.8 | 70.0%, 0.2 | 99% | 70% |

The ungated search at `k` = 4.5 is the one with a false count close to the
gated stage B's; the gate wins against it everywhere. The FDR routine
controls what it promises (its false discoveries are 2% of its discoveries
here) but admits three to four times the gate's false atoms at lower recall.
The gated stage B at `k` = 3.5 shows the gate can be pushed; its false count
doubles and it was not adopted.

Contamination of the region by an undetected weak-only domain, with the
level-1 atoms excluded by 1.5 `d_NN`: standard deviation measured over true,
1.054 at a mask radius of 21 A (337,000 voxels, 29 undetected atoms inside
the region), 1.356 at 16 A (30,500 voxels), 1.398 at 14 A (12,100 voxels);
robust spread 1.019, 1.221, 1.298. With the surface shell undetected
instead: 1.000 to 1.002, no weak atom inside the region.

What these checks do not change: they use the same idealised generator as
section 8, so they say nothing about reconstructions, backgrounds, the
five-Gaussian kernels or position refinement. That is what section 10.5 is
for.

## 9. Code map

The prototype touches `detect_atoms` and the modules under it, and nothing
else that runs today. Every new behaviour sits behind a condition that is
false for a run as it is given today (`element` present, `discover_species`
absent), so the present path is the same code executing the same branches.

New: `src/main/nano/simple_nano_species.f90`, array numerics with no image
dependence: the mixture fit, BIC, admissibility and the choice of `K`
(`fit_species_mixture(x, s_meas, nspecies_in, K, labels, post, mu, var)`; the
routine is not named after the key), the enclosed-fraction function `F` of
3.5, and the threshold calibration of 3.2 from an array of region maxima
counts. Tester
`simple_nano_species_tester.f90`, sub-suite `species` of `unit_single`
(table `suites_single` in `simple_commanders_test_class.f90`, list in
`simple_test_ui_class.f90`, checked by `scripts/check_test_registry.py`).

`src/defs/simple_defs_atoms.f90`: the pseudo-atom symbols `X1`, `X2`, `X3`
in `get_element_Z_and_radius`, each with its own sentinel `Z` in a reserved
range held in named constants (precedent `CDSE`, 999; the range must not
collide with it), radius 1.0 until the data give `0.4 d_NN`.

`src/main/nano/simple_atoms.f90`: `convolve` gets the pseudo-atom branch of
section 4 for `Z` in the reserved range, reading `beta` as `B` and
`occupancy` as `q`; `element_exists`, `set_element`, `guess_an_element` and
`Z_and_radius_from_name` accept the symbols. The per-element table, `epot`
and every existing branch stay as they are. `simple_atoms_tester.f90`
gains the closed-form tests of phase 1.

`src/main/nano/simple_nanoparticle_utils.f90`: `phasecorr_one_atom` gets an
optional kernel width (`B_ref`) and, when it is given, correlates with the
pseudo-atom instead of the element; `est_nn_dist(centers)`; `find_rMax` and
`calc_contact_scores` get an optional `d_NN` that replaces the element
lookup. Callers that do not pass the new arguments behave as before.

`src/main/nano/simple_nanoparticle.f90`:

- a logical `l_species_free`, set in `new` when `params%element` is empty.
  `new` then does no `Z` lookup; `theoretical_radius`, the neighbour cutoff
  and the contact-score ceiling are set from `d_NN` once the first
  binarisation has produced centres (section 2), in `identify_atomic_pos`,
  before `discard_small_ccs`. With an element, `new` is unchanged.
- `atom_stats`: `amp`, `bfac`, `bfac_stage1`, `aper_int`, `aper_snr`,
  `det_z`, `det_stage`, `det_level`, `gate_used`, `species`, `species_post`.
  `cendist`, `u_iso` and `isocorr` exist.
- `identify_atomic_pos`: the present sequence, with the pseudo-atom template
  and the `d_NN` length scales substituted at the element entry points when
  `l_species_free` (`phasecorr_one_atom` at `:587`, `t2c` through
  `write_centers` and `atoms%convolve`, `calc_contact_scores` at `:1045`,
  `:1121`, `:1156`, `get_lattice_params` at `:1032-1033`, `simulate_atoms`
  at `:1775-1776`). After `discard_atoms`, and only when
  `params%l_discover_species`: `fit_atoms_joint` on the level-1 atoms
  (stage 1, the provisional fit that the subtraction of 3.3 needs),
  `est_noise_region`, `calibrate_threshold`, `detect_residual_levels`,
  `prune_recovered`, `fit_atoms_joint` (both stages, merged set),
  `calc_aperture_int`, `assign_species`, a diagnostic render of the
  class-intensity, per-width simulation into a scratch image through the
  pseudo-atom `convolve` (for the recovered atoms' `valid_corr` column; not
  written to disk), `write_species_report`, `write_radial_profiles`; then
  the present output block, unchanged, from the level-1 atoms. The
  recovered atoms live in a second atom table (`recovered(:)` of the same
  type) so that nothing that writes the present products can see them until
  promotion. `l_species_free` belongs to the nanoparticle (set in `new`);
  `l_discover_species` and `l_nn_gate` belong to `parameters`, derived from
  the keys in `simple_parameters_phases.f90`.
- `simulate_atoms`, `write_centers`, `kill`: unchanged in the prototype,
  except that `kill` frees the new table.

`calc_isotropic_disp` stays for `atoms_stats`. Its log-linear estimator over
positive voxels is biased in noise (at noise 0.03 it returned amplitude 0.91
and sigma 0.371 A for true 1 and 0.350); `atoms_stats` can take `u_iso` from
the stage-2 fit once that exists.

Parameters and UI, prototype: `element` becomes optional
(`required_override=.false.`) in `detect_atoms` only
(`single_ui_atom.f90:226-227`); three new keys on `detect_atoms`:
`discover_species` (`yes|no`, default `no`: residual recovery and the species
call, diagnostics only), `nspecies` (integer, default 0 = automatic `K`),
`nn_gate` (`yes|no`, default `yes`); `vol_even` and `vol_odd` (existing
keys) accepted as the optional half maps of 3.2. The thresholds other than
the calibrated ones stay module constants. `exec_detect_atoms`
(`simple_commanders_atoms.f90:405-423`) validates the combination before
`params%new`: `nspecies` and `nn_gate` need `discover_species=yes`; both
half maps or neither. `atoms_stats`, `autorefine3D_nano`,
`conv_atom_denoise` and `analysis2D_nano` keep requiring `element` until
promotion.

Later, with promotion (section 10.5): `element` optional in the other four
programs; `simulate_atoms` sets `occupancy` and `beta` per atom and
`write_centers` writes symbol, chain, occupancy and `B`; data-derived
components for the recovered atoms (section 6); `autorefine3D_nano` passes
the half maps to `detect_atoms` and copies the three species files with the
other per-iteration products; `analysis2D_nano` without an element cannot
build the ideal-lattice start (Q4).

Seen on the way, not changed by this work: `nanoparticle%element` is
`character(len=2)`, so a five-character species string such as `CdSeW` is
cut to `Cd` before every lattice lookup, which then falls back to fcc with
`a` = 3.76 (`simple_nanoparticle.f90:189`, `:262-266`, `:1032-1033`,
`:1470-1471`; `simple_defs_atoms.f90:223-271`). Recorded for the maintainer.

## 10. Phases and the tests that gate them

Sections 10.1 to 10.4 are phases 0 to 3, the run; 10.5 is what comes after
it.

### 10.0 Rules of coexistence

The present detection is working production code with a regression test,
and it stays that way throughout. Concretely:

- No existing routine changes its behaviour for a call as it is made today.
  New behaviour is reached only through `l_species_free` (no `element`) or
  `l_discover_species` (`discover_species=yes`), both false for every
  existing caller and every existing test. Optional arguments added to
  utilities default to the present behaviour.
- `single_atoms_stats` is run on every phase's build and must pass with its
  present floors; its atom count, recall, precision, root-mean-square error
  and correlation are recorded in Progress each time and compared with
  phase 0. Any drift is a stop condition, not something to explain.
- The present products are written by the present code from the level-1
  atoms; the prototype writes only its own three files (section 6).
- No new CMake target, no per-file compile option, no new global flag
  (`doc/policies/compile_time_policy.md`); a new tester goes into the
  existing `unit_single` executable through the suite table; a new CTest
  entry raises `SIMPLE_CTEST_BUDGET` in `production/CMakeLists.txt` by one.
- Release 4 rules: nothing is aliased or kept for compatibility; the only
  on-disk format that must keep reading is the `.simple` project.
- The fast gate (`scripts/run_fast_gate.sh`: the test-registry, description
  and `scripts/check_flex_dag.py` checks, then the `fast` and `provisional`
  CTest labels) passes at the end of every phase, in Debug.

### 10.1 Phase 0: baseline

Build the copy as it was handed over (Debug with `BUILD_TESTS=ON`, the
machine's compiler and environment recorded), run the fast gate and
`single_atoms_stats`, and record its numbers. Then run the same simulation
and detection by hand (`simulate_nanoparticle element=Pt moldiam=20
box=160 smpd=0.358`, `detect_atoms element=Pt`) and keep the input volume,
`simatms.pdb` and `_ATMS.pdb` under the run's `scratch/keep/`: they are the
reference of phase 2's equivalence test. Record the atom count of
`simatms.pdb` and of `_ATMS.pdb`, the positions' root-mean-square
difference and the `_SIM.mrc` correlation. Exit: fast gate passes;
`single_atoms_stats` passes; the reference is kept and its numbers are in
Progress.

### 10.2 Phase 1: foundations that change nothing

The pseudo-atom kernel in `convolve` with its symbols and the PDB round
trip through `occupancy` and `beta`, so that `simulate_nanoparticle
pdbfile=` renders pseudo-atoms; and `simple_nano_species` with its tester.
Exit: unit tests on one and two pseudo-atoms against the closed forms (peak
`q (4 pi / B)^(3/2)`, voxel sum `q / smpd^3` within 1% at a cutoff of
`4 sigma`, two atoms adding linearly), and a round trip PDB to `atoms` to PDB
preserving symbol, occupancy and `B`; `discover_species` returns the right
`K`, labels and BIC ordering on fixed samples (one class; two classes at
`D` of 2.5 and 4, the first inadmissible; a class below the size floor;
three classes; `fit_species_mixture`) and the calibration returns 5.0
within 0.1 from the region counts of section 8.8 (white noise) scaled to
its geometry; the `species`
sub-suite is registered and `scripts/check_test_registry.py` is clean; fast
gate; `single_atoms_stats` unchanged against phase 0.

### 10.3 Phase 2: level 1 without an element

`d_NN`, the derived length scales, the Gaussian template, `element`
optional in `detect_atoms`. Exit: `detect_atoms` on the phase 0 input
volume without `element` gives the same atom count as the reference within
1%, positions matched within 0.5 A with a root-mean-square difference of
at most 0.1 A, and the same `_SIM.mrc` correlation within 0.01; the same
command with `element=Pt` gives a `_ATMS.pdb` identical to the reference
(same count, positions within 1e-3 A: the element path did not move); the
measured `d_NN` is within 2% of 2.766 A (Pt, `a` = 3.912); fast gate;
`single_atoms_stats` unchanged. If the species-free count differs by more
than 1% the phase stops and reports which step (template, threshold
ladder, splitting, pruning) lost or gained the atoms, with the counts after
each.

### 10.4 Phase 3: residual recovery and species call, diagnostics only

Noise region and half-map reference, calibrated thresholds, the two stages
with the gate switch, pruning of the recovered atoms on the atom table,
both fit stages, aperture intensities, species assignment, radial profiles,
the three files of section 6. The fixtures are rendered by
`simulate_nanoparticle pdbfile=` from pseudo-atom PDBs written by the test
itself (the generating models of cases 2, 9, R3 and R6 of section 8 on the
Pt lattice of phase 0's particle, so that level 1 is known to work on them),
with Gaussian noise added in process at a fixed seed and at the section 8
noise levels relative to the peak of a core atom of the strongest class.
Exit, as a new high-level test `single_species_discovery` (registered like
`single_atoms_stats`, budget raised by one): for the two-species fixtures,
recall of the weak class at least 0.90 where its predicted detection
signal-to-noise (section 7) is 6.5 or more, no more than one false atom per
particle, `K` = 2, every label of a found atom right, the intensity ratio
within 10%, the per-shell width of the strongest class within 3% and of the
light class within 10%; for the single-species fixtures, `K` = 1 and no
recovered atom; for every fixture, `_ATMS.pdb`, `_BIN.mrc`, `_CC.mrc`,
`_MSK.mrc` and `_SIM.mrc` identical to the run of the same command without
`discover_species` (the products did not move); fast gate;
`single_atoms_stats` unchanged. The floors are those of the emulation and
are set before the first run; a floor the code does not reach stops the
phase with the measured value, it is not lowered.

### 10.5 After the prototype (not part of the run)

Validation on an independent forward model: particles rendered with the
existing five-Gaussian element kernels of `convolve` (Pt and Ni, not the
Gaussian the detector assumes), passed through the reconstruction's
low-pass filtering and soft mask, with background and shell-dependent noise;
cross-half reproducibility of the labels on real half maps; an A/B of
`detect_atoms` with and without discovery on a real Pt/Ni map at matched
false-positive rate, judged by the maintainer. Then promotion: recovered
atoms into `_ATMS.pdb` with data-derived components in `_CC.mrc`, class
intensities and per-atom widths into `_SIM.mrc`, `element` optional in the
other four programs, `autorefine3D_nano` end to end on a simulated
two-species trajectory, `atoms_stats` without the lattice table.

## 11. Relation to `single_pcg_integration.md`

With pseudo-atoms the columns of `G` are the Gaussians of section 1, `Z`
drops out, and the least-squares amplitudes `c = (G^T G)^-1 G^T x` are the
per-atom amplitudes with neighbour overlap removed exactly. The fit of 3.4 is
a Gauss-Seidel approximation of that solve. Two points in that note need
changing for Z contrast. Its section 5, step 5 makes atoms with
`c_i < 0.2 median(c)` removal candidates, which would flag a whole light
species; the rule should be referenced to the noise. And its statement that
trailing reconstruction blends the simulated map away does not match the
accumulator-domain contract cited in section 0.

## 12. Deliberately absent

- **Element names from the data.** With two classes the intensity ratio
  constrains `(Z_2 / Z_1)^n` for an exponent between about 1.5 and 2, and the
  fitted lattice constant with the class fractions constrains the pair
  through Vegard's law. Both depend on calibrations (the exponent, the pixel
  size) that the map does not carry. Names come from the user (section 5) or
  not at all.
- **Spatial priors on the labels** (3.6).
- **A lattice model.** The gate and the length scales use `d_NN` only. Testing
  the empty sites of a fitted lattice at a lower threshold would reach weaker
  atoms, at the price of assuming a single crystal.
- **Anisotropic widths** in the fit. `calc_anisotropic_disp` stays in
  `atoms_stats`.
- **Per-atom amplitudes in `_SIM.mrc`** (section 4).
- **A bank of template widths or of element templates** in detection
  (sections 5 and 7).
- **Model-generated connected components** (section 6, withdrawn).
- **Fixed detection thresholds** (3.2): they are calibrated from the map.
- **Any change to the present products or the present pruning** before the
  validation of 10.5.

## 13. Open questions

- **Q1.** Is 3-6:1 the atom ratio or the intensity ratio? And what is the
  signal-to-noise of the weakest atoms you can see in a current map
  (section 7)? Together they say whether the light species is in the regime of
  cases 1-4 or of cases 6 and R2.
- **Q2.** Class symbols in the PDB element column (`X1`, `X2`, ...) with the
  class also in the chain: acceptable, or one symbol with the class in the
  chain only?
- **Q3.** Radius from the centre as the radial coordinate, or depth below the
  surface?
- **Q4.** Starting volume without an element: `solve3D_nano`, or the
  ideal-lattice start with a user-supplied lattice constant?
- **Q5.** Is the close-packed neighbour cutoff an acceptable default until the
  pair-distance form is in?
- **Q6.** Stage 2 assumes equal integrals within a class. Is that acceptable
  for surface atoms, or should the stage-1 width be the reported one where
  the two disagree?
- **Q7.** The target of 0.2 expected false atoms per particle behind the
  calibrated thresholds, and the gate on by default: both are prototype
  defaults to be judged on the Pt/Ni A/B of 10.5.

## 14. File table

The files the prototype may change, and in which phase. A phase edits only
its files; a file the work turns out to need is added here, with its phase
and a one-line reason, before it is edited. Testers of a listed Fortran file
(`*_tester.f90`) are covered by its row. Phase 0 edits nothing but this
note.

| File | What changes | Phases |
| --- | --- | --- |
| `src/main/nano/simple_atoms.f90` | pseudo-atom branch in `convolve`; the symbols accepted by `element_exists`, `set_element`, `guess_an_element`, `Z_and_radius_from_name`; closed-form tests in its tester | 1 |
| `src/defs/simple_defs_atoms.f90` | pseudo-atom symbols and their sentinel range in `get_element_Z_and_radius` | 1 |
| `src/main/nano/simple_nano_species.f90` | new: mixture fit, BIC, admissibility, `discover_species`, enclosed fraction, threshold calibration; its tester | 1, 3 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | `species` sub-suite in `suites_single` | 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | `unit_single` suite list | 1 |
| `src/main/nano/simple_nanoparticle_utils.f90` | optional kernel width in `phasecorr_one_atom`; `est_nn_dist`; optional `d_NN` in `find_rMax` and `calc_contact_scores` | 2 |
| `src/main/nano/simple_nanoparticle.f90` | `l_species_free` and the `d_NN` length scales (phase 2); the residual recovery, fits, species call and the three files behind `l_discover_species` (phase 3) | 2, 3 |
| `src/main/ui/single/single_ui_atom.f90` | `element` optional in `detect_atoms` (phase 2); `discover_species`, `nspecies`, `nn_gate`, `vol_even`, `vol_odd` on `detect_atoms` (phase 3) | 2, 3 |
| `src/main/ui/simple_ui_params_common.f90` | only if a new key has to be shared rather than declared on `detect_atoms` | 3 |
| `src/main/params/simple_parameters.f90` | `discover_species`, `nspecies`, `nn_gate`, `l_discover_species`, `l_nn_gate` | 3 |
| `src/main/params/simple_parameters_parse.f90` | registration of the new keys | 3 |
| `src/main/params/simple_parameters_phases.f90` | phase 2 only if the element validation at lines 937-941 needs it for an absent key; phase 3 the derived logicals and the validation of the new keys | 2, 3 |
| `src/main/commanders/simple/simple_commanders_atoms.f90` | `exec_detect_atoms`: command-line validation, half maps | 2, 3 |
| `src/main/commanders/test/simple_commanders_test_single.f90` | the high-level test `single_species_discovery` | 3 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | its program entry | 3 |
| `src/main/exec/simple_test_exec_single.f90` | its router case | 3 |
| `production/CMakeLists.txt` | its CTest entry; `SIMPLE_CTEST_BUDGET` raised by one | 3 |
| `doc/implementation_notes/planned/single_species_discovery.md` | this note: Progress, file-table rows, rulings | each |

## 15. Review of 2026-10-07 and what it changed

The first version of this note was reviewed on 2026-10-07 by a second
reviewer and by the maintainer, who declined to approve it as written and
approved a small, opt-in prototype instead. The findings, what was checked,
and what the note now says:

The pipeline order broke the pruning contract. Fitting and classifying
before `discard_atoms` would have lost every per-atom quantity, because
pruning rebuilds the atom table (section 0, last constraint but one), and
components generated from the model would have made the size test of the
pruning circular. Confirmed from the source. Pruning now comes first and
the recovered atoms are pruned on the atom table (3.9); model-generated
components are withdrawn in favour of data-derived ones at promotion
(section 6).

The aperture intensity lacked a factor `smpd^3` against the continuous
integral it is compared with. Confirmed; added (3.5).

A fixed threshold on local maxima of a z map is not statistically
controlled, the voxel standard deviation does not calibrate the maxima, and
undetected weak atoms can contaminate the noise region. Checked by
simulation (8.8): the voxel tail indeed overpredicts the number of noise
maxima six times, but the region's own maxima calibrate it to within 0.8 to
1.0, and the calibrated threshold equals the fixed one for the emulated
geometry (5.0) and differs for strongly correlated noise (4.7). The
threshold is now calibrated per map (3.2). Contamination is real for a
weak-only domain thicker than 1.5 `d_NN` inside a tight mask (36 to 40% on
the standard deviation) and absent for a surface shell; the half-map
difference is now the noise reference when the halves are given, and the
robust-to-plain spread ratio is reported otherwise. The existing
false-discovery-rate routine of `simple_segmentation.f90` was benchmarked on
the same maps: it admits three to four times the gate's false atoms at
lower recall, so it is not adopted.

The three-neighbour gate is a lattice prior that will disadvantage
surfaces, defects, interfaces and weak-only domains. Checked by simulation
(8.8): at the same threshold the gate costs about six points of recall at
low coordination (six to nine points); at the same false-positive rate it gains 7 to 15 points
everywhere, including the low-coordination bins and a weak-only domain. Its
assumptions and its blind spot are now stated (3.3), and it is a switch with
its use recorded per atom.

The evidence is model-matched: the generator produces the Gaussian atoms the
detector assumes, several production steps were not emulated, and some
threshold numbers came from an earlier variant of the algorithm. Agreed. The
numbers from the earlier variant are replaced (3.3), the independent check
is described for what it is (the preamble, 8.8), and validation on an
independent forward model with a real Pt/Ni comparison gates promotion
(10.5).

What the review asked for and the note now does: the existing path is
unchanged and remains the production path with its regression test; the
recovery runs behind an experimental key and produces diagnostics without
feeding `_SIM.mrc`; components stay data-derived; pruning precedes fitting
and classification; the forward-model validation and the Pt/Ni comparison
at matched false-positive rate come before promotion.

## 16. Progress

Each phase records here: the date, the files changed, the evidence with log
paths under the run directory, and each exit item of section 10 with how it
was met. The `single_atoms_stats` numbers are recorded for every phase next
to phase 0's.

(nothing yet)

## Appendix. Emulation scripts

`emul_detect.py`: `atoms%convolve`, `phase_corr`, `otsu`, `sortmeans`,
26-connected components, raw-weighted centroids and the bisection of
`binarize_and_find_centers`, written from the source. `emul_species.py`: the
generating model, sections 3.1 to 3.8 and 4 of this note, the three statistics
of 8.4, and the case list (`python3 emul_species.py`, or with case labels as
arguments; about 25 s per case). Seeds are fixed. The template-width bank of
section 7 and the "13 noise peaks returned as a species" observation of 3.3
come from an earlier variant of the residual levels (connected components,
exclusion at `0.5 d_NN`); the gate numbers of 3.3 are from the current
variant (8.8).

`codex_checks.py` (numpy only, kept with the other two): the independent
generator and detector of 8.8, with `python3 codex_checks.py null | contam |
gate | all` (about two minutes for `all`). It ports `detect_peak_thres_fdr`
and implements the threshold calibration of 3.2 as `calibrate_k`.
