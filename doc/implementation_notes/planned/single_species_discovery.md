# SINGLE without a species input: atom detection, Gaussian density model, species discovery and map simulation

Implementation note, 2026-10-07. Design for review; nothing here is
implemented. Written against master `108c8c9d`. Companion page:
`single_pcg_integration.md` (section 11 says how the two fit together).

Nothing in SIMPLE was compiled or run. The numbers in section 8 come from a
numpy emulation of the proposed pipeline on synthetic particles
(`emul_species.py`, which needs `emul_detect.py`; both are kept out of the
tree, like the other reference scripts). They show how the algorithm behaves
on idealised input and are not a measurement on data.

## 0. What is asked, and what constrains the design

The request: no `element` on the command line. The program finds the atoms,
finds how many kinds of atom there are and which atom is which, and simulates
a map from what it found. The first data are HAADF-STEM images of Pt/Ni
particles.

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
  system) becomes the median contact score of the inner 30% of the atoms.
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

Nine steps inside `detect_atoms`. Step 1 and step 9 are the present flow with
a different kernel. The rest is new.

### 3.1 Level 1: the present flow with a pseudo-atom

`phasecorr_one_atom`, the Otsu threshold ladder, the bisection,
`discard_small_ccs` and `split_atoms` run as now, with one change: the
template and the equal-atom simulation of `t2c` use the Gaussian of section 1
at a reference width `B_ref`, not an element. At 0.358 A the Pt template
correlates at 0.997 with a single Gaussian of sigma 0.42 A (`B` = 13.9 A^2),
which is the starting `B_ref`. After the first fit (3.4) `B_ref` becomes the
median `B` of the strong atoms, and if it moved by more than 10% the
detection is run once more. `detect_atoms` stays stateless.

Level 1 finds the strongest class and nothing else. That is the behaviour to
preserve for single-species data, and the reason more levels are needed: the
Otsu gate sits at about a quarter of the strong peak, and in every
multi-species case of section 8 level 1 found the strongest class (all of it
in all cases but one, 394 of 396 there) and none of the weakest.

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

Two stages. Stage A uses `k` = 5 and accepts any candidate farther than
`0.7 d_NN` from every accepted atom; it repeats until it adds nothing. Stage B
then uses `k` = 4 and accepts a candidate only if it passes a local lattice
gate: no accepted atom closer than `0.85 d_NN`, and at least three within
`1.15 d_NN`. It also repeats until it adds nothing. At most eight levels.

Three choices here were forced by the emulation.

- Local maxima, not connected components. At a fixed low threshold the
  components of neighbouring atoms merge when the signal is strong, and the
  merged centroid lands between two atoms.
- No Otsu on the residual. It has no noise reference and can fall to the
  noise floor.
- The gate. Lowering `k` to 4 without it admitted 6 to 13 noise peaks per
  particle, and in one single-species run 13 of them were returned as a
  second species holding 2.4% of the atoms. That is how a species gets
  invented. With the gate, `k` = 4 gave no false atom in the same four cases
  and raised the recall of weak, broad atoms (93 to 115 of 132, and 15 to 63).
  The gate uses only `d_NN`, so it does not assume one lattice across the
  particle.

The gate cannot find a weak atom with fewer than three accepted neighbours.
A region made only of the weak species is reached from its boundary inwards,
one shell per level.

This is the one step that matches against a shape, and section 7 treats what
that means at the surface.

### 3.4 Density fit, stage 1: free amplitude

Three sweeps over all atoms. For atom `i` the data are the map minus the
background minus the fitted density of every other atom, inside a sphere of
radius `d_NN / 2`. The amplitude is linear; `B_i` is found by a
one-dimensional search that minimises

    SSE(B) / sigma_n^2 + ((ln B - m_i) / tau_i)^2

`m_i` and `tau_i` are the median and the robust spread (floor 0.1) of `ln B`
over the strong atoms (amplitude above `10 s_A`) of the same radial shell,
within `d_NN` in radius; the global values are used when fewer than five
qualify. After each sweep the background is the residual smoothed with a
Gaussian of standard deviation `d_NN`.

A strong atom keeps its own width, since its data term dominates. A weak atom
cannot: with a free amplitude its width is not determined by its own density
(in the earlier estimator comparison the fitted sigma of an atom at 0.175 of
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

    I_i = sum over |r - r_i| < d_NN/2 of [ map - b - sum_{j /= i} rho_j ] / F(d_NN / 2 sigma_i)

`F(t) = erf(t / sqrt 2) - sqrt(2 / pi) t exp(-t^2 / 2)` is the fraction of a
Gaussian inside `t` standard deviations, and `sigma_i^2 = B_i / 8 pi^2`. The
correction is 0.2% for the narrowest atoms of the emulation and 6% at twice
their variance.

This is the classification statistic. It is a plain sum, so it involves no
matching; it measures the quantity that carries the species; and it depends
on `B_i` only through the small correction. Its cost is noise: about 2.3
times that of the matched amplitude for white noise, and more when the noise
is concentrated at low frequency (case 8b of section 8.1).

### 3.6 Species discovery

Input: the `I_i` and `s_I`. A one-dimensional Gaussian mixture is fitted by
EM for `K` = 1, 2, 3, with every component variance floored at `s_I^2`, from
two deterministic starts (sorted quantiles, as in `sortmeans`, and the `K-1`
largest gaps of the sorted sample). Each fit gets
`BIC = -2 ln L + (3K - 1) ln N`.

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

`discard_atoms` and `validate_atoms` run once, on the merged atom set, with
the length scales of section 2. `valid_corr` is measured against the
simulation of section 4, in which every atom has its class intensity and its
own width. Two penalties of the present equal-atom, equal-width simulation
disappear: the one on a weak atom next to strong neighbours, and the one on a
broad atom at the surface.

## 4. Simulation

`atoms%convolve` gets a pseudo-atom branch: for a sentinel `Z` the kernel is
the Gaussian of section 1 with `B` from `beta` and amplitude
`q (4 pi / B)^(3/2)`, `q` from `occupancy`. An atom with `q` = 1 has unit
integrated intensity whatever its width.

`_SIM.mrc` is rendered with the class intensity for `q` and the fitted `B_i`
of each atom. Against the noise-free generating density of the multi-species
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
`K` is the length of the list and the classes take its names. Today a
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

- `_ATMS.pdb`. Element column: the class symbol (`X1`, `X2`, `X3`, registered
  in `get_element_Z_and_radius` with sentinel `Z`, precedent `CDSE`, 999; or
  the element names when a species list is given). Chain: `A`, `B`, `C` by
  class, so a class is one selection in a viewer. Occupancy: `I_k / I_1`.
  B column: `B_i`. Today the B column holds `valid_corr`
  (`simple_nanoparticle.f90:1821`); within SIMPLE only the unit tester reads
  it back, and `valid_corr_in_bfac_field.pdb` keeps it for `atoms_stats`.
- `_species.csv`, one row per atom at full precision: position, radius,
  coordination number, level, detection z, amplitude, `B` of stage 1 and of
  stage 2, aperture intensity and its signal-to-noise, class, class
  probabilities.
- `_species_radial.csv`: the profiles of 3.8.
- `_species.txt`: `K`; class intensities, ratios and fractions; `D`; the BIC
  table with admissibility; the expected misclassification rate; `d_NN`,
  `B_ref`, the noise numbers; atoms added per level; the diagnostics of
  section 7.
- `_SIM.mrc` as in section 4.
- `_CC.mrc`, `_BIN.mrc`, `_MSK.mrc`. Atoms found at residual levels have no
  connected component of the level-1 threshold. In species-free mode the
  component of atom `i` is the set of voxels nearer to it than to any other
  atom and within `min(d_NN / 2, 2 sigma_i)`. `atoms_stats` reads `_CC.mrc`
  with `_ATMS.pdb`, so this keeps it working. Not emulated.

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
  z and level.
- **Half-map agreement.** With `vol_even` and `vol_odd` given, the aperture
  intensities are measured at the final positions in each half map, and the
  label agreement and the intensity correlation are reported. The halves
  share their alignment, so this measures the effect of noise on the labels
  and not model bias. Not emulated.
- **Model bias.** Labels are found from scratch at every call; nothing is
  carried between iterations. Since `_SIM.mrc` only aligns (section 0), a
  wrong label can reach the next map only through the orientations. A direct
  test is cheap and worth having during development: render a random 10% of
  the atoms at the pooled intensity and check that they separate as well as
  the rest in the next iteration. Not emulated.

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
statistic: it reads a broad heavy atom as a lighter species. The review of
2026-10-07 recommended it on the evidence of particles with uniform widths;
this table withdraws that recommendation.

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

`split_atoms`, `discard_atoms`, `validate_atoms`, the soft mask and the ICM
filter were not emulated. Atoms are exact Gaussians, the background is zero,
positions are not refined, the particle is a sphere, and the noise is
stationary. The noise levels are guesses; what transfers is the
signal-to-noise scale of section 7. The thresholds `k` = 5 and 4, `D` = 3,
the class-size floor and the 5% floor on `tau_k` are starting values.

## 9. Code map

New: `src/main/nano/simple_nano_species.f90`, array numerics with no image
dependence: the mixture fit, BIC, admissibility and the choice of `K`
(`discover_species(x, s_meas, nspecies_in, K, labels, post, mu, var)`), and the
enclosed-fraction function. Tester `simple_nano_species_tester.f90` in
`unit_single`.

`simple_atoms.f90`: pseudo-atom branch in `convolve`; the per-element table
and `epot` stay as they are.

`simple_defs_atoms.f90`: the pseudo-atom symbols.

`simple_nanoparticle_utils.f90`: `phasecorr_one_atom` takes a kernel width
when no element is given; `est_nn_dist(centers)`; `find_rMax` and its three
callers accept `d_NN`.

`simple_nanoparticle.f90`:

- `new`: no `Z` lookup without an element; keeps the volume as read next to
  the masked one.
- `atom_stats`: `amp`, `bfac`, `bfac_stage1`, `aper_int`, `aper_snr`, `det_z`,
  `det_level`, `species`, `species_post`. `cendist`, `u_iso` and `isocorr`
  exist.
- new procedures: `est_noise_region`, `detect_residual_levels`,
  `fit_atoms_joint` (both stages), `calc_aperture_int`, `assign_species`,
  `write_species_report`, `write_radial_profiles`, `make_cc_from_model`.
- `identify_atomic_pos`: sections 3.1 to 3.9 in order when species-free; the
  present sequence otherwise.
- `simulate_atoms`: sets `occupancy` and `beta` per atom; `write_centers`
  writes symbol, chain, occupancy and `B`.

`calc_isotropic_disp` stays for `atoms_stats`. Its log-linear estimator over
positive voxels is biased in noise (at noise 0.03 it returned amplitude 0.91
and sigma 0.371 A for true 1 and 0.350); `atoms_stats` can take `u_iso` from
the stage-2 fit once that exists.

Parameters and UI: `element` becomes optional (`required_override=.false.`) in
`detect_atoms`, `autorefine3D_nano`, `atoms_stats`, `conv_atom_denoise` and
`analysis2D_nano`; its absence selects discovery, a species list fixes `K` and
the names (section 5). One new key, `nspecies` (0 = automatic). The
thresholds stay module constants until data say otherwise.

`single_commanders_nano3D.f90`: `autorefine3D_nano` passes the half maps to
`detect_atoms` and copies the three species files with the other
per-iteration products.

`single_commanders_nano2D.f90`: without an element `analysis2D_nano` cannot
build the ideal-lattice start (Q4).

## 10. Order of work and the tests that gate each step

1. **Pseudo-atom kernel** in `convolve`, symbols, PDB round trip.
   `simulate_nanoparticle pdbfile=` renders pseudo-atoms. Gate: unit tests on
   one and two atoms against closed forms (peak `q (4 pi / B)^(3/2)`, sum over
   voxels `q / smpd^3`).
2. **`simple_nano_species`** with its tester. Gate: `K`, labels and BIC on
   fixed samples (one class; two classes at `D` of 2.5 and 4; a class below the
   size floor; three classes), expected values pinned from `emul_species.py`.
3. **Level 1 without an element**: `d_NN`, the derived length scales, the
   Gaussian template. Gate: the Pt case of `lib_single`, sub-suite
   `nanoparticle atoms`, gives the same atom count with and without `element`,
   positions within a stated tolerance.
4. **Noise region, residual levels, both fit stages, aperture intensity,
   species assignment, radial profiles, products.** Gate: a new sub-suite of
   `lib_single` on fixtures from step 1 with noise at a fixed seed: the
   generating models of cases 2, 9, R3 and R6, asserting recall and precision
   per class, `K`, label accuracy, the intensity ratio and the per-shell
   widths within stated tolerances, and for the single-species cases that no
   atom is added and `K` = 1.
5. **`autorefine3D_nano`** end to end on a simulated two-species trajectory
   (`single_workflow` with a species-free suite): final labels against the
   generating model, half-map agreement reported.
6. **`atoms_stats`** without the lattice table: neighbour cutoff, coordination
   and strain reference from the fitted lattice; `u_iso` from stage 2.

Steps 1 and 2 change nothing for existing runs. Step 3 changes nothing unless
`element` is left out.

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

## Appendix. Emulation scripts

`emul_detect.py`: `atoms%convolve`, `phase_corr`, `otsu`, `sortmeans`,
26-connected components, raw-weighted centroids and the bisection of
`binarize_and_find_centers`, written from the source. `emul_species.py`: the
generating model, sections 3.1 to 3.8 and 4 of this note, the three statistics
of 8.4, and the case list (`python3 emul_species.py`, or with case labels as
arguments; about 25 s per case). Seeds are fixed. The threshold scans behind
the gate numbers of 3.3 and the template-width bank of section 7 were run
with an earlier variant of the residual levels (connected components,
exclusion at `0.5 d_NN`).
