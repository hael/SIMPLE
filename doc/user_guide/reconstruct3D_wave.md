# reconstruct3D_wave — Reconstruct Under a Wave Model of the Image

`simple_exec prg=reconstruct3D_wave` reconstructs half maps from oriented particles with known CTF, like
`reconstruct3D`, but predicts each particle image with wave optics instead of SIMPLE's weak-phase model (image =
CTF times the central section of the map's transform). The wave models account for two effects the weak-phase
model leaves out:

- **Ewald-sphere curvature**: a thick specimen's image samples the map's transform on a curved surface (the
  Ewald paraboloid), not on the central plane. It limits large particles, low voltages and high resolution: beyond
  roughly `sqrt(2 lambda D)` (wavelength `lambda`, particle diameter `D`; 4.7 A for a 300 A particle at 100 kV) a
  flat reconstruction loses the particle's periphery.
- **Strong phase**: when the specimen shifts the electron wave's phase by much more than a few tenths of a radian
  (metal nanoparticles, very large complexes at low voltage), the image is no longer linear in the potential.

The models follow those of Ghostbuster, the reconstruction inside the specter package, against which they were
checked.

Poses, shifts, states, half sets and CTF come from the project and are not refined. The program runs in shared
memory on one machine.

## Models

`wave_model` selects the image formation; `ewald` the curvature sign of the paraboloid.

- `rytov` (default): the exponential of the wave's phase and amplitude change, on the paraboloid. Accounts for
  both effects; the cheapest model with curvature. Use it for cryo-EM data.
- `first_born`: the linear model on the paraboloid. Curvature only.
- `projection`: the exponential on the central section (`ewald=no` required). The strong phase only.
- `multislice`: slice-by-slice transmission and propagation through the rotated potential, multiple scattering
  included (`ewald=neg` or `pos` required: the sign is the beam direction). Some sixty times more expensive than
  Rytov; meant for nanoparticle boxes and small sets.

`ewald=neg` (default) is the sign of specter's default and of data simulated with it: the beam enters through the
map's first `z` section, and a feature at `+z` of the box centre images as if at larger underfocus. If the sign of your data is unknown,
reconstruct with both and compare their half-map FSC (Fourier shell correlation) beyond the curvature limit; the wrong sign is clearly worse.

## Running

```sh
simple_exec prg=reconstruct3D_wave projfile=my_proj.simple \
  wave_model=rytov ewald=neg total_dose=40 objfun=cc \
  mskdiam=320 pgrp=c1 nthr=24
```

Required: `projfile` (particles with poses, states, CTF parameters `dfx`, `dfy`, `angast`, `phshift`, `kv`, `cs`,
`fraca`), `total_dose`, `pgrp`, `mskdiam`.

- `total_dose`: the exposure in electrons per square angstrom. Its square root times the pixel size is the gain
  that links SIMPLE's noise-normalised particles to the relative intensity the wave models predict. The nonlinear
  models need the potential in absolute units (volts), so this is not optional.
- `objfun=euclid` (SIMPLE's default) weights the misfit by the project's noise spectra, which an earlier
  refinement or reconstruction registered; without them the program stops. `objfun=cc` weights every frequency
  alike and needs none, as in the example.
- `maxits` (default 4): nonlinear conjugate-gradient iterations per half, one pass over its particles each.
- `pcg_mskfile`: a real-space mask volume whose nonzero voxels are the solve's support; without it the
  `mskdiam` sphere is the support. `automsk` is refused.
- `box_crop`: reconstruct at a Fourier-cropped box, as in `reconstruct3D`.
- `maxits_pcg`, `maxits_ml`: the flat start solve, as for `reconstruct3D rec_backend=pcg`. `tau` (default 1),
  the strength of the PCG backend's prior, is accepted on the command line though the program's interface does
  not list it.
- `winsz` (default 1.5): the Kaiser-Bessel window of the operator's lattice; 3 gives seven taps and ten times
  smaller interpolation errors at about 1.4 times the run time. The default is far below a particle's noise.

`nparts` is refused: distributed execution comes when the wave models join `reconstruct3D` and `refine3D`.

## What it does

For every state:

1. The particles are read once into memory when they fit (the log reports `particle buffer: filled`); otherwise
   every pass reads them from the stacks. The even and odd halves are the project's.
2. Each half is reconstructed with SIMPLE's flat PCG (preconditioned conjugate gradients) backend; the FSC of
   the two maps sets the per-shell maximum-likelihood prior; each half's map is regularised
   with it. These flat maps are the starting point and are written as `recvol_stateNN_even_flat.mrc` and
   `_odd_flat.mrc`: the flat model under the same prior and support, the fair comparison for the wave maps.
3. Each half is solved under the wave model from its flat start, by nonlinear conjugate gradients on the
   potential in volts.
4. The maps are written at the scale of `reconstruct3D` maps, so the two compare directly.

## Outputs

In the new numbered execution directory (`mkdir=yes` by default):

| File | Contents |
| --- | --- |
| `recvol_stateNN.mrc` | the mean of the regularised halves |
| `recvol_stateNN_even.mrc`, `_odd.mrc` | the regularised half maps |
| `recvol_stateNN_even_unfil.mrc`, `_odd_unfil.mrc` | the half maps with the regularisation's Wiener filter undone, for FSC and postprocessing |
| `recvol_stateNN_even_flat.mrc`, `_odd_flat.mrc` | the regularised flat starts |
| `fsc_stateNN.bin`, `RESOLUTION_STATENN` | the FSC of the unfiltered pair |
| `reconstruct3D_wave_stateNN_<half>_diagnostics.txt` | model, the gain from the dose and the fitted gain, the solver's loss per iteration, particle passes |

The merged map and the FSC are registered in the project's output segment. The particles' records are not changed.

The log prints `GAIN FROM DOSE` and `GAIN FITTED` per half. The fitted gain is the scale the data prefer for the
nonlinear part of the image; it is reported, never used. Within a few percent of the dose's gain the scale link is
sound. Well below it, the particles' noise normalisation saw signal in the background (a box with little margin
at high defocus) or the dose is not the dose of the images (dose-weighted sums); the nonlinear terms are then
scaled wrongly by that factor. For first Born the fitted value is not a gain and can be ignored.

## What to expect

On simulated data with a known answer (the run that built the program, recorded in
`doc/implementation_notes/planned/ghostbuster_operators.md`):

- An ordinary protein at 300 kV, 1 A pixels: every model agrees with the flat map; Rytov is slightly ahead
  (truth FSC 0.762 against 0.733 between 3 and 2 A).
- A 300 A ribosome at 100 kV, where the curvature limits the map beyond 4.7 A: Rytov with the data's sign raises
  the truth FSC between 3 and 2.4 A from 0.66 (flat) to 0.81; the wrong sign drops it to 0.21.
- A Pt nanoparticle at 300 kV: Rytov and multislice place the atoms more accurately than the flat map (0.029 against
  0.035 A error) and partly repair the band near the particle's own curvature limit (about 1 A). The default prior
  is too strong in that band; `tau=10` (a ten times weaker prior) recovers most of it.

The wave solve costs a few times a PCG reconstruction of the same particles for Rytov and first Born (the ribosome:
about 9 min against 1.5 min at 24 threads), and far more for multislice (90 min for 2 000 particles at box 160 on
12 threads).

The conventions, the scale link and the tests are in `doc/policies/3D/wave_model_policy.md`.
