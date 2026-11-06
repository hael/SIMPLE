---
name: simple-main-ctf
description: Use when working in SIMPLE's src/main/ctf subsystem, including the microscope contrast transfer function model (the ctf type, defocus, astigmatism, amplitude contrast and phase-shift conventions and units, ctfparams versus ctfvars), micrograph CTF estimation (tiled power spectra, grid search, differential-evolution and continuous refinement, phase-shift fitting, patch defocus polynomials and the ctfdoc), its cost functions, and how fitted CTF parameters reach particle records.
---

# `src/main/ctf`: the CTF model and micrograph CTF estimation

Owns the `ctf` type (CTFFIND4's transfer function) and the fit of one micrograph:
`ctf_estimate_iter` (fit, diagnostics, the micrograph's row), `ctf_estimate_fit` (tiled spectra,
search, refinement, patch polynomials, the `ctfdoc`) and its costs. Applying a CTF to images is
`image/` (`simple_image_ctf`, `gen_fplane4rec`) and `utils/math/simple_math_ctf`; polar CTFs are
`pftc/`; the parameters live in `project/`; choosing micrographs and writing the project is the
ctf_estimate, preprocess and extract strategies in `strategies/parallelization` (the stream runs
`preprocess` jobs). `canonical_phshift` is in `src/defs/simple_type_defs.f90`.

## Rules

Silent ones first: nothing reports these.

- Units. Project rows and `ctfparams` (`ori%get_ctfvars`, `sp_project%get_ctfparams`): `dfx`, `dfy`
  in microns, underfocus positive; `angast` in degrees, the `dfx` axis from +x; `phshift` in radians,
  canonical in [0, 2pi); `kv` in kV, `cs` in mm, `fraca` a fraction, `smpd` in Angstrom. The `ctf`
  object holds defocus, Cs and wavelength in pixels and `angast` in radians, and `ctf%get_ctfvars`
  returns a `ctfvars` in those units: never mix it with a `ctfparams` (`ori%get_ctfvars` returns a
  `ctfparams` despite its name). Fit keys: `dfmin`, `dfmax`, `astigtol` in microns, `hp`, `lp` in
  Angstrom, `phshift_min/max/step` in degrees.
- Build the `ctf` at the sampling of the image it is applied to (a cropped particle needs
  `ctf(smpd_crop, ...)`), then `init(dfx, dfy, angast)` per record before evaluating: the
  constructor leaves defocus at zero. One object holds one record: one per thread in a parallel loop.
- The phase is canonical modulo 2pi, never pi (a pi fold negates the CTF). A spectrum fit fixes it
  only modulo pi; the search window picks the branch, and a window of 180 degrees or more only warns.
  Optimised routines take the phase as a required argument; zero is a real value.
- A fit leaves `dfx >= dfy` and `angast` in [0, 180] (`apply_convention`): compare fits up to the
  swap of `dfx` and `dfy` with `angast` + 90.
- The CTF formula has copies: `simple_ctf` (`eval_phase`, `eval_canonical`), `ft_map_ctf_kernel`
  (`utils/math/simple_math_ctf`), inline in `volume/simple_reconstructor_pcg` and
  `strategies/parallelization/simple_cls_expansion_strategy`, and the analytic-gradient cost of
  `ctf_estimate_cost4Dcont`. A change to the model changes all of them.
- `ctf_estimate_fit%new` memoises the process-wide tables of `image/simple_memoize_ft_maps` at its
  spectrum box (`THROW_HARD` inside an OpenMP region) and leaves them memoised: `kill` releases only
  the fit's own state, so another user of the maps (image CTF kernel, class averager) must memoise its
  own box before use, not rely on the fit's. `new` also reseeds the RNG (`seed_rnd`) per micrograph:
  intended.
- `ctf_estimate_fit` points at the caller's micrograph (`new` back-transforms it): keep it unchanged
  until `kill`. Order: `new` -> `fit` -> `fit_patches` (starts from the global fit; does nothing
  unless there are more than two tiles per axis) -> `write_doc`, `write_diagnostic` -> `kill`.
  `read_doc` + `pix2polyvals` evaluate a saved `ctfdoc`.
- `ctf_estimate_iter%iterate` fits over `params%hp`/`params%lp` (`preprocess` sets them from
  `hp_ctf_estimate`/`lp_ctf_estimate` first), resolves `fit_phshift` itself and fills the caller's
  micrograph row (`dfx`, `dfy`, `df`, `angast`, `phshift`, `ctfres`, `icefrac`, `astig`, ...;
  `ctfdoc` with `ctfpatch=yes`). It never reads or writes the project.
- A fit reaches particles only through extraction: `extract` copies the micrograph row's CTF to the
  stack and every particle row (`add_stk`, `smpd` reset to `params%smpd`), then sets each particle's
  `dfx`/`dfy` from the `ctfdoc` polynomial at its centre (1-based micrograph pixels). 2D/3D read
  `kv`, `cs`, `fraca`, `smpd`, `ctf` from `os_stk` and defocus and phase from the particle row,
  never `os_optics`. Refitting after extraction changes no particle: re-extract.
- `ctf=no` zeroes defocus and phase in `get_ctfparams`; `ctf=flip` marks phase-flipped stacks;
  `extractfrommov=yes` refuses `flip` (`THROW_HARD`).

## Where to start

- `simple_ctf_estimate_iter.f90` (the order); `simple_ctf_estimate_fit.f90`: `srch` (1D grid with
  the phase grid), `refine` (differential evolution, then continuous with gradients), `fit_patches`.
- Callers: `strategies/parallelization/simple_preprocess_strategy.f90` (motion correction, then CTF
  on the same `ctfparams`), `simple_ctf_estimate_strategy.f90`, `simple_extract_strategy.f90`.

## Tests

- `build/bin/simple_test_exec test=unit_image suite=ctf` (`unit_image` under 1 s): closed-form values,
  conventions, phase policy, extrema, the image kernel, restoration contracts. The tester is here,
  registered under `unit_image`.
- No unit test fits a spectrum. `test=stream_preproc` gates the fitted defocus of five simulated
  movies within 0.20 microns of the truth (`ctfpatch=no`). `test=simulated_workflow suite=6vxx`
  (or `1jxy`) fits with patches but judges only the final map, a robustness profile: one run shows
  breakage, not improvement.
- Phase across folders: `test=unit_project suite=project_records` and `suite=star_project`,
  `test=unit_ui suite=ui_visibility`.

## See also

`doc/policies/phase_shift_ctf_policy.md`, `doc/policies/microscope_parameters_policy.md`,
`doc/algorithms/ctf_estimation.md`; skills `simple-main-image`, `simple-main-project`, `simple-main-motion`.
