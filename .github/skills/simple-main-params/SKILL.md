---
name: simple-main-params
description: Use when working in SIMPLE's src/main/params subsystem, including the parameters type, command-line parameter binding, the generated list of accepted keys (simple_args.f90), typed parsing, derived parameter phases, dynamic defaults, filt_mode flags, validation of parsed workflow state, the process-wide side effects of params%new (seed, threads, execution directory), and adding or changing a command-line parameter.
---

# `src/main/params`: the parameters object and its parse, derive and validate phases

Owns `type(parameters)` and `params%new(cline)`: command-line values bound onto fields, values
derived from them, the program's UI record and the project's first records, and checks between keys.
Not here: the accepted-key list (generated), defaults a workflow sets on `cline` before `new`
(commanders, `solve/`, `strategies/`), the UI's display of defaults (`ui/`), `src/main/simple_cmdline.f90`.

## Rules

Silent ones first: nothing reports these.

- Every scalar declared in the type becomes an accepted key of every program, unless its comment says
  "not a command-line key" (a field `parameters` derives or sets itself): the build generates
  `simple_args.f90` from the declarations up to the type's `contains` and `cmdline` checks keys against
  it, not against the program. So no parameter for an internal mechanism or a test switch: decide
  from state the run already has (a key already given, the project, `which_iter`), force the other
  path of an A/B identity run with an existing key (`ncunits=1` runs a distributed iteration's parts
  in turn), and ask Hans before adding any parameter. The `internal:` keys in the type (controller-to-child
  settings) are ruled to be retired: add none, build nothing new on them.
- A declared key that `simple_parameters_parse.f90` does not bind would be accepted and ignored:
  `scripts/check_cli_chain.py` (fast gate) fails on it, on a key bound twice, on a UI input that is no
  accepted key and on a yes/no key not bound with `add_yesno`. A new field is bound, or marked "not a
  command-line key"; `nrots` is the one unbound accepted key, which a program reads from its command
  line itself (`KEYS_READ_OUTSIDE` in the check).
- `add_char` writes into the declared length and cuts a longer value (`len=3` makes `yess` `yes`).
- Bind files with `add_file`, directories with `add_dir`: they are made absolute and written back to
  `cline` before `mkdir=yes` changes directory (image files also set `ext`). A path bound with
  `add_char` stays relative to the launch directory.
- Presence is not value: phases test `cline%defined(key)` (`lp` sets `l_lpset`; `hp`, `box_crop`,
  `mskdiam`, `nthr`, `state`, `bfac`, `trs`, `nspace`; `part` makes the run a distributed worker).
  A child command line that restates a default, or keeps a stale key, changes the child. `smpd`,
  `box`, `nptcls`, `ctf` come from the project (or the stack/volume) unless given.
- Values are normalised without a word: `update_frac` above 0.99 becomes 1 and turns off
  `trail_rec`, `fillin`, `assign_all`; `lp`, `lpstart`, `lpstop` are raised to Nyquist.
- A key declared with the choices `(yes|no)` is bound with `add_yesno`, which stops on any other value
  (before `add_char` would cut `yess` to `yes`); the choices in the comment are what the gate check
  reads, and a yes/no (`binary`) UI input whose declaration lacks them fails it. A yes/no key with
  neither the choices nor a UI input is still compared with `'yes'` unchecked: give it the choices.
- `new` changes the process (layering policy section 4): `seed_rnd` reseeds (under `SIMPLE_SEED`
  each call takes the next seed of a count, so one `new` more or less shifts every later draw of a
  seeded run); it sets the OpenMP team, `nthr_glob`, `CWD_GLOB`, `part_glob`, `l_distr_worker_glob`;
  with `mkdir=yes` it creates and enters the execution directory and, when the program's UI says
  `sp_required`, copies the project there. It rewrites `cline` (`projfile`, absolute paths).
- Phase order is in `new` (`_phases.f90`). Validation also derives (`filt_mode` flags, `l_neigh`,
  `trs`, `l_shift_map`, `l_centring`, which clears `l_shift_prior`): set a flag after the last line
  that changes what it reads. A test's fixed seed goes after `new`.
- Loud: an unknown key stops in `cmdline`, a value a validation `select case` refuses stops
  (`THROW_HARD`); a too-large or zero `mskdiam` falls back with `THROW_WARN`.

## Adding or changing a key

1. Declare it with its default in `simple_parameters.f90` (a non-empty `type(string)` default in
   `init_dynamic_defaults`); the build regenerates `simple_args.f90` (never edit that file).
2. Bind it in the matching group of `simple_parameters_parse.f90` (`add_yesno` for a `(yes|no)` key);
   a field `parameters` only derives says "not a command-line key" in its comment instead.
3. Derive it in the phase that has its inputs; refuse bad values in validation.
4. Expose it in its programs (skill `simple-main-ui`); then read `params%<key>`, not `cline`.

Renamed keys get no alias. The whole chain and the steps that fail silently:
`.github/skills/simple-main-ui/references/cli_change_chain.md`.

## Where to start

`_registry.f90` (what each binding kind does), `_phases.f90` (`new`, the phases); callers of `new`:
commanders and `builder%init_params_and_build_*` (`src/main/simple_builder.f90`).

## Tests

- `build/bin/simple_test_exec test=unit_core suite=parameters` (a parsed key survives `new`, the accepted
  keys, a yes/no key refusing another value in a forked child) and `suite=command_line`; all of
  `unit_core` takes under a second. `python3 scripts/check_cli_chain.py` checks the chain statically.
- A changed derivation: `test=units` (all fifteen area suites), then the high-level entry of each
  workflow that reads the value. A change meant to alter nothing: two runs of that entry and
  `test=identity target=<entry> dir1=<run> dir2=<run>`.

## See also

`doc/policies/layering_and_project_ownership_policy.md` (sections 4, 6), `doc/policies/ui_layer_policy.md`
(defaults); skills `simple-main-ui`, `simple-main-root` (cmdline, builder).
