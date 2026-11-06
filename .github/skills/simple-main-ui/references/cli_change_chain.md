# Adding or changing a command-line key or a program: the chain

Each row is one step, where it is done, and what happens when it is missed. "Silent" means the build,
the fast gate and the run all pass.

## A key

Ask Hans before adding any key (skill `simple-main-params`: every field is a key of every program).

| step | where | when missed |
| --- | --- | --- |
| field and default | `src/main/params/simple_parameters.f90`; `simple_args.f90` is regenerated from it | the key is unknown: `cmdline` stops (`THROW_HARD`) |
| a key that is no field, read only by a front door or the queue | the generator's own list in `scripts/simple_args_generator.pl` | as above |
| binding | `src/main/params/simple_parameters_parse.f90` (`add_yesno` for a `(yes|no)` key); a field `parameters` only derives says "not a command-line key" in its comment instead | loud: `scripts/check_cli_chain.py` fails the fast gate |
| field length | the declaration | silent: a longer value is cut |
| non-empty string default | `init_dynamic_defaults` (`_core.f90`) | silent: the field is empty, the UI shows no default |
| derivation, validation | `_phases.f90` | silent: no derived flag, bad values pass |
| UI input | the program's `add_input` (`src/main/ui/simple/` or `single/`), shared metadata in `simple_ui_params_common.f90`; an input that is no accepted key fails the gate (`scripts/check_cli_chain.py`) | silent for most programs (the key is accepted anyway, without help or GUI field); `solve3D_addon` and `solve3D_chunks` refuse it (`accepts`), and `solve3D_chunks` drops it when it forwards to the add-on |
| UI default | generated (`scripts/default_audit.py`): the declaration, `init_dynamic_defaults`, or a literal `if( .not. cline%defined(k) ) call cline%set(k, v)` in the routed commander's `execute` before `params%new` (the rule for a new default), the front door, the router case or a helper listed in `TRACED_DEFAULT_HELPERS` / `DELEGATED_WORKFLOWS` | silent: a default set elsewhere (after `new`, `strategies/`, `solve/`) is not shown; the UI shows the baseline |
| stream launcher JSON | `print_stream_ui_json` in `simple_ui.f90` (`prg=print_ui_stream`), written by hand | silent: the nice_lite stream launcher does not offer it; a key it names that is no accepted key fails the gate, a user input that is no input of its process or has another default fails the build (`validate_complete_ui_json.py`), and `validate_stream_ui_contract` compares it with `production/stream_launcher_contract.json` |
| child command lines | the controller or strategy that builds them | silent: the child runs with the default |
| required in a public program whose name a worker shares (`refine2D`, `refine3D`, ...) | the input's `required` flag | loud: the workers' command lines stop in `parse_private` |

## A program

| step | where | when missed |
| --- | --- | --- |
| descriptor | `ui_program%new` and `add_input` in the category module, `add_ui_program(name, ..., UI_CATEGORY)`; a new category is a new module called from its group module (`simple_ui_simple_group.f90` and its siblings) | loud: the front door stops on the name; a bad descriptor stops every program in `make_ui` |
| `executable`, `sp_required` | `ui_program%new` | `executable`: loud; `sp_required`: silent, it decides project discovery and the project copy under `mkdir=yes` |
| route | a case in `src/main/exec/<exe>_<topic>.f90` (simple_exec, single_exec); `production/simple_stream.f90` for simple_stream; `src/main/exec/simple_test_exec_*.f90` for tests | loud: "not recognized" after the router chain; for tests `scripts/check_test_registry.py` fails the fast gate |
| new router module | the executable's `src/main/apis/<exe>_api.f90` and the call chain in `production/<exe>.f90` | loud: compile error or "not recognized" |
| worker of `simple_private_exec` | a case in `dispatch_private_prg` (`production/simple_private_exec_driver.f90`), the commander in `src/main/apis/simple_private_exec_api.f90`, required keys in `src/main/ui/simple_private_prgs.f90` (`add_private_prg`, `push_req`, `push_opt`) unless it has a public program's name | dispatch: loud; keys: loud, `scripts/check_cli_chain.py` fails the gate on a worker without an entry (`add_private_prg`; the count follows) |
| JSON and NICE | follows from the descriptor (`print_ui_json`, which NICE reads through `simple_private_exec`); every build validates it (`validate_complete_ui_json`) | loud at build |
| stream screens | `cmake --build build --target validate_stream_ui_contract`, then `update_stream_ui_contract` after review; `print_ui_stream` by hand | contract (descriptors and launcher): loud when run; `print_ui_stream` against the descriptors: loud at build |
| test program | test UI (`src/main/ui/simple_test/`), one router case; for a suite its table and `suite=` help; CTest registration for a suite or a registry entry | loud: `check_test_registry.py` in the fast gate |
| rename | no alias; NICE rewrites stored names by a migration (`nice/nice_lite/migrations/0006_rename_programs.py`) | NICE's stored jobs keep the old name |
