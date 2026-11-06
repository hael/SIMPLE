---
name: simple-dead-code
description: Use before calling SIMPLE code dead, unused, test-only or unreferenced, before making a module entity or type-bound procedure private, deleting or renaming it, and when asked who uses a module, procedure, type-bound procedure, constant or type. Covers what counts as test-only code, the reference forms grep misses or misreads in this code base, and the compiler probe (scripts/probe_entity_users.py) that proves who the users are.
---

# Finding the users of SIMPLE code

An entity has production users, only test users, or none. Grep finds candidates; only the compiler
proves an answer. Never report code as dead, unused or test-only from a search alone.

## Test-only code

What `src/CMakeLists.txt` drops without `BUILD_TESTS`: every `*_tester.f90`, `src/main/commanders/test/`
(the test commanders, the fixtures `simple_test_<x>_fixture`, the gate and registry),
`src/main/exec/simple_test_exec*.f90` and `src/main/apis/simple_test_exec_api.f90`. Everything else is
production, including `src/utils/simple_test_utils.f90`. A caller in the test front end is a test user,
not a production one.

## What a search misses or misreads

All seen in this code base:

- Component chains: `stage%sieve%get_n_chunks_coarse()` names no type. A `%name` search that also
  requires the type's name misses it; one that does not hits every type with a binding of that name
  (`%new`, `%kill`, `%init_params` of every stream stage).
- A binding's procedure called by its own name inside its module (`runs_to_final(self)`), not by `%`.
- Generics: a specific (`get_dp`, `getter_str`) is called through its generic (`generic :: get => ...`,
  an `interface` block, `operator(...)`, `assignment(=)`), never by name.
- Overrides: a binding that overrides a parent's is reached by calls on `class(parent)`.
- Imports: `use simple_core_module_api` and other unqualified uses re-export whole modules; renames
  (`local => remote`); names are case-insensitive; continuation lines; `!$` sentinel lines; macros
  (`THROW_HARD` is a call of `simple_exception`).
- Name-alikes: another module's procedure of the same name (`robust_spread` of `simple_test_profile_stats`).
- Procedure arguments and pointers (`add_suite(s, n, '...', run_all_x)`, `p => f`): a word search sees
  them, a `%` search does not.
- `bind(c)` procedures called from C/C++ (the starfile wrappers): invisible to the compiler as well.

## The probe

`scripts/probe_entity_users.py` renames entities (`<module>:<name>` or `<type>%<binding>`) in their
module family and in the files you accept as known users, backing every edited file up:

1. `python3 scripts/probe_entity_users.py apply . 'ptcl_sieve%get_n_chunks_coarse' 'simple_segmentation:canny'
   --accept src/main/sieve/simple_ptcl_sieve_tester.f90`. Read its warnings: a generic or a parent's
   class reaches the entity without naming it, so probe the generic or treat the entity as used.
2. Build with tests (the normal dev cycle). No edits to the tree until step 4.
3. `python3 scripts/probe_entity_users.py classify <build log>` lists each failing reference as
   `PRODUCTION` or `TEST`. A failing build shows some users, not all (make stops scheduling): accept or
   note them, `restore`, `apply` again, rebuild, until the build passes. An error on a component of the
   same name in another file is a name-alike, not a user.
4. `python3 scripts/probe_entity_users.py restore .` (before any other edit; `apply` refuses while a
   probe is pending).

A passing build proves that no file outside the family and the accepted files references the entities
(BUILD_TESTS=ON compiles production and tests). Report the run: "no reference outside its module and its
tester (probe, run N)"; a search result is a "candidate".

## What to do with the answer

- Production users: leave it public.
- Only its own tester: make it private and turn the tester into a submodule
  (`doc/policies/test_environment_policy.md` section 4.3).
- Only tests in several files: it stays public, or the tests move: say which.
- No production users: dead in production, whatever its tests. Policy deletes dead code with its tests;
  for the API of a general-purpose class (`image`, `oris`, `string`, the containers) the owner decides.
  List each entity with its probe run.

## See also

`doc/policies/test_environment_policy.md` (sections 4.1 and 4.3), `doc/policies/compile_time_policy.md`;
skills `simple-modern-fortran`, `simple-architecture`.
