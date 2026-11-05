# Adding a per-particle field (a named slot)

A key on a particle-typed record (`os_ptcl2D`, `os_ptcl3D`) is saved in the project only when it
is a named slot: `simple_binoris` writes particle segments as `ori2prec` (the slots) padded to
`N_PTCL_RECORD_REALS` reals, and drops the record's hash tables. A new per-particle field
therefore needs, all in `src/defs/simple_defs_ori.f90` unless stated:

1. An enumerator after the last one (`I_LAST_NAMED_ORIPARAM`); never renumber or reuse one. A
   retired field keeps its number and reports absent (`I_RETIRED_W`).
2. Its flag in `get_oriparam_ind` (flag -> slot) and `get_oriparam_flag` (slot -> flag).
3. `N_PTCL_ORIPARAMS` and `I_LAST_NAMED_ORIPARAM` raised. `N_PTCL_ORIPARAMS` must stay at most
   `N_PTCL_RECORD_REALS` (64): then the on-disk record width does not change and projects of
   earlier releases still read, with zeros in the new slot.
4. A case in `oriparam_isthere` only if a value of 0 means "not set" (as for class, proj, eo).
   Without one the slot always reports present.
5. The slot checks in `src/main/ori/simple_ori_tester.f90` (`test_cartesian_score_fields` pins
   each recent slot's number and flag); `src/main/project/simple_binoris_tester.f90` pins the
   record width.
6. If the field is exported or imported through STAR, the flag table in `src/main/star/`
   (skill `simple-main-star`).

Run `simple_test_exec test=unit_ori` and `test=unit_project suite=binoris`.
