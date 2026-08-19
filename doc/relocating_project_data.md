# Relocating project data

`update_project` can write a new project whose dataset paths point to files
that were moved under a different root directory:

```shell
simple_exec prg=update_project \
  projfile=project.simple \
  old_root=/old/data \
  new_root=/new/data
```

The global roots are used as a fallback for four path scopes:

| Scope | Project fields |
| --- | --- |
| Micrographs | `mic.movie`, `mic.intg`, `mic.boxfile` |
| Particles | `stk.stk`, `stk.stk_den`, `stk.boxfile` |
| Class averages | `out.stk`, `out.stkpath`, `out.frcs` (`frc2D`), `out.sigma2` |
| Volumes | `out.vol`, `out.fsc`, `out.frcs` (`frc3D`) |

`ptcl2D` and `ptcl3D` rows refer to raw or denoised images through the `stk`
segment, so particle remapping updates the stack paths used by both segments.

## Data stored under different roots

Each scope can have its own root pair:

```shell
simple_exec prg=update_project \
  projfile=project.simple \
  mic_old_root=/old/micrographs \
  mic_new_root=/new/micrographs \
  ptcl_old_root=/old/particles \
  ptcl_new_root=/new/particles \
  cavg_old_root=/old/class_averages \
  cavg_new_root=/new/class_averages \
  vol_old_root=/old/volumes \
  vol_new_root=/new/volumes
```

A scoped pair overrides the global `old_root` and `new_root` pair. This permits
one global mapping for most files and separate roots for exceptional scopes.
Every old-root argument must be accompanied by its matching new-root argument.
An explicitly supplied scope must match at least one stored project path.

Only paths equal to `old_root` or below it on a directory boundary are
rewritten. For example, `/old/data/movies/movie.eer` matches
`old_root=/old/data`, but `/old/data_backup/movie.eer` does not.

Before writing the output project, SIMPLE verifies that the new root and every
proposed file or directory exist. If validation fails, no output project is
written.

The input project is never modified. By default, the result is written beside
it as `<project>_remapped.simple`. Use `projfile_out` to choose another new
project filename:

```shell
simple_exec prg=update_project \
  projfile=project.simple \
  old_root=/old/data \
  new_root=/new/data \
  projfile_out=/projects/project_on_new_storage.simple
```

For safety, `projfile_out` must not already exist. If no root mapping is
supplied, `update_project` retains its original metadata-only behavior.
