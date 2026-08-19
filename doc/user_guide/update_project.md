# update_project

## NAME

`update_project` — update SIMPLE project metadata or write a new project with
relocated dataset paths

## SYNOPSIS

Metadata-only update:

```shell
simple_exec prg=update_project projfile=PROJECT.simple [COMPUTER_OPTION=value ...]
```

Global data-root relocation:

```shell
simple_exec prg=update_project \
  projfile=PROJECT.simple \
  old_root=OLD_ROOT \
  new_root=NEW_ROOT \
  [projfile_out=OUTPUT.simple] \
  [SCOPED_ROOT_PAIR ...] \
  [COMPUTER_OPTION=value ...]
```

Per-scope relocation:

```shell
simple_exec prg=update_project \
  projfile=PROJECT.simple \
  SCOPE_old_root=OLD_ROOT \
  SCOPE_new_root=NEW_ROOT \
  [projfile_out=OUTPUT.simple] \
  [COMPUTER_OPTION=value ...]
```

Valid scope prefixes are `mic`, `ptcl`, `cavg`, and `vol`.

## DESCRIPTION

`update_project` has two modes. The presence of any complete global or scoped
root pair selects relocation mode.

In metadata-only mode, the commander updates the project and computer
environment metadata in `projfile`. It does not load or rewrite the project
data segments.

In relocation mode, the commander reads the complete input project, replaces
supported dataset path prefixes, updates project and computer environment
metadata, and writes a new complete project. It does not modify the input
project.

The value of an old-root option is a prefix stored in the project. The old
directory does not need to remain accessible. The corresponding new root and
every remapped target must exist.

## REQUIRED OPTION

`projfile=PROJECT.simple`
: SIMPLE project to update. The filename must have the `.simple` extension.

## RELOCATION OPTIONS

`old_root=OLD_ROOT`
: Global old path prefix. It supplies the mapping for every scope that does not
  have an explicit scoped pair.

`new_root=NEW_ROOT`
: Existing global destination root. It must be supplied with `old_root`.

`projfile_out=OUTPUT.simple`
: Relocation output filename. Its parent directory must exist, its extension
  must be `.simple`, and it must not name an existing file. If omitted, the
  output is `<PROJECT>_remapped.simple` beside the input project. This option is
  invalid in metadata-only mode.

`mic_old_root=OLD_ROOT`, `mic_new_root=NEW_ROOT`
: Scoped mapping for movie, integrated-micrograph, and box paths.

`ptcl_old_root=OLD_ROOT`, `ptcl_new_root=NEW_ROOT`
: Scoped mapping for raw particle stacks, denoised particle stacks, and stack
  box paths.

`cavg_old_root=OLD_ROOT`, `cavg_new_root=NEW_ROOT`
: Scoped mapping for class-average stacks, stack directories, 2D FRC files,
  and sigma files.

`vol_old_root=OLD_ROOT`, `vol_new_root=NEW_ROOT`
: Scoped mapping for volumes, FSC files, and 3D FRC files.

Each old-root option requires its matching new-root option. A scoped pair
overrides the global pair for that scope. A scoped pair can also be used
without a global pair.

## PATH SCOPES

| Scope | Project segment and fields |
| --- | --- |
| Micrographs (`mic`) | `mic.movie`, `mic.intg`, `mic.boxfile` |
| Particles (`ptcl`) | `stk.stk`, `stk.stk_den`, `stk.boxfile` |
| Class averages (`cavg`) | `out.stk`, `out.stkpath`, `out.frcs` for `frc2D`, `out.sigma2` |
| Volumes (`vol`) | `out.vol`, `out.fsc`, `out.frcs` for `frc3D` |

Particle rows in the `ptcl2D` and `ptcl3D` segments refer to image data
through the `stk` segment. Relocating the particle scope therefore updates the
stack paths used by both particle segments.

## COMPUTER OPTIONS

These optional values update the project's computer environment metadata in
both operating modes.

| Option | Meaning |
| --- | --- |
| `user_email=ADDRESS` | Notification email address |
| `time_per_image=SECONDS` | Estimated processing time per image |
| `user_account=NAME` | Scheduler account name |
| `user_project=NAME` | Scheduler project name |
| `qsys_partition=NAME` | Scheduler partition or queue |
| `qsys_qos=NAME` | Scheduler quality-of-service or priority |
| `qsys_reservation=NAME` | Scheduler reservation name |
| `job_memory_per_task=MB` | Memory in MB per distributed part or computing node |
| `qsys_name=NAME` | Queue system: `local`, `coarray`, `slurm`, `pbs`, `lsf`, or `sge` |
| `walltime=SECONDS` | Maximum scheduler execution time |

## MATCHING AND VALIDATION

- Root matching occurs on a complete path-component boundary. For example,
  `/old/data/movie.mrc` matches `old_root=/old/data`, but
  `/old/data_backup/movie.mrc` does not.
- Empty roots, equal old and new roots, and filesystem roots used as an old
  root are rejected.
- The new root must be an existing directory.
- Every file or directory produced by a proposed mapping must exist.
- An explicitly supplied scoped pair must match at least one supported path in
  its scope.
- A global mapping can skip scopes that have no matching paths, but at least
  one supported path must match across all requested mappings.
- Validation failure prevents the output project from being written. The
  input project is not changed.

Path relocation applies only to the fields listed in [PATH SCOPES](#path-scopes).
Other strings in the project are not treated as dataset paths.

## EXAMPLES

Update scheduler metadata in the input project:

```shell
simple_exec prg=update_project \
  projfile=project.simple \
  qsys_name=slurm \
  qsys_partition=gpu \
  user_email=user@example.org
```

Relocate all supported dataset paths under one root:

```shell
simple_exec prg=update_project \
  projfile=project.simple \
  old_root=/old/storage/project \
  new_root=/new/storage/project
```

Write the relocated project to an explicit destination:

```shell
mkdir -p /work/relocated
simple_exec prg=update_project \
  projfile=/work/original/project.simple \
  old_root=/old/storage/project \
  new_root=/new/storage/project \
  projfile_out=/work/relocated/project.simple
```

Use independent roots for each data scope:

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

Use a global mapping with a micrograph-specific override:

```shell
simple_exec prg=update_project \
  projfile=project.simple \
  old_root=/old/project \
  new_root=/new/project \
  mic_old_root=/archive/micrographs \
  mic_new_root=/data/micrographs
```

Relocate the same scope from two old roots by chaining non-destructive updates:

```shell
simple_exec prg=update_project \
  projfile=project.simple \
  mic_old_root=/old/raw_movies \
  mic_new_root=/new/raw_movies \
  projfile_out=project_movies_remapped.simple

simple_exec prg=update_project \
  projfile=project_movies_remapped.simple \
  mic_old_root=/old/integrated_micrographs \
  mic_new_root=/new/integrated_micrographs \
  projfile_out=project_remapped.simple
```

## ENVIRONMENT

`SIMPLE_QSYS`
: Queue-system value used when `qsys_name` is not supplied. The commander
  reports an error if neither source provides a value.

`SIMPLE_EMAIL`
: Default email address when `user_email` is not supplied. If unset, SIMPLE
  uses its built-in placeholder address.

`SIMPLE_QSYS_PARTITION`
: Default scheduler partition when `qsys_partition` is not supplied.

## OUTPUT

In metadata-only mode, `projfile` is updated in place.

In relocation mode, the commander reports the number of remapped paths and
the output filename. The default output is `<PROJECT>_remapped.simple`; the
input project remains unchanged.

## EXIT STATUS

The command returns success after the requested metadata or relocation update
is written. Invalid option pairs, invalid output paths, unmatched explicit
scopes, missing targets, or an existing output file cause the command to
terminate with an error.

## SEE ALSO

- [SIMPLE User Guide](user_guide.md)
- [Relocating project data](../relocating_project_data.md)
