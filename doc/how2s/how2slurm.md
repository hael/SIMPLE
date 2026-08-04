# SLURM for beginners

SLURM is a cluster scheduler: you ask it for resources (CPUs, memory, time, and
optionally GPUs), and it starts your work when those resources are available. Do
not run compute-heavy work on a login node. Exact partitions, limits, modules,
and account names are site-specific; use your cluster's documentation or support
team for those values.

## Before your first submission

1. Log in to the cluster and find the local guidance (`module avail`, a cluster
   web page, or `sinfo`).
2. Run a tiny test first. Make sure its output, working directory, and software
   environment are correct before requesting long runs.
3. Choose the smallest realistic resource request. A job cannot start until all
   requested resources are free, and requests that are too small can be killed.

Useful discovery commands are:

```bash
sinfo                          # partitions and their current state
scontrol show partition        # limits and configuration (site-dependent output)
module avail                   # software modules, if your site uses Environment Modules
```

## Your first batch job

Save this as `hello.sbatch`, replacing the placeholder partition, account, and
module commands only when your site requires them.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=hello
#SBATCH --partition=<partition>       # omit if your site has a default
#SBATCH --account=<account>           # omit if not required
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0-00:05:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "Working directory: ${SLURM_SUBMIT_DIR:-$PWD}"

# module load <software>/<version>
# conda activate <environment>

srun hostname
echo "Finished: $(date)"
```

Create the `logs` directory **before** submitting—the scheduler opens the log
files before the script starts—then submit:

```bash
mkdir -p logs
sbatch hello.sbatch
```

`sbatch` returns a job ID, for example `Submitted batch job 123456`. Keep it:
it is the handle used to monitor, inspect, and cancel the job. `#SBATCH` lines
must appear before ordinary shell commands. Shell variables do **not** expand in
`#SBATCH` directives; use command-line options or a wrapper if you need dynamic
values.

## Requesting the right resources

For a normal multithreaded program, start with one task and request its thread
count with `--cpus-per-task`. Set the program's thread setting to the same value
(for example, `OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK` for OpenMP programs).

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
srun ./my_threaded_program input.dat
```

For an MPI program, request multiple tasks and launch it with `srun`; do not
assume `mpirun` has the right integration on a new cluster.

```bash
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1

srun ./my_mpi_program input.dat
```

For GPUs, consult local documentation for the required partition and GPU type.
A common request looks like this:

```bash
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
```

Memory is usually either a total job request (`--mem=32G`) or a per-CPU request
(`--mem-per-cpu=4G`). Use one deliberately: combining them can be confusing or
rejected by a site. Wall time accepts `days-hours:minutes:seconds`, such as
`1-12:00:00`.

## Where the job runs and where files go

Submit from the directory containing the input files, or explicitly change to it
in the script:

```bash
cd "$SLURM_SUBMIT_DIR"
```

`%j` in log paths becomes the job ID and `%x` becomes the job name. Prefer
separate output and error files while learning. Compute nodes may have fast,
temporary local storage; copy only disposable scratch data there and copy final
results back before the job exits. Never assume a compute-node filesystem is
preserved after a job ends.

## Monitor, investigate, and cancel

```bash
squeue -u "$USER"                         # queued and running jobs
squeue -j <JOB_ID>                        # one job
scontrol show job <JOB_ID>                # requested resources and pending reason
tail -f logs/<JOB_NAME>-<JOB_ID>.out      # follow standard output
tail -f logs/<JOB_NAME>-<JOB_ID>.err      # follow standard error
scancel <JOB_ID>                          # cancel one job
scancel -u "$USER"                        # cancel all of your jobs; use with care
```

Common `squeue` states: `PD` means pending (look at the reason column or
`scontrol show job`), `R` means running, `CG` means completing, and `F` means
failed. Once a job leaves the queue, review its accounting record:

```bash
sacct -j <JOB_ID> --format=JobID,JobName%30,State,ExitCode,Elapsed,AllocTRES,MaxRSS
```

`ExitCode=0:0` generally indicates success. `OUT_OF_MEMORY` means request more
memory or reduce use; `TIMEOUT` means request more time or make the work faster.
An exit code from the program itself calls for reading the job's error log.

## Useful variations

Submit with an override without editing the file:

```bash
sbatch --time=08:00:00 --job-name=test-long hello.sbatch
```

Run an interactive shell for short debugging (options vary by site):

```bash
salloc --time=00:30:00 --cpus-per-task=4 --mem=8G
srun --pty bash
```

Make a later job wait for a successful earlier job:

```bash
first=$(sbatch --parsable step1.sbatch)
sbatch --dependency=afterok:"$first" step2.sbatch
```

For many independent, similar jobs, use a job array; see
[the array guide](how2slurm_array.md). Avoid thousands of tiny jobs—bundle
short pieces of work into fewer tasks so the scheduler and your jobs spend less
time waiting.

## A practical checklist

- Verify input paths and write logs to a known location.
- Load software and activate environments inside the job script, not just in your
  login shell.
- Match program threads/ranks to the SLURM request and launch parallel work with
  `srun`.
- Start small, inspect `sacct` and logs, then scale up.
- Cancel jobs you no longer need and clean up only reproducible scratch files.

date
