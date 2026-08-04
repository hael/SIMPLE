# SLURM job arrays

Job arrays run many independent copies of one batch script. Start with the
[general SLURM guide](how2slurm.md) if you have not submitted an ordinary job
yet.

## Submit an array

A job array can be submitted simply by adding to the job script:

```bash
#SBATCH --array=x-y
```

where `x` and `y` are the inclusive array bounds. You can also provide a
comma-separated list of task numbers, for example, to rerun selected failed
tasks:

```bash
sbatch --array=4,8,15,16,23,42  job_script.sbatch
```

Limit the number of simultaneously active tasks with `%N`. This avoids consuming
all of your allocation or overloading a shared filesystem:

```bash
#SBATCH --array=1-200%10
```

## Naming output and error files

SLURM uses the `%A` and `%a` replacement strings for the master job ID and task ID, respectively.

```bash
#SBATCH --output=Array_test.%A_%a.out
#SBATCH --error=Array_test.%A_%a.error

! get the SLURM array task id from Fortran
CHARACTER(len=255) :: task_id
CALL get_environment_variable("SLURM_ARRAY_TASK_ID", task_id)
```

## Deleting job arrays and tasks

To delete all of the tasks of an array job, use scancel with the job ID:

```bash
scancel 292441
```

To delete a single task, add the task ID:

```bash
scancel 292441_5
```

Job arrays set `SLURM_ARRAY_JOB_ID` to the array's master job ID and
`SLURM_ARRAY_TASK_ID` to the current task index. They also set
`SLURM_ARRAY_TASK_COUNT`, `SLURM_ARRAY_TASK_MAX`, and `SLURM_ARRAY_TASK_MIN`.
Use `squeue -j <ARRAY_JOB_ID>` to see its tasks, `scontrol show job
<ARRAY_JOB_ID>_<TASK_ID>` for one task, and `sacct -j <ARRAY_JOB_ID>` after it
finishes. General monitoring and troubleshooting are in
[how2slurm.md](how2slurm.md).

## Batch many short runs per array task

```bash
#!/usr/bin/env bash
#SBATCH --job-name=mega_array       # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@example.org   # Where to send mail
#SBATCH --nodes=1                   # One node per array task
#SBATCH --ntasks=1                  # One task per array element
#SBATCH --mem-per-cpu=1G            # Memory per CPU
#SBATCH --time=00:10:00             # Wall time per array task
#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=1-5                 # 5 tasks × 1000 runs = 5000 runs
# This combines array tasks with a loop to process many short runs. This is
# friendlier to the scheduler than submitting one array task per very short run.
pwd; hostname; date

# Set the number of runs that each SLURM task should do.
PER_TASK=1000

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
  echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
  # Do your work for "$run" here.
done
```

date
