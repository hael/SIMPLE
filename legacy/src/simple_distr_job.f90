module simple_distr_job
use simple_defs ! singleton
implicit none

public :: distr_job
private

type distr_job
    private
end type distr_job

! ### USER PARAMETERS
! user_account
! user_email
! user_project

! ### QSYS PARAMETERS
! qsys_partition
! qsys_queue
! qsys_qos

! ### JOB PARAMETERS
! job_name
! job_exec_dir
! job_cmd_str
! job_outdoc
! job_output
! job_ntasks
! job_ntasks_per_socket
! job_cpus_per_task
! job_gpus_per_task
! job_memory
! job_time
! job_time_days
! job_time_hrs
! job_time_secs
! job_part
! job_fromp
! job_top

end module simple_distr_job
