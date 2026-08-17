!@descr: distributed-execution context for flex_pca stage fan-out
module simple_flex_pca_distr
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_qsys_env,   only: qsys_env
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_distr_init, flex_pca_distr_kill
public :: flex_pca_is_master, flex_pca_is_worker, flex_pca_nparts
public :: flex_pca_run_stage
public :: PCA_STAGE_SNR, PCA_STAGE_COLS, PCA_STAGE_SOLVE, PCA_STAGE_EMBED, PCA_STAGE_STATES
public :: PCA_STAGE_PROBE

! Stage selector carried to the worker in params%stage, exactly as flex_analysis does.
integer, parameter :: PCA_STAGE_SNR  = 1
integer, parameter :: PCA_STAGE_COLS = 2
integer, parameter :: PCA_STAGE_SOLVE = 3
integer, parameter :: PCA_STAGE_EMBED  = 4
integer, parameter :: PCA_STAGE_STATES = 5
! One qsys round per probe EM iteration: the basis changes every iteration, so workers must be
! re-launched against the master's refreshed flex_pca_pc*.mrc rather than looping locally.
integer, parameter :: PCA_STAGE_PROBE  = 6

! Master-side context. Distribution is injected at the level of the individual accumulate routines
! rather than by restructuring run_flex_pca into stage procedures: the master executes the ordinary
! single-node control flow and each distributable routine either loops in process (shared memory)
! or fans the particle range out and reduces the parts. That keeps ONE code path for the master-only
! algebra and the stage ordering, so the shared-memory result stays the reference.
type(qsys_env)   :: qenv
type(chash)      :: job_descr
logical          :: l_master = .false.
logical          :: l_worker = .false.
integer          :: nparts_run = 1

contains

    !> Called once by the commander. A run is a MASTER when nparts is given and part is not; a WORKER
    !! when part is given; otherwise plain shared memory.
    subroutine flex_pca_distr_init( params, cline )
        type(parameters), intent(inout) :: params
        class(cmdline),   intent(inout) :: cline
        integer :: fromp_glob, top_glob, nsel
        l_worker = cline%defined('part')
        l_master = cline%defined('nparts') .and. .not. l_worker
        if( l_master ) nparts_run = max(1, params%nparts)
        if( l_master .and. nparts_run > 1 )then
            ! qsys_env%new splits params%nptcls, i.e. the WHOLE project. flex_pca is routinely run on a
            ! window through fromp/top, and that window must be what gets partitioned -- otherwise the
            ! workers are handed the entire project and the master's selection is silently discarded.
            ! Split the window, then shift the ranges back into absolute project rows, which is what
            ! qsys_ctrl writes into each part's fromp/top and what sample4rec then selects within.
            fromp_glob = 1
            top_glob   = params%nptcls
            if( cline%defined('fromp') ) fromp_glob = params%fromp
            if( cline%defined('top')   ) top_glob   = params%top
            nsel = max(1, top_glob - fromp_glob + 1)
            nparts_run = min(nparts_run, nsel)
            call qenv%new(params, nparts_run, nptcls=nsel)
            if( fromp_glob > 1 ) qenv%parts = qenv%parts + (fromp_glob - 1)
            call cline%gen_job_descr(job_descr, prg=string('flex_pca'))
            call job_descr%set('mkdir',  'no')
            call job_descr%set('nparts', int2str(nparts_run))
            call job_descr%set('numlen', int2str(params%numlen))
            write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> DISTRIBUTED FLEX_PCA (MASTER), nparts=', &
                &nparts_run,' over project rows ',fromp_glob,'-',top_glob,' assigned through fromp/top'
            call flush(logfhandle)
        else
            l_master = .false.
            if( l_worker )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA (WORKER) part=',params%part,' stage=',params%stage
                call flush(logfhandle)
            endif
        endif
    end subroutine flex_pca_distr_init

    logical function flex_pca_is_master()
        flex_pca_is_master = l_master
    end function flex_pca_is_master

    logical function flex_pca_is_worker()
        flex_pca_is_worker = l_worker
    end function flex_pca_is_worker

    integer function flex_pca_nparts()
        flex_pca_nparts = nparts_run
    end function flex_pca_nparts

    !> One qsys round: every worker runs `stage` over its own fromp/top range and writes its part file.
    !! Returns when all parts are on disk. The caller then reduces them.
    subroutine flex_pca_run_stage( stage_id, label )
        integer,          intent(in) :: stage_id
        character(len=*), intent(in) :: label
        integer(timer_int_kind) :: t_round
        if( .not. l_master ) THROW_HARD('flex_pca_run_stage called outside the distributed master')
        t_round = tic()
        ! Clear the previous round's completion sentinels FIRST. qsys_ctrl decides a part is finished by
        ! the presence of JOB_FINISHED_<part>, so a second round started with the first round's sentinels
        ! still on disk concludes instantly and launches nothing -- the master then blocks on part files
        ! that were never written. Programs that schedule only one round never hit this.
        call del_files(JOB_FINISHED_FBODY, nparts_run)
        call del_files(SUBPROJECT_JOB_FINISHED_FBODY, nparts_run)
        call job_descr%set('stage', int2str(stage_id))
        write(logfhandle,'(A,A,A,I0,A)') '>>> FLEX_PCA distributing ',label,' over ',nparts_run,' parts'
        call flush(logfhandle)
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR)
        write(logfhandle,'(A,A,A,F8.1)') '>>> FLEX_PCA ',label,' qsys round seconds=',toc(t_round)
        call flush(logfhandle)
    end subroutine flex_pca_run_stage

    subroutine flex_pca_distr_kill( params )
        use simple_qsys_funs, only: qsys_cleanup
        type(parameters), intent(in) :: params
        if( l_master .and. nparts_run > 1 )then
            call qenv%kill
            call job_descr%kill
            call qsys_cleanup(params)
        endif
        l_master   = .false.
        l_worker   = .false.
        nparts_run = 1
    end subroutine flex_pca_distr_kill

end module simple_flex_pca_distr
