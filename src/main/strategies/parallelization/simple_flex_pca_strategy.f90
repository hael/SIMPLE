!@descr: flex_pca execution strategies: shared memory, distributed master, distributed worker
!!
!! The strategy owns the ROLE of a run and everything that follows from it: the qsys context,
!! the partition of the particle selection, the rounds and their control state, and the
!! worker's stage dispatch. Domain modules (model, em, rec3D) own the numerics; they receive the
!! strategy as a `flex_pca_rounds` (simple_flex_pca_rounds) and never query their role
!! otherwise. Same lifecycle as rec3D/refine3D: initialize -> execute -> finalize_run ->
!! cleanup. The flex_pca commander selects shared memory or the distributed master
!! (nparts>1); the worker commander runs the worker strategy. The master and shared-memory
!! roles settle the canonical sigma2 state before any part command line exists.
module simple_flex_pca_strategy
use simple_core_module_api, only: arr2txtfile, chash, del_file, int2str, int2str_pad, L_USE_SLURM_ARR, &
    &logfhandle, simple_exception, simple_mkdir, simple_rmdir, string, tic, timer_int_kind, toc, TXT_EXT
use simple_builder,         only: builder
use simple_cmdline,         only: cmdline
use simple_parameters,      only: parameters
use simple_qsys_env,        only: qsys_env
use simple_flex_pca_rounds,    only: flex_pca_rounds, flex_pca_rounds_shmem
use simple_flex_pca_stages,    only: flex_stage_request, FLEX_FIT_ALL, &
    &PCA_STAGE_PROBE, PCA_STAGE_POLISH, PCA_STAGE_EMBED
use simple_flex_pca_artifacts, only: flex_pca_local_part_dir
use simple_flex_pca_application, only: flex_pca_application
use simple_flex_pca_project_gateway, only: ensure_canonical_sigma_state
use simple_exec_helpers,         only: set_master_num_threads
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_strategy, flex_pca_worker_strategy, create_flex_pca_strategy

!> The lifecycle. Every strategy is also a flex_pca_rounds, which is what the domain code sees.
type, abstract :: flex_pca_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
end type flex_pca_strategy

type, extends(flex_pca_strategy) :: flex_pca_shmem_strategy
    type(flex_pca_rounds_shmem) :: rounds
contains
    procedure :: initialize   => shmem_initialize
    procedure :: execute      => shmem_execute
    procedure :: finalize_run => shmem_finalize_run
    procedure :: cleanup      => shmem_cleanup
end type flex_pca_shmem_strategy

type, extends(flex_pca_strategy) :: flex_pca_worker_strategy
    type(flex_pca_rounds_shmem) :: rounds
contains
    procedure :: initialize   => worker_initialize
    procedure :: execute      => worker_execute
    procedure :: finalize_run => worker_finalize_run
    procedure :: cleanup      => worker_cleanup
end type flex_pca_worker_strategy

!> The distributed master's qsys context: qenv, job_descr, one part_params per part (the
!! part's particle list as pindfile=), the master thread budget.
type, extends(flex_pca_rounds) :: flex_pca_master_rounds
    type(qsys_env)           :: qenv
    type(chash)              :: job_descr
    type(chash), allocatable :: part_params(:)
    integer                  :: nthr_master = 1   !< the master process's shared-memory budget (set_master_num_threads)
    integer                  :: nthr_worker = 1   !< the per-worker nthr of the command line, which the part headers request
contains
    procedure :: plan_partitions => master_plan_partitions
    procedure :: run_stage       => master_run_stage
end type flex_pca_master_rounds

type, extends(flex_pca_strategy) :: flex_pca_master_strategy
        character(len=:), allocatable :: part_dir
    type(flex_pca_master_rounds) :: rounds
contains
    procedure :: initialize   => master_initialize
    procedure :: execute      => master_execute
    procedure :: finalize_run => master_finalize_run
    procedure :: cleanup      => master_cleanup
end type flex_pca_master_strategy

abstract interface
    subroutine init_interface( self, params, build, cline )
        import :: flex_pca_strategy, parameters, builder, cmdline
        class(flex_pca_strategy), intent(inout) :: self
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
    end subroutine init_interface
    subroutine exec_interface( self, params, build, cline )
        import :: flex_pca_strategy, parameters, builder, cmdline
        class(flex_pca_strategy), intent(inout) :: self
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
    end subroutine exec_interface
    subroutine finalize_interface( self, params, build, cline )
        import :: flex_pca_strategy, parameters, builder, cmdline
        class(flex_pca_strategy), intent(inout) :: self
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
    end subroutine finalize_interface
    subroutine cleanup_interface( self, params, build, cline )
        import :: flex_pca_strategy, parameters, builder, cmdline
        class(flex_pca_strategy), intent(inout) :: self
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
    end subroutine cleanup_interface
end interface

contains

    !> The flex_pca commander's role from the command-line shape: nparts>1 makes the distributed
    !! master, anything else runs in shared memory (a part is the worker commander's)
    function create_flex_pca_strategy( cline ) result( strategy )
        class(cmdline), intent(in) :: cline
        class(flex_pca_strategy), allocatable :: strategy
        integer :: nparts
        if( cline%defined('part') ) THROW_HARD('a flex_pca part runs through the worker commander (simple_private_exec)')
        nparts = 1
        if( cline%defined('nparts') ) nparts = max(1, cline%get_iarg('nparts'))
        if( nparts > 1 )then
            allocate(flex_pca_master_strategy :: strategy)
        else
            allocate(flex_pca_shmem_strategy :: strategy)
        endif
    end function create_flex_pca_strategy

    ! ------------------------------------------------------------------ shared memory

    subroutine shmem_initialize( self, params, build, cline )
        class(flex_pca_shmem_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call ensure_canonical_sigma_state(params, build, cline)
        ! a command-line pcafit pins the whole job to one independent half (two-job halfset fits)
        self%rounds%fit_sel = params%pcafit
    end subroutine shmem_initialize

    subroutine shmem_execute( self, params, build, cline )
        class(flex_pca_shmem_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(flex_pca_application) :: app
        call app%run(params, build, cline, self%rounds)
    end subroutine shmem_execute

    subroutine shmem_finalize_run( self, params, build, cline )
        class(flex_pca_shmem_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
    end subroutine shmem_finalize_run

    subroutine shmem_cleanup( self, params, build, cline )
        class(flex_pca_shmem_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
    end subroutine shmem_cleanup

    ! ------------------------------------------------------------------ distributed master

    subroutine master_initialize( self, params, build, cline )
        class(flex_pca_master_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer :: fromp_glob, top_glob, nsel
        ! the master's budget is reported as every distributed master's; the master's own stages run
        ! at params%nthr, which FLEX sizes its per-thread scratch from
        call set_master_num_threads(self%rounds%nthr_master, string('FLEX_PCA'))
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        ! the canonical sigma2 state, and the grouping the workers load it under, before job_descr
        call ensure_canonical_sigma_state(params, build, cline)
        self%rounds%l_master   = .true.
        self%rounds%nparts_run = max(1, params%nparts)
        ! a command-line pcafit pins the whole job to one independent half; every scheduled
        ! round stamps the same selector
        if( params%pcafit /= FLEX_FIT_ALL ) self%rounds%fit_sel = params%pcafit
        ! the partition itself is planned once the master's particle selection exists
        ! (plan_partitions): each part receives its own index list, so the workers never
        ! re-derive the selection from fromp/top
        fromp_glob = 1
        top_glob   = params%nptcls
        if( cline%defined('fromp') ) fromp_glob = params%fromp
        if( cline%defined('top')   ) top_glob   = params%top
        nsel = max(1, top_glob - fromp_glob + 1)
        self%rounds%nparts_run = min(self%rounds%nparts_run, nsel)
        call cline%gen_job_descr(self%rounds%job_descr, prg=string('flex_pca'))
        call self%rounds%job_descr%set('mkdir',  'no')
        ! node-local part directory (local queue system + cache_dir): created here, adopted by
        ! every worker from the same derivation, removed by master_cleanup
        self%part_dir = flex_pca_local_part_dir(params, self%rounds%nparts_run)
        if( len_trim(self%part_dir) > 0 )then
            call simple_mkdir(self%part_dir)
            call self%rounds%set_part_dir(self%part_dir)
            write(logfhandle,'(A,A)') '>>> DISTRIBUTED FLEX_PCA: part files on local scratch ', trim(self%part_dir)
            call flush(logfhandle)
        endif
        call self%rounds%job_descr%set('nparts', int2str(self%rounds%nparts_run))
        call self%rounds%job_descr%set('numlen', int2str(params%numlen))
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> DISTRIBUTED FLEX_PCA (MASTER), nparts=', &
            &self%rounds%nparts_run,' over project rows ',fromp_glob,'-',top_glob
        call flush(logfhandle)
        self%rounds%nthr_worker = params%nthr
    end subroutine master_initialize

    subroutine master_execute( self, params, build, cline )
        class(flex_pca_master_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(flex_pca_application) :: app
        ! the master runs the shared-memory control flow; each distributable phase fans its
        ! particle partition out as one qsys round through self%rounds and reduces the parts
        call app%run(params, build, cline, self%rounds)
    end subroutine master_execute

    subroutine master_finalize_run( self, params, build, cline )
        class(flex_pca_master_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
    end subroutine master_finalize_run

    subroutine master_cleanup( self, params, build, cline )
        use simple_qsys_funs, only: qsys_cleanup
        class(flex_pca_master_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer :: ipart
        if( allocated(self%rounds%part_params) )then
            do ipart = 1, size(self%rounds%part_params)
                call del_file(self%rounds%part_params(ipart)%get('pindfile'))
                call self%rounds%part_params(ipart)%kill
            end do
            deallocate(self%rounds%part_params)
        endif
        call self%rounds%qenv%kill
        call self%rounds%job_descr%kill
        call qsys_cleanup(params)
        if( len_trim(self%part_dir) > 0 )then
            call simple_rmdir(self%part_dir)   ! consumers delete every part after the reduce; the directory is empty
            call self%rounds%kill_artifacts
        endif
    end subroutine master_cleanup

    !> Partition the master's particle selection: one contiguous slice of pinds per part, written
    !! to flex_pca_particles_part<NN>.txt and handed to that part as pindfile= through part_params.
    !! Called once, after validation.
    subroutine master_plan_partitions( self, params, pinds )
        class(flex_pca_master_rounds), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        integer,                       intent(in)    :: pinds(:)
        type(string) :: fname
        integer :: ipart, first, last, nsel, numlen
        nsel = size(pinds)
        if( nsel < 1 ) THROW_HARD('flex_pca plan_partitions: empty particle selection')
        self%nparts_run = min(self%nparts_run, nsel)
        ! qsys_nthr: the part scripts request the worker thread count, never the master's budget,
        ! which would inflate every part's CPU request against the scheduler's per-user cap
        call self%qenv%new(params, self%nparts_run, numlen=params%numlen, nptcls=nsel, qsys_nthr=self%nthr_worker)
        numlen = max(params%numlen, len(int2str(self%nparts_run)))
        if( allocated(self%part_params) ) deallocate(self%part_params)
        allocate(self%part_params(self%nparts_run))
        do ipart = 1, self%nparts_run
            first = self%qenv%parts(ipart,1)
            last  = self%qenv%parts(ipart,2)
            fname = string('flex_pca_particles_part')//int2str_pad(ipart,numlen)//TXT_EXT
            call arr2txtfile(pinds(first:last), fname)
            call self%part_params(ipart)%new(1)
            call self%part_params(ipart)%set('pindfile', fname%to_char())
            call fname%kill
        end do
        call self%job_descr%set('nparts', int2str(self%nparts_run))
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA partitioned ',nsel,' selected particles over ', &
            &self%nparts_run,' parts (one index list per part)'
        call flush(logfhandle)
    end subroutine master_plan_partitions

    !> One qsys round. The request's control state travels in job_descr under registered keys
    !! (the refine3D idiom): stage, pcafit, which_iter (the global EM iteration), maxits (the
    !! budget) and nfits (1 single fit, 2 paired halves). extra_params=params has the scheduler
    !! clear the previous round's sentinels itself.
    subroutine master_run_stage( self, params, req )
        class(flex_pca_master_rounds), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        type(flex_stage_request),      intent(in)    :: req
        integer(timer_int_kind) :: t_round
        if( .not. allocated(self%part_params) ) THROW_HARD('flex_pca run_stage before plan_partitions')
        t_round = tic()
        call self%job_descr%set('stage',      int2str(req%stage))
        call self%job_descr%set('pcafit',     int2str(self%fit_sel))
        call self%job_descr%set('which_iter', int2str(req%which_iter))
        call self%job_descr%set('maxits',     int2str(req%maxits))
        call self%job_descr%set('nfits',      int2str(req%nfits))
        write(logfhandle,'(A,A,A,I0,A)') '>>> FLEX_PCA distributing ',trim(req%label),' over ',self%nparts_run,' parts'
        call flush(logfhandle)
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, part_params=self%part_params, &
            &array=L_USE_SLURM_ARR, extra_params=params)
        write(logfhandle,'(A,A,A,F8.1)') '>>> FLEX_PCA ',trim(req%label),' qsys round seconds=',toc(t_round)
        call flush(logfhandle)
    end subroutine master_run_stage

    ! ------------------------------------------------------------------ distributed worker

    !> A worker runs in the master's directory (mkdir=no), takes its particle list (pindfile=) and
    !! its round state (stage, which_iter, maxits, nfits, pcafit) from the master's job_descr, and
    !! never re-derives a master decision.
    subroutine worker_initialize( self, params, build, cline )
        class(flex_pca_worker_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        if( .not. cline%defined('part') )     THROW_HARD('part= must be defined for flex_pca worker execution')
        if( .not. cline%defined('stage') )    THROW_HARD('stage= must be defined for flex_pca worker execution')
        if( .not. cline%defined('pindfile') ) THROW_HARD('pindfile= (the part''s particle list) must be defined for flex_pca worker execution')
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        select case(params%stage)
            case(PCA_STAGE_PROBE, PCA_STAGE_POLISH, PCA_STAGE_EMBED)
            case default
                THROW_HARD('invalid flex_pca worker stage')
        end select
        self%rounds%l_worker = .true.
        self%rounds%fit_sel  = params%pcafit
        ! adopt the master's node-local part directory (same derivation, same node, same cwd)
        call self%rounds%set_part_dir(flex_pca_local_part_dir(params, params%nparts))
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA (WORKER) part=',params%part, &
            &' stage=',params%stage,' fit=',self%rounds%fit_sel
        call flush(logfhandle)
    end subroutine worker_initialize

    subroutine worker_execute( self, params, build, cline )
        use simple_qsys_funs, only: qsys_declare_part_finished
        class(flex_pca_worker_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(flex_pca_application) :: app
        call app%run_worker(params, build, cline, self%rounds)
        call qsys_declare_part_finished(params, string('simple_flex_pca_strategy :: worker_execute'))
    end subroutine worker_execute

    subroutine worker_finalize_run( self, params, build, cline )
        class(flex_pca_worker_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
    end subroutine worker_finalize_run

    subroutine worker_cleanup( self, params, build, cline )
        class(flex_pca_worker_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
    end subroutine worker_cleanup

end module simple_flex_pca_strategy
