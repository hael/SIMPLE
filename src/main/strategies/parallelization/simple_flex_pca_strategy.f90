!@descr: flex_pca execution strategies: shared memory, distributed master, distributed worker
!!
!! The strategy owns the ROLE of a run and everything that follows from it: the qsys context,
!! the partition of the particle selection, the rounds and their control state, and the
!! worker's stage dispatch. Domain modules (model, em, rec3D) own the numerics; they receive the
!! strategy as a `flex_pca_rounds` (simple_flex_pca_rounds) and never query their role
!! otherwise. Same lifecycle as rec3D/refine3D: initialize -> execute -> finalize_run ->
!! cleanup, selected by the command-line shape (part= -> worker; nparts>1 -> master; else
!! shared memory).
module simple_flex_pca_strategy
use simple_core_module_api
use simple_builder,         only: builder
use simple_cmdline,         only: cmdline
use simple_parameters,      only: parameters
use simple_qsys_env,        only: qsys_env
use simple_flex_pca_rounds, only: flex_pca_rounds, flex_pca_rounds_shmem, FLEX_FIT_ALL, &
    &flex_pca_set_part_dir, flex_pca_local_part_dir, &
    &PCA_STAGE_PROBE, PCA_STAGE_POLISH, PCA_STAGE_EMBED, PCA_STAGE_STATES
use simple_flex_pca_model,  only: run_flex_pca, run_flex_pca_worker
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_strategy, create_flex_pca_strategy

!> The lifecycle. Every strategy is also a flex_pca_rounds, which is what the domain code sees.
type, abstract :: flex_pca_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
    !> a key the master decided after initialize (the sigma fallback) that every worker must carry;
    !! a no-op for the shared-memory and worker roles
    procedure :: set_worker_key => strategy_set_worker_key_noop
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
    integer                  :: nthr_master = 1
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
    procedure :: set_worker_key => master_set_worker_key
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

    !> Role from the command-line shape, as the July flex_analysis factory: part= makes a worker;
    !! nparts>1 without part= makes a master; anything else runs in shared memory.
    function create_flex_pca_strategy( cline ) result( strategy )
        class(cmdline), intent(in) :: cline
        class(flex_pca_strategy), allocatable :: strategy
        integer :: nparts
        logical :: is_worker, is_master
        nparts = 1
        if( cline%defined('nparts') ) nparts = max(1, cline%get_iarg('nparts'))
        is_worker = cline%defined('part')
        is_master = nparts > 1 .and. .not. is_worker
        if( is_master )then
            allocate(flex_pca_master_strategy :: strategy)
        else if( is_worker )then
            allocate(flex_pca_worker_strategy :: strategy)
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
        ! a command-line pcafit pins the whole job to one independent half (two-job halfset fits)
        self%rounds%fit_sel = params%pcafit
    end subroutine shmem_initialize

    subroutine shmem_execute( self, params, build, cline )
        class(flex_pca_shmem_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        call run_flex_pca(params, build, cline, self%rounds)
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

    subroutine strategy_set_worker_key_noop( self, key, val )
        class(flex_pca_strategy), intent(inout) :: self
        character(len=*),         intent(in)    :: key, val
    end subroutine strategy_set_worker_key_noop

    !> the part scripts are generated from job_descr, captured from the command line at initialize;
    !! a decision the master takes later (sigma_est=global by the fallback) has to be written here or
    !! the workers load the canonical state under the wrong grouping policy and die at once
    subroutine master_set_worker_key( self, key, val )
        class(flex_pca_master_strategy), intent(inout) :: self
        character(len=*),                intent(in)    :: key, val
        call self%rounds%job_descr%set(key, val)
    end subroutine master_set_worker_key

    subroutine master_initialize( self, params, build, cline )
        use omp_lib, only: omp_set_num_threads, omp_get_num_procs
        class(flex_pca_master_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer :: fromp_glob, top_glob, nsel, vovr, ncpu_own
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
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
            call flex_pca_set_part_dir(self%part_dir)
            write(logfhandle,'(A,A)') '>>> DISTRIBUTED FLEX_PCA: part files on local scratch ', trim(self%part_dir)
            call flush(logfhandle)
        endif
        call self%rounds%job_descr%set('nparts', int2str(self%rounds%nparts_run))
        call self%rounds%job_descr%set('numlen', int2str(params%numlen))
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> DISTRIBUTED FLEX_PCA (MASTER), nparts=', &
            &self%rounds%nparts_run,' over project rows ',fromp_glob,'-',top_glob
        call flush(logfhandle)
        ! ---- MASTER THREAD BOOST (profiled 2026-09-02: master-only phases own 81% of wall at
        ! 1.5 busy threads while every worker core idles). Master and worker phases NEVER overlap
        ! (one qsys round per iteration; the master blocks on part files), so the master uses
        ! the workers' full thread budget for its own stages. Expressed the refine3D way: the
        ! budget is a context component applied through omp_set_num_threads; job_descr and the
        ! cline already carry the PER-WORKER nthr. The master's OWN params must agree with its
        ! OpenMP budget because builder/matcher/reconstructor scratch is sized from params%nthr
        ! (nthr_glob) -- the same contract refine3D keeps by handing nthr_master to every
        ! sub-command line it runs in-process. SIMPLE_COV_MASTER_NTHR overrides.
        ! Capped at the cores this master process actually owns: under a scheduler that is the
        ! job's CPU allocation (SLURM_CPUS_PER_TASK), else the machine (omp_get_num_procs). Without
        ! the cap, params%nthr became nparts*nthr and the part scripts the master generates inherited
        ! it (--cpus-per-task=160 on a 96-core partition: unschedulable, five retries, master dead;
        ! verification 2026-09-16). No environment override: the cap is the whole policy.
        self%rounds%nthr_master = max(params%nthr, self%rounds%nparts_run*params%nthr)
        ncpu_own = omp_get_num_procs()
        vovr = 0
        call get_env_int_local('SLURM_CPUS_PER_TASK', vovr)
        if( vovr > 0 ) ncpu_own = vovr
        self%rounds%nthr_master = max(params%nthr, min(self%rounds%nthr_master, ncpu_own))
        call omp_set_num_threads(self%rounds%nthr_master)
        if( self%rounds%nthr_master > params%nthr )then
            params%nthr = self%rounds%nthr_master
            nthr_glob   = self%rounds%nthr_master
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA MASTER THREAD BOOST: master-only stages run at ', &
                &self%rounds%nthr_master, ' threads (worker phases never overlap)'
            call flush(logfhandle)
        endif
    end subroutine master_initialize

    subroutine master_execute( self, params, build, cline )
        class(flex_pca_master_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        ! the master runs the shared-memory control flow; each distributable phase fans its
        ! particle partition out as one qsys round through self%rounds and reduces the parts
        call run_flex_pca(params, build, cline, self%rounds)
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
            call flex_pca_set_part_dir('')
        endif
    end subroutine master_cleanup

    !> Partition the master's particle selection: one contiguous slice of pinds per part, written
    !! to flex_pca_particles_part<NN>.txt and handed to that part as pindfile= through part_params
    !! (the July flex_analysis prepare_particle_partitions). Called once, after validation.
    subroutine master_plan_partitions( self, params, pinds )
        class(flex_pca_master_rounds), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        integer,                       intent(in)    :: pinds(:)
        type(string) :: fname
        integer :: ipart, first, last, nsel, numlen
        nsel = size(pinds)
        if( nsel < 1 ) THROW_HARD('flex_pca plan_partitions: empty particle selection')
        self%nparts_run = min(self%nparts_run, nsel)
        call self%qenv%new(params, self%nparts_run, numlen=params%numlen, nptcls=nsel)
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

    !> One qsys round. Round state travels in job_descr under registered keys (the refine3D
    !! idiom): stage, pcafit, which_iter (the global EM iteration), maxits (the budget) and nfits
    !! (1 single fit, 2 paired halves). extra_params=params has the scheduler clear the previous
    !! round's sentinels itself.
    subroutine master_run_stage( self, params, stage_id, label, which_iter, maxits, nfits )
        class(flex_pca_master_rounds), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        integer,                       intent(in)    :: stage_id
        character(len=*),              intent(in)    :: label
        integer, optional,             intent(in)    :: which_iter, maxits, nfits
        integer(timer_int_kind) :: t_round
        integer :: it_here, maxits_here, nfits_here
        if( .not. allocated(self%part_params) ) THROW_HARD('flex_pca run_stage before plan_partitions')
        t_round = tic()
        it_here = 0;     if( present(which_iter) ) it_here     = which_iter
        maxits_here = 0; if( present(maxits) )     maxits_here = maxits
        nfits_here = 1;  if( present(nfits) )      nfits_here  = nfits
        call self%job_descr%set('stage',      int2str(stage_id))
        call self%job_descr%set('pcafit',     int2str(self%fit_sel))
        call self%job_descr%set('which_iter', int2str(it_here))
        call self%job_descr%set('maxits',     int2str(maxits_here))
        call self%job_descr%set('nfits',      int2str(nfits_here))
        write(logfhandle,'(A,A,A,I0,A)') '>>> FLEX_PCA distributing ',label,' over ',self%nparts_run,' parts'
        call flush(logfhandle)
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, part_params=self%part_params, &
            &array=L_USE_SLURM_ARR, extra_params=params)
        write(logfhandle,'(A,A,A,F8.1)') '>>> FLEX_PCA ',label,' qsys round seconds=',toc(t_round)
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
            case(PCA_STAGE_PROBE, PCA_STAGE_POLISH, PCA_STAGE_EMBED, PCA_STAGE_STATES)
            case default
                THROW_HARD('invalid flex_pca worker stage')
        end select
        self%rounds%l_worker = .true.
        self%rounds%fit_sel  = params%pcafit
        ! adopt the master's node-local part directory (same derivation, same node, same cwd)
        call flex_pca_set_part_dir(flex_pca_local_part_dir(params, params%nparts))
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA (WORKER) part=',params%part, &
            &' stage=',params%stage,' fit=',self%rounds%fit_sel
        call flush(logfhandle)
    end subroutine worker_initialize

    subroutine worker_execute( self, params, build, cline )
        use simple_qsys_funs, only: qsys_job_finished
        class(flex_pca_worker_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        call run_flex_pca_worker(params, build, cline, self%rounds)
        call qsys_job_finished(params, string('simple_flex_pca_strategy :: worker_execute'))
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

    subroutine get_env_int_local( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        character(len=32) :: envval
        integer :: stat, ln, ival
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) val = ival
    end subroutine get_env_int_local

end module simple_flex_pca_strategy
