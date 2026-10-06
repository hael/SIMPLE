!@descr: the stream's 2D pool as a type: its project and history, dimensions, mask and low-pass limits, iteration job, and what it writes for the GUI, snapshots, 3D and the final project
!==============================================================================
! MODULE: simple_stream_pool2D
!
! PURPOSE:
!   The 2D pool of stream p06 (doc/policies/stream/pool2D_policy.md). p06 holds
!   an empty pool from its start (new), hands it the sets it imports
!   (append_sets) and starts it at the first import (start). Each pass then
!   takes back a completed iteration (update_status, update), sets up and
!   submits the next (update_aln_params, iterate), and writes snapshots,
!   publications for 3D and the final project (write_snapshot, publish,
!   finalise). What the GUI shows of the pool is one record (stats); the
!   pool's state is private.
!
!   The pool's files have fixed names (POOL_DIR, POOL_EXIT_CODE, the iteration
!   files), so a folder holds one pool.
!
! LIFECYCLE:
!   new() -> start() at the first import -> { iterate() } -> finalise() -> kill()
!==============================================================================
module simple_stream_pool2D
use simple_core_module_api
use simple_defs_environment,      only: SIMPLE_STREAM_POOL_PARTITION
use simple_cmdline,               only: cmdline
use simple_parameters,            only: parameters
use simple_sp_project,            only: sp_project
use simple_class_frcs,            only: class_frcs
use simple_image,                 only: image
use simple_stack_io,              only: stack_io
use simple_guistats,              only: guistats
use simple_starproject,           only: starproject
use simple_starproject_stream,    only: starproject_stream
use simple_qsys_env,              only: qsys_env
use simple_qsys_funs,             only: qsys_cleanup
use simple_qsys_job_record,       only: cancel_queued_job
use simple_qsys_async_job,        only: qsys_async_job, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
use simple_imgarr_utils,          only: rank_cavgs_stk
use simple_optics_maps,           only: import_latest_optics_map
use simple_stream_refine2D_utils, only: tidy_2Dstream_iter, build_pool_publication, pool_publication_names,&
                                       &snapshot_cavgs_meta, log_rss, draw_new_classes, append_project_sets
implicit none

public :: stream_pool2D, stream_pool2D_stats
private
#include "simple_local_flags.inc"

integer,          parameter :: POOL_NHISTORY  = 5               ! the completed iterations kept for snapshots
character(len=*), parameter :: POOL_JOB_LABEL = 'refine2D_pool' ! names its script, log and exit status as POOL_* do
logical,          parameter :: DEBUG_HERE     = .false.

!> What the GUI shows of the pool (stats)
type :: stream_pool2D_stats
    integer      :: last_complete_iter = 0      ! the latest iteration whose results are back
    integer      :: nassigned          = 0      ! particles assigned
    integer      :: nrejected          = 0      ! particles rejected
    real         :: resolution         = 999.   ! the current estimated resolution (A)
    real         :: mskdiam            = 0.     ! the mask the next iteration uses: diameter (A)
    real         :: msk_crop           = 0.     ! and radius at the working dimensions (px)
    type(string) :: cavgs_jpeg                  ! the class averages' sprite sheet, absolute ('' before the first)
    type(string) :: cavgs_mrc                   ! the class averages, absolute ('' before the first)
    integer      :: jpeg_ntilesx       = 0
    integer      :: jpeg_ntilesy       = 0
    integer, allocatable :: jpeg_map(:)         ! per tile: class index, population and resolution;
    integer, allocatable :: jpeg_pop(:)         ! not allocated before the first
    real,    allocatable :: jpeg_res(:)
end type stream_pool2D_stats

type :: stream_pool2D
    private
    ! allocatable (compile-time policy), allocated in new and released in kill
    type(cmdline),            allocatable :: cline_refine2D  ! the iterations' command line
    type(sp_project),         allocatable :: proj            ! the pool: every imported particle
    type(sp_project),         allocatable :: history(:)      ! the last POOL_NHISTORY completed iterations, a ring (history_slot)
    type(qsys_env),           allocatable :: qenv            ! the iterations' queue environment
    type(starproject_stream), allocatable :: starproj_stream ! the final project's STAR files
    type(qsys_async_job)                  :: job             ! the current iteration's job (its files: POOL_EXIT_CODE...)
    integer :: history_iter(POOL_NHISTORY) = 0               ! the iteration each history slot holds (0: none)
    ! dimensions, mask and resolution; none is a parameter (they are known at the first import)
    type(scaled_dims) :: dims                                ! the working dimensions, downscaled from the native ones
    integer :: native_box  = 0                               ! particle box of the imported sets (px)
    real    :: native_smpd = 0.                              ! their pixel size (A)
    real    :: mskdiam     = 0.                              ! the mask diameter (A)
    real    :: user_lpstop = 0.                              ! the user's low-pass stop (A), 0 when none was given
    real    :: lpstop_hard = 0.                              ! the hard low-pass limit (A): the user's, never below Nyquist
    real    :: lpstart     = 0.                              ! the low-pass ramp, from (A)
    real    :: lpstop      = 0.                              ! to (A)
    real    :: lpcen       = 0.                              ! the centering low-pass limit (A)
    real    :: resolution  = 999.                            ! the current estimated resolution (A)
    real    :: resolutions(POOL_NPREV_RES) = 999.            ! its history
    logical :: l_scaling   = .false.
    ! the iterations
    integer :: iter               = 0                        ! the latest iteration dispatched
    integer :: last_complete_iter = 0
    integer :: ncls               = 0
    integer :: numlen             = 0
    integer :: nattempts          = 0                        ! submissions of the current iteration: a failed first is retried once
    integer :: lim_ufrac_nptcls   = 0                        ! selected particles beyond which an iteration samples stacks
    logical, allocatable :: stacks_mask(:)                   ! the stacks of the current iteration
    type(string) :: center                                   ! the centering (yes|no) on a full update
    type(string) :: refs                                     ! the latest class averages, relative to the stage's folder
    type(string) :: orig_projfile                            ! the stage's project file
    logical :: l_active            = .false.                 ! started
    logical :: l_available         = .false.                 ! no iteration running
    logical :: l_failed            = .false.                 ! the current iteration failed twice: the pool stops
    logical :: l_update_frac_given = .false.                 ! update_frac on the stage's command line
    logical :: l_cenlp_given       = .false.                 ! cenlp on the stage's command line
    ! counts and convergence
    integer :: nptcls          = 0                           ! particles in the pool
    integer :: nptcls_rejected = 0
    integer :: ncls_rejected   = 0
    real    :: conv_frac       = 0.
    real    :: conv_mi_class   = 0.
    real    :: conv_score      = 0.
    ! the GUI's class averages
    type(string) :: jpeg                                     ! the current class-average JPEG ('' before the first)
    integer :: jpeg_ntiles  = 0
    integer :: jpeg_ntilesx = 0
    integer :: jpeg_ntilesy = 0
    integer :: stats_iter   = 0                              ! the iteration write_stats wrote last
    integer, allocatable :: jpeg_map(:), jpeg_pop(:)
    real,    allocatable :: jpeg_res(:)
contains
    ! what p06 calls
    procedure          :: new
    procedure          :: start
    procedure          :: kill
    procedure          :: append_sets
    procedure          :: iterate
    procedure          :: update
    procedure          :: update_aln_params
    procedure          :: update_status
    procedure          :: set_mskdiam
    procedure, nopass  :: cancel
    procedure          :: write_stats
    procedure          :: write_snapshot
    procedure          :: publish
    procedure          :: finalise
    procedure          :: iteration
    procedure          :: available
    procedure          :: failed
    procedure          :: stats
    ! the state part of start, public for simple_stream_pool2D_tester (a pool without a queue)
    procedure          :: init_state
    ! private steps
    procedure, private :: init_queue
    procedure, private :: submit_iteration
    procedure, private :: job_failed
    procedure, private :: set_dimensions
    procedure, private :: set_mask
    procedure, private :: set_resolution_limits
    procedure, private :: update_dims
    procedure, private :: update_for_gui
    procedure, private :: write_jpeg
    procedure, private :: write_project
    procedure, private :: rescale_cavgs
    procedure, private :: rank_cavgs
    procedure, private :: apply_snapshot_selection
end type stream_pool2D

contains

    !---------------- lifecycle ----------------

    !> An empty pool, as p06 holds it before its first import: its project takes the imported sets,
    !! and the queries answer as before any iteration (iteration 0, not available, not failed).
    subroutine new( self )
        class(stream_pool2D), intent(inout) :: self
        call self%kill()
        allocate(self%cline_refine2D, self%proj, self%qenv, self%starproj_stream)
        allocate(self%history(POOL_NHISTORY))
    end subroutine new

    !> The pool, at its first import: its state (init_state), then its queue environment.
    subroutine start( self, params, cline, spproj, box, smpd, mskdiam )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(inout) :: params
        class(cmdline),       intent(inout) :: cline
        class(sp_project),    intent(inout) :: spproj
        integer,              intent(in)    :: box
        real,                 intent(in)    :: smpd, mskdiam
        call self%init_state(params, cline, spproj, box, smpd, mskdiam)
        call self%init_queue(params)
    end subroutine start

    !> The pool's state at its first import, without a queue: @p box and @p smpd are the imported
    !! sets' (px, A) and @p mskdiam (A) the pool's mask diameter. From the stage's command line
    !! @p cline it takes what the iterations read later (whether update_frac and cenlp were given).
    subroutine init_state( self, params, cline, spproj, box, smpd, mskdiam )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(inout) :: params
        class(cmdline),       intent(inout) :: cline
        class(sp_project),    intent(inout) :: spproj
        integer,              intent(in)    :: box
        real,                 intent(in)    :: smpd, mskdiam
        type(string) :: carg, pool_sigma_path
        call seed_rnd
        ! general parameters
        self%l_update_frac_given = cline%defined('update_frac')
        self%l_cenlp_given       = cline%defined('cenlp')
        self%native_box  = box
        self%native_smpd = smpd
        self%mskdiam     = mskdiam
        self%user_lpstop = 0.
        if( cline%defined('lpstop') ) self%user_lpstop = params%lpstop
        call mskdiam2lplimits(self%mskdiam, self%lpstart, self%lpstop, self%lpcen)
        self%l_scaling     = .true.
        self%ncls          = params%ncls
        self%ncls_rejected = 0
        self%orig_projfile = params%projfile
        params%nparts_pool = params%nparts ! backwards compatibility
        ! bookkeeping & directory structure
        self%numlen      = len(int2str(params%nparts))
        self%refs        = ''
        self%l_available = .true.
        self%iter        = 0
        call simple_mkdir(POOL_DIR, verbose=.false.)
        call simple_mkdir(POOL_DIR//STDERROUT_DIR)
        call simple_mkdir(DIR_SNAPSHOT)
        self%proj%projinfo = spproj%projinfo
        self%proj%compenv  = spproj%compenv
        call self%proj%projinfo%delete_entry('projname')
        call self%proj%projinfo%delete_entry('projfile')
        call self%proj%projinfo%delete_entry('sigma2_state')
        pool_sigma_path = string(POOL_DIR)//'sigma2_state.bin'
        call self%proj%set_sigma2_state_path(pool_sigma_path)
        call pool_sigma_path%kill
        ! computing environment of the pool
        if( cline%defined('walltime') ) call self%proj%compenv%set(1,'walltime', params%walltime)
        ! commit to disk
        call self%proj%write(string(POOL_DIR)//POOL_PROJFILE)
        ! Pool command line
        call self%cline_refine2D%set('prg',        'refine2D_distr')
        call self%cline_refine2D%set('oritype',    'ptcl2D')
        call self%cline_refine2D%set('trs',        MINSHIFT)
        call self%cline_refine2D%set('projfile',   POOL_PROJFILE)
        call self%cline_refine2D%set('projname',   get_fbody(POOL_PROJFILE,'simple'))
        call self%cline_refine2D%set('sigma_est', params%sigma_est)
        if( cline%defined('cls_init') )then
            call self%cline_refine2D%set('cls_init', params%cls_init)
        else
            call self%cline_refine2D%set('cls_init', 'rand')
        endif
        self%center = 'yes'
        if( cline%defined('center') )then
            carg        = cline%get_carg('center')
            self%center = carg
            call carg%kill
        endif
        call self%cline_refine2D%set('center', self%center)
        if( cline%defined('center_type') )then
            call self%cline_refine2D%set('center_type', params%center_type)
        else
            call self%cline_refine2D%set('center_type', 'seg')
        endif
        call self%cline_refine2D%set('extr_iter', 99999)
        call self%cline_refine2D%set('extr_lim',   MAX_EXTRLIM2D)
        call self%cline_refine2D%set('mkdir',      'no')
        call self%cline_refine2D%set('mskdiam',    self%mskdiam)
        call self%cline_refine2D%set('async',      'yes') ! to enable hard termination
        call self%cline_refine2D%set('stream2d',   'yes') ! the only place this flag should be turned on
        call self%cline_refine2D%set('nparts',     params%nparts)
        if( cline%defined('worker_server') ) call self%cline_refine2D%set('worker_server', cline%get_carg('worker_server'))
        if( cline%defined('worker_server_nthr') ) call self%cline_refine2D%set('worker_server_nthr', cline%get_iarg('worker_server_nthr'))
        call self%cline_refine2D%delete('autoscale')
        ! when the 2D analysis is started from raw particles
        ! set # of ptcls beyond which fractional updates will be used
        self%lim_ufrac_nptcls = STREAM_NPTCLS_MAX
        if( cline%defined('nsample_max') ) self%lim_ufrac_nptcls = params%nsample_max
        ! the iterations' threads (the master's resources table, SIMPLE_STREAM_POOL_NTHR over it)
        call self%cline_refine2D%set('nthr', params%nthr)
        ! objective function
        call self%cline_refine2D%set('objfun', 'euclid')
        call self%cline_refine2D%set('ml_reg', params%ml_reg)
        call self%cline_refine2D%set('tau',     params%tau)
        ! refinement
        select case(trim(params%refine))
            case('snhc','snhc_smpl')
                call self%cline_refine2D%set( 'refine', params%refine)
            case DEFAULT
                THROW_HARD('UNSUPPORTED REFINE PARAMETER!')
        end select
        ! Determines dimensions for downscaling
        call self%set_dimensions()
        ! updates command-lines with resolution limits
        call self%set_resolution_limits(params)
        self%l_active = .true.
    end subroutine init_state

    ! The iterations' queue environment: the local queue, on SIMPLE_STREAM_POOL_PARTITION when set.
    subroutine init_queue( self, params )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(inout) :: params
        character(len=STDLEN) :: pool_part_env
        integer               :: envlen
        call get_environment_variable(SIMPLE_STREAM_POOL_PARTITION, pool_part_env, envlen)
        if(envlen > 0) then
            call self%qenv%new(params, params%nparts, stream=.true., exec_bin=string('simple_private_exec'),qsys_name=string('local'),&
            &qsys_partition=string(trim(pool_part_env)))
        else
            call self%qenv%new(params, params%nparts, stream=.true., exec_bin=string('simple_private_exec'),qsys_name=string('local'))
        end if
    end subroutine init_queue

    !> Releases everything and resets every field: a killed pool is an empty one.
    subroutine kill( self )
        class(stream_pool2D), intent(inout) :: self
        integer :: i
        if( allocated(self%cline_refine2D) )then
            call self%cline_refine2D%kill
            deallocate(self%cline_refine2D)
        endif
        if( allocated(self%proj) )then
            call self%proj%kill
            deallocate(self%proj)
        endif
        if( allocated(self%history) )then
            do i = 1,size(self%history)
                call self%history(i)%kill
            enddo
            deallocate(self%history)
        endif
        if( allocated(self%qenv) )then
            call self%qenv%kill
            deallocate(self%qenv)
        endif
        if( allocated(self%starproj_stream) ) deallocate(self%starproj_stream)
        call self%job%kill
        if( allocated(self%stacks_mask) ) deallocate(self%stacks_mask)
        if( allocated(self%jpeg_map)    ) deallocate(self%jpeg_map)
        if( allocated(self%jpeg_pop)    ) deallocate(self%jpeg_pop)
        if( allocated(self%jpeg_res)    ) deallocate(self%jpeg_res)
        call self%center%kill
        call self%refs%kill
        call self%orig_projfile%kill
        call self%jpeg%kill
        self%history_iter        = 0
        self%dims                = scaled_dims()
        self%native_box          = 0
        self%native_smpd         = 0.
        self%mskdiam             = 0.
        self%user_lpstop         = 0.
        self%lpstop_hard         = 0.
        self%lpstart             = 0.
        self%lpstop              = 0.
        self%lpcen               = 0.
        self%resolution          = 999.
        self%resolutions         = 999.
        self%l_scaling           = .false.
        self%iter                = 0
        self%last_complete_iter  = 0
        self%ncls                = 0
        self%numlen              = 0
        self%nattempts           = 0
        self%lim_ufrac_nptcls    = 0
        self%l_active            = .false.
        self%l_available         = .false.
        self%l_failed            = .false.
        self%l_update_frac_given = .false.
        self%l_cenlp_given       = .false.
        self%nptcls              = 0
        self%nptcls_rejected     = 0
        self%ncls_rejected       = 0
        self%conv_frac           = 0.
        self%conv_mi_class       = 0.
        self%conv_score          = 0.
        self%jpeg_ntiles         = 0
        self%jpeg_ntilesx        = 0
        self%jpeg_ntilesy        = 0
        self%stats_iter          = 0
    end subroutine kill

    !---------------- iterations ----------------

    !> Appends the sieve sets @p sets of one import to the pool's project, in order
    !! (append_project_sets: micrographs, stacks renumbered after the pool's, particles as new rows);
    !! @p nmics and @p nsel are the pool's micrographs and selected particles after it.
    subroutine append_sets( self, sets, nmics, nsel )
        class(stream_pool2D), intent(inout) :: self
        class(sp_project),    intent(inout) :: sets(:)
        integer,              intent(out)   :: nmics, nsel
        call append_project_sets(self%proj, sets)
        nmics = self%proj%os_mic%get_noris()
        nsel  = self%proj%os_ptcl2D%count_state_gt_zero()
    end subroutine append_sets

    ! Performs one iteration:
    ! updates to command-line, particles sampling, temporary project & execution
    subroutine iterate( self, params )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(inout) :: params
        logical, parameter            :: L_BENCH = .false.
        type(sp_project)              :: spproj
        integer(timer_int_kind)       :: t_tot
        integer,          allocatable :: nptcls_per_stk(:)
        type(string) :: stkname
        real         :: frac_update, smpd
        integer      :: iptcl,i, nptcls_tot, fromp, top, nstks_tot, jptcl, islot
        integer      :: nptcls_sel, istk, nptcls2update, nstks2update, jjptcl, ncls
        if( .not. self%l_active    ) return
        if( .not. self%l_available ) return
        if( L_BENCH ) t_tot  = tic()
        nptcls_tot           = self%proj%os_ptcl2D%get_noris()
        self%nptcls          = nptcls_tot
        self%nptcls_rejected = 0
        if( nptcls_tot == 0 ) return
        ! the completed iteration goes into the history, replacing the iteration POOL_NHISTORY
        ! before it, whose files tidy_2Dstream_iter removes below: history and files keep the
        ! same iterations
        if(self%iter .gt. 0) then
            islot = history_slot(self%iter)
            call self%history(islot)%kill
            call self%history(islot)%copy(self%proj)
            call simple_copy_file(string(POOL_DIR)//FRCS_FILE, string(POOL_DIR)//swap_suffix(FRCS_FILE,"_iter"//int2str_pad(self%iter, 3)//".bin",".bin"))
            call self%history(islot)%get_cavgs_stk(stkname, ncls, smpd)
            call self%history(islot)%os_out%kill
            call self%history(islot)%add_cavgs2os_out(stkname, smpd, 'cavg')
            call self%history(islot)%add_frcs2os_out(string(POOL_DIR)//swap_suffix(FRCS_FILE,"_iter"//int2str_pad(self%iter, 3)//".bin",".bin"),'frc2D')
            self%history_iter(islot) = self%iter
        end if
        self%iter = self%iter + 1 ! Global iteration counter update
        call self%cline_refine2D%set('ncls',     self%ncls)
        call self%cline_refine2D%set('startit', self%iter)
        call self%cline_refine2D%set('maxits',   self%iter)
        call self%cline_refine2D%set('frcs',     FRCS_FILE)
        call self%cline_refine2D%set('refs', self%refs)
        if( self%iter==1 )then
            if( self%cline_refine2D%defined('cls_init') )then
                ! references taken care of by refine2D_distr
                call self%cline_refine2D%delete('frcs')
                call self%cline_refine2D%delete('refs')
            endif
        else
            call self%cline_refine2D%delete('cls_init')
        endif
        ! Project metadata update
        spproj%projinfo = self%proj%projinfo
        spproj%compenv  = self%proj%compenv
        call spproj%projinfo%delete_entry('projname')
        call spproj%projinfo%delete_entry('projfile')
        call spproj%update_projinfo( self%cline_refine2D )
        ! Sampling of stacks that will be used for this iteration
        ! counting number of stacks & selected particles
        nstks_tot  = self%proj%os_stk%get_noris()
        allocate(nptcls_per_stk(nstks_tot), source=0)
        !$omp parallel do schedule(static) proc_bind(close) private(istk,fromp,top,iptcl) default(shared)
        do istk = 1,nstks_tot
            fromp = self%proj%os_stk%get_fromp(istk)
            top   = self%proj%os_stk%get_top(istk)
            do iptcl = fromp,top
                if( self%proj%os_ptcl2D%get_state(iptcl) > 0 )then
                    nptcls_per_stk(istk)  = nptcls_per_stk(istk) + 1 ! # ptcls with state=1
                endif
            enddo
        enddo
        !$omp end parallel do
        self%nptcls_rejected = self%nptcls - sum(nptcls_per_stk)
        ! Update info for gui
        call spproj%projinfo%set(1,'nptcls_tot',     self%nptcls)
        call spproj%projinfo%set(1,'nptcls_rejected',self%nptcls_rejected)
        ! Uniformly sample stacks
        call uniform_stack_sampling
        nstks2update = count(self%stacks_mask)
        ! Transfer stacks and particles
        call spproj%os_stk%new(nstks2update, is_ptcl=.false.)
        call spproj%os_ptcl2D%new(nptcls2update, is_ptcl=.true.)
        i     = 0
        jptcl = 0
        do istk = 1,nstks_tot
            fromp = self%proj%os_stk%get_fromp(istk)
            top   = self%proj%os_stk%get_top(istk)
            if( self%stacks_mask(istk) )then
                ! transfer alignement parameters for selected particles
                i = i + 1 ! stack index in spproj
                call spproj%os_stk%transfer_ori(i, self%proj%os_stk, istk)
                call spproj%os_stk%set(i, 'fromp', jptcl+1)
                !$omp parallel do private(iptcl,jjptcl) proc_bind(close) default(shared)
                do iptcl = fromp,top
                    jjptcl = jptcl+iptcl-fromp+1
                    call spproj%os_ptcl2D%transfer_ori(jjptcl, self%proj%os_ptcl2D, iptcl)
                    call spproj%os_ptcl2D%set_stkind(jjptcl, i)
                enddo
                !$omp end parallel do
                jptcl = jptcl + (top-fromp+1)
                call spproj%os_stk%set(i, 'top', jptcl)
            endif
        enddo
        call spproj%os_ptcl3D%new(nptcls2update, is_ptcl=.true.)
        spproj%os_cls2D = self%proj%os_cls2D
        ! the new particles get a populated class, drawn reproducibly
        if( self%iter >= 2 ) call draw_new_classes(spproj, nptcls2update, self%iter, self%ncls)
        ! the sampled stacks are the update set (decision 5): every particle of the sample is
        ! updated, the others keep their parameters; only the user's update_frac thins the sample,
        ! and then the class averages are not centered
        call self%cline_refine2D%delete('update_frac')
        frac_update = 1.0
        if( self%l_update_frac_given ) frac_update = params%update_frac
        if( frac_update < 0.99999 )then
            call self%cline_refine2D%set('update_frac', frac_update)
            call self%cline_refine2D%set('center',      'no')
        else
            call self%cline_refine2D%set('center',      self%center)
        endif
        ! write project, and keep it as made for a retry
        call spproj%write(string(POOL_DIR)//POOL_PROJFILE)
        call spproj%kill
        call simple_copy_file(string(POOL_DIR)//POOL_PROJFILE, string(POOL_DIR)//POOL_INPUT_PROJFILE)
        ! pool stats
        call self%write_stats(params)
        ! execution
        self%nattempts = 0
        call self%submit_iteration(params)
        write(logfhandle,'(A,I6,A,I8,A3,I8,A)')'>>> POOL         INITIATED ITERATION ',self%iter,' WITH ',nptcls_sel,&
        &' / ', sum(nptcls_per_stk),' PARTICLES'
        if( L_BENCH ) print *,'timer analyze2D_pool tot : ',toc(t_tot)
        ! cleanup
        if( allocated(nptcls_per_stk) )      deallocate(nptcls_per_stk)
        ! the files of the iteration that has just left the history
        call tidy_2Dstream_iter(self%iter - 1 - POOL_NHISTORY)

      contains

        subroutine uniform_stack_sampling
            use simple_ran_tabu
            type(ran_tabu) :: random_generator
            integer        :: stk_order(nstks_tot)
            integer        :: i, j
            if( allocated(self%stacks_mask) ) deallocate(self%stacks_mask)
            allocate(self%stacks_mask(nstks_tot), source=.false.)
            stk_order        = (/(i,i=1,nstks_tot)/)
            random_generator = ran_tabu(nstks_tot)
            call random_generator%shuffle(stk_order)
            nptcls2update = 0 ! # of ptcls including state=0 within selected stacks
            nptcls_sel    = 0 ! # of ptcls excluding state=1 within selected stacks
            do i = 1,nstks_tot
                j = stk_order(i)
                if( nptcls_sel > self%lim_ufrac_nptcls ) cycle
                nptcls_sel    = nptcls_sel    + nptcls_per_stk(j)
                nptcls2update = nptcls2update + self%proj%os_stk%get_int(j, 'nptcls')
                self%stacks_mask(j) = .true.
            enddo
            call random_generator%kill
        end subroutine uniform_stack_sampling

    end subroutine iterate

    ! Submits the pool's current iteration as the pool's job (qsys_async_job, in the stage's folder:
    ! POOL_DIR is empty), which removes the status and job record of an earlier attempt; the attempt
    ! is counted.
    subroutine submit_iteration( self, params )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(in)    :: params
        type(cmdline), allocatable :: pool_clines(:)
        type(string)               :: cwd
        call simple_getcwd(cwd)
        self%nattempts = self%nattempts + 1
        if( params%cc_objfun == OBJFUN_EUCLID )then
            ! Stream pool membership changes invalidate row identity. Rebuild a
            ! complete canonical bootstrap for the exact pool layout, then run
            ! clustering in the same queued script so no consumer can observe
            ! a missing or stale state.
            allocate(pool_clines(2))
            pool_clines(1) = self%cline_refine2D
            call pool_clines(1)%set('prg', 'calc_pspec')
            call pool_clines(1)%set('mkdir', 'no')
            call pool_clines(1)%delete('stream2d')
            call pool_clines(1)%delete('update_frac')
            pool_clines(2) = self%cline_refine2D
            call self%job%start_seq(self%qenv, pool_clines, cwd, POOL_JOB_LABEL)
            call pool_clines(:)%kill
            deallocate(pool_clines)
        else
            call self%job%start(self%qenv, self%cline_refine2D, cwd, POOL_JOB_LABEL)
        endif
        self%l_available = .false.
    end subroutine submit_iteration

    ! .true. when the current iteration's job has ended without refine2D finishing: it exited,
    ! whatever its status, or vanished without one, which the job's liveness check finds
    ! (qsys_async_job: a walltime kill, a lost node).
    logical function job_failed( self )
        class(stream_pool2D), intent(inout) :: self
        integer :: job_status
        job_failed = .false.
        job_status = self%job%status()
        if( job_status /= ASYNC_JOB_DONE .and. job_status /= ASYNC_JOB_FAILED ) return
        ! refine2D touches its marker before it exits and the script writes the status after
        if( file_exists(POOL_DIR//REFINE2D_FINISHED) ) return
        job_failed = .true.
    end function job_failed

    ! Cancels the job of the pool's current iteration when it recorded itself and has not exited
    ! (simple_qsys_job_record): on stop, and on a restart for a job a crashed stage left running.
    ! It reads only the pool's files, so it works on a pool not started.
    subroutine cancel
        if( cancel_queued_job(string(POOL_DIR//POOL_EXIT_CODE)) )then
            write(logfhandle,'(A)') '>>> CANCELLED THE RUNNING POOL ITERATION'
        endif
    end subroutine cancel

    ! The pool's working dimensions from its native ones: downscaled to a pixel size of up to
    ! MAX_SMPD and never below a CHUNK_MINBOXSZ box (setup_downscaling's rule). The pool's command
    ! line carries the cropped dimensions only; the native ones come from its project.
    subroutine set_dimensions( self )
        class(stream_pool2D), intent(inout) :: self
        real    :: smpd, scale_factor
        integer :: box
        if( self%native_box == 0 ) THROW_HARD('the pool has no native box; set_dimensions')
        self%dims%smpd = self%native_smpd
        self%dims%box  = self%native_box
        if( self%l_scaling .and. self%native_box >= CHUNK_MINBOXSZ )then
            call autoscale(self%native_box, self%native_smpd, MAX_SMPD, box, smpd, scale_factor, minbox=CHUNK_MINBOXSZ)
            self%l_scaling = box < self%native_box
            if( self%l_scaling )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',self%native_box,'/',box
                self%dims%smpd = smpd
                self%dims%box  = box
            endif
        endif
        self%dims%boxpd = 2*round2even(KBALPHA*real(self%dims%box/2)) ! logics from parameters
        call self%set_mask()
        ! Scaling-related command lines update
        call self%cline_refine2D%set('smpd_crop',   self%dims%smpd)
        call self%cline_refine2D%set('box_crop',    self%dims%box)
    end subroutine set_dimensions

    ! The pool's mask radius at its working dimensions from its mask diameter, clamped (and logged)
    ! to (box - COSMSKHALFWIDTH)/2 pixels, so a diameter beyond the box never reaches the workers
    ! (D40); the diameter and the radius go on the pool's command line.
    subroutine set_mask( self )
        class(stream_pool2D), intent(inout) :: self
        real :: msk_max
        msk_max       = (real(self%dims%box) - COSMSKHALFWIDTH) / 2.
        self%dims%msk = round2even(self%mskdiam / self%dims%smpd / 2.)
        if( real(self%dims%msk) > msk_max )then
            write(logfhandle,'(A,F8.2,A,F8.2,A)') '>>> MASK DIAMETER ', self%mskdiam, ' A EXCEEDS THE POOL''S BOX; CLAMPED TO ',&
                &2. * msk_max * self%dims%smpd, ' A'
            self%dims%msk = floor(msk_max)
            self%mskdiam  = 2. * msk_max * self%dims%smpd
        endif
        call self%cline_refine2D%set('mskdiam',  self%mskdiam)
        call self%cline_refine2D%set('msk_crop', self%dims%msk)
    end subroutine set_mask

    ! The resolution limits on the pool's command line
    subroutine set_resolution_limits( self, params )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(inout) :: params
        self%lpstart     = max(self%lpstart, 2.0*self%dims%smpd)
        self%lpstop_hard = max(2.0*self%dims%smpd, self%user_lpstop)
        call self%cline_refine2D%set('lpstart',   self%lpstart)
        call self%cline_refine2D%set('lpstop',    self%lpstop_hard)
        if( .not. self%l_cenlp_given )then
            call self%cline_refine2D%set( 'cenlp', self%lpcen)
        else
            call self%cline_refine2D%set( 'cenlp', params%cenlp)
        endif
        write(logfhandle,'(A,F5.1)') '>>> STARTING LOW-PASS LIMIT  (IN A): ', self%lpstart
        write(logfhandle,'(A,F5.1)') '>>> HARD RESOLUTION LIMIT    (IN A): ', self%lpstop_hard
        write(logfhandle,'(A,F5.1)') '>>> CENTERING LOW-PASS LIMIT (IN A): ', self%lpcen
    end subroutine set_resolution_limits

    ! A new mask diameter (A) for the next pool iterations. The pool command line also carries
    ! the cropped mask radius (pixels), which the workers take over the one parameters would
    ! derive from mskdiam, so it is updated with it. Before the pool starts, the stage keeps the
    ! diameter and gives it to start.
    subroutine set_mskdiam( self, new_mskdiam )
        class(stream_pool2D), intent(inout) :: self
        integer,              intent(in)    :: new_mskdiam
        write(*,'(A,I4,A)')'>>> UPDATED MASK DIAMETER TO', new_mskdiam ,'Å'
        self%mskdiam = real(new_mskdiam)
        call self%cline_refine2D%set('mskdiam',   self%mskdiam)
        if( self%dims%smpd > 0. )then
            call self%set_mask()
            ! the low-pass ramp and the centering limit follow the mask
            call mskdiam2lplimits(self%mskdiam, self%lpstart, self%lpstop, self%lpcen)
            self%lpstart = max(self%lpstart, 2.0*self%dims%smpd)
            call self%cline_refine2D%set('lpstart', self%lpstart)
            if( .not. self%l_cenlp_given ) call self%cline_refine2D%set('cenlp', self%lpcen)
            write(logfhandle,'(A,F5.1,A,F5.1,A)') '>>> LOW-PASS RAMP FROM ', self%lpstart, ' A, CENTERING LOW-PASS ', self%lpcen, ' A'
        endif
    end subroutine set_mskdiam

    ! Reports alignment info from completed iteration of subset
    ! of particles back to the pool
    subroutine update( self, params )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(inout) :: params
        integer,      allocatable :: pops(:)
        type(sp_project) :: spproj
        type(oris)       :: os
        type(class_frcs) :: frcs
        type(string)     :: fname, cwd
        integer          :: i, it, jptcl, iptcl, istk
        if( .not. self%l_active    ) return
        if( .not. self%l_available ) return
        call del_file(POOL_DIR//REFINE2D_FINISHED)
        ! iteration info
        fname = POOL_DIR//STATS_FILE
        if( file_exists(fname) )then
            call os%new(1,is_ptcl=.false.)
            call os%read(fname)
            it = os%get_int(1,'ITERATION')
            if( it == self%iter )then
                self%conv_mi_class = os%get(1,'CLASS_OVERLAP')
                self%conv_frac     = os%get(1,'SEARCH_SPACE_SCANNED')
                self%conv_score    = os%get(1,'SCORE')
                ! new
                self%last_complete_iter = it
                call simple_getcwd(cwd)
                self%jpeg = cwd//'/'//CAVGS_ITER_FBODY//int2str_pad(it, 3)//'.jpg'
                ! its grid follows the classes, transferred below
                ! end new
                write(logfhandle,'(A,I6,A,F7.3,A,F7.3,A,F7.3)')'>>> POOL         ITERATION ',it,&
                    &'; CLASS OVERLAP: ',self%conv_mi_class,'; SEARCH SPACE SCANNED: ',self%conv_frac,'; SCORE: ',self%conv_score
            endif
            call os%kill
        endif
        ! transfer to pool
        call spproj%read_segment('cls2D', string(POOL_DIR)//POOL_PROJFILE)
        if( spproj%os_cls2D%get_noris() == 0 )then
            ! not executed yet, do nothing
        else
            if( .not.allocated(self%stacks_mask) )then
                THROW_HARD('Critical ERROR 0') ! first time
            endif
            ! transfer particles parameters
            call spproj%read_segment('stk',   string(POOL_DIR)//POOL_PROJFILE)
            call spproj%read_segment('ptcl2D',string(POOL_DIR)//POOL_PROJFILE)
            i = 0
            do istk = 1,size(self%stacks_mask)
                if( self%stacks_mask(istk) )then
                    i = i+1
                    iptcl = self%proj%os_stk%get_fromp(istk)
                    do jptcl = spproj%os_stk%get_fromp(i),spproj%os_stk%get_top(i)
                        if( spproj%os_ptcl2D%get_state(jptcl) > 0 )then
                            call self%proj%os_ptcl2D%transfer_2Dparams(iptcl, spproj%os_ptcl2D, jptcl)
                        endif
                        iptcl = iptcl+1
                    enddo
                endif
            enddo
            ! update classes info
            call self%proj%os_ptcl2D%get_pops(pops, 'class', maxn=self%ncls)
            self%proj%os_cls2D = spproj%os_cls2D
            call self%proj%os_cls2D%set_all('pop', real(pops))
            ! the sprite sheet's grid, as refine2D lays it out (mrc2jpeg_tiled: a tile per class,
            ! floor(sqrt(n)) across). Taken from the iteration's own classes: the pool has none
            ! before the first iteration's, which gave a 1 x 0 grid and blank tiles in the GUI.
            self%jpeg_ntiles  = self%proj%os_cls2D%get_noris()
            self%jpeg_ntilesx = max(1, floor(sqrt(real(self%jpeg_ntiles))))
            self%jpeg_ntilesy = ceiling(real(self%jpeg_ntiles)/real(self%jpeg_ntilesx))
            ! update thumbnail metadata
            if(allocated(self%jpeg_map)) deallocate(self%jpeg_map)
            if(allocated(self%jpeg_pop)) deallocate(self%jpeg_pop)
            if(allocated(self%jpeg_res)) deallocate(self%jpeg_res)
            self%jpeg_pop = self%proj%os_cls2D%get_all_asint('pop')
            self%jpeg_res = self%proj%os_cls2D%get_all('res')
            allocate(self%jpeg_map, mold=self%jpeg_pop)
            do i=1, size(self%jpeg_map)
                self%jpeg_map(i) = i
            end do
            ! estimate resolution
            call frcs%read(string(POOL_DIR)//FRCS_FILE)
            self%resolution = frcs%estimate_lp_for_align()
            write(logfhandle,'(A,F5.1)')'>>> CURRENT POOL RESOLUTION: ',self%resolution
            call frcs%kill
            ! deal with dimensions/resolution update
            call self%update_dims(params)
            ! for gui
            call self%update_for_gui(params)
        endif
        call spproj%kill
    end subroutine update

    ! This controls the evolution of the pool alignement parameters:
    ! lp, Gaussian filter, trs, extr_iter
    subroutine update_aln_params( self )
        class(stream_pool2D), intent(inout) :: self
        integer, parameter :: ITERLIM    = 20
        integer, parameter :: ITERSHIFT  = 5
        real :: lp, gamma
        if( .not. self%l_active    ) return
        if( .not. self%l_available ) return
        if( self%iter < ITERLIM )then
            gamma = min(1., max(0., real(ITERLIM-self%iter)/real(ITERLIM)))
            ! offset
            if( self%iter < ITERSHIFT )then
                call self%cline_refine2D%set('trs', 0.)
            else
                call self%cline_refine2D%set('trs', MINSHIFT)
            endif
            ! resolution limit
            lp = self%lpstop + (self%lpstart-self%lpstop) * gamma
            call self%cline_refine2D%set('lp', lp)
            ! Extremal iteration
            call self%cline_refine2D%set('extr_iter', self%iter+1)
            call self%cline_refine2D%set('extr_lim',   ITERLIM)
            ! Gaussian filter
            call self%cline_refine2D%set('gauref',   'yes')
            call self%cline_refine2D%set('gaufreq', lp)
        else
            call self%cline_refine2D%set('trs', MINSHIFT)
            call self%cline_refine2D%set('gauref', 'no')
            call self%cline_refine2D%delete('extr_iter')
            call self%cline_refine2D%delete('gaufreq')
            call self%cline_refine2D%delete('lp')
        endif
    end subroutine update_aln_params

    ! Deals with pool dimensions & resolution update
    subroutine update_dims( self, params )
        use simple_procimgstk,    only: scale_imgfile
        use simple_classaverager, only: cavger_pad_carried_sums
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(inout) :: params
        type(scaled_dims) :: new_dims, prev_dims
        type(oris)        :: os
        type(class_frcs)  :: frcs
        type(string)      :: str, str_tmp_mrc
        real              :: scale_factor
        integer           :: ldim(3)
        ! resolution book-keeping
        self%resolutions(1:POOL_NPREV_RES-1) = self%resolutions(2:POOL_NPREV_RES)
        self%resolutions(POOL_NPREV_RES)     = self%resolution
        ! optional
        if( trim(params%dynreslim).ne.'yes' ) return
        prev_dims = self%dims
        ! Auto-scaling?
        if( trim(params%autoscale) .ne. 'yes' ) return
        ! Hard limit reached?
        if( self%dims%smpd < POOL_SMPD_HARD_LIMIT ) return
        ! Too early?
        if( self%iter < 10 ) return
        ! Current resolution at Nyquist?
        if( abs(self%resolution-2.*self%dims%smpd) > 0.01 ) return
        ! When POOL_NPREV_RES iterations are at Nyquist the pool resolution may be updated
        if( any(abs(self%resolutions-self%resolution) > 0.01 ) ) return
        ! determines new dimensions
        new_dims%box   = find_larger_magic_box(self%dims%box+1)
        scale_factor   = real(new_dims%box) / real(self%native_box)
        if( scale_factor > 0.99 ) return ! safety
        new_dims%smpd  = self%native_smpd / scale_factor
        new_dims%boxpd = 2 * round2even(KBALPHA * real(new_dims%box/2)) ! logics from parameters
        ! New dimensions are accepted when new Nyquist is > 5/4 of original
        if( new_dims%smpd < 1.25*self%native_smpd ) return
        ! Update global variables
        self%l_scaling   = .true.
        self%dims        = new_dims
        call self%set_mask()
        self%lpstop_hard = max(2.0*self%dims%smpd, self%user_lpstop)
        call self%cline_refine2D%set('lpstop',     self%lpstop_hard)
        call self%cline_refine2D%set('smpd_crop', self%dims%smpd)
        call self%cline_refine2D%set('box_crop',   self%dims%box)
        write(logfhandle,'(A)')             '>>> UPDATING POOL DIMENSIONS '
        write(logfhandle,'(A,I5,A1,I5)')    '>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',self%native_box,'/',self%dims%box
        write(logfhandle,'(A,F5.2,A1,F5.2)')'>>> ORIGINAL/CROPPED PIXEL SIZE (Angs)  : ',self%native_smpd,'/',self%dims%smpd
        write(logfhandle,'(A,F5.1)')        '>>> POOL   HARD RESOLUTION LIMIT (Angs) : ',self%lpstop_hard
        ! upsample cavgs
        ldim = [self%dims%box,self%dims%box,1]
        str_tmp_mrc = TMP_STK_FNAME
        call scale_imgfile(self%refs, str_tmp_mrc, prev_dims%smpd, ldim, self%dims%smpd)
        call simple_rename(str_tmp_mrc,self%refs)
        str  = add2fbody(self%refs, MRC_EXT,'_even')
        call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, self%dims%smpd)
        call simple_rename(str_tmp_mrc,str)
        str  = add2fbody(self%refs, MRC_EXT,'_odd')
        call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, self%dims%smpd)
        call simple_rename(str_tmp_mrc,str)
        ! upsample the carried class sums, one set
        call cavger_pad_carried_sums(self%dims%box, self%dims%smpd)
        ! update cls2D field
        os = self%proj%os_cls2D
        call self%proj%os_out%kill
        call self%proj%add_cavgs2os_out(self%refs, self%dims%smpd, 'cavg', clspath=.true.)
        self%proj%os_cls2D = os
        call os%kill
        ! rescale frcs
        call frcs%read(string(FRCS_FILE))
        call frcs%pad(self%dims%smpd, self%dims%box)
        call frcs%write(string(FRCS_FILE))
        call frcs%kill
        call self%proj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
    end subroutine update_dims

    !> Points the pool project at the iteration's class averages and writes its class STAR file
    subroutine update_for_gui( self, params )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(in)    :: params
        type(oris)        :: os_backup
        type(starproject) :: starproj
        os_backup = self%proj%os_cls2D
        call self%proj%add_cavgs2os_out(string(POOL_DIR)//self%refs, self%dims%smpd, 'cavg')
        self%proj%os_cls2D = os_backup
        ! Write star file for iteration
        call starproj%export_cls2D(self%proj, self%iter)
        call self%proj%os_cls2D%delete_entry('stk')
        call os_backup%kill
        call starproj%kill
    end subroutine update_for_gui

    ! Flags pool availibility & updates the global name of references. An iteration whose job
    ! exited without finishing is submitted again once, from its project as made, after its log is
    ! kept aside and the files of its parts are removed; a second failure stops the pool (failed).
    subroutine update_status( self, params )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(in)    :: params
        type(string) :: failed_log
        if( .not. self%l_active ) return
        if( self%l_available .or. self%l_failed ) return
        self%l_available = file_exists(POOL_DIR//REFINE2D_FINISHED)
        if( self%l_available )then
            if( self%iter >= 1 ) self%refs = CAVGS_ITER_FBODY//int2str_pad(self%iter,3)//MRC_EXT
        else if( self%job_failed() )then
            failed_log = POOL_DIR//POOL_LOGFILE//'_failed_iter'//int2str_pad(self%iter,3)//'_attempt'//int2str(self%nattempts)
            if( file_exists(POOL_DIR//POOL_LOGFILE) ) call simple_rename(string(POOL_DIR//POOL_LOGFILE), failed_log)
            if( self%nattempts < 2 .and. file_exists(POOL_DIR//POOL_INPUT_PROJFILE) )then
                write(logfhandle,'(A,I6,A,A)') '>>> WARNING: POOL ITERATION ', self%iter, ' FAILED; RETRYING IT ONCE. LOG: ',&
                    &failed_log%to_char()
                call qsys_cleanup(params)
                call simple_copy_file(string(POOL_DIR)//POOL_INPUT_PROJFILE, string(POOL_DIR)//POOL_PROJFILE)
                call self%submit_iteration(params)
            else
                write(logfhandle,'(A,I6,A,A)') '>>> POOL ITERATION ', self%iter, ' FAILED AGAIN; THE POOL STOPS. LOG: ',&
                    &failed_log%to_char()
                self%l_failed = .true.
            endif
        endif
    end subroutine update_status

    !---------------- outputs ----------------

    ! write jpeg of the latest class averages
    subroutine write_jpeg( self, params, filename)
        class(stream_pool2D),    intent(inout) :: self
        class(parameters),       intent(in)    :: params
        class(string), optional, intent(in)    :: filename
        type(string)   :: jpeg_path, cwd
        type(image)    :: img, img_pad, img_jpeg
        type(stack_io) :: stkio_r
        integer        :: ldim_stk(3)
        integer        :: ncls_here, xtiles, ytiles, icls, ix, iy, ntiles
        call simple_getcwd(cwd)
        if(present(filename)) then
            jpeg_path = filename
        else
            jpeg_path = fname_new_ext(self%refs, "jpeg") ! temporarily jpeg so compatible with old pool_stats.
        end if
        if(.not. file_exists(self%refs)) return
        if(file_exists(jpeg_path))       return
        if(allocated(self%jpeg_map)) deallocate(self%jpeg_map)
        if(allocated(self%jpeg_pop)) deallocate(self%jpeg_pop)
        if(allocated(self%jpeg_res)) deallocate(self%jpeg_res)
        allocate(self%jpeg_map(0))
        allocate(self%jpeg_pop(0))
        allocate(self%jpeg_res(0))
        call find_ldim_nptcls(self%refs, ldim_stk, ncls_here)
        if(ncls_here .ne. self%proj%os_cls2D%get_noris()) THROW_HARD('ncls and n_noris mismatch')
        xtiles = floor(sqrt(real(self%ncls)))
        ytiles = ceiling(real(self%ncls) / real(xtiles))
        call img%new([ldim_stk(1), ldim_stk(1), 1], self%dims%smpd)
        call img_pad%new([JPEG_DIM, JPEG_DIM, 1], self%dims%smpd)
        call img_jpeg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], self%dims%smpd)
        call stkio_r%open(self%refs, self%dims%smpd, 'read', bufsz=ncls_here)
        call stkio_r%read_whole
        ix = 1
        iy = 1
        ntiles = 0
        ! mask memoization
        call img%memoize_mask_coords
        do icls=1, ncls_here
            if(self%proj%os_cls2D%get(icls,'state') < 0.5) cycle
            if(self%proj%os_cls2D%get(icls,'pop')   < 0.5) cycle
            self%jpeg_map = [self%jpeg_map, icls]
            self%jpeg_pop = [self%jpeg_pop, nint(self%proj%os_cls2D%get(icls,'pop'))]
            self%jpeg_res = [self%jpeg_res, self%proj%os_cls2D%get(icls,'res')]
            call img%zero_and_unflag_ft
            call stkio_r%get_image(icls, img)
            call img%mask2D_softavg(self%mskdiam / (2 * self%dims%smpd))
            call img%fft
            if(ldim_stk(1) > JPEG_DIM) then
                call img%clip(img_pad)
            else
                call img%pad(img_pad, backgr=0., antialiasing=.false.)
            end if
            call img_pad%ifft
            call img_jpeg%tile(img_pad, ix, iy)
            ntiles = ntiles + 1
            ix = ix + 1
            if(ix > xtiles) then
                ix = 1
                iy = iy + 1
            end if
        enddo
        call stkio_r%close()
        call img_jpeg%write_jpg(jpeg_path)
        self%jpeg = cwd%to_char() // '/' // jpeg_path%to_char()
        self%jpeg_ntiles  = ntiles
        self%jpeg_ntilesx = xtiles
        self%jpeg_ntilesy = ytiles
        call img%kill()
        call img_pad%kill()
        call img_jpeg%kill()
    end subroutine write_jpeg

    ! When refine2D has not written the sprite sheet of the latest completed iteration, writes the
    ! iteration's class averages as images (generate_2D_jpeg, write_jpeg)
    subroutine write_stats( self, params )
        class(stream_pool2D), intent(inout) :: self
        class(parameters),    intent(in)    :: params
        type(guistats) :: pool_stats
        type(string)   :: cwd
        call simple_getcwd(cwd)
        if(file_exists(cwd//'/'//CLS2D_STARFBODY//'_iter'//int2str_pad(self%iter,3)//STAR_EXT)) then
            self%stats_iter = self%iter
        else if(file_exists(cwd//'/'//CLS2D_STARFBODY//'_iter'//int2str_pad(self%iter - 1,3)//STAR_EXT)) then
            self%stats_iter = self%iter - 1
        endif
        if(.not. file_exists(cwd//'/'//CAVGS_ITER_FBODY//int2str_pad(self%stats_iter, 3)//'.jpg')) then
            call pool_stats%init
            call pool_stats%generate_2D_jpeg('latest', '', self%proj%os_cls2D, self%stats_iter, self%dims%smpd)
            call pool_stats%kill
            self%last_complete_iter = self%stats_iter
            call self%write_jpeg(params)
        endif
    end subroutine write_stats

    ! Deselects the particles of the classes @p selection does not name (all classes are kept
    ! when it names none in range).
    subroutine apply_snapshot_selection( self, snapshot_projfile, selection )
        class(stream_pool2D), intent(in)    :: self
        type(sp_project),     intent(inout) :: snapshot_projfile
        integer,              intent(in)    :: selection(:)
        logical, allocatable :: cls_mask(:)
        integer              :: nptcls_rejected, ncls_rejected, iptcl
        integer              :: icls, jcls, i
        if( snapshot_projfile%os_cls2D%get_noris() == 0 ) return
        allocate(cls_mask(self%ncls), source=.false.)
        do i = 1, size(selection)
            icls = selection(i)
            if( icls < 1 ) cycle
            if( icls > self%ncls ) cycle
            cls_mask(icls) = .true.
        enddo
        if( count(cls_mask) == 0 ) return
        ncls_rejected   = 0
        do icls = 1,self%ncls
            if(cls_mask(icls) ) cycle
            nptcls_rejected = 0
            !$omp parallel do private(iptcl,jcls) reduction(+:nptcls_rejected) proc_bind(close)
            do iptcl = 1,snapshot_projfile%os_ptcl2D%get_noris()
                if( snapshot_projfile%os_ptcl2D%get_state(iptcl) == 0 )cycle
                jcls = snapshot_projfile%os_ptcl2D%get_class(iptcl)
                if( jcls == icls )then
                    call snapshot_projfile%os_ptcl2D%reject(iptcl)
                    call snapshot_projfile%os_ptcl2D%delete_2Dclustering(iptcl)
                    nptcls_rejected = nptcls_rejected + 1
                endif
            enddo
            !$omp end parallel do
            if( nptcls_rejected > 0 )then
                ncls_rejected = ncls_rejected + 1
                call snapshot_projfile%os_cls2D%set_state(icls,0)
                call snapshot_projfile%os_cls2D%set(icls,'pop',0.)
                call snapshot_projfile%os_cls2D%set(icls,'corr',-1.)
                call snapshot_projfile%os_cls2D%set(icls,'prev_pop_even',0.)
                call snapshot_projfile%os_cls2D%set(icls,'prev_pop_odd', 0.)
                write(logfhandle,'(A,I6,A,I4)')'>>> USER REJECTED FROM SNAPSHOT: ',nptcls_rejected,' PARTICLE(S) IN CLASS ',icls
            endif
        enddo
    end subroutine apply_snapshot_selection

    ! The final class averages ranked by resolution, beside them (<refs>_ranked), over the classes
    ! of the stage's project; in this process (rank_cavgs_stk), as commander_rank_cavgs does
    subroutine rank_cavgs( self )
        class(stream_pool2D), intent(in) :: self
        type(sp_project) :: spproj
        type(string)     :: refs_ranked, stk
        refs_ranked = add2fbody(self%refs, MRC_EXT ,'_ranked')
        stk = string(POOL_DIR)//self%refs
        if( .not. file_exists(stk) ) return
        call spproj%read_segment('cls2D', self%orig_projfile)
        call rank_cavgs_stk(spproj%os_cls2D, stk, find_img_smpd(stk), 'res', refs_ranked, string('classdoc_ranked.txt'))
        call spproj%kill
    end subroutine rank_cavgs

    ! Final rescaling of references
    subroutine rescale_cavgs( self, src, dest )
        class(stream_pool2D), intent(in) :: self
        class(string),        intent(in) :: src, dest
        integer, allocatable :: cls_pop(:)
        type(image)    :: img, img_pad
        type(stack_io) :: stkio_r, stkio_w
        type(string)   :: dest_here
        integer        :: ldim(3),icls, ncls_here
        if(src == dest)then
            dest_here = 'tmp_cavgs.mrc'
        else
            dest_here = dest
        endif
        call img%new([self%dims%box,self%dims%box,1],self%dims%smpd)
        call img_pad%new([self%native_box,self%native_box,1],self%native_smpd)
        cls_pop = nint(self%proj%os_cls2D%get_all('pop'))
        call find_ldim_nptcls(src,ldim,ncls_here)
        call stkio_r%open(src, self%dims%smpd, 'read', bufsz=ncls_here)
        call stkio_r%read_whole
        call stkio_w%open(dest_here, self%native_smpd, 'write', box=self%native_box, bufsz=ncls_here)
        do icls = 1,ncls_here
            if( cls_pop(icls) > 0 )then
                call img%zero_and_unflag_ft
                call stkio_r%get_image(icls, img)
                call img%fft
                call img%pad(img_pad, backgr=0., antialiasing=.false.)
                call img_pad%ifft
            else
                img_pad = 0.
            endif
            call stkio_w%write(icls, img_pad)
        enddo
        call stkio_r%close
        call stkio_w%close
        if ( src == dest ) call simple_rename('tmp_cavgs.mrc',dest)
        call img%kill
        call img_pad%kill
    end subroutine rescale_cavgs

    !> Ends the pool and writes its final project (pool2D_policy.md section 10): from the last
    !! complete iteration, or, before the first, the raw project of every imported particle, both
    !! with the newest optics map's groups (with @p optics_dir); then cleans up.
    subroutine finalise( self, params, optics_dir )
        class(stream_pool2D),    intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(string), optional, intent(in)    :: optics_dir
        integer      :: ipart, lastmap
        if( self%iter <= 0 )then
            ! no 2D yet
            call write_raw_project
        else
            if( .not.self%l_available )then
                self%iter = self%iter-1 ! iteration self%iter not complete so fall back on previous iteration
                if( self%iter <= 0 )then
                    ! no 2D yet
                    call write_raw_project
                else
                    self%refs = CAVGS_ITER_FBODY//int2str_pad(self%iter,3)//MRC_EXT
                    ! tricking the asynchronous master process to come to a hard stop
                    call simple_touch(POOL_DIR//TERM_STREAM)
                    do ipart = 1,params%nparts_pool
                        call simple_touch(POOL_DIR//JOB_FINISHED_FBODY//int2str_pad(ipart,self%numlen))
                    enddo
                    call simple_touch(POOL_DIR//'CAVGASSEMBLE_FINISHED')
                endif
            endif
            if( self%iter >= 1 )then
                call self%write_project(params, write_star=.true., clspath=.true., optics_dir=optics_dir)
                call self%rank_cavgs()
            endif
        endif
        ! cleanup
        call del_file(POOL_DIR//POOL_PROJFILE)
        if( .not. DEBUG_HERE )then
            call qsys_cleanup(params)
        endif

        contains

            ! no complete pool iteration: the pool's imported micrographs, stacks and particles as they
            ! came (no 2D parameters but their shifts), with the newest optics map's groups, applied
            ! to the project as written (not before it is assembled, nor from a stub)
            subroutine write_raw_project
                type(string) :: projfile
                if( self%proj%os_ptcl2D%get_noris() == 0 ) return
                self%proj%os_ptcl3D = self%proj%os_ptcl2D
                call self%proj%os_ptcl3D%clean_entry('updatecnt', 'sampled')
                if( present(optics_dir) ) lastmap = import_latest_optics_map(self%proj, optics_dir)
                projfile = get_fbody(self%orig_projfile, METADATA_EXT, separator=.false.)//METADATA_EXT
                call self%proj%projinfo%set(1, 'projname', get_fbody(self%orig_projfile, METADATA_EXT, separator=.false.))
                call self%proj%projinfo%set(1, 'projfile', projfile)
                call self%starproj_stream%stream_export_micrographs(params, self%proj, params%cwd, optics_set=.true.)
                call self%starproj_stream%stream_export_particles_2D(params, self%proj, params%cwd, optics_set=.true.)
                call self%proj%write(projfile)
                write(logfhandle,'(A,A,A,I8,A)') '>>> NO COMPLETE POOL ITERATION; RAW PROJECT ', projfile%to_char(), ' WITH ',&
                    &self%proj%os_ptcl2D%count_state_gt_zero(), ' SELECTED PARTICLES'
                call self%proj%os_ptcl3D%kill
            end subroutine write_raw_project

    end subroutine finalise

    !> Writes the pool's project (the stage's final project): the class averages and FRCs at the
    !! native sampling, the newest optics map's groups (with @p optics_dir), ptcl3D prepared for
    !! 3D, the class STAR file and, with @p write_star, the micrograph and particle STAR files.
    subroutine write_project( self, params, write_star, clspath, optics_dir )
        class(stream_pool2D),    intent(inout) :: self
        class(parameters),       intent(inout) :: params
        logical,       optional, intent(in)    :: write_star
        logical,       optional, intent(in)    :: clspath
        class(string), optional, intent(in)    :: optics_dir
        type(class_frcs)        :: frcs
        type(oris)              :: os_backup
        type(starproject)       :: starproj
        type(string)            :: projfile, projfname, cavgsfname, frcsfname, pool_refs, frcs_src
        integer(timer_int_kind) :: t
        integer                 :: lastmap
        logical                 :: l_write_star, l_clspath
        l_write_star = .false.
        l_clspath    = .false.
        if(present(write_star)) l_write_star = write_star
        if(present(clspath))    l_clspath    = clspath
        ! file naming
        projfname  = get_fbody(self%orig_projfile, METADATA_EXT, separator=.false.)
        cavgsfname = get_fbody(self%refs, MRC_EXT, separator=.false.)
        frcsfname  = get_fbody(FRCS_FILE, BIN_EXT, separator=.false.)
        call self%proj%projinfo%set(1,'projname', projfname)
        projfile   = projfname//METADATA_EXT
        call self%proj%projinfo%set(1,'projfile', projfile)
        cavgsfname = cavgsfname//MRC_EXT
        frcsfname  = frcsfname//BIN_EXT
        pool_refs  = string(POOL_DIR)//self%refs
        ! the FRCs of the iteration written: its kept copy when there is one (a stop mid-iteration
        ! falls back on the previous iteration while frcs.bin is being rewritten), frcs.bin otherwise
        frcs_src   = string(POOL_DIR)//swap_suffix(FRCS_FILE, "_iter"//int2str_pad(self%iter, 3)//".bin", ".bin")
        if( .not. file_exists(frcs_src) ) frcs_src = string(POOL_DIR)//FRCS_FILE
        lastmap    = 0
        write(logfhandle,'(A,A,A,A)')'>>> WRITING PROJECT ', projfile%to_char(), ' AT: ',cast_time_char(simple_gettime())
        if( present(optics_dir) ) lastmap = import_latest_optics_map(self%proj, optics_dir)
        if( self%l_scaling )then
            os_backup = self%proj%os_cls2D
            call rescale_refs( cavgsfname )
            call self%proj%os_out%kill
            call self%proj%add_cavgs2os_out(cavgsfname, self%native_smpd, 'cavg', clspath=l_clspath)
            self%proj%os_cls2D = os_backup
            call os_backup%kill
            ! rescale frcs
            call frcs%read(frcs_src)
            call frcs%pad(self%native_smpd, self%native_box)
            call frcs%write(frcsfname)
            call frcs%kill
            call self%proj%add_frcs2os_out(frcsfname, 'frc2D')
        else
            call self%proj%os_out%kill
            call self%proj%add_cavgs2os_out(cavgsfname, self%native_smpd, 'cavg', clspath=l_clspath)
            if( .not. (frcs_src == frcsfname) ) call simple_copy_file(frcs_src, frcsfname)
            call self%proj%add_frcs2os_out(frcsfname, 'frc2D')
        endif
        ! the 3D field as the STAR files have it: 2D clustering removed, shifts kept
        self%proj%os_ptcl3D = self%proj%os_ptcl2D
        call self%proj%os_ptcl3D%delete_2Dclustering
        call self%proj%os_ptcl3D%clean_entry('updatecnt', 'sampled')
        call self%proj%write(projfile)
        ! write starfiles
        call starproj%export_cls2D(self%proj)
        call starproj%kill
        if(l_write_star) then
            if( DEBUG_HERE ) t = tic()
            call self%starproj_stream%stream_export_micrographs(params, self%proj, params%cwd, optics_set=.true.)
            if( DEBUG_HERE ) print *,'ms_export  : ', toc(t); call flush(6); t = tic()
            call self%starproj_stream%stream_export_particles_2D(params, self%proj, params%cwd, optics_set=.true.)
            if( DEBUG_HERE ) print *,'ptcl_export  : ', toc(t); call flush(6)
        end if
        call self%proj%os_ptcl3D%kill
        call self%proj%os_cls2D%delete_entry('stk')

        contains

        ! rescale classes to original scale
        subroutine rescale_refs( cavgs_fname )
            class(string), intent(in) :: cavgs_fname
            type(string) :: source, destination
            call self%rescale_cavgs(pool_refs, cavgs_fname)
            source  = add2fbody(pool_refs, MRC_EXT, '_even')
            destination = add2fbody(cavgs_fname, MRC_EXT, '_even')
            call self%rescale_cavgs(source, destination)
            source  = add2fbody(pool_refs, MRC_EXT,'_odd')
            destination = add2fbody(cavgs_fname, MRC_EXT,'_odd')
            call self%rescale_cavgs(source, destination)
        end subroutine rescale_refs

    end subroutine write_project

    !> Writes a snapshot of pool iteration @p iteration as @p projfile: the classes @p selection
    !! names (the others' particles deselected), with the class averages and FRCs of that
    !! iteration at the pool's sampling, the newest optics map's groups (when @p optics_dir is not
    !! empty) and micrograph and particle STAR files (@p starfile_base, optics group ids offset by
    !! @p optics_offset). The iteration is the current one or one of the history
    !! (POOL_NHISTORY). @p nptcls is the number of selected particles written; 0 when the
    !! iteration is no longer kept or its files are missing, and then nothing is written.
    !! @p cavgs_* report the selected class averages' sprite sheet for the GUI.
    subroutine write_snapshot( self, iteration, selection, projfile, starfile_base, optics_dir, optics_offset,&
            &nptcls, cavgs_jpeg, cavgs_mrc, cavgs_ntilesx, cavgs_ntilesy, cavgs_idx, cavgs_pop, cavgs_res )
        class(stream_pool2D), intent(inout) :: self
        integer,              intent(in)    :: iteration
        integer,              intent(in)    :: selection(:)
        class(string),        intent(in)    :: projfile, starfile_base, optics_dir
        integer,              intent(in)    :: optics_offset
        integer,              intent(out)   :: nptcls
        type(string),         intent(out)   :: cavgs_jpeg, cavgs_mrc
        integer,              intent(out)   :: cavgs_ntilesx, cavgs_ntilesy
        integer, allocatable, intent(out)   :: cavgs_idx(:), cavgs_pop(:)
        real,    allocatable, intent(out)   :: cavgs_res(:)
        type(sp_project) :: snapshot_proj
        type(string)     :: dir, stk, frcs, cavgsfname, frcsfname
        real             :: smpd
        integer          :: ncls, islot, lastmap
        logical          :: l_found
        nptcls        = 0
        cavgs_jpeg    = ''
        cavgs_mrc     = ''
        cavgs_ntilesx = 0
        cavgs_ntilesy = 0
        allocate(cavgs_idx(0), cavgs_pop(0), cavgs_res(0))
        write(logfhandle,'(A,I4,A,A,A,A)') '>>> WRITING SNAPSHOT FROM ITERATION ', iteration, ': ', projfile%to_char(),&
            &' AT: ', cast_time_char(simple_gettime())
        call log_rss('snapshot/before copy')
        ! the iteration: the current one, or one the history keeps
        l_found = .false.
        if( iteration == self%iter )then
            ! registered only when it exists (add_frcs2os_out requires the file)
            if( file_exists(string(POOL_DIR)//FRCS_FILE) )then
                call snapshot_proj%copy(self%proj)
                call snapshot_proj%add_frcs2os_out(string(POOL_DIR)//FRCS_FILE, 'frc2D')
                l_found = .true.
            endif
        else if( iteration >= 1 )then
            islot   = history_slot(iteration)
            l_found = self%history_iter(islot) == iteration
            if( l_found ) call snapshot_proj%copy(self%history(islot))
        endif
        ! and its files
        if( l_found )then
            call snapshot_proj%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
            call snapshot_proj%get_frcs(frcs, 'frc2D', fail=.false.)
            l_found = ncls > 0
            if( l_found ) l_found = file_exists(stk) .and. file_exists(frcs)
            if( l_found ) l_found = file_exists(add2fbody(stk, MRC_EXT, '_even'))
            if( l_found ) l_found = file_exists(add2fbody(stk, MRC_EXT, '_odd'))
        endif
        if( .not. l_found )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> WARNING: ITERATION ', iteration, ' IS NOT KEPT (THE POOL KEEPS THE LAST ',&
                &POOL_NHISTORY, '); SNAPSHOT NOT WRITTEN'
            call snapshot_proj%kill
            return
        endif
        call log_rss('snapshot/after copy')
        dir = stemname(projfile)
        if( .not. file_exists(stemname(dir)) ) call simple_mkdir(stemname(dir))
        if( .not. file_exists(dir) )           call simple_mkdir(dir)
        if( optics_dir%strlen() > 0 ) lastmap = import_latest_optics_map(snapshot_proj, optics_dir)
        call self%apply_snapshot_selection(snapshot_proj, selection)
        call log_rss('snapshot/after selection')
        cavgsfname = dir//'/cavgs'//STK_EXT
        frcsfname  = dir//'/'//FRCS_FILE
        call simple_copy_file(stk, cavgsfname)
        call simple_copy_file(add2fbody(stk, MRC_EXT, '_even'), dir//'/cavgs_even'//STK_EXT)
        call simple_copy_file(add2fbody(stk, MRC_EXT, '_odd'),  dir//'/cavgs_odd'//STK_EXT)
        call simple_copy_file(frcs, frcsfname)
        call snapshot_proj%os_out%kill
        ! copied as they are, at the pool's (possibly cropped) sampling
        call snapshot_proj%add_cavgs2os_out(cavgsfname, smpd, 'cavg')
        call snapshot_proj%add_frcs2os_out(frcsfname, 'frc2D')
        call snapshot_proj%set_cavgs_thumb(projfile)
        call snapshot_cavgs_meta(snapshot_proj, cavgsfname, cavgs_jpeg, cavgs_mrc, cavgs_ntilesx, cavgs_ntilesy,&
            &cavgs_idx, cavgs_pop, cavgs_res)
        snapshot_proj%os_ptcl3D = snapshot_proj%os_ptcl2D
        call snapshot_proj%write(projfile)
        call snapshot_proj%write_mics_star(starfile_base//"_micrographs.star", optics_offset=optics_offset)
        call snapshot_proj%write_ptcl2D_star(starfile_base//"_particles.star", optics_offset=optics_offset)
        nptcls = snapshot_proj%os_ptcl2D%count_state_gt_zero()
        call log_rss('snapshot/after write')
        call snapshot_proj%kill
    end subroutine write_snapshot

    !> Publishes the pool's classified state for 3D as @p projfile (build_pool_publication, with the
    !! newest optics map of @p optics_dir), with the pool's class averages (and their even and odd
    !! halves) and FRCs at the native sampling beside it, registered with the pool's mask diameter.
    !! The class averages and FRCs are written first and the project last, under a temporary name
    !! renamed into place, so a reader that finds the project finds it complete. @p l_final marks
    !! the pool's final publication (pool_final=yes in its out segment), which starts multistate
    !! 3D's final run. @p nstks is the number of stacks published (0: nothing was written).
    subroutine publish( self, projfile, nstks, optics_dir, l_final )
        class(stream_pool2D), intent(inout) :: self
        class(string),        intent(in)    :: projfile
        integer,              intent(out)   :: nstks
        class(string),        intent(in)    :: optics_dir
        logical,              intent(in)    :: l_final
        type(sp_project) :: pub
        type(class_frcs) :: frcs
        type(string)     :: pool_refs, cavgsfname, frcsfname
        call build_pool_publication(self%proj, pub, nstks, optics_dir)
        if( nstks == 0 ) return
        call pool_publication_names(projfile, cavgsfname, frcsfname)
        pool_refs = string(POOL_DIR)//self%refs
        if( self%l_scaling )then
            call self%rescale_cavgs(pool_refs, cavgsfname)
            call self%rescale_cavgs(add2fbody(pool_refs, MRC_EXT, '_even'), add2fbody(cavgsfname, MRC_EXT, '_even'))
            call self%rescale_cavgs(add2fbody(pool_refs, MRC_EXT, '_odd'),  add2fbody(cavgsfname, MRC_EXT, '_odd'))
            call frcs%read(string(POOL_DIR)//FRCS_FILE)
            call frcs%pad(self%native_smpd, self%native_box)
            call frcs%write(frcsfname)
            call frcs%kill
        else
            call simple_copy_file(pool_refs, cavgsfname)
            call simple_copy_file(add2fbody(pool_refs, MRC_EXT, '_even'), add2fbody(cavgsfname, MRC_EXT, '_even'))
            call simple_copy_file(add2fbody(pool_refs, MRC_EXT, '_odd'),  add2fbody(cavgsfname, MRC_EXT, '_odd'))
            call simple_copy_file(string(POOL_DIR)//FRCS_FILE, frcsfname)
        endif
        call pub%os_out%kill
        call pub%add_cavgs2os_out(cavgsfname, self%native_smpd, 'cavg', mskdiam=self%mskdiam)
        call pub%add_frcs2os_out(frcsfname, 'frc2D')
        if( l_final ) call pub%os_out%set(1, 'pool_final', 'yes')
        call pub%write(projfile, tempfile=.true.)
        write(logfhandle,'(A,A,A,I8,A,I8,A)') '>>> PUBLISHED THE POOL FOR 3D: ', projfile%to_char(), ', ',&
            &nstks, ' STACK(S), ', pub%os_ptcl2D%count_state_gt_zero(), ' PARTICLE(S)'
        if( l_final ) write(logfhandle,'(A)') '>>> THE POOL''S FINAL PUBLICATION'
        call pub%kill
    end subroutine publish

    !---------------- queries ----------------

    !> The latest iteration dispatched (0 before the first)
    integer function iteration( self )
        class(stream_pool2D), intent(in) :: self
        iteration = self%iter
    end function iteration

    !> Whether the pool is available for another iteration
    logical function available( self )
        class(stream_pool2D), intent(in) :: self
        available = self%l_available
    end function available

    !> .true. once an iteration has failed twice: the pool stops
    logical function failed( self )
        class(stream_pool2D), intent(in) :: self
        failed = self%l_failed
    end function failed

    !> What the GUI shows of the pool: the latest complete iteration, the particles assigned and
    !! rejected, the resolution, the mask, and the class averages with their sprite sheet. The class
    !! averages' path is absolute; the pool's reference stays relative to the stage's directory,
    !! where its workers find it.
    function stats( self ) result( st )
        class(stream_pool2D), intent(in) :: self
        type(stream_pool2D_stats) :: st
        type(string) :: cwd
        st%last_complete_iter = self%last_complete_iter
        st%nassigned          = self%nptcls - self%nptcls_rejected
        st%nrejected          = self%nptcls_rejected
        st%resolution         = self%resolution
        st%mskdiam            = self%mskdiam
        st%msk_crop           = self%dims%msk
        st%cavgs_jpeg         = self%jpeg
        st%cavgs_mrc          = ''
        if( self%refs%strlen() > 0 )then
            call simple_getcwd(cwd)
            st%cavgs_mrc = cwd//'/'//self%refs
        endif
        st%jpeg_ntilesx = self%jpeg_ntilesx
        st%jpeg_ntilesy = self%jpeg_ntilesy
        if( allocated(self%jpeg_map) )then
            st%jpeg_map = self%jpeg_map
            st%jpeg_pop = self%jpeg_pop
            st%jpeg_res = self%jpeg_res
        endif
    end function stats

    !---------------- helpers ----------------

    !> The history slot of pool iteration @p iter
    pure integer function history_slot( iter )
        integer, intent(in) :: iter
        history_slot = modulo(iter - 1, POOL_NHISTORY) + 1
    end function history_slot

end module simple_stream_pool2D
