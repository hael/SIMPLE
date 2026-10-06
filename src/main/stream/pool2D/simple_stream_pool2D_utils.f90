!@descr: utilities for running the pool 2D refinement
module simple_stream_pool2D_utils
use simple_stream_api
use simple_qsys_job_record, only: cancel_queued_job
use simple_qsys_async_job,  only: qsys_async_job, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
implicit none

! CALCULATORS
public :: init_pool_clustering
public :: iterate_pool
public :: draw_new_classes
! GETTERS
public :: get_pool_assigned
public :: get_pool_cavgs_jpeg
public :: get_pool_cavgs_jpeg_ntilesx
public :: get_pool_cavgs_jpeg_ntilesy
public :: get_pool_cavgs_mrc
public :: get_pool_iter
public :: get_pool_ptr
public :: get_pool_rejected
public :: get_pool_resolution
public :: is_pool_available
public :: is_pool_failed
! SETTERS
public :: set_pool_resolution_limits
! UPDATERS
public :: update_mskdiam
public :: update_pool
public :: update_pool_aln_params
public :: update_pool_status
public :: cancel_pool_job
! JPGS / GUI
public :: generate_pool_stats
private
#include "simple_local_flags.inc"

integer                 :: current_jpeg_ntiles              ! number of used tiles in current JPEG
integer                 :: current_jpeg_ntilesx             ! number of tiles in x
integer                 :: current_jpeg_ntilesy             ! number of tiles in y
integer                 :: lim_ufrac_nptcls  = 0            ! threshold for fractional updates
integer                 :: ncls_rejected_glob               ! number of rejected classes
integer                 :: nptcls_glob                      ! total particles in pool
integer                 :: nptcls_rejected_glob             ! rejected particles in pool
logical,    allocatable :: pool_stacks_mask(:)              ! subset of stacks undergoing 2D analysis
real                    :: current_resolution=999.          ! current estimated resolution
real                    :: resolutions(POOL_NPREV_RES)=999. ! pool resolution history (length POOL_NPREV_RES)
type(qsys_env)          :: pool_qenv                        ! qsys submission environment for pool
integer                 :: pool_nattempts    = 0            ! submissions of the current iteration: a failed first is retried once
logical                 :: l_pool_failed     = .false.      ! the current iteration failed twice: the pool stops
type(qsys_async_job)    :: pool_job                         ! the current iteration's job (its files: POOL_EXIT_CODE...)
character(len=*), parameter :: POOL_JOB_LABEL = 'refine2D_pool' ! names its script, log and exit status as POOL_* do
type(string)            :: pool_center                      ! the pool's centering (yes|no) on a full update
type(string)            :: current_jpeg                     ! filename of current pool JPEG (type(string))
! convergence
real                    :: conv_frac     = 0.0
real                    :: conv_mi_class = 0.0
real                    :: conv_score    = 0.0

contains

    ! CALCULATORS

    ! The pool, at its first import: @p box and @p smpd are the imported sets' (px, A) and
    ! @p mskdiam (A) the pool's mask diameter, kept as pool state (simple_stream2D_state).
    subroutine init_pool_clustering( params, cline, spproj, box, smpd, mskdiam )
        class(parameters), intent(inout) :: params
        class(cmdline),    target, intent(inout) :: cline
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: box
        real,              intent(in)    :: smpd, mskdiam
        type(string)          :: carg, pool_sigma_path
        character(len=STDLEN) :: pool_part_env
        integer               :: envlen
        call seed_rnd
        ! general parameters
        master_cline => cline
        pool_native_box  = box
        pool_native_smpd = smpd
        pool_mskdiam     = mskdiam
        pool_user_lpstop = 0.
        if( cline%defined('lpstop') ) pool_user_lpstop = params%lpstop
        call mskdiam2lplimits(pool_mskdiam, lpstart, lpstop, lpcen)
        l_scaling          = .true.
        ncls_glob          = params%ncls
        ncls_rejected_glob = 0
        orig_projfile      = params%projfile
        params%nparts_pool = params%nparts ! backwards compatibility
        ! bookkeeping & directory structure
        numlen             = len(int2str(params%nparts))
        refs_glob          = ''
        l_pool_available   = .true.
        pool_iter          = 0
        call simple_mkdir(POOL_DIR, verbose=.false.)
        call simple_mkdir(POOL_DIR//STDERROUT_DIR)
        call simple_mkdir(DIR_SNAPSHOT)
        pool_proj%projinfo = spproj%projinfo
        pool_proj%compenv  = spproj%compenv
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        call pool_proj%projinfo%delete_entry('sigma2_state')
        pool_sigma_path = string(POOL_DIR)//'sigma2_state.bin'
        call pool_proj%set_sigma2_state_path(pool_sigma_path)
        call pool_sigma_path%kill
        ! update to computational parameters to pool, will be transferred to chunks upon init
        if( cline%defined('walltime') ) call pool_proj%compenv%set(1,'walltime', params%walltime)
        ! commit to disk
        call pool_proj%write(string(POOL_DIR)//POOL_PROJFILE)
        ! Pool command line
        call cline_refine2D_pool%set('prg',        'refine2D_distr')
        call cline_refine2D_pool%set('oritype',    'ptcl2D')
        call cline_refine2D_pool%set('trs',        MINSHIFT)
        call cline_refine2D_pool%set('projfile',   POOL_PROJFILE)
        call cline_refine2D_pool%set('projname',   get_fbody(POOL_PROJFILE,'simple'))
        call cline_refine2D_pool%set('sigma_est', params%sigma_est)
        if( cline%defined('cls_init') )then
            call cline_refine2D_pool%set('cls_init', params%cls_init)
        else
            call cline_refine2D_pool%set('cls_init', 'rand')
        endif
        pool_center = 'yes'
        if( cline%defined('center') )then
            carg        = cline%get_carg('center')
            pool_center = carg
            call carg%kill
        endif
        call cline_refine2D_pool%set('center', pool_center)
        if( cline%defined('center_type') )then
            call cline_refine2D_pool%set('center_type', params%center_type)
        else
            call cline_refine2D_pool%set('center_type', 'seg')
        endif
        call cline_refine2D_pool%set('extr_iter', 99999)
        call cline_refine2D_pool%set('extr_lim',   MAX_EXTRLIM2D)
        call cline_refine2D_pool%set('mkdir',      'no')
        call cline_refine2D_pool%set('mskdiam',    pool_mskdiam)
        call cline_refine2D_pool%set('async',      'yes') ! to enable hard termination
        call cline_refine2D_pool%set('stream2d',   'yes') ! the only place this flag should be turned on
        call cline_refine2D_pool%set('nparts',     params%nparts)
        if( cline%defined('worker_server') ) call cline_refine2D_pool%set('worker_server', cline%get_carg('worker_server'))
        if( cline%defined('worker_server_nthr') ) call cline_refine2D_pool%set('worker_server_nthr', cline%get_iarg('worker_server_nthr'))
        call cline_refine2D_pool%delete('autoscale')
        ! when the 2D analysis is started from raw particles
        ! set # of ptcls beyond which fractional updates will be used
        lim_ufrac_nptcls = STREAM_NPTCLS_MAX
        if( master_cline%defined('nsample_max') ) lim_ufrac_nptcls = params%nsample_max
        ! the iterations' threads (the master's resources table, SIMPLE_STREAM_POOL_NTHR over it)
        call cline_refine2D_pool%set('nthr', params%nthr)
        call get_environment_variable(SIMPLE_STREAM_POOL_PARTITION, pool_part_env, envlen)
        if(envlen > 0) then
            call pool_qenv%new(params, params%nparts, stream=.true., exec_bin=string('simple_private_exec'),qsys_name=string('local'),&
            &qsys_partition=string(trim(pool_part_env)))
        else
            call pool_qenv%new(params, params%nparts, stream=.true., exec_bin=string('simple_private_exec'),qsys_name=string('local'))
        end if
        ! objective function
        call cline_refine2D_pool%set('objfun', 'euclid')
        call cline_refine2D_pool%set('ml_reg', params%ml_reg)
        call cline_refine2D_pool%set('tau',     params%tau)
        ! refinement
        select case(trim(params%refine))
            case('snhc','snhc_smpl')
                call cline_refine2D_pool%set( 'refine', params%refine)
            case DEFAULT
                THROW_HARD('UNSUPPORTED REFINE PARAMETER!')
        end select
        ! Determines dimensions for downscaling
        call set_pool_dimensions()
        ! updates command-lines with resolution limits
        call set_pool_resolution_limits(params)
        ! module variables
        l_stream2D_active = .true.
    end subroutine init_pool_clustering

    ! Performs one iteration:
    ! updates to command-line, particles sampling, temporary project & execution
    subroutine iterate_pool( params )
        class(parameters), intent(inout) :: params
        logical, parameter            :: L_BENCH = .false.
        type(sp_project)              :: spproj
        integer(timer_int_kind)       :: t_tot
        integer,          allocatable :: nptcls_per_stk(:)
        type(string) :: stkname
        real         :: frac_update, smpd
        integer      :: iptcl,i, nptcls_tot, fromp, top, nstks_tot, jptcl, islot
        integer      :: nptcls_sel, istk, nptcls2update, nstks2update, jjptcl, ncls
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available  ) return
        if( L_BENCH ) t_tot  = tic()
        nptcls_tot           = pool_proj%os_ptcl2D%get_noris()
        nptcls_glob          = nptcls_tot
        nptcls_rejected_glob = 0
        if( nptcls_tot == 0 ) return
        ! the completed iteration goes into the history, replacing the iteration POOL_NHISTORY
        ! before it, whose files tidy_2Dstream_iter removes below: history and files keep the
        ! same iterations
        if(pool_iter .gt. 0) then
            islot = pool_history_slot(pool_iter)
            call pool_proj_history(islot)%kill
            call pool_proj_history(islot)%copy(pool_proj)
            call simple_copy_file(string(POOL_DIR)//FRCS_FILE, string(POOL_DIR)//swap_suffix(FRCS_FILE,"_iter"//int2str_pad(pool_iter, 3)//".bin",".bin"))
            call pool_proj_history(islot)%get_cavgs_stk(stkname, ncls, smpd)
            call pool_proj_history(islot)%os_out%kill
            call pool_proj_history(islot)%add_cavgs2os_out(stkname, smpd, 'cavg')
            call pool_proj_history(islot)%add_frcs2os_out(string(POOL_DIR)//swap_suffix(FRCS_FILE,"_iter"//int2str_pad(pool_iter, 3)//".bin",".bin"),'frc2D')
            pool_history_iter(islot) = pool_iter
        end if
        pool_iter = pool_iter + 1 ! Global iteration counter update
        call cline_refine2D_pool%set('ncls',     ncls_glob)
        call cline_refine2D_pool%set('startit', pool_iter)
        call cline_refine2D_pool%set('maxits',   pool_iter)
        call cline_refine2D_pool%set('frcs',     FRCS_FILE)
        call cline_refine2D_pool%set('refs', refs_glob)
        if( pool_iter==1 )then
            if( cline_refine2D_pool%defined('cls_init') )then
                ! references taken care of by refine2D_distr
                call cline_refine2D_pool%delete('frcs')
                call cline_refine2D_pool%delete('refs')
            endif
        else
            call cline_refine2D_pool%delete('cls_init')
        endif
        ! Project metadata update
        spproj%projinfo = pool_proj%projinfo
        spproj%compenv  = pool_proj%compenv
        call spproj%projinfo%delete_entry('projname')
        call spproj%projinfo%delete_entry('projfile')
        call spproj%update_projinfo( cline_refine2D_pool )
        ! Sampling of stacks that will be used for this iteration
        ! counting number of stacks & selected particles
        nstks_tot  = pool_proj%os_stk%get_noris()
        allocate(nptcls_per_stk(nstks_tot), source=0)
        !$omp parallel do schedule(static) proc_bind(close) private(istk,fromp,top,iptcl) default(shared)
        do istk = 1,nstks_tot
            fromp = pool_proj%os_stk%get_fromp(istk)
            top   = pool_proj%os_stk%get_top(istk)
            do iptcl = fromp,top
                if( pool_proj%os_ptcl2D%get_state(iptcl) > 0 )then
                    nptcls_per_stk(istk)  = nptcls_per_stk(istk) + 1 ! # ptcls with state=1
                endif
            enddo
        enddo
        !$omp end parallel do
        nptcls_rejected_glob = nptcls_glob - sum(nptcls_per_stk)
        ! Update info for gui
        call spproj%projinfo%set(1,'nptcls_tot',     nptcls_glob)
        call spproj%projinfo%set(1,'nptcls_rejected',nptcls_rejected_glob)
        ! Uniformly sample stacks
        call uniform_stack_sampling
        nstks2update = count(pool_stacks_mask)
        ! Transfer stacks and particles
        call spproj%os_stk%new(nstks2update, is_ptcl=.false.)
        call spproj%os_ptcl2D%new(nptcls2update, is_ptcl=.true.)
        i     = 0
        jptcl = 0
        do istk = 1,nstks_tot
            fromp = pool_proj%os_stk%get_fromp(istk)
            top   = pool_proj%os_stk%get_top(istk)
            if( pool_stacks_mask(istk) )then
                ! transfer alignement parameters for selected particles
                i = i + 1 ! stack index in spproj
                call spproj%os_stk%transfer_ori(i, pool_proj%os_stk, istk)
                call spproj%os_stk%set(i, 'fromp', jptcl+1)
                !$omp parallel do private(iptcl,jjptcl) proc_bind(close) default(shared)
                do iptcl = fromp,top
                    jjptcl = jptcl+iptcl-fromp+1
                    call spproj%os_ptcl2D%transfer_ori(jjptcl, pool_proj%os_ptcl2D, iptcl)
                    call spproj%os_ptcl2D%set_stkind(jjptcl, i)
                enddo
                !$omp end parallel do
                jptcl = jptcl + (top-fromp+1)
                call spproj%os_stk%set(i, 'top', jptcl)
            endif
        enddo
        call spproj%os_ptcl3D%new(nptcls2update, is_ptcl=.true.)
        spproj%os_cls2D = pool_proj%os_cls2D
        ! the new particles get a populated class, drawn reproducibly
        if( pool_iter >= 2 ) call draw_new_classes(spproj, nptcls2update, pool_iter, ncls_glob)
        ! the sampled stacks are the update set (decision 5): every particle of the sample is
        ! updated, the others keep their parameters; only the user's update_frac thins the sample,
        ! and then the class averages are not centered
        call cline_refine2D_pool%delete('update_frac')
        frac_update = 1.0
        if( master_cline%defined('update_frac') ) frac_update = params%update_frac
        if( frac_update < 0.99999 )then
            call cline_refine2D_pool%set('update_frac', frac_update)
            call cline_refine2D_pool%set('center',      'no')
        else
            call cline_refine2D_pool%set('center',      pool_center)
        endif
        ! write project, and keep it as made for a retry
        call spproj%write(string(POOL_DIR)//POOL_PROJFILE)
        call spproj%kill
        call simple_copy_file(string(POOL_DIR)//POOL_PROJFILE, string(POOL_DIR)//POOL_INPUT_PROJFILE)
        ! pool stats
        call generate_pool_stats(params)
        ! execution
        pool_nattempts = 0
        call submit_pool_iteration(params)
        write(logfhandle,'(A,I6,A,I8,A3,I8,A)')'>>> POOL         INITIATED ITERATION ',pool_iter,' WITH ',nptcls_sel,&
        &' / ', sum(nptcls_per_stk),' PARTICLES'
        if( L_BENCH ) print *,'timer analyze2D_pool tot : ',toc(t_tot)
        ! cleanup
        if( allocated(nptcls_per_stk) )      deallocate(nptcls_per_stk)
        ! the files of the iteration that has just left the history
        call tidy_2Dstream_iter(pool_iter - 1 - POOL_NHISTORY)

      contains

        subroutine uniform_stack_sampling
            use simple_ran_tabu
            type(ran_tabu) :: random_generator
            integer        :: stk_order(nstks_tot)
            integer        :: i, j
            if( allocated(pool_stacks_mask) ) deallocate(pool_stacks_mask)
            allocate(pool_stacks_mask(nstks_tot), source=.false.)
            stk_order        = (/(i,i=1,nstks_tot)/)
            random_generator = ran_tabu(nstks_tot)
            call random_generator%shuffle(stk_order)
            nptcls2update = 0 ! # of ptcls including state=0 within selected stacks
            nptcls_sel    = 0 ! # of ptcls excluding state=1 within selected stacks
            do i = 1,nstks_tot
                j = stk_order(i)
                if( nptcls_sel > lim_ufrac_nptcls ) cycle
                nptcls_sel    = nptcls_sel    + nptcls_per_stk(j)
                nptcls2update = nptcls2update + pool_proj%os_stk%get_int(j, 'nptcls')
                pool_stacks_mask(j) = .true.
            enddo
            call random_generator%kill
        end subroutine uniform_stack_sampling

    end subroutine iterate_pool

    ! Classes for the never-updated particles among the first @p nptcls of @p spproj, drawn among
    ! its populated classes (of the first @p ncls) from a generator seeded with the iteration
    ! @p iter, in one thread, so a run is reproducible; the process's generator state is restored
    ! after. With no populated class they keep the class they have.
    subroutine draw_new_classes( spproj, nptcls, iter, ncls )
        use iso_fortran_env, only: int64
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: nptcls, iter, ncls
        integer, allocatable :: clspops(:), populated(:), saved_seed(:), seed(:)
        integer :: ncls_drawn, iptcl, nseed, i
        clspops    = spproj%os_cls2D%get_all_asint('pop')
        ncls_drawn = min(ncls, size(clspops))
        if( ncls_drawn < 1 ) return
        populated  = pack([(i, i=1,ncls_drawn)], clspops(:ncls_drawn) > 0)
        if( size(populated) == 0 ) return
        call random_seed(size=nseed)
        allocate(saved_seed(nseed), seed(nseed))
        call random_seed(get=saved_seed)
        do i = 1,nseed
            seed(i) = int(modulo(int(iter,int64) * 7919_int64 + 104729_int64 * int(i-1,int64), int(huge(0)-1,int64)) + 1_int64)
        enddo
        call random_seed(put=seed)
        do iptcl = 1,min(nptcls, spproj%os_ptcl2D%get_noris())
            if( spproj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
            if( spproj%os_ptcl2D%get_updatecnt(iptcl) /= 0 ) cycle
            call spproj%os_ptcl2D%set_class(iptcl, populated(irnd_uni(size(populated))))
        enddo
        call random_seed(put=saved_seed)
    end subroutine draw_new_classes

    ! Submits the pool's current iteration as the pool's job (qsys_async_job, in the stage's folder:
    ! POOL_DIR is empty), which removes the status and job record of an earlier attempt; the attempt
    ! is counted.
    subroutine submit_pool_iteration( params )
        class(parameters), intent(in) :: params
        type(cmdline), allocatable :: pool_clines(:)
        type(string)               :: cwd
        call simple_getcwd(cwd)
        pool_nattempts = pool_nattempts + 1
        if( params%cc_objfun == OBJFUN_EUCLID )then
            ! Stream pool membership changes invalidate row identity. Rebuild a
            ! complete canonical bootstrap for the exact pool layout, then run
            ! clustering in the same queued script so no consumer can observe
            ! a missing or stale state.
            allocate(pool_clines(2))
            pool_clines(1) = cline_refine2D_pool
            call pool_clines(1)%set('prg', 'calc_pspec')
            call pool_clines(1)%set('mkdir', 'no')
            call pool_clines(1)%delete('stream2d')
            call pool_clines(1)%delete('update_frac')
            pool_clines(2) = cline_refine2D_pool
            call pool_job%start_seq(pool_qenv, pool_clines, cwd, POOL_JOB_LABEL)
            call pool_clines(:)%kill
            deallocate(pool_clines)
        else
            call pool_job%start(pool_qenv, cline_refine2D_pool, cwd, POOL_JOB_LABEL)
        endif
        l_pool_available = .false.
    end subroutine submit_pool_iteration

    ! .true. when the current iteration's job has ended without refine2D finishing: it exited,
    ! whatever its status, or vanished without one, which the job's liveness check finds
    ! (qsys_async_job: a walltime kill, a lost node).
    logical function pool_job_failed()
        integer :: job_status
        pool_job_failed = .false.
        job_status = pool_job%status()
        if( job_status /= ASYNC_JOB_DONE .and. job_status /= ASYNC_JOB_FAILED ) return
        ! refine2D touches its marker before it exits and the script writes the status after
        if( file_exists(POOL_DIR//REFINE2D_FINISHED) ) return
        pool_job_failed = .true.
    end function pool_job_failed

    ! Cancels the job of the pool's current iteration when it recorded itself and has not exited
    ! (simple_qsys_job_record): on stop, and on a restart for a job a crashed stage left running.
    subroutine cancel_pool_job
        if( cancel_queued_job(string(POOL_DIR//POOL_EXIT_CODE)) )then
            write(logfhandle,'(A)') '>>> CANCELLED THE RUNNING POOL ITERATION'
        endif
    end subroutine cancel_pool_job

    ! GETTERS

    ! .true. once an iteration has failed twice: the pool stops
    logical function is_pool_failed()
        is_pool_failed = l_pool_failed
    end function is_pool_failed

    ! returns number currently assigned particles
    integer function get_pool_assigned()
        get_pool_assigned = nptcls_glob - nptcls_rejected_glob
    end function get_pool_assigned

    ! absolute path of the current pool class-average JPEG ('' before the first)
    type(string) function get_pool_cavgs_jpeg()
        get_pool_cavgs_jpeg = current_jpeg
    end function get_pool_cavgs_jpeg

    integer function get_pool_cavgs_jpeg_ntilesx()
        get_pool_cavgs_jpeg_ntilesx = current_jpeg_ntilesx
    end function get_pool_cavgs_jpeg_ntilesx

    integer function get_pool_cavgs_jpeg_ntilesy()
        get_pool_cavgs_jpeg_ntilesy = current_jpeg_ntilesy
    end function get_pool_cavgs_jpeg_ntilesy

    ! absolute path of the current pool class averages ('' before the first); refs_glob itself
    ! stays relative to the stage's directory, where the pool workers find it
    type(string) function get_pool_cavgs_mrc()
        type(string) :: cwd
        get_pool_cavgs_mrc = ''
        if( refs_glob%strlen() == 0 ) return
        call simple_getcwd(cwd)
        get_pool_cavgs_mrc = cwd//'/'//refs_glob
    end function get_pool_cavgs_mrc

    ! returns current pool iteration
    integer function get_pool_iter()
        get_pool_iter = pool_iter
    end function get_pool_iter

    ! returns pointer to pool project
    subroutine get_pool_ptr( ptr )
        class(sp_project), pointer, intent(out) :: ptr
        ptr => pool_proj
    end subroutine get_pool_ptr

    ! returns number currently rejected particles
    integer function get_pool_rejected()
        get_pool_rejected = nptcls_rejected_glob
    end function get_pool_rejected

    real function get_pool_resolution()
        get_pool_resolution = current_resolution
    end function get_pool_resolution

    ! whether the pool available for another iteration
    logical function is_pool_available()
        is_pool_available = l_pool_available
    end function is_pool_available

    ! SETTERS

    ! The pool's working dimensions from its native ones: downscaled to a pixel size of up to
    ! MAX_SMPD and never below a CHUNK_MINBOXSZ box (setup_downscaling's rule). The pool's command
    ! line carries the cropped dimensions only; the native ones come from its project.
    subroutine set_pool_dimensions
        real    :: smpd, scale_factor
        integer :: box
        if( pool_native_box == 0 ) THROW_HARD('the pool has no native box; set_pool_dimensions')
        pool_dims%smpd = pool_native_smpd
        pool_dims%box  = pool_native_box
        if( l_scaling .and. pool_native_box >= CHUNK_MINBOXSZ )then
            call autoscale(pool_native_box, pool_native_smpd, MAX_SMPD, box, smpd, scale_factor, minbox=CHUNK_MINBOXSZ)
            l_scaling = box < pool_native_box
            if( l_scaling )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',pool_native_box,'/',box
                pool_dims%smpd = smpd
                pool_dims%box  = box
            endif
        endif
        pool_dims%boxpd = 2*round2even(KBALPHA*real(pool_dims%box/2)) ! logics from parameters
        call set_pool_mask()
        ! chunk & pool have the same dimensions to start with (used for import)
        chunk_dims = pool_dims
        ! Scaling-related command lines update
        call cline_refine2D_pool%set('smpd_crop',   pool_dims%smpd)
        call cline_refine2D_pool%set('box_crop',    pool_dims%box)
    end subroutine set_pool_dimensions

    ! The pool's mask radius at its working dimensions from its mask diameter, clamped (and logged)
    ! to (box - COSMSKHALFWIDTH)/2 pixels, so a diameter beyond the box never reaches the workers
    ! (D40); the diameter and the radius go on the pool's command line.
    subroutine set_pool_mask
        real :: msk_max
        msk_max       = (real(pool_dims%box) - COSMSKHALFWIDTH) / 2.
        pool_dims%msk = round2even(pool_mskdiam / pool_dims%smpd / 2.)
        if( real(pool_dims%msk) > msk_max )then
            write(logfhandle,'(A,F8.2,A,F8.2,A)') '>>> MASK DIAMETER ', pool_mskdiam, ' A EXCEEDS THE POOL''S BOX; CLAMPED TO ',&
                &2. * msk_max * pool_dims%smpd, ' A'
            pool_dims%msk = floor(msk_max)
            pool_mskdiam  = 2. * msk_max * pool_dims%smpd
        endif
        call cline_refine2D_pool%set('mskdiam',  pool_mskdiam)
        call cline_refine2D_pool%set('msk_crop', pool_dims%msk)
    end subroutine set_pool_mask

    ! private routine for pool resolution-related updates to command-lines
    subroutine set_pool_resolution_limits( params )
        class(parameters), intent(inout) :: params
        lpstart     = max(lpstart, 2.0*pool_dims%smpd)
        pool_lpstop = max(2.0*pool_dims%smpd, pool_user_lpstop)
        call cline_refine2D_pool%set('lpstart',   lpstart)
        call cline_refine2D_pool%set('lpstop',    pool_lpstop)
        if( .not.master_cline%defined('cenlp') )then
            call cline_refine2D_pool%set( 'cenlp', lpcen)
        else
            call cline_refine2D_pool%set( 'cenlp', params%cenlp)
        endif
        write(logfhandle,'(A,F5.1)') '>>> STARTING LOW-PASS LIMIT  (IN A): ', lpstart
        write(logfhandle,'(A,F5.1)') '>>> HARD RESOLUTION LIMIT    (IN A): ', pool_lpstop
        write(logfhandle,'(A,F5.1)') '>>> CENTERING LOW-PASS LIMIT (IN A): ', lpcen
    end subroutine set_pool_resolution_limits

    ! UPDATERS

    ! A new mask diameter (A) for the next pool iterations. The pool command line also carries
    ! the cropped mask radius (pixels), which the workers take over the one parameters would
    ! derive from mskdiam, so it is updated with it. Before the pool starts, the stage keeps the
    ! diameter and gives it to init_pool_clustering.
    subroutine update_mskdiam( new_mskdiam )
        integer, intent(in) :: new_mskdiam
        write(*,'(A,I4,A)')'>>> UPDATED MASK DIAMETER TO', new_mskdiam ,'Å'
        pool_mskdiam = real(new_mskdiam)
        call cline_refine2D_pool%set('mskdiam',   pool_mskdiam)
        if( pool_dims%smpd > 0. )then
            call set_pool_mask()
            ! the low-pass ramp and the centering limit follow the mask
            call mskdiam2lplimits(pool_mskdiam, lpstart, lpstop, lpcen)
            lpstart = max(lpstart, 2.0*pool_dims%smpd)
            call cline_refine2D_pool%set('lpstart', lpstart)
            if( associated(master_cline) )then
                if( .not. master_cline%defined('cenlp') ) call cline_refine2D_pool%set('cenlp', lpcen)
            endif
            write(logfhandle,'(A,F5.1,A,F5.1,A)') '>>> LOW-PASS RAMP FROM ', lpstart, ' A, CENTERING LOW-PASS ', lpcen, ' A'
        endif
    end subroutine update_mskdiam

    ! Reports alignment info from completed iteration of subset
    ! of particles back to the pool
    subroutine update_pool( params )
        class(parameters), intent(inout) :: params
        integer,      allocatable :: pops(:)
        type(sp_project) :: spproj
        type(oris)       :: os
        type(class_frcs) :: frcs
        type(string)     :: fname, cwd
        integer          :: i, it, jptcl, iptcl, istk
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available  ) return
        call del_file(POOL_DIR//REFINE2D_FINISHED)
        ! iteration info
        fname = POOL_DIR//STATS_FILE
        if( file_exists(fname) )then
            call os%new(1,is_ptcl=.false.)
            call os%read(fname)
            it = os%get_int(1,'ITERATION')
            if( it == pool_iter )then
                conv_mi_class = os%get(1,'CLASS_OVERLAP')
                conv_frac     = os%get(1,'SEARCH_SPACE_SCANNED')
                conv_score    = os%get(1,'SCORE')
                ! new
                last_complete_iter = it
                call simple_getcwd(cwd)
                current_jpeg = cwd//'/'//CAVGS_ITER_FBODY//int2str_pad(it, 3)//'.jpg'
                current_jpeg_ntiles  = pool_proj%os_cls2D%get_noris()
                ! no classes yet after the first iteration (they are transferred below)
                current_jpeg_ntilesx = max(1, floor(sqrt(real(current_jpeg_ntiles))))
                current_jpeg_ntilesy = ceiling(real(current_jpeg_ntiles)/real(current_jpeg_ntilesx))
                ! end new
                write(logfhandle,'(A,I6,A,F7.3,A,F7.3,A,F7.3)')'>>> POOL         ITERATION ',it,&
                    &'; CLASS OVERLAP: ',conv_mi_class,'; SEARCH SPACE SCANNED: ',conv_frac,'; SCORE: ',conv_score
            endif
            call os%kill
        endif
        ! transfer to pool
        call spproj%read_segment('cls2D', string(POOL_DIR)//POOL_PROJFILE)
        if( spproj%os_cls2D%get_noris() == 0 )then
            ! not executed yet, do nothing
        else
            if( .not.allocated(pool_stacks_mask) )then
                THROW_HARD('Critical ERROR 0') ! first time
            endif
            ! transfer particles parameters
            call spproj%read_segment('stk',   string(POOL_DIR)//POOL_PROJFILE)
            call spproj%read_segment('ptcl2D',string(POOL_DIR)//POOL_PROJFILE)
            i = 0
            do istk = 1,size(pool_stacks_mask)
                if( pool_stacks_mask(istk) )then
                    i = i+1
                    iptcl = pool_proj%os_stk%get_fromp(istk)
                    do jptcl = spproj%os_stk%get_fromp(i),spproj%os_stk%get_top(i)
                        if( spproj%os_ptcl2D%get_state(jptcl) > 0 )then
                            call pool_proj%os_ptcl2D%transfer_2Dparams(iptcl, spproj%os_ptcl2D, jptcl)
                        endif
                        iptcl = iptcl+1
                    enddo
                endif
            enddo
            ! update classes info
            call pool_proj%os_ptcl2D%get_pops(pops, 'class', maxn=ncls_glob)
            pool_proj%os_cls2D = spproj%os_cls2D
            call pool_proj%os_cls2D%set_all('pop', real(pops))
            ! update thumbnail metadata
            if(allocated(pool_jpeg_map)) deallocate(pool_jpeg_map)
            if(allocated(pool_jpeg_pop)) deallocate(pool_jpeg_pop)
            if(allocated(pool_jpeg_res)) deallocate(pool_jpeg_res)
            pool_jpeg_pop = pool_proj%os_cls2D%get_all_asint('pop')
            pool_jpeg_res = pool_proj%os_cls2D%get_all('res')
            allocate(pool_jpeg_map, mold=pool_jpeg_pop)
            do i=1, size(pool_jpeg_map)
                pool_jpeg_map(i) = i
            end do
            ! estimate resolution
            call frcs%read(string(POOL_DIR)//FRCS_FILE)
            current_resolution = frcs%estimate_lp_for_align()
            write(logfhandle,'(A,F5.1)')'>>> CURRENT POOL RESOLUTION: ',current_resolution
            call frcs%kill
            ! deal with dimensions/resolution update
            call update_pool_dims(params)
            ! for gui
            call update_pool_for_gui(params)
        endif
        call spproj%kill
    end subroutine update_pool

    ! This controls the evolution of the pool alignement parameters:
    ! lp, Gaussian filter, trs, extr_iter
    subroutine update_pool_aln_params
        integer, parameter :: ITERLIM    = 20
        integer, parameter :: ITERSHIFT  = 5
        real :: lp, gamma
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available  ) return
        if( pool_iter < ITERLIM )then
            gamma = min(1., max(0., real(ITERLIM-pool_iter)/real(ITERLIM)))
            ! offset
            if( pool_iter < ITERSHIFT )then
                call cline_refine2D_pool%set('trs', 0.)
            else
                call cline_refine2D_pool%set('trs', MINSHIFT)
            endif
            ! resolution limit
            lp = lpstop + (lpstart-lpstop) * gamma
            call cline_refine2D_pool%set('lp', lp)
            ! Extremal iteration
            call cline_refine2D_pool%set('extr_iter', pool_iter+1)
            call cline_refine2D_pool%set('extr_lim',   ITERLIM)
            ! Gaussian filter
            call cline_refine2D_pool%set('gauref',   'yes')
            call cline_refine2D_pool%set('gaufreq', lp)
        else
            call cline_refine2D_pool%set('trs', MINSHIFT)
            call cline_refine2D_pool%set('gauref', 'no')
            call cline_refine2D_pool%delete('extr_iter')
            call cline_refine2D_pool%delete('gaufreq')
            call cline_refine2D_pool%delete('lp')
        endif
    end subroutine update_pool_aln_params

    ! Deals with pool dimensions & resolution update
    subroutine update_pool_dims( params )
        use simple_procimgstk,    only: scale_imgfile
        use simple_classaverager, only: cavger_pad_carried_sums
        class(parameters), intent(inout) :: params
        type(scaled_dims) :: new_dims, prev_dims
        type(oris)        :: os
        type(class_frcs)  :: frcs
        type(string)      :: str, str_tmp_mrc
        real              :: scale_factor
        integer           :: ldim(3)
        ! resolution book-keeping
        resolutions(1:POOL_NPREV_RES-1) = resolutions(2:POOL_NPREV_RES)
        resolutions(POOL_NPREV_RES)     = current_resolution
        ! optional
        if( trim(params%dynreslim).ne.'yes' ) return
        prev_dims = pool_dims
        ! Auto-scaling?
        if( trim(params%autoscale) .ne. 'yes' ) return
        ! Hard limit reached?
        if( pool_dims%smpd < POOL_SMPD_HARD_LIMIT ) return
        ! Too early?
        if( pool_iter < 10 ) return
        ! Current resolution at Nyquist?
        if( abs(current_resolution-2.*pool_dims%smpd) > 0.01 ) return
        ! When POOL_NPREV_RES iterations are at Nyquist the pool resolution may be updated
        if( any(abs(resolutions-current_resolution) > 0.01 ) ) return
        ! determines new dimensions
        new_dims%box   = find_larger_magic_box(pool_dims%box+1)
        scale_factor   = real(new_dims%box) / real(pool_native_box)
        if( scale_factor > 0.99 ) return ! safety
        new_dims%smpd  = pool_native_smpd / scale_factor
        new_dims%boxpd = 2 * round2even(KBALPHA * real(new_dims%box/2)) ! logics from parameters
        ! New dimensions are accepted when new Nyquist is > 5/4 of original
        if( new_dims%smpd < 1.25*pool_native_smpd ) return
        ! Update global variables
        l_scaling   = .true.
        pool_dims   = new_dims
        call set_pool_mask()
        pool_lpstop = max(2.0*pool_dims%smpd, pool_user_lpstop)
        call cline_refine2D_pool%set('lpstop',     pool_lpstop)
        call cline_refine2D_pool%set('smpd_crop', pool_dims%smpd)
        call cline_refine2D_pool%set('box_crop',   pool_dims%box)
        write(logfhandle,'(A)')             '>>> UPDATING POOL DIMENSIONS '
        write(logfhandle,'(A,I5,A1,I5)')    '>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',pool_native_box,'/',pool_dims%box
        write(logfhandle,'(A,F5.2,A1,F5.2)')'>>> ORIGINAL/CROPPED PIXEL SIZE (Angs)  : ',pool_native_smpd,'/',pool_dims%smpd
        write(logfhandle,'(A,F5.1)')        '>>> POOL   HARD RESOLUTION LIMIT (Angs) : ',pool_lpstop
        ! upsample cavgs
        ldim = [pool_dims%box,pool_dims%box,1]
        str_tmp_mrc = TMP_STK_FNAME
        call scale_imgfile(refs_glob, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
        call simple_rename(str_tmp_mrc,refs_glob)
        str  = add2fbody(refs_glob, MRC_EXT,'_even')
        call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
        call simple_rename(str_tmp_mrc,str)
        str  = add2fbody(refs_glob, MRC_EXT,'_odd')
        call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
        call simple_rename(str_tmp_mrc,str)
        ! upsample the carried class sums, one set
        call cavger_pad_carried_sums(pool_dims%box, pool_dims%smpd)
        ! update cls2D field
        os = pool_proj%os_cls2D
        call pool_proj%os_out%kill
        call pool_proj%add_cavgs2os_out(refs_glob, pool_dims%smpd, 'cavg', clspath=.true.)
        pool_proj%os_cls2D = os
        call os%kill
        ! rescale frcs
        call frcs%read(string(FRCS_FILE))
        call frcs%pad(pool_dims%smpd, pool_dims%box)
        call frcs%write(string(FRCS_FILE))
        call frcs%kill
        call pool_proj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
    end subroutine update_pool_dims

    !> Points the pool project at the iteration's class averages and writes its class STAR file
    subroutine update_pool_for_gui( params )
        class(parameters), intent(in) :: params
        type(oris)        :: os_backup
        type(starproject) :: starproj
        os_backup = pool_proj%os_cls2D
        call pool_proj%add_cavgs2os_out(string(POOL_DIR)//refs_glob, pool_dims%smpd, 'cavg')
        pool_proj%os_cls2D = os_backup
        ! Write star file for iteration
        call starproj%export_cls2D(pool_proj, pool_iter)
        call pool_proj%os_cls2D%delete_entry('stk')
        call os_backup%kill
        call starproj%kill
    end subroutine update_pool_for_gui

    ! Flags pool availibility & updates the global name of references. An iteration whose job
    ! exited without finishing is submitted again once, from its project as made, after its log is
    ! kept aside and the files of its parts are removed; a second failure stops the pool
    ! (is_pool_failed).
    subroutine update_pool_status( params )
        class(parameters), intent(in) :: params
        type(string) :: failed_log
        if( .not. l_stream2D_active ) return
        if( l_pool_available .or. l_pool_failed ) return
        l_pool_available = file_exists(POOL_DIR//REFINE2D_FINISHED)
        if( l_pool_available )then
            if( pool_iter >= 1 ) refs_glob = CAVGS_ITER_FBODY//int2str_pad(pool_iter,3)//MRC_EXT
        else if( pool_job_failed() )then
            failed_log = POOL_DIR//POOL_LOGFILE//'_failed_iter'//int2str_pad(pool_iter,3)//'_attempt'//int2str(pool_nattempts)
            if( file_exists(POOL_DIR//POOL_LOGFILE) ) call simple_rename(string(POOL_DIR//POOL_LOGFILE), failed_log)
            if( pool_nattempts < 2 .and. file_exists(POOL_DIR//POOL_INPUT_PROJFILE) )then
                write(logfhandle,'(A,I6,A,A)') '>>> WARNING: POOL ITERATION ', pool_iter, ' FAILED; RETRYING IT ONCE. LOG: ',&
                    &failed_log%to_char()
                call qsys_cleanup(params)
                call simple_copy_file(string(POOL_DIR)//POOL_INPUT_PROJFILE, string(POOL_DIR)//POOL_PROJFILE)
                call submit_pool_iteration(params)
            else
                write(logfhandle,'(A,I6,A,A)') '>>> POOL ITERATION ', pool_iter, ' FAILED AGAIN; THE POOL STOPS. LOG: ',&
                    &failed_log%to_char()
                l_pool_failed = .true.
            endif
        endif
    end subroutine update_pool_status

    ! JPGS / GUI

    ! write jpeg of refs_glob    
    subroutine generate_pool_jpeg( params, filename)
        class(parameters), intent(in) :: params
        class(string), optional, intent(in) :: filename
        type(string)   :: jpeg_path, cwd
        type(image)    :: img, img_pad, img_jpeg
        type(stack_io) :: stkio_r
        integer        :: ldim_stk(3)
        integer        :: ncls_here, xtiles, ytiles, icls, ix, iy, ntiles
        call simple_getcwd(cwd)
        if(present(filename)) then
            jpeg_path = filename
        else
            jpeg_path = fname_new_ext(refs_glob, "jpeg") ! temporarily jpeg so compatible with old pool_stats. 
        end if
        if(.not. file_exists(refs_glob)) return
        if(file_exists(jpeg_path))       return
        if(allocated(pool_jpeg_map)) deallocate(pool_jpeg_map)
        if(allocated(pool_jpeg_pop)) deallocate(pool_jpeg_pop)
        if(allocated(pool_jpeg_res)) deallocate(pool_jpeg_res)
        allocate(pool_jpeg_map(0))
        allocate(pool_jpeg_pop(0))
        allocate(pool_jpeg_res(0))
        call find_ldim_nptcls(refs_glob, ldim_stk, ncls_here)
        if(ncls_here .ne. pool_proj%os_cls2D%get_noris()) THROW_HARD('ncls and n_noris mismatch')
        xtiles = floor(sqrt(real(ncls_glob)))
        ytiles = ceiling(real(ncls_glob) / real(xtiles))
        call img%new([ldim_stk(1), ldim_stk(1), 1], pool_dims%smpd)
        call img_pad%new([JPEG_DIM, JPEG_DIM, 1], pool_dims%smpd)
        call img_jpeg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], pool_dims%smpd)
        call stkio_r%open(refs_glob, pool_dims%smpd, 'read', bufsz=ncls_here)
        call stkio_r%read_whole
        ix = 1
        iy = 1
        ntiles = 0
        ! mask memoization
        call img%memoize_mask_coords
        do icls=1, ncls_here
            if(pool_proj%os_cls2D%get(icls,'state') < 0.5) cycle
            if(pool_proj%os_cls2D%get(icls,'pop')   < 0.5) cycle
            pool_jpeg_map = [pool_jpeg_map, icls]
            pool_jpeg_pop = [pool_jpeg_pop, nint(pool_proj%os_cls2D%get(icls,'pop'))]
            pool_jpeg_res = [pool_jpeg_res, pool_proj%os_cls2D%get(icls,'res')]
            call img%zero_and_unflag_ft
            call stkio_r%get_image(icls, img)
            call img%mask2D_softavg(pool_mskdiam / (2 * pool_dims%smpd))
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
        current_jpeg = cwd%to_char() // '/' // jpeg_path%to_char()
        current_jpeg_ntiles  = ntiles
        current_jpeg_ntilesx = xtiles
        current_jpeg_ntilesy = ytiles
        call img%kill()
        call img_pad%kill()
        call img_jpeg%kill()
    end subroutine generate_pool_jpeg

    ! When refine2D has not written the sprite sheet of the latest completed iteration, writes the
    ! iteration's class averages as images (generate_2D_jpeg, generate_pool_jpeg)
    subroutine generate_pool_stats( params )
        class(parameters), intent(in) :: params
        type(guistats) :: pool_stats
        type(string)   :: cwd
        integer        :: iter_loc = 0
        call simple_getcwd(cwd)
        if(file_exists(cwd//'/'//CLS2D_STARFBODY//'_iter'//int2str_pad(pool_iter,3)//STAR_EXT)) then
            iter_loc = pool_iter
        else if(file_exists(cwd//'/'//CLS2D_STARFBODY//'_iter'//int2str_pad(pool_iter - 1,3)//STAR_EXT)) then
            iter_loc = pool_iter - 1
        endif
        if(.not. file_exists(cwd//'/'//CAVGS_ITER_FBODY//int2str_pad(iter_loc, 3)//'.jpg')) then
            call pool_stats%init
            call pool_stats%generate_2D_jpeg('latest', '', pool_proj%os_cls2D, iter_loc, pool_dims%smpd)
            call pool_stats%kill
            last_complete_iter = iter_loc
            call generate_pool_jpeg(params)
        endif
    end subroutine generate_pool_stats

end module simple_stream_pool2D_utils
