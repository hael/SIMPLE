module simple_stream_pool2D_utils
include 'simple_lib.f08'
use simple_cmdline,       only: cmdline
use simple_euclid_sigma2, only: sigma2_star_from_iter
use simple_parameters,    only: params_glob
use simple_sp_project,    only: sp_project
use simple_stream_chunk,  only: stream_chunk
use simple_image,         only: image
use simple_stack_io,      only: stack_io
use simple_guistats,      only: guistats
use simple_class_frcs,    only: class_frcs
use simple_stream_cluster2D_utils
use simple_stream2D_state
use simple_rec_list
implicit none

! LIFECYCLE
public :: import_records_into_pool
! CALCULATORS
public :: init_pool_clustering
public :: analyze2D_pool
public :: iterate_pool
! GETTERS
public :: get_pool_assigned
public :: get_pool_cavgs_jpeg
public :: get_pool_cavgs_jpeg_ntiles
public :: get_pool_cavgs_jpeg_ntilesx
public :: get_pool_cavgs_jpeg_ntilesy
public :: get_pool_cavgs_mrc
public :: get_pool_iter
public :: get_pool_ptr
public :: get_pool_rejected
public :: is_pool_available
! SETTERS
public :: set_lpthres_type
public :: set_pool_dimensions
public :: set_pool_resolution_limits
! UPDATERS
public :: update_mskdiam
public :: update_pool
public :: update_pool_aln_params
public :: update_pool_dims
public :: update_pool_for_gui
public :: update_pool_status
! JPGS / GUI
public :: generate_pool_jpeg
public :: generate_pool_stats
private
#include "simple_local_flags.inc"

character(16)           :: last_iteration_time = ""         ! string with last iteration timestamp
integer                 :: current_jpeg_ntiles              ! number of used tiles in current JPEG
integer                 :: current_jpeg_ntilesx             ! number of tiles in x
integer                 :: current_jpeg_ntilesy             ! number of tiles in y
integer                 :: iterswitch2euclid = -1           ! iteration index where pool switches to euclid objfun
integer                 :: lim_ufrac_nptcls  = 0            ! threshold for fractional updates
integer                 :: ncls_max                         ! maximum allowed classes
integer                 :: ncls_rejected_glob               ! number of rejected classes
integer                 :: nptcls_glob                      ! total particles in pool
integer                 :: nptcls_rejected_glob             ! rejected particles in pool
integer(8), allocatable :: pool_proj_history_timestamps(:)  ! timestamps (e.g. from time8()) for history entries
logical,    allocatable :: pool_stacks_mask(:)              ! subset of stacks undergoing 2D analysis
real                    :: current_jpeg_scale               ! tile scaling factor
real                    :: current_resolution=999.          ! current estimated resolution
real                    :: resolutions(POOL_NPREV_RES)=999. ! pool resolution history (length POOL_NPREV_RES)
type(qsys_env)          :: pool_qenv                        ! qsys submission environment for pool
type(string)            :: current_jpeg                     ! filename of current pool JPEG (type(string))
! convergence
real                    :: conv_frac     = 0.0
real                    :: conv_mi_class = 0.0
real                    :: conv_score    = 0.0

contains

    ! LIFECYCLE

    ! Appends new data for processing
    subroutine import_records_into_pool( project_list )
        class(rec_list), intent(inout) :: project_list
        type(sp_project)     :: spproj
        type(string)         :: projname
        type(project_rec)    :: prec
        type(rec_iterator)   :: it
        logical, allocatable :: incl_mask(:)
        integer :: nptcls2import, nmics2import, imic, nrecords
        integer :: fromp, i, nmics_imported, nptcls_imported, iptcl, irec
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available  ) return
        nrecords = project_list%size()
        if( nrecords == 0         ) return
        incl_mask = project_list%get_included_flags()
        if( all(incl_mask)        ) return
        nmics_imported  = pool_proj%os_mic%get_noris()
        nptcls_imported = pool_proj%os_ptcl2D%get_noris()
        nmics2import    = count(.not.incl_mask)
        nptcls2import   = project_list%get_nptcls_tot(l_not_included=.true.)
        ! reallocations
        nmics_imported  = pool_proj%os_mic%get_noris()
        nptcls_imported = pool_proj%os_ptcl2D%get_noris()
        if( nmics_imported == 0 )then
            call pool_proj%os_mic%new(nmics2import,     is_ptcl=.false.)
            call pool_proj%os_stk%new(nmics2import,     is_ptcl=.false.)
            call pool_proj%os_ptcl2D%new(nptcls2import, is_ptcl=.true. )
            fromp = 1
        else
            call pool_proj%os_mic%reallocate(nmics_imported+nmics2import)
            call pool_proj%os_stk%reallocate(nmics_imported+nmics2import)
            call pool_proj%os_ptcl2D%reallocate(nptcls_imported+nptcls2import)
            fromp = pool_proj%os_stk%get_top(nmics_imported)+1
        endif
        imic     = nmics_imported
        projname = ''
        it       = project_list%begin()
        do irec = 1,nrecords
            ! retrieve one record from the list with the iterator
            call it%get(prec)
            if( prec%included )then
                ! move the iterator
                call it%next()
                cycle
            endif
            if( projname /= prec%projname )then
                call spproj%read_mic_stk_ptcl2D_segments(prec%projname)
                projname = prec%projname
            endif
            ! mic & stack
            imic = imic + 1
            call pool_proj%os_mic%transfer_ori(imic, spproj%os_mic, prec%micind)
            call pool_proj%os_stk%transfer_ori(imic, spproj%os_stk, prec%micind)
            call pool_proj%os_stk%set(imic, 'fromp', fromp)
            call pool_proj%os_stk%set(imic, 'top',   fromp+prec%nptcls-1)
            ! particles
            !$omp parallel do private(i,iptcl) proc_bind(close) default(shared)
            do i = 1,prec%nptcls
                iptcl = fromp + i - 1
                call pool_proj%os_ptcl2D%transfer_ori(iptcl, spproj%os_ptcl2D, i)
                call pool_proj%os_ptcl2D%set_stkind(iptcl, imic)
                call pool_proj%os_ptcl2D%set(iptcl, 'updatecnt', 0)
                call pool_proj%os_ptcl2D%set(iptcl, 'frac',      0.0)
                call pool_proj%os_ptcl2D%set(iptcl, 'eo',        merge(0, 1, is_even(iptcl)))
                call pool_proj%os_ptcl2D%set(iptcl, 'w',         1.0)
                call pool_proj%os_ptcl2D%set_class(iptcl, irnd_uni(params_glob%ncls))
            enddo
            !$omp end parallel do
            fromp = fromp + prec%nptcls
            ! flag inclusion
            prec%included = .true.
            ! replace the node
            call project_list%replace_iterator(it, prec)
            ! move the iterator
            call it%next()
        enddo
        call spproj%kill
    end subroutine import_records_into_pool

    ! CALCULATORS

    subroutine init_pool_clustering( cline, spproj, projfilegui, reference_generation )
        class(cmdline),    target, intent(inout) :: cline
        class(sp_project), intent(inout) :: spproj
        class(string),     intent(in)    :: projfilegui
        logical, optional, intent(in)    :: reference_generation
        type(string)          :: carg
        character(len=STDLEN) :: pool_nthr_env, pool_part_env, refgen_nthr_env, refgen_part_env
        integer               :: envlen
        call seed_rnd
        nullify(snapshot_json)
        ! reference generation: used for generating references from raw particles
        l_no_chunks = .false.
        if( present(reference_generation) ) l_no_chunks = reference_generation
        ! general parameters
        master_cline => cline
        call mskdiam2lplimits(params_glob%mskdiam, lpstart, lpstop, lpcen)
        l_wfilt            = .false.
        l_scaling          = .true.
        ncls_max           = params_glob%ncls
        ncls_glob          = params_glob%ncls
        ncls_rejected_glob = 0
        orig_projfile      = params_glob%projfile
        projfile4gui       = projfilegui
        l_update_sigmas    = params_glob%l_needs_sigma
        l_abinitio2D       = cline%defined('algorithm')
        if( l_abinitio2D ) l_abinitio2D = str_has_substr(params_glob%algorithm,'abinitio')
        params_glob%nparts_pool = params_glob%nparts ! backwards compatibility
        ! bookkeeping & directory structure
        numlen             = len(int2str(params_glob%nparts))
        refs_glob          = ''
        l_pool_available   = .true.
        pool_iter          = 0
        call simple_mkdir(POOL_DIR, verbose=.false.)
        call simple_mkdir(POOL_DIR//STDERROUT_DIR)
        call simple_mkdir(DIR_SNAPSHOT)
        if( l_update_sigmas ) call simple_mkdir(SIGMAS_DIR)
        pool_proj%projinfo = spproj%projinfo
        pool_proj%compenv  = spproj%compenv
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        ! update to computational parameters to pool, will be transferred to chunks upon init
        if( cline%defined('walltime') ) call pool_proj%compenv%set(1,'walltime', params_glob%walltime)
        ! commit to disk
        call pool_proj%write(string(POOL_DIR)//POOL_PROJFILE)
        ! reference generation
        if( l_no_chunks )then
            iterswitch2euclid = 0
            ncls_glob         = params_glob%ncls
        else
            iterswitch2euclid = -1 ! because objfun=euclid always
        endif
        ! Pool command line
        call cline_cluster2D_pool%set('prg',       'cluster2D_distr')
        call cline_cluster2D_pool%set('oritype',   'ptcl2D')
        call cline_cluster2D_pool%set('trs',       MINSHIFT)
        call cline_cluster2D_pool%set('projfile',  POOL_PROJFILE)
        call cline_cluster2D_pool%set('projname',  get_fbody(POOL_PROJFILE,'simple'))
        call cline_cluster2D_pool%set('sigma_est', params_glob%sigma_est)
        if( cline%defined('cls_init') )then
            call cline_cluster2D_pool%set('cls_init', params_glob%cls_init)
        else
            call cline_cluster2D_pool%set('cls_init', 'rand')
        endif
        if( cline%defined('center') )then
            carg = cline%get_carg('center')
            call cline_cluster2D_pool%set('center',carg)
            call carg%kill
        else
            call cline_cluster2D_pool%set('center','yes')
        endif
        if( .not.cline%defined('center_type') )then
            call cline_cluster2D_pool%set('center_type', 'seg')
        endif
        call cline_cluster2D_pool%set('extr_iter', 99999)
        call cline_cluster2D_pool%set('extr_lim',  MAX_EXTRLIM2D)
        call cline_cluster2D_pool%set('mkdir',     'no')
        call cline_cluster2D_pool%set('mskdiam',   params_glob%mskdiam)
        call cline_cluster2D_pool%set('async',     'yes') ! to enable hard termination
        call cline_cluster2D_pool%set('stream',    'yes')
        call cline_cluster2D_pool%set('nparts',    params_glob%nparts)
        call cline_cluster2D_pool%delete('autoscale')
        if( l_update_sigmas ) call cline_cluster2D_pool%set('cc_iters', 0)
        ! when the 2D analysis is started from raw particles
        if( l_no_chunks ) l_update_sigmas = .false.
        ! set # of ptcls beyond which fractional updates will be used
        lim_ufrac_nptcls = MAX_STREAM_NPTCLS
        if( master_cline%defined('nsample_max') ) lim_ufrac_nptcls = params_glob%nsample_max
        ! EV override
        params_glob%nthr2D = params_glob%nthr ! will be deprecated
        if( l_no_chunks )then
            call get_environment_variable(SIMPLE_STREAM_REFGEN_NTHR, refgen_nthr_env, envlen)
            if(envlen > 0) then
                call cline_cluster2D_pool%set('nthr', str2int(refgen_nthr_env))
            else
                call cline_cluster2D_pool%set('nthr', params_glob%nthr)
            end if
            call get_environment_variable(SIMPLE_STREAM_REFGEN_PARTITION, refgen_part_env, envlen)
            if(envlen > 0) then
                call pool_qenv%new(params_glob%nparts,exec_bin=string('simple_private_exec'),qsys_name=string('local'),&
                &qsys_partition=string(trim(refgen_part_env)))
            else
                call pool_qenv%new(params_glob%nparts,exec_bin=string('simple_private_exec'),qsys_name=string('local'))
            end if
        else
            call get_environment_variable(SIMPLE_STREAM_POOL_NTHR, pool_nthr_env, envlen)
            if(envlen > 0) then
                call cline_cluster2D_pool%set('nthr',   str2int(pool_nthr_env))
                call cline_cluster2D_pool%set('nthr2D', str2int(pool_nthr_env))
            else
                call cline_cluster2D_pool%set('nthr', params_glob%nthr)
            end if
            call get_environment_variable(SIMPLE_STREAM_POOL_PARTITION, pool_part_env, envlen)
            if(envlen > 0) then
                call pool_qenv%new(params_glob%nparts,exec_bin=string('simple_private_exec'),qsys_name=string('local'),&
                &qsys_partition=string(trim(pool_part_env)))
            else
                call pool_qenv%new(params_glob%nparts,exec_bin=string('simple_private_exec'),qsys_name=string('local'))
            end if
        end if
        ! objective function
        call cline_cluster2D_pool%set('objfun', 'euclid')
        call cline_cluster2D_pool%set('ml_reg', params_glob%ml_reg)
        call cline_cluster2D_pool%set('tau',    params_glob%tau)
        ! refinement
        select case(trim(params_glob%refine))
            case('snhc','snhc_smpl','prob','prob_smpl')
                if( (.not.l_abinitio2D) .and. str_has_substr(params_glob%refine,'prob') )then
                    THROW_HARD('REFINE=PROBXX only compatible with algorithm=abinitio2D')
                endif
                call cline_cluster2D_pool%set( 'refine', params_glob%refine)
            case DEFAULT
                THROW_HARD('UNSUPPORTED REFINE PARAMETER!')
        end select
        ! polar representation
        if( master_cline%defined('polar') ) call cline_cluster2D_pool%set('polar', params_glob%polar)
        ! Determines dimensions for downscaling
        call set_pool_dimensions
        ! updates command-lines with resolution limits
        call set_pool_resolution_limits
        ! module variables
        l_stream2D_active = .true.
    end subroutine init_pool_clustering

    ! Performs one iteration:
    ! updates to command-line, particles sampling, temporary project & execution
    subroutine analyze2D_pool
        use simple_euclid_sigma2, only: consolidate_sigma2_groups, average_sigma2_groups
        use simple_ran_tabu
        logical,         parameter    :: L_BENCH = .false.
        type(ran_tabu)                :: random_generator
        type(sp_project)              :: spproj, spproj_history
        type(sp_project), allocatable :: pool_proj_history_tmp(:)
        type(cmdline),    allocatable :: clines(:)
        integer,          allocatable :: min_update_cnts_per_stk(:), nptcls_per_stk(:), stk_order(:)
        integer,          allocatable :: prev_eo_pops(:,:), prev_eo_pops_thread(:,:)
        type(string),     allocatable :: sigma_fnames(:)
        type(string)                  :: stack_fname, ext, fbody, stkname
        real                          :: frac_update, smpd
        integer                       :: iptcl,i, nptcls_tot, nptcls_old, fromp, top, nstks_tot, jptcl, nptcls_th, npool
        integer                       :: eo, icls, nptcls_sel, istk, nptcls2update, nstks2update, jjptcl, nsample, ncls
        integer(timer_int_kind)       :: t_tot
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available )  return
        if( L_BENCH ) t_tot  = tic()
        nptcls_tot           = pool_proj%os_ptcl2D%get_noris()
        nptcls_glob          = nptcls_tot
        nptcls_rejected_glob = 0
        if( nptcls_tot == 0 ) return
        ! save pool project to history
        if(pool_iter .gt. 0) then
            call spproj_history%copy(pool_proj)
            call simple_copy_file(string(POOL_DIR)//FRCS_FILE, string(POOL_DIR)//swap_suffix(string(FRCS_FILE), "_iter"//int2str_pad(pool_iter, 3)//".bin", ".bin"))
            call spproj_history%get_cavgs_stk(stkname, ncls, smpd)
            call spproj_history%os_out%kill
            call spproj_history%add_cavgs2os_out(stkname, smpd, 'cavg')
            call spproj_history%add_frcs2os_out(string(POOL_DIR)//swap_suffix(string(FRCS_FILE),"_iter"//int2str_pad(pool_iter, 3)//".bin",".bin"),'frc2D')
            if( allocated(pool_proj_history) )then
                ! copy history
                npool = size(pool_proj_history)
                allocate(pool_proj_history_tmp(npool + 1))
                do i = 1, npool
                    pool_proj_history_tmp(i) =  pool_proj_history(i)
                end do
                ! add new element
                pool_proj_history_tmp(npool + 1) = spproj_history
                ! swap arrays
                do i = 1, npool
                    call pool_proj_history(i)%kill
                end do
                deallocate(pool_proj_history)
                allocate(pool_proj_history(npool + 1))
                do i = 1, npool
                    pool_proj_history(i) =  pool_proj_history_tmp(i)
                end do
                 do i = 1, npool + 1
                    call pool_proj_history_tmp(i)%kill
                end do
                deallocate(pool_proj_history_tmp)
            else
                allocate(pool_proj_history(1), source=spproj_history)
            endif
            if(pool_iter .gt. 5) call pool_proj_history(pool_iter - 5)%kill()
            call spproj_history%kill
        end if
        pool_iter = pool_iter + 1 ! Global iteration counter update
        call cline_cluster2D_pool%set('ncls',    ncls_glob)
        call cline_cluster2D_pool%set('startit', pool_iter)
        call cline_cluster2D_pool%set('maxits',  pool_iter)
        call cline_cluster2D_pool%set('frcs',    FRCS_FILE)
        call cline_cluster2D_pool%set('refs', refs_glob)
        if( pool_iter==1 )then
            if( cline_cluster2D_pool%defined('cls_init') )then
                ! references taken care of by cluster2D_distr
                call cline_cluster2D_pool%delete('frcs')
                call cline_cluster2D_pool%delete('refs')
            endif
        else
            call cline_cluster2D_pool%delete('cls_init')
        endif
        if( l_no_chunks )then
            ! for reference generation everything is defined here
            call cline_cluster2D_pool%set('extr_iter', CHUNK_EXTR_ITER+pool_iter-1)
            if( pool_iter == 1 )then
                call cline_cluster2D_pool%delete('frcs')
                call cline_cluster2D_pool%delete('refs')
            endif
            call cline_cluster2D_pool%set('center','no')
            if( pool_iter < 5 )then
                call cline_cluster2D_pool%delete('lpstart')
                call cline_cluster2D_pool%delete('lpstop')
                call cline_cluster2D_pool%set('lp', lpstart)
                call cline_cluster2D_pool%set('trs', 0)
            else
                call cline_cluster2D_pool%set('trs',    MINSHIFT)
                call cline_cluster2D_pool%set('center', params_glob%center)
            endif
            call cline_cluster2D_pool%set('needs_sigma', 'no')
            call cline_cluster2D_pool%set('objfun',      'cc')
            call cline_cluster2D_pool%set('ml_reg',      'no')
            call cline_cluster2D_pool%delete('cc_iters')
            if( params_glob%cc_objfun .eq. OBJFUN_EUCLID )then
                if( iterswitch2euclid == 0 )then
                    ! switch to objfun=euclid when good enough resolution
                    if( current_resolution < lpstart .and. pool_iter > CHUNK_CC_ITERS )then
                        iterswitch2euclid = pool_iter+1
                    endif
                else if( iterswitch2euclid == pool_iter )then
                    ! switch
                    call cline_cluster2D_pool%set('needs_sigma','yes')
                    call cline_cluster2D_pool%set('objfun',     'euclid')
                    call cline_cluster2D_pool%set('cc_iters',    0)
                    call cline_cluster2D_pool%set('sigma_est',  'global')
                    call cline_cluster2D_pool%set('ml_reg',     params_glob%ml_reg_pool)
                    call cline_cluster2D_pool%set('lpstart',    lpstart)
                    call cline_cluster2D_pool%set('lpstop',     params_glob%lpstop)
                    call cline_cluster2D_pool%delete('lp')
                    allocate(clines(2))
                    ! sigmas2 are calculated first thing
                    call clines(1)%set('prg',        'calc_pspec_distr')
                    call clines(1)%set('oritype',    'ptcl2D')
                    call clines(1)%set('projfile',   POOL_PROJFILE)
                    call clines(1)%set('nthr',       cline_cluster2D_pool%get_iarg('nthr'))
                    call clines(1)%set('which_iter', pool_iter)
                    call clines(1)%set('mkdir',      'yes')
                    call clines(1)%set('sigma_est',  'global')
                    call clines(1)%set('nparts',     params_glob%nparts_pool)
                    clines(2) = cline_cluster2D_pool
                else
                    ! after switch
                    call cline_cluster2D_pool%set('needs_sigma','yes')
                    call cline_cluster2D_pool%set('objfun',     'euclid')
                    call cline_cluster2D_pool%set('cc_iters',   0)
                    call cline_cluster2D_pool%set('sigma_est',  'global')
                    call cline_cluster2D_pool%set('ml_reg',      params_glob%ml_reg_pool)
                    call cline_cluster2D_pool%set('lpstart',     lpstart)
                    call cline_cluster2D_pool%set('lpstop',      params_glob%lpstop)
                    do i = 1,params_glob%nparts_pool
                        call del_file(SIGMA2_FBODY//int2str_pad(i,numlen)//'.dat')
                    enddo
                endif
            else
                ! objfun = cc
                if( current_resolution < lpstart .and. pool_iter > CHUNK_CC_ITERS )then
                    call cline_cluster2D_pool%set('lpstart',lpstart)
                    call cline_cluster2D_pool%set('lpstop', params_glob%lpstop)
                    call cline_cluster2D_pool%delete('lp')
                endif
            endif
        endif
        ! project metadata update 
        spproj%projinfo = pool_proj%projinfo
        spproj%compenv  = pool_proj%compenv
        call spproj%projinfo%delete_entry('projname')
        call spproj%projinfo%delete_entry('projfile')
        call spproj%update_projinfo( cline_cluster2D_pool )
        ! Sampling of stacks that will be used for this iteration
        ! counting number of stacks & particles (old and newer)
        nstks_tot  = pool_proj%os_stk%get_noris()
        allocate(nptcls_per_stk(nstks_tot), min_update_cnts_per_stk(nstks_tot), source=0)
        nptcls_old = 0 ! Particles updated more than STREAM_SRCHLIM-1
        !$omp parallel do schedule(static) proc_bind(close) private(istk,fromp,top,iptcl)&
        !$omp default(shared) reduction(+:nptcls_old)
        do istk = 1,nstks_tot
            fromp = pool_proj%os_stk%get_fromp(istk)
            top   = pool_proj%os_stk%get_top(istk)
            min_update_cnts_per_stk(istk) = huge(istk)
            do iptcl = fromp,top
                if( pool_proj%os_ptcl2D%get_state(iptcl) > 0 )then
                    nptcls_per_stk(istk)          = nptcls_per_stk(istk) + 1 ! # ptcls with state=1
                    min_update_cnts_per_stk(istk) = min(min_update_cnts_per_stk(istk), pool_proj%os_ptcl2D%get_updatecnt(iptcl))
                endif
            enddo
            if( min_update_cnts_per_stk(istk) >= STREAM_SRCHLIM )then
                nptcls_old = nptcls_old + nptcls_per_stk(istk)
            endif
        enddo
        !$omp end parallel do
        nptcls_rejected_glob = nptcls_glob - sum(nptcls_per_stk)
        ! update info for gui
        call spproj%projinfo%set(1,'nptcls_tot',     nptcls_glob)
        call spproj%projinfo%set(1,'nptcls_rejected',nptcls_rejected_glob)
        ! flagging stacks to be skipped
        if( allocated(pool_stacks_mask) ) deallocate(pool_stacks_mask)
        allocate(pool_stacks_mask(nstks_tot), source=.false.)
        if( nptcls_old > params_glob%ncls_start*params_glob%nptcls_per_cls )then
            allocate(stk_order(nstks_tot), source=(/(istk,istk=1,nstks_tot)/))
            random_generator = ran_tabu(nstks_tot)
            call random_generator%shuffle(stk_order)
            nptcls2update = 0 ! # of ptcls including state=0 within selected stacks
            nptcls_sel    = 0 ! # of ptcls excluding state=0 within selected stacks
            do i = 1,nstks_tot
                istk = stk_order(i)
                if( (min_update_cnts_per_stk(istk) > STREAM_SRCHLIM) .and. (nptcls_sel > lim_ufrac_nptcls) ) cycle
                nptcls_sel    = nptcls_sel + nptcls_per_stk(istk)
                nptcls2update = nptcls2update + pool_proj%os_stk%get_int(istk, 'nptcls')
                pool_stacks_mask(istk) = .true.
            enddo
            call random_generator%kill
        else
            nptcls2update    = nptcls_tot
            nptcls_sel       = sum(nptcls_per_stk)
            pool_stacks_mask = nptcls_per_stk > 0
        endif
        nstks2update = count(pool_stacks_mask)
        ! transfer stacks and particles
        call spproj%os_stk%new(nstks2update, is_ptcl=.false.)
        call spproj%os_ptcl2D%new(nptcls2update, is_ptcl=.true.)
        allocate(prev_eo_pops(ncls_glob,2),prev_eo_pops_thread(ncls_glob,2),source=0)
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
            else
                ! keeps track of skipped particles
                prev_eo_pops_thread = 0
                !$omp parallel do private(iptcl,icls,eo) proc_bind(close) default(shared)&
                !$omp reduction(+:prev_eo_pops_thread)
                do iptcl = fromp,top
                    if( pool_proj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    icls = pool_proj%os_ptcl2D%get_class(iptcl)
                    eo   = pool_proj%os_ptcl2D%get_eo(iptcl) + 1
                    prev_eo_pops_thread(icls,eo) = prev_eo_pops_thread(icls,eo) + 1
                enddo
                !$omp end parallel do
                prev_eo_pops = prev_eo_pops + prev_eo_pops_thread
            endif
        enddo
        call spproj%os_ptcl3D%new(nptcls2update, is_ptcl=.true.)
        spproj%os_cls2D = pool_proj%os_cls2D
        ! consolidate sigmas doc
        if( l_update_sigmas )then
            if( trim(params_glob%sigma_est).eq.'group' )then
                allocate(sigma_fnames(nstks2update))
                do istk = 1,nstks2update
                    call spproj%os_stk%getter(istk,'stk',stack_fname)
                    stack_fname = basename(stack_fname)
                    ext         = fname2ext(stack_fname)
                    fbody       = get_fbody(stack_fname, ext)
                    sigma_fnames(istk) = SIGMAS_DIR//'/'//fbody%to_char()//STAR_EXT
                enddo
                call consolidate_sigma2_groups(sigma2_star_from_iter(pool_iter), sigma_fnames)
                deallocate(sigma_fnames)
            else
                ! sigma_est=global & first iteration
                if( pool_iter==1 )then
                    allocate(sigma_fnames(glob_chunk_id))
                    do i = 1,glob_chunk_id
                        sigma_fnames(i) = SIGMAS_DIR//'/chunk_'//int2str(i)//STAR_EXT
                    enddo
                    call average_sigma2_groups(sigma2_star_from_iter(pool_iter), sigma_fnames)
                    deallocate(sigma_fnames)
                endif
            endif
            do i = 1,params_glob%nparts_pool
                call del_file(SIGMA2_FBODY//int2str_pad(i,numlen)//'.dat')
            enddo
        endif
        ! update command line with fractional update parameters
        call cline_cluster2D_pool%delete('update_frac')
        call cline_cluster2D_pool%delete('maxpop')
        call cline_cluster2D_pool%delete('nsample')
        frac_update = 1.0
        if( trim(params_glob%autosample).eq.'yes' )then
            ! momentum & limits to # of particles per class
            nsample = MAXPOP_CLS
            if( master_cline%defined('nsample') ) nsample = params_glob%nsample
            nptcls_th = lim_ufrac_nptcls
            if( nptcls_sel > nptcls_th )then
                frac_update = real(nptcls_th) / real(nptcls_sel)
                call cline_cluster2D_pool%set('maxpop',    nsample)
                call cline_cluster2D_pool%set('sigma_est', 'global')
            endif
        else
            ! momentum only
            if( nptcls_sel > lim_ufrac_nptcls )then
                if( (sum(prev_eo_pops) > 0) .and. (nptcls_old > 0))then
                    frac_update = real(nptcls_old-sum(prev_eo_pops)) / real(nptcls_old)
                endif
            endif
        endif
        ! user override
        if( frac_update < 0.99999 )then
            if( master_cline%defined('update_frac') ) frac_update = params_glob%update_frac
            call cline_cluster2D_pool%set('update_frac', frac_update)
            call cline_cluster2D_pool%set('center',      'no')
            do icls = 1,ncls_glob
                call spproj%os_cls2D%set(icls,'prev_pop_even',prev_eo_pops(icls,1))
                call spproj%os_cls2D%set(icls,'prev_pop_odd', prev_eo_pops(icls,2))
            enddo
        endif
        ! write project
        call spproj%write(string(POOL_DIR)//POOL_PROJFILE)
        call spproj%kill
        ! pool stats
        call generate_pool_stats
        ! execution
        if( l_no_chunks .and. pool_iter == iterswitch2euclid )then
            write(logfhandle,'(A)')'>>> SWITCHING TO OBJFUN=EUCLID'
            call pool_qenv%exec_simple_prgs_in_queue_async(clines, string(POOL_DISTR_EXEC_FNAME), string(POOL_LOGFILE))
            call clines(:)%kill
            deallocate(clines)
        else
            call pool_qenv%exec_simple_prg_in_queue_async(cline_cluster2D_pool, string(POOL_DISTR_EXEC_FNAME), string(POOL_LOGFILE))
        endif
        l_pool_available = .false.
        write(logfhandle,'(A,I6,A,I8,A3,I8,A)')'>>> POOL         INITIATED ITERATION ',pool_iter,' WITH ',nptcls_sel,&
        &' / ', sum(nptcls_per_stk),' PARTICLES'
        if( L_BENCH ) print *,'timer analyze2D_pool tot : ',toc(t_tot)
        call tidy_2Dstream_iter
    end subroutine analyze2D_pool

    ! Performs one iteration:
    ! updates to command-line, particles sampling, temporary project & execution
    subroutine iterate_pool
        logical, parameter            :: L_BENCH = .false.
        type(sp_project)              :: spproj, spproj_history
        type(sp_project), allocatable :: pool_proj_history_tmp(:)
        integer(timer_int_kind)       :: t_tot
        integer,          allocatable :: nptcls_per_stk(:), prev_eo_pops(:,:), prev_eo_pops_thread(:,:), clspops(:)
        type(string) :: stkname
        real         :: frac_update, smpd
        integer      :: iptcl,i, nptcls_tot, nptcls_old, fromp, top, nstks_tot, jptcl, npool
        integer      :: eo, icls, nptcls_sel, istk, nptcls2update, nstks2update, jjptcl, ncls
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available  ) return
        if( l_no_chunks ) THROW_HARD('Designed for pre-clustered/matched particles!')
        if( L_BENCH ) t_tot  = tic()
        nptcls_tot           = pool_proj%os_ptcl2D%get_noris()
        nptcls_glob          = nptcls_tot
        nptcls_rejected_glob = 0
        if( nptcls_tot == 0 ) return
        ! save pool project to history
        if(pool_iter .gt. 0) then
            call spproj_history%copy(pool_proj)
            call simple_copy_file(string(POOL_DIR)//FRCS_FILE, string(POOL_DIR)//swap_suffix(FRCS_FILE,"_iter"//int2str_pad(pool_iter, 3)//".bin",".bin"))
            call spproj_history%get_cavgs_stk(stkname, ncls, smpd)
            call spproj_history%os_out%kill
            call spproj_history%add_cavgs2os_out(stkname, smpd, 'cavg')
            call spproj_history%add_frcs2os_out(string(POOL_DIR)//swap_suffix(FRCS_FILE,"_iter"//int2str_pad(pool_iter, 3)//".bin",".bin"),'frc2D')
            if( allocated(pool_proj_history) )then
                ! copy history
                npool = size(pool_proj_history)
                allocate(pool_proj_history_tmp(npool + 1))
                do i = 1, npool
                    pool_proj_history_tmp(i) =  pool_proj_history(i)
                end do
                ! add new element
                pool_proj_history_tmp(npool + 1) = spproj_history
                ! swap arrays
                do i = 1, npool
                    call pool_proj_history(i)%kill
                end do
                deallocate(pool_proj_history)
                allocate(pool_proj_history(npool + 1))
                do i = 1, npool
                    pool_proj_history(i) =  pool_proj_history_tmp(i)
                end do
                 do i = 1, npool + 1
                    call pool_proj_history_tmp(i)%kill
                end do
                deallocate(pool_proj_history_tmp)
            else
                allocate(pool_proj_history(1))
                pool_proj_history(1) = spproj_history
            endif
            if( allocated(pool_proj_history_timestamps) )then
                pool_proj_history_timestamps = [pool_proj_history_timestamps(:), time8()]
            else
                allocate(pool_proj_history_timestamps(1), source=time8())
            endif
            do i=1, size(pool_proj_history_timestamps)
                ! remove history older than 5 minutes
                if(pool_proj_history_timestamps(i) > 0 .and. pool_proj_history_timestamps(i) < time8() - 300) then
                    call pool_proj_history(i)%kill()
                    pool_proj_history_timestamps(i) = 0
                end if
            end do
            call spproj_history%kill
        end if
        pool_iter = pool_iter + 1 ! Global iteration counter update
        call cline_cluster2D_pool%set('ncls',    ncls_glob)
        call cline_cluster2D_pool%set('startit', pool_iter)
        call cline_cluster2D_pool%set('maxits',  pool_iter)
        call cline_cluster2D_pool%set('frcs',    FRCS_FILE)
        call cline_cluster2D_pool%set('refs', refs_glob)
        if( pool_iter==1 )then
            if( cline_cluster2D_pool%defined('cls_init') )then
                ! references taken care of by cluster2D_distr
                call cline_cluster2D_pool%delete('frcs')
                call cline_cluster2D_pool%delete('refs')
            endif
        else
            call cline_cluster2D_pool%delete('cls_init')
        endif
        ! Project metadata update
        spproj%projinfo = pool_proj%projinfo
        spproj%compenv  = pool_proj%compenv
        call spproj%projinfo%delete_entry('projname')
        call spproj%projinfo%delete_entry('projfile')
        call spproj%update_projinfo( cline_cluster2D_pool )
        ! Sampling of stacks that will be used for this iteration
        ! counting number of stacks & selected particles
        nstks_tot  = pool_proj%os_stk%get_noris()
        allocate(nptcls_per_stk(nstks_tot), source=0)
        nptcls_old = 0 ! Total # of particles with state=1
        !$omp parallel do schedule(static) proc_bind(close) private(istk,fromp,top,iptcl)&
        !$omp default(shared) reduction(+:nptcls_old)
        do istk = 1,nstks_tot
            fromp = pool_proj%os_stk%get_fromp(istk)
            top   = pool_proj%os_stk%get_top(istk)
            do iptcl = fromp,top
                if( pool_proj%os_ptcl2D%get_state(iptcl) > 0 )then
                    nptcls_per_stk(istk)  = nptcls_per_stk(istk) + 1 ! # ptcls with state=1
                endif
            enddo
            nptcls_old = nptcls_old + nptcls_per_stk(istk)
        enddo
        !$omp end parallel do
        nptcls_rejected_glob = nptcls_glob - sum(nptcls_per_stk)
        ! Update info for gui
        call spproj%projinfo%set(1,'nptcls_tot',     nptcls_glob)
        call spproj%projinfo%set(1,'nptcls_rejected',nptcls_rejected_glob)
        ! Uniformly sample stacks
        call uniform_stack_sampling
        ! call biased_stack_sampling
        nstks2update = count(pool_stacks_mask)
        ! Transfer stacks and particles
        call spproj%os_stk%new(nstks2update, is_ptcl=.false.)
        call spproj%os_ptcl2D%new(nptcls2update, is_ptcl=.true.)
        allocate(prev_eo_pops(ncls_glob,2),prev_eo_pops_thread(ncls_glob,2),source=0)
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
            else
                ! keeps track of skipped particles
                prev_eo_pops_thread = 0
                !$omp parallel do private(iptcl,icls,eo) proc_bind(close) default(shared)&
                !$omp reduction(+:prev_eo_pops_thread)
                do iptcl = fromp,top
                    if( pool_proj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    if( pool_proj%os_ptcl2D%get_updatecnt(iptcl) == 0 ) cycle
                    icls = pool_proj%os_ptcl2D%get_class(iptcl)
                    eo   = pool_proj%os_ptcl2D%get_eo(iptcl) + 1
                    prev_eo_pops_thread(icls,eo) = prev_eo_pops_thread(icls,eo) + 1
                enddo
                !$omp end parallel do
                prev_eo_pops = prev_eo_pops + prev_eo_pops_thread
            endif
        enddo
        call spproj%os_ptcl3D%new(nptcls2update, is_ptcl=.true.)
        spproj%os_cls2D = pool_proj%os_cls2D
        ! making sure the new particles are asigned a populated class
        if( pool_iter >= 2 )then
            clspops = spproj%os_cls2D%get_all_asint('pop')
            !$omp parallel do private(iptcl,icls) proc_bind(close) default(shared) schedule(static)
            do iptcl = 1,nptcls2update
                if( spproj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                if( spproj%os_ptcl2D%get_updatecnt(iptcl) == 0 )then
                    icls = irnd_uni(ncls_glob)
                    do while( clspops(icls) == 0 )
                        icls = irnd_uni(ncls_glob)
                    enddo
                    call spproj%os_ptcl2D%set_class(iptcl, icls)
                endif
            enddo
            !$omp end parallel do
        endif
        ! Consolidate sigmas doc
        call consolidate_sigmas(spproj, nstks2update)
        ! update command line with fractional update parameters
        call cline_cluster2D_pool%delete('update_frac')
        frac_update = 1.0
        if( nptcls_sel > lim_ufrac_nptcls )then
            if( (sum(prev_eo_pops) > 0) .and. (nptcls_old > 0))then
                frac_update = real(nptcls_old-sum(prev_eo_pops)) / real(nptcls_old)
            endif
        endif
        ! User override
        if( frac_update < 0.99999 )then
            if( master_cline%defined('update_frac') ) frac_update = params_glob%update_frac
            call cline_cluster2D_pool%set('update_frac', frac_update)
            call cline_cluster2D_pool%set('center',      'no')
            do icls = 1,ncls_glob
                call spproj%os_cls2D%set(icls,'prev_pop_even',prev_eo_pops(icls,1))
                call spproj%os_cls2D%set(icls,'prev_pop_odd', prev_eo_pops(icls,2))
            enddo
        endif
        ! write project
        call spproj%write(string(POOL_DIR)//POOL_PROJFILE)
        call spproj%kill
        ! pool stats
        call generate_pool_stats
        ! execution
        call pool_qenv%exec_simple_prg_in_queue_async(cline_cluster2D_pool, string(POOL_DISTR_EXEC_FNAME), string(POOL_LOGFILE))
        l_pool_available = .false.
        write(logfhandle,'(A,I6,A,I8,A3,I8,A)')'>>> POOL         INITIATED ITERATION ',pool_iter,' WITH ',nptcls_sel,&
        &' / ', sum(nptcls_per_stk),' PARTICLES'
        if( L_BENCH ) print *,'timer analyze2D_pool tot : ',toc(t_tot)
        ! cleanup
        if( allocated(nptcls_per_stk) )      deallocate(nptcls_per_stk)
        if( allocated(prev_eo_pops) )        deallocate(prev_eo_pops)
        if( allocated(prev_eo_pops_thread) ) deallocate(prev_eo_pops_thread)
        if( allocated(clspops) )             deallocate(clspops)
        call tidy_2Dstream_iter

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

        ! sampling biased towards older stacks
        subroutine biased_stack_sampling
            real, parameter :: B = 0.1
            real    :: vals(nstks_tot), diff
            integer :: counts(nstks_tot), inds(nstks_tot)
            integer :: fromp,top,iptcl,i,j,minc,maxc
            !$omp parallel do private(i,iptcl,fromp,top) proc_bind(close) default(shared)
            do i = 1,nstks_tot
                counts(i) = huge(i)
                fromp = pool_proj%os_stk%get_fromp(i)
                top   = pool_proj%os_stk%get_top(i)
                do iptcl = fromp,top
                    if( pool_proj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    counts(i) = min(counts(i), pool_proj%os_ptcl2D%get_updatecnt(iptcl))
                enddo

                call pool_proj%os_stk%set(i, 'cnt', counts(i))

            enddo
            !$omp end parallel do
            minc = minval(counts)
            maxc = maxval(counts)
            if( maxc == minc )then
                ! all stacks have been sampled uniformly, use uniform sampling
                call uniform_stack_sampling
            else
                if( allocated(pool_stacks_mask) ) deallocate(pool_stacks_mask)
                allocate(pool_stacks_mask(nstks_tot), source=.false.)
                ! generate biased order
                diff = real(1+maxc-minc)
                do i = 1,nstks_tot
                    inds(i) = i
                    vals(i) = real(1+counts(i)-minc) / diff
                    vals(i) = vals(i)**B - ran3()
                enddo
                ! draw
                call hpsort(vals, inds)
                nptcls2update = 0 ! # of ptcls including state=0 within selected stacks
                nptcls_sel    = 0 ! # of ptcls excluding state=0 within selected stacks
                do i = 1,nstks_tot
                    j = inds(i)
                    if( nptcls_sel > lim_ufrac_nptcls ) cycle
                    nptcls_sel    = nptcls_sel    + nptcls_per_stk(j)
                    nptcls2update = nptcls2update + pool_proj%os_stk%get_int(j, 'nptcls')
                    pool_stacks_mask(j) = .true.
                enddo
            endif
        end subroutine biased_stack_sampling

    end subroutine iterate_pool

    ! GETTERS

    ! returns number currently assigned particles
    integer function get_pool_assigned()
        get_pool_assigned = nptcls_glob - nptcls_rejected_glob
    end function get_pool_assigned

    type(string) function get_pool_cavgs_jpeg()
        get_pool_cavgs_jpeg = current_jpeg
    end function get_pool_cavgs_jpeg

    integer function get_pool_cavgs_jpeg_ntiles()
        get_pool_cavgs_jpeg_ntiles = current_jpeg_ntiles
    end function get_pool_cavgs_jpeg_ntiles

    integer function get_pool_cavgs_jpeg_ntilesx()
        get_pool_cavgs_jpeg_ntilesx = current_jpeg_ntilesx
    end function get_pool_cavgs_jpeg_ntilesx

    integer function get_pool_cavgs_jpeg_ntilesy()
        get_pool_cavgs_jpeg_ntilesy = current_jpeg_ntilesy
    end function get_pool_cavgs_jpeg_ntilesy

    type(string) function get_pool_cavgs_mrc()
        get_pool_cavgs_mrc = refs_glob
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

    ! whether the pool available for another iteration
    logical function is_pool_available()
        is_pool_available = l_pool_available
    end function is_pool_available

    ! SETTERS

    subroutine set_lpthres_type(type)
        character(len=*), intent(in) :: type
        lpthres_type = type
    end subroutine set_lpthres_type

    subroutine set_pool_dimensions
        call setup_downscaling
        pool_dims%smpd  = params_glob%smpd_crop
        pool_dims%box   = params_glob%box_crop
        pool_dims%boxpd = 2*round2even(params_glob%alpha*real(params_glob%box_crop/2)) ! logics from parameters
        pool_dims%msk   = params_glob%msk_crop
        ! chunk & pool have the same dimensions to start with (used for import)
        chunk_dims = pool_dims
        ! Scaling-related command lines update
        call cline_cluster2D_pool%set('smpd_crop',  pool_dims%smpd)
        call cline_cluster2D_pool%set('box_crop',   pool_dims%box)
        call cline_cluster2D_pool%set('msk_crop',   pool_dims%msk)
        call cline_cluster2D_pool%set('box',        params_glob%box)
        call cline_cluster2D_pool%set('smpd',       params_glob%smpd)
    end subroutine set_pool_dimensions

    ! private routine for pool resolution-related updates to command-lines
    subroutine set_pool_resolution_limits
        lpstart = max(lpstart, 2.0*params_glob%smpd_crop)
        if( l_no_chunks )then
            params_glob%lpstop = lpstop
        else
            if( master_cline%defined('lpstop') )then
                params_glob%lpstop = max(2.0*params_glob%smpd_crop,params_glob%lpstop)
            else
                params_glob%lpstop = 2.0*params_glob%smpd_crop
            endif
        endif
        call cline_cluster2D_pool%set('lpstart',  lpstart)
        call cline_cluster2D_pool%set('lpstop',   params_glob%lpstop)
        if( .not.master_cline%defined('cenlp') )then
            call cline_cluster2D_pool%set( 'cenlp', lpcen)
        else
            call cline_cluster2D_pool%set( 'cenlp', params_glob%cenlp)
        endif
        write(logfhandle,'(A,F5.1)') '>>> STARTING LOW-PASS LIMIT  (IN A): ', lpstart
        write(logfhandle,'(A,F5.1)') '>>> HARD RESOLUTION LIMIT    (IN A): ', params_glob%lpstop
        write(logfhandle,'(A,F5.1)') '>>> CENTERING LOW-PASS LIMIT (IN A): ', lpcen
    end subroutine set_pool_resolution_limits

    ! UPDATERS

    subroutine update_mskdiam(new_mskdiam)
        integer, intent(in) :: new_mskdiam
        write(*,'(A,I4,A)')'>>> UPDATED MASK DIAMETER TO',new_mskdiam ,'Å'
        params_glob%mskdiam = real(new_mskdiam)
        call cline_cluster2D_pool%set('mskdiam',   params_glob%mskdiam)
    end subroutine update_mskdiam

    ! Reports alignment info from completed iteration of subset
    ! of particles back to the pool
    subroutine update_pool
        use simple_euclid_sigma2, only: split_sigma2_into_groups
        type(string), allocatable :: sigma_fnames(:)
        integer,      allocatable :: pops(:)
        type(sp_project) :: spproj
        type(oris)       :: os
        type(class_frcs) :: frcs
        type(string)     :: stack_fname, ext, fbody, fname
        integer          :: i, it, jptcl, iptcl, istk, nstks
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available  ) return
        call del_file(POOL_DIR//CLUSTER2D_FINISHED)
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
            ! updates sigmas
            if( l_update_sigmas )then
                if( trim(params_glob%sigma_est).eq.'group' )then
                    ! propagate sigma2 changes back to the micrograph/stack document
                    nstks = spproj%os_stk%get_noris()
                    allocate(sigma_fnames(nstks))
                    do istk = 1,nstks
                        call spproj%os_stk%getter(istk,'stk',stack_fname)
                        stack_fname = basename(stack_fname)
                        ext         = fname2ext(stack_fname)
                        fbody       = get_fbody(stack_fname, ext)
                        sigma_fnames(istk) = SIGMAS_DIR//'/'//fbody%to_char()//STAR_EXT
                    enddo
                    call split_sigma2_into_groups(sigma2_star_from_iter(pool_iter+1), sigma_fnames)
                    deallocate(sigma_fnames)
                else
                    ! sigma_est=global, nothing to do
                endif
            endif
            call frcs%read(string(POOL_DIR)//FRCS_FILE)
            current_resolution = frcs%estimate_lp_for_align()
            write(logfhandle,'(A,F5.1)')'>>> CURRENT POOL RESOLUTION: ',current_resolution
            call frcs%kill
            ! deal with dimensions/resolution update
            call update_pool_dims
            ! for gui
            call update_pool_for_gui
        endif
        call spproj%kill
    end subroutine update_pool

    ! This controls the evolution of the pool alignement parameters:
    ! lp, ICM, trs, extr_iter
    subroutine update_pool_aln_params
        real,    parameter :: ICM_LAMBDA = 2.0
        integer, parameter :: ITERLIM    = 20
        integer, parameter :: ITERSHIFT  = 5
        real :: lp, lambda, gamma
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available  ) return
        if( pool_iter < ITERLIM )then
            gamma = min(1., max(0., real(ITERLIM-pool_iter)/real(ITERLIM)))
            ! offset
            if( pool_iter < ITERSHIFT )then
                call cline_cluster2D_pool%set('trs', 0.)
            else
                call cline_cluster2D_pool%set('trs', MINSHIFT)
            endif
            ! resolution limit
            lp = lpstop + (lpstart-lpstop) * gamma
            call cline_cluster2D_pool%set('lp', lp)
            ! Extremal iteration
            call cline_cluster2D_pool%set('extr_iter', pool_iter+1)
            call cline_cluster2D_pool%set('extr_lim',  ITERLIM)
            ! ICM
            lambda = ICM_LAMBDA * gamma
            call cline_cluster2D_pool%set('icm',    'yes')
            call cline_cluster2D_pool%set('lambda', lambda)
        else
            call cline_cluster2D_pool%set('trs', MINSHIFT)
            call cline_cluster2D_pool%set('icm', 'no')
            call cline_cluster2D_pool%delete('extr_iter')
            call cline_cluster2D_pool%delete('lambda')
            call cline_cluster2D_pool%delete('lp')
        endif
    end subroutine update_pool_aln_params

    ! Deals with pool dimensions & resolution update
    subroutine update_pool_dims
        use simple_procimgstk, only: pad_imgfile, scale_imgfile
        type(scaled_dims) :: new_dims, prev_dims
        type(oris)        :: os
        type(class_frcs)  :: frcs
        type(string)      :: str, str_tmp_mrc
        real              :: scale_factor
        integer           :: ldim(3), p
        ! resolution book-keeping
        resolutions(1:POOL_NPREV_RES-1) = resolutions(2:POOL_NPREV_RES)
        resolutions(POOL_NPREV_RES)     = current_resolution
        if( l_no_chunks ) return
        ! optional
        if( trim(params_glob%dynreslim).ne.'yes' ) return
        prev_dims = pool_dims
        ! Auto-scaling?
        if( trim(params_glob%autoscale) .ne. 'yes' ) return
        ! Hard limit reached?
        if( pool_dims%smpd < POOL_SMPD_HARD_LIMIT ) return
        ! Too early?
        if( pool_iter < 10 ) return
        if( ncls_glob < ncls_max ) return
        ! Current resolution at Nyquist?
        if( abs(current_resolution-2.*pool_dims%smpd) > 0.01 ) return
        ! When POOL_NPREV_RES iterations are at Nyquist the pool resolution may be updated
        if( any(abs(resolutions-current_resolution) > 0.01 ) ) return
        ! determines new dimensions
        new_dims%box   = find_larger_magic_box(pool_dims%box+1)
        scale_factor   = real(new_dims%box) / real(params_glob%box)
        if( scale_factor > 0.99 ) return ! safety
        new_dims%smpd  = params_glob%smpd / scale_factor
        new_dims%boxpd = 2 * round2even(params_glob%alpha * real(new_dims%box/2)) ! logics from parameters
        new_dims%msk   = round2even(params_glob%mskdiam / new_dims%smpd / 2.)
        ! New dimensions are accepted when new Nyquist is > 5/4 of original
        if( new_dims%smpd < 1.25*params_glob%smpd ) return
        ! Update global variables
        l_scaling = .true.
        pool_dims = new_dims
        if( master_cline%defined('lpstop') )then
            params_glob%lpstop = max(2.0*pool_dims%smpd, master_cline%get_rarg('lpstop'))
        else
            params_glob%lpstop = 2.0*pool_dims%smpd
        endif
        call cline_cluster2D_pool%set('lpstop',    params_glob%lpstop)
        call cline_cluster2D_pool%set('smpd_crop', pool_dims%smpd)
        call cline_cluster2D_pool%set('box_crop',  pool_dims%box)
        call cline_cluster2D_pool%set('msk_crop',  pool_dims%msk)
        write(logfhandle,'(A)')             '>>> UPDATING POOL DIMENSIONS '
        write(logfhandle,'(A,I5,A1,I5)')    '>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params_glob%box,'/',pool_dims%box
        write(logfhandle,'(A,F5.2,A1,F5.2)')'>>> ORIGINAL/CROPPED PIXEL SIZE (Angs)  : ',params_glob%smpd,'/',pool_dims%smpd
        write(logfhandle,'(A,F5.1)')        '>>> POOL   HARD RESOLUTION LIMIT (Angs) : ',params_glob%lpstop
        ! upsample cavgs
        ldim = [pool_dims%box,pool_dims%box,1]
        str_tmp_mrc = str_tmp_mrc
        if( l_wfilt )then
            str  = add2fbody(refs_glob, params_glob%ext, WFILT_SUFFIX)
            call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
            call simple_rename(str_tmp_mrc,str)
            str  = add2fbody(refs_glob, params_glob%ext,WFILT_SUFFIX//'_even')
            call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
            call simple_rename(str_tmp_mrc,str)
            str  = add2fbody(refs_glob, params_glob%ext,WFILT_SUFFIX//'_odd')
            call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
            call simple_rename(str_tmp_mrc,str)
        endif
        call scale_imgfile(refs_glob, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
        call simple_rename(str_tmp_mrc,refs_glob)
        str  = add2fbody(refs_glob, params_glob%ext,'_even')
        call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
        call simple_rename(str_tmp_mrc,str)
        str  = add2fbody(refs_glob, params_glob%ext,'_odd')
        call scale_imgfile(str, str_tmp_mrc, prev_dims%smpd, ldim, pool_dims%smpd)
        call simple_rename(str_tmp_mrc,str)
        ! upsample cavgs matrices
        ldim = [pool_dims%boxpd,pool_dims%boxpd,1]
        do p = 1,params_glob%nparts_pool
            str = POOL_DIR//'cavgs_even_part'//int2str_pad(p,numlen)//params_glob%ext%to_char()
            call pad_imgfile(str, str_tmp_mrc, ldim, pool_dims%smpd)
            call simple_rename(str_tmp_mrc,str)
            str = POOL_DIR//'cavgs_odd_part'//int2str_pad(p,numlen)//params_glob%ext%to_char()
            call pad_imgfile(str, str_tmp_mrc, ldim, pool_dims%smpd)
            call simple_rename(str_tmp_mrc,str)
            str = POOL_DIR//'ctfsqsums_even_part'//int2str_pad(p,numlen)//params_glob%ext%to_char()
            call pad_imgfile(str, str_tmp_mrc, ldim, pool_dims%smpd)
            call simple_rename(str_tmp_mrc,str)
            str = POOL_DIR//'ctfsqsums_odd_part'//int2str_pad(p,numlen)//params_glob%ext%to_char()
            call pad_imgfile(str, str_tmp_mrc, ldim, pool_dims%smpd)
            call simple_rename(str_tmp_mrc,str)
        enddo
        if( l_wfilt )then
            do p = 1,params_glob%nparts_pool
                str = POOL_DIR//'cavgs_even_wfilt_part'//int2str_pad(p,numlen)//params_glob%ext%to_char()
                call pad_imgfile(str, str_tmp_mrc, ldim, pool_dims%smpd)
                call simple_rename(str_tmp_mrc,str)
                str = POOL_DIR//'cavgs_odd_wfilt_part'//int2str_pad(p,numlen)//params_glob%ext%to_char()
                call pad_imgfile(str, str_tmp_mrc, ldim, pool_dims%smpd)
                call simple_rename(str_tmp_mrc,str)
                str = POOL_DIR//'ctfsqsums_even_wfilt_part'//int2str_pad(p,numlen)//params_glob%ext%to_char()
                call pad_imgfile(str, str_tmp_mrc, ldim, pool_dims%smpd)
                call simple_rename(str_tmp_mrc,str)
                str = POOL_DIR//'ctfsqsums_odd_wfilt_part'//int2str_pad(p,numlen)//params_glob%ext%to_char()
                call pad_imgfile(str, str_tmp_mrc, ldim, pool_dims%smpd)
                call simple_rename(str_tmp_mrc,str)
            enddo
        endif
        ! update cls2D field
        os = pool_proj%os_cls2D
        call pool_proj%os_out%kill
        call pool_proj%add_cavgs2os_out(refs_glob, pool_dims%smpd, 'cavg', clspath=.true.)
        if( l_wfilt )then
            str = add2fbody(refs_glob,params_glob%ext,WFILT_SUFFIX)
            call pool_proj%add_cavgs2os_out(str, pool_dims%smpd, 'cavg'//WFILT_SUFFIX, clspath=.true.)
        endif
        pool_proj%os_cls2D = os
        call os%kill
        ! rescale frcs
        call frcs%read(string(FRCS_FILE))
        call frcs%pad(pool_dims%smpd, pool_dims%box)
        call frcs%write(string(FRCS_FILE))
        call frcs%kill
        call pool_proj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
    end subroutine update_pool_dims

    !> Updates the project watched by the gui for display
    subroutine update_pool_for_gui
        type(oris)        :: os_backup
        type(starproject) :: starproj
        type(string)      :: src
        os_backup = pool_proj%os_cls2D
        call pool_proj%add_cavgs2os_out(string(POOL_DIR)//refs_glob, pool_dims%smpd, 'cavg')
        if( l_wfilt )then
            src = string(POOL_DIR)//add2fbody(refs_glob,params_glob%ext,WFILT_SUFFIX)
            call pool_proj%add_cavgs2os_out(src, pool_dims%smpd, 'cavg'//WFILT_SUFFIX)
        endif
        pool_proj%os_cls2D = os_backup
        call pool_proj%write_segment_inside('out',   projfile4gui)
        call pool_proj%write_segment_inside('cls2D', projfile4gui)
        ! Write star file for iteration
        call starproj%export_cls2D(pool_proj, pool_iter)
        call pool_proj%os_cls2D%delete_entry('stk')
        call os_backup%kill
        call starproj%kill
    end subroutine update_pool_for_gui

    ! Flags pool availibility & updates the global name of references
    subroutine update_pool_status
        if( .not. l_stream2D_active ) return
        if( .not. l_pool_available )then
            l_pool_available = file_exists(POOL_DIR//CLUSTER2D_FINISHED)
            if( l_pool_available .and. (pool_iter >= 1) )then
                refs_glob = CAVGS_ITER_FBODY//int2str_pad(pool_iter,3)//params_glob%ext%to_char()
            endif
        endif
    end subroutine update_pool_status

    ! JPGS / GUI

    ! write jpeg of refs_glob    
    subroutine generate_pool_jpeg(filename)
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
        do icls=1, ncls_here
            if(pool_proj%os_cls2D%get(icls,'state') < 0.5) cycle
            if(pool_proj%os_cls2D%get(icls,'pop')   < 0.5) cycle
            pool_jpeg_map = [pool_jpeg_map, icls]
            pool_jpeg_pop = [pool_jpeg_pop, nint(pool_proj%os_cls2D%get(icls,'pop'))]
            pool_jpeg_res = [pool_jpeg_res, pool_proj%os_cls2D%get(icls,'res')]
            call img%zero_and_unflag_ft
            call stkio_r%get_image(icls, img)
            call img%mask(params_glob%mskdiam / (2 * pool_dims%smpd), 'softavg')
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
        current_jpeg_scale = real(JPEG_DIM) / real(params_glob%box)
        call img%kill()
        call img_pad%kill()
        call img_jpeg%kill()
    end subroutine generate_pool_jpeg

    ! Reports stats to GUI
    subroutine generate_pool_stats
        type(guistats) :: pool_stats
        type(string)   :: cwd
        integer        :: iter_loc = 0
        real           :: lpthres_sugg
        call simple_getcwd(cwd)
        if(file_exists(cwd//'/'//CLS2D_STARFBODY//'_iter'//int2str_pad(pool_iter,3)//STAR_EXT)) then
            iter_loc = pool_iter
        else if(file_exists(cwd//'/'//CLS2D_STARFBODY//'_iter'//int2str_pad(pool_iter - 1,3)//STAR_EXT)) then
            iter_loc = pool_iter - 1
        endif
        call pool_stats%init
        call pool_stats%set('particles', 'particles_assigned',  int2str(nptcls_glob - nptcls_rejected_glob) // '_(' // int2str(ceiling(100.0 * real(nptcls_glob - nptcls_rejected_glob) / real(nptcls_glob))) // '%)')
        call pool_stats%set('particles', 'particles_rejected',  int2str(nptcls_rejected_glob) // '_(' // int2str(floor(100.0 * real(nptcls_rejected_glob) / real(nptcls_glob))) // '%)')
        call pool_stats%set('2D', 'iteration',                  iter_loc,             primary=.true.)
        call pool_stats%set('2D', 'number_classes',             ncls_glob,            primary=.true.)
        call pool_stats%set('2D', 'number_classes_rejected',    ncls_rejected_glob,   primary=.true.)
        call mskdiam2streamresthreshold(params_glob%mskdiam, lpthres_sugg)
        if(current_resolution < lpthres_sugg / 4 .and. params_glob%lpthres > lpthres_sugg) then
            call pool_stats%set('2D', 'maximum_resolution', current_resolution, primary=.true., alert=.true., alerttext='maximum_resolution_suggests_&
            &using_a_cutoff_of_' // int2str(int(lpthres_sugg)) // 'Å_or_better' , notify=.false.)
        else
            call pool_stats%set('2D', 'maximum_resolution', current_resolution, primary=.true., alert=.false., notify=.true.)
        endif
        if(.not. file_exists(cwd//'/'//CAVGS_ITER_FBODY//int2str_pad(iter_loc, 3)//'.jpg')) then
            if(iter_loc > 0) call pool_stats%set_now('2D', 'iteration_time')
            call pool_stats%generate_2D_thumbnail('2D', 'top_classes', pool_proj%os_cls2D, iter_loc)
            call pool_stats%generate_2D_jpeg('latest', '', pool_proj%os_cls2D, iter_loc, pool_dims%smpd)
            last_complete_iter = iter_loc
            !call pool_stats%generate_2D_jpeg('latest', '', pool_proj%os_cls2D, iter_loc, smpd)
            ! nice
            call set_iteration_time() ! called in wrong place
            ! current_jpeg = trim(adjustl(cwd)) // '/' // CAVGS_ITER_FBODY // int2str_pad(iter_loc, 3) // '.jpg'
            call generate_pool_jpeg()
        endif
        call pool_stats%write(string(POOLSTATS_FILE))
        call pool_stats%kill
 
        contains

        subroutine set_iteration_time()
            character(8)  :: date
            character(10) :: time
            character(5)  :: zone
            integer,dimension(8) :: values
            ! using keyword arguments
            call date_and_time(date,time,zone,values)
            call date_and_time(DATE=date,ZONE=zone)
            call date_and_time(TIME=time)
            call date_and_time(VALUES=values)
            write(last_iteration_time, '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') values(1), '/', values(2), '/', values(3), '_', values(5), ':', values(6)
        end subroutine set_iteration_time

    end subroutine generate_pool_stats

end module simple_stream_pool2D_utils