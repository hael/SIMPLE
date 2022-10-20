! concrete commander: cluster2D_stream for streaming 2D alignment and clustering of single-particle images
module simple_commander_cluster2D_stream_dev
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters, params_glob
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_image,          only: image
use simple_oris,           only: oris
use simple_stream_chunk,   only: stream_chunk
use simple_class_frcs,     only: class_frcs
use simple_stack_io,       only: stack_io
use simple_starproject,    only: starproject
use simple_qsys_funs
use simple_commander_cluster2D
use simple_timer
implicit none

public :: cluster2D_commander_subsets
public :: init_cluster2D_stream, update_projects_mask, write_project_stream2D, terminate_stream2D
public :: update_pool_status, update_pool, reject_from_pool, reject_from_pool_user, classify_pool
public :: update_chunks, classify_new_chunks, import_chunks_into_pool, update_user_params
public :: update_path

private
#include "simple_local_flags.inc"

type, extends(commander_base) :: cluster2D_commander_subsets
  contains
    procedure :: execute      => exec_cluster2D_subsets
end type cluster2D_commander_subsets

integer,               parameter   :: MINBOXSZ            = 128    ! minimum boxsize for scaling
real,                  parameter   :: GREEDY_TARGET_LP    = 15.0
! integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 600  ! dev settings
integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 7200  ! Frequency at which the original project file should be updated
integer,               parameter   :: FREQ_POOL_REJECTION = 5     !
character(len=STDLEN), parameter   :: USER_PARAMS         = 'stream2D_user_params.txt'
character(len=STDLEN), parameter   :: PROJFILE_POOL       = 'cluster2D.simple'
character(len=STDLEN), parameter   :: SCALE_DIR           = './scaled_stks/'
character(len=STDLEN), parameter   :: POOL_DIR            = '' ! should be './pool/' for tidyness but difficult with gui
character(len=STDLEN), parameter   :: SNAPSHOT_DIR        = './snapshot/'
character(len=STDLEN), parameter   :: DISTR_EXEC_FNAME    = './distr_cluster2D_pool'
logical,               parameter   :: DEBUG_HERE          = .false.
integer(timer_int_kind)            :: t

! Pool related
type(sp_project)                       :: pool_proj
type(qsys_env)                         :: qenv_pool
type(cmdline)                          :: cline_cluster2D_pool
type(starproject)                      :: starproj
logical,                   allocatable :: pool_stacks_mask(:)
integer                                :: pool_iter
logical                                :: pool_available
! Chunk related
type(stream_chunk),        allocatable :: chunks(:), converged_chunks(:)
type(cmdline)                          :: cline_cluster2D_chunk
integer                                :: glob_chunk_id
! Book-keeping
logical,                   allocatable :: spprojs_mask_glob(:) ! micrographs from preprocess to import
character(len=LONGSTRLEN), allocatable :: imported_stks(:)
character(len=LONGSTRLEN)              :: prev_snapshot_frcs, prev_snapshot_cavgs
character(len=:),          allocatable :: orig_projfile
integer                                :: origproj_time, n_spprojs_glob = 0
logical                                :: initiated = .false.
! Global parameters to avoid conflict with preprocess_stream
integer                                :: numlen
! GUI-related
character(len=:),          allocatable :: projfile4gui
! other
character(len=STDLEN) :: refs_glob
real                  :: orig_smpd, smpd, scale_factor, mskdiam     ! dimensions
integer               :: orig_box, box, boxpd
real                  :: lp_greedy, lpstart_stoch                   ! resolution limits
integer               :: max_ncls, nptcls_per_chunk, ncls_glob
logical               :: l_wfilt, do_autoscale, l_greedy

contains

    subroutine init_cluster2D_stream( cline, spproj, refs, projfilegui, do2D )
        class(cmdline),    intent(inout) :: cline
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: refs
        character(len=*),  intent(in)    :: projfilegui
        logical,           intent(inout) :: do2D
        character(len=:), allocatable :: carg
        real, parameter :: CENLP_2DSTREAM = 30.
        real    :: SMPD_TARGET = MAX_SMPD  ! target sampling distance
        real    :: rarg
        integer :: ldim(3), ichunk, maybe2D, ifoo
        ! check on strictly required parameters
        if( .not.cline%defined('nthr2D') )then
            THROW_HARD('Missing required agrument NTHR2D')
        endif
        ! check whether 2D classification will be performed based on 6 strictly required parameters
        maybe2D = merge(1,0,cline%defined('ncls'))
        maybe2D = maybe2D + merge(1,0,cline%defined('nptcls_per_cls'))
        maybe2D = maybe2D + merge(1,0,cline%defined('ncls_start'))
        maybe2D = maybe2D + merge(1,0,cline%defined('mskdiam'))
        if( maybe2D == 4 )then
            do2D = .true.
        else if( maybe2D > 0 )then
            THROW_HARD('Missing arguments for 2D classification')
        else
            do2D = .false.
        endif
        if( .not.do2D ) return
        call seed_rnd
        ! general parameters
        mskdiam             = cline%get_rarg('mskdiam')
        l_wfilt             = trim(params_glob%wiener) .eq. 'partial'
        do_autoscale        = trim(params_glob%autoscale) .eq. 'yes'
        max_ncls            = floor(cline%get_rarg('ncls')/real(params_glob%ncls_start))*params_glob%ncls_start ! effective maximum # of classes
        nptcls_per_chunk    = params_glob%nptcls_per_cls*params_glob%ncls_start         ! # of particles in each chunk
        ncls_glob           = 0
        l_greedy            = trim(params_glob%refine).eq.'greedy'
        lp_greedy           = GREEDY_TARGET_LP
        if( cline%defined('lp2D') ) lp_greedy = params_glob%lp2D
        lpstart_stoch       = lp_greedy
        prev_snapshot_cavgs = ''
        orig_projfile       = trim(params_glob%projfile)
        projfile4gui        = trim(projfilegui)
        if( cline%defined('box_extract') )then
            orig_box = nint(cline%get_rarg('box_extract'))
        else
            call find_ldim_nptcls(refs,ldim,ifoo)
            orig_box  = ldim(1)
        endif
        orig_smpd = params_glob%smpd / params_glob%scale
        call debug_print('cluster2D_stream orig_box: '//int2str(orig_box))
        call debug_print('cluster2D_stream orig_smpd: '//real2str(orig_smpd))
        call debug_print('cluster2D_stream mskdiam: '//real2str(mskdiam))
        numlen         = len(int2str(params_glob%nparts_pool))
        refs_glob      = 'start_cavgs'//params_glob%ext
        pool_available = .true.
        pool_iter      = 0
        call simple_mkdir(SCALE_DIR)
        call simple_mkdir(POOL_DIR)
        call simple_mkdir(trim(POOL_DIR)//trim(STDERROUT_DIR))
        call simple_mkdir(trim(SNAPSHOT_DIR))
        call simple_touch(trim(POOL_DIR)//trim(CLUSTER2D_FINISHED))
        call pool_proj%kill
        pool_proj%projinfo = spproj%projinfo
        pool_proj%compenv  = spproj%compenv
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        ! update to computational parameters to pool, will be transferred to chunks upon init
        if( cline%defined('job_memory_per_task2D') )then
            call pool_proj%compenv%set(1,'job_memory_per_task', real(params_glob%job_memory_per_task2D))
        endif
        if( cline%defined('qsys_partition2D') )then
            call pool_proj%compenv%set(1,'qsys_partition', params_glob%qsys_partition2D)
        endif
        call pool_proj%write(trim(POOL_DIR)//trim(PROJFILE_POOL))
        ! initialize chunks parameters and objects
        if( params_glob%nparts_chunk > 1 )then
            call cline_cluster2D_chunk%set('prg',    'cluster2D_distr')
            call cline_cluster2D_chunk%set('nparts', real(params_glob%nparts_chunk))
        else
            ! shared memory execution
            call cline_cluster2D_chunk%set('prg',       'cluster2D')
        endif
        call cline_cluster2D_chunk%set('oritype',   'ptcl2D')
        call cline_cluster2D_chunk%set('objfun',    'cc')
        call cline_cluster2D_chunk%set('center',    'no')
        call cline_cluster2D_chunk%set('match_filt','no')
        call cline_cluster2D_chunk%set('autoscale', 'no')
        call cline_cluster2D_chunk%set('ptclw',     'no')
        call cline_cluster2D_chunk%set('mkdir',     'no')
        call cline_cluster2D_chunk%set('stream',    'no')
        call cline_cluster2D_chunk%set('startit',   1.)
        call cline_cluster2D_chunk%set('mskdiam',   mskdiam)
        call cline_cluster2D_chunk%set('ncls',      real(params_glob%ncls_start))
        call cline_cluster2D_chunk%set('nthr',      real(params_glob%nthr2D))
        call cline_cluster2D_chunk%set('nonuniform',params_glob%nonuniform)
        call cline_cluster2D_chunk%set('nsearch',   real(params_glob%nsearch))
        call cline_cluster2D_chunk%set('smooth_ext',real(params_glob%smooth_ext))
        call cline_cluster2D_chunk%set('lp_lowres', real(params_glob%lp_lowres))
        if( l_wfilt ) call cline_cluster2D_chunk%set('wiener', 'partial')
        allocate(chunks(params_glob%nchunks))
        do ichunk = 1,params_glob%nchunks
            call chunks(ichunk)%init(ichunk, pool_proj)
        enddo
        glob_chunk_id = params_glob%nchunks
        ! initialize pool parameters and objects
        call cline_cluster2D_pool%set('prg',       'cluster2D_distr')
        call cline_cluster2D_pool%set('oritype',   'ptcl2D')
        call cline_cluster2D_pool%set('autoscale', 'no')
        call cline_cluster2D_pool%set('trs',       MINSHIFT)
        call cline_cluster2D_pool%set('projfile',  trim(PROJFILE_POOL))
        call cline_cluster2D_pool%set('projname',  trim(get_fbody(trim(PROJFILE_POOL),trim('simple'))))
        call cline_cluster2D_pool%set('objfun',    'cc')
        call cline_cluster2D_pool%set('ptclw',     'no')
        call cline_cluster2D_pool%set('nonuniform',params_glob%nonuniform)
        call cline_cluster2D_pool%set('nsearch',   real(params_glob%nsearch))
        call cline_cluster2D_pool%set('smooth_ext',real(params_glob%smooth_ext))
        call cline_cluster2D_pool%set('lp_lowres',real(params_glob%lp_lowres))
        if( cline%defined('match_filt') )then
            carg = cline%get_carg('match_filt')
            call cline_cluster2D_pool%set('match_filt',carg)
            deallocate(carg)
        else
            call cline_cluster2D_pool%set('match_filt','no')
        endif
        if( cline%defined('cenlp') )then
            rarg = cline%get_rarg('cenlp')
            call cline_cluster2D_pool%set('cenlp', rarg)
        else
            call cline_cluster2D_pool%set('cenlp', CENLP_2DSTREAM)
        endif
        if( cline%defined('center') )then
            carg = cline%get_carg('center')
            call cline_cluster2D_pool%set('center',carg)
            deallocate(carg)
        else
            call cline_cluster2D_pool%set('center','yes')
        endif
        if( l_wfilt ) call cline_cluster2D_pool%set('wiener', 'partial')
        call cline_cluster2D_pool%set('extr_iter', 100.)
        call cline_cluster2D_pool%set('mkdir',     'no')
        call cline_cluster2D_pool%set('mskdiam',   mskdiam)
        call cline_cluster2D_pool%set('async',     'yes') ! to enable hard termination
        call cline_cluster2D_pool%set('stream',    'yes') ! use for dual CTF treatment
        call cline_cluster2D_pool%set('nthr',     real(params_glob%nthr2D))
        call cline_cluster2D_pool%set('nparts',   real(params_glob%nparts_pool))
        call qenv_pool%new(params_glob%nparts_pool,exec_bin='simple_private_exec',qsys_name='local')
        ! auto-scaling
        if( orig_box == 0 ) THROW_HARD('FATAL ERROR')
        if( do_autoscale )then
            if( orig_box < MINBOXSZ )then
                do_autoscale = .false.
            else
                call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                if( box < MINBOXSZ )then
                    SMPD_TARGET = orig_smpd * real(orig_box) / real(MINBOXSZ)
                    call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                endif
                if( box >= orig_box )then
                    do_autoscale = .false.
                else
                    write(logfhandle,'(A,I5,F5.2)') '>>> 2D CLASSIFICATION DOWNSCALED IMAGE SIZE & PIXEL SIZE (IN A): ', box, smpd
                endif
            endif
        endif
        if( .not. do_autoscale )then
            smpd         = orig_smpd
            box          = orig_box
            scale_factor = 1.
        endif
        boxpd = 2 * round2even(params_glob%alpha * real(box/2)) ! logics from parameters
        call cline_cluster2D_chunk%set('box',  real(box))
        call cline_cluster2D_chunk%set('smpd', smpd)
        call cline_cluster2D_pool%set('box',   real(box))
        call cline_cluster2D_pool%set('smpd',  smpd)
        ! resolution-related updates to command-lines
        lp_greedy     = max(lp_greedy,    2.0*smpd)
        lpstart_stoch = max(lpstart_stoch,2.0*smpd)
        if( cline%defined('lpstop2D') )then
            params_glob%lpstop2D = max(2.0*smpd,params_glob%lpstop2D)
        else
            params_glob%lpstop2D = 2.0*smpd
        endif
        call cline_cluster2D_chunk%set('lp', lp_greedy)
        if( l_greedy )then
            call cline_cluster2D_chunk%set('maxits', 10.)
            call cline_cluster2D_chunk%set('refine', 'greedy')
        else
            call cline_cluster2D_chunk%set('refine',    'snhc')
            call cline_cluster2D_chunk%set('extr_iter', real(MAX_EXTRLIM2D-2))
            call cline_cluster2D_chunk%set('maxits',    12.)
        endif
        write(logfhandle,'(A,F5.1)')     '>>> CHUNK         LOW-PASS LIMIT (IN A) TO: ', lp_greedy
        if( l_greedy )then
            call cline_cluster2D_pool%set('refine', 'greedy')
            call cline_cluster2D_pool%set('lp',     lp_greedy)
            write(logfhandle,'(A,F5.1)') '>>> POOL          LOW-PASS LIMIT (IN A) TO: ', lp_greedy
        else
            call cline_cluster2D_pool%set('lpstart', lpstart_stoch)
            call cline_cluster2D_pool%set('lpstop',  params_glob%lpstop2D)
            write(logfhandle,'(A,F5.1)') '>>> POOL STARTING LOW-PASS LIMIT (IN A) TO: ', lpstart_stoch
        endif
        write(logfhandle,'(A,F5.1)')     '>>> POOL   HARD RESOLUTION LIMIT (IN A) TO: ', params_glob%lpstop2D
        ! module variables
        n_spprojs_glob = 0
        if( allocated(spprojs_mask_glob) ) deallocate(spprojs_mask_glob)
        origproj_time  = simple_gettime()
        initiated      = .true.
    end subroutine init_cluster2D_stream

    ! updates mask  of project to process
    subroutine update_projects_mask( spproj_fnames )
        character(len=LONGSTRLEN), allocatable, intent(in) :: spproj_fnames(:)
        logical, allocatable :: tmp_mask(:)
        integer :: n_spprojs_in
        if( .not.allocated(spproj_fnames) ) return
        if( .not.initiated ) return
        n_spprojs_in = size(spproj_fnames)
        if( n_spprojs_in == 0 ) return
        ! updates global mask of projects imported
        if( n_spprojs_glob == 0 )then
            ! first time
            n_spprojs_glob = n_spprojs_in
            allocate(spprojs_mask_glob(n_spprojs_glob),source=.false.)
        endif
        if( n_spprojs_in > n_spprojs_glob )then
            ! addition of new projects
            allocate(tmp_mask(n_spprojs_in),source=.false.)
            tmp_mask(1:n_spprojs_glob) = spprojs_mask_glob
            call move_alloc(tmp_mask, spprojs_mask_glob)
            n_spprojs_glob = n_spprojs_in
        else if( n_spprojs_in == n_spprojs_glob )then
            ! nothing to do
        else
            THROW_HARD('Error update_projects_mask 1')
        endif
        ! call debug_print('update_chunk_mask n_spprojs_glob: '//int2str(n_spprojs_glob))
    end subroutine update_projects_mask

    ! deals with chunk completion, rejection, reset
    subroutine update_chunks
        integer :: ichunk
        do ichunk = 1,params_glob%nchunks
            if( chunks(ichunk)%available ) cycle
            if( chunks(ichunk)%has_converged() )then
                call chunks(ichunk)%display_iter
                ! rejection
                call chunks(ichunk)%reject(params_glob%lpthres, params_glob%ndev2D, box)
                ! updates list of chunks to import
                if( allocated(converged_chunks) )then
                    converged_chunks = [converged_chunks(:), chunks(ichunk)]
                else
                    allocate(converged_chunks(1),source=[chunks(ichunk)])
                endif
                ! reinit
                glob_chunk_id = glob_chunk_id + 1
                call chunks(ichunk)%init(glob_chunk_id, pool_proj)
            endif
        enddo
    end subroutine update_chunks

    subroutine classify_new_chunks( spproj_fnames )
        character(len=LONGSTRLEN), allocatable, intent(in) :: spproj_fnames(:)
        type(sp_project)                       :: stream_proj
        character(len=LONGSTRLEN), allocatable :: spprojs_for_chunk(:)
        integer,                   allocatable :: spproj_nptcls(:)
        integer :: ichunk, n_avail_chunks, n_spprojs_in, iproj, nptcls, n2fill
        integer :: first2import, last2import, n2import, cnt
        integer :: inmics, instks, inptcls
        n_avail_chunks = count(chunks(:)%available)
        ! cannot import yet
        if( n_avail_chunks == 0 ) return
        if( .not.allocated(spprojs_mask_glob) )return
        n_spprojs_in = size(spprojs_mask_glob)
        if( n_spprojs_in == 0 ) return
        ! how many n2fill chunks to load
        allocate(spproj_nptcls(n_spprojs_in),source=0)
        n2fill       = 0
        nptcls       = 0
        first2import = 0
        do iproj = 1,n_spprojs_in
            if( spprojs_mask_glob(iproj) )cycle
            call stream_proj%read_data_info(spproj_fnames(iproj), inmics, instks, inptcls)
            spproj_nptcls(iproj) = inptcls
            if( spproj_nptcls(iproj) > 0 )then
                if( first2import == 0 ) first2import = iproj
                nptcls = nptcls + spproj_nptcls(iproj)
                if( nptcls > nptcls_per_chunk )then
                    n2fill = n2fill+1
                    if( n2fill >= n_avail_chunks )exit
                    nptcls = 0
                endif
            else
                spprojs_mask_glob(iproj) = .true. ! mask out empty projects
            endif
        enddo
        call stream_proj%kill
        if( n2fill == 0 ) return ! not enough particles
        do ichunk = 1,params_glob%nchunks
            if(.not.chunks(ichunk)%available) cycle
            if( n2fill == 0 ) exit
            n2fill   = n2fill - 1
            nptcls   = 0
            n2import = 0
            do iproj = first2import,n_spprojs_in
                if( spproj_nptcls(iproj) == 0 )cycle
                nptcls   = nptcls + spproj_nptcls(iproj)
                n2import = n2import + 1
                if( nptcls > nptcls_per_chunk )then
                    last2import = iproj
                    exit
                endif
            enddo
            if( nptcls > nptcls_per_chunk )then
                ! actual import
                allocate(spprojs_for_chunk(n2import))
                cnt = 0
                do iproj = first2import,last2import
                    if( spproj_nptcls(iproj) == 0 )cycle
                    cnt = cnt + 1
                    spprojs_for_chunk(cnt)   = trim( spproj_fnames(iproj) )
                    spprojs_mask_glob(iproj) = .true. ! mask out from future imports
                enddo
                call chunks(ichunk)%generate(spprojs_for_chunk, nptcls, 1)
                deallocate(spprojs_for_chunk)
                ! execution
                call chunks(ichunk)%exec_classify(cline_cluster2D_chunk, orig_smpd, orig_box, box)
                ! to avoid cycling through all projects
                first2import = last2import + 1
            endif
        enddo
    end subroutine classify_new_chunks

    subroutine import_chunks_into_pool( ignore_mic, nchunks_imported )
        logical, intent(in)  :: ignore_mic
        integer, intent(out) :: nchunks_imported
        type(class_frcs)                   :: frcs_glob, frcs_chunk, frcs_prev
        character(LONGSTRLEN), allocatable :: tmp(:)
        character(len=:),      allocatable :: cavgs_chunk, dir_chunk
        integer,               allocatable :: cls_pop(:), cls_chunk_pop(:), pinds(:), states(:)
        character(len=LONGSTRLEN) :: stk_relpath
        real    :: smpd_here
        integer :: ichunk, nchunks2import, nptcls2import, nmics2import, nptcls, imic, iproj
        integer :: ncls_tmp, fromp_prev, fromp, ii, jptcl, i, poolind, n_remap, pop, nptcls_sel
        integer :: nmics_imported, nptcls_imported, iptcl, ind, ncls_here, icls
        logical :: l_maxed
        nchunks_imported = 0
        if( .not.pool_available ) return
        if( .not.allocated(converged_chunks) ) return
        nchunks2import = size(converged_chunks)
        if( nchunks2import == 0 ) return
        nptcls2import   = 0
        nmics2import    = 0
        do ichunk = 1,nchunks2import
            nptcls2import = nptcls2import + converged_chunks(ichunk)%nptcls
            nmics2import  = nmics2import  + converged_chunks(ichunk)%nmics
        enddo
        if( nptcls2import == 0 ) return
        call debug_print('in import_chunks_into_pool 1 '//int2str(nptcls2import))
        ! reallocations
        nmics_imported  = pool_proj%os_mic%get_noris()
        nptcls_imported = pool_proj%os_ptcl2D%get_noris()
        if( nmics_imported == 0 )then
            call pool_proj%os_mic%new(nmics2import, is_ptcl=.false.)
            call pool_proj%os_stk%new(nmics2import, is_ptcl=.false.)
            call pool_proj%os_ptcl2D%new(nptcls2import, is_ptcl=.true.)
            allocate(imported_stks(nmics2import))
            fromp = 1
        else
            call pool_proj%os_mic%reallocate(nmics_imported+nmics2import)
            call pool_proj%os_stk%reallocate(nmics_imported+nmics2import)
            call pool_proj%os_ptcl2D%reallocate(nptcls_imported+nptcls2import)
            fromp = nint(pool_proj%os_stk%get(nmics_imported,'top'))+1
            call move_alloc(imported_stks, tmp)
            allocate(imported_stks(nmics_imported+nmics2import))
            imported_stks(1:nmics_imported) = tmp(:)
            deallocate(tmp)
        endif
        imic  = nmics_imported
        iptcl = nptcls_imported
        call debug_print('in import_chunks_into_pool 2 '//int2str(imic)//' '//int2str(iptcl))
        do ichunk = 1,nchunks2import
            fromp_prev = fromp
            call converged_chunks(ichunk)%read(boxpd)
            ! transfer micrographs, stacks & particles parameters
            jptcl = 0
            do iproj=1,converged_chunks(ichunk)%nmics
                imic  = imic+1
                ! micrographs & update paths so they are relative to root folder
                call pool_proj%os_mic%transfer_ori(imic, converged_chunks(ichunk)%spproj%os_mic, iproj)
                if( .not.ignore_mic )then
                    call update_path(pool_proj%os_mic, imic, 'mc_starfile')
                    call update_path(pool_proj%os_mic, imic, 'intg')
                    call update_path(pool_proj%os_mic, imic, 'forctf')
                    call update_path(pool_proj%os_mic, imic, 'thumb')
                    call update_path(pool_proj%os_mic, imic, 'mceps')
                    call update_path(pool_proj%os_mic, imic, 'ctfdoc')
                    call update_path(pool_proj%os_mic, imic, 'ctfjpg')
                    call update_path(pool_proj%os_mic, imic, 'boxfile')
                endif
                ! stacks & update path so they are relative to root folder
                call pool_proj%os_stk%transfer_ori(imic, converged_chunks(ichunk)%spproj%os_stk, iproj)
                nptcls = nint(converged_chunks(ichunk)%spproj%os_stk%get(iproj,'nptcls'))
                call pool_proj%os_stk%set(imic, 'fromp', real(fromp))
                call pool_proj%os_stk%set(imic, 'top',   real(fromp+nptcls-1))
                call make_relativepath(cwd_glob, converged_chunks(ichunk)%orig_stks(iproj), stk_relpath)
                imported_stks(imic) = trim(stk_relpath)
                ! particles
                do i = 1,nptcls
                    iptcl = iptcl + 1
                    jptcl = jptcl + 1
                    call pool_proj%os_ptcl2D%transfer_ori(iptcl, converged_chunks(ichunk)%spproj%os_ptcl2D, jptcl)
                    call pool_proj%os_ptcl2D%set(iptcl, 'stkind',    real(imic))
                    call pool_proj%os_ptcl2D%set(iptcl, 'updatecnt', 0.) ! new particle
                    call pool_proj%os_ptcl2D%set(iptcl, 'frac',      0.) ! new particle
                enddo
                fromp = fromp + nptcls
            enddo
            nptcls_sel = converged_chunks(ichunk)%spproj%os_ptcl2D%get_noris(consider_state=.true.)
            write(logfhandle,'(A,I6,A,I6,A)')'>>> TRANSFERRED ',nptcls_sel,' PARTICLES FROM CHUNK ',converged_chunks(ichunk)%id,' TO POOL'
            call flush(logfhandle)
            ! transfer classes
            l_maxed   = ncls_glob >= max_ncls ! max # of classes reached ?
            ncls_here = ncls_glob
            if( .not.l_maxed ) ncls_here = ncls_glob + params_glob%ncls_start
            states = nint(converged_chunks(ichunk)%spproj%os_ptcl2D%get_all('state'))
            call converged_chunks(ichunk)%spproj%get_cavgs_stk(cavgs_chunk, ncls_tmp, smpd_here)
            dir_chunk   = trim(converged_chunks(ichunk)%path)
            cavgs_chunk = trim(dir_chunk)//basename(cavgs_chunk)
            call frcs_chunk%new(params_glob%ncls_start, box, smpd, nstates=1)
            call frcs_chunk%read(trim(dir_chunk)//trim(FRCS_FILE))
            if( l_maxed )then
                call debug_print('in import_chunks_into_pool 3 '//int2str(ichunk))
                ! transfer all others
                cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
                n_remap = 0
                if( any(cls_pop==0) )then
                    if( all(cls_pop==0) ) THROW_HARD('Empty os_cls2D!')
                    ! remapping
                    cls_chunk_pop = nint(converged_chunks(ichunk)%spproj%os_cls2D%get_all('pop'))
                    call frcs_glob%new(ncls_glob, box, smpd, nstates=1)
                    call frcs_glob%read(trim(POOL_DIR)//trim(FRCS_FILE))
                    do icls=1,ncls_glob
                        if( cls_pop(icls)>0 ) cycle          ! class already filled
                        if( all(cls_chunk_pop == 0 ) ) exit  ! no more chunk class available
                        ind = irnd_uni(params_glob%ncls_start)    ! selects chunk class stochastically
                        do while( cls_chunk_pop(ind) == 0 )
                            ind = irnd_uni(params_glob%ncls_start)
                        enddo
                        cls_chunk_pop(ind) = 0             ! excludes from being picked again
                        n_remap = n_remap+1
                        ! class average
                        call transfer_cavg(cavgs_chunk, dir_chunk, ind, refs_glob, icls)
                        ! frcs
                        call frcs_glob%set_frc(icls,frcs_chunk%get_frc(ind, box, 1), 1)
                        ! class parameters transfer
                        call pool_proj%os_cls2D%transfer_ori(icls, converged_chunks(ichunk)%spproj%os_cls2D, ind)
                        call pool_proj%os_cls2D%set(icls, 'class', real(icls))
                        ! particles
                        call converged_chunks(ichunk)%spproj%os_ptcl2D%get_pinds(ind,'class',pinds,consider_w=.false.)
                        pop = size(pinds)
                        do i=1,pop
                            ii      = pinds(i)               ! in chunk
                            poolind = fromp_prev + ii - 1 ! in pool
                            call pool_proj%os_ptcl2D%set(poolind,'class',real(icls))
                        enddo
                        cls_pop(icls) = cls_pop(icls) + pop ! updates class populations
                    enddo
                    ! now transfer particles that were not remapped
                    do icls = 1,params_glob%ncls_start
                        if( cls_chunk_pop(icls) == 0 ) cycle
                        call converged_chunks(ichunk)%spproj%os_ptcl2D%get_pinds(icls,'class',pinds,consider_w=.false.)
                        pop = size(pinds)
                        do i = 1,pop
                            ii      = pinds(i)            ! in chunk
                            poolind = fromp_prev + ii - 1 ! in pool
                            ind     = irnd_uni(ncls_glob) ! stochastic labelling followed by greedy search
                            do while( cls_pop(ind) == 0 )
                                ind = irnd_uni(ncls_glob)
                            enddo
                            call pool_proj%os_ptcl2D%set(poolind,'class',real(ind))
                            cls_pop(ind) = cls_pop(ind) + 1 ! updates class populations
                        enddo
                    enddo
                    if( n_remap > 0 )then
                        call transfer_cavg(refs_glob,dir_chunk,ncls_glob,refs_glob,ncls_glob, self_transfer=.true.)
                        call frcs_glob%write(trim(POOL_DIR)//trim(FRCS_FILE))
                        write(logfhandle,'(A,I4)')'>>> # OF RE-MAPPED CLASS AVERAGES: ',n_remap
                    endif
                else
                    call debug_print('in import_chunks_into_pool 4 '//int2str(ichunk))
                    ! no remapping, just transfer particles & updates 2D population
                    do ii = 1,converged_chunks(ichunk)%nptcls
                        if( states(ii) /= 0 )then
                            icls          = irnd_uni(ncls_glob) ! stochastic labelling followed by greedy search
                            cls_pop(icls) = cls_pop(icls) + 1   ! updates populations
                            poolind       = fromp_prev + ii - 1
                            call pool_proj%os_ptcl2D%set(poolind, 'class', real(icls))
                        endif
                    enddo
                endif
                ! updates class populations
                call pool_proj%os_cls2D%set_all('pop',real(cls_pop))
                call debug_print('in import_chunks_into_pool 5 '//int2str(ichunk))
            else
                ! all new classes can be imported, no remapping
                if( ncls_glob == 0 )then
                    ! first transfer : copy classes, frcs & class parameters
                    refs_glob = 'start_cavgs'//params_glob%ext
                    do icls= 1,params_glob%ncls_start
                        call transfer_cavg(cavgs_chunk,dir_chunk,icls,refs_glob,icls)
                    enddo
                    call simple_copy_file(trim(dir_chunk)//trim(FRCS_FILE),trim(POOL_DIR)//trim(FRCS_FILE))
                    pool_proj%os_cls2D = converged_chunks(ichunk)%spproj%os_cls2D
                else
                    ! append new classes
                    do icls=1,params_glob%ncls_start
                        call transfer_cavg(cavgs_chunk, dir_chunk, icls, refs_glob, ncls_glob+icls)
                    enddo
                    ! FRCs
                    call frcs_glob%new(ncls_here, box, smpd, nstates=1)
                    call frcs_prev%new(ncls_glob, box, smpd, nstates=1)
                    call frcs_prev%read(trim(POOL_DIR)//trim(FRCS_FILE))
                    do icls=1,ncls_glob
                        call frcs_glob%set_frc(icls,frcs_prev%get_frc(icls, box, 1), 1)
                    enddo
                    do icls=1,params_glob%ncls_start
                        call frcs_glob%set_frc(ncls_glob+icls,frcs_chunk%get_frc(icls, box, 1), 1)
                    enddo
                    call frcs_glob%write(trim(POOL_DIR)//trim(FRCS_FILE))
                    ! class parameters
                    call pool_proj%os_cls2D%reallocate(ncls_here)
                    do icls = 1,params_glob%ncls_start
                        ind = ncls_glob+icls
                        call pool_proj%os_cls2D%transfer_ori(ind, converged_chunks(ichunk)%spproj%os_cls2D, icls)
                        call pool_proj%os_cls2D%set(ind, 'class', real(ind))
                    enddo
                endif
                call debug_print('in import_chunks_into_pool 7'//' '//int2str(ichunk))
                ! particles 2D
                do ii = 1,converged_chunks(ichunk)%nptcls
                    if( states(ii) /= 0 )then
                        poolind = fromp_prev+ii-1
                        icls    = ncls_glob+nint(converged_chunks(ichunk)%spproj%os_ptcl2D%get(ii,'class'))
                        call pool_proj%os_ptcl2D%set(poolind, 'class', real(icls))
                    endif
                enddo
                ! global # of classes
                ncls_glob = ncls_here
            endif
            ! tidy
            call converged_chunks(ichunk)%remove_folder
            call converged_chunks(ichunk)%kill
        enddo
        nchunks_imported = nchunks2import
        call update_pool_for_gui
        ! cleanup
        deallocate(converged_chunks)
        call frcs_prev%kill
        call frcs_glob%kill
        call frcs_chunk%kill
        call debug_print('end import_chunks_into_pool')
    end subroutine import_chunks_into_pool

    subroutine update_pool_status
        if( .not.pool_available )then
            pool_available = file_exists(trim(POOL_DIR)//trim(CLUSTER2D_FINISHED))
            if( pool_iter > 1 ) refs_glob = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter,3))//trim(params_glob%ext)
        endif
    end subroutine update_pool_status

    subroutine update_pool
        type(sp_project)      :: spproj
        type(oris)            :: os
        type(class_frcs)      :: frcs
        integer, allocatable  :: pops(:)
        character(len=STDLEN) :: fname
        real    :: corr, frac, mi_class
        integer :: i, it, jptcl, iptcl, istk
        if( .not.pool_available )return
        call del_file(trim(POOL_DIR)//trim(CLUSTER2D_FINISHED))
        ! iteration info
        fname = trim(POOL_DIR)//trim(STATS_FILE)
        if( file_exists(fname) )then
            call os%new(1,is_ptcl=.false.)
            call os%read(fname)
            it = nint(os%get(1,'ITERATION'))
            if( it == pool_iter )then
                mi_class = os%get(1,'CLASS_OVERLAP')
                frac     = os%get(1,'SEARCH_SPACE_SCANNED')
                corr     = os%get(1,'CORRELATION')
                write(logfhandle,'(A,I6,A,F7.3,A,F7.3,A,F7.3)')'>>> POOL         ITERATION ',it,&
                    &'; CLASS OVERLAP: ',mi_class,'; SEARCH SPACE SCANNED: ',frac,'; CORRELATION: ',corr
            endif
            call os%kill
        endif
        ! transfer to pool
        call spproj%read_segment('cls2D', trim(POOL_DIR)//trim(PROJFILE_POOL))
        if( spproj%os_cls2D%get_noris() == 0 )then
            ! not executed yet, do nothing
        else
            if( .not.allocated(pool_stacks_mask) )then
                THROW_HARD('Critical ERROR 0') ! first time
            endif
            call spproj%read_segment('stk',   trim(POOL_DIR)//trim(PROJFILE_POOL))
            call spproj%read_segment('ptcl2D',trim(POOL_DIR)//trim(PROJFILE_POOL))
            i = 0
            do istk = 1,size(pool_stacks_mask)
                if( pool_stacks_mask(istk) )then
                    i = i+1
                    iptcl = nint(pool_proj%os_stk%get(istk,'fromp'))
                    do jptcl = nint(spproj%os_stk%get(i,'fromp')),nint(spproj%os_stk%get(i,'top'))
                        if( spproj%os_ptcl2D%get_state(jptcl) > 0 )then
                            call pool_proj%os_ptcl2D%transfer_2Dparams(iptcl, spproj%os_ptcl2D, jptcl)
                        endif
                        iptcl = iptcl+1
                    enddo
                endif
            enddo
            call pool_proj%os_ptcl2D%get_pops(pops, 'class', consider_w=.false., maxn=ncls_glob)
            pool_proj%os_cls2D = spproj%os_cls2D
            call pool_proj%os_cls2D%set_all('pop', real(pops))
            if( .not.l_greedy )then
                call frcs%read(trim(POOL_DIR)//trim(FRCS_FILE))
                write(logfhandle,'(A,F5.1)')'>>> CURRENT POOOL RESOLUTION: ',frcs%estimate_lp_for_align()
                call frcs%kill
            endif
            ! for gui
            call update_pool_for_gui
        endif
        call spproj%kill
    end subroutine update_pool

    subroutine reject_from_pool
        type(image)          :: img
        logical, allocatable :: cls_mask(:)
        real                 :: ndev_here
        integer              :: nptcls_rejected, ncls_rejected, iptcl
        integer              :: icls, cnt
        if( .not.pool_available )return
        ! rejection frequency
        if( pool_iter <= 2*FREQ_POOL_REJECTION .or. mod(pool_iter,FREQ_POOL_REJECTION)/=0 ) return
        if( pool_proj%os_cls2D%get_noris() == 0 ) return
        ncls_rejected   = 0
        nptcls_rejected = 0
        allocate(cls_mask(ncls_glob), source=.true.)
        ! correlation & resolution
        ndev_here = 1.5*params_glob%ndev2D ! less stringent rejection than chunk
        call pool_proj%os_cls2D%find_best_classes(box,smpd,params_glob%lpthres,cls_mask,ndev_here)
        if( count(cls_mask) > 1 .and. count(cls_mask) < ncls_glob )then
            ncls_rejected = 0
            do iptcl=1,pool_proj%os_ptcl2D%get_noris()
                if( pool_proj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                icls = nint(pool_proj%os_ptcl2D%get(iptcl,'class'))
                if( cls_mask(icls) ) cycle
                nptcls_rejected = nptcls_rejected+1
                call pool_proj%os_ptcl2D%set_state(iptcl,0)
            enddo
            if( nptcls_rejected > 0 )then
                do icls=1,ncls_glob
                    if( .not.cls_mask(icls) )then
                        if( pool_proj%os_cls2D%get(icls,'pop') > 0.5 ) ncls_rejected = ncls_rejected+1
                        call pool_proj%os_cls2D%set(icls,'pop',0.)
                        call pool_proj%os_cls2D%set(icls,'corr',-1.)
                    endif
                enddo
                cnt = 0
                call img%new([box,box,1],smpd)
                do icls=1,ncls_glob
                    if( cls_mask(icls) ) cycle
                    cnt = cnt+1
                    if( debug_here )then
                        call img%read(trim(POOL_DIR)//trim(refs_glob),icls)
                        call img%write(trim(POOL_DIR)//'rejected_pool_'//int2str(pool_iter)//'.mrc',cnt)
                    endif
                    img = 0.
                    call img%write(trim(POOL_DIR)//trim(refs_glob),icls)
                    if( l_wfilt )then
                        call img%write(trim(POOL_DIR)//trim(add2fbody(refs_glob, params_glob%ext,trim(WFILT_SUFFIX))),icls)
                    endif
                enddo
                call img%read(trim(POOL_DIR)//trim(refs_glob), ncls_glob)
                call img%write(trim(POOL_DIR)//trim(refs_glob), ncls_glob)
                if( l_wfilt )then
                    call img%read(trim(POOL_DIR)//trim(add2fbody(refs_glob, params_glob%ext,trim(WFILT_SUFFIX))), ncls_glob)
                    call img%write(trim(POOL_DIR)//trim(add2fbody(refs_glob, params_glob%ext,trim(WFILT_SUFFIX))), ncls_glob)
                endif
                deallocate(cls_mask)
                write(logfhandle,'(A,I4,A,I6,A)')'>>> REJECTED FROM POOL: ',nptcls_rejected,' PARTICLES IN ',ncls_rejected,' CLUSTER(S)'
            endif
        else
            write(logfhandle,'(A,I4,A,I6,A)')'>>> NO PARTICLES FLAGGED FOR REJECTION FROM POOL'
        endif
        call img%kill
    end subroutine reject_from_pool

    subroutine reject_from_pool_user
        type(image)          :: img
        type(oris)           :: cls2reject
        logical, allocatable :: cls_mask(:)
        integer              :: nptcls_rejected, ncls_rejected, iptcl
        integer              :: icls, jcls, i, nl
        if( .not.pool_available )return
        if( pool_proj%os_cls2D%get_noris() == 0 ) return
        if( .not.file_exists(STREAM_REJECT_CLS) ) return
        nl = nlines(STREAM_REJECT_CLS)
        if( nl == 0 ) return
        call cls2reject%new(nl,is_ptcl=.false.)
        call cls2reject%read(STREAM_REJECT_CLS)
        call del_file(STREAM_REJECT_CLS)
        if( cls2reject%get_noris() == 0 ) return
        allocate(cls_mask(ncls_glob), source=.false.)
        do i = 1,cls2reject%get_noris()
            icls = nint(cls2reject%get(i,'class'))
            if( icls == 0 ) cycle
            if( icls > ncls_glob ) cycle
            if( nint(pool_proj%os_cls2D%get(icls,'pop')) == 0 ) cycle
            cls_mask(icls) = .true.
        enddo
        call cls2reject%kill
        if( count(cls_mask) == 0 ) return
        call img%new([box,box,1],smpd)
        img = 0.
        ncls_rejected   = 0
        do icls = 1,ncls_glob
            if( .not.cls_mask(icls) ) cycle
            nptcls_rejected = 0
            do iptcl=1,pool_proj%os_ptcl2D%get_noris()
                if( pool_proj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                jcls = nint(pool_proj%os_ptcl2D%get(iptcl,'class'))
                if( jcls == icls )then
                    call pool_proj%os_ptcl2D%reject(iptcl)
                    nptcls_rejected = nptcls_rejected + 1
                endif
            enddo
            if( nptcls_rejected > 0 )then
                ncls_rejected = ncls_rejected + 1
                call pool_proj%os_cls2D%set(icls,'pop',0.)
                call pool_proj%os_cls2D%set(icls,'corr',-1.)
                call img%write(trim(POOL_DIR)//trim(refs_glob),icls)
                if( l_wfilt )then
                    call img%write(trim(POOL_DIR)//trim(add2fbody(refs_glob, params_glob%ext,trim(WFILT_SUFFIX))), icls)
                endif
                write(logfhandle,'(A,I6,A,I4)')'>>> USER REJECTED FROM POOL: ',nptcls_rejected,' PARTICLE(S) IN CLASS ',icls
            endif

        enddo
        call img%read(trim(POOL_DIR)//trim(refs_glob), ncls_glob)
        call img%write(trim(POOL_DIR)//trim(refs_glob), ncls_glob)
        if( l_wfilt )then
            call img%read(trim(POOL_DIR)//trim(add2fbody(refs_glob, params_glob%ext,trim(WFILT_SUFFIX))), ncls_glob)
            call img%write(trim(POOL_DIR)//trim(add2fbody(refs_glob, params_glob%ext,trim(WFILT_SUFFIX))), ncls_glob)
        endif
        call img%kill
    end subroutine reject_from_pool_user

    !> updates current parameters with user input
    subroutine update_user_params
        type(oris) :: os
        real       :: lpthres, ndev
        call os%new(1, is_ptcl=.false.)
        if( file_exists(USER_PARAMS) )then
            call os%read(USER_PARAMS)
            if( os%isthere(1,'lpthres') )then
                lpthres = os%get(1,'lpthres')
                if( abs(lpthres-params_glob%lpthres) > 0.001 )then
                    params_glob%lpthres = lpthres
                    write(logfhandle,'(A,F8.2)')'>>> REJECTION lpthres UPDATED TO: ',params_glob%lpthres
                endif
            endif
            if( os%isthere(1,'ndev2D') )then
                ndev = os%get(1,'ndev2D')
                if( abs(ndev-params_glob%ndev2D) > 0.001 )then
                    params_glob%ndev2D = ndev
                    write(logfhandle,'(A,F8.2)')'>>> REJECTION NDEV2D   UPDATED TO: ',params_glob%ndev2D
                endif
            endif
            call del_file(USER_PARAMS)
        endif
        call os%kill
    end subroutine update_user_params

    subroutine classify_pool
        use simple_ran_tabu
        type(ran_tabu)          :: random_generator
        type(sp_project)        :: spproj
        logical, parameter      :: L_BENCH = .false.
        integer, allocatable    :: prev_eo_pops(:,:), min_update_cnts_per_stk(:), nptcls_per_stk(:), stk_order(:)
        ! character(len=XLONGSTRLEN) :: cwd
        real                    :: frac_update
        integer                 :: iptcl,i, nptcls_tot, nptcls_old, fromp, top, nstks_tot, jptcl
        integer                 :: eo, icls, nptcls_sel, istk, nptcls2update, nstks2update
        integer(timer_int_kind) :: t_tot
        if( .not.pool_available )return
        if( L_BENCH ) t_tot  = tic()
        nptcls_tot = pool_proj%os_ptcl2D%get_noris()
        if( nptcls_tot == 0 ) return
        pool_iter = pool_iter + 1
        call cline_cluster2D_pool%set('refs',    refs_glob)
        call cline_cluster2D_pool%set('ncls',    real(ncls_glob))
        call cline_cluster2D_pool%set('startit', real(pool_iter))
        call cline_cluster2D_pool%set('maxits',  real(pool_iter))
        call cline_cluster2D_pool%set('frcs',    trim(FRCS_FILE))
        spproj%projinfo = pool_proj%projinfo
        spproj%compenv  = pool_proj%compenv
        call spproj%projinfo%delete_entry('projname')
        call spproj%projinfo%delete_entry('projfile')
        call spproj%update_projinfo( cline_cluster2D_pool )
        ! counting number of stacks & particles (old and newer)
        nstks_tot  = pool_proj%os_stk%get_noris()
        allocate(nptcls_per_stk(nstks_tot), min_update_cnts_per_stk(nstks_tot), source=0)
        nptcls_old = 0
        do istk = 1,nstks_tot
            fromp = nint(pool_proj%os_stk%get(istk,'fromp'))
            top   = nint(pool_proj%os_stk%get(istk,'top'))
            min_update_cnts_per_stk(istk) = huge(istk)
            do iptcl = fromp,top
                if( pool_proj%os_ptcl2D%get_state(iptcl) > 0 )then
                    nptcls_per_stk(istk)          = nptcls_per_stk(istk) + 1 ! # ptcls with state=1
                    min_update_cnts_per_stk(istk) = min(min_update_cnts_per_stk(istk), nint(pool_proj%os_ptcl2D%get(iptcl,'updatecnt')))
                endif
            enddo
            if( min_update_cnts_per_stk(istk) >= STREAM_SRCHLIM )then
                nptcls_old = nptcls_old + nptcls_per_stk(istk)
            endif
        enddo
        ! flagging stacks to be skipped
        if( allocated(pool_stacks_mask) ) deallocate(pool_stacks_mask)
        allocate(pool_stacks_mask(nstks_tot), source=.false.)
        allocate(prev_eo_pops(ncls_glob,2),source=0)
        if( nptcls_old > params_glob%ncls_start*params_glob%nptcls_per_cls )then
            allocate(stk_order(nstks_tot))
            do istk = 1,nstks_tot
                stk_order(istk) = istk
            enddo
            random_generator = ran_tabu(nstks_tot)
            call random_generator%shuffle(stk_order)
            nptcls2update = 0 ! # of ptcls including state=0 within selected stacks
            nptcls_sel    = 0 ! # of ptcls excluding state=0 within selected stacks
            do i = 1,nstks_tot
                istk = stk_order(i)
                if( (min_update_cnts_per_stk(istk) > STREAM_SRCHLIM) .and. (nptcls_sel > MAX_STREAM_NPTCLS) ) cycle
                nptcls_sel    = nptcls_sel + nptcls_per_stk(istk)
                nptcls2update = nptcls2update + nint(pool_proj%os_stk%get(istk,'nptcls'))
                pool_stacks_mask(istk) = .true.
            enddo
            call random_generator%kill
        else
            nptcls2update = nptcls_tot
            nptcls_sel    = sum(nptcls_per_stk)
            do istk = 1,nstks_tot
                pool_stacks_mask(istk) = nptcls_per_stk(istk) > 0
            enddo
        endif
        nstks2update = count(pool_stacks_mask)
        ! transfer stacks and particles
        call spproj%os_stk%new(nstks2update, is_ptcl=.false.)
        call spproj%os_ptcl2D%new(nptcls2update, is_ptcl=.true.)
        i     = 0
        jptcl = 0
        do istk = 1,nstks_tot
            fromp = nint(pool_proj%os_stk%get(istk,'fromp'))
            top   = nint(pool_proj%os_stk%get(istk,'top'))
            if( pool_stacks_mask(istk) )then
                ! transfer alignement parameters for selected particles
                i = i + 1 ! stack index in spproj
                call spproj%os_stk%transfer_ori(i, pool_proj%os_stk, istk)
                call spproj%os_stk%set(i, 'fromp', real(jptcl+1))
                do iptcl = fromp,top
                    jptcl = jptcl+1
                    call spproj%os_ptcl2D%transfer_ori(jptcl, pool_proj%os_ptcl2D, iptcl)
                    call spproj%os_ptcl2D%set(jptcl, 'stkind', real(i))
                enddo
                call spproj%os_stk%set(i, 'top', real(jptcl))
            else
                ! keeps track of skipped particles
                do iptcl = fromp,top
                    icls = pool_proj%os_ptcl2D%get_class(iptcl)
                    eo   = nint(pool_proj%os_ptcl2D%get(iptcl,'eo')) + 1
                    prev_eo_pops(icls,eo) = prev_eo_pops(icls,eo) + 1
                enddo
            endif
        enddo
        call spproj%os_ptcl3D%new(nptcls2update, is_ptcl=.true.)
        spproj%os_cls2D = pool_proj%os_cls2D
        ! update command line and write project
        if( sum(prev_eo_pops) == 0 )then
            call cline_cluster2D_pool%delete('update_frac')
        else
            frac_update = real(nptcls_old-sum(prev_eo_pops)) / real(nptcls_old)
            call cline_cluster2D_pool%set('update_frac', frac_update)
            call cline_cluster2D_pool%set('center',      'no')
            do icls = 1,ncls_glob
                call spproj%os_cls2D%set(icls,'prev_pop_even',real(prev_eo_pops(icls,1)))
                call spproj%os_cls2D%set(icls,'prev_pop_odd', real(prev_eo_pops(icls,2)))
            enddo
        endif
        call spproj%write(trim(POOL_DIR)//trim(PROJFILE_POOL))
        call spproj%kill
        ! execution in correct directory
        call qenv_pool%exec_simple_prg_in_queue_async(cline_cluster2D_pool, DISTR_EXEC_FNAME, 'simple_log_cluster2D_pool')
        pool_available = .false.
        write(logfhandle,'(A,I6,A,I8,A3,I8,A)')'>>> POOL         INITIATED ITERATION ',pool_iter,' WITH ',nptcls_sel,&
        &' / ', sum(nptcls_per_stk),' PARTICLES'
        if( L_BENCH ) print *,'timer exec_classify_pool tot : ',toc(t_tot)
        call tidy_2Dstream_iter
    end subroutine classify_pool

    !> produces consolidated project at original scale
    subroutine write_project_stream2D( snapshot )
        logical,           intent(in) :: snapshot
        type(class_frcs)              :: frcs, frcs_sc
        type(oris)                    :: os_backup2, os_backup3
        character(len=:), allocatable :: projfile,projfname, cavgsfname, frcsfname, src, dest
        character(len=:), allocatable :: pool_refs
        integer :: istk
        logical :: do_write
        do_write = .not.snapshot
        if( snapshot )then
            ! writes snapshot every ORIGPROJ_WRITEFREQ seconds
            if( .not.pool_available )return
            if( (simple_gettime()-origproj_time) > ORIGPROJ_WRITEFREQ .and. pool_iter>1 )then
                do_write      = .true.
                origproj_time = simple_gettime()
            endif
        endif
        if( .not.do_write ) return
        ! file naming
        projfname  = get_fbody(orig_projfile, METADATA_EXT, separator=.false.)
        cavgsfname = get_fbody(refs_glob, params_glob%ext, separator=.false.)
        frcsfname  = get_fbody(FRCS_FILE, BIN_EXT, separator=.false.)
        call pool_proj%projinfo%set(1,'projname', projfname)
        projfile   = trim(projfname)//trim(METADATA_EXT)
        cavgsfname = trim(cavgsfname)//trim(params_glob%ext)
        frcsfname  = trim(frcsfname)//trim(BIN_EXT)
        if( snapshot )then
            cavgsfname = trim(SNAPSHOT_DIR)//trim(cavgsfname)
            frcsfname  = trim(SNAPSHOT_DIR)//trim(frcsfname)
        endif
        call pool_proj%projinfo%set(1,'projfile', projfile)
        ! removing previous snapshot
        if( (trim(prev_snapshot_cavgs) /= '') .and. snapshot)then
            call del_file(prev_snapshot_frcs)
            call del_file(prev_snapshot_cavgs)
            src = add2fbody(prev_snapshot_cavgs, params_glob%ext,'_even')
            call del_file(src)
            src = add2fbody(prev_snapshot_cavgs, params_glob%ext,'_odd')
            call del_file(src)
            if( l_wfilt )then
                src = add2fbody(prev_snapshot_cavgs, params_glob%ext, trim(WFILT_SUFFIX))
                call del_file(src)
                src = add2fbody(prev_snapshot_cavgs, params_glob%ext, trim(WFILT_SUFFIX)//'_even')
                call del_file(src)
                src = add2fbody(prev_snapshot_cavgs, params_glob%ext, trim(WFILT_SUFFIX)//'_odd')
                call del_file(src)
            endif
        endif
        write(logfhandle,'(A,A,A,A)')'>>> GENERATING PROJECT SNAPSHOT ',trim(projfile), ' AT: ',cast_time_char(simple_gettime())
        pool_refs = trim(POOL_DIR)//trim(refs_glob)
        if( do_autoscale )then
            os_backup3 = pool_proj%os_cls2D
            os_backup2 = pool_proj%os_stk
            ! rescale classes
            if( l_wfilt )then
                src  = add2fbody(pool_refs, params_glob%ext,trim(WFILT_SUFFIX))
                dest = add2fbody(cavgsfname,params_glob%ext,trim(WFILT_SUFFIX))
                call rescale_cavgs(src, dest)
                src  = add2fbody(pool_refs, params_glob%ext,trim(WFILT_SUFFIX)//'_even')
                dest = add2fbody(cavgsfname,params_glob%ext,trim(WFILT_SUFFIX)//'_even')
                call rescale_cavgs(src, dest)
                src  = add2fbody(pool_refs, params_glob%ext,trim(WFILT_SUFFIX)//'_odd')
                dest = add2fbody(cavgsfname,params_glob%ext,trim(WFILT_SUFFIX)//'_odd')
                call rescale_cavgs(src, dest)
            endif
            call rescale_cavgs(pool_refs, cavgsfname)
            src  = add2fbody(pool_refs, params_glob%ext,'_even')
            dest = add2fbody(cavgsfname,params_glob%ext,'_even')
            call rescale_cavgs(src, dest)
            src  = add2fbody(pool_refs, params_glob%ext,'_odd')
            dest = add2fbody(cavgsfname,params_glob%ext,'_odd')
            call rescale_cavgs(src, dest)
            call pool_proj%os_out%kill
            call pool_proj%add_cavgs2os_out(cavgsfname, orig_smpd, 'cavg')
            if( l_wfilt )then
                src = add2fbody(cavgsfname,params_glob%ext,trim(WFILT_SUFFIX))
                call pool_proj%add_cavgs2os_out(src, orig_smpd, 'cavg'//trim(WFILT_SUFFIX))
            endif
            pool_proj%os_cls2D = os_backup3
            call os_backup3%kill
            ! rescale frcs
            call frcs_sc%read(trim(POOL_DIR)//trim(FRCS_FILE))
            call frcs_sc%upsample(orig_smpd, orig_box, frcs)
            call frcs%write(frcsfname)
            call frcs%kill
            call frcs_sc%kill
            call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
            ! project updates to original scale
            call pool_proj%os_stk%set_all2single('box', real(orig_box))
            call pool_proj%os_stk%set_all2single('smpd',orig_smpd)
            call pool_proj%os_ptcl2D%mul_shifts( 1./scale_factor )
            do istk = 1,pool_proj%os_stk%get_noris()
                call pool_proj%os_stk%set(istk,'stk',imported_stks(istk))
            enddo
            ! write
            pool_proj%os_ptcl3D = pool_proj%os_ptcl2D
            call pool_proj%os_ptcl3D%delete_2Dclustering
            call pool_proj%write(projfile)
            call pool_proj%os_ptcl3D%kill
            ! preserve down-scaling
            call pool_proj%os_ptcl2D%mul_shifts( scale_factor )
            pool_proj%os_stk = os_backup2
            call os_backup2%kill
        else
            call pool_proj%os_out%kill
            call pool_proj%add_cavgs2os_out(cavgsfname, orig_smpd, 'cavg')
            if( l_wfilt )then
                src = add2fbody(cavgsfname,params_glob%ext,trim(WFILT_SUFFIX))
                call pool_proj%add_cavgs2os_out(src, orig_smpd, 'cavg'//trim(WFILT_SUFFIX))
            endif
            call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
            ! write
            pool_proj%os_ptcl3D = pool_proj%os_ptcl2D
            call pool_proj%os_ptcl3D%delete_2Dclustering
            call pool_proj%write(projfile)
            call pool_proj%os_ptcl3D%kill
        endif
        if( snapshot )then
            prev_snapshot_frcs  = trim(frcsfname)
            prev_snapshot_cavgs = trim(cavgsfname)
        endif
        ! classes export
        if(file_exists('clusters2D.star')) call del_file('clusters2D.star')
        call starproj%export_cls2D(pool_proj)
        call pool_proj%os_cls2D%delete_entry('stk')
    end subroutine write_project_stream2D

    subroutine terminate_stream2D( write_project )
        logical,        intent(in) :: write_project
        integer :: ichunk, ipart
        do ichunk = 1,params_glob%nchunks
            call chunks(ichunk)%terminate
        enddo
        if( .not.pool_available )then
            pool_iter = pool_iter-1
            refs_glob = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter,3))//trim(params_glob%ext)
            ! tricking the asynchronous master process to come to a hard stop
            call simple_touch(trim(POOL_DIR)//trim(TERM_STREAM))
            do ipart = 1,params_glob%nparts_pool
                call simple_touch(trim(POOL_DIR)//trim(JOB_FINISHED_FBODY)//int2str_pad(ipart,numlen))
            enddo
            call simple_touch(trim(POOL_DIR)//'CAVGASSEMBLE_FINISHED')
        endif
        if( pool_iter >= 1 .and. write_project )then
            ! updates project
            call write_project_stream2D(.false.)
            ! ranking
            call rank_cavgs
        endif
        ! cleanup
        call simple_rmdir(SCALE_DIR)
        if( .not.debug_here )then
            ! call qsys_cleanup
            ! call simple_rmdir(POOL_DIR)
        endif
    end subroutine terminate_stream2D

    ! Utilities

    ! For rankging class-averages
    subroutine rank_cavgs
        type(rank_cavgs_commander) :: xrank_cavgs
        type(cmdline)              :: cline_rank_cavgs
        character(len=STDLEN)      :: refs_ranked
        refs_ranked = add2fbody(refs_glob, params_glob%ext ,'_ranked')
        call cline_rank_cavgs%set('projfile', orig_projfile)
        call cline_rank_cavgs%set('stk',      trim(POOL_DIR)//refs_glob)
        call cline_rank_cavgs%set('outstk',   trim(refs_ranked))
        call xrank_cavgs%execute(cline_rank_cavgs)
        call cline_rank_cavgs%kill
    end subroutine rank_cavgs

    subroutine update_path(os, i, key)
        class(oris),      intent(inout) :: os
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        character(len=:), allocatable :: fname
        character(len=LONGSTRLEN)     :: newfname
        if( os%isthere(i,key) )then
            call os%getter(i,key,fname)
            if( len_trim(fname) > 3 )then
                if( fname(1:3) == '../' ) fname = trim(fname(4:))
            endif
            call make_relativepath(CWD_GLOB, fname, newfname)
            call os%set(i, key, newfname)
        endif
    end subroutine update_path

    !>  Updates the project watched by the gui for display
    subroutine update_pool_for_gui
        type(oris)                    :: os_backup
        type(starproject)             :: starproj
        character(len=:), allocatable :: src
        os_backup = pool_proj%os_cls2D
        call pool_proj%add_cavgs2os_out(trim(POOL_DIR)//trim(refs_glob), smpd, 'cavg')
        if( l_wfilt )then
            src = trim(POOL_DIR)//add2fbody(refs_glob,params_glob%ext,trim(WFILT_SUFFIX))
            call pool_proj%add_cavgs2os_out(src, smpd, 'cavg'//trim(WFILT_SUFFIX))
        endif
        pool_proj%os_cls2D = os_backup
        call pool_proj%write_segment_inside('out',   projfile4gui)
        call pool_proj%write_segment_inside('cls2D', projfile4gui)
        ! Write star file for iteration
        call starproj%export_cls2D(pool_proj, pool_iter)
        call pool_proj%os_cls2D%delete_entry('stk')
        call os_backup%kill
    end subroutine update_pool_for_gui

    !>  Convenience function
    subroutine transfer_cavg( refs_in, dir, indin, refs_out, indout, self_transfer )
        character(len=*),  intent(in) :: refs_in, dir, refs_out
        integer,           intent(in) :: indin, indout
        logical, optional, intent(in) :: self_transfer
        type(image)                   :: img
        character(len=:), allocatable :: stkout, stkin, refs_out_here, refs_in_here
        integer :: ipart
        logical :: l_self
        call debug_print('in transfer_cavg '//int2str(indin)//' '//int2str(indout))
        l_self = .false.
        if( present(self_transfer) ) l_self = self_transfer
        if( l_self )then
            call debug_print('in transfer_cavg self_transfer')
            refs_in_here = trim(POOL_DIR)//trim(refs_in)
        else
            refs_in_here = trim(refs_in)
        endif
        ! making sure we are writing to the correct folder
        refs_out_here = trim(POOL_DIR)//trim(refs_out)
        ! merged class
        call img%new([box,box,1],smpd)
        call img%read( refs_in_here, indin)
        call img%write(refs_out_here,indout)
        if( l_wfilt )then
            stkout = add2fbody(refs_out_here,params_glob%ext,trim(WFILT_SUFFIX))
            call img%write(stkout,indout)
        endif
        ! e/o
        stkin  = add2fbody(refs_in_here, params_glob%ext,'_even')
        stkout = add2fbody(refs_out_here,params_glob%ext,'_even')
        call img%read( stkin, indin)
        call img%write(stkout,indout)
        stkin  = add2fbody(refs_in_here,params_glob%ext,'_odd')
        stkout = add2fbody(refs_out_here,params_glob%ext,'_odd')
        call img%read( stkin, indin)
        call img%write(stkout,indout)
        ! temporary matrices, logics from chunk%read
        call img%new([boxpd,boxpd,1],smpd)
        call img%zero_and_flag_ft
        if( l_self )then
            do ipart = 1,params_glob%nparts_pool
                stkin = trim(POOL_DIR)//'cavgs_even_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                call img%read(stkin, indin)
                call img%write(stkin,indout)
                stkin = trim(POOL_DIR)//'cavgs_odd_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                call img%read(stkin, indin)
                call img%write(stkin,indout)
                stkin = trim(POOL_DIR)//'ctfsqsums_even_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                call img%read(stkin, indin)
                call img%write(stkin,indout)
                stkin = trim(POOL_DIR)//'ctfsqsums_odd_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                call img%read(stkin, indin)
                call img%write(stkin,indout)
                if( l_wfilt )then
                    stkin = trim(POOL_DIR)//'cavgs_even_wfilt_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                    call img%read(stkin, indin)
                    call img%write(stkin,indout)
                    stkin = trim(POOL_DIR)//'cavgs_odd_wfilt_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                    call img%read(stkin, indin)
                    call img%write(stkin,indout)
                    stkin = trim(POOL_DIR)//'ctfsqsums_even_wfilt_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                    call img%read(stkin, indin)
                    call img%write(stkin,indout)
                    stkin = trim(POOL_DIR)//'ctfsqsums_odd_wfilt_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                    call img%read(stkin, indin)
                    call img%write(stkin,indout)
                endif
            enddo
        else
            stkin  = trim(dir)//'/cavgs_even_part'//trim(params_glob%ext)
            call img%read(stkin, indin)
            do ipart = 1,params_glob%nparts_pool
                stkout = trim(POOL_DIR)//'cavgs_even_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                call img%write(stkout,indout)
                if( l_wfilt )then
                    stkout = trim(POOL_DIR)//'cavgs_even_wfilt_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                    call img%write(stkout,indout)
                endif
            enddo
            stkin  = trim(dir)//'/cavgs_odd_part'//trim(params_glob%ext)
            call img%read(stkin, indin)
            do ipart = 1,params_glob%nparts_pool
                stkout = trim(POOL_DIR)//'cavgs_odd_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                call img%write(stkout,indout)
                if( l_wfilt )then
                    stkout = trim(POOL_DIR)//'cavgs_odd_wfilt_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                    call img%write(stkout,indout)
                endif
            enddo
            stkin  = trim(dir)//'/ctfsqsums_even_part'//trim(params_glob%ext)
            call img%read(stkin, indin)
            do ipart = 1,params_glob%nparts_pool
                stkout = trim(POOL_DIR)//'ctfsqsums_even_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                call img%write(stkout,indout)
                if( l_wfilt )then
                    stkout = trim(POOL_DIR)//'ctfsqsums_even_wfilt_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                    call img%write(stkout,indout)
                endif
            enddo
            stkin  = trim(dir)//'/ctfsqsums_odd_part'//trim(params_glob%ext)
            call img%read(stkin, indin)
            do ipart = 1,params_glob%nparts_pool
                stkout = trim(POOL_DIR)//'ctfsqsums_odd_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                call img%write(stkout,indout)
                if( l_wfilt )then
                    stkout = trim(POOL_DIR)//'ctfsqsums_odd_wfilt_part'//int2str_pad(ipart,numlen)//trim(params_glob%ext)
                    call img%write(stkout,indout)
                endif
            enddo
        endif
        ! cleanup
        call img%kill
    end subroutine transfer_cavg

    subroutine rescale_cavgs( src, dest )
        character(len=*), intent(in) :: src, dest
        type(image)          :: img, img_pad
        type(stack_io)       :: stkio_r, stkio_w
        character(len=:), allocatable :: dest_here
        integer, allocatable :: cls_pop(:)
        integer              :: ldim(3),icls, ncls_here
        call debug_print('in rescale_cavgs '//trim(src))
        call debug_print('in rescale_cavgs '//trim(dest))
        if(trim(src) == trim(dest))then
            dest_here = 'tmp_cavgs.mrc'
        else
            dest_here = trim(dest)
        endif
        call img%new([box,box,1],smpd)
        call img_pad%new([orig_box,orig_box,1],orig_smpd)
        cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
        call find_ldim_nptcls(src,ldim,ncls_here)
        call stkio_r%open(trim(src), smpd, 'read', bufsz=ncls_here)
        call stkio_r%read_whole
        call stkio_w%open(dest_here, orig_smpd, 'write', box=orig_box, bufsz=ncls_here)
        do icls = 1,ncls_here
            if( cls_pop(icls) > 0 )then
                call img%zero_and_unflag_ft
                call stkio_r%get_image(icls, img)
                call img%fft
                call img%pad(img_pad, backgr=0.)
                call img_pad%ifft
            else
                img_pad = 0.
            endif
            call stkio_w%write(icls, img_pad)
        enddo
        call stkio_r%close
        call stkio_w%close
        if (trim(src) == trim(dest) ) call simple_rename('tmp_cavgs.mrc',dest)
        call img%kill
        call img_pad%kill
        call debug_print('end rescale_cavgs')
    end subroutine rescale_cavgs

    subroutine tidy_2Dstream_iter
        character(len=:), allocatable :: prefix
        if( pool_iter > 3 )then
            prefix = trim(POOL_DIR)//trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter-3,3))
            call del_file(prefix//'_even'//trim(params_glob%ext))
            call del_file(prefix//'_odd'//trim(params_glob%ext))
            if( l_wfilt )then
                call del_file(prefix//trim(WFILT_SUFFIX)//'_even'//trim(params_glob%ext))
                call del_file(prefix//trim(WFILT_SUFFIX)//'_odd'//trim(params_glob%ext))
            endif
        endif
    end subroutine tidy_2Dstream_iter

    subroutine debug_print( string )
        character(len=*), intent(in) :: string
        if( DEBUG_HERE )then
            write(logfhandle,*) trim(string)
            call flush(logfhandle)
        endif
    end subroutine debug_print

    !
    ! cluster2D_subsets (=cluster2d_stream offline)
    !
    subroutine exec_cluster2D_subsets( self, cline )
        use simple_class_frcs,        only: class_frcs
        class(cluster2D_commander_subsets), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        character(len=STDLEN),    parameter :: dir_preprocess = trim(PATH_HERE)//'spprojs/'
        integer,                  parameter :: WAITTIME       = 3
        real                                :: SMPD_TARGET    = MAX_SMPD  ! target sampling distance
        type(parameters)                       :: params
        type(sp_project)                       :: spproj
        type(class_frcs)                       :: frcs, frcs_sc
        type(oris)                             :: os_mic
        character(len=LONGSTRLEN), allocatable :: completed_fnames(:), tmp_fnames(:), pool_stacks(:)
        character(len=:),          allocatable :: one_projfile, one_projname, fname, stack_fname, orig_stack_fname, carg
        character(len=:),          allocatable :: frcsfname, src
        integer,                   allocatable :: stk_nptcls(:), stk_all_nptcls(:), chunks_map(:,:)
        logical,                   allocatable :: pool_stk_mask(:)
        integer :: ichunk, istk, imic, nstks, nptcls, nptcls_tot, ntot_chunks, cnt, ic, fromp, top, nsplit
        integer :: maxits, pool_nstks, iptcl, jptcl, jstk, nchunks_imported, tot_nchunks_imported
        logical :: all_chunks_submitted, all_chunks_imported, l_once
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')        ) call cline%set('cenlp',      30.0)
        if( .not. cline%defined('center')       ) call cline%set('center',     'yes')
        if( .not. cline%defined('autoscale')    ) call cline%set('autoscale',  'yes')
        if( .not. cline%defined('lp')           ) call cline%set('lp',          GREEDY_TARGET_LP)
        if( .not. cline%defined('lpthres')     ) call cline%set('lpthres',    30.0)
        if( .not. cline%defined('ndev')         ) call cline%set('ndev',        1.5)
        if( .not. cline%defined('oritype')      ) call cline%set('oritype',     'ptcl2D')
        if( .not. cline%defined('wiener')       ) call cline%set('wiener',      'partial')
        if( .not. cline%defined('walltime')     ) call cline%set('walltime',     29.0*60.0) ! 29 minutes
        if( .not. cline%defined('nparts_chunk') ) call cline%set('nparts_chunk', 1.0)
        if( .not. cline%defined('nchunks')      ) call cline%set('nchunks',      2.0)
        if( .not. cline%defined('numlen')       ) call cline%set('numlen',       5.0)
        if( .not. cline%defined('match_filt')   ) call cline%set('match_filt',    'no')
        if( .not. cline%defined('nonuniform')   ) call cline%set('nonuniform',    'no')
        call cline%set('ptclw',      'no')
        call seed_rnd
        call params%new(cline)
        l_wfilt = trim(params%wiener) .eq. 'partial'
        ! sanity
        if( .not.file_exists(params%projfile) )then
            THROW_HARD('project file: '//trim(params%projfile)//' does not exist!')
        endif
        orig_projfile = trim(params%projfile)
        call cline%set('mkdir','no')
        ! init
        do_autoscale      = params%autoscale .eq. 'yes'
        max_ncls          = floor(real(params%ncls)/real(params%ncls_start))*params%ncls_start ! effective maximum # of classes
        nptcls_per_chunk  = params%nptcls_per_cls*params%ncls_start         ! # of particles in each chunk
        ncls_glob         = 0
        l_greedy          = trim(params%refine).eq.'greedy'
        lp_greedy         = GREEDY_TARGET_LP
        if( cline%defined('lp') ) lp_greedy = params%lp
        lpstart_stoch     = lp_greedy
        prev_snapshot_cavgs = ''
        ! init command-lines
        call cline%delete('lp')
        call cline%delete('refine')
        call cline%delete('nptcls_per_cls')
        call cline%delete('ncls_start')
        call cline%delete('numlen')
        cline_cluster2D_pool   = cline
        cline_cluster2D_chunk  = cline
        ! chunk classification
        if( params%nparts_chunk > 1 )then
            call cline_cluster2D_chunk%set('prg',    'cluster2D_distr')
            call cline_cluster2D_chunk%set('nparts', real(params%nparts_chunk))
        else
            ! shared memory execution
            call cline_cluster2D_chunk%set('prg','cluster2D')
        endif
        call cline_cluster2D_chunk%delete('projfile')
        call cline_cluster2D_chunk%delete('projname')
        call cline_cluster2D_chunk%set('objfun',    'cc')
        call cline_cluster2D_chunk%set('center',    'no')
        call cline_cluster2D_chunk%set('match_filt','no')
        call cline_cluster2D_chunk%set('autoscale', 'no')
        call cline_cluster2D_chunk%set('ptclw',     'no')
        call cline_cluster2D_chunk%set('mkdir',     'no')
        call cline_cluster2D_chunk%set('stream',    'no')
        call cline_cluster2D_chunk%set('startit',   1.)
        call cline_cluster2D_chunk%set('ncls',      real(params%ncls_start))
        if( l_wfilt ) call cline_cluster2D_chunk%set('wiener', 'partial')
        call cline_cluster2D_chunk%delete('update_frac')
        call cline_cluster2D_chunk%delete('dir_target')
        call cline_cluster2D_chunk%delete('lpstop')
        ! pool classification: optional stochastic optimisation, optional match filter
        ! automated incremental learning, objective function is standard cross-correlation (cc)
        call cline_cluster2D_pool%set('prg',       'cluster2D_distr')
        call cline_cluster2D_pool%set('autoscale', 'no')
        call cline_cluster2D_pool%set('trs',       MINSHIFT)
        call cline_cluster2D_pool%set('projfile',  trim(PROJFILE_POOL))
        call cline_cluster2D_pool%set('projname',  trim(get_fbody(trim(PROJFILE_POOL),trim('simple'))))
        call cline_cluster2D_pool%set('objfun',    'cc')
        call cline_cluster2D_pool%set('ptclw',     'no')
        call cline_cluster2D_pool%set('extr_iter', 100.)
        call cline_cluster2D_pool%set('mkdir',     'no')
        call cline_cluster2D_pool%set('async',     'yes') ! to enable hard termination
        call cline_cluster2D_pool%set('stream',    'yes') ! use for dual CTF treatment
        call cline_cluster2D_pool%set('nparts',    real(params%nparts_pool))
        if( cline%defined('match_filt') )then
            carg = cline%get_carg('match_filt')
            call cline_cluster2D_pool%set('match_filt',carg)
            deallocate(carg)
        else
            call cline_cluster2D_pool%set('match_filt','no')
        endif
        if( l_wfilt ) call cline_cluster2D_pool%set('wiener', 'partial')
        call cline_cluster2D_pool%delete('lpstop')
        ! read strictly required fields project
        call spproj%read_non_data_segments(params%projfile)
        call spproj%read_segment('mic',   params%projfile)
        call spproj%read_segment('stk',   params%projfile)
        call spproj%read_segment('ptcl2D',params%projfile)
        nstks = spproj%os_stk%get_noris()
        ! some selection has been performed beforehand
        if( spproj%os_mic%get_noris() > nstks )then
            os_mic = spproj%os_mic
            call spproj%os_mic%kill
            call spproj%os_mic%new(nstks, is_ptcl=.false.)
            cnt = 0
            do imic = 1,os_mic%get_noris()
                if( (os_mic%get_state(imic) == 0) .or. (os_mic%get(imic,'nptcls') < 0.5) ) cycle
                cnt = cnt + 1
                call spproj%os_mic%transfer_ori(cnt,os_mic,imic)
            enddo
        endif
        ! sanity checks
        nptcls = spproj%get_nptcls()
        if( nptcls == 0 )then
            THROW_HARD('No particles found in project file: '//trim(params%projfile)//'; exec_cluster2d_subsets')
        endif
        if ( nptcls < 2*nptcls_per_chunk )then
            THROW_WARN('Not enough particles to classify more than one chunk')
            THROW_HARD('Review parameters or use cleanup2D/cluster2D instead')
        endif
        ! splitting
        if( spproj%os_stk%get_noris() == 1 )then
            spproj%os_ptcl3D = spproj%os_ptcl2D
            nsplit = floor(real(spproj%os_ptcl2D%get_noris())/real(nptcls_per_chunk))
            call spproj%split_stk(nsplit, dir=PATH_PARENT)
            call spproj%os_ptcl3D%kill
        endif
        ! directory structure & temporary project handling
        projfile4gui   = trim(orig_projfile)
        numlen         = len(int2str(params%nparts_pool))
        refs_glob      = 'start_cavgs'//params%ext
        pool_available = .true.
        pool_iter      = 0
        call simple_mkdir(SCALE_DIR)
        if( trim(POOL_DIR) /= '' ) call simple_mkdir(POOL_DIR)
        call simple_mkdir(trim(POOL_DIR)//trim(STDERROUT_DIR))
        call simple_mkdir(dir_preprocess)
        call simple_touch(trim(POOL_DIR)//trim(CLUSTER2D_FINISHED))
        call pool_proj%kill
        pool_proj%projinfo = spproj%projinfo
        pool_proj%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()>0 ) pool_proj%jobproc = spproj%jobproc
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        call pool_proj%write(trim(POOL_DIR)//PROJFILE_POOL)
        ! auto-scaling
        orig_box  = params%box
        orig_smpd = params%smpd
        if( orig_box == 0 ) THROW_HARD('FATAL ERROR')
        if( do_autoscale )then
            if( orig_box < MINBOXSZ )then
                do_autoscale = .false.
            else
                call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                if( box < MINBOXSZ )then
                    SMPD_TARGET = orig_smpd * real(orig_box) / real(MINBOXSZ)
                    call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                endif
                if( box >= orig_box )then
                    do_autoscale = .false.
                else
                    write(logfhandle,'(A,I5,F5.2)') '>>> 2D CLASSIFICATION DOWNSCALED IMAGE SIZE & PIXEL SIZE (IN A): ', box, smpd
                endif
            endif
        endif
        if( .not. do_autoscale )then
            smpd         = orig_smpd
            box          = orig_box
            scale_factor = 1.
        endif
        boxpd = 2 * round2even(params%alpha * real(box/2)) ! logics from parameters
        call cline_cluster2D_chunk%set('box',  real(box))
        call cline_cluster2D_chunk%set('smpd', smpd)
        call cline_cluster2D_pool%set('box',   real(box))
        call cline_cluster2D_pool%set('smpd',  smpd)
        ! resolution-related updates to command-lines
        lp_greedy     = max(lp_greedy,    2.0*smpd)
        lpstart_stoch = max(lpstart_stoch,2.0*smpd)
        if( cline%defined('lpstop2D') )then
            params%lpstop2D = max(2.0*smpd,params%lpstop2D)
        else
            params%lpstop2D = 2.0*smpd
        endif
        call cline_cluster2D_chunk%set('lp', lp_greedy)
        if( l_greedy )then
            call cline_cluster2D_chunk%set('maxits', 10.)
            call cline_cluster2D_chunk%set('refine', 'greedy')
        else
            call cline_cluster2D_chunk%set('refine',    'snhc')
            call cline_cluster2D_chunk%set('extr_iter', real(MAX_EXTRLIM2D-2))
            call cline_cluster2D_chunk%set('maxits',    12.)
        endif
        write(logfhandle,'(A,F5.1)')     '>>> CHUNK         LOW-PASS LIMIT (IN A) TO: ', lp_greedy
        if( l_greedy )then
            call cline_cluster2D_pool%set('refine', 'greedy')
            call cline_cluster2D_pool%set('lp',     lp_greedy)
            write(logfhandle,'(A,F5.1)') '>>> POOL          LOW-PASS LIMIT (IN A) TO: ', lp_greedy
        else
            call cline_cluster2D_pool%set('lpstart', lpstart_stoch)
            call cline_cluster2D_pool%set('lpstop',  params%lpstop2D)
            write(logfhandle,'(A,F5.1)') '>>> POOL STARTING LOW-PASS LIMIT (IN A) TO: ', lpstart_stoch
        endif
        write(logfhandle,'(A,F5.1)')     '>>> POOL   HARD RESOLUTION LIMIT (IN A) TO: ', params%lpstop2D
        ! initialize chunks
        allocate(chunks(params%nchunks))
        do ichunk = 1,params%nchunks
            call chunks(ichunk)%init(ichunk, pool_proj)
        enddo
        glob_chunk_id = params%nchunks
        ! general solution execution init
        call qenv_pool%new(params%nparts_pool,exec_bin='simple_private_exec',qsys_name='local')
        ! module variables
        n_spprojs_glob = 0
        if( allocated(spprojs_mask_glob) ) deallocate(spprojs_mask_glob)
        origproj_time  = simple_gettime()
        initiated      = .true.
        ! Particles/stacks map, number of chunks
        allocate(stk_all_nptcls(nstks),stk_nptcls(nstks),source=0)
        ntot_chunks     = 0
        cnt             = 0
        do istk = 1,nstks
            if( (spproj%os_stk%get_state(istk)==0) .or. (nint(spproj%os_stk%get(istk,'nptcls'))==0) )cycle
            fromp = nint(spproj%os_stk%get(istk,'fromp'))
            top   = nint(spproj%os_stk%get(istk,'top'))
            do iptcl = fromp,top
                stk_all_nptcls(istk) = stk_all_nptcls(istk) + 1 ! including state=0
                if( spproj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                stk_nptcls(istk) = stk_nptcls(istk) + 1 ! excluding state=0
            enddo
            cnt = cnt + stk_nptcls(istk)
            if( cnt > nptcls_per_chunk )then
                ntot_chunks = ntot_chunks + 1
                cnt = 0
            endif
        enddo
        nptcls_tot = sum(stk_nptcls)
        write(logfhandle,'(A,I6)')'>>> # OF STACKS   : ', nstks
        write(logfhandle,'(A,I6)')'>>> # OF PARTICLES: ', nptcls_tot
        write(logfhandle,'(A,I6)')'>>> # OF CHUNKS   : ', ntot_chunks
        ! chunks/stacks map
        allocate(chunks_map(ntot_chunks,2))
        cnt    = 0
        ichunk = 1
        do istk = 1,nstks
            if( ichunk > ntot_chunks ) exit
            if( cnt==0 ) chunks_map(ichunk,1) = istk
            cnt = cnt + stk_nptcls(istk)
            if( cnt > nptcls_per_chunk )then
                chunks_map(ichunk,2) = istk
                ichunk = ichunk + 1
                cnt = 0
            endif
        enddo
        chunks_map(ntot_chunks,2) = nstks ! border effect
        write(logfhandle,'(A)')'>>> CHUNKS MAP: '
        do ichunk = 1,ntot_chunks
            nptcls = sum(stk_nptcls(chunks_map(ichunk,1):chunks_map(ichunk,2)))
            write(logfhandle,'(A,I6,A,I6,A,I8)')'>>> CHUNK ID; # OF PARTICLES  : ',  ichunk, ' ; ',nptcls,' / ',sum(stk_nptcls)
        enddo
        ! Infinite loop
        allocate(completed_fnames(nstks))
        maxits = huge(maxits)     ! for convergence
        ichunk = 0                ! # of chunks that have been submitted
        tot_nchunks_imported = 0  ! Total # of chunks that are completed and imported into pool
        all_chunks_submitted = .false.
        all_chunks_imported  = .false.
        l_once = .true.
        do
            ! sequential chunk prep
            if( .not.all_chunks_submitted )then
                do ic = 1,params%nchunks
                    if( all_chunks_submitted ) exit
                    if( chunks(ic)%available )then
                        ichunk = ichunk + 1
                        ! setting nptcls_per_chunk with state=0 such that
                        ! sum(ptcls/state=1) = true nptcls_per_chunk to conform with generate_chunks
                        nptcls = sum(stk_all_nptcls(chunks_map(ichunk,1):chunks_map(ichunk,2)))
                        nptcls_per_chunk = nptcls-1
                        ! prep for current chunk
                        do istk = chunks_map(ichunk,1), chunks_map(ichunk,2)
                            one_projname = int2str_pad(istk, params%numlen)
                            one_projfile = trim(one_projname)//trim(METADATA_EXT)
                            call create_individual_project(spproj, istk, istk, dir_preprocess, one_projname, nptcls)
                            completed_fnames(istk) = trim(dir_preprocess)//trim(one_projfile)
                        enddo
                        ! submit chunk
                        tmp_fnames = completed_fnames(1:chunks_map(ichunk,2))
                        call update_projects_mask(tmp_fnames)
                        call classify_new_chunks(tmp_fnames)
                        if( ichunk == ntot_chunks )then
                            all_chunks_submitted = .true.
                            exit
                        endif
                    endif
                enddo
            endif
            ! streaming-like flow
            call update_chunks
            call update_pool_status
            call update_pool
            call reject_from_pool
            if( pool_iter /= maxits )then
                ! keep going
                if( all_chunks_imported )then
                    if( l_once )then
                        maxits = pool_iter + 2*STREAM_SRCHLIM ! for convergence
                        l_once = .false.
                        write(logfhandle,'(A)')'>>> ALL CHUNKS HAVE CONVERGED'
                        write(logfhandle,'(A,I6)')'>>> TERMINATING AT END OF ITERATION: ',maxits
                    endif
                else
                    call import_chunks_into_pool(.true., nchunks_imported)
                    tot_nchunks_imported = tot_nchunks_imported + nchunks_imported
                    all_chunks_imported  = tot_nchunks_imported == ntot_chunks
                endif
                if( pool_iter /= maxits ) call classify_pool
                if( .not.all_chunks_submitted ) call sleep(WAITTIME)
            else
                ! converged
                if( pool_available ) exit
            endif
        end do
        ! Termination
        call terminate_stream2D(.false.)
        deallocate(tmp_fnames,completed_fnames,stk_all_nptcls,stk_nptcls,chunks_map)
        ! maps stacks & gathering particles parameters
        ! only 2D parameters are transferred, everything else untouched
        pool_nstks = pool_proj%os_stk%get_noris()
        allocate(pool_stacks(pool_nstks), pool_stk_mask(pool_nstks))
        do istk = 1,pool_nstks
            call spproj%os_stk%getter(istk,'stk',fname)
            pool_stacks(istk) = basename(fname)  ! same as imported_stks ??
        enddo
        pool_stk_mask = .true.
        do istk = 1,nstks
            if( (spproj%os_stk%get_state(istk)==0) .or. (nint(spproj%os_stk%get(istk,'nptcls'))==0) )cycle
            call spproj%os_stk%getter(istk,'stk',fname)
            orig_stack_fname = basename(fname)
            stack_fname = trim(orig_stack_fname)
            if( do_autoscale ) stack_fname = add2fbody(orig_stack_fname,params%ext,SCALE_SUFFIX)
            do jstk = 1,pool_nstks
                if( pool_stk_mask(jstk) )then
                    if( trim(stack_fname) == pool_stacks(jstk) )then
                        fromp  = nint(spproj%os_stk%get(istk,'fromp'))
                        top    = nint(spproj%os_stk%get(istk,'top'))
                        jptcl  = nint(pool_proj%os_stk%get(jstk,'fromp'))
                        do iptcl = fromp,top
                            call spproj%os_ptcl2D%transfer_ori(iptcl, pool_proj%os_ptcl2D, jptcl)
                            jptcl = jptcl+1
                        enddo
                        pool_stk_mask(jstk) = .false. ! to be excluded from now
                    endif
                endif
            enddo
        enddo
        ! rescale class-averages & parameters
        refs_glob = trim(POOL_DIR)//trim(refs_glob)
        frcsfname = trim(POOL_DIR)//trim(FRCS_FILE)
        if( do_autoscale )then
            ! parameters
            call spproj%os_ptcl2D%mul_shifts(1.0/scale_factor)
            ! classes
            call rescale_cavgs(refs_glob, refs_glob)
            src  = add2fbody(refs_glob, params%ext,'_even')
            call rescale_cavgs(src, src)
            src  = add2fbody(refs_glob, params%ext,'_odd')
            call rescale_cavgs(src, src)
            call spproj%os_out%kill
            if( l_wfilt )then
                src  = add2fbody(refs_glob, params%ext,trim(WFILT_SUFFIX))
                call rescale_cavgs(src, src)
                src  = add2fbody(refs_glob, params%ext,trim(WFILT_SUFFIX)//'_even')
                call rescale_cavgs(src, src)
                src  = add2fbody(refs_glob, params%ext,trim(WFILT_SUFFIX)//'_odd')
                call rescale_cavgs(src, src)
            endif
            ! frcs
            call frcs_sc%read(trim(POOL_DIR)//trim(FRCS_FILE))
            call frcs_sc%upsample(orig_smpd, orig_box, frcs)
            call frcs%write(frcsfname)
            call frcs%kill
            call frcs_sc%kill
        endif
        call spproj%os_out%kill
        call spproj%add_cavgs2os_out(refs_glob, orig_smpd, 'cavg')
        if( l_wfilt )then
            src = add2fbody(refs_glob,params%ext,trim(WFILT_SUFFIX))
            call spproj%add_cavgs2os_out(src, orig_smpd, 'cavg'//trim(WFILT_SUFFIX))
        endif
        spproj%os_cls2D = pool_proj%os_cls2D
        call spproj%add_frcs2os_out(frcsfname, 'frc2D')
        call pool_proj%kill
        ! classes export
        if(file_exists('clusters2D.star')) call del_file('clusters2D.star')
        call starproj%export_cls2D(spproj)
        call spproj%os_cls2D%delete_entry('stk')
        ! write whole project
        if( os_mic%get_noris() > 0 )then
            spproj%os_mic = os_mic
            call os_mic%kill
        endif
        spproj%os_ptcl3D = spproj%os_ptcl2D
        call spproj%os_ptcl3D%delete_2Dclustering
        call spproj%write(orig_projfile)
        ! ranking
        call rank_cavgs
        ! cleanup
        call simple_rmdir(STDERROUT_DIR)
        call simple_rmdir(dir_preprocess)
        call del_file(trim(POOL_DIR)//trim(PROJFILE_POOL))
        call qsys_cleanup
        ! graceful end
        call simple_end('**** SIMPLE_CLUSTER2D_SUBSETS NORMAL STOP ****')
    end subroutine exec_cluster2D_subsets

    subroutine create_individual_project( sp_proj, stkind, proj_ind, output_dir, spprojname, nptcls )
        class(sp_project), intent(in)  :: sp_proj
        integer,           intent(in)  :: stkind, proj_ind
        character(len=*),  intent(in)  :: output_dir, spprojname
        integer,           intent(out) :: nptcls
        type(sp_project)              :: spproj_here
        type(cmdline)                 :: cline_here
        character(len=:), allocatable :: stack_fname
        character(len=LONGSTRLEN)     :: projfile, stack_relpath
        character(len=XLONGSTRLEN)    :: cwd, cwd_old
        integer :: fromp, top, i, iptcl
        cwd_old = trim(cwd_glob)
        call chdir(output_dir)
        call simple_getcwd(cwd)
        cwd_glob = trim(cwd)
        nptcls = nint(sp_proj%os_stk%get(stkind,'nptcls'))
        call spproj_here%os_mic%new(1, is_ptcl=.false.)
        if( nptcls==0 .or. sp_proj%os_stk%get_state(stkind)==0 )then
            ! we do not care about the mic_field here
            call spproj_here%os_mic%set(1,'nptcls',0.0)
            call spproj_here%os_mic%set_state(1,0)
        else
            fromp = nint(sp_proj%os_stk%get(stkind,'fromp'))
            top   = nint(sp_proj%os_stk%get(stkind,'top'))
            call spproj_here%os_stk%new(1, is_ptcl=.false.)
            call spproj_here%os_ptcl2D%new(nptcls, is_ptcl=.true.)
            if( sp_proj%os_mic%get_noris() >= proj_ind )then
                call spproj_here%os_mic%transfer_ori(1, sp_proj%os_mic, stkind)
            else
                ! in case mic field is not informed
                call spproj_here%os_mic%set(1,'nptcls',real(nptcls))
                call spproj_here%os_mic%set_state(1,1)
            endif
            call spproj_here%os_stk%transfer_ori(1,sp_proj%os_stk,stkind)
            call spproj_here%os_stk%set(1,'fromp',1.0)
            call spproj_here%os_stk%set(1,'top',  real(nptcls))
            call spproj_here%os_stk%getter(1,'stk',stack_fname)
            if( stack_fname(1:1) /= '/') stack_fname = trim('../')//trim(stack_fname)
            call make_relativepath(cwd_glob, stack_fname, stack_relpath)
            call spproj_here%os_stk%set(1,'stk',stack_relpath)
            iptcl = 0
            do i = fromp,top
                iptcl = iptcl + 1
                call spproj_here%os_ptcl2D%transfer_ori(iptcl, sp_proj%os_ptcl2D, i)
                call spproj_here%os_ptcl2D%set(iptcl,'stkind',1.0)
            enddo
            call spproj_here%os_ptcl2D%delete_2Dclustering
        endif
        projfile = trim(spprojname)//trim(METADATA_EXT)
        call cline_here%set('projname', spprojname)
        call cline_here%set('projfile', projfile)
        call spproj_here%update_projinfo(cline_here)
        spproj_here%compenv  = sp_proj%compenv
        spproj_here%jobproc  = sp_proj%jobproc
        call spproj_here%write
        call chdir(cwd_old)
        cwd_glob = trim(cwd_old)
        call spproj_here%kill
        call cline_here%kill
    end subroutine create_individual_project

end module simple_commander_cluster2D_stream_dev
