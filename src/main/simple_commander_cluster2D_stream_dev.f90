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
use simple_qsys_funs
use simple_commander_cluster2D
use simple_timer
implicit none

public :: init_cluster2D_stream
public :: update_chunk_mask, update_chunks, start_new_chunks

private
#include "simple_local_flags.inc"

integer,               parameter   :: MINBOXSZ            = 128    ! minimum boxsize for scaling
real,                  parameter   :: GREEDY_TARGET_LP    = 15.0
integer,               parameter   :: WAIT_WATCHER        = 5     ! seconds prior to new stack detection
! integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 600  ! dev settings
integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 7200  ! Frequency at which the original project file should be updated
integer,               parameter   :: FREQ_POOL_REJECTION = 5     !
character(len=STDLEN), parameter   :: USER_PARAMS         = 'stream2D_user_params.txt'
character(len=STDLEN), parameter   :: LAST_SNAPSHOT       = 'last_snapshot.txt'
character(len=STDLEN), parameter   :: SPPROJ_SNAPSHOT     = 'SIMPLE_PROJECT_SNAPSHOT'
character(len=STDLEN), parameter   :: PROJFILE_POOL       = 'cluster2D.simple'
character(len=STDLEN), parameter   :: SCALE_DIR           = './scaled_stks/'
character(len=STDLEN), parameter   :: POOL_DIR           = './pool/'
character(len=STDLEN), parameter   :: DISTR_EXEC_FNAME    = './distr_cluster2D_pool'
logical,               parameter   :: DEBUG_HERE          = .true.
integer(timer_int_kind)            :: t

type(stream_chunk), allocatable :: chunks(:)
type(sp_project)                :: pool_proj
type(qsys_env)                  :: qenv_pool
type(cmdline)                   :: cline_cluster2D_chunk, cline_cluster2D_pool
character(len=:),   allocatable :: orig_projfile
character(len=LONGSTRLEN)       :: prev_snapshot_frcs, prev_snapshot_cavgs
real                            :: orig_smpd, smpd, scale_factor, mskdiam
real                            :: lp_greedy, lpstart_stoch
integer                         :: orig_box, box, boxpd
integer                         :: glob_chunk_id, max_ncls, nptcls_per_chunk, ncls_glob
integer                         :: n_spprojs_glob = 0
logical                         :: l_wfilt, do_autoscale, l_greedy
logical                         :: initiated = .false.
logical, allocatable :: spprojs_mask_glob(:)

contains

    subroutine init_cluster2D_stream( cline, spproj, refs, do2D )
        class(cmdline),    intent(inout) :: cline
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: refs
        logical,           intent(inout) :: do2D
        real    :: SMPD_TARGET = MAX_SMPD  ! target sampling distance
        integer :: ldim(3), ichunk, maybe2D, ifoo
        ! check whether 2D classification will be performed
        maybe2D = merge(1,0,cline%defined('nchunks'))
        maybe2D = maybe2D + merge(1,0,cline%defined('nparts_chunk'))
        maybe2D = maybe2D + merge(1,0,cline%defined('ncls'))
        maybe2D = maybe2D + merge(1,0,cline%defined('nptcls_per_cls'))
        maybe2D = maybe2D + merge(1,0,cline%defined('ncls_start'))
        maybe2D = maybe2D + merge(1,0,cline%defined('mskdiam'))
        if( maybe2D == 6 )then
            do2D = .true.
        else if( maybe2D > 0 )then
            THROW_HARD('Missing arguments for 2D classification')
        else
            do2D = .false.
        endif
        if( .not.do2D ) return
        call seed_rnd
        ! general parameters
        mskdiam           = cline%get_rarg('mskdiam')
        l_wfilt           = trim(params_glob%wiener) .eq. 'partial'
        do_autoscale      = trim(params_glob%autoscale) .eq. 'yes'
        max_ncls          = floor(real(params_glob%ncls)/real(params_glob%ncls_start))*params_glob%ncls_start ! effective maximum # of classes
        nptcls_per_chunk  = params_glob%nptcls_per_cls*params_glob%ncls_start         ! # of particles in each chunk
        ncls_glob         = 0
        l_greedy          = trim(params_glob%refine).eq.'greedy'
        lp_greedy         = GREEDY_TARGET_LP
        if( cline%defined('lp2D') ) lp_greedy = params_glob%lp2D
        lpstart_stoch     = lp_greedy
        prev_snapshot_cavgs = ''
        orig_projfile = trim(params_glob%projfile)
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
        ! call write_user_params
        call simple_mkdir(SCALE_DIR)
        call simple_mkdir(POOL_DIR)
        call pool_proj%kill
        pool_proj%projinfo = spproj%projinfo
        pool_proj%compenv  = spproj%compenv
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        ! initialize chunks parameters and objects
        call cline_cluster2D_chunk%set('prg',       'cluster2D_distr')
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
        call cline_cluster2D_chunk%set('nparts',    real(params_glob%nparts_chunk))
        call cline_cluster2D_chunk%set('nthr',      real(params_glob%nthr))
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
        if( .not.cline%defined('match_filt') ) call cline_cluster2D_pool%set('match_filt','no')
        if( .not. cline%defined('cenlp')     ) call cline_cluster2D_pool%set('cenlp',     30.0)
        if( .not. cline%defined('center')    ) call cline_cluster2D_pool%set('center',    'yes')
        call cline_cluster2D_pool%set('extr_iter', 100.)
        call cline_cluster2D_pool%set('mkdir',     'no')
        call cline_cluster2D_pool%set('mskdiam',   mskdiam)
        call cline_cluster2D_pool%set('async',     'yes') ! to enable hard termination
        call cline_cluster2D_pool%set('stream',    'yes') ! use for dual CTF treatment
        call cline_cluster2D_pool%set('nthr',     real(params_glob%nthr))
        call cline_cluster2D_pool%set('nparts',   real(params_glob%nparts))
        call qenv_pool%new(params_glob%nparts,exec_bin='simple_private_exec',qsys_name='local')
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
        ! module varaiables
        n_spprojs_glob = 0
        if( allocated(spprojs_mask_glob) ) deallocate(spprojs_mask_glob)
        initiated      = .true.
    end subroutine init_cluster2D_stream

    ! updates mask  of project to process, check on chunks
    subroutine update_chunk_mask( spproj_fnames )
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
            THROW_HARD('Error update_chunks 1')
        endif
        call debug_print('update_chunk_mask n_spprojs_glob: '//int2str(n_spprojs_glob))
    end subroutine update_chunk_mask

    subroutine update_chunks
        integer :: ichunk
        ! deal with chunk completion, rejection, reset, import into pool
        do ichunk = 1,params_glob%nchunks
            if( chunks(ichunk)%available ) cycle
            if( chunks(ichunk)%has_converged() )then
                call chunks(ichunk)%display_iter
                ! rejection
                call chunks(ichunk)%reject(params_glob%lpthresh, params_glob%ndev, box)
                ! import into pool and folder removal happens here
                ! TODO
                ! free chunk
                glob_chunk_id = glob_chunk_id + 1
                call chunks(ichunk)%init(glob_chunk_id, pool_proj)
            endif
        enddo
    end subroutine update_chunks

    subroutine start_new_chunks( spproj_fnames )
        character(len=LONGSTRLEN), allocatable, intent(in) :: spproj_fnames(:)
        type(sp_project)                       :: stream_proj
        character(len=LONGSTRLEN), allocatable :: spprojs_for_chunk(:)
        integer,                   allocatable :: spproj_nptcls(:)
        integer :: ichunk, n_avail_chunks, n_spprojs_in, iproj, nptcls, n2fill
        integer :: first2import, last2import, n2import, cnt
        integer :: inmics, instks, inptcls
        n_avail_chunks = count(chunks(:)%available)
        call debug_print('start_new_chunks n_avail_chunks: '//int2str(n_avail_chunks))
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

            ! call stream_proj%read_segment('mic',spproj_fnames(iproj))
            ! spproj_nptcls(iproj) = nint(stream_proj%os_mic%get(1,'nptcls'))

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
    end subroutine

    ! Utilities

    subroutine debug_print( string )
        character(len=*), intent(in) :: string
        if( DEBUG_HERE )then
            write(logfhandle,*) trim(string)
            call flush(logfhandle)
        endif
    end subroutine debug_print

end module simple_commander_cluster2D_stream_dev
    