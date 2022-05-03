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
integer                         :: orig_box, box, boxpd
integer                         :: glob_chunk_id, max_ncls, nptcls_per_chunk, ncls_glob
real                            :: lp_greedy, lpstart_stoch
logical                         :: l_wfilt, do_autoscale, l_greedy

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
        ! if( cline%defined('lp2D') ) lp_greedy = params_glob%lp2D
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
    end subroutine init_cluster2D_stream

    ! Utilities

    subroutine debug_print( string )
        character(len=*), intent(in) :: string
        if( DEBUG_HERE )then
            write(logfhandle,*) trim(string)
            call flush(logfhandle)
        endif
    end subroutine debug_print

end module simple_commander_cluster2D_stream_dev
    