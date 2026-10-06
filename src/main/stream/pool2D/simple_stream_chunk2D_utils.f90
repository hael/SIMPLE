!@descr: utilities for chunk-based 2D clustering in stream
module simple_stream_chunk2D_utils
use simple_core_module_api
use simple_defs_environment
use simple_cmdline,                only: cmdline
use simple_parameters,             only: parameters
use simple_sp_project,             only: sp_project
use simple_stream_chunk,           only: stream_chunk
use simple_stream_refine2D_utils,  only: setup_downscaling
implicit none

! LIFECYCLE
public :: init_chunk_clustering
private
#include "simple_local_flags.inc"

contains

    ! LIFECYCLE

    ! The chunks of solve2D_chunks (@p chunks, params%nchunks of them) and their solve2D command
    ! line @p cline_refine2D_chunk, from the program's @p params and command line @p cline; the
    ! chunks take the computing environment of @p spproj.
    subroutine init_chunk_clustering( params, cline, spproj, cline_refine2D_chunk, chunks )
        class(parameters),               intent(inout) :: params
        class(cmdline),                  intent(inout) :: cline
        class(sp_project),               intent(inout) :: spproj
        type(cmdline),                   intent(inout) :: cline_refine2D_chunk
        type(stream_chunk), allocatable, intent(inout) :: chunks(:)
        type(sp_project)      :: template   ! the computing environment the chunks copy
        character(len=STDLEN) :: chunk_nthr_env
        integer               :: ichunk, envlen
        real                  :: lp_fixed, lpstart, lpstop, lpcen
        call seed_rnd
        ! general parameters
        params%ncls_start   = params%ncls ! backwards compatibility
        params%nparts_chunk = params%nparts ! required by chunk object, to remove
        ! the computing environment, used upon chunk generation
        template%projinfo = spproj%projinfo
        template%compenv  = spproj%compenv
        call template%projinfo%delete_entry('projname')
        call template%projinfo%delete_entry('projfile')
        if( cline%defined('walltime') ) call template%compenv%set(1,'walltime', params%walltime)
        ! chunk master command line
        call cline_refine2D_chunk%set('prg', 'solve2D')
        if( params%nparts > 1 )then
            call cline_refine2D_chunk%set('nparts',        params%nparts)
        endif
        if( cline%defined('cls_init') )then
            call cline_refine2D_chunk%set('cls_init',      params%cls_init)
        else
            call cline_refine2D_chunk%set('cls_init',      'rand')
        endif
        if( cline%defined('gaufreq') )then
            call cline_refine2D_chunk%set('gaufreq',       params%gaufreq)
        endif
        call cline_refine2D_chunk%set('oritype',    'ptcl2D')
        call cline_refine2D_chunk%set('center',     'no')
        call cline_refine2D_chunk%set('autoscale', 'no')
        call cline_refine2D_chunk%set('mkdir',      'no')
        call cline_refine2D_chunk%set('mskdiam',    params%mskdiam)
        call cline_refine2D_chunk%set('ncls',       params%ncls_start)
        call cline_refine2D_chunk%set('sigma_est', params%sigma_est)
        call cline_refine2D_chunk%set('rank_cavgs','yes')
        call cline_refine2D_chunk%set('chunk',      'yes')
        ! objective function
        call cline_refine2D_chunk%set('objfun', 'euclid')
        call cline_refine2D_chunk%set('ml_reg', params%ml_reg)
        call cline_refine2D_chunk%set('tau',     params%tau)
        ! refinement
        select case(trim(params%refine))
               case('snhc','snhc_smpl','prob','prob_snhc')
                call cline_refine2D_chunk%set('refine', params%refine)
            case DEFAULT
                THROW_HARD('UNSUPPORTED REFINE PARAMETER!')
        end select
        ! Determines dimensions for downscaling
        call set_chunk_dimensions( params )
        ! updates command-line with resolution limits, defaults are handled by solve2D
        if( cline%defined('lp') )then
            lp_fixed = max(params%lp, 2.0*params%smpd_crop)
            call cline_refine2D_chunk%set('lp', lp_fixed)
            write(logfhandle,'(A,F5.1)') '>>> FIXED RESOLUTION LIMIT    (IN A): ', lp_fixed
        else
            if( cline%defined('lpstart') )then
                lpstart = max(params%lpstart, 2.0*params%smpd_crop)
                call cline_refine2D_chunk%set('lpstart', lpstart)
                write(logfhandle,'(A,F5.1)') '>>> STARTING RESOLUTION LIMIT (IN A): ', lpstart
            endif
            if( cline%defined('lpstop') )then
                lpstop = max(params%lpstop, 2.0*params%smpd_crop)
                call cline_refine2D_chunk%set('lpstop', lpstop)
                write(logfhandle,'(A,F5.1)') '>>> HARD RESOLUTION LIMIT     (IN A): ', lpstop
            endif
        endif
        if( cline%defined('cenlp') )then
            lpcen = max(params%cenlp, 2.0*params%smpd_crop)
            call cline_refine2D_chunk%set('cenlp', lpcen)
            write(logfhandle,'(A,F5.1)') '>>> CENTERING LOW-PASS LIMIT  (IN A): ', lpcen
        endif
        ! EV override
        call get_environment_variable(SIMPLE_STREAM_CHUNK_NTHR, chunk_nthr_env, envlen)
        if(envlen > 0) then
            call cline_refine2D_chunk%set('nthr', str2int(chunk_nthr_env))
        else
            call cline_refine2D_chunk%set('nthr', params%nthr) ! cf comment just below about nthr2D
        end if
        ! Initialize subsets
        if( allocated(chunks) ) deallocate(chunks)
        allocate(chunks(params%nchunks))
        ! deal with nthr2d .ne. nthr
        ! Joe: the whole nthr/2d is confusing. Why not pass the number of threads to chunk%init?
        params%nthr2D = cline_refine2D_chunk%get_iarg('nthr') ! only used here   for backwards compatibility
        do ichunk = 1,params%nchunks
            call chunks(ichunk)%init_chunk(params, cline_refine2D_chunk, ichunk, template)
        enddo
        call template%kill

    contains

        subroutine set_chunk_dimensions( params )
            class(parameters), intent(inout) :: params
            type(scaled_dims) :: chunk_dims
            logical           :: l_scaling
            l_scaling = .true.
            call setup_downscaling(params, l_scaling)
            chunk_dims%smpd  = params%smpd_crop
            chunk_dims%box   = params%box_crop
            chunk_dims%boxpd = 2 * round2even(KBALPHA * real(params%box_crop/2)) ! logics from parameters
            chunk_dims%msk   = params%msk_crop
            ! Scaling-related command lines update
            call cline_refine2D_chunk%set('smpd_crop', chunk_dims%smpd)
            call cline_refine2D_chunk%set('box_crop',   chunk_dims%box)
            call cline_refine2D_chunk%set('msk_crop',   chunk_dims%msk)
            call cline_refine2D_chunk%set('box',        params%box)
            call cline_refine2D_chunk%set('smpd',       params%smpd)
        end subroutine set_chunk_dimensions

    end subroutine init_chunk_clustering

end module simple_stream_chunk2D_utils
