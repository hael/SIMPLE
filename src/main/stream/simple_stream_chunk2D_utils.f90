module simple_stream_chunk2D_utils
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_cmdline,      only: cmdline
use simple_stream_chunk, only: stream_chunk
use simple_sp_project,   only: sp_project
use simple_stream_cluster2D_utils
use simple_stream2D_state
use simple_gui_utils
use simple_rec_list
implicit none

! LIFECYCLE
public :: init_chunk_clustering
public :: analyze2D_new_chunks
public :: memoize_chunks
! GETTERS
public :: all_chunks_available
public :: get_chunk_rejected_jpeg
public :: get_chunk_rejected_jpeg_scale
public :: get_nchunks
! SETTERS
public :: set_chunk_dimensions
! UPDATERS
public :: update_chunks
private
#include "simple_local_flags.inc"

! Chunk rejection view
type(string) :: chunk_rejected_jpeg
integer      :: chunk_rejected_jpeg_ntiles  = 0
integer      :: chunk_rejected_thumbnail_id = 0
real         :: chunk_rejected_jpeg_scale   = 1.0

contains

    ! LIFECYCLE

    subroutine init_chunk_clustering( cline, spproj )
        class(cmdline),    target, intent(inout) :: cline
        class(sp_project),         intent(inout) :: spproj
        character(len=STDLEN) :: chunk_nthr_env
        integer               :: ichunk, envlen
        call seed_rnd
        ! general parameters
        master_cline => cline
        l_wfilt          = .false.
        l_scaling        = .true.
        params_glob%ncls_start = params_glob%ncls ! backwards compatibility
        nptcls_per_chunk = params_glob%nptcls_per_cls*params_glob%ncls_start
        ncls_glob        = 0
        l_update_sigmas  = params_glob%l_needs_sigma
        numlen           = len(int2str(params_glob%nparts))
        l_no_chunks      = .false. ! will be using chunk indeed
        l_abinitio2D     = cline%defined('algorithm')
        if( l_abinitio2D ) l_abinitio2D = str_has_substr(params_glob%algorithm,'abinitio')
        params_glob%nparts_chunk = params_glob%nparts ! required by chunk object, to remove
        ! bookkeeping & directory structure
        if( l_update_sigmas ) call simple_mkdir(SIGMAS_DIR)
        ! pool_proj is only iused as a placeholder for computational info here
        ! used upon chunk generation
        call pool_proj%kill
        pool_proj%projinfo = spproj%projinfo
        pool_proj%compenv  = spproj%compenv
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        if( cline%defined('walltime') ) call pool_proj%compenv%set(1,'walltime', params_glob%walltime)
        ! chunk master command line
        if( l_abinitio2D )then
            call cline_cluster2D_chunk%set('prg', 'abinitio2D')
            if( params_glob%nparts > 1 )then
                call cline_cluster2D_chunk%set('nparts',       params_glob%nparts)
            endif
            if( cline%defined('cls_init') )then
                call cline_cluster2D_chunk%set('cls_init',     params_glob%cls_init)
            else
                call cline_cluster2D_chunk%set('cls_init',     'rand')
            endif
            if( master_cline%defined('focusmskdiam') )then
                call cline_cluster2D_chunk%set('focusmskdiam', params_glob%focusmskdiam)
            endif
            if( cline%defined('gaufreq') )then
                call cline_cluster2D_chunk%set('gaufreq',      params_glob%gaufreq)
            endif
        else
            if( params_glob%nparts > 1 )then
                call cline_cluster2D_chunk%set('prg',    'cluster2D_distr')
                call cline_cluster2D_chunk%set('nparts', params_glob%nparts)
            else
                ! shared memory execution
                call cline_cluster2D_chunk%set('prg',    'cluster2D')
            endif
            call cline_cluster2D_chunk%set('minits',    CHUNK_MINITS)
            call cline_cluster2D_chunk%set('maxits',    CHUNK_MAXITS)
            call cline_cluster2D_chunk%set('extr_iter', CHUNK_EXTR_ITER)
            call cline_cluster2D_chunk%set('extr_lim',  MAX_EXTRLIM2D)
            call cline_cluster2D_chunk%set('startit',   1)
            if( l_update_sigmas ) call cline_cluster2D_chunk%set('cc_iters', CHUNK_CC_ITERS)
            if( cline%defined('cls_init') )then
                call cline_cluster2D_chunk%set('cls_init', params_glob%cls_init)
            else
                call cline_cluster2D_chunk%set('cls_init','ptcl')
            endif
        endif
        call cline_cluster2D_chunk%set('oritype',   'ptcl2D')
        call cline_cluster2D_chunk%set('center',    'no')
        call cline_cluster2D_chunk%set('autoscale', 'no')
        call cline_cluster2D_chunk%set('mkdir',     'no')
        call cline_cluster2D_chunk%set('stream',    'no')
        call cline_cluster2D_chunk%set('mskdiam',   params_glob%mskdiam)
        call cline_cluster2D_chunk%set('ncls',      params_glob%ncls_start)
        call cline_cluster2D_chunk%set('sigma_est', params_glob%sigma_est)
        call cline_cluster2D_chunk%set('rank_cavgs','no')
        call cline_cluster2D_chunk%set('chunk',     'yes')
        if( l_wfilt )then
            call cline_cluster2D_chunk%set('wiener',     'partial')
        endif
        ! objective function
        call cline_cluster2D_chunk%set('objfun', 'euclid')
        call cline_cluster2D_chunk%set('ml_reg', params_glob%ml_reg)
        call cline_cluster2D_chunk%set('tau',    params_glob%tau)
        ! refinement
        select case(trim(params_glob%refine))
            case('snhc','snhc_smpl','prob','prob_smpl')
                if( (.not.l_abinitio2D) .and. str_has_substr(params_glob%refine,'prob') )then
                    THROW_HARD('REFINE=PROBXX only compatible with algorithm=abinitio2D')
                endif
                call cline_cluster2D_chunk%set('refine', params_glob%refine)
            case DEFAULT
                THROW_HARD('UNSUPPORTED REFINE PARAMETER!')
        end select
        ! polar representation
        if( master_cline%defined('polar') ) call cline_cluster2D_chunk%set('polar', params_glob%polar)
        ! Determines dimensions for downscaling
        call set_chunk_dimensions
        ! updates command-line with resolution limits, defaults are handled by abinitio2D
        if( master_cline%defined('lpstart') )then
            lpstart = max(params_glob%lpstart, 2.0*params_glob%smpd_crop)
            call cline_cluster2D_chunk%set('lpstart', lpstart)
            write(logfhandle,'(A,F5.1)') '>>> STARTING RESOLUTION LIMIT (IN A): ', lpstart
        endif
        if( master_cline%defined('lpstop') )then
            lpstop = max(params_glob%lpstop, 2.0*params_glob%smpd_crop)
            call cline_cluster2D_chunk%set('lpstop', lpstop)
            write(logfhandle,'(A,F5.1)') '>>> HARD RESOLUTION LIMIT     (IN A): ', lpstop
        endif
        if( master_cline%defined('cenlp') )then
            lpcen = max(params_glob%cenlp, 2.0*params_glob%smpd_crop)
            call cline_cluster2D_chunk%set('cenlp', lpcen)
            write(logfhandle,'(A,F5.1)') '>>> CENTERING LOW-PASS LIMIT  (IN A): ', lpcen
        endif
        ! EV override
        call get_environment_variable(SIMPLE_STREAM_CHUNK_NTHR, chunk_nthr_env, envlen)
        if(envlen > 0) then
            call cline_cluster2D_chunk%set('nthr', str2int(chunk_nthr_env))
        else
            call cline_cluster2D_chunk%set('nthr', params_glob%nthr) ! cf comment just below about nthr2D
        end if
        ! Initialize subsets
        allocate(chunks(params_glob%nchunks))
        ! deal with nthr2d .ne. nthr
        ! Joe: the whole nthr/2d is confusing. Why not pass the number of threads to chunk%init?
        params_glob%nthr2D = cline_cluster2D_chunk%get_iarg('nthr') ! only used here  for backwards compatibility
        glob_chunk_id      = 0
        do ichunk = 1,params_glob%nchunks
            glob_chunk_id = glob_chunk_id + 1
            call chunks(ichunk)%init_chunk(ichunk, cline_cluster2D_chunk, pool_proj)
        enddo
        ! module variables
        l_stream2D_active = .true.
    end subroutine init_chunk_clustering

    ! Initiates analysis of all available chunks
    subroutine analyze2D_new_chunks( project_list, makecavgs )
        class(rec_list),   intent(inout) :: project_list
        logical, optional, intent(in)    :: makecavgs
        type(project_rec)  :: prec
        type(rec_iterator) :: it
        type(rec_list)     :: project_list_slice
        integer :: ichunk, n_avail_chunks, n_spprojs_in, iproj, nptcls, n2fill
        integer :: first2import, last2import, n2import
        if( .not. l_stream2D_active ) return
        n_avail_chunks = count(chunks(:)%is_available())
        ! cannot import yet
        if( n_avail_chunks == 0 ) return
        n_spprojs_in = project_list%size()
        if( n_spprojs_in == 0 ) return
        ! how many n2fill chunks to load
        n2fill       = 0
        nptcls       = 0
        first2import = 0
        it           = project_list%begin()
        do iproj = 1,n_spprojs_in
            ! retrieve one record from the list with the iterator
            call it%get(prec)
            if( prec%included )then
                ! move the iterator
                call it%next()
                cycle
            endif
            if( prec%nptcls_sel > 0 )then
                if( first2import == 0 ) first2import = iproj
                nptcls = nptcls + prec%nptcls_sel
                if( nptcls >= nptcls_per_chunk )then
                    n2fill = n2fill + 1
                    if( n2fill >= n_avail_chunks )exit
                    nptcls = 0
                endif
            else
                ! mask out empty stacks
                prec%included = .true. 
                ! replace the node
                call project_list%replace_iterator(it, prec)
            endif
            ! move the iterator
            call it%next()
        enddo
        if( n2fill == 0 ) return ! not enough particles
        do ichunk = 1,params_glob%nchunks
            if(.not.chunks(ichunk)%is_available()) cycle
            if( n2fill == 0 ) exit
            n2fill   = n2fill - 1
            nptcls   = 0
            n2import = 0
            it       = project_list%begin()
            do iproj = first2import,n_spprojs_in
                ! retrieve one record from the list with the iterator
                call it%get(prec)
                nptcls   = nptcls   + prec%nptcls_sel
                n2import = n2import + 1
                if( nptcls >= nptcls_per_chunk )then
                    last2import = iproj
                    exit
                endif
                ! move the iterator
                call it%next()
            enddo
            if( nptcls >= nptcls_per_chunk )then
                ! need a slice of project_list here
                call project_list%slice(first2import, last2import, project_list_slice)
                ! generate chunk from slice
                call chunks(ichunk)%generate(project_list_slice)
                ! flag inclusion in original list
                call project_list%set_included_flags([first2import,last2import])
                ! execution
                call chunks(ichunk)%analyze2D(l_update_sigmas, makecavgs)
                first2import = last2import + 1 ! to avoid cycling through all projects
                call project_list_slice%kill
            endif
        enddo
    end subroutine analyze2D_new_chunks

    ! Chunks Book-keeping
    subroutine memoize_chunks( list, nchunks_imported )
        class(rec_list), intent(inout) :: list
        integer,         intent(out)   :: nchunks_imported
        type(string) :: fname
        integer      :: i, id, nchunks2import
        nchunks_imported = 0
        if( .not.allocated(converged_chunks) ) return
        nchunks2import = size(converged_chunks)
        do i = 1,nchunks2import
            fname = converged_chunks(i)%get_projfile_fname()
            id    = converged_chunks(i)%get_id()
            ! append to list
            call list%push2chunk_list(fname, id, .false.)
            ! sigma2 book-keeping
            call converged_chunks(i)%split_sigmas_into(string(SIGMAS_DIR))
            ! destroy chunk
            call converged_chunks(i)%kill
        enddo
        nchunks_imported = nchunks2import
        deallocate(converged_chunks)
    end subroutine memoize_chunks

    ! GETTERS

    ! Are all chunks inactive
    logical function all_chunks_available()
        if( params_glob%nchunks == 0 )then
            all_chunks_available = .true.
        else
            all_chunks_available = all(chunks(:)%is_available())
        endif
    end function all_chunks_available

    type(string) function get_chunk_rejected_jpeg()
        get_chunk_rejected_jpeg = chunk_rejected_jpeg
    end function get_chunk_rejected_jpeg

    real function get_chunk_rejected_jpeg_scale()
        get_chunk_rejected_jpeg_scale = chunk_rejected_jpeg_scale
    end function get_chunk_rejected_jpeg_scale

    integer function get_nchunks()
        if(allocated(chunks)) then
            get_nchunks = size(chunks)
        else
            get_nchunks = 0
        end if
    end function get_nchunks

    ! SETTERS
    
    subroutine set_chunk_dimensions
        call setup_downscaling
        chunk_dims%smpd  = params_glob%smpd_crop
        chunk_dims%box   = params_glob%box_crop
        chunk_dims%boxpd = 2 * round2even(params_glob%alpha * real(params_glob%box_crop/2)) ! logics from parameters
        chunk_dims%msk   = params_glob%msk_crop
        ! Scaling-related command lines update
        call cline_cluster2D_chunk%set('smpd_crop', chunk_dims%smpd)
        call cline_cluster2D_chunk%set('box_crop',  chunk_dims%box)
        call cline_cluster2D_chunk%set('msk_crop',  chunk_dims%msk)
        call cline_cluster2D_chunk%set('box',       params_glob%box)
        call cline_cluster2D_chunk%set('smpd',      params_glob%smpd)
    end subroutine set_chunk_dimensions

    ! UPDATERS

    ! Deals with chunk completion, rejection, reset
    subroutine update_chunks
        type(stream_chunk), allocatable :: tmpchunks(:)
        integer :: ichunk, jchunk, nthr2D, n
        logical :: chunk_complete
        if( .not. l_stream2D_active ) return
        do ichunk = 1,params_glob%nchunks
            if( chunks(ichunk)%is_available() ) cycle
            chunk_complete = .false.
            if( chunks(ichunk)%to_analyze2D() )then
                ! chunk meant to be classified
                if( chunks(ichunk)%has_converged() )then
                    chunk_complete = .true.
                    call chunks(ichunk)%display_iter
                    ! rejection
                    if( trim(params_glob%reject_cls).ne.'no' )then
                        call chunks(ichunk)%reject(params_glob%lpthres, params_glob%ndev)
                        call mrc2jpeg_tiled(string('cls_rejected_chunks.mrc'), string('cls_rejected_chunks.jpeg'),&
                        &scale=chunk_rejected_jpeg_scale, ntiles=chunk_rejected_jpeg_ntiles)
                        chunk_rejected_jpeg = CWD_GLOB // '/' // 'cls_rejected_chunks.jpeg'
                        chunk_rejected_thumbnail_id = chunk_rejected_jpeg_ntiles
                    endif
                endif
            else
                ! placeholder chunk (no analysis performed, sigma2 only)
                if( chunks(ichunk)%has_converged() ) chunk_complete = .true.
            endif
            if( chunk_complete )then
                ! updates list of chunks to import
                if( allocated(converged_chunks) )then
                    ! append item
                    n = size(converged_chunks)
                    allocate(tmpchunks(n+1),source=[converged_chunks(:), chunks(ichunk)])
                    do jchunk = 1,n
                        call converged_chunks(jchunk)%kill
                    enddo
                    deallocate(converged_chunks)
                    allocate(converged_chunks(n+1),source=tmpchunks)
                    do jchunk = 1,n+1
                        call tmpchunks(jchunk)%kill
                    enddo
                    deallocate(tmpchunks)
                else
                    ! first item
                    allocate(converged_chunks(1),source=[chunks(ichunk)])
                endif
                ! reinit and deal with nthr2D != nthr
                glob_chunk_id = glob_chunk_id + 1
                ! deal with nthr2d .ne. nthr
                nthr2D = params_glob%nthr2D
                params_glob%nthr2D = cline_cluster2D_chunk%get_iarg('nthr')
                call chunks(ichunk)%init_chunk(glob_chunk_id, cline_cluster2D_chunk, pool_proj)
                params_glob%nthr2D = nthr2D
            endif
        enddo
    end subroutine update_chunks

end module simple_stream_chunk2D_utils
