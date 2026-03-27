!@descr: task 5 in the stream pipeline: multi-pass chunk-based 2D classification with automatic sieving of low-quality class averages
!==============================================================================
! MODULE: simple_stream_p05_sieve_cavgs_new
!
! PURPOSE:
!   Implements stream pipeline task 5: continuously ingests incoming project
!   files and drives a four-stage chunked2D classification pipeline that
!   produces progressively refined class averages, automatically rejecting
!   low-quality averages after each stage (sieving).
!
! TYPES:
!   stream_p05_sieve_cavgs - commander_base extension; entry point for the
!                            sieve-cavgs stream task.
!
! WORKFLOW:
!   1. Initialise parameters, queue environment, and chunked2D object.
!   2. Restore previously imported project history (if present).
!   3. Enter main loop (runs until termination signal):
!      a. Watch dir_target for newly completed project files.
!      b. Import new projects into the rec_list.
!      c. Call chunked_2D%cycle(), which performs per-cycle:
!           i.  collect_and_reject   — harvest completed chunks, sieve cavgs
!           ii. generate_microchunks_pass_1 — seed pass-1 chunks from new records
!           iii.generate_microchunks_pass_2 — promote pass-1 results to pass-2
!           iv. generate_refchunk          — build reference chunk from pass-2
!           v.  generate_microchunks_match — match-refine against reference
!           vi. submit                     — dispatch pending chunks to queue
!      d. Sleep for WAITTIME before the next cycle.
!
! PARAMETERS (hard-coded):
!   MAX_MOVIE_IMPORT    — maximum movies imported per loop cycle   (20)
!
! ENVIRONMENT:
!   SIMPLE_STREAM_CHUNK_PARTITION — queue partition for chunk jobs
!
! DEPENDENCIES:
!   simple_stream_api, simple_microchunked2D, simple_stream_pool2D_utils
!==============================================================================
module simple_stream_p05_sieve_cavgs_new
  use simple_stream_api
  use simple_fileio,               only: read_filetable
  use simple_microchunked2D,       only: microchunked2D
  use simple_stream_pool2D_utils,  only: set_lpthres_type

  implicit none
  public  :: stream_p05_sieve_cavgs
  private
#include "simple_local_flags.inc"

  type, extends(commander_base) :: stream_p05_sieve_cavgs
  contains
    procedure :: execute => exec_stream_p05_sieve_cavgs
  end type stream_p05_sieve_cavgs

contains

  ! Entry point for stream task 5. Continuously watches dir_target for
  ! completed project files, imports their micrographs and particles into a
  ! growing rec_list, and drives the four-stage chunked2D pipeline
  ! (collect/reject → pass-1 → pass-2 → refchunk → match → submit) on each
  ! loop cycle until a termination signal is detected.
  subroutine exec_stream_p05_sieve_cavgs( self, cline )
    class(stream_p05_sieve_cavgs), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline

    ! Hard-coded classification and import parameters
    integer, parameter :: MAX_MOVIE_IMPORT    = 20   ! max movies imported per cycle

    type(rec_list)               :: project_list
    type(parameters)             :: params
    type(qsys_env)               :: qenv
    type(sp_project)             :: spproj_glob
    type(microchunked2D)         :: chunked_2D
    type(stream_watcher)         :: project_buff
    type(string), allocatable    :: projects(:)
    integer                      :: nprojects, n_mics_imported, n_ptcls_imported, i

    ! Set fixed command-line defaults
    call cline%set('oritype',      'mic')
    call cline%set('mkdir',        'yes')
    call cline%set('autoscale',    'yes')
    call cline%set('reject_mics',  'no')
    call cline%set('prune',        'no')
    call cline%set('wiener',       'full')
    call cline%set('refine',       'snhc_smpl')
    call cline%set('ml_reg',       'no')
    call cline%set('objfun',       'euclid')
    call cline%set('sigma_est',    'global')
    call cline%set('numlen',       5)
    call cline%set('stream',       'yes')

    ! Apply user-overridable defaults
    if( .not. cline%defined('walltime')      ) call cline%set('walltime',      29 * 60)
    if( .not. cline%defined('nchunksperset') ) call cline%set('nchunksperset', 2)
    if( .not. cline%defined('remove_chunks') ) call cline%set('remove_chunks', 'yes')
    if( .not. cline%defined('center')        ) call cline%set('center',        'yes')
    if( .not. cline%defined('nmics')         ) call cline%set('nmics',         100)

    ! Create project file, folder structure, and initialise parameters
    call create_stream_project(spproj_glob, cline, string('sieve_cavgs'))
    call simple_mkdir(PATH_HERE // DIR_STREAM_COMPLETED)
    call params%new(cline)
    call init_stream_qenv(params, qenv, string(SIMPLE_STREAM_CHUNK_PARTITION))

    ! Sanity-check: must start from an empty project
    call spproj_glob%read(params%projfile)
    if( spproj_glob%os_mic%get_noris() /= 0 ) &
      THROW_HARD('stream_p05_sieve_cavgs must start from an empty project (e.g. from root project folder)')

    ! Initialise the project watcher and chunked2D pipeline
    project_buff = stream_watcher(LONGTIME, params%dir_target // '/' // DIR_STREAM_COMPLETED, &
                                  spproj=.true., nretries=10)

    ! Restore previously imported project history to avoid re-importing on restart
    if( file_exists('imported_projects.txt') ) then
      call read_filetable(string('imported_projects.txt'), projects)
      if( allocated(projects) ) then
        do i=1, size(projects)
          call project_buff%add2history(projects(i))
        end do
        deallocate(projects)
      endif
    endif
    call chunked_2D%new(params, string(PATH_HERE // DIR_STREAM_COMPLETED))

    n_mics_imported  = 0
    n_ptcls_imported = 0

    ! Main processing loop — runs until a termination signal is detected
    do
      ! Detect and import newly completed project files
      call project_buff%watch(nprojects, projects, max_nmovies=MAX_MOVIE_IMPORT)
      if( nprojects > 0 ) then
        call import_new_projects(project_list, projects, n_mics_imported, n_ptcls_imported)
        call project_buff%add2history(projects)
        write(logfhandle,'(A,I6,I9)') &
          '>>> # MICROGRAPHS / PARTICLES IMPORTED : ', n_mics_imported, n_ptcls_imported
      end if

      ! Drive the chunk pipeline for this cycle
      call chunked_2D%cycle(project_list)

      call sleep(WAITTIME)
    end do

    call chunked_2D%kill()

  end subroutine exec_stream_p05_sieve_cavgs

end module simple_stream_p05_sieve_cavgs_new