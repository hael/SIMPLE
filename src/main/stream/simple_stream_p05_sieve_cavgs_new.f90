!@descr: task 5 in the stream pipeline: chunk-based 2D clustering and automatic selection of high-quality class averages (sieving)
!==============================================================================
! MODULE: simple_stream_p05_sieve_cavgs_new
!
! PURPOSE:
!   Implements stream pipeline task 5: continuously ingests incoming project
!   files, partitions their particles into micro-, midi-, and maxichunks for
!   parallel ab initio 2D classification, and automatically rejects low-quality
!   class averages at each tier (sieving).
!
! TYPES:
!   stream_p05_sieve_cavgs - commander_base extension; entry point for the
!                            sieve-cavgs stream task.
!
! WORKFLOW:
!   1. Initialise parameters, queue environment, and chunked2D object.
!   2. Enter infinite loop:
!      a. Watch dir_target for newly completed project files.
!      b. Import any new projects into the rec_list.
!      c. Generate microchunks from un-included records.
!      d. Poll running chunks and reject low-quality class averages.
!      e. Generate midi- and maxichunks from rejection-complete tiers.
!      f. Submit pending chunks to the queue.
!      g. Sleep for WAITTIME before the next cycle.
!
! PARAMETERS (hard-coded):
!   MAX_MOVIE_IMPORT    — maximum movies imported per loop cycle   (20)
!
! ENVIRONMENT:
!   SIMPLE_STREAM_CHUNK_PARTITION — queue partition for chunk jobs
!
! DEPENDENCIES:
!   simple_stream_api, simple_chunked2D, simple_stream_pool2D_utils
!==============================================================================
module simple_stream_p05_sieve_cavgs_new
  use simple_stream_api
  use simple_chunked2D,            only: chunked2D
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
  ! growing rec_list, and drives the chunked2D pipeline (generate, submit,
  ! collect, reject) on each loop cycle until a termination signal is detected.
  subroutine exec_stream_p05_sieve_cavgs( self, cline )
    class(stream_p05_sieve_cavgs), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline

    ! Hard-coded classification and import parameters
    integer, parameter :: MAX_MOVIE_IMPORT    = 20   ! max movies imported per cycle

    type(rec_list)               :: project_list
    type(parameters)             :: params
    type(qsys_env)               :: qenv
    type(sp_project)             :: spproj_glob
    type(chunked2D)              :: chunked_2D
    type(stream_watcher)         :: project_buff
    type(string), allocatable    :: projects(:)
    integer                      :: nprojects, n_mics_imported, n_ptcls_imported

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

    ! Clean up any previous run artefacts before starting
    call cleanup4restart()

    ! Create project file, folder structure, and initialise parameters
    call create_stream_project(spproj_glob, cline, string('sieve_cavgs'))
    call simple_mkdir(PATH_HERE // DIR_STREAM_COMPLETED)
    call params%new(cline)
    call init_stream_qenv(params, qenv, string(SIMPLE_STREAM_CHUNK_PARTITION))

    ! Sanity-check: must start from an empty project
    call spproj_glob%read(params%projfile)
    if( spproj_glob%os_mic%get_noris() /= 0 ) &
      THROW_HARD('stream_cluster2D must start from an empty project (e.g. from root project folder)')

    ! Initialise the project watcher and chunked2D pipeline
    project_buff = stream_watcher(LONGTIME, params%dir_target // '/' // DIR_STREAM_COMPLETED, &
                                  spproj=.true., nretries=10)
    call chunked_2D%new(params)

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

  contains

    ! Removes artefacts from a previous run in the same directory so the
    ! pipeline can restart cleanly. Handles both outdir and dir_exec restart
    ! modes; deletes termination flags, sigma and completed-project directories,
    ! and any chunk or set subdirectories found in the working directory.
    subroutine cleanup4restart
      type(string), allocatable :: folders(:)
      type(string)              :: cwd_restart, str_dir
      logical                   :: l_restart
      integer                   :: i
      call simple_getcwd(cwd_restart)
      l_restart = .false.
      if( cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir')) ) then
        l_restart = .true.
        call simple_chdir(cline%get_carg('outdir'))
      end if
      if( cline%defined('dir_exec') ) then
        if( .not. file_exists(cline%get_carg('dir_exec')) ) then
          str_dir = cline%get_carg('dir_exec')
          THROW_HARD('Previous directory does not exist: ' // str_dir%to_char())
          call str_dir%kill()
        end if
        l_restart = .true.
      end if
      if( l_restart ) then
        write(logfhandle,'(A,A)') '>>> RESTARTING EXISTING JOB IN ', cwd_restart%to_char()
        if( cline%defined('dir_exec') ) call cline%delete('dir_exec')
        call del_file(TERM_STREAM)
        call del_file(USER_PARAMS2D)
        call simple_rmdir(SIGMAS_DIR)
        call simple_rmdir(DIR_STREAM_COMPLETED)
        folders = simple_list_dirs('.')
        if( allocated(folders) ) then
          do i = 1, size(folders)
            if( folders(i)%has_substr(DIR_CHUNK) .or. folders(i)%has_substr(DIR_SET) ) &
              call simple_rmdir(folders(i))
          end do
        end if
      end if
      call simple_chdir(cwd_restart)
    end subroutine cleanup4restart

  end subroutine exec_stream_p05_sieve_cavgs

end module simple_stream_p05_sieve_cavgs_new