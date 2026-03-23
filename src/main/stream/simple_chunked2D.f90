!@descr: Manages the creation, configuration, and submission of pass-1/pass-2 microchunks and a single reference chunk
!==============================================================================
! MODULE: simple_chunked2D
!
! PURPOSE:
!   Manages the creation, configuration, and submission of microchunks across
!   two refinement passes and a final reference chunk — sub-projects each
!   containing a bounded particle subset drawn from a larger project list —
!   for parallel ab initio 2D class averaging. Each pass merges surviving
!   particles from the previous pass for a progressively higher-resolution
!   result.
!
! TYPES:
!   chunk2D    - Holds the identity, particle counts (total and selected),
!                folder path, project file path, command line, and lifecycle
!                state flags for a single chunk.
!   chunked2D  - Owns arrays of pass-1 and pass-2 microchunks and a single
!                reference chunk, together with the queue environment, output
!                directories, and shared processing parameters (threads, mask
!                diameter, parallelism).
!
! WORKFLOW:
!   1. new()                            — initialise a chunked2D object, creating
!                                         output directories and the queue
!                                         environment from a parameters object.
!   2. cycle()                          — convenience wrapper: collects completed
!                                         jobs and runs rejection, generates
!                                         chunks at all three tiers, then
!                                         submits pending chunks to the queue.
!   3. generate_microchunks_pass_1()    — partition un-included records from a
!                                         rec_list into pass-1 microchunks, each
!                                         capped at MICROCHUNK_P1_THRESHOLD
!                                         particles; write a project file per
!                                         chunk.
!   4. generate_microchunks_pass_2()    — merge rejection-complete pass-1
!                                         microchunks into pass-2 microchunks,
!                                         each capped at MICROCHUNK_P2_THRESHOLD
!                                         selected particles; write a project
!                                         file per chunk.
!   5. generate_refchunk()              — once pass-2 selected particles exceed
!                                         REFCHUNK_THRESHOLD, merge all eligible
!                                         pass-2 microchunks into the single
!                                         reference chunk; no-op if the reference
!                                         chunk already exists or the threshold
!                                         is not yet met.
!   6. submit()                         — dispatch pending chunks to the queue up
!                                         to the nparallel concurrency limit;
!                                         the reference chunk takes priority over
!                                         pass-2, which takes priority over
!                                         pass-1.
!   7. collect_and_reject()             — poll running chunks for the
!                                         ABINITIO2D_FINISHED sentinel; transition
!                                         completed chunks and immediately run
!                                         class-average rejection on each newly
!                                         finished chunk.
!
! CONSTANTS:
!   MICROCHUNK_P1_THRESHOLD — maximum particles per pass-1 microchunk       (5000)
!   MICROCHUNK_P2_THRESHOLD — maximum selected particles per pass-2 chunk   (5000)
!   REFCHUNK_THRESHOLD      — minimum selected particles to form a ref chunk(20000)
!   DEFAULT_NCLS            — default number of 2D classes                    (100)
!   DEFAULT_MICRO_P1_LP     — pass-1 microchunk low-pass filter cutoff, Å   (15.0)
!   DEFAULT_MICRO_P2_LP     — pass-2 microchunk low-pass filter cutoff, Å   (10.0)
!   DEFAULT_REF_LP          — reference chunk low-pass filter cutoff, Å      (6.0)
!   DEFAULT_WALLTIME        — per-chunk job time limit, seconds             (1740)
!
! ENVIRONMENT:
!   SIMPLE_CHUNK_PARTITION — if set, overrides the queue partition used for
!                            all chunk job submission.
!
! DEPENDENCIES:
!   simple_defs, simple_image, simple_string, simple_fileio, simple_syslib,
!   simple_cmdline, simple_qsys_env, simple_rec_list, simple_stack_io,
!   simple_parameters, simple_sp_project, simple_defs_fname,
!   simple_string_utils, simple_imgarr_utils, simple_projfile_utils,
!   simple_stream_microchunk_utils
!==============================================================================
module simple_chunked2D
  use simple_defs,         only: logfhandle, STDLEN, CWD_GLOB, COSMSKHALFWIDTH
  use simple_image,        only: image
  use simple_string,       only: string
  use simple_fileio,       only: swap_suffix
  use simple_syslib,       only: simple_mkdir, simple_abspath, simple_chdir, &
                                 simple_getcwd, file_exists, del_file, dir_exists
  use simple_cmdline,      only: cmdline
  use simple_qsys_env,     only: qsys_env
  use simple_rec_list,     only: rec_list, rec_iterator, project_rec
  use simple_stack_io,     only: stack_io
  use simple_parameters,   only: parameters
  use simple_sp_project,   only: sp_project
  use simple_defs_fname,   only: METADATA_EXT, ABINITIO2D_FINISHED, FRCS_FILE
  use simple_string_utils, only: int2str
  use simple_imgarr_utils, only: read_cavgs_into_imgarr, dealloc_imgarr
  use simple_projfile_utils,          only: merge_chunk_projfiles
  use simple_stream_microchunk_utils, only: calc_rejection_score, reject_outliers, reject_auto, reject_basic

  implicit none
  public  :: chunked2D
  private
#include "simple_local_flags.inc"

  integer, parameter :: MICROCHUNK_P1_THRESHOLD = 5000
  integer, parameter :: MICROCHUNK_P2_THRESHOLD = 5000
  integer, parameter :: REFCHUNK_THRESHOLD      = 10000
  integer, parameter :: DEFAULT_NCLS            = 100
  integer, parameter :: DEFAULT_WALLTIME        = 29 * 60  ! 29 minutes in seconds
  real,    parameter :: DEFAULT_MICRO_P1_LP     = 15.0
  real,    parameter :: DEFAULT_MICRO_P2_LP     = 10.0
  real,    parameter :: DEFAULT_REF_LP          = 10.0

  ! Labels used to route rejection strategy inside reject_cavgs
  character(len=*), parameter :: LABEL_PASS_1 = 'MICROCHUNK PASS 1'
  character(len=*), parameter :: LABEL_PASS_2 = 'MICROCHUNK PASS 2'
  character(len=*), parameter :: LABEL_REF    = 'REFCHUNK'

  type :: chunk2D
    private
    type(string)  :: projfile
    type(string)  :: folder
    type(cmdline) :: cline
    integer :: id                  = 0
    integer :: nptcls              = 0
    integer :: nptcls_selected     = 0
    logical :: abinitio2D_running  = .false.
    logical :: abinitio2D_complete = .false.
    logical :: rejection_complete  = .false.
    logical :: complete            = .false.
    logical :: failed              = .false.
  end type chunk2D

  type :: chunked2D
    private
    type(chunk2D), allocatable :: microchunks_pass_1(:)
    type(chunk2D), allocatable :: microchunks_pass_2(:)
    type(chunk2D)              :: refchunk
    type(qsys_env)             :: qenv
    type(string)               :: outdir_microchunks_pass_1
    type(string)               :: outdir_microchunks_pass_2
    type(string)               :: outdir_refchunk
    integer                    :: nparallel = 1
    integer                    :: nthr      = 1
    real                       :: mskdiam   = 0.0
  contains
    procedure :: new
    procedure :: kill
    procedure :: cycle
    procedure :: import_existing_microchunks_pass_1
    procedure :: import_existing_microchunks_pass_2
    procedure :: get_n_microchunks_pass_1
    procedure :: get_n_microchunks_pass_2
    procedure :: get_n_chunks_running
    procedure :: get_n_pass_1_non_rejected_ptcls
    procedure :: get_n_pass_2_non_rejected_ptcls
    procedure :: append_microchunk_pass_1
    procedure :: append_microchunk_pass_2
    procedure :: generate_microchunks_pass_1
    procedure :: generate_microchunk_pass_1_cline
    procedure :: generate_microchunks_pass_2
    procedure :: generate_microchunk_pass_2_cline
    procedure :: generate_refchunk
    procedure :: generate_refchunk_cline
    procedure :: submit
    procedure :: collect_and_reject
    procedure :: reject_cavgs
  end type chunked2D

contains

  ! Initialises a chunked2D object from a parameters object: derives output
  ! directories from the current working directory, stores the concurrency
  ! limit, mask diameter, and thread count; creates all three output
  ! directories; allocates empty microchunk arrays; and initialises the queue
  ! environment.
  subroutine new( self, params )
    class(chunked2D), intent(inout) :: self
    type(parameters), intent(in)    :: params
    type(string)                    :: cwd
    call self%kill()
    call simple_getcwd(cwd)
    self%nparallel                 = params%nchunks
    self%mskdiam                   = params%mskdiam
    self%nthr                      = params%nthr
    self%outdir_microchunks_pass_1 = string(cwd%to_char() // '/microchunks_pass_1')
    self%outdir_microchunks_pass_2 = string(cwd%to_char() // '/microchunks_pass_2')
    self%outdir_refchunk           = string(cwd%to_char() // '/refchunk')
    allocate(self%microchunks_pass_1(0))
    allocate(self%microchunks_pass_2(0))
    if( dir_exists(self%outdir_microchunks_pass_1) ) call self%import_existing_microchunks_pass_1()
    if( dir_exists(self%outdir_microchunks_pass_2) ) call self%import_existing_microchunks_pass_2()
    call simple_mkdir(self%outdir_microchunks_pass_1)
    call simple_mkdir(self%outdir_microchunks_pass_2)
    call simple_mkdir(self%outdir_refchunk)
    call self%qenv%new(params, 1, exec_bin=string('simple_exec'))
  end subroutine new

  ! Kills all live command lines for pass-1, pass-2, and the reference chunk,
  ! then deallocates the two microchunk arrays, releasing all associated memory.
  subroutine kill( self )
    class(chunked2D), intent(inout) :: self
    integer :: i
    do i = 1, self%get_n_microchunks_pass_1()
      call self%microchunks_pass_1(i)%cline%kill()
    end do
    do i = 1, self%get_n_microchunks_pass_2()
      call self%microchunks_pass_2(i)%cline%kill()
    end do
    call self%refchunk%cline%kill()
    if( allocated(self%microchunks_pass_1) ) deallocate(self%microchunks_pass_1)
    if( allocated(self%microchunks_pass_2) ) deallocate(self%microchunks_pass_2)
  end subroutine kill

  subroutine import_existing_microchunks_pass_1( self )
    class(chunked2D), intent(inout) :: self
    type(sp_project)                :: chunk_project
    type(chunk2D)                   :: new_chunk
    type(string)                    :: chunk_folder
    integer                         :: chunk_id
    chunk_id = 0
    do
      chunk_id = chunk_id + 1
      chunk_folder = string(self%outdir_microchunks_pass_1%to_char() // '/microchunk_pass_1_' // int2str(chunk_id))
      if(.not. dir_exists(chunk_folder)) exit
      new_chunk%id              = chunk_id
      new_chunk%folder          = simple_abspath(chunk_folder)
      new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                   '/microchunk_pass_1_' // int2str(chunk_id) // METADATA_EXT)
      call chunk_project%read(new_chunk%projfile)
      new_chunk%nptcls          = chunk_project%os_ptcl2D%get_noris()
      new_chunk%nptcls_selected = chunk_project%os_ptcl2D%count_state_gt_zero()
      new_chunk%abinitio2D_running  = .false.
      new_chunk%abinitio2D_complete = .true.
      new_chunk%rejection_complete  = .true.
      new_chunk%complete            = .true.
      new_chunk%failed              = .false.
      call chunk_project%kill()
      call self%append_microchunk_pass_1(new_chunk)
      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> IMPORTED EXISTING MICROCHUNK PASS 1 # ', chunk_id, ' WITH ', new_chunk%nptcls, ' PARTICLES'
    end do
  end subroutine import_existing_microchunks_pass_1

  subroutine import_existing_microchunks_pass_2( self )
    class(chunked2D), intent(inout) :: self
    type(sp_project)                :: chunk_project
    type(chunk2D)                   :: new_chunk
    type(string)                    :: chunk_folder
    integer                         :: chunk_id
    chunk_id = 0
    do
      chunk_id = chunk_id + 1
      chunk_folder = string(self%outdir_microchunks_pass_2%to_char() // '/microchunk_pass_2_' // int2str(chunk_id))
      if(.not. dir_exists(chunk_folder)) exit
      new_chunk%id              = chunk_id
      new_chunk%folder          = simple_abspath(chunk_folder)
      new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                   '/microchunk_pass_2_' // int2str(chunk_id) // METADATA_EXT)
      call chunk_project%read(new_chunk%projfile)
      new_chunk%nptcls          = chunk_project%os_ptcl2D%get_noris()
      new_chunk%nptcls_selected = chunk_project%os_ptcl2D%count_state_gt_zero()
      new_chunk%abinitio2D_running  = .false.
      new_chunk%abinitio2D_complete = .true.
      new_chunk%rejection_complete  = .true.
      new_chunk%complete            = .false.
      new_chunk%failed              = .false.
      call chunk_project%kill()
      call self%append_microchunk_pass_2(new_chunk)
      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> IMPORTED EXISTING MICROCHUNK PASS 2 # ', chunk_id, ' WITH ', new_chunk%nptcls, ' PARTICLES'
    end do
  end subroutine import_existing_microchunks_pass_2


  ! Convenience wrapper for the main processing loop: collects completed jobs
  ! and runs rejection, generates new chunks at all three tiers, then submits
  ! any pending chunks to the queue.
  subroutine cycle( self, project_list )
    class(chunked2D), intent(inout) :: self
    type(rec_list),   intent(inout) :: project_list
    call self%collect_and_reject()
    call self%generate_microchunks_pass_1(project_list)
    call self%generate_microchunks_pass_2()
    call self%generate_refchunk()
    call self%submit()
  end subroutine cycle

  ! Returns the total number of pass-1 microchunks; zero if unallocated.
  pure integer function get_n_microchunks_pass_1( self )
    class(chunked2D), intent(in) :: self
    get_n_microchunks_pass_1 = 0
    if( allocated(self%microchunks_pass_1) ) get_n_microchunks_pass_1 = size(self%microchunks_pass_1)
  end function get_n_microchunks_pass_1

  ! Returns the total number of pass-2 microchunks; zero if unallocated.
  pure integer function get_n_microchunks_pass_2( self )
    class(chunked2D), intent(in) :: self
    get_n_microchunks_pass_2 = 0
    if( allocated(self%microchunks_pass_2) ) get_n_microchunks_pass_2 = size(self%microchunks_pass_2)
  end function get_n_microchunks_pass_2

  ! Returns the combined number of currently running pass-1, pass-2, and
  ! reference chunks.
  pure integer function get_n_chunks_running( self )
    class(chunked2D), intent(in) :: self
    integer :: i
    get_n_chunks_running = 0
    do i = 1, self%get_n_microchunks_pass_1()
      if( self%microchunks_pass_1(i)%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
    end do
    do i = 1, self%get_n_microchunks_pass_2()
      if( self%microchunks_pass_2(i)%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
    end do
    if( self%refchunk%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
  end function get_n_chunks_running

  ! Returns the total selected-particle count across all pass-1 microchunks
  ! that have completed rejection but have not yet been consumed into a
  ! pass-2 microchunk.
  pure integer function get_n_pass_1_non_rejected_ptcls( self )
    class(chunked2D), intent(in) :: self
    integer :: i
    get_n_pass_1_non_rejected_ptcls = 0
    do i = 1, self%get_n_microchunks_pass_1()
      associate( chunk => self%microchunks_pass_1(i) )
        if( chunk%rejection_complete .and. .not. chunk%complete ) &
          get_n_pass_1_non_rejected_ptcls = &
            get_n_pass_1_non_rejected_ptcls + chunk%nptcls_selected
      end associate
    end do
  end function get_n_pass_1_non_rejected_ptcls

  ! Returns the total selected-particle count across all pass-2 microchunks
  ! that have completed rejection but have not yet been consumed into the
  ! reference chunk.
  pure integer function get_n_pass_2_non_rejected_ptcls( self )
    class(chunked2D), intent(in) :: self
    integer :: i
    get_n_pass_2_non_rejected_ptcls = 0
    do i = 1, self%get_n_microchunks_pass_2()
      associate( chunk => self%microchunks_pass_2(i) )
        if( chunk%rejection_complete .and. .not. chunk%complete ) &
          get_n_pass_2_non_rejected_ptcls = &
            get_n_pass_2_non_rejected_ptcls + chunk%nptcls_selected
      end associate
    end do
  end function get_n_pass_2_non_rejected_ptcls

  ! Appends a single chunk2D to the end of the pass-1 microchunk array.
  subroutine append_microchunk_pass_1( self, new_chunk )
    class(chunked2D), intent(inout) :: self
    type(chunk2D),    intent(in)    :: new_chunk
    self%microchunks_pass_1 = [self%microchunks_pass_1, new_chunk]
  end subroutine append_microchunk_pass_1

  ! Appends a single chunk2D to the end of the pass-2 microchunk array.
  subroutine append_microchunk_pass_2( self, new_chunk )
    class(chunked2D), intent(inout) :: self
    type(chunk2D),    intent(in)    :: new_chunk
    self%microchunks_pass_2 = [self%microchunks_pass_2, new_chunk]
  end subroutine append_microchunk_pass_2

  ! Partitions un-included records from project_list into pass-1 microchunks,
  ! each holding up to MICROCHUNK_P1_THRESHOLD particles. For each chunk:
  ! creates its subdirectory, accumulates micrographs until the threshold is
  ! reached, builds and writes a project file, applies any
  ! SIMPLE_CHUNK_PARTITION environment override, configures the command line,
  ! and marks the consumed records as included in project_list.
  subroutine generate_microchunks_pass_1( self, project_list )
    class(chunked2D), intent(inout) :: self
    type(rec_list),   intent(inout) :: project_list

    logical,          allocatable :: included(:)
    integer,          allocatable :: ids(:), nptcls(:)
    type(chunk2D)                 :: new_chunk
    type(sp_project)              :: chunk_project
    type(rec_list)                :: chunk_project_list
    type(string)                  :: chunk_folder
    character(len=STDLEN)         :: chunk_part_env
    integer                       :: imic, chunk_nptcls, envlen, chunk_id

    do while( project_list%get_nptcls_tot(l_not_included=.true.) > MICROCHUNK_P1_THRESHOLD )
      included = project_list%get_included_flags()
      ids      = pack(project_list%get_ids(),    .not. included)
      nptcls   = pack(project_list%get_nptcls(), .not. included)

      ! Accumulate micrographs until the particle threshold is reached
      chunk_nptcls = 0
      do imic = 1, size(nptcls)
        chunk_nptcls = chunk_nptcls + nptcls(imic)
        if( chunk_nptcls > MICROCHUNK_P1_THRESHOLD ) exit
      end do

      ! Create the chunk folder and populate chunk fields
      chunk_id                  = self%get_n_microchunks_pass_1() + 1
      chunk_folder              = string(self%outdir_microchunks_pass_1%to_char() // &
                                         '/microchunk_pass_1_' // int2str(chunk_id))
      call simple_mkdir(chunk_folder)
      new_chunk%id              = chunk_id
      new_chunk%nptcls          = chunk_nptcls
      new_chunk%nptcls_selected = chunk_nptcls
      new_chunk%folder          = simple_abspath(chunk_folder)
      new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                         '/microchunk_pass_1_' // int2str(chunk_id) // METADATA_EXT)
      call self%generate_microchunk_pass_1_cline(new_chunk)

      ! Build the project file from the sliced record list
      call project_list%slice(ids(1), ids(imic), chunk_project_list)
      call chunk_project%projrecords2proj(chunk_project_list)
      call chunk_project%update_projinfo(new_chunk%cline)
      call chunk_project%update_compenv(new_chunk%cline)

      ! Apply queue partition override if the environment variable is set
      call get_environment_variable('SIMPLE_CHUNK_PARTITION', chunk_part_env, envlen)
      if( envlen > 0 ) call chunk_project%compenv%set(1, 'qsys_partition', trim(chunk_part_env))

      call chunk_project%write(new_chunk%projfile)
      call chunk_project%kill()

      call self%append_microchunk_pass_1(new_chunk)
      call project_list%set_included_flags([ids(1), ids(imic)])

      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> MICROCHUNK PASS 1 # ', chunk_id, ' GENERATED WITH ', chunk_nptcls, ' PARTICLES'
    end do
  end subroutine generate_microchunks_pass_1

  ! Merges rejection-complete pass-1 microchunks into pass-2 microchunks,
  ! each accumulating up to MICROCHUNK_P2_THRESHOLD selected particles. For
  ! each pass-2 chunk: creates its subdirectory, collects eligible pass-1
  ! project files and particle counts, merges the project files, clears stale
  ! 2D results and intermediate files, applies any SIMPLE_CHUNK_PARTITION
  ! environment override, configures the command line, and marks the consumed
  ! pass-1 chunks as complete.
  subroutine generate_microchunks_pass_2( self )
    class(chunked2D), intent(inout) :: self

    type(string),     allocatable :: projfiles(:)
    type(chunk2D)                 :: new_chunk
    type(sp_project)              :: chunk_project
    type(string)                  :: chunk_folder
    character(len=STDLEN)         :: chunk_part_env
    integer                       :: i, chunk_nptcls, chunk_nptcls_selected, envlen, chunk_id

    do while( self%get_n_pass_1_non_rejected_ptcls() > MICROCHUNK_P2_THRESHOLD )
      write(*,*) "NON REJECTED", self%get_n_pass_1_non_rejected_ptcls()
      allocate(projfiles(0))
      chunk_nptcls          = 0
      chunk_nptcls_selected = 0

      ! Accumulate rejection-complete pass-1 chunks until threshold is reached
      do i = 1, self%get_n_microchunks_pass_1()
        associate( src => self%microchunks_pass_1(i) )
          if( src%rejection_complete .and. .not. src%complete ) then
            chunk_nptcls_selected = chunk_nptcls_selected + src%nptcls_selected
            chunk_nptcls          = chunk_nptcls          + src%nptcls
            projfiles             = [projfiles, src%projfile]
            src%complete          = .true.
          end if
          if( chunk_nptcls_selected > MICROCHUNK_P2_THRESHOLD ) exit
        end associate
      end do

      ! Create the pass-2 chunk folder and populate chunk fields
      chunk_id                  = self%get_n_microchunks_pass_2() + 1
      chunk_folder              = string(self%outdir_microchunks_pass_2%to_char() // &
                                         '/microchunk_pass_2_' // int2str(chunk_id))
      call simple_mkdir(chunk_folder)
      new_chunk%id              = chunk_id
      new_chunk%nptcls          = chunk_nptcls
      new_chunk%nptcls_selected = chunk_nptcls_selected
      new_chunk%folder          = simple_abspath(chunk_folder)
      new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                         '/microchunk_pass_2_' // int2str(chunk_id) // METADATA_EXT)
      call self%generate_microchunk_pass_2_cline(new_chunk)

      ! Merge source project files and clear stale 2D results
      call merge_and_clear(projfiles, chunk_folder, chunk_project, new_chunk%cline)

      call chunk_project%write(new_chunk%projfile)
      call chunk_project%kill()

      call self%append_microchunk_pass_2(new_chunk)
      deallocate(projfiles)

      write(logfhandle,'(A,I6,A,I8,A,I8,A)') '>>> MICROCHUNK PASS 2 # ', chunk_id, &
        ' GENERATED WITH ', chunk_nptcls_selected, '/', chunk_nptcls, ' PARTICLES'
    end do
  end subroutine generate_microchunks_pass_2

  ! Merges all rejection-complete pass-2 microchunks into the single reference
  ! chunk once their combined selected-particle count exceeds REFCHUNK_THRESHOLD.
  ! No-op if the reference chunk already exists (nptcls > 0) or the threshold
  ! is not yet met. Populates reference chunk fields, merges project files,
  ! clears stale 2D results and intermediate files, applies any
  ! SIMPLE_CHUNK_PARTITION environment override, configures the command line,
  ! and marks all consumed pass-2 chunks as complete.
  subroutine generate_refchunk( self )
    class(chunked2D), intent(inout) :: self

    type(string),     allocatable :: projfiles(:)
    type(sp_project)              :: chunk_project
    character(len=STDLEN)         :: chunk_part_env
    integer                       :: i, chunk_nptcls, chunk_nptcls_selected, envlen

    ! Only generate once, and only when enough selected particles are available
    if( self%refchunk%nptcls > 0 )                                    return
    if( self%get_n_pass_2_non_rejected_ptcls() < REFCHUNK_THRESHOLD ) return

    allocate(projfiles(0))
    chunk_nptcls          = 0
    chunk_nptcls_selected = 0

    ! Consume all rejection-complete pass-2 chunks
    do i = 1, self%get_n_microchunks_pass_2()
      associate( src => self%microchunks_pass_2(i) )
        if( src%rejection_complete .and. .not. src%complete ) then
          chunk_nptcls_selected = chunk_nptcls_selected + src%nptcls_selected
          chunk_nptcls          = chunk_nptcls          + src%nptcls
          projfiles             = [projfiles, src%projfile]
          src%complete          = .true.
        end if
      end associate
    end do

    ! Populate reference chunk fields
    self%refchunk%id              = 1
    self%refchunk%nptcls          = chunk_nptcls
    self%refchunk%nptcls_selected = chunk_nptcls_selected
    self%refchunk%folder          = simple_abspath(self%outdir_refchunk)
    self%refchunk%projfile        = string(self%refchunk%folder%to_char() // &
                                           '/refchunk' // METADATA_EXT)
    call self%generate_refchunk_cline(self%refchunk)

    ! Merge source project files and clear stale 2D results
    call merge_and_clear(projfiles, self%outdir_refchunk, chunk_project, self%refchunk%cline)

    call chunk_project%write(self%refchunk%projfile)
    call chunk_project%kill()
    deallocate(projfiles)

    write(logfhandle,'(A,I8,A,I8,A)') '>>> REFCHUNK GENERATED WITH ', &
      chunk_nptcls_selected, '/', chunk_nptcls, ' PARTICLES'
  end subroutine generate_refchunk

  ! Populates the abinitio2D command line for a pass-1 microchunk with:
  ! program name, project file and name, no-mkdir flag, thread count, mask
  ! diameter, class count, low-pass filter cutoff (DEFAULT_MICRO_P1_LP), and
  ! wall-time limit.
  subroutine generate_microchunk_pass_1_cline( self, new_chunk )
    class(chunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: new_chunk
    associate( cline => new_chunk%cline )
      call cline%set('prg',      'abinitio2D')
      call cline%set('projfile', new_chunk%projfile)
      call cline%set('projname', 'microchunk_pass_1')
      call cline%set('mkdir',    'no')
      call cline%set('nthr',     self%nthr)
      call cline%set('mskdiam',  self%mskdiam)
      call cline%set('ncls',     DEFAULT_NCLS)
      call cline%set('lp',       DEFAULT_MICRO_P1_LP)
      call cline%set('walltime', DEFAULT_WALLTIME)
    end associate
  end subroutine generate_microchunk_pass_1_cline

  ! Populates the abinitio2D command line for a pass-2 microchunk with:
  ! program name, project file and name, no-mkdir flag, thread count, mask
  ! diameter, class count, tighter low-pass filter cutoff (DEFAULT_MICRO_P2_LP),
  ! and wall-time limit.
  subroutine generate_microchunk_pass_2_cline( self, new_chunk )
    class(chunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: new_chunk
    associate( cline => new_chunk%cline )
      call cline%set('prg',      'abinitio2D')
      call cline%set('projfile', new_chunk%projfile)
      call cline%set('projname', 'microchunk_pass_2')
      call cline%set('mkdir',    'no')
      call cline%set('nthr',     self%nthr)
      call cline%set('mskdiam',  self%mskdiam)
      call cline%set('ncls',     DEFAULT_NCLS)
      call cline%set('lp',       DEFAULT_MICRO_P2_LP)
      call cline%set('walltime', DEFAULT_WALLTIME)
    end associate
  end subroutine generate_microchunk_pass_2_cline

  ! Populates the abinitio2D command line for the reference chunk with:
  ! program name, project file and name, no-mkdir flag, thread count, mask
  ! diameter, class count, tightest low-pass filter cutoff (DEFAULT_REF_LP),
  ! and wall-time limit.
  subroutine generate_refchunk_cline( self, new_chunk )
    class(chunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: new_chunk
    associate( cline => new_chunk%cline )
      call cline%set('prg',      'abinitio2D')
      call cline%set('projfile', new_chunk%projfile)
      call cline%set('projname', 'refchunk')
      call cline%set('mkdir',    'no')
      call cline%set('nthr',     self%nthr)
      call cline%set('mskdiam',  self%mskdiam)
      call cline%set('ncls',     DEFAULT_NCLS)
      call cline%set('lp',       DEFAULT_REF_LP)
      call cline%set('walltime', DEFAULT_WALLTIME)
    end associate
  end subroutine generate_refchunk_cline

  ! Submits pending chunks to the queue asynchronously, up to the nparallel
  ! concurrency limit. The reference chunk takes priority, followed by pass-2,
  ! then pass-1. Skips chunks that are already running or complete. For each
  ! submitted chunk: changes into its working directory, dispatches the job,
  ! marks it as running, and releases its command line. Restores the original
  ! working directory on completion.
  subroutine submit( self )
    class(chunked2D), intent(inout) :: self
    type(string) :: cwd
    integer      :: i
    call simple_getcwd(cwd)

    ! Submit reference chunk first — highest priority
    if( self%refchunk%nptcls > 0 ) then
      if( self%get_n_chunks_running() < self%nparallel ) then
        if( .not. (self%refchunk%abinitio2D_running .or. self%refchunk%abinitio2D_complete) ) then
          call simple_chdir(self%refchunk%folder)
          CWD_GLOB = self%refchunk%folder%to_char()
          call self%qenv%exec_simple_prg_in_queue_async( &
            self%refchunk%cline, string('./distr_microchunk'), string('simple_log_refchunk'))
          self%refchunk%abinitio2D_running = .true.
          call self%refchunk%cline%kill()
          write(logfhandle,'(A)') '>>> INITIATED 2D ANALYSIS OF REFCHUNK'
        end if
      end if
    end if

    ! Submit pass-2 microchunks second
    do i = 1, self%get_n_microchunks_pass_2()
      if( self%get_n_chunks_running() >= self%nparallel ) exit
      associate( chunk => self%microchunks_pass_2(i) )
        if( chunk%abinitio2D_running .or. chunk%abinitio2D_complete ) cycle
        call simple_chdir(chunk%folder)
        CWD_GLOB = chunk%folder%to_char()
        call self%qenv%exec_simple_prg_in_queue_async( &
          chunk%cline, string('./distr_microchunk'), string('simple_log_microchunk_pass_2'))
        chunk%abinitio2D_running = .true.
        call chunk%cline%kill()
        write(logfhandle,'(A,I6)') '>>> INITIATED 2D ANALYSIS OF MICROCHUNK PASS 2 # ', chunk%id
      end associate
    end do

    ! Submit pass-1 microchunks last — lowest priority
    do i = 1, self%get_n_microchunks_pass_1()
      if( self%get_n_chunks_running() >= self%nparallel ) exit
      associate( chunk => self%microchunks_pass_1(i) )
        if( chunk%abinitio2D_running .or. chunk%abinitio2D_complete ) cycle
        call simple_chdir(chunk%folder)
        CWD_GLOB = chunk%folder%to_char()
        call self%qenv%exec_simple_prg_in_queue_async( &
          chunk%cline, string('./distr_microchunk'), string('simple_log_microchunk_pass_1'))
        chunk%abinitio2D_running = .true.
        call chunk%cline%kill()
        write(logfhandle,'(A,I6)') '>>> INITIATED 2D ANALYSIS OF MICROCHUNK PASS 1 # ', chunk%id
      end associate
    end do

    call simple_chdir(cwd)
    CWD_GLOB = cwd%to_char()
  end subroutine submit

  ! Polls all running pass-1, pass-2, and reference chunks for the
  ! ABINITIO2D_FINISHED sentinel file. For each newly completed chunk:
  ! transitions it from running to complete and immediately runs class-average
  ! rejection. The reference chunk rejection is only attempted if it has been
  ! generated (nptcls > 0).
  subroutine collect_and_reject( self )
    class(chunked2D), intent(inout) :: self
    integer :: i

    do i = 1, self%get_n_microchunks_pass_1()
      associate( chunk => self%microchunks_pass_1(i) )
        if( chunk%abinitio2D_running ) then
          if( file_exists(chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED) ) then
            chunk%abinitio2D_running  = .false.
            chunk%abinitio2D_complete = .true.
            write(logfhandle,'(A,I6)') '>>> COMPLETED 2D ANALYSIS OF MICROCHUNK PASS 1 # ', chunk%id
          end if
        end if
        call self%reject_cavgs(chunk, string(LABEL_PASS_1))
      end associate
    end do

    do i = 1, self%get_n_microchunks_pass_2()
      associate( chunk => self%microchunks_pass_2(i) )
        if( chunk%abinitio2D_running ) then
          if( file_exists(chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED) ) then
            chunk%abinitio2D_running  = .false.
            chunk%abinitio2D_complete = .true.
            write(logfhandle,'(A,I6)') '>>> COMPLETED 2D ANALYSIS OF MICROCHUNK PASS 2 # ', chunk%id
          end if
        end if
        call self%reject_cavgs(chunk, string(LABEL_PASS_2))
      end associate
    end do

    if( self%refchunk%nptcls > 0 ) then
      if( self%refchunk%abinitio2D_running ) then
        if( file_exists(self%refchunk%folder%to_char() // '/' // ABINITIO2D_FINISHED) ) then
          self%refchunk%abinitio2D_running  = .false.
          self%refchunk%abinitio2D_complete = .true.
          write(logfhandle,'(A)') '>>> COMPLETED 2D ANALYSIS OF REFCHUNK'
        end if
      end if
      call self%reject_cavgs(self%refchunk, string(LABEL_REF))
    end if
  end subroutine collect_and_reject

  ! Performs rejection on a single chunk. Reads the class-average stack,
  ! computes the mask radius from image dimensions and mskdiam, runs outlier
  ! then auto (pass-1) or basic (pass-2 / reference) rejection, writes three
  ! filtered stacks (outlier-rejected, strategy-rejected, selected), propagates
  ! rejection flags into project states via map2ptcls_state, records the final
  ! selected particle count, and writes the project. No-op if the chunk is not
  ! yet complete or has already been rejected. Cleans up all allocations on exit.
  subroutine reject_cavgs( self, chunk, label )
    class(chunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: chunk
    type(string),     intent(in)    :: label

    type(image), allocatable :: cavg_imgs(:)
    logical,     allocatable :: l_rejected(:)
    integer,     allocatable :: states(:)
    type(sp_project)         :: spproj
    type(string)             :: stkname
    integer                  :: ncls, ldim(3)
    real                     :: smpd, mskrad

    if( .not. chunk%abinitio2D_complete ) return
    if( chunk%rejection_complete )        return

    ! Read project and class averages
    call spproj%read(chunk%projfile)
    call spproj%get_cavgs_stk(stkname, ncls, smpd)
    cavg_imgs = read_cavgs_into_imgarr(spproj)
    smpd      = cavg_imgs(1)%get_smpd()
    ldim      = cavg_imgs(1)%get_ldim()
    mskrad    = min(real(ldim(1) / 2) - COSMSKHALFWIDTH - 1., 0.5 * self%mskdiam / smpd)
    allocate(l_rejected(ncls), source=.false.)

    ! Outlier rejection applied at every tier
    call reject_outliers(cavg_imgs, mskrad, l_rejected)
    call write_cavgs(stkname, string('_rejected_outliers.mrc'), selected=.false.)

    ! Pass-1 uses auto rejection; pass-2 and reference use basic rejection
    if( label == string(LABEL_PASS_1) ) then
      call reject_auto(cavg_imgs, l_rejected)
      call write_cavgs(stkname, string('_rejected_auto.mrc'),  selected=.false.)
    else
      call reject_basic(spproj%os_cls2D, l_rejected)
      call write_cavgs(stkname, string('_rejected_basic.mrc'), selected=.false.)
    end if
    call write_cavgs(stkname, string('_selected.mrc'), selected=.true.)

    ! Propagate rejection flags into project states
    states = spproj%os_cls2D%get_all_asint('state')
    where( l_rejected ) states = 0
    call spproj%os_cls2D%set_all('state', states)
    call spproj%os_cls3D%set_all('state', states)
    call spproj%map2ptcls_state()
    call spproj%write()

    chunk%nptcls_selected    = spproj%os_ptcl2D%count_state_gt_zero()
    chunk%rejection_complete = .true.
    write(logfhandle,'(A,A,A,I6,A,I8,A,I8,A)') '>>> COMPLETED REJECTION FOR ', &
      label%to_char(), ' # ', chunk%id, ' : ', &
      chunk%nptcls_selected, '/', chunk%nptcls, ' PARTICLES SELECTED'

    ! Cleanup
    call dealloc_imgarr(cavg_imgs)
    call spproj%kill()
    deallocate(l_rejected, states)

  contains

    ! Writes a filtered subset of class averages to a new stack. If selected
    ! is true, writes only un-rejected classes; if false, writes only rejected
    ! ones. The output filename is derived by swapping the .mrc suffix of
    ! stkpath with the provided suffix string.
    subroutine write_cavgs( stkpath, suffix, selected )
      type(string), intent(in) :: stkpath, suffix
      logical,      intent(in) :: selected
      type(stack_io) :: stkio_r
      type(string)   :: out_stkname
      type(image)    :: img
      integer        :: istk, icls
      out_stkname = swap_suffix(stkpath, suffix, string('.mrc'))
      call stkio_r%open(stkpath, smpd, 'read', bufsz=1)
      istk = 0
      do icls = 1, ncls
        if( selected .eqv. (.not. l_rejected(icls)) ) then
          istk = istk + 1
          call img%new([ldim(1), ldim(2), 1], smpd, wthreads=.false.)
          call stkio_r%read(icls, img)
          call img%write(out_stkname, istk)
          call img%kill()
        end if
      end do
      call stkio_r%close()
    end subroutine write_cavgs

  end subroutine reject_cavgs

  ! ============================================================================
  ! MODULE-LEVEL HELPERS
  ! ============================================================================

  ! Merges a list of source project files into chunk_project within the given
  ! output directory, updates projinfo and compenv from cline, kills stale 2D
  ! orientation sets and output oris, and removes intermediate files (FRCs,
  ! class-average stacks, sigma2 star) carried over from the merge. Applies any
  ! SIMPLE_CHUNK_PARTITION environment override to the queue partition.
  subroutine merge_and_clear( projfiles, outdir, chunk_project, cline )
    type(string),    intent(in)    :: projfiles(:)
    type(string),    intent(in)    :: outdir
    type(sp_project), intent(inout) :: chunk_project
    type(cmdline),    intent(inout) :: cline
    character(len=STDLEN) :: chunk_part_env
    integer               :: envlen, i
    do i=1, size(projfiles)
      write(*,*) projfiles(i)%to_char()
    end do
    call merge_chunk_projfiles(projfiles, outdir, chunk_project, write_proj=.false.)
    call chunk_project%update_projinfo(cline)
    call chunk_project%update_compenv(cline)
    call chunk_project%os_cls2D%kill()
    call chunk_project%os_cls3D%kill()
    call chunk_project%os_out%kill()
    call del_file(string(outdir%to_char() // '/' // FRCS_FILE))
    call del_file(string(outdir%to_char() // '/cavgs.mrc'))
    call del_file(string(outdir%to_char() // '/cavgs_odd.mrc'))
    call del_file(string(outdir%to_char() // '/cavgs_even.mrc'))
    call del_file(string(outdir%to_char() // '/sigma2_combined.star'))
    call get_environment_variable('SIMPLE_CHUNK_PARTITION', chunk_part_env, envlen)
    if( envlen > 0 ) call chunk_project%compenv%set(1, 'qsys_partition', trim(chunk_part_env))
  end subroutine merge_and_clear

end module simple_chunked2D