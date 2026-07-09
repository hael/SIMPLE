!@descr: multi-tier microchunk 2D classification driver (pass-1/pass-2)
!==============================================================================
! MODULE: simple_microchunked2D
!
! PURPOSE:
!   Manages creation, configuration, and submission of pass-1 and pass-2
!   microchunks (bounded particle subsets drawn from a larger project list)
!   for parallel ab initio 2D class averaging.
!
! TYPES:
!   chunk2D    - Holds the identity, particle counts (total and selected),
!                folder path, project file path, command line, and lifecycle
!                state flags (abinitio2D_running, abinitio2D_complete,
!                rejection_complete, complete) for a single chunk.
!   microchunked2D - Owns arrays of pass-1 and pass-2 microchunks, queue
!                    environment, output directories, and shared processing
!                    parameters.
!
! WORKFLOW:
!   1. new()                            — initialise a microchunked2D object: derives
!                                         output directories, imports any
!                                         existing chunks from previous runs,
!                                         creates directories, and initialises
!                                         the queue environment.
!   2. cycle()                          — convenience wrapper: collects completed
!                                         jobs and runs rejection, generates
!                                         chunks at all active tiers, then submits
!                                         pending chunks to the queue.
!   3. generate_chunks_coarse()    — partition un-included records from a
!                                         rec_list into pass-1 microchunks, each
!                                         capped at DEFAULT_COARSE_POP_THRESHOLD
!                                         particles; write a project file per
!                                         chunk and update the imported-projects
!                                         file table.
!   4. generate_chunks_fine()    — merge rejection-complete pass-1
!                                         microchunks into pass-2 microchunks.
!   5. submit()                         — dispatch pending chunks to the queue up
!                                         to the nparallel concurrency limit;
!                                         priority order: pass-2 > pass-1.
!   6. combine_completed_pass_1_chunks() — merge the project files of all
!                                         rejection-complete pass-1 chunks into
!                                         a single combined project file.
!                                         No-op if no pass-1 chunks are
!                                         rejection-complete or the combined
!                                         file already exists.
!   7. collect_and_reject()             — poll running chunks for the
!                                         ABINITIO2D_FINISHED sentinel; transition
!                                         completed chunks and immediately run
!                                         class-average rejection on each newly
!                                         finished chunk.
!
! SENTINEL FILES (written to each chunk folder):
!   ABINITIO2D_FINISHED  — written by the queue job on completion.
!   REJECTION_FINISHED   — written by reject_cavgs on successful rejection.
!   COMPLETE             — written when a chunk is consumed into the next tier.
!
! CONSTANTS:
!   DEFAULT_COARSE_POP_THRESHOLD — maximum particles per pass-1 microchunk       (5000)
!   DEFAULT_FINE_POP_THRESHOLD — maximum selected particles per pass-2 chunk  (10000)
!   DEFAULT_NCLS                          — default number of 2D classes         (100)
!   DEFAULT_LPSTART                       — shared low-pass start cutoff, Å    (15.0)
!   DEFAULT_MICRO_P1_LP                   — pass-1 low-pass stop cutoff, Å    (15.0)
!   DEFAULT_MICRO_P2_LP                   — pass-2 low-pass stop cutoff, Å    (10.0)
!   DEFAULT_WALLTIME                      — per-chunk job time limit, s  (1740 / 29 min)
!
! ENVIRONMENT:
!   SIMPLE_CHUNK_PARTITION — if set, overrides the queue partition used for
!                            all chunk job submission.
!
! DEPENDENCIES:
!   simple_defs, simple_image, simple_timer, simple_string, simple_fileio,
!   simple_syslib, simple_cmdline, simple_qsys_env, simple_rec_list,
!   simple_gui_utils, simple_parameters, simple_sp_project, simple_defs_fname,
!   simple_string_utils, simple_imgarr_utils, simple_projfile_utils,
!   simple_cluster2D_rejector
!==============================================================================
module simple_ptcl_sieve
  use unix,                               only: c_time, c_long
  use simple_defs,                        only: logfhandle, STDLEN, CWD_GLOB
  use simple_error,                       only: simple_exception
  use simple_image,                       only: image
  use simple_timer,                       only: timer_int_kind, tic, toc
  use simple_fileio,                      only: swap_suffix, simple_copy_file, write_filetable, simple_touch, basename
  use simple_string,                      only: string
  use simple_syslib,                      only: simple_mkdir, simple_abspath, simple_chdir, &
                                                simple_getcwd, file_exists, del_file, dir_exists
  use simple_cmdline,                     only: cmdline
  use simple_qsys_env,                    only: qsys_env
  use simple_rec_list,                    only: rec_list
  use simple_gui_utils,                   only: mrc2jpeg_tiled
  use simple_defs_fname,                  only: METADATA_EXT, ABINITIO2D_FINISHED, FRCS_FILE, JPG_EXT, MRC_EXT
  use simple_parameters,                  only: parameters
  use simple_sp_project,                  only: sp_project
  use simple_imgarr_utils,                only: read_cavgs_into_imgarr, dealloc_imgarr, write_imgarr
  use simple_string_utils,                only: int2str
  use simple_projfile_utils,              only: merge_chunk_projfiles
  use simple_cluster2D_rejector,          only: cluster2D_rejector
  use simple_cavg_quality_types,          only: cavg_quality_result, CAVG_QUALITY_CONTEXT_SIEVE
  use simple_class_compatibility,         only: class_compatibility, support_model_metrics
  use simple_cavg_quality_helpers,        only: cavg_rejection_reason_string
  use simple_cavg_quality_analysis,       only: evaluate_cavg_quality_hard_reject

  implicit none
  public  :: ptcl_sieve
  private
#include "simple_local_flags.inc"

  logical, parameter :: DEBUG                                 = .false.
  integer, parameter :: DEFAULT_COARSE_POP_THRESHOLD          = 5000
  integer, parameter :: DEFAULT_FINE_POP_THRESHOLD            = 10000
  integer, parameter :: DEFAULT_NCLS                          = 100
  integer, parameter :: DEFAULT_WALLTIME                      = 29 * 60  ! 29 minutes in seconds
  integer, parameter :: DEFAULT_COARSE_BOX                    = 128
  integer, parameter :: DEFAULT_FINE_BOX                      = 128
  integer, parameter :: DEFAULT_COARSE_NSAMPLE                = 2000
  integer, parameter :: DEFAULT_FINE_NSAMPLE                  = 2000
  integer, parameter :: COARSE_COMPATIBILITY_TRAINING_THRESHOLD = 10
  integer, parameter :: FINE_COMPATIBILITY_TRAINING_THRESHOLD = 10
  real,    parameter :: DEFAULT_LPSTART                       = 15.0
  real,    parameter :: DEFAULT_COARSE_LP                     = 15.0
  real,    parameter :: DEFAULT_FINE_LP                       = 10.0

  ! Labels used to route rejection strategy inside reject_cavgs
  character(len=*), parameter :: LABEL_COARSE     = 'COARSE CHUNK'
  character(len=*), parameter :: LABEL_FINE       = 'FINE CHUNK'
  character(len=*), parameter :: REJECTION_FAILED = 'REJECTION_FAILED'

  type :: chunk2D_state
    private
    type(string)  :: projfile
    type(string)  :: folder
    type(cmdline) :: cline
    integer       :: id                  = 0
    integer       :: nptcls              = 0
    integer       :: nptcls_selected     = 0
    logical       :: abinitio2D_running  = .false.
    logical       :: abinitio2D_complete = .false.
    logical       :: rejection_complete  = .false.
    logical       :: complete            = .false.
    logical       :: failed              = .false.
  end type chunk2D_state

  type :: chunk2D_coarse_defaults
    integer :: pop_threshold = DEFAULT_COARSE_POP_THRESHOLD
    integer :: ncls          = DEFAULT_NCLS
    integer :: nsample       = DEFAULT_COARSE_NSAMPLE
    integer :: box_crop      = DEFAULT_COARSE_BOX
    real    :: lpstart       = DEFAULT_LPSTART
    real    :: lpstop        = DEFAULT_COARSE_LP
  end type chunk2D_coarse_defaults

  type :: chunk2D_fine_defaults
    integer :: pop_threshold = DEFAULT_FINE_POP_THRESHOLD
    integer :: ncls          = DEFAULT_NCLS
    integer :: nsample       = DEFAULT_FINE_NSAMPLE
    integer :: box_crop      = DEFAULT_FINE_BOX
    real    :: lpstart       = DEFAULT_LPSTART
    real    :: lpstop        = DEFAULT_FINE_LP
  end type chunk2D_fine_defaults

  type :: ptcl_sieve
    private
    type(chunk2D_state), allocatable  :: chunks_coarse(:)
    type(chunk2D_state), allocatable  :: chunks_fine(:)
    integer,             allocatable  :: latest_jpeg_inds(:)
    integer,             allocatable  :: latest_jpeg_pops(:)
    integer,             allocatable  :: latest_jpeg_selection(:)
    real,                allocatable  :: latest_jpeg_res(:)
    type(class_compatibility)         :: coarse_compatibility_model
    type(class_compatibility)         :: fine_compatibility_model
    type(chunk2D_coarse_defaults)     :: coarse_defaults
    type(chunk2D_fine_defaults)       :: fine_defaults
    type(qsys_env)                    :: qenv
    type(string)                      :: outdir_chunks_coarse
    type(string)                      :: outdir_chunks_fine
    type(string)                      :: completedir
    type(string)                      :: latest_jpeg
    type(string)                      :: latest_stkname
    logical                           :: coarse_only        = .false.
    logical                           :: pre_chunked        = .false.
    logical                           :: final_ingestion    = .false.
    logical                           :: compatibility_training_pass_1_complete = .false.
    logical                           :: compatibility_training_pass_2_complete = .false.
    integer                           :: compatibility_training_pass_1_count    = 0
    integer                           :: compatibility_training_pass_2_count    = 0
    integer                           :: nparallel          = 1
    integer                           :: nthr               = 1
    integer                           :: nparts             = 1
    integer                           :: n_accepted_ptcls   = 0
    integer                           :: n_accepted_mics    = 0
    integer                           :: n_rejected_ptcls   = 0
    integer                           :: last_import        = 0
    integer                           :: latest_jpeg_xtiles = 0
    integer                           :: latest_jpeg_ytiles = 0
    real                              :: mskdiam            = 0.0
  contains
    procedure :: new
    procedure :: kill
    procedure :: cycle
    procedure :: import_existing_chunks_coarse
    procedure :: import_existing_chunks_fine
    procedure :: get_n_chunks_coarse
    procedure :: get_n_chunks_fine
    procedure :: get_n_chunks_running
    procedure :: get_n_pass_1_non_rejected_ptcls
    procedure :: get_n_pass_2_non_rejected_ptcls
    procedure :: get_n_accepted_ptcls
    procedure :: get_n_accepted_micrographs
    procedure :: get_n_rejected_ptcls
    procedure :: get_n_total_particles
    procedure :: get_latest
    procedure :: get_finished
    procedure :: set_final_ingestion
    procedure :: append_chunk_coarse
    procedure :: append_chunk_fine
    procedure :: generate_chunks_coarse
    procedure :: generate_chunk_coarse_cline
    procedure :: generate_chunks_fine
    procedure :: generate_chunk_fine_cline
    procedure :: combine_completed_chunks
    procedure :: submit
    procedure :: collect_and_reject
    procedure :: reject_cavgs
  end type ptcl_sieve

contains

  ! --------------------------------------------------------------------------
  ! LIFECYCLE
  ! --------------------------------------------------------------------------

  ! Initialises a microchunked2D object from a parameters object: derives output
  ! directories from the current working directory, stores the concurrency
  ! limit, mask diameter, and thread count; imports any existing chunks from a
  ! previous run; creates all four output directories; allocates empty chunk
  ! arrays; and initialises the queue environment.
  subroutine new( self, params, completedir, coarse_only, pre_chunked, compatibility_refs )
    class(ptcl_sieve), intent(inout) :: self
    type(parameters),           intent(in)    :: params
    type(string),               intent(in)    :: completedir
    type(string),     optional, intent(in)    :: compatibility_refs
    logical,          optional, intent(in)    :: coarse_only
    logical,          optional, intent(in)    :: pre_chunked
    type(string)                              :: cwd
    integer(timer_int_kind)                   :: t0
    t0 = timer_start()
    call self%kill()
    call simple_getcwd(cwd)
    self%completedir = completedir
    self%nparallel   = params%nchunks
    self%mskdiam     = params%mskdiam
    self%nthr        = params%nthr
    self%nparts      = params%nparts
    if( present(coarse_only) ) self%coarse_only = coarse_only
    if( present(pre_chunked) ) self%pre_chunked = pre_chunked
    self%outdir_chunks_coarse = string(cwd%to_char() // '/chunks_coarse')
    self%outdir_chunks_fine   = string(cwd%to_char() // '/chunks_fine')
    allocate(self%chunks_coarse(0))
    allocate(self%chunks_fine(0))
    ! Import existing chunks before creating directories.
    if( dir_exists(self%outdir_chunks_coarse)      ) call self%import_existing_chunks_coarse()
    if( dir_exists(self%outdir_chunks_fine) ) call self%import_existing_chunks_fine()
    call simple_mkdir(self%outdir_chunks_coarse)
    call simple_mkdir(self%outdir_chunks_fine)
    call self%qenv%new(params, 1, exec_bin=string('simple_exec'))
    call self%coarse_compatibility_model%new()
    call self%fine_compatibility_model%new()
    if( present(compatibility_refs) ) then
      if( file_exists(compatibility_refs) ) then
        call self%coarse_compatibility_model%train(compatibility_refs)
        call self%fine_compatibility_model%train(compatibility_refs)
      end if
    end if
    call timer_stop(t0, string('new'))
  end subroutine new

  ! Kills all live command lines for pass-1 and pass-2 chunks, deallocates
  ! all microchunk arrays,
  ! and resets all scalar state so the object can be safely reused via new().
  subroutine kill( self )
    class(ptcl_sieve), intent(inout) :: self
    integer(timer_int_kind) :: t0
    integer                 :: i
    t0 = timer_start()
    if( allocated(self%chunks_coarse) ) then
      do i = 1, self%get_n_chunks_coarse()
        call self%chunks_coarse(i)%cline%kill()
      end do
      deallocate(self%chunks_coarse)
    end if
    if( allocated(self%chunks_fine) ) then
      do i = 1, self%get_n_chunks_fine()
        call self%chunks_fine(i)%cline%kill()
      end do
      deallocate(self%chunks_fine)
    end if
    call self%coarse_compatibility_model%kill()
    call self%fine_compatibility_model%kill()
    self%compatibility_training_pass_1_complete = .false.
    self%compatibility_training_pass_2_complete = .false.
    self%compatibility_training_pass_1_count    = 0
    self%compatibility_training_pass_2_count    = 0
    self%n_accepted_ptcls    = 0
    self%n_accepted_mics     = 0
    self%n_rejected_ptcls    = 0
    self%coarse_only         = .false.
    self%pre_chunked         = .false.
    self%final_ingestion     = .false.
    call timer_stop(t0, string('kill'))
  end subroutine kill

  ! Sets the module-level final-ingestion flag. Once enabled, pass-2 final
  ! flushing logic is allowed to run without waiting on timeout.
  subroutine set_final_ingestion(self)
    class(ptcl_sieve), intent(inout) :: self
    self%final_ingestion = .true.
  end subroutine set_final_ingestion

  ! --------------------------------------------------------------------------
  ! IMPORT FROM PREVIOUS RUN
  ! --------------------------------------------------------------------------

  ! Scans the pass-1 output directory for existing chunk subdirectories and
  ! populates the pass-1 array with chunk records whose state is inferred from
  ! sentinel files (ABINITIO2D_FINISHED, REJECTION_FINISHED, COMPLETE).
  ! Regenerates the cline for any chunk that has not yet completed abinitio2D,
  ! so that interrupted chunks can be resubmitted after a restart.
  subroutine import_existing_chunks_coarse( self )
    class(ptcl_sieve), intent(inout) :: self
    type(sp_project)                     :: chunk_project
    type(chunk2D_state)                  :: new_chunk
    type(string)                         :: chunk_folder
    integer(timer_int_kind)              :: t0
    integer                              :: chunk_id
    t0       = timer_start()
    chunk_id = 0
    do
      chunk_id     = chunk_id + 1
      chunk_folder = string(self%outdir_chunks_coarse%to_char() // &
                            '/chunk_coarse_' // int2str(chunk_id))
      if( .not. dir_exists(chunk_folder) ) exit
      new_chunk%id       = chunk_id
      new_chunk%folder   = simple_abspath(chunk_folder)
      new_chunk%projfile = string(new_chunk%folder%to_char() // &
                                  '/chunk_coarse_' // int2str(chunk_id) // METADATA_EXT)
      if( .not. file_exists(new_chunk%projfile) ) then
        write(logfhandle,'(A,I6,A)') '>>> WARNING: missing project file for coarse chunk #', chunk_id, ', skipping'
        cycle
      end if
      call chunk_project%read(new_chunk%projfile)
      new_chunk%nptcls              = chunk_project%os_ptcl2D%get_noris()
      new_chunk%nptcls_selected     = chunk_project%os_ptcl2D%count_state_gt_zero()
      new_chunk%abinitio2D_running  = .false.
      new_chunk%abinitio2D_complete = file_exists(new_chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED)
      new_chunk%failed              = file_exists(new_chunk%folder%to_char() // '/' // REJECTION_FAILED)
      new_chunk%rejection_complete  = file_exists(new_chunk%folder%to_char() // '/REJECTION_FINISHED') .or. new_chunk%failed
      new_chunk%complete            = file_exists(new_chunk%folder%to_char() // '/COMPLETE') .or. new_chunk%failed
      call chunk_project%kill()
      if( .not. new_chunk%abinitio2D_complete .and. .not. new_chunk%failed ) call self%generate_chunk_coarse_cline(new_chunk, new_chunk%nptcls_selected)
      call self%append_chunk_coarse(new_chunk)
      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> IMPORTED EXISTING COARSE CHUNK # ', chunk_id, ' WITH ', new_chunk%nptcls, ' PARTICLES'
    end do
    call timer_stop(t0, string('import_existing_chunks_coarse'))
  end subroutine import_existing_chunks_coarse

  ! Scans the pass-2 output directory for existing chunk subdirectories and
  ! populates the pass-2 array with chunk records whose state is inferred from
  ! sentinel files. Regenerates the cline for any chunk that has not yet
  ! completed abinitio2D, so that interrupted chunks can be resubmitted.
  subroutine import_existing_chunks_fine( self )
    class(ptcl_sieve), intent(inout) :: self
    type(sp_project)                 :: chunk_project
    type(chunk2D_state)              :: new_chunk
    type(string)                     :: chunk_folder
    integer(timer_int_kind)          :: t0
    integer                          :: chunk_id
    t0       = timer_start()
    chunk_id = 0
    do
      chunk_id     = chunk_id + 1
      chunk_folder = string(self%outdir_chunks_fine%to_char() // &
                            '/chunk_fine_' // int2str(chunk_id))
      if( .not. dir_exists(chunk_folder) ) exit
      new_chunk%id       = chunk_id
      new_chunk%folder   = simple_abspath(chunk_folder)
      new_chunk%projfile = string(new_chunk%folder%to_char() // &
                                  '/chunk_fine_' // int2str(chunk_id) // METADATA_EXT)
      if( .not. file_exists(new_chunk%projfile) ) then
        write(logfhandle,'(A,I6,A)') '>>> WARNING: missing project file for fine chunk #', chunk_id, ', skipping'
        cycle
      end if
      call chunk_project%read(new_chunk%projfile)
      new_chunk%nptcls              = chunk_project%os_ptcl2D%get_noris()
      new_chunk%nptcls_selected     = chunk_project%os_ptcl2D%count_state_gt_zero()
      new_chunk%abinitio2D_running  = .false.
      new_chunk%abinitio2D_complete = file_exists(new_chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED)
      new_chunk%failed              = file_exists(new_chunk%folder%to_char() // '/' // REJECTION_FAILED)
      new_chunk%rejection_complete  = file_exists(new_chunk%folder%to_char() // '/REJECTION_FINISHED') .or. new_chunk%failed
      new_chunk%complete            = file_exists(new_chunk%folder%to_char() // '/COMPLETE') .or. new_chunk%failed
      call chunk_project%kill()
      if( .not. new_chunk%abinitio2D_complete .and. .not. new_chunk%failed ) call self%generate_chunk_fine_cline(new_chunk, new_chunk%nptcls_selected)
      call self%append_chunk_fine(new_chunk)
      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> IMPORTED EXISTING FINE CHUNK # ', chunk_id, ' WITH ', new_chunk%nptcls, ' PARTICLES'
    end do
    call timer_stop(t0, string('import_existing_chunks_fine'))
  end subroutine import_existing_chunks_fine

  ! --------------------------------------------------------------------------
  ! MAIN LOOP WRAPPER
  ! --------------------------------------------------------------------------

  ! Convenience wrapper for the main processing loop: collects completed jobs
  ! and runs rejection, generates new chunks at all tiers in tier order, then
  ! submits any pending chunks to the queue.
  subroutine cycle( self, project_list )
    class(ptcl_sieve), intent(inout) :: self
    type(rec_list),        intent(inout) :: project_list
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    call self%collect_and_reject()
    call self%generate_chunks_coarse(project_list)
    if( .not. self%coarse_only ) call self%generate_chunks_fine()
    call self%submit()
    call timer_stop(t0, string('cycle'))
  end subroutine cycle

  ! --------------------------------------------------------------------------
  ! QUERIES
  ! --------------------------------------------------------------------------

  ! Returns the total number of pass-1 microchunks; zero if unallocated.
  pure integer function get_n_chunks_coarse( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_chunks_coarse = 0
    if( allocated(self%chunks_coarse) ) get_n_chunks_coarse = size(self%chunks_coarse)
  end function get_n_chunks_coarse

  ! Returns the total number of pass-2 microchunks; zero if unallocated.
  pure integer function get_n_chunks_fine( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_chunks_fine = 0
    if( allocated(self%chunks_fine) ) get_n_chunks_fine = size(self%chunks_fine)
  end function get_n_chunks_fine

  ! Returns the combined number of currently running pass-1 and pass-2 chunks.
  pure integer function get_n_chunks_running( self )
    class(ptcl_sieve), intent(in) :: self
    integer :: i
    get_n_chunks_running = 0
    do i = 1, self%get_n_chunks_coarse()
      if( self%chunks_coarse(i)%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
    end do
    do i = 1, self%get_n_chunks_fine()
      if( self%chunks_fine(i)%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
    end do
  end function get_n_chunks_running

  ! Returns the total selected-particle count across all pass-1 microchunks
  ! that have completed rejection but have not yet been consumed into a
  ! pass-2 microchunk.
  pure integer function get_n_pass_1_non_rejected_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    integer :: i
    get_n_pass_1_non_rejected_ptcls = 0
    do i = 1, self%get_n_chunks_coarse()
      associate( chunk => self%chunks_coarse(i) )
        if( chunk%rejection_complete .and. .not. chunk%complete .and. .not. chunk%failed ) &
          get_n_pass_1_non_rejected_ptcls = &
            get_n_pass_1_non_rejected_ptcls + chunk%nptcls_selected
      end associate
    end do
  end function get_n_pass_1_non_rejected_ptcls

  ! Returns the total selected-particle count across all pass-2 microchunks
  ! that have completed rejection but are not yet finalised.
  pure integer function get_n_pass_2_non_rejected_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    integer :: i
    get_n_pass_2_non_rejected_ptcls = 0
    do i = 1, self%get_n_chunks_fine()
      associate( chunk => self%chunks_fine(i) )
        if( chunk%rejection_complete .and. .not. chunk%complete .and. .not. chunk%failed ) &
          get_n_pass_2_non_rejected_ptcls = &
            get_n_pass_2_non_rejected_ptcls + chunk%nptcls_selected
      end associate
    end do
  end function get_n_pass_2_non_rejected_ptcls

  ! Returns the cumulative number of particles accepted (state > 0) across all
  ! finalised pass-2 chunks. Updated by collect_and_reject when each pass-2 chunk
  ! is marked complete. Pass-1 rejections are not included; use
  ! get_n_pass_1/2_non_rejected_ptcls for unconsumed per-chunk selected counts.
  pure integer function get_n_accepted_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_accepted_ptcls = self%n_accepted_ptcls
  end function get_n_accepted_ptcls

  ! Returns the cumulative number of accepted micrographs across all
  ! finalised pass-2 chunks.
  pure integer function get_n_accepted_micrographs( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_accepted_micrographs = self%n_accepted_mics
  end function get_n_accepted_micrographs

  ! Returns the cumulative number of particles rejected (state = 0) across all
  ! finalised pass-2 chunks. Updated alongside n_accepted_ptcls by
  ! collect_and_reject when each pass-2 chunk is marked complete.
  pure integer function get_n_rejected_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_rejected_ptcls = self%n_rejected_ptcls
  end function get_n_rejected_ptcls

  ! Returns the cumulative total number of particles in finalised chunks,
  ! computed as accepted + rejected.
  pure integer function get_n_total_particles( self )
    class(ptcl_sieve), intent(in)  :: self
    get_n_total_particles = self%n_accepted_ptcls + self%n_rejected_ptcls
  end function get_n_total_particles

  logical function get_latest( self, jpeg_inds, jpeg_pops, jpeg_res, jpeg, stk, xtiles, ytiles, selection )
    class(ptcl_sieve), intent(in)    :: self
    integer,      allocatable,  intent(inout) :: jpeg_inds(:), jpeg_pops(:), selection(:)  ! class indices and populations in the JPEG tile order
    real,         allocatable,  intent(inout) :: jpeg_res(:)                               ! resolution (Angstrom) per class, same order
    type(string),               intent(out)   :: jpeg, stk                                 ! paths to the JPEG contact sheet and MRC reference stack
    integer,                    intent(out)   :: xtiles, ytiles                            ! tile grid dimensions of the JPEG contact sheet
    ! release any prior allocations so the caller gets a clean result on .false. return
    if( allocated(jpeg_inds)  ) deallocate(jpeg_inds)
    if( allocated(jpeg_pops)  ) deallocate(jpeg_pops)
    if( allocated(jpeg_res)   ) deallocate(jpeg_res)
    if( allocated(selection)  ) deallocate(selection) 
    jpeg   = ''
    stk    = ''
    xtiles = 0
    ytiles = 0
    if( self%latest_jpeg%strlen() == 0 ) then
      get_latest = .false.
      return
    end if
    ! guard against inconsistent state where latest_jpeg is set but the
    ! latest arrays have not yet been populated
    if( .not. (allocated(self%latest_jpeg_inds) .and. allocated(self%latest_jpeg_pops) .and. allocated(self%latest_jpeg_res) .and. allocated(self%latest_jpeg_selection)) ) then
      get_latest = .false.
      return
    end if
    allocate(jpeg_inds, source=self%latest_jpeg_inds)
    allocate(jpeg_pops, source=self%latest_jpeg_pops)
    allocate(jpeg_res,  source=self%latest_jpeg_res )
    allocate(selection, source=self%latest_jpeg_selection)
    stk    = self%latest_stkname
    jpeg   = self%latest_jpeg
    xtiles = self%latest_jpeg_xtiles
    ytiles = self%latest_jpeg_ytiles
    get_latest = .true.
  end function get_latest

  ! Returns true when all required tiers are complete.
  ! In coarse_only mode, only pass-1 chunks are required.
  ! In two-tier mode, pass-1 must be complete and pass-2 must be complete
  ! when present; if no pass-2 chunks were generated, pass-1 completion is
  ! treated as terminal.
  logical function get_finished( self )
    class(ptcl_sieve), intent(in) :: self
    ! pass-1: must have at least one chunk and all must be complete
    get_finished = self%get_n_chunks_coarse() > 0
    if( .not. get_finished ) return
    get_finished = all(self%chunks_coarse(:)%complete .or. self%chunks_coarse(:)%failed)
    if( .not. get_finished ) return
    ! pass-1-only runs terminate once all pass-1 chunks are done
    if( self%coarse_only ) return
    ! pass-2: when absent, pass-1 completion is terminal
    if( self%get_n_chunks_fine() == 0 )then
      get_finished = .true.
      return
    endif
    get_finished = all(self%chunks_fine(:)%complete .or. self%chunks_fine(:)%failed)
  end function get_finished

  ! --------------------------------------------------------------------------
  ! APPEND HELPERS
  ! --------------------------------------------------------------------------

  ! Appends a single chunk2D to the end of the pass-1 microchunk array.
  subroutine append_chunk_coarse( self, new_chunk )
    class(ptcl_sieve),   intent(inout) :: self
    type(chunk2D_state), intent(in)    :: new_chunk
    type(chunk2D_state), allocatable   :: tmp(:)
    integer :: n
    n = self%get_n_chunks_coarse()
    allocate(tmp(n + 1))
    if( n > 0 ) tmp(1:n) = self%chunks_coarse
    tmp(n + 1) = new_chunk
    call move_alloc(tmp, self%chunks_coarse)
  end subroutine append_chunk_coarse

  ! Appends a single chunk2D to the end of the pass-2 microchunk array.
  subroutine append_chunk_fine( self, new_chunk )
    class(ptcl_sieve),   intent(inout) :: self
    type(chunk2D_state), intent(in)    :: new_chunk
    type(chunk2D_state), allocatable   :: tmp(:)
    integer :: n
    n = self%get_n_chunks_fine()
    allocate(tmp(n + 1))
    if( n > 0 ) tmp(1:n) = self%chunks_fine
    tmp(n + 1) = new_chunk
    call move_alloc(tmp, self%chunks_fine)
  end subroutine append_chunk_fine


  ! --------------------------------------------------------------------------
  ! GENERATION
  ! --------------------------------------------------------------------------

  ! Partitions un-included records from project_list into pass-1 microchunks,
  ! each holding up to coarse_defaults%pop_threshold particles. For each chunk:
  ! creates its subdirectory, accumulates micrographs until the threshold is
  ! reached, builds and writes a project file, applies any
  ! SIMPLE_CHUNK_PARTITION environment override, configures the command line,
  ! marks the consumed records as included in project_list, and updates the
  ! imported_projects.txt file table with all currently included project files.
  subroutine generate_chunks_coarse( self, project_list )
    class(ptcl_sieve), intent(inout) :: self
    type(rec_list),        intent(inout) :: project_list

    type(string),     allocatable   :: projfiles(:), projfiles_all(:)
    type(string),     allocatable   :: unique_projfiles(:)
    logical,          allocatable   :: included(:)
    integer,          allocatable   :: ids(:), nptcls(:)
    type(chunk2D_state)             :: new_chunk
    type(sp_project)                :: chunk_project
    type(rec_list)                  :: chunk_project_list
    type(string)                    :: chunk_folder
    character(len=STDLEN)           :: chunk_part_env
    integer(timer_int_kind)         :: t0
    integer                         :: i, j, imic, chunk_nptcls, envlen, chunk_id

    t0 = timer_start()

    ! Pre-chunked mode: import existing per-chunk project files from
    ! project_list instead of creating new pass-1 chunks from particle slices.
    if( self%pre_chunked ) then
      if( self%get_n_chunks_coarse() > 0 ) then
        call timer_stop(t0, string('generate_chunks_coarse'))
        return
      end if

      included      = project_list%get_included_flags()
      ids           = pack(project_list%get_ids(),       .not. included)
      projfiles_all = project_list%get_projnames()
      projfiles     = pack(projfiles_all, .not. included)
      if( size(ids) == 0 ) then
        call timer_stop(t0, string('generate_chunks_coarse'))
        return
      end if
      
      unique_projfiles = remove_duplicates(projfiles)
      do i = 1, size(unique_projfiles)
        if( .not. file_exists(unique_projfiles(i)) ) then
          THROW_HARD('missing pre-chunked project file: '//unique_projfiles(i)%to_char())
        end if

        chunk_id                  = self%get_n_chunks_coarse() + 1
        chunk_folder              = string(self%outdir_chunks_coarse%to_char() // &
                                           '/chunk_coarse_' // int2str(chunk_id))
        call simple_mkdir(chunk_folder)
        new_chunk%id              = chunk_id
        new_chunk%folder          = simple_abspath(chunk_folder)
        new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                           '/chunk_coarse_' // int2str(chunk_id) // METADATA_EXT)
        call simple_copy_file(unique_projfiles(i), new_chunk%projfile)
        call chunk_project%read(new_chunk%projfile)
        new_chunk%nptcls          = chunk_project%os_ptcl2D%get_noris()
        new_chunk%nptcls_selected = chunk_project%os_ptcl2D%count_state_gt_zero()
        chunk_project%os_ptcl3D   = chunk_project%os_ptcl2D
        call chunk_project%write()
        call chunk_project%kill()
        call self%generate_chunk_coarse_cline(new_chunk, new_chunk%nptcls_selected)
        call self%append_chunk_coarse(new_chunk)

        do j = 1, size(ids)
          if( projfiles(j) == unique_projfiles(i) ) call project_list%set_included_flags([ids(j), ids(j)])
        end do

        self%last_import = c_time(0_c_long)
      end do

      included = project_list%get_included_flags()
      projfiles = project_list%get_projnames()
      projfiles = pack(projfiles, included)
      projfiles = remove_duplicates(projfiles)
      call write_filetable(string('imported_projects.txt'), projfiles)

      call timer_stop(t0, string('generate_chunks_coarse'))
      return
    end if

    ! NOTE: pass-1 uses a soft cap. The trigger micrograph that crosses the
    ! threshold is included in the current chunk to preserve contiguous slices.
    do while( project_list%get_nptcls_tot(l_not_included=.true.) > self%coarse_defaults%pop_threshold )
      included = project_list%get_included_flags()
      ids      = pack(project_list%get_ids(),    .not. included)
      nptcls   = pack(project_list%get_nptcls(), .not. included)
      if( size(ids) == 0 ) exit

      ! Accumulate micrographs until the soft cap is exceeded.
      chunk_nptcls = 0
      do imic = 1, size(nptcls)
        chunk_nptcls = chunk_nptcls + nptcls(imic)
        if( chunk_nptcls > self%coarse_defaults%pop_threshold ) exit
      end do
      if( chunk_nptcls == 0 ) exit

      ! Create the chunk folder and populate chunk fields
      chunk_id                  = self%get_n_chunks_coarse() + 1
      chunk_folder              = string(self%outdir_chunks_coarse%to_char() // &
                                         '/chunk_coarse_' // int2str(chunk_id))
      call simple_mkdir(chunk_folder)
      new_chunk%id              = chunk_id
      new_chunk%nptcls          = chunk_nptcls
      new_chunk%nptcls_selected = chunk_nptcls
      new_chunk%folder          = simple_abspath(chunk_folder)
      new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                         '/chunk_coarse_' // int2str(chunk_id) // METADATA_EXT)
      call self%generate_chunk_coarse_cline(new_chunk, new_chunk%nptcls_selected)

      ! Build the project file from the sliced record list
      call project_list%slice(ids(1), ids(imic), chunk_project_list)
      call chunk_project%projrecords2proj(chunk_project_list)
      call chunk_project_list%kill()
      call chunk_project%update_projinfo(new_chunk%cline)
      call chunk_project%update_compenv(new_chunk%cline)

      ! Apply queue partition override if the environment variable is set
      call get_environment_variable('SIMPLE_CHUNK_PARTITION', chunk_part_env, envlen)
      if( envlen > 0 ) call chunk_project%compenv%set(1, 'qsys_partition', trim(chunk_part_env))

      call chunk_project%write(new_chunk%projfile)
      call chunk_project%kill()

      call self%append_chunk_coarse(new_chunk)
      call project_list%set_included_flags([ids(1), ids(imic)])

      ! Update the imported-projects file table with all included project files
      included = project_list%get_included_flags()
      projfiles = project_list%get_projnames()
      projfiles = pack(projfiles, included)
      projfiles = remove_duplicates(projfiles)
      call write_filetable(string('imported_projects.txt'), projfiles)

      self%last_import = c_time(0_c_long)

      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> COARSE CHUNK # ', chunk_id, ' GENERATED WITH ', chunk_nptcls, ' PARTICLES'
    end do

    ! Final ingestion: sweep any remaining sub-threshold un-included particles
    ! into a pass-1 chunk flagged as rejection-complete so the pass-2 final
    ! flush consumes them without running 2D classification.
    if( self%final_ingestion ) then
      included = project_list%get_included_flags()
      ids      = pack(project_list%get_ids(),    .not. included)
      nptcls   = pack(project_list%get_nptcls(), .not. included)
      if( size(ids) > 0 ) then
        chunk_id                  = self%get_n_chunks_coarse() + 1
        chunk_folder              = string(self%outdir_chunks_coarse%to_char() // &
                                           '/chunk_coarse_' // int2str(chunk_id))
        call simple_mkdir(chunk_folder)
        new_chunk%id              = chunk_id
        new_chunk%nptcls          = sum(nptcls)
        new_chunk%nptcls_selected = sum(nptcls)
        new_chunk%folder          = simple_abspath(chunk_folder)
        new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                           '/chunk_coarse_' // int2str(chunk_id) // METADATA_EXT)
        call self%generate_chunk_coarse_cline(new_chunk, new_chunk%nptcls_selected)
        call project_list%slice(ids(1), ids(size(ids)), chunk_project_list)
        call chunk_project%projrecords2proj(chunk_project_list)
        call chunk_project_list%kill()
        call chunk_project%update_projinfo(new_chunk%cline)
        call chunk_project%update_compenv(new_chunk%cline)
        call chunk_project%write(new_chunk%projfile)
        call chunk_project%kill()
        new_chunk%abinitio2D_complete = .true.
        new_chunk%rejection_complete  = .true.
        call self%append_chunk_coarse(new_chunk)
        call project_list%set_included_flags([ids(1), ids(size(ids))])
        included  = project_list%get_included_flags()
        projfiles = project_list%get_projnames()
        projfiles = pack(projfiles, included)
        projfiles = remove_duplicates(projfiles)
        call write_filetable(string('imported_projects.txt'), projfiles)
        write(logfhandle,'(A,I8,A)') &
          '>>> FINAL INGESTION: STAGED ', new_chunk%nptcls, ' REMAINING PARTICLES FOR PASS 2'
      end if
    end if
    call timer_stop(t0, string('generate_chunks_coarse'))
  end subroutine generate_chunks_coarse

  ! Merges rejection-complete pass-1 microchunks into pass-2 microchunks,
  ! each accumulating up to fine_defaults%pop_threshold selected particles.
  subroutine generate_chunks_fine( self )
    class(ptcl_sieve), intent(inout) :: self
    type(string),      allocatable   :: projfiles(:)
    integer,           allocatable   :: consumed(:)
    type(chunk2D_state)              :: new_chunk
    type(sp_project)                 :: chunk_project
    type(string)                     :: chunk_folder
    integer(timer_int_kind)          :: t0
    integer                          :: i, chunk_nptcls, chunk_nptcls_selected, chunk_id, n_consumed

    t0 = timer_start()
    do while( self%get_n_pass_1_non_rejected_ptcls() > self%fine_defaults%pop_threshold )
      allocate(projfiles(self%get_n_chunks_coarse()), consumed(self%get_n_chunks_coarse()))
      chunk_nptcls          = 0
      chunk_nptcls_selected = 0
      n_consumed            = 0

      do i = 1, self%get_n_chunks_coarse()
        associate( src => self%chunks_coarse(i) )
          if( src%rejection_complete .and. .not. src%complete .and. .not. src%failed ) then
            n_consumed            = n_consumed            + 1
            chunk_nptcls_selected = chunk_nptcls_selected + src%nptcls_selected
            chunk_nptcls          = chunk_nptcls          + src%nptcls
            projfiles(n_consumed) = src%projfile
            consumed(n_consumed)  = i
          end if
          if( chunk_nptcls_selected > self%fine_defaults%pop_threshold ) exit
        end associate
      end do
      if( n_consumed == 0 ) then
        deallocate(projfiles, consumed)
        exit
      end if
      projfiles = projfiles(:n_consumed)
      consumed  = consumed(:n_consumed)

      chunk_id                  = self%get_n_chunks_fine() + 1
      chunk_folder              = string(self%outdir_chunks_fine%to_char() // &
                                         '/chunk_fine_' // int2str(chunk_id))
      call simple_mkdir(chunk_folder)
      new_chunk%id              = chunk_id
      new_chunk%nptcls          = chunk_nptcls
      new_chunk%nptcls_selected = chunk_nptcls_selected
      new_chunk%folder          = simple_abspath(chunk_folder)
      new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                         '/chunk_fine_' // int2str(chunk_id) // METADATA_EXT)
      call self%generate_chunk_fine_cline(new_chunk, new_chunk%nptcls_selected)

      call merge_and_clear(projfiles, chunk_folder, chunk_project, new_chunk%cline)

      do i = 1, size(consumed)
        self%chunks_coarse(consumed(i))%complete = .true.
        call simple_touch(self%chunks_coarse(consumed(i))%folder%to_char() // '/COMPLETE')
      end do

      ! flag this as the final sieve chunk once final ingestion is signaled and all pass-1 chunks are done
      if( self%final_ingestion .and. all(self%chunks_coarse(:)%complete .or. self%chunks_coarse(:)%failed) )then
        if( chunk_project%os_out%get_noris() < 1 ) call chunk_project%os_out%new(1, is_ptcl=.false.)
        call chunk_project%os_out%set(1, 'sieve_final', 'yes')
        write(logfhandle,'(A,I6,A,I8,A,I8,A)') '>>> FINAL FINE CHUNK # ', chunk_id, &
        ' GENERATED WITH ', chunk_nptcls_selected, '/', chunk_nptcls, ' PARTICLES'
      else
        write(logfhandle,'(A,I6,A,I8,A,I8,A)') '>>> FINE CHUNK # ', chunk_id, &
        ' GENERATED WITH ', chunk_nptcls_selected, '/', chunk_nptcls, ' PARTICLES'
      end if

      call chunk_project%write(new_chunk%projfile)
      call chunk_project%kill()

      call self%append_chunk_fine(new_chunk)
      deallocate(projfiles, consumed)
      
    end do

    ! Flush a final small pass-2 chunk when final ingestion is signaled.
    if( self%final_ingestion .and. &
        self%get_n_pass_1_non_rejected_ptcls() > 0 .and. all(self%chunks_coarse(:)%rejection_complete .or. self%chunks_coarse(:)%failed) ) then
      allocate(projfiles(self%get_n_chunks_coarse()), consumed(self%get_n_chunks_coarse()))
      chunk_nptcls          = 0
      chunk_nptcls_selected = 0
      n_consumed            = 0
      do i = 1, self%get_n_chunks_coarse()
        associate( src => self%chunks_coarse(i) )
          if( src%rejection_complete .and. .not. src%complete .and. .not. src%failed ) then
            n_consumed            = n_consumed + 1
            chunk_nptcls_selected = chunk_nptcls_selected + src%nptcls_selected
            chunk_nptcls          = chunk_nptcls + src%nptcls
            projfiles(n_consumed) = src%projfile
            consumed(n_consumed)  = i
          end if
        end associate
      end do
      if( n_consumed > 0 ) then
        projfiles = projfiles(:n_consumed)
        consumed  = consumed(:n_consumed)
        chunk_id                  = self%get_n_chunks_fine() + 1
        chunk_folder              = string(self%outdir_chunks_fine%to_char() // &
                                           '/chunk_fine_' // int2str(chunk_id))
        call simple_mkdir(chunk_folder)
        new_chunk%id              = chunk_id
        new_chunk%nptcls          = chunk_nptcls
        new_chunk%nptcls_selected = chunk_nptcls_selected
        new_chunk%folder          = simple_abspath(chunk_folder)
        new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                           '/chunk_fine_' // int2str(chunk_id) // METADATA_EXT)
        call self%generate_chunk_fine_cline(new_chunk, new_chunk%nptcls_selected)
        call merge_and_clear(projfiles, chunk_folder, chunk_project, new_chunk%cline)
        ! flag this as the final sieve chunk once final ingestion is signaled and all pass-1 chunks are done
        if( chunk_project%os_out%get_noris() < 1 ) call chunk_project%os_out%new(1, is_ptcl=.false.)
        call chunk_project%os_out%set(1, 'sieve_final', 'yes')
        call chunk_project%write(new_chunk%projfile)
        call chunk_project%kill()
        do i = 1, size(consumed)
          self%chunks_coarse(consumed(i))%complete = .true.
          call simple_touch(self%chunks_coarse(consumed(i))%folder%to_char() // '/COMPLETE')
        end do
        call self%append_chunk_fine(new_chunk)
        write(logfhandle,'(A,I6,A,I8,A,I8,A)') '>>> FINAL FINE CHUNK # ', chunk_id, &
          ' GENERATED WITH ', chunk_nptcls_selected, '/', chunk_nptcls, ' PARTICLES'
      end if
      deallocate(projfiles, consumed)
    end if

    call timer_stop(t0, string('generate_chunks_fine'))
  end subroutine generate_chunks_fine

  ! --------------------------------------------------------------------------
  ! COMMAND-LINE BUILDERS
  ! --------------------------------------------------------------------------

  ! Populates the abinitio2D command line for a pass-1 microchunk with:
  ! program name, project file and name, no-mkdir flag, thread count, mask
  ! diameter, class count, low-pass stop cutoff (DEFAULT_MICRO_P1_LP), and
  ! wall-time limit.
  subroutine generate_chunk_coarse_cline( self, new_chunk, nptcls )
    class(ptcl_sieve),   intent(inout) :: self
    type(chunk2D_state), intent(inout) :: new_chunk
    integer,             intent(in)    :: nptcls
    type(string)                       :: server_address
    server_address = self%qenv%get_persistent_worker_server_address()
    associate( cline => new_chunk%cline )
      call cline%set('prg',                               'abinitio2D')
      call cline%set('projfile',                    new_chunk%projfile)
      call cline%set('projname',                        'chunk_coarse')
      call cline%set('mkdir',                                     'no')
      call cline%set('nthr',                                 self%nthr)
      call cline%set('mskdiam',                           self%mskdiam)
      call cline%set('box_crop',         self%coarse_defaults%box_crop)
      call cline%set('nsample', min(self%coarse_defaults%nsample, nptcls))
      call cline%set('ncls',             self%coarse_defaults%ncls)
      call cline%set('lpstart',       self%coarse_defaults%lpstart)
      call cline%set('lpstop',         self%coarse_defaults%lpstop)
      call cline%set('walltime',                      DEFAULT_WALLTIME)
      if( self%nparts > 1             ) call cline%set('nparts',           self%nparts)
      if( server_address%strlen() > 0 ) call cline%set('worker_server', server_address)
    end associate
    call server_address%kill()
  end subroutine generate_chunk_coarse_cline

  ! Populates the abinitio2D command line for a pass-2 microchunk.
  subroutine generate_chunk_fine_cline( self, new_chunk, nptcls )
    class(ptcl_sieve),   intent(inout) :: self
    type(chunk2D_state), intent(inout) :: new_chunk
    integer,             intent(in)    :: nptcls
    type(string)                       :: server_address
    server_address = self%qenv%get_persistent_worker_server_address()
    associate( cline => new_chunk%cline )
      call cline%set('prg',                               'abinitio2D')
      call cline%set('projfile',                    new_chunk%projfile)
      call cline%set('projname',                          'chunk_fine')
      call cline%set('mkdir',                                     'no')
      call cline%set('nthr',                                 self%nthr)
      call cline%set('mskdiam',                           self%mskdiam)
      call cline%set('box_crop',            self%fine_defaults%box_crop)
      call cline%set('nsample',   min(self%fine_defaults%nsample, nptcls))
      call cline%set('ncls',                  self%fine_defaults%ncls)
      call cline%set('lpstart',            self%fine_defaults%lpstart)
      call cline%set('lpstop',              self%fine_defaults%lpstop)
      call cline%set('walltime',                      DEFAULT_WALLTIME)
      if( self%nparts > 1             ) call cline%set('nparts',           self%nparts)
      if( server_address%strlen() > 0 ) call cline%set('worker_server', server_address)
    end associate
    call server_address%kill()
  end subroutine generate_chunk_fine_cline

  ! --------------------------------------------------------------------------
  ! COMBINATION
  ! --------------------------------------------------------------------------

  ! Merges completed microchunk project files into a single combined project
  ! written to completedir. In coarse_only mode, combines rejection-complete
  ! pass-1 chunks; otherwise combines complete pass-2 chunks.
  ! No-op if no eligible chunks exist or if the combined file already exists.
  subroutine combine_completed_chunks( self, combined_projfile )
    class(ptcl_sieve), intent(inout) :: self
    type(string),          intent(in)    :: combined_projfile
    type(string),     allocatable :: projfiles(:)
    type(sp_project)              :: combined_project
    integer(timer_int_kind)       :: t0
    integer                       :: i, nptcls_tot, nptcls_sel_tot, n_complete, n_chunks
    logical                       :: use_pass_1

    use_pass_1 = self%coarse_only
    if( use_pass_1 ) then
      n_chunks = self%get_n_chunks_coarse()
    else
      n_chunks = self%get_n_chunks_fine()
    end if

    if( n_chunks == 0 ) return
    if( file_exists(combined_projfile) )       return

    t0 = timer_start()

    ! Collect projfiles from eligible completed chunks.
    allocate(projfiles(n_chunks))
    nptcls_tot     = 0
    nptcls_sel_tot = 0
    n_complete     = 0
    if( use_pass_1 ) then
      do i = 1, self%get_n_chunks_coarse()
        associate( chunk => self%chunks_coarse(i) )
          if( chunk%rejection_complete .and. .not. chunk%failed ) then
            n_complete            = n_complete + 1
            projfiles(n_complete) = chunk%projfile
            nptcls_tot            = nptcls_tot     + chunk%nptcls
            nptcls_sel_tot        = nptcls_sel_tot + chunk%nptcls_selected
          end if
        end associate
      end do
    else
      do i = 1, self%get_n_chunks_fine()
        associate( chunk => self%chunks_fine(i) )
          if( chunk%complete .and. .not. chunk%failed ) then
            n_complete            = n_complete + 1
            projfiles(n_complete) = chunk%projfile
            nptcls_tot            = nptcls_tot     + chunk%nptcls
            nptcls_sel_tot        = nptcls_sel_tot + chunk%nptcls_selected
          end if
        end associate
      end do
    end if
    if( n_complete == 0 ) then
      deallocate(projfiles)
      return
    end if
    projfiles = projfiles(:n_complete)

    call merge_chunk_projfiles(projfiles, simple_abspath(self%completedir), combined_project, write_proj=.false.)
    call combined_project%write(combined_projfile)
    call combined_project%kill()
    deallocate(projfiles)

    if( use_pass_1 ) then
      write(logfhandle,'(A,I6,A,I8,A,I8,A)') '>>> COMBINED ', n_complete, &
        ' COARSE CHUNKS INTO ', nptcls_sel_tot, '/', nptcls_tot, ' PARTICLES'
    else
      write(logfhandle,'(A,I6,A,I8,A,I8,A)') '>>> COMBINED ', n_complete, &
        ' FINE CHUNKS INTO ', nptcls_sel_tot, '/', nptcls_tot, ' PARTICLES'
    end if
    call timer_stop(t0, string('combine_completed_chunks'))
  end subroutine combine_completed_chunks

  ! --------------------------------------------------------------------------
  ! SUBMISSION AND COLLECTION
  ! --------------------------------------------------------------------------

  ! Submits pending chunks to the queue asynchronously, up to the nparallel
  ! concurrency limit. Priority order: pass-2 > pass-1.
  ! Skips chunks that are already running or complete. For each submitted
  ! chunk: changes into its working directory, dispatches the job, marks it
  ! as running, and releases its command line. Restores the original working
  ! directory on completion.
  subroutine submit( self )
    class(ptcl_sieve), intent(inout) :: self
    type(string)            :: cwd
    integer(timer_int_kind) :: t0
    integer                 :: i
    t0 = timer_start()
    call simple_getcwd(cwd)

    ! Submit pass-2 microchunks first
    do i = 1, self%get_n_chunks_fine()
      if( self%get_n_chunks_running() >= self%nparallel ) exit
      associate( chunk => self%chunks_fine(i) )
        if( chunk%failed ) cycle
        if( chunk%abinitio2D_running .or. chunk%abinitio2D_complete ) cycle
        call simple_chdir(chunk%folder)
        CWD_GLOB = chunk%folder%to_char()
        call self%qenv%exec_simple_prg_in_queue_async( &
          chunk%cline, string('./distr_ptcl_sieve'), string('simple_log_chunk_fine'))
        chunk%abinitio2D_running = .true.
        call chunk%cline%kill()
        write(logfhandle,'(A,I6)') '>>> INITIATED 2D ANALYSIS OF FINE CHUNK # ', chunk%id
      end associate
    end do

    ! Submit pass-1 microchunks last — lowest priority
    do i = 1, self%get_n_chunks_coarse()
      if( self%get_n_chunks_running() >= self%nparallel ) exit
      associate( chunk => self%chunks_coarse(i) )
        if( chunk%failed ) cycle
        if( chunk%abinitio2D_running .or. chunk%abinitio2D_complete ) cycle
        call simple_chdir(chunk%folder)
        CWD_GLOB = chunk%folder%to_char()
        call self%qenv%exec_simple_prg_in_queue_async( &
          chunk%cline, string('./distr_ptcl_sieve'), string('simple_log_chunk_coarse'))
        chunk%abinitio2D_running = .true.
        call chunk%cline%kill()
        write(logfhandle,'(A,I6)') '>>> INITIATED 2D ANALYSIS OF COARSE CHUNK # ', chunk%id
      end associate
    end do

    call simple_chdir(cwd)
    CWD_GLOB = cwd%to_char()
    call timer_stop(t0, string('submit'))
  end subroutine submit

  ! Polls all running pass-1 and pass-2 chunks for the ABINITIO2D_FINISHED
  ! sentinel file. For each newly completed chunk: transitions it from running
  ! to complete and immediately runs class-average rejection.
  subroutine collect_and_reject( self )
    class(ptcl_sieve), intent(inout) :: self
    logical,                      allocatable :: cls_msk(:)
    type(sp_project)        :: spproj
    integer(timer_int_kind) :: t0
    integer                 :: i, non_zero, ncls
    real                    :: smpd_dummy
    t0 = timer_start()

    do i = 1, self%get_n_chunks_coarse()
      associate( chunk => self%chunks_coarse(i) )
        if( chunk%failed ) cycle
        if( chunk%abinitio2D_running ) then
          if( file_exists(chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED) ) then
            chunk%abinitio2D_running  = .false.
            chunk%abinitio2D_complete = .true.
            write(logfhandle,'(A,I6)') '>>> COMPLETED 2D ANALYSIS OF COARSE CHUNK # ', chunk%id
          end if
        end if
        call self%reject_cavgs(chunk, string(LABEL_COARSE))
        if( self%coarse_only .and. chunk%rejection_complete .and. .not. chunk%complete ) then
          call spproj%read_segment('mic',    chunk%projfile)
          call spproj%read_segment('ptcl2D', chunk%projfile)
          call spproj%read_segment('cls2D',  chunk%projfile)
          call spproj%read_segment('out',    chunk%projfile)
          non_zero = spproj%os_ptcl2D%count_state_gt_zero()
          self%n_accepted_ptcls = self%n_accepted_ptcls + non_zero
          self%n_accepted_mics  = self%n_accepted_mics + spproj%os_mic%count_state_gt_zero()
          self%n_rejected_ptcls = self%n_rejected_ptcls + spproj%os_ptcl2D%get_noris() - non_zero
          call spproj%get_cavgs_stk(self%latest_stkname, ncls, smpd_dummy)
          self%latest_jpeg = swap_suffix(self%latest_stkname, JPG_EXT, MRC_EXT)
          call spproj%cavgs2jpg(self%latest_jpeg_inds, self%latest_jpeg, self%latest_jpeg_xtiles, self%latest_jpeg_ytiles, ignore_states=.true.)
          self%latest_jpeg_pops       = spproj%os_cls2D%get_all_asint('pop')
          self%latest_jpeg_res        = spproj%os_cls2D%get_all('res')
          self%latest_jpeg_selection  = spproj%os_cls2D%get_all_asint('state')
          allocate(cls_msk, source=self%latest_jpeg_inds /= 0)
          self%latest_jpeg_inds       = pack(self%latest_jpeg_inds, cls_msk)
          self%latest_jpeg_pops       = pack(self%latest_jpeg_pops, cls_msk)
          self%latest_jpeg_res        = pack(self%latest_jpeg_res,  cls_msk)
          self%latest_jpeg_selection  = pack(self%latest_jpeg_selection, cls_msk)
          deallocate(cls_msk)
          call spproj%kill()
          call simple_copy_file(chunk%projfile, self%completedir // '/' // basename(chunk%projfile))
          call simple_touch(chunk%folder%to_char() // '/COMPLETE')
          chunk%complete            = .true.
          write(logfhandle,'(A,I6)') '>>> FINALISED FINE CHUNK # ', chunk%id
        end if
        if( (.not. self%coarse_only) .and. chunk%rejection_complete .and. .not. chunk%complete ) then
          if( self%get_n_chunks_fine() == 0 ) then
              call spproj%read_segment('cls2D',  chunk%projfile)
              call spproj%read_segment('out',    chunk%projfile)
              call spproj%get_cavgs_stk(self%latest_stkname, ncls, smpd_dummy)
              self%latest_jpeg = swap_suffix(self%latest_stkname, JPG_EXT, MRC_EXT)
              call spproj%cavgs2jpg(self%latest_jpeg_inds, self%latest_jpeg, self%latest_jpeg_xtiles, self%latest_jpeg_ytiles, ignore_states=.true.)
              self%latest_jpeg_pops       = spproj%os_cls2D%get_all_asint('pop')
              self%latest_jpeg_res        = spproj%os_cls2D%get_all('res')
              self%latest_jpeg_selection  = spproj%os_cls2D%get_all_asint('state')
              allocate(cls_msk, source=self%latest_jpeg_inds /= 0)
              self%latest_jpeg_inds       = pack(self%latest_jpeg_inds, cls_msk)
              self%latest_jpeg_pops       = pack(self%latest_jpeg_pops, cls_msk)
              self%latest_jpeg_res        = pack(self%latest_jpeg_res,  cls_msk)
              self%latest_jpeg_selection  = pack(self%latest_jpeg_selection, cls_msk)
              deallocate(cls_msk)
              call spproj%kill()
          end if
          if( chunk%nptcls_selected == 0 ) then
            call simple_touch(chunk%folder%to_char() // '/COMPLETE')
            chunk%complete = .true.
            write(logfhandle,'(A,I6,A)') '>>> FINALISED EMPTY COARSE CHUNK # ', chunk%id, ' (NO SELECTED PARTICLES)'
          end if
        end if
      end associate
    end do

    do i = 1, self%get_n_chunks_fine()
      associate( chunk => self%chunks_fine(i) )
        if( chunk%failed ) cycle
        if( chunk%abinitio2D_running ) then
          if( file_exists(chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED) ) then
            chunk%abinitio2D_running  = .false.
            chunk%abinitio2D_complete = .true.
            write(logfhandle,'(A,I6)') '>>> COMPLETED 2D ANALYSIS OF FINE CHUNK # ', chunk%id
          end if
        end if
        call self%reject_cavgs(chunk, string(LABEL_FINE))
        if( chunk%rejection_complete .and. .not. chunk%complete ) then
          call spproj%read_segment('mic',    chunk%projfile)
          call spproj%read_segment('ptcl2D', chunk%projfile)
          call spproj%read_segment('cls2D',  chunk%projfile)
          call spproj%read_segment('out',    chunk%projfile)
          non_zero = spproj%os_ptcl2D%count_state_gt_zero()
          self%n_accepted_ptcls = self%n_accepted_ptcls + non_zero
          self%n_accepted_mics  = self%n_accepted_mics + spproj%os_mic%count_state_gt_zero()
          self%n_rejected_ptcls = self%n_rejected_ptcls + spproj%os_ptcl2D%get_noris() - non_zero
          call spproj%get_cavgs_stk(self%latest_stkname, ncls, smpd_dummy)
          self%latest_jpeg = swap_suffix(self%latest_stkname, JPG_EXT, MRC_EXT)
          call spproj%cavgs2jpg(self%latest_jpeg_inds, self%latest_jpeg, self%latest_jpeg_xtiles, self%latest_jpeg_ytiles, ignore_states=.true.)
          self%latest_jpeg_pops       = spproj%os_cls2D%get_all_asint('pop')
          self%latest_jpeg_res        = spproj%os_cls2D%get_all('res')
          self%latest_jpeg_selection  = spproj%os_cls2D%get_all_asint('state')
          allocate(cls_msk, source=self%latest_jpeg_inds /= 0)
          self%latest_jpeg_inds      = pack(self%latest_jpeg_inds, cls_msk)
          self%latest_jpeg_pops      = pack(self%latest_jpeg_pops, cls_msk)
          self%latest_jpeg_res       = pack(self%latest_jpeg_res,  cls_msk)
          self%latest_jpeg_selection = pack(self%latest_jpeg_selection, cls_msk)
          deallocate(cls_msk)
          call spproj%kill()
          call simple_copy_file(chunk%projfile, self%completedir // '/' // basename(chunk%projfile))
          call simple_touch(chunk%folder%to_char() // '/COMPLETE')
          chunk%complete            = .true.
          write(logfhandle,'(A,I6)') '>>> FINALISED FINE CHUNK # ', chunk%id
        end if
      end associate
    end do

    call timer_stop(t0, string('collect_and_reject'))
  end subroutine collect_and_reject

  ! Alternate rejection path mirroring run_cavg_quality_selection logic:
  ! class-average quality model scoring, state mapping via map_cavgs_selection,
  ! and selected/rejected stack export.
  subroutine reject_cavgs( self, chunk, label )
    class(ptcl_sieve),   intent(inout) :: self
    type(chunk2D_state), intent(inout) :: chunk
    type(string),        intent(in)    :: label
    type(image),         allocatable   :: cavg_imgs(:), out_imgs(:)
    integer,             allocatable   :: states(:)
    type(cavg_quality_result)          :: quality
    type(support_model_metrics)        :: compat_metrics
    type(sp_project)                   :: spproj
    type(string)                       :: stkname, jpgname
    integer(timer_int_kind)            :: t0
    integer                            :: ncls, iimg, nout
    real                               :: smpd_dummy

    if( .not. chunk%abinitio2D_complete ) return
    if( chunk%failed )                    return
    if( chunk%rejection_complete )        return

    t0 = timer_start()

    call spproj%read(chunk%projfile)
    cavg_imgs = read_cavgs_into_imgarr(spproj)
    if( label == LABEL_COARSE ) then
      call evaluate_cavg_quality_hard_reject(cavg_imgs, spproj%os_cls2D, 0.0, quality, CAVG_QUALITY_CONTEXT_SIEVE)
      do iimg = 1, size(cavg_imgs)
          if( quality%states(iimg) == 0 ) then
              call spproj%os_cls2D%set(iimg, 'state', 0)
              call spproj%os_cls2D%set(iimg, 'rejection_reason', string('coarse_reject: ')//cavg_rejection_reason_string(quality%reasons(iimg)))
          end if
      end do
      if( .not. self%coarse_compatibility_model%converged() ) then
        call self%coarse_compatibility_model%train(spproj)
      end if
      call self%coarse_compatibility_model%infer(spproj)
      call self%coarse_compatibility_model%get_support_model_metrics(compat_metrics)
      write(logfhandle,'(A,I6,A,3(F10.4,1X),A,3(F10.4,1X),A,L1,A,L1,A,L1)') &
        '>>> COARSE COMPAT METRICS CHUNK # ', chunk%id, ' a/b/c=', &
        compat_metrics%axis_a, compat_metrics%axis_b, compat_metrics%axis_c, ' da/db/dc=', &
        compat_metrics%delta_a, compat_metrics%delta_b, compat_metrics%delta_c, ' valid=', &
        compat_metrics%valid, ' delta_valid=', compat_metrics%delta_valid, ' converged=', compat_metrics%converged

    else
      call evaluate_cavg_quality_hard_reject(cavg_imgs, spproj%os_cls2D, 0.0, quality, CAVG_QUALITY_CONTEXT_SIEVE)
      do iimg = 1, size(cavg_imgs)
          if( quality%states(iimg) == 0 ) then
              call spproj%os_cls2D%set(iimg, 'state', 0)
              call spproj%os_cls2D%set(iimg, 'rejection_reason', string('coarse_reject: ')//cavg_rejection_reason_string(quality%reasons(iimg)))
          end if
      end do
      if( .not. self%fine_compatibility_model%converged() ) then
        call self%fine_compatibility_model%train(spproj)
      end if
      call self%fine_compatibility_model%infer(spproj)
      call self%fine_compatibility_model%get_support_model_metrics(compat_metrics)
      write(logfhandle,'(A,I6,A,3(F10.4,1X),A,3(F10.4,1X),A,L1,A,L1,A,L1)') &
        '>>> FINE COMPAT METRICS CHUNK # ', chunk%id, ' a/b/c=', &
        compat_metrics%axis_a, compat_metrics%axis_b, compat_metrics%axis_c, ' da/db/dc=', &
        compat_metrics%delta_a, compat_metrics%delta_b, compat_metrics%delta_c, ' valid=', &
        compat_metrics%valid, ' delta_valid=', compat_metrics%delta_valid, ' converged=', compat_metrics%converged
    end if
    states = spproj%os_cls2D%get_all_asint('state')
    call spproj%map_cavgs_selection(states)
    call spproj%write()

    ! write selected
    allocate(out_imgs(count(states > 0)))
    nout    = 0
    stkname = swap_suffix(chunk%projfile, '_selected'//MRC_EXT, METADATA_EXT)
    jpgname = swap_suffix(stkname, JPG_EXT, MRC_EXT)
    do iimg = 1, size(states)
      if( states(iimg) > 0 ) then
        nout = nout + 1
        call out_imgs(nout)%copy(cavg_imgs(iimg))
      end if
    end do
    call write_imgarr(out_imgs, stkname)
    call dealloc_imgarr(out_imgs)
    call mrc2jpeg_tiled(stkname, jpgname)
    write(logfhandle,'(A,A)') '>>> JPEG ', jpgname%to_char()

    ! write rejected
    allocate(out_imgs(count(states == 0)))
    nout    = 0
    stkname = swap_suffix(chunk%projfile, '_rejected'//MRC_EXT, METADATA_EXT)
    jpgname = swap_suffix(stkname, JPG_EXT, MRC_EXT)
    do iimg = 1, size(states)
      if( states(iimg) == 0 ) then
        nout = nout + 1
        call out_imgs(nout)%copy(cavg_imgs(iimg))
      end if
    end do
    call write_imgarr(out_imgs, stkname)
    call dealloc_imgarr(out_imgs)
    call mrc2jpeg_tiled(stkname, jpgname)
    write(logfhandle,'(A,A)') '>>> JPEG ', jpgname%to_char()


    call simple_touch(chunk%folder%to_char() // '/REJECTION_FINISHED')
    chunk%nptcls_selected    = spproj%os_ptcl2D%count_state_gt_zero()
    chunk%rejection_complete = .true.
    write(logfhandle,'(A,A,A,I6,A,I8,A,I8,A)') '>>> COMPLETED REJECTION FOR ', &
      label%to_char(), ' # ', chunk%id, ' : ', &
      chunk%nptcls_selected, '/', chunk%nptcls, ' PARTICLES SELECTED'

    call dealloc_imgarr(cavg_imgs)
    call spproj%kill()
    if( allocated(states) ) deallocate(states)
    call timer_stop(t0, string('reject_cavgs'))
  end subroutine reject_cavgs

  ! ============================================================================
  ! MODULE-LEVEL HELPERS
  ! ============================================================================

  ! Merges a list of source project files into chunk_project within the given
  ! output directory, updates projinfo and compenv from cline, kills stale 2D
  ! orientation sets and output oris, and removes intermediate files (FRCs,
  ! class-average stacks, sigma2 star file) carried over from the merge.
  ! Applies any SIMPLE_CHUNK_PARTITION environment override to the queue
  ! partition.
  subroutine merge_and_clear( projfiles, outdir, chunk_project, cline )
    type(string),     intent(in)    :: projfiles(:)
    type(string),     intent(in)    :: outdir
    type(sp_project), intent(inout) :: chunk_project
    type(cmdline),    intent(inout) :: cline
    character(len=STDLEN) :: chunk_part_env
    integer               :: envlen
    call merge_chunk_projfiles(projfiles, outdir, chunk_project, write_proj=.false.)
    call chunk_project%update_projinfo(cline)
    call chunk_project%update_compenv(cline)
    call chunk_project%os_cls2D%kill()
    call chunk_project%os_cls3D%kill()
    call chunk_project%os_out%kill()
    if( file_exists(string(outdir%to_char() // '/' // FRCS_FILE))        ) call del_file(string(outdir%to_char() // '/' // FRCS_FILE))
    if( file_exists(string(outdir%to_char() // '/cavgs.mrc'))            ) call del_file(string(outdir%to_char() // '/cavgs.mrc'))
    if( file_exists(string(outdir%to_char() // '/cavgs_odd.mrc'))        ) call del_file(string(outdir%to_char() // '/cavgs_odd.mrc'))
    if( file_exists(string(outdir%to_char() // '/cavgs_even.mrc'))       ) call del_file(string(outdir%to_char() // '/cavgs_even.mrc'))
    if( file_exists(string(outdir%to_char() // '/sigma2_combined.star')) ) call del_file(string(outdir%to_char() // '/sigma2_combined.star'))
    call get_environment_variable('SIMPLE_CHUNK_PARTITION', chunk_part_env, envlen)
    if( envlen > 0 ) call chunk_project%compenv%set(1, 'qsys_partition', trim(chunk_part_env))
  end subroutine merge_and_clear

  ! Returns a copy of arr with duplicate string entries removed, preserving
  ! the order of first occurrence.
  function remove_duplicates( arr ) result( unique )
    type(string), allocatable, intent(in) :: arr(:)
    type(string), allocatable             :: unique(:)
    integer :: i, j, n, new_size
    logical :: is_duplicate
    n        = size(arr)
    new_size = 0
    allocate(unique(n))
    do i = 1, n
      is_duplicate = .false.
      do j = 1, new_size
        if( arr(i) == unique(j) ) then
          is_duplicate = .true.
          exit
        end if
      end do
      if( .not. is_duplicate ) then
        new_size         = new_size + 1
        unique(new_size) = arr(i)
      end if
    end do
    unique = unique(1 : new_size)
  end function remove_duplicates

  ! Starts a timer and returns the tick value.
  integer(timer_int_kind) function timer_start()
    timer_start = tic()
  end function timer_start

  ! Stops the timer and, when DEBUG is true, logs the elapsed time with the
  ! provided subroutine name to logfhandle.
  subroutine timer_stop( t0, routine_name )
    integer(timer_int_kind), intent(in) :: t0
    type(string),            intent(in) :: routine_name
    real(timer_int_kind) :: elapsed
    if( .not. DEBUG ) return
    elapsed = toc(t0)
    write(logfhandle,'(A,A,A,F8.1)') 'simple_ptcl_sieve->', routine_name%to_char(), ' execution time:', elapsed
    call flush(logfhandle)
  end subroutine timer_stop

end module simple_ptcl_sieve
