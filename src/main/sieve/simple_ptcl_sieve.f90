!@descr: multi-tier particle sieve with coarse/fine 2D chunking and rejection
!==============================================================================
! MODULE: simple_ptcl_sieve
!
! PURPOSE:
!   Manages chunked ab initio 2D processing across two tiers:
!   coarse chunks (pass 1) and fine chunks (pass 2), including queue
!   submission, completion polling, rejection, compatibility filtering,
!   and final project combination.
!
! TYPES:
!   chunk2D_state - Per-chunk state container (identity, particle counts,
!                   paths, command line, and lifecycle flags).
!   ptcl_sieve    - Orchestrator owning coarse/fine chunk arrays, queue
!                   environment, defaults, compatibility models, and counters.
!
! WORKFLOW:
!   1. new()                     - initialize object state, output dirs,
!                                  queue environment, and optional
!                                  compatibility pretraining.
!   2. cycle()                   - collect completions and reject,
!                                  generate coarse/fine chunks, submit work.
!   3. generate_chunks_coarse()  - build coarse chunks from project list.
!   4. generate_chunks_fine()    - merge eligible coarse chunks into fine chunks.
!   5. submit()                  - submit pending chunks (fine prioritized).
!   6. collect_and_reject()      - detect finished jobs and apply rejection.
!   7. combine_completed_chunks()- combine eligible completed outputs.
!
! SENTINEL FILES (per chunk directory):
!   ABINITIO2D_FINISHED - queue job completion marker.
!   REJECTION_FINISHED  - rejection stage completion marker.
!   COMPLETE            - chunk consumed/finalized marker.
!   REJECTION_FAILED    - rejection failure marker.
!
! ENVIRONMENT:
!   SIMPLE_CHUNK_PARTITION - optional queue partition override.
!==============================================================================
module simple_ptcl_sieve
  use unix,                               only: c_time, c_long
  use simple_defs,                        only: logfhandle, STDLEN, CWD_GLOB, JPEG_DIM
  use simple_error,                       only: simple_exception
  use simple_image,                       only: image
  use simple_timer,                       only: timer_int_kind, tic, toc
  use simple_fileio,                      only: swap_suffix, simple_copy_file, write_filetable, simple_touch, basename, simple_list_files
  use simple_string,                      only: string
  use simple_syslib,                      only: simple_mkdir, simple_abspath, simple_chdir, &
                                                simple_getcwd, file_exists, del_file, dir_exists
  use simple_cmdline,                     only: cmdline
  use simple_qsys_env,                    only: qsys_env
  use simple_rec_list,                    only: rec_list
  use simple_image_bin,                   only: image_bin
  use simple_gui_utils,                   only: mrc2jpeg_tiled
  use simple_defs_fname,                  only: METADATA_EXT, ABINITIO2D_FINISHED, FRCS_FILE, JPG_EXT, MRC_EXT
  use simple_parameters,                  only: parameters
  use simple_sp_project,                  only: sp_project
  use simple_imgarr_utils,                only: read_cavgs_into_imgarr, dealloc_imgarr, write_imgarr
  use simple_string_utils,                only: int2str
  use simple_projfile_utils,              only: merge_chunk_projfiles
  use simple_commanders_cavgs,            only: commander_cluster_cavgs
  use simple_cavg_quality_model,          only: CAVG_QUALITY_MODEL_CHUNK_DEFAULT, cavg_quality_model
  use simple_cavg_quality_types,          only: cavg_quality_result, CAVG_QUALITY_CONTEXT_SIEVE, &
                                                CAVG_REJECT_REASON_POP, CAVG_REJECT_REASON_BAD_PIXELS, &
                                                CAVG_REJECT_REASON_NO_COMPONENT, CAVG_REJECT_REASON_MASK_GEOMETRY, &
                                                CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW, CAVG_REJECT_REASON_LOCVAR_FG_LOW, &
                                                CAVG_REJECT_REASON_FUZZY_BALL_SIGNAL_LOW
  use simple_class_compatibility,         only: class_compatibility, support_model_metrics, PREPROCESS_MORPH_SIZE
  use simple_cavg_quality_helpers,        only: cavg_rejection_reason_string
  use simple_cavg_quality_analysis,       only: evaluate_cavg_quality_hard_reject, evaluate_cavg_quality
  use simple_cavg_quality_feats,          only: I_BP40_100_CENTER_EDGE_VAR, SIEVE_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN
  use simple_segmentation,                only: otsu_img

  implicit none
  public :: ptcl_sieve
  public :: DEFAULT_COARSE_POP_THRESHOLD
  public :: DEFAULT_FINE_POP_THRESHOLD
  public :: DEFAULT_NCLS
  public :: DEFAULT_COARSE_BOX
  public :: DEFAULT_FINE_BOX
  public :: DEFAULT_COARSE_NSAMPLE
  public :: DEFAULT_FINE_NSAMPLE
  public :: DEFAULT_LPSTART
  public :: DEFAULT_COARSE_LP
  public :: DEFAULT_FINE_LP

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
  real,    parameter :: DEFAULT_LPSTART                       = 15.0
  real,    parameter :: DEFAULT_COARSE_LP                     = 15.0
  real,    parameter :: DEFAULT_FINE_LP                       = 10.0
  real,    parameter :: OVERFIT_CLUSTER_REJECT_FRAC           = 0.70

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
    logical                           :: coarse_only             = .false.
    logical                           :: pre_chunked             = .false.
    logical                           :: final_ingestion         = .false.
    logical                           :: model_rejection_enabled = .false.
    integer                           :: nparallel          = 1
    integer                           :: nthr               = 1
    integer                           :: nparts             = 1
    integer                           :: n_coarse_accepted_ptcls = 0
    integer                           :: n_coarse_rejected_ptcls = 0
    integer                           :: n_fine_accepted_ptcls   = 0
    integer                           :: n_fine_rejected_ptcls   = 0
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
    procedure :: get_n_coarse_accepted_ptcls
    procedure :: get_n_coarse_rejected_ptcls
    procedure :: get_n_fine_accepted_ptcls
    procedure :: get_n_fine_rejected_ptcls
    procedure :: get_n_accepted_ptcls
    procedure :: get_n_accepted_micrographs
    procedure :: get_n_rejected_ptcls
    procedure :: get_n_total_particles
    procedure :: get_latest
    procedure :: get_finished
    procedure :: set_final_ingestion
    procedure :: unset_final_ingestion
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
    procedure :: cleanup_chunk
  end type ptcl_sieve

contains

  ! --------------------------------------------------------------------------
  ! LIFECYCLE
  ! --------------------------------------------------------------------------

  ! Initializes a ptcl_sieve object from a parameters object: derives output
  ! directories from the current working directory, stores the concurrency
  ! limit, mask diameter, and thread count; imports any existing chunks from a
  ! previous run; creates all four output directories; allocates empty chunk
  ! arrays; and initialises the queue environment.
  subroutine new( self, params, completedir, pre_chunked )
    class(ptcl_sieve), intent(inout) :: self
    type(parameters),           intent(in)    :: params
    type(string),               intent(in)    :: completedir
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
    if( present(pre_chunked)        ) self%pre_chunked             = pre_chunked
    if( params%single_pass == 'yes' ) self%coarse_only             = .true.
    if( params%fine_model  == 'yes' ) self%model_rejection_enabled = .true.
    if( params%lpstart > 0.0 ) then
      if( self%coarse_defaults%lpstart /= params%lpstart ) self%coarse_defaults%lpstart = params%lpstart
      if( self%fine_defaults%lpstart   /= params%lpstart ) self%fine_defaults%lpstart   = params%lpstart
    end if
    if( params%lpstop_coarse  > 0.0 .and. self%coarse_defaults%lpstop   /= params%lpstop_coarse ) self%coarse_defaults%lpstop   = params%lpstop_coarse
    if( params%lpstop_fine    > 0.0 .and. self%fine_defaults%lpstop     /= params%lpstop_fine   ) self%fine_defaults%lpstop     = params%lpstop_fine
    if( params%box_coarse     > 0   .and. self%coarse_defaults%box_crop /= params%box_coarse    ) self%coarse_defaults%box_crop = params%box_coarse
    if( params%box_fine       > 0   .and. self%fine_defaults%box_crop   /= params%box_fine      ) self%fine_defaults%box_crop   = params%box_fine
    if( params%nsample_coarse > 0   .and. self%coarse_defaults%nsample  /= params%nsample_coarse) self%coarse_defaults%nsample  = params%nsample_coarse
    if( params%nsample_fine   > 0   .and. self%fine_defaults%nsample    /= params%nsample_fine  ) self%fine_defaults%nsample    = params%nsample_fine
    if( params%ncls_coarse    > 0   .and. self%coarse_defaults%ncls     /= params%ncls_coarse   ) self%coarse_defaults%ncls     = params%ncls_coarse
    if( params%ncls_fine      > 0   .and. self%fine_defaults%ncls       /= params%ncls_fine     ) self%fine_defaults%ncls       = params%ncls_fine
    self%outdir_chunks_coarse = string(cwd%to_char() // '/chunks_coarse')
    self%outdir_chunks_fine   = string(cwd%to_char() // '/chunks_fine')
    allocate(self%chunks_coarse(0))
    allocate(self%chunks_fine(0))
    ! Import existing chunks before creating directories.
    if( dir_exists(self%outdir_chunks_coarse) ) call self%import_existing_chunks_coarse()
    if( dir_exists(self%outdir_chunks_fine)   ) call self%import_existing_chunks_fine()
    call simple_mkdir(self%outdir_chunks_coarse)
    call simple_mkdir(self%outdir_chunks_fine)
    call self%qenv%new(params, 1, exec_bin=string('simple_exec'))
    call self%coarse_compatibility_model%new()
    call self%fine_compatibility_model%new()
    if( params%refs /= '' ) then
      if( file_exists(params%refs) ) then
        call self%coarse_compatibility_model%train(params%refs)
        call self%fine_compatibility_model%train(params%refs)
      else
        write(logfhandle,'(A,A)') '>>> WARNING: compatibility refs not found, skipping pretraining: ', params%refs%to_char()
      end if
    end if
    call timer_stop(t0, string('new'))
  end subroutine new

  ! Kills all live command lines for coarse and fine chunks, deallocates
  ! all chunk arrays,
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
    self%n_coarse_accepted_ptcls = 0
    self%n_coarse_rejected_ptcls = 0
    self%n_fine_accepted_ptcls   = 0
    self%n_fine_rejected_ptcls   = 0
    self%n_accepted_ptcls        = 0
    self%n_accepted_mics         = 0
    self%n_rejected_ptcls        = 0
    self%coarse_only             = .false.
    self%pre_chunked             = .false.
    self%final_ingestion         = .false.
    self%model_rejection_enabled = .false.
    call timer_stop(t0, string('kill'))
  end subroutine kill

  ! Sets the module-level final-ingestion flag. Once enabled, fine-tier final
  ! flushing logic is allowed to run without waiting on timeout.
  subroutine set_final_ingestion(self)
    class(ptcl_sieve), intent(inout) :: self
    self%final_ingestion = .true.
  end subroutine set_final_ingestion

  ! Clears the module-level final-ingestion flag.
  subroutine unset_final_ingestion(self)
    class(ptcl_sieve), intent(inout) :: self
    self%final_ingestion = .false.
  end subroutine unset_final_ingestion

  ! --------------------------------------------------------------------------
  ! IMPORT FROM PREVIOUS RUN
  ! --------------------------------------------------------------------------

  ! Scans the coarse output directory for existing chunk subdirectories and
  ! populates the coarse array with chunk records whose state is inferred from
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

  ! Scans the fine output directory for existing chunk subdirectories and
  ! populates the fine array with chunk records whose state is inferred from
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

  ! Returns the total number of coarse chunks; zero if unallocated.
  pure integer function get_n_chunks_coarse( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_chunks_coarse = 0
    if( allocated(self%chunks_coarse) ) get_n_chunks_coarse = size(self%chunks_coarse)
  end function get_n_chunks_coarse

  ! Returns the total number of fine chunks; zero if unallocated.
  pure integer function get_n_chunks_fine( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_chunks_fine = 0
    if( allocated(self%chunks_fine) ) get_n_chunks_fine = size(self%chunks_fine)
  end function get_n_chunks_fine

  ! Returns the combined number of currently running coarse and fine chunks.
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

  ! Returns the total selected-particle count across all coarse chunks
  ! that have completed rejection but have not yet been consumed into a
  ! fine chunk.
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

  ! Returns the total selected-particle count across all fine chunks
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

  ! Returns the cumulative number of coarse-tier accepted particles after
  ! coarse rejection and compatibility filtering.
  pure integer function get_n_coarse_accepted_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_coarse_accepted_ptcls = self%n_coarse_accepted_ptcls
  end function get_n_coarse_accepted_ptcls

  ! Returns the cumulative number of coarse-tier rejected particles after
  ! coarse rejection and compatibility filtering.
  pure integer function get_n_coarse_rejected_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_coarse_rejected_ptcls = self%n_coarse_rejected_ptcls
  end function get_n_coarse_rejected_ptcls

  ! Returns the cumulative number of fine-tier accepted particles after
  ! fine rejection and compatibility filtering.
  pure integer function get_n_fine_accepted_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_fine_accepted_ptcls = self%n_fine_accepted_ptcls
  end function get_n_fine_accepted_ptcls

  ! Returns the cumulative number of fine-tier rejected particles after
  ! fine rejection and compatibility filtering.
  pure integer function get_n_fine_rejected_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_fine_rejected_ptcls = self%n_fine_rejected_ptcls
  end function get_n_fine_rejected_ptcls

  ! Returns the cumulative number of particles accepted (state > 0) across all
  ! finalised fine chunks. Updated by collect_and_reject when each fine chunk
  ! is marked complete. Coarse rejections are not included; use
  ! get_n_pass_1/2_non_rejected_ptcls for unconsumed per-chunk selected counts.
  pure integer function get_n_accepted_ptcls( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_accepted_ptcls = self%n_accepted_ptcls
  end function get_n_accepted_ptcls

  ! Returns the cumulative number of accepted micrographs across all
  ! finalised fine chunks.
  pure integer function get_n_accepted_micrographs( self )
    class(ptcl_sieve), intent(in) :: self
    get_n_accepted_micrographs = self%n_accepted_mics
  end function get_n_accepted_micrographs

  ! Returns the cumulative number of particles rejected (state = 0) across all
  ! finalised fine chunks. Updated alongside n_accepted_ptcls by
  ! collect_and_reject when each fine chunk is marked complete.
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
  ! In coarse_only mode, only coarse chunks are required.
  ! In two-tier mode, coarse must be complete and fine must be complete
  ! when present; if no fine chunks were generated, coarse completion is
  ! treated as terminal.
  logical function get_finished( self )
    class(ptcl_sieve), intent(in) :: self
    ! Coarse tier: must have at least one chunk and all must be complete.
    get_finished = self%get_n_chunks_coarse() > 0
    if( .not. get_finished ) return
    get_finished = all(self%chunks_coarse(:)%complete .or. self%chunks_coarse(:)%failed)
    if( .not. get_finished ) return
    ! Coarse-only runs terminate once all coarse chunks are done.
    if( self%coarse_only ) return
    ! Fine tier: when absent, coarse completion is terminal.
    if( self%get_n_chunks_fine() == 0 )then
      get_finished = .true.
      return
    endif
    get_finished = all(self%chunks_fine(:)%complete .or. self%chunks_fine(:)%failed)
  end function get_finished

  ! --------------------------------------------------------------------------
  ! APPEND HELPERS
  ! --------------------------------------------------------------------------

  ! Appends a single chunk2D_state to the end of the coarse chunk array.
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

  ! Appends a single chunk2D_state to the end of the fine chunk array.
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

  ! Partitions un-included records from project_list into coarse chunks,
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
    ! project_list instead of creating new coarse chunks from particle slices.
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

    ! NOTE: coarse chunking uses a soft cap. The trigger micrograph that crosses the
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
    ! into a coarse chunk flagged as rejection-complete so the fine-tier final
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

  ! Merges rejection-complete coarse chunks into fine chunks,
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

      ! Flag this as the final sieve chunk once final ingestion is signaled and all coarse chunks are done.
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

    ! Flush a final small fine chunk when final ingestion is signaled.
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
        ! Flag this as the final sieve chunk once final ingestion is signaled and all coarse chunks are done.
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

  ! Populates the abinitio2D command line for a coarse chunk with:
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

  ! Populates the abinitio2D command line for a fine chunk.
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

  ! Merges completed chunk project files into a single combined project
  ! written to completedir. In coarse_only mode, combines rejection-complete
  ! coarse chunks; otherwise combines complete fine chunks.
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
  ! concurrency limit. Priority order: fine > coarse.
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

    ! Submit fine chunks first.
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

    ! Submit coarse chunks last (lowest priority).
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

  ! Polls all running coarse and fine chunks for ABINITIO2D_FINISHED,
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
          chunk%complete = .true.
          write(logfhandle,'(A,I6)') '>>> FINALISED COARSE CHUNK # ', chunk%id
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
          chunk%complete = .true.
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
    type(commander_cluster_cavgs)      :: cluster_cavg_commander
    type(support_model_metrics)        :: compat_metrics
    type(cavg_quality_result)          :: quality
    type(cavg_quality_model)           :: model
    type(sp_project)                   :: spproj
    type(sp_project)                   :: spproj_cluster
    type(cmdline)                      :: cline_cluster_cavgs
    type(string)                       :: stkname, jpgname, cwd_before_cluster, cluster_projfile, cluster_projbase
    integer(timer_int_kind)            :: t0
    integer                            :: ncls, iimg, nout, non_zero_ptcls, n_total_ptcls
    integer,             allocatable   :: clusters(:), uniq_clusters(:)
    type(string),        allocatable   :: overlay_reasons(:)
    logical,             allocatable   :: overfit_metric_fail(:)
    integer                            :: nuniq, ic, cid, n_in_cluster, n_fail_cluster
    real                               :: overfit_fail_frac
    real                               :: smpd_dummy
    real                               :: bp_score
    logical                            :: bp_fail
    real, parameter                    :: SIEVE_BP40_100_CENTER_EDGE_VAR_MIN_LOG = log(max(SIEVE_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN, tiny(1.0)))

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
      ! Map the coarse selection to the cluster-average level and write out the updated project file.
      states = spproj%os_cls2D%get_all_asint('state')
      call spproj%map_cavgs_selection(states)
      if( .true. ) then
        ! This branch is enabled for testing purposes. It would write out the updated project file after mapping the selection to cluster averages.
        call spproj%write()
        ! Run cluster-average quality evaluation
        call cline_cluster_cavgs%set('prg', 'cluster_cavgs')
        call cline_cluster_cavgs%set('projfile', chunk%projfile)
        call cline_cluster_cavgs%set('mskdiam', self%mskdiam)
        call cline_cluster_cavgs%set('skip_rejection', 'yes')
        call cline_cluster_cavgs%set('lp', self%coarse_defaults%lpstop)
        call cline_cluster_cavgs%set('nthr', 1)
        call cline_cluster_cavgs%set('mkdir', 'yes')
        call cline_cluster_cavgs%set('outdir', 'cluster_cavgs')
        call simple_getcwd(cwd_before_cluster)
        call simple_chdir(chunk%folder)
        CWD_GLOB = chunk%folder%to_char()
        call cluster_cavg_commander%execute(cline_cluster_cavgs)
        call simple_chdir(cwd_before_cluster)
        CWD_GLOB = cwd_before_cluster%to_char()
        call cwd_before_cluster%kill()
        call cline_cluster_cavgs%kill()

        cluster_projbase = basename(chunk%projfile)
        cluster_projfile = string(chunk%folder%to_char() // '/cluster_cavgs/' // cluster_projbase%to_char())
        if( file_exists(cluster_projfile) ) then
          call spproj_cluster%read(cluster_projfile)
          clusters = spproj_cluster%os_cls2D%get_all_asint('cluster')
          allocate(overfit_metric_fail(size(cavg_imgs)), source=.false.)
          do iimg = 1, size(cavg_imgs)
            if( quality%raw(iimg, I_BP40_100_CENTER_EDGE_VAR) < SIEVE_BP40_100_CENTER_EDGE_VAR_MIN_LOG) then
              overfit_metric_fail(iimg) = .true.
            end if
          end do
          
          nuniq = 0
          allocate(uniq_clusters(size(clusters)), source=0)
          do iimg = 1, size(clusters)
            if( clusters(iimg) <= 0 ) cycle
            if( .not. any(uniq_clusters(1:nuniq) == clusters(iimg)) ) then
              nuniq = nuniq + 1
              uniq_clusters(nuniq) = clusters(iimg)
            end if
          end do

          do ic = 1, nuniq
            cid = uniq_clusters(ic)
            n_in_cluster = count(clusters == cid)
            if( n_in_cluster == 0 ) cycle
            n_fail_cluster = count((clusters == cid) .and. overfit_metric_fail)
            overfit_fail_frac = real(n_fail_cluster) / real(n_in_cluster)
            if( n_in_cluster == 1 ) then
              do iimg = 1, size(clusters)
                if( clusters(iimg) == cid ) then
                  call spproj%os_cls2D%set(iimg, 'state', 1)
                  call spproj%os_cls2D%set(iimg, 'rejection_reason', '')
                end if
              end do
              write(logfhandle,'(A,I6,A,I6,A)') '>>> KEPT CLUSTER ', cid, ' IN COARSE CHUNK # ', chunk%id, &
                ' (SINGLE AVERAGE EXEMPT FROM OVERFIT CLUSTER GATE)'
            else if( overfit_fail_frac >= OVERFIT_CLUSTER_REJECT_FRAC ) then
              do iimg = 1, size(clusters)
                if( clusters(iimg) == cid ) then
                  call spproj%os_cls2D%set(iimg, 'state', 0)
                  call spproj%os_cls2D%set(iimg, 'rejection_reason', string('coarse_reject: overfit_cluster_frac_ge_0.70'))
                  bp_score    = quality%raw(iimg, I_BP40_100_CENTER_EDGE_VAR)
                  bp_fail     = bp_score    < SIEVE_BP40_100_CENTER_EDGE_VAR_MIN_LOG
                  write(logfhandle,'(A,I6,A,I6,A,F8.3,A,L1)') &
                    '>>> REJECTED CLUSTER CLASS chunk#', chunk%id, ' cls=', iimg, &
                    ' bp40_100=', bp_score, ' fail=', bp_fail
                end if
              end do
              write(logfhandle,'(A,I6,A,I6,A,I6,A,F6.3,A)') '>>> REJECTED CLUSTER ', cid, ' IN COARSE CHUNK # ', chunk%id, &
                ' DUE TO OVERFIT FAIL RATE ', n_fail_cluster, '/', overfit_fail_frac, ' (>= 0.70)'
                
            else
              do iimg = 1, size(clusters)
                if( clusters(iimg) == cid ) then
                  call spproj%os_cls2D%set(iimg, 'state', 1)
                  call spproj%os_cls2D%set(iimg, 'rejection_reason', '')
                end if
              end do
              write(logfhandle,'(A,I6,A,I6,A,F6.3,A)') '>>> KEPT CLUSTER ', cid, ' IN COARSE CHUNK # ', chunk%id, &
                ' OVERFIT FAIL FRAC ', overfit_fail_frac, ' (< 0.70)'
            end if
          end do

          if( allocated(clusters) ) deallocate(clusters)
          if( allocated(uniq_clusters) ) deallocate(uniq_clusters)
          if( allocated(overfit_metric_fail) ) deallocate(overfit_metric_fail)
          call spproj_cluster%kill()
        else
          write(logfhandle,'(A,A)') '>>> SKIPPING CLUSTER-LEVEL OVERFIT GATE; MISSING FILE: ', cluster_projfile%to_char()
        end if
      end if
      
      if( .not. self%coarse_compatibility_model%converged() ) then
        call self%coarse_compatibility_model%train(spproj)
        call self%coarse_compatibility_model%get_support_model_metrics(compat_metrics)
        write(logfhandle,'(A,I6,A,3(F10.4,1X),A,3(F10.4,1X),A,L1,A,L1,A,L1)') &
          '>>> COARSE COMPAT METRICS CHUNK # ', chunk%id, ' a/b/c=', &
          compat_metrics%axis_a, compat_metrics%axis_b, compat_metrics%axis_c, ' da/db/dc=', &
          compat_metrics%delta_a, compat_metrics%delta_b, compat_metrics%delta_c, ' valid=', &
          compat_metrics%valid, ' delta_valid=', compat_metrics%delta_valid, ' converged=', compat_metrics%converged
        if( compat_metrics%converged ) then
          write(logfhandle,'(A,I6)') '>>> COARSE COMPATIBILITY MODEL CONVERGED AT CHUNK # ', chunk%id
        end if
      end if
      call self%coarse_compatibility_model%infer(spproj)
    else
      if( self%model_rejection_enabled ) then
          call model%init_preset(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
          model%context = CAVG_QUALITY_CONTEXT_SIEVE
          call evaluate_cavg_quality_hard_reject(cavg_imgs, spproj%os_cls2D, 0.0, quality, CAVG_QUALITY_CONTEXT_SIEVE)
          call evaluate_cavg_quality(cavg_imgs, spproj%os_cls2D, 0.0, quality, model, CAVG_QUALITY_CONTEXT_SIEVE)
          call model%kill()
      else
          call evaluate_cavg_quality_hard_reject(cavg_imgs, spproj%os_cls2D, 0.0, quality, CAVG_QUALITY_CONTEXT_SIEVE)
      end if
      do iimg = 1, size(cavg_imgs)
          if( quality%states(iimg) == 0 ) then
              call spproj%os_cls2D%set(iimg, 'state', 0)
              call spproj%os_cls2D%set(iimg, 'rejection_reason', string('fine_reject: ')//cavg_rejection_reason_string(quality%reasons(iimg)))
          end if
      end do
      if( .not. self%fine_compatibility_model%converged() ) then
        call self%fine_compatibility_model%train(spproj)
        call self%fine_compatibility_model%get_support_model_metrics(compat_metrics)
        write(logfhandle,'(A,I6,A,3(F10.4,1X),A,3(F10.4,1X),A,L1,A,L1,A,L1)') &
          '>>> FINE COMPAT METRICS CHUNK # ', chunk%id, ' a/b/c=', &
          compat_metrics%axis_a, compat_metrics%axis_b, compat_metrics%axis_c, ' da/db/dc=', &
          compat_metrics%delta_a, compat_metrics%delta_b, compat_metrics%delta_c, ' valid=', &
          compat_metrics%valid, ' delta_valid=', compat_metrics%delta_valid, ' converged=', compat_metrics%converged
        if( compat_metrics%converged ) then
          write(logfhandle,'(A,I6)') '>>> FINE COMPATIBILITY MODEL CONVERGED AT CHUNK # ', chunk%id
        end if
      end if
      call self%fine_compatibility_model%infer(spproj)
    end if
    states = spproj%os_cls2D%get_all_asint('state')
    call spproj%map_cavgs_selection(states)
    non_zero_ptcls = spproj%os_ptcl2D%count_state_gt_zero()
    n_total_ptcls  = spproj%os_ptcl2D%get_noris()
    if( label == LABEL_COARSE ) then
      self%n_coarse_accepted_ptcls = self%n_coarse_accepted_ptcls + non_zero_ptcls
      self%n_coarse_rejected_ptcls = self%n_coarse_rejected_ptcls + n_total_ptcls - non_zero_ptcls
    else
      self%n_fine_accepted_ptcls   = self%n_fine_accepted_ptcls + non_zero_ptcls
      self%n_fine_rejected_ptcls   = self%n_fine_rejected_ptcls + n_total_ptcls - non_zero_ptcls
    end if
    call spproj%write(chunk%projfile)

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

    allocate(overlay_reasons(size(states)))
    do iimg = 1, size(states)
      overlay_reasons(iimg) = spproj%os_cls2D%get_str(iimg, 'rejection_reason')
    end do
    call write_rejection_reason_overlay_jpg(chunk, label, cavg_imgs, states, quality, overlay_reasons)
    deallocate(overlay_reasons)

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
      
    call self%cleanup_chunk(chunk, label)
    call dealloc_imgarr(cavg_imgs)
    call spproj%kill()
    if( allocated(states) ) deallocate(states)
    call timer_stop(t0, string('reject_cavgs'))
  end subroutine reject_cavgs

  subroutine cleanup_chunk( self, chunk, label )
    class(ptcl_sieve),   intent(inout) :: self
    type(chunk2D_state), intent(inout) :: chunk
    type(string),        intent(in)    :: label
    type(string), allocatable          :: files(:)
    type(string)                       :: fname_keep, fname_keep_last_iter_jpeg, fname_keep_last_iter_stk
    type(string)                       :: fname_keep_last_iter_even_stk, fname_keep_last_iter_odd_stk
    type(string)                       :: fname_keep_sigma_final, fname_keep_sigma_iter_1, fname_keep_sigma_iter_2, fname
    integer(timer_int_kind)            :: t0
    integer                            :: i, ndeleted, iter_idx, sigma_rank
    integer                            :: max_sigma_final_rank, max_sigma_iter_rank_1, max_sigma_iter_rank_2
    integer                            :: max_iter_even_stk, max_iter_odd_stk, max_iter_stk, max_iter_jpeg
    t0 = timer_start()

    call simple_list_files(chunk%folder%to_char() // '/*', files)
    if( .not. allocated(files) ) then
      call timer_stop(t0, string('cleanup_chunk'))
      return
    end if

    fname_keep                = basename(chunk%projfile)
    fname_keep_last_iter_jpeg = ''
    fname_keep_last_iter_stk  = ''
    fname_keep_last_iter_even_stk = ''
    fname_keep_last_iter_odd_stk  = ''
    fname_keep_sigma_final    = ''
    fname_keep_sigma_iter_1   = ''
    fname_keep_sigma_iter_2   = ''
    max_iter_jpeg             = -1
    max_iter_stk              = -1
    max_iter_even_stk         = -1
    max_iter_odd_stk          = -1
    max_sigma_final_rank      = -1
    max_sigma_iter_rank_1     = -1
    max_sigma_iter_rank_2     = -1
    do i = 1, size(files)
      fname = basename(files(i))

      iter_idx = iter_jpeg_index(fname)
      if( iter_idx > max_iter_jpeg ) then
        max_iter_jpeg             = iter_idx
        fname_keep_last_iter_jpeg = fname
      end if

      iter_idx = iter_stack_index(fname)
      if( .not. (fname%has_substr('_even') .or. fname%has_substr('_odd')) ) then
        if( iter_idx > max_iter_stk ) then
          max_iter_stk             = iter_idx
          fname_keep_last_iter_stk = fname
        end if
      end if
      if( fname%has_substr('_even') ) then
        if( iter_idx > max_iter_even_stk ) then
          max_iter_even_stk         = iter_idx
          fname_keep_last_iter_even_stk = fname
        end if
      end if
      if( fname%has_substr('_odd') ) then
        if( iter_idx > max_iter_odd_stk ) then
          max_iter_odd_stk         = iter_idx
          fname_keep_last_iter_odd_stk = fname
        end if
      end if

      sigma_rank = sigma_file_rank(fname)
      if( sigma_rank >= 1000000 ) then
        if( sigma_rank > max_sigma_final_rank ) then
          max_sigma_final_rank = sigma_rank
          fname_keep_sigma_final = fname
        end if
      else if( sigma_rank > 0 ) then
        if( sigma_rank > max_sigma_iter_rank_1 ) then
          max_sigma_iter_rank_2   = max_sigma_iter_rank_1
          fname_keep_sigma_iter_2 = fname_keep_sigma_iter_1
          max_sigma_iter_rank_1   = sigma_rank
          fname_keep_sigma_iter_1 = fname
        else if( sigma_rank > max_sigma_iter_rank_2 ) then
          max_sigma_iter_rank_2   = sigma_rank
          fname_keep_sigma_iter_2 = fname
        end if
      end if
    end do

    ndeleted   = 0
    do i = 1, size(files)
      fname = basename(files(i))
      if( fname == fname_keep ) cycle
      if( fname == string(ABINITIO2D_FINISHED) ) cycle
      if( fname == string('REJECTION_FINISHED') ) cycle
      if( fname == string('COMPLETE') ) cycle
      if( fname == string(REJECTION_FAILED) ) cycle
      if( fname == string(FRCS_FILE) ) cycle
      if( fname_keep_last_iter_jpeg%strlen() > 0 .and. fname == fname_keep_last_iter_jpeg ) cycle
      if( fname_keep_last_iter_stk%strlen()  > 0 .and. fname == fname_keep_last_iter_stk  ) cycle
      if( fname_keep_last_iter_even_stk%strlen() > 0 .and. fname == fname_keep_last_iter_even_stk ) cycle
      if( fname_keep_last_iter_odd_stk%strlen()  > 0 .and. fname == fname_keep_last_iter_odd_stk  ) cycle
      if( fname_keep_sigma_final%strlen()    > 0 .and. fname == fname_keep_sigma_final    ) cycle
      if( fname_keep_sigma_iter_1%strlen()   > 0 .and. fname == fname_keep_sigma_iter_1   ) cycle
      if( fname_keep_sigma_iter_2%strlen()   > 0 .and. fname == fname_keep_sigma_iter_2   ) cycle
      if( fname%has_substr('_selected.jpg') .or. fname%has_substr('_rejected.jpg') ) cycle
      if( fname%has_substr('_selected.jpeg') .or. fname%has_substr('_rejected.jpeg') ) cycle
      if( fname%has_substr('_all_reasons') ) cycle
      if( file_exists(files(i)) ) then
        call del_file(files(i))
        ndeleted = ndeleted + 1
      end if
    end do

    if( ndeleted > 0 ) then
      write(logfhandle,'(A,A,A,I6,A)') '>>> CLEANED ', label%to_char(), ' # ', chunk%id, ' (REMOVED NON-ESSENTIAL FILES)'
    end if

    call timer_stop(t0, string('cleanup_chunk'))
 
  contains

    ! Returns the iteration index for chunk-local JPEGs that follow the
    ! '*_iterNNN.jpg' or '*_iterNNN.jpeg' naming convention.
    ! Returns -1 when the file is not an iteration JPEG.
    pure integer function iter_jpeg_index( fname )
      type(string), intent(in) :: fname
      character(len=:), allocatable :: s
      integer :: n, p_iter, i, ch

      iter_jpeg_index = -1
      s = fname%to_char()
      n = len_trim(s)
      if( n < 4 ) return
      if( .not. ((n >= 4 .and. s(n-3:n) == '.jpg') .or. (n >= 5 .and. s(n-4:n) == '.jpeg')) ) return

      p_iter = index(s, '_iter', back=.true.)
      if( p_iter == 0 ) return
      i = p_iter + 5
      if( i > n ) return

      do while( i <= n )
        ch = iachar(s(i:i))
        if( ch < iachar('0') .or. ch > iachar('9') ) exit
        if( iter_jpeg_index < 0 ) iter_jpeg_index = 0
        iter_jpeg_index = 10 * iter_jpeg_index + ch - iachar('0')
        i = i + 1
      end do

      if( i == p_iter + 5 ) iter_jpeg_index = -1
    end function iter_jpeg_index

    ! Returns the iteration index for chunk-local stack files that follow the
    ! '*_iterNNN.mrc' or '*_iterNNN.mrcs' naming convention.
    ! Returns -1 when the file is not an iteration stack.
    pure integer function iter_stack_index( fname )
      type(string), intent(in) :: fname
      character(len=:), allocatable :: s
      integer :: n, p_iter, i, ch

      iter_stack_index = -1
      s = fname%to_char()
      n = len_trim(s)
      if( n < 4 ) return
      if( .not. ((n >= 4 .and. s(n-3:n) == '.mrc') .or. (n >= 5 .and. s(n-4:n) == '.mrcs')) ) return

      p_iter = index(s, '_iter', back=.true.)
      if( p_iter == 0 ) return
      i = p_iter + 5
      if( i > n ) return

      do while( i <= n )
        ch = iachar(s(i:i))
        if( ch < iachar('0') .or. ch > iachar('9') ) exit
        if( iter_stack_index < 0 ) iter_stack_index = 0
        iter_stack_index = 10 * iter_stack_index + ch - iachar('0')
        i = i + 1
      end do

      if( i == p_iter + 5 ) iter_stack_index = -1
    end function iter_stack_index

    ! Returns a rank for sigma candidate files.
    ! -1: not a sigma .star file
    ! >0 and <1000000: sigma .star with _iterNNN or _it_NNN (rank = NNN + 1)
    ! 1000000: sigma .star without an iteration suffix
    ! 2000000: sigma .star explicitly marked as final/combined
    pure integer function sigma_file_rank( fname )
      type(string), intent(in) :: fname
      character(len=:), allocatable :: s
      integer :: n, p_iter, i, ch, iter_val, ndigits

      sigma_file_rank = -1
      s = fname%to_char()
      n = len_trim(s)
      if( n < 5 ) return
      if( s(n-4:n) /= '.star' ) return
      if( index(s, 'sigma') == 0 ) return

      ! Explicit final/combined files should win over all iter files.
      if( index(s, 'final') > 0 .or. index(s, 'combined') > 0 ) then
        sigma_file_rank = 2000000
        return
      end if

      p_iter = index(s, '_iter', back=.true.)
      if( p_iter > 0 ) then
        i = p_iter + 5
      else
        p_iter = index(s, '_it_', back=.true.)
        if( p_iter > 0 ) then
          i = p_iter + 4
        else
          sigma_file_rank = 1000000
          return
        end if
      end if
      if( i > n ) then
        sigma_file_rank = 1000000
        return
      end if

      iter_val = 0
      ndigits = 0
      do while( i <= n )
        ch = iachar(s(i:i))
        if( ch < iachar('0') .or. ch > iachar('9') ) exit
        iter_val = 10 * iter_val + ch - iachar('0')
        ndigits = ndigits + 1
        i = i + 1
      end do

      if( ndigits > 0 ) then
        sigma_file_rank = iter_val + 1
      else
        sigma_file_rank = 1000000
      end if
    end function sigma_file_rank
    
  end subroutine cleanup_chunk

  ! Writes a third JPEG containing all class averages with reason-coded borders
  ! and emits a sidecar key file documenting the color mapping.
  subroutine write_rejection_reason_overlay_jpg( chunk, label, cavg_imgs, states, quality, rejection_reasons )
    type(chunk2D_state), intent(in)    :: chunk
    type(string),        intent(in)    :: label
    type(image),         intent(inout) :: cavg_imgs(:)
    integer,             intent(in)    :: states(:)
    type(cavg_quality_result), intent(in) :: quality
    type(string),        intent(in)    :: rejection_reasons(:)
    type(image)                      :: tile_img, work_img
    type(string)                     :: jpgname, keyname
    real, allocatable                :: rgb_map(:,:,:)
    integer                          :: ncls, xtiles, ytiles
    integer                          :: icls, ix, iy, ldim(3), borderw, g8
    integer                          :: x_start, x_end, y_start, y_end, iu
    integer                          :: reason_code

    ncls = size(cavg_imgs)
    if( ncls == 0 ) return

    xtiles = max(1, floor(sqrt(real(ncls))))
    ytiles = ceiling(real(ncls) / real(xtiles))
    call tile_img%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], cavg_imgs(1)%get_smpd())

    ldim = cavg_imgs(1)%get_ldim()
    call work_img%new(ldim, cavg_imgs(1)%get_smpd())

    do icls = 1, ncls
      call work_img%copy(cavg_imgs(icls))
      call work_img%fft
      if( ldim(1) > JPEG_DIM ) then
        call work_img%clip_inplace([JPEG_DIM, JPEG_DIM, 1])
      else
        call work_img%pad_inplace([JPEG_DIM, JPEG_DIM, 1], backgr=0., antialiasing=.false.)
      end if
      call work_img%ifft
      ix = mod(icls - 1, xtiles) + 1
      iy = (icls - 1) / xtiles + 1
      call tile_img%tile(work_img, ix, iy)
    end do

    rgb_map = tile_img%get_rmat()
    rgb_map(:,:,1) = max(0.0, min(255.0, rgb_map(:,:,1))) / 255.0

    ! Encode all class pixels as neutral RGB so the tile content stays grayscale
    ! while allowing colored overlays on selected border pixels.
    do iy = 1, size(rgb_map,2)
      do ix = 1, size(rgb_map,1)
        g8 = nint(max(0.0, min(1.0, rgb_map(ix,iy,1))) * 255.0)
        rgb_map(ix,iy,1) = rgb_code(g8, g8, g8)
      end do
    end do

    borderw = 3
    do icls = 1, ncls
      ix = mod(icls - 1, xtiles) + 1
      iy = (icls - 1) / xtiles + 1
      x_start = (ix - 1) * JPEG_DIM + 1
      x_end   = ix * JPEG_DIM
      y_start = (iy - 1) * JPEG_DIM + 1
      y_end   = iy * JPEG_DIM

      call paint_mask_outline(rgb_map, cavg_imgs(icls), x_start, y_start)

      if( states(icls) > 0 ) then
        reason_code = -1
      else if( rejection_reasons(icls)%strlen() > 0 ) then
        reason_code = reason_code_from_text(rejection_reasons(icls))
      else if( quality%states(icls) == 0 ) then
        reason_code = quality%reasons(icls)
      else
        reason_code = 100
      end if
      call paint_reason_border(rgb_map, x_start, x_end, y_start, y_end, borderw, reason_code)
    end do

    call draw_reason_key_swatches(rgb_map)
    ! Force deterministic [0,1] range for RGB packing in simple_jpg.
    rgb_map(1,1,1) = 0.0
    rgb_map(1,2,1) = 1.0

    call tile_img%set_rmat(rgb_map, .false.)
    jpgname = swap_suffix(chunk%projfile, '_all_reasons'//JPG_EXT, METADATA_EXT)
    call tile_img%write_jpg(jpgname, colorspec=3)
    write(logfhandle,'(A,A)') '>>> JPEG ', jpgname%to_char()

    keyname = string(jpgname%to_char() // '.key.txt')
    open(newunit=iu, file=keyname%to_char(), status='replace', action='write')
    write(iu,'(A)') 'Rejection Reason Color Key'
    write(iu,'(A)') 'Selected (state>0): green'
    write(iu,'(A)') 'Compatibility reject: orange'
    write(iu,'(A)') 'POP: red'
    write(iu,'(A)') 'BAD_PIXELS: magenta'
    write(iu,'(A)') 'NO_COMPONENT: cyan'
    write(iu,'(A)') 'MASK_GEOMETRY: yellow'
    write(iu,'(A)') 'BP_CENTER_EDGE_LOW: blue'
    write(iu,'(A)') 'LOCVAR_FG_LOW: purple'
    write(iu,'(A)') 'FUZZY_BALL_SIGNAL_NEG: brown'
    write(iu,'(A)') 'COARSE_OVERFIT_CLUSTER: dark orange'
    close(iu)
    write(logfhandle,'(A,A)') '>>> KEY  ', keyname%to_char()

    if( allocated(rgb_map) ) deallocate(rgb_map)
    call tile_img%kill()
    call work_img%kill()

  contains

    pure real function rgb_code(r, g, b)
      integer, intent(in) :: r, g, b
      integer :: packed
      packed = ishft(max(0, min(255, r)), 16) + ishft(max(0, min(255, g)), 8) + max(0, min(255, b))
      rgb_code = real(packed) / 16777215.0
    end function rgb_code

    pure real function reason_color(code)
      integer, intent(in) :: code
      select case(code)
      case(-1)
        reason_color = rgb_code( 46, 204, 113) ! selected
      case(100)
        reason_color = rgb_code(255, 165,   0) ! compatibility reject
      case(101)
        reason_color = rgb_code(204, 102,   0) ! coarse overfit cluster reject
      case(CAVG_REJECT_REASON_POP)
        reason_color = rgb_code(220,  20,  60)
      case(CAVG_REJECT_REASON_BAD_PIXELS)
        reason_color = rgb_code(255,   0, 255)
      case(CAVG_REJECT_REASON_NO_COMPONENT)
        reason_color = rgb_code(  0, 200, 200)
      case(CAVG_REJECT_REASON_MASK_GEOMETRY)
        reason_color = rgb_code(255, 220,   0)
      case(CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW)
        reason_color = rgb_code( 30, 144, 255)
      case(CAVG_REJECT_REASON_LOCVAR_FG_LOW)
        reason_color = rgb_code(153,  50, 204)
      case(CAVG_REJECT_REASON_FUZZY_BALL_SIGNAL_LOW)
        reason_color = rgb_code(139,  69,  19)
      case default
        reason_color = rgb_code(180, 180, 180)
      end select
    end function reason_color

    integer function reason_code_from_text(reason_text)
      type(string), intent(in) :: reason_text
      reason_code_from_text = 100
      if( reason_text%has_substr('overfit_cluster_frac_ge_0.70') ) then
        reason_code_from_text = 101
      else if( reason_text%has_substr('low population') ) then
        reason_code_from_text = CAVG_REJECT_REASON_POP
      else if( reason_text%has_substr('bad pixels') ) then
        reason_code_from_text = CAVG_REJECT_REASON_BAD_PIXELS
      else if( reason_text%has_substr('no mask component') ) then
        reason_code_from_text = CAVG_REJECT_REASON_NO_COMPONENT
      else if( reason_text%has_substr('mask geometry issues') ) then
        reason_code_from_text = CAVG_REJECT_REASON_MASK_GEOMETRY
      else if( reason_text%has_substr('low band-pass center-edge variance') ) then
        reason_code_from_text = CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW
      else if( reason_text%has_substr('low local variance in foreground') ) then
        reason_code_from_text = CAVG_REJECT_REASON_LOCVAR_FG_LOW
      else if( reason_text%has_substr('low fuzzy-ball signal') ) then
        reason_code_from_text = CAVG_REJECT_REASON_FUZZY_BALL_SIGNAL_LOW
      end if
    end function reason_code_from_text

    pure subroutine paint_reason_border(rmap, xs, xe, ys, ye, width, code)
      real,    intent(inout) :: rmap(:,:,:)
      integer, intent(in)    :: xs, xe, ys, ye, width, code
      integer :: i, j
      real    :: col
      col = reason_color(code)
      do j = ys, ye
        do i = xs, min(xe, xs + width - 1)
          rmap(i,j,1) = col
        end do
        do i = max(xs, xe - width + 1), xe
          rmap(i,j,1) = col
        end do
      end do
      do i = xs, xe
        do j = ys, min(ye, ys + width - 1)
          rmap(i,j,1) = col
        end do
        do j = max(ys, ye - width + 1), ye
          rmap(i,j,1) = col
        end do
      end do
    end subroutine paint_reason_border

    ! Paint a 1px contour of the Otsu binary mask in light blue on one tile.
    subroutine paint_mask_outline(rmap, src_img, x_start, y_start)
      real,        intent(inout) :: rmap(:,:,:)
      type(image), intent(inout) :: src_img
      integer,     intent(in)    :: x_start, y_start
      type(image)                :: local_img, mask_img
      type(image_bin)            :: mask_bin, cc_img
      real, allocatable          :: mask_vals(:,:,:)
      real, allocatable          :: weight_vals(:,:,:)
      integer, allocatable       :: cc_sizes(:)
      integer                    :: ldim_local(3), i, j, imorph
      integer                    :: nx, ny
      real                       :: cx, cy, r2, r2max, radial_w
      logical                    :: boundary

      ldim_local = src_img%get_ldim()
      call local_img%copy(src_img)
      call local_img%zero_edgeavg()
      call local_img%bp(0., 10.)
      call local_img%fft
      if( ldim_local(1) > JPEG_DIM ) then
        call local_img%clip_inplace([JPEG_DIM, JPEG_DIM, 1])
      else
        call local_img%pad_inplace([JPEG_DIM, JPEG_DIM, 1], backgr=0., antialiasing=.false.)
      end if
      call local_img%ifft

      ! Center-weight intensities before Otsu so dark/bright edge backgrounds
      ! are less likely to dominate the selected foreground component.
      weight_vals = local_img%get_rmat()
      nx = size(weight_vals, 1)
      ny = size(weight_vals, 2)
      cx = 0.5 * real(nx + 1)
      cy = 0.5 * real(ny + 1)
      r2max = max(1.0, (0.5 * real(min(nx, ny)))**2)
      do j = 1, ny
        do i = 1, nx
          r2 = (real(i) - cx)**2 + (real(j) - cy)**2
          radial_w = 0.35 + 0.65 * max(0.0, 1.0 - min(1.0, r2 / r2max))
          weight_vals(i,j,1) = weight_vals(i,j,1) * radial_w
        end do
      end do

      call mask_img%copy(local_img)
      call mask_img%set_rmat(weight_vals, .false.)
      call otsu_img(mask_img)
      ! Apply 5 px morphological closing to the binary Otsu mask.
      call mask_bin%transfer2bimg(mask_img)
      do imorph = 1, PREPROCESS_MORPH_SIZE
          call mask_bin%dilate()
      end do
      do imorph = 1, PREPROCESS_MORPH_SIZE
          call mask_bin%erode()
      end do

      ! Keep only the largest connected foreground component in the mask.
      call mask_bin%find_ccs(cc_img)
      cc_sizes = cc_img%size_ccs()
      if( size(cc_sizes) > 0 .and. maxval(cc_sizes) > 0 )then
          call cc_img%cc2bin(maxloc(cc_sizes, dim=1))
          call mask_bin%copy_bimg(cc_img)
      end if
      
      mask_vals = mask_bin%get_rmat()

      do j = 1, JPEG_DIM
        do i = 1, JPEG_DIM
          if( mask_vals(i,j,1) < 0.5 ) cycle
          boundary = (i == 1) .or. (i == JPEG_DIM) .or. (j == 1) .or. (j == JPEG_DIM)
          if( .not. boundary ) then
            boundary = mask_vals(i-1,j,1) < 0.5 .or. mask_vals(i+1,j,1) < 0.5 .or. &
                       mask_vals(i,j-1,1) < 0.5 .or. mask_vals(i,j+1,1) < 0.5
          end if
          if( boundary ) rmap(x_start + i - 1, y_start + j - 1, 1) = rgb_code(120, 205, 255)
        end do
      end do
      if( allocated(weight_vals) ) deallocate(weight_vals)
      if( allocated(cc_sizes)  ) deallocate(cc_sizes)
      if( allocated(mask_vals) ) deallocate(mask_vals)
      call mask_bin%kill_bimg()
      call mask_img%kill()
      call local_img%kill()
    end subroutine paint_mask_outline

    pure subroutine draw_reason_key_swatches(rmap)
      real, intent(inout) :: rmap(:,:,:)
      integer, parameter  :: ncolors = 10, box = 12, gap = 4, margin = 6
      integer             :: i, x0, y0, x1, y1, c
      integer             :: codes(ncolors)
      codes = [ -1, 100, CAVG_REJECT_REASON_POP, CAVG_REJECT_REASON_BAD_PIXELS, &
                CAVG_REJECT_REASON_NO_COMPONENT, CAVG_REJECT_REASON_MASK_GEOMETRY, &
                CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW, CAVG_REJECT_REASON_LOCVAR_FG_LOW, &
                CAVG_REJECT_REASON_FUZZY_BALL_SIGNAL_LOW, 101 ]
      y0 = margin
      do i = 1, ncolors
        x0 = margin + (i - 1) * (box + gap)
        x1 = min(size(rmap,1), x0 + box - 1)
        y1 = min(size(rmap,2), y0 + box - 1)
        do c = x0, x1
          rmap(c, y0:y1, 1) = reason_color(codes(i))
        end do
      end do
    end subroutine draw_reason_key_swatches

  end subroutine write_rejection_reason_overlay_jpg

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
