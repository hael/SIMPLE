!@descr: multi-tier microchunk 2D classification driver (pass-1/pass-2/refchunk/match)
!==============================================================================
! MODULE: simple_microchunked2D
!
! PURPOSE:
!   Manages the creation, configuration, and submission of microchunks across
!   three processing tiers and a final reference chunk — sub-projects each
!   containing a bounded particle subset drawn from a larger project list —
!   for parallel ab initio 2D class averaging. Each tier merges or matches
!   surviving particles from the previous tier for a progressively
!   higher-resolution result.
!
! TYPES:
!   chunk2D    - Holds the identity, particle counts (total and selected),
!                folder path, project file path, command line, and lifecycle
!                state flags for a single chunk.
!   chunked2D  - Owns arrays of pass-1, pass-2, and match microchunks and a
!                single reference chunk, together with the queue environment,
!                output directories, reference stack path, box size, shared
!                processing parameters (threads, mask diameter, parallelism),
!                and a timer for execution-time logging.
!
! WORKFLOW:
!   1. new()                            — initialise a chunked2D object: derives
!                                         output directories, imports any
!                                         existing chunks from previous runs,
!                                         creates directories, and initialises
!                                         the queue environment.
!   2. cycle()                          — convenience wrapper: collects completed
!                                         jobs and runs rejection, generates
!                                         chunks at all tiers, then submits
!                                         pending chunks to the queue.
!   3. generate_microchunks_pass_1()    — partition un-included records from a
!                                         rec_list into pass-1 microchunks, each
!                                         capped at MICROCHUNK_P1_THRESHOLD
!                                         particles; write a project file per
!                                         chunk and update the imported-projects
!                                         file table.
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
!   6. generate_microchunks_match()     — once the reference chunk rejection is
!                                         complete, create one match microchunk
!                                         per rejection-complete pass-2 chunk by
!                                         copying its project file and setting
!                                         refs/box on the command line for
!                                         template-guided classification.
!   7. submit()                         — dispatch pending chunks to the queue up
!                                         to the nparallel concurrency limit;
!                                         priority order: refchunk > match >
!                                         pass-2 > pass-1.
!   8. collect_and_reject()             — poll running chunks for the
!                                         ABINITIO2D_FINISHED sentinel; transition
!                                         completed chunks and immediately run
!                                         class-average rejection on each newly
!                                         finished chunk. Captures refs and box
!                                         from the refchunk after its rejection.
!
! SENTINEL FILES (written to each chunk folder):
!   ABINITIO2D_FINISHED  — written by the queue job on completion.
!   REJECTION_FINISHED   — written by reject_cavgs on successful rejection.
!   COMPLETE             — written when a chunk is consumed into the next tier.
!
! CONSTANTS:
!   MICROCHUNK_P1_THRESHOLD — maximum particles per pass-1 microchunk       (5000)
!   MICROCHUNK_P2_THRESHOLD — maximum selected particles per pass-2 chunk   (8000)
!   REFCHUNK_THRESHOLD      — minimum selected particles to form a ref chunk(10000)
!   DEFAULT_NCLS            — default number of 2D classes                    (100)
!   DEFAULT_MICRO_P1_LP     — pass-1 microchunk low-pass stop cutoff, Å    (15.0)
!   DEFAULT_MICRO_P2_LP     — pass-2 microchunk low-pass stop cutoff, Å    (10.0)
!   DEFAULT_REF_LP          — reference / match chunk low-pass stop cutoff, Å(8.0)
!   DEFAULT_WALLTIME        — per-chunk job time limit, seconds             (1740)
!
! ENVIRONMENT:
!   SIMPLE_CHUNK_PARTITION — if set, overrides the queue partition used for
!                            all chunk job submission.
!
! DEPENDENCIES:
!   simple_defs, simple_image, simple_timer, simple_string, simple_fileio,
!   simple_syslib, simple_cmdline, simple_qsys_env, simple_rec_list,
!   simple_stack_io, simple_parameters, simple_sp_project, simple_defs_fname,
!   simple_string_utils, simple_imgarr_utils, simple_projfile_utils,
!   simple_stream_microchunk_utils
!==============================================================================
module simple_microchunked2D
  use simple_defs,         only: logfhandle, STDLEN, CWD_GLOB, COSMSKHALFWIDTH
  use simple_image,        only: image
  use simple_timer,        only: timer_int_kind, tic, toc
  use simple_string,       only: string
  use simple_fileio,       only: swap_suffix, simple_copy_file, write_filetable, simple_touch, basename
  use simple_syslib,       only: simple_mkdir, simple_abspath, simple_chdir, &
                                 simple_getcwd, file_exists, del_file, dir_exists
  use simple_cmdline,      only: cmdline
  use simple_qsys_env,     only: qsys_env
  use simple_rec_list,     only: rec_list, rec_iterator, project_rec
  use simple_stack_io,     only: stack_io
  use simple_gui_utils,    only: mrc2jpeg_tiled
  use simple_parameters,   only: parameters
  use simple_sp_project,   only: sp_project
  use simple_defs_fname,   only: METADATA_EXT, ABINITIO2D_FINISHED, FRCS_FILE, JPG_EXT, MRC_EXT
  use simple_string_utils, only: int2str
  use simple_imgarr_utils, only: read_cavgs_into_imgarr, dealloc_imgarr
  use simple_projfile_utils,          only: merge_chunk_projfiles
  use simple_cluster2D_rejector,      only: cluster2D_rejector
  use simple_stream_microchunk_utils, only: calc_rejection_score, reject_outliers, reject_auto, reject_basic

  implicit none
  public  :: microchunked2D
  private
#include "simple_local_flags.inc"

  logical, parameter :: DEBUG                   = .true.
  integer, parameter :: MICROCHUNK_P1_THRESHOLD = 5000  
  integer, parameter :: MICROCHUNK_P2_THRESHOLD = 8000
  integer, parameter :: REFCHUNK_THRESHOLD      = 10000 
  integer, parameter :: DEFAULT_NCLS            = 100
  integer, parameter :: DEFAULT_WALLTIME        = 29 * 60  ! 29 minutes in seconds
  real,    parameter :: DEFAULT_LPSTART         = 15.0
  real,    parameter :: DEFAULT_MICRO_P1_LP     = 15.0
  real,    parameter :: DEFAULT_MICRO_P2_LP     = 10.0
  real,    parameter :: DEFAULT_REF_LP          =  8.0

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

  type :: microchunked2D
    private
    type(chunk2D), allocatable  :: microchunks_pass_1(:)
    type(chunk2D), allocatable  :: microchunks_pass_2(:)
    type(chunk2D), allocatable  :: microchunks_match(:)
    integer,       allocatable  :: refs_jpeg_inds(:), refs_jpeg_pops(:), ref_selection(:)
    integer,       allocatable  :: match_jpeg_inds(:), match_jpeg_pops(:)
    real,          allocatable  :: refs_jpeg_res(:), match_jpeg_res(:) 
    type(chunk2D)               :: refchunk
    type(qsys_env)              :: qenv
    type(string)                :: outdir_microchunks_pass_1
    type(string)                :: outdir_microchunks_pass_2
    type(string)                :: outdir_microchunks_match
    type(string)                :: outdir_refchunk
    type(string)                :: completedir
    type(string)                :: refs
    type(string)                :: refs_jpeg
    type(string)                :: match_stk
    type(string)                :: match_jpeg
    integer                     :: refs_jpeg_xtiles = 0
    integer                     :: refs_jpeg_ytiles = 0
    integer                     :: match_jpeg_xtiles = 0
    integer                     :: match_jpeg_ytiles = 0
    integer                     :: nparallel = 1
    integer                     :: nthr      = 1
    integer                     :: box       = 0
    integer                     :: n_accepted_ptcls = 0
    integer                     :: n_rejected_ptcls = 0
    real                        :: mskdiam   = 0.0
  contains
    procedure :: new
    procedure :: kill
    procedure :: cycle
    procedure :: import_existing_microchunks_pass_1
    procedure :: import_existing_microchunks_pass_2
    procedure :: import_existing_microchunks_match
    procedure :: import_existing_refchunk
    procedure :: get_n_microchunks_pass_1
    procedure :: get_n_microchunks_pass_2
    procedure :: get_n_microchunks_match
    procedure :: get_n_chunks_running
    procedure :: get_n_pass_1_non_rejected_ptcls
    procedure :: get_n_pass_2_non_rejected_ptcls
    procedure :: get_n_accepted_ptcls
    procedure :: get_n_rejected_ptcls
    procedure :: get_finished
    procedure :: get_references
    procedure :: get_latest_match
    procedure :: get_reference_selection
    procedure :: append_microchunk_pass_1
    procedure :: append_microchunk_pass_2
    procedure :: append_microchunk_match
    procedure :: generate_microchunks_pass_1
    procedure :: generate_microchunk_pass_1_cline
    procedure :: generate_microchunks_pass_2
    procedure :: generate_microchunk_pass_2_cline
    procedure :: generate_microchunks_match
    procedure :: generate_microchunk_match_cline
    procedure :: generate_refchunk
    procedure :: generate_refchunk_cline
    procedure :: combine_completed_match_chunks
    procedure :: submit
    procedure :: collect_and_reject
    procedure :: reject_cavgs
  end type microchunked2D

contains

  ! --------------------------------------------------------------------------
  ! LIFECYCLE
  ! --------------------------------------------------------------------------

  ! Initialises a chunked2D object from a parameters object: derives output
  ! directories from the current working directory, stores the concurrency
  ! limit, mask diameter, and thread count; imports any existing chunks from a
  ! previous run; creates all four output directories; allocates empty chunk
  ! arrays; and initialises the queue environment.
  subroutine new( self, params, completedir )
    class(microchunked2D), intent(inout) :: self
    type(parameters),      intent(in)    :: params
    type(string),          intent(in)    :: completedir
    type(string)                         :: cwd
    integer(timer_int_kind)              :: t0
    t0 = timer_start()
    call self%kill()
    call simple_getcwd(cwd)
    self%completedir               = completedir
    self%nparallel                 = params%nchunks
    self%mskdiam                   = params%mskdiam
    self%nthr                      = params%nthr
    self%outdir_microchunks_pass_1 = string(cwd%to_char() // '/microchunks_pass_1')
    self%outdir_microchunks_pass_2 = string(cwd%to_char() // '/microchunks_pass_2')
    self%outdir_microchunks_match  = string(cwd%to_char() // '/microchunks_match')
    self%outdir_refchunk           = string(cwd%to_char() // '/refchunk')
    allocate(self%microchunks_pass_1(0))
    allocate(self%microchunks_pass_2(0))
    allocate(self%microchunks_match(0))
    ! Import completed chunks from a previous run before creating directories
    if( dir_exists(self%outdir_microchunks_pass_1) ) call self%import_existing_microchunks_pass_1()
    if( dir_exists(self%outdir_microchunks_pass_2) ) call self%import_existing_microchunks_pass_2()
    if( dir_exists(self%outdir_microchunks_match)  ) call self%import_existing_microchunks_match()
    if( dir_exists(self%outdir_refchunk)           ) call self%import_existing_refchunk()
    call simple_mkdir(self%outdir_microchunks_pass_1)
    call simple_mkdir(self%outdir_microchunks_pass_2)
    call simple_mkdir(self%outdir_microchunks_match)
    call simple_mkdir(self%outdir_refchunk)
    call self%qenv%new(params, 1, exec_bin=string('simple_exec'))
    call timer_stop(t0, string('new'))
  end subroutine new

  ! Kills all live command lines for pass-1, pass-2, match, and the reference
  ! chunk, then deallocates all three microchunk arrays, releasing all
  ! associated memory.
  subroutine kill( self )
    class(microchunked2D), intent(inout) :: self
    integer(timer_int_kind) :: t0
    integer                 :: i
    t0 = timer_start()
    if( allocated(self%microchunks_pass_1) ) then
      do i = 1, self%get_n_microchunks_pass_1()
        call self%microchunks_pass_1(i)%cline%kill()
      end do
      deallocate(self%microchunks_pass_1)
    end if
    if( allocated(self%microchunks_pass_2) ) then
      do i = 1, self%get_n_microchunks_pass_2()
        call self%microchunks_pass_2(i)%cline%kill()
      end do
      deallocate(self%microchunks_pass_2)
    end if
    if( allocated(self%microchunks_match) ) then
      do i = 1, self%get_n_microchunks_match()
        call self%microchunks_match(i)%cline%kill()
      end do
      deallocate(self%microchunks_match)
    end if
    call self%refchunk%cline%kill()
    call timer_stop(t0, string('kill'))
  end subroutine kill

  ! --------------------------------------------------------------------------
  ! IMPORT FROM PREVIOUS RUN
  ! --------------------------------------------------------------------------

  ! Scans the pass-1 output directory for existing chunk subdirectories and
  ! populates the pass-1 array with chunk records whose state is inferred from
  ! sentinel files (ABINITIO2D_FINISHED, REJECTION_FINISHED, COMPLETE).
  subroutine import_existing_microchunks_pass_1( self )
    class(microchunked2D), intent(inout) :: self
    type(sp_project)                     :: chunk_project
    type(chunk2D)                        :: new_chunk
    type(string)                         :: chunk_folder
    integer(timer_int_kind)              :: t0
    integer                              :: chunk_id
    t0       = timer_start()
    chunk_id = 0
    do
      chunk_id     = chunk_id + 1
      chunk_folder = string(self%outdir_microchunks_pass_1%to_char() // &
                            '/microchunk_pass_1_' // int2str(chunk_id))
      if( .not. dir_exists(chunk_folder) ) exit
      new_chunk%id       = chunk_id
      new_chunk%folder   = simple_abspath(chunk_folder)
      new_chunk%projfile = string(new_chunk%folder%to_char() // &
                                  '/microchunk_pass_1_' // int2str(chunk_id) // METADATA_EXT)
      if( .not. file_exists(new_chunk%projfile) ) then
        write(logfhandle,'(A,I6,A)') '>>> WARNING: missing project file for microchunk pass 1 #', chunk_id, ', skipping'
        cycle
      end if
      call chunk_project%read(new_chunk%projfile)
      new_chunk%nptcls              = chunk_project%os_ptcl2D%get_noris()
      new_chunk%nptcls_selected     = chunk_project%os_ptcl2D%count_state_gt_zero()
      new_chunk%abinitio2D_running  = .false.
      new_chunk%abinitio2D_complete = file_exists(new_chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED)
      new_chunk%rejection_complete  = file_exists(new_chunk%folder%to_char() // '/REJECTION_FINISHED')
      new_chunk%complete            = file_exists(new_chunk%folder%to_char() // '/COMPLETE')
      new_chunk%failed              = .false.
      call chunk_project%kill()
      call self%append_microchunk_pass_1(new_chunk)
      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> IMPORTED EXISTING MICROCHUNK PASS 1 # ', chunk_id, ' WITH ', new_chunk%nptcls, ' PARTICLES'
    end do
    call timer_stop(t0, string('import_existing_microchunks_pass_1'))
  end subroutine import_existing_microchunks_pass_1

  ! Scans the pass-2 output directory for existing chunk subdirectories and
  ! populates the pass-2 array with chunk records whose state is inferred from
  ! sentinel files.
  subroutine import_existing_microchunks_pass_2( self )
    class(microchunked2D), intent(inout) :: self
    type(sp_project)                     :: chunk_project
    type(chunk2D)                        :: new_chunk
    type(string)                         :: chunk_folder
    integer(timer_int_kind)              :: t0
    integer                              :: chunk_id
    t0       = timer_start()
    chunk_id = 0
    do
      chunk_id     = chunk_id + 1
      chunk_folder = string(self%outdir_microchunks_pass_2%to_char() // &
                            '/microchunk_pass_2_' // int2str(chunk_id))
      if( .not. dir_exists(chunk_folder) ) exit
      new_chunk%id       = chunk_id
      new_chunk%folder   = simple_abspath(chunk_folder)
      new_chunk%projfile = string(new_chunk%folder%to_char() // &
                                  '/microchunk_pass_2_' // int2str(chunk_id) // METADATA_EXT)
      if( .not. file_exists(new_chunk%projfile) ) then
        write(logfhandle,'(A,I6,A)') '>>> WARNING: missing project file for microchunk pass 2 #', chunk_id, ', skipping'
        cycle
      end if
      call chunk_project%read(new_chunk%projfile)
      new_chunk%nptcls              = chunk_project%os_ptcl2D%get_noris()
      new_chunk%nptcls_selected     = chunk_project%os_ptcl2D%count_state_gt_zero()
      new_chunk%abinitio2D_running  = .false.
      new_chunk%abinitio2D_complete = file_exists(new_chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED)
      new_chunk%rejection_complete  = file_exists(new_chunk%folder%to_char() // '/REJECTION_FINISHED')
      new_chunk%complete            = file_exists(new_chunk%folder%to_char() // '/COMPLETE')
      new_chunk%failed              = .false.
      call chunk_project%kill()
      call self%append_microchunk_pass_2(new_chunk)
      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> IMPORTED EXISTING MICROCHUNK PASS 2 # ', chunk_id, ' WITH ', new_chunk%nptcls, ' PARTICLES'
    end do
    call timer_stop(t0, string('import_existing_microchunks_pass_2'))
  end subroutine import_existing_microchunks_pass_2

  ! Scans the match output directory for existing chunk subdirectories and
  ! populates the match array with chunk records whose state is inferred from
  ! sentinel files.
  subroutine import_existing_microchunks_match( self )
    class(microchunked2D), intent(inout) :: self
    logical,               allocatable   :: cls_msk(:)
    type(sp_project)                     :: chunk_project
    type(chunk2D)                        :: new_chunk
    type(string)                         :: chunk_folder, stkname
    integer(timer_int_kind)              :: t0
    integer                              :: chunk_id, ncls
    real                                 :: smpd
    t0       = timer_start()
    chunk_id = 0
    do
      chunk_id     = chunk_id + 1
      chunk_folder = string(self%outdir_microchunks_match%to_char() // &
                            '/microchunk_match_' // int2str(chunk_id))
      if( .not. dir_exists(chunk_folder) ) exit
      new_chunk%id       = chunk_id
      new_chunk%folder   = simple_abspath(chunk_folder)
      new_chunk%projfile = string(new_chunk%folder%to_char() // &
                                  '/microchunk_match_' // int2str(chunk_id) // METADATA_EXT)
      if( .not. file_exists(new_chunk%projfile) ) then
        write(logfhandle,'(A,I6,A)') '>>> WARNING: missing project file for microchunk match #', chunk_id, ', skipping'
        cycle
      end if
      call chunk_project%read(new_chunk%projfile)
      call chunk_project%get_cavgs_stk(stkname, ncls, smpd)
      new_chunk%nptcls              = chunk_project%os_ptcl2D%get_noris()
      new_chunk%nptcls_selected     = chunk_project%os_ptcl2D%count_state_gt_zero()
      new_chunk%abinitio2D_running  = .false.
      new_chunk%abinitio2D_complete = file_exists(new_chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED)
      new_chunk%rejection_complete  = file_exists(new_chunk%folder%to_char() // '/REJECTION_FINISHED')
      new_chunk%complete            = file_exists(new_chunk%folder%to_char() // '/COMPLETE')
      new_chunk%failed              = .false.
      self%match_stk                = stkname
      if( new_chunk%complete ) then
        self%match_jpeg = swap_suffix(self%match_stk, JPG_EXT, MRC_EXT)
        call chunk_project%cavgs2jpg(self%match_jpeg_inds, self%match_jpeg, self%match_jpeg_xtiles, self%match_jpeg_ytiles)
        self%match_jpeg_pops = chunk_project%os_cls2D%get_all_asint('pop')
        self%match_jpeg_res  = chunk_project%os_cls2D%get_all('res')
        allocate(cls_msk, source=self%match_jpeg_inds /= 0)
        self%match_jpeg_inds = pack(self%match_jpeg_inds, cls_msk)
        self%match_jpeg_pops = pack(self%match_jpeg_pops, cls_msk)
        self%match_jpeg_res  = pack(self%match_jpeg_res,  cls_msk)
        call chunk_project%kill()
        deallocate(cls_msk)
      end if
      call chunk_project%kill()
      call self%append_microchunk_match(new_chunk)
      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> IMPORTED EXISTING MICROCHUNK MATCH # ', chunk_id, ' WITH ', new_chunk%nptcls, ' PARTICLES'
    end do
    call timer_stop(t0, string('import_existing_microchunks_match'))
  end subroutine import_existing_microchunks_match

  ! Reads the reference chunk project file from the refchunk output directory
  ! if it exists and its project file is present, and restores the refchunk
  ! state from sentinel files. No-op if the directory or project file is absent.
  subroutine import_existing_refchunk( self )
    class(microchunked2D), intent(inout) :: self
    logical,               allocatable   :: cls_msk(:)
    integer,               allocatable   :: states(:)
    type(sp_project)                     :: chunk_project
    type(string)                         :: stkname
    integer(timer_int_kind)              :: t0
    integer                              :: ncls
    real                                 :: smpd
    t0 = timer_start()
    if( .not. dir_exists(self%outdir_refchunk) ) return
    self%refchunk%id      = 1
    self%refchunk%folder  = simple_abspath(self%outdir_refchunk)
    self%refchunk%projfile = string(self%refchunk%folder%to_char() // '/refchunk' // METADATA_EXT)
    if( .not. file_exists(self%refchunk%projfile) ) return
    call chunk_project%read(self%refchunk%projfile)
    call chunk_project%get_cavgs_stk(stkname, ncls, smpd)
    self%refchunk%nptcls              = chunk_project%os_ptcl2D%get_noris()
    self%refchunk%nptcls_selected     = chunk_project%os_ptcl2D%count_state_gt_zero()
    self%refchunk%abinitio2D_running  = .false.
    self%refchunk%abinitio2D_complete = file_exists(self%refchunk%folder%to_char() // '/' // ABINITIO2D_FINISHED)
    self%refchunk%rejection_complete  = file_exists(self%refchunk%folder%to_char() // '/REJECTION_FINISHED')
    self%refchunk%complete            = file_exists(self%refchunk%folder%to_char() // '/COMPLETE')
    self%refchunk%failed              = .false.
    self%refs                         = stkname
    if( self%refchunk%complete ) then
        self%refs_jpeg = swap_suffix(self%refs, JPG_EXT, MRC_EXT)
        call chunk_project%cavgs2jpg(self%refs_jpeg_inds, self%refs_jpeg, self%refs_jpeg_xtiles, self%refs_jpeg_ytiles)
        self%ref_selection  = self%refs_jpeg_inds
        self%refs_jpeg_pops = chunk_project%os_cls2D%get_all_asint('pop')
        self%refs_jpeg_res  = chunk_project%os_cls2D%get_all('res')
        allocate(cls_msk, source=self%refs_jpeg_inds /= 0)
        self%refs_jpeg_inds = pack(self%refs_jpeg_inds, cls_msk)
        self%refs_jpeg_pops = pack(self%refs_jpeg_pops, cls_msk)
        self%refs_jpeg_res  = pack(self%refs_jpeg_res,  cls_msk)
        deallocate(cls_msk)
        states = chunk_project%os_cls2D%get_all_asint('state')
        allocate(cls_msk, source=states /= 0)
        self%ref_selection = pack(self%ref_selection, cls_msk)
        call chunk_project%kill()
        deallocate(cls_msk, states)
    end if
    call chunk_project%kill()
    write(logfhandle,'(A,I8,A)') &
      '>>> IMPORTED EXISTING REFCHUNK WITH ', self%refchunk%nptcls, ' PARTICLES'
    call timer_stop(t0, string('import_existing_refchunk'))
  end subroutine import_existing_refchunk

  ! --------------------------------------------------------------------------
  ! MAIN LOOP WRAPPER
  ! --------------------------------------------------------------------------

  ! Convenience wrapper for the main processing loop: collects completed jobs
  ! and runs rejection, generates new chunks at all tiers in tier order, then
  ! submits any pending chunks to the queue.
  subroutine cycle( self, project_list )
    class(microchunked2D), intent(inout) :: self
    type(rec_list),        intent(inout) :: project_list
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    call self%collect_and_reject()
    call self%generate_microchunks_pass_1(project_list)
    call self%generate_microchunks_pass_2()
    call self%generate_refchunk()
    call self%generate_microchunks_match()
    call self%submit()
    call timer_stop(t0, string('cycle'))
  end subroutine cycle

  ! --------------------------------------------------------------------------
  ! QUERIES
  ! --------------------------------------------------------------------------

  ! Returns the total number of pass-1 microchunks; zero if unallocated.
  pure integer function get_n_microchunks_pass_1( self )
    class(microchunked2D), intent(in) :: self
    get_n_microchunks_pass_1 = 0
    if( allocated(self%microchunks_pass_1) ) get_n_microchunks_pass_1 = size(self%microchunks_pass_1)
  end function get_n_microchunks_pass_1

  ! Returns the total number of pass-2 microchunks; zero if unallocated.
  pure integer function get_n_microchunks_pass_2( self )
    class(microchunked2D), intent(in) :: self
    get_n_microchunks_pass_2 = 0
    if( allocated(self%microchunks_pass_2) ) get_n_microchunks_pass_2 = size(self%microchunks_pass_2)
  end function get_n_microchunks_pass_2

  ! Returns the total number of match microchunks; zero if unallocated.
  pure integer function get_n_microchunks_match( self )
    class(microchunked2D), intent(in) :: self
    get_n_microchunks_match = 0
    if( allocated(self%microchunks_match) ) get_n_microchunks_match = size(self%microchunks_match)
  end function get_n_microchunks_match

  ! Returns the combined number of currently running pass-1, pass-2, match,
  ! and reference chunks.
  pure integer function get_n_chunks_running( self )
    class(microchunked2D), intent(in) :: self
    integer :: i
    get_n_chunks_running = 0
    do i = 1, self%get_n_microchunks_pass_1()
      if( self%microchunks_pass_1(i)%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
    end do
    do i = 1, self%get_n_microchunks_pass_2()
      if( self%microchunks_pass_2(i)%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
    end do
    do i = 1, self%get_n_microchunks_match()
      if( self%microchunks_match(i)%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
    end do
    if( self%refchunk%abinitio2D_running ) get_n_chunks_running = get_n_chunks_running + 1
  end function get_n_chunks_running

  ! Returns the total selected-particle count across all pass-1 microchunks
  ! that have completed rejection but have not yet been consumed into a
  ! pass-2 microchunk.
  pure integer function get_n_pass_1_non_rejected_ptcls( self )
    class(microchunked2D), intent(in) :: self
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
    class(microchunked2D), intent(in) :: self
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

  ! Returns the cumulative number of particles accepted (not rejected) across
  ! all completed rejection passes. Updated by collect_and_reject as each chunk
  ! finishes; distinct from get_n_pass_1/2_non_rejected_ptcls, which count
  ! only unconsumed per-chunk selected particles computed on the fly.
  pure integer function get_n_accepted_ptcls( self )
    class(microchunked2D), intent(in) :: self
    get_n_accepted_ptcls = self%n_accepted_ptcls
  end function get_n_accepted_ptcls

  ! Returns the cumulative number of particles rejected across all completed
  ! rejection passes. Updated alongside n_accepted_ptcls by collect_and_reject.
  pure integer function get_n_rejected_ptcls( self )
    class(microchunked2D), intent(in) :: self
    get_n_rejected_ptcls = self%n_rejected_ptcls
  end function get_n_rejected_ptcls

  ! Returns true only when every chunk across all four tiers is marked complete.
  ! Checks proceed in pipeline order; an incomplete earlier tier returns early.
  logical function get_finished( self )
    class(microchunked2D), intent(in) :: self
    ! pass-1: all ab-initio microchunks must be rejected
    get_finished = self%get_n_microchunks_pass_1() == 0 .or. &
                   all(self%microchunks_pass_1(:)%rejection_complete)    
    if( .not. get_finished ) return
    ! pass-2: all merged microchunks must be done
    get_finished = self%get_n_microchunks_pass_2() == 0 .or. &
                   all(self%microchunks_pass_2(:)%complete)  
    if( .not. get_finished ) return
    ! refchunk: the reference class-average chunk must be done
    get_finished = self%refchunk%complete
    if( .not. get_finished ) return
    ! match: all template-match chunks must be done
    get_finished = self%get_n_microchunks_match() == 0 .or. &
                   all(self%microchunks_match(:)%complete)   
  end function get_finished

  ! Returns .true. and populates the reference class-average outputs when the
  ! refchunk is complete and all reference arrays are allocated; returns .false.
  ! otherwise. On a .false. return, all intent(out) and intent(inout) outputs
  ! are left in a defined, safe state (deallocated / zero / empty).
  logical function get_references( self, jpeg_inds, jpeg_pops, jpeg_res, jpeg, stk, xtiles, ytiles )
    class(microchunked2D), intent(in)    :: self
    integer, allocatable,  intent(inout) :: jpeg_inds(:), jpeg_pops(:) ! class indices and populations in the JPEG tile order
    real,    allocatable,  intent(inout) :: jpeg_res(:)                 ! resolution (Angstrom) per class, same order
    type(string),          intent(out)   :: jpeg, stk                   ! paths to the JPEG contact sheet and MRC reference stack
    integer,               intent(out)   :: xtiles, ytiles              ! tile grid dimensions of the JPEG contact sheet
    ! release any prior allocations so the caller gets a clean result on .false. return
    if( allocated(jpeg_inds) ) deallocate(jpeg_inds)
    if( allocated(jpeg_pops) ) deallocate(jpeg_pops)
    if( allocated(jpeg_res)  ) deallocate(jpeg_res)
    jpeg   = ''
    stk    = ''
    xtiles = 0
    ytiles = 0
    ! references are only available once the refchunk has finished
    if( .not. self%refchunk%complete ) then
        get_references = .false.
        return
    end if
    ! guard against inconsistent state where refchunk is marked complete but
    ! the reference arrays have not yet been populated
    if( .not. (allocated(self%refs_jpeg_inds) .and. allocated(self%refs_jpeg_pops) .and. allocated(self%refs_jpeg_res)) ) then
        get_references = .false.
        return
    end if
    allocate(jpeg_inds, source=self%refs_jpeg_inds)
    allocate(jpeg_pops, source=self%refs_jpeg_pops)
    allocate(jpeg_res,  source=self%refs_jpeg_res)
    stk    = self%refs
    jpeg   = self%refs_jpeg
    xtiles = self%refs_jpeg_xtiles
    ytiles = self%refs_jpeg_ytiles
    get_references = .true.
  end function get_references

  logical function get_latest_match( self, jpeg_inds, jpeg_pops, jpeg_res, jpeg, stk, xtiles, ytiles )
    class(microchunked2D), intent(in)    :: self
    integer, allocatable,  intent(inout) :: jpeg_inds(:), jpeg_pops(:) ! class indices and populations in the JPEG tile order
    real,    allocatable,  intent(inout) :: jpeg_res(:)                 ! resolution (Angstrom) per class, same order
    type(string),          intent(out)   :: jpeg, stk                   ! paths to the JPEG contact sheet and MRC reference stack
    integer,               intent(out)   :: xtiles, ytiles              ! tile grid dimensions of the JPEG contact sheet
    ! release any prior allocations so the caller gets a clean result on .false. return
    if( allocated(jpeg_inds) ) deallocate(jpeg_inds)
    if( allocated(jpeg_pops) ) deallocate(jpeg_pops)
    if( allocated(jpeg_res)  ) deallocate(jpeg_res)
    jpeg   = ''
    stk    = ''
    xtiles = 0
    ytiles = 0
    if( self%match_jpeg == '' ) then
        get_latest_match = .false.
        return
    end if
    ! guard against inconsistent state where match_jpeg is set but the
    ! match arrays have not yet been populated
    if( .not. (allocated(self%match_jpeg_inds) .and. allocated(self%match_jpeg_pops) .and. allocated(self%match_jpeg_res)) ) then
        get_latest_match = .false.
        return
    end if
    allocate(jpeg_inds, source=self%match_jpeg_inds)
    allocate(jpeg_pops, source=self%match_jpeg_pops)
    allocate(jpeg_res,  source=self%match_jpeg_res )
    stk    = self%match_stk
    jpeg   = self%match_jpeg
    xtiles = self%match_jpeg_xtiles
    ytiles = self%match_jpeg_ytiles
    get_latest_match = .true.
  end function get_latest_match

  logical function get_reference_selection( self, selection )
    class(microchunked2D), intent(in)    :: self
    integer, allocatable,  intent(out)   :: selection(:)
    if( .not. allocated(self%ref_selection) ) then
      get_reference_selection = .false.
      return
    end if
    allocate(selection, source=self%ref_selection)
    get_reference_selection = .true.
  end function get_reference_selection

  ! --------------------------------------------------------------------------
  ! APPEND HELPERS
  ! --------------------------------------------------------------------------

  ! Appends a single chunk2D to the end of the pass-1 microchunk array.
  subroutine append_microchunk_pass_1( self, new_chunk )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D),    intent(in)    :: new_chunk
    self%microchunks_pass_1 = [self%microchunks_pass_1, new_chunk]
  end subroutine append_microchunk_pass_1

  ! Appends a single chunk2D to the end of the pass-2 microchunk array.
  subroutine append_microchunk_pass_2( self, new_chunk )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D),    intent(in)    :: new_chunk
    self%microchunks_pass_2 = [self%microchunks_pass_2, new_chunk]
  end subroutine append_microchunk_pass_2

  ! Appends a single chunk2D to the end of the match microchunk array.
  subroutine append_microchunk_match( self, new_chunk )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D),    intent(in)    :: new_chunk
    self%microchunks_match = [self%microchunks_match, new_chunk]
  end subroutine append_microchunk_match

  ! --------------------------------------------------------------------------
  ! GENERATION
  ! --------------------------------------------------------------------------

  ! Partitions un-included records from project_list into pass-1 microchunks,
  ! each holding up to MICROCHUNK_P1_THRESHOLD particles. For each chunk:
  ! creates its subdirectory, accumulates micrographs until the threshold is
  ! reached, builds and writes a project file, applies any
  ! SIMPLE_CHUNK_PARTITION environment override, configures the command line,
  ! marks the consumed records as included in project_list, and updates the
  ! imported_projects.txt file table with all currently included project files.
  subroutine generate_microchunks_pass_1( self, project_list )
    class(microchunked2D), intent(inout) :: self
    type(rec_list),        intent(inout) :: project_list

    type(string),     allocatable   :: projfiles(:)
    logical,          allocatable   :: included(:)
    integer,          allocatable   :: ids(:), nptcls(:)
    type(chunk2D)                   :: new_chunk
    type(sp_project)                :: chunk_project
    type(rec_list)                  :: chunk_project_list
    type(string)                    :: chunk_folder
    character(len=STDLEN)           :: chunk_part_env
    integer(timer_int_kind)         :: t0
    integer                         :: imic, chunk_nptcls, envlen, chunk_id

    t0 = timer_start()
    do while( project_list%get_nptcls_tot(l_not_included=.true.) > MICROCHUNK_P1_THRESHOLD )
      included = project_list%get_included_flags()
      ids      = pack(project_list%get_ids(),    .not. included)
      nptcls   = pack(project_list%get_nptcls(), .not. included)
      if( size(ids) == 0 ) exit

      ! Accumulate micrographs until the particle threshold is exceeded
      chunk_nptcls = 0
      do imic = 1, size(nptcls)
        chunk_nptcls = chunk_nptcls + nptcls(imic)
        if( chunk_nptcls > MICROCHUNK_P1_THRESHOLD ) exit
      end do
      if( chunk_nptcls == 0 ) exit

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

      ! Update the imported-projects file table with all included project files
      included = project_list%get_included_flags()
      projfiles = project_list%get_projnames()
      projfiles = pack(projfiles, included)
      projfiles = remove_duplicates(projfiles)
      call write_filetable(string('imported_projects.txt'), projfiles)

      write(logfhandle,'(A,I6,A,I8,A)') &
        '>>> MICROCHUNK PASS 1 # ', chunk_id, ' GENERATED WITH ', chunk_nptcls, ' PARTICLES'
    end do
    call timer_stop(t0, string('generate_microchunks_pass_1'))
  end subroutine generate_microchunks_pass_1

  ! Merges rejection-complete pass-1 microchunks into pass-2 microchunks,
  ! each accumulating up to MICROCHUNK_P2_THRESHOLD selected particles. For
  ! each pass-2 chunk: creates its subdirectory, collects eligible pass-1
  ! project files and particle counts, merges the project files, clears stale
  ! 2D results and intermediate files, applies any SIMPLE_CHUNK_PARTITION
  ! environment override, configures the command line, and marks the consumed
  ! pass-1 chunks as complete (writing a COMPLETE sentinel file).
  subroutine generate_microchunks_pass_2( self )
    class(microchunked2D), intent(inout) :: self

    type(string),     allocatable :: projfiles(:)
    type(chunk2D)                 :: new_chunk
    type(sp_project)              :: chunk_project
    type(string)                  :: chunk_folder
    integer(timer_int_kind)       :: t0
    integer                       :: i, chunk_nptcls, chunk_nptcls_selected, chunk_id

    t0 = timer_start()
    do while( self%get_n_pass_1_non_rejected_ptcls() > MICROCHUNK_P2_THRESHOLD )
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
            call simple_touch(src%folder%to_char() // '/COMPLETE')
          end if
          if( chunk_nptcls_selected > MICROCHUNK_P2_THRESHOLD ) exit
        end associate
      end do
      if( size(projfiles) == 0 ) exit

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
    call timer_stop(t0, string('generate_microchunks_pass_2'))
  end subroutine generate_microchunks_pass_2

  ! Merges all rejection-complete pass-2 microchunks into the single reference
  ! chunk once their combined selected-particle count exceeds REFCHUNK_THRESHOLD.
  ! No-op if the reference chunk already exists (nptcls > 0) or the threshold
  ! is not yet met. Populates reference chunk fields, merges project files,
  ! clears stale 2D results and intermediate files, applies any
  ! SIMPLE_CHUNK_PARTITION environment override, and configures the command line.
  ! Note: pass-2 chunks consumed here are not marked complete at this stage;
  ! they are marked complete after match chunk generation.
  subroutine generate_refchunk( self )
    class(microchunked2D), intent(inout) :: self

    type(string),     allocatable :: projfiles(:)
    type(sp_project)              :: chunk_project
    integer(timer_int_kind)       :: t0
    integer                       :: i, chunk_nptcls, chunk_nptcls_selected

    t0 = timer_start()

    ! Only generate once, and only when enough selected particles are available
    if( self%refchunk%projfile%strlen() > 0 )                         return
    if( self%get_n_pass_2_non_rejected_ptcls() < REFCHUNK_THRESHOLD ) return

    allocate(projfiles(0))
    chunk_nptcls          = 0
    chunk_nptcls_selected = 0

    ! Consume all rejection-complete pass-2 chunks (marked complete after match)
    do i = 1, self%get_n_microchunks_pass_2()
      associate( src => self%microchunks_pass_2(i) )
        if( src%rejection_complete .and. .not. src%complete ) then
          chunk_nptcls_selected = chunk_nptcls_selected + src%nptcls_selected
          chunk_nptcls          = chunk_nptcls          + src%nptcls
          projfiles             = [projfiles, src%projfile]
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
    call timer_stop(t0, string('generate_refchunk'))
  end subroutine generate_refchunk

  ! Creates one match microchunk for each rejection-complete pass-2 microchunk
  ! that does not yet have a corresponding match chunk, provided the reference
  ! chunk rejection is complete (i.e. refs and box are available). Each match
  ! chunk reuses the pass-2 project file (copied into a new folder) and is
  ! configured for template-guided classification against the reference stack.
  ! Marks the source pass-2 chunk as complete and writes a COMPLETE sentinel.
  subroutine generate_microchunks_match( self )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D)                   :: new_chunk
    type(string)                    :: chunk_folder
    integer(timer_int_kind)         :: t0
    integer                         :: i, chunk_id

    t0 = timer_start()

    ! Match chunks require the reference stack — skip until refs is available
    if( self%box == 0 ) then
      call timer_stop(t0, string('generate_microchunks_match'))
      return
    end if

    do i = 1, self%get_n_microchunks_pass_2()
      associate( src => self%microchunks_pass_2(i) )
        if( .not. src%rejection_complete ) cycle
        if( src%complete )                 cycle
        chunk_id                  = src%id
        chunk_folder              = string(self%outdir_microchunks_match%to_char() // &
                                           '/microchunk_match_' // int2str(chunk_id))
        call simple_mkdir(chunk_folder)
        new_chunk%id              = chunk_id
        new_chunk%nptcls          = src%nptcls
        new_chunk%nptcls_selected = src%nptcls_selected
        new_chunk%folder          = simple_abspath(chunk_folder)
        new_chunk%projfile        = string(new_chunk%folder%to_char() // &
                                           '/microchunk_match_' // int2str(chunk_id) // METADATA_EXT)
        call self%generate_microchunk_match_cline(new_chunk)
        call simple_copy_file(src%projfile, new_chunk%projfile)
        src%complete = .true.
        call simple_touch(src%folder%to_char() // '/COMPLETE')
        call self%append_microchunk_match(new_chunk)
        write(logfhandle,'(A,I6,A,I8,A,I8,A)') '>>> MICROCHUNK MATCH # ', chunk_id, &
          ' GENERATED WITH ', src%nptcls_selected, '/', src%nptcls, ' PARTICLES'
      end associate
    end do
    call timer_stop(t0, string('generate_microchunks_match'))
  end subroutine generate_microchunks_match

  ! --------------------------------------------------------------------------
  ! COMMAND-LINE BUILDERS
  ! --------------------------------------------------------------------------

  ! Populates the abinitio2D command line for a pass-1 microchunk with:
  ! program name, project file and name, no-mkdir flag, thread count, mask
  ! diameter, class count, low-pass stop cutoff (DEFAULT_MICRO_P1_LP), and
  ! wall-time limit.
  subroutine generate_microchunk_pass_1_cline( self, new_chunk )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: new_chunk
    associate( cline => new_chunk%cline )
      call cline%set('prg',      'abinitio2D')
      call cline%set('projfile', new_chunk%projfile)
      call cline%set('projname', 'microchunk_pass_1')
      call cline%set('mkdir',    'no')
      call cline%set('nthr',     self%nthr)
      call cline%set('mskdiam',  self%mskdiam)
      call cline%set('ncls',     DEFAULT_NCLS)
      call cline%set('lpstart',  DEFAULT_LPSTART)
      call cline%set('lpstop',   DEFAULT_MICRO_P1_LP)
      call cline%set('walltime', DEFAULT_WALLTIME)
    end associate
  end subroutine generate_microchunk_pass_1_cline

  ! Populates the abinitio2D command line for a pass-2 microchunk with:
  ! program name, project file and name, no-mkdir flag, thread count, mask
  ! diameter, class count, low-pass stop cutoff (DEFAULT_MICRO_P2_LP), and
  ! wall-time limit.
  subroutine generate_microchunk_pass_2_cline( self, new_chunk )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: new_chunk
    associate( cline => new_chunk%cline )
      call cline%set('prg',      'abinitio2D')
      call cline%set('projfile', new_chunk%projfile)
      call cline%set('projname', 'microchunk_pass_2')
      call cline%set('mkdir',    'no')
      call cline%set('nthr',     self%nthr)
      call cline%set('mskdiam',  self%mskdiam)
      call cline%set('ncls',     DEFAULT_NCLS)
      call cline%set('lpstart',  DEFAULT_LPSTART)
      call cline%set('lpstop',   DEFAULT_MICRO_P2_LP)
      call cline%set('walltime', DEFAULT_WALLTIME)
    end associate
  end subroutine generate_microchunk_pass_2_cline

  ! Populates the abinitio2D command line for a match microchunk with:
  ! program name, project file and name, no-mkdir flag, thread count, mask
  ! diameter, class count, low-pass stop cutoff (DEFAULT_MICRO_P2_LP),
  ! reference stack path (refs), extremal limit, number of stages, and
  ! wall-time limit. The refs field enables template-guided classification
  ! against the reference chunk class averages.
  subroutine generate_microchunk_match_cline( self, new_chunk )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: new_chunk
    associate( cline => new_chunk%cline )
      call cline%set('prg',            'abinitio2D')
      call cline%set('projfile',       new_chunk%projfile)
      call cline%set('projname',       'microchunk_match')
      call cline%set('mkdir',          'no')
      call cline%set('nthr',           self%nthr)
      call cline%set('mskdiam',        self%mskdiam)
      call cline%set('ncls',           DEFAULT_NCLS)
      call cline%set('lpstart',        DEFAULT_LPSTART)
      call cline%set('lpstop',         DEFAULT_MICRO_P2_LP)
      call cline%set('refs',           self%refs)
      call cline%set('extr_lim',       14)
      call cline%set('eo_stage',       'no')
      call cline%set('nits_per_stage', 3)
      call cline%set('walltime',       DEFAULT_WALLTIME)
      call cline%set('refine',         'prob')
    end associate
  end subroutine generate_microchunk_match_cline

  ! Populates the abinitio2D command line for the reference chunk with:
  ! program name, project file and name, no-mkdir flag, thread count, mask
  ! diameter, class count, low-pass stop cutoff (DEFAULT_REF_LP), and
  ! wall-time limit.
  subroutine generate_refchunk_cline( self, new_chunk )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: new_chunk
    associate( cline => new_chunk%cline )
      call cline%set('prg',      'abinitio2D')
      call cline%set('projfile', new_chunk%projfile)
      call cline%set('projname', 'refchunk')
      call cline%set('mkdir',    'no')
      call cline%set('nthr',     self%nthr)
      call cline%set('mskdiam',  self%mskdiam)
      call cline%set('ncls',     DEFAULT_NCLS)
      call cline%set('lpstop',   DEFAULT_REF_LP)
      call cline%set('walltime', DEFAULT_WALLTIME)
      call cline%set('refine',   'prob')
    end associate
  end subroutine generate_refchunk_cline

  ! --------------------------------------------------------------------------
  ! COMBINATION
  ! --------------------------------------------------------------------------

  ! Merges the project files of all complete match chunks into a single combined
  ! project written to completedir/microchunks_match_combined.simple. Only the
  ! class assignment is retained in os_ptcl2D (shifts and other 2D parameters
  ! are stripped); stale os_cls2D, os_cls3D and os_out segments are killed.
  ! No-op if no match chunks are complete or if the combined file already exists.
  subroutine combine_completed_match_chunks( self, combined_projfile )
    class(microchunked2D), intent(inout) :: self
    type(string),          intent(in)    :: combined_projfile
    type(string),     allocatable :: projfiles(:)
    integer,          allocatable :: pops(:)
    type(sp_project)              :: combined_project, spproj_ref
    type(string)                  :: ref_stkname
    integer(timer_int_kind)       :: t0
    integer                       :: i, nptcls_tot, nptcls_sel_tot, n_complete, nrefs
    real                          :: smpd_refs

    if( self%get_n_microchunks_match() == 0 ) return
    if( file_exists(combined_projfile) )      return

    t0 = timer_start()

    ! Collect projfiles from all complete match chunks
    allocate(projfiles(self%get_n_microchunks_match()))
    nptcls_tot     = 0
    nptcls_sel_tot = 0
    n_complete     = 0
    do i = 1, self%get_n_microchunks_match()
      associate( chunk => self%microchunks_match(i) )
        if( chunk%complete ) then
          n_complete            = n_complete + 1
          projfiles(n_complete) = chunk%projfile
          nptcls_tot            = nptcls_tot     + chunk%nptcls
          nptcls_sel_tot        = nptcls_sel_tot + chunk%nptcls_selected
        end if
      end associate
    end do
    if( n_complete == 0 ) return
    projfiles = projfiles(:n_complete)

    call spproj_ref%read_segment('out',   self%refchunk%projfile)
    call spproj_ref%read_segment('cls2D', self%refchunk%projfile)
    call spproj_ref%get_cavgs_stk(ref_stkname, nrefs, smpd_refs)

    ! Merge, then strip all 2D parameters except class assignment
    call merge_chunk_projfiles(projfiles, simple_abspath(self%completedir), &
                               combined_project, write_proj=.false., update_classno=.false.)
    call combined_project%os_ptcl2D%delete_2Dclustering(keepshifts=.false., keepcls=.true.)
    call combined_project%os_cls2D%kill()
    call combined_project%os_cls3D%kill()
    call combined_project%os_out%kill()
    call combined_project%os_cls2D%copy(spproj_ref%os_cls2D)
    call combined_project%add_cavgs2os_out(ref_stkname, smpd_refs, 'cavg')
    call combined_project%os_ptcl2D%get_pops(pops, 'class')
    call combined_project%os_cls2D%set_all('pop', real(pops))
    call combined_project%map2ptcls_state()
    call combined_project%write(combined_projfile)
    call combined_project%kill()
    call spproj_ref%kill()
    deallocate(projfiles, pops)

    write(logfhandle,'(A,I6,A,I8,A,I8,A)') '>>> COMBINED ', n_complete, &
      ' MATCH CHUNKS INTO ', nptcls_sel_tot, '/', nptcls_tot, ' PARTICLES'
    call timer_stop(t0, string('combine_completed_match_chunks'))
  end subroutine combine_completed_match_chunks

  ! --------------------------------------------------------------------------
  ! SUBMISSION AND COLLECTION
  ! --------------------------------------------------------------------------

  ! Submits pending chunks to the queue asynchronously, up to the nparallel
  ! concurrency limit. Priority order: refchunk > match > pass-2 > pass-1.
  ! Skips chunks that are already running or complete. For each submitted
  ! chunk: changes into its working directory, dispatches the job, marks it
  ! as running, and releases its command line. Restores the original working
  ! directory on completion.
  subroutine submit( self )
    class(microchunked2D), intent(inout) :: self
    type(string)            :: cwd
    integer(timer_int_kind) :: t0
    integer                 :: i
    t0 = timer_start()
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

    ! Submit match microchunks second
    do i = 1, self%get_n_microchunks_match()
      if( self%get_n_chunks_running() >= self%nparallel ) exit
      associate( chunk => self%microchunks_match(i) )
        if( chunk%abinitio2D_running .or. chunk%abinitio2D_complete ) cycle
        call simple_chdir(chunk%folder)
        CWD_GLOB = chunk%folder%to_char()
        call self%qenv%exec_simple_prg_in_queue_async( &
          chunk%cline, string('./distr_microchunk'), string('simple_log_microchunk_match'))
        chunk%abinitio2D_running = .true.
        call chunk%cline%kill()
        write(logfhandle,'(A,I6)') '>>> INITIATED 2D ANALYSIS OF MICROCHUNK MATCH # ', chunk%id
      end associate
    end do

    ! Submit pass-2 microchunks third
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
    write(*,*) ' self%get_n_microchunks_pass_1()',  self%get_n_microchunks_pass_1()
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
    call timer_stop(t0, string('submit'))
  end subroutine submit

  ! Polls all running pass-1, pass-2, match, and reference chunks for the
  ! ABINITIO2D_FINISHED sentinel file. For each newly completed chunk:
  ! transitions it from running to complete and immediately runs class-average
  ! rejection. Capturing refs and box from the refchunk after its rejection
  ! enables subsequent match chunk generation. The refchunk COMPLETE sentinel
  ! is written after its rejection completes.
  subroutine collect_and_reject( self )
    class(microchunked2D), intent(inout) :: self
    logical,               allocatable   :: cls_msk(:)
    integer,               allocatable   :: states(:)
    type(sp_project)        :: spproj
    integer(timer_int_kind) :: t0
    integer                 :: i, j, non_zero, ncls
    real                    :: smpd
    t0 = timer_start()

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

    do i = 1, self%get_n_microchunks_match()
      associate( chunk => self%microchunks_match(i) )
        if( chunk%abinitio2D_running ) then
          if( file_exists(chunk%folder%to_char() // '/' // ABINITIO2D_FINISHED) ) then
            chunk%abinitio2D_running  = .false.
            chunk%abinitio2D_complete = .true.
            chunk%complete            = .true.
            call spproj%read_segment('out',    chunk%projfile)
            call spproj%read_segment('ptcl2D', chunk%projfile)
            call spproj%read_segment('cls2D',  chunk%projfile)
            non_zero = spproj%os_ptcl2D%count_state_gt_zero()
            self%n_accepted_ptcls = self%n_accepted_ptcls + non_zero
            self%n_rejected_ptcls = self%n_rejected_ptcls + spproj%os_ptcl2D%get_noris() - non_zero 
            call spproj%get_cavgs_stk(self%match_stk, ncls, smpd)
            self%match_jpeg = swap_suffix(self%match_stk, JPG_EXT, MRC_EXT)
            call spproj%cavgs2jpg(self%match_jpeg_inds, self%match_jpeg, self%match_jpeg_xtiles, self%match_jpeg_ytiles)
            self%match_jpeg_pops = spproj%os_cls2D%get_all_asint('pop')
            self%match_jpeg_res  = spproj%os_cls2D%get_all('res')
            allocate(cls_msk, source=self%match_jpeg_inds /= 0)
            self%match_jpeg_inds = pack(self%match_jpeg_inds, cls_msk)
            self%match_jpeg_pops = pack(self%match_jpeg_pops, cls_msk)
            self%match_jpeg_res  = pack(self%match_jpeg_res,  cls_msk)
            call spproj%kill()
            deallocate(cls_msk)
            call simple_copy_file(chunk%projfile, self%completedir // '/' // basename(chunk%projfile))
            call simple_touch(chunk%folder%to_char() // '/COMPLETE')
            write(logfhandle,'(A,I6)') '>>> COMPLETED 2D ANALYSIS OF MICROCHUNK MATCH # ', chunk%id
          end if
        end if
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
      if( self%refchunk%rejection_complete ) then
        call spproj%read_segment('cls2D', self%refchunk%projfile)
        call spproj%read_segment('out',   self%refchunk%projfile)
        self%refs_jpeg = swap_suffix(self%refs, JPG_EXT, MRC_EXT)
        call spproj%cavgs2jpg(self%refs_jpeg_inds, self%refs_jpeg, self%refs_jpeg_xtiles, self%refs_jpeg_ytiles)
        self%ref_selection  = self%refs_jpeg_inds
        self%refs_jpeg_pops = spproj%os_cls2D%get_all_asint('pop')
        self%refs_jpeg_res  = spproj%os_cls2D%get_all('res')
        allocate(cls_msk, source=self%refs_jpeg_inds /= 0)
        self%refs_jpeg_inds = pack(self%refs_jpeg_inds, cls_msk)
        self%refs_jpeg_pops = pack(self%refs_jpeg_pops, cls_msk)
        self%refs_jpeg_res  = pack(self%refs_jpeg_res,  cls_msk)
        deallocate(cls_msk)
        states = spproj%os_cls2D%get_all_asint('state')
        allocate(cls_msk, source=states /= 0)
        self%ref_selection = pack(self%ref_selection, cls_msk)
        call spproj%kill()
        deallocate(cls_msk, states)
        call simple_touch(self%refchunk%folder%to_char() // '/COMPLETE')
        self%refchunk%complete = .true.
      end if
    end if
    call timer_stop(t0, string('collect_and_reject'))
  end subroutine collect_and_reject

  ! Performs rejection on a single chunk. Reads the class-average stack,
  ! computes the mask radius from image dimensions and mskdiam, runs outlier
  ! then auto (pass-1) or basic (pass-2 / reference / match) rejection, writes
  ! three filtered stacks (outlier-rejected, strategy-rejected, selected),
  ! propagates rejection flags into project states via map2ptcls_state, records
  ! the final selected particle count, and writes the REJECTION_FINISHED
  ! sentinel. When called on the reference chunk, also captures the
  ! class-average stack path and box size into self%refs and self%box for use
  ! by match chunk generation.
  ! When DEBUG=.true., also writes a companion _deselected project containing
  ! only the rejected classes (state=1 for rejected, state=0 for selected) with
  ! particle states restored to their pre-rejection values, for inspection.
  ! No-op if the chunk is not yet complete or has already been rejected.
  ! Cleans up all allocations on exit.
  subroutine reject_cavgs( self, chunk, label )
    class(microchunked2D), intent(inout) :: self
    type(chunk2D),    intent(inout) :: chunk
    type(string),     intent(in)    :: label

    type(image), allocatable :: cavg_imgs(:)
    logical,     allocatable :: l_rejected(:)
    integer,     allocatable :: states(:), ptcl_states(:)
    type(cluster2D_rejector) :: rejector
    type(sp_project)         :: spproj
    type(string)             :: stkname
    integer(timer_int_kind)  :: t0
    integer                  :: ncls, ldim(3)
    real                     :: smpd, mskrad

    if( .not. chunk%abinitio2D_complete ) return
    if( chunk%rejection_complete )        return

    t0 = timer_start()

    ! Read project and class averages
    call spproj%read(chunk%projfile)
    call spproj%get_cavgs_stk(stkname, ncls, smpd)
    cavg_imgs = read_cavgs_into_imgarr(spproj)
    if( size(cavg_imgs) == 0 ) then
      write(logfhandle,'(A,A,A,I6)') '>>> WARNING: no class averages found for ', &
        label%to_char(), ' chunk #', chunk%id
      call spproj%kill()
      return
    end if
    smpd      = cavg_imgs(1)%get_smpd()
    ldim      = cavg_imgs(1)%get_ldim()
    mskrad    = min(real(ldim(1) / 2) - COSMSKHALFWIDTH - 1., 0.5 * self%mskdiam / smpd)
    allocate(l_rejected(ncls), source=.false.)

    ! ! Outlier rejection applied at every tier
    ! call reject_outliers(cavg_imgs, mskrad, l_rejected)
    ! call write_cavgs(stkname, string('_rejected_outliers.mrc'), selected=.false.)

    ! ! Pass-1 uses auto rejection; all other tiers use basic rejection
    ! if( label == string(LABEL_PASS_1) ) then
    !   call reject_auto(cavg_imgs, l_rejected)
    !   call write_cavgs(stkname, string('_rejected_auto.mrc'),  selected=.false.)
    ! else
    !   call reject_basic(spproj%os_cls2D, l_rejected)
    !   call write_cavgs(stkname, string('_rejected_basic.mrc'), selected=.false.)
    ! end if
    ! call write_cavgs(stkname, string('_selected.mrc'), selected=.true.)

    call rejector%new(cavg_imgs, self%mskdiam)
    call rejector%reject_pop(spproj%os_cls2D)
    call rejector%reject_res(spproj%os_cls2D)
    call rejector%reject_mask()
  !  call rejector%reject_brightness()! doesnt add anything
    call rejector%reject_local_variance()
    l_rejected = rejector%get_rejected()
    call write_cavgs(stkname, string('_rejected.mrc'), selected=.false.)
    call write_cavgs(stkname, string('_selected.mrc'),  selected=.true.)
    call rejector%kill()

    ! Capture refs and box from the refchunk for use in match chunk generation
    if( label == string(LABEL_REF) ) then
      self%refs = stkname
      self%box  = ldim(1)
    end if

    ! Propagate rejection flags into project states
    allocate(states(ncls))
    states = spproj%os_cls2D%get_all_asint('state')
    if( DEBUG ) ptcl_states = spproj%os_ptcl2D%get_all_asint('state')
    where( l_rejected ) states = 0
    call spproj%os_cls2D%set_all('state', states)
    call spproj%os_cls3D%set_all('state', states)
    call spproj%map2ptcls_state()
    call spproj%write()

    call simple_touch(chunk%folder%to_char() // '/REJECTION_FINISHED')
    chunk%nptcls_selected    = spproj%os_ptcl2D%count_state_gt_zero()
    chunk%rejection_complete = .true.
    write(logfhandle,'(A,A,A,I6,A,I8,A,I8,A)') '>>> COMPLETED REJECTION FOR ', &
      label%to_char(), ' # ', chunk%id, ' : ', &
      chunk%nptcls_selected, '/', chunk%nptcls, ' PARTICLES SELECTED'

    if( DEBUG ) then
      l_rejected = .not. l_rejected
      states     = 1
      where( l_rejected ) states = 0
      call spproj%os_cls2D%set_all('state', states)
      call spproj%os_cls3D%set_all('state', states)
      call spproj%os_ptcl2D%set_all('state', ptcl_states)
      call spproj%os_ptcl3D%set_all('state', ptcl_states)
      call spproj%map2ptcls_state()
      call spproj%write(swap_suffix(chunk%projfile, '_deselected'//METADATA_EXT, METADATA_EXT))
    end if

    ! Cleanup
    call dealloc_imgarr(cavg_imgs)
    call spproj%kill()
    deallocate(l_rejected, states)
    if( DEBUG ) deallocate(ptcl_states)
    call timer_stop(t0, string('reject_cavgs'))

  contains

    ! Writes a filtered subset of class averages to a new stack. If selected
    ! is true, writes only un-rejected classes; if false, writes only rejected
    ! ones. The output filename is derived by swapping the .mrc suffix of
    ! stkpath with the provided suffix string.
    subroutine write_cavgs( stkpath, suffix, selected )
      type(string), intent(in) :: stkpath, suffix
      logical,      intent(in) :: selected
      type(string)   :: out_stkname, out_jpgname
      integer        :: istk, icls
      out_stkname = swap_suffix(stkpath, suffix, string('.mrc'))
      out_jpgname = swap_suffix(out_stkname, JPG_EXT, MRC_EXT)
      istk = 0
      do icls = 1, ncls
        if( selected .eqv. (.not. l_rejected(icls)) ) then
          istk = istk + 1
          call cavg_imgs(icls)%write(out_stkname, istk)
        end if
      end do
      call mrc2jpeg_tiled(out_stkname, out_jpgname)
      write(logfhandle,'(A,A)') '>>> JPEG ', out_jpgname%to_char()
    end subroutine write_cavgs

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
    call del_file(string(outdir%to_char() // '/' // FRCS_FILE))
    call del_file(string(outdir%to_char() // '/cavgs.mrc'))
    call del_file(string(outdir%to_char() // '/cavgs_odd.mrc'))
    call del_file(string(outdir%to_char() // '/cavgs_even.mrc'))
    call del_file(string(outdir%to_char() // '/sigma2_combined.star'))
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
    write(logfhandle,'(A,A,A,F8.1)') 'simple_microchunked2D->', routine_name%to_char(), ' execution time:', elapsed
    call flush(logfhandle)
  end subroutine timer_stop

end module simple_microchunked2D