!@descr: GUI metadata for a stream quality update — thresholds and user selections broadcast from the GUI
!==============================================================================
! MODULE: simple_gui_metadata_stream_update
!
! PURPOSE:
!   Extends gui_metadata_base with the fields the GUI broadcasts to running
!   stream processes when the user adjusts settings mid-session.  Each field
!   may be set independently; only fields that changed need to be sent.
!
! FIELDS:
!   pickrefs_selection        — per-class selection array (1 = selected, 0 = not)
!   pickrefs_selection_length — number of valid entries in pickrefs_selection
!   sieverefs_selection        — per-class selection for sieve refchunk (1 = selected, 0 = not)
!   sieverefs_selection_length — number of valid entries in sieverefs_selection
!   increase_nmics            — additional micrographs to use before re-picking (0 = none)
!   ctfresupdate              — CTF resolution threshold (A); 0 = unset
!   astigmatismupdate         — astigmatism threshold (A);   0 = unset
!   icescoreupdate            — ice-contamination score;      0 = unset
!   mskdiam2D                    — mask diameter for 2D classification (pixels); 0 = unset
!   snapshot2D_id                — snapshot set ID; 0 = unset
!   snapshot2D_iteration         — 2D classification iteration to snapshot; 0 = unset
!   snapshot2D_selection         — class indices included in the snapshot
!   snapshot2D_selection_length  — number of valid entries in snapshot2D_selection
!   snapshot2D_filename          — project file name for the snapshot
!
! PROCEDURES:
!   set/get_ctfres_update, set/get_astigmatism_update, set/get_icescore_update
!   set/get_increase_nmics
!   set/get_pickrefs_selection, set/get_pickrefs_selection_length
!   set/get_sieverefs_selection, set/get_sieverefs_selection_length
!   set/get_mskdiam2D_update
!   set/get_snapshot2D_update, has_snapshot2D_update
!
! DEPENDENCIES:
!   json_kinds, json_module, simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_stream_update
use json_kinds
use json_module,              only: json_core, json_value
use simple_defs,              only: STDLEN
use simple_error,             only: simple_exception
use simple_string,            only: string
use simple_gui_metadata_base, only: gui_metadata_base

implicit none

public :: gui_metadata_stream_update
private
#include "simple_local_flags.inc"

type, extends(gui_metadata_base) :: gui_metadata_stream_update
  private
  integer(kind=2)       :: pickrefs_selection(1000)     = 0    ! per-class selection; 1 = selected, 0 = not selected
  integer               :: pickrefs_selection_length    = 0    ! number of classes in the selection
  integer(kind=2)       :: sieverefs_selection(1000)    = 0    ! per-class selection for sieve refchunk; 1 = selected, 0 = not selected
  integer               :: sieverefs_selection_length   = 0    ! number of sieve-ref classes in the selection
  integer               :: increase_nmics               = 0    ! additional micrographs requested before re-picking; 0 = no request
  real                  :: ctfresupdate                 = 0.0  ! CTF resolution threshold (A); 0 = unset
  real                  :: astigmatismupdate            = 0.0  ! astigmatism threshold (A);   0 = unset
  real                  :: icescoreupdate               = 0.0  ! ice-contamination score;      0 = unset
  real                  :: mskdiam2D                    = 0.0  ! mask diameter for 2D classification (pixels); 0 = unset
  integer               :: snapshot2D_id                = 0    ! snapshot set ID; 0 = unset
  integer               :: snapshot2D_iteration         = 0    ! 2D classification iteration to snapshot; 0 = unset
  integer(kind=2)       :: snapshot2D_selection(1000)   = 0    ! class indices included in the snapshot
  integer               :: snapshot2D_selection_length  = 0    ! number of valid entries in snapshot2D_selection
  character(len=STDLEN) :: snapshot2D_filename          = ''   ! project file name for the snapshot
contains
  procedure :: set_ctfres_update
  procedure :: get_ctfres_update
  procedure :: set_astigmatism_update
  procedure :: get_astigmatism_update
  procedure :: set_icescore_update
  procedure :: get_icescore_update
  procedure :: set_increase_nmics
  procedure :: get_increase_nmics
  procedure :: set_pickrefs_selection
  procedure :: get_pickrefs_selection
  procedure :: set_pickrefs_selection_length
  procedure :: get_pickrefs_selection_length
  procedure :: set_sieverefs_selection
  procedure :: get_sieverefs_selection
  procedure :: set_sieverefs_selection_length
  procedure :: get_sieverefs_selection_length
  procedure :: set_mskdiam2D_update
  procedure :: get_mskdiam2D_update
  procedure :: set_snapshot2D_update
  procedure :: get_snapshot2D_update
  procedure :: has_snapshot2D_update
end type gui_metadata_stream_update

contains

  ! Assign the updated CTF resolution threshold received from the GUI.
  subroutine set_ctfres_update( self, ctfresupdate )
    class(gui_metadata_stream_update), intent(inout) :: self
    real,                              intent(in)    :: ctfresupdate
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned   = .true.
    self%ctfresupdate = ctfresupdate
  end subroutine set_ctfres_update

  ! Retrieve the CTF resolution threshold.
  function get_ctfres_update( self ) result( ctfresupdate )
    class(gui_metadata_stream_update), intent(in) :: self
    real                                          :: ctfresupdate
    ctfresupdate = self%ctfresupdate
  end function get_ctfres_update

  ! Assign the updated astigmatism threshold received from the GUI.
  subroutine set_astigmatism_update( self, astigmatismupdate )
    class(gui_metadata_stream_update), intent(inout) :: self
    real,                              intent(in)    :: astigmatismupdate
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned        = .true.
    self%astigmatismupdate = astigmatismupdate
  end subroutine set_astigmatism_update

  ! Retrieve the astigmatism threshold.
  function get_astigmatism_update( self ) result( astigmatismupdate )
    class(gui_metadata_stream_update), intent(in) :: self
    real                                          :: astigmatismupdate
    astigmatismupdate = self%astigmatismupdate
  end function get_astigmatism_update

  ! Assign the updated ice score threshold received from the GUI.
  subroutine set_icescore_update( self, icescoreupdate )
    class(gui_metadata_stream_update), intent(inout) :: self
    real,                              intent(in)    :: icescoreupdate
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned     = .true.
    self%icescoreupdate = icescoreupdate
  end subroutine set_icescore_update

  ! Retrieve the ice score threshold.
  function get_icescore_update( self ) result( icescoreupdate )
    class(gui_metadata_stream_update), intent(in) :: self
    real                                          :: icescoreupdate
    icescoreupdate = self%icescoreupdate
  end function get_icescore_update

  ! Set the number of additional micrographs requested before re-picking.
  subroutine set_increase_nmics( self, increase_nmics )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: increase_nmics
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned     = .true.
    self%increase_nmics = increase_nmics
  end subroutine set_increase_nmics

  ! Retrieve the number of additional micrographs requested before re-picking.
  function get_increase_nmics( self ) result( increase_nmics )
    class(gui_metadata_stream_update), intent(in) :: self
    integer                                       :: increase_nmics
    increase_nmics = self%increase_nmics
  end function get_increase_nmics

  ! Store the user's class selection as an integer array (1 = selected, 0 = not selected).
  subroutine set_pickrefs_selection( self, selection )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: selection(:)
    integer :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    n = size(selection)
    if( n > size(self%pickrefs_selection) ) THROW_HARD('pickrefs_selection exceeds maximum size')
    self%l_assigned             = .true.
    self%pickrefs_selection_length = n
    self%pickrefs_selection(1:n)   = selection  ! only 1:n is ever read back
  end subroutine set_pickrefs_selection

  ! Retrieve the class selection as an integer array (1 = selected, 0 = not selected).
  function get_pickrefs_selection( self ) result( selection )
    class(gui_metadata_stream_update), intent(in)  :: self
    integer, allocatable                           :: selection(:)
    integer :: n
    n = self%pickrefs_selection_length
    allocate(selection(n))
    selection = self%pickrefs_selection(1:n)
  end function get_pickrefs_selection

  ! Set the number of classes in the selection.
  subroutine set_pickrefs_selection_length( self, n )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned             = .true.
    self%pickrefs_selection_length = n
  end subroutine set_pickrefs_selection_length

  ! Retrieve the number of classes in the selection.
  function get_pickrefs_selection_length( self ) result( n )
    class(gui_metadata_stream_update), intent(in) :: self
    integer                                       :: n
    n = self%pickrefs_selection_length
  end function get_pickrefs_selection_length

  ! Store the user's sieve-reference class selection as an integer array
  ! Called when the GUI returns a refs_selection array from the particle-sieving stage.
  subroutine set_sieverefs_selection( self, selection )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: selection(:)
    integer :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    n = size(selection)
    if( n > size(self%sieverefs_selection) ) THROW_HARD('sieverefs_selection exceeds maximum size')
    self%l_assigned              = .true.
    self%sieverefs_selection_length = n
    self%sieverefs_selection(1:n)   = selection  ! only 1:n is ever read back
  end subroutine set_sieverefs_selection

  ! Retrieve the sieve-reference class selection as an integer array
  function get_sieverefs_selection( self ) result( selection )
    class(gui_metadata_stream_update), intent(in)  :: self
    integer, allocatable                           :: selection(:)
    integer :: n
    n = self%sieverefs_selection_length
    allocate(selection(n))
    selection = self%sieverefs_selection(1:n)
  end function get_sieverefs_selection

  ! Set the number of sieve-reference classes in the selection.
  subroutine set_sieverefs_selection_length( self, n )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned              = .true.
    self%sieverefs_selection_length = n
  end subroutine set_sieverefs_selection_length

  ! Retrieve the number of sieve-reference classes in the selection.
  function get_sieverefs_selection_length( self ) result( n )
    class(gui_metadata_stream_update), intent(in) :: self
    integer                                       :: n
    n = self%sieverefs_selection_length
  end function get_sieverefs_selection_length

  ! Assign the mask diameter for 2D classification received from the GUI.
  subroutine set_mskdiam2D_update( self, mskdiam2D )
    class(gui_metadata_stream_update), intent(inout) :: self
    real,                              intent(in)    :: mskdiam2D
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%mskdiam2D  = mskdiam2D
  end subroutine set_mskdiam2D_update

  ! Retrieve the mask diameter for 2D classification.
  function get_mskdiam2D_update( self ) result( mskdiam2D )
    class(gui_metadata_stream_update), intent(in) :: self
    real                                          :: mskdiam2D
    mskdiam2D = self%mskdiam2D
  end function get_mskdiam2D_update

  ! Store a 2D-classification snapshot request from the GUI.
  subroutine set_snapshot2D_update( self, snapshot_id, iteration, selection, filename )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: snapshot_id, iteration
    integer,                           intent(in)    :: selection(:)
    type(string),                      intent(in)    :: filename
    integer :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    n = size(selection)
    if( n > size(self%snapshot2D_selection) ) THROW_HARD('snapshot2D_selection exceeds maximum size')
    self%l_assigned                   = .true.
    self%snapshot2D_id                = snapshot_id
    self%snapshot2D_iteration         = iteration
    self%snapshot2D_selection_length  = n
    self%snapshot2D_selection(1:n)    = selection
    self%snapshot2D_filename          = filename%to_char()
  end subroutine set_snapshot2D_update

  ! Retrieve the snapshot2D request fields.
  subroutine get_snapshot2D_update( self, snapshot_id, iteration, selection, filename )
    class(gui_metadata_stream_update), intent(in)  :: self
    integer,                           intent(out) :: snapshot_id, iteration
    integer,           allocatable,    intent(out) :: selection(:)
    type(string),                      intent(out) :: filename
    integer :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    snapshot_id = self%snapshot2D_id
    iteration   = self%snapshot2D_iteration
    n           = self%snapshot2D_selection_length
    allocate(selection(n))
    selection   = self%snapshot2D_selection(1:n)
    filename    = trim(self%snapshot2D_filename)
  end subroutine get_snapshot2D_update

  ! Returns .true. when a snapshot2D request is present (snapshot_id > 0).
  function has_snapshot2D_update( self ) result( l_has )
    class(gui_metadata_stream_update), intent(in) :: self
    logical :: l_has
    l_has = self%snapshot2D_id > 0
  end function has_snapshot2D_update

end module simple_gui_metadata_stream_update