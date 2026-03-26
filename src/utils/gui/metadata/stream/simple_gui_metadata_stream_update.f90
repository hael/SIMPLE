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
!   increase_nmics         — additional micrographs to use before re-picking (0 = none)
!   ctfresupdate           — CTF resolution threshold (A); 0 = unset
!   astigmatismupdate      — astigmatism threshold (A);   0 = unset
!   icescoreupdate         — ice-contamination score;      0 = unset
!
! PROCEDURES:
!   set/get_ctfres_update, set/get_astigmatism_update, set/get_icescore_update
!   set/get_increase_nmics
!   set/get_pickrefs_selection, set/get_pickrefs_selection_length
!
! DEPENDENCIES:
!   json_kinds, json_module, simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_stream_update
use json_kinds
use json_module,              only: json_core, json_value
use simple_error,             only: simple_exception
use simple_gui_metadata_base, only: gui_metadata_base

implicit none

public :: gui_metadata_stream_update
private
#include "simple_local_flags.inc"

type, extends(gui_metadata_base) :: gui_metadata_stream_update
  private
  integer :: pickrefs_selection(2000)        = 0       ! per-class selection; 1 = selected, 0 = not selected
  integer :: pickrefs_selection_length       = 0       ! number of classes in the selection
  integer :: increase_nmics               = 0       ! number of additional micrographs requested before re-picking; 0 = no request
  real    :: ctfresupdate                 = 0.0     ! CTF resolution threshold (A); 0 = unset
  real    :: astigmatismupdate            = 0.0     ! astigmatism threshold (A);   0 = unset
  real    :: icescoreupdate               = 0.0     ! ice-contamination score;      0 = unset
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

end module simple_gui_metadata_stream_update