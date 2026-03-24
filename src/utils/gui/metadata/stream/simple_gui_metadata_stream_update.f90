!@descr: GUI metadata for a stream quality update — CTF resolution, astigmatism, and ice score thresholds broadcast from the GUI
!==============================================================================
! MODULE: simple_gui_metadata_stream_update
!
! PURPOSE:
!   Extends gui_metadata_base with the three quality-threshold fields that the
!   GUI broadcasts to running stream processes when the user adjusts acceptance
!   cutoffs mid-session: CTF resolution limit (A), astigmatism limit (A), and
!   ice-contamination score limit.  Each field may be set independently when
!   only a subset of thresholds change.
!
! TYPES:
!   gui_metadata_stream_update — extends gui_metadata_base
!     set_ctfres_update()     — assign the CTF resolution threshold
!     set_astigmatism_update() — assign the astigmatism threshold
!     set_icescore_update()   — assign the ice score threshold
!     get_ctfres_update()     — retrieve the CTF resolution threshold
!     get_astigmatism_update() — retrieve the astigmatism threshold
!     get_icescore_update()   — retrieve the ice score threshold
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
  real :: ctfresupdate      = 0.0  ! CTF resolution threshold (A); 0 = unset
  real :: astigmatismupdate = 0.0  ! astigmatism threshold (A); 0 = unset
  real :: icescoreupdate    = 0.0  ! ice-contamination score threshold; 0 = unset
contains
  procedure :: set_ctfres_update
  procedure :: set_astigmatism_update
  procedure :: set_icescore_update
  procedure :: get_ctfres_update
  procedure :: get_astigmatism_update
  procedure :: get_icescore_update
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

  ! Assign the updated astigmatism threshold received from the GUI.
  subroutine set_astigmatism_update( self, astigmatismupdate )
    class(gui_metadata_stream_update), intent(inout) :: self
    real,                              intent(in)    :: astigmatismupdate
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned        = .true.
    self%astigmatismupdate = astigmatismupdate
  end subroutine set_astigmatism_update

  ! Assign the updated ice score threshold received from the GUI.
  subroutine set_icescore_update( self, icescoreupdate )
    class(gui_metadata_stream_update), intent(inout) :: self
    real,                              intent(in)    :: icescoreupdate
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned     = .true.
    self%icescoreupdate = icescoreupdate
  end subroutine set_icescore_update

  ! Retrieve the CTF resolution threshold.
  function get_ctfres_update( self ) result( ctfresupdate )
    class(gui_metadata_stream_update), intent(in) :: self
    real                                          :: ctfresupdate
    ctfresupdate = self%ctfresupdate
  end function get_ctfres_update

  ! Retrieve the astigmatism threshold.
  function get_astigmatism_update( self ) result( astigmatismupdate )
    class(gui_metadata_stream_update), intent(in) :: self
    real                                          :: astigmatismupdate
    astigmatismupdate = self%astigmatismupdate
  end function get_astigmatism_update

  ! Retrieve the ice score threshold.
  function get_icescore_update( self ) result( icescoreupdate )
    class(gui_metadata_stream_update), intent(in) :: self
    real                                          :: icescoreupdate
    icescoreupdate = self%icescoreupdate
  end function get_icescore_update

end module simple_gui_metadata_stream_update