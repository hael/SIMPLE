!@descr: GUI metadata for a stream quality update — thresholds and user selections broadcast from the GUI
! The master copies GUI response fields here and sends the whole object to the running p01, p03,
! p06 and p07. Unset values are 0; receivers ignore 0 and unchanged values. Readers: p01 thresholds,
! p03 pickrefs_selection+cycle, p06 mskdiam2D/snapshot2D, p07 snapshot3D.
module simple_gui_metadata_stream_update
use simple_error,             only: simple_exception
use simple_string,            only: string
use simple_gui_metadata_base, only: gui_metadata_base
use simple_gui_metadata_stream_solve3D_multistate, only: MAX_STATES_SOLVE3D_MULTISTATE

implicit none

public :: gui_metadata_stream_update
public :: MAX_PICKREFS_SELECTION, MAX_SNAPSHOT2D_SELECTION, MAX_SNAPSHOT_FNAME_LEN
private
#include "simple_local_flags.inc"

integer, parameter :: MAX_PICKREFS_SELECTION   = 500  ! classes a picking-reference selection holds
integer, parameter :: MAX_SNAPSHOT2D_SELECTION = 1000 ! classes a 2D snapshot selection holds
integer, parameter :: MAX_SNAPSHOT_FNAME_LEN   = 128  ! characters of a snapshot's file name (2D and 3D)

type, extends(gui_metadata_base) :: gui_metadata_stream_update
  private
  integer(kind=2)       :: pickrefs_selection(MAX_PICKREFS_SELECTION)     = 0 ! selected class indices
  integer               :: pickrefs_cycle               = 0    ! current pickrefs cycle
  integer               :: pickrefs_selection_length    = 0    ! number of classes in the selection
  real                  :: ctfresupdate                 = 0.0  ! CTF resolution threshold (A); 0 = unset
  real                  :: astigmatismupdate            = 0.0  ! astigmatism threshold (A);   0 = unset
  real                  :: icescoreupdate               = 0.0  ! ice-contamination score;      0 = unset
  real                  :: mskdiam2D                    = 0.0  ! mask diameter for 2D classification (A); 0 = unset
  integer               :: snapshot2D_id                = 0    ! snapshot set ID; 0 = unset
  integer               :: snapshot2D_iteration         = 0    ! 2D classification iteration to snapshot; 0 = unset
  integer(kind=2)       :: snapshot2D_selection(MAX_SNAPSHOT2D_SELECTION) = 0 ! class indices included in the snapshot
  integer               :: snapshot2D_selection_length  = 0    ! number of valid entries in snapshot2D_selection
  character(len=MAX_SNAPSHOT_FNAME_LEN) :: snapshot2D_filename = '' ! project file name for the snapshot
  integer               :: snapshot3D_id                = 0    ! 3D snapshot set ID; 0 = unset
  integer               :: snapshot3D_selection(MAX_STATES_SOLVE3D_MULTISTATE) = 0 ! states included in the 3D snapshot
  integer               :: snapshot3D_selection_length  = 0    ! number of valid entries in snapshot3D_selection
  character(len=MAX_SNAPSHOT_FNAME_LEN) :: snapshot3D_filename = '' ! project file name for the 3D snapshot
contains
  procedure :: kill => kill_override
  procedure :: set_ctfres_update
  procedure :: get_ctfres_update
  procedure :: set_astigmatism_update
  procedure :: get_astigmatism_update
  procedure :: set_icescore_update
  procedure :: get_icescore_update
  procedure :: set_pickrefs_selection
  procedure :: get_pickrefs_selection
  procedure :: set_pickrefs_cycle
  procedure :: get_pickrefs_cycle
  procedure :: get_pickrefs_selection_length
  procedure :: set_mskdiam2D_update
  procedure :: get_mskdiam2D_update
  procedure :: set_snapshot2D_update
  procedure :: get_snapshot2D_update
  procedure :: has_snapshot2D_update
  procedure :: set_snapshot3D_update
  procedure :: get_snapshot3D_update
  procedure :: has_snapshot3D_update
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

  ! Store the user's class selection as an integer array
  subroutine set_pickrefs_selection( self, selection )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: selection(:)
    integer :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    n = size(selection)
    if( n > size(self%pickrefs_selection) ) THROW_HARD('pickrefs_selection exceeds maximum size')
    self%l_assigned             = .true.
    self%pickrefs_selection_length = n
    self%pickrefs_selection(1:n)   = int(selection, kind=kind(self%pickrefs_selection))  ! only 1:n is ever read back
  end subroutine set_pickrefs_selection

  ! Retrieve the user's class selection as an integer array
  function get_pickrefs_selection( self ) result( selection )
    class(gui_metadata_stream_update), intent(in)  :: self
    integer, allocatable                           :: selection(:)
    integer :: n
    n = self%pickrefs_selection_length
    if( n < 0 ) THROW_HARD('pickrefs_selection_length is negative')
    if( n > size(self%pickrefs_selection) ) THROW_HARD('pickrefs_selection_length exceeds maximum size')
    allocate(selection(n))
    selection = self%pickrefs_selection(1:n)
  end function get_pickrefs_selection

  ! Set the p03 cycle whose class averages pickrefs_selection indexes.
  subroutine set_pickrefs_cycle( self, ncycle )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: ncycle
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned       = .true.
    self%pickrefs_cycle   = ncycle
  end subroutine set_pickrefs_cycle

  ! Retrieve the p03 cycle whose class averages pickrefs_selection indexes.
  function get_pickrefs_cycle( self ) result( ncycle )
    class(gui_metadata_stream_update), intent(in)  :: self
    integer                                        :: ncycle
    ncycle = self%pickrefs_cycle
  end function get_pickrefs_cycle

  ! Retrieve the number of classes in the selection.
  function get_pickrefs_selection_length( self ) result( n )
    class(gui_metadata_stream_update), intent(in) :: self
    integer                                       :: n
    n = self%pickrefs_selection_length
  end function get_pickrefs_selection_length

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
    if( filename%strlen() > MAX_SNAPSHOT_FNAME_LEN ) THROW_HARD('snapshot2D_filename exceeds maximum length')
    self%l_assigned                   = .true.
    self%snapshot2D_id                = snapshot_id
    self%snapshot2D_iteration         = iteration
    self%snapshot2D_selection_length  = n
    self%snapshot2D_selection(1:n)    = int(selection, kind=kind(self%snapshot2D_selection))
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
    if( n < 0 ) THROW_HARD('snapshot2D_selection_length is negative')
    if( n > size(self%snapshot2D_selection) ) THROW_HARD('snapshot2D_selection_length exceeds maximum size')
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

  ! Store a 3D snapshot request from the GUI: the states whose particles it holds, merged into one.
  subroutine set_snapshot3D_update( self, snapshot_id, selection, filename )
    class(gui_metadata_stream_update), intent(inout) :: self
    integer,                           intent(in)    :: snapshot_id
    integer,                           intent(in)    :: selection(:)
    type(string),                      intent(in)    :: filename
    integer :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    n = size(selection)
    if( n > size(self%snapshot3D_selection) ) THROW_HARD('snapshot3D_selection exceeds maximum size')
    if( filename%strlen() > MAX_SNAPSHOT_FNAME_LEN ) THROW_HARD('snapshot3D_filename exceeds maximum length')
    self%l_assigned                  = .true.
    self%snapshot3D_id               = snapshot_id
    self%snapshot3D_selection_length = n
    self%snapshot3D_selection(1:n)   = selection
    self%snapshot3D_filename         = filename%to_char()
  end subroutine set_snapshot3D_update

  ! Retrieve the snapshot3D request fields.
  subroutine get_snapshot3D_update( self, snapshot_id, selection, filename )
    class(gui_metadata_stream_update), intent(in)  :: self
    integer,                           intent(out) :: snapshot_id
    integer,           allocatable,    intent(out) :: selection(:)
    type(string),                      intent(out) :: filename
    integer :: n
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    n = self%snapshot3D_selection_length
    if( n < 0 ) THROW_HARD('snapshot3D_selection_length is negative')
    if( n > size(self%snapshot3D_selection) ) THROW_HARD('snapshot3D_selection_length exceeds maximum size')
    snapshot_id = self%snapshot3D_id
    selection   = self%snapshot3D_selection(1:n)
    filename    = trim(self%snapshot3D_filename)
  end subroutine get_snapshot3D_update

  ! Returns .true. when a snapshot3D request is present (snapshot_id > 0).
  function has_snapshot3D_update( self ) result( l_has )
    class(gui_metadata_stream_update), intent(in) :: self
    logical :: l_has
    l_has = self%snapshot3D_id > 0
  end function has_snapshot3D_update

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_stream_update), intent(inout) :: self
    select type( self )
      type is( gui_metadata_stream_update )
        self = gui_metadata_stream_update()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_stream_update
