!@descr: GUI metadata type for a single particle entry (sprite-sheet position + stats).
!==============================================================================
! MODULE: simple_gui_metadata_ptcl
!
! PURPOSE:
!   Extends gui_metadata_base with fields that describe one particle entry
!   in a sprite-sheet JPEG, matching the per-tile layout used for the
!   randomly-sampled particle montage produced by gui_metadata_project%set:
!     path            — absolute path to the particle-sample JPEG
!     pathlp          — absolute path to the source particle stack MRC
!     idx             — particle index within the ptcl2D segment
!     spritex/spritey — position of this tile within the sheet (percentage, 0–100)
!     spriteh/spritew — total sprite-sheet height/width (pixels)
!     df              — defocus estimate (microns, optional)
!     box             — particle box size in pixels (optional)
!   Provides set/get for all fields and a jsonise override that emits all
!   mandatory fields plus the optional df and box when they have been set.
!
! DEPENDENCIES:
!   json_module, simple_defs, simple_string, simple_error,
!   simple_gui_metadata_base, simple_gui_metadata_types, simple_gui_metadata_cavg2D
!==============================================================================
module simple_gui_metadata_ptcl
use json_module,                only: json_core, json_value
use simple_defs,                only: LONGSTRLEN
use simple_error,               only: simple_exception
use simple_string,              only: string
use simple_gui_metadata_base,   only: gui_metadata_base
use simple_gui_metadata_types,  only: GUI_METADATA_PTCL_TYPE
use simple_gui_metadata_cavg2D, only: sprite_sheet_pos

implicit none

public :: gui_metadata_ptcl
private
#include "simple_local_flags.inc"

type, extends( gui_metadata_base ) :: gui_metadata_ptcl
  private
  character(len=LONGSTRLEN) :: path       = ''                  ! absolute path to the particle-sample JPEG
  character(len=LONGSTRLEN) :: pathlp     = ''                  ! absolute path to the source particle stack MRC
  integer                   :: idx        = 0                   ! particle index within the ptcl2D segment
  type(sprite_sheet_pos)    :: sprite     = sprite_sheet_pos()  ! tile position and sheet dimensions
  real                      :: df         = 0.0   ! defocus estimate (microns); valid only when l_df
  integer                   :: box        = 0     ! particle box size in pixels; valid only when l_box
  logical                   :: l_df       = .false.
  logical                   :: l_box      = .false.
  integer                   :: i          = 1     ! index of this entry within the current batch (IPC routing)
  integer                   :: i_max      = 1     ! total entries in the current batch (IPC routing)
contains
  procedure :: set
  procedure :: get
  procedure :: get_idx
  procedure :: get_i
  procedure :: get_i_max
  procedure :: jsonise => jsonise_override
end type gui_metadata_ptcl

contains

  !---------------- setters ----------------

  ! Set all particle fields and mark the object as assigned.
  ! df and box are optional; omitting them clears the corresponding flag.
  ! i and i_max are required IPC routing fields (batch index / batch size).
  subroutine set( self, path, pathlp, idx, sprite, i, i_max, df, box )
    class(gui_metadata_ptcl),   intent(inout) :: self
    type(string),               intent(in)    :: path, pathlp
    integer,                    intent(in)    :: idx
    type(sprite_sheet_pos),     intent(in)    :: sprite
    integer,                    intent(in)    :: i, i_max
    real,    optional,          intent(in)    :: df
    integer, optional,          intent(in)    :: box
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%path       = path%to_char()
    self%pathlp     = pathlp%to_char()
    self%idx        = idx
    self%sprite     = sprite
    self%i          = i
    self%i_max      = i_max
    self%l_df       = present(df)
    if( self%l_df ) self%df = df
    self%l_box      = present(box)
    if( self%l_box ) self%box = box
  end subroutine set

  !---------------- getters ----------------

  ! Return all fields; result is .true. if the object has been assigned.
  ! df and box are set only when the corresponding optional was supplied to set().
  function get( self, path, pathlp, idx, sprite, df, box ) result( l_assigned )
    class(gui_metadata_ptcl),   intent(in)  :: self
    type(string),               intent(out) :: path, pathlp
    integer,                    intent(out) :: idx
    type(sprite_sheet_pos),     intent(out) :: sprite
    real,    optional,          intent(out) :: df
    integer, optional,          intent(out) :: box
    logical                                 :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
    path       = trim(self%path)
    pathlp     = trim(self%pathlp)
    idx        = self%idx
    sprite     = self%sprite
    if( present(df) .and. self%l_df ) df = self%df
    if( present(box) .and. self%l_box ) box = self%box
  end function get

  ! Return the particle index (tile position within the sprite sheet).
  function get_idx( self ) result( idx )
    class(gui_metadata_ptcl), intent(in) :: self
    integer                              :: idx
    idx = self%idx
  end function get_idx

  ! Return the index of this entry within the current IPC batch.
  function get_i( self ) result( i )
    class(gui_metadata_ptcl), intent(in) :: self
    integer                              :: i
    i = self%i
  end function get_i

  ! Return the total number of entries in the current IPC batch.
  function get_i_max( self ) result( i_max )
    class(gui_metadata_ptcl), intent(in) :: self
    integer                              :: i_max
    i_max = self%i_max
  end function get_i_max

  !---------------- serialisation ----------------

  ! Emit all mandatory fields plus optional df/box as a JSON object.
  ! Returns a null pointer when the object has not been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_ptcl), intent(inout) :: self
    type(json_core)                         :: json
    type(json_value),         pointer       :: json_ptr
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( .not.self%l_assigned ) then
      nullify(json_ptr)
      return
    endif
    call json%create_object(json_ptr, '')
    call json%add(json_ptr, "path",    trim(self%path))
    call json%add(json_ptr, "pathlp",  trim(self%pathlp))
    call json%add(json_ptr, "spritex", dble(self%sprite%x))
    call json%add(json_ptr, "spritey", dble(self%sprite%y))
    call json%add(json_ptr, "spriteh", self%sprite%h)
    call json%add(json_ptr, "spritew", self%sprite%w)
    call json%add(json_ptr, "idx",     self%idx)
    if( self%l_df ) call json%add(json_ptr, "df", dble(self%df))
    if( self%l_box ) call json%add(json_ptr, "box", self%box)
  end function jsonise_override

end module simple_gui_metadata_ptcl
