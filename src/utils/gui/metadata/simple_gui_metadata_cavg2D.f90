!@descr: GUI metadata type for a single 2D class average entry (sprite-sheet position + stats).
!==============================================================================
! MODULE: simple_gui_metadata_cavg2D
!
! PURPOSE:
!   Extends gui_metadata_base with fields that describe one class-average entry
!   in a sprite-sheet JPEG, matching the JSON produced by add_cls2D_accepted_to_json
!   in simple_stream_p05_sieve_cavgs:
!     path            — absolute path to the selection JPEG
!     idx             — class index within the JPEG map
!     spritex/spritey — position of this tile within the sheet (percentage, 0–100)
!     spriteh/spritew — total sprite-sheet height/width (pixels)
!     res             — resolution estimate (Angstroms, optional)
!     pop             — particle population count (optional)
!   Provides set/get for all fields and a jsonise override that emits all
!   mandatory fields plus the optional res and pop when they have been set.
!
! DEPENDENCIES:
!   json_module, simple_defs, simple_string, simple_error,
!   simple_gui_metadata_base, simple_gui_metadata_types
!==============================================================================
module simple_gui_metadata_cavg2D
use json_module,               only: json_core, json_value
use simple_defs,               only: LONGSTRLEN
use simple_error,              only: simple_exception
use simple_string,             only: string
use simple_gui_metadata_base,  only: gui_metadata_base
use simple_gui_metadata_types, only: GUI_METADATA_CAVG2D_TYPE

implicit none

public :: gui_metadata_cavg2D, sprite_sheet_pos
private
#include "simple_local_flags.inc"

! Position and dimensions of one tile within a sprite-sheet JPEG.
! x/y are the tile origin as a percentage of the sheet (0–100);
! h/w are the total sheet dimensions in pixels.
type :: sprite_sheet_pos
  real    :: x = 0.0   ! tile x-origin (%, 0–100)
  real    :: y = 0.0   ! tile y-origin (%, 0–100)
  integer :: h = 0     ! total sprite-sheet height (pixels)
  integer :: w = 0     ! total sprite-sheet width (pixels)
end type sprite_sheet_pos

type, extends( gui_metadata_base ) :: gui_metadata_cavg2D
  private
  character(len=LONGSTRLEN) :: path    = ''                  ! absolute path to the cavgs JPEG
  character(len=LONGSTRLEN) :: mrcpath = ''               ! absolute path to the cavgs MRC
  integer                   :: idx     = 0                   ! class index within the MRC
  type(sprite_sheet_pos)    :: sprite  = sprite_sheet_pos()  ! tile position and sheet dimensions
  real                      :: res     = 0.0   ! resolution estimate (Angstroms); valid only when l_res
  integer                   :: pop     = 0     ! particle population count; valid only when l_pop
  logical                   :: l_res   = .false.
  logical                   :: l_pop   = .false.
  integer                   :: i       = 1     ! index of this entry within the current batch (IPC routing)
  integer                   :: i_max   = 1     ! total entries in the current batch (IPC routing)
contains
  procedure :: set
  procedure :: get
  procedure :: get_idx
  procedure :: get_i
  procedure :: get_i_max
  procedure :: jsonise => jsonise_override
end type gui_metadata_cavg2D

contains

  !---------------- setters ----------------

  ! Set all class-average fields and mark the object as assigned.
  ! res and pop are optional; omitting them clears the corresponding flag.
  ! i and i_max are required IPC routing fields (batch index / batch size).
  subroutine set( self, path, mrcpath, idx, sprite, i, i_max, res, pop )
    class(gui_metadata_cavg2D), intent(inout) :: self
    type(string),               intent(in)    :: path, mrcpath
    integer,                    intent(in)    :: idx
    type(sprite_sheet_pos),     intent(in)    :: sprite
    integer,                    intent(in)    :: i, i_max
    real,    optional,          intent(in)    :: res
    integer, optional,          intent(in)    :: pop
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%path       = path%to_char()
    self%mrcpath    = mrcpath%to_char()
    self%idx        = idx
    self%sprite     = sprite
    self%i          = i
    self%i_max      = i_max
    self%l_res      = present(res)
    if( self%l_res ) self%res = res
    self%l_pop      = present(pop)
    if( self%l_pop ) self%pop = pop
  end subroutine set

  !---------------- getters ----------------

  ! Return all fields; result is .true. if the object has been assigned.
  ! res and pop are set only when the corresponding optional was supplied to set().
  function get( self, path, mrcpath, idx, sprite, res, pop ) result( l_assigned )
    class(gui_metadata_cavg2D), intent(in)  :: self
    type(string),               intent(out) :: path, mrcpath
    integer,                    intent(out) :: idx
    type(sprite_sheet_pos),     intent(out) :: sprite
    real,    optional,          intent(out) :: res
    integer, optional,          intent(out) :: pop
    logical                                 :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
    path       = trim(self%path)
    mrcpath    = trim(self%mrcpath)
    idx        = self%idx
    sprite     = self%sprite
    if( present(res) .and. self%l_res ) res = self%res
    if( present(pop) .and. self%l_pop ) pop = self%pop
  end function get

  ! Return the class index (tile position within the sprite sheet).
  function get_idx( self ) result( idx )
    class(gui_metadata_cavg2D), intent(in) :: self
    integer                                :: idx
    idx = self%idx
  end function get_idx

  ! Return the index of this entry within the current IPC batch.
  function get_i( self ) result( i )
    class(gui_metadata_cavg2D), intent(in) :: self
    integer                                :: i
    i = self%i
  end function get_i

  ! Return the total number of entries in the current IPC batch.
  function get_i_max( self ) result( i_max )
    class(gui_metadata_cavg2D), intent(in) :: self
    integer                                :: i_max
    i_max = self%i_max
  end function get_i_max

  !---------------- serialisation ----------------

  ! Emit all mandatory fields plus optional res/pop as a JSON object.
  ! Returns a null pointer when the object has not been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_cavg2D), intent(inout) :: self
    type(json_core)                           :: json
    type(json_value),           pointer       :: json_ptr
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( .not.self%l_assigned ) then
      nullify(json_ptr)
      return
    endif
    call json%create_object(json_ptr, '')
    call json%add(json_ptr, "path",    trim(self%path))
    call json%add(json_ptr, "mrcpath", trim(self%mrcpath))
    call json%add(json_ptr, "spritex", dble(self%sprite%x))
    call json%add(json_ptr, "spritey", dble(self%sprite%y))
    call json%add(json_ptr, "spriteh", self%sprite%h)
    call json%add(json_ptr, "spritew", self%sprite%w)
    call json%add(json_ptr, "idx",     self%idx)
    if( self%l_res ) call json%add(json_ptr, "res", dble(self%res))
    if( self%l_pop ) call json%add(json_ptr, "pop", self%pop)
  end function jsonise_override

end module simple_gui_metadata_cavg2D
