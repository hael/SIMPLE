!@descr: GUI metadata for the stream picking stages — micrograph counts, particle yield, and import timestamp
! Serves both initial picking (p03) and reference picking (p04). set() derives micrographs_rejected
! and particles_per_mic, and stamps last_micrograph_imported (Unix time) when the import count changes.
module simple_gui_metadata_stream_picking
  use unix,                     only: c_long, c_time
  use json_kinds
  use json_module,              only: json_core, json_value
  use simple_defs,              only: STDLEN
  use simple_string,            only: string
  use simple_error,             only: simple_exception
  use simple_gui_metadata_base, only: gui_metadata_base

  implicit none

  public :: gui_metadata_stream_picking
  private
#include "simple_local_flags.inc"

  type, extends(gui_metadata_base) :: gui_metadata_stream_picking
    private
    character(len=STDLEN) :: stage                    = 'unknown'
    integer               :: micrographs_imported     = 0  ! total micrographs received from import
    integer               :: micrographs_accepted     = 0  ! micrographs passing acceptance criteria
    integer               :: micrographs_rejected     = 0  ! micrographs_imported - micrographs_accepted
    integer               :: particles_extracted      = 0  ! total particles picked across accepted mics
    integer               :: particles_per_mic        = 0  ! average particles per accepted micrograph
    integer               :: last_micrograph_imported = 0  ! Unix timestamp of most recent import event
    integer               :: box_size                 = 0  ! particle box size in pixels
  contains
    procedure :: kill => kill_override
    procedure :: set
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_picking

contains

  ! Assign all fields. Derives micrographs_rejected and particles_per_mic
  ! automatically. Updates last_micrograph_imported only when the import
  ! count changes. Guards against division by zero when no mics are present.
  subroutine set( self, stage, micrographs_imported, micrographs_accepted, particles_extracted, box_size )
    class(gui_metadata_stream_picking), intent(inout) :: self
    type(string),                               intent(in)    :: stage
    integer,                                    intent(in)    :: micrographs_imported, micrographs_accepted
    integer,                                    intent(in)    :: particles_extracted, box_size
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned           = .true.
    self%stage                = stage%to_char()
    if( micrographs_imported /= self%micrographs_imported ) &
      self%last_micrograph_imported = int(c_time(0_c_long))
    self%micrographs_imported = micrographs_imported
    self%micrographs_accepted = micrographs_accepted
    self%particles_extracted  = particles_extracted
    self%box_size             = box_size
    self%micrographs_rejected = self%micrographs_imported - self%micrographs_accepted
    if( self%micrographs_accepted > 0 ) then
      self%particles_per_mic  = nint(real(self%particles_extracted) / real(self%micrographs_accepted))
    else
      self%particles_per_mic  = 0
    end if
  end subroutine set

  ! Retrieve all fields. Returns .true. if the object has been assigned.
  function get( self, stage, micrographs_imported, micrographs_accepted, micrographs_rejected, &
                particles_extracted, particles_per_mic, last_micrograph_imported, box_size ) result( l_assigned )
    class(gui_metadata_stream_picking), intent(inout) :: self
    type(string),                               intent(out)   :: stage
    integer,                                    intent(out)   :: micrographs_imported, micrographs_accepted
    integer,                                    intent(out)   :: micrographs_rejected, particles_extracted
    integer,                                    intent(out)   :: particles_per_mic, last_micrograph_imported
    integer,                                    intent(out)   :: box_size
    logical                                                   :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned               = self%l_assigned
    stage                    = trim(self%stage)
    micrographs_imported     = self%micrographs_imported
    micrographs_accepted     = self%micrographs_accepted
    micrographs_rejected     = self%micrographs_rejected
    particles_extracted      = self%particles_extracted
    particles_per_mic        = self%particles_per_mic
    last_micrograph_imported = self%last_micrograph_imported
    box_size                 = self%box_size
  end function get

  ! Serialise to a JSON object containing all fields. Returns a null pointer
  ! when the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_picking), intent(inout) :: self
    type(json_core)                                           :: json
    type(json_value),                           pointer       :: json_ptr
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'stage',                    trim(self%stage)              )
      call json%add(json_ptr, 'micrographs_imported',     self%micrographs_imported     )
      call json%add(json_ptr, 'micrographs_accepted',     self%micrographs_accepted     )
      call json%add(json_ptr, 'micrographs_rejected',     self%micrographs_rejected     )
      call json%add(json_ptr, 'particles_extracted',      self%particles_extracted      )
      call json%add(json_ptr, 'particles_per_mic',        self%particles_per_mic        )
      call json%add(json_ptr, 'last_micrograph_imported', self%last_micrograph_imported )
      call json%add(json_ptr, 'box_size',                 self%box_size                 )
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_stream_picking), intent(inout) :: self
    select type( self )
      type is( gui_metadata_stream_picking )
        self = gui_metadata_stream_picking()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_stream_picking
