!@descr: GUI metadata for the stream opening-2D stage — particle counts, masking parameters, and user-input flag
!==============================================================================
! MODULE: simple_gui_metadata_stream_opening2D
!
! PURPOSE:
!   Extends gui_metadata_base with fields specific to the opening 2-D
!   classification stage of the cryo-EM streaming pipeline.  Tracks the
!   number of particles imported and accepted, the rejection count, masking
!   geometry, and a flag indicating whether user input has been provided.
!   The Unix timestamp of the most recently imported particle batch is also
!   recorded and updated on each change.
!
! TYPES:
!   gui_metadata_stream_opening2D — extends gui_metadata_base
!     set()           — assign particle counts and masking fields
!     set_user_input() — set the user-input flag independently
!     get()           — retrieve particle counts and last-import timestamp;
!                       returns the l_assigned flag
!     jsonise()       — serialise all fields to a json_value tree (base override)
!
! DEPENDENCIES:
!   unix, json_kinds, json_module, simple_string, simple_defs,
!   simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_stream_opening2D
  use unix,                     only: c_long, c_time
  use json_kinds
  use json_module,              only: json_core, json_value
  use simple_defs,              only: STDLEN
  use simple_string,            only: string
  use simple_error,             only: simple_exception
  use simple_gui_metadata_base, only: gui_metadata_base

  implicit none

  public :: gui_metadata_stream_opening2D
  private
#include "simple_local_flags.inc"

  type, extends(gui_metadata_base) :: gui_metadata_stream_opening2D
    private
    character(len=STDLEN) :: stage                   = 'unknown'
    integer               :: particles_imported      = 0       ! total particles received from upstream
    integer               :: particles_accepted      = 0       ! particles passing 2-D selection criteria
    integer               :: particles_rejected      = 0       ! particles_imported - particles_accepted
    integer               :: mask_diam               = 0       ! circular mask diameter (pixels)
    integer               :: box_size                = 0       ! particle box size (pixels)
    integer               :: last_particles_imported = 0       ! Unix timestamp of most recent import event
    real                  :: mask_scale              = 0.0     ! fractional mask scale factor
    logical               :: user_input              = .false. ! .true. once the user has supplied input
  contains
    procedure :: set
    procedure :: set_user_input
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_opening2D

contains

  ! Assign particle counts and masking fields. Derives particles_rejected
  ! automatically. Updates last_particles_imported only when the import
  ! count changes.
  subroutine set( self, stage, particles_imported, particles_accepted, mask_diam, box_size, mask_scale )
    class(gui_metadata_stream_opening2D), intent(inout) :: self
    type(string),                         intent(in)    :: stage
    integer,                              intent(in)    :: particles_imported, particles_accepted
    integer,                              intent(in)    :: mask_diam, box_size
    real,                                 intent(in)    :: mask_scale
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned          = .true.
    self%stage               = stage%to_char()
    if( particles_imported /= self%particles_imported ) &
      self%last_particles_imported = int(c_time(0_c_long))
    self%particles_imported  = particles_imported
    self%particles_accepted  = particles_accepted
    self%particles_rejected  = self%particles_imported - self%particles_accepted
    self%mask_diam           = mask_diam
    self%box_size            = box_size
    self%mask_scale          = mask_scale
  end subroutine set

  ! Set the user-input flag. May be called independently of set().
  subroutine set_user_input( self, user_input )
    class(gui_metadata_stream_opening2D), intent(inout) :: self
    logical,                              intent(in)    :: user_input
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%user_input = user_input
  end subroutine set_user_input

  ! Retrieve particle counts and the last-import timestamp.
  ! Returns .true. if the object has been assigned.
  function get( self, stage, particles_imported, particles_accepted, last_particles_imported ) result( l_assigned )
    class(gui_metadata_stream_opening2D), intent(inout) :: self
    type(string),                         intent(out)   :: stage
    integer,                              intent(out)   :: particles_imported, particles_accepted
    integer,                              intent(out)   :: last_particles_imported
    logical                                             :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned              = self%l_assigned
    stage                   = trim(self%stage)
    particles_imported      = self%particles_imported
    particles_accepted      = self%particles_accepted
    last_particles_imported = self%last_particles_imported
  end function get

  ! Serialise all fields to a JSON object. Returns a null pointer when
  ! the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_opening2D), intent(inout) :: self
    type(json_core)                                     :: json
    type(json_value),                     pointer       :: json_ptr
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'stage',                    trim(self%stage)               )
      call json%add(json_ptr, 'particles_imported',       self%particles_imported        )
      call json%add(json_ptr, 'particles_accepted',       self%particles_accepted        )
      call json%add(json_ptr, 'particles_rejected',       self%particles_rejected        )
      call json%add(json_ptr, 'last_particles_imported',  self%last_particles_imported   )
      call json%add(json_ptr, 'mask_diam',                self%mask_diam                 )
      call json%add(json_ptr, 'box_size',                 self%box_size                  )
      call json%add(json_ptr, 'mskscale',                 dble(self%mask_scale)          )
      call json%add(json_ptr, 'user_input',               self%user_input                )
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

end module simple_gui_metadata_stream_opening2D
