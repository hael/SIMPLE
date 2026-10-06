!@descr: GUI metadata for the stream initial analysis stage — particle counts and masking parameters
! Filled by stream p03 (send_initial_analysis_status). particles_rejected = imported - accepted;
! last_particles_imported is a Unix timestamp, stamped when particles_imported changes.
! mask_diam is in A; mask_scale (JSON mskscale) is the box size in A, for overlay scaling.
! The assembler sends it under the JSON key opening2D, which NICE reads.
module simple_gui_metadata_stream_initial_analysis
  use unix,                     only: c_long, c_time
  use json_module,              only: json_core, json_value
  use simple_defs,              only: STDLEN
  use simple_string,            only: string
  use simple_error,             only: simple_exception
  use simple_gui_metadata_base, only: gui_metadata_base

  implicit none

  public :: gui_metadata_stream_initial_analysis
  private
#include "simple_local_flags.inc"

  type, extends(gui_metadata_base) :: gui_metadata_stream_initial_analysis
    private
    character(len=STDLEN) :: stage                        = 'unknown'
    integer               :: particles_imported           = 0       ! total particles received from upstream
    integer               :: particles_accepted           = 0       ! particles passing 2-D selection criteria
    integer               :: particles_rejected           = 0       ! particles_imported - particles_accepted
    integer               :: mask_diam                    = 0       ! mask diameter (A)
    integer               :: box_size                     = 0       ! particle box size (pixels)
    integer               :: cycle                        = 0       ! initial analysis plan cycle index
    integer               :: last_particles_imported      = 0       ! Unix timestamp of most recent import event
    real                  :: mask_scale                   = 0.0     ! box size in A (box_size * smpd)
  contains
    procedure :: kill => kill_override
    procedure :: set
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_initial_analysis

contains

  ! Assign particle counts and masking fields. Derives particles_rejected
  ! automatically. Updates last_particles_imported only when the import
  ! count changes.
  subroutine set( self, stage, particles_imported, particles_accepted, mask_diam, box_size, mask_scale, cycle )
    class(gui_metadata_stream_initial_analysis), intent(inout) :: self
    type(string),                                intent(in)    :: stage
    integer,                                     intent(in)    :: particles_imported, particles_accepted
    integer,                                     intent(in)    :: mask_diam, box_size
    real,                                        intent(in)    :: mask_scale
    integer,                                     intent(in)    :: cycle
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
    self%cycle               = cycle
  end subroutine set

  ! Retrieve particle counts and the last-import timestamp.
  ! Returns .true. if the object has been assigned.
  function get( self, stage, particles_imported, particles_accepted, last_particles_imported ) result( l_assigned )
    class(gui_metadata_stream_initial_analysis), intent(in)  :: self
    type(string),                                intent(out) :: stage
    integer,                                     intent(out) :: particles_imported, particles_accepted
    integer,                                     intent(out) :: last_particles_imported
    logical                                                  :: l_assigned
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
    class(gui_metadata_stream_initial_analysis), intent(in) :: self
    type(json_core)                                         :: json
    type(json_value),                            pointer    :: json_ptr
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
      call json%add(json_ptr, 'cycle',                    self%cycle                     )
      call json%add(json_ptr, 'mskscale',                 dble(self%mask_scale)          )
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_stream_initial_analysis), intent(inout) :: self
    select type( self )
      type is( gui_metadata_stream_initial_analysis )
        self = gui_metadata_stream_initial_analysis()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_stream_initial_analysis
