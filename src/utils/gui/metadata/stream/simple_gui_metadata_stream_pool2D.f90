!@descr: GUI metadata for the stream pool-2D stage — particle counts, mask geometry, initial reference selection, and user-input flag
!==============================================================================
! MODULE: simple_gui_metadata_stream_pool2D
!
! PURPOSE:
!   Extends gui_metadata_base with fields specific to the pool-2D
!   stage of the cryo-EM streaming pipeline.  Tracks the number of particles
!   imported and accepted, the rejection count, the initial reference
!   selection array, and a flag indicating whether user input has been
!   provided.  The Unix timestamp of the most recently imported particle
!   batch is recorded and updated on each change.
!
! TYPES:
!   gui_metadata_stream_pool2D — extends gui_metadata_base
!     set()                      — assign particle counts, mask geometry, and stage field
!     set_user_input()           — set the user-input flag independently
!     set_initial_ref_selection() — append an index to the ref-selection array
!     get()                      — retrieve particle counts, user-input flag,
!                                  and last-import timestamp; returns l_assigned
!     jsonise()                  — serialise all fields to a json_value tree
!                                  (base override)
!
! DEPENDENCIES:
!   unix, json_module, simple_string, simple_defs, simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_stream_pool2D
  use unix,                     only: c_long, c_time
  use simple_error,             only: simple_exception
  use json_module,              only: json_core, json_value
  use simple_defs,              only: STDLEN
  use simple_string,            only: string
  use simple_gui_metadata_base, only: gui_metadata_base


  implicit none

  public :: gui_metadata_stream_pool2D
  private
#include "simple_local_flags.inc"

  type, extends(gui_metadata_base) :: gui_metadata_stream_pool2D
    private
    
    character(len=STDLEN) :: stage                   = 'unknown'
    integer(kind=2)       :: initial_ref_selection(1500) = 0
    integer               :: n_initial_ref_selection = 0
    integer               :: iteration               = 0       ! current 2-D classification iteration
    integer               :: particles_imported      = 0       ! total particles received from upstream
    integer               :: particles_accepted      = 0       ! particles passing 2-D selection criteria
    integer               :: particles_rejected      = 0       ! particles_imported - particles_accepted
    integer               :: last_import_time        = 0       ! Unix timestamp of most recent import event
    integer               :: mskdiam                 = 0
    logical               :: user_input              = .false. ! .true. once the user has supplied input
    real                  :: mskscale                = 0.0
  contains
    procedure :: set
    procedure :: set_user_input
    procedure :: set_initial_ref_selection
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_pool2D

contains

  ! Assign particle counts and masking fields. Derives particles_rejected
  ! automatically. Updates last_particles_imported only when the import
  ! count changes.
  subroutine set( self, stage, iteration, particles_imported, particles_accepted, particles_rejected, mskdiam, mskscale )
    class(gui_metadata_stream_pool2D), intent(inout) :: self
    type(string), intent(in) :: stage
    integer,      intent(in) :: iteration
    integer,      intent(in) :: particles_imported, particles_accepted, particles_rejected
    integer,      intent(in) :: mskdiam
    real,         intent(in) :: mskscale
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned         = .true.
    self%stage              = stage%to_char()
    self%iteration          = iteration
    if( particles_imported /= self%particles_imported ) &
      self%last_import_time = int(c_time(0_c_long))
    self%particles_imported = particles_imported
    self%particles_accepted = particles_accepted
    self%particles_rejected = particles_rejected
    self%mskdiam            = mskdiam
    self%mskscale           = mskscale
  end subroutine set

  ! Set the user-input flag. May be called independently of set().
  subroutine set_user_input( self, user_input )
    class(gui_metadata_stream_pool2D), intent(inout) :: self
    logical,                                     intent(in)    :: user_input
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%user_input = user_input
  end subroutine set_user_input

  subroutine set_initial_ref_selection( self, idx )
    class(gui_metadata_stream_pool2D), intent(inout) :: self
    integer,                                     intent(in)    :: idx
    if( .not.self%l_initialized )                THROW_HARD('gui metadata object is uninitialised')
    if( self%n_initial_ref_selection >= size(self%initial_ref_selection) ) THROW_HARD('idx is out of range')
    self%l_assigned = .true.
    self%n_initial_ref_selection = self%n_initial_ref_selection + 1
    self%initial_ref_selection(self%n_initial_ref_selection) = int2(idx)
  end subroutine set_initial_ref_selection

  ! Retrieve particle counts, user-input flag, and the last-import timestamp.
  ! Returns .true. if the object has been assigned.
  function get( self, stage, iteration, particles_imported, particles_accepted, particles_rejected, last_import_time, user_input, mskdiam, mskscale ) result( l_assigned )
    class(gui_metadata_stream_pool2D), intent(in)  :: self
    type(string), intent(out) :: stage
    integer,      intent(out) :: iteration
    integer,      intent(out) :: particles_imported, particles_accepted, particles_rejected
    integer,      intent(out) :: last_import_time
    logical,      intent(out) :: user_input
    integer,      intent(out) :: mskdiam
    real,         intent(out) :: mskscale
    logical                   :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned         = self%l_assigned
    stage              = trim(self%stage)
    iteration          = self%iteration
    particles_imported = self%particles_imported
    particles_accepted = self%particles_accepted
    particles_rejected = self%particles_rejected
    last_import_time   = self%last_import_time
    user_input         = self%user_input
    mskdiam            = self%mskdiam
    mskscale           = self%mskscale
  end function get

  ! Serialise all fields to a JSON object. Returns a null pointer when
  ! the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_pool2D), intent(inout) :: self
    type(json_core)                                            :: json
    type(json_value),                             pointer      :: json_ptr, json_ref_selection_ptr => null()
    integer                                                    :: i_ref
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'stage',                   trim(self%stage)        )
      call json%add(json_ptr, 'iteration',               self%iteration          )
      call json%add(json_ptr, 'particles_imported',      self%particles_imported )
      call json%add(json_ptr, 'particles_accepted',      self%particles_accepted )
      call json%add(json_ptr, 'particles_rejected',      self%particles_rejected )
      call json%add(json_ptr, 'last_import_time',        self%last_import_time   )
      call json%add(json_ptr, 'user_input',              self%user_input         )
      call json%add(json_ptr, 'mskdiam',                 self%mskdiam            )
      call json%add(json_ptr, 'mskscale',                dble(self%mskscale)     )
      if( self%n_initial_ref_selection > 0 ) then
        call json%create_array(json_ref_selection_ptr, 'initial_ref_selection')
        do i_ref = 1, self%n_initial_ref_selection
          call json%add(json_ref_selection_ptr, '', int(self%initial_ref_selection(i_ref)))
        end do
        call json%add(json_ptr, json_ref_selection_ptr)
      end if
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

end module simple_gui_metadata_stream_pool2D
