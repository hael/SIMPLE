!@descr: GUI metadata for the stream preprocessing stage — movie counts, CTF quality averages, and acceptance cutoffs
! Filled by stream p01: movies_rate is the movie watcher's detection rate (movies per hour) and
! the cutoffs are the current ctfres/icefrac/astig thresholds.
module simple_gui_metadata_stream_preprocess
  use json_module,              only: json_core, json_value
  use simple_string,            only: string
  use simple_defs,              only: STDLEN
  use simple_error,             only: simple_exception
  use simple_gui_metadata_base, only: gui_metadata_base

  implicit none

  public :: gui_metadata_stream_preprocess
  private
#include "simple_local_flags.inc"

  type, extends(gui_metadata_base) :: gui_metadata_stream_preprocess
    private
    character(len=STDLEN) :: stage               = 'unknown'
    integer               :: movies_imported     = 0    ! total movies received from the file watcher
    integer               :: movies_processed    = 0    ! movies that completed motion-correction and CTF
    integer               :: movies_rejected     = 0    ! movies failing quality cutoffs
    integer               :: movies_rate         = 0    ! movie detection rate (movies per hour)
    real                  :: average_ctf_res     = 0.0  ! mean CTF resolution estimate (Å)
    real                  :: average_ice_score   = 0.0  ! mean ice-contamination score
    real                  :: average_astigmatism = 0.0  ! mean astigmatism magnitude (Å)
    real                  :: cutoff_ctf_res      = 0.0  ! rejection threshold: CTF resolution (Å)
    real                  :: cutoff_ice_score    = 0.0  ! rejection threshold: ice score
    real                  :: cutoff_astigmatism  = 0.0  ! rejection threshold: astigmatism (Å)
  contains
    procedure :: kill => kill_override
    procedure :: set
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_preprocess

contains

  ! Assign all fields from caller-supplied values.
  subroutine set( self, stage, movies_imported, movies_processed, movies_rejected, movies_rate, &
                  average_ctf_res, average_ice_score, average_astigmatism,                      &
                  cutoff_ctf_res, cutoff_ice_score, cutoff_astigmatism )
    class(gui_metadata_stream_preprocess), intent(inout) :: self
    type(string),                          intent(in)    :: stage
    integer,                               intent(in)    :: movies_imported, movies_processed
    integer,                               intent(in)    :: movies_rejected, movies_rate
    real,                                  intent(in)    :: average_ctf_res, average_ice_score
    real,                                  intent(in)    :: average_astigmatism, cutoff_ctf_res
    real,                                  intent(in)    :: cutoff_ice_score, cutoff_astigmatism
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned          = .true.
    self%stage               = stage%to_char()
    self%movies_imported     = movies_imported
    self%movies_processed    = movies_processed
    self%movies_rejected     = movies_rejected
    self%movies_rate         = movies_rate
    self%average_ctf_res     = average_ctf_res
    self%average_ice_score   = average_ice_score
    self%average_astigmatism = average_astigmatism
    self%cutoff_ctf_res      = cutoff_ctf_res
    self%cutoff_ice_score    = cutoff_ice_score
    self%cutoff_astigmatism  = cutoff_astigmatism
  end subroutine set

  ! Retrieve all fields. Returns .true. if the object has been assigned.
  function get( self, stage, movies_imported, movies_processed, movies_rejected, movies_rate, &
                average_ctf_res, average_ice_score, average_astigmatism,                      &
                cutoff_ctf_res, cutoff_ice_score, cutoff_astigmatism ) result( l_assigned )
    class(gui_metadata_stream_preprocess), intent(in)    :: self
    type(string),                          intent(out)   :: stage
    integer,                               intent(out)   :: movies_imported, movies_processed
    integer,                               intent(out)   :: movies_rejected, movies_rate
    real,                                  intent(out)   :: average_ctf_res, average_ice_score
    real,                                  intent(out)   :: average_astigmatism, cutoff_ctf_res
    real,                                  intent(out)   :: cutoff_ice_score, cutoff_astigmatism
    logical                                              :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned          = self%l_assigned
    stage               = trim(self%stage)
    movies_imported     = self%movies_imported
    movies_processed    = self%movies_processed
    movies_rejected     = self%movies_rejected
    movies_rate         = self%movies_rate
    average_ctf_res     = self%average_ctf_res
    average_ice_score   = self%average_ice_score
    average_astigmatism = self%average_astigmatism
    cutoff_ctf_res      = self%cutoff_ctf_res
    cutoff_ice_score    = self%cutoff_ice_score
    cutoff_astigmatism  = self%cutoff_astigmatism
  end function get

  ! Serialise all fields to a JSON object. Reals are promoted to double
  ! precision for JSON fidelity. Returns a null pointer when unassigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_preprocess), intent(in)    :: self
    type(json_core)                                      :: json
    type(json_value),                      pointer       :: json_ptr
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'stage',               trim(self%stage)               )
      call json%add(json_ptr, 'movies_imported',     self%movies_imported           )
      call json%add(json_ptr, 'movies_processed',    self%movies_processed          )
      call json%add(json_ptr, 'movies_rejected',     self%movies_rejected           )
      call json%add(json_ptr, 'movies_rate',         self%movies_rate               )
      call json%add(json_ptr, 'average_ctf_res',     dble(self%average_ctf_res)     )
      call json%add(json_ptr, 'average_ice_score',   dble(self%average_ice_score)   )
      call json%add(json_ptr, 'average_astigmatism', dble(self%average_astigmatism) )
      call json%add(json_ptr, 'cutoff_ctf_res',      dble(self%cutoff_ctf_res)      )
      call json%add(json_ptr, 'cutoff_ice_score',    dble(self%cutoff_ice_score)    )
      call json%add(json_ptr, 'cutoff_astigmatism',  dble(self%cutoff_astigmatism)  )
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_stream_preprocess), intent(inout) :: self
    select type( self )
      type is( gui_metadata_stream_preprocess )
        self = gui_metadata_stream_preprocess()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_stream_preprocess
