!@descr: gui metadata stream preprocess structure
module simple_gui_metadata_stream_preprocess
use json_kinds
use json_module
use simple_core_module_api
use simple_gui_metadata_base

implicit none

public :: gui_metadata_stream_preprocess
private
#include "simple_local_flags.inc"

type, extends( gui_metadata_base ) :: gui_metadata_stream_preprocess
  private
  character(len=STDLEN) :: stage               = 'unknown'
  integer               :: movies_imported     = 0
  integer               :: movies_processed    = 0
  integer               :: movies_rejected     = 0
  integer               :: movies_rate         = 0
  real                  :: average_ctf_res     = 0.0
  real                  :: average_ice_score   = 0.0
  real                  :: average_astigmatism = 0.0
  real                  :: cutoff_ctf_res      = 0.0
  real                  :: cutoff_ice_score    = 0.0
  real                  :: cutoff_astigmatism  = 0.0
  !last_movie_imported

  contains

  procedure :: set
  procedure :: get
  procedure :: jsonise => jsonise_override

end type gui_metadata_stream_preprocess

contains

  subroutine set( self, stage, movies_imported, movies_processed, movies_rejected, movies_rate,&
                  average_ctf_res, average_ice_score, average_astigmatism, cutoff_ctf_res,&
                  cutoff_ice_score, cutoff_astigmatism)
    class(gui_metadata_stream_preprocess), intent(inout) :: self
    type(string),                          intent(in)    :: stage 
    integer,                               intent(in)    :: movies_imported, movies_processed
    integer,                               intent(in)    :: movies_rejected, movies_rate
    real,                                  intent(in)    :: average_ctf_res, average_ice_score
    real,                                  intent(in)    :: average_astigmatism, cutoff_ctf_res
    real,                                  intent(in)    :: cutoff_ice_score, cutoff_astigmatism
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
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

  function get( self, stage, movies_imported, movies_processed, movies_rejected, movies_rate,&
                average_ctf_res, average_ice_score, average_astigmatism, cutoff_ctf_res,&
                cutoff_ice_score, cutoff_astigmatism ) result( l_assigned )
    class(gui_metadata_stream_preprocess), intent(inout) :: self
    type(string),                          intent(out)   :: stage                
    integer,                               intent(out)   :: movies_imported, movies_processed
    integer,                               intent(out)   :: movies_rejected, movies_rate
    real,                                  intent(out)   :: average_ctf_res, average_ice_score
    real,                                  intent(out)   :: average_astigmatism, cutoff_ctf_res
    real,                                  intent(out)   :: cutoff_ice_score, cutoff_astigmatism
    logical                                              :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
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

  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_preprocess), intent(inout) :: self
    type(json_core)                                      :: json
    type(json_value),                      pointer       :: json_ptr
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, "stage",               trim(self%stage)              )
      call json%add(json_ptr, "movies_imported",     self%movies_imported          )
      call json%add(json_ptr, "movies_processed",    self%movies_processed         )
      call json%add(json_ptr, "movies_rejected",     self%movies_rejected          )         
      call json%add(json_ptr, "movies_rate",         self%movies_rate              )             
      call json%add(json_ptr, "average_ctf_res",     dble(self%average_ctf_res)    )
      call json%add(json_ptr, "average_ice_score",   dble(self%average_ice_score)  )
      call json%add(json_ptr, "average_astigmatism", dble(self%average_astigmatism))
      call json%add(json_ptr, "cutoff_ctf_res",      dble(self%cutoff_ctf_res)     )
      call json%add(json_ptr, "cutoff_ice_score",    dble(self%cutoff_ice_score)   )
      call json%add(json_ptr, "cutoff_astigmatism",  dble(self%cutoff_astigmatism) )
    endif
  end function jsonise_override

end module simple_gui_metadata_stream_preprocess