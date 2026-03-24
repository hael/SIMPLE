!@descr: Utility functions for GUI metadata types.
!==============================================================================
! MODULE: simple_gui_metadata_utils
!
! PURPOSE:
!   Provides helpers that operate across all gui_metadata types.
!     max_metadata_size — return sizeof() the largest concrete metadata type,
!                         used to size fixed receive buffers for IPC transfer.
!
! DEPENDENCIES:
!   simple_gui_metadata_api
!==============================================================================
module simple_gui_metadata_utils
use simple_gui_metadata_api

implicit none

public :: max_metadata_size
private
#include "simple_local_flags.inc"

contains

  ! Return the size in bytes of the largest concrete gui_metadata type.
  function max_metadata_size() result( max_size )
    type(gui_metadata_base)                     :: meta_base
    type(gui_metadata_micrograph)               :: meta_micrograph
    type(gui_metadata_histogram)                :: meta_histogram
    type(gui_metadata_timeplot)                 :: meta_timeplot
    type(gui_metadata_optics_group)             :: meta_optics_group
    type(gui_metadata_stream_update)            :: meta_update
    type(gui_metadata_stream_preprocess)        :: meta_preprocess
    type(gui_metadata_stream_optics_assignment) :: meta_optics_assignment
    type(gui_metadata_stream_initial_picking)   :: meta_initial_picking
    type(gui_metadata_stream_opening2D)         :: meta_opening2D
    integer                                     :: max_size
    max_size = max(sizeof(meta_base),              &
                   sizeof(meta_micrograph),        &
                   sizeof(meta_histogram),         &
                   sizeof(meta_timeplot),          &
                   sizeof(meta_optics_group),      &
                   sizeof(meta_update),            &
                   sizeof(meta_preprocess),        &
                   sizeof(meta_optics_assignment), &
                   sizeof(meta_initial_picking),   &
                   sizeof(meta_opening2D))
  end function max_metadata_size

end module simple_gui_metadata_utils
