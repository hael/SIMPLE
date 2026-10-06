!@descr: Utility functions for GUI metadata types.
! max_metadata_size sizes the stream IPC receive buffers: a metadata type newly sent over the
! pipes must be added to its list.
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
    type(gui_metadata_cavg2D)                   :: meta_cavg2D
    type(gui_metadata_stream_update)            :: meta_update
    type(gui_metadata_stream_preprocess)        :: meta_preprocess
    type(gui_metadata_stream_optics_assignment) :: meta_optics_assignment
    type(gui_metadata_stream_picking)          :: meta_initial_picking
    type(gui_metadata_stream_initial_analysis)         :: meta_initial_analysis
    type(gui_metadata_stream_particle_sieving)  :: meta_particle_sieving
    type(gui_metadata_stream_pool2D)            :: meta_pool2D
    type(gui_metadata_stream_pool2D_snapshot)   :: meta_pool2D_snapshot
    type(gui_metadata_stream_solve3D_multistate) :: meta_solve3D_multistate
    type(gui_metadata_vol3D)                    :: meta_vol3D
    integer                                     :: max_size
    max_size = max(sizeof(meta_base),              &
                   sizeof(meta_micrograph),        &
                   sizeof(meta_histogram),         &
                   sizeof(meta_timeplot),          &
                   sizeof(meta_optics_group),      &
                   sizeof(meta_cavg2D),            &
                   sizeof(meta_update),            &
                   sizeof(meta_preprocess),        &
                   sizeof(meta_optics_assignment), &
                   sizeof(meta_initial_picking),   &
                   sizeof(meta_initial_analysis),         &
                   sizeof(meta_particle_sieving),  &
                   sizeof(meta_pool2D),            &
                   sizeof(meta_pool2D_snapshot),   &
                   sizeof(meta_solve3D_multistate), &
                   sizeof(meta_vol3D))
  end function max_metadata_size

end module simple_gui_metadata_utils
