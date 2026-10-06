!@descr: Utility functions for GUI metadata types.
! max_metadata_size sizes the stream IPC receive buffers: a metadata type newly sent over the
! pipes must be added to its list.
module simple_gui_metadata_utils
use simple_gui_metadata_base,                      only: gui_metadata_base
use simple_gui_metadata_micrograph,                only: gui_metadata_micrograph
use simple_gui_metadata_histogram,                 only: gui_metadata_histogram
use simple_gui_metadata_timeplot,                  only: gui_metadata_timeplot
use simple_gui_metadata_optics_group,              only: gui_metadata_optics_group
use simple_gui_metadata_cavg2D,                    only: gui_metadata_cavg2D
use simple_gui_metadata_vol3D,                     only: gui_metadata_vol3D
use simple_gui_metadata_stream_update,             only: gui_metadata_stream_update
use simple_gui_metadata_stream_preprocess,         only: gui_metadata_stream_preprocess
use simple_gui_metadata_stream_optics_assignment,  only: gui_metadata_stream_optics_assignment
use simple_gui_metadata_stream_picking,            only: gui_metadata_stream_picking
use simple_gui_metadata_stream_initial_analysis,   only: gui_metadata_stream_initial_analysis
use simple_gui_metadata_stream_particle_sieving,   only: gui_metadata_stream_particle_sieving
use simple_gui_metadata_stream_pool2D,             only: gui_metadata_stream_pool2D
use simple_gui_metadata_stream_pool2D_snapshot,    only: gui_metadata_stream_pool2D_snapshot
use simple_gui_metadata_stream_solve3D_multistate, only: gui_metadata_stream_solve3D_multistate

implicit none

public :: max_metadata_size
private

contains

  ! Return the size in bytes of the largest concrete gui_metadata type.
  function max_metadata_size() result( max_size )
    type(gui_metadata_base)                      :: meta_base
    type(gui_metadata_micrograph)                :: meta_micrograph
    type(gui_metadata_histogram)                 :: meta_histogram
    type(gui_metadata_timeplot)                  :: meta_timeplot
    type(gui_metadata_optics_group)              :: meta_optics_group
    type(gui_metadata_cavg2D)                    :: meta_cavg2D
    type(gui_metadata_stream_update)             :: meta_update
    type(gui_metadata_stream_preprocess)         :: meta_preprocess
    type(gui_metadata_stream_optics_assignment)  :: meta_optics_assignment
    type(gui_metadata_stream_picking)            :: meta_picking
    type(gui_metadata_stream_initial_analysis)   :: meta_initial_analysis
    type(gui_metadata_stream_particle_sieving)   :: meta_particle_sieving
    type(gui_metadata_stream_pool2D)             :: meta_pool2D
    type(gui_metadata_stream_pool2D_snapshot)    :: meta_pool2D_snapshot
    type(gui_metadata_stream_solve3D_multistate) :: meta_solve3D_multistate
    type(gui_metadata_vol3D)                     :: meta_vol3D
    integer                                      :: max_size
    max_size = max(sizeof(meta_base),               &
                   sizeof(meta_micrograph),         &
                   sizeof(meta_histogram),          &
                   sizeof(meta_timeplot),           &
                   sizeof(meta_optics_group),       &
                   sizeof(meta_cavg2D),             &
                   sizeof(meta_update),             &
                   sizeof(meta_preprocess),         &
                   sizeof(meta_optics_assignment),  &
                   sizeof(meta_picking),            &
                   sizeof(meta_initial_analysis),   &
                   sizeof(meta_particle_sieving),   &
                   sizeof(meta_pool2D),             &
                   sizeof(meta_pool2D_snapshot),    &
                   sizeof(meta_solve3D_multistate), &
                   sizeof(meta_vol3D))
  end function max_metadata_size

end module simple_gui_metadata_utils
