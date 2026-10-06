!@descr: Integer type-tag constants for all GUI metadata kinds.
! A tag names a message role, not a Fortran type (gui_metadata_micrograph alone serves several);
! the stream master's pipe listener dispatches on it. Values are enum-assigned and never
! persisted, so inserting a tag may renumber the later ones.
module simple_gui_metadata_types

implicit none

enum, bind(c)
  ! standalone display types
  enumerator :: GUI_METADATA_MICROGRAPH_TYPE   = 1
  enumerator :: GUI_METADATA_HISTOGRAM_TYPE        ! 2
  enumerator :: GUI_METADATA_TIMEPLOT_TYPE         ! 3
  enumerator :: GUI_METADATA_OPTICS_GROUP_TYPE     ! 4
  enumerator :: GUI_METADATA_CAVG2D_TYPE           ! 5
  enumerator :: GUI_METADATA_VOL3D_TYPE            ! 6
  enumerator :: GUI_METADATA_PTCL_TYPE             ! 7
  enumerator :: GUI_METADATA_PROJECT_TYPE          ! 8
  ! stream control
  enumerator :: GUI_METADATA_STREAM_UPDATE_TYPE    ! 9
  ! preprocess stage
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TYPE                   ! 10
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE        ! 11
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE  ! 12
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE ! 13
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE   ! 14
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE    ! 15
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE   ! 16
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE       ! 17
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE     ! 18
  ! optics assignment stage
  enumerator :: GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE              ! 19
  enumerator :: GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE ! 20
  ! initial picking stage
  enumerator :: GUI_METADATA_STREAM_INITIAL_PICKING_TYPE              ! 21
  enumerator :: GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE   ! 22
  ! initial analysis stage
  enumerator :: GUI_METADATA_STREAM_INITIAL_ANALYSIS_TYPE             ! 23
  enumerator :: GUI_METADATA_STREAM_INITIAL_ANALYSIS_CLS2D_TYPE       ! 24
  enumerator :: GUI_METADATA_STREAM_INITIAL_ANALYSIS_CLS2D_FINAL_TYPE ! 25
  enumerator :: GUI_METADATA_STREAM_INITIAL_ANALYSIS_VOL3D_TYPE       ! 26
  ! reference picking stage
  enumerator :: GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE            ! 27
  enumerator :: GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE ! 28
  enumerator :: GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE      ! 29
  ! particle sieving stage
  enumerator :: GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE           ! 30
  enumerator :: GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE     ! 31
  enumerator :: GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_REF_TYPE ! 32
  ! pool 2D stage
  enumerator :: GUI_METADATA_STREAM_POOL2D_TYPE               ! 33
  enumerator :: GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE         ! 34
  enumerator :: GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE      ! 35
  enumerator :: GUI_METADATA_STREAM_POOL2D_SNAPSHOT_CLS2D_TYPE ! 36
  ! multistate solve3D stage
  enumerator :: GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE          ! 37
  enumerator :: GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE   ! 38
  enumerator :: GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE            ! 39
end enum

end module simple_gui_metadata_types
