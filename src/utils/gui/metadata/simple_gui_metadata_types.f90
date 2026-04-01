!@descr: Integer type-tag constants for all GUI metadata kinds.
!==============================================================================
! MODULE: simple_gui_metadata_types
!
! PURPOSE:
!   Defines the integer parameters used as meta_type tags in gui_metadata_base.
!   Tags are assigned via enum, bind(c) so values are contiguous and
!   self-maintaining — inserting a new tag only renumbers later entries.
!   Each concrete gui_metadata subtype has exactly one tag.  Sub-message types
!   (e.g. per-micrograph rows or per-histogram variants within a stage) are
!   assigned sequential values immediately after their parent stage tag.
!
! USAGE:
!   Pass a tag to gui_metadata_base%new() to identify the subtype at runtime.
!   Tags are also used by IPC consumers to dispatch incoming buffers.
!==============================================================================
module simple_gui_metadata_types

implicit none

enum, bind(c)
  ! standalone display types
  enumerator :: GUI_METADATA_MICROGRAPH_TYPE   = 1
  enumerator :: GUI_METADATA_HISTOGRAM_TYPE        ! 2
  enumerator :: GUI_METADATA_TIMEPLOT_TYPE         ! 3
  enumerator :: GUI_METADATA_OPTICS_GROUP_TYPE     ! 4
  enumerator :: GUI_METADATA_CAVG2D_TYPE           ! 5
  ! stream control
  enumerator :: GUI_METADATA_STREAM_UPDATE_TYPE    ! 6
  ! preprocess stage
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TYPE                   ! 7
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE        ! 8
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE  ! 9
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE ! 10
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE   ! 11
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE    ! 12
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE   ! 13
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE       ! 14
  enumerator :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE     ! 15
  ! optics assignment stage
  enumerator :: GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE              ! 16
  enumerator :: GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE ! 17
  ! initial picking stage
  enumerator :: GUI_METADATA_STREAM_INITIAL_PICKING_TYPE              ! 18
  enumerator :: GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE   ! 19
  ! opening 2D stage
  enumerator :: GUI_METADATA_STREAM_OPENING2D_TYPE             ! 20
  enumerator :: GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE       ! 21
  enumerator :: GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE ! 22
  ! reference picking stage
  enumerator :: GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE            ! 23
  enumerator :: GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE ! 24
  enumerator :: GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE      ! 25
  ! particle sieving stage
  enumerator :: GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE           ! 26
  enumerator :: GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE     ! 27
  enumerator :: GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_REF_TYPE ! 28
  ! pool 2D stage
  enumerator :: GUI_METADATA_STREAM_POOL2D_TYPE           ! 29
  enumerator :: GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE     ! 30
end enum

end module simple_gui_metadata_types
