!@descr: Integer type-tag constants for all GUI metadata kinds.
!==============================================================================
! MODULE: simple_gui_metadata_types
!
! PURPOSE:
!   Defines the integer parameters used as meta_type tags in gui_metadata_base.
!   Tags are grouped by category:
!     1–99   — standalone display types (micrograph, histogram, timeplot, ...)
!     100+   — stream pipeline types, subdivided by stage
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

! --- standalone display types ---
integer, public, parameter :: GUI_METADATA_MICROGRAPH_TYPE   = 1
integer, public, parameter :: GUI_METADATA_HISTOGRAM_TYPE    = 2
integer, public, parameter :: GUI_METADATA_TIMEPLOT_TYPE     = 3
integer, public, parameter :: GUI_METADATA_OPTICS_GROUP_TYPE = 4

! --- stream pipeline types ---

! stream control
integer, public, parameter :: GUI_METADATA_STREAM_UPDATE_TYPE = 100

! preprocess stage (101–109)
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_TYPE                   = 101
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE        = 102
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE  = 103
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE = 104
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE   = 105
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE    = 106
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE   = 107
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE       = 108
integer, public, parameter :: GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE     = 109

! optics assignment stage (110–111)
integer, public, parameter :: GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE              = 110
integer, public, parameter :: GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE = 111

! initial picking stage (112–113)
integer, public, parameter :: GUI_METADATA_STREAM_INITIAL_PICKING_TYPE              = 112
integer, public, parameter :: GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE   = 113

! opening 2D stage (114)
integer, public, parameter :: GUI_METADATA_STREAM_OPENING2D_TYPE                    = 114

end module simple_gui_metadata_types
