!@descr: Aggregated public API for gui metadata modules
module simple_gui_metadata_api
  use json_kinds
  use json_module
  use simple_core_module_api
  use simple_gui_metadata_types
  use simple_gui_metadata_base,                     only: gui_metadata_base
  use simple_gui_metadata_micrograph,               only: gui_metadata_micrograph
  use simple_gui_metadata_histogram,                only: gui_metadata_histogram
  use simple_gui_metadata_timeplot,                 only: gui_metadata_timeplot
  use simple_gui_metadata_optics_group,             only: gui_metadata_optics_group
  use simple_gui_metadata_cavg2D,                   only: gui_metadata_cavg2D, sprite_sheet_pos
  use simple_gui_metadata_stream_update,            only: gui_metadata_stream_update
  use simple_gui_metadata_stream_preprocess,        only: gui_metadata_stream_preprocess
  use simple_gui_metadata_stream_optics_assignment, only: gui_metadata_stream_optics_assignment
  use simple_gui_metadata_stream_initial_picking,   only: gui_metadata_stream_initial_picking
  use simple_gui_metadata_stream_opening2D,         only: gui_metadata_stream_opening2D
end module simple_gui_metadata_api
