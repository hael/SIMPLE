!@descr: Aggregated public API for gui metadata modules
module simple_gui_metadata_api
  use json_kinds
  use json_module
  use simple_core_module_api
  use simple_gui_metadata_types
  use simple_gui_metadata_base,              only: gui_metadata_base
  use simple_gui_metadata_micrograph,        only: gui_metadata_micrograph
  use simple_gui_metadata_histogram,         only: gui_metadata_histogram
  use simple_gui_metadata_stream_preprocess, only: gui_metadata_stream_preprocess
end module simple_gui_metadata_api
