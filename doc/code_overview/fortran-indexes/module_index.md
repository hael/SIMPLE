# Module Index

## Module: CPlot2D_wrapper_module

Files:
- `extlibs/CPlot2D/simple_CPlot2D_wrapper.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `C_CDataPoint__delete` — subroutine
- `C_CDataPoint__new2` — function
- `C_CDataSet__AddDataPoint` — subroutine
- `C_CDataSet__delete` — subroutine
- `C_CDataSet__new` — function
- `C_CDataSet__SetDashedLine` — subroutine
- `C_CDataSet__SetDatasetColor` — subroutine
- `C_CDataSet__SetDrawLine` — subroutine
- `C_CDataSet__SetDrawMarker` — subroutine
- `C_CDataSet__SetFillArea` — subroutine
- `C_CDataSet__SetLineWidth` — subroutine
- `C_CDataSet__SetMarkerSize` — subroutine
- `C_CPlot2D__AddDataSet` — subroutine
- `C_CPlot2D__delete` — subroutine
- `C_CPlot2D__new` — function
- `C_CPlot2D__OutputPostScriptPlot` — subroutine
- `C_CPlot2D__SetDrawLegend` — subroutine
- `C_CPlot2D__SetDrawXAxisGridLines` — subroutine
- `C_CPlot2D__SetDrawYAxisGridLines` — subroutine
- `C_CPlot2D__SetFlipY` — subroutine
- `C_CPlot2D__SetXAxisSize` — subroutine
- `C_CPlot2D__SetXAxisTitle` — subroutine
- `C_CPlot2D__SetYAxisSize` — subroutine
- `C_CPlot2D__SetYAxisTitle` — subroutine
- `CDataPoint__delete` — subroutine
- `CDataPoint__new2` — subroutine
- `CDataPoint_type` — type
- `CDataSet__AddDataPoint` — subroutine
- `CDataSet__delete` — subroutine
- `CDataSet__new` — subroutine
- `CDataSet__SetDashedLine` — subroutine
- `CDataSet__SetDatasetColor` — subroutine
- `CDataSet__SetDrawLine` — subroutine
- `CDataSet__SetDrawMarker` — subroutine
- `CDataSet__SetFillArea` — subroutine
- `CDataSet__SetLineWidth` — subroutine
- `CDataSet__SetMarkerSize` — subroutine
- `CDataSet_addpoint` — subroutine
- `CDataSet_type` — type
- `CPlot2D__AddDataSet` — subroutine
- `CPlot2D__delete` — subroutine
- `CPlot2D__new` — subroutine
- `CPlot2D__OutputPostScriptPlot` — subroutine
- `CPlot2D__SetDrawLegend` — subroutine
- `CPlot2D__SetDrawXAxisGridLines` — subroutine
- `CPlot2D__SetDrawYAxisGridLines` — subroutine
- `CPlot2D__SetFlipY` — subroutine
- `CPlot2D__SetXAxisSize` — subroutine
- `CPlot2D__SetXAxisTitle` — subroutine
- `CPlot2D__SetYAxisSize` — subroutine
- `CPlot2D__SetYAxisTitle` — subroutine
- `CPlot2D_type` — type
- `plot2D` — subroutine

---
## Module: curl

Files:
- `extlibs/curl/src/curl.f90`

Uses:
- `curl_easy`
- `curl_multi`
- `curl_urlapi`
- `curl_util`

---
## Module: curl_easy

Files:
- `extlibs/curl/src/curl_easy.f90`

Uses:
- `curl_util`
- `iso_c_binding`
- `iso_fortran_env`

Public symbols:
- `curl_easy_cleanup` — subroutine
- `curl_easy_escape` — function
- `curl_easy_init` — function
- `curl_easy_pause` — function
- `curl_easy_perform` — function
- `curl_easy_setopt_char` — function
- `curl_easy_setopt_funptr` — function
- `curl_easy_setopt_int` — function
- `curl_easy_setopt_long` — function
- `curl_easy_setopt_ptr` — function
- `curl_easy_strerror` — function
- `curl_easy_unescape` — function
- `curl_escape` — function
- `curl_free` — subroutine
- `curl_global_cleanup` — subroutine
- `curl_global_init` — function
- `curl_mime_addpart` — function
- `curl_mime_data` — function
- `curl_mime_encoder` — function
- `curl_mime_filedata` — function
- `curl_mime_filename` — function
- `curl_mime_free` — subroutine
- `curl_mime_headers` — function
- `curl_mime_init` — function
- `curl_mime_name` — function
- `curl_mime_subparts` — function
- `curl_mime_type` — function
- `curl_slist` — type
- `curl_slist_append` — function
- `curl_slist_free_all` — subroutine
- `curl_unescape` — function
- `curl_version` — function
- `curl_version_info` — function
- `curl_version_info_data` — type
- `curl_version_now` — function

Private symbols:
- `curl_easy_escape_` — function
- `curl_easy_getinfo_c_double` — function
- `curl_easy_getinfo_c_long` — function
- `curl_easy_getinfo_c_ptr` — function
- `curl_easy_getinfo_char` — function
- `curl_easy_getinfo_double` — function
- `curl_easy_getinfo_int` — function
- `curl_easy_getinfo_long` — function
- `curl_easy_getinfo_ptr` — function
- `curl_easy_setopt_c_char` — function
- `curl_easy_setopt_c_funptr` — function
- `curl_easy_setopt_c_long` — function
- `curl_easy_setopt_c_ptr` — function
- `curl_easy_strerror_` — function
- `curl_easy_unescape_` — function
- `curl_escape_` — function
- `curl_global_init_` — function
- `curl_mime_encoder_` — function
- `curl_mime_filedata_` — function
- `curl_mime_filename_` — function
- `curl_mime_name_` — function
- `curl_mime_type_` — function
- `curl_slist_append_` — function
- `curl_unescape_` — function
- `curl_version_` — function
- `curl_version_info_` — function

---
## Module: curl_multi

Files:
- `extlibs/curl/src/curl_multi.f90`

Uses:
- `curl_easy`
- `curl_util`
- `iso_c_binding`

Public symbols:
- `curl_msg` — type
- `curl_multi_add_handle` — function
- `curl_multi_cleanup` — function
- `curl_multi_info_read` — function
- `curl_multi_init` — function
- `curl_multi_perform` — function
- `curl_multi_poll` — function
- `curl_multi_remove_handle` — function
- `curl_multi_strerror` — function
- `curl_multi_timeout` — function
- `curl_multi_wakeup` — function

Private symbols:
- `curl_multi_strerror_` — function

---
## Module: curl_urlapi

Files:
- `extlibs/curl/src/curl_urlapi.f90`

Uses:
- `curl_easy`
- `curl_util`
- `iso_c_binding`

Public symbols:
- `curl_url` — function
- `curl_url_cleanup` — subroutine
- `curl_url_dup` — function
- `curl_url_get` — function
- `curl_url_set` — function
- `curl_url_strerror` — function

Private symbols:
- `curl_url_get_` — function
- `curl_url_set_` — function
- `curl_url_strerror_` — function

---
## Module: curl_util

Files:
- `extlibs/curl/src/curl_util.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`

Public symbols:
- `c_f_str_ptr` — subroutine

Private symbols:
- `c_strlen` — function

---
## Module: gnufor2

Files:
- `extlibs/gnufor/gnufor2.f90`

Public symbols:
- `hist` — subroutine
- `image_1` — subroutine
- `image_2` — subroutine
- `image_3` — subroutine
- `image_4` — subroutine
- `image_5` — subroutine
- `image_6` — subroutine
- `my_date_and_time` — function
- `output_terminal` — function
- `plot3d` — subroutine
- `plot_1` — subroutine
- `plot_2` — subroutine
- `plot_3` — subroutine
- `plot_4` — subroutine
- `run_gnuplot` — subroutine
- `surf_1` — subroutine
- `surf_2` — subroutine
- `surf_3` — subroutine

---
## Module: json_file_module

Files:
- `extlibs/json/json_file_module.f90`

Uses:
- `iso_fortran_env`
- `json_kinds`
- `json_parameters`
- `json_string_utilities`
- `json_value_module`

Public symbols:
- `json_file` — type

Private symbols:
- `get_json_core_in_file` — subroutine
- `initialize_json_core_in_file` — subroutine
- `initialize_json_file` — function
- `initialize_json_file_v2` — function
- `json_file_add_double` — subroutine
- `json_file_add_double_vec` — subroutine
- `json_file_add_integer` — subroutine
- `json_file_add_integer_vec` — subroutine
- `json_file_add_logical` — subroutine
- `json_file_add_logical_vec` — subroutine
- `json_file_add_object` — subroutine
- `json_file_add_string` — subroutine
- `json_file_add_string_path_ascii` — subroutine
- `json_file_add_string_value_ascii` — subroutine
- `json_file_add_string_vec` — subroutine
- `json_file_add_string_vec_path_ascii` — subroutine
- `json_file_add_string_vec_vec_ascii` — subroutine
- `json_file_check_for_errors` — subroutine
- `json_file_clear_exceptions` — subroutine
- `json_file_destroy` — subroutine
- `json_file_failed` — function
- `json_file_get_alloc_string_vec` — subroutine
- `json_file_get_double` — subroutine
- `json_file_get_double_vec` — subroutine
- `json_file_get_integer` — subroutine
- `json_file_get_integer_vec` — subroutine
- `json_file_get_logical` — subroutine
- `json_file_get_logical_vec` — subroutine
- `json_file_get_object` — subroutine
- `json_file_get_root` — subroutine
- `json_file_get_string` — subroutine
- `json_file_get_string_vec` — subroutine
- `json_file_load` — subroutine
- `json_file_load_from_string` — subroutine
- `json_file_move_pointer` — subroutine
- `json_file_print_1` — subroutine
- `json_file_print_2` — subroutine
- `json_file_print_error_message` — subroutine
- `json_file_print_to_console` — subroutine
- `json_file_print_to_string` — subroutine
- `json_file_remove` — subroutine
- `json_file_rename` — subroutine
- `json_file_rename_name_ascii` — subroutine
- `json_file_rename_path_ascii` — subroutine
- `json_file_traverse` — subroutine
- `json_file_update_integer` — subroutine
- `json_file_update_logical` — subroutine
- `json_file_update_real` — subroutine
- `json_file_update_string` — subroutine
- `json_file_update_string_name_ascii` — subroutine
- `json_file_update_string_val_ascii` — subroutine
- `json_file_valid_path` — function
- `json_file_valid_path_op` — function
- `json_file_variable_info` — subroutine
- `json_file_variable_matrix_info` — subroutine
- `set_json_core_in_file` — subroutine
- `wrap_json_file_add_double` — subroutine
- `wrap_json_file_add_double_vec` — subroutine
- `wrap_json_file_add_integer` — subroutine
- `wrap_json_file_add_integer_vec` — subroutine
- `wrap_json_file_add_logical` — subroutine
- `wrap_json_file_add_logical_vec` — subroutine
- `wrap_json_file_add_object` — subroutine
- `wrap_json_file_add_string` — subroutine
- `wrap_json_file_add_string_vec` — subroutine
- `wrap_json_file_get_alloc_string_vec` — subroutine
- `wrap_json_file_get_double` — subroutine
- `wrap_json_file_get_double_vec` — subroutine
- `wrap_json_file_get_integer` — subroutine
- `wrap_json_file_get_integer_vec` — subroutine
- `wrap_json_file_get_logical` — subroutine
- `wrap_json_file_get_logical_vec` — subroutine
- `wrap_json_file_get_object` — subroutine
- `wrap_json_file_get_string` — subroutine
- `wrap_json_file_get_string_vec` — subroutine
- `wrap_json_file_load_from_string` — subroutine
- `wrap_json_file_remove` — subroutine
- `wrap_json_file_rename` — subroutine
- `wrap_json_file_update_integer` — subroutine
- `wrap_json_file_update_logical` — subroutine
- `wrap_json_file_update_real` — subroutine
- `wrap_json_file_update_string` — subroutine
- `wrap_json_file_valid_path` — function
- `wrap_json_file_valid_path_op` — function
- `wrap_json_file_variable_info` — subroutine
- `wrap_json_file_variable_matrix_info` — subroutine

---
## Module: json_kinds

Files:
- `extlibs/json/json_kinds.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`

---
## Module: json_module

Files:
- `extlibs/json/json_module.f90`

Uses:
- `json_file_module`
- `json_kinds`
- `json_parameters`
- `json_string_utilities`
- `json_value_module`

---
## Module: json_parameters

Files:
- `extlibs/json/json_parameters.f90`

Uses:
- `json_kinds`

---
## Module: json_string_utilities

Files:
- `extlibs/json/json_string_utilities.f90`

Uses:
- `json_kinds`
- `json_parameters`

Public symbols:
- `decode_rfc6901` — function
- `encode_rfc6901` — function
- `escape_string` — subroutine
- `integer_to_string` — subroutine
- `lowercase_string` — function
- `real_to_string` — subroutine
- `replace_string` — subroutine
- `string_to_integer` — subroutine
- `string_to_real` — subroutine
- `unescape_string` — subroutine
- `valid_json_hex` — function

Private symbols:
- `compact_real_string` — subroutine
- `default_comp_ucs4` — function
- `default_join_ucs4` — function
- `default_neq_ucs4` — function
- `to_uni` — function
- `to_uni_vec` — function
- `ucs4_comp_default` — function
- `ucs4_join_default` — function
- `ucs4_neq_default` — function

---
## Module: json_value_module

Files:
- `extlibs/json/json_value_module.f90`

Uses:
- `iso_fortran_env`
- `json_kinds`
- `json_parameters`
- `json_string_utilities`

Public symbols:
- `json_array_callback_func` — subroutine
- `json_core` — type
- `json_traverse_callback_func` — subroutine
- `json_value` — type

Private symbols:
- `annotate_invalid_json` — subroutine
- `convert` — subroutine
- `destroy_json_core` — subroutine
- `destroy_json_data` — subroutine
- `get_current_line_from_file_sequential` — subroutine
- `get_current_line_from_file_stream` — subroutine
- `initialize_json_core` — function
- `json_add_double_by_path` — subroutine
- `json_add_double_vec_by_path` — subroutine
- `json_add_integer_by_path` — subroutine
- `json_add_integer_vec_by_path` — subroutine
- `json_add_logical_by_path` — subroutine
- `json_add_logical_vec_by_path` — subroutine
- `json_add_member_by_path` — subroutine
- `json_add_string_by_path` — subroutine
- `json_add_string_by_path_path_ascii` — subroutine
- `json_add_string_by_path_value_ascii` — subroutine
- `json_add_string_vec_by_path` — subroutine
- `json_add_string_vec_by_path_path_ascii` — subroutine
- `json_add_string_vec_by_path_value_ascii` — subroutine
- `json_check_all_for_duplicate_keys` — subroutine
- `json_check_children_for_duplicate_keys` — subroutine
- `json_check_for_errors` — subroutine
- `json_clear_exceptions` — subroutine
- `json_clone` — subroutine
- `json_count` — function
- `json_create_by_path` — subroutine
- `json_failed` — function
- `json_get_alloc_string_vec` — subroutine
- `json_get_alloc_string_vec_by_path` — subroutine
- `json_get_array` — subroutine
- `json_get_array_by_path` — subroutine
- `json_get_by_path` — subroutine
- `json_get_by_path_default` — subroutine
- `json_get_by_path_jsonpath_bracket` — subroutine
- `json_get_by_path_rfc6901` — subroutine
- `json_get_double` — subroutine
- `json_get_double_by_path` — subroutine
- `json_get_double_vec` — subroutine
- `json_get_double_vec_by_path` — subroutine
- `json_get_integer` — subroutine
- `json_get_integer_by_path` — subroutine
- `json_get_integer_vec` — subroutine
- `json_get_integer_vec_by_path` — subroutine
- `json_get_logical` — subroutine
- `json_get_logical_by_path` — subroutine
- `json_get_logical_vec` — subroutine
- `json_get_logical_vec_by_path` — subroutine
- `json_get_next` — subroutine
- `json_get_parent` — subroutine
- `json_get_path` — subroutine
- `json_get_previous` — subroutine
- `json_get_string` — subroutine
- `json_get_string_by_path` — subroutine
- `json_get_string_vec` — subroutine
- `json_get_string_vec_by_path` — subroutine
- `json_get_tail` — subroutine
- `json_info` — subroutine
- `json_info_by_path` — subroutine
- `json_initialize` — subroutine
- `json_matrix_info` — subroutine
- `json_matrix_info_by_path` — subroutine
- `json_parse_file` — subroutine
- `json_parse_string` — subroutine
- `json_print_error_message` — subroutine
- `json_print_to_filename` — subroutine
- `json_print_to_unit` — subroutine
- `json_print_value_fast` — subroutine
- `json_print_value_fast_impl` — subroutine
- `json_rename_by_path` — subroutine
- `json_rename_by_path_name_ascii` — subroutine
- `json_rename_by_path_path_ascii` — subroutine
- `json_string_info` — subroutine
- `json_throw_exception` — subroutine
- `json_traverse` — subroutine
- `json_update_double` — subroutine
- `json_update_integer` — subroutine
- `json_update_logical` — subroutine
- `json_update_string` — subroutine
- `json_update_string_name_ascii` — subroutine
- `json_update_string_val_ascii` — subroutine
- `json_valid_path` — function
- `json_value_add_double` — subroutine
- `json_value_add_double_vec` — subroutine
- `json_value_add_integer` — subroutine
- `json_value_add_integer_vec` — subroutine
- `json_value_add_logical` — subroutine
- `json_value_add_logical_vec` — subroutine
- `json_value_add_member` — subroutine
- `json_value_add_null` — subroutine
- `json_value_add_string` — subroutine
- `json_value_add_string_name_ascii` — subroutine
- `json_value_add_string_val_ascii` — subroutine
- `json_value_add_string_vec` — subroutine
- `json_value_add_string_vec_name_ascii` — subroutine
- `json_value_add_string_vec_val_ascii` — subroutine
- `json_value_clone_func` — subroutine
- `json_value_create` — subroutine
- `json_value_create_array` — subroutine
- `json_value_create_double` — subroutine
- `json_value_create_integer` — subroutine
- `json_value_create_logical` — subroutine
- `json_value_create_null` — subroutine
- `json_value_create_object` — subroutine
- `json_value_create_string` — subroutine
- `json_value_destroy` — subroutine
- `json_value_get_child` — subroutine
- `json_value_get_child_by_index` — subroutine
- `json_value_get_child_by_name` — subroutine
- `json_value_insert_after` — subroutine
- `json_value_insert_after_child_by_index` — subroutine
- `json_value_is_child_of` — function
- `json_value_print` — subroutine
- `json_value_remove` — subroutine
- `json_value_remove_if_present` — subroutine
- `json_value_rename` — subroutine
- `json_value_replace` — subroutine
- `json_value_reverse` — subroutine
- `json_value_swap` — subroutine
- `json_value_to_string` — subroutine
- `json_value_to_string_fast` — subroutine
- `json_value_validate` — subroutine
- `name_equal` — function
- `name_strings_equal` — function
- `parse_array` — subroutine
- `parse_for_chars` — subroutine
- `parse_number` — subroutine
- `parse_object` — subroutine
- `parse_string` — subroutine
- `parse_value` — subroutine
- `pop_char` — subroutine
- `push_char` — subroutine
- `string_to_dble` — function
- `string_to_int` — function
- `to_array` — subroutine
- `to_double` — subroutine
- `to_integer` — subroutine
- `to_logical` — subroutine
- `to_null` — subroutine
- `to_object` — subroutine
- `to_string` — subroutine
- `wrap_json_add_double_by_path` — subroutine
- `wrap_json_add_double_vec_by_path` — subroutine
- `wrap_json_add_integer_by_path` — subroutine
- `wrap_json_add_integer_vec_by_path` — subroutine
- `wrap_json_add_logical_by_path` — subroutine
- `wrap_json_add_logical_vec_by_path` — subroutine
- `wrap_json_add_member_by_path` — subroutine
- `wrap_json_add_string_by_path` — subroutine
- `wrap_json_add_string_vec_by_path` — subroutine
- `wrap_json_create_by_path` — subroutine
- `wrap_json_get_alloc_string_vec_by_path` — subroutine
- `wrap_json_get_array_by_path` — subroutine
- `wrap_json_get_by_path` — subroutine
- `wrap_json_get_double_by_path` — subroutine
- `wrap_json_get_double_vec_by_path` — subroutine
- `wrap_json_get_integer_by_path` — subroutine
- `wrap_json_get_integer_vec_by_path` — subroutine
- `wrap_json_get_logical_by_path` — subroutine
- `wrap_json_get_logical_vec_by_path` — subroutine
- `wrap_json_get_path` — subroutine
- `wrap_json_get_string_by_path` — subroutine
- `wrap_json_get_string_vec_by_path` — subroutine
- `wrap_json_info_by_path` — subroutine
- `wrap_json_matrix_info_by_path` — subroutine
- `wrap_json_parse_string` — subroutine
- `wrap_json_rename_by_path` — subroutine
- `wrap_json_throw_exception` — subroutine
- `wrap_json_update_double` — subroutine
- `wrap_json_update_integer` — subroutine
- `wrap_json_update_logical` — subroutine
- `wrap_json_update_string` — subroutine
- `wrap_json_valid_path` — function
- `wrap_json_value_add_double` — subroutine
- `wrap_json_value_add_double_vec` — subroutine
- `wrap_json_value_add_integer` — subroutine
- `wrap_json_value_add_integer_vec` — subroutine
- `wrap_json_value_add_logical` — subroutine
- `wrap_json_value_add_logical_vec` — subroutine
- `wrap_json_value_add_null` — subroutine
- `wrap_json_value_add_string` — subroutine
- `wrap_json_value_add_string_vec` — subroutine
- `wrap_json_value_create_array` — subroutine
- `wrap_json_value_create_double` — subroutine
- `wrap_json_value_create_integer` — subroutine
- `wrap_json_value_create_logical` — subroutine
- `wrap_json_value_create_null` — subroutine
- `wrap_json_value_create_object` — subroutine
- `wrap_json_value_create_string` — subroutine
- `wrap_json_value_get_child_by_name` — subroutine
- `wrap_json_value_remove_if_present` — subroutine
- `wrap_json_value_rename` — subroutine

Local symbols (enclosing procedure scope):
- `json_check_all_for_duplicate_keys::duplicate_key_func` — subroutine
- `json_check_children_for_duplicate_keys::alloc_str` — type
- `json_get_alloc_string_vec::get_chars_from_array` — subroutine
- `json_get_double_vec::get_double_from_array` — subroutine
- `json_get_integer_vec::get_int_from_array` — subroutine
- `json_get_logical_vec::get_logical_from_array` — subroutine
- `json_get_path::add_to_path` — subroutine
- `json_get_string_vec::get_chars_from_array` — subroutine
- `json_print_value_fast_impl::append_text` — subroutine
- `json_print_value_fast_impl::ensure_capacity` — subroutine
- `json_print_value_fast_impl::write_it_fast` — subroutine
- `json_string_info::get_string_lengths` — subroutine
- `json_traverse::traverse` — subroutine
- `json_value_is_child_of::is_child_of_callback` — subroutine
- `json_value_print::write_it` — subroutine
- `json_value_swap::swap_pointers` — subroutine
- `json_value_validate::check_if_valid` — subroutine

---
## Module: simple_accum_blend_tester

Files:
- `main/image/simple_accum_blend_tester.f90`

Uses:
- `simple_image`
- `simple_oris`
- `simple_test_utils`

Public symbols:
- `run_all_accum_blend_tests` — subroutine

Private symbols:
- `test_trail_rec_blend` — subroutine
- `test_trail_rec_cohort` — subroutine
- `test_trail_rec_population` — subroutine

Local symbols (enclosing procedure scope):
- `test_trail_rec_blend::make_accum` — subroutine
- `test_trail_rec_blend::restored_at_origin` — function
- `test_trail_rec_blend::run_chain_start_iteration` — subroutine
- `test_trail_rec_blend::run_chain_start_then_update` — subroutine
- `test_trail_rec_blend::run_recurrence` — subroutine
- `test_trail_rec_blend::test_scale_mats_primitive` — subroutine
- `test_trail_rec_cohort::make_unit_density_accum` — subroutine
- `test_trail_rec_population::blend` — subroutine
- `test_trail_rec_population::make_unit_accum` — subroutine

---
## Module: simple_aff_prop

Files:
- `utils/clustering/simple_aff_prop.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `aff_prop` — type

Private symbols:
- `deterministic_tiebreak` — function
- `kill` — subroutine
- `new` — subroutine
- `propagate` — subroutine

Local symbols (enclosing procedure scope):
- `propagate::calc_exemplar_mask` — subroutine

---
## Module: simple_aff_prop_tester

Files:
- `utils/clustering/simple_aff_prop_tester.f90`

Uses:
- `simple_aff_prop`
- `simple_defs`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_aff_prop_tests` — subroutine

Private symbols:
- `make_data` — subroutine
- `medoids` — function
- `test_determinism` — subroutine
- `test_preference_limits` — subroutine
- `test_three_clusters` — subroutine

---
## Module: simple_ansi_ctrls

Files:
- `defs/simple_ansi_ctrls.f90`

Public symbols:
- `format_str` — function

---
## Module: simple_ansi_ctrls_tester

Files:
- `defs/simple_ansi_ctrls_tester.f90`

Uses:
- `simple_ansi_ctrls`
- `simple_test_utils`

Public symbols:
- `run_all_ansi_ctrls_tests` — subroutine

---
## Module: simple_atomic_species

Files:
- `main/nano/simple_atomic_species.f90`

Uses:
- `simple_atoms`
- `simple_core_module_api`
- `simple_defs`
- `simple_image`
- `simple_nano_species`

Public symbols:
- `call_atomic_species` — subroutine

Private symbols:
- `read_atoms_table` — subroutine
- `read_model` — subroutine
- `split_fields` — subroutine

Local symbols (enclosing procedure scope):
- `call_atomic_species::apertures` — function
- `call_atomic_species::bin_of` — function
- `call_atomic_species::class_shells` — subroutine
- `call_atomic_species::expected_int` — function
- `call_atomic_species::half_apertures` — subroutine
- `call_atomic_species::halfmap_agreement` — subroutine
- `call_atomic_species::median_of` — function
- `call_atomic_species::posterior_model_intensity` — function
- `call_atomic_species::rstr` — function
- `call_atomic_species::sdev` — function
- `call_atomic_species::species_occupancy` — function
- `call_atomic_species::tau_in_bin` — function
- `call_atomic_species::within_class_tau` — function
- `call_atomic_species::write_kv` — subroutine
- `call_atomic_species::write_radial_profiles` — subroutine
- `call_atomic_species::write_species_pdb` — subroutine
- `call_atomic_species::write_species_report` — subroutine
- `call_atomic_species::write_species_table` — subroutine
- `read_atoms_table::col` — function

---
## Module: simple_atoms

Files:
- `main/nano/simple_atoms.f90`

Uses:
- `simple_core_module_api`
- `simple_defs_atoms`
- `simple_image`
- `simple_molecule_data`
- `simple_string_utils`

Public symbols:
- `atoms` — type
- `parse_element_list` — function

Private symbols:
- `atom_validate` — subroutine
- `center_inbox` — subroutine
- `center_pdbcoord` — subroutine
- `check_center` — function
- `cif2mrc` — subroutine
- `cif2pdb` — subroutine
- `convolve` — subroutine
- `copy` — subroutine
- `element_exists` — function
- `extract_atom` — subroutine
- `get_atomicnumber` — function
- `get_beta` — function
- `get_coord` — function
- `get_element` — function
- `get_geom_center` — function
- `get_n` — function
- `get_name` — function
- `get_nres` — function
- `get_occupancy` — function
- `get_radius` — function
- `get_resnum` — function
- `guess_an_element` — subroutine
- `guess_element` — subroutine
- `kill` — subroutine
- `map_validate` — subroutine
- `model_validate` — subroutine
- `model_validate_eo` — subroutine
- `new_from_cif` — subroutine
- `new_from_file` — subroutine
- `new_from_molecule` — subroutine
- `new_from_pdb` — subroutine
- `new_instance` — subroutine
- `pdb2mrc` — subroutine
- `sample_atom_window` — subroutine
- `set_atom_corr` — subroutine
- `set_beta` — subroutine
- `set_chain` — subroutine
- `set_coord` — subroutine
- `set_element` — subroutine
- `set_name` — subroutine
- `set_num` — subroutine
- `set_occupancy` — subroutine
- `set_resname` — subroutine
- `set_resnum` — subroutine
- `translate` — subroutine
- `writepdb` — subroutine
- `writepdb_aniso` — subroutine
- `Z_and_radius_from_name` — subroutine

Local symbols (enclosing procedure scope):
- `convolve::epot` — function
- `element_exists::exists` — function
- `new_from_cif::element_from_name` — function
- `new_from_cif::find_col` — function
- `new_from_cif::scan_cif` — subroutine
- `new_from_cif::split_tokens` — subroutine
- `new_from_pdb::is_valid_entry` — function
- `new_from_pdb::strip_digits` — function
- `parse_element_list::capitalised` — function
- `parse_element_list::physical` — function

---
## Module: simple_atoms_tester

Files:
- `main/nano/simple_atoms_tester.f90`

Uses:
- `simple_atoms`
- `simple_cmdline`
- `simple_defs`
- `simple_defs_atoms`
- `simple_forked_process`
- `simple_image`
- `simple_nanoparticle_utils`
- `simple_parameters`
- `simple_rnd`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`
- `simple_tseries_graphene_subtr`

Public symbols:
- `run_all_atoms_tests` — subroutine

Private symbols:
- `lattice_lookup` — type
- `lookup_missing_lattice` — subroutine
- `make_model` — subroutine
- `make_pseudo` — subroutine
- `test_access` — subroutine
- `test_anisou` — subroutine
- `test_compound_selector` — subroutine
- `test_element_bfac` — subroutine
- `test_geometry` — subroutine
- `test_graphene_keeps_process` — subroutine
- `test_lattice_cutoff` — subroutine
- `test_lattice_params_missing` — subroutine
- `test_parse_element_list` — subroutine
- `test_pdb_hetatm` — subroutine
- `test_pdb_roundtrip` — subroutine
- `test_pseudo_density` — subroutine
- `test_pseudo_halo` — subroutine
- `test_pseudo_pdb_roundtrip` — subroutine
- `test_pseudo_symbols` — subroutine
- `test_validation` — subroutine

---
## Module: simple_bench_io

Files:
- `fileio/simple_bench_io.f90`

Uses:
- `simple_fileio`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `bench_report_fname` — function
- `open_bench_report` — subroutine

---
## Module: simple_binoris

Files:
- `fileio/simple_binoris.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_defs_ori`
- `simple_fileio`
- `simple_map_reduce`
- `simple_ori`
- `simple_oris`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `binoris` — type
- `binoris_seginfo` — type

Private symbols:
- `add_segment_oris` — subroutine
- `add_segment_strs` — subroutine
- `byte_manager4seg_inside_restore` — subroutine
- `byte_manager4seg_inside_save` — subroutine
- `clear_segments` — subroutine
- `close` — subroutine
- `empty_segment_inside` — subroutine
- `get_first_data_byte` — function
- `get_fromto` — function
- `get_n_bytes_per_record` — function
- `get_n_bytes_tot` — function
- `get_n_records` — function
- `get_n_segments` — function
- `get_segments_info` — subroutine
- `is_opened` — function
- `is_particle_seg` — function
- `open` — subroutine
- `read_first_segment_record` — subroutine
- `read_header` — subroutine
- `read_particle_record` — subroutine
- `read_record` — subroutine
- `read_segment_oris` — subroutine
- `read_segment_strs` — subroutine
- `update_byte_ranges` — subroutine
- `write_header` — subroutine
- `write_segment_inside_oris` — subroutine
- `write_segment_inside_strs` — subroutine
- `write_segment_oris` — subroutine
- `write_segment_strs` — subroutine

Local symbols (enclosing procedure scope):
- `open::open_local` — subroutine

---
## Module: simple_binoris_io

Files:
- `fileio/simple_binoris_io.f90`

Uses:
- `simple_core_module_api`
- `simple_sp_project`

Public symbols:
- `binread_ctfparams_state_eo` — subroutine
- `binread_nlines` — function
- `binread_oritab` — subroutine
- `binwrite_oritab` — subroutine

---
## Module: simple_binoris_tester

Files:
- `main/project/simple_binoris_tester.f90`

Uses:
- `simple_binoris`
- `simple_binoris_io`
- `simple_defs_ori`
- `simple_fileio`
- `simple_ori`
- `simple_oris`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_binoris_tests` — subroutine

Private symbols:
- `assert_ptcl_equal` — subroutine
- `assert_records_equal` — subroutine
- `cleanup` — subroutine
- `make_mic_oris` — subroutine
- `make_project` — subroutine
- `make_ptcl_oris` — subroutine
- `make_stk_oris` — subroutine
- `read_three` — subroutine
- `test_binoris_io_project_dispatch` — subroutine
- `test_binoris_io_text_dispatch` — subroutine
- `test_cartesian_record_slots` — subroutine
- `test_earlier_release_particle_records` — subroutine
- `test_open_close_and_empty_header` — subroutine
- `test_particle_segment_roundtrip` — subroutine
- `test_project_clear_segment_inside` — subroutine
- `test_project_write_segment_inside` — subroutine
- `test_string_segment_roundtrip` — subroutine
- `test_write_segment_inside` — subroutine
- `test_write_segment_inside_strings` — subroutine

---
## Module: simple_bspline_smoother

Files:
- `utils/filter/simple_bspline_smoother.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `bspline_smoother` — type

Private symbols:
- `a0xy` — function
- `a2xy` — function
- `b3` — function
- `fill_b` — subroutine
- `fill_b_3d` — subroutine
- `fill_r` — subroutine
- `fill_r_3d` — subroutine
- `kill_bspline_smoother` — subroutine
- `new_bspline_smoother` — subroutine
- `smooth` — subroutine
- `smooth_3d` — subroutine

Local symbols (enclosing procedure scope):
- `b3::b3_help` — function

---
## Module: simple_bspline_smoother_tester

Files:
- `utils/filter/simple_bspline_smoother_tester.f90`

Uses:
- `simple_bspline_smoother`
- `simple_defs`
- `simple_image`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_bspline_smoother_tests` — subroutine

Private symbols:
- `a0` — function
- `a2` — function
- `beta` — function
- `make_cosine` — subroutine
- `test_ft_state` — subroutine
- `test_transfer_2d` — subroutine
- `test_transfer_3d` — subroutine
- `transfer` — function

---
## Module: simple_builder

Files:
- `main/simple_builder.f90`

Uses:
- `simple_binoris_io`
- `simple_cartft_calc`
- `simple_class_frcs`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_euclid_sigma2`
- `simple_image`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_projector`
- `simple_sp_project`
- `simple_srchspace_map`
- `simple_ui`

Public symbols:
- `builder` — type

Private symbols:
- `build_eulspace` — subroutine
- `build_general_tbox` — subroutine
- `build_spproj` — subroutine
- `build_strategy2D_tbox` — subroutine
- `build_strategy3D_tbox` — subroutine
- `init_params_and_build_general_tbox` — subroutine
- `init_params_and_build_spproj` — subroutine
- `init_params_and_build_strategy2D_tbox` — subroutine
- `init_params_and_build_strategy3D_tbox` — subroutine
- `init_params_spproj_tbox2D` — subroutine
- `kill_general_tbox` — subroutine
- `kill_strategy2D_tbox` — subroutine
- `kill_strategy3D_tbox` — subroutine
- `new_vols` — subroutine
- `note_obs_mask` — subroutine
- `obs_on_box` — function

Local symbols (enclosing procedure scope):
- `build_spproj::set_field_ptr` — subroutine

---
## Module: simple_butterworth

Files:
- `utils/filter/simple_butterworth.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `butterworth` — function
- `butterworth_filter_band` — subroutine
- `butterworth_filter_img` — subroutine
- `butterworth_filter_krange` — subroutine
- `butterworth_filter_pair` — subroutine
- `butterworth_filter_vec` — subroutine

---
## Module: simple_calc_pspec_strategy

Files:
- `main/strategies/parallelization/simple_calc_pspec_strategy.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_cmdline`
- `simple_commanders_euclid_distr`
- `simple_core_module_api`
- `simple_eo_partition`
- `simple_image`
- `simple_matcher_ptcl_io`
- `simple_parameters`
- `simple_qsys_funs`
- `simple_sigma2_pspec_sums`

Public symbols:
- `calc_pspec_inmem_strategy` — type
- `calc_pspec_partitioned_strategy` — type
- `calc_pspec_strategy` — type
- `create_calc_pspec_strategy` — function

Private symbols:
- `cleanup_calc_pspec_outputs` — subroutine
- `cleanup_interface` — subroutine
- `compute_pspec_partitions` — subroutine
- `ensure_calc_pspec_eo_partition` — subroutine
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `partitioned_cleanup` — subroutine
- `partitioned_execute` — subroutine
- `partitioned_finalize_run` — subroutine
- `partitioned_initialize` — subroutine
- `read_sigma2_bootstrap_selection` — subroutine
- `sanity_check_calc_pspec_input` — subroutine
- `write_sigma2_bootstrap_selection` — subroutine

Local symbols (enclosing procedure scope):
- `compute_pspec_partitions::accumulate_run` — subroutine
- `compute_pspec_partitions::sanitize_computed_sigma2` — function
- `compute_pspec_partitions::slot_of` — function

---
## Module: simple_calpha_finder

Files:
- `main/nano/simple_calpha_finder.f90`

Uses:
- `simple_atoms`
- `simple_core_module_api`
- `simple_image`
- `simple_oris`

Public symbols:
- `calpha_finder` — type

Private symbols:
- `build_rotation_grid` — subroutine
- `clear_search_border` — subroutine
- `fill_rotated_target` — subroutine
- `interpolation_cell` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `sample_target` — subroutine
- `search` — subroutine
- `suppress_neighborhood` — subroutine
- `trilinear_value` — function
- `update_scores` — subroutine
- `write_candidate_files` — subroutine
- `write_candidates` — subroutine
- `zero_score_border` — subroutine

---
## Module: simple_calpha_finder_tester

Files:
- `main/nano/simple_calpha_finder_tester.f90`

Uses:
- `simple_atoms`
- `simple_calpha_finder`
- `simple_core_module_api`
- `simple_image`
- `simple_test_utils`

Public symbols:
- `run_all_calpha_finder_tests` — subroutine

Private symbols:
- `test_three_residues` — subroutine

---
## Module: simple_cartesian_fourier

Files:
- `main/interp/simple_cartesian_fourier.f90`

Uses:
- `simple_defs`
- `simple_image`

Public symbols:
- `bspline3_stencil_envelope_1d` — subroutine
- `center_crop_real3d` — function
- `center_embed_real3d` — function
- `extract_native_fourier_plane` — function
- `gather_packed_bspline3_window_grad` — subroutine
- `gather_packed_window` — function
- `gather_packed_window_grad` — subroutine

---
## Module: simple_cartesian_fourier_tester

Files:
- `main/interp/simple_cartesian_fourier_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_cartesian_fourier`
- `simple_core_module_api`
- `simple_defs`
- `simple_gridding`
- `simple_image`
- `simple_kbinterpol`
- `simple_ori`
- `simple_reconstructor_pcg`
- `simple_test_utils`

Public symbols:
- `run_all_cartesian_fourier_tests` — subroutine

Private symbols:
- `build_truth_volume` — subroutine
- `complex_is_finite` — function
- `complex_vector_is_finite` — function
- `gather_synthetic_field` — function
- `ideal_kb_reference` — subroutine
- `legacy_center_crop` — function
- `legacy_center_embed` — function
- `legacy_native_plane` — function
- `legacy_packed_gather` — subroutine
- `relative_complex_error` — function
- `run_kb_derivative` — subroutine
- `run_neutral_extract` — subroutine
- `run_packed_gather_derivative` — subroutine
- `sample_neutral` — subroutine
- `test_centered_crop_envelope` — subroutine
- `test_centered_embed_crop` — subroutine
- `test_fast_polynomial` — subroutine
- `test_fused_bspline3_gather_oracle` — subroutine
- `test_native_plane_extraction` — subroutine
- `test_normalized_stencil` — subroutine
- `test_packed_gathers` — subroutine
- `test_stencil_switch` — subroutine

---
## Module: simple_cartft_calc

Files:
- `main/cftc/simple_cartft_calc.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_cartesian_fourier`
- `simple_core_module_api`
- `simple_ctf`
- `simple_gridding`
- `simple_image`

Public symbols:
- `cartft_calc` — type

Private symbols:
- `add_dc_pixel` — subroutine
- `cart_refvols_header_compatible` — function
- `cartft_ptcl` — type
- `cartft_ref` — type
- `check_eval_args` — subroutine
- `ctf_value` — function
- `curvature_scale` — function
- `fill_ptcl` — subroutine
- `finalize_cc_normal_terms` — subroutine
- `finalize_euclid_normal_terms` — subroutine
- `gather` — subroutine
- `get_box` — function
- `get_nptcls` — function
- `get_ptcl_kfromto` — function
- `get_ptcl_npix` — function
- `get_ptcl_objfun` — function
- `get_ptcl_slot` — function
- `get_ptcl_taper` — function
- `get_ref_bytes` — function
- `has_refs` — function
- `kill` — subroutine
- `new` — subroutine
- `new_ptcls` — subroutine
- `objective_gradient` — subroutine
- `pose_normal_terms` — subroutine
- `pose_normal_terms_full_disk` — subroutine
- `predict` — subroutine
- `ptcl_is_valid` — function
- `read` — subroutine
- `ref_exists` — function
- `ref_index` — function
- `score` — function
- `set_ptcl_cc` — subroutine
- `set_ptcl_euclid` — subroutine
- `set_ptcl_inds` — subroutine
- `set_ref` — subroutine
- `set_refvol` — subroutine
- `set_sigma2` — subroutine
- `shift_normal_terms` — subroutine
- `shift_normal_terms_full_disk` — subroutine
- `shift_phase_tables` — subroutine
- `sigma_contribution` — subroutine
- `sigma_from_shell_residual` — subroutine
- `write` — subroutine

---
## Module: simple_cartft_calc_tester

Files:
- `main/cftc/simple_cartft_calc_tester.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_builder`
- `simple_cartft_calc`
- `simple_cartft_pose_opt`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_ctf`
- `simple_defs`
- `simple_euclid_sigma2`
- `simple_image`
- `simple_matcher_2dprep`
- `simple_matcher_ptcl_batch`
- `simple_matcher_refvol_utils`
- `simple_memoize_ft_maps`
- `simple_oris`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_projector`
- `simple_refine3d_fnames`
- `simple_sigma2_state`
- `simple_sigma2_state_file`
- `simple_sigma2_state_tester`
- `simple_sp_project`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_cartft_calc_tests` — subroutine

Private symbols:
- `build_batch_fixture` — subroutine
- `build_raw_particle` — subroutine
- `build_small_volume` — subroutine
- `build_test_volume` — subroutine
- `build_two_state_volumes` — subroutine
- `identity_rotation` — function
- `new_sigma2_fixture` — subroutine
- `new_test_calc` — subroutine
- `plane_to_image` — subroutine
- `remove_sigma2_files` — subroutine
- `same_predictions` — function
- `set_euclid_slot` — subroutine
- `set_unweighted_particle` — subroutine
- `slot_signature` — subroutine
- `slot_signatures` — subroutine
- `test_batch_preparation` — subroutine
- `test_canonical_sigma_round_trip` — subroutine
- `test_cc_objective_formula` — subroutine
- `test_concurrent_evaluation` — subroutine
- `test_curvature_scale` — subroutine
- `test_euclid_objective_formula` — subroutine
- `test_even_odd_references` — subroutine
- `test_five_parameter_gradient` — subroutine
- `test_gauss_newton_oracle` — subroutine
- `test_half_plane_kernels` — subroutine
- `test_matched_projector_boundary` — subroutine
- `test_observation_established_path` — subroutine
- `test_observation_preparation` — subroutine
- `test_particle_preparation_contract` — subroutine
- `test_polar_grid_pose_agreement` — subroutine
- `test_polar_sigma_agreement` — subroutine
- `test_polar_sigma_group` — subroutine
- `test_polar_sigma_signal` — subroutine
- `test_polish_sigma_fallback` — subroutine
- `test_reference_lifecycle` — subroutine
- `test_reference_volume_file` — subroutine
- `test_reference_volume_handoff` — subroutine
- `test_shift_phase_sign` — subroutine
- `test_sigma2_invalidates_slots` — subroutine
- `test_sigma_endpoint_contracts` — subroutine
- `test_sigma_from_shell_residual` — subroutine
- `test_sigma_owner` — subroutine
- `test_sigma_shell_contract` — subroutine
- `whitened_power` — function

Local symbols (enclosing procedure scope):
- `test_batch_preparation::new_batch_calc` — subroutine
- `test_concurrent_evaluation::evaluate_all` — subroutine
- `test_curvature_scale::normal_blocks` — subroutine
- `test_gauss_newton_oracle::objective_at` — function
- `test_matched_projector_boundary::compare` — subroutine
- `test_polar_grid_pose_agreement::grid_rotation` — function
- `test_polar_sigma_group::make_signal` — subroutine

---
## Module: simple_cartft_pose_opt

Files:
- `main/cftc/simple_cartft_pose_opt.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_cartft_calc`
- `simple_core_module_api`

Public symbols:
- `cartft_pose_opt` — type
- `cartft_pose_work` — type
- `right_increment_rotation` — function

Private symbols:
- `aborts_route` — function
- `add` — subroutine
- `add_solve` — subroutine
- `apply_pose_parameter_mask` — subroutine
- `build_pose_lm_system` — subroutine
- `cartft_stage` — type
- `get_motion` — subroutine
- `get_objectives` — subroutine
- `get_shell_residual` — subroutine
- `get_stage` — subroutine
- `get_status` — function
- `increase_lm_damping` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `positive_finite` — function
- `ratio_bin` — function
- `refine` — subroutine
- `refine_joint` — subroutine
- `refine_pose` — subroutine
- `refine_shift` — subroutine
- `rejection_terminal_status` — function
- `reset_record` — subroutine
- `solve_pose_cholesky` — subroutine
- `write_bench` — subroutine

Local symbols (enclosing procedure scope):
- `write_bench::put_count` — subroutine
- `write_bench::put_real` — subroutine

---
## Module: simple_cartft_pose_opt_tester

Files:
- `main/cftc/simple_cartft_pose_opt_tester.f90`

Uses:
- `simple_cartft_calc`
- `simple_cartft_pose_opt`
- `simple_core_module_api`
- `simple_defs`
- `simple_defs_ori`
- `simple_ori`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_cartft_pose_opt_tests` — subroutine

Private symbols:
- `build_small_volume` — subroutine
- `build_test_volume` — subroutine
- `determinant3` — function
- `fixture_record` — type
- `get_stage_counts` — subroutine
- `identity_rotation` — function
- `make_seed_ori` — subroutine
- `new_fixture` — subroutine
- `opt_stage_status` — function
- `rotation_distance` — function
- `set_centred_particle` — subroutine
- `test_committed_shell_residual` — subroutine
- `test_invalid_and_unobservable_inputs` — subroutine
- `test_joint_solver` — subroutine
- `test_ncc_solver` — subroutine
- `test_objective_recovery` — subroutine
- `test_ori_round_trip` — subroutine
- `test_parameter_bounds` — subroutine
- `test_rollback_transaction_contracts` — subroutine
- `test_rotation_increment` — subroutine
- `test_route_transaction_contracts` — subroutine
- `test_shift_solver` — subroutine
- `test_tiny_accepted_reduction_stop` — subroutine

---
## Module: simple_cavg_quality_analysis

Files:
- `main/cavg_quality/simple_cavg_quality_analysis.f90`

Uses:
- `simple_cavg_quality_feats`
- `simple_cavg_quality_lattice`
- `simple_cavg_quality_model`
- `simple_cavg_quality_relations`
- `simple_cavg_quality_stats`
- `simple_cavg_quality_types`
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_oris`
- `simple_parameters`
- `simple_stat`

Public symbols:
- `evaluate_cavg_quality` — subroutine
- `evaluate_cavg_quality_for_analysis` — subroutine
- `evaluate_cavg_quality_hard_reject` — subroutine
- `prepare_cavg_quality` — subroutine
- `write_cavg_quality_feature_table` — subroutine
- `write_cavg_quality_training_table` — subroutine

Private symbols:
- `calc_best_thresholds` — subroutine
- `evaluate_cavg_quality_in_context` — subroutine
- `write_analysis_class_header` — subroutine
- `write_analysis_class_row` — subroutine
- `write_feature_table_class_row` — subroutine
- `write_lattice_comments` — subroutine
- `write_model_as_analysis_comments` — subroutine
- `write_threshold_scan_comments` — subroutine

Local symbols (enclosing procedure scope):
- `calc_best_thresholds::eval_threshold` — subroutine
- `write_threshold_scan_comments::write_one_threshold` — subroutine

---
## Module: simple_cavg_quality_feats

Files:
- `main/cavg_quality/simple_cavg_quality_feats.f90`

Uses:
- `simple_cavg_quality_lattice`
- `simple_cavg_quality_types`
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_image_bin`
- `simple_oris`
- `simple_segmentation`
- `simple_stat`

Public symbols:
- `cavg_quality_bank_size` — function
- `cavg_quality_feature_description` — function
- `cavg_quality_feature_direction` — function
- `cavg_quality_feature_family` — function
- `cavg_quality_feature_name` — function
- `extract_cavg_quality_features` — subroutine
- `normalize_cavg_quality_features` — subroutine
- `validate_quality_context` — subroutine
- `write_cavg_quality_feature_inventory` — subroutine

Private symbols:
- `cavg_hard_reject_for_context` — function
- `feature_def` — function
- `measure_cavg_foreground_geometry` — subroutine
- `measure_cavg_image_metrics` — subroutine
- `measure_nano_locvar` — subroutine
- `nano_hard_reject` — function
- `normalize_robust` — subroutine
- `resolution_feature` — function

---
## Module: simple_cavg_quality_helpers

Files:
- `main/cavg_quality/simple_cavg_quality_helpers.f90`

Uses:
- `simple_cavg_quality_types`
- `simple_string`

Public symbols:
- `cavg_rejection_reason_string` — function

---
## Module: simple_cavg_quality_lattice

Files:
- `main/cavg_quality/simple_cavg_quality_lattice.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_math_ft`
- `simple_stat`

Public symbols:
- `cavg_quality_lattice` — type
- `model_gate_states` — subroutine

Private symbols:
- `extract` — subroutine
- `insert_desc` — subroutine
- `is_local_max` — function
- `kill` — subroutine
- `measure` — subroutine
- `new` — subroutine
- `soft_edge` — function
- `sphere2` — function
- `wrap` — function

---
## Module: simple_cavg_quality_lattice_tester

Files:
- `main/cavg_quality/simple_cavg_quality_lattice_tester.f90`

Uses:
- `simple_cavg_quality_analysis`
- `simple_cavg_quality_feats`
- `simple_cavg_quality_lattice`
- `simple_cavg_quality_learn`
- `simple_cavg_quality_model`
- `simple_cavg_quality_relations`
- `simple_cavg_quality_types`
- `simple_cmdline`
- `simple_defs`
- `simple_fileio`
- `simple_forked_process`
- `simple_image`
- `simple_oris`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_cavg_quality_lattice_tests` — subroutine

Private symbols:
- `class_oris` — subroutine
- `drop_last_value` — subroutine
- `gauss` — function
- `gauss_vec` — function
- `kill_imgs` — subroutine
- `make_cavg` — subroutine
- `model_reader` — type
- `next_uniform` — function
- `project_quality` — subroutine
- `read_short_model` — subroutine
- `test_band_presence_unmasked` — subroutine
- `test_bank_by_context` — subroutine
- `test_diameter_calibration` — subroutine
- `test_lattice_against_noise` — subroutine
- `test_model_file_round_trips` — subroutine
- `test_model_gate` — subroutine
- `test_nano_gates` — subroutine
- `test_nano_preset` — subroutine
- `test_normalization` — subroutine
- `test_spot_arc_straight_bent` — subroutine
- `test_table_round_trip_learn` — subroutine
- `test_team_identity` — subroutine
- `write_synthetic_table` — subroutine

Local symbols (enclosing procedure scope):
- `test_model_file_round_trips::check_short_list` — subroutine
- `test_nano_preset::gate_is_decision` — function
- `test_nano_preset::m_probability` — function
- `test_normalization::median_of` — function

---
## Module: simple_cavg_quality_learn

Files:
- `main/cavg_quality/simple_cavg_quality_learn.f90`

Uses:
- `simple_cavg_quality_feats`
- `simple_cavg_quality_model`
- `simple_cavg_quality_relations`
- `simple_cavg_quality_stats`
- `simple_cavg_quality_types`
- `simple_defs`
- `simple_error`
- `simple_quality_logistic`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `evaluate_cavg_quality_model` — subroutine
- `evaluate_cavg_quality_result` — subroutine
- `learn_cavg_quality_model` — subroutine

Private symbols:
- `append_feature_family` — subroutine
- `bad_overfit_signature` — function
- `balanced_fn_excess_rate` — function
- `build_logistic_problem` — subroutine
- `classify_training_dataset_detail` — subroutine
- `collect_learn_diagnostics` — subroutine
- `count_accepted_bad_overfit_signature` — function
- `count_hard_rejected_manual_good` — function
- `count_trainable_bad_overfit_signature` — function
- `count_trainable_classes` — function
- `count_trainable_manual_bad` — function
- `count_trainable_manual_good` — function
- `dataset_is_overfit_focus` — function
- `dataset_learn_role` — function
- `dataset_learn_role_name` — function
- `feature_policy_indices` — subroutine
- `feature_policy_mask` — subroutine
- `feature_policy_name` — function
- `fn_tolerant_specificity_score` — function
- `guarded_recall_score` — function
- `infer_training_context` — subroutine
- `is_analysis_data_line` — function
- `is_analysis_header_line` — function
- `kill_training_datasets` — subroutine
- `learn_balacc_from_confusion` — function
- `learn_cavg_quality_pairwise_logistic_model` — subroutine
- `load_quality_training_datasets` — subroutine
- `logistic_bad_example_loss_multiplier` — function
- `logistic_overfit_focus_bad_loss_multiplier` — function
- `logistic_solution_to_model` — subroutine
- `macro_balacc_for_model` — function
- `map_analysis_columns` — subroutine
- `n_feature_policies` — function
- `neutral_learn_base_spec` — function
- `overfit_false_positive_penalty` — function
- `overfit_signature_applies` — function
- `parse_training_metadata_line` — subroutine
- `read_quality_training_dataset` — subroutine
- `require_analysis_columns` — subroutine
- `require_relational_training_datasets` — subroutine
- `require_trainable_rows` — subroutine
- `validate_quality_context` — subroutine
- `validate_training_dataset_schemas` — subroutine
- `write_cavg_quality_evaluate_report` — subroutine
- `write_cavg_quality_logistic_learn_report` — subroutine
- `write_dataset_metric_table` — subroutine
- `write_evaluate_diagnostic` — subroutine
- `write_evaluate_diagnostics` — subroutine
- `write_feature_policy_grid` — subroutine
- `write_fixed_model_summary` — subroutine
- `write_logistic_coefficient_table` — subroutine
- `write_real_list` — subroutine

Local symbols (enclosing procedure scope):
- `build_logistic_problem::set_row` — subroutine

---
## Module: simple_cavg_quality_model

Files:
- `main/cavg_quality/simple_cavg_quality_model.f90`

Uses:
- `simple_cavg_quality_feats`
- `simple_cavg_quality_lattice`
- `simple_cavg_quality_types`
- `simple_defs`
- `simple_error`
- `simple_string_utils`

Public symbols:
- `cavg_quality_model` — type
- `write_cavg_quality_model_builtin_code` — subroutine

Private symbols:
- `apply_pairwise_logistic` — subroutine
- `builtin_names` — function
- `builtin_spec` — function
- `chunk100mics_model_spec` — function
- `classify` — subroutine
- `feature_index_from_token` — function
- `get_spec` — function
- `init_preset` — subroutine
- `init_spec` — subroutine
- `kill_model` — subroutine
- `nano_model_spec` — function
- `normalize` — subroutine
- `pairwise_logistic_probability` — function
- `parse_feature_weight_values` — subroutine
- `parse_model_key_value` — subroutine
- `parse_real_values` — subroutine
- `pool_model_spec` — function
- `read_feature_weights` — subroutine
- `read_interaction_terms` — subroutine
- `read_model` — subroutine
- `read_real_values_keyed` — subroutine
- `set_nano_threshold` — subroutine
- `set_pairwise_interactions_for_feature_count` — subroutine
- `sieve_model_spec` — function
- `supports_relational` — function
- `write_interaction_terms` — subroutine
- `write_logistic_spec_assignments` — subroutine
- `write_model` — subroutine
- `write_model_real_list` — subroutine
- `write_model_spec_function` — subroutine
- `write_real_array_assignment` — subroutine
- `write_weights_assignment` — subroutine

---
## Module: simple_cavg_quality_relations

Files:
- `main/cavg_quality/simple_cavg_quality_relations.f90`

Uses:
- `simple_cavg_quality_types`
- `simple_error`
- `simple_image`
- `simple_imgarr_utils`
- `simple_parameters`
- `simple_srch_sort_loc`
- `simple_stat`
- `simple_strategy2d_utils`
- `simple_type_defs`

Public symbols:
- `calculate_promoted_feature` — subroutine
- `cavg_quality_relation_analysis` — type

Private symbols:
- `calculate_relation_analysis` — subroutine
- `kill_relation_analysis` — subroutine
- `normalize_relation_feature` — subroutine
- `sorted_cc_neighbours` — subroutine

---
## Module: simple_cavg_quality_relations_tester

Files:
- `main/cavg_quality/simple_cavg_quality_relations_tester.f90`

Uses:
- `simple_cavg_quality_relations`
- `simple_cavg_quality_types`
- `simple_test_utils`

Public symbols:
- `run_all_cavg_quality_relations_tests` — subroutine

Private symbols:
- `test_promoted_feature` — subroutine

---
## Module: simple_cavg_quality_selection

Files:
- `main/cavg_quality/simple_cavg_quality_selection.f90`

Uses:
- `simple_cavg_quality_analysis`
- `simple_cavg_quality_model`
- `simple_cavg_quality_types`
- `simple_cmdline`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_image`
- `simple_imgarr_utils`
- `simple_parameters`
- `simple_rnd`
- `simple_sp_project`
- `simple_string`

Public symbols:
- `score_project_cavgs` — subroutine
- `write_cavg_selection_stacks` — subroutine
- `write_cavg_stack` — subroutine

---
## Module: simple_cavg_quality_selection_tester

Files:
- `main/cavg_quality/simple_cavg_quality_selection_tester.f90`

Uses:
- `simple_cavg_quality_selection`
- `simple_fileio`
- `simple_image`
- `simple_imghead`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_cavg_quality_selection_tests` — subroutine

Private symbols:
- `check_stack` — subroutine
- `cleanup` — subroutine
- `make_imgs` — subroutine
- `stack_names` — subroutine
- `test_selection_stacks_replace_existing` — subroutine
- `test_selection_stacks_split_and_order` — subroutine

---
## Module: simple_cavg_quality_stats

Files:
- `main/cavg_quality/simple_cavg_quality_stats.f90`

Uses:
- `simple_cavg_quality_types`
- `simple_defs`
- `simple_error`
- `simple_srch_sort_loc`
- `simple_stat`

Public symbols:
- `auc_for_values` — function
- `calc_binary_metrics` — subroutine
- `calc_confusion` — subroutine
- `mad_by_state` — function
- `median_by_state` — function
- `safe_div` — function

Local symbols (enclosing procedure scope):
- `auc_for_values::count_below` — function
- `auc_for_values::count_within_eps` — function

---
## Module: simple_cavg_quality_types

Files:
- `main/cavg_quality/simple_cavg_quality_types.f90`

Uses:
- `simple_defs`
- `simple_string`

Public symbols:
- `cavg_quality_feature_def` — type
- `cavg_quality_learn_diagnostics` — type
- `cavg_quality_model_spec` — type
- `cavg_quality_result` — type
- `cavg_quality_training_dataset` — type
- `reset_cavg_quality_result` — subroutine

---
## Module: simple_cavg_registration_tester

Files:
- `main/strategies/search/simple_cavg_registration_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_image`
- `simple_imgarr_utils`
- `simple_parameters`
- `simple_strategy2d_utils`
- `simple_test_utils`

Public symbols:
- `run_all_cavg_registration_tests` — subroutine

Private symbols:
- `test_rotated_copies` — subroutine

Local symbols (enclosing procedure scope):
- `test_rotated_copies::angdist` — function

---
## Module: simple_cavg_sums

Files:
- `main/class/simple_cavg_sums.f90`

Uses:
- `simple_core_module_api`
- `simple_ftiter`

Public symbols:
- `cavg_contrib_fname` — function
- `cavg_sums` — type

Private symbols:
- `accumulate` — subroutine
- `blend` — subroutine
- `check_shape` — subroutine
- `class_mass` — subroutine
- `get_eo_pops` — subroutine
- `get_l_frac` — function
- `get_mrep` — subroutine
- `get_offsets` — subroutine
- `get_part` — function
- `get_sums` — subroutine
- `header_bytes` — function
- `kill` — subroutine
- `matches` — function
- `new` — subroutine
- `pad_to` — subroutine
- `payload_bytes` — function
- `read` — subroutine
- `set_contrib_meta` — subroutine
- `set_mrep` — subroutine
- `set_sums` — subroutine
- `shift_class` — subroutine
- `write` — subroutine
- `zero` — subroutine

---
## Module: simple_cavg_sums_tester

Files:
- `main/class/simple_cavg_sums_tester.f90`

Uses:
- `simple_cavg_sums`
- `simple_core_module_api`
- `simple_test_utils`

Public symbols:
- `run_all_cavg_sums_tests` — subroutine

Private symbols:
- `fill` — subroutine
- `maxdiff` — function
- `test_cavg_sums_pad` — subroutine
- `test_cavg_sums_refusals` — subroutine
- `test_cavg_sums_roundtrip` — subroutine
- `test_owner_blend_split` — subroutine

Local symbols (enclosing procedure scope):
- `fill::next` — function
- `test_owner_blend_split::shifted_prev` — subroutine
- `test_owner_blend_split::zero_c` — function
- `test_owner_blend_split::zero_r` — function

---
## Module: simple_cavgassemble_strategy

Files:
- `main/strategies/parallelization/simple_cavgassemble_strategy.f90`

Uses:
- `simple_bench_io`
- `simple_builder`
- `simple_classaverager`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_parameters`
- `simple_refine2d_fnames`
- `simple_starproject`
- `simple_stream_utils`

Public symbols:
- `cavgassemble_exec` — subroutine
- `run_cavgassemble` — subroutine

---
## Module: simple_cavgs_state_match

Files:
- `main/class/simple_cavgs_state_match.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_imgarr_utils`
- `simple_math`
- `simple_oris`
- `simple_parameters`
- `simple_simple_volinterp`
- `simple_strategy2d_utils`
- `simple_type_defs`

Public symbols:
- `match_cavgs_to_states` — subroutine
- `matched_pose` — function
- `reduce_state_matches` — subroutine

Local symbols (enclosing procedure scope):
- `match_cavgs_to_states::resample` — subroutine

---
## Module: simple_cavgs_state_match_tester

Files:
- `main/class/simple_cavgs_state_match_tester.f90`

Uses:
- `simple_cavgs_state_match`
- `simple_image`
- `simple_imgarr_utils`
- `simple_linalg`
- `simple_ori_utils`
- `simple_oris`
- `simple_simple_volinterp`
- `simple_test_utils`

Public symbols:
- `run_all_cavgs_state_match_tests` — subroutine

Private symbols:
- `blob_volume` — subroutine
- `test_matched_pose` — subroutine
- `test_reduce_state_matches` — subroutine
- `test_two_volumes` — subroutine

---
## Module: simple_cavgs_work_project

Files:
- `main/simple_cavgs_work_project.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_linalg`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_sym`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `cavgs_work_project` — type

Private symbols:
- `cap_nparts` — subroutine
- `delete_files` — subroutine
- `eo_angular_distance` — function
- `eo_state_overlap` — function
- `even_row` — function
- `get_ncavgs` — function
- `get_projfile` — function
- `kill` — subroutine
- `label_rows` — subroutine
- `new` — subroutine
- `odd_row` — function

---
## Module: simple_cavgs_work_project_tester

Files:
- `main/simple_cavgs_work_project_tester.f90`

Uses:
- `simple_cavgs_work_project`
- `simple_cmdline`
- `simple_fileio`
- `simple_image`
- `simple_oris`
- `simple_parameters`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_cavgs_work_project_tests` — subroutine

Private symbols:
- `make_input` — subroutine
- `test_build` — subroutine
- `test_cap_nparts` — subroutine
- `test_eo_diagnostics` — subroutine
- `test_memberless_halves` — subroutine

---
## Module: simple_chash

Files:
- `utils/structs/simple_chash.f90`

Uses:
- `simple_ansi_ctrls`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `chash` — type

Private symbols:
- `alloc_chash` — subroutine
- `chash2str` — function
- `chash_strlen` — function
- `constructor_by_sz` — function
- `constructor_default` — function
- `copy` — subroutine
- `dealloc_chash` — subroutine
- `delete` — subroutine
- `gen_job_descr` — subroutine
- `get_by_ind` — function
- `get_by_key` — function
- `get_key` — function
- `get_static` — subroutine
- `isthere` — function
- `kill` — subroutine
- `lookup` — function
- `move_key_to_front_char` — subroutine
- `move_key_to_front_str` — subroutine
- `new_by_sz` — subroutine
- `new_default` — subroutine
- `parse_cmdline` — subroutine
- `print_key_val_pair_by_ind` — subroutine
- `print_key_val_pair_by_key` — subroutine
- `print_key_val_pairs` — subroutine
- `push_char` — subroutine
- `push_char_str` — subroutine
- `push_str` — subroutine
- `realloc_chash` — subroutine
- `reverselookup` — function
- `set_char` — subroutine
- `set_char_str` — subroutine
- `set_str` — subroutine
- `size_of` — function
- `sort` — subroutine

---
## Module: simple_chash_tester

Files:
- `utils/structs/simple_chash_tester.f90`

Uses:
- `simple_chash`
- `simple_defs`
- `simple_fileio`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_chash_tests` — subroutine

Private symbols:
- `test_chash2str_and_strlen` — subroutine
- `test_construction_and_kill` — subroutine
- `test_delete` — subroutine
- `test_get_by_key_and_index` — subroutine
- `test_isthere_lookup_reverse` — subroutine
- `test_move_key_to_front` — subroutine
- `test_push_and_set_basic` — subroutine
- `test_set_replaces_existing` — subroutine
- `test_sort` — subroutine

---
## Module: simple_cif2mrc_tester

Files:
- `main/nano/simple_cif2mrc_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_atoms`
- `simple_defs`
- `simple_image`
- `simple_imghead`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_cif2mrc_tests` — subroutine

Private symbols:
- `build_analytical_reference` — subroutine
- `expected_position` — function
- `fixture_position` — function
- `reference_coefficients` — subroutine
- `write_all_element_fixture` — subroutine

---
## Module: simple_class_compatibility

Files:
- `main/class/simple_class_compatibility.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_image_bin`
- `simple_imgarr_utils`
- `simple_segmentation`
- `simple_sp_project`
- `simple_srch_sort_loc`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `class_compatibility` — type
- `support_model_metrics` — type

Private symbols:
- `apply_support_model` — subroutine
- `converged` — function
- `get_support_model_metrics` — subroutine
- `infer` — subroutine
- `kill` — subroutine
- `kill_support_model` — subroutine
- `new` — subroutine
- `new_support_training_set` — subroutine
- `preprocess_img` — subroutine
- `support_model` — type
- `support_model_training_set` — type
- `train_1` — subroutine
- `train_2` — subroutine
- `update_support_model` — subroutine

Local symbols (enclosing procedure scope):
- `update_support_model::percentile_real` — function
- `update_support_model::percentile_weighted_real` — function

---
## Module: simple_class_compatibility_tester

Files:
- `main/class/simple_class_compatibility_tester.f90`

Uses:
- `simple_class_compatibility`
- `simple_core_module_api`
- `simple_image`
- `simple_imgarr_utils`
- `simple_sp_project`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_class_compatibility_tests` — subroutine

Private symbols:
- `init_project_for_cavgs` — subroutine
- `make_blob_stack` — subroutine
- `paint_square_blob` — subroutine
- `test_defaults_and_lifecycle` — subroutine
- `test_infer_rejects_size_outliers` — subroutine
- `test_infer_without_valid_model_preserves_selection` — subroutine
- `test_kill_support_model_resets_after_fit` — subroutine
- `test_train_1_from_project_sets_valid_fit` — subroutine
- `test_train_1_too_small_set_no_fit` — subroutine
- `test_train_2_retrain_sets_delta_and_convergence` — subroutine

---
## Module: simple_class_frcs

Files:
- `main/class/simple_class_frcs.f90`

Uses:
- `simple_core_module_api`
- `simple_fsc`

Public symbols:
- `class_frcs` — type
- `resample_filter` — function

Private symbols:
- `avg_frc_getter` — subroutine
- `bound_res` — subroutine
- `crop` — subroutine
- `estimate_find_for_eoavg` — function
- `estimate_lp_for_align` — function
- `estimate_res` — subroutine
- `frc_getter` — subroutine
- `get_box` — function
- `get_filtsz` — function
- `get_frc` — function
- `get_ncls` — function
- `getter` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `pad` — subroutine
- `plot_frcs` — subroutine
- `print_frcs` — subroutine
- `raise_exception` — subroutine
- `read` — subroutine
- `set_frc` — subroutine
- `write` — subroutine

---
## Module: simple_class_sample_io

Files:
- `fileio/simple_class_sample_io.f90`

Uses:
- `simple_error`
- `simple_fileio`
- `simple_string_utils`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `deallocate_class_samples` — subroutine
- `read_class_samples` — subroutine
- `write_class_samples` — subroutine

Private symbols:
- `serialize_class_sample` — function

---
## Module: simple_class_sample_io_tester

Files:
- `fileio/simple_class_sample_io_tester.f90`

Uses:
- `simple_class_sample_io`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_class_sample_io_tests` — subroutine

Private symbols:
- `assert_entry_equal` — subroutine
- `make_ragged` — subroutine
- `test_empty_class_roundtrip` — subroutine
- `test_ragged_roundtrip` — subroutine
- `test_read_replaces_previous_array` — subroutine

---
## Module: simple_classaverager

Files:
- `main/class/simple_classaverager.f90`
- `main/class/simple_classaverager_core.f90`
- `main/class/simple_classaverager_restore.f90`

Uses:
- `simple_bench_io`
- `simple_builder`
- `simple_cavg_sums`
- `simple_core_module_api`
- `simple_ctf`
- `simple_euclid_sigma2`
- `simple_fftw3`
- `simple_fplane_grid`
- `simple_ftiter`
- `simple_gridding`
- `simple_image`
- `simple_imgarr_utils`
- `simple_imgfile`
- `simple_matcher_ptcl_io`
- `simple_memoize_ft_maps`
- `simple_oris`
- `simple_parameters`
- `simple_sp_project`

Public symbols:
- `cavger_assemble_sums` — subroutine
- `cavger_assemble_sums` — subroutine
- `cavger_assemble_sums_from_parts` — subroutine
- `cavger_assemble_sums_from_parts` — subroutine
- `cavger_commit_carryover` — subroutine
- `cavger_commit_carryover` — subroutine
- `cavger_dealloc_online` — subroutine
- `cavger_dealloc_online` — subroutine
- `cavger_gen2Dclassdoc` — subroutine
- `cavger_gen2Dclassdoc` — subroutine
- `cavger_init_online` — subroutine
- `cavger_init_online` — subroutine
- `cavger_kill` — subroutine
- `cavger_kill` — subroutine
- `cavger_new` — subroutine
- `cavger_new` — subroutine
- `cavger_pad_carried_sums` — subroutine
- `cavger_pad_carried_sums` — subroutine
- `cavger_read_all` — subroutine
- `cavger_read_all` — subroutine
- `cavger_read_euclid_sigma2` — subroutine
- `cavger_read_euclid_sigma2` — subroutine
- `cavger_restore_cavgs` — subroutine
- `cavger_restore_cavgs` — subroutine
- `cavger_set_center_offset` — subroutine
- `cavger_set_center_offset` — subroutine
- `cavger_transf_oridat` — subroutine
- `cavger_transf_oridat` — subroutine
- `cavger_update_sums` — subroutine
- `cavger_write_all` — subroutine
- `cavger_write_all` — subroutine
- `cavger_write_contribution` — subroutine
- `cavger_write_contribution` — subroutine
- `cavger_write_eo` — subroutine
- `cavger_write_eo` — subroutine
- `cavger_write_merged` — subroutine
- `cavger_write_merged` — subroutine
- `cavger_write_sums_load` — subroutine
- `cavger_write_sums_load` — subroutine
- `transform_ptcls` — subroutine
- `transform_ptcls` — subroutine

Private symbols:
- `add_invnoisepower2rho` — subroutine
- `add_invnoisepower2rho` — subroutine
- `calc_cavgs_stats` — subroutine
- `calc_cavgs_stats` — subroutine
- `cavger_update_sums` — subroutine
- `cavgs2sums` — subroutine
- `cavgs_set` — type
- `commit_carryover` — subroutine
- `copy_fast` — subroutine
- `copy_fast` — subroutine
- `ctf_dens_correct` — subroutine
- `ctf_dens_correct` — subroutine
- `dealloc_cavgs` — subroutine
- `fft` — subroutine
- `fft` — subroutine
- `frc` — subroutine
- `frc` — subroutine
- `ifft` — subroutine
- `ifft` — subroutine
- `insert_lowres_serial` — subroutine
- `insert_lowres_serial` — subroutine
- `kill_set` — subroutine
- `kill_set` — subroutine
- `kill_stack` — subroutine
- `kill_stack` — subroutine
- `matrix_ptrs` — type
- `new_set` — subroutine
- `new_set` — subroutine
- `new_stack` — subroutine
- `new_stack` — subroutine
- `pad` — subroutine
- `pad` — subroutine
- `ptcl_record` — type
- `quadrant_swap` — subroutine
- `quadrant_swap` — subroutine
- `read_cavgs` — subroutine
- `read_cmat` — subroutine
- `read_cmat` — subroutine
- `read_ctfsq` — subroutine
- `read_ctfsq` — subroutine
- `softmask` — subroutine
- `softmask` — subroutine
- `stack` — type
- `stack_accumulate_fplane` — subroutine
- `stack_accumulate_fplane` — subroutine
- `sums2cavgs` — subroutine
- `write` — subroutine
- `write` — subroutine
- `write_assemble_sums_bench` — subroutine
- `write_cmat_debug` — subroutine
- `write_cmat_debug` — subroutine
- `write_ctfsq` — subroutine
- `write_ctfsq` — subroutine
- `write_ctfsq_debug` — subroutine
- `write_ctfsq_debug` — subroutine
- `zero` — subroutine
- `zero` — subroutine
- `zero_set` — subroutine
- `zero_set` — subroutine
- `zero_slice` — subroutine
- `zero_slice` — subroutine

Local symbols (enclosing procedure scope):
- `cavger_update_sums::bucket_of` — function
- `cavger_update_sums::order_tasks` — subroutine

---
## Module: simple_cls_expansion_strategy

Files:
- `main/strategies/parallelization/simple_cls_expansion_strategy.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_class_frcs`
- `simple_classaverager`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_ctf`
- `simple_estimate_ssnr`
- `simple_euclid_sigma2`
- `simple_exec_helpers`
- `simple_flex_cls_expansion`
- `simple_flex_cls_kernel`
- `simple_image`
- `simple_imgarr_utils`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_rnd`
- `simple_sigma2_files`
- `simple_sp_project`
- `simple_srch_sort_loc`
- `simple_view_partition_sampling`

Public symbols:
- `cls_expansion_master_strategy` — type
- `cls_expansion_shmem_strategy` — type
- `cls_expansion_strategy` — type
- `cls_expansion_worker_strategy` — type
- `create_cls_expansion_strategy` — function

Private symbols:
- `allocate_subclasses` — subroutine
- `apply_defaults` — subroutine
- `apply_split_project_updates` — subroutine
- `cleanup_interface` — subroutine
- `collect_split_classes` — subroutine
- `count_data_lines` — subroutine
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `flex_cavgs_part_fname` — function
- `flex_class_planes` — subroutine
- `flex_class_whole` — subroutine
- `flex_fit_band` — function
- `flex_frcs_part_fname` — function
- `flex_frcs_txt2bin` — subroutine
- `flex_planes_to_images` — subroutine
- `flex_res143` — function
- `flex_time` — subroutine
- `init_common` — subroutine
- `init_interface` — subroutine
- `master_cleanup` — subroutine
- `master_execute` — subroutine
- `master_finalize_run` — subroutine
- `master_initialize` — subroutine
- `merge_part_stacks` — subroutine
- `merge_worker_outputs` — subroutine
- `prepare_class_partitions` — subroutine
- `prepare_flex_sigma2` — subroutine
- `read_allocation` — subroutine
- `read_int_file` — subroutine
- `read_part_map` — subroutine
- `run_local_split` — subroutine
- `shmem_cleanup` — subroutine
- `shmem_execute` — subroutine
- `shmem_finalize_run` — subroutine
- `shmem_initialize` — subroutine
- `sort_combined_maps` — subroutine
- `sort_order_by_weight_desc` — subroutine
- `split_class_flex` — subroutine
- `split_label_inds` — function
- `validate_cls_expansion` — subroutine
- `worker_cleanup` — subroutine
- `worker_execute` — subroutine
- `worker_finalize_run` — subroutine
- `worker_initialize` — subroutine

Local symbols (enclosing procedure scope):
- `merge_worker_outputs::global_subclass` — function
- `run_local_split::record_class` — subroutine
- `split_class_flex::deliver_kernel` — subroutine
- `split_class_flex::fit_class_model` — subroutine
- `split_class_flex::half_mean` — subroutine

---
## Module: simple_cluster_seed

Files:
- `utils/clustering/simple_cluster_seed.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `gen_labelling` — subroutine

Private symbols:
- `draw_ranked` — subroutine
- `draw_ranked_uniform` — subroutine
- `draw_squared` — subroutine
- `draw_squared_uniform` — subroutine
- `draw_squared_uniform_projdir` — subroutine
- `draw_uniform` — subroutine
- `draw_uniform_corr` — subroutine

Local symbols (enclosing procedure scope):
- `draw_squared_uniform_projdir::bit_remove` — subroutine
- `draw_squared_uniform_projdir::bit_select` — function
- `draw_squared_uniform_projdir::draw` — subroutine

---
## Module: simple_clustering_utils

Files:
- `utils/clustering/simple_clustering_utils.f90`

Uses:
- `simple_aff_prop`
- `simple_core_module_api`
- `simple_hac`
- `simple_kmedoids`
- `simple_srch_sort_loc`
- `simple_stat`

Public symbols:
- `cluster_dmat` — subroutine
- `equal_mass_quantile_start` — subroutine
- `farthest_point_seeds` — subroutine

Local symbols (enclosing procedure scope):
- `cluster_dmat::merge_if_necessary` — subroutine

---
## Module: simple_cmdline

Files:
- `main/simple_cmdline.f90`

Uses:
- `simple_args`
- `simple_core_module_api`
- `simple_private_prgs`
- `simple_ui`
- `simple_ui_program`

Public symbols:
- `cmdline` — type
- `cmdline_err` — subroutine

Private symbols:
- `assign` — subroutine
- `check` — subroutine
- `checkvar` — subroutine
- `cmdarg` — type
- `copy` — subroutine
- `copy_arg` — subroutine
- `defined` — function
- `delete` — subroutine
- `gen_job_descr` — subroutine
- `get_argcnt` — function
- `get_carg` — function
- `get_iarg` — function
- `get_keys` — function
- `get_rarg` — function
- `kill` — subroutine
- `lookup` — function
- `parse` — subroutine
- `parse_command_line_value` — subroutine
- `parse_private` — subroutine
- `parse_private_line` — subroutine
- `printline` — subroutine
- `read` — subroutine
- `set_1` — subroutine
- `set_2` — subroutine
- `set_3` — subroutine
- `set_4` — subroutine
- `set_5` — subroutine
- `set_6` — subroutine
- `set_from_text` — subroutine
- `skip_mode_required_keys` — function

---
## Module: simple_cmdline_tester

Files:
- `main/simple_cmdline_tester.f90`

Uses:
- `simple_chash`
- `simple_cmdline`
- `simple_defs`
- `simple_string`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_cmdline_tests` — subroutine

Private symbols:
- `test_copy_and_assignment` — subroutine
- `test_defined_and_lookup` — subroutine
- `test_delete_behavior` — subroutine
- `test_gen_job_descr` — subroutine
- `test_get_keys` — subroutine
- `test_read_line_typed_values` — subroutine
- `test_read_parsing` — subroutine
- `test_set_and_get_char_and_string` — subroutine
- `test_set_and_get_numeric` — subroutine
- `test_set_from_text_and_copy_arg` — subroutine

---
## Module: simple_commander_base

Files:
- `main/commanders/simple/simple_commander_base.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs`

Public symbols:
- `commander_base` — type

Private symbols:
- `generic_execute` — subroutine

---
## Module: simple_commanders_api

Files:
- `main/apis/simple_commanders_api.f90`

Uses:
- `simple_binoris_io`
- `simple_builder`
- `simple_cmdline`
- `simple_commander_base`
- `simple_core_module_api`
- `simple_default_clines`
- `simple_exec_helpers`
- `simple_image`
- `simple_image_bin`
- `simple_image_msk`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_sp_project`
- `simple_stack_io`
- `simple_stream_utils`

---
## Module: simple_commanders_atoms

Files:
- `main/commanders/simple/simple_commanders_atoms.f90`

Uses:
- `simple_atomic_species`
- `simple_atoms`
- `simple_calpha_finder`
- `simple_commanders_api`
- `simple_commanders_resolest`
- `simple_commanders_sim`
- `simple_commanders_volops`
- `simple_nanoparticle`
- `simple_nanoparticle_utils`

Public symbols:
- `commander_atoms_register` — type
- `commander_atoms_rmsd` — type
- `commander_atoms_stats` — type
- `commander_cif2mrc` — type
- `commander_cif2pdb` — type
- `commander_conv_atom_denoise` — type
- `commander_core_atoms_analysis` — type
- `commander_crys_score` — type
- `commander_detect_atomic_species` — type
- `commander_detect_atoms` — type
- `commander_detect_calpha` — type
- `commander_map2model_fsc` — type
- `commander_map_validate` — type
- `commander_model_validate` — type
- `commander_model_validate_eo` — type
- `commander_pdb2mrc` — type
- `common_atoms` — type
- `exec_atoms_register` — subroutine
- `exec_atoms_rmsd` — subroutine
- `exec_atoms_stats` — subroutine
- `exec_cif2mrc` — subroutine
- `exec_cif2pdb` — subroutine
- `exec_conv_atom_denoise` — subroutine
- `exec_core_atoms_analysis` — subroutine
- `exec_crys_score` — subroutine
- `exec_detect_atomic_species` — subroutine
- `exec_detect_atoms` — subroutine
- `exec_detect_calpha` — subroutine
- `exec_map2model_fsc` — subroutine
- `exec_map_validate` — subroutine
- `exec_model_validate` — subroutine
- `exec_model_validate_eo` — subroutine
- `exec_pdb2mrc` — subroutine

Local symbols (enclosing procedure scope):
- `exec_crys_score::compute_stats` — subroutine
- `exec_crys_score::sort_cluster` — subroutine

---
## Module: simple_commanders_cavgs

Files:
- `main/commanders/simple/simple_commanders_cavgs.f90`

Uses:
- `simple_cavg_quality_analysis`
- `simple_cavg_quality_lattice`
- `simple_cavg_quality_learn`
- `simple_cavg_quality_model`
- `simple_cavg_quality_relations`
- `simple_cavg_quality_selection`
- `simple_cavg_quality_types`
- `simple_cavgs_state_match`
- `simple_clustering_utils`
- `simple_commanders_api`
- `simple_corrmat`
- `simple_imgarr_utils`
- `simple_projfile_utils`
- `simple_srch_sort_loc`
- `simple_strategy2d_utils`

Public symbols:
- `commander_assign_cavgs_states` — type
- `commander_cluster_cavgs` — type
- `commander_cluster_cavgs_selection` — type
- `commander_map_cavgs_selection` — type
- `commander_map_cavgs_states` — type
- `commander_match_cavgs` — type
- `commander_model_cavgs_rejection` — type
- `commander_rank_cavgs` — type
- `commander_select_clusters` — type
- `commander_shape_rank_cavgs` — type
- `exec_assign_cavgs_states` — subroutine
- `exec_cluster_cavgs` — subroutine
- `exec_cluster_cavgs_selection` — subroutine
- `exec_map_cavgs_selection` — subroutine
- `exec_map_cavgs_states` — subroutine
- `exec_match_cavgs` — subroutine
- `exec_model_cavgs_rejection` — subroutine
- `exec_rank_cavgs` — subroutine
- `exec_select_clusters` — subroutine
- `exec_shape_rank_cavgs` — subroutine

Local symbols (enclosing procedure scope):
- `exec_cluster_cavgs::build_temp_project_from_stack` — subroutine
- `exec_cluster_cavgs::prep_stack4clust` — subroutine
- `exec_cluster_cavgs::update_project` — subroutine
- `exec_match_cavgs::find_label_state` — function
- `exec_model_cavgs_rejection::annotate_nano_diameters` — subroutine
- `exec_model_cavgs_rejection::annotate_project` — subroutine
- `exec_model_cavgs_rejection::log_quality_model_selection` — subroutine
- `exec_model_cavgs_rejection::score_cls3d_state_groups` — subroutine
- `exec_model_cavgs_rejection::set_nano_hard_gate_model` — subroutine
- `exec_model_cavgs_rejection::write_hard_gate_stack` — subroutine
- `exec_model_cavgs_rejection::write_ranked_quality_stack` — subroutine
- `exec_shape_rank_cavgs::p1_lt_p2` — function

---
## Module: simple_commanders_checks

Files:
- `main/commanders/simple/simple_commanders_checks.f90`

Uses:
- `simple_commanders_api`
- `simple_decay_funs`
- `simple_ori`
- `simple_rnd`

Public symbols:
- `check_euler_shift` — subroutine
- `commander_check_box` — type
- `commander_check_nptcls` — type
- `commander_check_stoch_update` — type
- `commander_check_update_frac` — type
- `commander_info_image` — type
- `commander_info_stktab` — type
- `exec_check_box` — subroutine
- `exec_check_nptcls` — subroutine
- `exec_check_stoch_update` — subroutine
- `exec_check_update_frac` — subroutine
- `exec_info_image` — subroutine
- `exec_info_stktab` — subroutine

---
## Module: simple_commanders_denoise

Files:
- `main/commanders/simple/simple_commanders_denoise.f90`

Uses:
- `simple_cls_expansion_strategy`
- `simple_commanders_api`

Public symbols:
- `commander_cls_expansion` — type
- `exec_cls_expansion` — subroutine

---
## Module: simple_commanders_distr

Files:
- `main/commanders/simple/simple_commanders_distr.f90`

Uses:
- `simple_commanders_api`
- `simple_map_reduce`

Public symbols:
- `commander_split` — type
- `exec_split` — subroutine

---
## Module: simple_commanders_euclid

Files:
- `main/commanders/simple/simple_commanders_euclid.f90`

Uses:
- `simple_calc_pspec_strategy`
- `simple_commanders_api`
- `simple_strategy_sigma2`

Public symbols:
- `commander_calc_pspec` — type
- `commander_commit_sigma2` — type
- `exec_calc_pspec` — subroutine
- `exec_commit_sigma2` — subroutine

---
## Module: simple_commanders_euclid_distr

Files:
- `main/commanders/simple/simple_commanders_euclid_distr.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_commanders_api`
- `simple_ptcl_layout`
- `simple_sigma2_pspec_sums`
- `simple_sigma2_state`
- `simple_sigma2_state_file`

Public symbols:
- `commander_calc_pspec_assemble` — type
- `exec_calc_pspec_assemble` — subroutine

---
## Module: simple_commanders_flex_pca

Files:
- `main/commanders/simple/simple_commanders_flex_pca.f90`

Uses:
- `simple_commanders_api`
- `simple_flex_pca_strategy`
- `simple_magic_boxes`

Public symbols:
- `apply_flex_pca_defaults` — subroutine
- `apply_flex_pca_pcg_defaults` — subroutine
- `commander_flex_pca` — type
- `commander_flex_pca_worker` — type
- `derive_flex_pca_band` — subroutine
- `derive_flex_pca_sampling` — subroutine
- `exec_flex_pca` — subroutine
- `exec_flex_pca_worker` — subroutine
- `flex_pca_auto_states` — function
- `pickup_project_consensus_volume` — subroutine

---
## Module: simple_commanders_halfmap_enhance

Files:
- `main/commanders/simple/simple_commanders_halfmap_enhance.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_halfmap_enhancer`
- `simple_image`
- `simple_imghead`
- `simple_jiffys`
- `simple_parameters`
- `simple_string`

Public symbols:
- `commander_enhance` — type

Private symbols:
- `exec_enhance` — subroutine

---
## Module: simple_commanders_imgops

Files:
- `main/commanders/simple/simple_commanders_imgops.f90`

Uses:
- `simple_bspline_smoother`
- `simple_commanders_api`
- `simple_gui_communicator`
- `simple_procimgstk`
- `simple_segmentation`

Public symbols:
- `commander_binarize` — type
- `commander_filter` — type
- `commander_normalize` — type
- `commander_scale` — type
- `exec_binarize` — subroutine
- `exec_filter` — subroutine
- `exec_normalize` — subroutine
- `exec_scale` — subroutine

Local symbols (enclosing procedure scope):
- `exec_binarize::doit` — subroutine
- `exec_scale::parse_smpd_target` — subroutine

---
## Module: simple_commanders_imgproc

Files:
- `main/commanders/simple/simple_commanders_imgproc.f90`

Uses:
- `simple_commanders_api`
- `simple_ctf`
- `simple_estimate_ssnr`
- `simple_fileio`
- `simple_matcher_ptcl_io`
- `simple_memoize_ft_maps`
- `simple_ori`
- `simple_procimgstk`
- `simple_segmentation`

Public symbols:
- `commander_ctf_correct` — type
- `commander_ctfops` — type
- `commander_estimate_diam` — type
- `exec_ctf_correct` — subroutine
- `exec_ctfops` — subroutine
- `exec_estimate_diam` — subroutine

Local symbols (enclosing procedure scope):
- `exec_ctf_correct::correct_particles_batch` — subroutine
- `exec_ctf_correct::update_local_project_copy` — subroutine

---
## Module: simple_commanders_mask

Files:
- `main/commanders/simple/simple_commanders_mask.f90`

Uses:
- `simple_atoms`
- `simple_commanders_api`
- `simple_procimgstk`

Public symbols:
- `commander_auto_spher_mask` — type
- `commander_automask` — type
- `commander_automask2D` — type
- `commander_mask` — type
- `exec_auto_spher_mask` — subroutine
- `exec_automask` — subroutine
- `exec_automask2D` — subroutine
- `exec_mask` — subroutine

---
## Module: simple_commanders_mic_subset

Files:
- `main/commanders/simple/simple_commanders_mic_subset.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_jiffys`
- `simple_mic_subset`
- `simple_oris`
- `simple_parameters`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`

Public symbols:
- `commander_select_micrographs` — type

Private symbols:
- `absolute_data_paths` — subroutine
- `csv_field` — function
- `exec_select_micrographs` — subroutine

---
## Module: simple_commanders_misc

Files:
- `main/commanders/simple/simple_commanders_misc.f90`

Uses:
- `simple_class_frcs`
- `simple_commanders_api`
- `simple_fsc`
- `simple_simple_volinterp`

Public symbols:
- `commander_kstest` — type
- `commander_mkdir` — type
- `commander_pearsn` — type
- `commander_print_dose_weights` — type
- `commander_print_fsc` — type
- `commander_print_magic_boxes` — type
- `exec_kstest` — subroutine
- `exec_mkdir` — subroutine
- `exec_pearsn` — subroutine
- `exec_print_dose_weights` — subroutine
- `exec_print_fsc` — subroutine
- `exec_print_magic_boxes` — subroutine
- `read_nrs_dat` — subroutine

---
## Module: simple_commanders_mkcavgs

Files:
- `main/commanders/simple/simple_commanders_mkcavgs.f90`

Uses:
- `simple_cavgassemble_strategy`
- `simple_classaverager`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_make_cavgs_strategy`
- `simple_parameters`
- `simple_pftc_srch_api`
- `simple_refine2d_fnames`

Public symbols:
- `commander_bootstrap_cavgs` — type
- `commander_cavgassemble` — type
- `commander_make_cavgs` — type
- `commander_make_cavgs_distr` — type
- `commander_unbootstrap_cavgs` — type
- `commander_write_classes` — type
- `exec_bootstrap_cavgs` — subroutine
- `exec_cavgassemble` — subroutine
- `exec_make_cavgs` — subroutine
- `exec_make_cavgs_distr` — subroutine
- `exec_unbootstrap_cavgs` — subroutine
- `exec_write_classes` — subroutine
- `make_cavgs_exec_cavgassemble` — subroutine
- `run_make_cavgs_workflow` — subroutine

Local symbols (enclosing procedure scope):
- `exec_bootstrap_cavgs::append_bootstrap_row` — subroutine
- `exec_bootstrap_cavgs::append_skipped_class` — subroutine
- `exec_bootstrap_cavgs::bootstrap_parent_sample` — type
- `exec_bootstrap_cavgs::bootstrap_row` — type
- `exec_bootstrap_cavgs::bootstrap_skipped_class` — type
- `exec_bootstrap_cavgs::build_bootstrap_project` — subroutine
- `exec_bootstrap_cavgs::class_state_for_report` — function
- `exec_bootstrap_cavgs::clear_parent_sample` — subroutine
- `exec_bootstrap_cavgs::collect_bootstrap_classes` — subroutine
- `exec_bootstrap_cavgs::count_eo` — function
- `exec_bootstrap_cavgs::derive_bootstrap_output_names` — subroutine
- `exec_bootstrap_cavgs::rank_of` — function
- `exec_bootstrap_cavgs::report_bootstrap_class_accounting` — subroutine
- `exec_bootstrap_cavgs::select_bootstrap_anchor` — subroutine
- `exec_bootstrap_cavgs::set_bootstrap_class_row` — subroutine
- `exec_bootstrap_cavgs::shuffle_ints` — subroutine
- `exec_bootstrap_cavgs::sum_parent_pops` — function
- `exec_bootstrap_cavgs::validate_evenodd_present` — subroutine
- `exec_bootstrap_cavgs::write_bootstrap_outputs` — subroutine

---
## Module: simple_commanders_motion

Files:
- `main/commanders/simple/simple_commanders_motion.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_core_module_api`
- `simple_parameters`
- `simple_refine_motion_model_strategy`

Public symbols:
- `commander_refine_motion_model` — type

Private symbols:
- `exec_refine_motion_model` — subroutine

---
## Module: simple_commanders_ori

Files:
- `main/commanders/simple/simple_commanders_ori.f90`

Uses:
- `simple_commanders_api`

Public symbols:
- `check_states` — subroutine
- `commander_check_states` — type
- `commander_make_oris` — type
- `commander_measure_projspace_angres` — type
- `commander_oriconsensus` — type
- `commander_oriops` — type
- `commander_oristats` — type
- `commander_rotmats2oris` — type
- `commander_vizoris` — type
- `exec_make_oris` — subroutine
- `exec_measure_projspace_angres` — subroutine
- `exec_oriconsensus` — subroutine
- `exec_oriops` — subroutine
- `exec_oristats` — subroutine
- `exec_rotmats2oris` — subroutine
- `exec_vizoris` — subroutine

Local symbols (enclosing procedure scope):
- `check_states::generate_perm` — subroutine
- `check_states::read_state_order` — subroutine

---
## Module: simple_commanders_pick

Files:
- `main/commanders/simple/simple_commanders_pick.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_api`
- `simple_core_module_api`
- `simple_extract_strategy`
- `simple_gui_communicator`
- `simple_parameters`
- `simple_pick_strategy`
- `simple_picker_iter`
- `simple_reextract_strategy`

Public symbols:
- `commander_extract` — type
- `commander_make_pickrefs` — type
- `commander_pick` — type
- `commander_pick_extract` — type
- `commander_reextract` — type

Private symbols:
- `exec_extract` — subroutine
- `exec_make_pickrefs` — subroutine
- `exec_pick` — subroutine
- `exec_pick_extract` — subroutine
- `exec_reextract` — subroutine

---
## Module: simple_commanders_postprocess_nu

Files:
- `main/commanders/simple/simple_commanders_postprocess_nu.f90`

Uses:
- `simple_commanders_api`
- `simple_nu_filter`

Public symbols:
- `commander_postprocess_nu` — type
- `exec_postprocess_nu` — subroutine

---
## Module: simple_commanders_postprocess_pdb_tester

Files:
- `main/commanders/simple/simple_commanders_postprocess_pdb_tester.f90`

Uses:
- `simple_atoms`
- `simple_cmdline`
- `simple_exec_filter`
- `simple_fileio`
- `simple_image`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_postprocess_pdb_tests` — subroutine

---
## Module: simple_commanders_preprocess

Files:
- `main/commanders/simple/simple_commanders_preprocess.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_api`
- `simple_core_module_api`
- `simple_ctf_estimate_strategy`
- `simple_gen_pspecs_and_thumbs_strategy`
- `simple_gui_communicator`
- `simple_mini_stream_utils`
- `simple_motion_correct_strategy`
- `simple_motion_correct_utils`
- `simple_parameters`
- `simple_preprocess_strategy`

Public symbols:
- `commander_ctf_estimate` — type
- `commander_gen_pspecs_and_thumbs` — type
- `commander_motion_correct` — type
- `commander_preprocess` — type
- `exec_ctf_estimate` — subroutine
- `exec_gen_pspecs_and_thumbs` — subroutine
- `exec_motion_correct` — subroutine
- `exec_preprocess` — subroutine

---
## Module: simple_commanders_project_cls

Files:
- `main/commanders/simple/simple_commanders_project_cls.f90`

Uses:
- `simple_commanders_api`

Public symbols:
- `commander_export_cavgs` — type
- `commander_import_cavgs` — type
- `commander_sample_classes` — type
- `exec_export_cavgs` — subroutine
- `exec_import_cavgs` — subroutine
- `exec_sample_classes` — subroutine

---
## Module: simple_commanders_project_core

Files:
- `main/commanders/simple/simple_commanders_project_core.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_project_ptcl`
- `simple_gui_communicator`
- `simple_imgarr_utils`
- `simple_projfile_utils`
- `simple_srch_sort_loc`
- `simple_strategy2d_utils`

Public symbols:
- `commander_aggregate_chunks` — type
- `commander_concatenate_projects` — type
- `commander_extract_subproj` — type
- `commander_fix_projfile` — type
- `commander_merge_projects` — type
- `commander_new_project` — type
- `commander_print_project_field` — type
- `commander_print_project_info` — type
- `commander_print_project_vals` — type
- `commander_ptcl3D_state_consensus` — type
- `commander_replace_project_field` — type
- `commander_selection` — type
- `commander_update_project` — type
- `exec_aggregate_chunks` — subroutine
- `exec_concatenate_projects` — subroutine
- `exec_extract_subproj` — subroutine
- `exec_fix_projfile` — subroutine
- `exec_merge_projects` — subroutine
- `exec_new_project` — subroutine
- `exec_print_project_field` — subroutine
- `exec_print_project_info` — subroutine
- `exec_print_project_vals` — subroutine
- `exec_ptcl3D_state_consensus` — subroutine
- `exec_replace_project_field` — subroutine
- `exec_selection` — subroutine
- `exec_update_project` — subroutine

Local symbols (enclosing procedure scope):
- `exec_ptcl3D_state_consensus::apply_consensus` — subroutine
- `exec_ptcl3D_state_consensus::best_state_mapping` — subroutine
- `exec_ptcl3D_state_consensus::build_agreement` — subroutine
- `exec_ptcl3D_state_consensus::build_consensus` — subroutine
- `exec_ptcl3D_state_consensus::cleanup` — subroutine
- `exec_ptcl3D_state_consensus::generate_mapping` — subroutine
- `exec_ptcl3D_state_consensus::map_state_correspondence` — subroutine
- `exec_ptcl3D_state_consensus::sanitize_source_labels` — subroutine
- `exec_ptcl3D_state_consensus::write_consensus_report` — subroutine
- `exec_selection::perform_states_selection` — subroutine

---
## Module: simple_commanders_project_mov

Files:
- `main/commanders/simple/simple_commanders_project_mov.f90`

Uses:
- `simple_commanders_api`
- `simple_gui_communicator`
- `simple_stream_watcher`

Public symbols:
- `commander_import_movies` — type
- `commander_write_mic_filetab` — type
- `exec_import_movies` — subroutine
- `exec_write_mic_filetab` — subroutine

---
## Module: simple_commanders_project_ptcl

Files:
- `main/commanders/simple/simple_commanders_project_ptcl.f90`

Uses:
- `simple_commanders_api`
- `simple_gui_communicator`
- `simple_starproject`
- `simple_stream_watcher`

Public symbols:
- `commander_consolidate_particles` — type
- `commander_import_boxes` — type
- `commander_import_particles` — type
- `commander_prune_project` — type
- `commander_prune_project_distr` — type
- `commander_reimport_particles` — type
- `commander_scale_project` — type
- `commander_split_stack` — type
- `commander_zero_project_shifts` — type
- `exec_consolidate_particles` — subroutine
- `exec_import_boxes` — subroutine
- `exec_import_particles` — subroutine
- `exec_prune_project` — subroutine
- `exec_prune_project_distr` — subroutine
- `exec_reimport_particles` — subroutine
- `exec_scale_project` — subroutine
- `exec_split_stack` — subroutine
- `exec_zero_project_shifts` — subroutine

Local symbols (enclosing procedure scope):
- `exec_consolidate_particles::copy_run` — subroutine
- `exec_consolidate_particles::open_src` — subroutine
- `exec_consolidate_particles::report_progress` — subroutine

---
## Module: simple_commanders_ptcl_subset

Files:
- `main/commanders/simple/simple_commanders_ptcl_subset.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_error`
- `simple_jiffys`
- `simple_parameters`
- `simple_ptcl_subset`
- `simple_sp_project`
- `simple_string`
- `simple_sym`
- `simple_syslib`

Public symbols:
- `commander_select_particles` — type

Private symbols:
- `exec_select_particles` — subroutine

---
## Module: simple_commanders_rec

Files:
- `main/commanders/simple/simple_commanders_rec.f90`

Uses:
- `simple_builder`
- `simple_commanders_api`
- `simple_matcher_2dprep`
- `simple_matcher_3drec`
- `simple_parameters`
- `simple_rec3d_pcg_strategy`
- `simple_rec3d_strategy`
- `simple_rec3d_wave_strategy`
- `simple_sigma2_files`
- `simple_state_weight_set`

Public symbols:
- `commander_rec3D` — type
- `commander_rec3D_wave` — type
- `commander_rec3D_worker` — type
- `exec_random_rec` — subroutine
- `exec_rec3D` — subroutine
- `exec_rec3D_distr_worker` — subroutine
- `exec_rec3D_wave` — subroutine
- `random_rec_commander` — type
- `read_part_rows` — subroutine

---
## Module: simple_commanders_rec_distr

Files:
- `main/commanders/simple/simple_commanders_rec_distr.f90`

Uses:
- `simple_commanders_api`
- `simple_volassemble_strategy`

Public symbols:
- `commander_volassemble` — type

Private symbols:
- `exec_volassemble` — subroutine

---
## Module: simple_commanders_refine2D

Files:
- `main/commanders/simple/simple_commanders_refine2D.f90`

Uses:
- `simple_classaverager`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_imgops`
- `simple_commanders_mkcavgs`
- `simple_defs`
- `simple_eul_prob_tab_utils`
- `simple_gui_utils`
- `simple_matcher_smpl_and_lplims`
- `simple_pftc_srch_api`
- `simple_procimgstk`
- `simple_progress`
- `simple_refine2d_strategy`
- `simple_strategy2d_matcher`

Public symbols:
- `commander_refine2D` — type
- `commander_refine2D_distr_worker` — type
- `exec_refine2D` — subroutine
- `exec_refine2D_distr_worker` — subroutine

---
## Module: simple_commanders_refine3D

Files:
- `main/commanders/simple/simple_commanders_refine3D.f90`

Uses:
- `simple_cavgs_work_project`
- `simple_cmdline`
- `simple_commander_base`
- `simple_commanders_api`
- `simple_commanders_euclid`
- `simple_commanders_project_core`
- `simple_commanders_rec`
- `simple_core_module_api`
- `simple_defs`
- `simple_defs_conv`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_eul_prob_tab_utils`
- `simple_external_reference_pose_initialization`
- `simple_fileio`
- `simple_final_rec`
- `simple_flex_pca_stages`
- `simple_gui_communicator`
- `simple_halfmap_diagnostics`
- `simple_imghead`
- `simple_linalg`
- `simple_magic_boxes`
- `simple_math_ft`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_pftc_srch_api`
- `simple_refine3d_fnames`
- `simple_refine3d_stage_plan`
- `simple_refine3d_states_controller`
- `simple_refine3d_strategy`
- `simple_refine3d_workflow_controller`
- `simple_sigma2_bootstrap`
- `simple_solve3d_utils`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_strategy3d_cont`
- `simple_strategy3d_matcher`
- `simple_string_utils`
- `simple_syslib`
- `simple_view_partition_sampling`

Public symbols:
- `commander_bootstrap_rec3D` — type
- `commander_classify3D_refs` — type
- `commander_refine3D` — type
- `commander_refine3D_auto` — type
- `commander_refine3D_distr_worker` — type
- `commander_refine3D_local` — type
- `commander_refine3D_states` — type
- `commander_refine3D_states_cavgs` — type
- `exec_bootstrap_rec3D` — subroutine
- `exec_classify3D_refs` — subroutine
- `exec_refine3D` — subroutine
- `exec_refine3D_auto` — subroutine
- `exec_refine3D_distr_worker` — subroutine
- `exec_refine3D_local` — subroutine
- `exec_refine3D_states` — subroutine
- `exec_refine3D_states_cavgs` — subroutine
- `require_refine3D_local_value` — subroutine
- `run_refine3D_registration_pass` — subroutine

Private symbols:
- `configure_refine3D_grid` — subroutine
- `configure_refine3D_sampling` — subroutine
- `exec_refine3D_workflow` — subroutine
- `prepare_refine3D_state_work_project` — subroutine
- `seed_refine3D_auto_nonuniform_lpset` — subroutine

Local symbols (enclosing procedure scope):
- `exec_bootstrap_rec3D::prepare_bootstrap_rec_cline` — subroutine
- `exec_bootstrap_rec3D::register_bootstrap_rec_outputs` — subroutine
- `exec_bootstrap_rec3D::warn_for_forced_bootstrap_overrides` — subroutine
- `exec_classify3D_refs::ensure_all_active_particles_updated` — subroutine
- `exec_classify3D_refs::initialize_state_volumes` — subroutine
- `exec_classify3D_refs::prepare_classify3D_refs_class_sampling` — subroutine
- `exec_classify3D_refs::read_update_coverage` — subroutine
- `exec_classify3D_refs::reconstruct_all_particles_volumes` — subroutine
- `exec_classify3D_refs::run_classify3D_refs` — subroutine
- `exec_classify3D_refs::run_classify3D_refs_missing_update` — subroutine
- `exec_classify3D_refs::set_classify3D_refs_downscaling` — subroutine
- `exec_classify3D_refs::set_classify3D_refs_nsample` — subroutine
- `exec_classify3D_refs::set_classify3D_refs_nstates` — subroutine
- `exec_classify3D_refs::set_classify3D_refs_sampling` — subroutine
- `exec_classify3D_refs::validate_classify3D_refs_filtering` — subroutine
- `exec_classify3D_refs::validate_classify3D_refs_search_mode` — subroutine
- `exec_classify3D_refs::validate_input_volumes` — subroutine
- `exec_refine3D_states_cavgs::check_flex_fits` — subroutine
- `exec_refine3D_states_cavgs::commit` — subroutine
- `exec_refine3D_states_cavgs::copy_cls3D_labels` — subroutine
- `exec_refine3D_states_cavgs::half_active` — function
- `exec_refine3D_states_cavgs::pick_better_rows` — subroutine
- `exec_refine3D_states_cavgs::preflight` — subroutine
- `exec_refine3D_states_cavgs::refuse_inputs` — subroutine
- `exec_refine3D_states_cavgs::refuse_key` — subroutine
- `exec_refine3D_states_cavgs::renumber_states` — subroutine
- `exec_refine3D_states_cavgs::require_value` — subroutine
- `exec_refine3D_states_cavgs::seed_poses` — subroutine
- `exec_refine3D_states_cavgs::write_final_maps` — subroutine

---
## Module: simple_commanders_relion

Files:
- `main/commanders/simple/simple_commanders_relion.f90`

Uses:
- `simple_commanders_api`
- `simple_relion`

Public symbols:
- `commander_export_relion` — type
- `exec_export_relion` — subroutine

---
## Module: simple_commanders_reproject

Files:
- `main/commanders/simple/simple_commanders_reproject.f90`

Uses:
- `simple_commanders_api`
- `simple_gui_communicator`
- `simple_imgarr_utils`
- `simple_projector`
- `simple_simple_volinterp`

Public symbols:
- `commander_reproject` — type
- `exec_reproject` — subroutine

---
## Module: simple_commanders_resolest

Files:
- `main/commanders/simple/simple_commanders_resolest.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_euclid`
- `simple_fsc`
- `simple_nu_filter`
- `simple_opt_filter`
- `simple_pftc_srch_api`
- `simple_refine3d_fnames`

Public symbols:
- `commander_estimate_lpstages` — type
- `commander_fsc` — type
- `commander_fsc_area_score` — type
- `commander_nu_filt3D` — type
- `commander_uniform_filter2D` — type
- `commander_uniform_filter3D` — type
- `exec_estimate_lpstages` — subroutine
- `exec_fsc` — subroutine
- `exec_fsc_area_score` — subroutine
- `exec_nu_filt3D` — subroutine
- `exec_uniform_filter2D` — subroutine
- `exec_uniform_filter3D` — subroutine

Local symbols (enclosing procedure scope):
- `exec_estimate_lpstages::calc_lplim_final_stage` — function

---
## Module: simple_commanders_sieve

Files:
- `main/commanders/simple/simple_commanders_sieve.f90`

Uses:
- `simple_commanders_api`
- `simple_ptcl_sieve`
- `simple_ptcl_sieve_utils`
- `simple_rec_list`

Public symbols:
- `commander_sieve_ptcls` — type
- `exec_sieve_ptcls` — subroutine

---
## Module: simple_commanders_sim

Files:
- `main/commanders/simple/simple_commanders_sim.f90`

Uses:
- `simple_atoms`
- `simple_commanders_api`
- `simple_ctf`
- `simple_defs_atoms`
- `simple_matcher_ptcl_io`
- `simple_math_ctf`
- `simple_memoize_ft_maps`
- `simple_projector`
- `simple_rnd`
- `simple_simple_volinterp`
- `simple_simulator`
- `simple_string_utils`
- `simple_wave_imaging`
- `simple_wave_projector`

Public symbols:
- `commander_simulate_movie` — type
- `commander_simulate_nanoparticle` — type
- `commander_simulate_noise` — type
- `commander_simulate_particles` — type
- `commander_simulate_scattering` — type
- `commander_simulate_subtomogram` — type
- `exec_simulate_movie` — subroutine
- `exec_simulate_nanoparticle` — subroutine
- `exec_simulate_noise` — subroutine
- `exec_simulate_particles` — subroutine
- `exec_simulate_scattering` — subroutine
- `exec_simulate_subtomogram` — subroutine

Private symbols:
- `draw_states` — function
- `oritab_holds_states` — function
- `simulate_particles_wave` — subroutine
- `state_fractions` — function

Local symbols (enclosing procedure scope):
- `exec_simulate_movie::gen_ptcl_pos` — function

---
## Module: simple_commanders_solve2D

Files:
- `main/commanders/simple/simple_commanders_solve2D.f90`

Uses:
- `simple_bench_io`
- `simple_class_frcs`
- `simple_classaverager`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_euclid`
- `simple_commanders_imgops`
- `simple_commanders_refine2d`
- `simple_commanders_volops`
- `simple_defs_environment`
- `simple_eo_partition`
- `simple_gui_communicator`
- `simple_procimgstk`
- `simple_rec_list`
- `simple_refine2d_fnames`
- `simple_sigma2_bootstrap`
- `simple_solve2d_chunk`
- `simple_solve2d_controller`
- `simple_timer`

Public symbols:
- `commander_solve2D` — type
- `commander_solve2D_chunks` — type
- `execute_solve2D_staged` — subroutine

Private symbols:
- `exec_solve2D` — subroutine
- `exec_solve2D_chunks` — subroutine
- `exec_solve2D_workflow` — subroutine
- `init_chunk_clustering` — subroutine
- `setup_chunk_downscaling` — subroutine
- `terminal_cline_policy` — subroutine

Local symbols (enclosing procedure scope):
- `exec_solve2D_chunks::analyze_first_available_chunk` — subroutine
- `exec_solve2D_chunks::automatic_nchunks` — function
- `exec_solve2D_chunks::flag_converged_chunks` — subroutine
- `exec_solve2D_chunks::generate_chunk_projects` — subroutine
- `exec_solve2D_chunks::init_one_chunk` — subroutine
- `exec_solve2D_workflow::adopt_estimated_mskdiam` — subroutine
- `exec_solve2D_workflow::ensure_resume_sigma_state` — subroutine
- `exec_solve2D_workflow::ensure_seed_eo` — subroutine
- `exec_solve2D_workflow::execute_refine2D` — subroutine
- `exec_solve2D_workflow::execute_seed_pass` — subroutine
- `exec_solve2D_workflow::execute_terminal_pass` — subroutine
- `exec_solve2D_workflow::final_mskdiam` — function
- `exec_solve2D_workflow::gen_final_cavgs` — subroutine
- `exec_solve2D_workflow::inirefs` — subroutine
- `exec_solve2D_workflow::output_stats` — subroutine
- `exec_solve2D_workflow::prep_command_lines` — subroutine
- `exec_solve2D_workflow::record_stage_frc05` — subroutine
- `exec_solve2D_workflow::rederive_lplims` — subroutine
- `exec_solve2D_workflow::seed_from_previous_clustering` — subroutine
- `exec_solve2D_workflow::set_dims` — subroutine
- `exec_solve2D_workflow::set_lplims` — subroutine
- `exec_solve2D_workflow::set_sampling` — subroutine
- `exec_solve2D_workflow::write_seed_lineage` — subroutine
- `exec_solve2D_workflow::write_solve2D_benchmark` — subroutine
- `init_chunk_clustering::set_chunk_dimensions` — subroutine

---
## Module: simple_commanders_solve3D

Files:
- `main/commanders/simple/simple_commanders_solve3D.f90`

Uses:
- `simple_bench_io`
- `simple_cluster_seed`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_project_core`
- `simple_commanders_rec`
- `simple_commanders_refine3d`
- `simple_commanders_reproject`
- `simple_eo_partition`
- `simple_estimate_ssnr`
- `simple_external_reference_pose_initialization`
- `simple_final_rec`
- `simple_frozen_accum`
- `simple_gui_communicator`
- `simple_halfmap_diagnostics`
- `simple_procimgstk`
- `simple_project_superset`
- `simple_ptcl_layout`
- `simple_refine3d_fnames`
- `simple_sigma2_state`
- `simple_sigma2_state_file`
- `simple_solve3d_addon_report`
- `simple_solve3d_manifest`
- `simple_solve3d_utils`
- `simple_ui`
- `simple_ui_program`
- `simple_view_partition_sampling`

Public symbols:
- `commander_solve3D` — type
- `commander_solve3D_addon` — type
- `commander_solve3D_cavgs` — type
- `commander_solve3D_cavgs_conditional_restarts` — type

Private symbols:
- `exec_solve3D` — subroutine
- `exec_solve3D_addon` — subroutine
- `exec_solve3D_cavgs` — subroutine
- `exec_solve3D_cavgs_conditional_restarts` — subroutine
- `new_run_id` — function

Local symbols (enclosing procedure scope):
- `exec_solve3D::addon_boxes` — subroutine
- `exec_solve3D::addon_cohort_diagnostic` — subroutine
- `exec_solve3D::addon_epilogue` — subroutine
- `exec_solve3D::addon_final` — subroutine
- `exec_solve3D::addon_parse` — subroutine
- `exec_solve3D::addon_prologue` — subroutine
- `exec_solve3D::addon_register_pass` — subroutine
- `exec_solve3D::addon_stage_boundary` — subroutine
- `exec_solve3D::addon_starting_state` — subroutine
- `exec_solve3D::adopt_frozen_bundle` — subroutine
- `exec_solve3D::clean_ptcl3D_sampling` — subroutine
- `exec_solve3D::ensure_multistate_particle_assignments` — subroutine
- `exec_solve3D::ini3D_from_cavgs` — subroutine
- `exec_solve3D::next_refine3D_iteration` — function
- `exec_solve3D::open_addon_sigma2_state` — subroutine
- `exec_solve3D::prepare_state_continue_project` — subroutine
- `exec_solve3D::read_multistate_update_coverage` — subroutine
- `exec_solve3D::record_emitted_limits` — subroutine
- `exec_solve3D::report_frozen_provenance` — subroutine
- `exec_solve3D::reset_ptcl3D_from_ptcl2D_selection` — subroutine
- `exec_solve3D::run_multistate_missing_update` — subroutine
- `exec_solve3D::set_balance_default` — subroutine
- `exec_solve3D::start_from_frozen_artifact` — subroutine
- `exec_solve3D::switch_addon_membership` — subroutine
- `exec_solve3D::validate_cavg_ini_ext_states` — subroutine
- `exec_solve3D::verify_frozen_sigma` — subroutine
- `exec_solve3D::write_run_manifest` — subroutine
- `exec_solve3D::write_solve3D_benchmark` — subroutine
- `exec_solve3D_addon::log_override` — subroutine
- `exec_solve3D_addon::publish_result` — subroutine
- `exec_solve3D_cavgs::configure_cavgs_distributed_clines` — subroutine
- `exec_solve3D_cavgs::conv_eo` — subroutine
- `exec_solve3D_cavgs::conv_eo_states` — subroutine
- `exec_solve3D_cavgs::del_pproc_vols` — subroutine
- `exec_solve3D_cavgs::rank_cavgs` — subroutine
- `exec_solve3D_cavgs::reseed_weak_states` — subroutine
- `exec_solve3D_cavgs::rndstart` — subroutine
- `exec_solve3D_cavgs::strip_distributed_child` — subroutine
- `exec_solve3D_cavgs::sync_distributed_child` — subroutine

---
## Module: simple_commanders_solve3D_chunks

Files:
- `main/commanders/simple/simple_commanders_solve3D_chunks.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_solve3d`
- `simple_refine3d_fnames`
- `simple_solve3d_addon_report`
- `simple_solve3d_chunks_plan`
- `simple_solve3d_manifest`
- `simple_solve3d_utils`
- `simple_ui`
- `simple_ui_program`

Public symbols:
- `commander_solve3D_chunks` — type

Private symbols:
- `collect_solve3D_step_stats` — subroutine
- `copy_forwardable_addon_arg` — subroutine
- `exec_solve3D_chunks` — subroutine
- `publish_solve3D_result` — subroutine
- `write_solve3D_step_project` — subroutine

Local symbols (enclosing procedure scope):
- `exec_solve3D_chunks::forward_keys` — subroutine

---
## Module: simple_commanders_solve3D_pyramid

Files:
- `main/commanders/simple/simple_commanders_solve3D_pyramid.f90`

Uses:
- `simple_cavg_quality_model`
- `simple_chash`
- `simple_cmdline`
- `simple_commander_base`
- `simple_commanders_cavgs`
- `simple_commanders_denoise`
- `simple_commanders_refine3d`
- `simple_commanders_solve2d`
- `simple_commanders_solve3d`
- `simple_defs`
- `simple_defs_fname`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_fileio`
- `simple_jiffys`
- `simple_math_ft`
- `simple_parameters`
- `simple_sp_project`
- `simple_state_assignment`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `commander_solve3D_pyramid` — type

Private symbols:
- `command_text` — function
- `exec_solve3D_pyramid` — subroutine
- `pyramid_report` — type
- `write_report` — subroutine

Local symbols (enclosing procedure scope):
- `exec_solve3D_pyramid::fail` — subroutine
- `exec_solve3D_pyramid::run_job` — subroutine
- `exec_solve3D_pyramid::state_resolutions` — subroutine
- `exec_solve3D_pyramid::take_2D` — subroutine

---
## Module: simple_commanders_starproject

Files:
- `main/commanders/simple/simple_commanders_starproject.f90`

Uses:
- `simple_commanders_api`
- `simple_jiffys`
- `simple_starproject`

Public symbols:
- `commander_assign_optics_groups` — type
- `commander_export_starproject` — type
- `commander_import_starproject` — type
- `exec_assign_optics_groups` — subroutine
- `exec_export_starproject` — subroutine
- `exec_import_starproject` — subroutine

---
## Module: simple_commanders_stkops

Files:
- `main/commanders/simple/simple_commanders_stkops.f90`

Uses:
- `simple_clustering_utils`
- `simple_commanders_api`
- `simple_imgarr_utils`
- `simple_procimgstk`
- `simple_stackops`
- `simple_strategy2d_utils`

Public symbols:
- `commander_cluster_stack` — type
- `commander_convert` — type
- `commander_match_stacks` — type
- `commander_stack` — type
- `commander_stackops` — type
- `exec_cluster_stack` — subroutine
- `exec_convert` — subroutine
- `exec_match_stacks` — subroutine
- `exec_stack` — subroutine
- `exec_stackops` — subroutine

Local symbols (enclosing procedure scope):
- `exec_stackops::finish_cleanup` — subroutine
- `exec_stackops::handle_random_selection` — subroutine
- `exec_stackops::handle_state_class_selection` — subroutine

---
## Module: simple_commanders_stream_p00_master

Files:
- `main/commanders/stream/simple_commanders_stream_p00_master.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_commander_base`
- `simple_commanders_stream_p01_preprocess`
- `simple_commanders_stream_p02_assign_optics`
- `simple_commanders_stream_p03_initial_analysis`
- `simple_commanders_stream_p04_refpick_extract`
- `simple_commanders_stream_p05_sieve_cavgs`
- `simple_commanders_stream_p06_pool2d`
- `simple_commanders_stream_p07_solve3d_multistate`
- `simple_defs`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_error`
- `simple_fileio`
- `simple_gui_assembler`
- `simple_gui_metadata_utils`
- `simple_http_post`
- `simple_memory_monitor`
- `simple_parameters`
- `simple_qsys_env`
- `simple_stream_master_gui_commands`
- `simple_stream_master_meta_store`
- `simple_stream_master_resources`
- `simple_stream_master_stage`
- `simple_stream_master_stage_ids`
- `simple_stream_sigterm`
- `simple_string`
- `simple_syslib`
- `simple_timer`
- `unix`

Public symbols:
- `commander_stream_p00_master` — type

Private symbols:
- `exec_stream_p00_master` — subroutine
- `lock` — subroutine
- `master_shared` — type
- `metadata_listener` — subroutine
- `stage_commander` — function
- `unlock` — subroutine

Local symbols (enclosing procedure scope):
- `exec_stream_p00_master::abort_startup` — subroutine
- `exec_stream_p00_master::apply_commands` — subroutine
- `exec_stream_p00_master::ask_to_stop` — subroutine
- `exec_stream_p00_master::log_memory` — subroutine
- `exec_stream_p00_master::make_stage_clines` — subroutine
- `exec_stream_p00_master::send_heartbeat` — subroutine
- `exec_stream_p00_master::send_last_heartbeat` — subroutine
- `exec_stream_p00_master::stop_stages` — subroutine
- `exec_stream_p00_master::wait_for_stop` — subroutine

---
## Module: simple_commanders_stream_p01_preprocess

Files:
- `main/commanders/stream/simple_commanders_stream_p01_preprocess.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_defs_stream`
- `simple_error`
- `simple_jiffys`
- `simple_stream_sigterm`
- `simple_stream_stage_preprocess`
- `simple_string`

Public symbols:
- `commander_stream_p01_preprocess` — type
- `set_preprocess_stream_cline` — subroutine

Private symbols:
- `exec_stream_p01_preprocess` — subroutine
- `set_flipgain` — subroutine

---
## Module: simple_commanders_stream_p02_assign_optics

Files:
- `main/commanders/stream/simple_commanders_stream_p02_assign_optics.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_error`
- `simple_jiffys`
- `simple_stream_sigterm`
- `simple_stream_stage_optics`

Public symbols:
- `commander_stream_p02_assign_optics` — type

Private symbols:
- `exec_stream_p02_assign_optics` — subroutine

---
## Module: simple_commanders_stream_p03_initial_analysis

Files:
- `main/commanders/stream/simple_commanders_stream_p03_initial_analysis.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_error`
- `simple_jiffys`
- `simple_stream_sigterm`
- `simple_stream_stage_initial_analysis`

Public symbols:
- `commander_stream_p03_initial_analysis` — type
- `set_initial_analysis_cline` — subroutine

Private symbols:
- `check_pickrefs_settings` — subroutine
- `exec_stream_p03_initial_analysis` — subroutine

Local symbols (enclosing procedure scope):
- `set_initial_analysis_cline::default_unless_positive` — subroutine

---
## Module: simple_commanders_stream_p04_refpick_extract

Files:
- `main/commanders/stream/simple_commanders_stream_p04_refpick_extract.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_error`
- `simple_jiffys`
- `simple_stream_sigterm`
- `simple_stream_stage_refpick`

Public symbols:
- `commander_stream_p04_refpick_extract` — type
- `set_refpick_cline` — subroutine

Private symbols:
- `exec_stream_p04_refpick_extract` — subroutine

---
## Module: simple_commanders_stream_p05_sieve_cavgs

Files:
- `main/commanders/stream/simple_commanders_stream_p05_sieve_cavgs.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_jiffys`
- `simple_ptcl_sieve`
- `simple_stream_sigterm`
- `simple_stream_stage_sieve`

Public symbols:
- `commander_stream_p05_sieve_cavgs` — type
- `set_sieve_cline` — subroutine

Private symbols:
- `exec_stream_p05_sieve_cavgs` — subroutine

---
## Module: simple_commanders_stream_p06_pool2D

Files:
- `main/commanders/stream/simple_commanders_stream_p06_pool2D.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_defs_fname`
- `simple_jiffys`
- `simple_stream_sigterm`
- `simple_stream_stage_pool2d`

Public symbols:
- `commander_stream_p06_pool2D` — type
- `set_pool2D_cline` — subroutine

Private symbols:
- `exec_stream_p06_pool2D` — subroutine

---
## Module: simple_commanders_stream_p07_solve3D_multistate

Files:
- `main/commanders/stream/simple_commanders_stream_p07_solve3D_multistate.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_defs_fname`
- `simple_error`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_jiffys`
- `simple_stream_sigterm`
- `simple_stream_stage_solve3d`

Public symbols:
- `commander_stream_p07_solve3D_multistate` — type
- `set_solve3D_cline` — subroutine

Private symbols:
- `exec_stream_p07_solve3D_multistate` — subroutine

---
## Module: simple_commanders_test_class

Files:
- `main/commanders/test/simple_commanders_test_class.f90`

Uses:
- `simple_accum_blend_tester`
- `simple_aff_prop_tester`
- `simple_ansi_ctrls_tester`
- `simple_atoms_tester`
- `simple_binoris_tester`
- `simple_bspline_smoother_tester`
- `simple_calpha_finder_tester`
- `simple_cartesian_fourier_tester`
- `simple_cartft_calc_tester`
- `simple_cartft_pose_opt_tester`
- `simple_cavg_quality_lattice_tester`
- `simple_cavg_quality_relations_tester`
- `simple_cavg_quality_selection_tester`
- `simple_cavg_registration_tester`
- `simple_cavg_sums_tester`
- `simple_cavgs_state_match_tester`
- `simple_cavgs_work_project_tester`
- `simple_chash_tester`
- `simple_cif2mrc_tester`
- `simple_class_compatibility_tester`
- `simple_class_sample_io_tester`
- `simple_cmdline_tester`
- `simple_commanders_api`
- `simple_commanders_postprocess_pdb_tester`
- `simple_commanders_volops_tester`
- `simple_convergence_tester`
- `simple_ctf_tester`
- `simple_decay_funs_tester`
- `simple_dose_series_manifest_tester`
- `simple_eo_partition_tester`
- `simple_epu_xml_tester`
- `simple_eul_prob_tab2d_tester`
- `simple_eul_prob_tab_tester`
- `simple_fileio_tester`
- `simple_flex_cls_expansion_tester`
- `simple_flex_cls_kernel`
- `simple_flex_pca_application_tester`
- `simple_flex_pca_tester`
- `simple_flex_pcg_tester`
- `simple_forked_process_tester`
- `simple_frozen_accum_tester`
- `simple_ftexp_shsrch_tester`
- `simple_ftiter_tester`
- `simple_gauran_tester`
- `simple_gmm_tester`
- `simple_gui_assembler_tester`
- `simple_gui_metadata_tester`
- `simple_hac_tester`
- `simple_halfmap_enhancer_tester`
- `simple_hash_tester`
- `simple_http_post_tester`
- `simple_image_msk_tester`
- `simple_image_serialize_tester`
- `simple_image_tester`
- `simple_imghead_tester`
- `simple_ipc_tcp_socket_tester`
- `simple_jpg_tester`
- `simple_kbinterpol_tester`
- `simple_kcenter_tester`
- `simple_kmeans_tester`
- `simple_linalg_tester`
- `simple_linked_list_tester`
- `simple_lpstages_tester`
- `simple_mic_import_tester`
- `simple_mic_selection_tester`
- `simple_mic_subset_tester`
- `simple_motion_gain_tester`
- `simple_motion_model_tester`
- `simple_mrc_validate_tester`
- `simple_nano_density_tester`
- `simple_nano_species_tester`
- `simple_nlcg_solver_tester`
- `simple_nu_filter_tester`
- `simple_online_var_tester`
- `simple_openmp_offload_tester`
- `simple_opt_tester`
- `simple_optics_groups_tester`
- `simple_optics_maps_tester`
- `simple_ori_tester`
- `simple_oris_tester`
- `simple_parameters_tester`
- `simple_particle_extractor_tester`
- `simple_pcg_atomic_prior_tester`
- `simple_pcg_halfset_tester`
- `simple_pcg_lattice_tester`
- `simple_pcg_sparse_tester`
- `simple_pdb2mrc_tester`
- `simple_persistent_worker_message_tester`
- `simple_persistent_worker_server_tester`
- `simple_pftc_inplane_tester`
- `simple_polar_shiftmap_tester`
- `simple_polarft_corr_tester`
- `simple_project_merge_tester`
- `simple_project_superset_tester`
- `simple_ptcl_sieve_tester`
- `simple_ptcl_subset_tester`
- `simple_qsys_ctrl_tester`
- `simple_qsys_env_tester`
- `simple_quality_logistic_tester`
- `simple_rec3d_service_tester`
- `simple_rec3d_strategy_tester`
- `simple_rec3d_wave_strategy_tester`
- `simple_rec_list_tester`
- `simple_reconstructor_tester`
- `simple_refine3d_local_tester`
- `simple_refine3d_provenance_tester`
- `simple_refine3d_startup_tester`
- `simple_refine3d_strategy_tester`
- `simple_refine_motion_model_tester`
- `simple_rnd_tester`
- `simple_segmentation_tester`
- `simple_sigma2_bootstrap_tester`
- `simple_sigma2_state_tester`
- `simple_simulator_tester`
- `simple_solve2d_controller_tester`
- `simple_solve3d_addon_dock_tester`
- `simple_solve3d_addon_report_tester`
- `simple_solve3d_chunks_plan_tester`
- `simple_solve3d_manifest_tester`
- `simple_solve3d_reseed_tester`
- `simple_sp_project_tester`
- `simple_srch_sort_loc_tester`
- `simple_stack_io_tester`
- `simple_starfile_tester`
- `simple_starproject_stream_tester`
- `simple_starproject_tester`
- `simple_stat_tester`
- `simple_state_assignment_tester`
- `simple_state_weight_set_tester`
- `simple_strategy3d_cont_1jyx_tester`
- `simple_strategy3d_cont_tester`
- `simple_strategy3d_inplane_tester`
- `simple_stream_chain_tester`
- `simple_stream_job_sets_tester`
- `simple_stream_master_tester`
- `simple_stream_meta_plots_tester`
- `simple_stream_pipe_tester`
- `simple_stream_pool2d_tester`
- `simple_stream_stage_initial_analysis_tester`
- `simple_stream_stage_optics_tester`
- `simple_stream_stage_pool2d_tester`
- `simple_stream_stage_preprocess_tester`
- `simple_stream_stage_refpick_tester`
- `simple_stream_stage_sieve_tester`
- `simple_stream_stage_solve3d_tester`
- `simple_stream_tester`
- `simple_stream_watcher_tester`
- `simple_string_tester`
- `simple_sym_tester`
- `simple_syslib_tester`
- `simple_test_dispatch_tester`
- `simple_test_tools_tester`
- `simple_test_truth_metrics_tester`
- `simple_test_utils`
- `simple_trail_chain_manifest_tester`
- `simple_ui`
- `simple_ui_hash_tester`
- `simple_ui_visibility_tester`
- `simple_volpair_metrics_tester`
- `simple_volume_shape_tester`
- `simple_vrefhash_tester`
- `simple_wave_field_tester`
- `simple_wave_imaging_tester`
- `simple_wave_multislice_tester`
- `simple_wave_projector_tester`
- `simple_xd_gmm_tester`

Public symbols:
- `add_suite` — subroutine
- `commander_test_forked_process` — type
- `commander_test_lib_cart_align3D` — type
- `commander_test_lib_heterogeneity` — type
- `commander_test_lib_reconstruction` — type
- `commander_test_lib_single` — type
- `commander_test_lib_stream` — type
- `commander_test_openmp_offload` — type
- `commander_test_stream_chains` — type
- `commander_test_unit_cart_align3D` — type
- `commander_test_unit_core` — type
- `commander_test_unit_heterogeneity` — type
- `commander_test_unit_image` — type
- `commander_test_unit_ipc` — type
- `commander_test_unit_motion` — type
- `commander_test_unit_numerics` — type
- `commander_test_unit_ori` — type
- `commander_test_unit_parallel` — type
- `commander_test_unit_pftc_align2D3D` — type
- `commander_test_unit_project` — type
- `commander_test_unit_reconstruction` — type
- `commander_test_unit_single` — type
- `commander_test_unit_stream` — type
- `commander_test_unit_ui` — type
- `commander_test_units` — type
- `exec_test_forked_process` — subroutine
- `exec_test_lib_cart_align3D` — subroutine
- `exec_test_lib_heterogeneity` — subroutine
- `exec_test_lib_reconstruction` — subroutine
- `exec_test_lib_single` — subroutine
- `exec_test_lib_stream` — subroutine
- `exec_test_openmp_offload` — subroutine
- `exec_test_stream_chains` — subroutine
- `exec_test_unit_cart_align3D` — subroutine
- `exec_test_unit_core` — subroutine
- `exec_test_unit_heterogeneity` — subroutine
- `exec_test_unit_image` — subroutine
- `exec_test_unit_ipc` — subroutine
- `exec_test_unit_motion` — subroutine
- `exec_test_unit_numerics` — subroutine
- `exec_test_unit_ori` — subroutine
- `exec_test_unit_parallel` — subroutine
- `exec_test_unit_pftc_align2D3D` — subroutine
- `exec_test_unit_project` — subroutine
- `exec_test_unit_reconstruction` — subroutine
- `exec_test_unit_single` — subroutine
- `exec_test_unit_stream` — subroutine
- `exec_test_unit_ui` — subroutine
- `exec_test_units` — subroutine
- `no_arg_test` — subroutine
- `run_unit_suites` — subroutine
- `suite_id` — function
- `suite_ui_json` — subroutine
- `suites_cart_align3D` — subroutine
- `suites_core` — subroutine
- `suites_heterogeneity` — subroutine
- `suites_image` — subroutine
- `suites_ipc` — subroutine
- `suites_lib_cart_align3D` — subroutine
- `suites_lib_heterogeneity` — subroutine
- `suites_lib_reconstruction` — subroutine
- `suites_lib_single` — subroutine
- `suites_lib_stream` — subroutine
- `suites_motion` — subroutine
- `suites_numerics` — subroutine
- `suites_ori` — subroutine
- `suites_parallel` — subroutine
- `suites_pftc_align2D3D` — subroutine
- `suites_project` — subroutine
- `suites_reconstruction` — subroutine
- `suites_single` — subroutine
- `suites_stream` — subroutine
- `suites_stream_chains` — subroutine
- `suites_ui` — subroutine
- `test_euler_shift` — subroutine
- `unit_suite` — type

---
## Module: simple_commanders_test_highlevel

Files:
- `main/commanders/test/simple_commanders_test_highlevel.f90`
- `main/commanders/test/simple_commanders_test_highlevel_addon.f90`
- `main/commanders/test/simple_commanders_test_highlevel_backends.f90`
- `main/commanders/test/simple_commanders_test_highlevel_cont.f90`
- `main/commanders/test/simple_commanders_test_highlevel_fixture.f90`
- `main/commanders/test/simple_commanders_test_highlevel_flex.f90`
- `main/commanders/test/simple_commanders_test_highlevel_frac.f90`
- `main/commanders/test/simple_commanders_test_highlevel_matchbench.f90`
- `main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90`
- `main/commanders/test/simple_commanders_test_highlevel_pcg.f90`
- `main/commanders/test/simple_commanders_test_highlevel_stream.f90`
- `main/commanders/test/simple_commanders_test_highlevel_suite.f90`
- `main/commanders/test/simple_commanders_test_highlevel_workflow.f90`

Uses:
- `ieee_arithmetic`
- `iso_c_binding`
- `iso_fortran_env`
- `simple_atoms`
- `simple_bench_io`
- `simple_builder`
- `simple_cartft_pose_opt`
- `simple_class_frcs`
- `simple_commanders_api`
- `simple_commanders_denoise`
- `simple_commanders_flex_pca`
- `simple_commanders_pick`
- `simple_commanders_preprocess`
- `simple_commanders_project_cls`
- `simple_commanders_project_core`
- `simple_commanders_project_mov`
- `simple_commanders_project_ptcl`
- `simple_commanders_rec`
- `simple_commanders_refine3d`
- `simple_commanders_reproject`
- `simple_commanders_sim`
- `simple_commanders_solve2d`
- `simple_commanders_solve3d`
- `simple_commanders_validate`
- `simple_defs_fname`
- `simple_eo_partition`
- `simple_euclid_sigma2`
- `simple_fileio`
- `simple_flex_pca_application_tester`
- `simple_flex_pca_embedding_io`
- `simple_frozen_accum`
- `simple_image`
- `simple_imghead`
- `simple_matcher_ptcl_io`
- `simple_molecule_data`
- `simple_oris`
- `simple_parameters`
- `simple_polar_matching_bench`
- `simple_rec3d_pcg_strategy`
- `simple_reconstructor_pcg`
- `simple_refine3d_fnames`
- `simple_refine3d_provenance`
- `simple_rnd`
- `simple_sigma2_files`
- `simple_sigma2_state_file`
- `simple_simple_volinterp`
- `simple_solve3d_addon_report`
- `simple_solve3d_manifest`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_string_utils`
- `simple_sym`
- `simple_syslib`
- `simple_test_dispatch`
- `simple_test_exec_class`
- `simple_test_exec_highlevel`
- `simple_test_exec_single`
- `simple_test_exec_stream`
- `simple_test_gate`
- `simple_test_identity`
- `simple_test_profile_stats`
- `simple_test_registry`
- `simple_test_state_fixture`
- `simple_test_truth_metrics`
- `simple_test_utils`
- `simple_trail_chain_manifest`
- `simple_ui`

Public symbols:
- `commander_generate_solve3D_addon_snapshots` — type
- `commander_test_cont_refine3D_1jxy` — type
- `commander_test_cont_refine3D_bench` — type
- `commander_test_flex_pca_blobs` — type
- `commander_test_identity` — type
- `commander_test_mini_stream` — type
- `commander_test_multistate_workflows` — type
- `commander_test_nightly` — type
- `commander_test_pcg_frac_update` — type
- `commander_test_pcg_recon` — type
- `commander_test_polar_matching_bench` — type
- `commander_test_profile` — type
- `commander_test_rec3D_backends` — type
- `commander_test_simulate_particles` — type
- `commander_test_simulated_workflow` — type
- `commander_test_solve3D_addon` — type
- `commander_test_state_fixture` — type
- `commander_test_write_state_weights_labels` — type
- `commander_test_write_state_weights_mixed` — type
- `exec_generate_solve3D_addon_snapshots` — subroutine
- `exec_generate_solve3D_addon_snapshots` — subroutine
- `exec_test_cont_refine3D_1jxy` — subroutine
- `exec_test_cont_refine3D_1jxy` — subroutine
- `exec_test_cont_refine3D_bench` — subroutine
- `exec_test_cont_refine3D_bench` — subroutine
- `exec_test_flex_pca_blobs` — subroutine
- `exec_test_flex_pca_blobs` — subroutine
- `exec_test_identity` — subroutine
- `exec_test_identity` — subroutine
- `exec_test_mini_stream_quantitative` — subroutine
- `exec_test_mini_stream_quantitative` — subroutine
- `exec_test_multistate_workflows` — subroutine
- `exec_test_multistate_workflows` — subroutine
- `exec_test_nightly` — subroutine
- `exec_test_nightly` — subroutine
- `exec_test_pcg_frac_update` — subroutine
- `exec_test_pcg_frac_update` — subroutine
- `exec_test_pcg_recon` — subroutine
- `exec_test_pcg_recon` — subroutine
- `exec_test_polar_matching_bench` — subroutine
- `exec_test_polar_matching_bench` — subroutine
- `exec_test_profile` — subroutine
- `exec_test_profile` — subroutine
- `exec_test_rec3D_backends` — subroutine
- `exec_test_rec3D_backends` — subroutine
- `exec_test_simulate_particles` — subroutine
- `exec_test_simulate_particles` — subroutine
- `exec_test_simulated_workflow` — subroutine
- `exec_test_simulated_workflow` — subroutine
- `exec_test_solve3D_addon` — subroutine
- `exec_test_solve3D_addon` — subroutine
- `exec_test_state_fixture` — subroutine
- `exec_test_state_fixture` — subroutine
- `exec_test_write_state_weights_labels` — subroutine
- `exec_test_write_state_weights_labels` — subroutine
- `exec_test_write_state_weights_mixed` — subroutine
- `exec_test_write_state_weights_mixed` — subroutine

Private symbols:
- `backends_run_summary` — type
- `best_latent_correlation` — function
- `check_betagal` — subroutine
- `check_blob_modes` — subroutine
- `check_case` — subroutine
- `check_memberless_halves` — subroutine
- `check_path` — subroutine
- `core_budget` — function
- `evaluate` — subroutine
- `fields_equal` — function
- `fills_case` — subroutine
- `fixed` — function
- `flex_gate_result` — type
- `identity_case` — subroutine
- `inspect_eigenvolume` — subroutine
- `inspect_flex_result` — subroutine
- `inspect_states_cavgs` — subroutine
- `key_value` — function
- `kill_result` — subroutine
- `layout_name` — function
- `load_average` — function
- `log_start` — subroutine
- `map_pair_correlation` — subroutine
- `match_truth_maps` — subroutine
- `output_root` — function
- `prepare_states_cavgs_project` — subroutine
- `read_truth_states` — function
- `refine_flex_state` — subroutine
- `report_case` — subroutine
- `route_test` — subroutine
- `run_case` — subroutine
- `run_cavgs_route` — subroutine
- `run_cls_expansion_stage` — subroutine
- `run_cont_refine3D_bench` — subroutine
- `run_cont_refine3D_gate` — subroutine
- `run_flex_mode` — subroutine
- `run_flex_pca_blobs_gate` — subroutine
- `run_rec3D_backends_single` — subroutine
- `run_refine3D_states` — subroutine
- `run_solve3D` — subroutine
- `run_solve3D_addon_gate` — subroutine
- `run_state_refinements` — subroutine
- `run_states_cavgs_gate` — subroutine
- `run_states_cavgs_stage` — subroutine
- `run_states_route` — subroutine
- `scaled_ncpu` — function
- `selected_entries` — function
- `thread_args` — function
- `thread_list` — function
- `time_case` — subroutine
- `time_one` — function
- `validate_rec3D_pcg_fractional_updates` — subroutine
- `vector_correlation` — function
- `with_key` — function
- `write_mean_and_difference` — subroutine
- `write_report` — subroutine
- `write_state_weights_from_labels` — subroutine

Local symbols (enclosing procedure scope):
- `check_case::check_run` — subroutine
- `evaluate::docked_fname` — function
- `exec_test_mini_stream_quantitative::find_simulation` — function
- `exec_test_mini_stream_quantitative::nrecords` — function
- `exec_test_mini_stream_quantitative::require` — subroutine
- `exec_test_pcg_recon::band_rms` — subroutine
- `exec_test_pcg_recon::corr_of` — function
- `exec_test_profile::compare_profile_dirs` — subroutine
- `exec_test_simulate_particles::check_images_against_projections` — subroutine
- `exec_test_simulate_particles::check_mixture` — subroutine
- `exec_test_simulate_particles::check_oris` — subroutine
- `exec_test_simulate_particles::check_stack` — subroutine
- `exec_test_simulate_particles::gate_check_oris` — subroutine
- `exec_test_simulate_particles::gate_check_stack` — subroutine
- `exec_test_simulated_workflow::require` — subroutine
- `exec_test_simulated_workflow::return_to_stage_root` — subroutine
- `exec_test_simulated_workflow::run_solve3D_step` — subroutine
- `exec_test_simulated_workflow::update_project_path` — subroutine
- `exec_test_simulated_workflow::wrong_truth_control` — subroutine
- `load_average::getloadavg` — function
- `run_cavgs_route::simulate_half` — subroutine
- `run_cavgs_route::solve_halves_as_particles` — subroutine
- `run_cont_refine3D_bench::perturb` — subroutine
- `run_cont_refine3D_bench::pose_errors` — subroutine
- `run_cont_refine3D_bench::run_case` — subroutine
- `run_cont_refine3D_bench::run_regime` — subroutine
- `run_cont_refine3D_gate::check_no_polar` — subroutine
- `run_cont_refine3D_gate::enter_stage` — subroutine
- `run_cont_refine3D_gate::finish_gate` — subroutine
- `run_cont_refine3D_gate::gate_e25` — subroutine
- `run_cont_refine3D_gate::geodesic_angle` — function
- `run_cont_refine3D_gate::perturb` — subroutine
- `run_cont_refine3D_gate::pose_errors` — subroutine
- `run_cont_refine3D_gate::run_nu_handoff` — subroutine
- `run_cont_refine3D_gate::run_refine3D_cont` — subroutine
- `run_cont_refine3D_gate::run_refine3D_discrete` — subroutine
- `run_cont_refine3D_gate::sym_err` — function
- `run_flex_mode::job_exited_cleanly` — function
- `run_rec3D_backends_single::background_mean` — function
- `run_rec3D_backends_single::cline_tok` — function
- `run_rec3D_backends_single::gate_fail` — subroutine
- `run_rec3D_backends_single::ls_scale_profile` — subroutine
- `run_rec3D_backends_single::real2str_trim` — function
- `run_rec3D_backends_single::real_tok` — function
- `run_rec3D_backends_single::report_pair_fsc` — subroutine
- `run_rec3D_backends_single::safe_ratio` — function
- `run_solve3D_addon_gate::dock_and_compare` — subroutine
- `run_solve3D_addon_gate::shuffled_cohort_control` — subroutine
- `validate_rec3D_pcg_fractional_updates::accumulate_raw` — subroutine
- `validate_rec3D_pcg_fractional_updates::add_weighted` — subroutine
- `validate_rec3D_pcg_fractional_updates::build_blend` — subroutine
- `validate_rec3D_pcg_fractional_updates::collect_state` — subroutine
- `validate_rec3D_pcg_fractional_updates::collect_state_half` — subroutine
- `validate_rec3D_pcg_fractional_updates::finalize_and_solve` — subroutine
- `validate_rec3D_pcg_fractional_updates::frac_fname` — function
- `validate_rec3D_pcg_fractional_updates::load_weighted` — subroutine
- `validate_rec3D_pcg_fractional_updates::new_reduction` — subroutine
- `validate_rec3D_pcg_fractional_updates::require_raw` — subroutine
- `validate_rec3D_pcg_fractional_updates::split_complementary` — subroutine
- `validate_rec3D_pcg_fractional_updates::validate_half` — subroutine
- `validate_rec3D_pcg_fractional_updates::volume_rel_error` — function

---
## Module: simple_commanders_test_single

Files:
- `main/commanders/test/simple_commanders_test_single.f90`

Uses:
- `ieee_arithmetic`
- `simple_atoms`
- `simple_commanders_api`
- `simple_commanders_atoms`
- `simple_commanders_project_core`
- `simple_commanders_project_ptcl`
- `simple_commanders_reproject`
- `simple_commanders_sim`
- `simple_commanders_stkops`
- `simple_defs_atoms`
- `simple_imghead`
- `simple_nano_density`
- `simple_nano_density_tester`
- `simple_nano_species`
- `simple_nanoparticle_utils`
- `simple_refine3d_fnames`
- `simple_rnd`
- `simple_srch_sort_loc`
- `simple_test_gate`
- `simple_test_truth_metrics`
- `simple_test_utils`
- `single_commanders_nano2d`
- `single_commanders_nano3d`
- `single_commanders_trajectory`
- `single_commanders_trajectory_solve3d`

Public symbols:
- `commander_test_single_atoms_stats` — type
- `commander_test_single_workflow` — type
- `commander_test_species_discovery` — type
- `exec_test_single_atoms_stats` — subroutine
- `exec_test_single_workflow` — subroutine
- `exec_test_species_discovery` — subroutine
- `make_glc_trajectory_oris` — subroutine

Local symbols (enclosing procedure scope):
- `exec_test_single_atoms_stats::check_binary_lattice` — subroutine
- `exec_test_single_atoms_stats::check_coordination_statistics` — subroutine
- `exec_test_single_atoms_stats::check_line_count` — subroutine
- `exec_test_single_atoms_stats::check_model_products` — subroutine
- `exec_test_single_atoms_stats::check_statistics_files` — subroutine
- `exec_test_single_atoms_stats::check_volume` — subroutine
- `exec_test_single_atoms_stats::column_of` — subroutine
- `exec_test_single_atoms_stats::compare_atom_coordinates` — subroutine
- `exec_test_single_atoms_stats::compare_density_volumes` — subroutine
- `exec_test_single_atoms_stats::fail` — subroutine
- `exec_test_single_atoms_stats::shifted_control` — subroutine
- `exec_test_single_workflow::enter_workflow_stage` — subroutine
- `exec_test_single_workflow::group_pose_error` — function
- `exec_test_single_workflow::return_to_project_dir` — subroutine
- `exec_test_single_workflow::run_nano_refine` — subroutine
- `exec_test_single_workflow::run_trajectory_solve3D` — subroutine
- `exec_test_species_discovery::aperture_intensity` — function
- `exec_test_species_discovery::check_false_weak` — subroutine
- `exec_test_species_discovery::compare_products` — subroutine
- `exec_test_species_discovery::density_case` — subroutine
- `exec_test_species_discovery::discovery` — function
- `exec_test_species_discovery::elem_index` — function
- `exec_test_species_discovery::evaluate` — subroutine
- `exec_test_species_discovery::ground_truth` — subroutine
- `exec_test_species_discovery::interp` — function
- `exec_test_species_discovery::key_value` — function
- `exec_test_species_discovery::light_plain_recall` — subroutine
- `exec_test_species_discovery::make_fixture` — subroutine
- `exec_test_species_discovery::median_nn` — function
- `exec_test_species_discovery::outside_sdev` — function
- `exec_test_species_discovery::own_peak` — function
- `exec_test_species_discovery::profile` — function
- `exec_test_species_discovery::prune_eligibility` — subroutine
- `exec_test_species_discovery::radial_filter` — subroutine
- `exec_test_species_discovery::read_table` — subroutine
- `exec_test_species_discovery::read_transfer_smoothed` — subroutine
- `exec_test_species_discovery::recall_of` — function
- `exec_test_species_discovery::report_string` — function
- `exec_test_species_discovery::run_detect` — subroutine
- `exec_test_species_discovery::run_species` — subroutine
- `exec_test_species_discovery::species_run` — function
- `exec_test_species_discovery::width_rise` — function
- `exec_test_species_discovery::write_vol` — subroutine

---
## Module: simple_commanders_test_stream

Files:
- `main/commanders/test/simple_commanders_test_stream.f90`

Uses:
- `ieee_arithmetic`
- `iso_c_binding`
- `simple_commanders_api`
- `simple_commanders_sim`
- `simple_commanders_stream_p01_preprocess`
- `simple_oris`
- `simple_starfile_wrappers`
- `simple_syslib`
- `simple_test_gate`
- `simple_test_truth_metrics`
- `simple_ui`

Public symbols:
- `commander_test_stream_preproc` — type
- `exec_test_stream_preproc` — subroutine

Local symbols (enclosing procedure scope):
- `exec_test_stream_preproc::assert_output_file` — subroutine
- `exec_test_stream_preproc::check_worker_jobs` — subroutine
- `exec_test_stream_preproc::find_simulation` — function
- `exec_test_stream_preproc::max_motion_error` — function
- `exec_test_stream_preproc::require` — subroutine

---
## Module: simple_commanders_validate

Files:
- `main/commanders/simple/simple_commanders_validate.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_pick`
- `simple_commanders_preprocess`
- `simple_commanders_project_core`
- `simple_commanders_project_mov`
- `simple_commanders_solve2d`
- `simple_mini_stream_utils`

Public symbols:
- `commander_check_refpick` — type
- `commander_mini_stream` — type
- `exec_check_refpick` — subroutine
- `exec_mini_stream` — subroutine

Private symbols:
- `transfer_phase_fit_args` — subroutine

---
## Module: simple_commanders_volops

Files:
- `main/commanders/simple/simple_commanders_volops.f90`

Uses:
- `simple_atoms`
- `simple_butterworth`
- `simple_clustering_utils`
- `simple_commanders_api`
- `simple_commanders_reproject`
- `simple_dock_vols`
- `simple_halfmap_diagnostics`
- `simple_image_msk`
- `simple_imghead`
- `simple_procimgstk`
- `simple_projector`
- `simple_segmentation`
- `simple_simple_volinterp`
- `simple_strategy2d_utils`
- `simple_symanalyzer`
- `simple_vol_pproc_policy`
- `simple_volanalyzer`
- `simple_volcluster`
- `simple_volpft_symsrch`

Public symbols:
- `commander_centervol` — type
- `commander_dock_volpair` — type
- `commander_noisevol` — type
- `commander_postprocess` — type
- `commander_postprocess_pdb` — type
- `commander_symaxis_search` — type
- `commander_symmetrize_map` — type
- `commander_symmetry_test` — type
- `commander_volanalyze` — type
- `commander_volcluster` — type
- `commander_volops` — type
- `exec_centervol` — subroutine
- `exec_dock_volpair` — subroutine
- `exec_noisevol` — subroutine
- `exec_postprocess` — subroutine
- `exec_postprocess_pdb` — subroutine
- `exec_symaxis_search` — subroutine
- `exec_symmetrize_map` — subroutine
- `exec_symmetry_test` — subroutine
- `exec_volanalyze` — subroutine
- `exec_volcluster` — subroutine
- `exec_volops` — subroutine
- `postprocess_volume_from_files` — subroutine

---
## Module: simple_commanders_volops_tester

Files:
- `main/commanders/simple/simple_commanders_volops_tester.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_volops`
- `simple_defs`
- `simple_image`
- `simple_string`
- `simple_test_state_fixture`
- `simple_test_truth_metrics`
- `simple_test_utils`
- `simple_ui`

Public symbols:
- `run_all_volops_lib_tests` — subroutine

Private symbols:
- `map_corr` — function
- `mass_centre` — function
- `test_centervol` — subroutine
- `test_dock_volpair` — subroutine
- `test_symaxis_search_d2` — subroutine
- `write_blobs` — subroutine

---
## Module: simple_convergence

Files:
- `main/simple_convergence.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_cartft_pose_opt`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs`
- `simple_parameters`
- `simple_progress`

Public symbols:
- `cont_particle_stable` — function
- `convergence` — type

Private symbols:
- `anneal_tail_start` — function
- `append_stats` — subroutine
- `calc_continuous_inplane_stats` — subroutine
- `calc_pass_scores` — subroutine
- `calc_polish_stats` — subroutine
- `check_cont_conv` — function
- `check_conv2D` — function
- `check_conv3D` — function
- `get` — function
- `plot_projdirs` — subroutine
- `prev_iter_res` — function
- `print_iteration` — subroutine
- `print_verdict` — subroutine
- `record_res` — subroutine
- `res_two_back` — subroutine

---
## Module: simple_convergence_tester

Files:
- `main/simple_convergence_tester.f90`

Uses:
- `simple_cartft_pose_opt`
- `simple_cmdline`
- `simple_convergence`
- `simple_core_module_api`
- `simple_defs`
- `simple_defs_conv`
- `simple_defs_fname`
- `simple_parameters`
- `simple_test_utils`

Public symbols:
- `run_all_convergence_tests` — subroutine

Private symbols:
- `clear_samples` — subroutine
- `cont_converged` — function
- `fresh_cont_converged` — function
- `frozen_converged` — function
- `frozen_field` — subroutine
- `lpexit_converged` — function
- `make_cont_field` — subroutine
- `make_field` — subroutine
- `make_multistate_field` — subroutine
- `multistate_converged` — function
- `sample_half` — subroutine
- `test_cartesian_pass_scores` — subroutine
- `test_cartesian_pass_statistics` — subroutine
- `test_cont_particle_stable` — subroutine
- `test_cont_rule_full_sample` — subroutine
- `test_cont_rule_streak` — subroutine
- `test_cont_rule_tunable` — subroutine
- `test_lpexit_gain_guard` — subroutine
- `test_multistate_scan_exit` — subroutine
- `test_polar_pass_scores` — subroutine
- `test_posterior_anneal_tail` — subroutine
- `test_posterior_anneal_tail2D` — subroutine
- `test_posterior_assignment_bypass` — subroutine
- `test_posterior_exits_weighted_insertion` — subroutine

Local symbols (enclosing procedure scope):
- `test_posterior_anneal_tail2D::converged_field` — subroutine
- `test_posterior_exits_weighted_insertion::run_stage` — subroutine

---
## Module: simple_core_module_api

Files:
- `main/apis/simple_core_api.f90`

Uses:
- `iso_c_binding`
- `simple_binoris`
- `simple_chash`
- `simple_class_sample_io`
- `simple_defs`
- `simple_defs_conv`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_edges_sqwins`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_fileio`
- `simple_hash`
- `simple_imghead`
- `simple_is_check_assert`
- `simple_jiffys`
- `simple_kbinterpol`
- `simple_linalg`
- `simple_magic_boxes`
- `simple_map_reduce`
- `simple_math`
- `simple_math_ft`
- `simple_nrtxtfile`
- `simple_ori`
- `simple_ori_utils`
- `simple_oris`
- `simple_ran_tabu`
- `simple_rnd`
- `simple_srch_sort_loc`
- `simple_stat`
- `simple_string`
- `simple_string_utils`
- `simple_sym`
- `simple_syslib`
- `simple_timer`
- `simple_type_defs`

---
## Module: simple_corrmat

Files:
- `main/strategies/search/simple_corrmat.f90`

Uses:
- `simple_builder`
- `simple_defs`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`

Public symbols:
- `calc_inpl_invariant_cc_nomirr` — function

Private symbols:
- `calc_cartesian_corrmat_1` — subroutine
- `calc_cartesian_corrmat_2` — subroutine

---
## Module: simple_ctf

Files:
- `main/ctf/simple_ctf.f90`

Uses:
- `simple_core_module_api`
- `simple_math_ctf`

Public symbols:
- `ctf` — type

Private symbols:
- `apply_convention` — subroutine
- `constructor` — function
- `eval_canonical` — function
- `eval_df` — function
- `eval_phase` — function
- `get_ctfvars` — function
- `init` — subroutine
- `kV2wl` — function
- `nextrema` — function
- `SpaFreqSqAtNthZero` — function

---
## Module: simple_ctf_estimate_cost

Files:
- `main/ctf/simple_ctf_estimate_cost.f90`

Uses:
- `simple_core_module_api`
- `simple_ctf`
- `simple_image`
- `simple_memoize_ft_maps`
- `simple_opt_de`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `ctf_estimate_cost1D` — type
- `ctf_estimate_cost2D` — type
- `ctf_estimate_cost4Dcont` — type

Private symbols:
- `calc_cost2D` — function
- `cost1D` — function
- `cost2D` — function
- `cost2D_wrapper` — function
- `df_costfun` — subroutine
- `f_costfun` — function
- `fdf` — subroutine
- `fdf_costfun` — subroutine
- `init1D` — subroutine
- `init2D` — subroutine
- `init4Dcont` — subroutine
- `kill1D` — subroutine
- `kill2D` — subroutine
- `kill4Dcont` — subroutine
- `minimize` — subroutine
- `minimize4Dcont` — subroutine

---
## Module: simple_ctf_estimate_fit

Files:
- `main/ctf/simple_ctf_estimate_fit.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_core_module_api`
- `simple_ctf`
- `simple_ctf_estimate_cost`
- `simple_image`
- `simple_memoize_ft_maps`
- `simple_starfile_wrappers`
- `simple_timer`

Public symbols:
- `ctf_estimate_fit` — type

Private symbols:
- `calc_ctfres` — subroutine
- `calc_icefrac` — subroutine
- `ctf2pspecimg` — subroutine
- `dampen_pspec4srch` — subroutine
- `fit` — subroutine
- `fit_patches` — subroutine
- `fit_polynomial` — subroutine
- `ft2img` — subroutine
- `gen_resmsk` — subroutine
- `gen_roavspec1d` — subroutine
- `gen_tiles` — subroutine
- `get_astig` — function
- `get_ccfit` — function
- `get_ctfres` — function
- `get_icefrac` — function
- `get_parms` — subroutine
- `get_pspec` — subroutine
- `kill` — subroutine
- `mic2spec` — subroutine
- `new` — subroutine
- `norm_pspec` — subroutine
- `pix2poly` — subroutine
- `pix2polyvals` — subroutine
- `plot_ctf` — subroutine
- `poly` — function
- `read_doc` — subroutine
- `refine` — subroutine
- `srch` — subroutine
- `subtr_backgr` — subroutine
- `write_diagnostic` — subroutine
- `write_doc` — subroutine
- `write_star` — subroutine

Local symbols (enclosing procedure scope):
- `calc_ctfres::calc_frc` — subroutine
- `calc_ctfres::calc_pspec_stats` — subroutine
- `calc_ctfres::ctfres_shell` — function
- `calc_ctfres::gen_profiles` — subroutine
- `calc_ctfres::get_shell` — function
- `calc_ctfres::mask_central_disc` — subroutine
- `calc_ctfres::pix2logical` — function
- `calc_ctfres::rank_spec` — subroutine
- `calc_icefrac::mean_positive_amps` — function
- `pix2polyvals::poly2val` — function

---
## Module: simple_ctf_estimate_iter

Files:
- `main/ctf/simple_ctf_estimate_iter.f90`

Uses:
- `simple_core_module_api`
- `simple_ctf_estimate_fit`
- `simple_image`
- `simple_parameters`

Public symbols:
- `ctf_estimate_iter` — type

Private symbols:
- `gen_thumbnail` — subroutine
- `iterate` — subroutine
- `kill` — subroutine

---
## Module: simple_ctf_estimate_strategy

Files:
- `main/strategies/parallelization/simple_ctf_estimate_strategy.f90`

Uses:
- `simple_binoris_io`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_ctf_estimate_iter`
- `simple_parameters`
- `simple_qsys_env`
- `simple_sp_project`

Public symbols:
- `create_ctf_estimate_strategy` — function
- `ctf_estimate_distr_strategy` — type
- `ctf_estimate_inmem_strategy` — type
- `ctf_estimate_strategy` — type

Private symbols:
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_end_message` — function
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `endmsg_interface` — function
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_end_message` — function
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `set_ctf_estimate_defaults` — subroutine

---
## Module: simple_ctf_tester

Files:
- `main/ctf/simple_ctf_tester.f90`

Uses:
- `simple_ctf`
- `simple_ctf_estimate_fit`
- `simple_defs`
- `simple_fplane_grid`
- `simple_image`
- `simple_memoize_ft_maps`
- `simple_string_utils`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_ctf_tests` — subroutine

Private symbols:
- `angle` — function
- `closed_form` — function
- `ctf_at` — function
- `make_tfun` — function
- `spafreqsq` — function
- `test_closed_form_values` — subroutine
- `test_constructor_and_conventions` — subroutine
- `test_extrema` — subroutine
- `test_fast_kernel_and_ft2img` — subroutine
- `test_fit_kill_keeps_ft_maps` — subroutine
- `test_phase_shift_policy` — subroutine
- `test_restoration_contracts` — subroutine

---
## Module: simple_decay_funs

Files:
- `utils/math/simple_decay_funs.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `calc_nsampl_fromto` — function
- `calc_update_frac` — function
- `calc_update_frac_dyn` — function
- `cos_decay` — function
- `extremal_decay` — function
- `extremal_decay2D` — function
- `inv_cos_decay` — function
- `inv_nsampl_decay` — function

---
## Module: simple_decay_funs_tester

Files:
- `utils/math/simple_decay_funs_tester.f90`

Uses:
- `simple_decay_funs`
- `simple_defs`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_decay_funs_tests` — subroutine

Private symbols:
- `test_calc_nsampl_fromto` — subroutine
- `test_cos_decays` — subroutine
- `test_extremal_decays` — subroutine
- `test_inv_nsampl_decay` — subroutine
- `test_update_frac` — subroutine

---
## Module: simple_default_clines

Files:
- `defs/simple_default_clines.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_estimate_ssnr`

Public symbols:
- `set_automask2D_defaults` — subroutine
- `set_refine2D_defaults` — subroutine

---
## Module: simple_defs

Files:
- `defs/simple_defs.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`

---
## Module: simple_defs_atoms

Files:
- `defs/simple_defs_atoms.f90`

Uses:
- `simple_error`

Public symbols:
- `get_element_Z_and_radius` — subroutine
- `get_lattice_params` — subroutine

---
## Module: simple_defs_conv

Files:
- `defs/simple_defs_conv.f90`

---
## Module: simple_defs_environment

Files:
- `defs/simple_defs_environment.f90`

Uses:
- `iso_c_binding`

---
## Module: simple_defs_flex

Files:
- `defs/simple_defs_flex.f90`

Uses:
- `simple_defs`

---
## Module: simple_defs_fname

Files:
- `defs/simple_defs_fname.f90`

Uses:
- `iso_c_binding`

---
## Module: simple_defs_ori

Files:
- `defs/simple_defs_ori.f90`

Public symbols:
- `get_oriparam_flag` — function
- `get_oriparam_ind` — function
- `oriparam_is_spare` — function
- `oriparam_isthere` — function

---
## Module: simple_defs_stream

Files:
- `defs/simple_defs_stream.f90`

Public symbols:
- `scaled_dims` — type

---
## Module: simple_defs_string

Files:
- `defs/simple_defs_string.f90`

Uses:
- `iso_c_binding`

---
## Module: simple_discrete_stack_io

Files:
- `fileio/simple_discrete_stack_io.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_imgfile`
- `simple_imghead`

Public symbols:
- `dstack_io` — type

Private symbols:
- `check_index` — subroutine
- `close_stack` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `open_1` — subroutine
- `read` — subroutine

---
## Module: simple_dock_vols

Files:
- `main/volume/simple_dock_vols.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_projector`
- `simple_simple_volinterp`
- `simple_volpft_corrcalc`

Public symbols:
- `dock_vols` — type

Private symbols:
- `get_dock_info` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `rotate_target` — subroutine
- `rotpeak_interp` — subroutine
- `set_dock_info` — subroutine
- `set_ref` — subroutine
- `set_target` — subroutine
- `setup_srch_spaces` — subroutine
- `srch` — subroutine
- `srch_rots` — subroutine
- `srch_shift` — subroutine

---
## Module: simple_dose_series_manifest

Files:
- `main/motion/simple_dose_series_manifest.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_defs`
- `simple_halfmap_diagnostics`
- `simple_imghead`
- `simple_string`
- `simple_syslib`

Public symbols:
- `dose_series_manifest` — type

Private symbols:
- `add_window` — subroutine
- `check_header` — function
- `check_support` — function
- `close_value` — function
- `dose_window` — type
- `get_geometry` — subroutine
- `get_nwindows` — function
- `get_window` — subroutine
- `header_slot` — function
- `kill` — subroutine
- `new` — subroutine
- `overlap_predecessor` — function
- `parent_directory` — function
- `parse_int` — function
- `parse_real` — function
- `portable_path` — function
- `read` — function
- `read_line` — subroutine
- `remove_temporary` — subroutine
- `resolve_path` — function
- `select_windows` — function
- `split_tokens` — subroutine
- `valid_path` — function
- `validate` — function
- `validate_files` — function
- `write` — function

---
## Module: simple_dose_series_manifest_tester

Files:
- `main/motion/simple_dose_series_manifest_tester.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_defs`
- `simple_dose_series_manifest`
- `simple_image`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_dose_series_manifest_tests` — subroutine

Private symbols:
- `make_maps` — subroutine
- `make_series` — subroutine
- `test_bad_headers` — subroutine
- `test_bad_metadata` — subroutine
- `test_final_line` — subroutine
- `test_publication` — subroutine
- `test_roundtrip` — subroutine
- `write_case` — subroutine

---
## Module: simple_edges_sqwins

Files:
- `main/interp/simple_edges_sqwins.f90`

Uses:
- `simple_defs`
- `simple_error`

Public symbols:
- `cosedge_1` — function
- `cosedge_2` — function
- `cosedge_inner_1` — function
- `cosedge_inner_2` — function
- `cosedge_r2_2d` — function
- `cosedge_r2_3d` — function
- `hardedge_1` — function
- `hardedge_2` — function
- `hardedge_3` — function
- `hardedge_4` — function
- `hardedge_inner_1` — function
- `hardedge_inner_2` — function
- `hardedge_inner_3` — function
- `hardedge_inner_4` — function
- `hardedge_r2_2d` — function
- `hardedge_r2_3d` — function
- `sqwin_1d_1` — subroutine
- `sqwin_1d_2` — subroutine
- `sqwin_2d_1` — subroutine
- `sqwin_2d_2` — subroutine
- `sqwin_3d_1` — subroutine
- `sqwin_3d_2` — subroutine

---
## Module: simple_eer_factory

Files:
- `fileio/simple_eer_factory.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_tifflib`

Public symbols:
- `eer_decoder` — type

Private symbols:
- `decode` — subroutine
- `decode_frame_7bit_16k` — subroutine
- `decode_frame_7bit_4k` — subroutine
- `decode_frame_7bit_8k` — subroutine
- `decode_frame_8bit_16k` — subroutine
- `decode_frame_8bit_4k` — subroutine
- `decode_frame_8bit_8k` — subroutine
- `decode_frames` — subroutine
- `get_ldim` — function
- `get_nframes` — function
- `get_smpd` — function
- `get_smpd_out` — function
- `kill` — subroutine
- `new` — subroutine
- `prep_gainref` — subroutine
- `read` — subroutine
- `set_dims` — subroutine

---
## Module: simple_eo_partition

Files:
- `main/project/simple_eo_partition.f90`

Uses:
- `iso_fortran_env`
- `simple_core_module_api`
- `simple_oris`
- `simple_sp_project`
- `simple_type_defs`

Public symbols:
- `ensure_field_halves` — subroutine
- `eo_partition_status` — function
- `eo_transaction` — subroutine
- `repartition_selected_halves` — subroutine

Private symbols:
- `adopt_legacy` — subroutine
- `all_marked` — function
- `any_marked` — function
- `check_marked` — subroutine
- `copy_halves` — subroutine
- `eo_shuffle` — subroutine
- `field_layout` — subroutine
- `index_half` — function
- `label_by_index_all` — subroutine
- `next_random` — function
- `rows_match` — function
- `run_transaction` — subroutine
- `set_label` — subroutine
- `stable_sort_desc` — subroutine
- `stack_has_rows` — function

Local symbols (enclosing procedure scope):
- `repartition_selected_halves::assign_half` — subroutine

---
## Module: simple_eo_partition_tester

Files:
- `main/project/simple_eo_partition_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_defs`
- `simple_eo_partition`
- `simple_forked_process`
- `simple_oris`
- `simple_projfile_utils`
- `simple_rnd`
- `simple_sp_project`
- `simple_string`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_eo_partition_tests` — subroutine

Private symbols:
- `active_halves` — subroutine
- `append_rows` — subroutine
- `fill_particle_row` — subroutine
- `fill_stack_row` — subroutine
- `groups_whole` — function
- `make_groups_project` — subroutine
- `refuse_partition` — subroutine
- `refusing_partition` — type
- `test_appended_stacks` — subroutine
- `test_class_average_route` — subroutine
- `test_earlier_release` — subroutine
- `test_inconsistent` — subroutine
- `test_isolated_stream_and_seed` — subroutine
- `test_micrograph_groups` — subroutine
- `test_no_split_by_size` — subroutine
- `test_per_index_route` — subroutine
- `test_project_file_and_merge` — subroutine
- `test_refusals` — subroutine
- `test_selected_score_pairs` — subroutine
- `test_too_few_groups_status` — subroutine

---
## Module: simple_epu_xml

Files:
- `fileio/simple_epu_xml.f90`

Public symbols:
- `read_epu_beamshift` — subroutine

Private symbols:
- `element_content` — function
- `element_real` — subroutine
- `read_text` — function

---
## Module: simple_epu_xml_tester

Files:
- `fileio/simple_epu_xml_tester.f90`

Uses:
- `simple_epu_xml`
- `simple_test_utils`

Public symbols:
- `run_all_epu_xml_tests` — subroutine

Private symbols:
- `delete_doc` — subroutine
- `test_indented_document` — subroutine
- `test_longer_name_first` — subroutine
- `test_missing_beamshift` — subroutine
- `test_missing_component` — subroutine
- `test_missing_file` — subroutine
- `test_single_line_document` — subroutine
- `write_doc` — subroutine

---
## Module: simple_error

Files:
- `fileio/simple_error.f90`

Uses:
- `iso_fortran_env`
- `simple_defs`

Public symbols:
- `simple_error_check` — subroutine
- `simple_exception_hard` — subroutine
- `simple_exception_opt` — subroutine

---
## Module: simple_estimate_ssnr

Files:
- `utils/filter/simple_estimate_ssnr.f90`

Uses:
- `simple_defs`
- `simple_defs_conv`
- `simple_defs_stream`
- `simple_fileio`
- `simple_magic_boxes`
- `simple_math_ft`
- `simple_srch_sort_loc`
- `simple_string_utils`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `bfac_cap_filter` — function
- `calc_dose_weights` — subroutine
- `fsc2cref` — function
- `fsc2optlp` — function
- `fsc2optlp_sub` — subroutine
- `fsc2ssnr` — function
- `fsc2wiener_regularizer` — function
- `gaussian_filter` — subroutine
- `get_resolution` — subroutine
- `get_resolution_at_fsc` — subroutine
- `lowpass_from_klim` — subroutine
- `lpstages` — subroutine
- `lpstages_aggressive` — subroutine
- `lpstages_fast` — subroutine
- `lpstages_setlims` — subroutine
- `mskdiam2lplimits` — subroutine
- `mskdiam2streamresthreshold` — subroutine
- `ssnr2fsc` — function
- `ssnr2optlp` — function

Local symbols (enclosing procedure scope):
- `lpstages::apply_trslim_override` — subroutine
- `lpstages::calc_lpinfo` — subroutine
- `lpstages::calc_scaleinfo` — subroutine
- `lpstages::print_lpinfo` — subroutine
- `lpstages::print_scaleinfo` — subroutine
- `lpstages_fast::calc_scaleinfo` — subroutine

---
## Module: simple_euclid_sigma2

Files:
- `main/sigma2/simple_euclid_sigma2.f90`

Uses:
- `iso_fortran_env`
- `simple_cartft_calc`
- `simple_core_module_api`
- `simple_fplane_grid`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_sigma2_state`
- `simple_sigma2_state_file`
- `simple_syslib`

Public symbols:
- `euclid_sigma2` — type
- `sigma2_commit_iter` — function
- `sigma2_half` — function

Private symbols:
- `calc_sigma2_cart` — subroutine
- `calc_sigma2_cart_residual` — subroutine
- `calc_sigma2_pftc` — subroutine
- `check_record` — subroutine
- `get_kfromto` — function
- `get_sigma2_part` — function
- `kill` — subroutine
- `new_cart` — subroutine
- `new_pftc` — subroutine
- `prep_fplane_grid` — subroutine
- `read_global` — subroutine
- `read_part` — subroutine
- `read_pending_range` — subroutine
- `set_kfromto` — subroutine
- `write_sigma2` — subroutine

---
## Module: simple_eul_prob_tab

Files:
- `main/strategies/search/probabilistic/simple_eul_prob_tab.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `omp_lib`
- `simple_builder`
- `simple_defs`
- `simple_eul_prob_tab_utils`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_polarft_calc`
- `simple_rnd`
- `simple_state_weight_set`
- `simple_type_defs`

Public symbols:
- `eul_prob_tab` — type

Private symbols:
- `assign_batch` — subroutine
- `assign_candidate` — subroutine
- `assign_rows` — subroutine
- `close_assignment` — subroutine
- `entry_candidate` — function
- `fill_tab_range` — subroutine
- `full_to_compact_ref` — function
- `get_fill_profile` — function
- `get_row` — subroutine
- `initialize_storage` — subroutine
- `keep_batch_rows` — subroutine
- `kill` — subroutine
- `max_refined` — function
- `new_common` — subroutine
- `new_pipeline` — subroutine
- `open_assignment` — subroutine
- `ref_full` — function
- `ref_proj` — function
- `ref_state` — function
- `set_insert_cone` — subroutine
- `set_iteration_model` — subroutine
- `set_model_from_field` — subroutine
- `write_posterior_diag` — subroutine

Local symbols (enclosing procedure scope):
- `assign_rows::count_row_validity` — subroutine
- `assign_rows::select_extras` — subroutine
- `assign_rows::valid_score` — function
- `fill_tab_range::add_counters` — subroutine
- `fill_tab_range::deterministic_best` — subroutine
- `fill_tab_range::estimate_shift_seed` — subroutine
- `fill_tab_range::get_particle_context` — subroutine
- `fill_tab_range::posterior_diag` — subroutine
- `fill_tab_range::process_particle_with_map` — subroutine
- `fill_tab_range::process_particle_with_shift` — subroutine
- `fill_tab_range::process_particle_without_shift` — subroutine
- `fill_tab_range::rank_in_row` — function
- `fill_tab_range::refine_best_refs` — subroutine
- `fill_tab_range::score_ref` — subroutine
- `fill_tab_range::score_ref_bank` — subroutine
- `fill_tab_range::spread_about` — function
- `new_common::weighted_states_exist` — subroutine

---
## Module: simple_eul_prob_tab2D

Files:
- `main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_builder`
- `simple_decay_funs`
- `simple_eul_prob_tab_utils`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_polarft_calc`
- `simple_rnd`
- `simple_srch_sort_loc`
- `simple_type_defs`

Public symbols:
- `eul_prob_tab2D` — type

Private symbols:
- `alloc_seed_bank2D` — subroutine
- `assign_batch` — subroutine
- `assign_rows_dense` — subroutine
- `assign_rows_sparse` — subroutine
- `batch_rows_from_buffers` — subroutine
- `batch_store_from_buffers` — subroutine
- `close_assignment` — subroutine
- `dealloc_eval2D_sparse_ws` — subroutine
- `eval2D_sparse_ws` — type
- `fill_tab_prob_neigh_range` — subroutine
- `fill_tab_range` — subroutine
- `finalize_candidate_assignment` — subroutine
- `get_fill_profile` — function
- `init_eval2D_sparse_ws` — subroutine
- `initialize_storage` — subroutine
- `kill` — subroutine
- `memoize_seed_bank2D` — subroutine
- `new_common` — subroutine
- `new_disk_objects2D` — subroutine
- `new_pipeline` — subroutine
- `open_assignment` — subroutine
- `record_fill_profile2D` — subroutine
- `record_sparse_eval` — subroutine
- `row_candidate` — function
- `score_seeded_class` — subroutine
- `seed_start` — function
- `set_iteration_model` — subroutine
- `set_model_from_field` — subroutine
- `setup_shiftmap2D` — subroutine

Local symbols (enclosing procedure scope):
- `assign_rows_sparse::count_npeaks` — function
- `assign_rows_sparse::gather_row` — subroutine
- `assign_rows_sparse::valid_sparse_candidate` — function
- `fill_tab_prob_neigh_range::process_particle` — subroutine
- `fill_tab_prob_neigh_range::refine_best_classes` — subroutine
- `fill_tab_prob_neigh_range::score_class` — subroutine
- `record_fill_profile2D::add_counters` — subroutine

---
## Module: simple_eul_prob_tab2D_tester

Files:
- `main/strategies/search/probabilistic/simple_eul_prob_tab2D_tester.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab2d`
- `simple_eul_prob_tab_utils`
- `simple_parameters`
- `simple_rnd`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_eul_prob_tab2D_tests` — subroutine

Private symbols:
- `assert_assignment` — subroutine
- `assert_shifts` — subroutine
- `assign_batch_of` — subroutine
- `set_candidate` — subroutine
- `test_dense_batch` — subroutine
- `test_dense_posterior_assignment` — subroutine
- `test_sparse_batch` — subroutine
- `test_sparse_posterior_assignment` — subroutine
- `test_sparse_shared_best_class` — subroutine

Local symbols (enclosing procedure scope):
- `test_dense_posterior_assignment::assign_dense` — subroutine
- `test_sparse_posterior_assignment::assign_sparse` — subroutine

---
## Module: simple_eul_prob_tab_neigh

Files:
- `main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90`

Uses:
- `iso_fortran_env`
- `omp_lib`
- `simple_builder`
- `simple_decay_funs`
- `simple_defs`
- `simple_eul_prob_tab`
- `simple_eul_prob_tab_utils`
- `simple_ori`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_polarft_calc`
- `simple_rnd`

Public symbols:
- `eul_prob_tab_neigh` — type

Private symbols:
- `alloc_coarse_ws` — subroutine
- `alloc_eval_stoch_ws` — subroutine
- `alloc_eval_subspace_ws` — subroutine
- `alloc_seed_bank` — subroutine
- `assign_batch_neigh` — subroutine
- `assign_rows_neigh` — subroutine
- `batch_store_from_buffers` — subroutine
- `close_assignment_neigh` — subroutine
- `coarse_search_ws` — type
- `dealloc_coarse_ws` — subroutine
- `dealloc_eval_ws` — subroutine
- `estimate_shift_seed` — subroutine
- `eval_ws` — type
- `fill_tab_neigh_range` — subroutine
- `fill_tab_stoch_range` — subroutine
- `fill_tab_subspace_range` — subroutine
- `get_particle_context` — subroutine
- `init_eval_stoch_ws` — subroutine
- `init_eval_subspace_ws` — subroutine
- `kill_neigh` — subroutine
- `memoize_seed_bank` — subroutine
- `new_neigh_pipeline` — subroutine
- `open_assignment_neigh` — subroutine
- `record_neigh_fill_profile` — subroutine
- `record_sparse_eval` — subroutine
- `refine_best_neighbors` — subroutine
- `score_seeded_ref` — subroutine
- `set_model_from_field_neigh` — subroutine
- `setup_shiftmap_neigh` — subroutine
- `validate_neigh_mode` — subroutine

Local symbols (enclosing procedure scope):
- `assign_rows_neigh::assign_particles_posterior` — subroutine
- `assign_rows_neigh::find_candidate_position` — function
- `assign_rows_neigh::gather_state_row` — subroutine
- `assign_rows_neigh::has_valid_candidate` — function
- `assign_rows_neigh::particle_evaluated_count` — function
- `assign_rows_neigh::search_frac` — function
- `assign_rows_neigh::seed_fallback_if_empty` — subroutine
- `assign_rows_neigh::select_extras` — subroutine
- `assign_rows_neigh::valid_compact_ref` — function
- `fill_tab_stoch_range::calc_previous_corr` — subroutine
- `fill_tab_stoch_range::evaluate_direct_stochastic_refs` — subroutine
- `fill_tab_stoch_range::process_particle_stoch` — subroutine
- `fill_tab_stoch_range::score_direct_ref` — subroutine
- `fill_tab_subspace_range::add_pooled_subspace` — subroutine
- `fill_tab_subspace_range::build_geometric_neighborhood` — subroutine
- `fill_tab_subspace_range::build_pooled_neighborhood` — subroutine
- `fill_tab_subspace_range::consider_peak_subspace` — subroutine
- `fill_tab_subspace_range::evaluate_neighborhood` — subroutine
- `fill_tab_subspace_range::find_peak_subspaces` — subroutine
- `fill_tab_subspace_range::process_particle_subspace` — subroutine
- `fill_tab_subspace_range::score_subspace_ref` — subroutine

---
## Module: simple_eul_prob_tab_tester

Files:
- `main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab`
- `simple_eul_prob_tab_neigh`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_rnd`
- `simple_sym`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_eul_prob_tab_tests` — subroutine

Private symbols:
- `assigned_ref` — function
- `check_posterior_assignment` — subroutine
- `make_dense_rows` — subroutine
- `test_bounded_and_power_samplers` — subroutine
- `test_draw_posterior_row` — subroutine
- `test_insert_cone` — subroutine
- `test_insert_extras_multistate` — subroutine
- `test_insert_pose_loss` — subroutine
- `test_insert_table` — subroutine
- `test_neff_readout` — subroutine
- `test_neff_report` — subroutine
- `test_posterior_row_stats` — subroutine
- `test_prob_beta_at` — subroutine
- `test_record_seed_shift` — subroutine
- `test_ref_assign_neigh_posterior` — subroutine
- `test_ref_assign_posterior` — subroutine
- `test_screen_width_shortlist` — subroutine
- `test_top_posterior_extras` — subroutine
- `unrefined` — function

Local symbols (enclosing procedure scope):
- `check_posterior_assignment::assign_rows_all` — subroutine
- `test_bounded_and_power_samplers::dist_at` — function
- `test_bounded_and_power_samplers::heap_sift_down` — subroutine
- `test_bounded_and_power_samplers::reference_bounded` — subroutine
- `test_bounded_and_power_samplers::reference_power` — subroutine
- `test_bounded_and_power_samplers::reference_select` — subroutine
- `test_draw_posterior_row::oracle_cdf` — subroutine
- `test_draw_posterior_row::valid` — function
- `test_insert_cone::geodesic_deg` — function
- `test_insert_extras_multistate::assign_master` — subroutine
- `test_ref_assign_neigh_posterior::assign_store` — subroutine
- `test_ref_assign_neigh_posterior::set_store` — subroutine

---
## Module: simple_eul_prob_tab_utils

Files:
- `main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_defs`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_fileio`
- `simple_linalg`
- `simple_math`
- `simple_math_ft`
- `simple_ori_utils`
- `simple_oris`
- `simple_refine3d_fnames`
- `simple_rnd`
- `simple_string`
- `simple_string_utils`
- `simple_sym`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `add_prob_readout` — subroutine
- `calc_athres` — function
- `calc_num2sample` — subroutine
- `disk_radius_px` — function
- `draw_posterior_row` — subroutine
- `eulprob_corr_switch` — function
- `eulprob_dist_switch` — function
- `insert_cone_from_fsc` — function
- `insert_pose_loss` — function
- `keyed_candidate_order` — subroutine
- `map_prior_weight` — function
- `merge_prob_readouts` — subroutine
- `next_shift_width` — function
- `posterior_row_stats` — subroutine
- `posterior_temperature` — function
- `print_disk_set` — subroutine
- `prob_beta_at` — function
- `prob_candidate` — type
- `prob_candidate_buffer` — type
- `prob_candidate_store` — type
- `prob_dense_entry` — type
- `prob_fill_profile` — type
- `prob_insert_cone` — type
- `prob_insert_pose` — type
- `prob_insert_poses` — type
- `prob_insert_table` — type
- `prob_neff_readout` — type
- `prob_posterior_stat` — type
- `prob_readout` — type
- `prob_readout_fname` — function
- `prob_refined_entry` — type
- `prob_rnd_scope` — function
- `read_prob_readout` — subroutine
- `record_seed_shift` — subroutine
- `record_shift_state` — subroutine
- `report_prob_readout` — subroutine
- `sample_bounded_dists` — subroutine
- `sample_power_dists` — subroutine
- `screen_pool_size` — function
- `shift_box_radius` — function
- `shift_needs_disk` — function
- `shift_prior_sigma2` — function
- `temperature_or_argmin` — function
- `top_posterior_extras` — subroutine
- `update_screen_best` — subroutine
- `write_fill_profile_bench` — subroutine
- `write_prob_readout` — subroutine

Private symbols:
- `angle_sampling_1` — function
- `angle_sampling_2` — function
- `append_candidate` — subroutine
- `append_or_replace_candidate` — subroutine
- `get_buffer_candidate_inpl` — function
- `insert_cone_get_angle` — function
- `insert_cone_has_space` — function
- `insert_cone_is_open` — function
- `insert_cone_kill` — subroutine
- `insert_cone_new` — subroutine
- `insert_cone_set_space` — subroutine
- `insert_cone_within` — function
- `insert_poses_has` — function
- `insert_poses_kill` — subroutine
- `insert_poses_setup_range` — subroutine
- `insert_report_counts` — subroutine
- `insert_table_add` — subroutine
- `insert_table_counts` — subroutine
- `insert_table_get_committed_weight` — function
- `insert_table_get_extra` — function
- `insert_table_get_frac_capped` — function
- `insert_table_get_kmax` — function
- `insert_table_get_mean_eligible` — function
- `insert_table_get_nextra` — function
- `insert_table_get_nptcls` — function
- `insert_table_get_pind` — function
- `insert_table_kill` — subroutine
- `insert_table_new` — subroutine
- `insert_table_report` — subroutine
- `kill_candidate_buffer` — subroutine
- `kill_candidate_store` — subroutine
- `maxheap_sift_down` — subroutine
- `neff_readout_kill` — subroutine
- `neff_readout_mean_neff` — function
- `neff_readout_new` — subroutine
- `neff_readout_nkept` — function
- `neff_readout_offer` — subroutine
- `neff_readout_report` — subroutine
- `neff_readout_sum_neff` — function
- `neff_readout_wanted` — function
- `neff_report_sums` — subroutine
- `new_fixed_candidate_store` — subroutine
- `new_ragged_candidate_store` — subroutine
- `select_lowest_dists` — subroutine

Local symbols (enclosing procedure scope):
- `read_prob_readout::read_dp` — subroutine
- `read_prob_readout::read_i8` — subroutine
- `read_prob_readout::read_int` — subroutine
- `read_prob_readout::read_key` — subroutine
- `read_prob_readout::read_log` — subroutine
- `write_fill_profile_bench::key` — function

---
## Module: simple_exec_api

Files:
- `main/apis/simple_exec_api.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_exec_cavgproc`
- `simple_exec_denoise`
- `simple_exec_dock`
- `simple_exec_filter`
- `simple_exec_helpers`
- `simple_exec_image`
- `simple_exec_mask`
- `simple_exec_ori`
- `simple_exec_other`
- `simple_exec_preproc`
- `simple_exec_print`
- `simple_exec_project`
- `simple_exec_refine2d`
- `simple_exec_refine3d`
- `simple_exec_res`
- `simple_exec_sieve`
- `simple_exec_sim`
- `simple_exec_solve3d`
- `simple_exec_sym`
- `simple_exec_validate`
- `simple_exec_volume`
- `simple_jiffys`
- `simple_ui`
- `simple_ui_program`

---
## Module: simple_exec_cavgproc

Files:
- `main/exec/simple_exec_cavgproc.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_cavgs`
- `simple_commanders_stkops`

Public symbols:
- `exec_cavgproc_commander` — subroutine

---
## Module: simple_exec_denoise

Files:
- `main/exec/simple_exec_denoise.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_denoise`
- `simple_commanders_flex_pca`

Public symbols:
- `exec_denoise_commander` — subroutine

---
## Module: simple_exec_dock

Files:
- `main/exec/simple_exec_dock.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_volops`

Public symbols:
- `exec_dock_commander` — subroutine

---
## Module: simple_exec_filter

Files:
- `main/exec/simple_exec_filter.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_halfmap_enhance`
- `simple_commanders_imgops`
- `simple_commanders_postprocess_nu`
- `simple_commanders_resolest`
- `simple_commanders_volops`

Public symbols:
- `exec_filter_commander` — subroutine

---
## Module: simple_exec_helpers

Files:
- `main/exec/simple_exec_helpers.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_parameters`
- `simple_qsys_env`
- `simple_sp_project`

Public symbols:
- `async_exec` — subroutine
- `copy_project_file_to_root_dir` — subroutine
- `exec_screen` — subroutine
- `gen_exec_cmd` — function
- `restarted_exec` — subroutine
- `script_exec` — subroutine
- `set_master_num_threads` — subroutine
- `set_shmem_flag` — function
- `update_job_descriptions_in_project` — subroutine

---
## Module: simple_exec_image

Files:
- `main/exec/simple_exec_image.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_imgops`
- `simple_commanders_imgproc`
- `simple_commanders_stkops`

Public symbols:
- `exec_image_commander` — subroutine

---
## Module: simple_exec_mask

Files:
- `main/exec/simple_exec_mask.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_mask`

Public symbols:
- `exec_mask_commander` — subroutine

---
## Module: simple_exec_ori

Files:
- `main/exec/simple_exec_ori.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_ori`

Public symbols:
- `exec_ori_commander` — subroutine

---
## Module: simple_exec_other

Files:
- `main/exec/simple_exec_other.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_atoms`
- `simple_commanders_distr`
- `simple_commanders_project_ptcl`

Public symbols:
- `exec_other_commander` — subroutine

---
## Module: simple_exec_preproc

Files:
- `main/exec/simple_exec_preproc.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_motion`
- `simple_commanders_pick`
- `simple_commanders_preprocess`
- `simple_commanders_starproject`

Public symbols:
- `exec_preproc_commander` — subroutine

---
## Module: simple_exec_print

Files:
- `main/exec/simple_exec_print.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_checks`
- `simple_commanders_misc`

Public symbols:
- `exec_print_commander` — subroutine

---
## Module: simple_exec_project

Files:
- `main/exec/simple_exec_project.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_mic_subset`
- `simple_commanders_project_cls`
- `simple_commanders_project_core`
- `simple_commanders_project_mov`
- `simple_commanders_project_ptcl`
- `simple_commanders_ptcl_subset`
- `simple_commanders_relion`
- `simple_commanders_starproject`

Public symbols:
- `exec_project_commander` — subroutine

---
## Module: simple_exec_refine2D

Files:
- `main/exec/simple_exec_refine2D.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_cavgs`
- `simple_commanders_mkcavgs`
- `simple_commanders_project_cls`
- `simple_commanders_refine2d`
- `simple_commanders_solve2d`
- `simple_exec_helpers`
- `simple_string`

Public symbols:
- `exec_refine2D_commander` — subroutine

---
## Module: simple_exec_refine3D

Files:
- `main/exec/simple_exec_refine3D.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_mask`
- `simple_commanders_rec`
- `simple_commanders_refine3d`
- `simple_commanders_volops`
- `simple_exec_helpers`
- `simple_string`

Public symbols:
- `exec_refine3D_commander` — subroutine

---
## Module: simple_exec_res

Files:
- `main/exec/simple_exec_res.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_resolest`

Public symbols:
- `exec_res_commander` — subroutine

---
## Module: simple_exec_sieve

Files:
- `main/exec/simple_exec_sieve.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_sieve`

Public symbols:
- `exec_sieve_commander` — subroutine

---
## Module: simple_exec_sim

Files:
- `main/exec/simple_exec_sim.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_atoms`
- `simple_commanders_sim`

Public symbols:
- `exec_sim_commander` — subroutine

---
## Module: simple_exec_solve3D

Files:
- `main/exec/simple_exec_solve3D.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_resolest`
- `simple_commanders_solve3d`
- `simple_commanders_solve3d_chunks`
- `simple_commanders_solve3d_pyramid`
- `simple_commanders_volops`
- `simple_exec_helpers`
- `simple_string`

Public symbols:
- `exec_solve3D_commander` — subroutine

---
## Module: simple_exec_sym

Files:
- `main/exec/simple_exec_sym.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_volops`

Public symbols:
- `exec_sym_commander` — subroutine

---
## Module: simple_exec_validate

Files:
- `main/exec/simple_exec_validate.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_atoms`
- `simple_commanders_validate`

Public symbols:
- `exec_validate_commander` — subroutine

---
## Module: simple_exec_volume

Files:
- `main/exec/simple_exec_volume.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_atoms`
- `simple_commanders_reproject`
- `simple_commanders_volops`

Public symbols:
- `exec_volume_commander` — subroutine

---
## Module: simple_external_reference_pose_initialization

Files:
- `main/simple_external_reference_pose_initialization.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_euclid`
- `simple_estimate_ssnr`
- `simple_refine3d_fnames`

Public symbols:
- `initialize_poses_against_external_references` — subroutine

Private symbols:
- `prepare_pose_initialization_coverage` — subroutine
- `validate_pose_initialized_states` — subroutine

---
## Module: simple_extract_strategy

Files:
- `main/strategies/parallelization/simple_extract_strategy.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_api`
- `simple_ctf_estimate_fit`
- `simple_parameters`
- `simple_particle_extractor`
- `simple_qsys_env`
- `simple_sp_project`

Public symbols:
- `create_extract_strategy` — function
- `extract_distr_strategy` — type
- `extract_inmem_strategy` — type
- `extract_strategy` — type

Private symbols:
- `apply_defaults_interface` — subroutine
- `boxfile_from_mic_impl` — function
- `cleanup_interface` — subroutine
- `distr_apply_defaults` — subroutine
- `distr_cleanup` — subroutine
- `distr_end_message` — function
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `endmsg_interface` — function
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_apply_defaults` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_end_message` — function
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `set_extract_defaults` — subroutine

Local symbols (enclosing procedure scope):
- `inmem_execute::box_inside` — function
- `inmem_execute::killimgbatch` — subroutine
- `inmem_execute::prepimgbatch` — subroutine

---
## Module: simple_fftw3

Files:
- `defs/simple_fftw3.f90`

Uses:
- `iso_c_binding`

Public symbols:
- `fftw_alignment_of` — function
- `fftw_alloc_complex` — function
- `fftw_alloc_real` — function
- `fftw_cleanup` — subroutine
- `fftw_cleanup_threads` — subroutine
- `fftw_cost` — function
- `fftw_destroy_plan` — subroutine
- `fftw_estimate_cost` — function
- `fftw_execute_dft` — subroutine
- `fftw_execute_dft_c2r` — subroutine
- `fftw_execute_dft_r2c` — subroutine
- `fftw_execute_r2r` — subroutine
- `fftw_execute_split_dft` — subroutine
- `fftw_execute_split_dft_c2r` — subroutine
- `fftw_execute_split_dft_r2c` — subroutine
- `fftw_export_wisdom` — subroutine
- `fftw_export_wisdom_to_file` — subroutine
- `fftw_export_wisdom_to_filename` — function
- `fftw_export_wisdom_to_string` — function
- `fftw_flops` — subroutine
- `fftw_forget_wisdom` — subroutine
- `fftw_fprint_plan` — subroutine
- `fftw_free` — subroutine
- `fftw_free_str` — subroutine
- `fftw_import_system_wisdom` — function
- `fftw_import_wisdom` — function
- `fftw_import_wisdom_from_file` — function
- `fftw_import_wisdom_from_filename` — function
- `fftw_import_wisdom_from_string` — function
- `fftw_init_threads` — function
- `fftw_iodim` — type
- `fftw_iodim64` — type
- `fftw_malloc` — function
- `fftw_plan_dft` — function
- `fftw_plan_dft_1d` — function
- `fftw_plan_dft_2d` — function
- `fftw_plan_dft_3d` — function
- `fftw_plan_dft_c2r` — function
- `fftw_plan_dft_c2r_1d` — function
- `fftw_plan_dft_c2r_2d` — function
- `fftw_plan_dft_c2r_3d` — function
- `fftw_plan_dft_r2c` — function
- `fftw_plan_dft_r2c_1d` — function
- `fftw_plan_dft_r2c_2d` — function
- `fftw_plan_dft_r2c_3d` — function
- `fftw_plan_guru64_dft` — function
- `fftw_plan_guru64_dft_c2r` — function
- `fftw_plan_guru64_dft_r2c` — function
- `fftw_plan_guru64_r2r` — function
- `fftw_plan_guru64_split_dft` — function
- `fftw_plan_guru64_split_dft_c2r` — function
- `fftw_plan_guru64_split_dft_r2c` — function
- `fftw_plan_guru_dft` — function
- `fftw_plan_guru_dft_c2r` — function
- `fftw_plan_guru_dft_r2c` — function
- `fftw_plan_guru_r2r` — function
- `fftw_plan_guru_split_dft` — function
- `fftw_plan_guru_split_dft_c2r` — function
- `fftw_plan_guru_split_dft_r2c` — function
- `fftw_plan_many_dft` — function
- `fftw_plan_many_dft_c2r` — function
- `fftw_plan_many_dft_r2c` — function
- `fftw_plan_many_r2r` — function
- `fftw_plan_r2r` — function
- `fftw_plan_r2r_1d` — function
- `fftw_plan_r2r_2d` — function
- `fftw_plan_r2r_3d` — function
- `fftw_plan_with_nthreads` — subroutine
- `fftw_print_plan` — subroutine
- `fftw_set_timelimit` — subroutine
- `fftw_sprint_plan` — function
- `fftw_strlen` — function
- `fftwf_alignment_of` — function
- `fftwf_alloc_complex` — function
- `fftwf_alloc_real` — function
- `fftwf_cleanup` — subroutine
- `fftwf_cleanup_threads` — subroutine
- `fftwf_cost` — function
- `fftwf_destroy_plan` — subroutine
- `fftwf_estimate_cost` — function
- `fftwf_execute_dft` — subroutine
- `fftwf_execute_dft_c2r` — subroutine
- `fftwf_execute_dft_r2c` — subroutine
- `fftwf_execute_r2r` — subroutine
- `fftwf_execute_split_dft` — subroutine
- `fftwf_execute_split_dft_c2r` — subroutine
- `fftwf_execute_split_dft_r2c` — subroutine
- `fftwf_export_wisdom` — subroutine
- `fftwf_export_wisdom_to_file` — subroutine
- `fftwf_export_wisdom_to_filename` — function
- `fftwf_export_wisdom_to_string` — function
- `fftwf_flops` — subroutine
- `fftwf_forget_wisdom` — subroutine
- `fftwf_fprint_plan` — subroutine
- `fftwf_free` — subroutine
- `fftwf_import_system_wisdom` — function
- `fftwf_import_wisdom` — function
- `fftwf_import_wisdom_from_file` — function
- `fftwf_import_wisdom_from_filename` — function
- `fftwf_import_wisdom_from_string` — function
- `fftwf_init_threads` — function
- `fftwf_iodim` — type
- `fftwf_iodim64` — type
- `fftwf_malloc` — function
- `fftwf_plan_dft` — function
- `fftwf_plan_dft_1d` — function
- `fftwf_plan_dft_2d` — function
- `fftwf_plan_dft_3d` — function
- `fftwf_plan_dft_c2r` — function
- `fftwf_plan_dft_c2r_1d` — function
- `fftwf_plan_dft_c2r_2d` — function
- `fftwf_plan_dft_c2r_3d` — function
- `fftwf_plan_dft_r2c` — function
- `fftwf_plan_dft_r2c_1d` — function
- `fftwf_plan_dft_r2c_2d` — function
- `fftwf_plan_dft_r2c_3d` — function
- `fftwf_plan_guru64_dft` — function
- `fftwf_plan_guru64_dft_c2r` — function
- `fftwf_plan_guru64_dft_r2c` — function
- `fftwf_plan_guru64_r2r` — function
- `fftwf_plan_guru64_split_dft` — function
- `fftwf_plan_guru64_split_dft_c2r` — function
- `fftwf_plan_guru64_split_dft_r2c` — function
- `fftwf_plan_guru_dft` — function
- `fftwf_plan_guru_dft_c2r` — function
- `fftwf_plan_guru_dft_r2c` — function
- `fftwf_plan_guru_r2r` — function
- `fftwf_plan_guru_split_dft` — function
- `fftwf_plan_guru_split_dft_c2r` — function
- `fftwf_plan_guru_split_dft_r2c` — function
- `fftwf_plan_many_dft` — function
- `fftwf_plan_many_dft_c2r` — function
- `fftwf_plan_many_dft_r2c` — function
- `fftwf_plan_many_r2r` — function
- `fftwf_plan_r2r` — function
- `fftwf_plan_r2r_1d` — function
- `fftwf_plan_r2r_2d` — function
- `fftwf_plan_r2r_3d` — function
- `fftwf_plan_with_nthreads` — subroutine
- `fftwf_print_plan` — subroutine
- `fftwf_set_timelimit` — subroutine
- `fftwf_sprint_plan` — function

---
## Module: simple_fileio

Files:
- `fileio/simple_fileio.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `add2fbody_char` — function
- `add2fbody_str` — function
- `add2fbody_str_char` — function
- `append2basename_char` — function
- `append2basename_str` — function
- `arr2file_dp` — subroutine
- `arr2file_sp` — subroutine
- `arr2txtfile_int` — subroutine
- `arr2txtfile_real` — subroutine
- `basename` — function
- `del_files_list` — subroutine
- `del_files_numbered` — subroutine
- `fclose` — subroutine
- `file2drarr` — function
- `file2imat` — subroutine
- `file2lmat` — subroutine
- `file2rarr` — function
- `file2rmat` — subroutine
- `fileiochk` — subroutine
- `filelength` — function
- `filepath_char` — function
- `filepath_char_str` — function
- `filepath_str` — function
- `filepath_str_char` — function
- `files_identical` — function
- `fname2ext` — function
- `fname2format` — function
- `fname2iter` — function
- `fname_new_ext_char` — function
- `fname_new_ext_str` — function
- `fopen` — subroutine
- `funit_size` — function
- `get_fbody_char` — function
- `get_fbody_str` — function
- `get_fpath` — function
- `get_open_funits` — subroutine
- `get_relative_path` — function
- `imat2file` — subroutine
- `is_funit_open` — function
- `lmat2file` — subroutine
- `make_dirnames` — function
- `make_filenames` — function
- `move_files2dir` — subroutine
- `move_files_in_cwd` — subroutine
- `nlines` — function
- `read_exit_code` — subroutine
- `read_filetable` — subroutine
- `rm_from_fbody` — function
- `rmat2file` — subroutine
- `simple_copy_file` — subroutine
- `stemname` — function
- `swap_suffix_char` — function
- `swap_suffix_str` — function
- `wait_for_closure` — subroutine
- `write_filetable` — subroutine
- `write_singlelineoftext` — subroutine

---
## Module: simple_fileio_tester

Files:
- `fileio/simple_fileio_tester.f90`

Uses:
- `simple_defs`
- `simple_fileio`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_fileio_tests` — subroutine

Private symbols:
- `test_add2fbody_rm_swap_getfbody` — subroutine
- `test_arr2file_sp_dp_roundtrip` — subroutine
- `test_arr2txtfile_roundtrip` — subroutine
- `test_filepath_overloads` — subroutine
- `test_files_identical` — subroutine
- `test_fname2ext_and_basename` — subroutine
- `test_fname2format` — subroutine
- `test_fname2iter` — subroutine
- `test_fopen_and_fclose_basic` — subroutine
- `test_get_fpath_and_stemname` — subroutine
- `test_get_relative_path` — subroutine
- `test_make_dir_and_file_names` — subroutine
- `test_nlines_and_filelength` — subroutine
- `test_read_write_filetable` — subroutine
- `test_rmat_lmat_roundtrip` — subroutine
- `test_simple_copy_file` — subroutine
- `test_simple_list_files_and_regexp` — subroutine
- `test_simple_mkdir_dir_exists_chdir_getcwd_rmdir` — subroutine
- `test_simple_touch_rename_abspath` — subroutine
- `test_write_singlelineoftext_and_read_exit_code` — subroutine

Local symbols (enclosing procedure scope):
- `test_files_identical::write_bytes` — subroutine

---
## Module: simple_final_rec

Files:
- `main/simple_final_rec.f90`

Uses:
- `iso_fortran_env`
- `simple_commanders_api`
- `simple_parameters`
- `simple_reconstructor_pcg`
- `simple_refine3d_fnames`
- `simple_sigma2_bootstrap`
- `simple_sigma2_files`
- `simple_solve3d_utils`
- `simple_state_weight_set`
- `simple_vol_pproc_policy`

Public symbols:
- `calc_final_rec` — subroutine
- `calc_final_rec_addon` — subroutine

Private symbols:
- `addon_membership_switch` — subroutine
- `final_rec_bootstrap_sigma_iter` — function
- `final_stage_uses_ml_reg` — function
- `prep_final_rec_cline` — subroutine
- `register_state_map` — subroutine

Local symbols (enclosing procedure scope):
- `calc_final_rec::canonical_final_rec_needs_bootstrap` — function
- `calc_final_rec::final_rec_box_changed` — function
- `calc_final_rec_addon::set_addon_rec_cline` — subroutine

---
## Module: simple_flex_cls_expansion

Files:
- `main/flex/simple_flex_cls_expansion.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_flex_pca_posterior`
- `simple_rnd`

Public symbols:
- `flex_cls_allocate_subclasses` — subroutine
- `flex_cls_ctf_mean` — subroutine
- `flex_cls_extend_basis` — subroutine
- `flex_cls_fit` — subroutine
- `flex_cls_frc_cutoff` — function
- `flex_cls_half_frc` — subroutine
- `flex_cls_model` — type
- `flex_cls_pose_tangents` — subroutine
- `flex_cls_shell_frc` — subroutine
- `flex_cls_shell_noise` — subroutine
- `flex_cls_split_floor` — function
- `flex_cls_weighted_mean` — subroutine

Private symbols:
- `denuisanced` — function
- `estep` — subroutine
- `init_basis` — subroutine
- `member_residual` — function
- `model_kill` — subroutine
- `model_new` — subroutine
- `mstep` — subroutine
- `orthonormalise` — subroutine
- `spd_inverse` — subroutine
- `spd_solve_inverse` — subroutine
- `weighted_average` — subroutine
- `weighted_residual` — function

Local symbols (enclosing procedure scope):
- `flex_cls_extend_basis::frc_weights` — subroutine
- `flex_cls_extend_basis::solve3` — subroutine
- `flex_cls_extend_basis::unpack_lhs` — subroutine
- `flex_cls_pose_tangents::value_at` — function

---
## Module: simple_flex_cls_expansion_tester

Files:
- `main/flex/simple_flex_cls_expansion_tester.f90`

Uses:
- `simple_flex_cls_expansion`
- `simple_test_utils`

Public symbols:
- `run_all_flex_cls_expansion_tests` — subroutine

Private symbols:
- `test_allocation_capacity` — subroutine
- `test_allocation_equal_views` — subroutine
- `test_allocation_group_labels` — subroutine
- `test_allocation_one_per_class` — subroutine
- `test_allocation_views` — subroutine
- `test_split_floor` — subroutine

---
## Module: simple_flex_cls_kernel

Files:
- `main/flex/simple_flex_cls_kernel.f90`
- `main/flex/simple_flex_cls_kernel_tester.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_flex_cls_expansion`
- `simple_linalg`
- `simple_rnd`
- `simple_stat`
- `simple_test_utils`

Public symbols:
- `flex_cls_kernel` — type
- `run_all_flex_cls_kernel_tests` — subroutine
- `run_all_flex_cls_kernel_tests` — subroutine

Private symbols:
- `average` — subroutine
- `cgauss` — function
- `coeff_sets` — subroutine
- `cross_validate` — subroutine
- `cv_errors` — subroutine
- `get_member` — function
- `get_nrep` — function
- `get_zscore` — function
- `hand_model` — subroutine
- `kernel_mean` — subroutine
- `kill` — subroutine
- `labels` — function
- `lattice` — subroutine
- `logweights` — subroutine
- `make_class` — subroutine
- `mean_cov` — function
- `ml_latents` — subroutine
- `neff` — function
- `new` — subroutine
- `precisions` — subroutine
- `preimage` — subroutine
- `project` — subroutine
- `report` — subroutine
- `select` — subroutine
- `signal_space` — subroutine
- `split_class` — subroutine
- `test_coeff_sets` — subroutine
- `test_kernel_budget` — subroutine
- `test_kernel_null` — subroutine
- `test_kernel_two_states` — subroutine
- `test_ml_latents_exact` — subroutine
- `test_representatives` — subroutine
- `test_stein_correction` — subroutine
- `two_states` — subroutine
- `whole` — subroutine

Local symbols (enclosing procedure scope):
- `test_stein_correction::latent_dev` — function

---
## Module: simple_flex_pca_application

Files:
- `main/flex/run/simple_flex_pca_application.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_flex_pca_basis`
- `simple_flex_pca_delivery_3d`
- `simple_flex_pca_embed`
- `simple_flex_pca_embedding_io`
- `simple_flex_pca_fit_driver`
- `simple_flex_pca_merge`
- `simple_flex_pca_mstep`
- `simple_flex_pca_pcg`
- `simple_flex_pca_planes`
- `simple_flex_pca_project_gateway`
- `simple_flex_pca_records`
- `simple_flex_pca_rounds`
- `simple_flex_pca_run_types`
- `simple_flex_pca_stages`
- `simple_flex_pca_state_service`
- `simple_flex_pca_targets`
- `simple_flex_pca_weights`
- `simple_image`
- `simple_parameters`
- `simple_reconstructor`

Public symbols:
- `flex_pca_application` — type

Private symbols:
- `app_kill` — subroutine
- `app_run` — subroutine
- `embed_all` — subroutine
- `execute_worker_stage` — subroutine
- `infer_states` — subroutine
- `obtain_embedding` — subroutine
- `prepare` — subroutine
- `publish` — subroutine
- `reconstruct_states` — subroutine
- `resume_from_cache` — subroutine

---
## Module: simple_flex_pca_application_tester

Files:
- `main/flex/run/simple_flex_pca_application_tester.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_cmdline`
- `simple_commanders_flex_pca`
- `simple_commanders_sim`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_image`
- `simple_imghead`
- `simple_oris`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_string`
- `simple_syslib`
- `simple_test_truth_metrics`
- `simple_test_utils`
- `simple_timer`
- `simple_type_defs`
- `simple_ui`

Public symbols:
- `build_flex_pca_phantom_project` — subroutine
- `configure_flex_pca_phantom` — subroutine
- `create_flex_pca_phantom_fixture` — subroutine
- `read_state_weight_table` — subroutine
- `run_all_flex_pca_application_tests` — subroutine
- `write_phantom_truths` — subroutine

Private symbols:
- `ab_dump` — subroutine
- `arm_result` — type
- `interleave_stacks` — function
- `leading_pc_correlation` — function
- `phantom_layout` — subroutine
- `run_application_arm` — subroutine
- `simulate_stack` — function
- `state_truth_correlations` — subroutine
- `test_two_state_backend_matrix` — subroutine
- `validate_arm` — subroutine
- `write_uniform_orientations` — subroutine

---
## Module: simple_flex_pca_artifacts

Files:
- `main/flex/run/simple_flex_pca_artifacts.f90`

Uses:
- `simple_core_module_api`
- `simple_parameters`

Public symbols:
- `flex_pca_artifact_catalog` — type
- `flex_pca_local_part_dir` — function

Private symbols:
- `artifact_catalog_kill` — subroutine
- `artifact_catalog_new` — subroutine
- `artifact_part_fname` — function
- `artifact_part_path` — function

---
## Module: simple_flex_pca_basis

Files:
- `main/flex/fit/simple_flex_pca_basis.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_defs_flex`
- `simple_flex_pca_fit_types`
- `simple_flex_pca_mstep`
- `simple_flex_pca_pcg`
- `simple_flex_pca_planes`
- `simple_flex_pca_records`
- `simple_flex_pca_rounds`
- `simple_flex_pca_util`
- `simple_flex_reconstructor_latent_ops`
- `simple_gridding`
- `simple_image`
- `simple_linalg`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_math`
- `simple_ori`
- `simple_parameters`
- `simple_reconstructor`

Public symbols:
- `align_basis_to_reference` — subroutine
- `apply_cached_mean_scale` — subroutine
- `basis_recs_from_images` — subroutine
- `cov_herm_inner` — function
- `covariance_kfromto` — function
- `cross_half_subspace_angles` — subroutine
- `deflate_against_basis` — subroutine
- `estimate_covariance_mean` — subroutine
- `estimate_mean_scale` — subroutine
- `init_basis_datafree` — subroutine
- `init_mean_reconstructor` — subroutine
- `load_probe_basis` — subroutine
- `load_probe_state` — subroutine
- `orthonormalize_representatives` — subroutine
- `save_probe_state` — subroutine

Private symbols:
- `basis_to_real_representatives` — subroutine
- `cov_half_parity` — function
- `cov_herm_selfpower` — subroutine
- `em_calibrate_noise_prior` — subroutine
- `estimate_mean_from_data` — subroutine
- `form_reconstruction_plane` — subroutine
- `mean_scale_cache_fname` — function
- `particle_contrast` — function
- `pick_next_lowfreq` — subroutine
- `plane_hf_power` — subroutine
- `plane_shell_cross_accum` — subroutine
- `realize_hermitian_volume` — subroutine
- `select_frequencies_lowfreq` — subroutine
- `write_mean_scale` — subroutine

Local symbols (enclosing procedure scope):
- `cross_half_subspace_angles::inv_sqrt_sym` — subroutine

---
## Module: simple_flex_pca_crossfsc

Files:
- `main/flex/fit/simple_flex_pca_crossfsc.f90`

Uses:
- `simple_core_module_api`
- `simple_flex_reconstructor_latent_ops`

Public symbols:
- `crossfsc_file` — type
- `crossfsc_harvest_h` — subroutine
- `crossfsc_inband_mean` — function
- `crossfsc_record` — type
- `crossfsc_to_invtau2` — subroutine

Private symbols:
- `crossfsc_append` — subroutine
- `crossfsc_kill` — subroutine
- `crossfsc_kill_record` — subroutine
- `crossfsc_latest_upto` — function
- `crossfsc_load` — subroutine
- `crossfsc_write` — subroutine
- `read_record` — subroutine
- `write_record` — subroutine

---
## Module: simple_flex_pca_deconv

Files:
- `main/flex/states/simple_flex_pca_deconv.f90`

Uses:
- `simple_core_module_api`
- `simple_linalg`
- `simple_xd_gmm`

Public symbols:
- `calibrate_noise_scale` — subroutine
- `deconvolve_latent` — subroutine

Private symbols:
- `noise_and_projection` — subroutine
- `spd_inverse_or_diag` — subroutine

---
## Module: simple_flex_pca_delivery_3d

Files:
- `main/flex/run/simple_flex_pca_delivery_3d.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_flex_pca_plot`
- `simple_flex_pca_records`
- `simple_parameters`
- `simple_umap`

Public symbols:
- `deliver_latent_readouts` — subroutine
- `write_covariance_eigenvolumes` — subroutine
- `write_covariance_manifest` — subroutine
- `write_covariance_tables` — subroutine
- `write_states_umap_figure` — subroutine

---
## Module: simple_flex_pca_embed

Files:
- `main/flex/fit/simple_flex_pca_embed.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_flex_pca_artifacts`
- `simple_flex_pca_basis`
- `simple_flex_pca_fit_types`
- `simple_flex_pca_planes`
- `simple_flex_pca_posterior`
- `simple_flex_pca_records`
- `simple_flex_pca_rounds`
- `simple_flex_pca_stages`
- `simple_flex_reconstructor_latent_ops`
- `simple_linalg`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_ori`
- `simple_parameters`
- `simple_srch_sort_loc`
- `simple_stat`

Public symbols:
- `embed_latents_with_contrast` — subroutine

Private symbols:
- `read_embed_stats_part` — subroutine
- `reduce_embed_zhalf_parts` — subroutine
- `write_embed_stats_part` — subroutine

---
## Module: simple_flex_pca_embedding_io

Files:
- `main/flex/run/simple_flex_pca_embedding_io.f90`

Uses:
- `simple_core_module_api`
- `simple_flex_pca_records`
- `simple_syslib`

Public symbols:
- `read_deconv_block` — subroutine
- `read_embedding_cache` — subroutine
- `read_noise_scale` — subroutine
- `write_deconv_block` — subroutine
- `write_embedding_cache` — subroutine
- `write_noise_scale` — subroutine

Private symbols:
- `open_trailing_block` — subroutine
- `raw_block_end` — function
- `write_trailing_block` — subroutine

---
## Module: simple_flex_pca_fit_driver

Files:
- `main/flex/fit/simple_flex_pca_fit_driver.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_defs_flex`
- `simple_flex_pca_basis`
- `simple_flex_pca_embed`
- `simple_flex_pca_pairmerge`
- `simple_flex_pca_pcg`
- `simple_flex_pca_planes`
- `simple_flex_pca_records`
- `simple_flex_pca_rounds`
- `simple_flex_pca_stages`
- `simple_flex_pca_util`
- `simple_flex_probe_fit`
- `simple_image`
- `simple_parameters`
- `simple_reconstructor`

Public symbols:
- `embed_worker_pass` — subroutine
- `probe_worker_pass` — subroutine
- `run_flex_pca_paired` — subroutine

Private symbols:
- `prepare_fit_environments` — subroutine
- `run_flex_pca_paired_worker` — subroutine
- `write_paired_eigen_table` — subroutine

---
## Module: simple_flex_pca_fit_types

Files:
- `main/flex/fit/simple_flex_pca_fit_types.f90`

Uses:
- `simple_core_module_api`
- `simple_flex_pca_crossfsc`
- `simple_flex_pca_mstep`
- `simple_flex_pca_polar`
- `simple_flex_pca_records`
- `simple_image`
- `simple_reconstructor`

Public symbols:
- `cleanup_plane` — subroutine
- `flex_fit` — type
- `flex_fit_diag` — type
- `flex_fit_history` — type
- `flex_fit_iter` — type
- `flex_fit_spec` — type
- `flex_probe_part` — type
- `probe_part_borrow` — subroutine
- `probe_part_restore` — subroutine
- `xfsc_ctx_t` — type

Private symbols:
- `flex_fit_diag_kill` — subroutine
- `flex_fit_history_kill` — subroutine
- `flex_fit_iter_kill` — subroutine
- `flex_fit_kill` — subroutine
- `flex_fit_new` — subroutine
- `flex_fit_spec_kill` — subroutine
- `probe_part_kill` — subroutine
- `xfsc_ctx_kill` — subroutine
- `xfsc_ctx_new` — subroutine

---
## Module: simple_flex_pca_gmm

Files:
- `main/flex/states/simple_flex_pca_gmm.f90`

Uses:
- `simple_clustering_utils`
- `simple_core_module_api`
- `simple_flex_pca_targets`
- `simple_flex_pca_util`
- `simple_gmm`
- `simple_stat`

Public symbols:
- `gmm_auto_state_weights` — subroutine
- `gmm_state_weights` — subroutine

Private symbols:
- `continuum_region_weights` — subroutine

---
## Module: simple_flex_pca_merge

Files:
- `main/flex/states/simple_flex_pca_merge.f90`

Uses:
- `simple_core_module_api`
- `simple_defs_flex`
- `simple_flex_pca_pcg`
- `simple_flex_pca_records`
- `simple_hac`
- `simple_image`
- `simple_parameters`
- `simple_stat`

Public symbols:
- `two_gate_state_merge` — subroutine

Private symbols:
- `pair_map_ratio` — subroutine
- `view_coverage_chi2` — subroutine

---
## Module: simple_flex_pca_mstep

Files:
- `main/flex/fit/simple_flex_pca_mstep.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_flex_pca_pcg`
- `simple_flex_reconstructor_latent_ops`
- `simple_image`
- `simple_ori`
- `simple_parameters`
- `simple_reconstructor`

Public symbols:
- `flex_fit_mstep` — type
- `init_basis_reconstructor` — subroutine

Private symbols:
- `mstep_accumulate_batch` — subroutine
- `mstep_apply_ridge` — subroutine
- `mstep_begin_iteration` — subroutine
- `mstep_kill` — subroutine
- `mstep_kill_iteration` — subroutine
- `mstep_reduce_local` — subroutine
- `mstep_snapshot_for_pair_merge` — subroutine
- `mstep_solve_halves` — subroutine

Local symbols (enclosing procedure scope):
- `mstep_solve_halves::log_pcg_outcome` — subroutine

---
## Module: simple_flex_pca_pairmerge

Files:
- `main/flex/fit/simple_flex_pca_pairmerge.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_flex_pca_artifacts`
- `simple_flex_pca_basis`
- `simple_flex_pca_crossfsc`
- `simple_flex_pca_mstep`
- `simple_flex_pca_pcg`
- `simple_flex_pca_records`
- `simple_flex_pca_util`
- `simple_flex_probe_fit`
- `simple_flex_reconstructor_latent_ops`
- `simple_gridding`
- `simple_image`
- `simple_linalg`
- `simple_parameters`
- `simple_reconstructor`

Public symbols:
- `probe_paired_merge` — subroutine

Private symbols:
- `gauge_fix_to_frame` — subroutine
- `init_deflated_matchcos` — subroutine
- `matched_abs_cos` — subroutine
- `merge_deflate_mean_shaped` — subroutine
- `merge_pair_algebra_selfcheck` — subroutine
- `merged_delivery` — subroutine
- `polar_factor` — subroutine
- `rotate_pair_voxel` — subroutine
- `rotate_rho4_packed_add` — subroutine
- `rotate_rho_packed_add` — subroutine
- `rotate_rhs_packed_add` — subroutine

---
## Module: simple_flex_pca_pcg

Files:
- `main/flex/fit/simple_flex_pca_pcg.f90`
- `main/flex/fit/simple_flex_pca_pcg_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_core_module_api`
- `simple_flex_reconstructor_latent_ops`
- `simple_gridding`
- `simple_image`
- `simple_linalg`
- `simple_math`
- `simple_ori_utils`
- `simple_pcg_lattice`
- `simple_pcg_solver`
- `simple_reconstructor`
- `simple_rnd`
- `simple_test_utils`

Public symbols:
- `flex_pcg_environment` — type
- `flex_pcg_t` — type
- `test_flex_pcg_operator` — subroutine
- `test_flex_pcg_operator` — subroutine

Private symbols:
- `accumulate` — subroutine
- `accumulate_rhs` — subroutine
- `alloc_accum` — subroutine
- `alloc_packed` — subroutine
- `alloc_rhs_accum` — subroutine
- `alloc_rhs_packed` — subroutine
- `apply_operator` — subroutine
- `apply_precond` — subroutine
- `bandlimit_img` — subroutine
- `build_lists` — subroutine
- `bytes_accum` — function
- `bytes_packed` — function
- `bytes_rhs_accum` — function
- `bytes_rhs_packed` — function
- `cg_core` — subroutine
- `clear_ridge` — subroutine
- `dot_all` — function
- `ensure_pool` — subroutine
- `finalize` — subroutine
- `finalize_rhs` — subroutine
- `flex_env_active` — function
- `flex_env_copy_from` — subroutine
- `flex_env_kill` — subroutine
- `flex_env_new` — subroutine
- `flex_mskfile_set` — function
- `flex_pcg_adapter` — type
- `flex_pcg_install_window` — subroutine
- `flex_pcg_support_volume` — subroutine
- `flex_set_window_sphere` — subroutine
- `flex_vector_apply` — subroutine
- `flex_vector_dot` — function
- `flex_vector_precond` — subroutine
- `flex_vector_size` — function
- `flex_window_apply` — subroutine
- `flex_window_apply_rec` — subroutine
- `fold_accum` — subroutine
- `fold_accum_dense` — subroutine
- `fold_rhs` — subroutine
- `fold_rhs_dense` — subroutine
- `get_nexp` — function
- `get_npairs` — function
- `get_npk` — function
- `is_ready` — function
- `kill` — subroutine
- `new` — subroutine
- `prep_floor` — subroutine
- `require_lists` — subroutine
- `scatter_pairs_nowrap` — subroutine
- `scatter_pairs_nowrap_dense` — subroutine
- `scatter_pairs_wrap` — subroutine
- `scatter_pairs_wrap_dense` — subroutine
- `scatter_rhs_nowrap` — subroutine
- `scatter_rhs_nowrap_dense` — subroutine
- `scatter_rhs_wrap` — subroutine
- `scatter_rhs_wrap_dense` — subroutine
- `set_band` — subroutine
- `set_lambda_relative` — subroutine
- `set_ridge` — subroutine
- `set_verbose` — subroutine
- `solve` — subroutine
- `win_wraps` — function
- `window_product` — subroutine

Local symbols (enclosing procedure scope):
- `cg_core::audit_terms` — subroutine
- `flex_pcg_support_volume::occupancy` — subroutine
- `solve::put_back` — subroutine
- `test_flex_pcg_operator::adjoint_add` — subroutine
- `test_flex_pcg_operator::compare3` — subroutine
- `test_flex_pcg_operator::deposit_density` — subroutine
- `test_flex_pcg_operator::exps` — subroutine
- `test_flex_pcg_operator::lowpass` — subroutine
- `test_flex_pcg_operator::sample_exact` — subroutine

---
## Module: simple_flex_pca_plane_cache

Files:
- `main/flex/fit/simple_flex_pca_plane_cache.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_image`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_parameters`
- `simple_syslib`

Public symbols:
- `flex_pca_plane_cache` — type
- `plane_cache_contract_header` — function
- `plane_cache_master_matches` — function
- `plane_cache_worker_matches` — function

Private symbols:
- `cache_available` — function
- `cache_delete` — subroutine
- `cache_kill` — subroutine
- `cache_new` — subroutine
- `cache_pristine` — function
- `fold_str` — subroutine
- `header_key` — function
- `open_for_read` — subroutine
- `path_hash` — function
- `plane_cache_adopt` — subroutine
- `plane_cache_applicable` — function
- `plane_cache_dir` — function
- `plane_cache_ensure` — subroutine
- `plane_cache_fill` — subroutine
- `plane_cache_fname` — function
- `plane_cache_probe` — function
- `plane_cache_read_batch` — subroutine
- `plane_cache_run_token` — function
- `plane_cache_space_ok` — function
- `project_provenance` — subroutine
- `selection_hash` — function
- `set_geometry` — subroutine

---
## Module: simple_flex_pca_planes

Files:
- `main/flex/fit/simple_flex_pca_planes.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_flex_pca_plane_cache`
- `simple_image`
- `simple_parameters`

Public symbols:
- `flex_plane_store` — type

Private symbols:
- `all_held` — function
- `mem_available_gb` — function
- `plane_gb` — function
- `plane_store_adopt_cache` — subroutine
- `plane_store_cache_in_use` — function
- `plane_store_enable` — subroutine
- `plane_store_ensure_cache` — subroutine
- `plane_store_fetch` — subroutine
- `plane_store_fill_cached_image` — subroutine
- `plane_store_kill` — subroutine
- `plane_store_read_cache_batch` — subroutine
- `plane_store_release` — subroutine
- `plane_store_store` — subroutine

---
## Module: simple_flex_pca_plot

Files:
- `main/flex/simple_flex_pca_plot.f90`

Uses:
- `simple_core_module_api`
- `simple_jpg`
- `simple_math`
- `simple_srch_sort_loc`

Public symbols:
- `flex_plot_latent_jpg` — subroutine

Private symbols:
- `density_panel` — subroutine
- `fill` — subroutine
- `frame` — subroutine
- `inside` — function
- `limits` — subroutine
- `scatter_panel` — subroutine

---
## Module: simple_flex_pca_polar

Files:
- `main/flex/fit/simple_flex_pca_polar.f90`

Uses:
- `simple_core_module_api`
- `simple_flex_reconstructor_latent_ops`
- `simple_reconstructor`

Public symbols:
- `flex_polar_bank` — type
- `polar_assign_directions` — subroutine
- `polar_grid_build` — subroutine
- `polar_grid_kill` — subroutine
- `polar_grid_t` — type
- `polar_project_recs` — subroutine
- `polar_relative_inplane` — subroutine
- `polar_sample_particle_fused` — subroutine

Private symbols:
- `flex_polar_bank_kill` — subroutine
- `polar_nang` — function

---
## Module: simple_flex_pca_posterior

Files:
- `main/flex/fit/simple_flex_pca_posterior.f90`

Uses:
- `simple_core_module_api`
- `simple_flex_pca_fit_types`
- `simple_kmeans`
- `simple_linalg`

Public symbols:
- `map_sampling_precision` — subroutine
- `mcfa_condition` — subroutine
- `mcfa_init` — subroutine
- `mcfa_mstep` — subroutine
- `probe_solve_mix` — subroutine
- `probe_solve_plain` — subroutine
- `quad_form` — function
- `spd_inv_dp` — subroutine
- `spd_solve_dp` — subroutine

Private symbols:
- `mean_abs_diag` — function

---
## Module: simple_flex_pca_project_gateway

Files:
- `main/flex/run/simple_flex_pca_project_gateway.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_cmdline`
- `simple_commanders_rec`
- `simple_core_module_api`
- `simple_eo_partition`
- `simple_flex_pca_rounds`
- `simple_parameters`
- `simple_refine3d_fnames`
- `simple_sigma2_bootstrap`
- `simple_sigma2_files`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_syslib`

Public symbols:
- `deliver_state_maps` — subroutine
- `ensure_canonical_sigma_state` — subroutine
- `load_and_validate_sigma` — subroutine
- `register_embedding_artifact` — subroutine
- `validate_covariance_inputs` — subroutine
- `write_discrete_state_project` — subroutine
- `write_state_weight_set` — subroutine

Private symbols:
- `check_sigma_state` — subroutine
- `read_pind_list` — subroutine
- `write_sigma_state` — subroutine

---
## Module: simple_flex_pca_records

Files:
- `main/flex/run/simple_flex_pca_records.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_reconstructor`

Public symbols:
- `flex_fit_model` — type
- `flex_latent` — type
- `flex_latent_readout` — type
- `flex_selection` — type
- `flex_state_set` — type

Private symbols:
- `flex_fit_model_kill` — subroutine
- `latent_kill` — subroutine
- `readout_kill` — subroutine
- `selection_kill` — subroutine
- `state_set_kill` — subroutine

---
## Module: simple_flex_pca_rounds

Files:
- `main/flex/run/simple_flex_pca_rounds.f90`

Uses:
- `simple_core_module_api`
- `simple_flex_pca_artifacts`
- `simple_flex_pca_stages`
- `simple_parameters`

Public symbols:
- `flex_pca_rounds` — type
- `flex_pca_rounds_shmem` — type

Private symbols:
- `plan_iface` — subroutine
- `rounds_distributed` — function
- `rounds_is_master` — function
- `rounds_is_worker` — function
- `rounds_kill_artifacts` — subroutine
- `rounds_nparts` — function
- `rounds_part_fname` — function
- `rounds_part_path` — function
- `rounds_set_part_dir` — subroutine
- `run_iface` — subroutine
- `shmem_plan_partitions` — subroutine
- `shmem_run_stage` — subroutine

---
## Module: simple_flex_pca_run_types

Files:
- `main/flex/run/simple_flex_pca_run_types.f90`

Uses:
- `simple_core_module_api`
- `simple_flex_pca_records`
- `simple_reconstructor`

Public symbols:
- `flex_run_session` — type

Private symbols:
- `session_clamp_state_axis` — subroutine
- `session_kill` — subroutine

---
## Module: simple_flex_pca_stages

Files:
- `main/flex/run/simple_flex_pca_stages.f90`

Public symbols:
- `flex_pca_half_of` — function
- `flex_stage_request` — type

---
## Module: simple_flex_pca_state_service

Files:
- `main/flex/run/simple_flex_pca_state_service.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_flex_pca_deconv`
- `simple_flex_pca_embedding_io`
- `simple_flex_pca_records`
- `simple_flex_pca_rounds`
- `simple_flex_pca_weights`
- `simple_image`
- `simple_parameters`
- `simple_rec3d_service`
- `simple_rnd`
- `simple_srch_sort_loc`

Public symbols:
- `apply_latent_deconvolution` — subroutine
- `auto_box_crop` — function
- `auto_min_neff` — function
- `auto_state_count` — function
- `cv_select_bandwidths` — subroutine
- `delete_state_halves` — subroutine
- `infile_path` — function
- `place_states_with_population_floor` — subroutine
- `prune_underpopulated_states` — subroutine
- `reconstruct_state_halves` — subroutine
- `state_half_fname` — function

---
## Module: simple_flex_pca_strategy

Files:
- `main/strategies/parallelization/simple_flex_pca_strategy.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_exec_helpers`
- `simple_flex_pca_application`
- `simple_flex_pca_artifacts`
- `simple_flex_pca_project_gateway`
- `simple_flex_pca_rounds`
- `simple_flex_pca_stages`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`

Public symbols:
- `create_flex_pca_strategy` — function
- `flex_pca_strategy` — type
- `flex_pca_worker_strategy` — type

Private symbols:
- `cleanup_interface` — subroutine
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `flex_pca_master_rounds` — type
- `flex_pca_master_strategy` — type
- `flex_pca_shmem_strategy` — type
- `init_interface` — subroutine
- `master_cleanup` — subroutine
- `master_execute` — subroutine
- `master_finalize_run` — subroutine
- `master_initialize` — subroutine
- `master_plan_partitions` — subroutine
- `master_run_stage` — subroutine
- `shmem_cleanup` — subroutine
- `shmem_execute` — subroutine
- `shmem_finalize_run` — subroutine
- `shmem_initialize` — subroutine
- `worker_cleanup` — subroutine
- `worker_execute` — subroutine
- `worker_finalize_run` — subroutine
- `worker_initialize` — subroutine

---
## Module: simple_flex_pca_targets

Files:
- `main/flex/states/simple_flex_pca_targets.f90`

Uses:
- `simple_core_module_api`
- `simple_kcenter`
- `simple_kd_tree`
- `simple_kmeans`
- `simple_linalg`
- `simple_srch_sort_loc`

Public symbols:
- `component_reliability_proxy` — subroutine
- `diffusion_kcenter_targets` — subroutine
- `kmeans_latent_targets` — subroutine
- `orth_block` — subroutine
- `path_latent_targets` — subroutine
- `reliability_path_targets` — subroutine
- `sort_block_desc` — subroutine

---
## Module: simple_flex_pca_tester

Files:
- `main/flex/simple_flex_pca_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_defs_flex`
- `simple_flex_pca_basis`
- `simple_flex_pca_crossfsc`
- `simple_flex_pca_deconv`
- `simple_flex_pca_embedding_io`
- `simple_flex_pca_fit_types`
- `simple_flex_pca_pcg`
- `simple_flex_pca_plane_cache`
- `simple_flex_pca_polar`
- `simple_flex_pca_posterior`
- `simple_flex_pca_records`
- `simple_flex_pca_rounds`
- `simple_flex_pca_state_service`
- `simple_flex_pca_util`
- `simple_flex_pca_weights`
- `simple_flex_probe_fit`
- `simple_flex_reconstructor_latent_ops`
- `simple_gridding`
- `simple_image`
- `simple_kbinterpol`
- `simple_linalg`
- `simple_ori`
- `simple_parameters`
- `simple_reconstructor`
- `simple_rnd`
- `simple_sp_project`
- `simple_sym`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_flex_pca_lib_tests` — subroutine
- `run_all_flex_pca_tests` — subroutine

Private symbols:
- `diagonal` — function
- `invert_gauss_jordan` — subroutine
- `test_auto_settings` — subroutine
- `test_crossfit_fsc` — subroutine
- `test_deconvolution` — subroutine
- `test_embedding_cache_io` — subroutine
- `test_estep_former` — subroutine
- `test_flex_crossfsc` — subroutine
- `test_flex_estep_former` — subroutine
- `test_kernel_bandwidth` — subroutine
- `test_mstep_toy` — subroutine
- `test_plane_cache_contract` — subroutine
- `test_population_floor` — subroutine
- `test_posterior_solve` — subroutine
- `test_probe_codec` — subroutine
- `test_probe_part_codec` — subroutine
- `test_prune_states` — subroutine
- `test_state_weights` — subroutine

Local symbols (enclosing procedure scope):
- `test_deconvolution::gauss` — function
- `test_flex_crossfsc::init_context` — subroutine
- `test_flex_crossfsc::install_bases` — subroutine
- `test_flex_estep_former::blob` — function
- `test_flex_estep_former::build_exact_positions` — subroutine
- `test_flex_estep_former::count_halfplane` — function
- `test_flex_estep_former::direct_cartesian_stats` — subroutine
- `test_flex_estep_former::direct_polar_stats` — subroutine
- `test_flex_estep_former::linearity_relerr` — function
- `test_flex_estep_former::oracle_gather` — function
- `test_flex_estep_former::oracle_weights` — subroutine
- `test_flex_estep_former::prepare_reconstructor` — subroutine
- `test_flex_estep_former::stats_relerr` — function
- `test_mstep_toy::assemble_dense_system` — subroutine
- `test_mstep_toy::dense_cholesky` — subroutine
- `test_mstep_toy::finalize_reconstructors` — subroutine
- `test_mstep_toy::init_reference_plane` — subroutine
- `test_mstep_toy::new_empty_reconstructor` — subroutine

---
## Module: simple_flex_pca_util

Files:
- `main/flex/simple_flex_pca_util.f90`

Uses:
- `simple_core_module_api`
- `simple_defs_flex`
- `simple_image`
- `simple_oris`
- `simple_stat`

Public symbols:
- `chi2_median` — function
- `cov_accum_bytes` — function
- `cov_dim_budget` — function
- `cov_signal_rank` — function
- `cov_stage_subsample` — subroutine
- `dilation_template` — subroutine
- `kernel_weights_at_bandwidth` — subroutine
- `project_onto_target_polyline` — subroutine
- `two_gauss_unimodal` — function

---
## Module: simple_flex_pca_weights

Files:
- `main/flex/states/simple_flex_pca_weights.f90`

Uses:
- `simple_core_module_api`
- `simple_defs_flex`
- `simple_flex_pca_gmm`
- `simple_flex_pca_records`
- `simple_flex_pca_targets`
- `simple_flex_pca_util`
- `simple_linalg`
- `simple_srch_sort_loc`
- `simple_stat`

Public symbols:
- `bandwidth_cv_adopt` — subroutine
- `bandwidth_cv_error` — function
- `bandwidth_cv_grid` — subroutine
- `build_covariance_state_weights` — subroutine

---
## Module: simple_flex_pcg_tester

Files:
- `main/flex/simple_flex_pcg_tester.f90`

Uses:
- `simple_flex_pca_pcg`
- `simple_test_utils`

Public symbols:
- `run_all_flex_pcg_lib_tests` — subroutine
- `run_all_flex_pcg_sweep_tests` — subroutine
- `run_all_flex_pcg_tests` — subroutine

Private symbols:
- `check_operator` — subroutine

---
## Module: simple_flex_probe_fit

Files:
- `main/flex/fit/simple_flex_probe_fit.f90`
- `main/flex/fit/simple_flex_probe_fit_crossfsc.f90`
- `main/flex/fit/simple_flex_probe_fit_engine.f90`
- `main/flex/fit/simple_flex_probe_fit_estep.f90`
- `main/flex/fit/simple_flex_probe_fit_update.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_core_module_api`
- `simple_defs_flex`
- `simple_flex_pca_artifacts`
- `simple_flex_pca_basis`
- `simple_flex_pca_crossfsc`
- `simple_flex_pca_fit_types`
- `simple_flex_pca_mstep`
- `simple_flex_pca_pcg`
- `simple_flex_pca_planes`
- `simple_flex_pca_polar`
- `simple_flex_pca_posterior`
- `simple_flex_pca_records`
- `simple_flex_pca_rounds`
- `simple_flex_pca_stages`
- `simple_flex_pca_util`
- `simple_flex_reconstructor_latent_ops`
- `simple_gridding`
- `simple_image`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_math`
- `simple_ori`
- `simple_parameters`
- `simple_reconstructor`

Public symbols:
- `fit_engine_iterate` — subroutine
- `fit_engine_iterate` — subroutine
- `fit_estep_former_polar` — subroutine
- `fit_estep_former_polar` — subroutine
- `flex_mean_ref` — type
- `flex_probe_fit` — type
- `polar_ring_gram` — subroutine
- `polar_ring_gram` — subroutine
- `polar_ring_selfpower` — function
- `polar_ring_selfpower` — function
- `probe_subspace_iteration` — subroutine
- `probe_subspace_iteration` — subroutine
- `reduce_probe_parts` — subroutine
- `reduce_probe_parts` — subroutine
- `write_probe_part` — subroutine
- `write_probe_part` — subroutine
- `xfsc_paired_record` — subroutine
- `xfsc_paired_record` — subroutine

Private symbols:
- `close_probe_part_read` — subroutine
- `close_probe_part_read` — subroutine
- `close_probe_part_write` — subroutine
- `close_probe_part_write` — subroutine
- `cov_polar_enabled` — function
- `cov_polar_enabled` — function
- `cov_polar_ndir` — function
- `cov_polar_ndir` — function
- `dispatch_round` — subroutine
- `fit_batch_insert` — subroutine
- `fit_batch_insert` — subroutine
- `fit_estep_bank_prepare` — subroutine
- `fit_estep_bank_prepare` — subroutine
- `fit_estep_batch_begin` — subroutine
- `fit_estep_batch_begin` — subroutine
- `fit_estep_begin_stage` — subroutine
- `fit_estep_begin_stage` — subroutine
- `fit_estep_particle` — subroutine
- `fit_estep_particle` — subroutine
- `fit_estep_pass` — subroutine
- `fit_estep_pass` — subroutine
- `fit_estep_solve_stats` — subroutine
- `fit_estep_solve_stats` — subroutine
- `fit_iter_begin` — subroutine
- `fit_iter_begin` — subroutine
- `fit_iter_finish` — subroutine
- `fit_iter_finish` — subroutine
- `fit_iter_reduce` — subroutine
- `fit_iter_reduce` — subroutine
- `fit_polar_bank_build` — subroutine
- `fit_polar_bank_build` — subroutine
- `fold_probe_part_fit` — subroutine
- `fold_probe_part_fit` — subroutine
- `fold_probe_part_kernels` — subroutine
- `open_probe_part_read` — subroutine
- `open_probe_part_read` — subroutine
- `open_probe_part_write` — subroutine
- `open_probe_part_write` — subroutine
- `paired_reduce_parts` — subroutine
- `paired_reduce_parts` — subroutine
- `pk_mask_c1` — subroutine
- `pk_mask_r1` — subroutine
- `pk_mask_r2` — subroutine
- `pk_mask_to_idx` — subroutine
- `polar_hybrid_exact_accum` — subroutine
- `polar_hybrid_exact_accum` — subroutine
- `probe_fit_merge_stash` — subroutine
- `probe_fit_merge_stash` — subroutine
- `probe_part_mixture_from_threads` — subroutine
- `reduce_single_fit_parts` — subroutine
- `subtract_mean_banded` — subroutine
- `subtract_mean_banded` — subroutine
- `write_probe_part_fit` — subroutine
- `write_probe_part_fit` — subroutine
- `write_probe_part_kernels` — subroutine
- `write_worker_parts` — subroutine
- `xfsc_build_invtau2` — subroutine
- `xfsc_prep_iter` — subroutine
- `xfsc_prep_iter` — subroutine
- `xfsc_setup` — subroutine
- `xfsc_setup` — subroutine
- `xfsc_teardown` — subroutine
- `xfsc_teardown` — subroutine

---
## Module: simple_flex_reconstructor_latent_ops

Files:
- `main/flex/fit/simple_flex_reconstructor_latent_ops.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_flex_pca_planes`
- `simple_fplane_grid`
- `simple_image`
- `simple_linalg`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_reconstructor`

Public symbols:
- `add_invtausq2rho_coupled` — subroutine
- `pair_index` — function
- `planes_batch_load` — subroutine
- `prep_imgs4projected_model` — subroutine
- `project_fplanes_mean_basis` — subroutine
- `projected_model_kfromto` — function
- `solve_coupled_basis_exp` — subroutine

---
## Module: simple_flex_traj_denoise

Files:
- `main/flex/simple_flex_traj_denoise.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_flex_cls_expansion`
- `simple_linalg`
- `simple_rnd`
- `simple_stat`
- `simple_timer`

Public symbols:
- `flex_traj_denoise` — subroutine
- `flex_traj_fit_band` — function
- `flex_traj_noise_spectrum` — subroutine

Private symbols:
- `coeff_sets` — subroutine
- `cross_validate` — subroutine
- `cv_direction` — subroutine
- `denoise_set` — subroutine
- `gemm_nn` — subroutine
- `kernel_means` — subroutine
- `latents` — subroutine
- `log_time` — subroutine
- `pick` — function
- `preimages` — subroutine
- `probe_blas` — subroutine
- `whitened_latents` — type

Local symbols (enclosing procedure scope):
- `cv_direction::cv_block` — subroutine
- `cv_direction::kernel_col` — subroutine
- `kernel_means::mean_chunk` — subroutine
- `preimages::kernel_row` — subroutine

---
## Module: simple_forked_process

Files:
- `utils/simple_forked_process.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `simple_cmdline`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_memory_monitor`
- `simple_string`
- `simple_syslib`
- `unix`

Public symbols:
- `forked_process` — type
- `is_sig_ign` — function
- `sig_ign` — function

Private symbols:
- `await_final_status` — subroutine
- `destroy` — subroutine
- `execute_test` — subroutine
- `get_exit_code` — function
- `get_failtime` — function
- `get_pid` — function
- `get_queuetime` — function
- `get_starttime` — function
- `get_stoptime` — function
- `has_exited` — function
- `kill_forked_process` — subroutine
- `set_cline` — subroutine
- `signal_group` — subroutine
- `skip` — subroutine
- `start` — subroutine
- `status` — function
- `terminate` — subroutine

Local symbols (enclosing procedure scope):
- `execute_test::sigterm_handler` — subroutine
- `has_exited::c_waitid` — function

---
## Module: simple_forked_process_tester

Files:
- `utils/simple_forked_process_tester.f90`

Uses:
- `simple_cmdline`
- `simple_forked_process`
- `simple_memory_monitor`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_forked_process_tests` — subroutine

Private symbols:
- `test_destroy` — subroutine
- `test_fail_timestamps` — subroutine
- `test_fork_with_running_monitor` — subroutine
- `test_kill` — subroutine
- `test_logfile_redirection` — subroutine
- `test_new_session` — subroutine
- `test_sigint_ignored` — subroutine
- `test_start` — subroutine
- `test_terminate` — subroutine
- `test_timestamps` — subroutine

---
## Module: simple_fplane_grid

Files:
- `main/image/simple_fplane_grid.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_math_ft`

Public symbols:
- `fplane_grid` — type

Private symbols:
- `has_sigma2` — function
- `kill` — subroutine
- `matches` — function
- `new` — subroutine
- `set_sigma2` — subroutine

---
## Module: simple_frozen_accum

Files:
- `main/volume/simple_frozen_accum.f90`

Uses:
- `simple_core_module_api`
- `simple_imgfile`
- `simple_reconstructor`
- `simple_refine3d_fnames`

Public symbols:
- `carried_bundle` — type
- `frozen_accum` — type
- `read_carried_membership` — subroutine
- `write_carried_membership` — subroutine

Private symbols:
- `add_gridding_set` — subroutine
- `bundle_store` — function
- `component` — function
- `file_size` — function
- `get_nfrozen_state` — function
- `gridding_set_status` — subroutine
- `import_gridding_set` — subroutine
- `kill` — subroutine
- `kill_bundle` — subroutine
- `load` — subroutine
- `new` — subroutine
- `read` — subroutine
- `read_bundle` — subroutine
- `set_dir` — subroutine
- `set_manifest` — function
- `set_weighting` — subroutine
- `validate` — subroutine
- `validate_bundle` — subroutine
- `write` — subroutine
- `write_bundle` — subroutine
- `write_gridding_set` — subroutine
- `write_set_manifest` — subroutine

Local symbols (enclosing procedure scope):
- `add_gridding_set::read_one` — subroutine

---
## Module: simple_frozen_accum_tester

Files:
- `main/volume/simple_frozen_accum_tester.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_fileio`
- `simple_forked_process`
- `simple_fplane_grid`
- `simple_frozen_accum`
- `simple_image`
- `simple_memoize_ft_maps`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_reconstructor`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_string`
- `simple_sym`
- `simple_syslib`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_frozen_accum_tests` — subroutine

Private symbols:
- `accumulator_params` — subroutine
- `delete_gridding_set` — subroutine
- `make_context` — subroutine
- `make_images` — subroutine
- `make_orientations` — subroutine
- `max_rel_cdiff` — function
- `max_rel_rdiff` — function
- `read_rho_into` — subroutine
- `read_text` — subroutine
- `rho_reader` — type
- `same_text` — function
- `test_carried_chain` — subroutine
- `test_context_refusals` — subroutine
- `test_context_round_trip` — subroutine
- `test_gridding_union_equals_sum` — subroutine
- `test_rho_header` — subroutine
- `truncate_file` — subroutine
- `write_text` — subroutine

Local symbols (enclosing procedure scope):
- `test_gridding_union_equals_sum::round_trip` — subroutine

---
## Module: simple_fsc

Files:
- `utils/filter/simple_fsc.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `fsc_area_score_result` — type
- `phase_rand_fsc` — subroutine
- `plot_fsc` — subroutine
- `plot_fsc2` — subroutine
- `plot_fsc_area_score` — subroutine
- `plot_phrand_fsc` — subroutine
- `write_fsc_area_score_outputs` — subroutine

Private symbols:
- `calc_fsc_area_score` — subroutine
- `fibonacci_sphere_dirs` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `score_cfsc_curve` — subroutine

Local symbols (enclosing procedure scope):
- `plot_fsc_area_score::add_crossing_bars` — subroutine
- `plot_fsc_area_score::add_curve` — subroutine
- `plot_fsc_area_score::add_curve_segment` — subroutine
- `plot_fsc_area_score::add_filled_band` — subroutine
- `plot_fsc_area_score::add_filled_band_segment` — subroutine
- `plot_fsc_area_score::add_filled_rect` — subroutine
- `plot_fsc_area_score::add_invisible_bounds` — subroutine
- `plot_fsc_area_score::add_segment` — subroutine
- `plot_fsc_area_score::bar_edges` — subroutine

---
## Module: simple_ft_expanded

Files:
- `main/image/simple_ft_expanded.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `ft_expanded` — type
- `ftexp_transfmat_init` — subroutine
- `ftexp_transfmat_kill` — subroutine

Private symbols:
- `add` — subroutine
- `add2cmat` — subroutine
- `add_uncond` — subroutine
- `alloc_and_calc_tmp_cmat12` — subroutine
- `copy` — subroutine
- `corr` — function
- `corr_fdfshifted_cost_8` — subroutine
- `corr_gshifted_cost_8` — subroutine
- `corr_normalize_dp` — subroutine
- `corr_normalize_sp` — subroutine
- `corr_shifted_cost_8` — function
- `corr_unnorm` — function
- `corr_unnorm_serial` — function
- `dealloc_tmp_cmat12` — subroutine
- `div` — subroutine
- `get_flims` — function
- `get_hp_lp` — subroutine
- `get_sumsq` — function
- `kill` — subroutine
- `new` — subroutine
- `normalize_mat` — subroutine
- `set_cmat` — subroutine
- `shift` — subroutine
- `shift_and_add` — subroutine
- `subtr` — subroutine
- `zero` — subroutine

---
## Module: simple_ftexp_shsrch

Files:
- `main/image/simple_ftexp_shsrch.f90`

Uses:
- `simple_core_module_api`
- `simple_ft_expanded`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `ftexp_shsrch` — type

Private symbols:
- `ftexp_shsrch_corr_shifted_8` — function
- `ftexp_shsrch_cost_8` — function
- `ftexp_shsrch_fdfcost_8` — subroutine
- `ftexp_shsrch_gcost_8` — subroutine
- `ftexp_shsrch_kill` — subroutine
- `ftexp_shsrch_minimize` — function
- `ftexp_shsrch_new` — subroutine
- `set_factr_pgtol` — subroutine
- `set_shsrch_tol` — subroutine

---
## Module: simple_ftexp_shsrch_tester

Files:
- `main/image/simple_ftexp_shsrch_tester.f90`

Uses:
- `simple_defs`
- `simple_ft_expanded`
- `simple_ftexp_shsrch`
- `simple_image`
- `simple_rnd`
- `simple_test_utils`

Public symbols:
- `run_all_ftexp_shsrch_tests` — subroutine

Private symbols:
- `test_correlator` — subroutine
- `test_optimiser` — subroutine

---
## Module: simple_ftiter

Files:
- `main/image/simple_ftiter.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `ftiter` — type

Private symbols:
- `comp_addr_logi` — function
- `comp_addr_phys1` — function
- `comp_addr_phys2` — function
- `comp_addr_phys3` — function
- `constructor` — function
- `dynfind` — function
- `get_find` — function
- `get_ldim` — function
- `get_lfny` — function
- `get_lhp` — function
- `get_lp` — function
- `get_smpd` — function
- `get_spat_freq` — function
- `loop_lims` — function
- `new` — subroutine
- `reset` — subroutine

---
## Module: simple_ftiter_tester

Files:
- `main/image/simple_ftiter_tester.f90`

Uses:
- `simple_defs`
- `simple_ftiter`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_ftiter_tests` — subroutine

Private symbols:
- `box_tag` — function
- `test_address_maps` — subroutine
- `test_logical_limits` — subroutine
- `test_lowpass_limits` — subroutine
- `test_reset` — subroutine
- `test_resolution` — subroutine

Local symbols (enclosing procedure scope):
- `test_logical_limits::check` — subroutine

---
## Module: simple_gauran_tester

Files:
- `main/image/simple_gauran_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_image`
- `simple_test_utils`

Public symbols:
- `run_all_gauran_tests` — subroutine

Private symbols:
- `make_signal` — subroutine
- `pcorr` — function
- `pmean` — function
- `pvar` — function
- `test_add_gauran_snr` — subroutine
- `test_gauran_moments` — subroutine
- `test_replay_and_independence` — subroutine

Local symbols (enclosing procedure scope):
- `test_add_gauran_snr::snr_of` — function

---
## Module: simple_gen_pspecs_and_thumbs_strategy

Files:
- `main/strategies/parallelization/simple_gen_pspecs_and_thumbs_strategy.f90`

Uses:
- `simple_binoris_io`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_parameters`
- `simple_pspec_thumb_iter`
- `simple_qsys_env`
- `simple_sp_project`

Public symbols:
- `create_gen_pspecs_and_thumbs_strategy` — function
- `gen_pspecs_and_thumbs_distr_strategy` — type
- `gen_pspecs_and_thumbs_inmem_strategy` — type
- `gen_pspecs_and_thumbs_strategy` — type

Private symbols:
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_end_message` — function
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `endmsg_interface` — function
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_end_message` — function
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `set_gen_pspecs_and_thumbs_defaults` — subroutine

---
## Module: simple_gmm

Files:
- `utils/clustering/simple_gmm.f90`

Uses:
- `simple_core_module_api`
- `simple_linalg`
- `simple_stat`

Public symbols:
- `gmm` — type

Private symbols:
- `fit` — subroutine
- `get_bic` — function
- `get_cov` — subroutine
- `get_icl` — function
- `get_labels` — subroutine
- `get_loglik` — function
- `get_mass` — subroutine
- `get_means` — subroutine
- `get_niters` — function
- `get_pairsep` — subroutine
- `get_pi` — subroutine
- `get_resp` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `order_by_population` — subroutine

---
## Module: simple_gmm_tester

Files:
- `utils/clustering/simple_gmm_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_gmm`
- `simple_rnd`
- `simple_test_utils`

Public symbols:
- `run_all_gmm_tests` — subroutine

Private symbols:
- `test_gmm_pi_floor` — subroutine
- `test_gmm_recovery` — subroutine
- `two_blobs` — subroutine

---
## Module: simple_gpu_utils

Files:
- `utils/simple_gpu_utils.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_core_module_api`

Public symbols:
- `cublasCreate` — function
- `cublasDestroy` — function
- `cublasSgemm` — function
- `cuda_device_prop` — type
- `cudaDeviceSynchronize` — function
- `cudaGetDeviceProperties` — function
- `cudaMemGetInfo` — function
- `cudaSetDevice` — function
- `cufftDestroy` — function
- `cufftExecC2R` — function
- `cufftExecR2C` — function
- `cufftPlan1d` — function
- `cufftPlan2d` — function
- `cufftPlan3d` — function
- `cufftPlanMany` — function
- `print_gpu_specs` — subroutine
- `set_offload_device_cline` — subroutine
- `set_offload_device_id` — subroutine

---
## Module: simple_gridding

Files:
- `main/interp/simple_gridding.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_projector`

Public symbols:
- `deapodize3D_inplace` — subroutine
- `kb_stencil_centered_crop_inv_envelope_1d` — subroutine
- `kb_stencil_envelope_1d` — subroutine
- `kb_stencil_inv_envelope_1d` — subroutine
- `prep2D_inv_instrfun4mul` — function
- `prep3D_inv_kbenvelope4mul` — function

---
## Module: simple_gui_assembler

Files:
- `utils/gui/simple_gui_assembler.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_gui_metadata_base`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_histogram`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_optics_group`
- `simple_gui_metadata_stream_initial_analysis`
- `simple_gui_metadata_stream_optics_assignment`
- `simple_gui_metadata_stream_particle_sieving`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_stream_pool2d`
- `simple_gui_metadata_stream_preprocess`
- `simple_gui_metadata_stream_snapshot`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_gui_metadata_timeplot`
- `simple_gui_metadata_vol3d`
- `simple_string`
- `unix`

Public symbols:
- `gui_assembler` — type
- `gui_stage_status` — type

Private symbols:
- `add_assigned_cavgs2D` — subroutine
- `add_assigned_histograms` — subroutine
- `add_assigned_micrographs` — subroutine
- `add_assigned_optics_groups` — subroutine
- `add_assigned_timeplots` — subroutine
- `add_item` — subroutine
- `assemble_batch_heartbeat` — subroutine
- `assemble_batch_metadata` — subroutine
- `assemble_stream_heartbeat` — subroutine
- `assemble_stream_initial_analysis` — subroutine
- `assemble_stream_initial_picking` — subroutine
- `assemble_stream_optics_assignment` — subroutine
- `assemble_stream_particle_sieving` — subroutine
- `assemble_stream_pool2D` — subroutine
- `assemble_stream_preprocess` — subroutine
- `assemble_stream_reference_picking` — subroutine
- `assemble_stream_solve3D_multistate` — subroutine
- `clear_hashes` — subroutine
- `close_list` — subroutine
- `commit_section` — subroutine
- `is_associated` — function
- `kill` — subroutine
- `new` — subroutine
- `open_list` — subroutine
- `open_section` — function
- `set_stoptime` — subroutine
- `to_string` — function

Local symbols (enclosing procedure scope):
- `assemble_stream_heartbeat::add_stage_status` — subroutine
- `assemble_stream_solve3D_multistate::add_state_tiles` — subroutine

---
## Module: simple_gui_assembler_tester

Files:
- `utils/gui/simple_gui_assembler_tester.f90`

Uses:
- `simple_gui_assembler`
- `simple_gui_metadata_api`
- `simple_gui_project_builder`
- `simple_sp_project`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_gui_assembler_tests` — subroutine

Private symbols:
- `heartbeat_mismatch` — function
- `json_count` — function
- `json_int` — function
- `json_str_field` — function
- `parse_payload` — subroutine
- `test_batch_heartbeat` — subroutine
- `test_clear_hashes` — subroutine
- `test_initial_analysis` — subroutine
- `test_initial_picking` — subroutine
- `test_new_kill` — subroutine
- `test_new_kill_reuse` — subroutine
- `test_optics_assignment` — subroutine
- `test_particle_sieving` — subroutine
- `test_pool2D` — subroutine
- `test_preprocess` — subroutine
- `test_project` — subroutine
- `test_reference_picking` — subroutine
- `test_set_stoptime` — subroutine
- `test_solve3D_multistate` — subroutine
- `test_stream_heartbeat` — subroutine

---
## Module: simple_gui_communicator

Files:
- `utils/gui/simple_gui_communicator.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_core_module_api`
- `simple_gui_assembler`
- `simple_gui_metadata_project`
- `simple_gui_metadata_types`
- `simple_gui_project_builder`
- `simple_http_post`
- `simple_parameters`
- `simple_sp_project`
- `unix`

Public symbols:
- `gui_communicator` — type

Private symbols:
- `add_metadata_projfile` — subroutine
- `add_metadata_spproj` — subroutine
- `communication_worker` — subroutine
- `gui_comm_args` — type
- `kill_gui_communicator` — subroutine
- `new` — subroutine

Local symbols (enclosing procedure scope):
- `communication_worker::safe_destroy_json_ptr` — subroutine

---
## Module: simple_gui_metadata_api

Files:
- `main/apis/simple_gui_metadata_api.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_core_module_api`
- `simple_gui_metadata_base`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_histogram`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_optics_group`
- `simple_gui_metadata_project`
- `simple_gui_metadata_ptcl`
- `simple_gui_metadata_stream_initial_analysis`
- `simple_gui_metadata_stream_optics_assignment`
- `simple_gui_metadata_stream_particle_sieving`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_stream_pool2d`
- `simple_gui_metadata_stream_preprocess`
- `simple_gui_metadata_stream_snapshot`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_timeplot`
- `simple_gui_metadata_types`
- `simple_gui_metadata_vol3d`

---
## Module: simple_gui_metadata_base

Files:
- `utils/gui/metadata/simple_gui_metadata_base.f90`

Uses:
- `json_module`
- `simple_error`

Public symbols:
- `gui_metadata_base` — type

Private symbols:
- `assigned` — function
- `initialized` — function
- `jsonise` — function
- `kill` — subroutine
- `new` — subroutine
- `serialise` — subroutine
- `set_assigned` — subroutine
- `type` — function

---
## Module: simple_gui_metadata_cavg2D

Files:
- `utils/gui/metadata/simple_gui_metadata_cavg2D.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`

Public symbols:
- `gui_metadata_cavg2D` — type
- `sprite_sheet_pos` — type

Private symbols:
- `get` — function
- `get_i` — function
- `get_i_max` — function
- `get_idx` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_histogram

Files:
- `utils/gui/metadata/simple_gui_metadata_histogram.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`

Public symbols:
- `gui_metadata_histogram` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_micrograph

Files:
- `utils/gui/metadata/simple_gui_metadata_micrograph.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`

Public symbols:
- `gui_metadata_micrograph` — type

Private symbols:
- `clear_coordinates` — subroutine
- `get` — function
- `get_i` — function
- `get_i_max` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine
- `set_coordinate` — subroutine

---
## Module: simple_gui_metadata_optics_group

Files:
- `utils/gui/metadata/simple_gui_metadata_optics_group.f90`

Uses:
- `json_module`
- `simple_error`
- `simple_gui_metadata_base`

Public symbols:
- `gui_metadata_optics_group` — type

Private symbols:
- `get` — function
- `get_i` — function
- `get_i_max` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_project

Files:
- `utils/gui/metadata/simple_gui_metadata_project.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_ptcl`
- `simple_gui_metadata_vol3d`
- `simple_string`
- `simple_string_utils`
- `unix`

Public symbols:
- `gui_metadata_project` — type

Private symbols:
- `get` — function
- `gui_metadata_cavg2D_stage` — type
- `gui_metadata_vol3D_stage` — type
- `jsonise_override` — function
- `kill_override` — subroutine
- `serialise_override` — subroutine
- `set_cavgs2D` — subroutine
- `set_micrographs` — subroutine
- `set_movies` — subroutine
- `set_particles` — subroutine
- `set_project` — subroutine
- `set_summary` — subroutine
- `set_vols3D` — subroutine
- `stage_slot` — subroutine
- `start_update` — subroutine

---
## Module: simple_gui_metadata_ptcl

Files:
- `utils/gui/metadata/simple_gui_metadata_ptcl.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_gui_metadata_cavg2d`
- `simple_string`

Public symbols:
- `gui_metadata_ptcl` — type

Private symbols:
- `get` — function
- `get_i` — function
- `get_i_max` — function
- `get_idx` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_initial_analysis

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_initial_analysis.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_initial_analysis` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_optics_assignment

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_optics_assignment.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_optics_assignment` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_particle_sieving

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_particle_sieving.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_particle_sieving` — type

Private symbols:
- `clear_selection` — subroutine
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine
- `set_selection` — subroutine

---
## Module: simple_gui_metadata_stream_picking

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_picking.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_picking` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_pool2D

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_pool2D.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_pool2D` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine
- `set_user_input` — subroutine

---
## Module: simple_gui_metadata_stream_preprocess

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_preprocess.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`

Public symbols:
- `gui_metadata_stream_preprocess` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_snapshot

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_snapshot.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_snapshot` — type

Private symbols:
- `get` — function
- `get_states` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_solve3D_multistate

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_solve3D_multistate.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_solve3D_multistate` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine
- `set_state_stats` — subroutine

---
## Module: simple_gui_metadata_stream_update

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_update.f90`

Uses:
- `simple_error`
- `simple_gui_metadata_base`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_string`

Public symbols:
- `gui_metadata_stream_update` — type

Private symbols:
- `get_astigmatism_update` — function
- `get_ctfres_update` — function
- `get_icescore_update` — function
- `get_mskdiam2D_update` — function
- `get_pickrefs_cycle` — function
- `get_pickrefs_selection` — function
- `get_pickrefs_selection_length` — function
- `get_snapshot2D_update` — subroutine
- `get_snapshot3D_update` — subroutine
- `has_snapshot2D_update` — function
- `has_snapshot3D_update` — function
- `kill_override` — subroutine
- `set_astigmatism_update` — subroutine
- `set_ctfres_update` — subroutine
- `set_icescore_update` — subroutine
- `set_mskdiam2D_update` — subroutine
- `set_pickrefs_cycle` — subroutine
- `set_pickrefs_selection` — subroutine
- `set_snapshot2D_update` — subroutine
- `set_snapshot3D_update` — subroutine

---
## Module: simple_gui_metadata_tester

Files:
- `utils/gui/metadata/simple_gui_metadata_tester.f90`

Uses:
- `simple_gui_metadata_api`
- `simple_test_utils`

Public symbols:
- `run_all_gui_metadata_tests` — subroutine

Private symbols:
- `json_text` — function
- `test_jsonise_cavg2D` — subroutine
- `test_jsonise_histogram` — subroutine
- `test_jsonise_micrograph` — subroutine
- `test_jsonise_optics_group` — subroutine
- `test_jsonise_ptcl` — subroutine
- `test_jsonise_stream_initial_analysis` — subroutine
- `test_jsonise_stream_initial_picking` — subroutine
- `test_jsonise_stream_optics_assignment` — subroutine
- `test_jsonise_stream_particle_sieving` — subroutine
- `test_jsonise_stream_pool2D` — subroutine
- `test_jsonise_stream_pool2D_snapshot` — subroutine
- `test_jsonise_stream_preprocess` — subroutine
- `test_jsonise_stream_reference_picking` — subroutine
- `test_jsonise_stream_solve3D_multistate` — subroutine
- `test_jsonise_stream_update` — subroutine
- `test_jsonise_timeplot` — subroutine
- `test_jsonise_vol3D` — subroutine
- `test_kill_resets_fields` — subroutine
- `test_new_kill` — subroutine
- `test_serialise` — subroutine
- `test_serialise_cavg2D` — subroutine
- `test_serialise_histogram` — subroutine
- `test_serialise_micrograph` — subroutine
- `test_serialise_optics_group` — subroutine
- `test_serialise_ptcl` — subroutine
- `test_serialise_stream_initial_analysis` — subroutine
- `test_serialise_stream_initial_picking` — subroutine
- `test_serialise_stream_optics_assignment` — subroutine
- `test_serialise_stream_particle_sieving` — subroutine
- `test_serialise_stream_pool2D` — subroutine
- `test_serialise_stream_pool2D_snapshot` — subroutine
- `test_serialise_stream_preprocess` — subroutine
- `test_serialise_stream_reference_picking` — subroutine
- `test_serialise_stream_solve3D_multistate` — subroutine
- `test_serialise_stream_update` — subroutine
- `test_serialise_timeplot` — subroutine
- `test_serialise_vol3D` — subroutine
- `test_serialise_vol3D_tiles` — subroutine
- `test_set_get_cavg2D` — subroutine
- `test_set_get_histogram` — subroutine
- `test_set_get_micrograph` — subroutine
- `test_set_get_optics_group` — subroutine
- `test_set_get_ptcl` — subroutine
- `test_set_get_stream_initial_analysis` — subroutine
- `test_set_get_stream_initial_picking` — subroutine
- `test_set_get_stream_optics_assignment` — subroutine
- `test_set_get_stream_particle_sieving` — subroutine
- `test_set_get_stream_pool2D` — subroutine
- `test_set_get_stream_pool2D_snapshot` — subroutine
- `test_set_get_stream_preprocess` — subroutine
- `test_set_get_stream_reference_picking` — subroutine
- `test_set_get_stream_solve3D_multistate` — subroutine
- `test_set_get_stream_update` — subroutine
- `test_set_get_timeplot` — subroutine
- `test_set_get_vol3D` — subroutine
- `test_stream_snapshot3D` — subroutine

---
## Module: simple_gui_metadata_timeplot

Files:
- `utils/gui/metadata/simple_gui_metadata_timeplot.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`

Public symbols:
- `gui_metadata_timeplot` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `set` — subroutine

---
## Module: simple_gui_metadata_types

Files:
- `utils/gui/metadata/simple_gui_metadata_types.f90`

---
## Module: simple_gui_metadata_utils

Files:
- `utils/gui/metadata/simple_gui_metadata_utils.f90`

Uses:
- `simple_gui_metadata_base`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_histogram`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_optics_group`
- `simple_gui_metadata_stream_initial_analysis`
- `simple_gui_metadata_stream_optics_assignment`
- `simple_gui_metadata_stream_particle_sieving`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_stream_pool2d`
- `simple_gui_metadata_stream_preprocess`
- `simple_gui_metadata_stream_snapshot`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_timeplot`
- `simple_gui_metadata_vol3d`

Public symbols:
- `max_metadata_size` — function

---
## Module: simple_gui_metadata_vol3D

Files:
- `utils/gui/metadata/simple_gui_metadata_vol3D.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_gui_metadata_cavg2d`
- `simple_string`

Public symbols:
- `gui_metadata_vol3D` — type

Private symbols:
- `get` — function
- `get_fsc` — function
- `get_i` — function
- `get_i_max` — function
- `get_oridist` — function
- `get_state` — function
- `jsonise_override` — function
- `kill_override` — subroutine
- `serialise_override` — subroutine
- `set` — subroutine
- `set_fsc` — subroutine
- `set_minmax` — subroutine
- `set_oridist` — subroutine
- `set_reprojtiles` — subroutine

---
## Module: simple_gui_project_builder

Files:
- `utils/gui/simple_gui_project_builder.f90`

Uses:
- `simple_defs`
- `simple_defs_fname`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_fileio`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_project`
- `simple_gui_metadata_ptcl`
- `simple_gui_metadata_types`
- `simple_gui_metadata_vol3d`
- `simple_gui_utils`
- `simple_image`
- `simple_imghead`
- `simple_math`
- `simple_math_ft`
- `simple_motion_gain_helpers`
- `simple_nrtxtfile`
- `simple_oris_utils`
- `simple_procimgstk`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `gui_project_builder` — type

Private symbols:
- `add_cavgs2D` — subroutine
- `add_micrographs` — subroutine
- `add_movies` — subroutine
- `add_particles` — subroutine
- `add_vols3D` — subroutine
- `apply` — subroutine
- `build` — subroutine
- `kill` — subroutine
- `stage_jpeg` — function

---
## Module: simple_gui_utils

Files:
- `utils/gui/simple_gui_utils.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_stack_io`

Public symbols:
- `mic2thumb` — subroutine
- `mrc2jpeg_tiled` — subroutine

---
## Module: simple_hac

Files:
- `utils/clustering/simple_hac.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `hac` — type

Private symbols:
- `cluster` — subroutine
- `find_medoids` — subroutine
- `get_history` — subroutine
- `get_nclust` — function
- `kill` — subroutine
- `label_clusters` — subroutine
- `merge_clusters` — subroutine
- `new` — subroutine

---
## Module: simple_hac_tester

Files:
- `utils/clustering/simple_hac_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_hac`
- `simple_test_utils`

Public symbols:
- `run_all_hac_tests` — subroutine

Private symbols:
- `line_dmat` — function
- `run_line` — subroutine
- `test_hac_complete_linkage` — subroutine
- `test_hac_dense_group_stays_whole` — subroutine
- `test_hac_edges` — subroutine
- `test_hac_mask` — subroutine
- `test_hac_population_weighted_update` — subroutine
- `test_hac_separated_groups` — subroutine
- `test_hac_threshold_stop` — subroutine

---
## Module: simple_halfmap_diagnostics

Files:
- `main/volume/simple_halfmap_diagnostics.f90`

Uses:
- `simple_core_module_api`
- `simple_fileio`
- `simple_fsc`
- `simple_image`
- `simple_image_msk`
- `simple_parameters`

Public symbols:
- `copy_support_provenance` — subroutine
- `evaluate_halfmap_pair` — subroutine
- `halfmap_diagnostics_result` — type
- `read_support_provenance` — subroutine
- `remove_support_provenance` — subroutine
- `rename_support_provenance` — subroutine
- `support_provenance_fname` — function
- `write_base_halfmap_fsc` — subroutine
- `write_halfmap_diagnostics` — subroutine
- `write_support_provenance` — subroutine

Private symbols:
- `kill_halfmap_diagnostics_result` — subroutine

---
## Module: simple_halfmap_enhancer

Files:
- `main/volume/simple_halfmap_enhancer.f90`

Uses:
- `ieee_arithmetic`
- `simple_butterworth`
- `simple_defs`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_image`
- `simple_math_ft`
- `simple_stat`

Public symbols:
- `halfmap_enhancer` — type

Private symbols:
- `count_outer` — function
- `enhance` — subroutine
- `estimate_noise` — subroutine
- `get_noise_stats` — subroutine
- `get_window` — function
- `kill` — subroutine
- `local_envelope_confidence` — function
- `new` — subroutine
- `regularize_fsc` — subroutine
- `window_starts` — function

---
## Module: simple_halfmap_enhancer_tester

Files:
- `main/volume/simple_halfmap_enhancer_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_halfmap_enhancer`
- `simple_image`
- `simple_test_utils`

Public symbols:
- `run_all_halfmap_enhancer_tests` — subroutine

Private symbols:
- `array_correlation` — function
- `blob_map` — function
- `radius_sq` — function
- `run_enhancer` — subroutine
- `test_constant_pair` — subroutine
- `test_correlated_background_pair` — subroutine
- `test_noisy_blob_pair` — subroutine
- `test_unsupported_pair` — subroutine
- `white_noise` — function

---
## Module: simple_hash

Files:
- `utils/structs/simple_hash.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `hash` — type

Private symbols:
- `alloc_hash` — subroutine
- `constructor_by_sz` — function
- `constructor_default` — function
- `copy` — subroutine
- `dealloc_hash` — subroutine
- `delete` — subroutine
- `get` — function
- `get_dp` — function
- `get_key` — function
- `get_keys` — function
- `get_value_at` — function
- `getter_dp` — subroutine
- `getter_int` — subroutine
- `getter_real` — subroutine
- `hash2str` — function
- `hash_strlen` — function
- `isthere` — function
- `kill` — subroutine
- `lookup` — function
- `new_by_sz` — subroutine
- `new_default` — subroutine
- `print` — subroutine
- `push_dp` — subroutine
- `push_int` — subroutine
- `realloc_hash` — subroutine
- `set_dp` — subroutine
- `set_int` — subroutine
- `set_real` — subroutine
- `size_of` — function

---
## Module: simple_hash_tester

Files:
- `utils/structs/simple_hash_tester.f90`

Uses:
- `simple_defs`
- `simple_hash`
- `simple_string`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_hash_tests` — subroutine

Private symbols:
- `test_constructor_and_kill` — subroutine
- `test_copy_and_realloc` — subroutine
- `test_delete` — subroutine
- `test_get_keys_and_get_str` — subroutine
- `test_getters_numeric` — subroutine
- `test_hash2str_and_strlen` — subroutine
- `test_lookup_and_isthere` — subroutine
- `test_push_and_set` — subroutine

---
## Module: simple_histogram

Files:
- `utils/math/simple_histogram.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `histogram` — type

Private symbols:
- `add` — subroutine
- `div` — subroutine
- `entropy` — function
- `find_hill` — function
- `find_next_valley` — function
- `get` — function
- `get_bin` — function
- `get_nbins` — function
- `get_x` — function
- `HD` — function
- `hmode` — function
- `JSD` — function
- `kill` — subroutine
- `KLD` — function
- `mean` — function
- `moments` — function
- `new_by_nbins` — subroutine
- `new_from_hist` — subroutine
- `new_from_img` — subroutine
- `new_from_vec` — subroutine
- `new_from_vec_nbins` — subroutine
- `npeaks` — function
- `plot` — subroutine
- `quantize` — subroutine
- `reset` — subroutine
- `skew` — function
- `smooth` — subroutine
- `topdf` — subroutine
- `TVD` — function
- `update` — subroutine
- `variance` — function
- `zero` — subroutine

---
## Module: simple_http_post

Files:
- `utils/comm/simple_http_post.f90`

Uses:
- `curl`
- `iso_c_binding`
- `simple_error`
- `simple_string`

Public symbols:
- `http_post` — type
- `http_response` — type

Private symbols:
- `initialised` — function
- `kill` — subroutine
- `new` — subroutine
- `request` — function
- `response_callback` — function

---
## Module: simple_http_post_tester

Files:
- `utils/comm/simple_http_post_tester.f90`

Uses:
- `iso_c_binding`
- `simple_http_post`
- `simple_ipc_tcp_socket_server`
- `simple_string`
- `simple_string_utils`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_http_post_tests` — subroutine

Private symbols:
- `content_length` — function
- `destroy_listener_args` — subroutine
- `http_listener` — subroutine
- `init_listener_args` — subroutine
- `loopback_url` — function
- `test_create_and_kill` — subroutine
- `test_request_no_body` — subroutine
- `test_request_with_body` — subroutine
- `test_request_without_listener` — subroutine
- `test_response_reset` — subroutine

---
## Module: simple_image

Files:
- `main/image/simple_image.f90`
- `main/image/simple_image_access.f90`
- `main/image/simple_image_arith.f90`
- `main/image/simple_image_calc.f90`
- `main/image/simple_image_checks.f90`
- `main/image/simple_image_core.f90`
- `main/image/simple_image_ctf.f90`
- `main/image/simple_image_fft.f90`
- `main/image/simple_image_filt.f90`
- `main/image/simple_image_freq_anal.f90`
- `main/image/simple_image_geom.f90`
- `main/image/simple_image_io.f90`
- `main/image/simple_image_norm.f90`
- `main/image/simple_image_ops.f90`
- `main/image/simple_image_polar.f90`
- `main/image/simple_image_seg.f90`
- `main/image/simple_image_vis.f90`

Uses:
- `gnufor2`
- `ieee_arithmetic`
- `simple_core_module_api`
- `simple_ctf`
- `simple_fftw3`
- `simple_fplane_grid`
- `simple_ftiter`
- `simple_imgfile`
- `simple_imghead`
- `simple_jpg`
- `simple_math`
- `simple_math_ctf`
- `simple_memoize_ft_maps`
- `simple_neighs`
- `simple_online_var`
- `simple_parameters`
- `simple_tent_smooth`
- `simple_winfuns`

Public symbols:
- `image` — type
- `image_stack` — type
- `unmemoize_mask_coords` — subroutine
- `unmemoize_powspec_coords` — subroutine

Private symbols:
- `acf` — subroutine
- `acf` — subroutine
- `add_cmat_at` — subroutine
- `add_cmat_at` — subroutine
- `add_cmats_to_cmat` — subroutine
- `add_cmats_to_cmat` — subroutine
- `add_cmplx_at` — subroutine
- `add_cmplx_at` — subroutine
- `add_cmplx_at_kdiv` — subroutine
- `add_cmplx_at_kdiv` — subroutine
- `add_dble_cmat2mat` — subroutine
- `add_dble_cmat2mat` — subroutine
- `add_gauran` — subroutine
- `add_gauran` — subroutine
- `add_img` — subroutine
- `add_img` — subroutine
- `add_mat2cmat_cmplx` — subroutine
- `add_mat2cmat_cmplx` — subroutine
- `add_mat2cmat_real` — subroutine
- `add_mat2cmat_real` — subroutine
- `add_real_at` — subroutine
- `add_real_at` — subroutine
- `add_rmat2mat_workshare` — subroutine
- `add_rmat2mat_workshare` — subroutine
- `add_scalar` — subroutine
- `add_scalar` — subroutine
- `add_window` — subroutine
- `add_window` — subroutine
- `add_workshare_img` — subroutine
- `add_workshare_img` — subroutine
- `add_workshare_weighted` — subroutine
- `add_workshare_weighted` — subroutine
- `add_zero2mask` — subroutine
- `add_zero2mask` — subroutine
- `addition` — function
- `addition` — function
- `addition_const_real` — function
- `addition_const_real` — function
- `allocate_cmat` — function
- `allocate_cmat` — function
- `apply_bfac` — subroutine
- `apply_bfac` — subroutine
- `apply_ctf` — subroutine
- `apply_ctf` — subroutine
- `apply_ctf_wiener` — subroutine
- `apply_ctf_wiener` — subroutine
- `apply_ctf_wiener_wpad` — subroutine
- `apply_ctf_wiener_wpad` — subroutine
- `apply_ctf_wpad` — subroutine
- `apply_ctf_wpad` — subroutine
- `apply_dose_weighing` — subroutine
- `apply_dose_weighing` — subroutine
- `apply_filter_pair` — subroutine
- `apply_filter_pair` — subroutine
- `apply_filter_serial` — subroutine
- `apply_filter_serial` — subroutine
- `apply_filter_single` — subroutine
- `apply_filter_single` — subroutine
- `assign` — subroutine
- `assign` — subroutine
- `assign_c2img` — subroutine
- `assign_c2img` — subroutine
- `assign_r2img` — subroutine
- `assign_r2img` — subroutine
- `avg_loc_sdev` — function
- `avg_loc_sdev` — function
- `bartlett_reg_3D` — subroutine
- `bartlett_reg_3D` — subroutine
- `before_after` — subroutine
- `before_after` — subroutine
- `bin2logical` — function
- `bin2logical` — function
- `bin_inv` — subroutine
- `bin_inv` — subroutine
- `binarize_npix` — subroutine
- `binarize_npix` — subroutine
- `binarize_thres` — subroutine
- `binarize_thres` — subroutine
- `binarize_thres_mask` — subroutine
- `binarize_thres_mask` — subroutine
- `bp` — subroutine
- `bp` — subroutine
- `bpgau2D` — subroutine
- `bpgau2D` — subroutine
- `bpgau3D` — subroutine
- `bpgau3D` — subroutine
- `bs_smooth` — subroutine
- `bs_smooth` — subroutine
- `bwd_ft` — subroutine
- `bwd_ft` — subroutine
- `calc_3D_shape_descriptors` — subroutine
- `calc_3D_shape_descriptors` — subroutine
- `calc_bin_thres` — subroutine
- `calc_bin_thres` — subroutine
- `calc_gradient` — subroutine
- `calc_gradient` — subroutine
- `calc_ice_frac` — subroutine
- `calc_ice_frac` — subroutine
- `calc_ice_score` — subroutine
- `calc_ice_score` — subroutine
- `calc_principal_axes` — subroutine
- `calc_principal_axes` — subroutine
- `calc_shiftcen` — function
- `calc_shiftcen` — function
- `calc_shiftcen_serial` — function
- `calc_shiftcen_serial` — function
- `ccf` — function
- `ccf` — function
- `ccf_into` — subroutine
- `ccf_into` — subroutine
- `cendist` — subroutine
- `cendist` — subroutine
- `center_edge_snr` — function
- `center_edge_snr` — function
- `checkimg4nans` — subroutine
- `checkimg4nans` — subroutine
- `clear_polar_memo` — subroutine
- `clip` — subroutine
- `clip` — subroutine
- `clip_inplace` — subroutine
- `clip_inplace` — subroutine
- `collage` — subroutine
- `collage` — subroutine
- `comp_addr_logi` — function
- `comp_addr_logi` — function
- `comp_addr_phys1` — function
- `comp_addr_phys1` — function
- `comp_addr_phys2` — function
- `comp_addr_phys2` — function
- `comp_addr_phys3` — function
- `comp_addr_phys3` — function
- `conical_fsc` — subroutine
- `conical_fsc` — subroutine
- `conjugate` — subroutine
- `conjugate` — subroutine
- `construct_thread_safe_tmp_imgs` — subroutine
- `construct_thread_safe_tmp_imgs` — subroutine
- `constructor` — function
- `contains_nans` — function
- `contains_nans` — function
- `contrast` — function
- `contrast` — function
- `copy` — subroutine
- `copy` — subroutine
- `copy_fast` — subroutine
- `copy_fast` — subroutine
- `corner` — subroutine
- `corner` — subroutine
- `corners` — subroutine
- `corners` — subroutine
- `corr` — function
- `corr` — function
- `corr_shifted` — function
- `corr_shifted` — function
- `ctf2img` — subroutine
- `ctf2img` — subroutine
- `ctf_dens_correct` — subroutine
- `ctf_dens_correct` — subroutine
- `cure` — subroutine
- `cure` — subroutine
- `cure_outliers` — subroutine
- `cure_outliers` — subroutine
- `cylinder` — subroutine
- `cylinder` — subroutine
- `dampen_pspec_central_cross` — subroutine
- `dampen_pspec_central_cross` — subroutine
- `dead_hot_positions` — function
- `dead_hot_positions` — function
- `density_inoutside` — subroutine
- `density_inoutside` — subroutine
- `disc_img` — subroutine
- `disc_img` — subroutine
- `disc_lmsk` — subroutine
- `disc_lmsk` — subroutine
- `disc_sideview` — subroutine
- `disc_sideview` — subroutine
- `distance_transform` — subroutine
- `distance_transform` — subroutine
- `div_below` — subroutine
- `div_below` — subroutine
- `div_cmat_at_cmplx` — subroutine
- `div_cmat_at_cmplx` — subroutine
- `div_cmat_at_real` — subroutine
- `div_cmat_at_real` — subroutine
- `div_img` — subroutine
- `div_img` — subroutine
- `div_kernel_at` — subroutine
- `div_kernel_at` — subroutine
- `div_real_at` — subroutine
- `div_real_at` — subroutine
- `div_rmat_at` — subroutine
- `div_rmat_at` — subroutine
- `div_scalar` — subroutine
- `div_scalar` — subroutine
- `division` — function
- `division` — function
- `downsample_gain` — subroutine
- `downsample_gain` — subroutine
- `ellipse` — subroutine
- `estimate_background` — subroutine
- `estimate_background` — subroutine
- `euclid_norm` — function
- `euclid_norm` — function
- `even_dims` — function
- `even_dims` — function
- `exists` — function
- `exists` — function
- `expand_ft` — function
- `expand_ft` — function
- `fcomps_below_noise_power_stats` — subroutine
- `fcomps_below_noise_power_stats` — subroutine
- `fcorr_shift` — subroutine
- `fcorr_shift` — subroutine
- `fcorr_shift3D` — subroutine
- `fcorr_shift3D` — subroutine
- `fft_clip_shift` — subroutine
- `fft_clip_shift` — subroutine
- `fft_clip_shift_ctf_flip` — subroutine
- `fft_clip_shift_ctf_flip` — subroutine
- `fft_noshift` — subroutine
- `fft_noshift` — subroutine
- `flip` — subroutine
- `flip` — subroutine
- `frame_ref_phase_correlation` — subroutine
- `frame_ref_phase_correlation` — subroutine
- `frc_pspec` — subroutine
- `frc_pspec` — subroutine
- `fsc` — subroutine
- `fsc` — subroutine
- `fsvar` — subroutine
- `fsvar` — subroutine
- `ft2img` — subroutine
- `ft2img` — subroutine
- `fwd_ft` — subroutine
- `fwd_ft` — subroutine
- `gabor_filter2D` — subroutine
- `gabor_filter2D` — subroutine
- `gauimg2D` — subroutine
- `gauimg2D` — subroutine
- `gauimg3D` — subroutine
- `gauimg3D` — subroutine
- `gauimg_center` — subroutine
- `gauimg_center` — subroutine
- `gauimg_offset` — subroutine
- `gauimg_offset` — subroutine
- `gauran` — subroutine
- `gauran` — subroutine
- `gen_argtransf_comp` — function
- `gen_argtransf_comp` — function
- `gen_fplane4rec` — subroutine
- `gen_fplane4rec` — subroutine
- `generate_orthogonal_reprojs` — subroutine
- `generate_orthogonal_reprojs` — subroutine
- `get` — function
- `get` — function
- `get_array_shape` — function
- `get_array_shape` — function
- `get_box` — function
- `get_box` — function
- `get_cmat` — function
- `get_cmat` — function
- `get_cmat_at` — function
- `get_cmat_at` — function
- `get_cmat_ptr` — subroutine
- `get_cmat_ptr` — subroutine
- `get_cmat_sub` — subroutine
- `get_cmat_sub` — subroutine
- `get_fcomp` — function
- `get_fcomp` — function
- `get_fcomp2D` — function
- `get_fcomp2D` — function
- `get_filtsz` — function
- `get_filtsz` — function
- `get_find` — function
- `get_find` — function
- `get_ldim` — function
- `get_ldim` — function
- `get_lfny` — function
- `get_lfny` — function
- `get_lhp` — function
- `get_lhp` — function
- `get_lp` — function
- `get_lp` — function
- `get_nyq` — function
- `get_nyq` — function
- `get_res` — function
- `get_res` — function
- `get_rmat` — function
- `get_rmat` — function
- `get_rmat_at` — function
- `get_rmat_at` — function
- `get_rmat_ptr` — subroutine
- `get_rmat_ptr` — subroutine
- `get_rmat_ptr_padded` — subroutine
- `get_rmat_ptr_padded` — subroutine
- `get_rmat_sub` — subroutine
- `get_rmat_sub` — subroutine
- `get_shconst` — function
- `get_shconst` — function
- `get_slice` — subroutine
- `get_slice` — subroutine
- `get_smpd` — function
- `get_smpd` — function
- `get_spat_freq` — function
- `get_spat_freq` — function
- `get_subimg` — subroutine
- `get_subimg` — subroutine
- `get_sum_int` — function
- `get_sum_int` — function
- `GLCM` — subroutine
- `GLCM` — subroutine
- `gradient` — subroutine
- `gradient` — subroutine
- `gradients_magnitude` — subroutine
- `gradients_magnitude` — subroutine
- `guinier` — function
- `guinier` — function
- `guinier_bfac` — function
- `guinier_bfac` — function
- `hannw` — function
- `hannw` — function
- `ICM2D` — subroutine
- `ICM2D` — subroutine
- `ifft_mask_fft` — subroutine
- `ifft_mask_fft` — subroutine
- `ifft_mask_pad_fft` — subroutine
- `ifft_mask_pad_fft` — subroutine
- `img2ft` — subroutine
- `img2ft` — subroutine
- `img2spec` — subroutine
- `img2spec` — subroutine
- `init_gridcorr_mats` — subroutine
- `init_gridcorr_mats` — subroutine
- `insert` — subroutine
- `insert` — subroutine
- `insert_lowres` — subroutine
- `insert_lowres` — subroutine
- `insert_lowres_serial` — subroutine
- `insert_lowres_serial` — subroutine
- `inv` — subroutine
- `inv` — subroutine
- `inv_gain` — subroutine
- `inv_gain` — subroutine
- `is_2d` — function
- `is_2d` — function
- `is_3d` — function
- `is_3d` — function
- `is_empty` — function
- `is_empty` — function
- `is_ft` — function
- `is_ft` — function
- `is_wthreads` — function
- `is_wthreads` — function
- `iter_gridcorr` — subroutine
- `iter_gridcorr` — subroutine
- `kill` — subroutine
- `kill` — subroutine
- `kill_thread_safe_tmp_imgs` — subroutine
- `kill_thread_safe_tmp_imgs` — subroutine
- `kurt` — function
- `kurt` — function
- `loc_sdev` — subroutine
- `loc_sdev` — subroutine
- `loc_var` — subroutine
- `loc_var` — subroutine
- `loc_var3D` — subroutine
- `loc_var3D` — subroutine
- `loc_var_masked` — subroutine
- `loc_var_masked` — subroutine
- `logical2bin` — subroutine
- `logical2bin` — subroutine
- `loop_lims` — function
- `loop_lims` — function
- `lp_background` — subroutine
- `lp_background` — subroutine
- `mask2D_hard` — subroutine
- `mask2D_hard` — subroutine
- `mask2D_soft` — subroutine
- `mask2D_soft` — subroutine
- `mask2D_softavg` — subroutine
- `mask2D_softavg` — subroutine
- `mask3D_hard` — subroutine
- `mask3D_hard` — subroutine
- `mask3D_soft` — subroutine
- `mask3D_soft` — subroutine
- `mask3D_softavg` — subroutine
- `mask3D_softavg` — subroutine
- `mask_fft` — subroutine
- `mask_fft` — subroutine
- `mask_pad_fft` — subroutine
- `mask_pad_fft` — subroutine
- `masked_sqsum` — function
- `masked_sqsum` — function
- `masked_subtr_corr` — subroutine
- `masked_subtr_corr` — subroutine
- `masked_subtr_sqsum` — function
- `masked_subtr_sqsum` — function
- `masscen` — subroutine
- `masscen` — subroutine
- `max_abs_padding` — function
- `max_abs_padding` — function
- `mean` — function
- `mean` — function
- `memoize4polarize` — subroutine
- `memoize4polarize` — subroutine
- `memoize4polarize_oversamp` — subroutine
- `memoize4polarize_oversamp` — subroutine
- `memoize_mask_coords` — subroutine
- `memoize_mask_coords` — subroutine
- `memoize_polar_point` — subroutine
- `memoize_powspec_coords` — subroutine
- `memoize_powspec_coords` — subroutine
- `mic2spec` — subroutine
- `mic2spec` — subroutine
- `minmax` — function
- `minmax` — function
- `mirror` — subroutine
- `mirror` — subroutine
- `mul_cmat_at_cmplx` — subroutine
- `mul_cmat_at_cmplx` — subroutine
- `mul_cmat_at_real` — subroutine
- `mul_cmat_at_real` — subroutine
- `mul_cmat_masked` — subroutine
- `mul_cmat_masked` — subroutine
- `mul_cmat_unmasked` — subroutine
- `mul_cmat_unmasked` — subroutine
- `mul_cmplx_at` — subroutine
- `mul_cmplx_at` — subroutine
- `mul_img` — subroutine
- `mul_img` — subroutine
- `mul_img_lp` — subroutine
- `mul_img_lp` — subroutine
- `mul_real_at` — subroutine
- `mul_real_at` — subroutine
- `mul_rmat` — subroutine
- `mul_rmat` — subroutine
- `mul_rmat_at` — subroutine
- `mul_rmat_at` — subroutine
- `mul_scalar` — subroutine
- `mul_scalar` — subroutine
- `multiplication` — function
- `multiplication` — function
- `multiplication_const_int` — function
- `multiplication_const_int` — function
- `multiplication_const_real` — function
- `multiplication_const_real` — function
- `nbackground` — function
- `nbackground` — function
- `neg` — subroutine
- `neg` — subroutine
- `new` — subroutine
- `new` — subroutine
- `nforeground` — function
- `nforeground` — function
- `NLmean2D` — subroutine
- `NLmean2D` — subroutine
- `NLmean2D_eo` — subroutine
- `NLmean2D_eo` — subroutine
- `NLmean3D` — subroutine
- `NLmean3D` — subroutine
- `NLmean3D_eo` — subroutine
- `NLmean3D_eo` — subroutine
- `noisesdev` — function
- `noisesdev` — function
- `norm` — subroutine
- `norm` — subroutine
- `norm4viz` — subroutine
- `norm4viz` — subroutine
- `norm_ext` — subroutine
- `norm_ext` — subroutine
- `norm_minmax` — subroutine
- `norm_minmax` — subroutine
- `norm_noise` — subroutine
- `norm_noise` — subroutine
- `norm_noise_fft` — subroutine
- `norm_noise_fft` — subroutine
- `norm_noise_fft_clip_shift` — subroutine
- `norm_noise_fft_clip_shift` — subroutine
- `norm_noise_fft_clip_shift_ctf_flip` — subroutine
- `norm_noise_fft_clip_shift_ctf_flip` — subroutine
- `norm_noise_mask_fft_powspec` — subroutine
- `norm_noise_mask_fft_powspec` — subroutine
- `norm_noise_taper_edge_pad_fft` — subroutine
- `norm_noise_taper_edge_pad_fft` — subroutine
- `norm_within` — subroutine
- `norm_within` — subroutine
- `nu_objective` — subroutine
- `nu_objective` — subroutine
- `nu_objective_noise_scale` — function
- `nu_objective_noise_scale` — function
- `one_at_edge` — subroutine
- `one_at_edge` — subroutine
- `open` — subroutine
- `open` — subroutine
- `oshift_int` — function
- `oshift_int` — function
- `oshift_real` — function
- `oshift_real` — function
- `pad` — subroutine
- `pad` — subroutine
- `pad_fft` — subroutine
- `pad_fft` — subroutine
- `pad_inplace` — subroutine
- `pad_inplace` — subroutine
- `pad_mats` — subroutine
- `pad_mats` — subroutine
- `pad_mirr_inplace` — subroutine
- `pad_mirr_inplace` — subroutine
- `pad_mirr_into` — subroutine
- `pad_mirr_into` — subroutine
- `phase_corr` — subroutine
- `phase_corr` — subroutine
- `phase_rand` — subroutine
- `phase_rand` — subroutine
- `polar_memo_valid` — function
- `polarize` — subroutine
- `polarize` — subroutine
- `polarize_oversamp` — subroutine
- `polarize_oversamp` — subroutine
- `power_spectrum` — subroutine
- `power_spectrum` — subroutine
- `prenorm4real_corr_inplace` — subroutine
- `prenorm4real_corr_inplace` — subroutine
- `prenorm4real_corr_sxx` — subroutine
- `prenorm4real_corr_sxx` — subroutine
- `prenorm4real_corr_sxx_masked` — subroutine
- `prenorm4real_corr_sxx_masked` — subroutine
- `presence` — function
- `presence` — function
- `print_cmat` — subroutine
- `print_cmat` — subroutine
- `print_rmat` — subroutine
- `print_rmat` — subroutine
- `quantize_bwd` — subroutine
- `quantize_bwd` — subroutine
- `quantize_fwd` — subroutine
- `quantize_fwd` — subroutine
- `radial_cc` — subroutine
- `radial_cc` — subroutine
- `ran` — subroutine
- `ran` — subroutine
- `read` — subroutine
- `read` — subroutine
- `read_and_crop` — subroutine
- `read_and_crop` — subroutine
- `read_raw_mrc` — subroutine
- `read_raw_mrc` — subroutine
- `read_single_mrc_image` — subroutine
- `read_single_mrc_image` — subroutine
- `real_corr_masked` — function
- `real_corr_masked` — function
- `real_corr_prenorm_inplace` — function
- `real_corr_prenorm_inplace` — function
- `real_corr_prenorm_sxx` — function
- `real_corr_prenorm_sxx` — function
- `real_corr_prenorm_sxx_masked` — function
- `real_corr_prenorm_sxx_masked` — function
- `real_corr_unmasked` — function
- `real_corr_unmasked` — function
- `real_space_filter` — subroutine
- `real_space_filter` — subroutine
- `remove_edge` — subroutine
- `remove_edge` — subroutine
- `remove_neg` — subroutine
- `remove_neg` — subroutine
- `reset_mats` — subroutine
- `reset_mats` — subroutine
- `resmsk` — subroutine
- `resmsk` — subroutine
- `ring` — subroutine
- `ring` — subroutine
- `rmsd` — subroutine
- `rmsd` — subroutine
- `roavg` — subroutine
- `roavg` — subroutine
- `rtsq` — subroutine
- `rtsq` — subroutine
- `rtsq_serial` — subroutine
- `rtsq_serial` — subroutine
- `salt_n_pepper` — subroutine
- `salt_n_pepper` — subroutine
- `same_dims` — function
- `same_dims` — function
- `same_dims_img` — function
- `same_dims_img` — function
- `same_smpd` — function
- `same_smpd` — function
- `scale_mats` — subroutine
- `scale_mats` — subroutine
- `scale_pixels` — subroutine
- `scale_pixels` — subroutine
- `scale_pspec4viz` — subroutine
- `scale_pspec4viz` — subroutine
- `serialize_masked` — function
- `serialize_masked` — function
- `serialize_thres` — function
- `serialize_thres` — function
- `set_cmat` — subroutine
- `set_cmat` — subroutine
- `set_cmat_2d` — subroutine
- `set_cmat_2d` — subroutine
- `set_cmat_3d` — subroutine
- `set_cmat_3d` — subroutine
- `set_cmat_at` — subroutine
- `set_cmat_at` — subroutine
- `set_cmat_img` — subroutine
- `set_cmat_img` — subroutine
- `set_cmat_scalar` — subroutine
- `set_cmat_scalar` — subroutine
- `set_fcomp` — subroutine
- `set_fcomp` — subroutine
- `set_ft` — subroutine
- `set_ft` — subroutine
- `set_img` — subroutine
- `set_img` — subroutine
- `set_real_at` — subroutine
- `set_real_at` — subroutine
- `set_rmat` — subroutine
- `set_rmat` — subroutine
- `set_rmat_at` — subroutine
- `set_rmat_at` — subroutine
- `set_slice` — subroutine
- `set_slice` — subroutine
- `set_smpd` — subroutine
- `set_smpd` — subroutine
- `set_within` — subroutine
- `set_within` — subroutine
- `set_wthreads` — subroutine
- `set_wthreads` — subroutine
- `shift` — subroutine
- `shift` — subroutine
- `shift2D_demoivre` — subroutine
- `shift2D_demoivre` — subroutine
- `shift2Dserial_inplace` — subroutine
- `shift2Dserial_inplace` — subroutine
- `shift2Dserial_into` — subroutine
- `shift2Dserial_into` — subroutine
- `shift_phorig` — subroutine
- `shift_phorig` — subroutine
- `skew` — function
- `skew` — function
- `snr` — function
- `snr` — function
- `soft_ring` — subroutine
- `soft_ring` — subroutine
- `spectrum` — subroutine
- `spectrum` — subroutine
- `sqeuclid` — function
- `sqeuclid` — function
- `square` — subroutine
- `square` — subroutine
- `square_dims` — function
- `square_dims` — function
- `stats_mskimg` — subroutine
- `stats_mskimg` — subroutine
- `stats_region` — subroutine
- `stats_region` — subroutine
- `subtr_backgr` — subroutine
- `subtr_backgr` — subroutine
- `subtr_backgr_ramp` — subroutine
- `subtr_backgr_ramp` — subroutine
- `subtr_cmplx_at` — subroutine
- `subtr_cmplx_at` — subroutine
- `subtr_cmplx_at_kdiv` — subroutine
- `subtr_cmplx_at_kdiv` — subroutine
- `subtr_img` — subroutine
- `subtr_img` — subroutine
- `subtr_scalar` — subroutine
- `subtr_scalar` — subroutine
- `subtract_background` — subroutine
- `subtract_background` — subroutine
- `subtraction` — function
- `subtraction` — function
- `sum_reduce_mats` — subroutine
- `sum_reduce_mats` — subroutine
- `taper_edges` — subroutine
- `taper_edges` — subroutine
- `taper_edges_hann` — subroutine
- `taper_edges_hann` — subroutine
- `taper_edges_particle` — subroutine
- `taper_edges_particle` — subroutine
- `taper_edges_vol` — subroutine
- `taper_edges_vol` — subroutine
- `tile` — subroutine
- `tile` — subroutine
- `tophat` — subroutine
- `tophat` — subroutine
- `update_header_stats` — subroutine
- `update_header_stats` — subroutine
- `upsample_gain` — subroutine
- `upsample_gain` — subroutine
- `upsample_square_background` — subroutine
- `upsample_square_background` — subroutine
- `variance` — function
- `variance` — function
- `vis` — subroutine
- `vis` — subroutine
- `warp_frames_and_sum` — subroutine
- `warp_frames_and_sum` — subroutine
- `weighted_sqsum` — function
- `weighted_sqsum` — function
- `weighted_subtr_corr` — subroutine
- `weighted_subtr_corr` — subroutine
- `weighted_subtr_sqsum` — function
- `weighted_subtr_sqsum` — function
- `whiten_noise_power` — subroutine
- `whiten_noise_power` — subroutine
- `win2arr` — function
- `win2arr` — function
- `win2arr_rad` — subroutine
- `win2arr_rad` — subroutine
- `window` — subroutine
- `window` — subroutine
- `window_slim` — subroutine
- `window_slim` — subroutine
- `write` — subroutine
- `write` — subroutine
- `write_jpg` — subroutine
- `write_jpg` — subroutine
- `zero` — subroutine
- `zero` — subroutine
- `zero2one` — subroutine
- `zero2one` — subroutine
- `zero_and_flag_ft` — subroutine
- `zero_and_flag_ft` — subroutine
- `zero_and_unflag_ft` — subroutine
- `zero_and_unflag_ft` — subroutine
- `zero_background` — subroutine
- `zero_background` — subroutine
- `zero_below` — subroutine
- `zero_below` — subroutine
- `zero_edgeavg` — subroutine
- `zero_edgeavg` — subroutine
- `zero_env_background` — subroutine
- `zero_env_background` — subroutine
- `zero_neg` — subroutine
- `zero_neg` — subroutine
- `zero_padding` — subroutine
- `zero_padding` — subroutine

Local symbols (enclosing procedure scope):
- `distance_transform::dt_1d` — subroutine
- `frame_ref_phase_correlation::interp_peak` — function
- `NLmean2D::reflect_ind` — function
- `nu_objective::squared_loss_scaled` — function
- `warp_frames_and_sum::interp_bilin` — function
- `zero_background::add_edge` — subroutine
- `zero_background::edge_median` — function
- `zero_background::set_face_indices` — subroutine

---
## Module: simple_image_bin

Files:
- `main/image/simple_image_bin.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_neighs`

Public symbols:
- `image_bin` — type

Private symbols:
- `apply_mask` — subroutine
- `border_mask` — subroutine
- `cc2bin` — subroutine
- `copy_bimg` — subroutine
- `cos_edge` — subroutine
- `diameter_cc` — subroutine
- `dilate` — subroutine
- `elim_cc` — subroutine
- `elim_ccs` — subroutine
- `erode` — subroutine
- `feret_minmax` — subroutine
- `find_ccs` — subroutine
- `get_cc_voxels` — subroutine
- `get_imat` — subroutine
- `get_nccs` — subroutine
- `grow_bins` — subroutine
- `inv_bimg` — subroutine
- `kill_bimg` — subroutine
- `masscen_cc` — subroutine
- `max_dist` — subroutine
- `new_bimg` — subroutine
- `order_ccs` — subroutine
- `read_bimg` — subroutine
- `set_edgecc2background` — subroutine
- `set_imat` — subroutine
- `set_largestcc2background` — subroutine
- `size_ccs` — function
- `transfer2bimg` — subroutine
- `update_img_rmat` — subroutine
- `vol_shape_descr` — subroutine
- `write_bimg` — subroutine

Local symbols (enclosing procedure scope):
- `cos_edge::local_versine` — function
- `cos_edge::update_mask_2d` — subroutine
- `cos_edge::update_mask_3d` — subroutine
- `find_ccs::dequeue` — subroutine
- `find_ccs::empty_queue` — subroutine
- `find_ccs::enqueue` — subroutine
- `find_ccs::flood_fill` — subroutine
- `find_ccs::is_empty` — function

---
## Module: simple_image_msk

Files:
- `main/image/simple_image_msk.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_image_bin`
- `simple_parameters`
- `simple_segmentation`

Public symbols:
- `automask2D_mskdiam` — subroutine
- `automask2D_settings` — type
- `automask2D_settings_from` — function
- `automask2D_support_pix` — subroutine
- `density_inoutside_mask` — subroutine
- `image_msk` — type

Private symbols:
- `apply_3Dmask` — subroutine
- `automask2D_binary_one` — subroutine
- `automask2D_mskrad` — function
- `automask2D_params` — subroutine
- `automask2D_set` — subroutine
- `automask3D` — subroutine
- `automask3D_binarize` — subroutine
- `binary_imat_to_pix` — subroutine
- `env_rproject` — subroutine
- `envmask3D_from_lmask` — subroutine
- `estimate_spher_mask_diam` — subroutine

---
## Module: simple_image_msk_tester

Files:
- `main/image/simple_image_msk_tester.f90`

Uses:
- `simple_defs`
- `simple_image`
- `simple_image_bin`
- `simple_image_msk`
- `simple_math`
- `simple_math_ft`
- `simple_parameters`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_image_bin_tests` — subroutine
- `run_all_mask_tests` — subroutine
- `run_all_nano_mask_tests` — subroutine

Private symbols:
- `cleanup_nano_mask_files` — subroutine
- `symmetric_2d` — function
- `symmetric_3d` — function
- `test_automask2D_mskdiam` — subroutine
- `test_bounds_from_mask3D` — subroutine
- `test_cc_bookkeeping` — subroutine
- `test_cc_voxels` — subroutine
- `test_ccs_connectivity_2D` — subroutine
- `test_ccs_connectivity_3D` — subroutine
- `test_ccs_degenerate` — subroutine
- `test_density_support_envelope` — subroutine
- `test_diameter_cc_3D` — subroutine
- `test_disc_and_cos_edge` — subroutine
- `test_estimate_spher_mask_diam` — subroutine
- `test_graphene_mask` — subroutine
- `test_holes_and_feret` — subroutine
- `test_image_bin_examples` — subroutine
- `test_mask2D_semantics` — subroutine
- `test_mask3D_semantics` — subroutine
- `test_masks_parallel_equals_serial` — subroutine
- `test_morphology` — subroutine

Local symbols (enclosing procedure scope):
- `test_bounds_from_mask3D::brute_force_bounds` — subroutine
- `test_density_support_envelope::blobs` — function
- `test_graphene_mask::check_bands` — subroutine
- `test_masks_parallel_equals_serial::check_2d` — subroutine
- `test_masks_parallel_equals_serial::check_3d` — subroutine
- `test_masks_parallel_equals_serial::fill` — subroutine
- `test_masks_parallel_equals_serial::mask_2d` — subroutine
- `test_masks_parallel_equals_serial::mask_3d` — subroutine

---
## Module: simple_image_serialize_tester

Files:
- `main/image/simple_image_serialize_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_test_utils`

Public symbols:
- `run_all_image_serialize_tests` — subroutine

Private symbols:
- `make_distinct_image` — subroutine
- `test_masked_serialize` — subroutine

---
## Module: simple_image_slab

Files:
- `fileio/simple_image_slab.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `image_slab` — type

Private symbols:
- `bytes` — function
- `copy_slot` — subroutine
- `exists` — function
- `get` — subroutine
- `get_many` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `put` — subroutine
- `seal` — subroutine

---
## Module: simple_image_tester

Files:
- `main/image/simple_image_tester.f90`

Uses:
- `ieee_arithmetic`
- `ieee_exceptions`
- `simple_cmdline`
- `simple_ctf`
- `simple_defs`
- `simple_forked_process`
- `simple_image`
- `simple_linalg`
- `simple_memoize_ft_maps`
- `simple_ori`
- `simple_projector`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_image_tests` — subroutine

Private symbols:
- `apply_ctf_other_box` — subroutine
- `band_equals` — function
- `ctf_box_mismatch` — type
- `flat_spectrum` — subroutine
- `pad_dims` — function
- `pad_tag` — function
- `pearson` — function
- `smooth_noise` — subroutine
- `test_acf` — subroutine
- `test_apply_filter` — subroutine
- `test_bandpass` — subroutine
- `test_checkups_and_stats` — subroutine
- `test_construct_and_access` — subroutine
- `test_copy_and_inplace_threading` — subroutine
- `test_corr` — subroutine
- `test_ctf_maps_box_check` — subroutine
- `test_dose_weighing_modes` — subroutine
- `test_fft_roundtrip_and_nyquist` — subroutine
- `test_file_roundtrip` — subroutine
- `test_file_roundtrip_sizes` — subroutine
- `test_flip` — subroutine
- `test_fproject` — subroutine
- `test_masscen` — subroutine
- `test_odd_box_shift_constant` — subroutine
- `test_padding_after_arithmetic` — subroutine
- `test_padding_after_binary_ops` — subroutine
- `test_padding_after_entering_real_space` — subroutine
- `test_padding_after_pad_and_collage` — subroutine
- `test_padding_after_transforms` — subroutine
- `test_padding_query` — subroutine
- `test_power_spectrum` — subroutine
- `test_rmat_ptr_is_box` — subroutine
- `test_rotations_and_roavg` — subroutine
- `test_shifts` — subroutine
- `test_win2arr_rad` — subroutine

Local symbols (enclosing procedure scope):
- `test_file_roundtrip::del_all` — subroutine
- `test_file_roundtrip_sizes::pattern` — function
- `test_rotations_and_roavg::in_disc` — function

---
## Module: simple_imgarr_utils

Files:
- `main/image_processing/simple_imgarr_utils.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_sp_project`
- `simple_stack_io`

Public symbols:
- `alloc_imgarr` — subroutine
- `dealloc_imgarr` — subroutine
- `extract_imgarr` — function
- `join_imgarrs` — function
- `pack_imgarr` — function
- `rank_cavgs_stk` — subroutine
- `read_cavgs_into_imgarr` — function
- `read_stk_into_imgarr` — function
- `write_imgarr_1` — subroutine
- `write_imgarr_2` — subroutine
- `write_imgarr_3` — subroutine
- `write_junk_cavgs` — subroutine
- `write_selected_cavgs` — subroutine

---
## Module: simple_imgfile

Files:
- `fileio/simple_imgfile.f90`

Uses:
- `gnufor2`
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_core_module_api`
- `simple_imghead`
- `simple_string_utils`
- `simple_tifflib`

Public symbols:
- `imgfile` — type
- `simple_float16_to_float32` — subroutine

Private symbols:
- `close` — subroutine
- `getDim` — function
- `getDims` — function
- `getIform` — function
- `getMode` — function
- `open` — subroutine
- `open_local` — subroutine
- `prepare_float16_write_buffer` — subroutine
- `real16_array_to_real32` — subroutine
- `real32_to_real16` — function
- `rSlices` — subroutine
- `rTiffSlices` — subroutine
- `setDims` — subroutine
- `setIform` — subroutine
- `setMean` — subroutine
- `setMinMax` — subroutine
- `setMode` — subroutine
- `setPixSz` — subroutine
- `setRMSD` — subroutine
- `simple_float32_to_float16_f16c` — function
- `slice2bytepos` — subroutine
- `slice2recpos` — subroutine
- `update_MRC_stats` — subroutine
- `wmrcSlices` — subroutine
- `write_float16_slices` — subroutine
- `wSlices` — subroutine

---
## Module: simple_imghead

Files:
- `fileio/simple_imghead.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_tifflib`

Public symbols:
- `find_img_smpd` — function
- `find_ldim_nptcls` — subroutine
- `get_mrc_minmax` — subroutine
- `get_mrcfile_info` — subroutine
- `has_ldim_nptcls` — function
- `ImgHead` — type
- `MrcImgHead` — type
- `SpiImgHead` — type
- `TiffImgHead` — type
- `update_stack_nimgs` — subroutine

Private symbols:
- `bytesPerPix` — function
- `CloseTiff` — subroutine
- `decode` — subroutine
- `firstDataByte` — function
- `get_eerfile_info` — subroutine
- `get_spifile_info` — subroutine
- `get_tiffile_info` — subroutine
- `getDim` — function
- `getDims` — function
- `getIform` — function
- `getLabbyt` — function
- `getLabrec` — function
- `getLenbyt` — function
- `getMaxim` — function
- `getMaxPixVal` — function
- `getMean` — function
- `getMinPixVal` — function
- `getMode` — function
- `getPixSz` — function
- `getRMSD` — function
- `getStackSz` — function
- `isTiff` — function
- `kill` — subroutine
- `new` — subroutine
- `pixIsComplex` — function
- `pixIsSigned` — function
- `print_imghead` — subroutine
- `read` — subroutine
- `read_tiff` — subroutine
- `reset2default` — subroutine
- `setDim` — subroutine
- `setDims` — subroutine
- `setIform` — subroutine
- `setMaxim` — subroutine
- `setMaxPixVal` — subroutine
- `setMean` — subroutine
- `setMinimal` — subroutine
- `setMinPixVal` — subroutine
- `setMode` — subroutine
- `setPixSz` — subroutine
- `setRMSD` — subroutine
- `transfer_byte_array2obj` — subroutine
- `transfer_obj2byte_array` — subroutine
- `write` — subroutine

Local symbols (enclosing procedure scope):
- `setMinimal::is_even` — function

---
## Module: simple_imghead_tester

Files:
- `fileio/simple_imghead_tester.f90`

Uses:
- `simple_defs`
- `simple_imghead`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_imghead_tests` — subroutine

Private symbols:
- `test_mrc_roundtrip` — subroutine
- `test_spider_geometry` — subroutine
- `test_spider_probe` — subroutine
- `test_spider_roundtrip` — subroutine

---
## Module: simple_ipc_tcp_socket_client

Files:
- `utils/comm/simple_ipc_tcp_socket_client.f90`

Uses:
- `iso_c_binding`
- `simple_core_module_api`
- `simple_ipc_tcp_socket_helpers`
- `unix`

Public symbols:
- `ipc_tcp_socket_client` — type
- `tcp_sockaddr_in` — type

Private symbols:
- `c_timeval` — type
- `connect` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `send_recv_msg` — subroutine
- `set_retry_backoff_ms` — subroutine

---
## Module: simple_ipc_tcp_socket_helpers

Files:
- `utils/comm/simple_ipc_tcp_socket_helpers.f90`

Uses:
- `iso_c_binding`
- `unix`

Public symbols:
- `accept_connection` — function
- `c_pollfd` — type
- `fd_is_healthy` — function
- `poll_fds` — subroutine

Private symbols:
- `c_poll` — function
- `c_recv` — function

---
## Module: simple_ipc_tcp_socket_server

Files:
- `utils/comm/simple_ipc_tcp_socket_server.f90`

Uses:
- `iso_c_binding`
- `simple_core_module_api`
- `simple_error`
- `simple_ipc_tcp_socket_helpers`
- `unix`

Public symbols:
- `ipc_tcp_socket_server` — type
- `listener_args` — type
- `recv_msg` — subroutine
- `repl_msg` — subroutine
- `tcp_sockaddr_in` — type

Private symbols:
- `find_server_ips` — subroutine
- `forget_inherited` — subroutine
- `get_port` — function
- `get_server_ips` — function
- `is_listening` — function
- `join_listener` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `start_listener` — subroutine

---
## Module: simple_ipc_tcp_socket_tester

Files:
- `utils/comm/simple_ipc_tcp_socket_tester.f90`

Uses:
- `iso_c_binding`
- `simple_core_module_api`
- `simple_ipc_tcp_socket_client`
- `simple_ipc_tcp_socket_helpers`
- `simple_ipc_tcp_socket_server`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_ipc_tcp_socket_tests` — subroutine

Private symbols:
- `assert_port_in_range` — subroutine
- `destroy_listener_args` — subroutine
- `init_listener_args` — subroutine
- `listener_echo_ok_once` — subroutine
- `listener_ready_once` — subroutine
- `poke_server_port` — subroutine
- `test_client_send_recv_without_server_fails` — subroutine
- `test_client_server_roundtrip` — subroutine
- `test_helpers_accept_connection_invalid_fd` — subroutine
- `test_helpers_constants` — subroutine
- `test_helpers_fd_is_healthy_invalid_fd` — subroutine
- `test_helpers_poll_fds_bounds_and_zero` — subroutine
- `test_server_find_ips_and_port_range` — subroutine
- `test_server_new_kill_lifecycle` — subroutine
- `test_server_start_join_lifecycle` — subroutine

---
## Module: simple_is_check_assert

Files:
- `utils/simple_is_check_assert.f90`

Uses:
- `simple_defs`
- `simple_error`

Public symbols:
- `assert_eq2` — function
- `assert_eq3` — function
- `assert_eq4` — function
- `assert_eqn` — function
- `check4nans2D_cmplx` — subroutine
- `check4nans2D_real` — subroutine
- `check4nans3D_cmplx` — subroutine
- `check4nans3D_real` — subroutine
- `check4nans_cmplx` — subroutine
- `check4nans_real` — subroutine
- `is_a_number_cmplx` — function
- `is_a_number_cmplx_dp` — function
- `is_a_number_dp` — function
- `is_a_number_real` — function
- `is_equal_cmplx` — function
- `is_equal_dp` — function
- `is_equal_real` — function
- `is_even_array` — function
- `is_even_scalar` — function
- `is_gt_zero_dp` — function
- `is_gt_zero_int` — function
- `is_gt_zero_real` — function
- `is_odd` — function
- `is_zero_dp` — function
- `is_zero_int` — function
- `is_zero_real` — function

---
## Module: simple_jiffys

Files:
- `utils/simple_jiffys.f90`

Uses:
- `simple_defs`
- `simple_defs_fname`
- `simple_string`
- `simple_syslib`
- `simple_timer`

Public symbols:
- `haloween_end` — subroutine
- `masked_swap_rm` — subroutine
- `masked_swap_rs` — subroutine
- `masked_swap_rv` — subroutine
- `progress` — subroutine
- `progress_gfortran` — subroutine
- `simple_end` — subroutine
- `simple_print_git_version` — subroutine
- `simple_print_timer` — subroutine
- `swap_c` — subroutine
- `swap_cm` — subroutine
- `swap_cv` — subroutine
- `swap_i` — subroutine
- `swap_r` — subroutine
- `swap_rv` — subroutine

---
## Module: simple_jpg

Files:
- `extlibs/jpg/simple_jpg.f90`

Uses:
- `iso_c_binding`
- `simple_core_module_api`

Public symbols:
- `jpg_img` — type
- `write_rgb_jpeg` — subroutine

Private symbols:
- `save_jpeg_r4` — function
- `save_jpeg_r4_3D` — function
- `stbi_write_jpg` — function

---
## Module: simple_jpg_tester

Files:
- `main/image/simple_jpg_tester.f90`

Uses:
- `iso_c_binding`
- `simple_image`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_mrc2jpeg_tests` — subroutine

Private symbols:
- `cleanup_files` — subroutine
- `jpeg_name` — function
- `section_name` — function
- `source_value` — function
- `stbi_image_free` — subroutine
- `stbi_load` — function
- `verify_jpeg` — subroutine
- `write_source_stack` — subroutine

---
## Module: simple_kbinterpol

Files:
- `main/interp/simple_kbinterpol.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_edges_sqwins`

Public symbols:
- `apod_device` — function
- `apod_fast_device` — function
- `apod_kb15_a2` — function
- `kb_windim` — function
- `kbinterpol` — type

Private symbols:
- `apod` — function
- `apod_fast` — function
- `apod_fast_value_deriv` — subroutine
- `apod_mat_2d` — subroutine
- `apod_mat_2d_fast` — subroutine
- `apod_mat_3d` — subroutine
- `apod_mat_3d_fast` — subroutine
- `apod_mat_3d_fast_grad` — subroutine
- `bessi0` — function
- `bessi0_dp` — function
- `bessi0f` — function
- `bessi1` — function
- `constructor` — function
- `copy` — subroutine
- `get_alpha` — function
- `get_wdim` — function
- `get_winsz` — function
- `instr` — function
- `new` — subroutine
- `sinhc` — function

---
## Module: simple_kbinterpol_tester

Files:
- `main/interp/simple_kbinterpol_tester.f90`

Uses:
- `simple_defs`
- `simple_kbinterpol`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_kbinterpol_tests` — subroutine

Private symbols:
- `apod_ref` — function
- `bessel_i0` — function
- `reference_2d` — subroutine
- `reference_3d` — subroutine
- `standard_kernel` — function
- `test_apod_closed_form` — subroutine
- `test_apod_fast_polynomial` — subroutine
- `test_apod_fast_value_deriv` — subroutine
- `test_apod_mat_2d` — subroutine
- `test_apod_mat_3d` — subroutine
- `test_apod_mat_3d_fast_grad` — subroutine
- `test_device_forms` — subroutine
- `test_fast_stencils` — subroutine
- `test_instr` — subroutine
- `test_window_geometry` — subroutine

---
## Module: simple_kcenter

Files:
- `utils/clustering/simple_kcenter.f90`

Uses:
- `simple_clustering_utils`
- `simple_core_module_api`

Public symbols:
- `kcenter` — type

Private symbols:
- `cluster` — subroutine
- `kill` — subroutine
- `new` — subroutine

---
## Module: simple_kcenter_tester

Files:
- `utils/clustering/simple_kcenter_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_kcenter`
- `simple_test_utils`

Public symbols:
- `run_all_kcenter_tests` — subroutine

Private symbols:
- `fixture` — subroutine
- `test_kcenter_coverage` — subroutine
- `test_kcenter_edges` — subroutine

---
## Module: simple_kd_tree

Files:
- `utils/structs/simple_kd_tree.f90`

Uses:
- `ieee_arithmetic`
- `simple_core_module_api`

Public symbols:
- `kd_tree` — type
- `knn_table` — type

Private symbols:
- `build_kd_tree` — subroutine
- `build_node` — subroutine
- `get_kd_tree_height` — function
- `get_kd_tree_nnodes` — function
- `grow_nodes` — subroutine
- `heap_consider` — subroutine
- `heap_sift_down` — subroutine
- `kill_kd_tree` — subroutine
- `kill_knn_table` — subroutine
- `median_of_three_position` — function
- `new_node` — function
- `node_height` — function
- `pair_better` — function
- `pair_worse` — function
- `point_less` — function
- `query_all_kd_tree` — subroutine
- `query_kd_tree` — subroutine
- `query_node` — subroutine
- `reserve_nodes` — subroutine
- `select_kth` — subroutine
- `sort_neighbor_pairs` — subroutine
- `swap_int` — subroutine
- `swap_pair` — subroutine
- `validate_query` — subroutine
- `widest_dimension` — function

---
## Module: simple_kmeans

Files:
- `utils/clustering/simple_kmeans.f90`

Uses:
- `simple_clustering_utils`
- `simple_core_module_api`

Public symbols:
- `kmeans` — type

Private symbols:
- `cluster` — subroutine
- `get_nchanged` — function
- `get_niters` — function
- `kill` — subroutine
- `new` — subroutine

---
## Module: simple_kmeans_tester

Files:
- `utils/clustering/simple_kmeans_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_kmeans`
- `simple_test_utils`

Public symbols:
- `run_all_kmeans_tests` — subroutine

Private symbols:
- `ring` — subroutine
- `test_kmeans_edges` — subroutine
- `test_kmeans_separated_groups` — subroutine
- `test_kmeans_weighted_metric` — subroutine

---
## Module: simple_kmedoids

Files:
- `utils/clustering/simple_kmedoids.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `kmedoids` — type

Private symbols:
- `assign_labels` — subroutine
- `cluster` — subroutine
- `find_medoid` — subroutine
- `find_medoids` — subroutine
- `get_labels` — subroutine
- `get_medoids` — subroutine
- `init` — subroutine
- `kill` — subroutine
- `merge` — subroutine
- `merge_ranked` — subroutine
- `new_by_ncls` — subroutine
- `new_from_labels` — subroutine

---
## Module: simple_linalg

Files:
- `utils/math/simple_linalg.f90`

Uses:
- `simple_defs`
- `simple_is_check_assert`

Public symbols:
- `ang2vox` — function
- `arg` — function
- `chol_backward` — subroutine
- `chol_forward` — subroutine
- `cholesky` — subroutine
- `fit_lsq_plane` — subroutine
- `fit_straight_line` — subroutine
- `gemm_tn` — subroutine
- `plane_from_points` — function
- `projz` — subroutine
- `qr_solve` — subroutine
- `solve_real_spd_complex` — subroutine
- `spd_inverse` — subroutine
- `spd_logdet` — subroutine
- `trace` — function
- `vector_angle_norm` — function
- `vox2ang` — function

Private symbols:
- `deg2rad_dp` — function
- `deg2rad_sp` — function
- `dgelsy` — subroutine
- `dgesvd` — subroutine
- `dgetrf` — subroutine
- `dgetri` — subroutine
- `dnrm2` — function
- `dsyev` — subroutine
- `eigh_sp` — subroutine
- `eigsrt_dp` — subroutine
- `eigsrt_sp` — subroutine
- `euclid_dp` — function
- `euclid_sp_1d` — function
- `euclid_sp_2d` — function
- `hyp_int_xy` — function
- `hyp_int_xyz` — function
- `hyp_real_xy` — function
- `hyp_real_xyz` — function
- `jacobi_dp` — subroutine
- `jacobi_sp` — subroutine
- `lapack_stop` — subroutine
- `matinv_dp` — subroutine
- `matinv_sp` — subroutine
- `myacos_dp` — function
- `myacos_sp` — function
- `norm_2_dp` — function
- `norm_2_sp` — function
- `pythag_dp` — function
- `pythag_sp` — function
- `rad2deg_dp` — function
- `rad2deg_sp` — function
- `sgemm` — subroutine
- `sgesvd` — subroutine
- `sgetrf` — subroutine
- `sgetri` — subroutine
- `ssyev` — subroutine
- `ssyevr` — subroutine
- `svbksb_dp` — subroutine
- `svbksb_sp` — subroutine
- `svd_multifit_dp` — subroutine
- `svd_multifit_sp` — subroutine
- `svdcmp_dp` — subroutine
- `svdcmp_sp` — subroutine
- `svdfit_dp` — subroutine
- `svdfit_sp` — subroutine
- `vabs_dp` — function
- `vabs_sp` — function

Local symbols (enclosing procedure scope):
- `svd_multifit_dp::funcs` — function
- `svd_multifit_sp::funcs` — function
- `svdfit_dp::funcs` — function
- `svdfit_sp::funcs` — function

---
## Module: simple_linalg_tester

Files:
- `utils/math/simple_linalg_tester.f90`

Uses:
- `simple_defs`
- `simple_linalg`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_linalg_tests` — subroutine

Private symbols:
- `affine2` — function
- `assert_eigenpairs` — subroutine
- `example_inverse` — function
- `example_matrix` — function
- `identity` — function
- `poly3` — function
- `sort_desc` — subroutine
- `test_eigh_largest` — subroutine
- `test_eigh_smallest` — subroutine
- `test_fit_straight_line` — subroutine
- `test_fit_straight_line_recovery` — subroutine
- `test_gemm_tn` — subroutine
- `test_jacobi_eigsrt` — subroutine
- `test_matinv` — subroutine
- `test_plane_fits` — subroutine
- `test_real_spd_complex_solve` — subroutine
- `test_spd_cholesky` — subroutine
- `test_svd_multifit_plane` — subroutine
- `test_svdcmp` — subroutine
- `test_svdfit_polynomial` — subroutine
- `test_vector_helpers` — subroutine

---
## Module: simple_linked_list

Files:
- `utils/structs/simple_linked_list.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `linked_list` — type
- `list_iterator` — type

Private symbols:
- `advance` — subroutine
- `append` — function
- `assign` — subroutine
- `at` — subroutine
- `at_char` — subroutine
- `at_complex` — subroutine
- `at_int` — subroutine
- `at_logical` — subroutine
- `at_real` — subroutine
- `back` — subroutine
- `begin` — function
- `end` — function
- `end_iter` — function
- `equals` — function
- `front` — subroutine
- `getter` — subroutine
- `has_value` — function
- `index` — function
- `is_empty` — function
- `kill` — subroutine
- `kill_node` — subroutine
- `next` — subroutine
- `node` — type
- `pop_front` — subroutine
- `push_back` — subroutine
- `push_back_char` — subroutine
- `push_back_complex` — subroutine
- `push_back_int` — subroutine
- `push_back_logical` — subroutine
- `push_back_real` — subroutine
- `push_front` — subroutine
- `replace_at` — subroutine
- `replace_iterator` — subroutine
- `replace_with` — subroutine
- `size` — function
- `slice` — subroutine

---
## Module: simple_linked_list_tester

Files:
- `utils/structs/simple_linked_list_tester.f90`

Uses:
- `simple_linked_list`
- `simple_test_utils`

Public symbols:
- `run_all_list_tests` — subroutine

Private symbols:
- `particle` — type
- `test_append` — subroutine
- `test_assign_and_copy_semantics` — subroutine
- `test_complex` — subroutine
- `test_derived_type_independence` — subroutine
- `test_derived_type_values` — subroutine
- `test_end_iterator_behavior` — subroutine
- `test_front_back_at` — subroutine
- `test_intrinsic_integer` — subroutine
- `test_iteration` — subroutine
- `test_iterator_index_and_advance` — subroutine
- `test_kill` — subroutine
- `test_logical` — subroutine
- `test_move_semantics` — subroutine
- `test_nested_allocations` — subroutine
- `test_pop_front` — subroutine
- `test_push_and_size` — subroutine
- `test_real` — subroutine
- `test_replace_iterator_stability` — subroutine
- `test_slice` — subroutine
- `test_string_values` — subroutine
- `test_type_safe_wrappers` — subroutine

---
## Module: simple_lpstages_tester

Files:
- `utils/filter/simple_lpstages_tester.f90`

Uses:
- `simple_butterworth`
- `simple_defs`
- `simple_estimate_ssnr`
- `simple_image`
- `simple_string_utils`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_lpstages_tests` — subroutine

Private symbols:
- `assert_crop_consistent` — subroutine
- `frc_curve` — function
- `test_bfac_cap_filter` — subroutine
- `test_butterworth` — subroutine
- `test_filter_nyquist_support` — subroutine
- `test_fsc_weightings` — subroutine
- `test_lpstages_aggressive` — subroutine
- `test_lpstages_aggressive_from_stage` — subroutine
- `test_lpstages_cavgs_vs_ptcls` — subroutine
- `test_lpstages_fast` — subroutine
- `test_lpstages_flat_frc_fallback` — subroutine
- `test_lpstages_from_frc` — subroutine
- `test_lpstages_setlims` — subroutine
- `test_lpstages_single_stage` — subroutine
- `test_mskdiam2lplimits` — subroutine
- `test_trslim_override` — subroutine

Local symbols (enclosing procedure scope):
- `test_lpstages_aggressive_from_stage::set_prefix` — subroutine
- `test_lpstages_setlims::count_setlims_diagnostics` — subroutine

---
## Module: simple_magic_boxes

Files:
- `utils/simple_magic_boxes.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_math`

Public symbols:
- `autoscale` — subroutine
- `find_larger_magic_box` — function
- `find_magic_box` — function
- `find_magic_boxes4scale` — function
- `print_magic_box_range` — subroutine

Private symbols:
- `magic_pftsz_crop` — function
- `magic_pftsz_msk` — function

---
## Module: simple_make_cavgs_strategy

Files:
- `main/strategies/parallelization/simple_make_cavgs_strategy.f90`

Uses:
- `simple_builder`
- `simple_classaverager`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_eo_partition`
- `simple_exec_helpers`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_sigma2_bootstrap`
- `simple_sp_project`

Public symbols:
- `create_make_cavgs_strategy` — function
- `make_cavgs_hooks` — type
- `make_cavgs_master_strategy` — type
- `make_cavgs_shmem_strategy` — type
- `make_cavgs_strategy` — type
- `make_cavgs_worker_strategy` — type

Private symbols:
- `apply_defaults_interface` — subroutine
- `apply_distr_entry_defaults` — subroutine
- `base_after_end` — subroutine
- `cavgassemble_hook` — subroutine
- `cavger_prepare_and_assemble` — subroutine
- `cleanup_interface` — subroutine
- `endmsg_interface` — function
- `ensure_canonical_sigma2` — subroutine
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `generate_cluster_centers` — subroutine
- `init_interface` — subroutine
- `master_apply_defaults` — subroutine
- `master_cleanup` — subroutine
- `master_end_message` — function
- `master_execute` — subroutine
- `master_finalize_run` — subroutine
- `master_initialize` — subroutine
- `print_raw_source_verification` — subroutine
- `setup_evenodd` — subroutine
- `shmem_apply_defaults` — subroutine
- `shmem_bookkeeping` — subroutine
- `shmem_cleanup` — subroutine
- `shmem_end_message` — function
- `shmem_execute` — subroutine
- `shmem_finalize_run` — subroutine
- `shmem_initialize` — subroutine
- `worker_apply_defaults` — subroutine
- `worker_cleanup` — subroutine
- `worker_end_message` — function
- `worker_execute` — subroutine
- `worker_finalize_run` — subroutine
- `worker_initialize` — subroutine
- `write_oritype_segment_for_compute` — subroutine

---
## Module: simple_map_reduce

Files:
- `utils/simple_map_reduce.f90`

Uses:
- `simple_defs`
- `simple_fileio`
- `simple_jiffys`
- `simple_srch_sort_loc`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `split_nobjs_active` — function
- `split_nobjs_even` — function
- `split_pairs_in_parts` — subroutine

---
## Module: simple_matcher_2Dprep

Files:
- `main/strategies/search/simple_matcher_2Dprep.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_butterworth`
- `simple_eul_prob_tab_utils`
- `simple_matcher_ptcl_io`
- `simple_matcher_smpl_and_lplims`
- `simple_memoize_ft_maps`
- `simple_pftc_srch_api`
- `simple_projector`
- `simple_strategy2d_utils`
- `simple_timer`

Public symbols:
- `calc_2Dref_offset` — subroutine
- `prep2Dref` — subroutine
- `prepimg4align` — subroutine
- `prepimg4align_cart` — subroutine
- `prepimg4align_cropped` — subroutine
- `ptcl_on_disk` — function

---
## Module: simple_matcher_3Drec

Files:
- `main/strategies/search/simple_matcher_3Drec.f90`

Uses:
- `simple_bench_io`
- `simple_builder`
- `simple_core_module_api`
- `simple_euclid_sigma2`
- `simple_eul_prob_tab_utils`
- `simple_fplane_grid`
- `simple_image`
- `simple_imgarr_utils`
- `simple_matcher_ptcl_io`
- `simple_math_ft`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_reconstructor`
- `simple_refine3d_fnames`
- `simple_state_weight_set`
- `simple_timer`

Public symbols:
- `calc_3Drec` — subroutine
- `cleanup_rec_buffers` — subroutine
- `gen_rec_plane` — subroutine
- `init_rec` — subroutine
- `prep_imgs4rec` — subroutine
- `write_rec3D_bench` — subroutine
- `write_state_half_partial` — subroutine

Private symbols:
- `group_pinds_by_state_eo` — subroutine
- `group_pinds_by_weights` — subroutine
- `init_state_half_rec` — subroutine
- `insert_extra_poses` — subroutine
- `kill_state_half_rec` — subroutine
- `mark_empty_state` — subroutine
- `state_eo_group` — function
- `update_state_half_rec` — subroutine

Local symbols (enclosing procedure scope):
- `group_pinds_by_weights::state_members` — type

---
## Module: simple_matcher_pftc_prep

Files:
- `main/strategies/search/simple_matcher_pftc_prep.f90`

Uses:
- `simple_builder`
- `simple_classaverager`
- `simple_core_module_api`
- `simple_image`
- `simple_matcher_2dprep`
- `simple_matcher_ptcl_batch`
- `simple_matcher_ptcl_io`
- `simple_parameters`
- `simple_strategy2d_alloc`

Public symbols:
- `prep_pftc4align2D` — subroutine

---
## Module: simple_matcher_ptcl_batch

Files:
- `main/strategies/search/simple_matcher_ptcl_batch.f90`

Uses:
- `iso_fortran_env`
- `simple_builder`
- `simple_cartft_calc`
- `simple_euclid_sigma2`
- `simple_imgarr_utils`
- `simple_matcher_2dprep`
- `simple_matcher_ptcl_io`
- `simple_pftc_srch_api`
- `simple_sigma2_state_file`

Public symbols:
- `alloc_ptcl_imgs` — subroutine
- `build_batch_particles2D` — subroutine
- `build_batch_particles3D` — subroutine
- `check_pass_team` — subroutine
- `clean_batch_particles2D` — subroutine
- `clean_batch_particles3D` — subroutine
- `crop_particles2D` — function
- `prep_cart_batch` — subroutine
- `prep_sigmas_objfun` — subroutine
- `write_sigma2_bench_context` — subroutine
- `write_sigma2_bench_timings` — subroutine

Private symbols:
- `polarize_batch_particles3D` — subroutine
- `prep_cart_batch_of_build` — subroutine

---
## Module: simple_matcher_ptcl_io

Files:
- `main/strategies/search/simple_matcher_ptcl_io.f90`

Uses:
- `iso_fortran_env`
- `simple_builder`
- `simple_cartesian_fourier`
- `simple_euclid_sigma2`
- `simple_image`
- `simple_image_slab`
- `simple_pftc_srch_api`
- `simple_stack_read_session`
- `simple_syslib`

Public symbols:
- `fill_pass_buffer` — function
- `get_rec_observation` — subroutine
- `get_records` — subroutine
- `killimgbatch` — subroutine
- `open_record_buffer` — function
- `pass_buffer_bytes` — function
- `pass_buffer_status_name` — function
- `pass_team` — function
- `prep_rec_observation` — subroutine
- `prep_rec_observation_ptcl` — subroutine
- `prepimgbatch` — subroutine
- `ptcl_buffer_budget` — function
- `ptcl_read_counts` — type
- `ptcl_read_counts_now` — function
- `ptcl_read_counts_since` — function
- `ptcl_read_session` — type
- `put_record` — subroutine
- `read_rec_obs_batch` — subroutine
- `rec_obs_inv_sig2` — subroutine
- `rec_obs_oris` — subroutine
- `rec_state_half` — subroutine
- `record_slot` — function
- `records_open` — function
- `release_pass_buffer` — subroutine
- `seal_record_buffer` — subroutine
- `take_records` — function
- `write_mem_bench_context` — subroutine
- `write_ptcl_read_bench_context` — subroutine

Private symbols:
- `buffer_covers` — function
- `expect_ptcl_pass` — subroutine
- `kill_ptcl_read_session` — subroutine
- `new_ptcl_read_session` — subroutine
- `read_imgbatch_1` — subroutine
- `read_imgbatch_2` — subroutine
- `read_imgbatch_3` — subroutine
- `read_ptcl_batch` — subroutine

Local symbols (enclosing procedure scope):
- `read_ptcl_batch::read_particles` — subroutine

---
## Module: simple_matcher_refvol_utils

Files:
- `main/strategies/search/simple_matcher_refvol_utils.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_image`
- `simple_opt_filter`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_refine3d_fnames`
- `simple_refine3d_provenance`

Public symbols:
- `adopt_reprojection_model_range` — subroutine
- `any_volume_source_defined` — function
- `complete_volume_source_available` — function
- `complete_volume_source_defined` — function
- `estimate_lp_from_refs` — subroutine
- `materialize_reprojection_model_from_volumes` — subroutine
- `pick_lp_est_state` — subroutine
- `read_cart_refvols` — subroutine
- `read_mask_filter_reproject_refvols` — subroutine
- `read_reprojection_model` — subroutine
- `remove_ref_section_files` — subroutine
- `reprojection_model_available` — function

Private symbols:
- `calcrefvolshift_and_mapshifts2ptcls` — subroutine
- `cart_refs_needed` — function
- `read_mask_filter_refvols` — subroutine
- `read_reprojection_model_file_header` — subroutine
- `read_reprojection_model_header` — subroutine
- `reprojection_model_diag_available` — function
- `reprojection_model_header_compatible` — function

Local symbols (enclosing procedure scope):
- `read_mask_filter_refvols::apply_automatic_envelope_to_refs` — subroutine
- `read_mask_filter_refvols::mask_matching_reference` — subroutine
- `read_mask_filter_refvols::regularize_ref_with_noise` — subroutine
- `read_mask_filter_reproject_refvols::blend_lowres_eo_for_registration` — subroutine

---
## Module: simple_matcher_smpl_and_lplims

Files:
- `main/strategies/search/simple_matcher_smpl_and_lplims.f90`

Uses:
- `iso_fortran_env`
- `simple_builder`
- `simple_eul_prob_tab_utils`
- `simple_pftc_srch_api`
- `simple_refine3d_fnames`
- `simple_rnd`
- `simple_solve2d_controller`

Public symbols:
- `adopt_run_seed` — subroutine
- `all_active_ptcls_2D_assigned` — function
- `master_draws_sample2D` — function
- `master_draws_sample3D` — function
- `refine2D_blends_carryover` — function
- `refine2D_requires_full_assignment` — function
- `sample_ptcls4fillin` — subroutine
- `sample_ptcls4iter3D` — subroutine
- `sample_ptcls4missing3D` — subroutine
- `sample_ptcls4prob2D` — subroutine
- `sample_ptcls4prob3D` — subroutine
- `sample_ptcls4update2D` — subroutine
- `sample_ptcls4update3D` — subroutine
- `set_bp_range2D` — subroutine
- `set_bp_range3D` — subroutine

Private symbols:
- `sample_rnd_key` — function
- `try_nu_project_working_lp` — function

---
## Module: simple_math

Files:
- `utils/math/simple_math.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_is_check_assert`
- `simple_linalg`
- `simple_srch_sort_loc`

Public symbols:
- `bounds_from_mask3D` — subroutine
- `calc_score_thres` — function
- `ceil_div` — function
- `create_hist_vector` — subroutine
- `cross` — function
- `elim_dup` — subroutine
- `equispaced_vals` — function
- `floor_div` — function
- `gaussian1D` — function
- `gaussian2D` — function
- `gaussian3D` — function
- `gauwfun` — function
- `gcd` — function
- `get_pixel_pos` — subroutine
- `great_circle_samples` — function
- `hac_1d` — subroutine
- `hac_1d_fast` — subroutine
- `hac_med_thres` — subroutine
- `lcm` — function
- `mode` — subroutine
- `nextPow2` — function
- `otsu_binarize` — subroutine
- `otsu_mask` — subroutine
- `otsu_thresh` — subroutine
- `pixels_dist_int` — function
- `pixels_dist_real` — function
- `put_last` — subroutine
- `quadri` — function
- `quantize_vec` — subroutine
- `quantize_vec_serial` — subroutine
- `rotate_vec` — function
- `rotmat2d` — subroutine
- `round2even` — function
- `round2odd` — function
- `SavitzkyGolay_filter` — subroutine
- `shft` — subroutine
- `sinc` — function
- `sortmeans` — subroutine
- `vis_2Dinteger_mat` — subroutine
- `vis_2Dreal_mat` — subroutine
- `vis_3Dinteger_mat` — subroutine
- `vis_3Dreal_mat` — subroutine

Local symbols (enclosing procedure scope):
- `hac_med_thres::find_medoids` — subroutine
- `SavitzkyGolay_filter::savgol` — subroutine

---
## Module: simple_math_ctf

Files:
- `utils/math/simple_math_ctf.f90`

Uses:
- `simple_defs`
- `simple_memoize_ft_maps`

Public symbols:
- `electron_interaction_constant` — function
- `electron_scale_link` — function
- `electron_wavelength` — function
- `ft_map_ctf_kernel` — function

---
## Module: simple_math_ft

Files:
- `utils/math/simple_math_ft.f90`

Uses:
- `ieee_arithmetic`
- `iso_c_binding`
- `simple_defs`
- `simple_error`
- `simple_is_check_assert`
- `simple_srch_sort_loc`
- `simple_type_defs`

Public symbols:
- `calc_fourier_index` — function
- `calc_graphene_mask` — function
- `calc_lowpass_lim` — function
- `csq_dp` — function
- `csq_fast_dp` — function
- `csq_fast_sp` — function
- `csq_sp` — function
- `cyci_1d` — function
- `cyci_1d_static` — function
- `fdim` — function
- `fplane_get_cmplx` — function
- `fplane_get_ctfsq` — function
- `fplane_shift` — subroutine
- `get_find_at_crit` — function
- `get_find_at_res` — function
- `get_resarr` — function
- `mycabs` — function
- `phase_angle` — function
- `resample_sigma2` — subroutine
- `resang` — function

---
## Module: simple_mem_estimator

Files:
- `utils/qsys/simple_mem_estimator.f90`

Uses:
- `simple_core_module_api`
- `simple_parameters`

Public symbols:
- `estimate_mem_usage` — subroutine
- `estimate_mem_usage_2D` — function
- `estimate_mem_usage_ctf_estimation` — function
- `estimate_mem_usage_extract` — function
- `estimate_mem_usage_motion_corr` — function
- `estimate_mem_usage_pick` — function
- `estimate_mem_usage_preprocess` — function
- `estimate_mem_usage_pspec` — function
- `estimate_mem_usage_scale` — function

---
## Module: simple_memoize_ft_maps

Files:
- `main/image/simple_memoize_ft_maps.f90`

Uses:
- `simple_core_module_api`
- `simple_ftiter`

Public symbols:
- `forget_ft_maps` — subroutine
- `ft_map_get_farray_shape` — function
- `ft_map_get_lfny` — function
- `ft_maps_match` — function
- `memoize_ft_maps` — subroutine

---
## Module: simple_memory_monitor

Files:
- `utils/simple_memory_monitor.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_defs`
- `simple_string`
- `simple_syslib`

Public symbols:
- `mem_monitor_finish` — subroutine
- `mem_monitor_init` — subroutine
- `mem_monitor_is_enabled` — function
- `mem_report` — subroutine

Private symbols:
- `log_sample` — subroutine
- `start_monitor` — subroutine

---
## Module: simple_mic_import

Files:
- `main/stream/shared/simple_mic_import.f90`

Uses:
- `simple_oris`
- `simple_sp_project`
- `simple_string`

Public symbols:
- `append_mics_from_projects` — subroutine

Local symbols (enclosing procedure scope):
- `append_mics_from_projects::selected` — function

---
## Module: simple_mic_import_tester

Files:
- `main/stream/shared/simple_mic_import_tester.f90`

Uses:
- `simple_fileio`
- `simple_mic_import`
- `simple_oris`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_mic_import_tests` — subroutine

Private symbols:
- `delete_fixtures` — subroutine
- `test_append_accepted_only` — subroutine
- `test_append_all_in_file_order` — subroutine
- `test_append_to_existing_and_empty_list` — subroutine
- `write_fixtures` — subroutine

---
## Module: simple_mic_selection

Files:
- `main/stream/shared/simple_mic_selection.f90`

Uses:
- `simple_fileio`
- `simple_oris`
- `simple_string`

Public symbols:
- `reject_mics_by_thresholds` — subroutine
- `reject_mics_without_particles` — subroutine

Local symbols (enclosing procedure scope):
- `reject_mics_by_thresholds::exceeds` — function

---
## Module: simple_mic_selection_tester

Files:
- `main/stream/shared/simple_mic_selection_tester.f90`

Uses:
- `simple_fileio`
- `simple_mic_selection`
- `simple_oris`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_mic_selection_tests` — subroutine

Private symbols:
- `make_mics` — subroutine
- `test_absent_thresholds_and_missing_keys` — subroutine
- `test_each_threshold_rejects` — subroutine
- `test_reject_mics_without_particles` — subroutine
- `test_rejected_mics_stay_rejected_and_uncounted` — subroutine

---
## Module: simple_mic_subset

Files:
- `main/project/simple_mic_subset.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_eo_partition`
- `simple_ori`
- `simple_oris`
- `simple_ptcl_subset`
- `simple_sp_project`
- `simple_srch_sort_loc`
- `simple_sym`

Public symbols:
- `apply_mic_subset` — subroutine
- `build_mic_support` — subroutine
- `rank_mic_support` — subroutine
- `select_mic_support` — subroutine

Private symbols:
- `mic_subset_objective` — function

Local symbols (enclosing procedure scope):
- `select_mic_support::preferred` — function
- `select_mic_support::search` — subroutine

---
## Module: simple_mic_subset_tester

Files:
- `main/project/simple_mic_subset_tester.f90`

Uses:
- `simple_defs`
- `simple_eo_partition`
- `simple_mic_subset`
- `simple_oris`
- `simple_ptcl_subset`
- `simple_sp_project`
- `simple_sym`
- `simple_test_utils`

Public symbols:
- `run_all_mic_subset_tests` — subroutine

Private symbols:
- `fixture` — subroutine
- `test_counts_and_refusals` — subroutine
- `test_exhaustive_oracle` — subroutine
- `test_global_search` — subroutine
- `test_project_mask` — subroutine
- `test_shared_coverage` — subroutine

---
## Module: simple_micproc

Files:
- `main/image_processing/simple_micproc.f90`

Uses:
- `simple_bspline_smoother`
- `simple_core_module_api`
- `simple_histogram`
- `simple_image`
- `simple_image_bin`

Public symbols:
- `binarize_mic_den` — subroutine
- `bs_smoother` — subroutine
- `bs_smoother_biomol` — subroutine
- `cascade_filter` — subroutine
- `cascade_filter_biomol` — subroutine
- `flag_amorphous_carbon` — subroutine
- `flag_ice` — subroutine
- `read_mic` — subroutine
- `read_mic_subtr_backgr` — subroutine
- `read_mic_subtr_backgr_shrink` — subroutine
- `sample_filetab` — function

Local symbols (enclosing procedure scope):
- `flag_amorphous_carbon::cluster` — subroutine
- `flag_amorphous_carbon::clustering_rejection` — subroutine

---
## Module: simple_mini_stream_utils

Files:
- `main/stream/simple_mini_stream_utils.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_gui_utils`
- `simple_image`
- `simple_image_bin`
- `simple_micproc`
- `simple_nrtxtfile`
- `simple_parameters`
- `simple_picksegdiam`
- `simple_segdiam_bin_picker`
- `simple_sp_project`

Public symbols:
- `segdiampick_mics` — subroutine
- `segdiampick_preprocess` — subroutine

---
## Module: simple_molecule_data

Files:
- `main/nano/simple_molecule_data.f90`

Uses:
- `simple_defs`
- `simple_error`

Public symbols:
- `betagal_1jyx` — function
- `molecule_data` — type
- `sars_cov2_spkgp_6vxx` — function

Private symbols:
- `is_atom` — function
- `read_molecule_pdb` — function

---
## Module: simple_motion_align_hybrid

Files:
- `main/motion/simple_motion_align_hybrid.f90`

Uses:
- `simple_core_module_api`
- `simple_ft_expanded`
- `simple_ftexp_shsrch`
- `simple_image`
- `simple_parameters`

Public symbols:
- `motion_align_hybrid` — type

Private symbols:
- `align` — subroutine
- `align_corr` — subroutine
- `align_dcorr` — subroutine
- `align_discrete` — subroutine
- `calc_rmsd` — function
- `calc_shifts` — subroutine
- `dealloc_ftexps` — subroutine
- `dealloc_images` — subroutine
- `fit_polynomial` — subroutine
- `gen_weights` — subroutine
- `get_coords` — subroutine
- `get_corr` — function
- `get_corrs` — subroutine
- `get_opt_shifts` — subroutine
- `get_shifts_toplot` — subroutine
- `get_weights` — subroutine
- `init_ftexps` — subroutine
- `init_images` — subroutine
- `is_fitshifts` — function
- `kill` — subroutine
- `new` — subroutine
- `poly` — function
- `polynomial2shift` — subroutine
- `recenter_shifts` — subroutine
- `refine` — subroutine
- `set_bfactor` — subroutine
- `set_coords` — subroutine
- `set_fitshifts` — subroutine
- `set_fixed_frame` — subroutine
- `set_rand_init_shifts` — subroutine
- `set_reslims` — subroutine
- `set_shsrch_tol` — subroutine
- `set_trs` — subroutine
- `set_weights` — subroutine

Local symbols (enclosing procedure scope):
- `calc_shifts::interp_peak` — function
- `polynomial2shift::poly2val` — function

---
## Module: simple_motion_align_nano

Files:
- `main/motion/simple_motion_align_nano.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_parameters`

Public symbols:
- `motion_align_nano` — type

Private symbols:
- `align` — subroutine
- `align_dcorr` — subroutine
- `calc_rmsd` — function
- `calc_shifts` — subroutine
- `dealloc_images` — subroutine
- `gen_weights` — subroutine
- `get_opt_shifts` — subroutine
- `get_reference` — subroutine
- `init_images` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `set_bfactor` — subroutine
- `set_reslims` — subroutine
- `set_trs` — subroutine
- `shift_frames_gen_ref` — subroutine

Local symbols (enclosing procedure scope):
- `calc_shifts::interp_peak` — function

---
## Module: simple_motion_correct

Files:
- `main/motion/simple_motion_correct.f90`

Uses:
- `simple_core_module_api`
- `simple_eer_factory`
- `simple_ft_expanded`
- `simple_image`
- `simple_motion_align_hybrid`
- `simple_motion_correct_utils`
- `simple_motion_model`
- `simple_motion_patched`
- `simple_parameters`

Public symbols:
- `motion_correct_calc_bid` — subroutine
- `motion_correct_iso` — subroutine
- `motion_correct_iso_calc_sums` — subroutine
- `motion_correct_iso_kill` — subroutine
- `motion_correct_iso_shift_frames` — subroutine
- `motion_correct_kill_common` — subroutine
- `motion_correct_mic2spec` — subroutine
- `motion_correct_patched` — subroutine
- `motion_correct_patched_calc_sums` — subroutine
- `motion_correct_patched_kill` — subroutine
- `motion_correct_write_docs` — subroutine

Private symbols:
- `cure_outliers` — subroutine
- `motion_correct_init` — subroutine

---
## Module: simple_motion_correct_iter

Files:
- `main/motion/simple_motion_correct_iter.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_image`
- `simple_motion_correct`
- `simple_parameters`
- `simple_stackops`

Public symbols:
- `motion_correct_iter` — type

Private symbols:
- `calc_npatches` — subroutine
- `get_moviename` — function
- `iterate` — subroutine

---
## Module: simple_motion_correct_strategy

Files:
- `main/strategies/parallelization/simple_motion_correct_strategy.f90`

Uses:
- `simple_binoris_io`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_motion_correct_iter`
- `simple_motion_correct_utils`
- `simple_motion_gain_helpers`
- `simple_parameters`
- `simple_qsys_env`
- `simple_sp_project`

Public symbols:
- `create_motion_correct_strategy` — function
- `motion_correct_distr_strategy` — type
- `motion_correct_inmem_strategy` — type
- `motion_correct_strategy` — type

Private symbols:
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_end_message` — function
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `endmsg_interface` — function
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_end_message` — function
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `set_motion_correct_defaults` — subroutine

---
## Module: simple_motion_correct_utils

Files:
- `main/motion/simple_motion_correct_utils.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_eer_factory`
- `simple_image`

Public symbols:
- `apply_patch_poly` — function
- `calc_eer_fraction` — subroutine
- `correct_gain` — subroutine
- `eer_scale_movie_convention` — subroutine
- `flip_gain` — subroutine
- `patch_poly` — function
- `pix2polycoords` — function

---
## Module: simple_motion_gain_analysis

Files:
- `main/motion/simple_motion_gain_analysis.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `gain_flip_analyzer` — type

Private symbols:
- `analyze_if_due` — subroutine
- `analyzer_kill` — subroutine
- `analyzer_new` — subroutine
- `get_converged` — function
- `get_flip_mode` — function

---
## Module: simple_motion_gain_helpers

Files:
- `main/motion/simple_motion_gain_helpers.f90`

Uses:
- `simple_core_module_api`
- `simple_eer_factory`
- `simple_image`

Public symbols:
- `add_movies_to_gain_sum` — subroutine
- `gainref_to_jpg` — subroutine
- `normalized_inverse_average_intensity` — subroutine
- `read_movies_and_sum_frames` — subroutine
- `write_gain_from_sum` — subroutine

---
## Module: simple_motion_gain_tester

Files:
- `main/motion/simple_motion_gain_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_motion_correct_utils`
- `simple_motion_gain_analysis`
- `simple_motion_gain_helpers`
- `simple_test_utils`

Public symbols:
- `run_all_motion_gain_tests` — subroutine

Private symbols:
- `create_constant_image` — subroutine
- `create_gainref` — subroutine
- `create_movie_stack` — subroutine
- `test_add_movies_to_gain_sum_accumulates` — subroutine
- `test_correct_gain_flipping` — subroutine
- `test_gain_flip_analyzer_batch_updates` — subroutine
- `test_gain_flip_mode_mapping` — subroutine
- `test_normalized_inverse_average_intensity` — subroutine
- `test_read_movies_and_sum_frames_counts` — subroutine
- `test_write_gain_from_sum` — subroutine

---
## Module: simple_motion_model

Files:
- `main/motion/simple_motion_model.f90`

Uses:
- `iso_fortran_env`
- `simple_core_module_api`
- `simple_image`
- `simple_motion_correct_utils`
- `simple_parameters`
- `simple_starfile_wrappers`

Public symbols:
- `motion_model` — type

Private symbols:
- `add_patches` — subroutine
- `determine_patch_coordinates` — subroutine
- `display` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `read` — subroutine
- `read_string` — subroutine
- `refit_polynomial` — subroutine
- `refit_weighted_polynomial` — subroutine
- `set_drift_offsets` — subroutine
- `set_frameweights` — subroutine
- `set_local_offsets` — subroutine
- `set_model_coeffs` — subroutine
- `set_outlier_coords` — subroutine
- `write` — subroutine
- `write_bin` — subroutine
- `write_doc` — subroutine
- `write_string` — subroutine

---
## Module: simple_motion_model_tester

Files:
- `main/motion/simple_motion_model_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_forked_process`
- `simple_image`
- `simple_motion_correct_iter`
- `simple_motion_correct_utils`
- `simple_motion_model`
- `simple_parameters`
- `simple_starfile_wrappers`
- `simple_test_utils`

Public symbols:
- `run_all_motion_model_tests` — subroutine

Private symbols:
- `binary_files_equal` — function
- `correct_grouped_movie` — subroutine
- `grouped_movie` — type
- `test_binary_roundtrip` — subroutine
- `test_binary_roundtrip_with_optional_arrays` — subroutine
- `test_binary_roundtrip_with_rejected_patch` — subroutine
- `test_binary_roundtrip_without_optional_arrays` — subroutine
- `test_flipgain_lifecycle` — subroutine
- `test_movie_geometry` — subroutine
- `test_nframesgrp_outside_tseries` — subroutine
- `test_refit_polynomial_known_coefficients` — subroutine
- `write_reference_binary` — subroutine
- `write_reference_string` — subroutine

---
## Module: simple_motion_patched

Files:
- `main/motion/simple_motion_patched.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_core_module_api`
- `simple_ft_expanded`
- `simple_image`
- `simple_motion_align_hybrid`
- `simple_motion_correct_utils`
- `simple_opt_factory`
- `simple_opt_lbfgsb`
- `simple_opt_spec`
- `simple_optimizer`
- `simple_parameters`

Public symbols:
- `motion_patched` — type

Private symbols:
- `allocate_patches` — subroutine
- `correct` — subroutine
- `deallocate_patches` — subroutine
- `det_shifts` — subroutine
- `fit_polynomial` — subroutine
- `gen_micrograph` — subroutine
- `gen_patch` — subroutine
- `get_local_shift` — subroutine
- `get_poly4star` — subroutine
- `get_poly_coeffs` — subroutine
- `get_polyfit_rmsd` — function
- `kill` — subroutine
- `new` — subroutine
- `plot_shifts` — subroutine
- `set_bfactor` — subroutine
- `set_fitshifts` — subroutine
- `set_fixed_frame` — subroutine
- `set_frameweights` — subroutine
- `set_interp_fixed_frame` — subroutine
- `set_nframes` — subroutine
- `set_poly_coeffs` — subroutine
- `set_size_frames_ref` — subroutine

---
## Module: simple_mrc_validate_tester

Files:
- `main/image/simple_mrc_validate_tester.f90`

Uses:
- `simple_defs`
- `simple_image`
- `simple_imghead`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_mrc_validate_tests` — subroutine

Private symbols:
- `cleanup_files` — subroutine
- `verify_density` — subroutine
- `verify_header` — subroutine
- `voxel_value` — function

---
## Module: simple_nano_density

Files:
- `main/nano/simple_nano_density.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_stat`

Public symbols:
- `fit_atom` — subroutine
- `gauss_filter3D` — subroutine
- `gauss_kernel` — function
- `halo_kernel` — function
- `nano_density` — type
- `regress_profile` — subroutine
- `render_atoms` — subroutine
- `robust_spread` — function
- `transfer_kernels` — type

Private symbols:
- `add_atom` — subroutine
- `add_atoms` — subroutine
- `add_halo` — subroutine
- `add_model` — subroutine
- `aperture_intensities` — subroutine
- `centroid` — function
- `check_residual` — subroutine
- `fit_model` — subroutine
- `fit_one` — subroutine
- `fit_profiles` — subroutine
- `fit_transfer_rounds` — subroutine
- `gather` — subroutine
- `get_atoms` — subroutine
- `get_coordination` — function
- `get_gauge` — subroutine
- `get_intensity_profile` — subroutine
- `get_natoms` — function
- `get_noise` — subroutine
- `get_residual` — subroutine
- `get_rounds` — subroutine
- `get_transfer` — subroutine
- `get_width_profile` — subroutine
- `halo_at` — function
- `has_profiles` — function
- `is_plain` — function
- `keep_atoms` — subroutine
- `kernel_at` — function
- `kernel_values` — function
- `kill` — subroutine
- `kill_kernels` — subroutine
- `measure_transfer` — subroutine
- `new` — subroutine
- `new_kernels` — subroutine
- `noise_from_region` — subroutine
- `noise_from_residual` — subroutine
- `ntab` — function
- `phantom_noise` — subroutine
- `prof_at` — function
- `profile` — subroutine
- `refresh_halo` — subroutine
- `render_fit` — subroutine
- `render_narrow` — subroutine
- `render_sim` — subroutine
- `rerender` — subroutine
- `resid_over_noise` — function
- `residual_correction` — function
- `set_kernels` — subroutine
- `set_transfer` — subroutine
- `set_variances` — subroutine
- `sweep` — subroutine
- `t_of` — function
- `update_background` — subroutine
- `update_coordination` — subroutine
- `width_prior` — subroutine

Local symbols (enclosing procedure scope):
- `fit_atom::cost` — function
- `fit_atom::int_at` — function
- `fit_atom::kernel_on` — subroutine
- `measure_transfer::gauge` — subroutine
- `measure_transfer::med3` — function
- `new_kernels::t_at` — function
- `regress_profile::spread2` — function
- `sweep::cell_colour` — function

---
## Module: simple_nano_density_tester

Files:
- `main/nano/simple_nano_density_tester.f90`

Uses:
- `simple_atoms`
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_nano_density`
- `simple_rnd`
- `simple_test_utils`

Public symbols:
- `rendered_transfer` — subroutine
- `run_all_nano_density_lib_tests` — subroutine
- `run_all_nano_density_tests` — subroutine

Private symbols:
- `check_artefact` — subroutine
- `check_profiles` — subroutine
- `check_transfer` — subroutine
- `fit_measured` — subroutine
- `make_particle` — subroutine
- `notch_table` — function
- `report_rounds` — subroutine
- `sum_g2` — function
- `test_atom_set` — subroutine
- `test_convolution_identity` — subroutine
- `test_fit_atom` — subroutine
- `test_gauss_filter` — subroutine
- `test_kernels` — subroutine
- `test_noise` — subroutine
- `test_particle_gradient` — subroutine
- `test_particle_uniform` — subroutine
- `test_regress_profile` — subroutine
- `test_render` — subroutine
- `test_robust_spread` — subroutine
- `test_transfer_kernels` — subroutine
- `width_rise` — function

Local symbols (enclosing procedure scope):
- `test_transfer_kernels::table_value` — function

---
## Module: simple_nano_species

Files:
- `main/nano/simple_nano_species.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_srch_sort_loc`
- `simple_xd_gmm`

Public symbols:
- `calibrate_threshold` — subroutine
- `class_posterior` — subroutine
- `class_separation` — function
- `enclosed_fraction` — function
- `expected_false_count` — function
- `fit_gauss_width` — subroutine
- `fit_species_mixture` — subroutine
- `local_maxima` — subroutine
- `maxima_shape` — function

Private symbols:
- `fit_k` — subroutine
- `fit_k_cov` — subroutine
- `is_admissible` — function
- `median_dp` — function
- `sort_desc` — function

Local symbols (enclosing procedure scope):
- `fit_gauss_width::amp_at` — function
- `fit_gauss_width::cost` — function
- `fit_k::fit_from` — subroutine
- `fit_k_cov::alternate` — subroutine
- `fit_k_cov::shifted` — function
- `local_maxima::neighbour_mask` — function

---
## Module: simple_nano_species_tester

Files:
- `main/nano/simple_nano_species_tester.f90`

Uses:
- `simple_defs`
- `simple_nano_species`
- `simple_rnd`
- `simple_test_utils`

Public symbols:
- `run_all_nano_species_tests` — subroutine

Private symbols:
- `draw` — subroutine
- `gauss_samples` — subroutine
- `test_calibration` — subroutine
- `test_coordination_term` — subroutine
- `test_enclosed_fraction` — subroutine
- `test_four_classes` — subroutine
- `test_gauss_width` — subroutine
- `test_local_maxima` — subroutine
- `test_multiplicative_classes` — subroutine
- `test_one_class` — subroutine
- `test_posterior_error` — subroutine
- `test_separation` — subroutine
- `test_size_floor` — subroutine
- `test_three_classes` — subroutine

---
## Module: simple_nanoparticle

Files:
- `main/nano/simple_nanoparticle.f90`

Uses:
- `simple_atoms`
- `simple_bench_io`
- `simple_core_module_api`
- `simple_defs_atoms`
- `simple_image`
- `simple_image_bin`
- `simple_image_msk`
- `simple_nano_density`
- `simple_nano_species`
- `simple_nanoparticle_utils`
- `simple_parameters`

Public symbols:
- `nanoparticle` — type

Private symbols:
- `atom_stats` — type
- `atominfo2centers` — function
- `atominfo2centers_A` — function
- `binarize_and_find_centers` — subroutine
- `calc_anisotropic_disp` — subroutine
- `calc_isotropic_disp` — subroutine
- `calc_longest_atm_dist` — subroutine
- `center_on_atom` — subroutine
- `check_neighbors_cn` — subroutine
- `contact_scores` — subroutine
- `conv_denoise` — subroutine
- `discard_atoms` — subroutine
- `discard_small_ccs` — subroutine
- `empty_zone` — subroutine
- `fillin_atominfo` — subroutine
- `find_centers` — subroutine
- `fit_density_model` — subroutine
- `get_img` — subroutine
- `get_img_raw` — subroutine
- `get_ldim` — subroutine
- `get_natoms` — function
- `get_valid_corrs` — function
- `identify_atomic_pos` — subroutine
- `identify_lattice_params` — subroutine
- `kill` — subroutine
- `masscen` — function
- `new` — subroutine
- `np_binary_lattice` — function
- `np_lattice_bond` — function
- `pack_instance4stats` — subroutine
- `set_atomic_coords_from_pdb` — subroutine
- `set_atomic_coords_from_xyz` — subroutine
- `set_coords4stats` — subroutine
- `set_img` — subroutine
- `set_ncc` — subroutine
- `set_nn_length_scales` — subroutine
- `simulate_atoms` — subroutine
- `split_atoms` — subroutine
- `update_ncc` — subroutine
- `validate_atoms` — subroutine
- `write_atominfo` — subroutine
- `write_centers_1` — subroutine
- `write_centers_2` — subroutine
- `write_centers_aniso` — subroutine
- `write_cn_stats` — subroutine
- `write_csv_files` — subroutine
- `write_individual_atoms` — subroutine
- `write_np_stats` — subroutine

Local symbols (enclosing procedure scope):
- `binarize_and_find_centers::set_pseudo_bfacs` — subroutine
- `binarize_and_find_centers::t2c` — function
- `discard_atoms::relabel` — subroutine
- `discard_atoms::remove_lowly_contacted` — subroutine
- `discard_atoms::remove_small_and_lowly_contacted` — subroutine
- `discard_atoms::zero_component` — subroutine
- `fillin_atominfo::calc_cn_stats` — subroutine
- `fillin_atominfo::calc_zscore` — subroutine
- `fillin_atominfo::write_2D_slice` — subroutine
- `fillin_atominfo::write_cn_atoms` — subroutine
- `fit_density_model::band_passed` — function
- `fit_density_model::build_masks` — subroutine
- `fit_density_model::calibrate` — subroutine
- `fit_density_model::centroid` — function
- `fit_density_model::compute_zmap` — subroutine
- `fit_density_model::corr_in` — function
- `fit_density_model::count_extrema` — subroutine
- `fit_density_model::detect_residual_level` — subroutine
- `fit_density_model::kv` — subroutine
- `fit_density_model::map_corrs` — subroutine
- `fit_density_model::mark_ball` — subroutine
- `fit_density_model::prune_weak` — subroutine
- `fit_density_model::read_transfer_tab` — subroutine
- `fit_density_model::sdev_of` — function
- `fit_density_model::set_noise` — subroutine
- `fit_density_model::stage_b_mask` — subroutine
- `fit_density_model::weak_valid_corr` — subroutine
- `fit_density_model::write_atoms_csv` — subroutine
- `fit_density_model::write_model_txt` — subroutine
- `fit_density_model::write_transfer_tsv` — subroutine
- `fit_density_model::write_vol` — subroutine
- `simulate_atoms::set_atom` — subroutine
- `split_atoms::brightest` — function
- `split_atoms::claim` — subroutine
- `split_atoms::release_unclaimed` — subroutine
- `split_atoms::split_atom` — subroutine
- `write_centers_1::set_atm_info` — subroutine
- `write_centers_1::set_model_info` — subroutine

---
## Module: simple_nanoparticle_utils

Files:
- `main/nano/simple_nanoparticle_utils.f90`

Uses:
- `simple_atoms`
- `simple_core_module_api`
- `simple_defs_atoms`
- `simple_image`
- `simple_linalg`

Public symbols:
- `atm_rmsd_stats` — function
- `atoms_register` — subroutine
- `binary_lattice` — function
- `calc_contact_scores` — subroutine
- `contact_cutoff` — function
- `dists_btw_common` — function
- `est_nn_dist` — function
- `find_atoms_subset` — subroutine
- `find_couples` — subroutine
- `find_rMax` — function
- `fit_lattice` — subroutine
- `Kabsch_algo` — subroutine
- `lattice_bond` — function
- `lattice_cutoff` — function
- `phasecorr_one_atom` — subroutine
- `read_pdb2matrix` — subroutine
- `remove_atoms` — subroutine
- `run_cn_analysis` — subroutine
- `strain_analysis` — subroutine
- `thres_detect_conv_atom_denoised` — subroutine
- `write_matrix2pdb` — subroutine

Private symbols:
- `el_symbol` — function
- `write_ideal_lattice_pdb` — subroutine

Local symbols (enclosing procedure scope):
- `fit_lattice::nnbdl_lat` — subroutine
- `fit_lattice::qr` — subroutine
- `Kabsch_algo::determinant` — function
- `strain_analysis::gaussian_kernel` — subroutine

---
## Module: simple_neighs

Files:
- `utils/math/simple_neighs.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `neigh_4_3D_coords` — subroutine
- `neigh_4_3D_imat` — subroutine
- `neigh_8_3D_coords` — subroutine
- `neigh_8_3D_imat` — subroutine
- `neigh_8_3D_rmat` — subroutine
- `neigh_8_coords` — subroutine
- `neigh_8_imat` — subroutine
- `neigh_8_rmat` — subroutine

---
## Module: simple_nlcg_solver

Files:
- `main/opt/simple_nlcg_solver.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_error`

Public symbols:
- `nlcg_objective` — type
- `nlcg_solve` — subroutine
- `nlcg_solver_options` — type
- `nlcg_solver_outcome` — type

Private symbols:
- `nlcg_curvature_i` — function
- `nlcg_dot_i` — function
- `nlcg_eval_i` — subroutine
- `nlcg_line_losses_i` — subroutine
- `nlcg_precond_i` — subroutine
- `nlcg_project_i` — subroutine
- `nlcg_size_i` — function
- `nlcg_solver_outcome_kill` — subroutine

Local symbols (enclosing procedure scope):
- `nlcg_solve::finish` — subroutine
- `nlcg_solve::projected_norm` — function
- `nlcg_solve::rises` — function
- `nlcg_solve::try_step` — subroutine

---
## Module: simple_nlcg_solver_tester

Files:
- `main/opt/simple_nlcg_solver_tester.f90`

Uses:
- `simple_defs`
- `simple_nlcg_solver`
- `simple_pcg_solver`
- `simple_test_utils`

Public symbols:
- `run_all_nlcg_solver_tests` — subroutine

Private symbols:
- `lsq_curvature` — function
- `lsq_dot` — function
- `lsq_eval` — subroutine
- `lsq_line_losses` — subroutine
- `lsq_obj` — type
- `lsq_precond` — subroutine
- `lsq_project` — subroutine
- `lsq_resid` — function
- `lsq_size` — function
- `qop_apply` — subroutine
- `qop_dot` — function
- `qop_precond` — subroutine
- `qop_size` — function
- `quad_curvature` — function
- `quad_dot` — function
- `quad_eval` — subroutine
- `quad_line_losses` — subroutine
- `quad_obj` — type
- `quad_op` — type
- `quad_precond` — subroutine
- `quad_project` — subroutine
- `quad_size` — function
- `test_nonconvex_least_squares` — subroutine
- `test_quadratic_is_linear_pcg` — subroutine

---
## Module: simple_nrtxtfile

Files:
- `fileio/simple_nrtxtfile.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `nrtxtfile` — type

Private symbols:
- `get_ndatalines` — function
- `get_nrecs_per_line` — function
- `kill` — subroutine
- `new` — subroutine
- `readNextDataLine` — subroutine
- `writeCommentLine` — subroutine
- `writeDataLineInt` — subroutine
- `writeDataLineReal` — subroutine

---
## Module: simple_nu_filter

Files:
- `main/nu_filt/simple_nu_filter.f90`
- `main/nu_filt/simple_nu_filter_apply.f90`
- `main/nu_filt/simple_nu_filter_bank.f90`
- `main/nu_filt/simple_nu_filter_envmask.f90`
- `main/nu_filt/simple_nu_filter_evidence.f90`
- `main/nu_filt/simple_nu_filter_potts.f90`
- `main/nu_filt/simple_nu_filter_sharpen.f90`
- `main/nu_filt/simple_nu_filter_state.f90`
- `main/nu_filt/simple_nu_filter_stats.f90`

Uses:
- `ieee_arithmetic`
- `simple_butterworth`
- `simple_core_module_api`
- `simple_image`
- `simple_image_msk`
- `simple_neighs`
- `simple_nu_filter_vars`
- `simple_tent_smooth`

Public symbols:
- `analyze_filtmap_neighbor_continuity` — subroutine
- `analyze_filtmap_neighbor_continuity` — subroutine
- `assert_nu_evidence_replay_ready` — subroutine
- `assert_nu_evidence_replay_ready` — subroutine
- `build_nu_evidence_state` — subroutine
- `build_nu_evidence_state` — subroutine
- `calc_filtmap_lowpass_histogram` — subroutine
- `calc_filtmap_lowpass_histogram` — subroutine
- `calc_filtmap_lowpass_stats` — subroutine
- `calc_filtmap_lowpass_stats` — subroutine
- `calc_nu_evidence_margin` — subroutine
- `calc_nu_evidence_margin` — subroutine
- `cleanup_nu_filter` — subroutine
- `cleanup_nu_filter` — subroutine
- `clear_nu_evidence_null_shell` — subroutine
- `clear_nu_evidence_null_shell` — subroutine
- `clear_nu_solvent_envelope` — subroutine
- `clear_nu_solvent_envelope` — subroutine
- `get_nu_evidence_summary` — subroutine
- `get_nu_evidence_summary` — subroutine
- `get_nu_filter_bank_finest_lp` — function
- `get_nu_filter_bank_finest_lp` — function
- `get_nu_filter_timings` — subroutine
- `get_nu_filter_timings` — subroutine
- `get_nu_filtmap_finest_selected_lp` — function
- `get_nu_filtmap_finest_selected_lp` — function
- `nu_envmask_params` — type
- `nu_envmask_stats` — type
- `nu_evidence_envelope` — subroutine
- `nu_evidence_envelope` — subroutine
- `nu_evidence_sharpen_vol` — subroutine
- `nu_evidence_sharpen_vol` — subroutine
- `nu_evidence_state` — type
- `nu_evidence_state_is_valid` — function
- `nu_evidence_state_is_valid` — function
- `nu_evidence_summary` — type
- `nu_filter_vols` — subroutine
- `nu_filter_vols` — subroutine
- `optimize_nu_cutoff_finds` — subroutine
- `optimize_nu_cutoff_finds` — subroutine
- `pack_filtmap_lowpass_limits` — subroutine
- `pack_filtmap_lowpass_limits` — subroutine
- `print_filtmap_lowpass_histogram` — subroutine
- `print_filtmap_lowpass_histogram` — subroutine
- `print_nu_envmask_stats` — subroutine
- `print_nu_envmask_stats` — subroutine
- `print_nu_evidence_summary` — subroutine
- `print_nu_evidence_summary` — subroutine
- `print_nu_filtmap_lowpass_stats` — subroutine
- `print_nu_filtmap_lowpass_stats` — subroutine
- `set_nu_evidence_null_shell` — subroutine
- `set_nu_evidence_null_shell` — subroutine
- `set_nu_filter_report` — subroutine
- `set_nu_filter_report` — subroutine
- `set_nu_solvent_envelope` — subroutine
- `set_nu_solvent_envelope` — subroutine
- `setup_nu_dmats` — subroutine
- `setup_nu_dmats` — subroutine
- `unpack_nu_evidence_state` — subroutine
- `unpack_nu_evidence_state` — subroutine
- `write_nu_evidence_envmask` — subroutine
- `write_nu_evidence_envmask` — subroutine
- `write_nu_evidence_map` — subroutine
- `write_nu_evidence_map` — subroutine
- `write_nu_local_resolution_map` — subroutine
- `write_nu_local_resolution_map` — subroutine

Private symbols:
- `accumulate_nu_evidence_raw` — subroutine
- `accumulate_nu_evidence_raw` — subroutine
- `active_nu_mask_at` — function
- `add_nu_evidence_density` — subroutine
- `add_nu_evidence_density` — subroutine
- `apply_nu_solvent_clamp` — subroutine
- `apply_nu_solvent_clamp` — subroutine
- `calc_nu_evidence_score` — subroutine
- `calc_nu_evidence_score` — subroutine
- `calc_nu_label_smooth_site_energy` — function
- `calc_nu_label_smooth_site_energy` — function
- `calculate_nu_source_fingerprint` — subroutine
- `calculate_nu_source_fingerprint` — subroutine
- `clamp_nu_filtmap_labels` — subroutine
- `clamp_nu_filtmap_labels` — subroutine
- `cleanup_aux_bank` — subroutine
- `cleanup_aux_bank` — subroutine
- `compose_nu_labels_from_volume` — subroutine
- `count_active_nu_mask` — function
- `count_nu_solvent_clamped` — function
- `cutoff_find_to_lowpass_limit` — function
- `cutoff_find_to_lowpass_limit` — function
- `estimate_evidence_beta` — function
- `estimate_nu_label_smooth_beta` — function
- `estimate_nu_label_smooth_beta` — function
- `evidence_site_energy` — function
- `evidence_unary` — function
- `filter_nu_candidate` — subroutine
- `filter_nu_candidate` — subroutine
- `init_nu_filter` — subroutine
- `init_nu_filter` — subroutine
- `load_nu_smooth_norm` — subroutine
- `load_nu_smooth_norm` — subroutine
- `log_nu_aux_replacement_margin_stats` — subroutine
- `log_nu_candidate_coords` — subroutine
- `log_nu_candidate_coords` — subroutine
- `log_nu_candidate_selection_counts` — subroutine
- `log_nu_candidate_selection_counts` — subroutine
- `log_nu_objective_smoothing_bank` — subroutine
- `lowpass_histogram_value_at` — function
- `nu_candidate_coord_for_label` — function
- `nu_candidate_coord_for_label` — function
- `nu_evidence_baseline_floor` — function
- `nu_evidence_baseline_floor` — function
- `nu_evidence_calibration_mask` — subroutine
- `nu_evidence_null_mask` — subroutine
- `nu_label_is_aux_replacement` — function
- `nu_label_is_aux_replacement` — function
- `nu_label_lowpass_limit` — function
- `nu_label_lowpass_limit` — function
- `nu_label_smooth_color` — function
- `nu_label_smooth_color` — function
- `nu_label_smooth_coord_pair_cost` — function
- `nu_label_smooth_coord_pair_cost` — function
- `nu_label_smooth_is_better` — function
- `nu_label_smooth_is_better` — function
- `nu_label_smooth_neighborhood_cost` — function
- `nu_label_smooth_neighborhood_cost` — function
- `nu_label_smooth_pair_cost` — function
- `nu_label_smooth_pair_cost` — function
- `nu_objective_smooth_radius_angstrom` — function
- `nu_objective_smooth_radius_angstrom` — function
- `nu_objective_smooth_radius_pixels` — function
- `nu_objective_smooth_radius_pixels` — function
- `nu_potts_coord_for_label` — function
- `nu_potts_coord_for_label` — function
- `nu_static_ladder_count` — function
- `nu_static_ladder_count` — function
- `pack_nu_dmat_candidate` — subroutine
- `pack_nu_dmat_candidate` — subroutine
- `pack_nu_full_to_mask` — subroutine
- `pack_nu_full_to_mask` — subroutine
- `pack_nu_raw_candidate` — subroutine
- `pack_nu_raw_candidate` — subroutine
- `prepare_nu_candidate_fts` — subroutine
- `prepare_nu_candidate_fts` — subroutine
- `prepare_nu_smooth_norm` — subroutine
- `refine_nu_candidate_map_ordered_labels` — subroutine
- `refine_nu_candidate_map_ordered_labels` — subroutine
- `regularize_evidence_labels` — subroutine
- `release_nu_filter_unary_storage` — subroutine
- `release_nu_filter_unary_storage` — subroutine
- `release_nu_smooth_norm` — subroutine
- `release_nu_smooth_norm` — subroutine
- `require_valid_stats_mask` — subroutine
- `segment_nu_evidence` — subroutine
- `segment_nu_evidence` — subroutine
- `setup_evidence_candidate_geometry` — subroutine
- `setup_nu_candidate_coords` — subroutine
- `setup_nu_candidate_coords` — subroutine
- `setup_nu_mask_voxels` — subroutine
- `setup_nu_mask_voxels` — subroutine
- `setup_nu_observed_mask` — subroutine
- `setup_nu_observed_mask` — subroutine
- `smooth_nu_objective` — subroutine
- `smooth_nu_objective` — subroutine
- `stash_aux_volumes` — subroutine
- `stash_aux_volumes` — subroutine
- `unpack_nu_raw_candidate` — subroutine
- `unpack_nu_raw_candidate` — subroutine
- `validate_aux_volumes` — subroutine
- `validate_aux_volumes` — subroutine
- `validate_compact_evidence_state` — subroutine

Local symbols (enclosing procedure scope):
- `nu_filter_vols::overlay_aux_label` — subroutine
- `print_nu_envmask_stats::warn_if_invalid` — subroutine
- `refine_nu_candidate_map_ordered_labels::pack_color_voxels` — subroutine
- `regularize_evidence_labels::clamp_evidence_background` — subroutine
- `set_nu_evidence_null_shell::check_input` — subroutine

---
## Module: simple_nu_filter_tester

Files:
- `main/nu_filt/simple_nu_filter_tester.f90`

Uses:
- `simple_butterworth`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_forked_process`
- `simple_image`
- `simple_nu_filter`
- `simple_nu_filter_vars`
- `simple_rnd`
- `simple_tent_smooth`
- `simple_test_utils`

Public symbols:
- `run_all_nu_filter_tests` — subroutine

Private symbols:
- `clamped_setup` — type
- `convolve_axis` — subroutine
- `fill_gaussian` — subroutine
- `make_pair` — subroutine
- `reference_tent_smooth` — subroutine
- `setup_with_clamp` — subroutine
- `smoothed_field` — function
- `support_sphere` — subroutine
- `test_clamp_between_pairs` — subroutine
- `test_coarser_region` — subroutine
- `test_composition_per_label` — subroutine
- `test_evidence_envelope` — subroutine
- `test_label_field_threads` — subroutine
- `test_noise_scale` — subroutine
- `test_objective_closed_form` — subroutine
- `test_ranking_scale_free` — subroutine
- `test_tent_smooth_exact` — subroutine

Local symbols (enclosing procedure scope):
- `reference_tent_smooth::filter_1d` — subroutine
- `reference_tent_smooth::pass` — subroutine
- `test_composition_per_label::count_label_diffs` — function
- `test_evidence_envelope::lcg_uniform` — function
- `test_objective_closed_form::new_from` — subroutine

---
## Module: simple_nu_filter_vars

Files:
- `main/nu_filt/simple_nu_filter_vars.f90`

Uses:
- `simple_image`
- `simple_timer`

---
## Module: simple_nu_state_filter

Files:
- `main/volume/simple_nu_state_filter.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_image_msk`
- `simple_nu_filter`
- `simple_parameters`

Public symbols:
- `nonuniform_filter_state` — subroutine
- `nu_aux_member` — function
- `nu_state_filter_timings` — type

Local symbols (enclosing procedure scope):
- `nonuniform_filter_state::add_filter_split` — subroutine
- `nonuniform_filter_state::cleanup_nu_aux_images` — subroutine
- `nonuniform_filter_state::record_nu_alignment_lowpass_limit` — subroutine
- `nonuniform_filter_state::write_nonuniform_outputs` — subroutine

---
## Module: simple_online_var

Files:
- `utils/math/simple_online_var.f90`

Uses:
- `simple_defs`

Public symbols:
- `online_var` — type

Private symbols:
- `add` — subroutine
- `get_mean` — function
- `get_var` — function

---
## Module: simple_online_var_tester

Files:
- `utils/math/simple_online_var_tester.f90`

Uses:
- `simple_defs`
- `simple_online_var`
- `simple_rnd`
- `simple_stat`
- `simple_test_utils`

Public symbols:
- `run_all_online_var_tests` — subroutine

Private symbols:
- `test_against_moment` — subroutine
- `test_closed_forms` — subroutine
- `test_few_samples` — subroutine
- `test_large_offset` — subroutine

---
## Module: simple_openmp_offload_tester

Files:
- `utils/simple_openmp_offload_tester.f90`

Uses:
- `iso_c_binding`
- `omp_lib`
- `omp_lib_kinds`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_fftw3`
- `simple_gpu_utils`
- `simple_kbinterpol`
- `simple_parameters`

Public symbols:
- `run_openmp_offload_tests` — subroutine

Private symbols:
- `dummytype` — type

Local symbols (enclosing procedure scope):
- `run_openmp_offload_tests::async_test` — subroutine
- `run_openmp_offload_tests::banner` — subroutine
- `run_openmp_offload_tests::init_dummy_1d` — subroutine
- `run_openmp_offload_tests::init_dummy_2d` — subroutine
- `run_openmp_offload_tests::init_dummy_3d` — subroutine
- `run_openmp_offload_tests::numerical_test` — subroutine
- `run_openmp_offload_tests::persistence_test` — subroutine
- `run_openmp_offload_tests::persistence_test2` — subroutine
- `run_openmp_offload_tests::pointer_kernel` — subroutine
- `run_openmp_offload_tests::print_cpair` — subroutine
- `run_openmp_offload_tests::relerr_c1d` — function
- `run_openmp_offload_tests::relerr_c2d` — function
- `run_openmp_offload_tests::report_fftw` — subroutine
- `run_openmp_offload_tests::test_cublas_sgemm` — subroutine
- `run_openmp_offload_tests::test_fftw_vs_cufft_1d_in_place` — subroutine
- `run_openmp_offload_tests::test_fftw_vs_cufft_1d_many` — subroutine
- `run_openmp_offload_tests::test_fftw_vs_cufft_1d_many_c2r` — subroutine
- `run_openmp_offload_tests::test_fftw_vs_cufft_2d_in_place` — subroutine
- `run_openmp_offload_tests::test_fftw_vs_cufft_2d_in_place_roundtrip_timing` — subroutine
- `run_openmp_offload_tests::test_fftw_vs_cufft_2d_many` — subroutine
- `run_openmp_offload_tests::test_fftw_vs_cufft_3d_in_place_roundtrip_timing` — subroutine
- `run_openmp_offload_tests::test_fftw_vs_cufft_movie` — subroutine
- `run_openmp_offload_tests::test_kb` — subroutine
- `run_openmp_offload_tests::test_pointers` — subroutine

---
## Module: simple_opt_de

Files:
- `main/opt/simple_opt_de.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `opt_de` — type

Private symbols:
- `de_minimize` — subroutine
- `kill_de` — subroutine
- `new_de` — subroutine

---
## Module: simple_opt_factory

Files:
- `main/opt/simple_opt_factory.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_de`
- `simple_opt_lbfgsb`
- `simple_opt_simplex`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `opt_factory` — type

Private symbols:
- `new` — subroutine

---
## Module: simple_opt_filter

Files:
- `utils/filter/simple_opt_filter.f90`

Uses:
- `simple_butterworth`
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `estimate_lplim3D` — subroutine
- `estimate_lplims2D` — subroutine

Private symbols:
- `estimate_lplim_find` — subroutine
- `estimate_lplim_lp` — subroutine

---
## Module: simple_opt_lbfgsb

Files:
- `main/opt/simple_opt_lbfgsb.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `opt_lbfgsb` — type

Private symbols:
- `daxpy` — subroutine
- `dcopy` — subroutine
- `ddot` — function
- `dscal` — subroutine
- `kill_opt_lbfgsb` — subroutine
- `lbfgsb_minimize` — subroutine
- `new_opt_lbfgsb` — subroutine

Local symbols (enclosing procedure scope):
- `lbfgsb_minimize::active` — subroutine
- `lbfgsb_minimize::bmv` — subroutine
- `lbfgsb_minimize::cauchy` — subroutine
- `lbfgsb_minimize::cmprlb` — subroutine
- `lbfgsb_minimize::dcsrch` — subroutine
- `lbfgsb_minimize::dcstep` — subroutine
- `lbfgsb_minimize::dpofa` — subroutine
- `lbfgsb_minimize::dtrsl` — subroutine
- `lbfgsb_minimize::errclb` — subroutine
- `lbfgsb_minimize::formk` — subroutine
- `lbfgsb_minimize::formt` — subroutine
- `lbfgsb_minimize::freev` — subroutine
- `lbfgsb_minimize::hpsolb` — subroutine
- `lbfgsb_minimize::lbfgsbmin` — subroutine
- `lbfgsb_minimize::lnsrlb` — subroutine
- `lbfgsb_minimize::mainlb` — subroutine
- `lbfgsb_minimize::matupd` — subroutine
- `lbfgsb_minimize::projgr` — subroutine
- `lbfgsb_minimize::setulb` — subroutine
- `lbfgsb_minimize::subsm` — subroutine

---
## Module: simple_opt_mask

Files:
- `main/image_processing/simple_opt_mask.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `estimate_spher_mask` — subroutine

---
## Module: simple_opt_simplex

Files:
- `main/opt/simple_opt_simplex.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`
- `simple_opt_subs`
- `simple_optimizer`

Public symbols:
- `opt_simplex` — type

Private symbols:
- `kill_opt_simplex` — subroutine
- `new_opt_simplex` — subroutine
- `simplex_minimize` — subroutine

Local symbols (enclosing procedure scope):
- `simplex_minimize::init` — subroutine

---
## Module: simple_opt_spec

Files:
- `main/opt/simple_opt_spec.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `costfun` — function
- `opt_spec` — type

Private symbols:
- `change_opt` — subroutine
- `costfun_8` — function
- `eval_df_4` — subroutine
- `eval_df_8` — subroutine
- `eval_f_4` — function
- `eval_f_8` — function
- `eval_fdf_4` — subroutine
- `eval_fdf_8` — subroutine
- `fdfcostfun` — subroutine
- `fdfcostfun_8` — subroutine
- `gcostfun` — subroutine
- `gcostfun_8` — subroutine
- `kill` — subroutine
- `opt_callback` — subroutine
- `set_costfun` — subroutine
- `set_costfun_8` — subroutine
- `set_fdfcostfun` — subroutine
- `set_fdfcostfun_8` — subroutine
- `set_gcostfun` — subroutine
- `set_gcostfun_8` — subroutine
- `set_inipop` — subroutine
- `set_limits` — subroutine
- `set_limits_init` — subroutine
- `specify` — subroutine

Local symbols (enclosing procedure scope):
- `set_costfun::fun` — function
- `set_costfun_8::fun` — function
- `set_fdfcostfun::fun` — subroutine
- `set_fdfcostfun_8::fun` — subroutine
- `set_gcostfun::fun` — subroutine
- `set_gcostfun_8::fun` — subroutine

---
## Module: simple_opt_subs

Files:
- `main/opt/simple_opt_subs.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`

Public symbols:
- `amoeba` — subroutine

Local symbols (enclosing procedure scope):
- `amoeba::amoeba_private` — subroutine
- `amoeba::amotry` — function

---
## Module: simple_opt_tester

Files:
- `main/opt/simple_opt_tester.f90`

Uses:
- `simple_defs`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`
- `simple_pcg_solver`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_opt_tests` — subroutine

Private symbols:
- `dense_pcg_apply` — subroutine
- `dense_pcg_dot` — function
- `dense_pcg_operator` — type
- `dense_pcg_precond` — subroutine
- `dense_pcg_size` — function
- `dircost` — function
- `dircost_grad` — subroutine
- `eval_recorder` — type
- `quad1d` — function
- `quad1d_grad` — subroutine
- `quad2d` — function
- `quad2d_grad` — subroutine
- `quad2d_rec` — function
- `rosenbrock` — function
- `rosenbrock_grad` — subroutine
- `set_spd_operator` — subroutine
- `test_de_npop_inipop` — subroutine
- `test_de_quadratic` — subroutine
- `test_lbfgsb_bounded` — subroutine
- `test_lbfgsb_direction` — subroutine
- `test_lbfgsb_quadratic` — subroutine
- `test_lbfgsb_rosenbrock` — subroutine
- `test_pcg_cold_restart` — subroutine
- `test_pcg_dense_spd` — subroutine
- `test_pcg_diagonal_one_step` — subroutine
- `test_simplex_quadratic` — subroutine
- `test_spec_bookkeeping` — subroutine
- `test_start_at_origin` — subroutine

---
## Module: simple_optics_groups

Files:
- `main/stream/shared/simple_optics_groups.f90`

Uses:
- `simple_defs`
- `simple_math`
- `simple_ori`
- `simple_sp_project`
- `simple_srch_sort_loc`
- `simple_string_utils`

Public symbols:
- `assign_optics_groups` — subroutine

Private symbols:
- `kept_ids` — function
- `link_shifts` — subroutine

Local symbols (enclosing procedure scope):
- `assign_optics_groups::id_lt` — function
- `kept_ids::pair_first` — function
- `kept_ids::pair_lt` — function
- `link_shifts::cell_lt` — function
- `link_shifts::find_cell` — function
- `link_shifts::link_cells` — subroutine
- `link_shifts::root` — function
- `link_shifts::unite` — subroutine

---
## Module: simple_optics_groups_tester

Files:
- `main/stream/shared/simple_optics_groups_tester.f90`

Uses:
- `simple_optics_groups`
- `simple_sp_project`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_optics_groups_tests` — subroutine

Private symbols:
- `add_mic` — subroutine
- `check_two_groups` — subroutine
- `make_line` — subroutine
- `make_mics` — subroutine
- `set_mic` — subroutine
- `test_chain_is_one_group` — subroutine
- `test_ctf_constants_from_first_accepted` — subroutine
- `test_equal_shifts_one_group` — subroutine
- `test_group_offset` — subroutine
- `test_ids_kept_across_regroupings` — subroutine
- `test_two_shift_clusters_with_beamtilt` — subroutine
- `test_two_shift_clusters_without_beamtilt` — subroutine

---
## Module: simple_optics_maps

Files:
- `main/stream/shared/simple_optics_maps.f90`

Uses:
- `simple_defs_fname`
- `simple_fileio`
- `simple_nrtxtfile`
- `simple_sp_project`
- `simple_stream_utils`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `copy_project_with_optics_map` — subroutine
- `import_latest_optics_map` — function
- `latest_optics_map_id` — function
- `latest_optics_map_table` — function
- `publish_optics_map` — subroutine

Private symbols:
- `map_prefix` — function

---
## Module: simple_optics_maps_tester

Files:
- `main/stream/shared/simple_optics_maps_tester.f90`

Uses:
- `simple_defs_fname`
- `simple_fileio`
- `simple_optics_maps`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_optics_maps_tests` — subroutine

Private symbols:
- `make_grouped_mics` — subroutine
- `map_exists` — function
- `new_test_dir` — function
- `test_copy_with_optics_map` — subroutine
- `test_empty_directory` — subroutine
- `test_import_latest` — subroutine
- `test_publish_keeps_newest` — subroutine

---
## Module: simple_optimizer

Files:
- `main/opt/simple_optimizer.f90`

Uses:
- `simple_opt_spec`

Public symbols:
- `optimizer` — type

Private symbols:
- `generic_kill` — subroutine
- `generic_minimize` — subroutine
- `generic_new` — subroutine

---
## Module: simple_ori

Files:
- `main/ori/simple_ori.f90`

Uses:
- `simple_ori_api`
- `simple_sauron`

Public symbols:
- `ori` — type

Private symbols:
- `append_ori` — subroutine
- `calc_offset2D` — subroutine
- `chash2ori` — subroutine
- `compeuler` — subroutine
- `compose2dshift3d` — subroutine
- `compose3d2d` — subroutine
- `constructor` — function
- `copy` — subroutine
- `delete_2Dclustering` — subroutine
- `delete_3Dalignment` — subroutine
- `delete_entry` — subroutine
- `e1get` — function
- `e1set` — subroutine
- `e2get` — function
- `e2set` — subroutine
- `e3get` — function
- `e3set` — subroutine
- `euldist` — function
- `exists` — function
- `gau_rnd_shift` — subroutine
- `geodesic_dist_frobdev` — function
- `geodesic_dist_trace` — function
- `get` — function
- `get_2Dshift` — function
- `get_3Dshift` — function
- `get_class` — function
- `get_ctfvars` — function
- `get_dfx` — function
- `get_dfy` — function
- `get_eo` — function
- `get_euler` — function
- `get_fromp` — function
- `get_int` — function
- `get_keys` — function
- `get_mat` — function
- `get_normal` — function
- `get_proj` — function
- `get_sampled` — function
- `get_state` — function
- `get_static` — subroutine
- `get_str` — function
- `get_top` — function
- `get_updatecnt` — function
- `getter_int` — subroutine
- `getter_real` — subroutine
- `getter_str` — subroutine
- `has_been_searched` — function
- `inplrotdist` — function
- `is_particle` — function
- `ischar` — function
- `isstate` — function
- `isstatezero` — function
- `isthere` — function
- `kill` — subroutine
- `kill_chash` — subroutine
- `kill_hash` — subroutine
- `map3dshift22d` — subroutine
- `mirror2d` — subroutine
- `mirror3d` — subroutine
- `new_ori` — subroutine
- `ori2chash` — function
- `ori2json` — subroutine
- `ori2prec` — subroutine
- `ori2str` — function
- `ori_from_rotmat` — subroutine
- `ori_strlen_trim` — function
- `pparms2str` — function
- `pparms_strlen` — function
- `prec2ori` — subroutine
- `print_mat` — subroutine
- `print_ori` — subroutine
- `read` — subroutine
- `reject` — subroutine
- `reset_pparms` — subroutine
- `rnd_euler_full` — subroutine
- `rnd_euler_lims` — subroutine
- `rnd_euler_near` — subroutine
- `rnd_euler_near_projinpl` — subroutine
- `rnd_euler_offset` — subroutine
- `rnd_inpl` — subroutine
- `rnd_ori` — subroutine
- `rnd_shift` — subroutine
- `round_shifts` — subroutine
- `set_char` — subroutine
- `set_class` — subroutine
- `set_ctfvars` — subroutine
- `set_dfx` — subroutine
- `set_dfy` — subroutine
- `set_dp` — subroutine
- `set_euler` — subroutine
- `set_int` — subroutine
- `set_ogid` — subroutine
- `set_real` — subroutine
- `set_shift` — subroutine
- `set_state` — subroutine
- `set_stkind` — subroutine
- `set_str` — subroutine
- `str2ori` — subroutine
- `swape1e3` — subroutine
- `transfer_2Dparams` — subroutine
- `transfer_3Dparams` — subroutine
- `transp` — subroutine
- `write` — subroutine
- `write2bild` — subroutine

Local symbols (enclosing procedure scope):
- `compose3d2d::transfmat2inpls` — subroutine

---
## Module: simple_ori_api

Files:
- `main/apis/simple_ori_api.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_chash`
- `simple_defs`
- `simple_defs_ori`
- `simple_error`
- `simple_fileio`
- `simple_hash`
- `simple_is_check_assert`
- `simple_linalg`
- `simple_math`
- `simple_math_ft`
- `simple_nrtxtfile`
- `simple_ori_utils`
- `simple_ran_tabu`
- `simple_rnd`
- `simple_srch_sort_loc`
- `simple_stat`
- `simple_string`
- `simple_string_utils`
- `simple_type_defs`

---
## Module: simple_ori_tester

Files:
- `main/ori/simple_ori_tester.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_chash`
- `simple_defs`
- `simple_defs_ori`
- `simple_linalg`
- `simple_ori`
- `simple_ori_utils`
- `simple_string`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_ori_tests` — subroutine

Private symbols:
- `has_key` — function
- `test_2Dshift_set_get_and_rounding` — subroutine
- `test_append_ori` — subroutine
- `test_cartesian_score_fields` — subroutine
- `test_copy_nonptcl_to_ptcl` — subroutine
- `test_copy_ptcl_to_nonptcl` — subroutine
- `test_copy_ptcl_to_ptcl` — subroutine
- `test_ctfvars_roundtrip` — subroutine
- `test_ctor_and_flags` — subroutine
- `test_delete_2Dclustering` — subroutine
- `test_delete_3Dalignment` — subroutine
- `test_delete_entry` — subroutine
- `test_euler_basic_set_get` — subroutine
- `test_euler_compose_vs_compeuler` — subroutine
- `test_euler_matrix_roundtrip` — subroutine
- `test_genfire_conversion` — subroutine
- `test_geodesic_metrics` — subroutine
- `test_geodesic_trace_self_distance` — subroutine
- `test_get_dfx_dfy_and_setters` — subroutine
- `test_get_keys` — subroutine
- `test_isthere_and_ischar` — subroutine
- `test_mirror2d_involution` — subroutine
- `test_mirror3d_involution` — subroutine
- `test_normal_and_mat_consistency` — subroutine
- `test_ori2chash_chash2ori_roundtrip` — subroutine
- `test_ori2json` — subroutine
- `test_ori2str_and_str2ori_roundtrip` — subroutine
- `test_ori_from_rotmat` — subroutine
- `test_ori_strlen_trim` — subroutine
- `test_reject` — subroutine
- `test_state_class_proj_eo_sampled_updatecnt` — subroutine
- `test_transfer_2Dparams` — subroutine
- `test_transfer_3Dparams` — subroutine
- `test_transp_involution` — subroutine

---
## Module: simple_ori_utils

Files:
- `main/ori/simple_ori_utils.f90`

Uses:
- `simple_defs`
- `simple_linalg`
- `simple_math`
- `simple_rnd`

Public symbols:
- `dm2euler` — function
- `euler2genfire` — function
- `euler2m` — function
- `euler_compose` — subroutine
- `euler_dist` — function
- `euler_inplrotdist` — function
- `euler_mirror` — subroutine
- `genfire2euler` — function
- `geodesic_frobdev` — function
- `m2euler` — function
- `make_transfmat` — function
- `rnd_romat` — subroutine

Private symbols:
- `cross3` — function
- `euler_normal` — function
- `genfire_rotmat` — function
- `rotmat` — function

---
## Module: simple_oris

Files:
- `main/ori/simple_oris.f90`
- `main/ori/simple_oris_dists.f90`
- `main/ori/simple_oris_getters.f90`
- `main/ori/simple_oris_io.f90`
- `main/ori/simple_oris_life.f90`
- `main/ori/simple_oris_neig.f90`
- `main/ori/simple_oris_reshape.f90`
- `main/ori/simple_oris_sampling.f90`
- `main/ori/simple_oris_setters.f90`
- `main/ori/simple_oris_stats.f90`
- `main/ori/simple_oris_transform.f90`

Uses:
- `simple_jpg`
- `simple_linalg`
- `simple_map_reduce`
- `simple_ori`
- `simple_ori_api`
- `simple_ran_tabu`
- `simple_rnd`

Public symbols:
- `class_sample_quotas` — subroutine
- `class_sample_quotas` — subroutine
- `class_sample_sweep` — function
- `class_sample_sweep` — function
- `oris` — type

Private symbols:
- `add_shift2class` — subroutine
- `add_shift2class` — subroutine
- `alloc_unit_quotas` — subroutine
- `any_active_updated` — function
- `any_state_zero` — function
- `any_state_zero` — function
- `append_at` — subroutine
- `append_at` — subroutine
- `append_oris` — subroutine
- `append_oris` — subroutine
- `calc_avg_offset2D` — subroutine
- `calc_avg_offset2D` — subroutine
- `calc_avg_offset3D` — subroutine
- `calc_avg_offset3D` — subroutine
- `calc_sum` — subroutine
- `calc_sum` — subroutine
- `clean_entry` — subroutine
- `clean_entry` — subroutine
- `compress` — subroutine
- `compress` — subroutine
- `constructor` — function
- `constructor` — function
- `copy_as_kind` — subroutine
- `copy_as_kind` — subroutine
- `copy_same_kind` — subroutine
- `copy_same_kind` — subroutine
- `count_state_gt_zero` — function
- `count_state_gt_zero` — function
- `counts2fracs` — subroutine
- `delete` — subroutine
- `delete` — subroutine
- `delete_2Dclustering_all` — subroutine
- `delete_2Dclustering_all` — subroutine
- `delete_2Dclustering_by_ind` — subroutine
- `delete_2Dclustering_by_ind` — subroutine
- `delete_3Dalignment` — subroutine
- `delete_3Dalignment` — subroutine
- `delete_entry_all` — subroutine
- `delete_entry_all` — subroutine
- `delete_entry_by_ind` — subroutine
- `delete_entry_by_ind` — subroutine
- `detect_peaks` — subroutine
- `detect_peaks` — subroutine
- `discretize` — subroutine
- `discretize` — subroutine
- `diststat_between` — subroutine
- `diststat_between` — subroutine
- `diststat_within` — subroutine
- `diststat_within` — subroutine
- `e1get` — function
- `e1get` — function
- `e1set` — subroutine
- `e1set` — subroutine
- `e2get` — function
- `e2get` — function
- `e2set` — subroutine
- `e2set` — subroutine
- `e3get` — function
- `e3get` — function
- `e3set` — subroutine
- `e3set` — subroutine
- `e3swapsgn` — subroutine
- `e3swapsgn` — subroutine
- `empty_allowed` — function
- `euldist_inds` — function
- `euldist_inds` — function
- `euldist_ori` — function
- `euldist_ori` — function
- `exists` — function
- `exists` — function
- `expand_classes` — subroutine
- `expand_classes` — subroutine
- `extract_subset_inds` — function
- `extract_subset_inds` — function
- `extract_subset_range` — function
- `extract_subset_range` — function
- `extract_subspace` — subroutine
- `extract_subspace` — subroutine
- `extremal_bound` — function
- `extremal_bound` — function
- `find_angres` — function
- `find_angres` — function
- `find_closest_proj` — function
- `find_closest_proj` — function
- `gau_rnd_shifts` — subroutine
- `gau_rnd_shifts` — subroutine
- `gen_balanced_partitions` — subroutine
- `gen_balanced_partitions` — subroutine
- `gen_mask` — subroutine
- `gen_mask` — subroutine
- `gen_ptcl_mask` — function
- `gen_ptcl_mask` — function
- `get` — function
- `get` — function
- `get_2Dshift` — function
- `get_2Dshift` — function
- `get_all` — function
- `get_all` — function
- `get_all_asint` — function
- `get_all_asint` — function
- `get_all_normals` — function
- `get_all_normals` — function
- `get_all_rmats` — function
- `get_all_rmats` — function
- `get_all_sampled` — function
- `get_all_sampled` — function
- `get_arr` — function
- `get_arr` — function
- `get_avg` — function
- `get_avg` — function
- `get_class` — function
- `get_class` — function
- `get_class_sample_stats` — subroutine
- `get_class_sample_stats` — subroutine
- `get_ctfvars` — function
- `get_ctfvars` — function
- `get_dfx` — function
- `get_dfx` — function
- `get_dfy` — function
- `get_dfy` — function
- `get_eo` — function
- `get_eo` — function
- `get_euler` — function
- `get_euler` — function
- `get_fromp` — function
- `get_fromp` — function
- `get_group_update_counts` — subroutine
- `get_group_update_counts` — subroutine
- `get_half` — function
- `get_half` — function
- `get_int` — function
- `get_int` — function
- `get_mat` — function
- `get_mat` — function
- `get_n` — function
- `get_n` — function
- `get_neven` — function
- `get_neven` — function
- `get_nevenodd` — function
- `get_nevenodd` — function
- `get_nodd` — function
- `get_nodd` — function
- `get_noris` — function
- `get_noris` — function
- `get_normal` — function
- `get_normal` — function
- `get_ori` — subroutine
- `get_ori` — subroutine
- `get_pinds` — subroutine
- `get_pinds` — subroutine
- `get_pop_array` — function
- `get_pop_array` — function
- `get_pop_scalar` — function
- `get_pop_scalar` — function
- `get_pops` — subroutine
- `get_pops` — subroutine
- `get_proj` — function
- `get_proj` — function
- `get_sample_ind` — function
- `get_sample_ind` — function
- `get_sampled` — function
- `get_sampled` — function
- `get_state` — function
- `get_state` — function
- `get_state_rec_pops` — subroutine
- `get_state_rec_pops` — subroutine
- `get_state_update_fracs` — subroutine
- `get_state_update_fracs` — subroutine
- `get_static` — subroutine
- `get_static` — subroutine
- `get_str` — function
- `get_str` — function
- `get_sum` — function
- `get_sum` — function
- `get_top` — function
- `get_top` — function
- `get_tseries_neighs` — subroutine
- `get_tseries_neighs` — subroutine
- `get_update_frac` — function
- `get_update_frac` — function
- `get_update_rows` — subroutine
- `get_update_rows` — subroutine
- `get_updatecnt` — function
- `get_updatecnt` — function
- `getter_int` — subroutine
- `getter_int` — subroutine
- `getter_real` — subroutine
- `getter_real` — subroutine
- `getter_str` — subroutine
- `getter_str` — subroutine
- `group_quotas` — subroutine
- `has_been_sampled` — function
- `has_been_sampled` — function
- `has_been_searched` — function
- `has_been_searched` — function
- `included` — function
- `included` — function
- `incr_sampled_updatecnt` — subroutine
- `incr_sampled_updatecnt` — subroutine
- `ini_tseries` — subroutine
- `ini_tseries` — subroutine
- `introd_alig_err` — subroutine
- `introd_alig_err` — subroutine
- `introd_ctf_err` — subroutine
- `introd_ctf_err` — subroutine
- `is_particle` — function
- `is_particle` — function
- `ischar` — function
- `ischar` — function
- `isthere_any` — function
- `isthere_any` — function
- `isthere_by_ind` — function
- `isthere_by_ind` — function
- `keyed_or_intrinsic_partial_shuffle` — subroutine
- `keyed_or_intrinsic_shuffle` — subroutine
- `kill` — subroutine
- `kill` — subroutine
- `kill_chash` — subroutine
- `kill_chash` — subroutine
- `map3dshift22d_all` — subroutine
- `map3dshift22d_all` — subroutine
- `map3dshift22d_by_ind` — subroutine
- `map3dshift22d_by_ind` — subroutine
- `mask_from_state` — subroutine
- `mask_from_state` — subroutine
- `max_ori_strlen_trim` — function
- `max_ori_strlen_trim` — function
- `median_key` — function
- `median_key` — function
- `merge` — subroutine
- `merge` — subroutine
- `merge_classes` — subroutine
- `merge_classes` — subroutine
- `min_euldist` — subroutine
- `min_euldist` — subroutine
- `minmax` — subroutine
- `minmax` — subroutine
- `mirror2d` — subroutine
- `mirror2d` — subroutine
- `mirror3d` — subroutine
- `mirror3d` — subroutine
- `mul_shifts` — subroutine
- `mul_shifts` — subroutine
- `nearest_proj_neighbors_k` — subroutine
- `nearest_proj_neighbors_k` — subroutine
- `nearest_proj_neighbors_mat` — subroutine
- `nearest_proj_neighbors_mat` — subroutine
- `nearest_proj_neighbors_thres` — subroutine
- `nearest_proj_neighbors_thres` — subroutine
- `new_by_n` — subroutine
- `new_by_n` — subroutine
- `new_from_arr` — subroutine
- `new_from_arr` — subroutine
- `order` — function
- `order` — function
- `order_cls` — function
- `order_cls` — function
- `ori2json` — subroutine
- `ori2json` — subroutine
- `ori2prec` — subroutine
- `ori2prec` — subroutine
- `ori2str` — function
- `ori2str` — function
- `overlap` — function
- `overlap` — function
- `partition_eo` — subroutine
- `partition_eo` — subroutine
- `population_blend_weights_counts` — subroutine
- `population_blend_weights_counts` — subroutine
- `population_blend_weights_mass` — subroutine
- `population_blend_weights_mass` — subroutine
- `prec2ori` — subroutine
- `prec2ori` — subroutine
- `print` — subroutine
- `print` — subroutine
- `print_matrices` — subroutine
- `print_matrices` — subroutine
- `proj2class` — subroutine
- `proj2class` — subroutine
- `projs_exist` — function
- `projs_exist` — function
- `read` — subroutine
- `read` — subroutine
- `read_ctfparams_state_eo` — subroutine
- `read_ctfparams_state_eo` — subroutine
- `reallocate` — subroutine
- `reallocate` — subroutine
- `reject` — subroutine
- `reject` — subroutine
- `remap_cls` — subroutine
- `remap_cls` — subroutine
- `remap_projs` — subroutine
- `remap_projs` — subroutine
- `replace_with_closest` — subroutine
- `replace_with_closest` — subroutine
- `reseed_classes` — subroutine
- `reseed_classes` — subroutine
- `revorisgn` — subroutine
- `revorisgn` — subroutine
- `revshsgn` — subroutine
- `revshsgn` — subroutine
- `rnd_cls` — subroutine
- `rnd_cls` — subroutine
- `rnd_corrs` — subroutine
- `rnd_corrs` — subroutine
- `rnd_ctf` — subroutine
- `rnd_ctf` — subroutine
- `rnd_inpls` — subroutine
- `rnd_inpls` — subroutine
- `rnd_lps` — subroutine
- `rnd_lps` — subroutine
- `rnd_ori` — subroutine
- `rnd_ori` — subroutine
- `rnd_oris` — subroutine
- `rnd_oris` — subroutine
- `rnd_proj_space` — subroutine
- `rnd_proj_space` — subroutine
- `rnd_states` — subroutine
- `rnd_states` — subroutine
- `rot_all` — subroutine
- `rot_all` — subroutine
- `rot_by_ind` — subroutine
- `rot_by_ind` — subroutine
- `rot_transp_all` — subroutine
- `rot_transp_all` — subroutine
- `rot_transp_by_ind` — subroutine
- `rot_transp_by_ind` — subroutine
- `round_shifts` — subroutine
- `round_shifts` — subroutine
- `sample4rec` — subroutine
- `sample4rec` — subroutine
- `sample4update_all` — subroutine
- `sample4update_all` — subroutine
- `sample4update_class` — subroutine
- `sample4update_class` — subroutine
- `sample4update_cnt` — subroutine
- `sample4update_cnt` — subroutine
- `sample4update_fillin` — subroutine
- `sample4update_fillin` — subroutine
- `sample4update_missing` — subroutine
- `sample4update_missing` — subroutine
- `sample4update_reprod` — subroutine
- `sample4update_reprod` — subroutine
- `sample4update_rescore` — subroutine
- `sample4update_rescore` — subroutine
- `sample4update_rnd` — subroutine
- `sample4update_rnd` — subroutine
- `sample4update_updated` — subroutine
- `sample4update_updated` — subroutine
- `sample_balanced_frac_best` — subroutine
- `sample_balanced_frac_best` — subroutine
- `sample_balanced_inv` — subroutine
- `sample_balanced_inv` — subroutine
- `sample_balanced_mode` — subroutine
- `sample_balanced_mode` — subroutine
- `sample_balanced_parts` — subroutine
- `sample_balanced_parts` — subroutine
- `sample_ranked_parts` — subroutine
- `sample_ranked_parts` — subroutine
- `set_all2single_char` — subroutine
- `set_all2single_char` — subroutine
- `set_all2single_int` — subroutine
- `set_all2single_int` — subroutine
- `set_all2single_real` — subroutine
- `set_all2single_real` — subroutine
- `set_all_char` — subroutine
- `set_all_char` — subroutine
- `set_all_int` — subroutine
- `set_all_int` — subroutine
- `set_all_real` — subroutine
- `set_all_real` — subroutine
- `set_char` — subroutine
- `set_char` — subroutine
- `set_class` — subroutine
- `set_class` — subroutine
- `set_ctfvars` — subroutine
- `set_ctfvars` — subroutine
- `set_dfx` — subroutine
- `set_dfx` — subroutine
- `set_dfy` — subroutine
- `set_dfy` — subroutine
- `set_dp` — subroutine
- `set_dp` — subroutine
- `set_euler` — subroutine
- `set_euler` — subroutine
- `set_extremal_vars` — subroutine
- `set_extremal_vars` — subroutine
- `set_field2single_char` — subroutine
- `set_field2single_char` — subroutine
- `set_field2single_int` — subroutine
- `set_field2single_int` — subroutine
- `set_field2single_real` — subroutine
- `set_field2single_real` — subroutine
- `set_int` — subroutine
- `set_int` — subroutine
- `set_nonzero_updatecnt` — subroutine
- `set_nonzero_updatecnt` — subroutine
- `set_ogid` — subroutine
- `set_ogid` — subroutine
- `set_ori` — subroutine
- `set_ori` — subroutine
- `set_projs` — subroutine
- `set_projs` — subroutine
- `set_real` — subroutine
- `set_real` — subroutine
- `set_shift` — subroutine
- `set_shift` — subroutine
- `set_state` — subroutine
- `set_state` — subroutine
- `set_stkind` — subroutine
- `set_stkind` — subroutine
- `set_str` — subroutine
- `set_str` — subroutine
- `set_updatecnt` — subroutine
- `set_updatecnt` — subroutine
- `spiral_full` — subroutine
- `spiral_full` — subroutine
- `spiral_lims` — subroutine
- `spiral_lims` — subroutine
- `split_class` — subroutine
- `split_class` — subroutine
- `split_state` — subroutine
- `split_state` — subroutine
- `states_exist` — function
- `states_exist` — function
- `stats_all` — subroutine
- `stats_all` — subroutine
- `stats_masked` — subroutine
- `stats_masked` — subroutine
- `stats_range` — subroutine
- `stats_range` — subroutine
- `str2ori` — subroutine
- `str2ori` — subroutine
- `str2ori_ctfparams_state_eo` — subroutine
- `str2ori_ctfparams_state_eo` — subroutine
- `swape1e3` — subroutine
- `swape1e3` — subroutine
- `symmetrize` — subroutine
- `symmetrize` — subroutine
- `transfer_2Dparams_ori` — subroutine
- `transfer_2Dparams_ori` — subroutine
- `transfer_2Dparams_oris` — subroutine
- `transfer_2Dparams_oris` — subroutine
- `transfer_2Dshifts` — subroutine
- `transfer_2Dshifts` — subroutine
- `transfer_3Dparams_ori` — subroutine
- `transfer_3Dparams_ori` — subroutine
- `transfer_3Dparams_oris` — subroutine
- `transfer_3Dparams_oris` — subroutine
- `transfer_class_assignment` — subroutine
- `transfer_class_assignment` — subroutine
- `transfer_ori` — subroutine
- `transfer_ori` — subroutine
- `unit_groups` — subroutine
- `write2bild` — subroutine
- `write2bild` — subroutine
- `write_by_ind` — subroutine
- `write_by_ind` — subroutine
- `write_projdir_heatmap` — subroutine
- `write_projdir_heatmap` — subroutine
- `write_range` — subroutine
- `write_range` — subroutine
- `zero` — subroutine
- `zero` — subroutine
- `zero_projs` — subroutine
- `zero_projs` — subroutine
- `zero_shifts` — subroutine
- `zero_shifts` — subroutine

Local symbols (enclosing procedure scope):
- `alloc_unit_quotas::takes_remainder` — function
- `spiral_lims::gen_c1` — subroutine
- `write_projdir_heatmap::calc_heatmap` — subroutine
- `write_projdir_heatmap::draw_hline` — subroutine
- `write_projdir_heatmap::draw_smooth_angle_glyph` — subroutine
- `write_projdir_heatmap::draw_smooth_tick_glyph` — subroutine
- `write_projdir_heatmap::draw_smooth_tick_text` — subroutine
- `write_projdir_heatmap::draw_text` — subroutine
- `write_projdir_heatmap::draw_vline` — subroutine
- `write_projdir_heatmap::fill_rect` — subroutine
- `write_projdir_heatmap::glyph_rows` — function
- `write_projdir_heatmap::greek_token_glyph` — function
- `write_projdir_heatmap::greek_token_length` — function
- `write_projdir_heatmap::heat_color` — subroutine
- `write_projdir_heatmap::rgb_code` — function
- `write_projdir_heatmap::segment_distance_sq` — function
- `write_projdir_heatmap::smooth_tick_text_width` — function
- `write_projdir_heatmap::text_width` — function
- `write_projdir_heatmap::token_matches` — function
- `write_projdir_heatmap::upper_character` — function

---
## Module: simple_oris_tester

Files:
- `main/ori/simple_oris_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_forked_process`
- `simple_oris`
- `simple_oris_utils`
- `simple_rnd`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_oris_tests` — subroutine

Private symbols:
- `free_units` — subroutine
- `half_of_unhalved` — type
- `make_units` — subroutine
- `take_unhalved` — subroutine
- `test_assignment` — subroutine
- `test_blend_mass_scripted` — subroutine
- `test_blend_weight_cases` — subroutine
- `test_class_sample_sweep` — subroutine
- `test_cohort_rescore` — subroutine
- `test_compress_and_masks` — subroutine
- `test_constructors_and_basic_props` — subroutine
- `test_empty_partition_sampling` — subroutine
- `test_extract_and_copy` — subroutine
- `test_get_half` — subroutine
- `test_get_pops_halves` — subroutine
- `test_getters_setters` — subroutine
- `test_group_update_counts` — subroutine
- `test_keyed_sampling` — subroutine
- `test_misc_flags` — subroutine
- `test_nested_class_quota` — subroutine
- `test_oridist_from_oris` — subroutine
- `test_proj_space_and_remap` — subroutine
- `test_randomization_and_symmetry` — subroutine
- `test_ranked_balanced_sampling` — subroutine
- `test_reallocate` — subroutine
- `test_reseed_classes` — subroutine
- `test_rnd_oris_bounds` — subroutine
- `test_rotations_and_errors` — subroutine
- `test_sample4rec_global_coverage` — subroutine
- `test_sample4update_cnt_large` — subroutine
- `test_sampling_and_updatecnt` — subroutine
- `test_stats_and_ordering` — subroutine
- `test_write_read_roundtrip` — subroutine
- `unit_counts` — function

Local symbols (enclosing procedure scope):
- `test_blend_mass_scripted::sample` — subroutine
- `test_keyed_sampling::fresh_field` — subroutine

---
## Module: simple_oris_utils

Files:
- `main/ori/simple_oris_utils.f90`

Uses:
- `simple_linalg`
- `simple_oris`

Public symbols:
- `oridist_from_oris` — subroutine

---
## Module: simple_parameters

Files:
- `main/params/simple_parameters.f90`
- `main/params/simple_parameters_core.f90`
- `main/params/simple_parameters_parse.f90`
- `main/params/simple_parameters_phases.f90`

Uses:
- `ieee_arithmetic`
- `simple_atoms`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_decay_funs`
- `simple_gpu_utils`
- `simple_parameters_registry`
- `simple_sp_project`
- `simple_syslib`
- `simple_ui`
- `simple_ui_program`

Public symbols:
- `parameters` — type

Private symbols:
- `bind_input_carg` — subroutine
- `bind_input_carg` — subroutine
- `bind_input_dirs` — subroutine
- `bind_input_dirs` — subroutine
- `bind_input_files` — subroutine
- `bind_input_files` — subroutine
- `bind_input_iarg` — subroutine
- `bind_input_iarg` — subroutine
- `bind_input_rarg` — subroutine
- `bind_input_rarg` — subroutine
- `derive_image_settings` — subroutine
- `derive_image_settings` — subroutine
- `derive_parallel_settings` — subroutine
- `derive_parallel_settings` — subroutine
- `derive_sampling_settings` — subroutine
- `derive_sampling_settings` — subroutine
- `init_dynamic_defaults` — subroutine
- `init_dynamic_defaults` — subroutine
- `is_final_planned_iter` — function
- `is_final_planned_iter` — function
- `new` — subroutine
- `new` — subroutine
- `parse_inputs` — subroutine
- `parse_inputs` — subroutine
- `resolve_parameter_sources` — subroutine
- `resolve_parameter_sources` — subroutine
- `set_img_format` — subroutine
- `set_img_format` — subroutine
- `setup_execution_context` — subroutine
- `setup_execution_context` — subroutine
- `validate_parameter_consistency` — subroutine
- `validate_parameter_consistency` — subroutine

Local symbols (enclosing procedure scope):
- `resolve_parameter_sources::check_file_formats` — subroutine
- `resolve_parameter_sources::check_vol` — subroutine
- `resolve_parameter_sources::double_check_file_formats` — subroutine
- `resolve_parameter_sources::mkfnames` — subroutine
- `resolve_parameter_sources::read_masks` — subroutine
- `resolve_parameter_sources::set_ldim_box_from_stk` — subroutine

---
## Module: simple_parameters_registry

Files:
- `main/params/simple_parameters_registry.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`

Public symbols:
- `param_registry` — type

Private symbols:
- `add_char_fixed` — subroutine
- `add_char_string` — subroutine
- `add_dir` — subroutine
- `add_file` — subroutine
- `add_int` — subroutine
- `add_real` — subroutine
- `add_yesno` — subroutine
- `apply_cmdline` — subroutine
- `char_binding` — type
- `clear` — subroutine
- `dir_binding` — type
- `ensure_char_capacity` — subroutine
- `ensure_dir_capacity` — subroutine
- `ensure_file_capacity` — subroutine
- `ensure_int_capacity` — subroutine
- `ensure_real_capacity` — subroutine
- `ensure_string_capacity` — subroutine
- `file_binding` — type
- `int_binding` — type
- `real_binding` — type
- `string_binding` — type

---
## Module: simple_parameters_tester

Files:
- `main/params/simple_parameters_tester.f90`

Uses:
- `simple_args`
- `simple_cmdline`
- `simple_defs`
- `simple_forked_process`
- `simple_parameters`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_parameters_tests` — subroutine

Private symbols:
- `base_cline` — subroutine
- `build_parameters` — subroutine
- `parameters_child` — type
- `test_accepted_keys` — subroutine
- `test_suite_binding` — subroutine
- `test_ui_keys_bound` — subroutine
- `test_worker_project_recorded` — subroutine
- `test_yes_no_refused` — subroutine

---
## Module: simple_particle_extractor

Files:
- `main/preprocess/simple_particle_extractor.f90`

Uses:
- `simple_core_module_api`
- `simple_eer_factory`
- `simple_image`
- `simple_motion_correct_utils`
- `simple_motion_model`
- `simple_ori`
- `simple_parameters`

Public symbols:
- `ptcl_extractor` — type

Private symbols:
- `add_eer_gain_defects` — subroutine
- `apply_gain_correction` — subroutine
- `cure_movie_outliers` — subroutine
- `cure_outlier_coords` — subroutine
- `cure_outliers` — subroutine
- `cure_outliers_from_coords` — subroutine
- `dealloc_particles_frames` — subroutine
- `display` — subroutine
- `does_exist` — function
- `evaluate_local_shift_1` — subroutine
- `evaluate_local_shift_2` — subroutine
- `extract_all_particles_frames` — subroutine
- `extract_particles` — subroutine
- `extract_particles_from_mic` — subroutine
- `extract_ptcl` — subroutine
- `extract_single_particle_frame` — subroutine
- `from_model` — function
- `from_mov` — function
- `generate_cumulative_particles_stack` — subroutine
- `get_box` — function
- `get_l_neg` — function
- `get_nframes` — function
- `get_particle_mask` — subroutine
- `has_frames` — function
- `has_hotpix_coords` — function
- `has_isoshifts` — function
- `has_model_frameweights` — function
- `has_particle_frames` — function
- `has_weights` — function
- `init_mask` — subroutine
- `init_mic` — subroutine
- `init_model` — subroutine
- `init_mov` — subroutine
- `init_particle_buffers` — subroutine
- `kill` — subroutine
- `load_model_metadata` — subroutine
- `post_process` — subroutine
- `prepare_frames` — subroutine
- `prepare_motion` — subroutine
- `read_movie_frames` — subroutine
- `test_movie_outlier_curation` — subroutine

---
## Module: simple_particle_extractor_tester

Files:
- `main/preprocess/simple_particle_extractor_tester.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_forked_process`
- `simple_image`
- `simple_motion_model`
- `simple_ori`
- `simple_parameters`
- `simple_particle_extractor`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_particle_extractor_tests` — subroutine

Private symbols:
- `assert_dose_weighting_report` — subroutine
- `assert_micrograph_state` — subroutine
- `convention_extract` — type
- `extract_other_convention` — subroutine
- `test_binary_model_constructors` — subroutine
- `test_init_mic_selects_micrograph_mode` — subroutine
- `test_missing_motion_model` — subroutine
- `test_movie_outlier_curation` — subroutine
- `test_particle_frame_cache` — subroutine

Local symbols (enclosing procedure scope):
- `test_binary_model_constructors::check_convention_refused` — subroutine
- `test_particle_frame_cache::check_cached_sum` — subroutine

---
## Module: simple_pcg_atomic_prior

Files:
- `main/volume/simple_pcg_atomic_prior.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_pcg_solver`

Public symbols:
- `pcg_atomic_prior` — type
- `pcg_atomic_prior_stats` — type

Private symbols:
- `amplitudes` — subroutine
- `apply` — subroutine
- `build` — subroutine
- `build_gram` — subroutine
- `estimate_lambda` — function
- `fit_width` — subroutine
- `gather` — function
- `get_bfac` — function
- `get_bfac_min` — function
- `get_box` — function
- `get_natoms` — function
- `get_stats` — function
- `gram_matvec` — subroutine
- `is_built` — function
- `kill` — subroutine
- `model_map` — subroutine
- `new` — subroutine
- `overlap` — function
- `project` — function
- `scatter` — function
- `scatter_add` — subroutine
- `solve_gram` — subroutine
- `uses_fitted_model` — function

Local symbols (enclosing procedure scope):
- `build_gram::cell_of` — function
- `new::drop_duplicates` — subroutine

---
## Module: simple_pcg_atomic_prior_tester

Files:
- `main/volume/simple_pcg_atomic_prior_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_image`
- `simple_oris`
- `simple_pcg_atomic_prior`
- `simple_reconstructor_pcg`
- `simple_test_utils`

Public symbols:
- `run_all_pcg_atomic_prior_tests` — subroutine

Private symbols:
- `bfac_of_sigma` — function
- `build_atoms` — subroutine
- `dotv` — function
- `new_fixture_prior` — subroutine
- `new_fixture_reconstructor` — subroutine
- `new_fixture_reconstructor_prior` — subroutine
- `normv` — function
- `probe_map` — function
- `project_atoms` — subroutine
- `sphere_support` — function
- `test_apply_symmetric` — subroutine
- `test_atom_widths_projector` — subroutine
- `test_column_scaling` — subroutine
- `test_fitted_model_projector` — subroutine
- `test_gram_solve` — subroutine
- `test_projector` — subroutine
- `test_reconstructor_amplitudes` — subroutine
- `test_reconstructor_operator_modes` — subroutine
- `test_reconstructor_zero_strength` — subroutine
- `test_width_fit` — subroutine
- `varied` — function

---
## Module: simple_pcg_halfset_tester

Files:
- `main/volume/simple_pcg_halfset_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_fplane_grid`
- `simple_image`
- `simple_memoize_ft_maps`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_projector`
- `simple_reconstructor`
- `simple_reconstructor_pcg`
- `simple_sp_project`
- `simple_sym`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `build_truth_volume` — subroutine
- `run_all_pcg_halfset_tests` — subroutine

Private symbols:
- `apply_support` — subroutine
- `build_disjoint_halves` — subroutine
- `build_observations` — subroutine
- `calc_fsc` — subroutine
- `corr3` — function
- `forward_residual` — subroutine
- `fsc_low_high` — subroutine
- `interior_minimum` — function
- `reconstruct_half` — subroutine
- `reconstruct_half_gridding` — subroutine
- `reconstruct_half_pair` — subroutine
- `reconstruct_half_trajectory` — subroutine
- `rel_error` — function
- `rel_l2_error` — function
- `test_half_spectrum_weights` — subroutine
- `test_halfset_fsc_96` — subroutine
- `test_halfset_matrix_48` — subroutine
- `test_independent_observations` — subroutine
- `test_truth_volume_fixture` — subroutine

---
## Module: simple_pcg_lattice

Files:
- `main/volume/simple_pcg_lattice.f90`

Uses:
- `iso_c_binding`
- `simple_cartesian_fourier`
- `simple_core_module_api`
- `simple_ftiter`
- `simple_gridding`
- `simple_image`

Public symbols:
- `pcg_lattice` — type
- `plane_loc` — function
- `rim_free` — function
- `scatter_window` — subroutine
- `scatter_window_cmplx_nowrap` — subroutine
- `scatter_window_nowrap` — subroutine
- `scatter_window_pair` — subroutine
- `scatter_window_pair_nowrap` — subroutine
- `scatter_window_real` — subroutine
- `slab_count` — function
- `slab_k_range` — subroutine
- `slab_planes` — subroutine
- `win_below_fold` — function
- `win_wraps` — function

Private symbols:
- `divide_deposition_plane` — subroutine
- `ewald_curv` — function
- `fold_lattice_accum` — function
- `gather_cube` — subroutine
- `gather_plane` — subroutine
- `install_support` — subroutine
- `isqrt` — function
- `kill_lattice` — subroutine
- `mask_mul` — subroutine
- `native_deposition_envelope` — function
- `new_lattice` — subroutine
- `rim_row_start` — function
- `scatter_cube` — subroutine
- `scatter_plane` — subroutine
- `scatter_rim` — subroutine
- `set_window_sphere` — subroutine
- `set_window_volume` — subroutine

---
## Module: simple_pcg_lattice_tester

Files:
- `main/volume/simple_pcg_lattice_tester.f90`

Uses:
- `simple_cartesian_fourier`
- `simple_defs`
- `simple_image`
- `simple_pcg_lattice`
- `simple_test_utils`

Public symbols:
- `run_all_pcg_lattice_tests` — subroutine

Private symbols:
- `lattice_dot` — function
- `plane_dot` — function
- `point_lattice` — function
- `random_lattice` — function
- `random_plane` — function
- `random_rotation` — function
- `test_adjoint` — subroutine
- `test_cube` — subroutine
- `test_curvature_phase` — subroutine
- `test_thread_identity` — subroutine

---
## Module: simple_pcg_solvent_sidecar

Files:
- `main/volume/simple_pcg_solvent_sidecar.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_image_msk`
- `simple_parameters`
- `simple_reconstructor_pcg`
- `simple_refine3d_fnames`

Public symbols:
- `build_solvent_check_support` — subroutine
- `build_solvent_prior_weight` — subroutine
- `estimate_solvent_prior_lambda` — subroutine
- `pcg_solvent_lambda_stats` — type
- `pcg_solvent_stats` — type
- `prepare_solvent_prior_on_pair` — subroutine
- `solvent_prior_cross_half_objective` — function
- `solvent_prior_provenance` — function
- `write_pcg_solvent_pair` — subroutine

Private symbols:
- `build_pcg_solvent_prior_weight` — subroutine

Local symbols (enclosing procedure scope):
- `estimate_solvent_prior_lambda::shrink_half` — subroutine

---
## Module: simple_pcg_solver

Files:
- `main/opt/simple_pcg_solver.f90`

Uses:
- `ieee_arithmetic`
- `simple_core_module_api`

Public symbols:
- `pcg_operator` — type
- `pcg_solve` — subroutine
- `pcg_solver_options` — type
- `pcg_solver_outcome` — type

Private symbols:
- `pcg_apply_i` — subroutine
- `pcg_dot_i` — function
- `pcg_precond_i` — subroutine
- `pcg_size_i` — function
- `pcg_solver_outcome_kill` — subroutine
- `pcg_tic` — function
- `pcg_toc` — function

---
## Module: simple_pcg_sparse_tester

Files:
- `main/volume/simple_pcg_sparse_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_error`
- `simple_fplane_grid`
- `simple_image`
- `simple_memoize_ft_maps`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_reconstructor`
- `simple_reconstructor_pcg`
- `simple_sp_project`
- `simple_sym`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_pcg_sparse_tests` — subroutine

Private symbols:
- `beyond_nyquist_fraction` — function
- `build_atoms` — subroutine
- `build_truth` — subroutine
- `build_view_set` — subroutine
- `centre_spike` — function
- `corr_sphere` — function
- `project_atoms` — subroutine
- `reconstruct_gridding` — subroutine
- `reconstruct_pcg` — subroutine
- `test_atom_fixture` — subroutine
- `test_view_set` — subroutine

---
## Module: simple_pdb2mrc_tester

Files:
- `main/nano/simple_pdb2mrc_tester.f90`

Uses:
- `simple_atoms`
- `simple_fileio`
- `simple_imghead`
- `simple_molecule_data`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_pdb2mrc_tests` — subroutine

Private symbols:
- `check_map` — subroutine
- `test_pdb2mrc_model` — subroutine

---
## Module: simple_persistent_worker_message_base

Files:
- `utils/persistent_worker/message/simple_persistent_worker_message_base.f90`

Public symbols:
- `qsys_persistent_worker_message_base` — type

Private symbols:
- `kill_qsys_persistent_worker_message_base` — subroutine
- `new_qsys_persistent_worker_message_base` — subroutine
- `serialise` — subroutine

---
## Module: simple_persistent_worker_message_heartbeat

Files:
- `utils/persistent_worker/message/simple_persistent_worker_message_heartbeat.f90`

Uses:
- `simple_persistent_worker_message_base`
- `simple_persistent_worker_message_types`

Public symbols:
- `qsys_persistent_worker_message_heartbeat` — type

Private symbols:
- `kill_qsys_persistent_worker_message_heartbeat` — subroutine
- `new_qsys_persistent_worker_message_heartbeat` — subroutine
- `serialise_qsys_persistent_worker_message_heartbeat` — subroutine

---
## Module: simple_persistent_worker_message_status

Files:
- `utils/persistent_worker/message/simple_persistent_worker_message_status.f90`

Uses:
- `simple_defs`
- `simple_persistent_worker_message_base`
- `simple_persistent_worker_message_types`

Public symbols:
- `qsys_persistent_worker_message_status` — type

Private symbols:
- `kill_qsys_persistent_worker_message_status` — subroutine
- `new_qsys_persistent_worker_message_status` — subroutine
- `serialise_qsys_persistent_worker_message_status` — subroutine

---
## Module: simple_persistent_worker_message_task

Files:
- `utils/persistent_worker/message/simple_persistent_worker_message_task.f90`

Uses:
- `simple_defs`
- `simple_persistent_worker_message_base`
- `simple_persistent_worker_message_types`

Public symbols:
- `qsys_persistent_worker_message_task` — type

Private symbols:
- `kill_qsys_persistent_worker_message_task` — subroutine
- `new_qsys_persistent_worker_message_task` — subroutine
- `serialise_qsys_persistent_worker_message_task` — subroutine

---
## Module: simple_persistent_worker_message_terminate

Files:
- `utils/persistent_worker/message/simple_persistent_worker_message_terminate.f90`

Uses:
- `simple_defs`
- `simple_persistent_worker_message_base`
- `simple_persistent_worker_message_types`

Public symbols:
- `qsys_persistent_worker_message_terminate` — type

Private symbols:
- `kill_qsys_persistent_worker_message_terminate` — subroutine
- `new_qsys_persistent_worker_message_terminate` — subroutine
- `serialise_qsys_persistent_worker_message_terminate` — subroutine

---
## Module: simple_persistent_worker_message_tester

Files:
- `utils/persistent_worker/message/simple_persistent_worker_message_tester.f90`

Uses:
- `simple_persistent_worker_message_base`
- `simple_persistent_worker_message_heartbeat`
- `simple_persistent_worker_message_status`
- `simple_persistent_worker_message_task`
- `simple_persistent_worker_message_terminate`
- `simple_persistent_worker_message_types`
- `simple_test_utils`

Public symbols:
- `run_all_persistent_worker_message_tests` — subroutine

Private symbols:
- `test_base_message_new_kill_and_serialise` — subroutine
- `test_heartbeat_new_kill_and_serialise` — subroutine
- `test_message_type_enum_values` — subroutine
- `test_status_new_kill_and_serialise` — subroutine
- `test_task_new_kill_and_serialise` — subroutine
- `test_terminate_new_kill_and_serialise` — subroutine

---
## Module: simple_persistent_worker_message_types

Files:
- `utils/persistent_worker/message/simple_persistent_worker_message_types.f90`

---
## Module: simple_persistent_worker_server

Files:
- `utils/persistent_worker/simple_persistent_worker_server.f90`

Uses:
- `iso_c_binding`
- `simple_core_module_api`
- `simple_ipc_tcp_socket_client`
- `simple_ipc_tcp_socket_helpers`
- `simple_ipc_tcp_socket_server`
- `simple_persistent_worker_message_heartbeat`
- `simple_persistent_worker_message_status`
- `simple_persistent_worker_message_task`
- `simple_persistent_worker_message_terminate`
- `simple_persistent_worker_message_types`
- `simple_string`
- `unix`

Public symbols:
- `forget_inherited_persistent_worker` — subroutine
- `persistent_worker_server` — type

Private symbols:
- `claim_warmup_worker_ids` — subroutine
- `get_host_ips` — function
- `get_n_active_workers` — function
- `get_port` — function
- `get_queue_pressure_workers_required` — function
- `is_running` — function
- `kill` — subroutine
- `mark_worker_slots_launch_pending` — subroutine
- `new` — subroutine
- `persistent_worker_data` — type
- `persistent_worker_runtime` — type
- `queue_task` — function
- `set_warmup_cooldown_enabled` — subroutine
- `worker_listener_thread` — subroutine

Local symbols (enclosing procedure scope):
- `worker_listener_thread::clear_registry_entry_by_fd` — subroutine
- `worker_listener_thread::clear_stale_registry_entry` — subroutine
- `worker_listener_thread::count_active_workers` — function
- `worker_listener_thread::count_effective_workers` — function
- `worker_listener_thread::find_first_dispatchable` — function
- `worker_listener_thread::handle_heartbeat_msg` — subroutine
- `worker_listener_thread::handle_new_task_msg` — subroutine
- `worker_listener_thread::handle_terminate_msg` — subroutine
- `worker_listener_thread::handle_unknown_msg` — subroutine
- `worker_listener_thread::pack_task_queue` — subroutine
- `worker_listener_thread::send_idle_reply` — subroutine
- `worker_listener_thread::set_warmup_request` — subroutine
- `worker_listener_thread::signal_ready` — subroutine
- `worker_listener_thread::update_queue_pressure` — subroutine
- `worker_listener_thread::update_worker_registry` — subroutine

---
## Module: simple_persistent_worker_server_tester

Files:
- `utils/persistent_worker/simple_persistent_worker_server_tester.f90`

Uses:
- `simple_persistent_worker_message_task`
- `simple_persistent_worker_server`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_persistent_worker_server_tests` — subroutine

Private symbols:
- `build_test_task` — subroutine
- `test_claim_warmup_marks_launch_pending_and_clears_requests` — subroutine
- `test_forget_inherited_persistent_worker` — subroutine
- `test_get_host_ips_default_empty` — subroutine
- `test_mark_worker_slots_launch_pending_sets_flags` — subroutine
- `test_module_defaults_and_constants` — subroutine
- `test_new_and_kill_lifecycle` — subroutine
- `test_new_client_only_invalid_address_rejected` — subroutine
- `test_new_client_only_valid_address_parse` — subroutine
- `test_new_invalid_args_do_not_start` — subroutine
- `test_new_is_idempotent_while_running` — subroutine
- `test_queue_task_eventually_rejects_when_full` — subroutine
- `test_queue_task_priority_requests` — subroutine
- `test_server_default_state` — subroutine
- `test_startup_pending_prevents_duplicate_warmup_claim` — subroutine

---
## Module: simple_pftc_api

Files:
- `main/apis/simple_pftc_api.f90`

Uses:
- `gnufor2`
- `simple_class_frcs`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_ctf`
- `simple_fftw3`
- `simple_image`
- `simple_imgarr_utils`
- `simple_parameters`
- `simple_sp_project`

---
## Module: simple_pftc_inplane_tester

Files:
- `main/pftc/simple_pftc_inplane_tester.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_builder`
- `simple_eul_prob_tab_utils`
- `simple_image`
- `simple_matcher_smpl_and_lplims`
- `simple_oris`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_polarft_calc`
- `simple_rnd`
- `simple_strategy2d_srch`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_pftc_inplane_tests` — subroutine

Private symbols:
- `angle_of` — function
- `angular_error` — function
- `build_fixture` — subroutine
- `evaluate` — subroutine
- `fd_check` — subroutine
- `gradient_probe_set` — subroutine
- `grid_identity` — subroutine
- `kill_fixture` — subroutine
- `series_floor_scan` — subroutine
- `test_cc_evaluator` — subroutine
- `test_euclid_band_edges` — subroutine
- `test_euclid_routes_and_recovery` — subroutine
- `test_grid_fit_search` — subroutine
- `test_grid_fit_step` — subroutine
- `test_grid_ptcl_search` — subroutine
- `test_loaded_reference` — subroutine
- `test_particle_bank` — subroutine
- `test_posterior_threads` — subroutine
- `test_search_threads` — subroutine
- `test_shift_frames` — subroutine
- `test_shift_prior` — subroutine
- `test_shift_prior_frame` — subroutine
- `test_shift_prior_grid_fit` — subroutine
- `test_shmat_searched_shells` — subroutine
- `test_strategy2D_route_flags` — subroutine
- `write_phantom` — subroutine

Local symbols (enclosing procedure scope):
- `test_grid_fit_search::grid_fit_case` — subroutine
- `test_grid_fit_search::same_result` — function
- `test_grid_fit_step::fill_quadratic` — subroutine
- `test_grid_ptcl_search::grid_ptcl_case` — subroutine
- `test_particle_bank::direct_reference` — function
- `test_posterior_threads::draws_off_row` — function
- `test_posterior_threads::score_refs` — subroutine
- `test_search_threads::run_searches` — subroutine
- `test_shmat_searched_shells::reset_particle` — subroutine

---
## Module: simple_pftc_shsrch_grad

Files:
- `main/pftc/simple_pftc_shsrch_grad.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_builder`
- `simple_core_module_api`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `fit_quadratic_step` — subroutine
- `grid_stencil_exists` — function
- `pftc_shsrch_grad` — type
- `shift_ptcl2ref` — function
- `shift_ref2ptcl` — function

Private symbols:
- `build_grid_bank` — subroutine
- `build_rot_table` — subroutine
- `clear_shift_prior` — subroutine
- `coarse_search` — subroutine
- `coarse_search_opt_angle` — subroutine
- `does_opt_angle` — function
- `grad_shsrch_costfun` — function
- `grad_shsrch_fdfcostfun` — subroutine
- `grad_shsrch_gcostfun` — subroutine
- `grad_shsrch_get_angle` — function
- `grad_shsrch_get_bank_builds` — function
- `grad_shsrch_get_fit_profile` — subroutine
- `grad_shsrch_get_profile` — subroutine
- `grad_shsrch_kill` — subroutine
- `grad_shsrch_minimize` — function
- `grad_shsrch_minimize_joint` — function
- `grad_shsrch_new_alternating` — subroutine
- `grad_shsrch_new_fixed` — subroutine
- `grad_shsrch_new_joint` — subroutine
- `grad_shsrch_new_mode` — subroutine
- `grad_shsrch_reset_profile` — subroutine
- `grad_shsrch_set_indices` — subroutine
- `grad_shsrch_update_discrete_angle` — subroutine
- `grad_shsrch_update_discrete_angle_wrapper` — subroutine
- `has_shift_prior` — function
- `in_limits` — function
- `joint_cost_to_score` — function
- `joint_grad_at_angle` — subroutine
- `minimize_in_cell` — subroutine
- `objfun_vals` — subroutine
- `prior_centre_angle` — subroutine
- `prior_centre_rot` — function
- `prior_factor` — function
- `prior_factor_ptcl` — function
- `prior_factors` — subroutine
- `select_best_discrete_angle` — subroutine
- `set_limits` — subroutine
- `set_prior_centre` — subroutine
- `set_shift_prior` — subroutine
- `shift_prior_penalty` — function
- `shift_prior_penalty_ptcl` — function
- `try_grid_fit` — subroutine
- `uses_joint_inplane` — function

---
## Module: simple_pftc_srch_api

Files:
- `main/apis/simple_pftc_srch_api.f90`

Uses:
- `simple_class_frcs`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_ctf`
- `simple_image`
- `simple_imgarr_utils`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_sp_project`
- `simple_stack_io`
- `simple_starproject`
- `simple_stat`

---
## Module: simple_pick_strategy

Files:
- `main/strategies/parallelization/simple_pick_strategy.f90`

Uses:
- `simple_binoris_io`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_default_clines`
- `simple_gui_utils`
- `simple_image`
- `simple_image_msk`
- `simple_imghead`
- `simple_parameters`
- `simple_picker_iter`
- `simple_qsys_env`
- `simple_sp_project`
- `simple_stack_io`
- `simple_syslib`

Public symbols:
- `create_pick_strategy` — function
- `make_pickrefs_impl` — subroutine
- `pick_distr_strategy` — type
- `pick_inmem_strategy` — type
- `pick_strategy` — type
- `template_lowpass` — function

Private symbols:
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_end_message` — function
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `endmsg_interface` — function
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_end_message` — function
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `set_pick_defaults` — subroutine
- `validate_pick_mode` — subroutine

---
## Module: simple_picker_iter

Files:
- `main/pick/simple_picker_iter.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_image`
- `simple_parameters`
- `simple_picker_utils`

Public symbols:
- `picker_iter` — type

Private symbols:
- `iterate` — subroutine
- `kill` — subroutine
- `read_pickrefs` — subroutine

---
## Module: simple_picker_utils

Files:
- `main/pick/simple_picker_utils.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_imgarr_utils`
- `simple_parameters`
- `simple_pickref`
- `simple_picksegdiam`
- `simple_string_utils`

Public symbols:
- `exec_gaupick` — subroutine
- `exec_refpick` — subroutine
- `exec_segdiampick` — subroutine

---
## Module: simple_pickref

Files:
- `main/pick/simple_pickref.f90`

Uses:
- `simple_core_module_api`
- `simple_gui_utils`
- `simple_image`
- `simple_micproc`
- `simple_pickref_corr_batch`
- `simple_segmentation`

Public symbols:
- `multiref_merge` — subroutine
- `pickref` — type
- `read_mic_raw_pickref` — subroutine

Private symbols:
- `detect_peaks` — subroutine
- `distance_filter` — subroutine
- `distance_filter_grid` — subroutine
- `get_box` — function
- `get_maxdiam` — function
- `get_nboxes` — function
- `get_positions` — subroutine
- `get_scores` — subroutine
- `kill` — subroutine
- `limit_nboxes` — subroutine
- `match_boximgs` — subroutine
- `new` — subroutine
- `peak_vs_nonpeak_stats` — subroutine
- `refine_upscaled` — subroutine
- `refine_upscaled_batch` — subroutine
- `refpick` — subroutine
- `report_boxfile` — subroutine
- `report_thumb_den` — subroutine
- `setup_iterators` — subroutine
- `write_boxfile` — subroutine

---
## Module: simple_pickref_corr_batch

Files:
- `main/pick/simple_pickref_corr_batch.f90`

Uses:
- `simple_defs`
- `simple_image`
- `simple_linalg`
- `simple_timer_omp`

Public symbols:
- `pickref_corr_batch` — type

Private symbols:
- `kill` — subroutine
- `new` — subroutine
- `report_stats` — subroutine
- `score_positions` — subroutine

---
## Module: simple_picksegdiam

Files:
- `main/pick/simple_picksegdiam.f90`

Uses:
- `simple_core_module_api`
- `simple_gui_utils`
- `simple_image`
- `simple_image_bin`
- `simple_linked_list`
- `simple_micproc`
- `simple_nrtxtfile`

Public symbols:
- `picksegdiam` — type

Private symbols:
- `get_diameters` — subroutine
- `get_nboxes` — function
- `kill` — subroutine
- `pick_1` — subroutine
- `pick_2` — subroutine
- `pick_3` — subroutine
- `write_diameters` — subroutine
- `write_pos_and_diams` — subroutine

Local symbols (enclosing procedure scope):
- `pick_1::dealloc_objs` — subroutine
- `pick_3::dealloc_objs` — subroutine

---
## Module: simple_polar_matching_bench

Files:
- `main/commanders/test/simple_polar_matching_bench.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `simple_builder`
- `simple_cmdline`
- `simple_ctf`
- `simple_defs`
- `simple_error`
- `simple_eul_prob_tab`
- `simple_eul_prob_tab2d`
- `simple_eul_prob_tab_neigh`
- `simple_eul_prob_tab_utils`
- `simple_fftw3`
- `simple_map_reduce`
- `simple_math_ft`
- `simple_oris`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_rnd`
- `simple_sp_project`
- `simple_test_utils`

Public symbols:
- `compare_fill_results` — subroutine
- `kill_fill_result` — subroutine
- `pmb_fill_result` — type
- `pmb_route_name` — function
- `polar_matching_bench` — type

Private symbols:
- `cgemm` — subroutine
- `get_gmat_build_time` — function
- `get_parts` — function
- `get_sizes` — subroutine
- `get_vecs` — subroutine
- `kill` — subroutine
- `kill_bufs` — subroutine
- `make_blocks` — subroutine
- `make_bufs` — subroutine
- `memo_gbytes` — function
- `new` — subroutine
- `pair_crvec` — subroutine
- `pair_finish` — subroutine
- `pair_sums` — subroutine
- `pmb_buf` — type
- `prod_kernels` — subroutine
- `ref_bytes_per_pair` — function
- `run_fill` — subroutine
- `score` — function
- `score_blocked` — function
- `score_matrix` — function
- `score_per_pair` — function
- `score_replica` — function
- `score_states_inner` — function

Local symbols (enclosing procedure scope):
- `new::build_fill_toolbox` — subroutine
- `new::fill_band` — subroutine
- `new::memoize_ptcls_here` — subroutine
- `new::memoize_ref` — subroutine
- `run_fill::alloc_result` — subroutine
- `run_fill::dense_projs_ns` — function
- `run_fill::gather_buffers` — subroutine
- `run_fill::gather_dense_batch3D` — subroutine
- `run_fill::gather_tab2D` — subroutine
- `run_fill::gather_tab3D` — subroutine
- `run_fill::iteration_shift_prior` — function
- `run_fill::perturb_committed` — subroutine

---
## Module: simple_polar_shiftmap

Files:
- `main/pftc/simple_polar_shiftmap.f90`

Uses:
- `simple_defs`
- `simple_error`

Public symbols:
- `polar_shiftmap` — type
- `shiftmap_ndirections` — function

Private symbols:
- `cell_areas` — subroutine
- `coeffs` — subroutine
- `cross_spectrum` — subroutine
- `eval` — function
- `eval_rotcoeffs` — function
- `exists` — function
- `get_alpha` — function
- `get_cell_area` — function
- `get_kfromto` — function
- `get_log_area` — function
- `get_nalpha` — function
- `get_nharm` — function
- `get_nharm_max` — function
- `get_nrho` — function
- `get_rho` — function
- `get_shift` — function
- `harmonic_order` — function
- `kill` — subroutine
- `map` — subroutine
- `matches` — function
- `new` — subroutine
- `peak` — subroutine
- `shell_weights` — subroutine

---
## Module: simple_polar_shiftmap_tester

Files:
- `main/pftc/simple_polar_shiftmap_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_eul_prob_tab_utils`
- `simple_parameters`
- `simple_polar_shiftmap`
- `simple_polarft_calc`
- `simple_rnd`
- `simple_test_utils`
- `simple_timer`

Public symbols:
- `run_all_polar_shiftmap_tests` — subroutine

Private symbols:
- `check_map_objective` — subroutine
- `direct_shift_width` — function
- `kweighted_sqsum` — function
- `magnitudes` — subroutine
- `make_pftc` — subroutine
- `make_pftc_kernel` — subroutine
- `rotate_half_ring` — subroutine
- `shift_args` — subroutine
- `test_against_pftc_cc` — subroutine
- `test_angular_nyquist` — subroutine
- `test_cell_areas` — subroutine
- `test_closed_form_ramp` — subroutine
- `test_disk_within_radius` — subroutine
- `test_kernel_draw` — subroutine
- `test_kernel_identity` — subroutine
- `test_kernel_threads` — subroutine
- `test_kernel_timing` — subroutine
- `test_map_cache_lifecycle` — subroutine
- `test_orders_beyond_ring` — subroutine
- `test_peak_recovery` — subroutine
- `test_physical_width` — subroutine
- `test_screen_width_identity` — subroutine
- `test_zero_prior_width_rounding` — subroutine

---
## Module: simple_polarft_calc

Files:
- `main/pftc/simple_polarft_access.f90`
- `main/pftc/simple_polarft_calc.f90`
- `main/pftc/simple_polarft_core.f90`
- `main/pftc/simple_polarft_corr.f90`
- `main/pftc/simple_polarft_ctf.f90`
- `main/pftc/simple_polarft_geom.f90`
- `main/pftc/simple_polarft_memo.f90`
- `main/pftc/simple_polarft_ops_io.f90`
- `main/pftc/simple_polarft_shiftmap.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_eul_prob_tab_utils`
- `simple_gpu_utils`
- `simple_math`
- `simple_parameters`
- `simple_pftc_api`
- `simple_polar_shiftmap`
- `simple_projector`
- `simple_projector_pft_batch`
- `simple_rnd`
- `simple_type_defs`

Public symbols:
- `polarft_calc` — type
- `vol_pad2ref_pfts` — subroutine
- `vol_pad2ref_pfts` — subroutine
- `vol_pad2ref_pfts_opt` — subroutine
- `vol_pad2ref_pfts_opt` — subroutine

Private symbols:
- `alloc_memo_ptcls` — subroutine
- `alloc_memo_ptcls` — subroutine
- `alloc_memo_refs` — subroutine
- `alloc_memo_refs` — subroutine
- `allocate_memo_workspace` — subroutine
- `allocate_memo_workspace` — subroutine
- `allocate_pft` — function
- `allocate_pft` — function
- `allocate_ptcl_pft` — function
- `allocate_ptcl_pft` — function
- `assert_polarize_pdim` — subroutine
- `assign_sigma2_noise` — subroutine
- `assign_sigma2_noise` — subroutine
- `bump_ptcl_gen` — subroutine
- `bump_ptcl_gen` — subroutine
- `calc_frc` — subroutine
- `calc_frc` — subroutine
- `cc_shells_sum` — subroutine
- `cc_sums` — subroutine
- `check_euclid_weights` — subroutine
- `check_euclid_weights` — subroutine
- `check_ptcl_memo` — subroutine
- `check_ptcl_memo` — subroutine
- `check_ref_aln_8` — subroutine
- `cp_even2odd_ref` — subroutine
- `cp_even2odd_ref` — subroutine
- `cp_even_ref2ptcl` — subroutine
- `cp_even_ref2ptcl` — subroutine
- `cp_odd2even_ref` — subroutine
- `cp_odd2even_ref` — subroutine
- `cp_refs` — subroutine
- `cp_refs` — subroutine
- `create_polar_absctfmats` — subroutine
- `create_polar_absctfmats` — subroutine
- `draw_posterior_crvec` — subroutine
- `euclid_dists_from_crvec` — subroutine
- `euclid_shell_sum` — subroutine
- `euclid_sums` — subroutine
- `eval_joint_coeffs_at_rotind` — subroutine
- `eval_series_at_rotind` — subroutine
- `exists` — function
- `exists` — function
- `fftw_cmat` — type
- `fftw_crmat` — type
- `fftw_crvec` — type
- `fftw_drvec` — type
- `gen_best_euclid_val` — subroutine
- `gen_best_euclid_val` — subroutine
- `gen_best_euclid_val_bank` — subroutine
- `gen_best_euclid_val_bank` — subroutine
- `gen_best_euclid_val_cross` — subroutine
- `gen_best_euclid_val_cross` — subroutine
- `gen_best_objfun_val` — subroutine
- `gen_best_objfun_val` — subroutine
- `gen_corr_cc_for_rot_8_noshift` — function
- `gen_corr_cc_for_rot_8_noshift` — function
- `gen_corr_cc_for_rot_8_shifted` — function
- `gen_corr_cc_for_rot_8_shifted` — function
- `gen_corr_cc_grad_for_rot_8` — subroutine
- `gen_corr_cc_grad_for_rot_8` — subroutine
- `gen_corr_for_rot_8_aln` — function
- `gen_corr_for_rot_8_aln` — function
- `gen_corr_for_rot_8_noshift` — function
- `gen_corr_for_rot_8_noshift` — function
- `gen_corr_for_rot_8_shifted` — function
- `gen_corr_for_rot_8_shifted` — function
- `gen_corr_grad_at_angle` — subroutine
- `gen_corr_grad_at_angle` — subroutine
- `gen_corr_grad_for_rot_8_aln` — subroutine
- `gen_corr_grad_for_rot_8_aln` — subroutine
- `gen_corrs` — subroutine
- `gen_corrs` — subroutine
- `gen_corrs_bank` — subroutine
- `gen_euclid_crvec` — subroutine
- `gen_euclid_crvec_bank` — subroutine
- `gen_euclid_dists` — subroutine
- `gen_euclid_dists` — subroutine
- `gen_euclid_dists_bank` — subroutine
- `gen_euclid_dists_bank` — subroutine
- `gen_euclid_for_rot_8_noshift` — function
- `gen_euclid_for_rot_8_noshift` — function
- `gen_euclid_for_rot_8_shifted` — function
- `gen_euclid_for_rot_8_shifted` — function
- `gen_euclid_grad_for_rot_8` — subroutine
- `gen_euclid_grad_for_rot_8` — subroutine
- `gen_euclid_residual_grad` — subroutine
- `gen_euclid_val_at` — subroutine
- `gen_euclid_val_at` — subroutine
- `gen_euclid_val_at_bank` — subroutine
- `gen_euclid_val_at_bank` — subroutine
- `gen_euclids` — subroutine
- `gen_euclids` — subroutine
- `gen_many2many_euclids_cufft` — subroutine
- `gen_many2many_euclids_cufft` — subroutine
- `gen_many2many_euclids_cufft_kernel` — subroutine
- `gen_many_euclids` — subroutine
- `gen_many_euclids` — subroutine
- `gen_many_euclids_gpu` — subroutine
- `gen_many_euclids_gpu` — subroutine
- `gen_many_objfun_vals` — subroutine
- `gen_many_objfun_vals` — subroutine
- `gen_normalized_corr_grad_at_angle` — subroutine
- `gen_objfun_vals` — subroutine
- `gen_objfun_vals` — subroutine
- `gen_objfun_vals_bank` — subroutine
- `gen_objfun_vals_bank` — subroutine
- `gen_posterior_val` — subroutine
- `gen_posterior_val` — subroutine
- `gen_posterior_val_bank` — subroutine
- `gen_posterior_val_bank` — subroutine
- `gen_prob_euclid_val` — subroutine
- `gen_prob_euclid_val` — subroutine
- `gen_prob_objfun_val` — subroutine
- `gen_prob_objfun_val` — subroutine
- `gen_prob_power_euclid_val` — subroutine
- `gen_prob_power_euclid_val` — subroutine
- `gen_prob_power_objfun_val` — subroutine
- `gen_prob_power_objfun_val` — subroutine
- `gen_raw_euclid_grad_at_angle` — subroutine
- `gen_raw_euclid_grad_at_angle` — subroutine
- `gen_raw_euclid_grad_for_rot_8` — subroutine
- `gen_raw_euclid_grad_for_rot_8` — subroutine
- `gen_raw_euclid_vals` — subroutine
- `gen_raw_euclid_vals` — subroutine
- `gen_raw_euclid_vals_impl` — subroutine
- `gen_shiftmap` — subroutine
- `gen_shiftmap` — subroutine
- `gen_shiftmap_kernel` — subroutine
- `gen_shmat` — subroutine
- `gen_shmat` — subroutine
- `gen_shmat4aln` — subroutine
- `gen_shmat4aln` — subroutine
- `gen_shmat4aln_8` — subroutine
- `gen_shmat4aln_8` — subroutine
- `gen_shmat_8` — subroutine
- `gen_shmat_8` — subroutine
- `gen_sigma_contrib` — subroutine
- `gen_sigma_contrib` — subroutine
- `get_coord` — function
- `get_coord` — function
- `get_dang` — function
- `get_dang` — function
- `get_kfromto` — function
- `get_kfromto` — function
- `get_nptcls` — function
- `get_nptcls` — function
- `get_nrefs` — function
- `get_nrefs` — function
- `get_nrots` — function
- `get_nrots` — function
- `get_pdim_interp` — function
- `get_pdim_interp` — function
- `get_pdim_srch` — function
- `get_pdim_srch` — function
- `get_pftsz` — function
- `get_pftsz` — function
- `get_pinds` — subroutine
- `get_pinds` — subroutine
- `get_polar_coords` — subroutine
- `get_polar_coords` — subroutine
- `get_precalc_objfun_vals` — subroutine
- `get_precalc_objfun_vals` — subroutine
- `get_ptcl_line` — subroutine
- `get_ptcl_line` — subroutine
- `get_ptcl_pft` — subroutine
- `get_ptcl_pft` — subroutine
- `get_ref_pft` — subroutine
- `get_ref_pft` — subroutine
- `get_roind` — function
- `get_roind` — function
- `get_roind_fast` — function
- `get_roind_fast` — function
- `get_rot` — function
- `get_rot` — function
- `get_shiftmap_profile` — subroutine
- `get_shiftmap_profile` — subroutine
- `has_shiftmap` — function
- `has_shiftmap` — function
- `heap_vars` — type
- `is_cc_objfun` — function
- `is_cc_objfun` — function
- `is_euclid_objfun` — function
- `is_euclid_objfun` — function
- `is_joint_grad_objfun` — function
- `is_joint_grad_objfun` — function
- `is_with_ctf` — function
- `is_with_ctf` — function
- `kill` — subroutine
- `kill` — subroutine
- `kill_memo_ptcls` — subroutine
- `kill_memo_ptcls` — subroutine
- `kill_memo_refs` — subroutine
- `kill_memo_refs` — subroutine
- `kill_memo_workspace` — subroutine
- `kill_memo_workspace` — subroutine
- `kill_one_shiftmap` — subroutine
- `kill_one_shiftmap` — subroutine
- `kill_shiftmap` — subroutine
- `kill_shiftmap` — subroutine
- `load_ref_aln_8` — subroutine
- `load_ref_aln_8` — subroutine
- `memoize_ptcls` — subroutine
- `memoize_ptcls` — subroutine
- `memoize_refs` — subroutine
- `memoize_refs` — subroutine
- `memoize_shifted_ptcl` — subroutine
- `memoize_shifted_ptcl` — subroutine
- `memoize_sigma2_weights` — subroutine
- `memoize_sigma2_weights` — subroutine
- `memoize_sqsum_ptcl` — subroutine
- `memoize_sqsum_ptcl` — subroutine
- `mirror_ref_pft_by_ind` — subroutine
- `mirror_ref_pft_by_ind` — subroutine
- `mirror_ref_pft_by_pft` — subroutine
- `mirror_ref_pft_by_pft` — subroutine
- `new` — subroutine
- `new` — subroutine
- `new_shiftmap` — subroutine
- `new_shiftmap` — subroutine
- `open_pft_or_ctf2_array_for_write` — subroutine
- `polarize_ptcl_pft` — subroutine
- `polarize_ptcl_pft` — subroutine
- `polarize_ref_pft` — subroutine
- `polarize_ref_pft` — subroutine
- `ptcl_iseven` — function
- `ptcl_iseven` — function
- `ptcl_posterior_temp` — function
- `ptcl_posterior_temp` — function
- `read_ref_pfts` — subroutine
- `read_ref_pfts` — subroutine
- `reallocate_ptcls` — subroutine
- `reallocate_ptcls` — subroutine
- `ref_slots_identical` — function
- `ref_slots_identical` — function
- `rotate_ctf` — subroutine
- `rotate_ctf` — subroutine
- `rotate_ptcl` — subroutine
- `rotate_ptcl` — subroutine
- `rotate_ref_8` — subroutine
- `rotate_ref_8` — subroutine
- `set_eo` — subroutine
- `set_eo` — subroutine
- `set_ptcl_ctfmat` — subroutine
- `set_ptcl_ctfmat` — subroutine
- `set_ptcl_fcomp` — subroutine
- `set_ptcl_fcomp` — subroutine
- `set_ptcl_pft` — subroutine
- `set_ptcl_pft` — subroutine
- `set_ref_fcomp` — subroutine
- `set_ref_fcomp` — subroutine
- `set_ref_pft` — subroutine
- `set_ref_pft` — subroutine
- `set_rnd_scope` — subroutine
- `set_rnd_scope` — subroutine
- `set_with_ctf` — subroutine
- `set_with_ctf` — subroutine
- `setup_shiftmaps` — subroutine
- `setup_shiftmaps` — subroutine
- `shift_ptcl` — subroutine
- `shift_ptcl` — subroutine
- `shiftmap_cache` — type
- `shiftmap_cell_area` — function
- `shiftmap_cell_area` — function
- `shiftmap_cross` — function
- `shiftmap_cross` — function
- `shiftmap_dims` — subroutine
- `shiftmap_dims` — subroutine
- `shiftmap_draw` — subroutine
- `shiftmap_draw` — subroutine
- `shiftmap_mass` — function
- `shiftmap_mode` — subroutine
- `shiftmap_mode` — subroutine
- `shiftmap_objfun_grid` — subroutine
- `shiftmap_objfun_grid` — subroutine
- `shiftmap_peak` — subroutine
- `shiftmap_peak` — subroutine
- `shiftmap_prior` — subroutine
- `shiftmap_prior` — subroutine
- `shiftmap_shift` — function
- `shiftmap_shift` — function
- `shiftmap_work` — type
- `write_ptcl_pft_range` — subroutine
- `write_ptcl_pft_range` — subroutine
- `write_ref_pfts` — subroutine
- `write_ref_pfts` — subroutine
- `wsqsum_ptcl_band` — function
- `wsqsum_ptcl_band` — function
- `wsqsum_ptcl_now` — function
- `wsqsum_ptcl_now` — function

Local symbols (enclosing procedure scope):
- `alloc_memo_ptcls::memo_ptcls_alloc_error` — subroutine
- `setup_shiftmaps::build_map` — subroutine

---
## Module: simple_polarft_corr_tester

Files:
- `main/pftc/simple_polarft_corr_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_eul_prob_tab_utils`
- `simple_forked_process`
- `simple_image`
- `simple_parameters`
- `simple_pftc_shsrch_grad`
- `simple_polarft_calc`
- `simple_rnd`
- `simple_test_utils`

Public symbols:
- `run_all_polarft_corr_tests` — subroutine

Private symbols:
- `order_violation` — type
- `score_out_of_order` — subroutine
- `test_bank_and_at_pose_kernels` — subroutine
- `test_call_order_refusals` — subroutine
- `test_cross_half_euclid` — subroutine
- `test_generated_images` — subroutine
- `test_posterior_inpl_draw` — subroutine
- `test_ref_pfts_offset_io` — subroutine

Local symbols (enclosing procedure scope):
- `test_generated_images::check_frc_identity` — subroutine
- `test_generated_images::check_frc_peak` — subroutine
- `test_generated_images::cyclic_dist` — function
- `test_posterior_inpl_draw::draws_against_posterior` — subroutine

---
## Module: simple_preprocess_strategy

Files:
- `main/strategies/parallelization/simple_preprocess_strategy.f90`

Uses:
- `simple_binoris_io`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_ctf_estimate_iter`
- `simple_epu_xml`
- `simple_mini_stream_utils`
- `simple_motion_correct_iter`
- `simple_motion_correct_utils`
- `simple_motion_gain_helpers`
- `simple_parameters`
- `simple_qsys_env`
- `simple_sp_project`

Public symbols:
- `create_preprocess_strategy` — function
- `preprocess_distr_strategy` — type
- `preprocess_inmem_strategy` — type
- `preprocess_strategy` — type

Private symbols:
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_end_message` — function
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `endmsg_interface` — function
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_end_message` — function
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `set_preprocess_defaults` — subroutine

---
## Module: simple_private_exec_api

Files:
- `main/apis/simple_private_exec_api.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_cavgs`
- `simple_commanders_checks`
- `simple_commanders_denoise`
- `simple_commanders_distr`
- `simple_commanders_euclid`
- `simple_commanders_euclid_distr`
- `simple_commanders_flex_pca`
- `simple_commanders_imgops`
- `simple_commanders_mask`
- `simple_commanders_misc`
- `simple_commanders_motion`
- `simple_commanders_ori`
- `simple_commanders_pick`
- `simple_commanders_preprocess`
- `simple_commanders_project_cls`
- `simple_commanders_project_core`
- `simple_commanders_project_ptcl`
- `simple_commanders_rec`
- `simple_commanders_rec_distr`
- `simple_commanders_refine2d`
- `simple_commanders_refine3d`
- `simple_commanders_volops`
- `simple_core_module_api`
- `simple_jiffys`
- `simple_private_prgs`
- `simple_symanalyzer`
- `simple_syslib`
- `simple_ui`
- `single_commanders_trajectory`
- `single_commanders_tseries`

---
## Module: simple_private_prgs

Files:
- `main/ui/simple_private_prgs.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `get_n_private_keys_required` — function
- `get_private_keys_required` — function
- `make_private_ui` — subroutine
- `print_private_cmdline` — subroutine

Private symbols:
- `add_private_prg` — subroutine
- `init_cmd_dict` — subroutine
- `new_private_prgs` — subroutine
- `push_opt` — subroutine
- `push_opt_key` — subroutine
- `push_req` — subroutine
- `push_req_key` — subroutine
- `set_name` — subroutine
- `simple_private_prg` — type

---
## Module: simple_procimgstk

Files:
- `main/image_processing/simple_procimgstk.f90`

Uses:
- `simple_bspline_smoother`
- `simple_core_module_api`
- `simple_ctf`
- `simple_image`
- `simple_memoize_ft_maps`
- `simple_sp_project`
- `simple_stack_io`

Public symbols:
- `add_noise_imgfile` — subroutine
- `apply_ctf_imgfile` — subroutine
- `bp_imgfile` — subroutine
- `bs_smoother_imgfile` — subroutine
- `clip_imgfile` — subroutine
- `copy_imgfile` — subroutine
- `icm_imgfile` — subroutine
- `loc_sdev_imgfile` — subroutine
- `mask_imgfile` — subroutine
- `mirror_imgfile` — subroutine
- `neg_imgfile` — subroutine
- `nlmean_imgfile` — subroutine
- `noise_imgfile` — subroutine
- `noise_norm_imgfile` — subroutine
- `norm_imgfile` — subroutine
- `pad_imgfile` — subroutine
- `phase_rand_imgfile` — subroutine
- `quantize_imgfile` — subroutine
- `raise_exception_imgfile` — subroutine
- `random_cls_from_imgfile` — subroutine
- `random_selection_from_imgfile` — subroutine
- `real_filter_imgfile` — subroutine
- `roavg_imgfile` — subroutine
- `scale_and_clip_imgfile` — subroutine
- `scale_imgfile` — subroutine
- `selection_from_tseries_imgfile` — subroutine
- `shift_imgfile` — subroutine
- `taper_edges_imgfile` — subroutine

---
## Module: simple_progress

Files:
- `utils/simple_progress.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `lastfoundfile_update` — subroutine
- `progress_estimate_2D` — function
- `progress_estimate_preprocess_stream` — function
- `progressfile_complete_parts` — subroutine
- `progressfile_init` — subroutine
- `progressfile_init_part` — subroutine
- `progressfile_init_parts` — subroutine
- `progressfile_update` — subroutine
- `progressfile_update_part` — subroutine

---
## Module: simple_project_merge_tester

Files:
- `main/project/simple_project_merge_tester.f90`

Uses:
- `simple_fileio`
- `simple_image`
- `simple_oris`
- `simple_projfile_utils`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_project_merge_tests` — subroutine

Private symbols:
- `cleanup_files` — subroutine
- `create_empty_file` — subroutine
- `fill_class_row` — subroutine
- `fill_particle_row` — subroutine
- `make_project` — subroutine
- `test_fix_projfile_refuses_particles_without_stacks` — subroutine
- `test_fix_projfile_repairs_unpruned_legacy_stack` — subroutine
- `test_fix_projfile_repairs_wrong_segment_stkind` — subroutine
- `test_merge_pruned_stack_indexing_heterogeneous_ctf` — subroutine
- `test_remap_project_paths` — subroutine
- `test_remap_project_paths_allows_unmatched_scope` — subroutine
- `test_remap_project_paths_scoped_roots` — subroutine
- `test_sigma2_state_path_registration` — subroutine

---
## Module: simple_project_superset

Files:
- `main/project/simple_project_superset.f90`

Uses:
- `simple_core_module_api`
- `simple_sp_project`

Public symbols:
- `project_superset` — type

Private symbols:
- `get_masks` — subroutine
- `get_ncohort` — function
- `get_nfrozen` — function
- `get_nfrozen_state` — function
- `get_nnever_updated` — function
- `get_nretired` — function
- `image_id` — function
- `is_frozen_row` — function
- `is_small_cohort` — function
- `kill` — subroutine
- `mask` — subroutine
- `new` — subroutine
- `restore` — subroutine
- `retire_from_frozen` — subroutine
- `same_ctf` — function
- `same_image` — function
- `validate_cohort_states` — subroutine
- `validate_identity` — subroutine

Local symbols (enclosing procedure scope):
- `validate_identity::check_appended_stacks` — subroutine
- `validate_identity::name_particle` — subroutine

---
## Module: simple_project_superset_tester

Files:
- `main/project/simple_project_superset_tester.f90`

Uses:
- `simple_defs`
- `simple_project_superset`
- `simple_sp_project`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_project_superset_tests` — subroutine

Private symbols:
- `expect_refusal` — subroutine
- `make_current` — subroutine
- `make_frozen` — subroutine
- `names_particle` — function
- `test_carried_membership` — subroutine
- `test_identity_accepts_superset` — subroutine
- `test_identity_refusals` — subroutine
- `test_mask_and_restore` — subroutine
- `test_membership` — subroutine

---
## Module: simple_projector

Files:
- `main/image/simple_projector.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `projector` — type

Private symbols:
- `expand_cmat` — subroutine
- `fproject` — subroutine
- `fproject_serial` — subroutine
- `get_kbwin` — function
- `interp_fcomp` — function
- `interp_fcomp_oversamp` — function
- `is_expanded` — function
- `kill_expanded` — subroutine
- `reset_expanded` — subroutine

---
## Module: simple_projector_pft_batch

Files:
- `main/image/simple_projector_pft_batch.f90`

Uses:
- `simple_core_module_api`
- `simple_kbinterpol`
- `simple_polarft_calc`
- `simple_projector`

Public symbols:
- `fproject_polar_batch` — subroutine
- `fproject_polar_batch_mirr` — subroutine
- `fproject_polar_batch_opt` — subroutine

Local symbols (enclosing procedure scope):
- `fproject_polar_batch_opt::fproject_polar_batch_opt_kernel` — subroutine

---
## Module: simple_projfile_utils

Files:
- `fileio/simple_projfile_utils.f90`

Uses:
- `iso_fortran_env`
- `simple_class_frcs`
- `simple_core_module_api`
- `simple_image`
- `simple_oris`
- `simple_ptcl_layout`
- `simple_sigma2_state`
- `simple_sigma2_state_file`
- `simple_sp_project`

Public symbols:
- `absolutize_project_stack_paths` — subroutine
- `fix_project_file` — subroutine
- `merge_chunk_projfiles` — subroutine
- `merge_chunk_projfiles_without_sigma2` — subroutine
- `merge_selected_project_files` — subroutine
- `remap_project_paths` — subroutine

Local symbols (enclosing procedure scope):
- `fix_project_file::err` — subroutine
- `fix_project_file::fix_particle_indstk` — subroutine
- `fix_project_file::fix_stack_nptcls_stk` — subroutine
- `fix_project_file::repair` — subroutine
- `fix_project_file::repair_stack_counts_without_particles` — subroutine
- `fix_project_file::repair_stack_ranges` — subroutine
- `fix_project_file::set_stack_int_if_changed` — subroutine
- `fix_project_file::stkind_for_project_row` — function
- `fix_project_file::warn` — subroutine
- `merge_chunk_projfiles::concatenate_canonical_states` — subroutine
- `merge_chunk_projfiles::validate_chunk_canonical_identity` — subroutine
- `merge_selected_project_files::copy_particle_row` — subroutine
- `merge_selected_project_files::count_stack_particles` — function
- `merge_selected_project_files::ctf_enabled_for_row` — function
- `merge_selected_project_files::data_segments_present` — function
- `merge_selected_project_files::make_prefix_offsets` — subroutine
- `merge_selected_project_files::max_ogid_in_oris` — function
- `merge_selected_project_files::merge_selected_canonical_field` — subroutine
- `merge_selected_project_files::merge_selected_canonical_states` — subroutine
- `merge_selected_project_files::remap_row_ogid` — subroutine
- `merge_selected_project_files::require_row_field` — subroutine
- `merge_selected_project_files::resolved_particle_indstk` — function
- `merge_selected_project_files::source_max_ogid` — function
- `merge_selected_project_files::validate_field_presence` — subroutine
- `merge_selected_project_files::validate_mic_sampling` — subroutine
- `merge_selected_project_files::validate_particle_field` — subroutine
- `merge_selected_project_files::validate_selected_source_sigma` — subroutine
- `merge_selected_project_files::validate_source_project` — subroutine
- `merge_selected_project_files::validate_stack_dimensions` — subroutine
- `merge_selected_project_files::validate_stack_particle_ranges` — subroutine
- `remap_project_paths::is_filesystem_root` — function
- `remap_project_paths::normalize_root` — function
- `remap_project_paths::normalize_separators` — function
- `remap_project_paths::process_cavg_fields` — subroutine
- `remap_project_paths::process_field` — subroutine
- `remap_project_paths::process_mic_fields` — subroutine
- `remap_project_paths::process_ptcl_fields` — subroutine
- `remap_project_paths::process_scope` — subroutine
- `remap_project_paths::process_vol_fields` — subroutine
- `remap_project_paths::remap_one_path` — subroutine

---
## Module: simple_pspec_thumb_iter

Files:
- `main/preprocess/simple_pspec_thumb_iter.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_image`
- `simple_parameters`

Public symbols:
- `pspec_thumb_iter` — type

Private symbols:
- `iterate` — subroutine

---
## Module: simple_pspecs

Files:
- `main/strategies/search/simple_pspecs.f90`

Uses:
- `simple_core_module_api`
- `simple_fsc`
- `simple_image`
- `simple_oris`

Public symbols:
- `pspecs` — type

Private symbols:
- `calc_distmat_1` — function
- `calc_distmat_2` — function
- `get_nspecs` — function
- `kill` — subroutine
- `new` — subroutine

---
## Module: simple_ptcl_layout

Files:
- `main/project/simple_ptcl_layout.f90`

Uses:
- `iso_fortran_env`
- `simple_oris`
- `simple_sigma2_state_file`
- `simple_sp_project`
- `simple_string`

Private symbols:
- `ptcl_layout_digest_arrays` — function
- `ptcl_layout_digest_project` — function

---
## Module: simple_ptcl_sieve

Files:
- `main/sieve/simple_ptcl_sieve.f90`

Uses:
- `simple_cavg_quality_helpers`
- `simple_cavg_quality_model`
- `simple_cavg_quality_selection`
- `simple_cavg_quality_types`
- `simple_class_compatibility`
- `simple_cmdline`
- `simple_defs`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_error`
- `simple_fileio`
- `simple_gui_utils`
- `simple_image`
- `simple_image_bin`
- `simple_imgarr_utils`
- `simple_imghead`
- `simple_optics_maps`
- `simple_parameters`
- `simple_projfile_utils`
- `simple_qsys_async_job`
- `simple_qsys_env`
- `simple_qsys_job_record`
- `simple_rec_list`
- `simple_segmentation`
- `simple_sigma2_state_file`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_timer`

Public symbols:
- `ptcl_sieve` — type
- `ptcl_sieve_settings` — type
- `read_chunked_mics` — subroutine
- `rebuild_chunked_mics` — subroutine
- `sieve_lpstart` — function
- `sieve_settings` — function

Private symbols:
- `append_chunk_coarse` — subroutine
- `append_chunk_fine` — subroutine
- `blanks_for_underscores` — function
- `cancel` — subroutine
- `carry_chunk_lists` — subroutine
- `chunk2D_coarse_defaults` — type
- `chunk2D_fine_defaults` — type
- `chunk2D_state` — type
- `chunk_folder_ids` — function
- `chunk_has_cavgs` — function
- `chunk_job_failed` — function
- `cleanup_chunk` — subroutine
- `collect_and_reject` — subroutine
- `combine_completed_chunks` — subroutine
- `cycle` — subroutine
- `fail_chunk` — subroutine
- `generate_chunk_coarse_cline` — subroutine
- `generate_chunk_fine_cline` — subroutine
- `generate_chunks_coarse` — subroutine
- `generate_chunks_fine` — subroutine
- `get_finished` — function
- `get_latest` — function
- `get_n_accepted_micrographs` — function
- `get_n_accepted_ptcls` — function
- `get_n_chunks_coarse` — function
- `get_n_chunks_fine` — function
- `get_n_chunks_running` — function
- `get_n_coarse_accepted_ptcls` — function
- `get_n_coarse_rejected_ptcls` — function
- `get_n_failed_chunks` — function
- `get_n_failed_ptcls` — function
- `get_n_fine_accepted_ptcls` — function
- `get_n_pass_1_non_rejected_ptcls` — function
- `get_n_pass_2_non_rejected_ptcls` — function
- `get_n_rejected_ptcls` — function
- `get_n_total_particles` — function
- `hand_off` — subroutine
- `hand_off_final_set` — subroutine
- `import_existing_chunks_coarse` — subroutine
- `import_existing_chunks_fine` — subroutine
- `kill` — subroutine
- `merge_and_clear` — subroutine
- `new` — subroutine
- `read_consumed_coarse` — subroutine
- `reason_token` — function
- `reconcile_consumed_coarse` — subroutine
- `reject_cavgs` — subroutine
- `remove_duplicates` — function
- `renew_unfinished_chunk` — subroutine
- `restore_counters` — subroutine
- `set_final_ingestion` — subroutine
- `submit` — subroutine
- `timer_start` — function
- `timer_stop` — subroutine
- `underscores_for_blanks` — function
- `unset_final_ingestion` — subroutine
- `write_chunk_mics` — subroutine
- `write_chunked_mics` — subroutine
- `write_consumed_coarse` — subroutine
- `write_rejection_reason_overlay_jpg` — subroutine

Local symbols (enclosing procedure scope):
- `cleanup_chunk::iter_jpeg_index` — function
- `cleanup_chunk::iter_stack_index` — function
- `cleanup_chunk::sigma_file_rank` — function
- `generate_chunks_fine::all_done_but` — function
- `rebuild_chunked_mics::chunk_paths` — subroutine
- `restore_counters::count_chunk` — subroutine
- `submit::prepare_chunk_attempt` — subroutine
- `write_rejection_reason_overlay_jpg::draw_reason_key_swatches` — subroutine
- `write_rejection_reason_overlay_jpg::paint_mask_outline` — subroutine
- `write_rejection_reason_overlay_jpg::paint_reason_border` — subroutine
- `write_rejection_reason_overlay_jpg::reason_code_from_text` — function
- `write_rejection_reason_overlay_jpg::reason_color` — function
- `write_rejection_reason_overlay_jpg::rgb_code` — function

---
## Module: simple_ptcl_sieve_tester

Files:
- `main/sieve/simple_ptcl_sieve_tester.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_fname`
- `simple_image`
- `simple_optics_maps`
- `simple_parameters`
- `simple_projfile_utils`
- `simple_ptcl_sieve`
- `simple_rec_list`
- `simple_sp_project`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_ptcl_sieve_tests` — subroutine

Private symbols:
- `init_test_params` — subroutine
- `make_chunk_project` — subroutine
- `make_completed_coarse_chunk` — subroutine
- `setup_workspace` — subroutine
- `teardown_workspace` — subroutine
- `test_collect_and_reject_hard_gates` — subroutine
- `test_cycle_empty_project_list` — subroutine
- `test_failed_chunk_dropped_after_retry` — subroutine
- `test_failed_chunk_retried_once` — subroutine
- `test_finished_semantics` — subroutine
- `test_hand_off_applies_optics_map` — subroutine
- `test_import_existing_chunks_and_counts` — subroutine
- `test_import_restores_final_ingestion_chunk` — subroutine
- `test_lpstart_given_or_derived` — subroutine
- `test_merge_staged_chunk_only` — subroutine
- `test_new_accepts_tuning_overrides` — subroutine
- `test_new_kill_and_empty_queries` — subroutine
- `test_rebuild_chunked_mics` — subroutine
- `test_restart_reconciles_and_counts` — subroutine
- `test_settings_from_params` — subroutine
- `test_single_pass_ignores_incomplete_fine` — subroutine
- `write_text` — subroutine

---
## Module: simple_ptcl_sieve_utils

Files:
- `main/sieve/simple_ptcl_sieve_utils.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_parameters`
- `simple_sp_project`
- `simple_stream_api`
- `simple_string`
- `simple_string_utils`
- `simple_type_defs`

Public symbols:
- `generate_sieve_projects` — subroutine

---
## Module: simple_ptcl_subset

Files:
- `main/project/simple_ptcl_subset.f90`

Uses:
- `ieee_arithmetic`
- `simple_eo_partition`
- `simple_error`
- `simple_fileio`
- `simple_ori`
- `simple_oris`
- `simple_sp_project`
- `simple_string`
- `simple_sym`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `select_ptcl_subset` — subroutine
- `write_ptcl_subset` — subroutine

---
## Module: simple_ptcl_subset_tester

Files:
- `main/project/simple_ptcl_subset_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_eo_partition`
- `simple_image`
- `simple_imghead`
- `simple_oris`
- `simple_parameters`
- `simple_ptcl_subset`
- `simple_sp_project`
- `simple_string`
- `simple_sym`
- `simple_syslib`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_ptcl_subset_tests` — subroutine

Private symbols:
- `cleanup` — subroutine
- `fixture` — subroutine
- `test_export` — subroutine
- `test_selection` — subroutine

---
## Module: simple_qsys_async_job

Files:
- `utils/qsys/simple_qsys_async_job.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_defs_fname`
- `simple_fileio`
- `simple_qsys_env`
- `simple_qsys_job_record`
- `simple_string`
- `simple_syslib`

Public symbols:
- `qsys_async_job` — type

Private symbols:
- `cancel` — subroutine
- `get_dir` — function
- `get_log` — function
- `kill` — subroutine
- `start` — subroutine
- `start_seq` — subroutine
- `status` — function

---
## Module: simple_qsys_base

Files:
- `utils/qsys/simple_qsys_base.f90`

Uses:
- `simple_chash`
- `simple_string`

Public symbols:
- `qsys_base` — type

Private symbols:
- `generic_kill` — subroutine
- `generic_new` — subroutine
- `generic_submit_cmd` — function
- `generic_write_array_instr` — subroutine
- `generic_write_instr` — subroutine

---
## Module: simple_qsys_coarray

Files:
- `utils/qsys/simple_qsys_coarray.f90`

Uses:
- `simple_core_module_api`
- `simple_qsys_base`

Public symbols:
- `qsys_coarray` — type

Private symbols:
- `get_coarray_submit_cmd` — function
- `kill_coarray_env` — subroutine
- `new_coarray_env` — subroutine
- `write_coarray_array_header` — subroutine
- `write_coarray_header` — subroutine

---
## Module: simple_qsys_ctrl

Files:
- `utils/qsys/simple_qsys_ctrl.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_mem_estimator`
- `simple_parameters`
- `simple_persistent_worker_message_task`
- `simple_persistent_worker_server`
- `simple_qsys_base`
- `simple_qsys_coarray`
- `simple_qsys_job_record`
- `simple_qsys_local`
- `simple_qsys_lsf`
- `simple_qsys_persistent_worker`
- `simple_qsys_slurm`
- `simple_rnd`
- `simple_syslib`
- `unix`

Public symbols:
- `qsys_ctrl` — type

Private symbols:
- `add_to_stream_stack` — subroutine
- `add_to_streaming` — subroutine
- `augment_partition_job_descr` — subroutine
- `cancel_streaming` — subroutine
- `clear_job_status` — subroutine
- `clear_partition_job_descr` — subroutine
- `clear_stack` — subroutine
- `constructor` — function
- `declare_coarray_jobs_finished` — subroutine
- `dispatch_task_to_persistent_worker` — subroutine
- `exists` — function
- `failed_parts` — function
- `free_all_cunits` — subroutine
- `generate_array_script` — subroutine
- `generate_script_cline` — subroutine
- `generate_script_job` — subroutine
- `generate_script_jobs` — subroutine
- `generate_script_part` — subroutine
- `get_done_stacksz` — function
- `get_failed_stacksz` — function
- `get_jobs_status` — subroutine
- `get_ncomputing_units_avail` — function
- `get_stack_range` — function
- `get_stacksz` — function
- `get_stream_done_stack` — subroutine
- `get_stream_fail_stack` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `poll_pause` — subroutine
- `prep_part_jobs` — subroutine
- `print_jobs_status` — subroutine
- `report_phase_completion` — subroutine
- `schedule_array_jobs` — subroutine
- `schedule_jobs` — subroutine
- `schedule_streaming` — subroutine
- `service_persistent_worker_warmup` — subroutine
- `set_jobs_status` — subroutine
- `shell_quote` — function
- `stage_coarray_jobs` — subroutine
- `stop_on_failed_parts` — subroutine
- `submit_coarray_jobs` — subroutine
- `submit_script` — subroutine
- `submit_scripts` — subroutine
- `update_queue` — subroutine
- `worker_warmup_callback` — subroutine
- `write_job_record` — subroutine
- `write_seed_export` — subroutine

Local symbols (enclosing procedure scope):
- `schedule_streaming::updatestack` — subroutine

---
## Module: simple_qsys_ctrl_tester

Files:
- `utils/qsys/simple_qsys_ctrl_tester.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_qsys_ctrl`
- `simple_qsys_job_record`
- `simple_qsys_local`
- `simple_test_utils`

Public symbols:
- `run_all_qsys_ctrl_tests` — subroutine

Private symbols:
- `file_text` — function
- `make_ctrl` — subroutine
- `test_async_job_record_and_fresh_dir` — subroutine
- `test_cancel_unfinished_jobs` — subroutine
- `test_computing_units` — subroutine
- `test_constructor_and_kill` — subroutine
- `test_failed_local_parts` — subroutine
- `test_failed_parts` — subroutine
- `test_job_liveness` — subroutine
- `test_jobs_status_round_trip` — subroutine
- `test_multi_job_script` — subroutine
- `test_prep_part_jobs` — subroutine
- `test_single_job_scripts` — subroutine
- `test_streaming_stack` — subroutine
- `write_lines` — subroutine

---
## Module: simple_qsys_env

Files:
- `utils/qsys/simple_qsys_env.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_environment`
- `simple_parameters`
- `simple_persistent_worker_server`
- `simple_qsys_base`
- `simple_qsys_coarray`
- `simple_qsys_ctrl`
- `simple_qsys_factory`
- `simple_qsys_funs`
- `simple_qsys_local`
- `simple_qsys_lsf`
- `simple_qsys_pbs`
- `simple_qsys_persistent_worker`
- `simple_qsys_sge`
- `simple_qsys_slurm`
- `simple_sp_project`

Public symbols:
- `forget_warmup_envs` — subroutine
- `is_warmup_owner` — function
- `n_warmup_envs` — function
- `qsys_env` — type
- `register_warmup_env` — subroutine
- `unregister_warmup_env` — subroutine

Private symbols:
- `diagnose_parallelism` — subroutine
- `exec_simple_prg_in_queue` — subroutine
- `exec_simple_prg_in_queue_async` — subroutine
- `exec_simple_prgs_in_queue_async` — subroutine
- `exists` — function
- `gen_script` — subroutine
- `gen_scripts_and_schedule_jobs` — subroutine
- `get_exec_bin` — function
- `get_n_active_persistent_workers` — function
- `get_navail_computing_units` — function
- `get_persistent_worker_nthr` — function
- `get_persistent_worker_server_address` — function
- `get_qsys` — function
- `init_qdescr_from_runtime` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `persistent_worker_warmup_callback` — subroutine
- `qsys_env_ref` — type
- `service_persistent_worker_warmup` — subroutine
- `start_persistent_workers` — subroutine

Local symbols (enclosing procedure scope):
- `new::check_worker_claim` — subroutine
- `new::standard_exec_path` — subroutine
- `new::validate_requested_qsys` — subroutine

---
## Module: simple_qsys_env_tester

Files:
- `utils/qsys/simple_qsys_env_tester.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_parameters`
- `simple_qsys_env`
- `simple_sp_project`
- `simple_test_utils`

Public symbols:
- `run_all_qsys_env_tests` — subroutine

Private symbols:
- `active_counts` — function
- `test_active_balanced_qsys_parts` — subroutine
- `test_active_balanced_split` — subroutine
- `test_forget_warmup_envs` — subroutine
- `test_installation_path_policy` — subroutine
- `test_missing_project_scheduler` — subroutine
- `test_warmup_owner_lifetime` — subroutine

---
## Module: simple_qsys_factory

Files:
- `utils/qsys/simple_qsys_factory.f90`

Uses:
- `simple_core_module_api`
- `simple_qsys_base`
- `simple_qsys_coarray`
- `simple_qsys_local`
- `simple_qsys_lsf`
- `simple_qsys_pbs`
- `simple_qsys_persistent_worker`
- `simple_qsys_sge`
- `simple_qsys_slurm`

Public symbols:
- `qsys_factory` — type

Private symbols:
- `kill` — subroutine
- `new` — subroutine

---
## Module: simple_qsys_funs

Files:
- `utils/qsys/simple_qsys_funs.f90`

Uses:
- `simple_core_module_api`
- `simple_parameters`

Public symbols:
- `qsys_cleanup` — subroutine
- `qsys_declare_part_finished` — subroutine
- `qsys_watcher_diag` — subroutine
- `qsys_watcher_file` — subroutine
- `qsys_watcher_files` — subroutine

---
## Module: simple_qsys_job_record

Files:
- `utils/qsys/simple_qsys_job_record.f90`

Uses:
- `simple_defs`
- `simple_defs_fname`
- `simple_fileio`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_timer`

Public symbols:
- `cancel_queued_job` — function
- `cancel_unfinished_jobs` — subroutine
- `check_job_lost` — function
- `fresh_job_dir` — subroutine
- `job_left_unfinished` — function
- `query_job` — function
- `set_aside_dir` — subroutine

---
## Module: simple_qsys_local

Files:
- `utils/qsys/simple_qsys_local.f90`

Uses:
- `simple_core_module_api`
- `simple_qsys_base`

Public symbols:
- `qsys_local` — type

Private symbols:
- `get_local_submit_cmd` — function
- `kill_local_env` — subroutine
- `new_local_env` — subroutine
- `write_local_array_header` — subroutine
- `write_local_header` — subroutine

---
## Module: simple_qsys_lsf

Files:
- `utils/qsys/simple_qsys_lsf.f90`

Uses:
- `simple_core_module_api`
- `simple_qsys_base`

Public symbols:
- `qsys_lsf` — type

Private symbols:
- `get_lsf_submit_cmd` — function
- `kill_lsf_env` — subroutine
- `new_lsf_env` — subroutine
- `write_lsf_array_header` — subroutine
- `write_lsf_header` — subroutine

---
## Module: simple_qsys_pbs

Files:
- `utils/qsys/simple_qsys_pbs.f90`

Uses:
- `simple_core_module_api`
- `simple_qsys_base`

Public symbols:
- `qsys_pbs` — type

Private symbols:
- `get_pbs_submit_cmd` — function
- `kill_pbs_env` — subroutine
- `new_pbs_env` — subroutine
- `write_pbs_array_header` — subroutine
- `write_pbs_header` — subroutine

Local symbols (enclosing procedure scope):
- `write_pbs_header::write_formatted` — subroutine

---
## Module: simple_qsys_persistent_worker

Files:
- `utils/qsys/simple_qsys_persistent_worker.f90`

Uses:
- `simple_core_module_api`
- `simple_persistent_worker_server`
- `simple_qsys_base`

Public symbols:
- `qsys_persistent_worker` — type

Private symbols:
- `get_worker_submit_cmd` — function
- `kill_worker_env` — subroutine
- `new_worker_env` — subroutine
- `write_worker_array_header` — subroutine
- `write_worker_header` — subroutine

---
## Module: simple_qsys_sge

Files:
- `utils/qsys/simple_qsys_sge.f90`

Uses:
- `simple_core_module_api`
- `simple_qsys_base`

Public symbols:
- `qsys_sge` — type

Private symbols:
- `get_sge_submit_cmd` — function
- `kill_sge_env` — subroutine
- `new_sge_env` — subroutine
- `write_sge_array_header` — subroutine
- `write_sge_header` — subroutine

---
## Module: simple_qsys_slurm

Files:
- `utils/qsys/simple_qsys_slurm.f90`

Uses:
- `simple_core_module_api`
- `simple_qsys_base`

Public symbols:
- `qsys_slurm` — type

Private symbols:
- `get_slurm_submit_cmd` — function
- `kill_slurm_env` — subroutine
- `new_slurm_env` — subroutine
- `write_slurm_array_header` — subroutine
- `write_slurm_header` — subroutine

---
## Module: simple_quality_logistic

Files:
- `main/quality/simple_quality_logistic.f90`

Uses:
- `simple_error`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `dataset_learn_role` — function
- `logistic_dataset_class_weights` — subroutine
- `quality_logistic_problem` — type
- `robust_macro_score` — function
- `stable_sigmoid_8` — function

Private symbols:
- `fdf` — subroutine
- `finalize` — subroutine
- `fit` — subroutine
- `kill` — subroutine
- `logistic_cost_wrapper` — function
- `logistic_fdf_wrapper` — subroutine
- `logistic_gradient_wrapper` — subroutine
- `new` — subroutine
- `set_row` — subroutine

---
## Module: simple_quality_logistic_tester

Files:
- `main/quality/simple_quality_logistic_tester.f90`

Uses:
- `simple_cavg_quality_stats`
- `simple_defs`
- `simple_quality_logistic`
- `simple_test_utils`

Public symbols:
- `run_all_quality_logistic_tests` — subroutine

Private symbols:
- `test_auc_eps_ties` — subroutine
- `test_auc_matches_pairwise` — subroutine
- `test_auc_trajectory_scale` — subroutine
- `test_fit_repeatable` — subroutine
- `test_roles_weights_macro` — subroutine
- `test_toy_plane_recovered` — subroutine
- `toy_problem` — subroutine

---
## Module: simple_ran_tabu

Files:
- `utils/math/simple_ran_tabu.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_rnd`

Public symbols:
- `ran_tabu` — type

Private symbols:
- `balanced` — subroutine
- `constructor` — function
- `insert` — subroutine
- `irnd` — function
- `irnd_gau` — function
- `irnd_pair` — subroutine
- `is` — function
- `kill` — subroutine
- `mnomal` — function
- `ne_mnomal_iarr` — subroutine
- `ne_ran_iarr` — subroutine
- `remove` — subroutine
- `reset` — subroutine
- `shuffle_int` — subroutine
- `shuffle_real` — subroutine

---
## Module: simple_rec3D_pcg_strategy

Files:
- `main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`

Uses:
- `ieee_arithmetic`
- `simple_atoms`
- `simple_bench_io`
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_eul_prob_tab_utils`
- `simple_halfmap_diagnostics`
- `simple_image`
- `simple_image_msk`
- `simple_imgarr_utils`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_nu_state_filter`
- `simple_oris`
- `simple_parameters`
- `simple_pcg_atomic_prior`
- `simple_pcg_solvent_sidecar`
- `simple_reconstructor_pcg`
- `simple_refine3d_fnames`
- `simple_refine3d_provenance`
- `simple_sigma2_files`
- `simple_state_weight_set`
- `simple_trail_chain_manifest`

Public symbols:
- `accumulate_pcg_extras` — subroutine
- `accumulate_pcg_state_half` — subroutine
- `calculate_pcg_state_diagnostics` — subroutine
- `execute_rec3D_pcg_distributed_master` — subroutine
- `execute_rec3D_pcg_worker` — subroutine
- `pcg_assemble_bench` — type
- `pcg_extra_entries` — subroutine
- `rec3D_master_nthr` — function
- `set_pcg_solve_support` — subroutine
- `validate_pcg_common` — subroutine
- `write_pcg_assemble_bench` — subroutine

Private symbols:
- `build_pcg_density_support` — subroutine
- `build_pcg_state_support` — subroutine
- `matched_band_kstop` — function
- `pcg_chain_provenance` — function
- `pcg_density_support_enabled` — function
- `pcg_raw_provenance` — function
- `pcg_trail_seed_requested` — function
- `report_beyond_band_excess` — subroutine
- `resolve_pcg_fsc_txt_fname` — function
- `write_output_diagnostics` — subroutine

Local symbols (enclosing procedure scope):
- `build_pcg_state_support::handle_missing_reference` — subroutine
- `execute_rec3D_pcg_distributed_master::add_raw_read_size` — subroutine
- `execute_rec3D_pcg_distributed_master::atomic_prior_provenance` — function
- `execute_rec3D_pcg_distributed_master::chain_manifest_valid` — function
- `execute_rec3D_pcg_distributed_master::check_solvent_lambda_by_resolve_distributed` — subroutine
- `execute_rec3D_pcg_distributed_master::column_of` — function
- `execute_rec3D_pcg_distributed_master::count_first_time` — function
- `execute_rec3D_pcg_distributed_master::count_full_state_half` — function
- `execute_rec3D_pcg_distributed_master::discard_stale_trail_chain_pair` — subroutine
- `execute_rec3D_pcg_distributed_master::distributed_half_job` — type
- `execute_rec3D_pcg_distributed_master::field_value` — function
- `execute_rec3D_pcg_distributed_master::find_empty_states` — subroutine
- `execute_rec3D_pcg_distributed_master::finish_distributed_half_job` — subroutine
- `execute_rec3D_pcg_distributed_master::fit_amplitudes_pair` — subroutine
- `execute_rec3D_pcg_distributed_master::full_state_half` — subroutine
- `execute_rec3D_pcg_distributed_master::nfields` — function
- `execute_rec3D_pcg_distributed_master::prepare_atomic_prior_on_pair` — subroutine
- `execute_rec3D_pcg_distributed_master::prepare_distributed_half_job` — subroutine
- `execute_rec3D_pcg_distributed_master::prepare_regularized_half_job` — subroutine
- `execute_rec3D_pcg_distributed_master::previous_state_maps_exist` — function
- `execute_rec3D_pcg_distributed_master::read_species_prior_fields` — subroutine
- `execute_rec3D_pcg_distributed_master::reduce_solve_state_pair` — subroutine
- `execute_rec3D_pcg_distributed_master::set_chain_blend_weights` — subroutine
- `execute_rec3D_pcg_distributed_master::solve_distributed_half_pair` — subroutine
- `execute_rec3D_pcg_distributed_master::solve_prepared_half_job` — subroutine
- `execute_rec3D_pcg_distributed_master::write_atomic_prior_products` — subroutine
- `execute_rec3D_pcg_distributed_master::write_chain_half` — subroutine
- `execute_rec3D_pcg_distributed_master::write_chain_manifest` — subroutine
- `execute_rec3D_pcg_distributed_master::write_distributed_diagnostics` — subroutine
- `execute_rec3D_pcg_worker::accumulate_worker_state_half` — subroutine
- `execute_rec3D_pcg_worker::collect_weighted_state_half` — subroutine

---
## Module: simple_rec3D_service

Files:
- `main/strategies/parallelization/simple_rec3D_service.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_commanders_volops`
- `simple_core_module_api`
- `simple_matcher_3drec`
- `simple_parameters`
- `simple_qsys_env`
- `simple_rec3d_pcg_strategy`
- `simple_refine3d_fnames`
- `simple_sigma2_files`
- `simple_state_weight_set`
- `simple_volassemble_strategy`

Public symbols:
- `rec3D_backend_id` — function
- `rec3D_backend_is_wired` — function
- `rec3D_request` — type
- `rec3D_service` — type

Private symbols:
- `apply_request` — subroutine
- `assemble_gridding` — subroutine
- `assemble_pcg` — subroutine
- `execute` — subroutine
- `execute_inproc` — subroutine
- `execute_queue` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `postprocess_states` — subroutine
- `rec3D_saved_params` — type
- `register_outputs` — subroutine
- `remove_pcg_raw_files` — subroutine
- `rename_outputs` — subroutine
- `report_weight_set` — subroutine
- `restore_params` — subroutine
- `write_pinds` — subroutine

---
## Module: simple_rec3D_service_tester

Files:
- `main/strategies/parallelization/simple_rec3D_service_tester.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_commanders_rec`
- `simple_defs`
- `simple_eo_partition`
- `simple_error`
- `simple_fileio`
- `simple_flex_pca_application_tester`
- `simple_halfmap_diagnostics`
- `simple_image`
- `simple_parameters`
- `simple_rec3d_service`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_truth_metrics`
- `simple_test_utils`
- `simple_trail_chain_manifest`

Public symbols:
- `run_all_rec3D_service_lib_tests` — subroutine

Private symbols:
- `conformer` — function
- `fractional_update_then_new_generation` — subroutine
- `fractional_weights` — subroutine
- `make_project` — subroutine
- `publish_weights` — subroutine
- `run_rec3D` — subroutine
- `state_map_fname` — function
- `test_empty_state_carried` — subroutine
- `test_fractional_phantom` — subroutine
- `test_labels_equal_hard` — subroutine
- `test_request_settings` — subroutine
- `test_trailing_weight_identity` — subroutine
- `test_unmasked_gridding_export` — subroutine
- `write_mixture` — subroutine

---
## Module: simple_rec3D_strategy

Files:
- `main/strategies/parallelization/simple_rec3D_strategy.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_eo_partition`
- `simple_exec_helpers`
- `simple_parameters`
- `simple_qsys_funs`
- `simple_rec3d_service`
- `simple_sigma2_bootstrap`

Public symbols:
- `create_rec3D_strategy` — function
- `rec3D_distr_strategy` — type
- `rec3D_inmem_strategy` — type
- `rec3D_pcg_inmem_strategy` — type
- `rec3D_strategy` — type
- `run_rec3D` — subroutine

Private symbols:
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_execute` — subroutine
- `distr_initialize` — subroutine
- `ensure_canonical_sigma_state` — subroutine
- `exec_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_execute` — subroutine
- `inmem_initialize` — subroutine
- `pcg_inmem_execute` — subroutine
- `sample_rows` — subroutine
- `sync_resolved_rec_params` — subroutine

---
## Module: simple_rec3D_strategy_tester

Files:
- `main/strategies/parallelization/simple_rec3D_strategy_tester.f90`

Uses:
- `simple_cmdline`
- `simple_eul_prob_tab_utils`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_rec3d_pcg_strategy`
- `simple_rec3d_service`
- `simple_rec3d_strategy`
- `simple_reconstructor_pcg`
- `simple_rnd`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_rec3D_strategy_tests` — subroutine

Private symbols:
- `strategy_kind` — function
- `test_backend_resolution` — subroutine
- `test_defaults` — subroutine
- `test_factory_branches` — subroutine
- `test_pcg_extra_entries` — subroutine

Local symbols (enclosing procedure scope):
- `test_pcg_extra_entries::accumulate_entries` — subroutine
- `test_pcg_extra_entries::accumulate_plain` — subroutine

---
## Module: simple_rec3D_wave_strategy

Files:
- `main/strategies/parallelization/simple_rec3D_wave_strategy.f90`

Uses:
- `simple_builder`
- `simple_defs`
- `simple_eo_partition`
- `simple_error`
- `simple_fileio`
- `simple_halfmap_diagnostics`
- `simple_image`
- `simple_imgarr_utils`
- `simple_matcher_ptcl_io`
- `simple_math_ctf`
- `simple_nlcg_solver`
- `simple_oris`
- `simple_parameters`
- `simple_rec3d_pcg_strategy`
- `simple_reconstructor_pcg`
- `simple_refine3d_fnames`
- `simple_sigma2_files`
- `simple_string`
- `simple_string_utils`
- `simple_timer`
- `simple_type_defs`
- `simple_wave_imaging`
- `simple_wave_projector`

Public symbols:
- `exec_rec3D_wave` — subroutine
- `rec3D_wave_objective` — type

Private symbols:
- `objective_curvature` — function
- `objective_dot` — function
- `objective_eval` — subroutine
- `objective_fitted_gain` — function
- `objective_kill` — subroutine
- `objective_line_losses` — subroutine
- `objective_new` — subroutine
- `objective_precond` — subroutine
- `objective_project` — subroutine
- `objective_size` — function
- `raw_density` — type

Local symbols (enclosing procedure scope):
- `exec_rec3D_wave::base_solve` — subroutine
- `exec_rec3D_wave::regularize_start` — subroutine
- `exec_rec3D_wave::state_scale` — subroutine
- `exec_rec3D_wave::write_wave_diagnostics` — subroutine

---
## Module: simple_rec3D_wave_strategy_tester

Files:
- `main/strategies/parallelization/simple_rec3D_wave_strategy_tester.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_rec`
- `simple_commanders_sim`
- `simple_defs`
- `simple_fileio`
- `simple_image`
- `simple_ori`
- `simple_oris`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_rec3D_wave_lib_tests` — subroutine

Private symbols:
- `band_fsc` — function
- `make_data` — subroutine
- `reconstruct` — subroutine

---
## Module: simple_rec_list

Files:
- `utils/structs/simple_rec_list.f90`

Uses:
- `simple_core_module_api`
- `simple_linked_list`

Public symbols:
- `chunk_rec` — type
- `process_rec` — type
- `project_rec` — type
- `rec` — type
- `rec_iterator` — type
- `rec_list` — type

Private symbols:
- `append` — function
- `assign` — subroutine
- `at_chunk_rec` — subroutine
- `at_process_rec` — subroutine
- `at_project_rec` — subroutine
- `at_rec` — subroutine
- `begin` — function
- `end` — function
- `end_iter` — function
- `get_busy_flags` — function
- `get_ids` — function
- `get_included_flags` — function
- `get_nptcls` — function
- `get_nptcls_sel_tot` — function
- `get_nptcls_tot` — function
- `get_processed_flags` — function
- `get_projfiles` — function
- `get_projnames` — function
- `iter_advance` — subroutine
- `iter_equals` — function
- `iter_get_chunk_rec` — subroutine
- `iter_get_process_rec` — subroutine
- `iter_get_project_rec` — subroutine
- `iter_get_rec` — subroutine
- `iter_next` — subroutine
- `iter_replace` — subroutine
- `iter_valid` — function
- `kill` — subroutine
- `print_file` — subroutine
- `print_proc` — subroutine
- `print_proj` — subroutine
- `push2chunk_list` — subroutine
- `push2project_list` — subroutine
- `push_back` — subroutine
- `rec_print` — subroutine
- `replace_at` — subroutine
- `replace_iterator` — subroutine
- `replace_with` — subroutine
- `set_included_flags` — subroutine
- `size` — function
- `slice` — subroutine
- `subset_to` — subroutine

---
## Module: simple_rec_list_tester

Files:
- `utils/structs/simple_rec_list_tester.f90`

Uses:
- `simple_defs`
- `simple_rec_list`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_rec_list_tests` — subroutine

Private symbols:
- `test_append_operator` — subroutine
- `test_assign_semantics` — subroutine
- `test_default_and_size` — subroutine
- `test_flag_included` — subroutine
- `test_included_mask` — subroutine
- `test_iterator_replace` — subroutine
- `test_iterator_traversal` — subroutine
- `test_kill_behavior` — subroutine
- `test_particle_sums` — subroutine
- `test_push_back_and_at` — subroutine
- `test_replace_at` — subroutine
- `test_slice_subset` — subroutine

---
## Module: simple_reconstructor

Files:
- `main/volume/simple_reconstructor.f90`

Uses:
- `simple_core_module_api`
- `simple_fftw3`
- `simple_gridding`
- `simple_image`
- `simple_math`
- `simple_parameters`
- `simple_sp_project`

Public symbols:
- `exp_samples` — type
- `gridding_half_restore` — type
- `insert_planes_multi` — subroutine
- `read_rho_header` — subroutine
- `reconstructor` — type
- `wiener_ridge` — function

Private symbols:
- `add_invtausq2rho` — subroutine
- `alloc_exp_samples` — subroutine
- `alloc_rho` — subroutine
- `apply_weight` — subroutine
- `apply_weight_sums` — subroutine
- `check_rho_header` — subroutine
- `compress_exp` — subroutine
- `dealloc_exp` — subroutine
- `dealloc_rho` — subroutine
- `deapodize_self` — subroutine
- `deapodize_volume` — subroutine
- `exp_samples_gather` — subroutine
- `exp_samples_gather_vbuf` — subroutine
- `exp_samples_get_n` — function
- `exp_samples_kill` — subroutine
- `exp_samples_new` — subroutine
- `exp_samples_new_plane` — subroutine
- `exp_samples_put_plane` — subroutine
- `exp_samples_set` — subroutine
- `expand_exp` — subroutine
- `finalize_gridding_half_restore` — subroutine
- `floor_rho_shellwise` — subroutine
- `gather_samples` — subroutine
- `get_kbwin` — function
- `get_rho_copy` — subroutine
- `insert_plane_oversamp` — subroutine
- `kill_gridding_half_restore` — subroutine
- `kill_reconstructor` — subroutine
- `new_accumulator` — subroutine
- `new_gridding_half_restore` — subroutine
- `pad_with_zeros` — subroutine
- `prepare_gridding_half_final` — subroutine
- `project_fplane` — subroutine
- `read_raw_rho` — subroutine
- `read_rho` — subroutine
- `reset` — subroutine
- `reset_exp` — subroutine
- `restore_base` — subroutine
- `restore_final` — subroutine
- `sampl_dens_correct` — subroutine
- `set_sh_lim` — subroutine
- `sum_reduce` — subroutine
- `validate_restore_source` — subroutine
- `validate_restore_target` — subroutine
- `write_absfc_as_mrc` — subroutine
- `write_raw_accum` — subroutine
- `write_rho` — subroutine
- `write_rho_as_mrc` — subroutine

Local symbols (enclosing procedure scope):
- `insert_plane_oversamp::kb_apod_vecs_3d_fast` — subroutine
- `insert_planes_multi::kb_weights` — subroutine

---
## Module: simple_reconstructor_openmpoffload

Files:
- `main/volume/simple_reconstructor_openmpoffload.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_kbinterpol`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_math`
- `simple_parameters`
- `simple_reconstructor`

Public symbols:
- `calc_3Drec_gpu` — subroutine

Private symbols:
- `insert_all_slices` — subroutine
- `insert_slices_batch` — subroutine
- `update_ptcls_arrays` — subroutine

Local symbols (enclosing procedure scope):
- `insert_all_slices::prep_batch` — subroutine

---
## Module: simple_reconstructor_pcg

Files:
- `main/volume/simple_reconstructor_pcg.f90`

Uses:
- `ieee_arithmetic`
- `iso_c_binding`
- `iso_fortran_env`
- `simple_cartesian_fourier`
- `simple_core_module_api`
- `simple_ctf`
- `simple_gridding`
- `simple_image`
- `simple_pcg_atomic_prior`
- `simple_pcg_lattice`
- `simple_pcg_solver`
- `simple_rnd`

Public symbols:
- `handle_cold_restart_outcome` — subroutine
- `measure_closed_form_agreement` — subroutine
- `pcg_raw_accum_compatible` — function
- `read_pcg_raw_accum_header` — subroutine
- `reconstructor_pcg` — type
- `report_closed_form_agreement` — subroutine
- `report_pcg_solve` — subroutine
- `validate_solved_map` — subroutine
- `write_closed_form_diagnostics` — subroutine

Private symbols:
- `absT2_plane` — subroutine
- `accumulate_absT2` — subroutine
- `accumulate_batch` — subroutine
- `accumulate_rhs_density` — subroutine
- `add_half_slice` — subroutine
- `add_raw_accum` — subroutine
- `add_raw_accum_weighted` — subroutine
- `add_ridge` — subroutine
- `adjoint_plane_add` — subroutine
- `apply_adjoint_all` — function
- `apply_fourier_diagonal` — function
- `apply_ml_prior` — function
- `apply_normal` — function
- `apply_normal_kernel` — function
- `apply_normal_matrixfree` — function
- `apply_precond` — function
- `apply_preconditioner` — function
- `atomic_model_map` — function
- `begin_accum` — subroutine
- `begin_reduction` — subroutine
- `build_env` — subroutine
- `build_hk_luts` — subroutine
- `build_kernel` — subroutine
- `build_ml_prior_from_density` — subroutine
- `build_operators` — subroutine
- `build_precond` — subroutine
- `build_transfer` — function
- `calibrate_kernel` — subroutine
- `compare_raw_accum` — subroutine
- `deapod_mul` — subroutine
- `dot_real_volume` — function
- `embed_fft` — subroutine
- `end_accum` — subroutine
- `ensure_wimg` — subroutine
- `extract_native_plane` — function
- `finalize_density_accum` — subroutine
- `finalize_khat` — subroutine
- `fit_atom_amplitudes` — subroutine
- `fold_and_ifft` — function
- `fold_ridges_into_precond` — subroutine
- `forward_plane` — subroutine
- `fourier_dot` — function
- `get_atomic_lambda` — function
- `get_data_scale` — function
- `get_effective_lambda` — function
- `get_env` — function
- `get_invenv` — function
- `get_lims2` — function
- `get_lims3` — function
- `get_ml_prior` — subroutine
- `get_ml_prior_stats` — subroutine
- `get_nptcls` — function
- `get_raw_accum` — subroutine
- `get_realspace_diagonal` — function
- `get_rhs` — subroutine
- `get_solve_support` — subroutine
- `get_solvent_lambda` — function
- `ifft_crop` — subroutine
- `kill` — subroutine
- `measure_kernel_scale` — function
- `ml_shrink_factor` — function
- `new` — subroutine
- `pack_half_slice` — subroutine
- `pcg_tic` — function
- `pcg_toc` — function
- `prep_particles` — subroutine
- `prepare_fused_planes` — subroutine
- `reconstruction_vector_apply` — subroutine
- `reconstruction_vector_dot` — function
- `reconstruction_vector_precond` — subroutine
- `reconstruction_vector_size` — function
- `reconstructor_pcg_adapter` — type
- `regularize_finalized` — subroutine
- `report_finalize_profile` — subroutine
- `report_profile` — subroutine
- `reset_finalize_profile` — subroutine
- `reset_profile` — subroutine
- `scale_raw_accum` — subroutine
- `set_atomic_prior` — subroutine
- `set_deapod` — subroutine
- `set_lambda_relative` — subroutine
- `set_ml_prior` — subroutine
- `set_op_mode` — subroutine
- `set_solvent_prior` — subroutine
- `set_sym` — subroutine
- `set_volume` — subroutine
- `shrink_by_ml_prior` — subroutine
- `solve` — subroutine
- `solve_accum` — subroutine
- `solve_core` — subroutine
- `solve_regularized` — subroutine
- `solve_with_cold_restart` — subroutine
- `sym_rots` — function
- `transfer_plane_cmplx` — subroutine
- `unshrink_by_ml_prior` — subroutine
- `update_lambda_from_density` — subroutine
- `window_div` — subroutine
- `window_mul` — subroutine
- `write_raw_accum` — subroutine

---
## Module: simple_reconstructor_tester

Files:
- `main/volume/simple_reconstructor_tester.f90`

Uses:
- `simple_defs`
- `simple_fplane_grid`
- `simple_image`
- `simple_math_ft`
- `simple_memoize_ft_maps`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_reconstructor`
- `simple_sp_project`
- `simple_sym`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_reconstructor_tests` — subroutine

Private symbols:
- `test_fplane_shift` — subroutine
- `test_weighted_insertion_linearity` — subroutine
- `test_wiener_ridge_scale_free` — subroutine

---
## Module: simple_reextract_strategy

Files:
- `main/strategies/parallelization/simple_reextract_strategy.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_parameters`
- `simple_particle_extractor`
- `simple_qsys_env`
- `simple_sp_project`

Public symbols:
- `create_reextract_strategy` — function
- `reextract_distr_strategy` — type
- `reextract_inmem_strategy` — type
- `reextract_strategy` — type

Private symbols:
- `box_inside` — function
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_end_message` — function
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `endmsg_interface` — function
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_end_message` — function
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `set_reextract_defaults` — subroutine
- `validate_reextract_cline` — subroutine

Local symbols (enclosing procedure scope):
- `inmem_execute::killimgbatch_local` — subroutine
- `inmem_execute::prepimgbatch_local` — subroutine

---
## Module: simple_refine2D_fnames

Files:
- `defs/simple_refine2D_fnames.f90`

Uses:
- `simple_defs_fname`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `refine2D_cavgs_fbody` — function
- `refine2D_cavgs_fname` — function
- `refine2D_cavgs_half_fname` — function
- `refine2D_cavgs_iter` — function
- `refine2D_cavgs_ranked_fname` — function

---
## Module: simple_refine2D_strategy

Files:
- `main/strategies/parallelization/simple_refine2D_strategy.f90`

Uses:
- `iso_fortran_env`
- `simple_bench_io`
- `simple_builder`
- `simple_cavg_sums`
- `simple_cavgassemble_strategy`
- `simple_classaverager`
- `simple_cmdline`
- `simple_commanders_euclid`
- `simple_commanders_imgops`
- `simple_commanders_mkcavgs`
- `simple_convergence`
- `simple_core_module_api`
- `simple_eo_partition`
- `simple_euclid_sigma2`
- `simple_eul_prob_tab_utils`
- `simple_exec_helpers`
- `simple_gui_utils`
- `simple_matcher_ptcl_io`
- `simple_matcher_smpl_and_lplims`
- `simple_parameters`
- `simple_procimgstk`
- `simple_progress`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_refine2d_fnames`
- `simple_sigma2_bootstrap`
- `simple_starproject`
- `simple_strategy2d_matcher`
- `simple_strategy_sigma2`
- `simple_stream_utils`

Public symbols:
- `create_refine2D_strategy` — function
- `refine2D_distr_strategy` — type
- `refine2D_inmem_strategy` — type
- `refine2D_strategy` — type

Private symbols:
- `carryover_needs_full_update` — function
- `check_refine2D_convergence` — subroutine
- `cleanup_distributed_iteration_artifacts` — subroutine
- `cleanup_interface` — subroutine
- `distr_anneal_tail_start` — function
- `distr_cleanup` — subroutine
- `distr_execute_iteration` — subroutine
- `distr_finalize_iteration` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `ensure_canonical_sigma2` — subroutine
- `exec_iter_interface` — subroutine
- `execute_make_cavgs` — subroutine
- `finalize_iter_interface` — subroutine
- `finalize_run_interface` — subroutine
- `gen_jpeg` — subroutine
- `init_interface` — subroutine
- `init_refine2D_refs` — subroutine
- `init_standard_refs` — subroutine
- `init_tseries_refs` — subroutine
- `inmem_anneal_tail_start` — function
- `inmem_cleanup` — subroutine
- `inmem_execute_iteration` — subroutine
- `inmem_finalize_iteration` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `refine2D_bench_state` — type
- `restore_cavgs_inmem` — subroutine
- `run_distributed_cavg_assembly` — subroutine
- `tail_start_interface` — function
- `write_strategy_bench_report` — subroutine

---
## Module: simple_refine3D_fnames

Files:
- `defs/simple_refine3D_fnames.f90`

Uses:
- `simple_defs_fname`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `refine3D_carried_bundle_fname` — function
- `refine3D_carried_context_fname` — function
- `refine3D_carried_membership_fname` — function
- `refine3D_carried_rec_fbody` — function
- `refine3D_cart_refvols_fname` — function
- `refine3D_frozen_context_fname` — function
- `refine3D_frozen_manifest_fname` — function
- `refine3D_frozen_rec_fbody` — function
- `refine3D_fsc_fbody` — function
- `refine3D_fsc_fname` — function
- `refine3D_fsc_plot_fbody` — function
- `refine3D_iter_vol_fname` — function
- `refine3D_oris_heatmap_fname` — function
- `refine3D_partial_rec_fbody` — function
- `refine3D_partial_rec_fname` — function
- `refine3D_partial_rho_fname` — function
- `refine3D_pcg_raw_accum_fname` — function
- `refine3D_pcg_trail_accum_fname` — function
- `refine3D_pcg_trail_manifest_fname` — function
- `refine3D_reproj_model_diag_fname` — function
- `refine3D_reproj_model_fname` — function
- `refine3D_reprojs_fname` — function
- `refine3D_resolution_txt_fbody` — function
- `refine3D_startvol_fbody` — function
- `refine3D_startvol_fname` — function
- `refine3D_startvol_half_fname` — function
- `refine3D_state_halfvol_fname` — function
- `refine3D_state_vol_fbody` — function
- `refine3D_state_vol_fname` — function
- `refine3D_state_vol_suffix_fname` — function
- `refine3D_trail_manifest_fname` — function
- `refine3D_trail_rec_fbody` — function
- `refine3D_trail_rec_fname` — function
- `refine3D_trail_rho_fname` — function

Private symbols:
- `half_suffix` — function
- `iter_tag` — function
- `part_tag` — function
- `state_tag` — function

---
## Module: simple_refine3D_local_tester

Files:
- `main/commanders/simple/simple_refine3D_local_tester.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_cmdline`
- `simple_commanders_refine3d`
- `simple_fileio`
- `simple_image`
- `simple_parameters`
- `simple_refine3d_provenance`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_dispatch`
- `simple_test_utils`

Public symbols:
- `run_all_refine3D_local_tests` — subroutine

Private symbols:
- `expected_reason` — function
- `file_snapshot` — type
- `log_has` — function
- `rejection_route` — subroutine
- `snapshot` — subroutine
- `write_inputs` — subroutine

---
## Module: simple_refine3D_provenance

Files:
- `main/volume/simple_refine3D_provenance.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_defs`
- `simple_defs_fname`
- `simple_error`
- `simple_fileio`
- `simple_halfmap_diagnostics`
- `simple_oris`
- `simple_parameters`
- `simple_ptcl_layout`
- `simple_sigma2_state_file`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `begin_reference_set` — subroutine
- `local_import_marker` — function
- `provenance_file` — function
- `provenance_pending_file` — function
- `publish_reference_set` — subroutine
- `publish_reuse_provenance` — subroutine
- `refine3D_provenance` — type
- `reuse_provenance_matches` — function

Private symbols:
- `contains_file` — function
- `get_anchor` — function
- `hash_observations` — subroutine
- `hash_reals` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `particle_identity` — subroutine
- `read_record` — function
- `reference_derivatives_match` — function
- `reference_matches_project` — function
- `reuse_matches` — function
- `settings_identity` — function
- `validate_reference_set` — function
- `write_record` — subroutine

Local symbols (enclosing procedure scope):
- `publish_reference_set::append` — subroutine

---
## Module: simple_refine3D_provenance_tester

Files:
- `main/volume/simple_refine3D_provenance_tester.f90`

Uses:
- `simple_fileio`
- `simple_halfmap_diagnostics`
- `simple_image`
- `simple_parameters`
- `simple_refine3d_provenance`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_refine3D_provenance_tests` — subroutine

Private symbols:
- `test_provenance` — subroutine
- `test_support_provenance` — subroutine

---
## Module: simple_refine3D_stage_plan

Files:
- `main/simple_refine3D_stage_plan.f90`

Uses:
- `simple_error`
- `simple_estimate_ssnr`
- `simple_math_ft`
- `simple_type_defs`

Public symbols:
- `plan_refine3D_frequency_stages` — subroutine
- `refine3D_stage_plan_entry` — type

---
## Module: simple_refine3D_startup

Files:
- `main/simple_refine3D_startup.f90`

Uses:
- `simple_cmdline`
- `simple_parameters`
- `simple_string`

Public symbols:
- `refine3D_startup` — type

Private symbols:
- `begin_iteration` — function
- `finish_iteration` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `preserves_history` — function
- `read_cline` — subroutine
- `write_cline` — subroutine

---
## Module: simple_refine3D_startup_tester

Files:
- `main/simple_refine3D_startup_tester.f90`

Uses:
- `simple_cmdline`
- `simple_parameters`
- `simple_refine3d_startup`
- `simple_test_utils`

Public symbols:
- `run_all_refine3D_startup_tests` — subroutine

---
## Module: simple_refine3D_states_controller

Files:
- `main/simple_refine3D_states_controller.f90`

Uses:
- `simple_class_sample_io`
- `simple_cmdline`
- `simple_commander_base`
- `simple_commanders_flex_pca`
- `simple_commanders_rec`
- `simple_defs`
- `simple_defs_fname`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_fileio`
- `simple_gui_communicator`
- `simple_linalg`
- `simple_math`
- `simple_oris`
- `simple_parameters`
- `simple_refine3d_fnames`
- `simple_refine3d_stage_plan`
- `simple_solve3d_utils`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_type_defs`
- `simple_view_partition_sampling`

Public symbols:
- `refine3D_states_controller` — type
- `refine3D_states_policy` — type

Private symbols:
- `default_nsample` — function
- `ensure_all_updated` — subroutine
- `initialize_state_volumes` — subroutine
- `initialize_states` — subroutine
- `kill` — subroutine
- `label` — function
- `new` — subroutine
- `next_iteration` — function
- `normalize_inputs` — subroutine
- `plan` — subroutine
- `prepare_sampling_units` — function
- `project_has_multistate_labels` — function
- `project_state_volumes_compatible` — function
- `read_maxits` — subroutine
- `read_state_overlap` — function
- `read_update_coverage` — subroutine
- `reject_input_volumes` — subroutine
- `report_coverage` — subroutine
- `resolve_sampling` — subroutine
- `resolve_states` — subroutine
- `run_flex_pca` — subroutine
- `run_frequency_march` — subroutine
- `run_march` — subroutine
- `run_missing_update` — subroutine
- `run_stage` — subroutine
- `set_downscaling` — subroutine
- `set_nstates` — subroutine
- `set_sampling` — subroutine
- `startup_maps` — subroutine
- `validate` — subroutine

---
## Module: simple_refine3D_strategy

Files:
- `main/strategies/parallelization/simple_refine3D_strategy.f90`

Uses:
- `iso_fortran_env`
- `simple_bench_io`
- `simple_builder`
- `simple_cluster_seed`
- `simple_cmdline`
- `simple_commanders_euclid`
- `simple_commanders_rec`
- `simple_commanders_volops`
- `simple_convergence`
- `simple_core_module_api`
- `simple_decay_funs`
- `simple_defs_fname`
- `simple_eo_partition`
- `simple_error`
- `simple_eul_prob_tab_utils`
- `simple_exec_helpers`
- `simple_fileio`
- `simple_fsc`
- `simple_halfmap_diagnostics`
- `simple_imghead`
- `simple_matcher_ptcl_io`
- `simple_matcher_refvol_utils`
- `simple_matcher_smpl_and_lplims`
- `simple_oris`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_rec3d_pcg_strategy`
- `simple_refine3d_fnames`
- `simple_refine3d_provenance`
- `simple_refine3d_startup`
- `simple_sigma2_bootstrap`
- `simple_sigma2_files`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_strategy3d_cont`
- `simple_strategy3d_matcher`
- `simple_strategy_sigma2`
- `simple_string`
- `simple_syslib`
- `simple_trail_chain_manifest`
- `simple_volassemble_strategy`

Public symbols:
- `adopt_refine3D_reference` — subroutine
- `create_refine3D_strategy` — function
- `refine3D_distr_strategy` — type
- `refine3D_inmem_strategy` — type
- `refine3D_strategy` — type
- `select_refine3D_state_rows` — subroutine
- `stage_reference_derivatives` — function
- `stage_refine3D_local_chain` — function
- `stage_refine3D_local_references` — subroutine
- `stage_refine3D_local_sigma2` — subroutine
- `strip_refine3D_search_only_args` — subroutine

Private symbols:
- `activate_ptcl3D_states_from_selection` — subroutine
- `assemble_refine3D_pcg` — subroutine
- `assert_multistate_populations` — subroutine
- `canonical_sigma_commit_deferred` — function
- `carry_over_trail_rec_chains` — subroutine
- `check_cls3D_field` — subroutine
- `check_polish_request` — subroutine
- `cleanup_interface` — subroutine
- `distr_anneal_tail_start` — function
- `distr_cleanup` — subroutine
- `distr_execute_iteration` — subroutine
- `distr_finalize_iteration` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `ensure_canonical_sigma2` — subroutine
- `exec_iter_interface` — subroutine
- `finalize_iter_interface` — subroutine
- `finalize_run_interface` — subroutine
- `get_state_pops` — subroutine
- `init_interface` — subroutine
- `inmem_anneal_tail_start` — function
- `inmem_cleanup` — subroutine
- `inmem_execute_iteration` — subroutine
- `inmem_finalize_iteration` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `invalidate_fresh_start_refs_from_volumes` — subroutine
- `labels_may_change` — function
- `materialize_reprojection_model` — subroutine
- `materialize_reprojection_model_master` — subroutine
- `polish_follows` — function
- `prepare_assembly_cline` — subroutine
- `publish_before_file_assembly` — function
- `refine3D_bench_state` — type
- `refresh_matching_lp_from_project` — subroutine
- `refresh_resolution_fields_from_fsc` — subroutine
- `remove_partial_rec_files` — subroutine
- `remove_pcg_raw_files` — subroutine
- `report_rec_master_phase` — subroutine
- `reset_refine3D_bench` — subroutine
- `seed_multistate_startup_labels` — subroutine
- `sigma_update_enabled` — function
- `tail_start_interface` — function
- `write_iteration_states` — subroutine
- `write_strategy_bench_report` — subroutine

---
## Module: simple_refine3D_strategy_tester

Files:
- `main/strategies/parallelization/simple_refine3D_strategy_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_builder`
- `simple_cmdline`
- `simple_defs_fname`
- `simple_fileio`
- `simple_halfmap_diagnostics`
- `simple_image`
- `simple_matcher_refvol_utils`
- `simple_parameters`
- `simple_refine3d_fnames`
- `simple_refine3d_provenance`
- `simple_refine3d_strategy`
- `simple_sigma2_state_file`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`
- `simple_trail_chain_manifest`

Public symbols:
- `run_all_refine3D_staging_tests` — subroutine

Private symbols:
- `test_chain_staging` — subroutine
- `test_reference_derivatives` — subroutine
- `test_reference_transition` — subroutine
- `test_startup_staging` — subroutine

Local symbols (enclosing procedure scope):
- `test_reference_transition::write_triplet` — subroutine

---
## Module: simple_refine3D_workflow_controller

Files:
- `main/simple_refine3D_workflow_controller.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_defs_fname`
- `simple_error`
- `simple_fileio`
- `simple_halfmap_diagnostics`
- `simple_image`
- `simple_imghead`
- `simple_oris`
- `simple_parameters`
- `simple_refine3d_fnames`
- `simple_refine3d_provenance`
- `simple_refine3d_strategy`
- `simple_sigma2_files`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_strategy3d_cont`
- `simple_string`
- `simple_syslib`

Public symbols:
- `refine3D_workflow_controller` — type

Private symbols:
- `kill` — subroutine
- `load_project` — subroutine
- `new` — subroutine
- `preflight` — subroutine
- `preflight_cont_startup` — subroutine
- `prepare_external_init_vol` — subroutine
- `prepare_startup` — subroutine
- `publish_project` — subroutine
- `reference_grid_matches` — function
- `stage_local_trailing_chain` — subroutine

---
## Module: simple_refine_motion_model_strategy

Files:
- `main/strategies/parallelization/simple_refine_motion_model_strategy.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_dose_series_manifest`
- `simple_eo_partition`
- `simple_matcher_ptcl_io`
- `simple_motion_model`
- `simple_oris`
- `simple_parameters`
- `simple_particle_extractor`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_rec3d_strategy`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_stack_io`
- `simple_syslib`

Public symbols:
- `create_refine_motion_model_strategy` — function
- `prepare_fraction_dose_manifest` — function
- `reconstruct_motion_fractions` — subroutine
- `refine_motion_model_distr_strategy` — type
- `refine_motion_model_inmem_strategy` — type
- `refine_motion_model_strategy` — type

Private symbols:
- `assemble_motion_fractions` — subroutine
- `cleanup_interface` — subroutine
- `count_fractioned_stacks` — function
- `distr_cleanup` — subroutine
- `distr_end_message` — function
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `endmsg_interface` — function
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `fraction_frame_range` — function
- `fraction_output_dir` — function
- `fraction_state1_rows` — subroutine
- `gnu_quote` — function
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_end_message` — function
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `number` — function
- `plot_fraction_fsc` — subroutine
- `remove_fraction_partials` — subroutine
- `shell_quote` — function
- `validate_fraction_range` — subroutine
- `write_fraction_fsc_matrix` — subroutine

Local symbols (enclosing procedure scope):
- `inmem_execute::validate_input` — subroutine
- `plot_fraction_fsc::gp` — subroutine

---
## Module: simple_refine_motion_model_tester

Files:
- `main/strategies/parallelization/simple_refine_motion_model_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_cmdline`
- `simple_defs`
- `simple_dose_series_manifest`
- `simple_fileio`
- `simple_flex_pca_application_tester`
- `simple_image`
- `simple_motion_model`
- `simple_parameters`
- `simple_refine_motion_model_strategy`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `check_fraction_dose_export` — subroutine
- `run_all_refine_motion_model_lib_tests` — subroutine
- `run_all_refine_motion_model_tests` — subroutine

Private symbols:
- `check_manifest_partition` — subroutine
- `write_saved_model` — subroutine

---
## Module: simple_relion

Files:
- `main/star/simple_relion.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_epu_xml`
- `simple_sp_project`
- `simple_starfile_wrappers`
- `simple_syslib`

Public symbols:
- `relion_project` — type

Private symbols:
- `allocate_opticsgroups` — subroutine
- `create` — subroutine
- `find_movienames` — subroutine
- `generate_epu_tiltgroups` — subroutine
- `generate_single_tiltgroup` — subroutine
- `generate_xml_tiltgroups` — subroutine
- `h_clust` — subroutine
- `write_corrected_micrographs_star` — subroutine
- `write_micrographs_star` — subroutine
- `write_particles2D_star` — subroutine

Local symbols (enclosing procedure scope):
- `write_particles2D_star::write_pose` — subroutine

---
## Module: simple_rnd

Files:
- `utils/math/simple_rnd.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `omp_lib`
- `simple_defs`
- `simple_error`
- `simple_math`
- `simple_syslib`

Public symbols:
- `ensure_run_seed` — function
- `get_rnd_state` — subroutine
- `has_global_seed` — function
- `irnd_gasdev` — function
- `irnd_uni` — function
- `irnd_uni_pair` — function
- `mnorm_smp` — function
- `multinomal` — function
- `poisson_rnd` — function
- `put_rnd_state` — subroutine
- `r8po_fa` — subroutine
- `ran3` — function
- `rnd_4dim_sphere_pnt` — function
- `rnd_bits_to_unit` — function
- `rnd_draw_key` — function
- `rnd_inds` — function
- `rnd_key` — type
- `rnd_key_lanes` — function
- `rnd_keyed` — function
- `rnd_keyed_bits` — function
- `rnd_keyed_int` — function
- `rnd_keyed_partial_shuffle` — subroutine
- `rnd_run_seed` — function
- `rnd_scope` — function
- `rnd_scope_key` — function
- `seed_rnd` — subroutine
- `set_run_seed` — subroutine
- `shcloc` — function

Private symbols:
- `gasdev_bounded` — function
- `gasdev_mean_sdev` — function
- `gasdev_std` — function
- `greedy_sampling_alloc` — function
- `greedy_sampling_prealloc` — function
- `partial_shuffle_int` — subroutine
- `partial_shuffle_real` — subroutine
- `put_fixed_seed` — subroutine
- `put_simple_seed` — subroutine
- `ran3arr_1d` — subroutine
- `ran3arr_2d` — subroutine
- `randn_1d` — function
- `randn_2d` — function
- `rnd_mix32` — function
- `seed_rnd_fixed_i4` — subroutine
- `seed_rnd_fixed_i8` — subroutine
- `shuffle_int` — subroutine
- `shuffle_real` — subroutine
- `simple_seed_base` — function

Local symbols (enclosing procedure scope):
- `put_simple_seed::c_setenv` — function
- `seed_rnd::lcg` — function

---
## Module: simple_rnd_tester

Files:
- `utils/math/simple_rnd_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_defs`
- `simple_error`
- `simple_math`
- `simple_rnd`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_rnd_tests` — subroutine

Private symbols:
- `test_global_seed` — subroutine
- `test_irnd_uni` — subroutine
- `test_keyed_int` — subroutine
- `test_keyed_mapping` — subroutine
- `test_keyed_order` — subroutine
- `test_keyed_partial_shuffle` — subroutine
- `test_keyed_scope` — subroutine
- `test_keyed_sensitivity` — subroutine
- `test_keyed_vectors` — subroutine
- `test_multinomal` — subroutine
- `test_partial_shuffle` — subroutine
- `test_poisson` — subroutine
- `test_rnd_state_restore` — subroutine
- `test_shuffle` — subroutine
- `vec_key` — function

Local symbols (enclosing procedure scope):
- `test_keyed_partial_shuffle::perm_code` — function
- `test_keyed_partial_shuffle::subset_code` — function

---
## Module: simple_sauron

Files:
- `main/ori/simple_sauron.f90`

Uses:
- `simple_chash`
- `simple_defs`
- `simple_defs_ori`
- `simple_hash`
- `simple_string_utils`

Private symbols:
- `sauron_ori_parser_hash` — subroutine
- `sauron_ori_parser_pparms` — subroutine
- `str2realdp_local` — subroutine

---
## Module: simple_segdiam_bin_picker

Files:
- `main/pick/simple_segdiam_bin_picker.f90`

Uses:
- `simple_core_module_api`
- `simple_gui_utils`
- `simple_image`
- `simple_micproc`
- `simple_nrtxtfile`
- `simple_picksegdiam`
- `simple_sp_project`

Public symbols:
- `segdiam_bin_picker` — type

Private symbols:
- `bin_index_from_bounds` — function
- `bins_set` — function
- `ensure_real_capacity` — subroutine
- `get_accepted_bins` — function
- `get_box` — function
- `get_mskdiam` — function
- `kill` — subroutine
- `new` — subroutine
- `pick` — subroutine
- `print_diam_stats` — subroutine

---
## Module: simple_segmentation

Files:
- `main/image_processing/simple_segmentation.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_image_bin`
- `simple_neighs`

Public symbols:
- `automatic_thresh_sobel` — subroutine
- `canny` — subroutine
- `canny_edge` — subroutine
- `detect_peak_thres_1` — subroutine
- `detect_peak_thres_2` — subroutine
- `detect_peak_thres_fdr` — subroutine
- `detect_peak_thres_for_npeaks` — subroutine
- `double_thresh` — subroutine
- `non_max_supp` — subroutine
- `otsu_img` — subroutine
- `otsu_robust_fast` — subroutine
- `refine_peak_thres_sortmeans` — subroutine
- `sauvola` — subroutine
- `sobel` — subroutine

Local symbols (enclosing procedure scope):
- `detect_peak_thres_fdr::threshold_for_no_peaks` — subroutine
- `detect_peak_thres_fdr::threshold_for_npeaks` — subroutine

---
## Module: simple_segmentation_tester

Files:
- `main/image_processing/simple_segmentation_tester.f90`

Uses:
- `simple_defs`
- `simple_image`
- `simple_math`
- `simple_rnd`
- `simple_segmentation`
- `simple_test_utils`

Public symbols:
- `run_all_segmentation_tests` — subroutine

Private symbols:
- `all_binary` — function
- `four_level` — subroutine
- `test_binarize` — subroutine
- `test_canny` — subroutine
- `test_detect_peak_thres` — subroutine
- `test_otsu` — subroutine
- `test_otsu_img` — subroutine
- `test_otsu_robust_fast` — subroutine
- `test_peak_thres_fdr` — subroutine
- `test_sauvola` — subroutine
- `test_sobel` — subroutine
- `three_level` — subroutine
- `two_level` — subroutine

---
## Module: simple_shiftmap_srch

Files:
- `main/strategies/search/simple_shiftmap_srch.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_defs`
- `simple_error`
- `simple_eul_prob_tab_utils`
- `simple_parameters`
- `simple_pftc_shsrch_grad`
- `simple_polarft_calc`
- `simple_srch_sort_loc`

Public symbols:
- `setup_shiftmap_srch` — subroutine
- `shiftmap_srch` — type

Private symbols:
- `get_result` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `on_disk` — function
- `prep` — subroutine
- `srch` — subroutine

Local symbols (enclosing procedure scope):
- `srch::alloc_results` — subroutine

---
## Module: simple_sigma2_bootstrap

Files:
- `main/sigma2/simple_sigma2_bootstrap.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_euclid`
- `simple_sigma2_files`

Public symbols:
- `ensure_sigma2_for_iteration` — subroutine
- `prepare_pspec_cline` — subroutine
- `prepare_residual_sigma2_pass_cline` — subroutine
- `sigma2_estimate_available` — function

---
## Module: simple_sigma2_bootstrap_tester

Files:
- `main/sigma2/simple_sigma2_bootstrap_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_defs`
- `simple_fileio`
- `simple_ptcl_layout`
- `simple_sigma2_bootstrap`
- `simple_sigma2_pspec_sums`
- `simple_sigma2_state`
- `simple_sigma2_state_file`
- `simple_sigma2_state_tester`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_sigma2_bootstrap_tests` — subroutine

Private symbols:
- `check_refused_sums` — subroutine
- `make_project` — subroutine
- `make_template` — subroutine
- `test_new_exec_dir` — subroutine
- `test_pspec_cline` — subroutine
- `test_pspec_sums` — subroutine
- `test_rebuild_entry` — subroutine
- `test_residual_pass_cline` — subroutine
- `write_store` — subroutine

---
## Module: simple_sigma2_files

Files:
- `fileio/simple_sigma2_files.f90`

Uses:
- `iso_fortran_env`
- `simple_core_module_api`
- `simple_euclid_sigma2`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_ptcl_layout`
- `simple_sigma2_state`
- `simple_sigma2_state_file`
- `simple_sp_project`

Public symbols:
- `canonical_sigma2_consumable` — function
- `load_sigma2` — subroutine

---
## Module: simple_sigma2_pspec_sums

Files:
- `main/sigma2/simple_sigma2_pspec_sums.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_defs`
- `simple_error`
- `simple_string_utils`

Public symbols:
- `bootstrap_seed_spectra` — subroutine
- `pspec_part_sums_fname` — function
- `read_pspec_part_sums` — subroutine
- `sum_valid_records` — subroutine
- `write_pspec_part_sums` — subroutine

Private symbols:
- `remove_negative_sigmas` — subroutine

---
## Module: simple_sigma2_state

Files:
- `main/sigma2/simple_sigma2_state.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_defs`
- `simple_sigma2_state_file`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `sigma2_state_assemble_candidate` — subroutine
- `sigma2_state_candidate_path` — function
- `sigma2_state_concatenate` — subroutine
- `sigma2_state_extend` — subroutine
- `sigma2_state_next_generation` — subroutine
- `sigma2_state_prepare_update` — subroutine
- `sigma2_state_range_path` — function
- `sigma2_state_reduce` — subroutine
- `sigma2_state_union_means` — subroutine
- `sigma2_state_validate_identity` — subroutine

Private symbols:
- `add_records` — subroutine
- `form_global` — subroutine
- `record_is_valid` — function
- `transaction_stem` — function

Local symbols (enclosing procedure scope):
- `sigma2_state_assemble_candidate::close_units` — subroutine
- `sigma2_state_concatenate::close_units` — subroutine
- `sigma2_state_extend::close_units` — subroutine

---
## Module: simple_sigma2_state_file

Files:
- `fileio/simple_sigma2_state_file.f90`

Uses:
- `iso_fortran_env`
- `simple_fileio`
- `simple_string`
- `simple_syslib`

Public symbols:
- `sigma2_state_artifact_signature` — function
- `sigma2_state_close_candidate` — subroutine
- `sigma2_state_create_candidate` — subroutine
- `sigma2_state_digest_begin` — function
- `sigma2_state_digest_file` — function
- `sigma2_state_digest_integer` — subroutine
- `sigma2_state_digest_text` — subroutine
- `sigma2_state_header` — type
- `sigma2_state_init_header` — subroutine
- `sigma2_state_open_candidate` — subroutine
- `sigma2_state_open_local_range` — subroutine
- `sigma2_state_open_store` — subroutine
- `sigma2_state_publish` — subroutine
- `sigma2_state_read_global` — subroutine
- `sigma2_state_read_header` — subroutine
- `sigma2_state_read_local_range` — subroutine
- `sigma2_state_read_local_range_rows` — subroutine
- `sigma2_state_read_particles` — subroutine
- `sigma2_state_read_particles_into` — subroutine
- `sigma2_state_read_particles_unit` — subroutine
- `sigma2_state_validate_file` — subroutine
- `sigma2_state_write_global` — subroutine
- `sigma2_state_write_local_range` — subroutine
- `sigma2_state_write_particles` — subroutine
- `sigma2_state_write_particles_unit` — subroutine

Private symbols:
- `digest_word` — subroutine
- `fnv_integer` — function
- `pack_state_header` — subroutine
- `set_state_layout` — subroutine
- `unpack_state_header` — subroutine
- `validate_header_values` — subroutine
- `write_header_unit` — subroutine

Local symbols (enclosing procedure scope):
- `sigma2_state_open_local_range::close_range` — subroutine

---
## Module: simple_sigma2_state_tester

Files:
- `main/sigma2/simple_sigma2_state_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_fileio`
- `simple_ptcl_layout`
- `simple_sigma2_bootstrap`
- `simple_sigma2_state`
- `simple_sigma2_state_file`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `patch_sigma2_header_word` — subroutine
- `run_all_sigma2_state_tests` — subroutine
- `write_test_sigma2_store` — subroutine

Private symbols:
- `assert_committed_generation` — subroutine
- `assert_file_unchanged` — subroutine
- `cleanup` — subroutine
- `read_file_bytes` — subroutine
- `require` — subroutine
- `require_full_reduction` — subroutine
- `require_ok` — subroutine
- `test_addon_copy_on_write` — subroutine
- `test_concatenation` — subroutine
- `test_estimate_available_on_disk` — subroutine
- `test_failed_publication` — subroutine
- `test_global_transaction` — subroutine
- `test_invalid_record_skip` — subroutine
- `test_particle_io_layout` — subroutine
- `test_recovery_guards` — subroutine
- `test_retired_header_words` — subroutine
- `test_unscheduled_rows_from_committed` — subroutine
- `test_update_preparation` — subroutine
- `write_committed_store` — subroutine

Local symbols (enclosing procedure scope):
- `test_retired_header_words::write_valid_store` — subroutine

---
## Module: simple_simple_volinterp

Files:
- `main/volume/simple_volinterp.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_projector`

Public symbols:
- `reproject` — function
- `rotvol` — function
- `rotvol_slim` — subroutine

---
## Module: simple_simulator

Files:
- `main/simple_simulator.f90`

Uses:
- `simple_core_module_api`
- `simple_ctf`
- `simple_image`
- `simple_wave_imaging`

Public symbols:
- `simimg` — subroutine
- `simimg_wave` — subroutine

---
## Module: simple_simulator_tester

Files:
- `main/simple_simulator_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_ctf`
- `simple_gridding`
- `simple_image`
- `simple_kbinterpol`
- `simple_math_ctf`
- `simple_memoize_ft_maps`
- `simple_projector`
- `simple_simulator`
- `simple_test_utils`
- `simple_wave_imaging`
- `simple_wave_projector`

Public symbols:
- `run_all_simulator_tests` — subroutine

Private symbols:
- `blob_map` — subroutine
- `deapod_envelope` — function
- `flat_particle` — subroutine
- `max_phase` — function
- `odd_even` — subroutine
- `test_bfac_per_particle` — subroutine
- `test_detector_noise` — subroutine
- `test_particle` — subroutine
- `test_weak_limit` — subroutine
- `wave_particle` — subroutine

---
## Module: simple_solve2D_chunk

Files:
- `main/solve/simple_solve2D_chunk.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_environment`
- `simple_parameters`
- `simple_qsys_env`
- `simple_rec_list`
- `simple_sp_project`

Public symbols:
- `solve2D_chunk` — type

Private symbols:
- `analyze2D` — subroutine
- `assign` — subroutine
- `calc_sigma2` — subroutine
- `copy` — subroutine
- `debug_print` — subroutine
- `display_iter` — subroutine
- `gen_final_cavgs` — subroutine
- `generate` — subroutine
- `get_projfile_fname` — function
- `has_converged` — function
- `init_chunk` — subroutine
- `is_available` — function
- `is_finished` — function
- `kill` — subroutine
- `print_info` — subroutine

---
## Module: simple_solve2D_controller

Files:
- `main/solve/simple_solve2D_controller.f90`

Uses:
- `simple_class_frcs`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_default_clines`
- `simple_image`
- `simple_image_msk`
- `simple_parameters`
- `simple_refine2d_fnames`
- `simple_sp_project`

Public symbols:
- `determine_solve2D_stages` — subroutine
- `frc05_ptcl_quantile` — subroutine
- `mskdiam2lplimits_refine2D` — subroutine
- `set_cline_refine2D_seed_pass` — subroutine
- `set_cline_refine2D_stage` — subroutine
- `set_solve2D_sampling_policy` — subroutine
- `solve2D_class_frc05_quantile` — subroutine
- `solve2D_first_shift_stage` — function
- `solve2D_mskdiam_estimate` — subroutine
- `solve2D_seed_pass_iter` — function
- `solve2D_stage_behind_frc05` — function
- `stage_params` — type

Private symbols:
- `build_refine2D_stage_cfg` — subroutine
- `emit_refine2D_stage_cfg` — subroutine
- `phase1_stage_endit` — function
- `refine2D_stage_cfg` — type
- `refine2D_stage_requires_full_assignment` — function
- `sample_coverage_nits` — function
- `sampled_stage_nits` — function
- `set_refine2D_stage_iteration_policy` — subroutine
- `set_refine2D_stage_mask_policy` — subroutine
- `set_refine2D_stage_objfun_policy` — subroutine
- `set_refine2D_stage_phase_policy` — subroutine
- `set_refine2D_stage_reference_policy` — subroutine
- `set_refine2D_stage_regular_refs` — subroutine
- `set_refine2D_stage_search_policy` — subroutine

---
## Module: simple_solve2D_controller_tester

Files:
- `main/solve/simple_solve2D_controller_tester.f90`

Uses:
- `simple_solve2d_controller`
- `simple_test_utils`

Public symbols:
- `run_all_solve2D_stage_rule_tests` — subroutine

Private symbols:
- `make_ladder` — subroutine
- `test_omission_never` — subroutine
- `test_omission_rule` — subroutine
- `test_quantile_empty` — subroutine
- `test_quantile_junk_classes` — subroutine
- `test_quantile_median` — subroutine

---
## Module: simple_solve3D_addon_dock_tester

Files:
- `main/solve/simple_solve3D_addon_dock_tester.f90`

Uses:
- `simple_defs`
- `simple_image`
- `simple_solve3d_addon_report`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_solve3D_addon_dock_tests` — subroutine

Private symbols:
- `test_rotated_frame` — subroutine
- `write_fixtures` — subroutine
- `write_map` — subroutine

---
## Module: simple_solve3D_addon_report

Files:
- `main/solve/simple_solve3D_addon_report.f90`

Uses:
- `simple_core_module_api`
- `simple_dock_vols`
- `simple_ori`
- `simple_oris`
- `simple_volpair_metrics`

Public symbols:
- `solve3D_addon_report` — type

Private symbols:
- `any_regressed` — function
- `check_state` — subroutine
- `compare_cohort` — subroutine
- `compare_fsc` — subroutine
- `compare_maps` — subroutine
- `corr_lp` — function
- `dock` — subroutine
- `get_cohort_res0143` — function
- `get_corr` — function
- `get_dock_angle` — function
- `get_dshell` — function
- `get_nstates` — function
- `get_verdict` — function
- `kill` — subroutine
- `new` — subroutine
- `print` — subroutine
- `read` — subroutine
- `set_populations` — subroutine
- `set_stage` — subroutine
- `stage_record` — type
- `state_record` — type
- `write` — subroutine

---
## Module: simple_solve3D_addon_report_tester

Files:
- `main/solve/simple_solve3D_addon_report_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_solve3d_addon_report`
- `simple_test_utils`

Public symbols:
- `run_all_solve3D_addon_report_tests` — subroutine

Private symbols:
- `fsc_curve` — function
- `test_fsc_verdicts` — subroutine
- `test_map_comparison` — subroutine
- `test_round_trip` — subroutine
- `write_fixtures` — subroutine
- `write_map` — subroutine

---
## Module: simple_solve3D_chunks_plan

Files:
- `main/solve/simple_solve3D_chunks_plan.f90`

Uses:
- `simple_core_module_api`
- `simple_sp_project`

Public symbols:
- `chunk_nsample` — function
- `chunks_report` — type
- `chunks_step` — type
- `derive_chunk_run_count` — subroutine
- `partition_selected_balanced` — subroutine
- `prepare_chunks_source` — subroutine
- `select_chunks_step` — subroutine

Private symbols:
- `add_step` — subroutine
- `get_nsteps` — function
- `kill` — subroutine
- `new_chunks` — subroutine
- `write` — subroutine

---
## Module: simple_solve3D_chunks_plan_tester

Files:
- `main/solve/simple_solve3D_chunks_plan_tester.f90`

Uses:
- `simple_defs`
- `simple_solve3d_chunks_plan`
- `simple_sp_project`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_solve3D_chunks_plan_tests` — subroutine

Private symbols:
- `chunk_size` — function
- `make_project` — subroutine
- `make_stacks` — subroutine
- `test_balanced_partition` — subroutine
- `test_chunk_nsample` — subroutine
- `test_chunk_run_count` — subroutine
- `test_chunks_report` — subroutine
- `test_source_and_step_projects` — subroutine

---
## Module: simple_solve3D_manifest

Files:
- `main/solve/simple_solve3D_manifest.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_ptcl_layout`
- `simple_sigma2_state_file`
- `simple_sp_project`

Public symbols:
- `manifest_overridable_input` — function
- `manifest_records_input` — function
- `solve3D_manifest` — type
- `solve3D_stage_record` — type

Private symbols:
- `add_artifact` — subroutine
- `add_input` — subroutine
- `build_lines` — subroutine
- `find_artifact` — function
- `fmt_real` — function
- `get_artifact` — subroutine
- `get_box` — function
- `get_fname` — function
- `get_input` — subroutine
- `get_last_stage` — function
- `get_nstages` — function
- `get_nstates` — function
- `get_rec_backend` — function
- `get_run_id` — function
- `get_smpd` — function
- `get_stage` — function
- `int64_str` — function
- `kill` — subroutine
- `logical_str` — function
- `matches_artifact` — function
- `new` — subroutine
- `optics_ctf_digest` — function
- `read` — subroutine
- `read_registered` — subroutine
- `record_artifacts` — subroutine
- `record_inputs` — subroutine
- `register` — subroutine
- `replay` — subroutine
- `set_ladder` — subroutine
- `set_sampling` — subroutine
- `set_solution` — subroutine
- `set_stage_line` — subroutine
- `stack_table_digest` — function
- `validate_frozen` — subroutine
- `write` — subroutine

Local symbols (enclosing procedure scope):
- `build_lines::push` — subroutine

---
## Module: simple_solve3D_manifest_tester

Files:
- `main/solve/simple_solve3D_manifest_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_defs`
- `simple_fileio`
- `simple_image`
- `simple_ptcl_layout`
- `simple_sigma2_state_file`
- `simple_solve3d_manifest`
- `simple_sp_project`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_solve3D_manifest_tests` — subroutine

Private symbols:
- `file_contains` — function
- `ladder` — function
- `make_manifest` — subroutine
- `make_project` — subroutine
- `overwrite_sigma_value` — subroutine
- `read_lines` — subroutine
- `same_text` — function
- `test_carried_bundle_artifact` — subroutine
- `test_file_refusals` — subroutine
- `test_frozen_validation` — subroutine
- `test_long_paths` — subroutine
- `test_overridable_inputs` — subroutine
- `test_registration` — subroutine
- `test_replay` — subroutine
- `test_round_trip` — subroutine
- `write_lines` — subroutine
- `write_sigma_state` — subroutine
- `write_volume` — subroutine
- `write_with_checksum` — subroutine

Local symbols (enclosing procedure scope):
- `test_frozen_validation::make_foreign` — subroutine

---
## Module: simple_solve3D_reseed_tester

Files:
- `main/solve/simple_solve3D_reseed_tester.f90`

Uses:
- `simple_ori`
- `simple_solve3d_utils`
- `simple_sp_project`
- `simple_string`
- `simple_sym`
- `simple_test_utils`

Public symbols:
- `run_all_solve3D_reseed_tests` — subroutine

Private symbols:
- `test_all_in_one_state` — subroutine
- `test_no_weak_state` — subroutine
- `test_one_empty_state` — subroutine
- `test_stage_skipping` — subroutine
- `test_symmetry_randomize` — subroutine
- `test_symsrch_stage` — subroutine
- `test_too_few_to_reseed` — subroutine
- `test_weak_by_fraction` — subroutine

---
## Module: simple_solve3D_utils

Files:
- `main/solve/simple_solve3D_controller.f90`
- `main/solve/simple_solve3D_utils.f90`

Uses:
- `simple_class_frcs`
- `simple_commanders_api`
- `simple_commanders_euclid`
- `simple_commanders_volops`
- `simple_halfmap_diagnostics`
- `simple_matcher_refvol_utils`
- `simple_parameters`
- `simple_refine3d_fnames`
- `simple_sigma2_bootstrap`
- `simple_solve3d_manifest`

Public symbols:
- `active_lp_schedule_nstages` — function
- `apply_refine3D_reconstruction_controls` — subroutine
- `calc_frozen_rec` — subroutine
- `calc_rec` — subroutine
- `commit_deferred_sigma_update` — subroutine
- `configure_final_pcg_solve_budget` — subroutine
- `exec_refine3D` — subroutine
- `gen_ortho_reprojs4viz` — subroutine
- `generate_random_volumes` — subroutine
- `inject_refine3D_volume` — subroutine
- `normalize_input_volumes` — subroutine
- `prep_class_command_lines` — subroutine
- `print_states` — subroutine
- `register_stage_volume` — subroutine
- `reseed_state_labels` — subroutine
- `set_cline_refine3D` — subroutine
- `set_cline_refine3D` — subroutine
- `set_lplims_from_frcs` — subroutine
- `set_lplims_from_input` — subroutine
- `set_lplims_from_manifest` — subroutine
- `set_symmetry_class_vars` — subroutine
- `solve3D_addon_ctx` — type
- `solve3D_aggressive_march_stage` — function
- `solve3D_aggressive_march_stage` — function
- `solve3D_cavgs_early_nstages` — function
- `solve3D_cavgs_early_nstages` — function
- `solve3D_cenlp_default` — function
- `solve3D_cenlp_default` — function
- `solve3D_full_sample_switch_frac` — function
- `solve3D_full_sample_switch_frac` — function
- `solve3D_independent_lpstop_default` — function
- `solve3D_independent_lpstop_default` — function
- `solve3D_independent_nstages_default` — function
- `solve3D_independent_nstages_default` — function
- `solve3D_lpstart_aggressive` — function
- `solve3D_lpstart_aggressive` — function
- `solve3D_lpstart_bounds` — function
- `solve3D_lpstart_bounds` — function
- `solve3D_lpstart_ini3D` — function
- `solve3D_lpstart_ini3D` — function
- `solve3D_lpstop_bounds` — function
- `solve3D_lpstop_bounds` — function
- `solve3D_lpstop_ini3D` — function
- `solve3D_lpstop_ini3D` — function
- `solve3D_lpsymsrch_lb` — function
- `solve3D_lpsymsrch_lb` — function
- `solve3D_ml_reg_start_stage` — function
- `solve3D_ml_reg_start_stage` — function
- `solve3D_nsample_default` — function
- `solve3D_nsample_default` — function
- `solve3D_nstages` — function
- `solve3D_nstages` — function
- `solve3D_nstages_ini3D` — function
- `solve3D_nstages_ini3D` — function
- `solve3D_nstages_ini3D_max` — function
- `solve3D_nstages_ini3D_max` — function
- `solve3D_rec_fbody` — function
- `solve3D_rec_fbody` — function
- `solve3D_remaining_niters` — function
- `solve3D_remaining_niters` — function
- `solve3D_stage_behind_fsc05` — function
- `solve3D_stage_behind_fsc05` — function
- `solve3D_stage_box_crop` — function
- `solve3D_stage_skipping` — function
- `solve3D_stage_skipping` — function
- `solve3D_stage_smpd_crop` — function
- `solve3D_state_fsc_lowpass` — function
- `solve3D_stoch_sampl_stage` — function
- `solve3D_stoch_sampl_stage` — function
- `solve3D_symsrch_stage` — function
- `solve3D_symsrch_stage` — function
- `solve3D_update_frac_max` — function
- `solve3D_update_frac_max` — function
- `solve3D_uses_delayed_input_lp` — function
- `solve3D_uses_delayed_input_lp` — function
- `solve3D_vol_supplied_stage` — function
- `solve3D_vol_supplied_stage` — function
- `stage_rec_is_euclid` — function
- `strip_pcg_backend_keys` — subroutine
- `strip_refine3D_planning_keys` — subroutine
- `strip_sampling_keys` — subroutine
- `symmetrize` — subroutine
- `symmetry_randomize_poses` — subroutine
- `write_final_rec_outputs` — subroutine
- `write_solve3D_lowpass_snapshot` — subroutine

Private symbols:
- `active_refine3D_nstages` — function
- `apply_refine3D_search_overrides` — subroutine
- `build_refine3D_stage_cfg` — subroutine
- `emit_refine3D_stage_cfg` — subroutine
- `force_full_sampling_mode` — function
- `full_budget_last_stage` — function
- `init_refine3D_iteration` — subroutine
- `nu_start_stage` — function
- `project_best_fsc05_resolution` — function
- `promote_stage_lp_from_fsc05` — subroutine
- `refine3D_stage_cfg` — type
- `set_refine3D_automsk_policy` — subroutine
- `set_refine3D_backend_policy` — subroutine
- `set_refine3D_balance_policy` — subroutine
- `set_refine3D_envfsc_policy` — subroutine
- `set_refine3D_filtering_policy` — subroutine
- `set_refine3D_gauref_policy` — subroutine
- `set_refine3D_mode_policy` — subroutine
- `set_refine3D_register_policy` — subroutine
- `set_refine3D_stage_controls` — subroutine
- `set_refine3D_symmetry_policy` — subroutine
- `set_refine3D_trailrec_policy` — subroutine
- `set_refine3D_update_policy` — subroutine
- `stage_matching_lp` — function

Local symbols (enclosing procedure scope):
- `set_lplims_from_frcs::calc_lplim_final_stage` — function

---
## Module: simple_sp_project

Files:
- `main/project/simple_sp_project.f90`
- `main/project/simple_sp_project_cls.f90`
- `main/project/simple_sp_project_core.f90`
- `main/project/simple_sp_project_io.f90`
- `main/project/simple_sp_project_mic.f90`
- `main/project/simple_sp_project_optics.f90`
- `main/project/simple_sp_project_out.f90`
- `main/project/simple_sp_project_ptcl.f90`
- `main/project/simple_sp_project_stk.f90`

Uses:
- `iso_fortran_env`
- `json_kinds`
- `json_module`
- `simple_ansi_ctrls`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs`
- `simple_defs_ori`
- `simple_discrete_stack_io`
- `simple_gui_utils`
- `simple_histogram`
- `simple_image`
- `simple_map_reduce`
- `simple_rec_list`
- `simple_stack_io`
- `simple_starfile`

Public symbols:
- `oritype2segment` — function
- `oritype2segment` — function
- `sp_project` — type

Private symbols:
- `add_cavgs2os_out` — subroutine
- `add_cavgs2os_out` — subroutine
- `add_entry2os_out` — subroutine
- `add_entry2os_out` — subroutine
- `add_frcs2os_out` — subroutine
- `add_frcs2os_out` — subroutine
- `add_fsc2os_out` — subroutine
- `add_fsc2os_out` — subroutine
- `add_intgs` — subroutine
- `add_intgs` — subroutine
- `add_movies` — subroutine
- `add_movies` — subroutine
- `add_scale_tag` — subroutine
- `add_scale_tag` — subroutine
- `add_single_movie` — subroutine
- `add_single_movie` — subroutine
- `add_single_stk` — subroutine
- `add_single_stk` — subroutine
- `add_state_weights2os_out` — subroutine
- `add_state_weights2os_out` — subroutine
- `add_stk` — subroutine
- `add_stk` — subroutine
- `add_stktab_1` — subroutine
- `add_stktab_1` — subroutine
- `add_stktab_2` — subroutine
- `add_stktab_2` — subroutine
- `add_vol2os_out` — subroutine
- `add_vol2os_out` — subroutine
- `append_job_descr2jobproc` — subroutine
- `append_job_descr2jobproc` — subroutine
- `append_project` — subroutine
- `append_project` — subroutine
- `assign` — subroutine
- `assign` — subroutine
- `cavgs2jpg` — subroutine
- `cavgs2jpg` — subroutine
- `cavgs2mrc` — subroutine
- `cavgs2mrc` — subroutine
- `clear_segment_inside` — subroutine
- `clear_segment_inside` — subroutine
- `copy` — subroutine
- `copy` — subroutine
- `count_state_gt_zero` — function
- `count_state_gt_zero` — function
- `ensure_phase_shift_fields` — subroutine
- `get_all_fscs` — subroutine
- `get_all_fscs` — subroutine
- `get_all_vols` — subroutine
- `get_all_vols` — subroutine
- `get_box` — function
- `get_box` — function
- `get_boxcoords` — subroutine
- `get_boxcoords` — subroutine
- `get_cavgs_stk` — subroutine
- `get_cavgs_stk` — subroutine
- `get_ctfflag` — function
- `get_ctfflag` — function
- `get_ctfflag_type` — function
- `get_ctfflag_type` — function
- `get_ctfparams` — function
- `get_ctfparams` — function
- `get_frcs` — subroutine
- `get_frcs` — subroutine
- `get_fsc` — subroutine
- `get_fsc` — subroutine
- `get_imgdims_from_osout` — subroutine
- `get_imgdims_from_osout` — subroutine
- `get_imginfo_from_osout` — subroutine
- `get_imginfo_from_osout` — subroutine
- `get_mic2stk_inds` — subroutine
- `get_mic2stk_inds` — subroutine
- `get_mic_kind` — function
- `get_mic_kind` — function
- `get_micname` — function
- `get_micname` — function
- `get_micparams` — function
- `get_micparams` — function
- `get_mics_table` — subroutine
- `get_mics_table` — subroutine
- `get_movies_table` — subroutine
- `get_movies_table` — subroutine
- `get_mskdiam` — subroutine
- `get_mskdiam` — subroutine
- `get_n_insegment` — function
- `get_n_insegment` — function
- `get_n_insegment_state` — function
- `get_n_insegment_state` — function
- `get_nframes` — function
- `get_nframes` — function
- `get_nintgs` — function
- `get_nintgs` — function
- `get_nmovies` — function
- `get_nmovies` — function
- `get_nptcls` — function
- `get_nptcls` — function
- `get_nstks` — function
- `get_nstks` — function
- `get_os_out_entry_index` — function
- `get_os_out_entry_index` — function
- `get_rnd_seed` — subroutine
- `get_rnd_seed` — subroutine
- `get_selected_clsinds` — function
- `get_selected_clsinds` — function
- `get_sigma2_state_path` — subroutine
- `get_sigma2_state_path` — subroutine
- `get_smpd` — function
- `get_smpd` — function
- `get_sp_oris` — subroutine
- `get_sp_oris` — subroutine
- `get_state_weights` — subroutine
- `get_state_weights` — subroutine
- `get_stkname` — function
- `get_stkname` — function
- `get_stkname_and_ind` — subroutine
- `get_stkname_and_ind` — subroutine
- `get_vol` — subroutine
- `get_vol` — subroutine
- `has_boxfile` — function
- `has_boxfile` — function
- `has_valid_ptcl3D_alignment` — function
- `has_valid_ptcl3D_alignment` — function
- `import_optics_map` — subroutine
- `import_optics_map` — subroutine
- `is_virgin_field` — function
- `is_virgin_field` — function
- `isthere_in_osout` — function
- `isthere_in_osout` — function
- `kill` — subroutine
- `kill` — subroutine
- `map2Dshifts23D` — subroutine
- `map2Dshifts23D` — subroutine
- `map2ptcls` — subroutine
- `map2ptcls` — subroutine
- `map2ptcls_state` — subroutine
- `map2ptcls_state` — subroutine
- `map_cavgs_selection` — subroutine
- `map_cavgs_selection` — subroutine
- `map_cls2D_flag_to_ptcls` — subroutine
- `map_cls2D_flag_to_ptcls` — subroutine
- `map_cls3D_states2ptcls` — subroutine
- `map_cls3D_states2ptcls` — subroutine
- `map_ptcl_ind2stk_ind` — subroutine
- `map_ptcl_ind2stk_ind` — subroutine
- `map_ptcls_state_to_cls` — subroutine
- `map_ptcls_state_to_cls` — subroutine
- `merge_algndocs` — subroutine
- `merge_algndocs` — subroutine
- `new_seg_with_ptr` — subroutine
- `new_seg_with_ptr` — subroutine
- `print_info` — subroutine
- `print_info` — subroutine
- `print_info_json` — subroutine
- `print_info_json` — subroutine
- `print_segment` — subroutine
- `print_segment` — subroutine
- `print_segment_json` — subroutine
- `print_segment_json` — subroutine
- `projrecords2proj` — subroutine
- `projrecords2proj` — subroutine
- `prune_particles` — subroutine
- `prune_particles` — subroutine
- `ptr2oritype` — subroutine
- `ptr2oritype` — subroutine
- `read` — subroutine
- `read` — subroutine
- `read_ctfparams_state_eo` — subroutine
- `read_ctfparams_state_eo` — subroutine
- `read_data_info` — subroutine
- `read_data_info` — subroutine
- `read_mic_stk_ptcl2D_segments` — subroutine
- `read_mic_stk_ptcl2D_segments` — subroutine
- `read_non_data_segments` — subroutine
- `read_non_data_segments` — subroutine
- `read_segment` — subroutine
- `read_segment` — subroutine
- `read_segments_info` — subroutine
- `read_segments_info` — subroutine
- `refuse_worker_write` — subroutine
- `remove_entry_from_osout` — subroutine
- `remove_entry_from_osout` — subroutine
- `remove_state_artifacts_from_osout` — subroutine
- `remove_state_artifacts_from_osout` — subroutine
- `remove_state_weights_from_osout` — subroutine
- `remove_state_weights_from_osout` — subroutine
- `replace_project` — subroutine
- `replace_project` — subroutine
- `report_state2mic` — subroutine
- `report_state2mic` — subroutine
- `report_state2stk` — subroutine
- `report_state2stk` — subroutine
- `scale_projfile` — subroutine
- `scale_projfile` — subroutine
- `segment2info` — function
- `segment2oritype` — function
- `segreader` — subroutine
- `segreader` — subroutine
- `segwriter` — subroutine
- `segwriter` — subroutine
- `segwriter_inside` — subroutine
- `segwriter_inside` — subroutine
- `set_boxcoords` — subroutine
- `set_boxcoords` — subroutine
- `set_boxfile` — subroutine
- `set_boxfile` — subroutine
- `set_cavgs_thumb` — subroutine
- `set_cavgs_thumb` — subroutine
- `set_ptcl2D_thumb` — subroutine
- `set_ptcl2D_thumb` — subroutine
- `set_rnd_seed` — subroutine
- `set_rnd_seed` — subroutine
- `set_sigma2_state_path` — subroutine
- `set_sigma2_state_path` — subroutine
- `set_sp_oris` — subroutine
- `set_sp_oris` — subroutine
- `shape_ranked_cavgs2jpg` — subroutine
- `shape_ranked_cavgs2jpg` — subroutine
- `sigma2_state_project_dir` — function
- `split_stk` — subroutine
- `split_stk` — subroutine
- `update_compenv` — subroutine
- `update_compenv` — subroutine
- `update_projinfo_1` — subroutine
- `update_projinfo_1` — subroutine
- `update_projinfo_2` — subroutine
- `update_projinfo_2` — subroutine
- `write` — subroutine
- `write` — subroutine
- `write_mics_star` — subroutine
- `write_mics_star` — subroutine
- `write_non_data_segments` — subroutine
- `write_non_data_segments` — subroutine
- `write_optics_map` — subroutine
- `write_optics_map` — subroutine
- `write_ptcl2D_star` — subroutine
- `write_ptcl2D_star` — subroutine
- `write_segment_inside` — subroutine
- `write_segment_inside` — subroutine
- `write_substk_1` — subroutine
- `write_substk_1` — subroutine
- `write_substk_2` — subroutine
- `write_substk_2` — subroutine

Local symbols (enclosing procedure scope):
- `merge_algndocs::register_doc_range` — subroutine
- `print_segment_json::add_fsc` — subroutine
- `print_segment_json::add_oriplot` — subroutine
- `print_segment_json::calculate_histogram` — subroutine
- `print_segment_json::calculate_indices` — subroutine
- `print_segment_json::calculate_optics_plot` — subroutine
- `print_segment_json::calculate_plot` — subroutine
- `print_segment_json::get_projections` — subroutine
- `print_segment_json::sort_oris` — function

---
## Module: simple_sp_project_tester

Files:
- `main/project/simple_sp_project_tester.f90`

Uses:
- `iso_fortran_env`
- `json_kinds`
- `json_module`
- `simple_binoris_io`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_forked_process`
- `simple_image`
- `simple_sp_project`
- `simple_test_utils`

Public symbols:
- `run_all_sp_project_tests` — subroutine

Private symbols:
- `assert_phase_at` — subroutine
- `assert_segments_equal` — subroutine
- `cleanup` — subroutine
- `make_project` — subroutine
- `range_reader` — type
- `read_mic_range` — subroutine
- `test_map_cls3D_states2ptcls` — subroutine
- `test_partial_docs_merge` — subroutine
- `test_partial_docs_merge_uneven` — subroutine
- `test_phase_records` — subroutine
- `test_phase_serialisation` — subroutine
- `test_print_segment_json` — subroutine
- `test_prune_particles` — subroutine
- `test_read_does_not_mutate_projinfo` — subroutine
- `test_read_segment_range` — subroutine
- `test_rnd_seed_roundtrip` — subroutine
- `test_worker_writes_no_project` — subroutine
- `test_write_read_roundtrip` — subroutine
- `worker_writer` — type
- `worker_writes` — subroutine
- `write_dummy_stack` — subroutine
- `zero_phase_ctf` — function

Local symbols (enclosing procedure scope):
- `test_print_segment_json::dump_json` — subroutine

---
## Module: simple_srch_sort_loc

Files:
- `utils/math/simple_srch_sort_loc.f90`

Uses:
- `simple_defs`
- `simple_error`

Public symbols:
- `dists2order` — function
- `find_int` — subroutine
- `find_real` — subroutine
- `hpsort_int` — subroutine
- `hpsort_int_cmp` — subroutine
- `hpsort_int_inds` — subroutine
- `hpsort_real` — subroutine
- `hpsort_real_inds` — subroutine
- `hpsort_real_real` — subroutine
- `locate_int` — function
- `locate_real` — function
- `mask2inds` — function
- `maxnloc` — function
- `min3` — function
- `minnloc` — function
- `peakfinder` — function
- `reorder_int` — subroutine
- `reorder_real` — subroutine
- `reverse_drarr` — subroutine
- `reverse_f` — subroutine
- `reverse_iarr` — subroutine
- `reverse_rarr` — subroutine
- `scores2order` — function
- `selec_dp` — function
- `selec_sp` — function
- `unique` — subroutine

Local symbols (enclosing procedure scope):
- `hpsort_int_cmp::p1_lt_p2` — function

---
## Module: simple_srch_sort_loc_tester

Files:
- `utils/math/simple_srch_sort_loc_tester.f90`

Uses:
- `simple_defs`
- `simple_srch_sort_loc`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_srch_sort_loc_tests` — subroutine

Private symbols:
- `descending` — function
- `is_permutation` — function
- `sorted_ints` — function
- `sorted_reals` — function
- `test_find` — subroutine
- `test_hpsort_edge_sizes` — subroutine
- `test_hpsort_integers` — subroutine
- `test_hpsort_reals` — subroutine
- `test_hpsort_with_comparator` — subroutine
- `test_locate` — subroutine
- `test_mask2inds_and_reorder` — subroutine
- `test_maxnloc_minnloc` — subroutine
- `test_min3` — subroutine
- `test_orderings` — subroutine
- `test_peakfinder` — subroutine
- `test_reverse` — subroutine
- `test_selec` — subroutine
- `test_unique` — subroutine

---
## Module: simple_srchspace_map

Files:
- `main/strategies/search/simple_srchspace_map.f90`

Uses:
- `simple_error`
- `simple_srch_sort_loc`

Public symbols:
- `srchspace_map` — type

Private symbols:
- `full2sub` — function
- `get_full2sub_map` — function
- `get_inds_in_full` — function
- `get_nspace` — function
- `get_nsub` — function
- `get_sub2full_map` — function
- `has_data` — function
- `kill` — subroutine
- `new_from_distmat` — subroutine
- `new_from_maps` — subroutine
- `sub2full` — function

---
## Module: simple_stack_io

Files:
- `fileio/simple_stack_io.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_imgfile`
- `simple_imghead`

Public symbols:
- `stack_io` — type

Private symbols:
- `close` — subroutine
- `get_image` — subroutine
- `get_ldim` — function
- `get_nptcls` — function
- `open_box` — subroutine
- `open_ldim` — subroutine
- `read` — subroutine
- `read_whole` — subroutine
- `same_stk` — function
- `stk_is_open` — function
- `write` — subroutine
- `write_buffer` — subroutine

---
## Module: simple_stack_io_tester

Files:
- `fileio/simple_stack_io_tester.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `simple_cmdline`
- `simple_discrete_stack_io`
- `simple_forked_process`
- `simple_image`
- `simple_image_slab`
- `simple_imghead`
- `simple_stack_io`
- `simple_stack_read_session`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_stack_io_tests` — subroutine

Private symbols:
- `assert_f16_payload` — subroutine
- `assert_mrc_header` — subroutine
- `assert_pattern` — subroutine
- `check_pixel_value` — subroutine
- `check_session_against_single_reads` — subroutine
- `check_span_counts` — subroutine
- `cleanup` — subroutine
- `create_synthetic_stack` — subroutine
- `delete_stacks` — subroutine
- `expect_child_stops` — subroutine
- `failing_reader` — type
- `fill_chunk_pattern` — subroutine
- `fill_pattern` — subroutine
- `fill_with_junk` — subroutine
- `read_and_fail` — subroutine
- `read_f16_first_row` — subroutine
- `read_stacks_concurrently` — subroutine
- `session_pixel` — function
- `session_stack_name` — function
- `set_f16_rounding_pattern` — subroutine
- `set_pattern` — subroutine
- `test_buffer_sizes_and_forward_skips` — subroutine
- `test_dstack_parallel_read_float32` — subroutine
- `test_dstack_parallel_read_int16` — subroutine
- `test_float16_chunked_roundtrip` — subroutine
- `test_float16_encoder_boundaries` — subroutine
- `test_float16_header_and_roundtrip` — subroutine
- `test_float16_image_roundtrip` — subroutine
- `test_float16_rounding_and_payload` — subroutine
- `test_float32_header` — subroutine
- `test_image_reader_to_stack_io_writer` — subroutine
- `test_open_state` — subroutine
- `test_posix_descriptor_exhaustion_stops` — subroutine
- `test_posix_extended_header` — subroutine
- `test_posix_fallback_mode` — subroutine
- `test_posix_missing_and_truncated_stop` — subroutine
- `test_posix_offsets_beyond_4GiB` — subroutine
- `test_posix_parallel_opens_no_leak` — subroutine
- `test_posix_short_reads_eintr_and_advice` — subroutine
- `test_read_image_written_stack` — subroutine
- `test_session_bitexact` — subroutine
- `test_session_dimension_mismatch_stops` — subroutine
- `test_session_gap_budget` — subroutine
- `test_session_many_stacks` — subroutine
- `test_session_one_open_per_stack` — subroutine
- `test_session_prefetch` — subroutine
- `test_session_spans` — subroutine
- `test_slab_budget_and_failed_allocation` — subroutine
- `test_slab_float32_nan_payload` — subroutine
- `test_slab_get_many_threads` — subroutine
- `test_slab_roundtrip_reuse_and_immutable` — subroutine
- `test_stack_io_copy` — subroutine
- `verify_stack` — subroutine
- `write_int16_stack` — subroutine
- `write_raw_float_stack` — subroutine
- `write_raw_integer_stack` — subroutine
- `write_session_stacks` — subroutine
- `wrong_guess` — function

---
## Module: simple_stack_read_session

Files:
- `fileio/simple_stack_read_session.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `simple_core_module_api`
- `simple_discrete_stack_io`
- `simple_image`
- `simple_imgfile`
- `simple_imghead`
- `simple_syslib`

Public symbols:
- `plan_spans` — subroutine
- `stack_read_counts` — type
- `stack_read_session` — type

Private symbols:
- `decode_image` — subroutine
- `find_held` — function
- `get_counts` — function
- `group_requests` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `open_stack` — type
- `open_stack_file` — subroutine
- `prefetch` — subroutine
- `raw_scratch` — type
- `read` — subroutine
- `reads_by_descriptor` — function
- `release_held` — subroutine
- `setup_streams` — subroutine
- `span_cap` — function
- `window_span` — type

Local symbols (enclosing procedure scope):
- `plan_spans::open_span` — subroutine

---
## Module: simple_stackops

Files:
- `main/image_processing/simple_stackops.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `acf_stack` — subroutine
- `frameavg_stack` — subroutine
- `make_avg_stack` — subroutine
- `make_pcavec_stack` — subroutine
- `prep_imgfile4movie` — subroutine
- `stats_imgfile` — subroutine

Private symbols:
- `raise_exception` — subroutine

---
## Module: simple_starfile

Files:
- `main/star/simple_starfile.f90`

Uses:
- `simple_core_module_api`
- `simple_starfile_wrappers`

Public symbols:
- `starfile` — type

Private symbols:
- `calc_part_boundaries` — subroutine
- `check_origin_smpd` — subroutine
- `complete` — subroutine
- `has_origin` — function
- `init` — subroutine
- `write_mics_table` — subroutine
- `write_optics_table` — subroutine
- `write_ptcl2D_table` — subroutine
- `write_ptcl2D_table_parallel` — subroutine

Local symbols (enclosing procedure scope):
- `write_mics_table::get_relative_path_here` — function
- `write_ptcl2D_table::get_relative_path_here` — function
- `write_ptcl2D_table_parallel::get_relative_path_here` — function

---
## Module: simple_starfile_tester

Files:
- `main/star/simple_starfile_tester.f90`

Uses:
- `omp_lib`
- `simple_core_module_api`
- `simple_starfile`
- `simple_starfile_wrappers`
- `simple_test_utils`

Public symbols:
- `run_all_starfile_tests` — subroutine

Private symbols:
- `count_objects` — function
- `fill_ptcl_stk` — subroutine
- `make_oris_n` — subroutine
- `readback_double` — function
- `readback_int` — function
- `readback_origins` — subroutine
- `setup_tmpdir` — subroutine
- `test_complete_noop_without_tmp` — subroutine
- `test_complete_overwrites_existing` — subroutine
- `test_complete_renames_tmp` — subroutine
- `test_init_creates_no_files` — subroutine
- `test_init_removes_stale_tmp` — subroutine
- `test_table_wrappers_roundtrip` — subroutine
- `test_verbose_flag_no_content_change` — subroutine
- `test_write_mics_table_all_fields` — subroutine
- `test_write_mics_table_basic` — subroutine
- `test_write_mics_table_empty` — subroutine
- `test_write_mics_table_state_filter` — subroutine
- `test_write_optics_table_all_fields` — subroutine
- `test_write_optics_table_basic` — subroutine
- `test_write_optics_table_empty` — subroutine
- `test_write_optics_table_state_filter` — subroutine
- `test_write_ptcl2D_origins_angstrom` — subroutine
- `test_write_ptcl2D_parallel_basic` — subroutine
- `test_write_ptcl2D_parallel_count` — subroutine
- `test_write_ptcl2D_parallel_empty` — subroutine
- `test_write_ptcl2D_parallel_vs_serial` — subroutine
- `test_write_ptcl2D_table_bad_stkind` — subroutine
- `test_write_ptcl2D_table_basic` — subroutine
- `test_write_ptcl2D_table_indstk_explicit` — subroutine
- `test_write_ptcl2D_table_state_filter` — subroutine
- `test_write_ptcl2D_table_with_mics` — subroutine
- `tmp` — function
- `write_dummy_file` — subroutine

---
## Module: simple_starfile_wrappers

Files:
- `main/star/simple_starfile_wrapper.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `C_dealloc_str` — subroutine
- `C_print_pointer` — subroutine
- `C_starfile_table__addObject` — subroutine
- `C_starfile_table__append` — subroutine
- `C_starfile_table__clear` — subroutine
- `C_starfile_table__close_ofile` — subroutine
- `C_starfile_table__delete` — subroutine
- `C_starfile_table__firstobject` — function
- `C_starfile_table__getcomment` — subroutine
- `C_starfile_table__getnames_cnt` — subroutine
- `C_starfile_table__getnames_nr` — subroutine
- `C_starfile_table__getValue_bool` — function
- `C_starfile_table__getValue_double` — function
- `C_starfile_table__getValue_float` — function
- `C_starfile_table__getValue_int` — function
- `C_starfile_table__getValue_string` — subroutine
- `C_starfile_table__hascomment` — function
- `C_starfile_table__hasLabel` — function
- `C_starfile_table__new` — function
- `C_starfile_table__nextobject` — function
- `C_starfile_table__numberofobjects` — function
- `C_starfile_table__open_ofile` — subroutine
- `C_starfile_table__read` — subroutine
- `C_starfile_table__setcomment` — subroutine
- `C_starfile_table__setIsList` — subroutine
- `C_starfile_table__setname` — subroutine
- `C_starfile_table__setValue_bool` — subroutine
- `C_starfile_table__setValue_double` — subroutine
- `C_starfile_table__setValue_float` — subroutine
- `C_starfile_table__setValue_int` — subroutine
- `C_starfile_table__setValue_string` — subroutine
- `C_starfile_table__write_ofile` — subroutine
- `C_starfile_table__write_omem` — subroutine
- `starfile_table__addObject` — subroutine
- `starfile_table__append` — subroutine
- `starfile_table__clear` — subroutine
- `starfile_table__close_ofile` — subroutine
- `starfile_table__delete` — subroutine
- `starfile_table__firstobject` — function
- `starfile_table__getcomment` — subroutine
- `starfile_table__getnames` — subroutine
- `starfile_table__getValue_bool` — function
- `starfile_table__getValue_double` — function
- `starfile_table__getValue_float` — function
- `starfile_table__getValue_int` — function
- `starfile_table__getValue_string` — function
- `starfile_table__hascomment` — function
- `starfile_table__haslabel` — function
- `starfile_table__new` — subroutine
- `starfile_table__nextobject` — function
- `starfile_table__numberofobjects` — function
- `starfile_table__open_ofile` — subroutine
- `starfile_table__read` — subroutine
- `starfile_table__setcomment` — subroutine
- `starfile_table__setIsList` — subroutine
- `starfile_table__setname` — subroutine
- `starfile_table__setValue_bool` — subroutine
- `starfile_table__setValue_double` — subroutine
- `starfile_table__setValue_float` — subroutine
- `starfile_table__setValue_int` — subroutine
- `starfile_table__setValue_string` — subroutine
- `starfile_table__write_ofile` — subroutine
- `starfile_table__write_omem` — subroutine
- `starfile_table_type` — type

---
## Module: simple_starproject

Files:
- `main/star/simple_starproject.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_epu_xml`
- `simple_refine2d_fnames`
- `simple_rnd`
- `simple_sp_project`
- `simple_starproject_utils`

Public symbols:
- `starproject` — type

Private symbols:
- `apply_optics_offset` — subroutine
- `assign_initial_tiltgroups` — subroutine
- `assign_optics` — subroutine
- `assign_xml_tiltinfo` — subroutine
- `cluster_tiltinfo` — subroutine
- `coords_to_corners` — subroutine
- `export_cls2D` — subroutine
- `export_mics` — subroutine
- `export_opticsgroups` — subroutine
- `export_stardata` — subroutine
- `finalise_ptcl_import` — subroutine
- `get_image_basename` — subroutine
- `import_cls2D` — subroutine
- `import_mics` — subroutine
- `import_ptcls` — subroutine
- `import_ptcls2D` — subroutine
- `import_ptcls3D` — subroutine
- `import_stardata` — subroutine
- `initialise` — subroutine
- `kill` — subroutine
- `origins_to_pixels` — subroutine
- `plot_opticsgroups` — subroutine
- `populate_opticsmap` — subroutine
- `populate_stkmap` — subroutine
- `read_layout` — subroutine
- `read_starheaders` — subroutine
- `set_verbose` — subroutine
- `sort_optics_maxpop` — subroutine

Local symbols (enclosing procedure scope):
- `coords_to_corners::field_to_corners` — subroutine
- `origins_to_pixels::field_to_pixels` — subroutine

---
## Module: simple_starproject_stream

Files:
- `main/star/simple_starproject_stream.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_parameters`
- `simple_sp_project`
- `simple_starfile_wrappers`

Public symbols:
- `starproject_stream` — type

Private symbols:
- `assign_optics_single` — subroutine
- `copy_optics` — subroutine
- `starfile_deinit` — subroutine
- `starfile_init` — subroutine
- `starfile_set_micrographs_table` — subroutine
- `starfile_set_optics_group_table` — subroutine
- `starfile_set_optics_table` — subroutine
- `starfile_set_particles2D_table` — subroutine
- `starfile_write_table` — subroutine
- `stream_export_micrographs` — subroutine
- `stream_export_particles_2D` — subroutine
- `stream_write_optics` — subroutine

Local symbols (enclosing procedure scope):
- `starfile_set_micrographs_table::get_relative_path_here` — function
- `starfile_set_particles2D_table::get_relative_path_here` — function

---
## Module: simple_starproject_stream_tester

Files:
- `main/star/simple_starproject_stream_tester.f90`

Uses:
- `simple_defs`
- `simple_fileio`
- `simple_parameters`
- `simple_sp_project`
- `simple_starfile_wrappers`
- `simple_starproject_stream`
- `simple_string`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_starproject_stream_tests` — subroutine

Private symbols:
- `ends_with` — function
- `field_with` — function
- `is_blank` — function
- `is_image` — function
- `make_two_stack_project` — subroutine
- `read_particle_rows` — subroutine
- `test_particles2D_names` — subroutine
- `test_particles2D_origins` — subroutine

---
## Module: simple_starproject_tester

Files:
- `main/star/simple_starproject_tester.f90`

Uses:
- `omp_lib`
- `simple_cmdline`
- `simple_defs`
- `simple_fileio`
- `simple_forked_process`
- `simple_image`
- `simple_oris`
- `simple_relion`
- `simple_sp_project`
- `simple_starfile`
- `simple_starfile_wrappers`
- `simple_starproject`
- `simple_starproject_utils`
- `simple_string`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_starproject_tests` — subroutine

Private symbols:
- `check_origins_px` — subroutine
- `create_directory` — subroutine
- `delete_directory` — subroutine
- `force_openmp_threads` — subroutine
- `make_origins_stacks` — subroutine
- `make_tiltinfo` — function
- `read_particle_doubles` — subroutine
- `run_star_import_child` — subroutine
- `setup_tmpdir` — subroutine
- `star_import_child` — type
- `test_export_cls2D_keeps_project` — subroutine
- `test_finalise_once` — subroutine
- `test_import_origins_angstrom` — subroutine
- `test_large_particle_table_export` — subroutine
- `test_optics_clustering_basic` — subroutine
- `test_particle_export_readback` — subroutine
- `test_relion30_layout_refused` — subroutine
- `test_relion_allocate_opticsgroups` — subroutine
- `test_relion_epu_tiltgroups` — subroutine
- `test_relion_export_poses` — subroutine
- `test_relion_find_movienames` — subroutine
- `test_relion_phase_shift` — subroutine
- `test_relion_single_tiltgroup` — subroutine
- `test_relion_write_corrected_micrographs_star` — subroutine
- `test_relion_write_micrographs_star` — subroutine
- `test_relion_write_particles2D_star` — subroutine
- `test_relion_writer_micrographs` — subroutine
- `test_resolve_stack_path` — subroutine
- `test_roundtrip_micrographs` — subroutine
- `test_roundtrip_origins` — subroutine
- `test_split_dataline_terminal_field` — subroutine
- `test_star_export_micrographs` — subroutine
- `test_star_micrographs_readback` — subroutine
- `test_starfile_append_mode` — subroutine
- `test_starfile_basic_init` — subroutine
- `test_starfile_multiple_tables` — subroutine
- `test_starfile_write_and_readback` — subroutine
- `test_write_omem_parallel_stability` — subroutine
- `test_xml_tiltinfo` — subroutine
- `tmpfile` — function
- `write_origins_star` — subroutine
- `write_textfile` — subroutine

---
## Module: simple_starproject_utils

Files:
- `main/star/simple_starproject_utils.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `enable_rlnflag` — subroutine
- `enable_splflag` — subroutine
- `enable_splflags` — subroutine
- `find_separators` — subroutine
- `get_value_from_ptcls` — function
- `h_clust` — subroutine
- `resolve_stack_path` — subroutine
- `split_dataline` — subroutine
- `star_data` — type
- `star_file` — type
- `star_flag` — type
- `tilt_info` — type

---
## Module: simple_stat

Files:
- `utils/math/simple_stat.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_is_check_assert`
- `simple_rnd`
- `simple_srch_sort_loc`
- `simple_type_defs`

Public symbols:
- `analyze_smat` — subroutine
- `avg_frac_smallest` — function
- `avg_sdev_1d` — subroutine
- `avg_sdev_1d_dp` — subroutine
- `avg_sdev_2d` — subroutine
- `avg_sdev_3d` — subroutine
- `calc_ap_pref` — function
- `calc_stats` — subroutine
- `conv2rank_weights` — subroutine
- `corr2distweight` — function
- `corrs2weights` — function
- `dists2scores_percen` — subroutine
- `dmat2smat` — function
- `kish_ess` — function
- `kstwo` — subroutine
- `kurtosis_1d` — function
- `kurtosis_2d` — function
- `logsumexp` — function
- `mad_dp` — function
- `mad_gau` — function
- `mad_sp` — function
- `median_dp` — function
- `median_nocopy` — function
- `median_sp` — function
- `medoid_from_dmat` — subroutine
- `medoid_from_smat` — subroutine
- `medoid_ranking_from_smat` — subroutine
- `merge_dmats_pair` — function
- `merge_dmats_triple` — function
- `merge_smats` — function
- `moment_1d` — subroutine
- `moment_1d_masked` — subroutine
- `moment_2d` — subroutine
- `moment_3d` — subroutine
- `moment_serial` — subroutine
- `neighfrac2nsmpl` — function
- `norm_corr` — function
- `norm_corr_8` — function
- `normalize_1d` — subroutine
- `normalize_1d_masked` — subroutine
- `normalize_2d` — subroutine
- `normalize_3d` — subroutine
- `normalize_minmax_1d` — subroutine
- `normalize_minmax_2d` — subroutine
- `normalize_sigm_1d` — subroutine
- `normalize_sigm_2d` — subroutine
- `normalize_sigm_3d` — subroutine
- `pearsn_1d` — function
- `pearsn_2d` — function
- `pearsn_3d` — function
- `pearsn_serial` — function
- `pearsn_serial_8` — function
- `power_sampling` — subroutine
- `prepare_tree_sub_distmat` — subroutine
- `rank_centroid_weights` — subroutine
- `rank_exponent_weights` — subroutine
- `rank_inverse_weights` — subroutine
- `rank_sum_weights` — subroutine
- `robust_normalization` — subroutine
- `robust_normalize_minmax` — subroutine
- `robust_scaling` — subroutine
- `robust_sigma_thres` — function
- `robust_z_scores` — function
- `scores2scores_percen` — subroutine
- `skewness_1d` — function
- `skewness_2d` — function
- `smat2dmat` — function
- `spear` — function
- `squared_sampling` — subroutine
- `std_mean_diff` — function
- `z_scores` — function

Local symbols (enclosing procedure scope):
- `kstwo::probks` — function

---
## Module: simple_stat_tester

Files:
- `utils/math/simple_stat_tester.f90`

Uses:
- `simple_defs`
- `simple_stat`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_stat_tests` — subroutine

Private symbols:
- `test_calc_stats` — subroutine
- `test_conv2rank_weights` — subroutine
- `test_corrs2weights_corrw` — subroutine
- `test_corrs2weights_other_criteria` — subroutine
- `test_logsumexp` — subroutine
- `test_median` — subroutine
- `test_rank_weight_kernels` — subroutine

Local symbols (enclosing procedure scope):
- `test_rank_weight_kernels::check_kernel` — subroutine

---
## Module: simple_state_assignment

Files:
- `main/solve/simple_state_assignment.f90`

Uses:
- `ieee_arithmetic`
- `simple_cavg_quality_model`
- `simple_cavg_quality_selection`
- `simple_cavg_quality_types`
- `simple_cmdline`
- `simple_defs`
- `simple_defs_fname`
- `simple_dock_vols`
- `simple_eo_partition`
- `simple_error`
- `simple_fileio`
- `simple_gui_utils`
- `simple_image`
- `simple_imgarr_utils`
- `simple_imghead`
- `simple_math`
- `simple_oris`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_sym`
- `simple_syslib`

Public symbols:
- `assign_stats` — type
- `cavgs_fbody` — function
- `check_state_refs` — subroutine
- `compact_states` — subroutine
- `count_active_classes` — function
- `cut_whole_stacks` — function
- `docked_cc_common_grid` — function
- `docked_pair_cc` — function
- `expansion_budget` — function
- `fit_mskdiam` — function
- `flag_ambiguous` — function
- `log_jpeg` — subroutine
- `map_assignment` — function
- `map_refine_result` — function
- `map_state_poses` — function
- `median_of` — function
- `model_rejected_rows` — function
- `read_refs_record` — function
- `refs_check` — type
- `refs_usable` — function
- `reject_by_fit` — subroutine
- `score_cavgs_with_model` — function
- `select_route_classes` — subroutine
- `select_state_classes` — function
- `set_assign_cline` — subroutine
- `set_cavgs3D_cline` — subroutine
- `set_clsexp_cline` — subroutine
- `set_refine_cont_cline` — subroutine
- `set_refine_seeds` — subroutine
- `set_solve2D_route_cline` — subroutine
- `set_state_cavgs3D_cline` — subroutine
- `usable_job_2D` — function
- `write_matched_cavgs` — subroutine
- `write_refine_project` — subroutine
- `write_state_cavgs` — subroutine

Local symbols (enclosing procedure scope):
- `check_state_refs::write_check_cavgs` — subroutine
- `docked_cc_common_grid::fourier_crop` — subroutine
- `map_assignment::log_fit` — subroutine
- `write_state_cavgs::write_cavgs_jpeg` — subroutine

---
## Module: simple_state_assignment_tester

Files:
- `main/solve/simple_state_assignment_tester.f90`

Uses:
- `simple_cmdline`
- `simple_fileio`
- `simple_image`
- `simple_imghead`
- `simple_sp_project`
- `simple_state_assignment`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_state_assignment_tests` — subroutine

Private symbols:
- `test_check_rules` — subroutine
- `test_cut_whole_stacks` — subroutine
- `test_fit_rules` — subroutine
- `test_job_clines` — subroutine
- `test_state_cavgs` — subroutine
- `test_state_runs` — subroutine

---
## Module: simple_state_weight_set

Files:
- `main/project/simple_state_weight_set.f90`

Uses:
- `ieee_arithmetic`
- `iso_fortran_env`
- `simple_error`
- `simple_fileio`
- `simple_oris`
- `simple_ptcl_layout`
- `simple_sp_project`
- `simple_state_weights_file`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `discard_state_weight_set` — subroutine
- `state_weight_set` — type

Private symbols:
- `applied_ess` — function
- `applied_mass` — function
- `check_state` — subroutine
- `close_enough` — function
- `column_populations` — subroutine
- `get_applied_weights` — subroutine
- `get_ess` — function
- `get_generation` — function
- `get_kind` — function
- `get_labels` — subroutine
- `get_layout_digest` — function
- `get_manifest_path` — function
- `get_mass` — function
- `get_members` — subroutine
- `get_nptcls` — function
- `get_nstates` — function
- `get_parent_generation` — function
- `get_parent_layout_digest` — function
- `get_parent_state` — function
- `get_pop` — function
- `get_producer` — function
- `get_selection` — subroutine
- `get_weights` — subroutine
- `is_fractional` — function
- `kill` — subroutine
- `load_state` — subroutine
- `new` — subroutine
- `open_set` — subroutine
- `publish` — subroutine
- `publish_work_state` — subroutine
- `row_populations` — subroutine

---
## Module: simple_state_weight_set_tester

Files:
- `main/project/simple_state_weight_set_tester.f90`

Uses:
- `iso_fortran_env`
- `simple_oris`
- `simple_sp_project`
- `simple_state_weight_set`
- `simple_state_weights_file`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_state_weight_set_tests` — subroutine

Private symbols:
- `make_project` — subroutine
- `partition_table` — subroutine
- `test_crash_between_files_and_manifest` — subroutine
- `test_discard` — subroutine
- `test_earlier_release_entries` — subroutine
- `test_refusals` — subroutine
- `test_round_trip_kernel` — subroutine
- `test_round_trip_partition` — subroutine
- `test_work_state` — subroutine

---
## Module: simple_state_weights_file

Files:
- `fileio/simple_state_weights_file.f90`

Uses:
- `iso_fortran_env`
- `simple_fileio`
- `simple_sigma2_state_file`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `state_weights_file_checksum` — function
- `state_weights_file_header` — type
- `state_weights_fname` — function
- `state_weights_kind_id` — function
- `state_weights_kind_name` — function
- `state_weights_manifest` — type
- `state_weights_read_manifest` — subroutine
- `state_weights_read_state_file` — subroutine
- `state_weights_write_manifest` — subroutine
- `state_weights_write_state_file` — subroutine

Private symbols:
- `header_from_words` — function
- `header_words` — function

Local symbols (enclosing procedure scope):
- `state_weights_read_manifest::read_keyed` — subroutine

---
## Module: simple_strategy2D

Files:
- `main/strategies/search/simple_strategy2D.f90`

Uses:
- `simple_builder`
- `simple_oris`
- `simple_parameters`
- `simple_strategy2d_srch`

Public symbols:
- `strategy2D` — type
- `strategy2D_per_ptcl` — type

Private symbols:
- `generic_kill` — subroutine
- `generic_new` — subroutine
- `generic_srch` — subroutine

---
## Module: simple_strategy2D_alloc

Files:
- `main/strategies/search/simple_strategy2D_alloc.f90`

Uses:
- `simple_eul_prob_tab_utils`
- `simple_pftc_srch_api`

Public symbols:
- `clean_strategy2D` — subroutine
- `is_fresh_2D_start` — function
- `prep_strategy2D_batch` — subroutine
- `prep_strategy2D_glob` — subroutine
- `prep_strategy2D_thread` — subroutine
- `set_strategy2D_stoch_bound` — subroutine

Private symbols:
- `strategy2D_alloc` — type

---
## Module: simple_strategy2D_greedy

Files:
- `main/strategies/search/simple_strategy2D_greedy.f90`

Uses:
- `simple_builder`
- `simple_oris`
- `simple_parameters`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_alloc`
- `simple_strategy2d_srch`

Public symbols:
- `strategy2D_greedy` — type

Private symbols:
- `kill_greedy` — subroutine
- `new_greedy` — subroutine
- `srch_greedy` — subroutine

---
## Module: simple_strategy2D_greedy_smpl

Files:
- `main/strategies/search/simple_strategy2D_greedy_smpl.f90`

Uses:
- `simple_builder`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_alloc`
- `simple_strategy2d_srch`

Public symbols:
- `strategy2D_greedy_smpl` — type

Private symbols:
- `kill_greedy_smpl` — subroutine
- `new_greedy_smpl` — subroutine
- `srch_greedy_smpl` — subroutine

---
## Module: simple_strategy2D_inpl

Files:
- `main/strategies/search/simple_strategy2D_inpl.f90`

Uses:
- `simple_builder`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_srch`

Public symbols:
- `strategy2D_inpl` — type

Private symbols:
- `kill_inpl` — subroutine
- `new_inpl` — subroutine
- `srch_inpl` — subroutine

---
## Module: simple_strategy2D_inpl_smpl

Files:
- `main/strategies/search/simple_strategy2D_inpl_smpl.f90`

Uses:
- `simple_builder`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_alloc`
- `simple_strategy2d_srch`
- `simple_type_defs`

Public symbols:
- `strategy2D_inpl_smpl` — type

Private symbols:
- `kill_inpl_smpl` — subroutine
- `new_inpl_smpl` — subroutine
- `srch_inpl_smpl` — subroutine

---
## Module: simple_strategy2D_matcher

Files:
- `main/strategies/search/simple_strategy2D_matcher.f90`

Uses:
- `iso_fortran_env`
- `simple_bench_io`
- `simple_binoris_io`
- `simple_builder`
- `simple_classaverager`
- `simple_decay_funs`
- `simple_eo_partition`
- `simple_eul_prob_tab2d`
- `simple_eul_prob_tab_utils`
- `simple_imgarr_utils`
- `simple_matcher_2dprep`
- `simple_matcher_pftc_prep`
- `simple_matcher_ptcl_batch`
- `simple_matcher_ptcl_io`
- `simple_matcher_smpl_and_lplims`
- `simple_pftc_srch_api`
- `simple_qsys_funs`
- `simple_shiftmap_srch`
- `simple_strategy2d`
- `simple_strategy2d_alloc`
- `simple_strategy2d_greedy`
- `simple_strategy2d_greedy_smpl`
- `simple_strategy2d_inpl`
- `simple_strategy2d_inpl_smpl`
- `simple_strategy2d_prob`
- `simple_strategy2d_snhc`
- `simple_strategy2d_snhc_smpl`
- `simple_strategy2d_snhc_smpl_many`
- `simple_strategy2d_srch`
- `simple_strategy2d_tseries`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `refine2D_exec` — subroutine

Private symbols:
- `refine2D_ctrl` — type

Local symbols (enclosing procedure scope):
- `refine2D_exec::allocate_strategy_for_particle` — subroutine
- `refine2D_exec::begin_pipeline` — subroutine
- `refine2D_exec::build_batch_particles_local` — subroutine
- `refine2D_exec::cleanup_search_state` — subroutine
- `refine2D_exec::compute_neigh_frac` — subroutine
- `refine2D_exec::ensure_even_odd_partition` — subroutine
- `refine2D_exec::init_ctrl` — subroutine
- `refine2D_exec::maybe_write_bench` — subroutine
- `refine2D_exec::pipeline_batch` — subroutine
- `refine2D_exec::pipeline_context_bytes` — function
- `refine2D_exec::prepare_alignment_references` — subroutine
- `refine2D_exec::prepare_batches` — subroutine
- `refine2D_exec::prepare_class_averages_and_restoration` — subroutine
- `refine2D_exec::restore_class_averages_for_batch` — subroutine
- `refine2D_exec::sample_particles_for_update` — subroutine
- `refine2D_exec::write_empty_readout` — subroutine
- `refine2D_exec::write_orientations` — subroutine

---
## Module: simple_strategy2D_prob

Files:
- `main/strategies/search/probabilistic/simple_strategy2D_prob.f90`

Uses:
- `simple_builder`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_srch`

Public symbols:
- `strategy2D_prob` — type

Private symbols:
- `kill_prob` — subroutine
- `new_prob` — subroutine
- `srch_prob` — subroutine

---
## Module: simple_strategy2D_snhc

Files:
- `main/strategies/search/simple_strategy2D_snhc.f90`

Uses:
- `simple_builder`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_alloc`
- `simple_strategy2d_srch`

Public symbols:
- `strategy2D_snhc` — type

Private symbols:
- `kill_snhc` — subroutine
- `new_snhc` — subroutine
- `srch_snhc` — subroutine

---
## Module: simple_strategy2D_snhc_smpl

Files:
- `main/strategies/search/simple_strategy2D_snhc_smpl.f90`

Uses:
- `simple_builder`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_alloc`
- `simple_strategy2d_srch`
- `simple_type_defs`

Public symbols:
- `strategy2D_snhc_smpl` — type

Private symbols:
- `kill_snhc_smpl` — subroutine
- `new_snhc_smpl` — subroutine
- `srch_snhc_smpl` — subroutine

---
## Module: simple_strategy2D_snhc_smpl_many

Files:
- `main/strategies/search/simple_strategy2D_snhc_smpl_many.f90`

Uses:
- `simple_builder`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_alloc`
- `simple_strategy2d_srch`
- `simple_type_defs`

Public symbols:
- `strategy2D_snhc_smpl_many` — type

Private symbols:
- `kill_snhc_smpl_many` — subroutine
- `new_snhc_smpl_many` — subroutine
- `srch_snhc_smpl_many` — subroutine

---
## Module: simple_strategy2D_srch

Files:
- `main/strategies/search/simple_strategy2D_srch.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab2d`
- `simple_eul_prob_tab_utils`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_rnd`
- `simple_shiftmap_srch`
- `simple_strategy2d_alloc`

Public symbols:
- `strategy2D_spec` — type
- `strategy2D_srch` — type

Private symbols:
- `assign_ori` — subroutine
- `continuous_pose_committed` — function
- `get_continuous_route_status` — subroutine
- `inpl_srch` — subroutine
- `inpl_srch_first` — subroutine
- `inpl_srch_peaks` — subroutine
- `kill` — subroutine
- `map_srch` — subroutine
- `new` — subroutine
- `prep4srch` — subroutine
- `refine_selected_continuously` — subroutine
- `store_continuous_e3` — subroutine
- `store_solution` — subroutine
- `uses_continuous_refinement` — function

Local symbols (enclosing procedure scope):
- `inpl_srch::commit_mapped_class` — subroutine
- `prep4srch::pick_populated_class` — subroutine

---
## Module: simple_strategy2D_tseries

Files:
- `main/strategies/search/simple_strategy2D_tseries.f90`

Uses:
- `simple_builder`
- `simple_pftc_srch_api`
- `simple_strategy2d`
- `simple_strategy2d_alloc`
- `simple_strategy2d_srch`

Public symbols:
- `strategy2D_tseries` — type

Private symbols:
- `kill_tseries` — subroutine
- `new_tseries` — subroutine
- `srch_tseries` — subroutine

Local symbols (enclosing procedure scope):
- `srch_tseries::per_ref_srch` — subroutine

---
## Module: simple_strategy2D_utils

Files:
- `main/strategies/search/simple_strategy2D_utils.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_histogram`
- `simple_image_bin`
- `simple_image_msk`
- `simple_parameters`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_pspecs`
- `simple_segmentation`
- `simple_sp_project`

Public symbols:
- `align_and_score_cavg_clusters` — function
- `calc_cavg_offset` — subroutine
- `calc_cavg_pairwise_algninfo` — function
- `calc_cavg_sigstats_components` — subroutine
- `calc_cluster_cavgs_dmat` — function
- `calc_match_cavgs_dmat` — function
- `cavg_sigstats_matrices` — type
- `flag_non_junk_cavgs` — subroutine
- `id_junk` — subroutine
- `id_junk_and_prep_cavgs4clust` — subroutine
- `match_imgs` — function
- `match_imgs2ref` — function
- `prep_cavgs4clust` — subroutine
- `rtsq_imgs` — subroutine
- `write_aligned_cavgs` — subroutine

Private symbols:
- `align_clusters2medoids` — function
- `calc_cc_and_res_dmats` — subroutine
- `calc_cc_and_res_dmats_ref` — subroutine
- `calc_sigstats_dmats` — subroutine
- `calc_sigstats_dmats_ref` — subroutine
- `kill_cavg_sigstats_matrices` — subroutine

Local symbols (enclosing procedure scope):
- `align_and_score_cavg_clusters::calc_scores` — subroutine
- `align_and_score_cavg_clusters::identify_good_bad_clusters` — subroutine
- `align_and_score_cavg_clusters::rank_clusters` — subroutine

---
## Module: simple_strategy3D

Files:
- `main/strategies/search/simple_strategy3D.f90`

Uses:
- `simple_builder`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d_srch`

Public symbols:
- `strategy3D` — type

Private symbols:
- `generic_kill` — subroutine
- `generic_new` — subroutine
- `generic_oris_assign` — subroutine
- `generic_srch` — subroutine

---
## Module: simple_strategy3D_alloc

Files:
- `main/strategies/search/simple_strategy3D_alloc.f90`

Uses:
- `simple_builder`
- `simple_eul_prob_tab_utils`
- `simple_pftc_srch_api`

Public symbols:
- `clean_strategy3D` — subroutine
- `prep_strategy3D` — subroutine
- `prep_strategy3D_thread` — subroutine
- `ref_proj_from_index` — function
- `ref_state_from_index` — function
- `seed_continuous_inplane_candidate` — subroutine

Private symbols:
- `strategy3D_alloc` — type

---
## Module: simple_strategy3D_cont

Files:
- `main/strategies/search/simple_strategy3D_cont.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_cartft_pose_opt`
- `simple_core_module_api`
- `simple_linalg`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d`
- `simple_strategy3d_srch`
- `simple_type_defs`

Public symbols:
- `check_cont_seeds` — subroutine
- `cont_reference_geometry_is_valid` — function
- `cont_seed_is_valid` — function
- `strategy3D_cont` — type

Private symbols:
- `add_work_cont` — subroutine
- `bind` — subroutine
- `kill_cont` — subroutine
- `new_cont` — subroutine
- `oris_assign_cont` — subroutine
- `srch_cont` — subroutine

---
## Module: simple_strategy3D_cont_1jyx_tester

Files:
- `main/strategies/search/simple_strategy3D_cont_1jyx_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_atoms`
- `simple_builder`
- `simple_cartft_calc`
- `simple_cartft_pose_opt`
- `simple_cmdline`
- `simple_commanders_sim`
- `simple_defs`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_fplane_grid`
- `simple_image`
- `simple_imghead`
- `simple_matcher_2dprep`
- `simple_memoize_ft_maps`
- `simple_molecule_data`
- `simple_ori`
- `simple_ori_utils`
- `simple_oris`
- `simple_parameters`
- `simple_reconstructor`
- `simple_sigma2_state_tester`
- `simple_sp_project`
- `simple_stat`
- `simple_strategy3d`
- `simple_strategy3d_cont`
- `simple_strategy3d_srch`
- `simple_string`
- `simple_sym`
- `simple_syslib`
- `simple_test_utils`
- `simple_type_defs`
- `simple_ui`

Public symbols:
- `run_all_strategy3D_cont_1jyx_tests` — subroutine

Private symbols:
- `active_upper_shell` — function
- `assert_pose_improvement` — subroutine
- `case_label` — function
- `centered_array_correlation` — function
- `cleanup_1jyx_fixture` — subroutine
- `estimate_noise_model` — subroutine
- `kill_work_images` — subroutine
- `make_perturbed_pose` — subroutine
- `memoize_mask_coords_of_box` — subroutine
- `new_work_images` — subroutine
- `prepare_particle` — subroutine
- `reconstruct_and_score` — subroutine
- `refine_all_particles` — subroutine
- `rotation_distance` — function
- `run_pose_cont_1jyx_reconstruction` — subroutine

Local symbols (enclosing procedure scope):
- `reconstruct_and_score::insert` — subroutine

---
## Module: simple_strategy3D_cont_tester

Files:
- `main/strategies/search/simple_strategy3D_cont_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_cartft_calc`
- `simple_cartft_pose_opt`
- `simple_core_module_api`
- `simple_defs_ori`
- `simple_ori`
- `simple_oris`
- `simple_parameters`
- `simple_sigma2_state_tester`
- `simple_strategy3d`
- `simple_strategy3d_cont`
- `simple_strategy3d_srch`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `run_all_strategy3D_cont_tests` — subroutine

Private symbols:
- `build_fixture` — subroutine
- `build_volume` — subroutine
- `kill_fixture` — subroutine
- `rotation_distance` — function
- `run_particle` — subroutine
- `test_commit_contract` — subroutine
- `test_lifecycle_and_rejection` — subroutine
- `test_outside_sample_untouched` — subroutine
- `test_particle_seed_contract` — subroutine
- `test_polish_keeps_convergence_fields` — subroutine
- `test_reference_geometry` — subroutine
- `test_rotation_bound` — subroutine
- `test_seed_contract` — subroutine
- `test_selected_seeds` — subroutine
- `test_threaded_batch` — subroutine
- `truth_e3` — function

Local symbols (enclosing procedure scope):
- `test_threaded_batch::record` — subroutine
- `test_threaded_batch::run_batch` — subroutine

---
## Module: simple_strategy3D_eval

Files:
- `main/strategies/search/simple_strategy3D_eval.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d_pftc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`

Public symbols:
- `strategy3D_eval` — type

Private symbols:
- `kill_eval` — subroutine
- `new_eval` — subroutine
- `oris_assign_eval` — subroutine
- `srch_eval` — subroutine

---
## Module: simple_strategy3D_greedy

Files:
- `main/strategies/search/simple_strategy3D_greedy.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d_alloc`
- `simple_strategy3d_pftc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`

Public symbols:
- `strategy3D_greedy` — type

Private symbols:
- `kill_greedy` — subroutine
- `new_greedy` — subroutine
- `oris_assign_greedy` — subroutine
- `srch_greedy` — subroutine

---
## Module: simple_strategy3D_greedy_inpl

Files:
- `main/strategies/search/simple_strategy3D_greedy_inpl.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d_alloc`
- `simple_strategy3d_pftc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`

Public symbols:
- `strategy3D_greedy_inpl` — type

Private symbols:
- `kill_greedy_inpl` — subroutine
- `new_greedy_inpl` — subroutine
- `oris_assign_greedy_inpl` — subroutine
- `srch_greedy_inpl` — subroutine

---
## Module: simple_strategy3D_greedy_sub

Files:
- `main/strategies/search/simple_strategy3D_greedy_sub.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d_alloc`
- `simple_strategy3d_pftc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`
- `simple_type_defs`

Public symbols:
- `strategy3D_greedy_sub` — type

Private symbols:
- `kill_greedy_sub` — subroutine
- `new_greedy_sub` — subroutine
- `oris_assign_greedy_sub` — subroutine
- `srch_greedy_sub` — subroutine

---
## Module: simple_strategy3D_inplane_tester

Files:
- `main/strategies/search/simple_strategy3D_inplane_tester.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs`
- `simple_forked_process`
- `simple_image`
- `simple_linked_list`
- `simple_matcher_2dprep`
- `simple_matcher_ptcl_batch`
- `simple_matcher_ptcl_io`
- `simple_matcher_smpl_and_lplims`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_pftc_shsrch_grad`
- `simple_polarft_calc`
- `simple_refine3d_strategy`
- `simple_rnd`
- `simple_strategy3d_alloc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`
- `simple_ui`
- `simple_ui_program`

Public symbols:
- `run_all_strategy3D_inplane_tests` — subroutine

Private symbols:
- `alloc_state` — subroutine
- `has_search_input` — function
- `one_particle_fixture` — subroutine
- `oversized_team_pass` — type
- `rotate_half_ring` — subroutine
- `run_oversized_team_pass` — subroutine
- `test_continuous_seed_beyond_window` — subroutine
- `test_discrete_seed_storage_and_invalid_predicate` — subroutine
- `test_grid_seed` — subroutine
- `test_inpl_cont_derivation` — subroutine
- `test_inpl_cont_policy` — subroutine
- `test_joint_storage` — subroutine
- `test_key_derivation` — subroutine
- `test_map_initialised_polish` — subroutine
- `test_pass_state_team_sized` — subroutine
- `test_pass_team_check` — subroutine
- `test_previous_shift_is_stored_shift` — subroutine
- `test_rec_observation_mask` — subroutine
- `test_resolve_inplane_e3` — subroutine
- `test_storage_with_continuous_arrays` — subroutine
- `test_storage_without_continuous_arrays` — subroutine
- `test_table_increment_axes` — subroutine

Local symbols (enclosing procedure scope):
- `test_rec_observation_mask::record_and_reread` — subroutine

---
## Module: simple_strategy3D_matcher

Files:
- `main/strategies/search/simple_strategy3D_matcher.f90`

Uses:
- `iso_fortran_env`
- `simple_bench_io`
- `simple_binoris_io`
- `simple_builder`
- `simple_cartft_pose_opt`
- `simple_eo_partition`
- `simple_eul_prob_tab`
- `simple_eul_prob_tab_neigh`
- `simple_eul_prob_tab_utils`
- `simple_matcher_2dprep`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_batch`
- `simple_matcher_ptcl_io`
- `simple_matcher_refvol_utils`
- `simple_matcher_smpl_and_lplims`
- `simple_pftc_srch_api`
- `simple_qsys_funs`
- `simple_rec3d_pcg_strategy`
- `simple_shiftmap_srch`
- `simple_state_weight_set`
- `simple_strategy3d`
- `simple_strategy3d_alloc`
- `simple_strategy3d_cont`
- `simple_strategy3d_eval`
- `simple_strategy3d_greedy`
- `simple_strategy3d_greedy_inpl`
- `simple_strategy3d_greedy_sub`
- `simple_strategy3d_prob`
- `simple_strategy3d_shc`
- `simple_strategy3d_srch`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `refine3D_exec` — subroutine

Private symbols:
- `print_flags` — subroutine
- `refine3D_ctrl` — type
- `write_partial_recs` — subroutine

Local symbols (enclosing procedure scope):
- `refine3D_exec::begin_pipeline` — subroutine
- `refine3D_exec::build_batch_particles_local` — subroutine
- `refine3D_exec::choose_and_run_strategy` — subroutine
- `refine3D_exec::ensure_even_odd_partition` — subroutine
- `refine3D_exec::init_ctrl` — subroutine
- `refine3D_exec::maybe_write_orientations` — subroutine
- `refine3D_exec::on_disk_classic` — function
- `refine3D_exec::pipeline_batch` — subroutine
- `refine3D_exec::pipeline_context_bytes` — function
- `refine3D_exec::prepare_particles_batches` — subroutine
- `refine3D_exec::prepare_refs_sigmas_and_pftc` — subroutine
- `refine3D_exec::sample_particles_for_update` — subroutine
- `refine3D_exec::strategy3D_per_ptcl` — type
- `refine3D_exec::write_bench_file` — subroutine
- `refine3D_exec::write_empty_readout` — subroutine

---
## Module: simple_strategy3D_pftc

Files:
- `main/strategies/search/simple_strategy3D_pftc.f90`

Uses:
- `simple_strategy3d`
- `simple_strategy3d_srch`

Public symbols:
- `strategy3D_pftc` — type

---
## Module: simple_strategy3D_prob

Files:
- `main/strategies/search/probabilistic/simple_strategy3D_prob.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d_pftc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`

Public symbols:
- `strategy3D_prob` — type

Private symbols:
- `has_valid_assignment` — function
- `kill_prob` — subroutine
- `new_prob` — subroutine
- `oris_assign_prob` — subroutine
- `resolve_extra_poses` — subroutine
- `srch_prob` — subroutine

---
## Module: simple_strategy3D_shc

Files:
- `main/strategies/search/simple_strategy3D_shc.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d_alloc`
- `simple_strategy3d_pftc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`

Public symbols:
- `strategy3D_shc` — type

Private symbols:
- `kill_shc` — subroutine
- `new_shc` — subroutine
- `oris_assign_shc` — subroutine
- `srch_shc` — subroutine

---
## Module: simple_strategy3D_srch

Files:
- `main/strategies/search/simple_strategy3D_srch.f90`

Uses:
- `ieee_arithmetic`
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab`
- `simple_eul_prob_tab_utils`
- `simple_parameters`
- `simple_pftc_shsrch_grad`
- `simple_shiftmap_srch`
- `simple_strategy3d_alloc`

Public symbols:
- `strategy3D_spec` — type
- `strategy3D_srch` — type

Private symbols:
- `get_continuous_route_status` — subroutine
- `inpl_srch` — subroutine
- `inpl_srch_first` — subroutine
- `inpl_srch_peaks` — subroutine
- `joint_evaluation_invalid` — function
- `kill` — subroutine
- `map_srch` — subroutine
- `new` — subroutine
- `prep4prob` — subroutine
- `prep4srch` — subroutine
- `record_shift` — subroutine
- `refine_assignment_continuously` — subroutine
- `refine_assignment_continuously_impl` — subroutine
- `refine_extra_pose_continuously` — subroutine
- `refine_selected_continuously` — subroutine
- `refine_selected_continuously_impl` — subroutine
- `store_continuous_solution` — subroutine
- `store_discrete_seed_solution` — subroutine
- `store_solution` — subroutine
- `uses_continuous_refinement` — function
- `write_continuous_route_status` — subroutine

---
## Module: simple_strategy3D_utils

Files:
- `main/strategies/search/simple_strategy3D_utils.f90`

Uses:
- `simple_core_module_api`
- `simple_pftc_shsrch_grad`
- `simple_strategy3d_alloc`
- `simple_strategy3d_srch`

Public symbols:
- `assign_ori` — subroutine
- `candidate_pose` — subroutine
- `extract_peak_ori` — subroutine
- `extract_peak_oris` — subroutine
- `resolve_inplane_e3` — function

---
## Module: simple_strategy_sigma2

Files:
- `main/strategies/parallelization/simple_strategy_sigma2.f90`

Uses:
- `iso_fortran_env`
- `simple_bench_io`
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs`
- `simple_error`
- `simple_parameters`
- `simple_refine3d_provenance`
- `simple_sigma2_state`
- `simple_sigma2_state_file`
- `simple_string`

Public symbols:
- `commit_sigma2_exec` — subroutine
- `prepare_canonical_sigma_update` — subroutine
- `read_sigma2_bench_counts` — subroutine
- `run_commit_sigma2` — subroutine

Local symbols (enclosing procedure scope):
- `commit_sigma2_exec::consolidate_canonical` — subroutine
- `commit_sigma2_exec::write_commit_bench` — subroutine

---
## Module: simple_stream_api

Files:
- `main/apis/simple_stream_api.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_class_frcs`
- `simple_cmdline`
- `simple_commander_base`
- `simple_core_module_api`
- `simple_defs_environment`
- `simple_gui_utils`
- `simple_image`
- `simple_parameters`
- `simple_progress`
- `simple_projfile_utils`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_rec_list`
- `simple_sp_project`
- `simple_stack_io`
- `simple_starproject_stream`
- `simple_stream_refine2d_utils`
- `simple_stream_utils`
- `simple_stream_watcher`

---
## Module: simple_stream_chain_tester

Files:
- `main/stream/simple_stream_chain_tester.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_reproject`
- `simple_commanders_sim`
- `simple_commanders_stream_p01_preprocess`
- `simple_commanders_stream_p03_initial_analysis`
- `simple_commanders_stream_p04_refpick_extract`
- `simple_commanders_stream_p05_sieve_cavgs`
- `simple_commanders_stream_p06_pool2d`
- `simple_commanders_stream_p07_solve3d_multistate`
- `simple_defs`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_image`
- `simple_ori`
- `simple_oris`
- `simple_sp_project`
- `simple_stream_stage_initial_analysis`
- `simple_stream_stage_optics`
- `simple_stream_stage_pool2d`
- `simple_stream_stage_preprocess`
- `simple_stream_stage_refpick`
- `simple_stream_stage_sieve`
- `simple_stream_stage_solve3d`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_truth_metrics`
- `simple_test_utils`
- `simple_timer`
- `simple_type_defs`
- `simple_ui`

Public symbols:
- `run_all_stream_chain_movie_tests` — subroutine
- `run_all_stream_chain_sieve_tests` — subroutine

Private symbols:
- `assert_multistate_truth` — subroutine
- `base_cline` — subroutine
- `enter` — subroutine
- `initial_analysis_cline` — subroutine
- `interleave_stacks` — function
- `no_folder` — subroutine
- `optics_cline` — subroutine
- `pool_cline` — subroutine
- `preprocess_cline` — subroutine
- `refpick_cline` — subroutine
- `reproject_truth` — function
- `sieve_cline` — subroutine
- `simulate_movies` — subroutine
- `simulate_particle_stack` — function
- `solve3D_cline` — subroutine
- `stage_dir` — function
- `test_movies_to_solve3D` — subroutine
- `test_sieve_to_solve3D` — subroutine
- `write_moldiam` — subroutine
- `write_picking_sets` — subroutine
- `write_repeated_orientations` — subroutine
- `write_truth_volume` — subroutine

---
## Module: simple_stream_communicator

Files:
- `utils/comm/simple_stream_communicator.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_core_module_api`
- `simple_fileio`
- `unix`

Public symbols:
- `stream_http_communicator` — type

Private symbols:
- `active_status` — function
- `add_to_json_char` — subroutine
- `add_to_json_dp` — subroutine
- `add_to_json_int` — subroutine
- `add_to_json_jval` — subroutine
- `add_to_json_logical` — subroutine
- `arg_associated` — function
- `background_heartbeat` — subroutine
- `create` — subroutine
- `curl_request` — subroutine
- `destroy_arg` — subroutine
- `evaluate_checksum` — function
- `exit_status` — function
- `get_json_arg_char` — subroutine
- `get_json_arg_dp` — subroutine
- `get_json_arg_int` — subroutine
- `get_json_arg_int_arr` — subroutine
- `get_json_arg_logical` — subroutine
- `join_background_heartbeat` — subroutine
- `remove_from_json_if_present` — subroutine
- `send_heartbeat` — subroutine
- `send_jobstats` — subroutine
- `stop_status` — function
- `term` — subroutine
- `update_json_char` — subroutine
- `update_json_dp` — subroutine
- `update_json_int` — subroutine
- `update_json_logical` — subroutine

---
## Module: simple_stream_gui_senders

Files:
- `main/stream/shared/simple_stream_gui_senders.f90`

Uses:
- `simple_defs`
- `simple_fileio`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_micrograph`
- `simple_nrtxtfile`
- `simple_oris`
- `simple_stream_pipe`
- `simple_string`

Public symbols:
- `send_cavgs` — subroutine
- `send_recent_micrographs` — subroutine
- `send_reproj_tiles` — subroutine

Local symbols (enclosing procedure scope):
- `send_cavgs::tile_pos` — function

---
## Module: simple_stream_job_sets

Files:
- `main/stream/shared/simple_stream_job_sets.f90`

Uses:
- `simple_cmdline`
- `simple_defs_fname`
- `simple_fileio`
- `simple_qsys_env`
- `simple_qsys_job_record`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `stream_job_sets` — type

Private symbols:
- `cancel` — subroutine
- `collect` — subroutine
- `complete` — subroutine
- `get_completed_dir` — function
- `get_counter` — function
- `get_job_dir` — function
- `kill` — subroutine
- `make_job_dir` — subroutine
- `new` — subroutine
- `restore` — subroutine
- `schedule` — subroutine
- `set_number` — function
- `submit` — subroutine
- `write_set` — subroutine

---
## Module: simple_stream_job_sets_tester

Files:
- `main/stream/shared/simple_stream_job_sets_tester.f90`

Uses:
- `simple_cmdline`
- `simple_defs_fname`
- `simple_fileio`
- `simple_sp_project`
- `simple_stream_job_sets`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_stream_job_sets_tests` — subroutine

Private symbols:
- `test_complete` — subroutine
- `test_new_makes_folders` — subroutine
- `test_restore` — subroutine
- `test_restore_past_unfinished` — subroutine
- `test_write_set` — subroutine

---
## Module: simple_stream_master_gui_commands

Files:
- `main/stream/master/simple_stream_master_gui_commands.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_defs`
- `simple_defs_fname`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_types`
- `simple_stream_master_stage_ids`
- `simple_string`

Public symbols:
- `stream_master_gui_commands` — type

Private symbols:
- `kill` — subroutine
- `parse` — function

Local symbols (enclosing procedure scope):
- `parse::flag` — function
- `parse::snapshot_fname_ok` — function
- `parse::warn_dropped` — subroutine

---
## Module: simple_stream_master_meta_store

Files:
- `main/stream/master/simple_stream_master_meta_store.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_error`
- `simple_gui_assembler`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_histogram`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_optics_group`
- `simple_gui_metadata_stream_initial_analysis`
- `simple_gui_metadata_stream_optics_assignment`
- `simple_gui_metadata_stream_particle_sieving`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_stream_pool2d`
- `simple_gui_metadata_stream_preprocess`
- `simple_gui_metadata_stream_snapshot`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_gui_metadata_timeplot`
- `simple_gui_metadata_types`
- `simple_gui_metadata_vol3d`
- `simple_stream_master_stage_ids`
- `simple_string_utils`

Public symbols:
- `stream_master_meta_store` — type

Private symbols:
- `assemble` — subroutine
- `clear_stage` — subroutine
- `fits` — function
- `frame_fits` — function
- `kill` — subroutine
- `log_state` — subroutine
- `new` — subroutine
- `place_cavg2D` — subroutine
- `place_micrograph` — subroutine
- `place_optics_group` — subroutine
- `place_vol3D` — subroutine
- `store` — subroutine

Local symbols (enclosing procedure scope):
- `log_state::ncavg` — function
- `log_state::nmic` — function
- `log_state::nog` — function
- `log_state::nvol` — function

---
## Module: simple_stream_master_resources

Files:
- `main/stream/master/simple_stream_master_resources.f90`

Uses:
- `simple_defs`
- `simple_defs_environment`
- `simple_string_utils`

Public symbols:
- `stream_resources` — type
- `stream_resources_from_env` — function

Private symbols:
- `log_resources` — subroutine

Local symbols (enclosing procedure scope):
- `stream_resources_from_env::env_override` — subroutine

---
## Module: simple_stream_master_stage

Files:
- `main/stream/master/simple_stream_master_stage.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_defs`
- `simple_error`
- `simple_forked_process`
- `simple_gui_assembler`
- `simple_persistent_worker_server`
- `simple_qsys_env`
- `simple_qsys_job_record`
- `simple_stream_master_stage_ids`
- `simple_stream_pipe`
- `simple_string`

Public symbols:
- `fork_gui_status` — function
- `stream_master_stage` — type
- `stream_master_stage_fork` — type

Private symbols:
- `discard_pipes` — subroutine
- `execute_stage` — subroutine
- `force_stop` — subroutine
- `get_label` — function
- `is_new_update` — function
- `is_running` — function
- `is_skipped` — function
- `kill` — subroutine
- `new` — subroutine
- `request_stop` — subroutine
- `send_update` — subroutine
- `skip` — subroutine
- `start` — subroutine

---
## Module: simple_stream_master_stage_ids

Files:
- `main/stream/master/simple_stream_master_stage_ids.f90`

Uses:
- `simple_defs_stream`
- `simple_error`
- `simple_stream_state`
- `unix`

Public symbols:
- `close_other_pipe_ends` — subroutine
- `close_stage_pipes` — subroutine
- `master_fds` — subroutine
- `open_stage_pipes` — subroutine
- `stage_fds` — subroutine
- `stage_gui_key` — function
- `stage_job_name` — function
- `stage_label` — function

Private symbols:
- `apply` — subroutine
- `pipe_op` — subroutine
- `stage_pipes` — subroutine

---
## Module: simple_stream_master_tester

Files:
- `main/stream/master/simple_stream_master_tester.f90`

Uses:
- `simple_cmdline`
- `simple_commander_base`
- `simple_forked_process`
- `simple_gui_assembler`
- `simple_gui_metadata_api`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_utils`
- `simple_stream_master_gui_commands`
- `simple_stream_master_meta_store`
- `simple_stream_master_stage`
- `simple_stream_master_stage_ids`
- `simple_stream_pipe`
- `simple_string`
- `simple_string_utils`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_master_tests` — subroutine
- `run_stream_heartbeat_tests` — subroutine

Private symbols:
- `cavg_idx` — function
- `cavg_message` — function
- `execute_noop` — subroutine
- `int_list` — function
- `noop_commander` — type
- `test_fork_gui_status` — subroutine
- `test_gui_commands` — subroutine
- `test_gui_commands_fresh_each_answer` — subroutine
- `test_gui_commands_invalid` — subroutine
- `test_gui_commands_oversized` — subroutine
- `test_gui_commands_snapshot3D` — subroutine
- `test_gui_commands_snapshot_name` — subroutine
- `test_skipped_stays_skipped` — subroutine
- `test_stage_names` — subroutine
- `test_stage_pipes` — subroutine
- `test_store_clear_stage` — subroutine
- `test_store_list` — subroutine
- `test_store_snapshot3D` — subroutine
- `test_store_status` — subroutine
- `test_store_volume` — subroutine
- `test_update_dedupe` — subroutine

Local symbols (enclosing procedure scope):
- `test_gui_commands_snapshot3D::snapshot3D_kept` — function
- `test_gui_commands_snapshot_name::snapshot_kept` — function

---
## Module: simple_stream_meta_plots

Files:
- `main/stream/shared/simple_stream_meta_plots.f90`

Uses:
- `simple_gui_metadata_histogram`
- `simple_gui_metadata_timeplot`
- `simple_histogram`
- `simple_oris`
- `simple_string`

Public symbols:
- `recent_shifts_by_optics_group` — subroutine
- `set_histogram_from_oris` — subroutine
- `set_rate_timeplot` — subroutine
- `set_timeplot_from_oris` — subroutine

---
## Module: simple_stream_meta_plots_tester

Files:
- `main/stream/shared/simple_stream_meta_plots_tester.f90`

Uses:
- `simple_gui_metadata_timeplot`
- `simple_gui_metadata_types`
- `simple_oris`
- `simple_stream_meta_plots`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_stream_meta_plots_tests` — subroutine

Private symbols:
- `get_plot` — subroutine
- `make_oris` — subroutine
- `test_rate_timeplot` — subroutine
- `test_rate_timeplot_newest` — subroutine
- `test_timeplot_widens` — subroutine
- `test_timeplot_windows` — subroutine

---
## Module: simple_stream_pipe

Files:
- `main/stream/shared/simple_stream_pipe.f90`

Uses:
- `iso_c_binding`
- `simple_error`
- `simple_gui_metadata_base`
- `unix`

Public symbols:
- `stream_pipe` — type

Private symbols:
- `consume` — subroutine
- `discard` — subroutine
- `extract_frame` — function
- `get_pending_bytes` — function
- `is_broken` — function
- `kill` — subroutine
- `limit_partial_frames` — subroutine
- `new` — subroutine
- `receive` — function
- `send` — subroutine
- `send_meta` — subroutine

---
## Module: simple_stream_pipe_tester

Files:
- `main/stream/shared/simple_stream_pipe_tester.f90`

Uses:
- `iso_c_binding`
- `simple_stream_pipe`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_pipe_tests` — subroutine

Private symbols:
- `close_loopback` — subroutine
- `open_loopback` — subroutine
- `test_bad_length_resyncs` — subroutine
- `test_discard` — subroutine
- `test_empty_and_closed_ends` — subroutine
- `test_frame_split_across_writes` — subroutine
- `test_partial_frame_abandoned` — subroutine
- `test_queued_frames_all_delivered` — subroutine
- `test_single_frame_round_trip` — subroutine

---
## Module: simple_stream_pool2D

Files:
- `main/stream/pool2D/simple_stream_pool2D.f90`

Uses:
- `simple_class_frcs`
- `simple_classaverager`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_environment`
- `simple_image`
- `simple_imgarr_utils`
- `simple_optics_maps`
- `simple_parameters`
- `simple_procimgstk`
- `simple_qsys_async_job`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_qsys_job_record`
- `simple_ran_tabu`
- `simple_refine2d_fnames`
- `simple_sp_project`
- `simple_stack_io`
- `simple_starproject`
- `simple_starproject_stream`
- `simple_stream_refine2d_utils`
- `simple_stream_utils`

Public symbols:
- `stream_pool2D` — type
- `stream_pool2D_stats` — type

Private symbols:
- `append_sets` — subroutine
- `apply_snapshot_selection` — subroutine
- `available` — function
- `cancel` — subroutine
- `failed` — function
- `finalise` — subroutine
- `history_slot` — function
- `init_queue` — subroutine
- `init_state` — subroutine
- `iterate` — subroutine
- `iteration` — function
- `job_failed` — function
- `kill` — subroutine
- `new` — subroutine
- `npublishable` — function
- `publish` — subroutine
- `rank_cavgs` — subroutine
- `rescale_cavgs` — subroutine
- `set_dimensions` — subroutine
- `set_mask` — subroutine
- `set_mskdiam` — subroutine
- `set_resolution_limits` — subroutine
- `start` — subroutine
- `stats` — function
- `submit_iteration` — subroutine
- `update` — subroutine
- `update_aln_params` — subroutine
- `update_dims` — subroutine
- `update_for_gui` — subroutine
- `update_status` — subroutine
- `write_jpeg` — subroutine
- `write_project` — subroutine
- `write_snapshot` — subroutine
- `write_stats` — subroutine

Local symbols (enclosing procedure scope):
- `finalise::write_raw_project` — subroutine
- `iterate::uniform_stack_sampling` — subroutine
- `write_project::rescale_refs` — subroutine

---
## Module: simple_stream_pool2D_tester

Files:
- `main/stream/pool2D/simple_stream_pool2D_tester.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_fileio`
- `simple_parameters`
- `simple_sp_project`
- `simple_stream_pool2d`
- `simple_stream_refine2d_utils`
- `simple_stream_utils`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_stream_pool2D_tests` — subroutine

Private symbols:
- `make_set` — subroutine
- `make_started_pool` — subroutine
- `set_test_cline` — subroutine
- `test_append_project_sets` — subroutine
- `test_append_sets` — subroutine
- `test_draw_new_classes` — subroutine
- `test_mask_clamped_to_box` — subroutine
- `test_new_kill` — subroutine
- `test_restart_after_kill` — subroutine
- `test_snapshot_without_frcs` — subroutine
- `test_unclassified_publishes_nothing` — subroutine

Local symbols (enclosing procedure scope):
- `test_draw_new_classes::make_draw_project` — subroutine

---
## Module: simple_stream_refine2D_utils

Files:
- `main/stream/pool2D/simple_stream_refine2D_utils.f90`

Uses:
- `iso_fortran_env`
- `simple_class_frcs`
- `simple_core_module_api`
- `simple_image`
- `simple_optics_maps`
- `simple_refine2d_fnames`
- `simple_sp_project`

Public symbols:
- `append_project_sets` — subroutine
- `build_pool_publication` — subroutine
- `build_sieve_publication` — subroutine
- `cleanup_root_folder` — subroutine
- `combine_sieve_classes` — function
- `delete_pool_publication` — subroutine
- `draw_new_classes` — subroutine
- `pool_publication_names` — subroutine
- `pool_publication_nselected` — function
- `snapshot_cavgs_meta` — subroutine
- `tidy_2Dstream_iter` — subroutine

---
## Module: simple_stream_sigterm

Files:
- `main/stream/shared/simple_stream_sigterm.f90`

Uses:
- `iso_c_binding`
- `unix`

Public symbols:
- `install_sigterm_handler` — subroutine
- `restore_sigterm_handler` — subroutine
- `sigterm_received` — function

Private symbols:
- `on_sigterm` — subroutine

---
## Module: simple_stream_stage_initial_analysis

Files:
- `main/stream/stages/simple_stream_stage_initial_analysis.f90`

Uses:
- `simple_cavg_quality_model`
- `simple_cavg_quality_selection`
- `simple_cavg_quality_types`
- `simple_class_compatibility`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_environment`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_stream_initial_analysis`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_gui_metadata_vol3d`
- `simple_gui_utils`
- `simple_image`
- `simple_image_bin`
- `simple_image_msk`
- `simple_imgarr_utils`
- `simple_imghead`
- `simple_mic_import`
- `simple_mic_selection`
- `simple_parameters`
- `simple_procimgstk`
- `simple_ptcl_sieve`
- `simple_qsys_async_job`
- `simple_qsys_env`
- `simple_rec_list`
- `simple_refine3d_fnames`
- `simple_segdiam_bin_picker`
- `simple_sp_project`
- `simple_stream_gui_senders`
- `simple_stream_pipe`
- `simple_stream_sigterm`
- `simple_stream_state`
- `simple_stream_utils`
- `simple_stream_watcher`

Public symbols:
- `stream_stage_initial_analysis` — type

Private symbols:
- `all_projfile` — function
- `apply_gui_updates` — subroutine
- `attach_upstream` — subroutine
- `balance_classes` — subroutine
- `balanced_projfile` — function
- `choose_state` — function
- `clear_previous_run` — subroutine
- `collect_extractions` — subroutine
- `cycle_projfile` — function
- `duplicate_balanced_stack` — subroutine
- `estimate_mskdiam` — function
- `finalize` — subroutine
- `find_final_solve3D_cavgs_dir` — subroutine
- `finish_extract` — subroutine
- `finish_solve3D` — subroutine
- `finished` — function
- `import_projects` — subroutine
- `init_gui` — subroutine
- `init_params` — subroutine
- `init_queue` — subroutine
- `iterate` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `pick_extract_all` — subroutine
- `publish_pickrefs` — subroutine
- `publish_reprojections` — subroutine
- `rebuild_init_mics` — subroutine
- `restore_cavgs_selection` — subroutine
- `restore_pickrefs` — subroutine
- `run_cycle1` — subroutine
- `run_cycle2` — subroutine
- `save_pickrefs_selection` — subroutine
- `select_and_send` — subroutine
- `select_project_cavgs` — subroutine
- `send_initial_analysis_status` — subroutine
- `send_picking_status` — subroutine
- `send_pickrefs` — subroutine
- `start_extract` — subroutine
- `start_reproject` — subroutine
- `start_sieve` — subroutine
- `start_solve2D` — subroutine
- `start_solve3D` — subroutine
- `state_vetoes` — function
- `unique_projnames` — function
- `update_os_out_stk` — subroutine
- `update_upstream_state` — subroutine

Local symbols (enclosing procedure scope):
- `finish_solve3D::state_fsc_res` — function
- `finish_solve3D::vol_mskdiam` — function
- `kill::release_heavy` — subroutine

---
## Module: simple_stream_stage_initial_analysis_tester

Files:
- `main/stream/stages/simple_stream_stage_initial_analysis_tester.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_commanders_stream_p03_initial_analysis`
- `simple_default_clines`
- `simple_defs`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_stream_initial_analysis`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_image`
- `simple_image_msk`
- `simple_imghead`
- `simple_pick_strategy`
- `simple_qsys_async_job`
- `simple_segdiam_bin_picker`
- `simple_sp_project`
- `simple_stream_pipe`
- `simple_stream_stage_initial_analysis`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_stage_initial_analysis_tests` — subroutine

Private symbols:
- `check_stack_size` — subroutine
- `close_loopback` — subroutine
- `count_class` — function
- `make_disc_cavg` — subroutine
- `make_test_stage` — subroutine
- `make_upstream` — subroutine
- `meta_type_of` — function
- `open_loopback` — subroutine
- `set_test_cline` — subroutine
- `swap_ext` — function
- `test_attach_upstream_waits` — subroutine
- `test_balance_classes` — subroutine
- `test_choose_state` — subroutine
- `test_clear_previous_run` — subroutine
- `test_commander_defaults` — subroutine
- `test_cycle1_setup_waits` — subroutine
- `test_estimate_mskdiam` — subroutine
- `test_find_final_solve3D_dir` — subroutine
- `test_finished` — subroutine
- `test_gui_selection_ends_stage` — subroutine
- `test_import_projects` — subroutine
- `test_init_params` — subroutine
- `test_iterate_passes` — subroutine
- `test_paths` — subroutine
- `test_published_pickrefs_are_final` — subroutine
- `test_rebuild_init_mics` — subroutine
- `test_restore_cavgs_selection` — subroutine
- `test_segdiam_bins_count_erosion` — subroutine
- `test_status_messages` — subroutine
- `write_cavgs_project` — subroutine
- `write_class_stack` — subroutine
- `write_completed_project` — subroutine

---
## Module: simple_stream_stage_optics

Files:
- `main/stream/stages/simple_stream_stage_optics.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_gui_metadata_optics_group`
- `simple_gui_metadata_stream_optics_assignment`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_mic_import`
- `simple_optics_groups`
- `simple_optics_maps`
- `simple_parameters`
- `simple_sp_project`
- `simple_starproject_stream`
- `simple_stream_meta_plots`
- `simple_stream_pipe`
- `simple_stream_state`
- `simple_stream_watcher`
- `simple_string`
- `simple_syslib`
- `simple_timer`

Public symbols:
- `stream_stage_optics` — type

Private symbols:
- `assign_and_publish` — subroutine
- `attach_upstream` — subroutine
- `end_restore` — subroutine
- `finalize` — subroutine
- `finished` — function
- `import_new_projects` — subroutine
- `iterate` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `restore_groups` — subroutine
- `send_group_shifts` — subroutine
- `send_status` — subroutine

---
## Module: simple_stream_stage_optics_tester

Files:
- `main/stream/stages/simple_stream_stage_optics_tester.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_gui_metadata_optics_group`
- `simple_gui_metadata_stream_optics_assignment`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_sp_project`
- `simple_stream_pipe`
- `simple_stream_stage_optics`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_stage_optics_tests` — subroutine

Private symbols:
- `check_two_groups` — subroutine
- `close_loopback` — subroutine
- `make_test_stage` — subroutine
- `make_upstream` — subroutine
- `map_file` — function
- `open_loopback` — subroutine
- `set_test_cline` — subroutine
- `test_assign_and_publish` — subroutine
- `test_attach_upstream_waits` — subroutine
- `test_beamtilt_from_command_line` — subroutine
- `test_finished` — subroutine
- `test_import_new_projects` — subroutine
- `test_iterate_passes` — subroutine
- `test_map_ids_continue_after_restart` — subroutine
- `test_new_starts_empty` — subroutine
- `test_restart_removes_termination_file` — subroutine
- `test_restart_restores_groups` — subroutine
- `test_send_group_shifts` — subroutine
- `test_send_status` — subroutine
- `write_completed_project` — subroutine
- `write_mic_project` — subroutine

---
## Module: simple_stream_stage_pool2D

Files:
- `main/stream/stages/simple_stream_stage_pool2D.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_error`
- `simple_fileio`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_stream_pool2d`
- `simple_gui_metadata_stream_snapshot`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_optics_maps`
- `simple_parameters`
- `simple_rec_list`
- `simple_sp_project`
- `simple_srch_sort_loc`
- `simple_stream_gui_senders`
- `simple_stream_pipe`
- `simple_stream_pool2d`
- `simple_stream_refine2d_utils`
- `simple_stream_state`
- `simple_stream_utils`
- `simple_stream_watcher`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `stream_stage_pool2D` — type

Private symbols:
- `apply_final_mskdiam` — subroutine
- `apply_gui_updates` — subroutine
- `apply_pause_policy` — subroutine
- `attach_upstream` — subroutine
- `clean_previous_run` — subroutine
- `collect_sieve_set` — subroutine
- `collects_sieve_sets` — function
- `default_mskdiam` — function
- `export_pool_state` — subroutine
- `exports_after` — function
- `finalize` — subroutine
- `finished` — function
- `first_set_due` — function
- `import_sets` — subroutine
- `init_gui` — subroutine
- `init_params` — subroutine
- `is_final_set` — function
- `iterate` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `owes_final_run` — function
- `pause_rate_factor` — function
- `publication_fname` — function
- `publish_first_from_sieve` — subroutine
- `publish_first_set` — subroutine
- `publish_sieve_set` — subroutine
- `publishes_final` — function
- `read_final_mskdiam` — subroutine
- `restore_export_id` — subroutine
- `run_iteration` — subroutine
- `runs_to_final` — function
- `send_pool_cavgs` — subroutine
- `send_snapshot` — subroutine
- `send_status` — subroutine
- `set_mskdiam` — subroutine
- `sort_sets` — subroutine
- `start_pool` — subroutine
- `target_nptcls` — function
- `transfer_sets` — subroutine
- `unpause` — subroutine
- `update_pool_progress` — subroutine
- `watch_sets` — subroutine
- `write_snapshot` — subroutine

Local symbols (enclosing procedure scope):
- `sort_sets::handed_off_before` — function
- `sort_sets::trailing_number` — function

---
## Module: simple_stream_stage_pool2D_tester

Files:
- `main/stream/stages/simple_stream_stage_pool2D_tester.f90`

Uses:
- `iso_c_binding`
- `simple_class_frcs`
- `simple_cmdline`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_gui_metadata_stream_pool2d`
- `simple_gui_metadata_stream_snapshot`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_image`
- `simple_imghead`
- `simple_optics_maps`
- `simple_rec_list`
- `simple_sp_project`
- `simple_stream_pipe`
- `simple_stream_refine2d_utils`
- `simple_stream_stage_pool2d`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_stage_pool2D_tests` — subroutine

Private symbols:
- `close_loopback` — subroutine
- `make_optics_map` — subroutine
- `make_test_stage` — subroutine
- `make_upstream` — subroutine
- `open_loopback` — subroutine
- `set_test_cline` — subroutine
- `test_attach_and_watch` — subroutine
- `test_export_numbering` — subroutine
- `test_finished` — subroutine
- `test_given_mskdiam_kept` — subroutine
- `test_gui_mskdiam_update` — subroutine
- `test_init_params` — subroutine
- `test_iterate_waits` — subroutine
- `test_pause_rules` — subroutine
- `test_publication_holds_classified_stacks` — subroutine
- `test_restart_cleans` — subroutine
- `test_send_snapshot` — subroutine
- `test_send_status` — subroutine
- `test_sets_taken_in_order` — subroutine
- `test_sieve_final_set` — subroutine
- `test_sieve_publication` — subroutine
- `test_transfer_sets` — subroutine
- `test_transfer_stepwise` — subroutine
- `write_empty_final_set` — subroutine
- `write_sieved_set` — function

Local symbols (enclosing procedure scope):
- `test_sieve_publication::make_sieve_set` — subroutine
- `test_sieve_publication::write_class_files` — subroutine

---
## Module: simple_stream_stage_preprocess

Files:
- `main/stream/stages/simple_stream_stage_preprocess.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_defs_environment`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_error`
- `simple_fileio`
- `simple_gui_metadata_histogram`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_stream_preprocess`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_timeplot`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_image`
- `simple_mic_import`
- `simple_mic_selection`
- `simple_motion_correct_utils`
- `simple_motion_gain_analysis`
- `simple_motion_gain_helpers`
- `simple_oris`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_sp_project`
- `simple_starproject_stream`
- `simple_stream_job_sets`
- `simple_stream_meta_plots`
- `simple_stream_pipe`
- `simple_stream_sigterm`
- `simple_stream_state`
- `simple_stream_watcher`
- `simple_string`
- `simple_syslib`
- `simple_timer`
- `simple_type_defs`

Public symbols:
- `stream_stage_preprocess` — type

Private symbols:
- `apply_gui_updates` — subroutine
- `apply_thresholds` — subroutine
- `build_worker_cline` — subroutine
- `collect_jobs` — subroutine
- `create_movies_set_project` — subroutine
- `detect_gain_flip` — function
- `finalize` — subroutine
- `finished` — function
- `generate_gain_from_movies` — function
- `idle` — subroutine
- `import_completed` — subroutine
- `import_previous_projects` — subroutine
- `init_gui` — subroutine
- `init_job_dirs` — subroutine
- `init_movie_watcher` — subroutine
- `init_params` — subroutine
- `init_queue` — subroutine
- `iterate` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `new_movie_watcher` — subroutine
- `next_movie_batch` — subroutine
- `process_imports` — subroutine
- `resolve_gain` — subroutine
- `restore_thresholds` — subroutine
- `resume_previous_run` — subroutine
- `save_thresholds` — subroutine
- `schedule_jobs` — subroutine
- `send_plots` — subroutine
- `send_status` — subroutine
- `set_threshold` — subroutine
- `submit_new_movies` — subroutine
- `update_idle_marker` — subroutine
- `write_mic_star_and_field` — subroutine

Local symbols (enclosing procedure scope):
- `kill::release_heavy` — subroutine

---
## Module: simple_stream_stage_preprocess_tester

Files:
- `main/stream/stages/simple_stream_stage_preprocess_tester.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_image`
- `simple_oris`
- `simple_sp_project`
- `simple_stream_pipe`
- `simple_stream_stage_preprocess`
- `simple_stream_watcher`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_stage_preprocess_tests` — subroutine

Private symbols:
- `close_loopback` — subroutine
- `make_test_stage` — subroutine
- `open_loopback` — subroutine
- `set_mic` — subroutine
- `set_test_cline` — subroutine
- `test_apply_gui_updates` — subroutine
- `test_build_worker_cline` — subroutine
- `test_create_movies_set_project` — subroutine
- `test_finished` — subroutine
- `test_import_completed` — subroutine
- `test_import_previous_projects` — subroutine
- `test_init_params_split_mode` — subroutine
- `test_process_imports` — subroutine
- `test_resolve_gain_static_flip` — subroutine
- `test_send_status` — subroutine
- `write_gainref` — subroutine
- `write_mic_project` — subroutine
- `write_movie` — subroutine

---
## Module: simple_stream_stage_refpick

Files:
- `main/stream/stages/simple_stream_stage_refpick.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_defs_environment`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_error`
- `simple_fileio`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_micrograph`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_gui_utils`
- `simple_mic_import`
- `simple_optics_maps`
- `simple_oris`
- `simple_parameters`
- `simple_ptcl_sieve`
- `simple_qsys_async_job`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_sp_project`
- `simple_starproject_stream`
- `simple_stream_gui_senders`
- `simple_stream_job_sets`
- `simple_stream_pipe`
- `simple_stream_state`
- `simple_stream_utils`
- `simple_stream_watcher`
- `simple_string`
- `simple_syslib`
- `simple_timer`

Public symbols:
- `stream_stage_refpick` — type

Private symbols:
- `apply_optics_map` — function
- `attach_upstream` — subroutine
- `build_worker_cline` — subroutine
- `collect_jobs` — subroutine
- `create_set_project` — subroutine
- `finalize` — subroutine
- `finished` — function
- `first_mics_smpd` — function
- `idle` — subroutine
- `import_finished_sets` — subroutine
- `import_previous_sets` — subroutine
- `import_sets` — subroutine
- `init_gui` — subroutine
- `init_job_dirs` — subroutine
- `init_params` — subroutine
- `init_queue` — subroutine
- `install_pickrefs` — subroutine
- `iterate` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `pickrefs_available` — function
- `prepare_pickrefs` — subroutine
- `process_imports` — subroutine
- `reset_dir` — subroutine
- `resume_previous_run` — subroutine
- `schedule_jobs` — subroutine
- `send_status` — subroutine
- `submit_new_projects` — subroutine
- `update_idle_marker` — subroutine
- `write_mic_star` — subroutine
- `write_project` — subroutine

Local symbols (enclosing procedure scope):
- `kill::release_heavy` — subroutine

---
## Module: simple_stream_stage_refpick_tester

Files:
- `main/stream/stages/simple_stream_stage_refpick_tester.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_optics_maps`
- `simple_sp_project`
- `simple_stream_pipe`
- `simple_stream_stage_refpick`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_stage_refpick_tests` — subroutine

Private symbols:
- `close_loopback` — subroutine
- `make_extracted_set` — subroutine
- `make_test_stage` — subroutine
- `make_upstream` — subroutine
- `open_loopback` — subroutine
- `publish_test_optics_map` — subroutine
- `set_test_cline` — subroutine
- `test_attach_upstream` — subroutine
- `test_create_set_project` — subroutine
- `test_finished` — subroutine
- `test_idle_marker_follows_upstream` — subroutine
- `test_import_finished_sets` — subroutine
- `test_init_params` — subroutine
- `test_iterate_waits` — subroutine
- `test_optics_map_applied` — subroutine
- `test_pickrefs_available` — subroutine
- `test_restart_clear` — subroutine
- `test_restart_history` — subroutine
- `test_send_status` — subroutine
- `test_without_optics_map` — subroutine
- `test_write_project` — subroutine
- `write_upstream_project` — function

---
## Module: simple_stream_stage_sieve

Files:
- `main/stream/stages/simple_stream_stage_sieve.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_defs_environment`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_error`
- `simple_fileio`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_stream_particle_sieving`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_oris`
- `simple_parameters`
- `simple_ptcl_sieve`
- `simple_qsys_env`
- `simple_rec_list`
- `simple_sp_project`
- `simple_stream_gui_senders`
- `simple_stream_pipe`
- `simple_stream_state`
- `simple_stream_utils`
- `simple_stream_watcher`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_timer`

Public symbols:
- `stream_stage_sieve` — type

Private symbols:
- `attach_upstream` — subroutine
- `finalize` — subroutine
- `finished` — function
- `import_projects` — subroutine
- `init_gui` — subroutine
- `init_params` — subroutine
- `init_queue` — subroutine
- `iterate` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `read_mask_diameter` — subroutine
- `restore_imports` — subroutine
- `resumable` — function
- `send_latest_cavgs` — subroutine
- `send_status` — subroutine
- `start_sieve` — subroutine
- `update_final_ingestion` — subroutine

Local symbols (enclosing procedure scope):
- `kill::release_heavy` — subroutine
- `restore_imports::any_set` — function

---
## Module: simple_stream_stage_sieve_tester

Files:
- `main/stream/stages/simple_stream_stage_sieve_tester.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_gui_metadata_stream_particle_sieving`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_oris`
- `simple_ptcl_sieve`
- `simple_sp_project`
- `simple_stream_pipe`
- `simple_stream_stage_sieve`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_stage_sieve_tests` — subroutine

Private symbols:
- `close_loopback` — subroutine
- `make_test_stage` — subroutine
- `make_upstream` — subroutine
- `open_loopback` — subroutine
- `set_test_cline` — subroutine
- `test_final_ingestion_follows_marker` — subroutine
- `test_finished` — subroutine
- `test_import_projects` — subroutine
- `test_init_params` — subroutine
- `test_iterate_waits` — subroutine
- `test_read_mask_diameter` — subroutine
- `test_restart_removes_term_stream` — subroutine
- `test_restart_resumes_sieve` — subroutine
- `test_restore_and_attach` — subroutine
- `test_send_status` — subroutine
- `test_start_sieve` — subroutine
- `write_chunked_mics` — subroutine
- `write_completed_set` — function
- `write_moldiam` — subroutine

---
## Module: simple_stream_stage_solve3D

Files:
- `main/stream/stages/simple_stream_stage_solve3D.f90`

Uses:
- `simple_cavg_quality_model`
- `simple_cmdline`
- `simple_defs`
- `simple_defs_environment`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_fileio`
- `simple_gui_metadata_cavg2d`
- `simple_gui_metadata_stream_snapshot`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_gui_metadata_vol3d`
- `simple_imghead`
- `simple_math_ft`
- `simple_oris_utils`
- `simple_parameters`
- `simple_qsys_async_job`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_qsys_job_record`
- `simple_rec_list`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_state_assignment`
- `simple_stream_gui_senders`
- `simple_stream_pipe`
- `simple_stream_state`
- `simple_stream_utils`
- `simple_stream_watcher`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `stream_stage_solve3D` — type

Private symbols:
- `accept_best` — subroutine
- `accept_refs` — subroutine
- `adopt_job_2D` — function
- `advance_jobs` — subroutine
- `apply_gui_updates` — subroutine
- `attach_upstream` — subroutine
- `best_attempt` — function
- `finalize` — subroutine
- `finish_assign` — subroutine
- `finish_cavgs3D` — subroutine
- `finish_clsexp` — subroutine
- `finish_expand` — subroutine
- `finish_refine` — subroutine
- `finish_solve2D_first` — subroutine
- `finished` — function
- `fit_deselected` — function
- `flag_model_rows` — subroutine
- `import_sets` — subroutine
- `in_first_set` — function
- `init_gui` — subroutine
- `init_params` — subroutine
- `init_queue` — subroutine
- `is_final_publication` — function
- `iterate` — subroutine
- `job_done` — subroutine
- `job_failed` — subroutine
- `kill` — subroutine
- `merge_publication` — subroutine
- `n_refs` — function
- `new` — subroutine
- `next_attempt` — subroutine
- `next_job` — function
- `prepare_refine_seeds` — subroutine
- `prune_quality_dirs` — subroutine
- `prune_run_dirs` — subroutine
- `publication_problem` — function
- `read_state_fsc` — function
- `record_run_publication` — subroutine
- `register_attempt` — function
- `restore_refs` — function
- `rows_problem` — function
- `select_first_set_classes` — subroutine
- `send_refs` — subroutine
- `send_state_volume` — subroutine
- `send_status` — subroutine
- `send_volumes` — subroutine
- `set_job_resources` — subroutine
- `stack_index` — function
- `start_assign` — subroutine
- `start_cavgs3D` — subroutine
- `start_clsexp` — subroutine
- `start_expand` — subroutine
- `start_fallback` — subroutine
- `start_refine` — subroutine
- `start_solve2D_first` — subroutine
- `take_assignment` — function
- `take_cavgs` — subroutine
- `take_first_publication` — subroutine
- `take_mskdiam` — subroutine
- `take_refine_result` — function
- `watch_sets` — subroutine
- `write_assignment_cavgs` — subroutine
- `write_snapshot` — subroutine
- `write_stage_project` — subroutine

Local symbols (enclosing procedure scope):
- `import_sets::add_rejected` — subroutine
- `import_sets::is_rejected` — function

---
## Module: simple_stream_stage_solve3D_tester

Files:
- `main/stream/stages/simple_stream_stage_solve3D_tester.f90`

Uses:
- `iso_c_binding`
- `json_module`
- `simple_cmdline`
- `simple_defs_fname`
- `simple_defs_stream`
- `simple_fileio`
- `simple_gui_metadata_stream_snapshot`
- `simple_gui_metadata_stream_solve3d_multistate`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_types`
- `simple_gui_metadata_utils`
- `simple_gui_metadata_vol3d`
- `simple_image`
- `simple_math_ft`
- `simple_ori`
- `simple_rec_list`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_strategy3d_cont`
- `simple_stream_pipe`
- `simple_stream_stage_solve3d`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_stream_stage_solve3D_lib_tests` — subroutine
- `run_all_stream_stage_solve3D_tests` — subroutine

Private symbols:
- `add_cavgs_entry` — subroutine
- `check_final_publication` — subroutine
- `close_loopback` — subroutine
- `make_cavgs3D_result` — subroutine
- `make_set` — subroutine
- `make_test_stage` — subroutine
- `open_loopback` — subroutine
- `set_test_cline` — subroutine
- `test_assignment` — subroutine
- `test_failed_attempts` — subroutine
- `test_finished` — subroutine
- `test_first_publication` — subroutine
- `test_init_params` — subroutine
- `test_iterate_waits` — subroutine
- `test_local_refinement` — subroutine
- `test_merge_publications` — subroutine
- `test_mskdiam_from_each_publication` — subroutine
- `test_reference_docking` — subroutine
- `test_reference_volumes` — subroutine
- `test_restart_removes_term_stream` — subroutine
- `test_retention` — subroutine
- `test_rows_problem` — subroutine
- `test_rules` — subroutine
- `test_send_refs` — subroutine
- `test_send_status` — subroutine
- `test_send_volumes` — subroutine
- `test_snapshot3D` — subroutine
- `test_take_cavgs` — subroutine
- `test_watch_order_and_mskdiam` — subroutine
- `write_export` — subroutine

Local symbols (enclosing procedure scope):
- `test_assignment::make_matching` — subroutine
- `test_reference_docking::make_docking_result` — subroutine
- `test_reference_docking::write_blob_volume` — subroutine
- `test_snapshot3D::receive_report` — function
- `test_snapshot3D::send_request` — subroutine
- `test_take_cavgs::write_stack` — subroutine

---
## Module: simple_stream_state

Files:
- `main/stream/shared/simple_stream_state.f90`

---
## Module: simple_stream_tester

Files:
- `main/stream/simple_stream_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_cmdline`
- `simple_commanders_pick`
- `simple_commanders_stream_p02_assign_optics`
- `simple_core_module_api`
- `simple_image`
- `simple_sp_project`
- `simple_test_utils`

Public symbols:
- `run_all_stream_optics_tests` — subroutine
- `run_all_stream_pick_extract_tests` — subroutine
- `run_all_stream_pickrefs_tests` — subroutine

Private symbols:
- `distinct_position_match` — function
- `test_assign_optics_two_shift_clusters` — subroutine
- `test_make_pickrefs_expansion` — subroutine
- `test_pick_extract_three_particles` — subroutine

---
## Module: simple_stream_utils

Files:
- `main/stream/simple_stream_utils.f90`

Uses:
- `iso_c_binding`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_parameters`
- `simple_qsys_env`
- `simple_rec_list`
- `simple_sp_project`
- `simple_syslib`

Public symbols:
- `create_stream_project` — subroutine
- `get_latest_optics_map_id` — function
- `import_new_projects` — subroutine
- `init_stream_qenv` — subroutine
- `log_rss` — subroutine
- `terminate_stream` — subroutine
- `upstream_done` — function

---
## Module: simple_stream_watcher

Files:
- `main/stream/simple_stream_watcher.f90`

Uses:
- `simple_core_module_api`
- `simple_progress`

Public symbols:
- `sniff_folders_SJ` — subroutine
- `stream_watcher` — type
- `workout_directory_structure` — subroutine

Private symbols:
- `add2history_1` — subroutine
- `add2history_2` — subroutine
- `add2watchdirs` — subroutine
- `clear_history` — subroutine
- `constructor` — function
- `detect_and_add_dirs` — subroutine
- `history_pos` — function
- `is_past` — function
- `kill` — subroutine
- `watch` — subroutine
- `watchdirs` — subroutine

---
## Module: simple_stream_watcher_tester

Files:
- `main/stream/simple_stream_watcher_tester.f90`

Uses:
- `simple_fileio`
- `simple_stream_watcher`
- `simple_string`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_stream_watcher_tests` — subroutine

Private symbols:
- `test_history_growth` — subroutine
- `test_history_lookup` — subroutine
- `test_rate_after_restored_history` — subroutine

---
## Module: simple_string

Files:
- `utils/text/simple_string.f90`

Uses:
- `simple_defs`
- `simple_defs_string`
- `simple_error`

Public symbols:
- `string` — type

Private symbols:
- `append` — function
- `assign` — subroutine
- `char2str` — function
- `dble2str` — function
- `ends_with_substr_char` — function
- `ends_with_substr_str` — function
- `has_substr_char` — function
- `has_substr_str` — function
- `int2str` — function
- `is_allocated` — function
- `is_blank` — function
- `is_equal_char` — function
- `is_equal_str` — function
- `is_not_equal_char` — function
- `is_not_equal_str` — function
- `kill` — subroutine
- `print` — subroutine
- `raw` — function
- `readfile` — subroutine
- `readline` — subroutine
- `real2str` — function
- `set_dble` — subroutine
- `set_int` — subroutine
- `set_real` — subroutine
- `strlen` — function
- `strlen_trim` — function
- `substr_ind_char` — function
- `substr_ind_str` — function
- `substr_remove` — function
- `substr_replace` — subroutine
- `to_char_full` — function
- `to_char_range` — function
- `to_dble` — function
- `to_fnv1a_hash64` — function
- `to_int` — function
- `to_real` — function
- `to_static_full` — subroutine
- `to_static_range` — subroutine
- `writeline` — subroutine

Local symbols (enclosing procedure scope):
- `substr_remove::rm_substr` — subroutine
- `substr_replace::rpl_substr` — subroutine

---
## Module: simple_string_tester

Files:
- `utils/text/simple_string_tester.f90`

Uses:
- `simple_ansi_ctrls`
- `simple_defs`
- `simple_string`
- `simple_string_utils`
- `simple_test_utils`

Public symbols:
- `run_all_string_tests` — subroutine

Private symbols:
- `test_ansi_format_str` — subroutine
- `test_append_operator` — subroutine
- `test_assign_and_ctor` — subroutine
- `test_default_and_kill` — subroutine
- `test_ends_with` — subroutine
- `test_equality_operators` — subroutine
- `test_has_substr` — subroutine
- `test_is_blank_and_is_allocated` — subroutine
- `test_list_of_ints2arr` — subroutine
- `test_numeric_conversion` — subroutine
- `test_readfile` — subroutine
- `test_readline_empty_line` — subroutine
- `test_readline_eof_behavior` — subroutine
- `test_readline_long_line` — subroutine
- `test_readline_multiple_lines` — subroutine
- `test_readline_writeline` — subroutine
- `test_strlen` — subroutine
- `test_substr_ind` — subroutine
- `test_substr_remove` — subroutine
- `test_substr_replace` — subroutine
- `test_to_char` — subroutine
- `test_to_fnv1a_hash64` — subroutine
- `test_writeline_unallocated` — subroutine

---
## Module: simple_string_utils

Files:
- `utils/text/simple_string_utils.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_defs_string`
- `simple_error`
- `simple_string`

Public symbols:
- `append_limited_char` — subroutine
- `cast_str_types` — function
- `char_is_a_letter` — function
- `char_is_a_number` — function
- `cnt_recs_per_line` — function
- `compact` — subroutine
- `csv_field` — function
- `findloc_str_char` — function
- `findloc_str_str` — function
- `first_blank` — function
- `first_non_blank` — function
- `fortran_logical` — function
- `fortran_quote` — function
- `fortran_symbol_from_string` — function
- `int2str` — function
- `int2str_pad` — function
- `is_fortran_symbol_char` — function
- `is_fortran_symbol_start` — function
- `lex_sort_char` — subroutine
- `lex_sort_str` — subroutine
- `list_of_ints2arr` — function
- `list_of_strs2arr` — function
- `lowercase` — function
- `map_str_nrs` — function
- `parsestr` — subroutine
- `realdp2str` — function
- `realsp2str` — function
- `removepunct` — subroutine
- `removesp` — subroutine
- `spaces` — function
- `split_char` — subroutine
- `split_str` — subroutine
- `split_string` — subroutine
- `str2format_dp` — subroutine
- `str2format_real` — subroutine
- `str2int_char` — function
- `str2int_char_stat` — function
- `str2int_str` — function
- `str2int_str_stat` — function
- `str2real_char` — function
- `str2real_str` — function
- `str_has_substr` — function
- `str_is_blank` — function
- `str_is_comment` — function
- `str_is_true` — function
- `str_pad` — function
- `strings_are_equal` — function
- `to_cstring` — function
- `uppercase` — function

Local symbols (enclosing procedure scope):
- `lex_sort_char::key_less` — function
- `lex_sort_char::merge_idx` — subroutine
- `lex_sort_char::mergesort_idx` — subroutine
- `lex_sort_char::norm_key` — function
- `lex_sort_str::key_less` — function
- `lex_sort_str::merge_idx` — subroutine
- `lex_sort_str::mergesort_idx` — subroutine
- `lex_sort_str::norm_less` — function
- `list_of_ints2arr::add_entry` — subroutine
- `list_of_strs2arr::add_entry` — subroutine

---
## Module: simple_sym

Files:
- `main/ori/simple_sym.f90`

Uses:
- `simple_math`
- `simple_ori`
- `simple_ori_api`
- `simple_oris`

Public symbols:
- `is_valid_pointgroup` — function
- `sym` — type

Private symbols:
- `apply2all` — subroutine
- `apply_euls` — subroutine
- `apply_ori` — subroutine
- `apply_sym_with_shift` — subroutine
- `build_eullims` — subroutine
- `build_refspiral` — subroutine
- `constructor` — function
- `find_closest_proj` — function
- `get_all_subgrps_descr` — function
- `get_eullims` — function
- `get_nsubgrp` — function
- `get_nsym` — function
- `get_pgrp` — function
- `get_subgrp` — function
- `get_subgrp_descr` — function
- `get_sym_rmat` — subroutine
- `get_symori` — subroutine
- `has_subgrp` — function
- `is_circular` — function
- `is_dihedral` — function
- `is_icosahedral` — function
- `is_platonic` — function
- `kill` — subroutine
- `make_c_and_d` — subroutine
- `make_i_relion` — subroutine
- `make_i_spider` — subroutine
- `make_o` — subroutine
- `make_t` — subroutine
- `matcreate` — function
- `matextract` — function
- `nearest_proj_neighbors_k` — subroutine
- `nearest_proj_neighbors_thres` — subroutine
- `nearest_sym_neighbors` — subroutine
- `new` — subroutine
- `rnd_euler_full` — subroutine
- `rnd_euler_near` — subroutine
- `rnd_euler_near_projinpl` — subroutine
- `rot_to_asym` — subroutine
- `rotall_to_asym` — subroutine
- `set_subgrps` — subroutine
- `sym_dist_table` — subroutine
- `sym_dists_euls` — subroutine
- `sym_dists_ori` — subroutine
- `symrandomize` — subroutine
- `within_asymunit` — function
- `within_platonic_asymunit` — function

Local symbols (enclosing procedure scope):
- `build_refspiral::gen_c1` — subroutine
- `set_subgrps::fmtsymstr` — function
- `set_subgrps::getevensym` — subroutine
- `within_platonic_asymunit::e2_bound` — function
- `within_platonic_asymunit::e2_upper_bound` — function

---
## Module: simple_sym_tester

Files:
- `main/ori/simple_sym_tester.f90`

Uses:
- `simple_defs`
- `simple_linalg`
- `simple_ori`
- `simple_oris`
- `simple_sym`
- `simple_test_utils`

Public symbols:
- `run_all_sym_tests` — subroutine

Private symbols:
- `det3` — function
- `matches_some_op` — function
- `related_by_symop` — function
- `test_apply` — subroutine
- `test_build_refspiral` — subroutine
- `test_construction_and_order` — subroutine
- `test_eullims` — subroutine
- `test_is_valid_pointgroup` — subroutine
- `test_neighbour_searches` — subroutine
- `test_rnd_euler_within_limits` — subroutine
- `test_rot_to_asym` — subroutine
- `test_subgroups` — subroutine
- `test_sym_dist_table` — subroutine
- `test_symops_form_a_group` — subroutine
- `test_symrandomize` — subroutine
- `within_lims` — function

Local symbols (enclosing procedure scope):
- `test_build_refspiral::check_spiral` — subroutine
- `test_eullims::check_lims` — subroutine
- `test_neighbour_searches::insertion_sort` — subroutine
- `test_rnd_euler_within_limits::check_rnd` — subroutine
- `test_rot_to_asym::check_c1_south_pole` — subroutine
- `test_rot_to_asym::check_rot` — subroutine
- `test_symops_form_a_group::check_group` — subroutine

---
## Module: simple_symanalyzer

Files:
- `main/volume/simple_symanalyzer.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_parameters`
- `simple_projector`
- `simple_simple_volinterp`
- `simple_volpft_symsrch`

Public symbols:
- `print_subgroups` — subroutine
- `symmetrize_map` — subroutine
- `symmetry_tester` — subroutine

Private symbols:
- `eval_point_groups` — subroutine
- `sym_stats` — type

Local symbols (enclosing procedure scope):
- `eval_point_groups::find_symaxis` — subroutine
- `eval_point_groups::symaverage` — subroutine

---
## Module: simple_syslib

Files:
- `fileio/simple_syslib.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `simple_defs`
- `simple_defs_environment`
- `simple_defs_fname`
- `simple_error`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `c_atomic_replace_path` — function
- `c_sync_directory_path` — function
- `c_sync_file_path` — function
- `del_file` — subroutine
- `dir_exists` — function
- `exec_cmdline` — subroutine
- `fadvise_willneed` — function
- `file_exists` — function
- `find_next_int_dir_prefix` — function
- `fs_avail_bytes` — function
- `fs_avail_bytes_c` — function
- `fs_is_rotational` — function
- `fs_is_rotational_c` — function
- `fsync` — function
- `get_absolute_pathname` — function
- `get_current_rss_bytes` — function
- `get_node_mem_in_use_bytes` — function
- `get_node_mem_total_bytes` — function
- `get_peak_rss_bytes` — function
- `get_process_id` — function
- `get_sysinfo` — function
- `io_name_of` — function
- `io_nstreams_for_device` — function
- `io_read_nstreams` — function
- `io_set_faults` — subroutine
- `is_file_open` — function
- `is_io` — function
- `is_open` — function
- `isdir` — function
- `list_dirs` — function
- `makedir` — function
- `memory_monitor_report_c` — function
- `memory_monitor_start_c` — function
- `memory_monitor_stop_c` — function
- `mkdir` — function
- `open_fd_count` — function
- `pread_full` — function
- `prefetch_drain` — subroutine
- `prefetch_free` — subroutine
- `prefetch_new` — function
- `prefetch_submit` — function
- `prefetch_wait` — function
- `prefetch_window` — function
- `print_slurm_env` — subroutine
- `redirect_stdout_stderr` — subroutine
- `regexp_match` — function
- `removedir` — function
- `resolve_project_files` — subroutine
- `restore_stdout_stderr` — subroutine
- `rmdir` — function
- `ro_close` — function
- `ro_open` — function
- `ro_open_c` — function
- `safe_flush` — subroutine
- `set_open_files_limit` — function
- `simple_abspath` — function
- `simple_atomic_replace` — subroutine
- `simple_chdir` — subroutine
- `simple_chmod` — function
- `simple_file_stat` — subroutine
- `simple_getcwd` — subroutine
- `simple_getenv` — function
- `simple_list_dirs` — function
- `simple_list_files` — subroutine
- `simple_list_files_regexp` — subroutine
- `simple_mkdir` — subroutine
- `simple_redirect_output_c` — function
- `simple_rename` — subroutine
- `simple_restore_output_c` — function
- `simple_rmdir` — subroutine
- `simple_rmfile` — subroutine
- `simple_sync_dir` — subroutine
- `simple_sync_file` — subroutine
- `simple_touch` — subroutine
- `sscanf` — function
- `symlink` — function
- `syslib_c2fortran_string` — subroutine
- `syslib_symlink` — subroutine
- `touch` — function
- `unlink` — function
- `wait_pid` — function

Private symbols:
- `syslib_c_path` — subroutine

Local symbols (enclosing procedure scope):
- `exec_cmdline::raise_sys_error` — subroutine

---
## Module: simple_syslib_tester

Files:
- `fileio/simple_syslib_tester.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`

Public symbols:
- `run_all_syslib_tests` — subroutine

Private symbols:
- `setup_test_env` — subroutine
- `teardown_test_env` — subroutine
- `test_exec_cmdline` — subroutine
- `test_find_next_int_dir_prefix` — subroutine
- `test_get_current_rss_bytes` — subroutine
- `test_get_peak_rss_bytes` — subroutine
- `test_get_process_id` — subroutine
- `test_getenv_wrapper` — subroutine
- `test_is_file_open` — subroutine
- `test_is_io_and_is_open` — subroutine
- `test_print_slurm_env_smoke` — subroutine
- `test_project_discovery` — subroutine
- `test_simple_abspath` — subroutine
- `test_simple_atomic_replace_and_sync` — subroutine
- `test_simple_file_stat` — subroutine
- `test_simple_getcwd_and_chdir` — subroutine
- `test_simple_list_dirs` — subroutine
- `test_simple_list_files` — subroutine
- `test_simple_list_files_regexp` — subroutine
- `test_simple_mkdir_dir_exists_rmdir` — subroutine
- `test_simple_rename_and_del_file` — subroutine
- `test_simple_rmfile` — subroutine
- `test_simple_touch_and_file_exists` — subroutine
- `test_syslib_c2fortran_string` — subroutine
- `test_syslib_symlink` — subroutine

---
## Module: simple_tent_smooth

Files:
- `utils/filter/simple_tent_smooth.f90`

Public symbols:
- `box_filter_1d` — subroutine
- `box_filter_x` — subroutine
- `box_filter_y` — subroutine
- `box_filter_z` — subroutine
- `clamp_index` — function
- `is_zero` — function
- `is_zero_slab` — function
- `tent_smooth_3d` — subroutine

---
## Module: simple_test_dispatch

Files:
- `main/commanders/test/simple_test_dispatch.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `simple_cmdline`
- `simple_defs`
- `simple_error`
- `simple_forked_process`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `exec_state_name` — function
- `next_to_start` — function
- `proc_available` — function
- `set_test_environment` — subroutine
- `test_dispatcher` — type
- `test_execution` — type
- `test_route` — subroutine

Private symbols:
- `add` — subroutine
- `clock_ticks` — function
- `command_line_of` — subroutine
- `dispatch_child` — type
- `execute_child` — subroutine
- `get_execution` — function
- `get_nexecs` — function
- `kill` — subroutine
- `member` — function
- `new` — subroutine
- `on_child_term` — subroutine
- `on_interrupt` — subroutine
- `page_bytes` — function
- `processes_left` — function
- `put_env` — subroutine
- `read_mark` — function
- `read_stat` — function
- `readlink_self` — function
- `resource_budgets` — subroutine
- `run` — subroutine
- `scan_proc` — subroutine

Local symbols (enclosing procedure scope):
- `clock_ticks::c_sysconf` — function
- `execute_child::c_prctl` — function
- `on_child_term::c_exit_now` — subroutine
- `page_bytes::c_sysconf` — function
- `put_env::c_setenv` — function
- `readlink_self::c_readlink` — function
- `run::delete_outputs` — subroutine
- `run::finish_execution` — subroutine
- `run::sample_sessions` — subroutine
- `run::signal_execution` — subroutine
- `run::start_execution` — subroutine

---
## Module: simple_test_dispatch_tester

Files:
- `main/commanders/test/simple_test_dispatch_tester.f90`

Uses:
- `iso_c_binding`
- `iso_fortran_env`
- `simple_cmdline`
- `simple_string`
- `simple_syslib`
- `simple_test_dispatch`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_test_dispatch_tests` — subroutine

Private symbols:
- `group_gone` — function
- `hang` — subroutine
- `left_after_wait` — function
- `probe_route` — subroutine
- `test_interrupt` — subroutine
- `test_outcomes` — subroutine

---
## Module: simple_test_exec_api

Files:
- `main/apis/simple_test_exec_api.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_exec_helpers`
- `simple_jiffys`
- `simple_test_exec_class`
- `simple_test_exec_highlevel`
- `simple_test_exec_single`
- `simple_test_exec_stream`
- `simple_ui`
- `simple_ui_program`

---
## Module: simple_test_exec_class

Files:
- `main/exec/simple_test_exec_class.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_class`

Public symbols:
- `exec_test_class_commander` — subroutine

---
## Module: simple_test_exec_highlevel

Files:
- `main/exec/simple_test_exec_highlevel.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_highlevel`

Public symbols:
- `exec_test_highlevel_commander` — subroutine

---
## Module: simple_test_exec_single

Files:
- `main/exec/simple_test_exec_single.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_single`

Public symbols:
- `exec_test_single_commander` — subroutine

---
## Module: simple_test_exec_stream

Files:
- `main/exec/simple_test_exec_stream.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_stream`

Public symbols:
- `exec_test_stream_commander` — subroutine

---
## Module: simple_test_gate

Files:
- `main/commanders/test/simple_test_gate.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_string`

Public symbols:
- `test_gate` — type

Private symbols:
- `check` — subroutine
- `kill` — subroutine
- `metric` — subroutine
- `new` — subroutine
- `passed` — function
- `report` — subroutine
- `success` — subroutine
- `write_row` — subroutine

---
## Module: simple_test_identity

Files:
- `main/commanders/test/simple_test_identity.f90`

Uses:
- `iso_fortran_env`
- `simple_fileio`
- `simple_oris`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `compare_run_dirs` — subroutine
- `manifest_line_kept` — function
- `same_bytes` — function
- `same_manifest` — function

Private symbols:
- `has_sigma2_magic` — function
- `same_file` — function
- `same_project_records` — function
- `same_text_lines` — function

Local symbols (enclosing procedure scope):
- `compare_run_dirs::has` — function
- `compare_run_dirs::record` — subroutine
- `compare_run_dirs::relative` — function
- `same_manifest::next_kept` — subroutine
- `same_project_records::same_oris` — function
- `same_project_records::segment` — function

---
## Module: simple_test_profile_stats

Files:
- `main/commanders/test/simple_test_profile_stats.f90`

Uses:
- `ieee_arithmetic`
- `simple_stat`

Public symbols:
- `compare_profiles` — subroutine
- `gate_row` — type
- `meets_success_floor` — function
- `metric_dist` — type
- `profile_table` — subroutine
- `read_gate_table` — subroutine
- `read_profile` — subroutine
- `robust_spread` — function
- `summarize_metric` — function
- `write_profile` — subroutine

Private symbols:
- `split_tabs` — subroutine

---
## Module: simple_test_registry

Files:
- `main/commanders/test/simple_test_registry.f90`

Public symbols:
- `find_test_entry` — function
- `get_test_registry` — function
- `test_entry` — type

Private symbols:
- `outputs_of` — function

---
## Module: simple_test_state_fixture

Files:
- `main/commanders/test/simple_test_state_fixture.f90`

Uses:
- `simple_atoms`
- `simple_cmdline`
- `simple_commanders_sim`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_flex_pca_application_tester`
- `simple_image`
- `simple_imghead`
- `simple_molecule_data`
- `simple_ori`
- `simple_ori_utils`
- `simple_oris`
- `simple_simple_volinterp`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_truth_metrics`
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `betagal_frame` — subroutine
- `point_group_check` — subroutine
- `state_fixture` — type

Private symbols:
- `box` — function
- `build` — subroutine
- `key_text` — function
- `kill` — subroutine
- `make_betagal_truths` — subroutine
- `make_truths` — subroutine
- `mskdiam` — function
- `new` — subroutine
- `project_fname` — function
- `read_pgrps` — subroutine
- `simulate` — subroutine
- `smpd` — function
- `stack_fname` — function
- `subset` — function
- `truth_doc` — function
- `truth_fname` — function
- `write_lowpassed` — subroutine
- `write_pgrps` — subroutine
- `write_project` — subroutine
- `write_text` — subroutine

Local symbols (enclosing procedure scope):
- `betagal_frame::unit` — function
- `make_betagal_truths::assert_in_frame` — subroutine

---
## Module: simple_test_tools_tester

Files:
- `main/commanders/test/simple_test_tools_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_string`
- `simple_syslib`
- `simple_test_dispatch`
- `simple_test_identity`
- `simple_test_profile_stats`
- `simple_test_registry`
- `simple_test_utils`

Public symbols:
- `run_all_test_tools_tests` — subroutine

Private symbols:
- `test_gate_tables` — subroutine
- `test_manifest_rules` — subroutine
- `test_masked_bytes` — subroutine
- `test_packing` — subroutine
- `test_registry` — subroutine
- `test_robust_spread` — subroutine
- `test_run_dir_identity` — subroutine
- `test_summarize_metric` — subroutine
- `write_bytes` — subroutine
- `write_table` — subroutine

---
## Module: simple_test_truth_metrics

Files:
- `main/commanders/test/simple_test_truth_metrics.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_dock_vols`
- `simple_error`
- `simple_estimate_ssnr`
- `simple_fileio`
- `simple_image`
- `simple_imghead`
- `simple_linalg`
- `simple_ori`
- `simple_oris`
- `simple_rnd`
- `simple_srch_sort_loc`
- `simple_stat`
- `simple_string`
- `simple_sym`
- `simple_syslib`
- `simple_test_gate`
- `simple_test_utils`
- `simple_volpair_metrics`

Public symbols:
- `accuracy_null` — subroutine
- `add_gaussian_blob` — subroutine
- `adjusted_rand_index` — function
- `axis_operator` — function
- `carry_pose` — subroutine
- `collapsed_accuracy` — function
- `compare_to_truth` — subroutine
- `discrimination` — subroutine
- `dock_both_hands` — subroutine
- `farthest_truth_error` — function
- `finite_min` — function
- `group_errors` — function
- `group_operators` — function
- `mapped_accuracy` — function
- `maps_distinct` — function
- `match_maps` — function
- `match_states` — function
- `offdiag` — function
- `pair_pose_error` — function
- `population_errors` — function
- `rot_angle` — function
- `state_recalls` — function
- `validate_reconstructed_volume` — subroutine

Private symbols:
- `confusion` — function
- `masked_corr` — function

Local symbols (enclosing procedure scope):
- `match_maps::try` — subroutine
- `match_states::try` — subroutine
- `validate_reconstructed_volume::geometry_matches` — function

---
## Module: simple_test_truth_metrics_tester

Files:
- `main/commanders/test/simple_test_truth_metrics_tester.f90`

Uses:
- `ieee_arithmetic`
- `simple_defs`
- `simple_image`
- `simple_ori`
- `simple_ori_utils`
- `simple_oris`
- `simple_projector`
- `simple_simple_volinterp`
- `simple_stat`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_test_gate`
- `simple_test_truth_metrics`
- `simple_test_utils`

Public symbols:
- `run_all_truth_docking_tests` — subroutine
- `run_all_truth_metrics_tests` — subroutine

Private symbols:
- `blob_volume` — subroutine
- `proj_corr` — function
- `project` — subroutine
- `test_ari` — subroutine
- `test_discrimination` — subroutine
- `test_docking_frame` — subroutine
- `test_frame_rule` — subroutine
- `test_group_error` — subroutine
- `test_map_mapping` — subroutine
- `test_matching` — subroutine
- `test_nonfinite` — subroutine
- `test_null` — subroutine

---
## Module: simple_test_ui_class

Files:
- `main/ui/simple_test/simple_test_ui_class.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_class_programs` — subroutine
- `new_forked_process` — subroutine
- `new_lib_cart_align3D` — subroutine
- `new_lib_heterogeneity` — subroutine
- `new_lib_reconstruction` — subroutine
- `new_lib_single` — subroutine
- `new_lib_stream` — subroutine
- `new_openmp_offload` — subroutine
- `new_unit_cart_align3D` — subroutine
- `new_unit_core` — subroutine
- `new_unit_heterogeneity` — subroutine
- `new_unit_image` — subroutine
- `new_unit_ipc` — subroutine
- `new_unit_motion` — subroutine
- `new_unit_numerics` — subroutine
- `new_unit_ori` — subroutine
- `new_unit_parallel` — subroutine
- `new_unit_pftc_align2D3D` — subroutine
- `new_unit_project` — subroutine
- `new_unit_reconstruction` — subroutine
- `new_unit_single` — subroutine
- `new_unit_stream` — subroutine
- `new_unit_ui` — subroutine
- `new_units` — subroutine

---
## Module: simple_test_ui_highlevel

Files:
- `main/ui/simple_test/simple_test_ui_highlevel.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_highlevel_programs` — subroutine
- `new_cont_refine3D_1jxy` — subroutine
- `new_cont_refine3D_bench` — subroutine
- `new_flex_pca_blobs` — subroutine
- `new_identity` — subroutine
- `new_mini_stream` — subroutine
- `new_multistate_workflows` — subroutine
- `new_nightly` — subroutine
- `new_pcg_frac_update` — subroutine
- `new_pcg_recon` — subroutine
- `new_polar_matching_bench` — subroutine
- `new_profile` — subroutine
- `new_rec3D_backends` — subroutine
- `new_simulate_particles` — subroutine
- `new_simulated_workflow` — subroutine
- `new_single_atoms_stats` — subroutine
- `new_single_workflow` — subroutine
- `new_solve3D_addon` — subroutine
- `new_solve3D_addon_snapshots` — subroutine
- `new_species_discovery` — subroutine
- `new_state_fixture` — subroutine
- `new_stream_chains` — subroutine
- `new_stream_preproc` — subroutine
- `new_write_state_weights` — subroutine
- `print_test_highlevel_programs` — subroutine

---
## Module: simple_test_utils

Files:
- `utils/simple_test_utils.f90`

Uses:
- `iso_fortran_env`
- `simple_defs`
- `simple_fileio`
- `simple_rnd`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`

Public symbols:
- `assert_char` — subroutine
- `assert_double` — subroutine
- `assert_false` — subroutine
- `assert_int` — subroutine
- `assert_real` — subroutine
- `assert_string_eq` — subroutine
- `assert_true` — subroutine
- `begin_test_suite` — subroutine
- `end_test_suite` — subroutine
- `enter_fixture` — subroutine
- `fix_test_seed` — subroutine
- `leave_fixture` — subroutine
- `offset_seed` — function
- `report_summary` — subroutine
- `reset_test_report` — subroutine
- `set_fixed_seed` — subroutine
- `set_test_seed_offset` — subroutine

Private symbols:
- `record_assertion` — subroutine
- `test_failure` — type
- `test_suite_result` — type
- `write_report_file` — subroutine
- `write_summary` — subroutine

---
## Module: simple_testfuns

Files:
- `utils/math/simple_testfuns.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `get_testfun` — subroutine
- `testfun` — function
- `testfun1` — function
- `testfun10` — function
- `testfun11` — function
- `testfun12` — function
- `testfun13` — function
- `testfun14` — function
- `testfun15` — function
- `testfun16` — function
- `testfun17` — function
- `testfun18` — function
- `testfun19` — function
- `testfun2` — function
- `testfun20` — function
- `testfun21` — function
- `testfun22` — function
- `testfun23` — function
- `testfun24` — function
- `testfun25` — function
- `testfun26` — function
- `testfun27` — function
- `testfun28` — function
- `testfun29` — function
- `testfun3` — function
- `testfun30` — function
- `testfun31` — function
- `testfun32` — function
- `testfun4` — function
- `testfun5` — function
- `testfun6` — function
- `testfun7` — function
- `testfun8` — function
- `testfun9` — function

---
## Module: simple_tifflib

Files:
- `fileio/simple_tifflib.f90`

Uses:
- `iso_c_binding`

Public symbols:
- `TIFFAllocateStripBuffer` — function
- `TIFFClose` — subroutine
- `TIFFCurrentDirectory` — function
- `TIFFfree` — function
- `TIFFGetBitsPerSample` — function
- `TIFFGetCompression` — function
- `TIFFGetField` — function
- `TIFFGetLength` — function
- `TIFFGetMaxVal` — function
- `TIFFGetMinVal` — function
- `TIFFGetRowsPerStrip` — function
- `TIFFGetSampleFormat` — function
- `TIFFGetSamplesPerPixel` — function
- `TIFFGetVersion` — function
- `TIFFGetWidth` — function
- `TIFFIsTiled` — function
- `TIFFLastDirectory` — function
- `TIFFmalloc` — function
- `TIFFMuteWarnings` — subroutine
- `TIFFNumberOfStrips` — function
- `TIFFNumDirectories` — function
- `TIFFOpen` — function
- `TIFFPrintInfo` — subroutine
- `TIFFRawStripSizer` — function
- `TIFFReadDirectory` — function
- `TIFFReadEncodedStrip` — function
- `TIFFReadRawStrip` — function
- `TIFFSetDirectory` — function
- `TIFFStripSize` — function
- `TIFFUnMuteWarnings` — subroutine

---
## Module: simple_timer

Files:
- `utils/simple_timer.f90`

Uses:
- `ifport`
- `simple_defs`
- `simple_error`
- `simple_is_check_assert`
- `simple_string_utils`

Public symbols:
- `cast_time_char` — function
- `in_timer_loop` — function
- `now` — subroutine
- `reset_timer` — subroutine
- `simple_gettime` — function
- `tdiff` — function
- `tdiff_i4` — function
- `tic` — function
- `tic_i4` — function
- `tickrate` — function
- `tickrate_i4` — function
- `timer_loop_end` — subroutine
- `timer_loop_start` — subroutine
- `timer_profile_break` — subroutine
- `timer_profile_report` — subroutine
- `timer_profile_setup` — subroutine
- `timer_profile_start` — subroutine
- `toc` — function
- `toc_i4` — function
- `tocprint` — subroutine

---
## Module: simple_timer_omp

Files:
- `utils/simple_timer_omp.f90`

Uses:
- `simple_defs`

Public symbols:
- `now_omp` — subroutine
- `reset_timer_omp` — subroutine
- `tdiff_omp` — function
- `tic_omp` — function
- `tickrate_omp` — function
- `toc_omp` — function
- `tocprint_omp` — function

---
## Module: simple_trail_chain_manifest

Files:
- `main/volume/simple_trail_chain_manifest.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `trail_chain_manifest` — type

Private symbols:
- `compatible` — function
- `get_box` — function
- `get_gen` — function
- `get_mrep` — function
- `get_ncontrib` — function
- `get_nptcls` — function
- `get_nstates` — function
- `get_size` — function
- `get_smpd` — function
- `get_state` — function
- `get_wset_id` — function
- `kill` — subroutine
- `new` — subroutine
- `read` — subroutine
- `write` — subroutine

---
## Module: simple_trail_chain_manifest_tester

Files:
- `main/volume/simple_trail_chain_manifest_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_reconstructor_pcg`
- `simple_test_utils`
- `simple_trail_chain_manifest`

Public symbols:
- `run_all_trail_chain_manifest_tests` — subroutine

Private symbols:
- `test_compatibility` — subroutine
- `test_gridding_manifest` — subroutine
- `test_pcg_chain_header` — subroutine

---
## Module: simple_tseries_graphene_subtr

Files:
- `main/nano/simple_tseries_graphene_subtr.f90`

Uses:
- `simple_atoms`
- `simple_pftc_srch_api`

Public symbols:
- `calc_peaks` — subroutine
- `init_graphene_subtr` — subroutine
- `kill_graphene_subtr` — subroutine
- `remove_lattices` — subroutine

Local symbols (enclosing procedure scope):
- `calc_peaks::interp_peak` — subroutine
- `calc_peaks::range_convention` — subroutine
- `remove_lattices::obscure_peak` — subroutine
- `remove_lattices::remove_lattice` — subroutine

---
## Module: simple_type_defs

Files:
- `defs/simple_type_defs.f90`

Uses:
- `simple_defs`
- `simple_string`

Public symbols:
- `canonical_phshift_dp` — function
- `canonical_phshift_sp` — function
- `class_sample` — type
- `clust_info` — type
- `clust_inpl` — type
- `ctfparams` — type
- `ctfvars` — type
- `fplane_type` — type
- `get_centype` — function
- `inpl_struct` — type
- `lp_crop_inf` — type
- `ptcl_ref` — type
- `sigma_array` — type
- `stats_struct` — type

---
## Module: simple_ui

Files:
- `main/ui/simple_ui.f90`

Uses:
- `json_file_module`
- `json_kinds`
- `json_module`
- `simple_ansi_ctrls`
- `simple_core_module_api`
- `simple_linked_list`
- `simple_ui_hash`
- `simple_ui_param`
- `simple_ui_params_common`
- `simple_ui_program`
- `simple_ui_simple_group`
- `simple_ui_single_group`
- `simple_ui_stream_group`
- `simple_ui_test_group`

Public symbols:
- `count_prgs_in_category` — function
- `get_prg_ptr` — subroutine
- `get_test_prg_ptr` — subroutine
- `list_simple_prgs_in_ui` — subroutine
- `list_simple_test_prgs_in_ui` — subroutine
- `list_single_prgs_in_ui` — subroutine
- `list_stream_prgs_in_ui` — subroutine
- `make_test_ui` — subroutine
- `make_ui` — subroutine
- `print_stream_ui_json` — subroutine
- `print_ui_json` — subroutine
- `validate_ui_json` — subroutine
- `validate_ui_presentation` — subroutine
- `write_ui_json` — subroutine

Private symbols:
- `print_registered_programs` — subroutine
- `validate_category_metadata` — subroutine

Local symbols (enclosing procedure scope):
- `validate_ui_presentation::validate_input_list` — subroutine

---
## Module: simple_ui_cavgproc

Files:
- `main/ui/simple/simple_ui_cavgproc.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_cavgproc_programs` — subroutine
- `new_cluster_cavgs` — subroutine
- `new_cluster_cavgs_selection` — subroutine
- `new_cluster_stack` — subroutine
- `new_match_cavgs` — subroutine
- `new_match_stacks` — subroutine
- `new_model_cavgs_rejection` — subroutine
- `new_select_clusters` — subroutine

---
## Module: simple_ui_descriptor_types

Files:
- `main/ui/simple_ui_descriptor_types.f90`

Uses:
- `simple_string`

Public symbols:
- `ui_choice` — type
- `ui_choices` — function
- `ui_choices_are_valid` — function

---
## Module: simple_ui_dock

Files:
- `main/ui/simple/simple_ui_dock.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_dock_programs` — subroutine
- `new_dock_volpair` — subroutine
- `new_volanalyze` — subroutine
- `new_volcluster` — subroutine

---
## Module: simple_ui_filter

Files:
- `main/ui/simple/simple_ui_filter.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_filter_programs` — subroutine
- `new_filter` — subroutine
- `new_nu_filt3D` — subroutine
- `new_uniform_filter2D` — subroutine
- `new_uniform_filter3D` — subroutine

---
## Module: simple_ui_hash

Files:
- `main/ui/simple_ui_hash.f90`

Uses:
- `simple_string`
- `simple_ui_program`
- `simple_vrefhash`

Public symbols:
- `ui_hash` — type

Private symbols:
- `get_ui_program` — subroutine
- `set_ui_program` — subroutine

---
## Module: simple_ui_hash_tester

Files:
- `main/ui/simple_ui_hash_tester.f90`

Uses:
- `simple_string`
- `simple_test_utils`
- `simple_ui_hash`
- `simple_ui_program`

Public symbols:
- `run_all_ui_hash_tests` — subroutine

Private symbols:
- `make_prg` — subroutine
- `test_found_is_optional` — subroutine
- `test_keys_are_trimmed` — subroutine
- `test_missing_key_is_a_typed_miss` — subroutine
- `test_overwrite_retargets_pointer` — subroutine
- `test_set_get_by_key` — subroutine
- `test_wrong_dynamic_type_is_a_typed_miss` — subroutine

---
## Module: simple_ui_heterogeneity

Files:
- `main/ui/simple/simple_ui_heterogeneity.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_heterogeneity_programs` — subroutine
- `new_classify3D_refs` — subroutine
- `new_flex_pca` — subroutine
- `new_ptcl3D_state_consensus` — subroutine
- `new_refine3D_states` — subroutine
- `new_refine3D_states_cavgs` — subroutine

---
## Module: simple_ui_image

Files:
- `main/ui/simple/simple_ui_image.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_image_programs` — subroutine
- `new_binarize` — subroutine
- `new_convert` — subroutine
- `new_ctf_correct` — subroutine
- `new_ctfops` — subroutine
- `new_normalize` — subroutine
- `new_scale` — subroutine
- `new_split_` — subroutine
- `new_split_stack` — subroutine
- `new_stack` — subroutine
- `new_stackops` — subroutine

---
## Module: simple_ui_mask

Files:
- `main/ui/simple/simple_ui_mask.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_mask_programs` — subroutine
- `new_auto_spher_mask` — subroutine
- `new_automask` — subroutine
- `new_automask2D` — subroutine
- `new_mask` — subroutine

---
## Module: simple_ui_modules

Files:
- `main/ui/simple_ui_modules.f90`

Uses:
- `simple_ansi_ctrls`
- `simple_core_module_api`
- `simple_ui_descriptor_types`
- `simple_ui_hash`
- `simple_ui_params_common`
- `simple_ui_program`
- `simple_ui_utils`
- `simple_ui_visibility`

---
## Module: simple_ui_ori

Files:
- `main/ui/simple/simple_ui_ori.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_ori_programs` — subroutine
- `new_make_oris` — subroutine
- `new_measure_projspace_angres` — subroutine
- `new_oriops` — subroutine
- `new_oristats` — subroutine
- `new_vizoris` — subroutine

---
## Module: simple_ui_other

Files:
- `main/ui/simple/simple_ui_other.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_other_programs` — subroutine
- `new_cif2pdb` — subroutine

---
## Module: simple_ui_param

Files:
- `main/ui/simple_ui_param.f90`

Uses:
- `simple_error`
- `simple_string`
- `simple_string_utils`
- `simple_ui_descriptor_types`

Public symbols:
- `choice_is_available` — function
- `finalize` — subroutine
- `infer_units` — function
- `normalize_numeric_example` — function
- `set_choices` — subroutine
- `set_generated_default` — subroutine
- `set_param_1` — subroutine
- `set_param_2` — subroutine
- `set_placeholder` — subroutine
- `standardize_placeholder` — function
- `ui_param` — type
- `ui_placeholder_is_standard` — function

---
## Module: simple_ui_params_common

Files:
- `main/ui/simple_ui_params_common.f90`

Uses:
- `simple_core_module_api`
- `simple_ui_descriptor_types`
- `simple_ui_param`

Public symbols:
- `set_ui_params` — subroutine

---
## Module: simple_ui_postprocess

Files:
- `main/ui/simple/simple_ui_postprocess.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_postprocess_programs` — subroutine
- `new_enhance` — subroutine
- `new_postprocess` — subroutine
- `new_postprocess_nu` — subroutine
- `new_postprocess_pdb` — subroutine

---
## Module: simple_ui_preproc

Files:
- `main/ui/simple/simple_ui_preproc.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_preproc_programs` — subroutine
- `new_assign_optics_groups` — subroutine
- `new_ctf_estimate` — subroutine
- `new_extract` — subroutine
- `new_gen_pspecs_and_thumbs` — subroutine
- `new_motion_correct` — subroutine
- `new_particle_sieving` — subroutine
- `new_pick` — subroutine
- `new_preprocess` — subroutine
- `new_reextract` — subroutine
- `new_refine_motion_model` — subroutine

---
## Module: simple_ui_print

Files:
- `main/ui/simple/simple_ui_print.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_print_programs` — subroutine
- `new_info_image` — subroutine
- `new_info_stktab` — subroutine
- `new_print_dose_weights` — subroutine
- `new_print_fsc` — subroutine
- `new_print_magic_boxes` — subroutine

---
## Module: simple_ui_program

Files:
- `main/ui/simple_ui_program.f90`

Uses:
- `json_module`
- `simple_ansi_ctrls`
- `simple_core_module_api`
- `simple_linked_list`
- `simple_ui_default_values`
- `simple_ui_descriptor_types`
- `simple_ui_param`
- `simple_ui_visibility`

Public symbols:
- `accepts` — function
- `add_input_num` — subroutine
- `add_input_param` — subroutine
- `add_input_str` — subroutine
- `add_program_groups_json` — subroutine
- `add_program_input_json` — subroutine
- `add_program_requirements_json` — subroutine
- `add_requirement` — subroutine
- `append_input` — subroutine
- `append_required_keys_from_list` — subroutine
- `apply_generated_default` — subroutine
- `bind_input` — subroutine
- `cardinality_text` — function
- `category_descriptor` — type
- `count_required_in_list` — function
- `create_json_entry` — subroutine
- `create_section_from_list` — subroutine
- `find_input_text` — subroutine
- `get_category` — function
- `get_display_name` — function
- `get_executable` — function
- `get_name` — function
- `get_nrequired_keys` — function
- `get_required_keys` — function
- `group_label` — function
- `has_input` — function
- `has_input_key` — function
- `input_text_for_key` — function
- `key_is_defined` — function
- `kill` — subroutine
- `list_has_input_key` — function
- `list_has_key` — function
- `new` — subroutine
- `print_cmdline` — subroutine
- `print_help` — subroutine
- `print_param_hash` — subroutine
- `print_param_list` — subroutine
- `print_requirement_groups` — subroutine
- `print_requirement_inputs` — subroutine
- `print_ui` — subroutine
- `register_group` — subroutine
- `requirements_satisfied` — function
- `requires_sp_project` — function
- `set_category` — subroutine
- `ui_activation` — type
- `ui_activation_equals_any` — function
- `ui_cli_param_choices` — function
- `ui_cli_param_summary` — function
- `ui_input_group` — type
- `ui_program` — type
- `ui_program_input` — type
- `ui_requirement_group` — type
- `write2json` — subroutine

---
## Module: simple_ui_project

Files:
- `main/ui/simple/simple_ui_project.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_project_programs` — subroutine
- `new_consolidate_particles` — subroutine
- `new_export_relion` — subroutine
- `new_export_starproject` — subroutine
- `new_extract_subproj` — subroutine
- `new_fix_projfile` — subroutine
- `new_import_boxes` — subroutine
- `new_import_cavgs` — subroutine
- `new_import_movies` — subroutine
- `new_import_particles` — subroutine
- `new_import_starproject` — subroutine
- `new_merge_projects` — subroutine
- `new_new_project` — subroutine
- `new_print_project_field` — subroutine
- `new_print_project_info` — subroutine
- `new_prune_project` — subroutine
- `new_reimport_particles` — subroutine
- `new_replace_project_field` — subroutine
- `new_select_micrographs` — subroutine
- `new_select_particles` — subroutine
- `new_selection` — subroutine
- `new_update_project` — subroutine
- `new_write_mic_filetab` — subroutine
- `new_zero_project_shifts` — subroutine

---
## Module: simple_ui_reconstruct3D

Files:
- `main/ui/simple/simple_ui_reconstruct3D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_reconstruct3D_programs` — subroutine
- `new_bootstrap_rec3D` — subroutine
- `new_reconstruct3D` — subroutine
- `new_reconstruct3D_wave` — subroutine

---
## Module: simple_ui_refine2D

Files:
- `main/ui/simple/simple_ui_refine2D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_refine2D_programs` — subroutine
- `new_bootstrap_cavgs` — subroutine
- `new_cls_expansion` — subroutine
- `new_make_cavgs` — subroutine
- `new_map_cavgs_selection` — subroutine
- `new_refine2D` — subroutine
- `new_sample_classes` — subroutine
- `new_solve2D` — subroutine
- `new_solve2D_chunks` — subroutine
- `new_solve2D_descriptor` — subroutine
- `new_unbootstrap_cavgs` — subroutine
- `new_write_classes` — subroutine

---
## Module: simple_ui_refine3D

Files:
- `main/ui/simple/simple_ui_refine3D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `add_refine3D_workflow_inputs` — subroutine
- `construct_refine3D_programs` — subroutine
- `new_refine3D` — subroutine
- `new_refine3D_auto` — subroutine
- `new_refine3D_local` — subroutine

---
## Module: simple_ui_res

Files:
- `main/ui/simple/simple_ui_res.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_res_programs` — subroutine
- `new_fsc` — subroutine
- `new_fsc_area_score` — subroutine

---
## Module: simple_ui_sim

Files:
- `main/ui/simple/simple_ui_sim.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_sim_programs` — subroutine
- `new_cif2mrc` — subroutine
- `new_pdb2mrc` — subroutine
- `new_simulate_movie` — subroutine
- `new_simulate_noise` — subroutine
- `new_simulate_particles` — subroutine
- `new_simulate_scattering` — subroutine

---
## Module: simple_ui_simple_group

Files:
- `main/ui/simple_ui_simple_group.f90`

Uses:
- `simple_ui_cavgproc`
- `simple_ui_dock`
- `simple_ui_filter`
- `simple_ui_hash`
- `simple_ui_heterogeneity`
- `simple_ui_image`
- `simple_ui_mask`
- `simple_ui_ori`
- `simple_ui_other`
- `simple_ui_postprocess`
- `simple_ui_preproc`
- `simple_ui_print`
- `simple_ui_project`
- `simple_ui_reconstruct3d`
- `simple_ui_refine2d`
- `simple_ui_refine3d`
- `simple_ui_res`
- `simple_ui_sim`
- `simple_ui_solve3d`
- `simple_ui_sym`
- `simple_ui_validate`
- `simple_ui_volume`

Public symbols:
- `add_simple_programs` — subroutine

---
## Module: simple_ui_single_group

Files:
- `main/ui/simple_ui_single_group.f90`

Uses:
- `simple_ui_hash`
- `single_ui_atom`
- `single_ui_map`
- `single_ui_nano2d`
- `single_ui_nano3d`
- `single_ui_trajectory`
- `single_ui_tseries`
- `single_ui_validate`

Public symbols:
- `add_single_programs` — subroutine

---
## Module: simple_ui_solve3D

Files:
- `main/ui/simple/simple_ui_solve3D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `add_solve3D_base_inputs` — subroutine
- `construct_solve3D_programs` — subroutine
- `new_estimate_lpstages` — subroutine
- `new_noisevol` — subroutine
- `new_solve3D` — subroutine
- `new_solve3D_addon` — subroutine
- `new_solve3D_cavgs` — subroutine
- `new_solve3D_chunks` — subroutine
- `new_solve3D_pyramid` — subroutine

---
## Module: simple_ui_stream

Files:
- `main/ui/simple/simple_ui_stream.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_stream_programs` — subroutine
- `new_assign_optics` — subroutine
- `new_gen_pickrefs` — subroutine
- `new_master` — subroutine
- `new_pick_extract` — subroutine
- `new_pool2D` — subroutine
- `new_preproc` — subroutine
- `new_sieve_cavgs` — subroutine
- `new_solve3D_stream` — subroutine

---
## Module: simple_ui_stream_group

Files:
- `main/ui/simple_ui_stream_group.f90`

Uses:
- `simple_ui_hash`
- `simple_ui_stream`

Public symbols:
- `add_stream_programs` — subroutine

---
## Module: simple_ui_sym

Files:
- `main/ui/simple/simple_ui_sym.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_symmetry_programs` — subroutine
- `new_symaxis_search` — subroutine
- `new_symmetrize_map` — subroutine
- `new_symmetry_test` — subroutine

---
## Module: simple_ui_test_group

Files:
- `main/ui/simple_ui_test_group.f90`

Uses:
- `simple_test_ui_class`
- `simple_test_ui_highlevel`
- `simple_ui_hash`

Public symbols:
- `add_test_programs` — subroutine

---
## Module: simple_ui_utils

Files:
- `main/ui/simple_ui_utils.f90`

Uses:
- `simple_error`
- `simple_ui_hash`
- `simple_ui_program`

Public symbols:
- `add_ui_program` — subroutine

---
## Module: simple_ui_validate

Files:
- `main/ui/simple/simple_ui_validation.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_validate_programs` — subroutine
- `new_check_refpick` — subroutine
- `new_mini_stream` — subroutine
- `new_model_validate` — subroutine

---
## Module: simple_ui_visibility

Files:
- `main/ui/simple_ui_visibility.f90`

Public symbols:
- `ui_visibility_is_valid` — function
- `ui_visibility_name` — function

---
## Module: simple_ui_visibility_tester

Files:
- `main/ui/simple_ui_visibility_tester.f90`

Uses:
- `simple_linked_list`
- `simple_string`
- `simple_test_utils`
- `simple_ui`
- `simple_ui_default_values`
- `simple_ui_descriptor_types`
- `simple_ui_param`
- `simple_ui_program`
- `simple_ui_visibility`

Public symbols:
- `run_all_ui_visibility_tests` — subroutine

Private symbols:
- `assert_input_binding` — subroutine
- `assert_input_choice` — subroutine
- `assert_input_placeholder` — subroutine
- `assert_input_visibility` — subroutine
- `assert_program_not_registered` — subroutine
- `assert_registered_category` — subroutine
- `assert_registered_cli_summary` — subroutine
- `assert_registered_required_keys` — subroutine
- `assert_registered_requirement` — subroutine
- `assert_registered_search_cli_summary` — subroutine
- `assert_registered_test_category` — subroutine
- `assert_true_all_valid` — subroutine
- `assert_ui_param` — subroutine
- `find_ui_param` — subroutine
- `program_has_key` — function
- `test_cont_surface_contract` — subroutine
- `test_param_descriptors` — subroutine
- `test_particle_storage_contract` — subroutine
- `test_phshift_contract` — subroutine
- `test_program_descriptor` — subroutine
- `test_refine_motion_model_contract` — subroutine
- `test_registered_programs` — subroutine
- `test_registered_test_programs` — subroutine
- `test_solve3D_trs_contract` — subroutine
- `test_visibility_levels` — subroutine

---
## Module: simple_ui_volume

Files:
- `main/ui/simple/simple_ui_volume.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_volume_programs` — subroutine
- `new_center` — subroutine
- `new_detect_calpha` — subroutine
- `new_reproject` — subroutine
- `new_volops` — subroutine

---
## Module: simple_umap

Files:
- `main/flex/simple_umap.f90`

Uses:
- `simple_core_module_api`
- `simple_kd_tree`

Public symbols:
- `umap_embed` — subroutine
- `umap_subsample` — subroutine

Private symbols:
- `fuzzy_graph` — subroutine
- `nbr_pos` — function
- `pca_init` — subroutine
- `rnd_below` — function
- `rnd_uni` — function
- `seed_rng` — function

---
## Module: simple_view_partition_sampling

Files:
- `main/strategies/search/simple_view_partition_sampling.f90`

Uses:
- `simple_clustering_utils`
- `simple_oris`
- `simple_pftc_srch_api`
- `simple_strategy2d_utils`

Public symbols:
- `cluster_cavgs_into_groups` — subroutine
- `has_clusterable_cavgs` — function
- `has_selected_cavgs` — function
- `make_class_samples` — subroutine
- `report_class_sample_coverage` — subroutine

Private symbols:
- `cluster_selected_cavgs` — subroutine
- `drop_inactive_rows` — subroutine

Local symbols (enclosing procedure scope):
- `cluster_cavgs_into_groups::write_group_table` — subroutine

---
## Module: simple_vol_pproc_policy

Files:
- `main/volume/simple_vol_pproc_policy.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `state_mask_is_compatible` — subroutine

---
## Module: simple_volanalyzer

Files:
- `main/volume/simple_volanalyzer.f90`

Uses:
- `simple_core_module_api`
- `simple_dock_vols`
- `simple_image`
- `simple_ori`
- `simple_parameters`
- `simple_simple_volinterp`

Public symbols:
- `calc_volpair_corrmat` — subroutine
- `dock_compare_volumes` — subroutine
- `init_volanalyzer` — subroutine

---
## Module: simple_volassemble_strategy

Files:
- `main/strategies/parallelization/simple_volassemble_strategy.f90`

Uses:
- `simple_bench_io`
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_frozen_accum`
- `simple_halfmap_diagnostics`
- `simple_image`
- `simple_imgfile`
- `simple_nu_filter`
- `simple_nu_state_filter`
- `simple_oris`
- `simple_parameters`
- `simple_reconstructor`
- `simple_reconstructor_pcg`
- `simple_refine3d_fnames`
- `simple_refine3d_provenance`
- `simple_state_weight_set`
- `simple_trail_chain_manifest`
- `simple_vol_pproc_policy`

Public symbols:
- `filter_pcg_nonuniform_maps` — subroutine
- `run_volassemble` — subroutine
- `volassemble_exec` — subroutine

Private symbols:
- `restore_state_from_parts` — subroutine
- `restore_timings_t` — type

Local symbols (enclosing procedure scope):
- `restore_state_from_parts::add_frozen_accumulators` — subroutine
- `restore_state_from_parts::blend_trailing_accumulators` — subroutine
- `restore_state_from_parts::calc_gridding_pair_diagnostics` — subroutine
- `restore_state_from_parts::cleanup_restore_state` — subroutine
- `restore_state_from_parts::discard_trail_chain_set` — subroutine
- `restore_state_from_parts::keep_half` — subroutine
- `restore_state_from_parts::previous_state_volume_exists` — function
- `restore_state_from_parts::read_gridding_pair_accumulators` — subroutine
- `restore_state_from_parts::reduce_partials` — subroutine
- `restore_state_from_parts::resolution_outfile_fbody` — function
- `restore_state_from_parts::resolve_fsc_txt_fname` — function
- `restore_state_from_parts::restore_eos_and_write_fsc` — subroutine
- `restore_state_from_parts::restore_gridding_pair` — subroutine
- `restore_state_from_parts::restore_merged_volume` — subroutine
- `restore_state_from_parts::set_state_filenames` — subroutine
- `restore_state_from_parts::sum_eos_after_density_correction_if_needed` — subroutine
- `restore_state_from_parts::sum_eos_before_density_correction_if_needed` — subroutine
- `restore_state_from_parts::sum_pair_into_sum_rec` — subroutine
- `restore_state_from_parts::trail_chain_component` — function
- `restore_state_from_parts::trail_chain_file_size` — function
- `restore_state_from_parts::use_static_nu_aux_replacement` — function
- `restore_state_from_parts::validate_trail_chain` — function
- `restore_state_from_parts::write_gridding_pair_accumulators` — subroutine
- `restore_state_from_parts::write_trail_chain_set` — subroutine
- `volassemble_exec::assemble_state` — subroutine
- `volassemble_exec::carry_forward_dropped_state` — subroutine
- `volassemble_exec::cleanup_context` — subroutine
- `volassemble_exec::collect_restore_timings` — subroutine
- `volassemble_exec::determine_dropped_states` — subroutine
- `volassemble_exec::determine_trailing_update_fraction` — subroutine
- `volassemble_exec::initialize_bench_timers` — subroutine
- `volassemble_exec::initialize_context` — subroutine
- `volassemble_exec::log_nu_alignment_lowpass_summary` — subroutine
- `volassemble_exec::postprocess_state` — subroutine
- `volassemble_exec::refresh_state_populations` — subroutine
- `volassemble_exec::run_state_nonuniform_filter` — subroutine
- `volassemble_exec::state_has_partials` — function
- `volassemble_exec::update_project_nu_alignment_lowpass` — subroutine
- `volassemble_exec::update_project_resolution_metadata` — subroutine
- `volassemble_exec::use_static_nu_aux_replacement` — function
- `volassemble_exec::write_benchmark` — subroutine

---
## Module: simple_volcluster

Files:
- `main/volume/simple_volcluster.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_parameters`

Public symbols:
- `calc_volcluster_dmat` — function
- `read_volcluster_volumes` — subroutine
- `write_volcluster_report` — subroutine

Private symbols:
- `kill_volumes` — subroutine

---
## Module: simple_volpair_metrics

Files:
- `main/volume/simple_volpair_metrics.f90`

Uses:
- `ieee_arithmetic`
- `simple_core_module_api`
- `simple_image`
- `simple_imghead`

Public symbols:
- `compare_volpair` — subroutine

---
## Module: simple_volpair_metrics_tester

Files:
- `main/volume/simple_volpair_metrics_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_test_utils`
- `simple_volpair_metrics`

Public symbols:
- `run_all_volpair_metrics_tests` — subroutine

Private symbols:
- `blob_map` — function
- `test_band_correlation` — subroutine
- `test_identical_maps` — subroutine
- `test_independent_maps` — subroutine
- `test_mismatched_maps` — subroutine
- `white_noise` — function
- `write_fixtures` — subroutine
- `write_map` — subroutine

---
## Module: simple_volpft_corrcalc

Files:
- `main/volume/simple_volpft_corrcalc.f90`

Uses:
- `simple_core_module_api`
- `simple_projector`

Public symbols:
- `volpft_corrcalc` — type

Private symbols:
- `corr_1` — function
- `corr_2` — function
- `corr_mag` — function
- `extract_ref` — subroutine
- `extract_target_1` — subroutine
- `extract_target_2` — subroutine
- `get_kfromto` — function
- `get_nspace` — function
- `get_nspace_nonred` — function
- `kill` — subroutine
- `new` — subroutine

---
## Module: simple_volpft_symsrch

Files:
- `main/volume/simple_volpft_symsrch.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`
- `simple_parameters`
- `simple_projector`
- `simple_volpft_corrcalc`

Public symbols:
- `volpft_srch4symaxis` — subroutine
- `volpft_symsrch_init` — subroutine

Private symbols:
- `opt4openMP` — type
- `volpft_symsrch_costfun` — function
- `volpft_symsrch_kill` — subroutine
- `volpft_symsrch_scorefun` — function

---
## Module: simple_volume_shape_tester

Files:
- `main/image/simple_volume_shape_tester.f90`

Uses:
- `simple_image`
- `simple_image_bin`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_volume_shape_tests` — subroutine

Private symbols:
- `test_vol_shape_descr_components` — subroutine
- `verify_shape` — subroutine

---
## Module: simple_vrefhash

Files:
- `utils/structs/simple_vrefhash.f90`

Uses:
- `simple_error`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `bucket_index` — function
- `clear` — subroutine
- `count` — function
- `del_char` — subroutine
- `del_str` — subroutine
- `destroy` — subroutine
- `finalize_node` — subroutine
- `find_node` — subroutine
- `fnv1a_hash_bytes` — function
- `get_ref_char` — subroutine
- `get_ref_str` — subroutine
- `has_key_char` — function
- `has_key_str` — function
- `init` — subroutine
- `keys` — function
- `keys_sorted` — function
- `set_ref_char` — subroutine
- `set_ref_str` — subroutine
- `vrefhash` — type
- `vrefhash_bucket` — type
- `vrefhash_node` — type

---
## Module: simple_vrefhash_tester

Files:
- `utils/structs/simple_vrefhash_tester.f90`

Uses:
- `simple_string`
- `simple_test_utils`
- `simple_vrefhash`

Public symbols:
- `run_all_vrefhash_tests` — subroutine

Private symbols:
- `inc` — subroutine
- `my_cfg` — type
- `my_obj` — type
- `test_collision_sanity_many_keys` — subroutine
- `test_delete` — subroutine
- `test_init_clear_destroy_and_count` — subroutine
- `test_keys_returns_all_keys` — subroutine
- `test_keys_sorted_returns_sorted_keys` — subroutine
- `test_missing_key_get_ref_returns_null` — subroutine
- `test_reference_semantics_updates_visible` — subroutine
- `test_replace_pointer_on_same_key` — subroutine
- `test_set_get_has_key_with_char_and_string_keys` — subroutine
- `test_store_multiple_dynamic_types` — subroutine

---
## Module: simple_wave_field

Files:
- `main/image/simple_wave_field.f90`

Uses:
- `iso_c_binding`
- `simple_error`
- `simple_fftw3`

Public symbols:
- `wave_field` — type
- `wave_freq` — function

Private symbols:
- `bwd` — subroutine
- `crop` — subroutine
- `crop_t` — subroutine
- `crop_t_real` — subroutine
- `embed` — subroutine
- `embed_t` — subroutine
- `fwd` — subroutine
- `fwd_real_crop` — subroutine
- `get_disk` — subroutine
- `get_n` — function
- `kill` — subroutine
- `new` — subroutine
- `ptr` — function
- `put_disk` — subroutine
- `rptr` — function
- `zero` — subroutine

---
## Module: simple_wave_field_tester

Files:
- `main/image/simple_wave_field_tester.f90`

Uses:
- `simple_defs`
- `simple_test_utils`
- `simple_wave_field`

Public symbols:
- `run_all_wave_field_tests` — subroutine

Private symbols:
- `random_disk` — function
- `random_field` — function
- `redot` — function
- `redot_disk` — function
- `test_dft_and_round_trip` — subroutine
- `test_disk_transposes` — subroutine
- `test_embedded_amplitude` — subroutine
- `test_freq` — subroutine
- `test_native_grid` — subroutine
- `test_real_transforms` — subroutine
- `test_transform_adjoints` — subroutine

---
## Module: simple_wave_imaging

Files:
- `main/image/simple_wave_imaging.f90`

Uses:
- `simple_ctf`
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_math_ctf`
- `simple_rnd`
- `simple_type_defs`
- `simple_wave_field`

Public symbols:
- `detect_intensity` — subroutine
- `wave_imaging` — type

Private symbols:
- `adjoint` — subroutine
- `adjoint_wave` — subroutine
- `forward` — subroutine
- `forward_wave` — subroutine
- `get_c` — function
- `get_exit_wave` — subroutine
- `get_n` — function
- `get_wl` — function
- `kill` — subroutine
- `lens_adjoint` — subroutine
- `lens_forward` — subroutine
- `lens_linearize` — subroutine
- `linearize` — subroutine
- `linearize_wave` — subroutine
- `new` — subroutine
- `set_optics` — subroutine

---
## Module: simple_wave_imaging_tester

Files:
- `main/image/simple_wave_imaging_tester.f90`

Uses:
- `simple_ctf`
- `simple_defs`
- `simple_image`
- `simple_math_ctf`
- `simple_test_utils`
- `simple_type_defs`
- `simple_wave_imaging`

Public symbols:
- `run_all_wave_imaging_tests` — subroutine

Private symbols:
- `count_disk` — function
- `hermitian` — function
- `misfit2d` — function
- `optics` — function
- `random_disk` — function
- `redot` — function
- `spectrum_scaled` — function
- `test_adjoint_and_gradient` — subroutine
- `test_born_is_rytov_linearised` — subroutine
- `test_detector` — subroutine
- `test_exit_wave` — subroutine
- `test_weak_phase` — subroutine

---
## Module: simple_wave_multislice

Files:
- `main/volume/simple_wave_multislice.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_error`
- `simple_fftw3`
- `simple_type_defs`
- `simple_wave_field`
- `simple_wave_imaging`

Public symbols:
- `wave_multislice` — type

Private symbols:
- `adjoint` — subroutine
- `cube_to_sections` — subroutine
- `forward` — subroutine
- `get_exit_wave` — subroutine
- `get_n` — function
- `kill` — subroutine
- `linearize` — subroutine
- `new` — subroutine
- `propagate` — subroutine
- `section_field` — subroutine
- `sections_to_cube_t` — subroutine
- `set_optics` — subroutine

---
## Module: simple_wave_multislice_tester

Files:
- `main/volume/simple_wave_multislice_tester.f90`

Uses:
- `simple_ctf`
- `simple_defs`
- `simple_math_ctf`
- `simple_test_utils`
- `simple_type_defs`
- `simple_wave_multislice`

Public symbols:
- `run_all_wave_multislice_tests` — subroutine

Private symbols:
- `count_sphere` — function
- `misfit` — function
- `optics` — function
- `random_cube` — function
- `random_disk` — function
- `redot2` — function
- `redot3` — function
- `test_adjoint_and_gradient` — subroutine
- `test_analytic_limits` — subroutine
- `test_weak_sheet` — subroutine

Local symbols (enclosing procedure scope):
- `test_analytic_limits::vacuum_dev` — function
- `test_weak_sheet::peak_of` — function
- `test_weak_sheet::rel_err` — function
- `test_weak_sheet::sheet` — subroutine
- `test_weak_sheet::weak_phase` — function

---
## Module: simple_wave_projector

Files:
- `main/volume/simple_wave_projector.f90`

Uses:
- `iso_c_binding`
- `simple_ctf`
- `simple_defs`
- `simple_error`
- `simple_gridding`
- `simple_image`
- `simple_ori`
- `simple_oris`
- `simple_pcg_lattice`
- `simple_rnd`
- `simple_sym`
- `simple_type_defs`
- `simple_wave_imaging`
- `simple_wave_multislice`

Public symbols:
- `wave_projector` — type

Private symbols:
- `adjoint_planes` — subroutine
- `begin_grad` — subroutine
- `curvature` — function
- `ensure_cubes` — subroutine
- `ensure_imgs` — subroutine
- `ensure_work` — subroutine
- `gain_terms` — subroutine
- `get_exit_wave` — subroutine
- `get_grad` — subroutine
- `get_model` — function
- `get_nptcls` — function
- `kill` — subroutine
- `line_losses` — subroutine
- `linearize_planes` — subroutine
- `loss_grad` — function
- `misfit` — function
- `new` — subroutine
- `pass` — subroutine
- `pass_work` — type
- `plane_to_image` — subroutine
- `predict` — subroutine
- `prep_optics` — subroutine
- `prep_particles` — subroutine
- `set_direction` — subroutine
- `set_sym` — subroutine
- `set_volume` — subroutine
- `shift_vectors` — subroutine
- `spectrum` — subroutine
- `spectrum_cube` — subroutine
- `transform` — subroutine
- `weight` — function

Local symbols (enclosing procedure scope):
- `ensure_cubes::grow` — subroutine

---
## Module: simple_wave_projector_golden_tester

Files:
- `main/volume/simple_wave_projector_golden_tester.f90`

---
## Module: simple_wave_projector_tester

Files:
- `main/volume/simple_wave_projector_tester.f90`

Uses:
- `simple_defs`
- `simple_gridding`
- `simple_image`
- `simple_math_ctf`
- `simple_ori`
- `simple_oris`
- `simple_reconstructor_pcg`
- `simple_sym`
- `simple_test_utils`
- `simple_type_defs`
- `simple_wave_imaging`
- `simple_wave_projector`
- `simple_wave_projector_golden_tester`

Public symbols:
- `run_all_wave_projector_tests` — subroutine

Private symbols:
- `blobs` — function
- `golden_case` — subroutine
- `misfit3d` — function
- `particles` — function
- `random_planes` — function
- `set_particle` — subroutine
- `test_adjoint_and_gradient` — subroutine
- `test_golden` — subroutine
- `test_line_losses` — subroutine
- `test_model_hierarchy` — subroutine
- `test_multislice_orientation` — subroutine
- `test_multislice_weak_limit` — subroutine
- `test_oversampling_one` — subroutine
- `test_pcg_scale` — subroutine
- `test_point_scatterer` — subroutine
- `test_thread_identity` — subroutine
- `test_weak_phase_identity` — subroutine

Local symbols (enclosing procedure scope):
- `test_multislice_orientation::interior` — subroutine

---
## Module: simple_winfuns

Files:
- `main/interp/simple_winfuns.f90`

Uses:
- `simple_defs`

Public symbols:
- `winfuns` — type

Private symbols:
- `bessi0` — function
- `bman_apod` — function
- `bman_instr` — function
- `constructor` — function
- `eval_apod` — function
- `eval_instr` — function
- `get_Whalf` — function
- `hann_apod` — function
- `hann_instr` — function
- `ifun` — function
- `kb_apod` — function
- `kb_instr` — function
- `lin_apod` — function
- `lin_instr` — function
- `nn_apod` — function
- `nn_instr` — function
- `sinc_apod` — function
- `sinc_instr` — function
- `wfun` — function
- `which` — function

---
## Module: simple_xd_gmm

Files:
- `utils/clustering/simple_xd_gmm.f90`

Uses:
- `simple_clustering_utils`
- `simple_core_module_api`
- `simple_linalg`
- `simple_stat`

Public symbols:
- `xd_gmm` — type
- `xd_select_k` — subroutine

Private symbols:
- `component_factor` — subroutine
- `component_moments` — subroutine
- `fit` — subroutine
- `get_covs` — subroutine
- `get_loglik` — function
- `get_means` — subroutine
- `get_pi` — subroutine
- `init` — subroutine
- `kill` — subroutine
- `loglik` — function
- `new` — subroutine
- `order_by_weight` — subroutine
- `outer` — function
- `posterior` — subroutine
- `psd_clip` — subroutine

Local symbols (enclosing procedure scope):
- `xd_select_k::fit_half` — function

---
## Module: simple_xd_gmm_tester

Files:
- `utils/clustering/simple_xd_gmm_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_rnd`
- `simple_test_utils`
- `simple_xd_gmm`

Public symbols:
- `run_all_xd_gmm_tests` — subroutine

Private symbols:
- `fixture` — subroutine
- `test_xd_gmm_recovery` — subroutine
- `test_xd_gmm_select_k` — subroutine

---
## Module: single_commanders_experimental

Files:
- `main/commanders/single/single_commanders_experimental.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_reproject`
- `simple_matcher_2dprep`
- `simple_matcher_ptcl_io`
- `simple_opt_mask`
- `single_commanders_nano3d`

Public symbols:
- `commander_cavgseoproc_nano` — type
- `commander_cavgsproc_nano` — type
- `commander_ptclsproc_nano` — type
- `commander_trajectory_make_projavgs` — type
- `commander_tsegmaps_core_finder` — type
- `commander_validate_cavgs_vs_model` — type
- `exec_cavgseoproc_nano` — subroutine
- `exec_cavgsproc_nano` — subroutine
- `exec_ptclsproc_nano` — subroutine
- `exec_trajectory_make_projavgs` — subroutine
- `exec_tsegmaps_core_finder` — subroutine
- `exec_validate_cavgs_vs_model` — subroutine

---
## Module: single_commanders_nano2D

Files:
- `main/commanders/single/single_commanders_nano2D.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_imgproc`
- `simple_commanders_mkcavgs`
- `simple_commanders_refine2d`
- `simple_commanders_sim`
- `simple_refine2d_fnames`

Public symbols:
- `commander_analysis2D_nano` — type
- `commander_center2D_nano` — type
- `commander_refine2D_nano` — type
- `exec_analysis2D_nano` — subroutine
- `exec_center2D_nano` — subroutine
- `exec_refine2D_nano` — subroutine

---
## Module: single_commanders_nano3D

Files:
- `main/commanders/single/single_commanders_nano3D.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_atoms`
- `simple_commanders_ori`
- `simple_commanders_rec`
- `simple_commanders_refine2d`
- `simple_commanders_refine3d`
- `simple_commanders_reproject`
- `simple_commanders_solve3d`
- `simple_nanoparticle`
- `simple_refine3d_fnames`
- `simple_solve3d_utils`

Public symbols:
- `commander_autorefine3D_nano` — type
- `commander_refine3D_nano` — type
- `commander_solve3D_nano` — type
- `commander_trajectory_reconstruct3D_distr` — type
- `exec_autorefine3D_nano` — subroutine
- `exec_commander_trajectory_reconstruct3D_distr` — subroutine
- `exec_refine3D_nano` — subroutine
- `exec_solve3D_nano` — subroutine

Local symbols (enclosing procedure scope):
- `exec_autorefine3D_nano::copy_atomic_prior_products` — subroutine
- `exec_autorefine3D_nano::copy_model_products` — subroutine
- `exec_autorefine3D_nano::delete_model_products` — subroutine
- `exec_autorefine3D_nano::normalize_detection_products` — subroutine
- `exec_autorefine3D_nano::set_detection_inputs` — subroutine

---
## Module: single_commanders_tomo

Files:
- `main/commanders/single/single_commanders_tomo.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_rec`
- `simple_commanders_refine3d`
- `simple_ori_utils`
- `simple_refine3d_fnames`

Public symbols:
- `commander_import_tilt_series` — type
- `commander_refine_tomo` — type

Private symbols:
- `exec_import_tilt_series` — subroutine
- `exec_refine_tomo` — subroutine
- `read_tilt_angles` — subroutine

Local symbols (enclosing procedure scope):
- `exec_refine_tomo::make_mask` — subroutine
- `exec_refine_tomo::map_res0143` — function
- `exec_refine_tomo::reconstruct` — subroutine
- `exec_refine_tomo::scheduled_lp` — function
- `exec_refine_tomo::write_poses` — subroutine

---
## Module: single_commanders_trajectory

Files:
- `main/commanders/single/single_commanders_trajectory.f90`

Uses:
- `simple_commanders_api`
- `simple_flex_traj_denoise`
- `simple_memoize_ft_maps`
- `simple_tseries_graphene_subtr`
- `single_tseries_tracker`

Public symbols:
- `commander_trajectory_backgr_subtr` — type
- `commander_trajectory_denoise` — type
- `commander_trajectory_extract_substk` — type
- `commander_trajectory_graphene_subtr` — type
- `commander_trajectory_import` — type
- `commander_trajectory_swap_stack` — type
- `commander_tseries_track_particles` — type
- `commander_tseries_track_particles_distr` — type
- `exec_trajectory_backgr_subtr` — subroutine
- `exec_trajectory_denoise` — subroutine
- `exec_trajectory_extract_substk` — subroutine
- `exec_trajectory_graphene_subtr` — subroutine
- `exec_trajectory_import` — subroutine
- `exec_trajectory_swap_stack` — subroutine
- `exec_tseries_track_particles` — subroutine
- `exec_tseries_track_particles_distr` — subroutine

---
## Module: single_commanders_trajectory_solve3D

Files:
- `main/commanders/single/single_commanders_trajectory_solve3D.f90`

Uses:
- `simple_cavg_quality_lattice`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_imgproc`
- `simple_commanders_reproject`
- `simple_commanders_sim`
- `simple_refine3d_fnames`
- `single_commanders_nano2d`
- `single_commanders_nano3d`
- `single_commanders_trajectory`

Public symbols:
- `commander_trajectory_solve3D` — type
- `exec_trajectory_solve3D` — subroutine

Local symbols (enclosing procedure scope):
- `exec_trajectory_solve3D::begin_stage` — function
- `exec_trajectory_solve3D::chdir_checked` — subroutine
- `exec_trajectory_solve3D::end_stage` — subroutine
- `exec_trajectory_solve3D::estimate_diameters` — subroutine
- `exec_trajectory_solve3D::final_map` — function
- `exec_trajectory_solve3D::model_gate` — subroutine
- `exec_trajectory_solve3D::read_records` — subroutine
- `exec_trajectory_solve3D::record_2D` — subroutine
- `exec_trajectory_solve3D::record_frames` — subroutine
- `exec_trajectory_solve3D::record_gate1` — subroutine
- `exec_trajectory_solve3D::record_shifts` — subroutine
- `exec_trajectory_solve3D::require_file` — subroutine
- `exec_trajectory_solve3D::simple_rmdir_checked` — subroutine
- `exec_trajectory_solve3D::simulate_startvol` — subroutine
- `exec_trajectory_solve3D::sort_order` — function
- `exec_trajectory_solve3D::write_gate_pairs` — subroutine
- `exec_trajectory_solve3D::write_pair` — subroutine

---
## Module: single_commanders_tseries

Files:
- `main/commanders/single/single_commanders_tseries.f90`

Uses:
- `simple_bspline_smoother`
- `simple_commanders_api`
- `simple_motion_correct_iter`

Public symbols:
- `commander_tseries_import` — type
- `commander_tseries_make_pickavg` — type
- `commander_tseries_motion_correct` — type
- `commander_tseries_motion_correct_distr` — type
- `exec_tseries_import` — subroutine
- `exec_tseries_make_pickavg` — subroutine
- `exec_tseries_motion_correct` — subroutine
- `exec_tseries_motion_correct_distr` — subroutine

---
## Module: single_exec_api

Files:
- `main/apis/single_exec_api.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_exec_helpers`
- `simple_jiffys`
- `simple_ui`
- `simple_ui_program`
- `single_exec_atom`
- `single_exec_map`
- `single_exec_nano2d`
- `single_exec_nano3d`
- `single_exec_trajectory`
- `single_exec_tseries`
- `single_exec_validate`

---
## Module: single_exec_atom

Files:
- `main/exec/single_exec_atom.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_atoms`
- `simple_commanders_sim`

Public symbols:
- `exec_atom_commander` — subroutine

---
## Module: single_exec_map

Files:
- `main/exec/single_exec_map.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_atoms`
- `single_commanders_experimental`

Public symbols:
- `exec_map_commander` — subroutine

---
## Module: single_exec_nano2D

Files:
- `main/exec/single_exec_nano2D.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_imgproc`
- `single_commanders_nano2d`

Public symbols:
- `exec_nano2D_commander` — subroutine

---
## Module: single_exec_nano3D

Files:
- `main/exec/single_exec_nano3D.f90`

Uses:
- `simple_cmdline`
- `simple_exec_helpers`
- `simple_string`
- `single_commanders_nano3d`
- `single_commanders_tomo`

Public symbols:
- `exec_nano3D_commander` — subroutine

---
## Module: single_exec_trajectory

Files:
- `main/exec/single_exec_trajectory.f90`

Uses:
- `simple_cmdline`
- `single_commanders_experimental`
- `single_commanders_nano3d`
- `single_commanders_trajectory`
- `single_commanders_trajectory_solve3d`

Public symbols:
- `exec_trajectory_commander` — subroutine

---
## Module: single_exec_tseries

Files:
- `main/exec/single_exec_tseries.f90`

Uses:
- `simple_cmdline`
- `single_commanders_trajectory`
- `single_commanders_tseries`

Public symbols:
- `exec_tseries_commander` — subroutine

---
## Module: single_exec_validate

Files:
- `main/exec/single_exec_validate.f90`

Uses:
- `simple_cmdline`
- `single_commanders_experimental`

Public symbols:
- `exec_validate_commander` — subroutine

---
## Module: single_tseries_tracker

Files:
- `main/nano/single_tseries_tracker.f90`

Uses:
- `simple_bspline_smoother`
- `simple_core_module_api`
- `simple_image`
- `simple_image_bin`
- `simple_motion_align_nano`
- `simple_parameters`
- `simple_segmentation`

Public symbols:
- `init_tracker` — subroutine
- `kill_tracker` — subroutine
- `track_particle` — subroutine
- `tracker_get_nnfname` — function
- `tracker_get_nnn` — function
- `tracker_get_stkname` — function

Private symbols:
- `center_reference` — function
- `update_background_pspec` — subroutine
- `write_tester_dump` — subroutine
- `write_trajectory` — subroutine

Local symbols (enclosing procedure scope):
- `center_reference::cleanup_centering_images` — subroutine
- `update_background_pspec::identify_neighbours` — subroutine

---
## Module: single_ui_atom

Files:
- `main/ui/single/single_ui_atom.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_atom_programs` — subroutine
- `new_atoms_register` — subroutine
- `new_atoms_rmsd` — subroutine
- `new_atoms_stats` — subroutine
- `new_core_atoms_analysis` — subroutine
- `new_crys_score` — subroutine
- `new_detect_atomic_species` — subroutine
- `new_detect_atoms` — subroutine
- `new_simulate_nanoparticle` — subroutine

---
## Module: single_ui_map

Files:
- `main/ui/single/single_ui_map.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_map_programs` — subroutine
- `new_conv_atom_denoise` — subroutine
- `new_tsegmaps_core_finder` — subroutine

---
## Module: single_ui_nano2D

Files:
- `main/ui/single/single_ui_nano2D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_nano2D_programs` — subroutine
- `new_analysis2D_nano` — subroutine
- `new_center2D_nano` — subroutine
- `new_estimate_diam` — subroutine
- `new_refine2D_nano` — subroutine

---
## Module: single_ui_nano3D

Files:
- `main/ui/single/single_ui_nano3D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_nano3D_programs` — subroutine
- `new_autorefine3D_nano` — subroutine
- `new_import_tilt_series` — subroutine
- `new_refine3D_nano` — subroutine
- `new_refine_tomo` — subroutine
- `new_solve3D_nano` — subroutine

---
## Module: single_ui_trajectory

Files:
- `main/ui/single/single_ui_trajectory.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_trajectory_programs` — subroutine
- `new_trajectory_denoise` — subroutine
- `new_trajectory_extract_substk` — subroutine
- `new_trajectory_graphene_subtr` — subroutine
- `new_trajectory_import` — subroutine
- `new_trajectory_make_projavgs` — subroutine
- `new_trajectory_reconstruct3D` — subroutine
- `new_trajectory_solve3D` — subroutine
- `new_trajectory_swap_stack` — subroutine

---
## Module: single_ui_tseries

Files:
- `main/ui/single/single_ui_tseries.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_tseries_programs` — subroutine
- `new_tseries_import` — subroutine
- `new_tseries_make_pickavg` — subroutine
- `new_tseries_motion_correct` — subroutine
- `new_tseries_track_particles` — subroutine

---
## Module: single_ui_validate

Files:
- `main/ui/single/single_ui_validate.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_validate_programs` — subroutine
- `new_cavgseoproc_nano` — subroutine
- `new_cavgsproc_nano` — subroutine
- `new_ptclsproc_nano` — subroutine
- `new_validate_cavgs_vs_model` — subroutine

---
## Module: unix

Files:
- `extlibs/unix/src/unix.f90`

Uses:
- `iso_c_binding`
- `unix_dirent`
- `unix_errno`
- `unix_fcntl`
- `unix_inet`
- `unix_ioctl`
- `unix_ipc`
- `unix_mqueue`
- `unix_msg`
- `unix_netdb`
- `unix_pthread`
- `unix_regex`
- `unix_semaphore`
- `unix_signal`
- `unix_socket`
- `unix_stat`
- `unix_stdio`
- `unix_stdlib`
- `unix_string`
- `unix_syslog`
- `unix_termios`
- `unix_time`
- `unix_types`
- `unix_unistd`
- `unix_utsname`
- `unix_wait`

Public symbols:
- `c_f_str_chars` — subroutine
- `c_f_str_ptr` — subroutine
- `c_int32_to_uint16` — function
- `c_int64_to_uint32` — function
- `c_uint16_to_int32` — function
- `c_uint32_to_int64` — function
- `f_c_str_chars` — subroutine
- `f_readdir` — function
- `f_strerror` — function

---
## Module: unix_dirent

Files:
- `extlibs/unix/src/unix_dirent.F90`

Uses:
- `iso_c_binding`

Public symbols:
- `c_closedir` — function
- `c_dirent` — type
- `c_dirent` — type
- `c_opendir` — function
- `c_readdir` — function

---
## Module: unix_errno

Files:
- `extlibs/unix/src/unix_errno.F90`

Uses:
- `iso_c_binding`

Public symbols:
- `c_errno` — function

---
## Module: unix_fcntl

Files:
- `extlibs/unix/src/unix_fcntl.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_fcntl` — function
- `c_open` — function

---
## Module: unix_inet

Files:
- `extlibs/unix/src/unix_inet.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_htonl` — function
- `c_htons` — function
- `c_inet_addr` — function

---
## Module: unix_ioctl

Files:
- `extlibs/unix/src/unix_ioctl.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_ioctl` — function

---
## Module: unix_ipc

Files:
- `extlibs/unix/src/unix_ipc.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_ftok` — function

---
## Module: unix_mqueue

Files:
- `extlibs/unix/src/unix_mqueue.F90`

Uses:
- `iso_c_binding`
- `unix_time`
- `unix_types`

Public symbols:
- `c_mq_attr` — type
- `c_mq_close` — function
- `c_mq_getattr` — function
- `c_mq_open` — function
- `c_mq_receive` — function
- `c_mq_send` — function
- `c_mq_setattr` — function
- `c_mq_timedreceive` — function
- `c_mq_timedsend` — function
- `c_mq_unlink` — function

---
## Module: unix_msg

Files:
- `extlibs/unix/src/unix_msg.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_msgctl` — function
- `c_msgget` — function
- `c_msgrcv` — function
- `c_msgsnd` — function

---
## Module: unix_netdb

Files:
- `extlibs/unix/src/unix_netdb.F90`

Uses:
- `iso_c_binding`
- `unix_fcntl`
- `unix_types`

Public symbols:
- `c_addrinfo` — type
- `c_addrinfo` — type
- `c_freeaddrinfo` — subroutine
- `c_gai_strerror` — function
- `c_getaddrinfo` — function
- `c_in_addr` — type
- `c_sockaddr` — type
- `c_sockaddr` — type
- `c_sockaddr_in` — type
- `c_sockaddr_in` — type

---
## Module: unix_pthread

Files:
- `extlibs/unix/src/unix_pthread.F90`

Uses:
- `iso_c_binding`

Public symbols:
- `c_pthread_cancel` — function
- `c_pthread_create` — function
- `c_pthread_detach` — function
- `c_pthread_exit` — subroutine
- `c_pthread_join` — function
- `c_pthread_mutex_destroy` — function
- `c_pthread_mutex_init` — function
- `c_pthread_mutex_lock` — function
- `c_pthread_mutex_t` — type
- `c_pthread_mutex_trylock` — function
- `c_pthread_mutex_unlock` — function
- `c_pthread_setcancelstate` — function
- `c_pthread_setcanceltype` — function
- `c_pthread_t` — type

---
## Module: unix_regex

Files:
- `extlibs/unix/src/unix_regex.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_regcomp` — function
- `c_regerror` — function
- `c_regex_t` — type
- `c_regex_t` — type
- `c_regexec` — function
- `c_regfree` — subroutine
- `c_regmatch_t` — type

---
## Module: unix_semaphore

Files:
- `extlibs/unix/src/unix_semaphore.F90`

Uses:
- `iso_c_binding`
- `unix_time`
- `unix_types`

Public symbols:
- `c_sem_close` — function
- `c_sem_destroy` — function
- `c_sem_getvalue` — function
- `c_sem_init` — function
- `c_sem_open` — function
- `c_sem_post` — function
- `c_sem_t` — type
- `c_sem_timedwait` — function
- `c_sem_trywait` — function
- `c_sem_unlink` — function
- `c_sem_wait` — function

---
## Module: unix_signal

Files:
- `extlibs/unix/src/unix_signal.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_kill` — function
- `c_signal` — function

---
## Module: unix_socket

Files:
- `extlibs/unix/src/unix_socket.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_accept` — function
- `c_bind` — function
- `c_connect` — function
- `c_listen` — function
- `c_send` — function
- `c_setsockopt` — function
- `c_socket` — function

---
## Module: unix_stat

Files:
- `extlibs/unix/src/unix_stat.F90`

Uses:
- `iso_c_binding`
- `unix_time`
- `unix_types`

Public symbols:
- `c_fstat` — function
- `c_lstat` — function
- `c_mkdir` — function
- `c_mkfifo` — function
- `c_stat` — function
- `c_stat_type` — type
- `c_stat_type` — type
- `c_stat_type` — type
- `c_umask` — function

---
## Module: unix_stdio

Files:
- `extlibs/unix/src/unix_stdio.F90`

Uses:
- `iso_c_binding`

Public symbols:
- `c_fclose` — function
- `c_fdopen` — function
- `c_fflush` — function
- `c_fgetc` — function
- `c_fgets` — function
- `c_fopen` — function
- `c_fprintf` — function
- `c_fputs` — function
- `c_fread` — function
- `c_fwrite` — function
- `c_getchar` — function
- `c_getline` — function
- `c_pclose` — function
- `c_perror` — subroutine
- `c_popen` — function
- `c_putchar` — function
- `c_scanf` — function
- `c_setbuf` — subroutine
- `c_setvbuf` — function

---
## Module: unix_stdlib

Files:
- `extlibs/unix/src/unix_stdlib.F90`

Uses:
- `iso_c_binding`

Public symbols:
- `c_exit` — subroutine
- `c_free` — subroutine

---
## Module: unix_string

Files:
- `extlibs/unix/src/unix_string.F90`

Uses:
- `iso_c_binding`

Public symbols:
- `c_memcpy` — function
- `c_memset` — function
- `c_strerror` — function
- `c_strlen` — function

---
## Module: unix_syslog

Files:
- `extlibs/unix/src/unix_syslog.F90`

Uses:
- `iso_c_binding`

Public symbols:
- `c_closelog` — subroutine
- `c_openlog` — subroutine
- `c_syslog` — subroutine

---
## Module: unix_termios

Files:
- `extlibs/unix/src/unix_termios.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_cfgetispeed` — function
- `c_cfgetospeed` — function
- `c_cfmakeraw` — subroutine
- `c_cfsetispeed` — function
- `c_cfsetospeed` — function
- `c_cfsetspeed` — function
- `c_tcdrain` — function
- `c_tcflow` — function
- `c_tcflush` — function
- `c_tcgetattr` — function
- `c_tcsendbreak` — function
- `c_tcsetattr` — function
- `c_termios` — type
- `c_termios` — type

---
## Module: unix_time

Files:
- `extlibs/unix/src/unix_time.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_asctime` — function
- `c_clock_gettime` — function
- `c_ctime` — function
- `c_gettimeofday` — function
- `c_gmtime` — function
- `c_gmtime_r` — function
- `c_localtime` — function
- `c_localtime_r` — function
- `c_mktime` — function
- `c_strftime` — function
- `c_time` — function
- `c_timegm` — function
- `c_timespec` — type
- `c_timeval` — type
- `c_timezone` — type
- `c_tm` — type

---
## Module: unix_types

Files:
- `extlibs/unix/src/unix_types.F90`

Uses:
- `iso_c_binding`

---
## Module: unix_unistd

Files:
- `extlibs/unix/src/unix_unistd.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_chdir` — function
- `c_close` — function
- `c_dup` — function
- `c_dup2` — function
- `c_execl` — function
- `c_fork` — function
- `c_getpid` — function
- `c_pipe` — function
- `c_read` — function
- `c_setsid` — function
- `c_unlink` — function
- `c_usleep` — function
- `c_write` — function

---
## Module: unix_utsname

Files:
- `extlibs/unix/src/unix_utsname.F90`

Uses:
- `iso_c_binding`

Public symbols:
- `c_uname` — function
- `c_utsname` — type
- `c_utsname` — type

---
## Module: unix_wait

Files:
- `extlibs/unix/src/unix_wait.F90`

Uses:
- `iso_c_binding`
- `unix_types`

Public symbols:
- `c_wait` — function
- `c_waitpid` — function

---
