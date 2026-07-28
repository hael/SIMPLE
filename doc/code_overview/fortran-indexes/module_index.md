# Module Index

## Module: CPlot2D_wrapper_module

Files:
- `extlibs/CPlot2D/simple_CPlot2D_wrapper.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `CDataPoint_type` — type

Private symbols:
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
- `test_CPlot2D` — subroutine

---
## Module: curl

Files:
- `extlibs/curl/src/curl.f90`

---
## Module: curl_easy

Files:
- `extlibs/curl/src/curl_easy.f90`

Public symbols:
- `curl_easy_cleanup` — subroutine
- `curl_easy_init` — function
- `curl_easy_pause` — function
- `curl_easy_perform` — function
- `curl_free` — subroutine
- `curl_global_cleanup` — subroutine
- `curl_mime_addpart` — function
- `curl_mime_data` — function
- `curl_mime_free` — subroutine
- `curl_mime_headers` — function
- `curl_mime_init` — function
- `curl_mime_subparts` — function
- `curl_slist_free_all` — subroutine
- `curl_version_now` — function

Private symbols:
- `curl_easy_escape_` — function
- `curl_easy_getinfo_c_double` — function
- `curl_easy_getinfo_c_long` — function
- `curl_easy_getinfo_c_ptr` — function
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
- `curl_slist` — type
- `curl_slist_append_` — function
- `curl_unescape_` — function
- `curl_version_` — function
- `curl_version_info_` — function
- `curl_version_info_data` — type

---
## Module: curl_multi

Files:
- `extlibs/curl/src/curl_multi.f90`

Public symbols:
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
- `curl_msg` — type
- `curl_multi_strerror_` — function

---
## Module: curl_urlapi

Files:
- `extlibs/curl/src/curl_urlapi.f90`

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

Public symbols:
- `c_f_str_ptr` — subroutine

Private symbols:
- `c_strlen` — function

---
## Module: elemental

Files:
- `main/image/simple_image.f90`
- `main/image/simple_image_access.f90`
- `main/image/simple_image_arith.f90`

---
## Module: FoX_common

Files:
- `extlibs/xml/common/FoX_common.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_count_parse_input`
- `fox_m_fsys_format`
- `fox_m_fsys_parse_input`
- `m_common_attrs`
- `m_common_error`

---
## Module: FoX_dom

Files:
- `extlibs/xml/dom/FoX_dom.f90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_format`
- `m_dom_dom`
- `m_dom_error`
- `m_dom_extras`
- `m_dom_parse`
- `m_dom_utils`

---
## Module: fox_m_fsys_abort_flush

Files:
- `extlibs/xml/fsys/fox_m_fsys_abort_flush.F90`

Uses:
- `f90_unix_io`
- `f90_unix_proc`

Public symbols:
- `abort` — subroutine
- `pure_pxfabort` — function
- `pxfabort` — subroutine
- `pxfflush` — subroutine

---
## Module: fox_m_fsys_array_str

Files:
- `extlibs/xml/fsys/fox_m_fsys_array_str.F90`

---
## Module: fox_m_fsys_count_parse_input

Files:
- `extlibs/xml/fsys/fox_m_fsys_count_parse_input.F90`

Uses:
- `fox_m_fsys_realtypes`

---
## Module: fox_m_fsys_format

Files:
- `extlibs/xml/fsys/fox_m_fsys_format.F90`

Uses:
- `fox_m_fsys_abort_flush`
- `fox_m_fsys_realtypes`

---
## Module: fox_m_fsys_parse_input

Files:
- `extlibs/xml/fsys/fox_m_fsys_parse_input.F90`

Uses:
- `fox_m_fsys_realtypes`

---
## Module: fox_m_fsys_realtypes

Files:
- `extlibs/xml/fsys/fox_m_fsys_realtypes.f90`

---
## Module: fox_m_fsys_string

Files:
- `extlibs/xml/fsys/fox_m_fsys_string.F90`

Public symbols:
- `toLower` — function

---
## Module: fox_m_fsys_string_list

Files:
- `extlibs/xml/fsys/fox_m_fsys_string_list.F90`

Uses:
- `fox_m_fsys_array_str`

---
## Module: fox_m_utils_mtprng

Files:
- `extlibs/xml/utils/fox_m_utils_mtprng.F90`

Private symbols:
- `mtprng_init` — subroutine
- `mtprng_init_by_array` — subroutine
- `mtprng_rand` — function
- `mtprng_rand64` — function
- `mtprng_rand_range` — function
- `mtprng_rand_real1` — function
- `mtprng_rand_real2` — function
- `mtprng_rand_real3` — function

---
## Module: fox_m_utils_uri

Files:
- `extlibs/xml/utils/fox_m_utils_uri.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_format`
- `fox_m_fsys_string`

Public symbols:
- `copyURI` — function
- `destroyURI` — subroutine
- `dumpURI` — subroutine
- `expressURI` — function
- `getAuthority` — function
- `getFragment` — function
- `getHost` — function
- `getPath` — function
- `getPort` — function
- `getQuery` — function
- `getScheme` — function
- `getUserinfo` — function
- `hasAuthority` — function
- `hasFragment` — function
- `hasHost` — function
- `hasPort` — function
- `hasQuery` — function
- `hasScheme` — function
- `hasUserinfo` — function
- `isAbsoluteURI` — function
- `parseURI` — function
- `rebaseURI` — function

Private symbols:
- `appendPaths` — function
- `checkAuthority` — function
- `checkFragment` — function
- `checkHost` — function
- `checkIpvX` — function
- `checkNonOpaquePath` — function
- `checkOpaquePart` — function
- `checkPath` — function
- `checkPathSegment` — function
- `checkQuery` — function
- `checkScheme` — function
- `cleanUp` — subroutine
- `expressSegments` — function
- `expressURI_len` — function
- `normalizepath` — function
- `pctEncode` — function
- `pctEncode_len` — function
- `produceResult` — subroutine
- `unEscape_alloc` — function
- `verifyWithPctEncoding` — function

---
## Module: fox_m_utils_uuid

Files:
- `extlibs/xml/utils/fox_m_utils_uuid.F90`

Uses:
- `fox_m_utils_mtprng`

Public symbols:
- `generate_uuid` — function

Private symbols:
- `get_utc_since_1582` — function
- `int32ToHexOctets` — function
- `int64ToHexOctets` — function
- `isLeapYear` — function

---
## Module: FoX_sax

Files:
- `extlibs/xml/sax/FoX_sax.f90`

Uses:
- `fox_common`
- `m_sax_operate`

---
## Module: FoX_utils

Files:
- `extlibs/xml/utils/FoX_utils.f90`

Uses:
- `fox_m_utils_uri`
- `fox_m_utils_uuid`

---
## Module: FoX_wcml

Files:
- `extlibs/xml/wcml/FoX_wcml.f90`

Uses:
- `fox_wxml`
- `m_wcml_coma`
- `m_wcml_core`
- `m_wcml_geometry`
- `m_wcml_lattice`
- `m_wcml_lists`
- `m_wcml_metadata`
- `m_wcml_molecule`
- `m_wcml_parameter`
- `m_wcml_property`
- `m_wcml_stml`

---
## Module: FoX_wkml

Files:
- `extlibs/xml/wkml/FoX_wkml.f90`

Uses:
- `fox_wxml`
- `m_wkml_chart`
- `m_wkml_color`
- `m_wkml_contours`
- `m_wkml_core`
- `m_wkml_coverage`
- `m_wkml_features`
- `m_wkml_lowlevel`
- `m_wkml_styling`

---
## Module: FoX_wxml

Files:
- `extlibs/xml/wxml/FoX_wxml.f90`

Uses:
- `m_wxml_core`
- `m_wxml_overloads`

---
## Module: function

Files:
- `main/image/simple_image.f90`
- `main/image/simple_image_access.f90`
- `main/image/simple_image_arith.f90`
- `main/image/simple_image_calc.f90`
- `main/image/simple_image_checks.f90`
- `main/image/simple_image_fft.f90`
- `main/image/simple_image_filt.f90`
- `main/image/simple_image_freq_anal.f90`
- `main/image/simple_image_geom.f90`
- `main/image/simple_image_seg.f90`
- `main/nu_filt/simple_nu_filter.f90`
- `main/nu_filt/simple_nu_filter_state.f90`
- `main/ori/simple_oris.f90`
- `main/ori/simple_oris_dists.f90`
- `main/ori/simple_oris_getters.f90`
- `main/ori/simple_oris_life.f90`
- `main/ori/simple_oris_neig.f90`
- `main/ori/simple_oris_sampling.f90`
- `main/ori/simple_oris_stats.f90`
- `main/params/simple_parameters.f90`
- `main/params/simple_parameters_core.f90`
- `main/pftc/simple_polarft_access.f90`
- `main/pftc/simple_polarft_calc.f90`
- `main/project/simple_sp_project.f90`
- `main/project/simple_sp_project_cls.f90`
- `main/project/simple_sp_project_core.f90`
- `main/project/simple_sp_project_io.f90`
- `main/project/simple_sp_project_mic.f90`
- `main/project/simple_sp_project_stk.f90`
- `main/simple_abinitio_controller.f90`
- `main/simple_abinitio_utils.f90`

Uses:
- `simple_online_var`

---
## Module: gnufor2

Files:
- `extlibs/gnufor/gnufor2.f90`

---
## Module: integer

Files:
- `main/nu_filt/simple_nu_filter.f90`
- `main/nu_filt/simple_nu_filter_bank.f90`
- `main/nu_filt/simple_nu_filter_potts.f90`
- `main/nu_filt/simple_nu_filter_state.f90`
- `main/pftc/simple_polarft_access.f90`
- `main/pftc/simple_polarft_calc.f90`
- `main/project/simple_sp_project.f90`
- `main/project/simple_sp_project_core.f90`
- `main/project/simple_sp_project_mic.f90`
- `main/project/simple_sp_project_ptcl.f90`
- `main/project/simple_sp_project_stk.f90`

---
## Module: json_file_module

Files:
- `extlibs/json/json_file_module.f90`

Uses:
- `json_kinds`
- `json_parameters`
- `json_string_utilities`
- `json_value_module`

Private symbols:
- `json_file` — type

---
## Module: json_kinds

Files:
- `extlibs/json/json_kinds.f90`

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

---
## Module: json_value_module

Files:
- `extlibs/json/json_value_module.f90`

Uses:
- `json_kinds`
- `json_parameters`
- `json_string_utilities`

Private symbols:
- `json_core` — type
- `json_value` — type

---
## Module: logical

Files:
- `main/image/simple_image.f90`
- `main/image/simple_image_calc.f90`
- `main/image/simple_image_checks.f90`
- `main/nu_filt/simple_nu_filter.f90`
- `main/nu_filt/simple_nu_filter_potts.f90`
- `main/nu_filt/simple_nu_filter_state.f90`
- `main/pftc/simple_polarft_access.f90`
- `main/pftc/simple_polarft_calc.f90`
- `main/project/simple_sp_project.f90`
- `main/project/simple_sp_project_mic.f90`
- `main/project/simple_sp_project_out.f90`
- `main/project/simple_sp_project_ptcl.f90`

---
## Module: m_common_attrs

Files:
- `extlibs/xml/common/m_common_attrs.F90`

Uses:
- `fox_m_fsys_array_str`
- `m_common_element`
- `m_common_error`

---
## Module: m_common_buffer

Files:
- `extlibs/xml/common/m_common_buffer.F90`

Uses:
- `fox_m_fsys_format`
- `m_common_charset`
- `m_common_error`

---
## Module: m_common_charset

Files:
- `extlibs/xml/common/m_common_charset.F90`

Uses:
- `fox_m_fsys_string`

Public symbols:
- `allowed_encoding` — function
- `checkChars` — function
- `isInitialNameChar` — function
- `isInitialNCNameChar` — function
- `isLegalChar` — function
- `isLegalCharRef` — function
- `isNameChar` — function
- `isNCNameChar` — function
- `isRepCharRef` — function
- `isUSASCII` — function
- `isXML1_0_NameChar` — function
- `isXML1_1_NameChar` — function

---
## Module: m_common_content_model

Files:
- `extlibs/xml/common/m_common_content_model.F90`

Uses:
- `fox_m_fsys_array_str`

Public symbols:
- `checkCP` — function
- `checkCPToEnd` — function
- `destroyCPtree` — subroutine
- `dumpCPtree` — subroutine
- `elementContentCP` — function
- `emptyContentCP` — function
- `newCP` — function
- `transformCPPlus` — subroutine

Private symbols:
- `copyCP` — function
- `copyCPtree` — function
- `destroyCP` — subroutine
- `dumpCP` — subroutine
- `nextCPafterfail` — function
- `nextCPaftermatch` — function
- `nextCPMustMatch` — function

---
## Module: m_common_element

Files:
- `extlibs/xml/common/m_common_element.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_string_list`
- `m_common_charset`
- `m_common_content_model`
- `m_common_error`
- `m_common_namecheck`

---
## Module: m_common_elstack

Files:
- `extlibs/xml/common/m_common_elstack.F90`

Uses:
- `fox_m_fsys_array_str`
- `m_common_content_model`
- `m_common_error`

Private symbols:
- `elstack_item` — type
- `elstack_t` — type

---
## Module: m_common_entities

Files:
- `extlibs/xml/common/m_common_entities.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_format`
- `fox_m_utils_uri`
- `m_common_charset`
- `m_common_error`

---
## Module: m_common_entity_expand

Files:
- `extlibs/xml/common/m_common_entity_expand.F90`

Uses:
- `fox_m_fsys_array_str`
- `m_common_entities`
- `m_common_error`
- `m_common_namecheck`
- `m_common_struct`

Public symbols:
- `expand_entity_value_alloc` — function

---
## Module: m_common_error

Files:
- `extlibs/xml/common/m_common_error.F90`

Uses:
- `fox_m_fsys_abort_flush`
- `fox_m_fsys_array_str`

---
## Module: m_common_io

Files:
- `extlibs/xml/common/m_common_io.F90`

Uses:
- `f90_iostat`
- `m_common_error`

Public symbols:
- `get_unit` — subroutine
- `setup_io` — subroutine

Private symbols:
- `find_eor_eof` — subroutine

---
## Module: m_common_namecheck

Files:
- `extlibs/xml/common/m_common_namecheck.F90`

Uses:
- `fox_m_fsys_format`
- `fox_m_fsys_string`
- `m_common_charset`

Public symbols:
- `checkAttValue` — function
- `checkCharacterEntityReference` — function
- `checkEncName` — function
- `checkName` — function
- `checkNames` — function
- `checkNCName` — function
- `checkNCNames` — function
- `checkNmtoken` — function
- `checkNmtokens` — function
- `checkPEDef` — function
- `checkPITarget` — function
- `checkPseudoAttValue` — function
- `checkPublicId` — function
- `checkQName` — function
- `checkQNames` — function
- `checkRepCharEntityReference` — function
- `likeCharacterEntityReference` — function
- `localpartOfQname` — function
- `prefixOfQName` — function

---
## Module: m_common_namespaces

Files:
- `extlibs/xml/common/m_common_namespaces.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_utils_uri`
- `m_common_attrs`
- `m_common_charset`
- `m_common_error`
- `m_common_namecheck`
- `m_common_struct`

---
## Module: m_common_notations

Files:
- `extlibs/xml/common/m_common_notations.F90`

Uses:
- `fox_m_fsys_array_str`
- `m_common_error`

Public symbols:
- `add_notation` — subroutine
- `destroy_notation_list` — subroutine
- `init_notation_list` — subroutine
- `notation_exists` — function

---
## Module: m_common_struct

Files:
- `extlibs/xml/common/m_common_struct.F90`

Uses:
- `fox_utils`
- `m_common_charset`
- `m_common_element`
- `m_common_entities`
- `m_common_notations`

Public symbols:
- `destroy_xml_doc_state` — subroutine
- `init_xml_doc_state` — subroutine
- `register_external_GE` — subroutine
- `register_external_PE` — subroutine
- `register_internal_GE` — subroutine
- `register_internal_PE` — subroutine

---
## Module: m_contours

Files:
- `extlibs/xml/wkml/m_contours.F90`

Uses:
- `fox_common`
- `m_common_error`

Public symbols:
- `destroy_contourObject` — subroutine
- `make_contours_on_a_complex_grid` — function
- `make_contours_on_a_simple_grid` — function
- `make_contours_on_simplest_grid` — function

Private symbols:
- `add_line` — subroutine
- `add_point` — subroutine
- `add_polygon` — subroutine
- `addLinesToContourObject` — subroutine
- `addLineToLineList` — subroutine
- `addPolygonsToContourObject` — subroutine
- `addPolygonToPolygonList` — subroutine
- `checkPointOnPolygon` — function
- `complexScalePoints` — subroutine
- `concatenate_lines` — subroutine
- `destroy_contourLines` — subroutine
- `draw` — subroutine
- `GCONTR` — subroutine
- `init_contourObject` — subroutine
- `init_line` — subroutine
- `joinUpLines` — subroutine
- `nextLineOnEdge` — function
- `pointInPolygon` — function
- `polygonsInPolygons` — subroutine
- `scale_contours_complex` — subroutine
- `scale_contours_simple` — subroutine
- `setpnew` — subroutine
- `shave_polygon` — subroutine
- `shave_polygons` — subroutine
- `simpleScalePoints` — subroutine
- `whichLevel` — function

---
## Module: m_dom_dom

Files:
- `extlibs/xml/dom/m_dom_dom.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_format`
- `fox_m_fsys_string`
- `fox_m_utils_uri`
- `m_common_charset`
- `m_common_element`
- `m_common_namecheck`
- `m_common_struct`
- `m_dom_error`

---
## Module: m_dom_error

Files:
- `extlibs/xml/dom/m_dom_error.f90`

Uses:
- `fox_m_fsys_abort_flush`
- `m_common_error`

Public symbols:
- `dom_error` — subroutine
- `getExceptionCode` — function
- `inException` — function
- `internal_error` — subroutine
- `throw_exception` — subroutine

Private symbols:
- `destroyDOMException` — subroutine
- `errorString` — function

---
## Module: m_dom_extras

Files:
- `extlibs/xml/dom/m_dom_extras.F90`

Uses:
- `fox_m_fsys_parse_input`
- `fox_m_fsys_realtypes`
- `m_dom_dom`
- `m_dom_error`

---
## Module: m_dom_parse

Files:
- `extlibs/xml/dom/m_dom_parse.f90`

Uses:
- `fox_common`
- `fox_m_fsys_array_str`
- `fox_m_utils_uri`
- `fox_sax`
- `m_common_attrs`
- `m_common_entities`
- `m_common_error`
- `m_common_struct`
- `m_dom_dom`
- `m_dom_error`
- `m_sax_parser`

Public symbols:
- `parsefile` — function
- `parsestring` — function

Private symbols:
- `characters_handler` — subroutine
- `comment_handler` — subroutine
- `endCdata_handler` — subroutine
- `endDocument_Handler` — subroutine
- `endDTD_handler` — subroutine
- `endElement_handler` — subroutine
- `endEntity_handler` — subroutine
- `entityErrorHandler` — subroutine
- `externalEntityDecl_handler` — subroutine
- `FoX_endDTD_handler` — subroutine
- `ignorableWhitespace_handler` — subroutine
- `internalEntityDecl_handler` — subroutine
- `normalErrorHandler` — subroutine
- `notationDecl_handler` — subroutine
- `processingInstruction_handler` — subroutine
- `runParser` — subroutine
- `skippedEntity_handler` — subroutine
- `startCdata_handler` — subroutine
- `startDocument_handler` — subroutine
- `startDTD_handler` — subroutine
- `startElement_handler` — subroutine
- `startEntity_handler` — subroutine
- `unparsedEntityDecl_handler` — subroutine

---
## Module: m_dom_utils

Files:
- `extlibs/xml/dom/m_dom_utils.f90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_format`
- `fox_wxml`
- `m_common_attrs`
- `m_common_element`
- `m_common_struct`
- `m_dom_dom`
- `m_dom_error`

Public symbols:
- `dumpTree` — subroutine
- `serialize` — subroutine

Private symbols:
- `dump2` — subroutine
- `iter_dmp_xml` — subroutine

---
## Module: m_ieee

Files:
- `extlibs/xml/fsys/m_ieee.f90`

Public symbols:
- `generate_nan` — function

---
## Module: m_sax_operate

Files:
- `extlibs/xml/sax/m_sax_operate.F90`

Uses:
- `fox_common`
- `m_common_error`
- `m_sax_parser`
- `m_sax_reader`
- `m_sax_types`

Public symbols:
- `close_xml_t` — subroutine
- `open_xml_file` — subroutine
- `open_xml_string` — subroutine
- `parse` — subroutine
- `stop_parser` — subroutine

Private symbols:
- `attributeDecl_handler` — subroutine
- `characters_handler` — subroutine
- `comment_handler` — subroutine
- `elementDecl_handler` — subroutine
- `endCdata_handler` — subroutine
- `endDocument_handler` — subroutine
- `endDTD_handler` — subroutine
- `endElement_handler` — subroutine
- `endEntity_handler` — subroutine
- `endPrefixMapping_handler` — subroutine
- `error_handler` — subroutine
- `externalEntityDecl_handler` — subroutine
- `fatalError_handler` — subroutine
- `ignorableWhitespace_handler` — subroutine
- `internalEntityDecl_handler` — subroutine
- `notationDecl_handler` — subroutine
- `processingInstruction_handler` — subroutine
- `skippedEntity_handler` — subroutine
- `startCdata_handler` — subroutine
- `startDocument_handler` — subroutine
- `startDTD_handler` — subroutine
- `startElement_handler` — subroutine
- `startEntity_handler` — subroutine
- `startPrefixMapping_handler` — subroutine
- `unparsedEntityDecl_handler` — subroutine
- `warning_handler` — subroutine

---
## Module: m_sax_parser

Files:
- `extlibs/xml/sax/m_sax_parser.F90`

Uses:
- `fox_common`
- `fox_m_fsys_array_str`
- `fox_m_fsys_string_list`
- `fox_utils`
- `m_common_attrs`
- `m_common_charset`
- `m_common_element`
- `m_common_elstack`
- `m_common_entities`
- `m_common_entity_expand`
- `m_common_error`
- `m_common_namecheck`
- `m_common_namespaces`
- `m_common_notations`
- `m_common_struct`
- `m_sax_reader`
- `m_sax_tokenizer`
- `m_sax_types`

Public symbols:
- `getNSDict` — function
- `sax_parse` — subroutine
- `sax_parser_destroy` — subroutine
- `sax_parser_init` — subroutine

Private symbols:
- `add_entity` — subroutine
- `attributeDecl_handler` — subroutine
- `characters_handler` — subroutine
- `checkAttributes` — subroutine
- `checkIdRefs` — subroutine
- `checkXMLAttributes` — subroutine
- `close_tag` — subroutine
- `comment_handler` — subroutine
- `elementDecl_handler` — subroutine
- `endCdata_handler` — subroutine
- `endDocument_handler` — subroutine
- `endDTD_handler` — subroutine
- `endDTDchecks` — subroutine
- `endElement_handler` — subroutine
- `endEntity_handler` — subroutine
- `endPrefixMapping_handler` — subroutine
- `error_handler` — subroutine
- `error_handler` — subroutine
- `externalEntityDecl_handler` — subroutine
- `fatalError_handler` — subroutine
- `FoX_endDTD_handler` — subroutine
- `getDefaultAttributes` — subroutine
- `getLocalNameofQName` — function
- `getURIofQName` — function
- `ignorableWhitespace_handler` — subroutine
- `internalEntityDecl_handler` — subroutine
- `notationDecl_handler` — subroutine
- `open_tag` — subroutine
- `parseDTD` — subroutine
- `processingInstruction_handler` — subroutine
- `sax_error` — subroutine
- `skippedEntity_handler` — subroutine
- `startCdata_handler` — subroutine
- `startDocument_handler` — subroutine
- `startDTD_handler` — subroutine
- `startElement_handler` — subroutine
- `startEntity_handler` — subroutine
- `startPrefixMapping_handler` — subroutine
- `unparsedEntityDecl_handler` — subroutine
- `URIlength` — function
- `warning_handler` — subroutine

---
## Module: m_sax_reader

Files:
- `extlibs/xml/sax/m_sax_reader.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_format`
- `fox_utils`
- `m_common_charset`
- `m_common_error`
- `m_common_io`
- `m_sax_xml_source`

Public symbols:
- `close_file` — subroutine
- `column` — function
- `get_all_characters` — function
- `get_character` — function
- `line` — function
- `open_file` — subroutine
- `open_new_file` — subroutine
- `open_new_string` — subroutine
- `parse_text_declaration` — subroutine
- `parse_xml_declaration` — subroutine
- `pop_buffer_stack` — subroutine
- `push_chars` — subroutine
- `reading_first_entity` — function
- `reading_main_file` — function

Private symbols:
- `close_actual_file` — subroutine
- `open_actual_file` — subroutine

---
## Module: m_sax_tokenizer

Files:
- `extlibs/xml/sax/m_sax_tokenizer.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_utils_uri`
- `m_common_charset`
- `m_common_entities`
- `m_common_error`
- `m_common_namecheck`
- `m_sax_reader`
- `m_sax_types`

Public symbols:
- `expand_pe_text` — function
- `normalize_attribute_text` — function
- `sax_tokenize` — subroutine

Private symbols:
- `tokenizeDTD` — subroutine

---
## Module: m_sax_types

Files:
- `extlibs/xml/sax/m_sax_types.F90`

Uses:
- `m_common_attrs`
- `m_common_elstack`
- `m_common_entities`
- `m_common_error`
- `m_common_namespaces`
- `m_common_notations`
- `m_common_struct`
- `m_sax_reader`

---
## Module: m_sax_xml_source

Files:
- `extlibs/xml/sax/m_sax_xml_source.F90`

Uses:
- `fox_m_fsys_array_str`
- `fox_m_fsys_format`
- `fox_utils`
- `m_common_charset`
- `m_common_error`
- `m_common_io`

Public symbols:
- `get_char_from_file` — function
- `parse_declaration` — subroutine
- `push_file_chars` — subroutine

Private symbols:
- `read_single_char` — function
- `rewind_source` — subroutine

---
## Module: m_wcml_coma

Files:
- `extlibs/xml/wcml/m_wcml_coma.F90`

Uses:
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_wcml_stml`
- `m_wxml_overloads`

---
## Module: m_wcml_core

Files:
- `extlibs/xml/wcml/m_wcml_core.F90`

Uses:
- `fox_common`
- `fox_utils`
- `fox_wxml`
- `m_common_error`
- `m_wcml_metadata`

Public symbols:
- `cmlAddNamespace` — subroutine
- `cmlBeginFile` — subroutine
- `cmlEndCml` — subroutine
- `cmlFinishFile` — subroutine
- `cmlStartCml` — subroutine

---
## Module: m_wcml_geometry

Files:
- `extlibs/xml/wcml/m_wcml_geometry.F90`

Uses:
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_wxml_overloads`

---
## Module: m_wcml_lattice

Files:
- `extlibs/xml/wcml/m_wcml_lattice.F90`

Uses:
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_wxml_overloads`

---
## Module: m_wcml_lists

Files:
- `extlibs/xml/wcml/m_wcml_lists.F90`

Uses:
- `fox_common`
- `fox_wxml`

Public symbols:
- `cmlEndbandList` — subroutine
- `cmlEndkpointList` — subroutine
- `cmlEndmetadataList` — subroutine
- `cmlEndmodule` — subroutine
- `cmlEndparameterList` — subroutine
- `cmlEndpropertyList` — subroutine
- `cmlEndStep` — subroutine
- `cmlStartbandList` — subroutine
- `cmlStartkpointList` — subroutine
- `cmlStartmetadataList` — subroutine
- `cmlStartmodule` — subroutine
- `cmlStartparameterList` — subroutine
- `cmlStartpropertyList` — subroutine
- `cmlStartStep` — subroutine

---
## Module: m_wcml_metadata

Files:
- `extlibs/xml/wcml/m_wcml_metadata.F90`

Uses:
- `fox_wxml`

Public symbols:
- `cmlAddMetadata` — subroutine

---
## Module: m_wcml_molecule

Files:
- `extlibs/xml/wcml/m_wcml_molecule.F90`

Uses:
- `fox_m_fsys_format`
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_common_error`
- `m_wxml_overloads`

---
## Module: m_wcml_parameter

Files:
- `extlibs/xml/wcml/m_wcml_parameter.F90`

Uses:
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_wcml_stml`
- `m_wxml_overloads`

---
## Module: m_wcml_property

Files:
- `extlibs/xml/wcml/m_wcml_property.F90`

Uses:
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_wcml_stml`
- `m_wxml_overloads`

---
## Module: m_wcml_stml

Files:
- `extlibs/xml/wcml/m_wcml_stml.F90`

Uses:
- `fox_wxml`
- `m_wxml_overloads`

---
## Module: m_wkml_chart

Files:
- `extlibs/xml/wkml/m_wkml_chart.F90`

Uses:
- `fox_common`
- `fox_wxml`
- `m_common_error`
- `m_wkml_lowlevel`

---
## Module: m_wkml_color

Files:
- `extlibs/xml/wkml/m_wkml_color.F90`

Uses:
- `fox_common`
- `fox_wxml`
- `m_common_error`
- `m_wkml_color_def`

---
## Module: m_wkml_color_def

Files:
- `extlibs/xml/wkml/m_wkml_color_def.F90`

---
## Module: m_wkml_contours

Files:
- `extlibs/xml/wkml/m_wkml_contours.F90`

Uses:
- `fox_common`
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_common_error`
- `m_contours`
- `m_wkml_color`
- `m_wkml_features`
- `m_wkml_lowlevel`
- `m_wkml_styling`

---
## Module: m_wkml_core

Files:
- `extlibs/xml/wkml/m_wkml_core.F90`

Uses:
- `fox_wxml`
- `m_common_error`
- `m_wkml_lowlevel`

Public symbols:
- `kmlAddNamespace` — subroutine
- `kmlBeginFile` — subroutine
- `kmlFinishFile` — subroutine

---
## Module: m_wkml_coverage

Files:
- `extlibs/xml/wkml/m_wkml_coverage.F90`

Uses:
- `fox_common`
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_common_error`
- `m_wkml_chart`
- `m_wkml_color`
- `m_wkml_features`
- `m_wkml_lowlevel`
- `m_wkml_styling`

---
## Module: m_wkml_features

Files:
- `extlibs/xml/wkml/m_wkml_features.F90`

Uses:
- `fox_common`
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_common_error`
- `m_wkml_chart`
- `m_wkml_color`
- `m_wkml_lowlevel`
- `m_wkml_styling`

---
## Module: m_wkml_lowlevel

Files:
- `extlibs/xml/wkml/m_wkml_lowlevel.F90`

Uses:
- `fox_common`
- `fox_m_fsys_realtypes`
- `fox_utils`
- `fox_wxml`
- `m_common_error`

---
## Module: m_wkml_styling

Files:
- `extlibs/xml/wkml/m_wkml_styling.F90`

Uses:
- `fox_m_fsys_realtypes`
- `fox_wxml`
- `m_common_error`
- `m_wkml_color`
- `m_wkml_lowlevel`

---
## Module: m_wxml_core

Files:
- `extlibs/xml/wxml/m_wxml_core.F90`

Uses:
- `fox_m_fsys_abort_flush`
- `fox_m_fsys_array_str`
- `fox_m_fsys_string`
- `fox_m_utils_uri`
- `m_common_attrs`
- `m_common_buffer`
- `m_common_charset`
- `m_common_element`
- `m_common_elstack`
- `m_common_entities`
- `m_common_error`
- `m_common_io`
- `m_common_namecheck`
- `m_common_namespaces`
- `m_common_notations`
- `m_common_struct`
- `m_wxml_escape`

---
## Module: m_wxml_escape

Files:
- `extlibs/xml/wxml/m_wxml_escape.F90`

Uses:
- `fox_m_fsys_format`
- `m_common_charset`
- `m_common_error`

Public symbols:
- `escape_string` — function
- `escape_string_len` — function

---
## Module: m_wxml_overloads

Files:
- `extlibs/xml/wxml/m_wxml_overloads.F90`

Uses:
- `fox_m_fsys_format`
- `fox_m_fsys_realtypes`
- `m_wxml_core`

---
## Module: mod_phasecorr

Files:
- `main/commanders/test/simple_test_mod_phasecorr.f90`

Uses:
- `simple_image`

Public symbols:
- `alloc` — subroutine
- `gen_images` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `run` — subroutine
- `t_phasecorr` — type
- `vis_imgs` — subroutine

---
## Module: procedure

Files:
- `defs/simple_type_defs.f90`
- `extlibs/curl/src/curl_easy.f90`
- `extlibs/gnufor/gnufor2.f90`
- `extlibs/json/json_file_module.f90`
- `extlibs/json/json_string_utilities.f90`
- `extlibs/json/json_value_module.f90`
- `extlibs/unix/src/unix.f90`
- `extlibs/xml/common/m_common_attrs.F90`
- `extlibs/xml/common/m_common_buffer.F90`
- `extlibs/xml/common/m_common_element.F90`
- `extlibs/xml/common/m_common_elstack.F90`
- `extlibs/xml/common/m_common_entities.F90`
- `extlibs/xml/common/m_common_error.F90`
- `extlibs/xml/common/m_common_namespaces.F90`
- `extlibs/xml/dom/m_dom_dom.F90`
- `extlibs/xml/dom/m_dom_extras.F90`
- `extlibs/xml/fsys/fox_m_fsys_array_str.F90`
- `extlibs/xml/fsys/fox_m_fsys_count_parse_input.F90`
- `extlibs/xml/fsys/fox_m_fsys_format.F90`
- `extlibs/xml/fsys/fox_m_fsys_parse_input.F90`
- `extlibs/xml/fsys/fox_m_fsys_string_list.F90`
- `extlibs/xml/wcml/m_wcml_coma.F90`
- `extlibs/xml/wcml/m_wcml_geometry.F90`
- `extlibs/xml/wcml/m_wcml_lattice.F90`
- `extlibs/xml/wcml/m_wcml_molecule.F90`
- `extlibs/xml/wcml/m_wcml_parameter.F90`
- `extlibs/xml/wcml/m_wcml_property.F90`
- `extlibs/xml/wcml/m_wcml_stml.F90`
- `extlibs/xml/wkml/m_wkml_chart.F90`
- `extlibs/xml/wkml/m_wkml_color.F90`
- `extlibs/xml/wkml/m_wkml_contours.F90`
- `extlibs/xml/wkml/m_wkml_coverage.F90`
- `extlibs/xml/wkml/m_wkml_features.F90`
- `extlibs/xml/wkml/m_wkml_lowlevel.F90`
- `extlibs/xml/wkml/m_wkml_styling.F90`
- `extlibs/xml/wxml/m_wxml_core.F90`
- `extlibs/xml/wxml/m_wxml_overloads.F90`
- `fileio/simple_fileio.f90`
- `main/ctf/simple_ctf.f90`
- `main/image/simple_ftiter.f90`
- `main/image/simple_image.f90`
- `main/interp/simple_edges_sqwins.f90`
- `main/interp/simple_kbinterpol.f90`
- `main/interp/simple_winfuns.f90`
- `main/ori/simple_ori.f90`
- `main/ori/simple_oris.f90`
- `main/simple_abinitio_controller.f90`
- `main/simple_eul_prob_tab_utils.f90`
- `main/simple_sym.f90`
- `main/strategies/search/simple_matcher_ptcl_io.f90`
- `main/stream/simple_stream_watcher.f90`
- `utils/filter/simple_butterworth.f90`
- `utils/filter/simple_opt_filter.f90`
- `utils/math/simple_linalg.f90`
- `utils/math/simple_math.f90`
- `utils/math/simple_math_ft.f90`
- `utils/math/simple_neighs.f90`
- `utils/math/simple_ran_tabu.f90`
- `utils/math/simple_rnd.f90`
- `utils/math/simple_srch_sort_loc.f90`
- `utils/math/simple_stat.f90`
- `utils/qsys/simple_qsys_ctrl.f90`
- `utils/qsys/simple_qsys_funs.f90`
- `utils/simple_corrmat.f90`
- `utils/simple_gauss2Dfit.f90`
- `utils/simple_gpu_utils.f90`
- `utils/simple_imgarr_utils.f90`
- `utils/simple_is_check_assert.f90`
- `utils/simple_jiffys.f90`
- `utils/simple_magic_boxes.f90`
- `utils/simple_sauron.f90`
- `utils/simple_segmentation.f90`
- `utils/simple_string.f90`
- `utils/simple_string_utils.f90`
- `utils/structs/simple_chash.f90`
- `utils/structs/simple_hash.f90`

Uses:
- `fox_common`
- `fox_wxml`
- `iso_fortran_env`
- `m_wkml_color`
- `simple_cmdline`
- `simple_image`
- `simple_math`
- `simple_qsys_local`
- `simple_rnd`
- `simple_sauron`
- `simple_string`
- `simple_syslib`
- `simple_type_defs`

Public symbols:
- `add_attribute` — function
- `add_element` — function
- `add_entity` — subroutine
- `add_error` — subroutine
- `add_external_entity` — subroutine
- `add_internal_entity` — subroutine
- `add_item_to_dict` — subroutine
- `add_to_buffer` — subroutine
- `add_to_path` — subroutine
- `addDefaultNS` — subroutine
- `addPrefix` — subroutine
- `addPrefixedNS` — subroutine
- `addPrefixedURI` — subroutine
- `adoptNode` — function
- `alloc_str` — type
- `annotate_invalid_json` — subroutine
- `append` — function
- `append_text` — subroutine
- `appendChild` — function
- `appendData` — subroutine
- `att_value_normalize` — function
- `att_value_normalize_len` — function
- `attribute_has_default` — function
- `attributeDecl_handler` — subroutine
- `buffer_length` — function
- `buffer_to_chararray` — function
- `buffer_to_str` — function
- `c_f_str_chars` — subroutine
- `c_f_str_ptr` — subroutine
- `c_int32_to_uint16` — function
- `c_int64_to_uint32` — function
- `c_uint16_to_int32` — function
- `c_uint32_to_int64` — function
- `check_buffer` — subroutine
- `check_duplicates` — function
- `check_if_valid` — subroutine
- `checkContentModel` — function
- `checkContentModelToEnd` — function
- `checkEndNamespaces` — subroutine
- `checkNamespaces` — subroutine
- `checkNamespacesWriting` — subroutine
- `cloneNode` — function
- `cmlEndMolecule` — subroutine
- `cmlStartMolecule` — subroutine
- `compact_real_string` — subroutine
- `convert` — subroutine
- `copyDOMConfig` — subroutine
- `copyURIMapping` — subroutine
- `createAttribute` — function
- `createAttributeNS` — function
- `createCdataSection` — function
- `createComment` — function
- `createDocument` — function
- `createDocumentFragment` — function
- `createDocumentType` — function
- `createElement` — function
- `createElementNS` — function
- `createEmptyDocument` — function
- `createEmptyEntityReference` — function
- `createEntity` — function
- `createEntityReference` — function
- `createNamespaceNode` — function
- `createNotation` — function
- `createProcessingInstruction` — function
- `createTextNode` — function
- `curl_easy_escape` — function
- `curl_easy_getinfo_char` — function
- `curl_easy_getinfo_double` — function
- `curl_easy_getinfo_int` — function
- `curl_easy_getinfo_long` — function
- `curl_easy_getinfo_ptr` — function
- `curl_easy_setopt_char` — function
- `curl_easy_setopt_funptr` — function
- `curl_easy_setopt_int` — function
- `curl_easy_setopt_long` — function
- `curl_easy_setopt_ptr` — function
- `curl_easy_strerror` — function
- `curl_easy_unescape` — function
- `curl_escape` — function
- `curl_global_init` — function
- `curl_mime_encoder` — function
- `curl_mime_filedata` — function
- `curl_mime_filename` — function
- `curl_mime_name` — function
- `curl_mime_type` — function
- `curl_slist_append` — function
- `curl_unescape` — function
- `curl_version` — function
- `curl_version_info` — function
- `declared_element` — function
- `decode_rfc6901` — function
- `default_comp_ucs4` — function
- `default_join_ucs4` — function
- `default_neq_ucs4` — function
- `deleteData` — subroutine
- `destroy_attribute_list` — subroutine
- `destroy_attribute_t` — subroutine
- `destroy_dict` — subroutine
- `destroy_dict_item` — subroutine
- `destroy_element_list` — subroutine
- `destroy_elstack` — subroutine
- `destroy_entity` — subroutine
- `destroy_entity_list` — subroutine
- `destroy_error_stack` — subroutine
- `destroy_json_core` — subroutine
- `destroy_json_data` — subroutine
- `destroyAllNodesRecursively` — subroutine
- `destroyDocument` — subroutine
- `destroyNamedNodeMap` — subroutine
- `destroyNamespaceDictionary` — subroutine
- `destroyNodeList` — subroutine
- `dump_buffer` — subroutine
- `dumpnsdict` — subroutine
- `duplicate_key_func` — subroutine
- `elementContent` — function
- `emptyContent` — function
- `encode_rfc6901` — function
- `end_prefix_handler` — subroutine
- `end_prefix_handler` — subroutine
- `ensure_capacity` — subroutine
- `escape_string` — subroutine
- `existing_attribute` — function
- `existing_element` — function
- `existing_entity` — function
- `expand_char_entity` — function
- `expand_char_entity_len` — function
- `expand_entity` — function
- `expand_entity_len` — function
- `expand_entity_text` — function
- `expand_entity_text_len` — function
- `express_att_decl_len` — function
- `express_attribute_declaration` — function
- `f_c_str_chars` — subroutine
- `f_readdir` — function
- `f_strerror` — function
- `FoX_error` — subroutine
- `FoX_error_base` — subroutine
- `FoX_fatal_base` — subroutine
- `FoX_get_fatal_errors` — function
- `FoX_get_fatal_warnings` — function
- `FoX_set_fatal_errors` — subroutine
- `FoX_set_fatal_warnings` — subroutine
- `FoX_warning_base` — subroutine
- `get_att_index_pointer` — subroutine
- `get_att_type_enum` — function
- `get_attdecl_by_index` — function
- `get_attdecl_by_name` — function
- `get_attlist_size` — function
- `get_attribute` — function
- `get_chars_from_array` — subroutine
- `get_chars_from_array` — subroutine
- `get_current_line_from_file_sequential` — subroutine
- `get_current_line_from_file_stream` — subroutine
- `get_double_from_array` — subroutine
- `get_element` — function
- `get_int_from_array` — subroutine
- `get_json_core_in_file` — subroutine
- `get_key` — function
- `get_key_index` — function
- `get_key_index_ns` — function
- `get_key_len` — function
- `get_localName_by_index` — function
- `get_localName_by_keyname` — function
- `get_localname_by_keyname_len` — function
- `get_logical_from_array` — subroutine
- `get_nsURI_by_index` — function
- `get_nsURI_by_keyname` — function
- `get_nsURI_by_keyname_len` — function
- `get_prefix_by_index` — function
- `get_prefix_by_keyname` — function
- `get_prefix_by_keyname_len` — function
- `get_string_lengths` — subroutine
- `get_top_elstack` — function
- `get_value_by_index` — function
- `get_value_by_index_len` — function
- `get_value_by_key` — function
- `get_value_by_key_len` — function
- `get_value_by_key_ns` — function
- `get_value_by_key_ns_len` — function
- `getAttribute` — function
- `getAttributeNode` — function
- `getAttributeNodeNS` — function
- `getAttributeNS` — function
- `getAttributes` — function
- `getBase` — function
- `getBase_len` — function
- `getBaseURI` — function
- `getchildNodes` — function
- `getdata` — function
- `getdocType` — function
- `getdocumentElement` — function
- `getdocumentURI` — function
- `getdomConfig` — function
- `getElementById` — function
- `getElementsByTagName` — function
- `getElementsByTagNameNS` — function
- `getEntities` — function
- `getEntityByIndex` — function
- `getEntityByName` — function
- `getEntityNameByIndex` — function
- `getfirstChild` — function
- `getFoX_checks` — function
- `getillFormed` — function
- `getImplementation` — function
- `getInputEncoding` — function
- `getInternalSubset` — function
- `getisElementContentWhitespace` — function
- `getIsId_by_index` — function
- `getlastChild` — function
- `getLength` — function
- `getliveNodeLists` — function
- `getLocalName` — function
- `getname` — function
- `getNamedItem` — function
- `getNamedItemNS` — function
- `getnamespaceNodes` — function
- `getNamespaceURI` — function
- `getnextSibling` — function
- `getnodeName` — function
- `getNodePath` — function
- `getnodeType` — function
- `getNodeValue` — function
- `getnotationName` — function
- `getNotations` — function
- `getNumberOfPrefixes` — function
- `getOwnerDocument` — function
- `getownerElement` — function
- `getParameter` — function
- `getParameterNames` — function
- `getparentNode` — function
- `getPrefix` — function
- `getPrefixByIndex` — function
- `getPrefixIndex` — function
- `getpreviousSibling` — function
- `getpublicId` — function
- `getreadonly` — function
- `getspecified` — function
- `getstrictErrorChecking` — function
- `getstringValue` — function
- `getsystemId` — function
- `gettagName` — function
- `getTarget` — function
- `getTextContent` — function
- `getType_by_index` — function
- `getType_by_index_len` — function
- `getType_by_keyname` — function
- `getType_by_keyname_len` — function
- `getURIofDefaultNS` — function
- `getURIofPrefixedNS` — function
- `getWhitespaceHandling` — function
- `getXds` — function
- `getXmlEncoding` — function
- `getxmlStandalone` — function
- `getXmlVersion` — function
- `has_key` — function
- `has_key_ns` — function
- `hasAttribute` — function
- `hasAttributeNS` — function
- `hasAttributes` — function
- `hasChildNodes` — function
- `hasFeature` — function
- `hist` — subroutine
- `image_1` — subroutine
- `image_2` — subroutine
- `image_3` — subroutine
- `image_4` — subroutine
- `image_5` — subroutine
- `image_6` — subroutine
- `importNode` — function
- `in_error` — function
- `init_attribute_list` — subroutine
- `init_dict` — subroutine
- `init_element_list` — subroutine
- `init_elstack` — subroutine
- `init_entity_list` — subroutine
- `init_error_stack` — subroutine
- `initialize_json_core` — function
- `initialize_json_core_in_file` — subroutine
- `initialize_json_file` — function
- `initialize_json_file_v2` — function
- `initNamespaceDictionary` — subroutine
- `insertBefore` — function
- `insertData` — subroutine
- `integer_to_string` — subroutine
- `is_child_of_callback` — subroutine
- `is_empty_elstack` — function
- `is_external_entity` — function
- `is_unparsed_entity_` — function
- `is_unparsed_entity_from_list` — function
- `isDeclared_by_index` — function
- `isDeclared_by_key` — function
- `isDeclared_by_keyNS` — function
- `isDefaultNamespace` — function
- `isDefaultNSInForce` — function
- `isEqualNode` — function
- `isPrefixInForce` — function
- `isSameNode` — function
- `isSpecified_by_index` — function
- `isSpecified_by_key` — function
- `isSpecified_by_keyNS` — function
- `isSupported` — function
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
- `json_array_callback_func` — subroutine
- `json_check_all_for_duplicate_keys` — subroutine
- `json_check_children_for_duplicate_keys` — subroutine
- `json_check_for_errors` — subroutine
- `json_clear_exceptions` — subroutine
- `json_clone` — subroutine
- `json_count` — function
- `json_create_by_path` — subroutine
- `json_failed` — function
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
- `json_traverse_callback_func` — subroutine
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
- `kmlAddLegend` — subroutine
- `kmlCloseInnerBoundaryIs` — subroutine
- `kmlCloseLinearRing` — subroutine
- `kmlCloseLineString` — subroutine
- `kmlCloseouterBoundaryIs` — subroutine
- `kmlClosePolygon` — subroutine
- `kmlCloseStyle` — subroutine
- `kmlCreateLineStyle` — subroutine
- `kmlCreatePolygonStyle` — subroutine
- `kmlEndRegion` — subroutine
- `kmlGetColorHex` — function
- `kmlMakeColorMap` — function
- `kmlOpenInnerBoundaryIs` — subroutine
- `kmlOpenLinearRing` — subroutine
- `kmlOpenLineString` — subroutine
- `kmlOpenOuterBoundaryIs` — subroutine
- `kmlOpenPolygon` — subroutine
- `kmlOpenStyle` — subroutine
- `kmlSetCustomColor` — subroutine
- `lookupNamespaceURI` — function
- `lookupPrefix` — function
- `lowercase_string` — function
- `make_token_group` — function
- `make_token_group_len` — function
- `my_date_and_time` — function
- `name_equal` — function
- `name_strings_equal` — function
- `namespaceFixup` — subroutine
- `newDOMConfig` — function
- `normalize` — subroutine
- `normalizeDocument` — subroutine
- `number_of_items` — function
- `output_terminal` — function
- `parse_array` — subroutine
- `parse_dtd_attlist` — subroutine
- `parse_dtd_element` — subroutine
- `parse_for_chars` — subroutine
- `parse_number` — subroutine
- `parse_object` — subroutine
- `parse_string` — subroutine
- `parse_value` — subroutine
- `plot3d` — subroutine
- `plot_1` — subroutine
- `plot_2` — subroutine
- `plot_3` — subroutine
- `plot_4` — subroutine
- `pop_char` — subroutine
- `pop_elstack` — function
- `pop_entity_list` — function
- `pop_nl` — function
- `print_buffer` — subroutine
- `print_dict` — subroutine
- `print_elstack` — subroutine
- `print_entity_list` — subroutine
- `push_char` — subroutine
- `push_elstack` — subroutine
- `real_to_string` — subroutine
- `remove_key_by_index` — subroutine
- `remove_nl` — function
- `removeAttribute` — subroutine
- `removeAttributeNode` — function
- `removeAttributeNodeNS` — function
- `removeAttributeNS` — subroutine
- `removeChild` — function
- `removeDefaultNS` — subroutine
- `removeNamedItem` — function
- `removeNamedItemNS` — function
- `removePrefix` — subroutine
- `removePrefixedNS` — subroutine
- `removePrefixedURI` — subroutine
- `renameNode` — function
- `replace_string` — subroutine
- `replaceChild` — function
- `replaceData` — subroutine
- `report_declarations` — subroutine
- `reset_buffer` — subroutine
- `reset_dict` — subroutine
- `reset_elstack` — subroutine
- `reset_entity_list` — subroutine
- `resize_elstack` — subroutine
- `run_gnuplot` — subroutine
- `set_json_core_in_file` — subroutine
- `set_localName_by_index_s` — subroutine
- `set_localName_by_index_vs` — subroutine
- `set_nsURI_by_index` — subroutine
- `set_prefix_by_index` — subroutine
- `setAttribute` — subroutine
- `setAttributeNode` — function
- `setAttributeNodeNS` — function
- `setAttributeNS` — subroutine
- `setBase` — subroutine
- `setData` — subroutine
- `setDeclared` — subroutine
- `setDocType` — subroutine
- `setDocumentElement` — subroutine
- `setdocumentURI` — subroutine
- `setdomConfig` — subroutine
- `setFoX_checks` — subroutine
- `setGCstate` — subroutine
- `setIdAttribute` — subroutine
- `setIdAttributeNode` — subroutine
- `setIdAttributeNS` — subroutine
- `setillFormed` — subroutine
- `setIsElementContentWhitespace` — subroutine
- `setIsId_by_index` — subroutine
- `setliveNodeLists` — subroutine
- `setNamedItem` — function
- `setNamedItemNS` — function
- `setNodeValue` — subroutine
- `setParameter` — subroutine
- `setPrefix` — subroutine
- `setReadOnlyMap` — subroutine
- `setReadOnlyNode` — subroutine
- `setSpecified` — subroutine
- `setspecified` — subroutine
- `setstrictErrorChecking` — subroutine
- `setstringValue` — subroutine
- `setTextContent` — subroutine
- `setValue` — subroutine
- `setXds` — subroutine
- `setxmlStandalone` — subroutine
- `setXmlVersion` — subroutine
- `shallow_copy_entity` — function
- `sigma_array` — type
- `size_el` — function
- `sortAttrs` — subroutine
- `splitText` — function
- `start_prefix_handler` — subroutine
- `str_to_int_10` — function
- `str_to_int_16` — function
- `str_vs` — function
- `string_to_dble` — function
- `string_to_int` — function
- `string_to_integer` — subroutine
- `string_to_real` — subroutine
- `strip_spaces` — function
- `subStringData` — function
- `surf_1` — subroutine
- `surf_2` — subroutine
- `surf_3` — subroutine
- `swap_pointers` — subroutine
- `to_array` — subroutine
- `to_double` — subroutine
- `to_integer` — subroutine
- `to_logical` — subroutine
- `to_null` — subroutine
- `to_object` — subroutine
- `to_string` — subroutine
- `to_uni` — function
- `to_uni_vec` — function
- `traverse` — subroutine
- `ucs4_comp_default` — function
- `ucs4_join_default` — function
- `ucs4_neq_default` — function
- `unescape_string` — subroutine
- `valid_json_hex` — function
- `vs_str` — function
- `vs_str_alloc` — function
- `vs_vs_alloc` — function
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
- `write_it` — subroutine
- `write_it_fast` — subroutine

Private symbols:
- `add2fbody_1` — function
- `add2fbody_2` — function
- `add2fbody_3` — function
- `add2history_1` — subroutine
- `add2history_2` — subroutine
- `add2watchdirs` — subroutine
- `add_eol` — subroutine
- `add_string` — subroutine
- `add_to_stream_stack` — subroutine
- `add_to_streaming` — subroutine
- `addBondArray` — subroutine
- `addcoords_x3_dp` — subroutine
- `addcoords_x3_sp` — subroutine
- `addcoords_xfrac_dp` — subroutine
- `addcoords_xfrac_sp` — subroutine
- `addcoords_xyz3_dp` — subroutine
- `addcoords_xyz3_sp` — subroutine
- `addcoords_xyzfrac_dp` — subroutine
- `addcoords_xyzfrac_sp` — subroutine
- `addDlpolyMatrix_3_dp` — subroutine
- `addDlpolyMatrix_3_sp` — subroutine
- `addDlpolyMatrix_dp` — subroutine
- `addDlpolyMatrix_sp` — subroutine
- `alloc_chash` — subroutine
- `alloc_hash` — subroutine
- `alloc_imgarr` — subroutine
- `analyze_smat` — subroutine
- `ang2vox` — function
- `angle_sampling_1` — function
- `angle_sampling_2` — function
- `apod` — function
- `apod_fast` — function
- `apod_fast_device` — function
- `apod_kb15_a2` — function
- `apod_mat_2d` — subroutine
- `apod_mat_2d_fast` — subroutine
- `apod_mat_3d` — subroutine
- `apod_mat_3d_fast` — subroutine
- `append2basename_1` — function
- `append2basename_2` — function
- `append_candidate` — subroutine
- `append_limited_char` — subroutine
- `append_nl` — subroutine
- `append_nnm` — subroutine
- `append_or_replace_candidate` — subroutine
- `append_ori` — subroutine
- `appendNSNode` — subroutine
- `apply2all` — subroutine
- `apply_1` — subroutine
- `apply_2` — subroutine
- `apply_convention` — subroutine
- `apply_refine3D_search_overrides` — subroutine
- `apply_sym_with_shift` — subroutine
- `arg` — function
- `arpack_stop` — subroutine
- `arr2file_dp` — subroutine
- `arr2file_sp` — subroutine
- `arr2txtfile_1` — subroutine
- `arr2txtfile_2` — subroutine
- `arraytocomplexdp` — subroutine
- `arraytocomplexsp` — subroutine
- `arraytointeger` — subroutine
- `arraytological` — subroutine
- `arraytorealdp` — subroutine
- `arraytorealsp` — subroutine
- `arraytostring` — subroutine
- `assert_eq2` — function
- `assert_eq3` — function
- `assert_eq4` — function
- `assert_eqn` — function
- `assign` — subroutine
- `AttributeArrayCh` — subroutine
- `AttributeArrayCmplxDp` — subroutine
- `AttributeArrayCmplxSp` — subroutine
- `AttributeArrayInt` — subroutine
- `AttributeArrayLg` — subroutine
- `AttributeArrayRealDp` — subroutine
- `AttributeArrayRealSp` — subroutine
- `AttributeMatrixCh` — subroutine
- `AttributeMatrixCmplxDp` — subroutine
- `AttributeMatrixCmplxSp` — subroutine
- `AttributeMatrixInt` — subroutine
- `AttributeMatrixLg` — subroutine
- `AttributeMatrixRealDp` — subroutine
- `AttributeMatrixRealSp` — subroutine
- `AttributeScalarCmplxDp` — subroutine
- `AttributeScalarCmplxSp` — subroutine
- `AttributeScalarInt` — subroutine
- `AttributeScalarLg` — subroutine
- `AttributeScalarRealDp` — subroutine
- `AttributeScalarRealSp` — subroutine
- `augment_partition_job_descr` — subroutine
- `automatic_thresh_sobel` — subroutine
- `autoscale` — subroutine
- `avg_frac_smallest` — function
- `avg_sdev_1` — subroutine
- `avg_sdev_2` — subroutine
- `avg_sdev_3` — subroutine
- `avg_sdev_4` — subroutine
- `balanced` — subroutine
- `basename` — function
- `batch_gauss2Dfit_1` — subroutine
- `batch_gauss2Dfit_2` — subroutine
- `bessi0` — function
- `bman_apod` — function
- `bman_instr` — function
- `bounds_from_mask3D` — subroutine
- `build_eullims` — subroutine
- `build_pind_lookup` — subroutine
- `build_refine3D_stage_cfg` — subroutine
- `build_refspiral` — subroutine
- `butterworth` — function
- `butterworth_filter_1` — subroutine
- `butterworth_filter_2` — subroutine
- `butterworth_filter_3` — subroutine
- `butterworth_filter_4` — subroutine
- `butterworth_filter_5` — subroutine
- `calc_ap_pref` — function
- `calc_athres` — function
- `calc_cartesian_corrmat_1` — subroutine
- `calc_cartesian_corrmat_2` — subroutine
- `calc_graphene_mask` — function
- `calc_inpl_invariant_cc_nomirr` — function
- `calc_num2sample` — subroutine
- `calc_offset2D` — subroutine
- `calc_score_thres` — function
- `calc_stats` — subroutine
- `canny` — subroutine
- `canny_edge` — subroutine
- `canSetParameter_ch` — function
- `canSetParameter_log` — function
- `cast_str_types` — function
- `char2str` — function
- `CharactersArrayCh` — subroutine
- `CharactersArrayCmplxDp` — subroutine
- `CharactersArrayCmplxSp` — subroutine
- `CharactersArrayInt` — subroutine
- `CharactersArrayLg` — subroutine
- `CharactersArrayRealDp` — subroutine
- `CharactersArrayRealSp` — subroutine
- `CharactersMatrixCh` — subroutine
- `CharactersMatrixCmplxDp` — subroutine
- `CharactersMatrixCmplxSp` — subroutine
- `CharactersMatrixInt` — subroutine
- `CharactersMatrixLg` — subroutine
- `CharactersMatrixRealDp` — subroutine
- `CharactersMatrixRealSp` — subroutine
- `CharactersScalarCmplxDp` — subroutine
- `CharactersScalarCmplxSp` — subroutine
- `CharactersScalarInt` — subroutine
- `CharactersScalarLg` — subroutine
- `CharactersScalarRealDp` — subroutine
- `CharactersScalarRealSp` — subroutine
- `chash2ori` — subroutine
- `chash2str` — function
- `check4nans2D_1` — subroutine
- `check4nans2D_2` — subroutine
- `check4nans3D_1` — subroutine
- `check4nans3D_2` — subroutine
- `check4nans_1` — subroutine
- `check4nans_2` — subroutine
- `check_xf` — subroutine
- `checkBondIdRefs` — subroutine
- `checkColorHex` — function
- `checkExistingRefs` — function
- `checkExistingRefsInAttValue` — function
- `checkFmt` — function
- `checkParsedRefsInAttValue` — function
- `clear_history` — subroutine
- `clear_partition_job_descr` — subroutine
- `clear_stack` — subroutine
- `close_start_tag` — subroutine
- `cmlAddAngle_dp` — subroutine
- `cmlAddAngle_sp` — subroutine
- `cmlAddAtoms_3_dp` — subroutine
- `cmlAddAtoms_3_dp_sh` — subroutine
- `cmlAddAtoms_3_sp` — subroutine
- `cmlAddAtoms_3_sp_sh` — subroutine
- `cmlAddAtomsdp` — subroutine
- `cmlAddAtomsdp_sh` — subroutine
- `cmlAddAtomssp` — subroutine
- `cmlAddAtomssp_sh` — subroutine
- `cmlAddBandListdp` — subroutine
- `cmlAddBandListsp` — subroutine
- `cmlAddCoords_dp` — subroutine
- `cmlAddCoords_sp` — subroutine
- `cmlAddCrystaldp` — subroutine
- `cmlAddCrystalsp` — subroutine
- `cmlAddEigenValuedp` — subroutine
- `cmlAddEigenValuesp` — subroutine
- `cmlAddEigenValueVectorCmplxdp` — subroutine
- `cmlAddEigenValueVectorCmplxsp` — subroutine
- `cmlAddEigenValueVectordp` — subroutine
- `cmlAddEigenValueVectorsp` — subroutine
- `cmlAddKPointdp` — subroutine
- `cmlAddKPointsp` — subroutine
- `cmlAddLatticedp` — subroutine
- `cmlAddLatticesp` — subroutine
- `cmlAddLength_dp` — subroutine
- `cmlAddLength_sp` — subroutine
- `cmlAddMolecule_3_dp` — subroutine
- `cmlAddMolecule_3_dp_sh` — subroutine
- `cmlAddMolecule_3_sp` — subroutine
- `cmlAddMolecule_3_sp_sh` — subroutine
- `cmlAddMoleculedp` — subroutine
- `cmlAddMoleculedp_sh` — subroutine
- `cmlAddMoleculesp` — subroutine
- `cmlAddMoleculesp_sh` — subroutine
- `cmlAddParticles_3_dp` — subroutine
- `cmlAddParticles_3_dp_sh` — subroutine
- `cmlAddParticles_3_sp` — subroutine
- `cmlAddParticles_3_sp_sh` — subroutine
- `cmlAddParticlesdp` — subroutine
- `cmlAddParticlesdp_sh` — subroutine
- `cmlAddParticlessp` — subroutine
- `cmlAddParticlessp_sh` — subroutine
- `cmlAddSymmetrydp` — subroutine
- `cmlAddSymmetryNoOps` — subroutine
- `cmlAddSymmetrysp` — subroutine
- `cmlAddTorsion_dp` — subroutine
- `cmlAddTorsion_sp` — subroutine
- `cmlEndBand` — subroutine
- `cmlEndKpoint` — subroutine
- `cmlStartBand` — subroutine
- `cmlStartKPointdp` — subroutine
- `cmlStartKPointsp` — subroutine
- `cnt_recs_per_line` — function
- `comp_addr_logi` — function
- `comp_addr_phys1` — function
- `comp_addr_phys2` — function
- `comp_addr_phys3` — function
- `comp_addr_phys_orig` — function
- `compact` — subroutine
- `compeuler` — subroutine
- `compose2dshift3d` — subroutine
- `compose3d2d` — subroutine
- `concat_complex_dp_str` — function
- `concat_complex_sp_str` — function
- `concat_int_str` — function
- `concat_logical_str` — function
- `concat_real_dp_str` — function
- `concat_real_sp_str` — function
- `concat_str_complex_dp` — function
- `concat_str_complex_sp` — function
- `concat_str_int` — function
- `concat_str_logical` — function
- `concat_str_real_dp` — function
- `concat_str_real_sp` — function
- `constructor` — function
- `constructor` — function
- `constructor` — function
- `constructor` — function
- `constructor` — function
- `constructor` — function
- `constructor` — function
- `constructor` — function
- `constructor` — function
- `constructor` — function
- `constructor_1` — function
- `constructor_1` — function
- `constructor_2` — function
- `constructor_2` — function
- `conv2rank_weights` — subroutine
- `copy` — subroutine
- `copy` — subroutine
- `copy` — subroutine
- `copy` — subroutine
- `copy_imgarr` — function
- `corrs2weights` — function
- `cosedge_1` — function
- `cosedge_2` — function
- `cosedge_inner_1` — function
- `cosedge_inner_2` — function
- `cosedge_r2_2d` — function
- `cosedge_r2_3d` — function
- `countcomplexdp` — function
- `countcomplexsp` — function
- `countinteger` — function
- `countlogical` — function
- `countrealdp` — function
- `countrealsp` — function
- `countstring` — function
- `create_hist_vector` — subroutine
- `createEmptyElement` — function
- `createEmptyElementNS` — function
- `createNode` — function
- `cross` — function
- `csq_1` — function
- `csq_2` — function
- `csv_field` — function
- `cyci_1d` — function
- `cyci_1d_static` — function
- `dble2str` — function
- `dealloc_chash` — subroutine
- `dealloc_hash` — subroutine
- `dealloc_imgarr` — subroutine
- `declare_coarray_jobs_finished` — subroutine
- `deg2rad_dp` — function
- `deg2rad_sp` — function
- `del_files_1` — subroutine
- `del_files_2` — subroutine
- `delete` — subroutine
- `delete` — subroutine
- `delete_2Dclustering` — subroutine
- `delete_3Dalignment` — subroutine
- `delete_entry` — subroutine
- `destroy_string_list` — subroutine
- `destroy_vs` — subroutine
- `destroyDocumentType` — subroutine
- `destroyDOMConfig` — subroutine
- `destroyElementOrAttribute` — subroutine
- `destroyEntityOrNotation` — subroutine
- `destroyNode` — subroutine
- `destroyNodeContents` — subroutine
- `detect_and_add_dirs` — subroutine
- `detect_peak_thres_1` — subroutine
- `detect_peak_thres_2` — subroutine
- `detect_peak_thres_fdr` — subroutine
- `detect_peak_thres_for_npeaks` — subroutine
- `detect_peak_thres_sortmeans` — subroutine
- `dgelss` — subroutine
- `dgelsy` — subroutine
- `dgesvd` — subroutine
- `dgetrf` — subroutine
- `dgetri` — subroutine
- `discrete_read_imgbatch` — subroutine
- `discrete_read_imgbatch_source` — subroutine
- `dispatch_task_to_persistent_worker` — subroutine
- `dist_avg` — function
- `dist_btw_farthest` — function
- `dist_btw_nearest` — function
- `dist_eval_fun` — function
- `dists2order` — function
- `dists2scores_percen` — subroutine
- `dmat2smat` — function
- `dnrm2` — function
- `double_thresh` — subroutine
- `dposv` — subroutine
- `dsyev` — subroutine
- `dynfind` — function
- `e1get` — function
- `e1set` — subroutine
- `e2get` — function
- `e2set` — subroutine
- `e3get` — function
- `e3set` — subroutine
- `eigh_sp` — subroutine
- `eigsrt_dp` — subroutine
- `eigsrt_sp` — subroutine
- `elim_dup` — subroutine
- `ellipse` — subroutine
- `emit_refine3D_stage_cfg` — subroutine
- `equispaced_vals` — function
- `estimate_lplim3D` — subroutine
- `estimate_lplim_1` — subroutine
- `estimate_lplim_2` — subroutine
- `estimate_lplims2D` — subroutine
- `euclid_dp` — function
- `euclid_sp_1` — function
- `euclid_sp_2` — function
- `euldist` — function
- `eulprob_corr_switch` — function
- `eulprob_dist_switch` — function
- `eval_apod` — function
- `eval_instr` — function
- `exists` — function
- `extract_imgarr` — function
- `extractDataAttNSChArr` — subroutine
- `extractDataAttNSChMat` — subroutine
- `extractDataAttNSChSca` — subroutine
- `extractDataAttNSCmplxDpArr` — subroutine
- `extractDataAttNSCmplxDpMat` — subroutine
- `extractDataAttNSCmplxDpSca` — subroutine
- `extractDataAttNSCmplxSpArr` — subroutine
- `extractDataAttNSCmplxSpMat` — subroutine
- `extractDataAttNSCmplxSpSca` — subroutine
- `extractDataAttNSIntArr` — subroutine
- `extractDataAttNSIntMat` — subroutine
- `extractDataAttNSIntSca` — subroutine
- `extractDataAttNSLgArr` — subroutine
- `extractDataAttNSLgMat` — subroutine
- `extractDataAttNSLgSca` — subroutine
- `extractDataAttNSRealDpArr` — subroutine
- `extractDataAttNSRealDpMat` — subroutine
- `extractDataAttNSRealDpSca` — subroutine
- `extractDataAttNSRealSpArr` — subroutine
- `extractDataAttNSRealSpMat` — subroutine
- `extractDataAttNSRealSpSca` — subroutine
- `extractDataAttributeChArr` — subroutine
- `extractDataAttributeChMat` — subroutine
- `extractDataAttributeChSca` — subroutine
- `extractDataAttributeCmplxDpArr` — subroutine
- `extractDataAttributeCmplxDpMat` — subroutine
- `extractDataAttributeCmplxDpSca` — subroutine
- `extractDataAttributeCmplxSpArr` — subroutine
- `extractDataAttributeCmplxSpMat` — subroutine
- `extractDataAttributeCmplxSpSca` — subroutine
- `extractDataAttributeIntArr` — subroutine
- `extractDataAttributeIntMat` — subroutine
- `extractDataAttributeIntSca` — subroutine
- `extractDataAttributeLgArr` — subroutine
- `extractDataAttributeLgMat` — subroutine
- `extractDataAttributeLgSca` — subroutine
- `extractDataAttributeRealDpArr` — subroutine
- `extractDataAttributeRealDpMat` — subroutine
- `extractDataAttributeRealDpSca` — subroutine
- `extractDataAttributeRealSpArr` — subroutine
- `extractDataAttributeRealSpMat` — subroutine
- `extractDataAttributeRealSpSca` — subroutine
- `extractDataContentChArr` — subroutine
- `extractDataContentChMat` — subroutine
- `extractDataContentChSca` — subroutine
- `extractDataContentCmplxDpArr` — subroutine
- `extractDataContentCmplxDpMat` — subroutine
- `extractDataContentCmplxDpSca` — subroutine
- `extractDataContentCmplxSpArr` — subroutine
- `extractDataContentCmplxSpMat` — subroutine
- `extractDataContentCmplxSpSca` — subroutine
- `extractDataContentIntArr` — subroutine
- `extractDataContentIntMat` — subroutine
- `extractDataContentIntSca` — subroutine
- `extractDataContentLgArr` — subroutine
- `extractDataContentLgMat` — subroutine
- `extractDataContentLgSca` — subroutine
- `extractDataContentRealDpArr` — subroutine
- `extractDataContentRealDpMat` — subroutine
- `extractDataContentRealDpSca` — subroutine
- `extractDataContentRealSpArr` — subroutine
- `extractDataContentRealSpMat` — subroutine
- `extractDataContentRealSpSca` — subroutine
- `fclose` — subroutine
- `fdim` — function
- `file2drarr` — function
- `file2imat` — subroutine
- `file2lmat` — subroutine
- `file2rarr` — function
- `file2rmat` — subroutine
- `fileiochk` — subroutine
- `filelength` — function
- `filepath_1` — function
- `filepath_2` — function
- `filepath_3` — function
- `filepath_4` — function
- `find_1` — subroutine
- `find_2` — subroutine
- `find_closest_proj` — function
- `find_larger_magic_box` — function
- `find_magic_box` — function
- `find_magic_boxes4scale` — function
- `find_medoids` — subroutine
- `findloc_str_1` — function
- `findloc_str_2` — function
- `fit_lsq_plane` — subroutine
- `fit_straight_line` — subroutine
- `fmtsymstr` — function
- `fname2ext` — function
- `fname2format` — function
- `fname_new_ext_1` — function
- `fname_new_ext_2` — function
- `fopen` — subroutine
- `fortran_logical` — function
- `fortran_quote` — function
- `fortran_symbol_from_string` — function
- `fplane_get_cmplx` — function
- `fplane_get_ctfsq` — function
- `free_all_cunits` — subroutine
- `funcs` — function
- `funcs` — function
- `funcs` — function
- `funcs` — function
- `funit_size` — function
- `gasdev_2` — function
- `gasdev_3` — function
- `gau_rnd_shift` — subroutine
- `gauss2Dfit` — subroutine
- `gauwfun` — function
- `gcd` — function
- `gen_c1` — subroutine
- `gen_job_descr` — subroutine
- `generate_array_script` — subroutine
- `generate_script_1` — subroutine
- `generate_script_2` — subroutine
- `generate_script_3` — subroutine
- `generate_script_4` — subroutine
- `generate_scripts_subprojects` — subroutine
- `geodesic_dist_trace` — function
- `get_1` — function
- `get_2` — function
- `get_2Dshift` — function
- `get_3Dshift` — function
- `get_all_subgrps_descr` — function
- `get_axis_angle` — subroutine
- `get_ctfvars` — function
- `get_ctfvars` — function
- `get_euler` — function
- `get_eullims` — function
- `get_fbody_1` — function
- `get_fbody_2` — function
- `get_find` — function
- `get_find_at_crit` — function
- `get_find_at_res` — function
- `get_fpath` — function
- `get_jobs_status` — subroutine
- `get_key` — function
- `get_key` — function
- `get_keys` — function
- `get_keys` — function
- `get_last_string` — function
- `get_ldim` — function
- `get_lfny` — function
- `get_lhp` — function
- `get_llp` — function
- `get_lp` — function
- `get_mat` — function
- `get_normal` — function
- `get_nsubgrp` — function
- `get_open_funits` — subroutine
- `get_pgrp` — function
- `get_pixel_pos` — subroutine
- `get_relative_path` — function
- `get_resarr` — function
- `get_spat_freq` — function
- `get_static` — subroutine
- `get_static` — subroutine
- `get_str` — function
- `get_stream_done_stack` — subroutine
- `get_stream_fail_stack` — subroutine
- `get_subgrp` — function
- `get_subgrp_descr` — function
- `get_sym_rmat` — subroutine
- `get_symori` — subroutine
- `get_Whalf` — function
- `getAttribute_len` — function
- `getAttributesNS_len` — function
- `getdata_len` — function
- `getdocumentURI_len` — function
- `getevensym` — subroutine
- `getGCstate` — function
- `getInputEncoding_len` — function
- `getInternalSubset_len` — function
- `getisId_DOM` — function
- `getLength_characterdata` — function
- `getLength_nl` — function
- `getLength_nnm` — function
- `getLocalName_len` — function
- `getname_len` — function
- `getNamespaceURI_len` — function
- `getnodeName_len` — function
- `getNodeValue_len` — function
- `getnotationName_len` — function
- `getPrefix_len` — function
- `getpublicId_len` — function
- `getstringValue_len` — function
- `getsystemId_len` — function
- `gettagName_len` — function
- `getTarget_len` — function
- `getter_1` — subroutine
- `getter_1` — subroutine
- `getter_2` — subroutine
- `getter_2` — subroutine
- `getter_3` — subroutine
- `getter_3` — subroutine
- `getTextContent_len` — function
- `getValue_DOM` — function
- `getXmlEncoding_len` — function
- `getXmlVersionEnum` — function
- `great_circle_samples` — function
- `greedy_sampling_1` — function
- `greedy_sampling_2` — function
- `hac_1d` — subroutine
- `hac_1d_fast` — subroutine
- `hac_med_thres` — subroutine
- `haloween_end` — subroutine
- `hann_apod` — function
- `hann_instr` — function
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
- `has_subgrp` — function
- `hash2str` — function
- `hermitian_eigh_z` — subroutine
- `hermitian_invert_dp` — subroutine
- `hermitian_invert_z` — subroutine
- `hermitian_solve_dp` — subroutine
- `hermitian_solve_z` — subroutine
- `hough_line` — subroutine
- `hpsort_1` — subroutine
- `hpsort_2` — subroutine
- `hpsort_3` — subroutine
- `hpsort_4` — subroutine
- `hpsort_5` — subroutine
- `hpsort_6` — subroutine
- `image_stack` — type
- `imat2file` — subroutine
- `init` — subroutine
- `init_refine3D_iteration` — subroutine
- `init_string_list` — subroutine
- `inplrotdist` — function
- `insert` — subroutine
- `instr` — function
- `int2str` — function
- `int2str` — function
- `int2str_pad` — function
- `irnd` — function
- `irnd_gasdev` — function
- `irnd_gau` — function
- `irnd_pair` — subroutine
- `irnd_uni` — function
- `irnd_uni_pair` — function
- `is` — function
- `is_even_2` — function
- `is_particle` — function
- `isCharData` — function
- `isthere` — function
- `isthere` — function
- `item_nl` — function
- `item_nnm` — function
- `jacobi_dp` — subroutine
- `jacobi_sp` — subroutine
- `join_imgarrs` — function
- `kb_apod` — function
- `kb_instr` — function
- `kill` — subroutine
- `kill` — subroutine
- `kill` — subroutine
- `kill` — subroutine
- `kill` — subroutine
- `kill` — subroutine
- `kill` — subroutine
- `kill` — subroutine
- `kill_candidate_buffer` — subroutine
- `kill_candidate_store` — subroutine
- `kill_chash` — subroutine
- `kill_hash` — subroutine
- `killimgbatch` — subroutine
- `kmlAddaddress` — subroutine
- `kmlAddaltitude` — subroutine
- `kmlAddaltitudeMode` — subroutine
- `kmlAddBalloonStyle` — subroutine
- `kmlAddbegin` — subroutine
- `kmlAddBgColor_c` — subroutine
- `kmlAddBgColor_h` — subroutine
- `kmlAddChart_dp` — subroutine
- `kmlAddChart_sp` — subroutine
- `kmlAddColor_c` — subroutine
- `kmlAddColor_h` — subroutine
- `kmlAddcolorMode` — subroutine
- `kmlAddcookie` — subroutine
- `kmlAddCoordinates_array_dp` — subroutine
- `kmlAddCoordinates_array_sp` — subroutine
- `kmlAddcoordinates_dp` — subroutine
- `kmlAddcoordinates_sp` — subroutine
- `kmlAdddescription_ch` — subroutine
- `kmlAdddescription_dp` — subroutine
- `kmlAdddescription_sp` — subroutine
- `kmlAdddrawOrder` — subroutine
- `kmlAddeast` — subroutine
- `kmlAddend` — subroutine
- `kmlAddextrude` — subroutine
- `kmlAddfill` — subroutine
- `kmlAddflyToView` — subroutine
- `kmlAddheading` — subroutine
- `kmlAddhref` — subroutine
- `kmlAddIcon_href` — subroutine
- `kmlAddIcon_refresh` — subroutine
- `kmlAddIcon_view` — subroutine
- `kmlAddInnerBoundary_1d_dp` — subroutine
- `kmlAddInnerBoundary_1d_sp` — subroutine
- `kmlAddInnerBoundary_2d_dp` — subroutine
- `kmlAddInnerBoundary_2d_sp` — subroutine
- `kmlAddkey` — subroutine
- `kmlAddLabelStyle_color` — subroutine
- `kmlAddLabelStyle_scale` — subroutine
- `kmlAddlatitude_dp` — subroutine
- `kmlAddlatitude_sp` — subroutine
- `kmlAddLineStyle` — subroutine
- `kmlAddlinkDescription` — subroutine
- `kmlAddlinkName` — subroutine
- `kmlAddlistItemType` — subroutine
- `kmlAddListStyle_bgcolor` — subroutine
- `kmlAddlongitude_dp` — subroutine
- `kmlAddlongitude_sp` — subroutine
- `kmlAddmaxAltitude` — subroutine
- `kmlAddmaxFadeExtent` — subroutine
- `kmlAddmaxLodPixels` — subroutine
- `kmlAddmessage` — subroutine
- `kmlAddminAltitude` — subroutine
- `kmlAddminFadeExtent` — subroutine
- `kmlAddminLodPixels` — subroutine
- `kmlAddminRefreshPeriod` — subroutine
- `kmlAddname` — subroutine
- `kmlAddNamespace` — subroutine
- `kmlAddnorth` — subroutine
- `kmlAddopen` — subroutine
- `kmlAddoutline` — subroutine
- `kmlAddphoneNumber` — subroutine
- `kmlAddPoint_dp` — subroutine
- `kmlAddPoint_sp` — subroutine
- `kmlAddPointStyle_s_dp` — subroutine
- `kmlAddPointStyle_s_int` — subroutine
- `kmlAddPointStyle_s_none_hd_dp` — subroutine
- `kmlAddPointStyle_s_none_hd_int` — subroutine
- `kmlAddPointStyle_s_none_hd_none` — subroutine
- `kmlAddPointStyle_s_none_hd_sp` — subroutine
- `kmlAddPointStyle_s_sp` — subroutine
- `kmlAddrange` — subroutine
- `kmlAddrefreshInterval` — subroutine
- `kmlAddrefreshMode` — subroutine
- `kmlAddrefreshVisibility` — subroutine
- `kmlAddrequest` — subroutine
- `kmlAddroll` — subroutine
- `kmlAddrotation` — subroutine
- `kmlAddScale_dp` — subroutine
- `kmlAddScale_int` — subroutine
- `kmlAddScale_sp` — subroutine
- `kmlAddsouth` — subroutine
- `kmlAddstate` — subroutine
- `kmlAddstyleUrl` — subroutine
- `kmlAddtargetHref` — subroutine
- `kmlAddtessellate` — subroutine
- `kmlAddTextColor_c` — subroutine
- `kmlAddTextColor_h` — subroutine
- `kmlAddtilt` — subroutine
- `kmlAddviewBoundScale` — subroutine
- `kmlAddviewrefreshMode` — subroutine
- `kmlAddviewRefreshTime` — subroutine
- `kmlAddvisibility` — subroutine
- `kmlAddwest` — subroutine
- `kmlAddwhen` — subroutine
- `kmlAddwidth` — subroutine
- `kmlCloseAddressDetail` — subroutine
- `kmlCloseBalloonStyle` — subroutine
- `kmlCloseChange` — subroutine
- `kmlCloseColorStyle` — subroutine
- `kmlCloseContainer` — subroutine
- `kmlClosecoordinates` — subroutine
- `kmlCloseCreate` — subroutine
- `kmlCloseDelete` — subroutine
- `kmlCloseDocument` — subroutine
- `kmlCloseFeature` — subroutine
- `kmlCloseFolder` — subroutine
- `kmlCloseGeometry` — subroutine
- `kmlCloseGeometryCollection` — subroutine
- `kmlCloseGroundOverlay` — subroutine
- `kmlCloseIcon` — subroutine
- `kmlCloseIconStyle` — subroutine
- `kmlCloseItemIcon` — subroutine
- `kmlCloseLabelStyle` — subroutine
- `kmlCloseLatLonAltBox` — subroutine
- `kmlCloseLatLonBox` — subroutine
- `kmlCloseLineStyle` — subroutine
- `kmlCloseLink` — subroutine
- `kmlCloseListStyle` — subroutine
- `kmlCloseLocation` — subroutine
- `kmlCloseLod` — subroutine
- `kmlCloseLookAt` — subroutine
- `kmlCloseModel` — subroutine
- `kmlCloseMultiGeometry` — subroutine
- `kmlCloseNetworkLink` — subroutine
- `kmlCloseNetworkLinkControl` — subroutine
- `kmlCloseObjArrayField` — subroutine
- `kmlCloseObject` — subroutine
- `kmlCloseObjField` — subroutine
- `kmlCloseOrientation` — subroutine
- `kmlCloseOverlay` — subroutine
- `kmlClosePair` — subroutine
- `kmlClosePlacemark` — subroutine
- `kmlClosePoint` — subroutine
- `kmlClosePolyStyle` — subroutine
- `kmlCloseRegion` — subroutine
- `kmlCloseResponse` — subroutine
- `kmlCloseScale` — subroutine
- `kmlCloseSchema` — subroutine
- `kmlCloseSchemaField` — subroutine
- `kmlCloseScreenOverlay` — subroutine
- `kmlCloseSimpleArrayField` — subroutine
- `kmlCloseSimpleField` — subroutine
- `kmlCloseSnippet` — subroutine
- `kmlCloseStatus` — subroutine
- `kmlCloseStyleMap` — subroutine
- `kmlCloseStyleSelector` — subroutine
- `kmlCloseTimePrimitive` — subroutine
- `kmlCloseTimeSpan` — subroutine
- `kmlCloseTimeStamp` — subroutine
- `kmlCloseUpdate` — subroutine
- `kmlCloseUrl` — subroutine
- `kmlcreateCells3_dp` — subroutine
- `kmlcreateCells3_sp` — subroutine
- `kmlCreateCells_dp` — subroutine
- `kmlCreateCells_longlat2_dp` — subroutine
- `kmlCreateCells_longlat2_sp` — subroutine
- `kmlCreateCells_longlat_dp` — subroutine
- `kmlCreateCells_longlat_sp` — subroutine
- `kmlCreateCells_sp` — subroutine
- `kmlCreateContours_longlat2_sp` — subroutine
- `kmlCreateContours_longlat_sp` — subroutine
- `kmlCreateContours_old` — subroutine
- `kmlCreateContours_sp` — subroutine
- `kmlCreateLine_1d_dp` — subroutine
- `kmlCreateLine_1d_sp` — subroutine
- `kmlCreateLine_2d_dp` — subroutine
- `kmlCreateLine_2d_sp` — subroutine
- `kmlCreateLine_seg_sh_dp` — subroutine
- `kmlCreateLineStyle_old` — subroutine
- `kmlCreatePoints_0d_dp` — subroutine
- `kmlCreatePoints_0d_sp` — subroutine
- `kmlCreatePoints_1d_dp` — subroutine
- `kmlCreatePoints_1d_sp` — subroutine
- `kmlCreatePoints_2d_dp` — subroutine
- `kmlCreatePoints_2d_sp` — subroutine
- `kmlCreatePointStyleOld` — subroutine
- `kmlCreateRGBCells_dp` — subroutine
- `kmlCreateRGBCells_sp` — subroutine
- `kmlGetColor_byName` — function
- `kmlGetColor_index` — function
- `kmlOpenAddressDetail` — subroutine
- `kmlOpenBalloonStyle` — subroutine
- `kmlOpenChange` — subroutine
- `kmlOpenColorStyle` — subroutine
- `kmlOpenContainer` — subroutine
- `kmlOpencoordinates` — subroutine
- `kmlOpenCreate` — subroutine
- `kmlOpenDelete` — subroutine
- `kmlOpenDocument` — subroutine
- `kmlOpenFeature` — subroutine
- `kmlOpenFolder` — subroutine
- `kmlOpenGeometry` — subroutine
- `kmlOpenGeometryCollection` — subroutine
- `kmlOpenGroundOverlay` — subroutine
- `kmlOpenIcon` — subroutine
- `kmlOpenIconStyle` — subroutine
- `kmlOpenItemIcon` — subroutine
- `kmlOpenLabelStyle` — subroutine
- `kmlOpenLatLonAltBox` — subroutine
- `kmlOpenLatLonBox` — subroutine
- `kmlOpenLineStyle` — subroutine
- `kmlOpenLink` — subroutine
- `kmlOpenListStyle` — subroutine
- `kmlOpenLocation` — subroutine
- `kmlOpenLod` — subroutine
- `kmlOpenLookAt` — subroutine
- `kmlOpenModel` — subroutine
- `kmlOpenMultiGeometry` — subroutine
- `kmlOpenNetworkLink` — subroutine
- `kmlOpenNetworkLinkControl` — subroutine
- `kmlOpenObjArrayField` — subroutine
- `kmlOpenObject` — subroutine
- `kmlOpenObjField` — subroutine
- `kmlOpenOrientation` — subroutine
- `kmlOpenOverlay` — subroutine
- `kmlOpenPair` — subroutine
- `kmlOpenPlacemark` — subroutine
- `kmlOpenPoint` — subroutine
- `kmlOpenPolyStyle` — subroutine
- `kmlOpenRegion` — subroutine
- `kmlOpenResponse` — subroutine
- `kmlOpenScale` — subroutine
- `kmlOpenSchema` — subroutine
- `kmlOpenSchemaField` — subroutine
- `kmlOpenScreenOverlay` — subroutine
- `kmlOpenSimpleArrayField` — subroutine
- `kmlOpenSimpleField` — subroutine
- `kmlOpenSnippet` — subroutine
- `kmlOpenStatus` — subroutine
- `kmlOpenStyleMap` — subroutine
- `kmlOpenStyleSelector` — subroutine
- `kmlOpenTimePrimitive` — subroutine
- `kmlOpenTimeSpan` — subroutine
- `kmlOpenTimeStamp` — subroutine
- `kmlOpenUpdate` — subroutine
- `kmlOpenUrl` — subroutine
- `kmlOutputContours` — subroutine
- `kmlStartPolygon_1d_dp` — subroutine
- `kmlStartPolygon_1d_sp` — subroutine
- `kmlStartPolygon_2d_dp` — subroutine
- `kmlStartPolygon_2d_sp` — subroutine
- `kstwo` — subroutine
- `l1dist_dp` — function
- `l1dist_sp` — function
- `lapack_stop` — subroutine
- `lcg` — function
- `lex_sort_1` — subroutine
- `lex_sort_2` — subroutine
- `lin_apod` — function
- `lin_instr` — function
- `list_of_ints2arr` — function
- `lmat2file` — subroutine
- `locate_1` — function
- `locate_2` — function
- `lookup` — function
- `lookupNamespaceURI_len` — function
- `lookupPrefix_len` — function
- `loop_lims` — function
- `lowercase` — function
- `make_c_and_d` — subroutine
- `make_dirnames` — function
- `make_filenames` — function
- `make_i_relion` — subroutine
- `make_i_spider` — subroutine
- `make_o` — subroutine
- `make_t` — subroutine
- `map3dshift22d` — subroutine
- `map_str_nrs` — function
- `mask2inds` — function
- `masked_swap_rm` — subroutine
- `masked_swap_rs` — subroutine
- `masked_swap_rv` — subroutine
- `matcreate` — function
- `materialize_seed_shift` — subroutine
- `matextract` — function
- `matinv_dp` — subroutine
- `matinv_sp` — subroutine
- `matrixtocomplexdp` — subroutine
- `matrixtocomplexsp` — subroutine
- `matrixtointeger` — subroutine
- `matrixtological` — subroutine
- `matrixtorealdp` — subroutine
- `matrixtorealsp` — subroutine
- `matrixtostring` — subroutine
- `maxheap_sift_down` — subroutine
- `maxnloc` — function
- `median` — function
- `median_nocopy` — function
- `medoid_from_dmat` — subroutine
- `medoid_from_smat` — subroutine
- `medoid_ranking_from_smat` — subroutine
- `merge_dmats_1` — function
- `merge_dmats_2` — function
- `merge_idx` — subroutine
- `merge_idx` — subroutine
- `merge_smats` — function
- `mergesort_idx` — subroutine
- `mergesort_idx` — subroutine
- `min3` — function
- `minnloc` — function
- `mirror2d` — subroutine
- `mirror3d` — subroutine
- `mnomal` — function
- `mnorm_smp` — function
- `mode` — subroutine
- `moment_1` — subroutine
- `moment_2` — subroutine
- `moment_3` — subroutine
- `moment_4` — subroutine
- `moment_serial` — subroutine
- `mostOfLineStyle` — subroutine
- `mostOfPointStyle` — subroutine
- `move_files2dir` — subroutine
- `move_files_in_cwd` — subroutine
- `move_key_to_front_1` — subroutine
- `move_key_to_front_2` — subroutine
- `multinomal` — function
- `myacos_dp` — function
- `myacos_sp` — function
- `mycabs` — function
- `ne_mnomal_iarr` — subroutine
- `ne_ran_iarr` — subroutine
- `nearest_proj_neighbors_1` — subroutine
- `nearest_proj_neighbors_2` — subroutine
- `nearest_sym_neighbors` — subroutine
- `neigh_4_3D_1` — subroutine
- `neigh_4_3D_2` — subroutine
- `neigh_8_1` — subroutine
- `neigh_8_2` — subroutine
- `neigh_8_3` — subroutine
- `neigh_8_3D_0` — subroutine
- `neigh_8_3D_1` — subroutine
- `neigh_8_3D_2` — subroutine
- `new` — subroutine
- `new` — subroutine
- `new` — subroutine
- `new` — subroutine
- `new_1` — subroutine
- `new_1` — subroutine
- `new_2` — subroutine
- `new_2` — subroutine
- `new_fixed_candidate_store` — subroutine
- `new_ori` — subroutine
- `new_ragged_candidate_store` — subroutine
- `nextPow2` — function
- `nlines` — function
- `nn_apod` — function
- `nn_instr` — function
- `non_max_supp` — subroutine
- `norm_2_dp` — function
- `norm_2_sp` — function
- `norm_corr` — function
- `norm_corr_8` — function
- `norm_key` — function
- `normal_solve` — subroutine
- `normalize_1` — subroutine
- `normalize_2` — subroutine
- `normalize_3` — subroutine
- `normalize_4` — subroutine
- `normalize_minmax_1` — subroutine
- `normalize_minmax_2` — subroutine
- `normalize_sigm_1` — subroutine
- `normalize_sigm_2` — subroutine
- `normalize_sigm_3` — subroutine
- `ori2chash` — function
- `ori2json` — subroutine
- `ori2prec` — subroutine
- `ori2str` — function
- `ori_from_rotmat` — subroutine
- `otsu_1` — subroutine
- `otsu_2` — subroutine
- `otsu_3` — subroutine
- `otsu_img` — subroutine
- `otsu_robust_fast` — subroutine
- `outerprod_d` — function
- `outerprod_r` — function
- `outputContourLines` — subroutine
- `outputContourRegions` — subroutine
- `p1_lt_p2` — function
- `pack_imgarr` — function
- `parameterChArrSh` — subroutine
- `parameterChArrSi` — subroutine
- `parameterChMatSh` — subroutine
- `parameterChMatSi` — subroutine
- `parameterChSca` — subroutine
- `parameterCmplxDpArrSh` — subroutine
- `parameterCmplxDpArrSi` — subroutine
- `parameterCmplxDpMatSh` — subroutine
- `parameterCmplxDpMatSi` — subroutine
- `parameterCmplxDpSca` — subroutine
- `parameterCmplxSpArrSh` — subroutine
- `parameterCmplxSpArrSi` — subroutine
- `parameterCmplxSpMatSh` — subroutine
- `parameterCmplxSpMatSi` — subroutine
- `parameterCmplxSpSca` — subroutine
- `parameterIntArrSh` — subroutine
- `parameterIntArrSi` — subroutine
- `parameterIntMatSh` — subroutine
- `parameterIntMatSi` — subroutine
- `parameterIntSca` — subroutine
- `parameterLgArrSh` — subroutine
- `parameterLgArrSi` — subroutine
- `parameterLgMatSh` — subroutine
- `parameterLgMatSi` — subroutine
- `parameterLgSca` — subroutine
- `parameterRealDpArrSh` — subroutine
- `parameterRealDpArrSi` — subroutine
- `parameterRealDpMatSh` — subroutine
- `parameterRealDpMatSi` — subroutine
- `parameterRealDpSca` — subroutine
- `parameterRealSpArrSh` — subroutine
- `parameterRealSpArrSi` — subroutine
- `parameterRealSpMatSh` — subroutine
- `parameterRealSpMatSi` — subroutine
- `parameterRealSpSca` — subroutine
- `parse_cmdline` — subroutine
- `parsestr` — subroutine
- `partial_shuffle_1` — subroutine
- `partial_shuffle_2` — subroutine
- `peakfinder` — function
- `pearsn_1` — function
- `pearsn_2` — function
- `pearsn_3` — function
- `pearsn_serial` — function
- `pearsn_serial_8` — function
- `phase_angle` — function
- `pixels_dist_1` — function
- `pixels_dist_2` — function
- `plane_from_points` — function
- `power_sampling` — subroutine
- `pparms2str` — function
- `prec2ori` — subroutine
- `prep_part_jobs` — subroutine
- `prepare_tree_sub_distmat` — subroutine
- `prepimgbatch` — subroutine
- `print` — subroutine
- `print` — subroutine
- `print_gpu_specs` — subroutine
- `print_jobs_status` — subroutine
- `print_key_val_pair_1` — subroutine
- `print_key_val_pair_2` — subroutine
- `print_key_val_pairs` — subroutine
- `print_magic_box_range` — subroutine
- `print_mat` — subroutine
- `print_ori` — subroutine
- `probks` — function
- `progress` — subroutine
- `progress_gfortran` — subroutine
- `projz` — subroutine
- `propertyChArrSh` — subroutine
- `propertyChArrSi` — subroutine
- `propertyChMatSh` — subroutine
- `propertyChMatSi` — subroutine
- `propertyChSca` — subroutine
- `propertyCmplxDpArrSh` — subroutine
- `propertyCmplxDpArrSi` — subroutine
- `propertyCmplxDpMatSh` — subroutine
- `propertyCmplxDpMatSi` — subroutine
- `propertyCmplxDpSca` — subroutine
- `propertyCmplxSpArrSh` — subroutine
- `propertyCmplxSpArrSi` — subroutine
- `propertyCmplxSpMatSh` — subroutine
- `propertyCmplxSpMatSi` — subroutine
- `propertyCmplxSpSca` — subroutine
- `propertyIntArrSh` — subroutine
- `propertyIntArrSi` — subroutine
- `propertyIntMatSh` — subroutine
- `propertyIntMatSi` — subroutine
- `propertyIntSca` — subroutine
- `propertyLgArrSh` — subroutine
- `propertyLgArrSi` — subroutine
- `propertyLgMatSh` — subroutine
- `propertyLgMatSi` — subroutine
- `propertyLgSca` — subroutine
- `propertyRealDpArrSh` — subroutine
- `propertyRealDpArrSi` — subroutine
- `propertyRealDpMatSh` — subroutine
- `propertyRealDpMatSi` — subroutine
- `propertyRealDpSca` — subroutine
- `propertyRealSpArrSh` — subroutine
- `propertyRealSpArrSi` — subroutine
- `propertyRealSpMatSh` — subroutine
- `propertyRealSpMatSi` — subroutine
- `propertyRealSpSca` — subroutine
- `PseudoAttributeArrayCh` — subroutine
- `PseudoAttributeArrayCmplxDp` — subroutine
- `PseudoAttributeArrayCmplxSp` — subroutine
- `PseudoAttributeArrayInt` — subroutine
- `PseudoAttributeArrayLg` — subroutine
- `PseudoAttributeArrayRealDp` — subroutine
- `PseudoAttributeArrayRealSp` — subroutine
- `PseudoAttributeMatrixCh` — subroutine
- `PseudoAttributeMatrixCmplxDp` — subroutine
- `PseudoAttributeMatrixCmplxSp` — subroutine
- `PseudoAttributeMatrixInt` — subroutine
- `PseudoAttributeMatrixLg` — subroutine
- `PseudoAttributeMatrixRealDp` — subroutine
- `PseudoAttributeMatrixRealSp` — subroutine
- `PseudoAttributeScalarCmplxDp` — subroutine
- `PseudoAttributeScalarCmplxSp` — subroutine
- `PseudoAttributeScalarInt` — subroutine
- `PseudoAttributeScalarLg` — subroutine
- `PseudoAttributeScalarRealDp` — subroutine
- `PseudoAttributeScalarRealSp` — subroutine
- `push_1` — subroutine
- `push_1` — subroutine
- `push_2` — subroutine
- `push_2` — subroutine
- `push_3` — subroutine
- `put_last` — subroutine
- `putNodesInDocument` — subroutine
- `pythag_dp` — function
- `pythag_sp` — function
- `qr_solve` — subroutine
- `qsys_cleanup` — subroutine
- `qsys_declare_part_finished` — subroutine
- `qsys_job_finished` — subroutine
- `qsys_subproject_job_finished` — subroutine
- `qsys_watcher_1` — subroutine
- `qsys_watcher_2` — subroutine
- `qsys_watcher_diag` — subroutine
- `quadri` — function
- `quantize_vec` — subroutine
- `quantize_vec_serial` — subroutine
- `r8po_fa` — subroutine
- `rad2deg_1` — function
- `rad2deg_2` — function
- `ran3` — function
- `ran3arr_1` — subroutine
- `ran3arr_2` — subroutine
- `randn_1` — function
- `randn_2` — function
- `rank_centroid_weights` — subroutine
- `rank_exponent_weights` — subroutine
- `rank_inverse_weights` — subroutine
- `rank_sum_weights` — subroutine
- `raw` — function
- `read` — subroutine
- `read_cavgs_into_imgarr` — function
- `read_exit_code` — subroutine
- `read_filetable` — subroutine
- `read_imgbatch_1` — subroutine
- `read_imgbatch_2` — subroutine
- `read_imgbatch_3` — subroutine
- `read_seed_shift_table` — subroutine
- `read_stk_into_imgarr` — function
- `readfile` — subroutine
- `readline` — subroutine
- `real2str` — function
- `real_dp_str` — function
- `real_sp_str` — function
- `realdp2str` — function
- `realloc_chash` — subroutine
- `realloc_hash` — subroutine
- `realsp2str` — function
- `refine_peak_thres_sortmeans` — subroutine
- `registered_string` — function
- `reject` — subroutine
- `remove` — subroutine
- `remove_last_string` — subroutine
- `remove_node_nl` — subroutine
- `removeNodesFromDocument` — subroutine
- `removepunct` — subroutine
- `removesp` — subroutine
- `reorder_1` — subroutine
- `reorder_2` — subroutine
- `resang` — function
- `reset` — subroutine
- `reset` — subroutine
- `reset_pparms` — subroutine
- `resetParameter` — subroutine
- `reverse_drarr` — subroutine
- `reverse_f` — subroutine
- `reverse_iarr` — subroutine
- `reverse_rarr` — subroutine
- `reverselookup` — function
- `ring_stats` — subroutine
- `rm_from_fbody` — function
- `rm_substr` — subroutine
- `rmat2file` — subroutine
- `rnd_4dim_sphere_pnt` — function
- `rnd_euler_1` — subroutine
- `rnd_euler_1` — subroutine
- `rnd_euler_2` — subroutine
- `rnd_euler_2` — subroutine
- `rnd_euler_3` — subroutine
- `rnd_euler_3` — subroutine
- `rnd_euler_4` — subroutine
- `rnd_euler_5` — subroutine
- `rnd_inds` — function
- `rnd_inpl` — subroutine
- `rnd_ori` — subroutine
- `rnd_shift` — subroutine
- `robust_normalization` — subroutine
- `robust_normalize_minmax` — subroutine
- `robust_scaling` — subroutine
- `robust_sigma_thres` — function
- `robust_z_scores` — function
- `rot_to_asym` — subroutine
- `rotall_to_asym` — subroutine
- `rotate_vec` — function
- `rotmat2d` — subroutine
- `round2even` — function
- `round2odd` — function
- `round_shifts` — subroutine
- `rpl_substr` — subroutine
- `same_energy_euclid` — function
- `sample_bounded_dist` — subroutine
- `sample_likelihood_dist` — subroutine
- `sample_likelihood_index` — subroutine
- `sample_power_dist` — subroutine
- `sauron_ori_parser_1` — subroutine
- `sauron_ori_parser_2` — subroutine
- `sauvola` — subroutine
- `savgol` — subroutine
- `SavitzkyGolay_filter` — subroutine
- `scalartocomplexdp` — subroutine
- `scalartocomplexsp` — subroutine
- `scalartointeger` — subroutine
- `scalartological` — subroutine
- `scalartorealdp` — subroutine
- `scalartorealsp` — subroutine
- `scalartostring` — subroutine
- `schedule_array_jobs` — subroutine
- `schedule_jobs` — subroutine
- `schedule_streaming` — subroutine
- `schedule_subproject_jobs` — subroutine
- `scores2order` — function
- `scores2scores_percen` — subroutine
- `seed_rnd` — subroutine
- `service_persistent_worker_warmup` — subroutine
- `set_1` — subroutine
- `set_1` — subroutine
- `set_1` — subroutine
- `set_2` — subroutine
- `set_2` — subroutine
- `set_2` — subroutine
- `set_3` — subroutine
- `set_3` — subroutine
- `set_3` — subroutine
- `set_4` — subroutine
- `set_5` — subroutine
- `set_class` — subroutine
- `set_ctfvars` — subroutine
- `set_dble` — subroutine
- `set_dfx` — subroutine
- `set_dfy` — subroutine
- `set_euler` — subroutine
- `set_hp` — subroutine
- `set_int` — subroutine
- `set_jobs_status` — subroutine
- `set_offload_device_1` — subroutine
- `set_offload_device_2` — subroutine
- `set_ogid` — subroutine
- `set_real` — subroutine
- `set_refine3D_automsk_policy` — subroutine
- `set_refine3D_balance_policy` — subroutine
- `set_refine3D_filtering_policy` — subroutine
- `set_refine3D_gauref_policy` — subroutine
- `set_refine3D_mode_policy` — subroutine
- `set_refine3D_stage_controls` — subroutine
- `set_refine3D_symmetry_policy` — subroutine
- `set_refine3D_trailrec_policy` — subroutine
- `set_refine3D_update_policy` — subroutine
- `set_shift` — subroutine
- `set_state` — subroutine
- `set_stkind` — subroutine
- `set_subgrps` — subroutine
- `setisId_DOM` — subroutine
- `setnamespaceURI` — subroutine
- `setownerDocument` — subroutine
- `sgelsy` — subroutine
- `sgesvd` — subroutine
- `sgetrf` — subroutine
- `sgetri` — subroutine
- `shcloc` — function
- `shell_quote` — function
- `shft` — subroutine
- `shift` — subroutine
- `shuffle_1` — subroutine
- `shuffle_1` — subroutine
- `shuffle_2` — subroutine
- `shuffle_2` — subroutine
- `simple_copy_file` — subroutine
- `simple_end` — subroutine
- `simple_print_git_version` — subroutine
- `simple_print_timer` — subroutine
- `sinc` — function
- `sinc_apod` — function
- `sinc_instr` — function
- `sinhc` — function
- `smat2dmat` — function
- `sniff_folders_SJ` — subroutine
- `snrm2` — function
- `sobel` — subroutine
- `sort` — subroutine
- `sortmeans` — subroutine
- `spaces` — function
- `sparse_eigh` — subroutine
- `spear` — function
- `split_1` — subroutine
- `split_2` — subroutine
- `split_str` — subroutine
- `squared_sampling` — subroutine
- `sqwin_1d_1` — subroutine
- `sqwin_1d_2` — subroutine
- `sqwin_2d_1` — subroutine
- `sqwin_2d_2` — subroutine
- `sqwin_3d_1` — subroutine
- `sqwin_3d_2` — subroutine
- `ssaupd` — subroutine
- `sseupd` — subroutine
- `ssyev` — subroutine
- `ssyevr` — subroutine
- `stage_coarray_jobs` — subroutine
- `stemname` — function
- `stmAddChArr` — subroutine
- `stmAddChMat` — subroutine
- `stmAddChSca` — subroutine
- `stmAddCmplxDpArr` — subroutine
- `stmAddCmplxDpMat` — subroutine
- `stmAddCmplxDpSca` — subroutine
- `stmAddCmplxSpArr` — subroutine
- `stmAddCmplxSpMat` — subroutine
- `stmAddCmplxSpSca` — subroutine
- `stmAddIntArr` — subroutine
- `stmAddIntMat` — subroutine
- `stmAddIntSca` — subroutine
- `stmAddLgArr` — subroutine
- `stmAddLgMat` — subroutine
- `stmAddLgSca` — subroutine
- `stmAddRealDpArr` — subroutine
- `stmAddRealDpMat` — subroutine
- `stmAddRealDpSca` — subroutine
- `stmAddRealSpArr` — subroutine
- `stmAddRealSpMat` — subroutine
- `stmAddRealSpSca` — subroutine
- `str2format_1` — subroutine
- `str2format_2` — subroutine
- `str2ori` — subroutine
- `str2real_1` — function
- `str2real_2` — function
- `str2realdp_local` — subroutine
- `str_complex_dp` — function
- `str_complex_dp_array` — function
- `str_complex_dp_array_fmt` — function
- `str_complex_dp_array_fmt_chk` — function
- `str_complex_dp_array_fmt_len` — function
- `str_complex_dp_array_len` — function
- `str_complex_dp_fmt` — function
- `str_complex_dp_fmt_chk` — function
- `str_complex_dp_fmt_len` — function
- `str_complex_dp_len` — function
- `str_complex_dp_matrix` — function
- `str_complex_dp_matrix_fmt` — function
- `str_complex_dp_matrix_fmt_chk` — function
- `str_complex_dp_matrix_fmt_len` — function
- `str_complex_dp_matrix_len` — function
- `str_complex_sp` — function
- `str_complex_sp_array` — function
- `str_complex_sp_array_fmt` — function
- `str_complex_sp_array_fmt_chk` — function
- `str_complex_sp_array_fmt_len` — function
- `str_complex_sp_array_len` — function
- `str_complex_sp_fmt` — function
- `str_complex_sp_fmt_chk` — function
- `str_complex_sp_fmt_len` — function
- `str_complex_sp_len` — function
- `str_complex_sp_matrix` — function
- `str_complex_sp_matrix_fmt` — function
- `str_complex_sp_matrix_fmt_chk` — function
- `str_complex_sp_matrix_fmt_len` — function
- `str_complex_sp_matrix_len` — function
- `str_integer` — function
- `str_integer_array` — function
- `str_integer_array_fmt` — function
- `str_integer_array_fmt_len` — function
- `str_integer_array_len` — function
- `str_integer_base_len` — function
- `str_integer_fmt` — function
- `str_integer_fmt_len` — function
- `str_integer_len` — function
- `str_integer_matrix` — function
- `str_integer_matrix_fmt` — function
- `str_integer_matrix_fmt_len` — function
- `str_integer_matrix_len` — function
- `str_logical` — function
- `str_logical_array` — function
- `str_logical_array_len` — function
- `str_logical_len` — function
- `str_logical_matrix` — function
- `str_logical_matrix_len` — function
- `str_pad` — function
- `str_real_dp` — function
- `str_real_dp_array` — function
- `str_real_dp_array_fmt` — function
- `str_real_dp_array_fmt_chk` — function
- `str_real_dp_array_fmt_len` — function
- `str_real_dp_array_len` — function
- `str_real_dp_fmt` — function
- `str_real_dp_fmt_chk` — function
- `str_real_dp_fmt_len` — function
- `str_real_dp_len` — function
- `str_real_dp_matrix` — function
- `str_real_dp_matrix_fmt` — function
- `str_real_dp_matrix_fmt_chk` — function
- `str_real_dp_matrix_fmt_len` — function
- `str_real_dp_matrix_len` — function
- `str_real_sp` — function
- `str_real_sp_array` — function
- `str_real_sp_array_fmt` — function
- `str_real_sp_array_fmt_chk` — function
- `str_real_sp_array_fmt_len` — function
- `str_real_sp_array_len` — function
- `str_real_sp_fmt` — function
- `str_real_sp_fmt_chk` — function
- `str_real_sp_fmt_len` — function
- `str_real_sp_len` — function
- `str_real_sp_matrix` — function
- `str_real_sp_matrix_fmt` — function
- `str_real_sp_matrix_fmt_chk` — function
- `str_real_sp_matrix_fmt_len` — function
- `str_real_sp_matrix_len` — function
- `str_string` — function
- `str_string_array` — function
- `str_string_array_len` — function
- `str_string_matrix` — function
- `str_string_matrix_len` — function
- `string` — type
- `submit_coarray_jobs` — subroutine
- `submit_script` — subroutine
- `submit_scripts` — subroutine
- `substr_remove` — function
- `substr_replace` — subroutine
- `svbksb_dp` — subroutine
- `svbksb_sp` — subroutine
- `svd_multifit_dp` — subroutine
- `svd_multifit_sp` — subroutine
- `svd_solve` — subroutine
- `svdcmp_dp` — subroutine
- `svdcmp_sp` — subroutine
- `svdfit_dp` — subroutine
- `svdfit_sp` — subroutine
- `svdvar` — subroutine
- `swap_c` — subroutine
- `swap_cm` — subroutine
- `swap_cv` — subroutine
- `swap_i` — subroutine
- `swap_r` — subroutine
- `swap_rv` — subroutine
- `swap_suffix_1` — function
- `swap_suffix_2` — function
- `swape1e3` — subroutine
- `sym_dists_1` — subroutine
- `sym_dists_2` — subroutine
- `sym_tester` — subroutine
- `symrandomize` — subroutine
- `test_addr` — subroutine
- `test_eigh` — subroutine
- `test_ftiter` — subroutine
- `test_image` — subroutine
- `test_image_local` — subroutine
- `threshold_for_no_peaks` — subroutine
- `threshold_for_npeaks` — subroutine
- `to_char_1` — function
- `to_char_2` — function
- `to_cstring` — function
- `to_fnv1a_hash64` — function
- `to_lower` — function
- `to_static_1` — subroutine
- `to_static_2` — subroutine
- `tokenize_and_add_strings` — subroutine
- `tokenize_to_string_list` — function
- `trace` — function
- `transfer_2Dparams` — subroutine
- `transfer_3Dparams` — subroutine
- `transfmat2inpls` — subroutine
- `transp` — subroutine
- `unique` — subroutine
- `unmemoize_mask_coords` — subroutine
- `unmemoize_powspec_coords` — subroutine
- `update_queue` — subroutine
- `updateNodeLists` — subroutine
- `updatestack` — subroutine
- `updateTextContentLength` — subroutine
- `uppercase` — function
- `upsample_sigma2` — subroutine
- `vabs_dp` — function
- `vabs_sp` — function
- `vector_angle_norm` — function
- `vis_2Dinteger_mat` — subroutine
- `vis_2Dreal_mat` — subroutine
- `vis_3Dinteger_mat` — subroutine
- `vis_3Dreal_mat` — subroutine
- `vox2ang` — function
- `vs_s_concat` — function
- `wait_for_closure` — subroutine
- `watch` — subroutine
- `watchdirs` — subroutine
- `which` — function
- `workout_directory_structure` — subroutine
- `write` — subroutine
- `write2bild` — subroutine
- `write_attributes` — subroutine
- `write_checkpoint` — subroutine
- `write_filetable` — subroutine
- `write_imgarr_1` — subroutine
- `write_imgarr_2` — subroutine
- `write_imgarr_3` — subroutine
- `write_junk_cavgs` — subroutine
- `write_seed_shift_table` — subroutine
- `write_selected_cavgs` — subroutine
- `write_singlelineoftext` — subroutine
- `writeline` — subroutine
- `wxml_error_xf` — subroutine
- `wxml_fatal_xf` — subroutine
- `wxml_warning_xf` — subroutine
- `xml_AddAttlistToDTD` — subroutine
- `xml_AddAttribute_Ch` — subroutine
- `xml_AddCharacters_ch` — subroutine
- `xml_AddComment` — subroutine
- `xml_AddDOCTYPE` — subroutine
- `xml_AddElementToDTD` — subroutine
- `xml_AddEntityReference` — subroutine
- `xml_AddExternalEntity` — subroutine
- `xml_AddInternalEntity` — subroutine
- `xml_AddNewline` — subroutine
- `xml_AddNotation` — subroutine
- `xml_AddParameterEntity` — subroutine
- `xml_AddPEReferenceToDTD` — subroutine
- `xml_AddPseudoAttribute_Ch` — subroutine
- `xml_AddXMLDeclaration` — subroutine
- `xml_AddXMLPI` — subroutine
- `xml_AddXMLStylesheet` — subroutine
- `xml_Close` — subroutine
- `xml_DeclareNamespace` — subroutine
- `xml_EndElement` — subroutine
- `xml_NewElement` — subroutine
- `xml_OpenFile` — subroutine
- `xml_UndeclareNamespace` — subroutine
- `xmlf_GetPretty_print` — function
- `xmlf_name` — function
- `xmlf_opentag` — function
- `xmlf_opentag_len` — function
- `xmlf_SetPretty_print` — subroutine
- `z_scores` — function
- `zheev` — subroutine
- `zposv` — subroutine

---
## Module: pure

Files:
- `main/image/simple_image.f90`
- `main/image/simple_image_access.f90`
- `main/image/simple_image_arith.f90`
- `main/image/simple_image_calc.f90`
- `main/image/simple_image_checks.f90`
- `main/ori/simple_oris.f90`
- `main/ori/simple_oris_dists.f90`
- `main/ori/simple_oris_getters.f90`
- `main/pftc/simple_polarft_access.f90`
- `main/pftc/simple_polarft_calc.f90`
- `main/pftc/simple_polarft_core.f90`
- `main/pftc/simple_polarft_geom.f90`

---
## Module: real

Files:
- `main/image/simple_image.f90`
- `main/image/simple_image_calc.f90`
- `main/nu_filt/simple_nu_filter.f90`
- `main/nu_filt/simple_nu_filter_bank.f90`
- `main/nu_filt/simple_nu_filter_potts.f90`
- `main/nu_filt/simple_nu_filter_state.f90`
- `main/nu_filt/simple_nu_filter_stats.f90`
- `main/ori/simple_oris.f90`
- `main/ori/simple_oris_neig.f90`
- `main/pftc/simple_polarft_calc.f90`
- `main/pftc/simple_polarft_corr.f90`
- `main/project/simple_sp_project.f90`
- `main/project/simple_sp_project_core.f90`

Public symbols:
- `test_oris` — subroutine

---
## Module: recursive

Files:
- `main/image/simple_image.f90`
- `main/image/simple_image_geom.f90`

---
## Module: simple_abinitio2D_controller

Files:
- `main/simple_abinitio2D_controller.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_parameters`

Public symbols:
- `determine_abinitio2D_stages` — subroutine
- `mskdiam2lplimits_cluster2D` — subroutine
- `set_abinitio2D_sampling_policy` — subroutine
- `set_cline_cluster2D_stage` — subroutine

Private symbols:
- `build_cluster2D_stage_cfg` — subroutine
- `emit_cluster2D_stage_cfg` — subroutine
- `set_cluster2D_stage_iteration_policy` — subroutine
- `set_cluster2D_stage_objfun_policy` — subroutine
- `set_cluster2D_stage_phase_policy` — subroutine
- `set_cluster2D_stage_reference_policy` — subroutine
- `set_cluster2D_stage_regular_refs` — subroutine
- `set_cluster2D_stage_search_policy` — subroutine

---
## Module: simple_abinitio_utils

Files:
- `main/simple_abinitio_utils.f90`

Uses:
- `simple_class_frcs`
- `simple_cluster_seed`
- `simple_commanders_api`
- `simple_commanders_rec`
- `simple_commanders_volops`
- `simple_euclid_sigma2`
- `simple_matcher_refvol_utils`
- `simple_parameters`
- `simple_refine3d_fnames`

---
## Module: simple_aff_prop

Files:
- `utils/clustering/simple_aff_prop.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `test_aff_prop` — subroutine

Private symbols:
- `calc_exemplar_mask` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `propagate` — subroutine

---
## Module: simple_ansi_ctrls

Files:
- `defs/simple_ansi_ctrls.f90`

Public symbols:
- `format_str` — function

---
## Module: simple_atoms

Files:
- `main/nano/simple_atoms.f90`

Uses:
- `simple_core_module_api`
- `simple_defs_atoms`
- `simple_image`
- `simple_molecule_data`

Public symbols:
- `atoms` — type
- `test_atoms` — subroutine

Private symbols:
- `assert_close` — subroutine
- `assert_int_eq` — subroutine
- `assert_true` — subroutine
- `assert_vec3_close` — subroutine
- `atom_validate` — subroutine
- `cc_res` — function
- `center_inbox` — subroutine
- `center_pdbcoord` — subroutine
- `cif2mrc` — subroutine
- `cif2pdb` — subroutine
- `convolve` — subroutine
- `copy` — subroutine
- `extract_atom` — subroutine
- `find_masscen` — function
- `fit_bfactors` — subroutine
- `geometry_analysis_pdb` — subroutine
- `get_beta` — function
- `get_coord` — function
- `get_geom_center` — function
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
- `print_atom` — subroutine
- `rotate` — subroutine
- `scan_cif` — subroutine
- `set_atom_corr` — subroutine
- `set_beta` — subroutine
- `set_chain` — subroutine
- `set_coord` — subroutine
- `set_element` — subroutine
- `set_name` — subroutine
- `set_num` — subroutine
- `set_occupancy` — subroutine
- `set_resnum` — subroutine
- `split_tokens` — subroutine
- `strip_digits` — function
- `translate` — subroutine
- `writepdb` — subroutine
- `writepdb_aniso` — subroutine
- `Z_and_radius_from_name` — subroutine

---
## Module: simple_binary_tree

Files:
- `utils/structs/simple_binary_tree.f90`

Uses:
- `simple_error`

Public symbols:
- `binary_tree` — type
- `bt_node` — type
- `deserialize_tree` — subroutine
- `serialize_tree` — subroutine

Private symbols:
- `build_balanced_index_tree` — subroutine
- `build_from_hclust` — subroutine
- `build_range` — function
- `get_node` — function
- `get_root_node` — function
- `int_list` — type
- `kill` — subroutine

---
## Module: simple_binary_tree_tester

Files:
- `utils/structs/simple_binary_tree_tester.f90`

Uses:
- `simple_binary_tree`
- `simple_test_utils`

Public symbols:
- `run_all_tree_tests` — subroutine

Private symbols:
- `assert_idx_valid` — subroutine
- `build_test_tree` — subroutine
- `test_build_balanced_index_tree_deterministic` — subroutine
- `test_build_from_hclust_basic_invariants` — subroutine
- `test_get_node_copy_semantics` — subroutine
- `test_initial_state` — subroutine
- `test_kill_and_rebuild` — subroutine
- `to_str` — function

---
## Module: simple_binoris

Files:
- `fileio/simple_binoris.f90`

Uses:
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

Private symbols:
- `add_segment_1` — subroutine
- `add_segment_2` — subroutine
- `byte_manager4seg_inside_1` — subroutine
- `byte_manager4seg_inside_2` — subroutine
- `clear_segments` — subroutine
- `close` — subroutine
- `get_fromto` — function
- `get_segments_info` — subroutine
- `open` — subroutine
- `open_local` — subroutine
- `read_first_segment_record` — subroutine
- `read_header` — subroutine
- `read_record` — subroutine
- `read_segment_1` — subroutine
- `read_segment_2` — subroutine
- `update_byte_ranges` — subroutine
- `write_header` — subroutine
- `write_segment_1` — subroutine
- `write_segment_2` — subroutine
- `write_segment_inside_1` — subroutine
- `write_segment_inside_2` — subroutine

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
## Module: simple_block_tree

Files:
- `utils/simple_block_tree.f90`

Uses:
- `simple_binary_tree`
- `simple_core_module_api`
- `simple_eulspace_neigh_map`
- `simple_multi_dendro`
- `simple_parameters`

Public symbols:
- `gen_block_index_tree` — function
- `gen_eulspace_block_tree` — function
- `gen_eulspace_block_tree_map` — function
- `gen_multi_block_index_tree` — function
- `gen_single_block_index_tree` — function
- `srch_eul_bl_tree` — subroutine
- `srch_eul_bl_tree_exhaustive` — subroutine
- `srch_eul_bl_tree_prob` — subroutine

Private symbols:
- `push_stack` — subroutine
- `sample_two` — function

---
## Module: simple_block_tree_corr

Files:
- `utils/simple_block_tree_corr.f90`

Uses:
- `simple_aff_prop`
- `simple_core_module_api`
- `simple_corrmat`
- `simple_eulspace_neigh_map`
- `simple_hclust`
- `simple_image`
- `simple_multi_dendro`
- `simple_parameters`

Public symbols:
- `gen_corr_block_tree_aff_prop` — function
- `gen_corr_block_tree_hac` — function
- `gen_eulspace_block_tree_corr` — function

---
## Module: simple_block_tree_io

Files:
- `fileio/simple_block_tree_io.f90`

Uses:
- `simple_error`
- `simple_fileio`
- `simple_multi_dendro`
- `simple_string_utils`

Public symbols:
- `read_block_tree` — subroutine
- `test_io` — subroutine
- `write_block_tree` — subroutine

---
## Module: simple_bspline_smoother

Files:
- `utils/filter/simple_bspline_smoother.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `bspline_smoother` — type
- `test_bspline_smoother` — subroutine
- `test_bspline_smoother_3d` — subroutine

Private symbols:
- `a0xy` — function
- `a2xy` — function
- `b3` — function
- `b3_help` — function
- `fill_b` — subroutine
- `fill_b_3d` — subroutine
- `fill_r` — subroutine
- `fill_r_3d` — subroutine
- `kill_bspline_smoother` — subroutine
- `new_bspline_smoother` — subroutine
- `new_bspline_smoother_img` — subroutine
- `smooth` — subroutine
- `smooth_3d` — subroutine

---
## Module: simple_builder

Files:
- `main/simple_builder.f90`

Uses:
- `simple_binoris_io`
- `simple_class_frcs`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_euclid_sigma2`
- `simple_image`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_projector`
- `simple_reconstructor_eo`
- `simple_sp_project`
- `simple_srchspace_map`
- `simple_ui`

Public symbols:
- `builder` — type

Private symbols:
- `build_general_tbox` — subroutine
- `build_rec_eo_tbox` — subroutine
- `build_spproj` — subroutine
- `build_strategy2D_tbox` — subroutine
- `build_strategy3D_tbox` — subroutine
- `init_params_and_build_general_tbox` — subroutine
- `init_params_and_build_spproj` — subroutine
- `init_params_and_build_strategy2D_tbox` — subroutine
- `init_params_and_build_strategy3D_tbox` — subroutine
- `init_params_spproj_tbox2D` — subroutine
- `kill_general_tbox` — subroutine
- `kill_rec_eo_tbox` — subroutine
- `kill_strategy2D_tbox` — subroutine
- `kill_strategy3D_tbox` — subroutine
- `set_field_ptr` — subroutine

---
## Module: simple_butterworth

Files:
- `utils/filter/simple_butterworth.f90`

Uses:
- `simple_core_module_api`

---
## Module: simple_calc_pspec_strategy

Files:
- `main/strategies/parallelization/simple_calc_pspec_strategy.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_commanders_euclid_distr`
- `simple_core_module_api`
- `simple_image`
- `simple_matcher_ptcl_io`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_ran_tabu`
- `simple_sigma2_binfile`

Public symbols:
- `calc_pspec_distr_strategy` — type
- `calc_pspec_inmem_strategy` — type
- `calc_pspec_strategy` — type
- `create_calc_pspec_strategy` — function

Private symbols:
- `cleanup_calc_pspec_outputs` — subroutine
- `cleanup_interface` — subroutine
- `compute_pspec_channel` — subroutine
- `distr_cleanup` — subroutine
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `ensure_calc_pspec_eo_partition` — subroutine
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `read_sigma2_bootstrap_selection` — subroutine
- `sanity_check_calc_pspec_input` — subroutine
- `write_sigma2_bootstrap_selection` — subroutine

---
## Module: simple_cavg_quality_analysis

Files:
- `main/cavg_quality/simple_cavg_quality_analysis.f90`

Uses:
- `simple_cavg_quality_feats`
- `simple_cavg_quality_model`
- `simple_cavg_quality_relations`
- `simple_cavg_quality_stats`
- `simple_cavg_quality_types`
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_oris`
- `simple_parameters`

Public symbols:
- `evaluate_cavg_quality` — subroutine
- `evaluate_cavg_quality_hard_reject` — subroutine
- `prepare_cavg_quality` — subroutine
- `write_cavg_quality_feature_table` — subroutine
- `write_cavg_quality_training_table` — subroutine

Private symbols:
- `calc_best_thresholds` — subroutine
- `eval_threshold` — subroutine
- `write_analysis_class_header` — subroutine
- `write_analysis_class_row` — subroutine
- `write_feature_table_class_row` — subroutine
- `write_model_as_analysis_comments` — subroutine
- `write_one_threshold` — subroutine
- `write_threshold_scan_comments` — subroutine

---
## Module: simple_cavg_quality_feats

Files:
- `main/cavg_quality/simple_cavg_quality_feats.f90`

Uses:
- `simple_cavg_quality_types`
- `simple_defs`
- `simple_error`
- `simple_image`
- `simple_image_bin`
- `simple_oris`
- `simple_segmentation`
- `simple_stat`

Public symbols:
- `cavg_quality_feature_description` — function
- `cavg_quality_feature_direction` — function
- `cavg_quality_feature_family` — function
- `cavg_quality_feature_name` — function
- `extract_cavg_quality_features` — subroutine
- `normalize_cavg_quality_features` — subroutine
- `write_cavg_quality_feature_inventory` — subroutine

Private symbols:
- `measure_cavg_foreground_geometry` — subroutine
- `measure_cavg_image_metrics` — subroutine
- `validate_quality_context` — subroutine

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
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`
- `simple_srch_sort_loc`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `evaluate_cavg_quality_model` — subroutine
- `evaluate_cavg_quality_result` — subroutine
- `learn_cavg_quality_model` — subroutine

Private symbols:
- `abinitio_learn_base_spec` — function
- `append_feature_family` — subroutine
- `apply_feature_policy` — subroutine
- `best_score_threshold_balacc` — subroutine
- `best_tie_minsep_span` — subroutine
- `build_logistic_problem` — subroutine
- `build_policy_caches` — subroutine
- `calc_suggested_training_weights` — subroutine
- `cavg_quality_logistic_problem` — type
- `classify_training_dataset` — subroutine
- `classify_training_dataset_cached_detail` — subroutine
- `classify_training_dataset_detail` — subroutine
- `collect_learn_diagnostics` — subroutine
- `consider_model_candidate` — subroutine
- `dataset_learn_role_name` — function
- `dataset_short_name` — function
- `evaluate_candidate_spec` — subroutine
- `evaluate_policy_grid` — subroutine
- `feature_auc_summary` — subroutine
- `feature_policy_indices` — subroutine
- `feature_policy_mask` — subroutine
- `feature_policy_name` — function
- `fit_logistic_problem` — subroutine
- `kill_logistic_problem` — subroutine
- `kill_policy_caches` — subroutine
- `kill_training_datasets` — subroutine
- `learn_cavg_quality_pairwise_logistic_model` — subroutine
- `load_quality_training_datasets` — subroutine
- `lodo_feature_policy_weights` — subroutine
- `logistic_cost_wrapper` — function
- `logistic_dataset_class_weights` — subroutine
- `logistic_fdf_wrapper` — subroutine
- `logistic_gradient_wrapper` — subroutine
- `logistic_problem_fdf` — subroutine
- `logistic_solution_to_model` — subroutine
- `map_analysis_columns` — subroutine
- `min_accept_policy_level` — function
- `parse_training_metadata_line` — subroutine
- `policy_level` — function
- `read_quality_training_dataset` — subroutine
- `record_top_model_candidate` — subroutine
- `require_analysis_columns` — subroutine
- `require_relational_training_datasets` — subroutine
- `rescue_policy_level` — function
- `score_dataset_feature_policy` — subroutine
- `select_preferred_best_tie` — subroutine
- `sort_real_prefix_ascending` — subroutine
- `validate_quality_context` — subroutine
- `validate_training_dataset_schemas` — subroutine
- `write_bool_grid_diagnostic` — subroutine
- `write_candidate_row` — subroutine
- `write_candidate_table_header` — subroutine
- `write_cavg_quality_evaluate_report` — subroutine
- `write_cavg_quality_learn_report` — subroutine
- `write_cavg_quality_logistic_learn_report` — subroutine
- `write_dataset_metric_table` — subroutine
- `write_evaluate_diagnostic` — subroutine
- `write_evaluate_diagnostics` — subroutine
- `write_feature_drop_diagnostics` — subroutine
- `write_feature_policy_grid` — subroutine
- `write_feature_policy_lodo_row` — subroutine
- `write_feature_policy_screen` — subroutine
- `write_feature_screen_diagnostics` — subroutine
- `write_feature_signal_diagnostics` — subroutine
- `write_fixed_model_summary` — subroutine
- `write_grid_position_diagnostic` — subroutine
- `write_learn_search_diagnostics` — subroutine
- `write_logical_list` — subroutine
- `write_logistic_coefficient_table` — subroutine
- `write_minsep_diagnostic` — subroutine
- `write_otsu_ablation_diagnostics` — subroutine
- `write_policy_parameter_diagnostics` — subroutine
- `write_real_list` — subroutine
- `write_search_diagnostic` — subroutine

---
## Module: simple_cavg_quality_model

Files:
- `main/cavg_quality/simple_cavg_quality_model.f90`

Uses:
- `simple_cavg_quality_feats`
- `simple_cavg_quality_stats`
- `simple_cavg_quality_types`
- `simple_clustering_utils`
- `simple_defs`
- `simple_error`
- `simple_srch_sort_loc`
- `simple_string_utils`

Public symbols:
- `apply_cached_decision_to_quality` — subroutine
- `build_classify_cache` — subroutine
- `cached_decision_confusion` — subroutine
- `cavg_quality_classify_cache` — type
- `cavg_quality_model` — type
- `kill_classify_cache` — subroutine
- `write_cavg_quality_model_builtin_code` — subroutine

Private symbols:
- `accept_fit_as_single_cluster` — subroutine
- `accumulate_accept_all_confusion` — subroutine
- `apply_pairwise_logistic` — subroutine
- `builtin_names` — function
- `builtin_spec` — function
- `cavg_quality_cached_decision` — type
- `choose_tied_good_label` — subroutine
- `chunk100mics_model_spec` — function
- `classify` — subroutine
- `forced_reason_name` — function
- `get_spec` — function
- `init_preset` — subroutine
- `init_spec` — subroutine
- `kill_model` — subroutine
- `mark_threshold_accepts_all` — subroutine
- `normalize` — subroutine
- `otsu_score_threshold` — subroutine
- `parse_feature_weight_values` — subroutine
- `parse_model_key_value` — subroutine
- `parse_real_values` — subroutine
- `prepare_cached_decision` — subroutine
- `read_feature_weights` — subroutine
- `read_interaction_terms` — subroutine
- `read_model` — subroutine
- `read_real_values_keyed` — subroutine
- `set_pairwise_interactions_for_feature_count` — subroutine
- `sieve_model_spec` — function
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
- `simple_test_utils`
- `simple_type_defs`

Public symbols:
- `cavg_quality_relation_analysis` — type
- `test_cavg_quality_relations` — subroutine

Private symbols:
- `calculate_promoted_feature` — subroutine
- `calculate_relation_analysis` — subroutine
- `kill_relation_analysis` — subroutine
- `normalize_relation_feature` — subroutine
- `sorted_cc_neighbours` — subroutine

---
## Module: simple_cavg_quality_stats

Files:
- `main/cavg_quality/simple_cavg_quality_stats.f90`

Uses:
- `simple_cavg_quality_types`
- `simple_error`
- `simple_stat`

Public symbols:
- `calc_binary_metrics` — subroutine
- `calc_confusion` — subroutine
- `normalize_quality_dmat` — subroutine

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
- `resample_filter` — function

Private symbols:
- `avg_frc_getter` — subroutine
- `bound_res` — subroutine
- `crop` — subroutine
- `estimate_find_for_eoavg` — function
- `estimate_lp_for_align` — function
- `estimate_res` — subroutine
- `frc_getter` — subroutine
- `get_frc` — function
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
- `class_samples_same` — function
- `deallocate_class_samples` — subroutine
- `print_class_sample` — subroutine
- `read_class_samples` — subroutine
- `write_class_samples` — subroutine

Private symbols:
- `serialize_class_sample` — function
- `unserialize_class_sample` — function

---
## Module: simple_classaverager

Files:
- `main/class/simple_classaverager.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_ctf`
- `simple_euclid_sigma2`
- `simple_fftw3`
- `simple_ftiter`
- `simple_image`
- `simple_imgfile`
- `simple_memoize_ft_maps`
- `simple_parameters`

Public symbols:
- `fourier_2d_accumulator` — type

Private symbols:
- `matrix_ptrs` — type
- `stack` — type

---
## Module: simple_cls_split_strategy

Files:
- `main/strategies/parallelization/simple_cls_split_strategy.f90`

Uses:
- `simple_builder`
- `simple_classaverager`
- `simple_clustering_utils`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_default_clines`
- `simple_defs_fname`
- `simple_diff_map_graphs`
- `simple_diffusion_maps`
- `simple_exec_helpers`
- `simple_image`
- `simple_image_msk`
- `simple_imgarr_utils`
- `simple_imgproc`
- `simple_kpca_svd`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_sp_project`
- `simple_srch_sort_loc`
- `simple_srchspace_map2d_io`

Public symbols:
- `cls_split_master_strategy` — type
- `cls_split_shmem_strategy` — type
- `cls_split_strategy` — type
- `cls_split_worker_strategy` — type
- `create_cls_split_strategy` — function

Private symbols:
- `apply_defaults` — subroutine
- `apply_split_project_updates` — subroutine
- `cleanup_interface` — subroutine
- `collect_split_classes` — subroutine
- `count_data_lines` — subroutine
- `determine_phase_flip` — subroutine
- `determine_split_label` — subroutine
- `estimate_icm_rank_stats` — subroutine
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_common` — subroutine
- `init_interface` — subroutine
- `make_split_embedding` — subroutine
- `master_cleanup` — subroutine
- `master_execute` — subroutine
- `master_finalize_run` — subroutine
- `master_initialize` — subroutine
- `merge_worker_outputs` — subroutine
- `prepare_class_partitions` — subroutine
- `read_int_file` — subroutine
- `read_part_map` — subroutine
- `run_icm_rank_seed` — subroutine
- `run_local_split` — subroutine
- `sanitize_distance_matrix` — subroutine
- `sanitize_embedding_coords` — subroutine
- `select_kmedoids_by_silhouette` — subroutine
- `shmem_cleanup` — subroutine
- `shmem_execute` — subroutine
- `shmem_finalize_run` — subroutine
- `shmem_initialize` — subroutine
- `sort_combined_maps` — subroutine
- `sort_order_by_weight_desc` — subroutine
- `split_one_parent_class` — subroutine
- `trim_split_embedding` — subroutine
- `update_icm_rank_site` — subroutine
- `validate_cls_split` — subroutine
- `worker_cleanup` — subroutine
- `worker_execute` — subroutine
- `worker_finalize_run` — subroutine
- `worker_initialize` — subroutine

---
## Module: simple_cluster2D_strategy

Files:
- `main/strategies/parallelization/simple_cluster2D_strategy.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_commanders_euclid`
- `simple_commanders_imgops`
- `simple_commanders_mkcavgs`
- `simple_commanders_prob`
- `simple_convergence`
- `simple_core_module_api`
- `simple_euclid_sigma2`
- `simple_exec_helpers`
- `simple_gui_utils`
- `simple_matcher_smpl_and_lplims`
- `simple_parameters`
- `simple_procimgstk`
- `simple_progress`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_starproject`
- `simple_strategy2d_matcher`
- `simple_stream_utils`

Public symbols:
- `cluster2D_distr_strategy` — type
- `cluster2D_inmem_strategy` — type
- `cluster2D_strategy` — type
- `create_cluster2D_strategy` — function

Private symbols:
- `cleanup_distributed_iteration_artifacts` — subroutine
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_execute_iteration` — subroutine
- `distr_finalize_iteration` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `exec_iter_interface` — subroutine
- `execute_make_cavgs` — subroutine
- `finalize_iter_interface` — subroutine
- `finalize_run_interface` — subroutine
- `gen_jpeg` — subroutine
- `init_cluster2D_refs` — subroutine
- `init_interface` — subroutine
- `init_standard_refs` — subroutine
- `init_tseries_refs` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_execute_iteration` — subroutine
- `inmem_finalize_iteration` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `run_distributed_cavg_assembly` — subroutine

---
## Module: simple_cluster_seed

Files:
- `utils/clustering/simple_cluster_seed.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `gen_labelling` — subroutine

Private symbols:
- `draw` — subroutine
- `draw_ranked` — subroutine
- `draw_ranked_uniform` — subroutine
- `draw_squared` — subroutine
- `draw_squared_uniform` — subroutine
- `draw_squared_uniform_projdir` — subroutine
- `draw_uniform` — subroutine
- `draw_uniform_corr` — subroutine

---
## Module: simple_clustering_utils

Files:
- `utils/clustering/simple_clustering_utils.f90`

Uses:
- `simple_aff_prop`
- `simple_core_module_api`
- `simple_hclust`
- `simple_kmedoids`
- `simple_stat`

Public symbols:
- `aggregate` — subroutine
- `cluster_dmat` — subroutine
- `labels2smat` — function

Private symbols:
- `merge_if_necessary` — subroutine

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
- `cmdline_err` — subroutine

Private symbols:
- `assign` — subroutine
- `check` — subroutine
- `checkvar` — subroutine
- `copy` — subroutine
- `defined` — function
- `delete` — subroutine
- `gen_job_descr` — subroutine
- `get_carg` — function
- `get_keys` — function
- `get_rarg` — function
- `kill` — subroutine
- `lookup` — function
- `parse` — subroutine
- `parse_command_line_value` — subroutine
- `parse_oldschool` — subroutine
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
- `writeline` — subroutine

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
- `test_read_parsing` — subroutine
- `test_set_and_get_char_and_string` — subroutine
- `test_set_and_get_numeric` — subroutine

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
## Module: simple_commanders_abinitio

Files:
- `main/commanders/simple/simple_commanders_abinitio.f90`

Uses:
- `simple_abinitio_utils`
- `simple_cavg_quality_analysis`
- `simple_cavg_quality_model`
- `simple_cavg_quality_types`
- `simple_cluster_seed`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_project_core`
- `simple_commanders_rec`
- `simple_commanders_refine3d`
- `simple_commanders_reproject`
- `simple_estimate_ssnr`
- `simple_imgarr_utils`
- `simple_procimgstk`
- `simple_qsys_funs`
- `simple_refine3d_fnames`

Public symbols:
- `commander_abinitio3D` — type
- `commander_abinitio3D_cavgs` — type
- `commander_abinitio3D_cavgs_reject` — type

Private symbols:
- `annotate_project` — subroutine
- `assign_good_bad_states` — subroutine
- `average_consensus_volume` — subroutine
- `best_state_mapping` — subroutine
- `build_consensus` — subroutine
- `build_consensus_volume` — subroutine
- `build_selection_states` — subroutine
- `clean_ptcl3D_sampling` — subroutine
- `cleanup` — subroutine
- `configure_cavgs_distributed_clines` — subroutine
- `conv_eo` — subroutine
- `conv_eo_states` — subroutine
- `dock_consensus_volumes` — subroutine
- `ensure_docked_multistate_particle_assignments` — subroutine
- `ensure_multistate_particle_assignments` — subroutine
- `evaluate_quality` — subroutine
- `exec_abinitio3D` — subroutine
- `exec_abinitio3D_cavgs` — subroutine
- `exec_abinitio3D_cavgs_reject` — subroutine
- `ini3D_from_cavgs` — subroutine
- `init_quality_model` — subroutine
- `map_state_correspondence` — subroutine
- `prepare_state_continue_project` — subroutine
- `rank_cavgs` — subroutine
- `read_multistate_assignment_coverage` — subroutine
- `read_restart_labels` — subroutine
- `reset_ptcl3D_from_ptcl2D_selection` — subroutine
- `rndstart` — subroutine
- `run_multistate_missing_update` — subroutine
- `sanitize_restart_labels` — subroutine
- `select_restart_consensus_volumes` — subroutine
- `strip_distributed_child` — subroutine
- `submit_restarts` — subroutine
- `sync_distributed_child` — subroutine
- `validate_cavg_ini_ext_states` — subroutine
- `write_consensus_report` — subroutine
- `write_consensus_volume_report` — subroutine
- `write_rejection_stack` — subroutine
- `write_selection_outputs` — subroutine

---
## Module: simple_commanders_abinitio2D

Files:
- `main/commanders/simple/simple_commanders_abinitio2D.f90`

Uses:
- `simple_abinitio2d_controller`
- `simple_class_frcs`
- `simple_classaverager`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_cluster2d`
- `simple_commanders_imgops`
- `simple_commanders_volops`
- `simple_procimgstk`
- `simple_timer`

Public symbols:
- `commander_abinitio2D` — type

Private symbols:
- `exec_abinitio2D` — subroutine
- `execute_cluster2D` — subroutine
- `execute_terminal_prob_pass` — subroutine
- `gen_final_cavgs` — subroutine
- `inirefs` — subroutine
- `output_stats` — subroutine
- `prep_command_lines` — subroutine
- `set_dims` — subroutine
- `set_lplims` — subroutine
- `set_sampling` — subroutine
- `write_abinitio_benchmark` — subroutine

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
- `simple_euclid_sigma2`
- `simple_exec_helpers`
- `simple_image`
- `simple_image_bin`
- `simple_image_msk`
- `simple_nice`
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
- `simple_atoms`
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
- `commander_detect_atoms` — type
- `commander_map2model_fsc` — type
- `commander_map_validate` — type
- `commander_model_validate` — type
- `commander_model_validate_eo` — type
- `commander_pdb2mrc` — type
- `common_atoms` — type
- `compute_stats` — subroutine
- `exec_atoms_register` — subroutine
- `exec_atoms_rmsd` — subroutine
- `exec_atoms_stats` — subroutine
- `exec_cif2mrc` — subroutine
- `exec_cif2pdb` — subroutine
- `exec_conv_atom_denoise` — subroutine
- `exec_core_atoms_analysis` — subroutine
- `exec_crys_score` — subroutine
- `exec_detect_atoms` — subroutine
- `exec_map2model_fsc` — subroutine
- `exec_map_validate` — subroutine
- `exec_model_validate` — subroutine
- `exec_model_validate_eo` — subroutine
- `exec_pdb2mrc` — subroutine
- `sort_cluster` — subroutine

---
## Module: simple_commanders_cavgs

Files:
- `main/commanders/simple/simple_commanders_cavgs.f90`

Uses:
- `simple_cavg_quality_analysis`
- `simple_cavg_quality_learn`
- `simple_cavg_quality_model`
- `simple_cavg_quality_relations`
- `simple_cavg_quality_types`
- `simple_clustering_utils`
- `simple_commanders_api`
- `simple_corrmat`
- `simple_imgarr_utils`
- `simple_projfile_utils`
- `simple_srch_sort_loc`
- `simple_strategy2d_utils`
- `simple_string_utils`

Public symbols:
- `annotate_project` — subroutine
- `build_temp_project_from_stack` — subroutine
- `commander_cluster_cavgs` — type
- `commander_cluster_cavgs_selection` — type
- `commander_map_cavgs_selection` — type
- `commander_map_cavgs_states` — type
- `commander_match_cavgs` — type
- `commander_model_cavgs_rejection` — type
- `commander_rank_cavgs` — type
- `commander_select_clusters` — type
- `commander_shape_rank_cavgs` — type
- `exec_cluster_cavgs` — subroutine
- `exec_cluster_cavgs_selection` — subroutine
- `exec_map_cavgs_selection` — subroutine
- `exec_map_cavgs_states` — subroutine
- `exec_match_cavgs` — subroutine
- `exec_model_cavgs_rejection` — subroutine
- `exec_rank_cavgs` — subroutine
- `exec_select_clusters` — subroutine
- `exec_shape_rank_cavgs` — subroutine
- `find_label_state` — function
- `p1_lt_p2` — function
- `prep_stack4clust` — subroutine
- `resolve_quality_context` — function
- `update_project` — subroutine
- `write_hard_gate_stack` — subroutine
- `write_quality_stack` — subroutine
- `write_ranked_quality_stack` — subroutine

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
## Module: simple_commanders_cluster2D

Files:
- `main/commanders/simple/simple_commanders_cluster2D.f90`

Uses:
- `simple_classaverager`
- `simple_cluster2d_strategy`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_euclid`
- `simple_commanders_imgops`
- `simple_commanders_mkcavgs`
- `simple_gui_utils`
- `simple_imgproc`
- `simple_kpca_svd`
- `simple_matcher_smpl_and_lplims`
- `simple_pca`
- `simple_pca_svd`
- `simple_pftc_srch_api`
- `simple_ppca`
- `simple_procimgstk`
- `simple_progress`
- `simple_strategy2d_matcher`

Public symbols:
- `commander_cluster2D` — type
- `commander_cluster2D_distr_worker` — type
- `commander_ppca_denoise_classes` — type
- `exec_cluster2D` — subroutine
- `exec_cluster2D_distr_worker` — subroutine
- `exec_ppca_denoise_classes` — subroutine
- `log_ppca_rank_scan` — subroutine

---
## Module: simple_commanders_denoise

Files:
- `main/commanders/simple/simple_commanders_denoise.f90`

Uses:
- `simple_cls_split_strategy`
- `simple_commanders_api`
- `simple_commanders_mkcavgs`
- `simple_denoise_project_strategy`
- `simple_ori`
- `simple_oris`

Public symbols:
- `build_den2raw_map` — subroutine
- `commander_cls_split` — type
- `commander_denoise_project` — type
- `commander_map_params_from_den` — type
- `copy_assignment_keys` — subroutine
- `exec_cls_split` — subroutine
- `exec_denoise_project` — subroutine
- `exec_map_params_from_den` — subroutine
- `validate_den2raw_map` — subroutine

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
- `simple_sigma2_binfile`

Public symbols:
- `commander_calc_group_sigmas` — type
- `commander_calc_pspec` — type
- `consolidate_channel` — subroutine
- `exec_calc_group_sigmas` — subroutine
- `exec_calc_pspec` — subroutine

---
## Module: simple_commanders_euclid_distr

Files:
- `main/commanders/simple/simple_commanders_euclid_distr.f90`

Uses:
- `simple_commanders_api`
- `simple_sigma2_binfile`

Public symbols:
- `commander_calc_pspec_assemble` — type
- `exec_calc_pspec_assemble` — subroutine
- `remove_negative_sigmas` — subroutine

---
## Module: simple_commanders_flex_analysis

Files:
- `main/commanders/simple/simple_commanders_flex_analysis.f90`

Uses:
- `simple_commanders_api`
- `simple_flex_analysis_strategy`

Public symbols:
- `commander_flex_analysis` — type
- `exec_flex_analysis` — subroutine

---
## Module: simple_commanders_imgops

Files:
- `main/commanders/simple/simple_commanders_imgops.f90`

Uses:
- `simple_bspline_smoother`
- `simple_commanders_api`
- `simple_diff_map_denoise`
- `simple_diff_map_graphs`
- `simple_imgarr_utils`
- `simple_imgproc`
- `simple_kpca_svd`
- `simple_pca`
- `simple_pca_svd`
- `simple_ppca`
- `simple_procimgstk`
- `simple_segmentation`

Public symbols:
- `commander_binarize` — type
- `commander_edge_detect` — type
- `commander_filter` — type
- `commander_normalize` — type
- `commander_ppca_denoise` — type
- `commander_scale` — type
- `doit` — subroutine
- `doit` — subroutine
- `exec_binarize` — subroutine
- `exec_edge_detect` — subroutine
- `exec_filter` — subroutine
- `exec_normalize` — subroutine
- `exec_ppca_denoise` — subroutine
- `exec_scale` — subroutine

---
## Module: simple_commanders_imgproc

Files:
- `main/commanders/simple/simple_commanders_imgproc.f90`

Uses:
- `simple_commanders_api`
- `simple_ctf`
- `simple_matcher_ptcl_io`
- `simple_memoize_ft_maps`
- `simple_ori`
- `simple_procimgstk`
- `simple_segmentation`

Public symbols:
- `commander_ctf_phaseflip` — type
- `commander_ctfops` — type
- `commander_estimate_diam` — type
- `exec_ctf_phaseflip` — subroutine
- `exec_ctfops` — subroutine
- `exec_estimate_diam` — subroutine

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
## Module: simple_commanders_misc

Files:
- `main/commanders/simple/simple_commanders_misc.f90`

Uses:
- `simple_class_frcs`
- `simple_commanders_api`
- `simple_fsc`
- `simple_micrograph_generator`
- `simple_motion_correct_utils`
- `simple_simple_volinterp`
- `simple_starproject`

Public symbols:
- `commander_fractionate_movies` — type
- `commander_fractionate_movies_distr` — type
- `commander_kstest` — type
- `commander_mkdir` — type
- `commander_pearsn` — type
- `commander_print_dose_weights` — type
- `commander_print_fsc` — type
- `commander_print_magic_boxes` — type
- `exec_fractionate_movies` — subroutine
- `exec_fractionate_movies_distr` — subroutine
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
- `simple_classaverager`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_make_cavgs_strategy`
- `simple_parameters`
- `simple_pftc_srch_api`
- `simple_timer`

Public symbols:
- `append_bootstrap_row` — subroutine
- `append_skipped_class` — subroutine
- `build_bootstrap_project` — subroutine
- `clear_parent_sample` — subroutine
- `collect_bootstrap_classes` — subroutine
- `commander_bootstrap_cavgs` — type
- `commander_cavgassemble` — type
- `commander_make_cavgs` — type
- `commander_make_cavgs_distr` — type
- `commander_unbootstrap_cavgs` — type
- `commander_write_classes` — type
- `derive_bootstrap_output_names` — subroutine
- `exec_bootstrap_cavgs` — subroutine
- `exec_cavgassemble` — subroutine
- `exec_make_cavgs` — subroutine
- `exec_make_cavgs_distr` — subroutine
- `exec_unbootstrap_cavgs` — subroutine
- `exec_write_classes` — subroutine
- `make_cavgs_exec_cavgassemble` — subroutine
- `report_bootstrap_class_accounting` — subroutine
- `run_make_cavgs_workflow` — subroutine
- `select_bootstrap_anchor` — subroutine
- `set_bootstrap_class_row` — subroutine
- `shuffle_ints` — subroutine
- `validate_evenodd_present` — subroutine
- `write_bootstrap_outputs` — subroutine

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
- `commander_oriconsensus` — type
- `commander_oriops` — type
- `commander_oristats` — type
- `commander_rotmats2oris` — type
- `commander_vizoris` — type
- `exec_make_oris` — subroutine
- `exec_oriconsensus` — subroutine
- `exec_oriops` — subroutine
- `exec_oristats` — subroutine
- `exec_rotmats2oris` — subroutine
- `exec_vizoris` — subroutine
- `generate_perm` — subroutine
- `read_state_order` — subroutine

---
## Module: simple_commanders_pcg_recon

Files:
- `main/commanders/simple/simple_commanders_pcg_recon.f90`

Uses:
- `simple_commanders_api`
- `simple_matcher_ptcl_io`
- `simple_math_ft`
- `simple_pcg_reconstruction`
- `simple_sigma2_files`

Public symbols:
- `commander_reconstruct3D_pcg` — type

Private symbols:
- `exec_reconstruct3D_pcg` — subroutine

---
## Module: simple_commanders_pick

Files:
- `main/commanders/simple/simple_commanders_pick.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_api`
- `simple_core_module_api`
- `simple_extract_strategy`
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
## Module: simple_commanders_preprocess

Files:
- `main/commanders/simple/simple_commanders_preprocess.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_api`
- `simple_core_module_api`
- `simple_ctf_estimate_strategy`
- `simple_gen_pspecs_and_thumbs_strategy`
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
## Module: simple_commanders_prob

Files:
- `main/commanders/simple/simple_commanders_prob.f90`

Uses:
- `simple_builder`
- `simple_classaverager`
- `simple_commanders_api`
- `simple_eul_prob_tab`
- `simple_eul_prob_tab2d`
- `simple_eul_prob_tab_neigh`
- `simple_imgarr_utils`
- `simple_matcher_2dprep`
- `simple_matcher_pftc_prep`
- `simple_matcher_ptcl_batch`
- `simple_matcher_refvol_utils`
- `simple_matcher_smpl_and_lplims`
- `simple_pftc_srch_api`
- `simple_strategy2d_matcher`

Public symbols:
- `cleanup_prob_align_outputs` — subroutine
- `cleanup_prob_neigh_workspace` — subroutine
- `commander_prob_align` — type
- `commander_prob_align2D` — type
- `commander_prob_align_neigh` — type
- `commander_prob_tab` — type
- `commander_prob_tab2D` — type
- `commander_prob_tab_neigh` — type
- `exec_prob_align` — subroutine
- `exec_prob_align2D` — subroutine
- `exec_prob_align_neigh` — subroutine
- `exec_prob_tab` — subroutine
- `exec_prob_tab2D` — subroutine
- `exec_prob_tab_neigh` — subroutine
- `prepare_prob_neigh_workspace` — subroutine
- `run_prob_tab_neigh_batch` — subroutine

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
- `simple_imgarr_utils`
- `simple_projfile_utils`
- `simple_strategy2d_utils`
- `simple_stream_communicator`

Public symbols:
- `apply_consensus` — subroutine
- `best_state_mapping` — subroutine
- `build_agreement` — subroutine
- `build_consensus` — subroutine
- `cleanup` — subroutine
- `commander_aggregate_chunks` — type
- `commander_concatenate_projects` — type
- `commander_extract_subproj` — type
- `commander_merge_projects` — type
- `commander_new_project` — type
- `commander_print_project_field` — type
- `commander_print_project_info` — type
- `commander_print_project_vals` — type
- `commander_ptcl3D_state_consensus` — type
- `commander_replace_project_field` — type
- `commander_selection` — type
- `commander_update_project` — type
- `commander_validate_projfile` — type
- `exec_aggregate_chunks` — subroutine
- `exec_concatenate_projects` — subroutine
- `exec_extract_subproj` — subroutine
- `exec_merge_projects` — subroutine
- `exec_new_project` — subroutine
- `exec_print_project_field` — subroutine
- `exec_print_project_info` — subroutine
- `exec_print_project_vals` — subroutine
- `exec_ptcl3D_state_consensus` — subroutine
- `exec_replace_project_field` — subroutine
- `exec_selection` — subroutine
- `exec_update_project` — subroutine
- `exec_validate_projfile` — subroutine
- `generate_mapping` — subroutine
- `map_state_correspondence` — subroutine
- `sanitize_source_labels` — subroutine
- `write_consensus_report` — subroutine

---
## Module: simple_commanders_project_mov

Files:
- `main/commanders/simple/simple_commanders_project_mov.f90`

Uses:
- `simple_commanders_api`
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
- `simple_starproject`
- `simple_stream_watcher`

Public symbols:
- `commander_import_boxes` — type
- `commander_import_particles` — type
- `commander_prune_project` — type
- `commander_prune_project_distr` — type
- `commander_reimport_particles` — type
- `commander_scale_project` — type
- `commander_split_stack` — type
- `commander_zero_project_shifts` — type
- `exec_import_boxes` — subroutine
- `exec_import_particles` — subroutine
- `exec_prune_project` — subroutine
- `exec_prune_project_distr` — subroutine
- `exec_reimport_particles` — subroutine
- `exec_scale_project` — subroutine
- `exec_split_stack` — subroutine
- `exec_zero_project_shifts` — subroutine
- `validate_denoised_stack_pair` — subroutine

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
- `simple_rec3d_strategy`
- `simple_refine3d_fnames`

Public symbols:
- `accumulate_bootstrap_eo_counts` — subroutine
- `calc_bootstrap_state_scales` — subroutine
- `commander_bootstrap_rec3D` — type
- `commander_rec3D` — type
- `commander_rec3D_worker` — type
- `condition_bootstrap_sigma2` — subroutine
- `exec_bootstrap_rec3D` — subroutine
- `exec_random_rec` — subroutine
- `exec_rec3D` — subroutine
- `exec_rec3D_distr_worker` — subroutine
- `prepare_bootstrap_rec_cline` — subroutine
- `random_rec_commander` — type
- `register_bootstrap_rec_outputs` — subroutine
- `warn_for_forced_bootstrap_overrides` — subroutine
- `write_bootstrap_sigma2_from_halfmaps` — subroutine

---
## Module: simple_commanders_rec_distr

Files:
- `main/commanders/simple/simple_commanders_rec_distr.f90`

Uses:
- `simple_commanders_api`
- `simple_fsc`
- `simple_gridding`
- `simple_nu_filter`
- `simple_reconstructor_eo`
- `simple_refine3d_fnames`
- `simple_vol_pproc_policy`

Public symbols:
- `commander_volassemble` — type

Private symbols:
- `assemble_state` — subroutine
- `build_nonuniform_mask` — subroutine
- `capture_nonuniform_source_halves` — subroutine
- `cleanup_context` — subroutine
- `cleanup_nonuniform_state` — subroutine
- `cleanup_nu_aux_images` — subroutine
- `cleanup_restore_state` — subroutine
- `collect_restore_timings` — subroutine
- `determine_trailing_update_fraction` — subroutine
- `exec_volassemble` — subroutine
- `initialize_bench_timers` — subroutine
- `initialize_context` — subroutine
- `log_nonuniform_filter_stats` — subroutine
- `log_nu_alignment_lowpass_summary` — subroutine
- `nu_highres_steps_fname` — function
- `postprocess_state` — subroutine
- `prepare_grid_correction` — subroutine
- `read_previous_halfmaps` — subroutine
- `record_nu_alignment_lowpass_limit` — subroutine
- `reduce_partials` — subroutine
- `refine_nonuniform_filter_bank` — subroutine
- `refresh_state_populations` — subroutine
- `release_nonuniform_aux_inputs` — subroutine
- `release_nonuniform_base_inputs` — subroutine
- `resolution_outfile_fbody` — function
- `resolve_fsc_txt_fname` — function
- `restore_eos_and_write_fsc` — subroutine
- `restore_merged_volume` — subroutine
- `restore_state_from_parts` — subroutine
- `restore_timings_t` — type
- `run_state_automask` — subroutine
- `run_state_nonuniform_filter` — subroutine
- `set_state_filenames` — subroutine
- `setup_nonuniform_filter` — subroutine
- `sum_eos_after_density_correction_if_needed` — subroutine
- `sum_eos_before_density_correction_if_needed` — subroutine
- `trail_restored_halves_if_needed` — subroutine
- `update_project_nu_alignment_lowpass` — subroutine
- `update_project_resolution_metadata` — subroutine
- `write_benchmark` — subroutine
- `write_nonuniform_outputs` — subroutine
- `write_nu_highres_steps_for_state` — subroutine

---
## Module: simple_commanders_refine3D

Files:
- `main/commanders/simple/simple_commanders_refine3D.f90`

Uses:
- `simple_abinitio_utils`
- `simple_commanders_api`
- `simple_commanders_rec`
- `simple_core_module_api`
- `simple_estimate_ssnr`
- `simple_nu_filter`
- `simple_pftc_srch_api`
- `simple_refine3d_fnames`
- `simple_refine3d_strategy`
- `simple_strategy3d_matcher`

Public symbols:
- `cleanup_init_vols` — subroutine
- `commander_refine3D` — type
- `commander_refine3D_auto` — type
- `commander_refine3D_distr_worker` — type
- `commander_refine3D_multi` — type
- `configure_refine3D_multi_stages` — subroutine
- `ensure_all_active_particles_updated` — subroutine
- `exec_nspace` — subroutine
- `exec_refine3D` — subroutine
- `exec_refine3D_auto` — subroutine
- `exec_refine3D_distr_worker` — subroutine
- `exec_refine3D_multi` — subroutine
- `initialize_state_volumes` — subroutine
- `make_projdir_class_samples` — subroutine
- `nspace_commander` — type
- `prepare_external_init_vol` — subroutine
- `prepare_nu_bootstrap_refs_from_raw_halves` — subroutine
- `prepare_refine3D_multi_class_sampling` — subroutine
- `prepare_startup_reconstruct3D_cline` — subroutine
- `read_update_coverage` — subroutine
- `run_refine3D_multi_missing_update` — subroutine
- `run_refine3D_multi_stage` — subroutine
- `seed_refine3D_auto_nonuniform_lpset` — subroutine
- `set_refine3D_auto_sampling` — subroutine
- `set_refine3D_multi_downscaling` — subroutine
- `set_refine3D_multi_nstates` — subroutine
- `set_refine3D_multi_sampling` — subroutine
- `validate_input_volumes` — subroutine
- `validate_refine3D_multi_combine_eo` — subroutine
- `validate_refine3D_multi_filtering` — subroutine
- `validate_refine3D_multi_mode` — subroutine
- `validate_refine3D_multi_prob_neigh_mode` — subroutine

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
- `calc_lplim_final_stage` — function
- `commander_estimate_lpstages` — type
- `commander_fsc` — type
- `commander_fsc_area_score` — type
- `commander_icm2D` — type
- `commander_icm3D` — type
- `commander_nu_filt3D` — type
- `commander_uniform_filter2D` — type
- `commander_uniform_filter3D` — type
- `exec_estimate_lpstages` — subroutine
- `exec_fsc` — subroutine
- `exec_fsc_area_score` — subroutine
- `exec_icm2D` — subroutine
- `exec_icm3D` — subroutine
- `exec_nu_filt3D` — subroutine
- `exec_uniform_filter2D` — subroutine
- `exec_uniform_filter3D` — subroutine

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
- `simple_memoize_ft_maps`
- `simple_projector`
- `simple_simple_volinterp`
- `simple_simulator`

Public symbols:
- `commander_simulate_movie` — type
- `commander_simulate_nanoparticle` — type
- `commander_simulate_noise` — type
- `commander_simulate_particles` — type
- `commander_simulate_subtomogram` — type
- `exec_simulate_movie` — subroutine
- `exec_simulate_nanoparticle` — subroutine
- `exec_simulate_noise` — subroutine
- `exec_simulate_particles` — subroutine
- `exec_simulate_subtomogram` — subroutine
- `gen_ptcl_pos` — function

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
- `commander_export_manifoldem_starproject` — type
- `commander_export_starproject` — type
- `commander_import_starproject` — type
- `exec_assign_optics_groups` — subroutine
- `exec_export_manifoldem_starproject` — subroutine
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
- `finish_cleanup` — subroutine
- `handle_random_selection` — subroutine
- `handle_state_class_selection` — subroutine

---
## Module: simple_commanders_test_class

Files:
- `main/commanders/test/simple_commanders_test_class.f90`

Uses:
- `simple_aff_prop`
- `simple_atoms`
- `simple_binary_tree_tester`
- `simple_bspline_smoother`
- `simple_chash_tester`
- `simple_cmdline_tester`
- `simple_commanders_api`
- `simple_fileio_tester`
- `simple_forked_process_tester`
- `simple_ftexp_shsrch`
- `simple_ftiter`
- `simple_gui_assembler_tester`
- `simple_gui_metadata_tester`
- `simple_hash_tester`
- `simple_hclust`
- `simple_http_post_tester`
- `simple_image`
- `simple_imghead`
- `simple_ipc_tcp_socket_tester`
- `simple_linked_list_tester`
- `simple_motion_gain_tester`
- `simple_multi_dendro_tester`
- `simple_online_var`
- `simple_ori_tester`
- `simple_oris`
- `simple_oris_tester`
- `simple_persistent_worker_message_tester`
- `simple_persistent_worker_server_tester`
- `simple_project_merge_tester`
- `simple_rec_list_tester`
- `simple_starfile_tester`
- `simple_strategy2d_utils`
- `simple_stream_api`
- `simple_string_tester`
- `simple_syslib_tester`
- `simple_test_utils`
- `simple_ui`
- `simple_ui_hash`
- `simple_vrefhash_tester`

Public symbols:
- `commander_test_strategy2D` — type
- `commander_test_ui_hash_test` — type
- `commander_test_units` — type
- `exec_test_strategy2D` — subroutine
- `exec_test_ui_hash_test` — subroutine
- `exec_test_units` — subroutine
- `simple_test_fit_line` — subroutine
- `test_euler_shift` — subroutine
- `test_multinomal` — subroutine

---
## Module: simple_commanders_test_fft

Files:
- `main/commanders/test/simple_commanders_test_fft.f90`

Uses:
- `gnufor2`
- `mod_phasecorr`
- `simple_builder`
- `simple_commanders_api`
- `simple_core_module_api`
- `simple_ftexp_shsrch`
- `simple_matcher_smpl_and_lplims`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_projector_pft`
- `simple_timer`

Public symbols:
- `commander_test_corrs2weights_test` — type
- `commander_test_eval_polarftcc` — type
- `commander_test_ft_expanded` — type
- `commander_test_gencorrs_fft` — type
- `commander_test_order_corr` — type
- `commander_test_phasecorr` — type
- `commander_test_rank_weights` — type
- `commander_test_rotate_ref` — type
- `exec_test_corrs2weights_test` — subroutine
- `exec_test_eval_polarftcc` — subroutine
- `exec_test_ft_expanded` — subroutine
- `exec_test_gencorrs_fft` — subroutine
- `exec_test_order_corr` — subroutine
- `exec_test_phasecorr` — subroutine
- `exec_test_rank_weights` — subroutine
- `exec_test_rotate_ref` — subroutine
- `fast_rotate_ref` — subroutine
- `rotate_ref` — subroutine
- `rotate_ref_2` — subroutine
- `rotate_ref_t` — subroutine

---
## Module: simple_commanders_test_geometry

Files:
- `main/commanders/test/simple_commanders_test_geometry.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_commanders_project_core`
- `simple_core_module_api`
- `simple_sp_project`
- `simple_sym`

Public symbols:
- `commander_test_angres` — type
- `commander_test_ori_test` — type
- `commander_test_oris_test` — type
- `commander_test_sym_test` — type
- `commander_test_uniform_euler` — type
- `commander_test_uniform_rot` — type
- `exec_test_angres` — subroutine
- `exec_test_ori_test` — subroutine
- `exec_test_oris_test` — subroutine
- `exec_test_sym_test` — subroutine
- `exec_test_uniform_euler` — subroutine
- `exec_test_uniform_rot` — subroutine
- `rgauss` — subroutine

---
## Module: simple_commanders_test_highlevel

Files:
- `main/commanders/test/simple_commanders_test_highlevel.f90`

Uses:
- `simple_atoms`
- `simple_builder`
- `simple_commanders_abinitio`
- `simple_commanders_abinitio2d`
- `simple_commanders_api`
- `simple_commanders_pick`
- `simple_commanders_preprocess`
- `simple_commanders_project_core`
- `simple_commanders_project_mov`
- `simple_commanders_project_ptcl`
- `simple_commanders_reproject`
- `simple_commanders_sim`
- `simple_commanders_stkops`
- `simple_commanders_validate`
- `simple_flex_analysis_strategy`
- `simple_flex_diffmap_rec3d`
- `simple_imghead`
- `simple_micproc`
- `simple_molecule_data`
- `simple_parameters`
- `simple_pcg_reconstruction`
- `simple_projfile_utils`
- `simple_qsys_env`
- `simple_stream_api`
- `simple_string_utils`
- `simple_ui`

Public symbols:
- `commander_test_flex_preimage_basis_ab` — type
- `commander_test_flex_preimage_identity` — type
- `commander_test_mini_stream` — type
- `commander_test_pcg_recon_ctf_free` — type
- `commander_test_pcg_recon_ctf_hetero` — type
- `commander_test_pcg_recon_deapod` — type
- `commander_test_pcg_recon_kernel` — type
- `commander_test_ptcls_ppca_subproject_distr` — type
- `commander_test_reproject` — type
- `commander_test_simulate_particles` — type
- `commander_test_simulated_workflow` — type
- `commander_test_subproject_distr` — type
- `exec_test_flex_preimage_basis_ab` — subroutine
- `exec_test_flex_preimage_identity` — subroutine
- `exec_test_mini_stream` — subroutine
- `exec_test_pcg_recon_ctf_free` — subroutine
- `exec_test_pcg_recon_ctf_hetero` — subroutine
- `exec_test_pcg_recon_deapod` — subroutine
- `exec_test_pcg_recon_kernel` — subroutine
- `exec_test_ptcls_ppca_subproject_distr` — subroutine
- `exec_test_reproject` — subroutine
- `exec_test_simulate_particles` — subroutine
- `exec_test_simulated_workflow` — subroutine
- `exec_test_subproject_distr` — subroutine
- `return_to_stage_root` — subroutine
- `update_project_path` — subroutine

---
## Module: simple_commanders_test_io

Files:
- `main/commanders/test/simple_commanders_test_io.f90`

Uses:
- `simple_atoms`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_image`
- `simple_imgfile`
- `simple_jpg`
- `simple_parameters`
- `simple_sp_project`
- `simple_stack_io`
- `simple_starfile`
- `simple_starfile_wrappers`
- `simple_starproject_tester`

Public symbols:
- `commander_test_imgfile` — type
- `commander_test_inside_write` — type
- `commander_test_io` — type
- `commander_test_io_parallel` — type
- `commander_test_mrc2jpeg` — type
- `commander_test_mrc_validation` — type
- `commander_test_stack_io` — type
- `commander_test_star_export` — type
- `commander_test_starfile_test` — type
- `exec_test_imgfile` — subroutine
- `exec_test_inside_write` — subroutine
- `exec_test_io` — subroutine
- `exec_test_io_parallel` — subroutine
- `exec_test_mrc2jpeg` — subroutine
- `exec_test_mrc_validation` — subroutine
- `exec_test_stack_io` — subroutine
- `exec_test_star_export` — subroutine
- `exec_test_starfile_test` — subroutine

---
## Module: simple_commanders_test_masks

Files:
- `main/commanders/test/simple_commanders_test_masks.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_api`
- `simple_commanders_atoms`
- `simple_image`
- `simple_image_msk`
- `simple_imgarr_utils`
- `simple_math`
- `simple_parameters`
- `simple_projector`

Public symbols:
- `commander_test_bounds_from_mask3D_test` — type
- `commander_test_graphene_mask` — type
- `commander_test_image_bin` — type
- `commander_test_mask` — type
- `commander_test_msk_routines` — type
- `commander_test_nano_mask` — type
- `commander_test_otsu_test` — type
- `commander_test_ptcl_center` — type
- `exec_test_bounds_from_mask3D_test` — subroutine
- `exec_test_graphene_mask` — subroutine
- `exec_test_image_bin` — subroutine
- `exec_test_mask` — subroutine
- `exec_test_msk_routines` — subroutine
- `exec_test_nano_mask` — subroutine
- `exec_test_otsu_test` — subroutine
- `exec_test_ptcl_center` — subroutine
- `test_center` — subroutine

---
## Module: simple_commanders_test_network

Files:
- `main/commanders/test/simple_commanders_test_network.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_commanders_api`
- `simple_distr_comm`
- `simple_socket_comm`
- `simple_string_utils`
- `unix`

Public symbols:
- `commander_test_socket_client` — type
- `commander_test_socket_comm_distr` — type
- `commander_test_socket_io` — type
- `commander_test_socket_server` — type
- `exec_test_socket_client` — subroutine
- `exec_test_socket_comm_distr` — subroutine
- `exec_test_socket_io` — subroutine
- `exec_test_socket_server` — subroutine
- `socket_server` — subroutine

---
## Module: simple_commanders_test_numerics

Files:
- `main/commanders/test/simple_commanders_test_numerics.f90`

Uses:
- `simple_commanders_api`
- `simple_kbinterpol`
- `simple_linalg`
- `simple_test_utils`
- `simple_timer`

Public symbols:
- `commander_test_eigh_test` — type
- `commander_test_kbinterpol_fast` — type
- `commander_test_maxnloc_test` — type
- `commander_test_neigh` — type
- `exec_test_eigh_test` — subroutine
- `exec_test_kbinterpol_fast` — subroutine
- `exec_test_maxnloc_test` — subroutine
- `exec_test_neigh` — subroutine
- `rnd_romat` — subroutine

---
## Module: simple_commanders_test_optimize

Files:
- `main/commanders/test/simple_commanders_test_optimize.f90`

Uses:
- `simple_builder`
- `simple_butterworth`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_commanders_atoms`
- `simple_commanders_reproject`
- `simple_image`
- `simple_math`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`
- `simple_parameters`

Public symbols:
- `commander_test_lbfgsb` — type
- `commander_test_lbfgsb_cosine` — type
- `commander_test_lplims` — type
- `commander_test_lpstages_test` — type
- `commander_test_opt_lp` — type
- `costfct` — function
- `costfct` — function
- `exec_test_lbfgsb` — subroutine
- `exec_test_lbfgsb_cosine` — subroutine
- `exec_test_lplims` — subroutine
- `exec_test_lpstages_test` — subroutine
- `exec_test_opt_lp` — subroutine
- `gradfct` — subroutine
- `gradfct` — subroutine

---
## Module: simple_commanders_test_parallel

Files:
- `main/commanders/test/simple_commanders_test_parallel.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_sim`
- `simple_image`
- `simple_timer`

Public symbols:
- `commander_test_coarrays` — type
- `commander_test_openacc` — type
- `commander_test_openmp` — type
- `commander_test_simd` — type
- `exec_test_coarrays` — subroutine
- `exec_test_openacc` — subroutine
- `exec_test_openmp` — subroutine
- `exec_test_simd` — subroutine

---
## Module: simple_commanders_test_single

Files:
- `main/commanders/test/simple_commanders_test_single.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_atoms`
- `simple_commanders_project_core`
- `simple_commanders_reproject`
- `simple_commanders_sim`
- `single_commanders_nano3d`
- `single_commanders_trajectory`

Public symbols:
- `commander_test_atoms_stats` — type
- `commander_test_detect_atoms` — type
- `commander_test_simulate_nanoparticle` — type
- `commander_test_single_workflow` — type
- `exec_test_atoms_stats` — subroutine
- `exec_test_detect_atoms` — subroutine
- `exec_test_simulate_nanoparticle` — subroutine
- `exec_test_single_workflow` — subroutine

---
## Module: simple_commanders_test_stats

Files:
- `main/commanders/test/simple_commanders_test_stats.f90`

Uses:
- `simple_aff_prop`
- `simple_binoris_io`
- `simple_class_sample_io`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_ctf`
- `simple_decay_funs`
- `simple_image`
- `simple_kpca_svd`
- `simple_linalg`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_pca_svd`
- `simple_ppca`
- `simple_refine3d_fnames`
- `simple_rnd`
- `simple_sp_project`

Public symbols:
- `commander_test_class_sample_test` — type
- `commander_test_clustering` — type
- `commander_test_ctf_test` — type
- `commander_test_eo_diff` — type
- `commander_test_extr_frac` — type
- `commander_test_multinomal_test` — type
- `commander_test_pca_all` — type
- `commander_test_pca_imgvar` — type
- `commander_test_sp_project` — type
- `exec_test_class_sample_test` — subroutine
- `exec_test_clustering` — subroutine
- `exec_test_ctf_test` — subroutine
- `exec_test_eo_diff` — subroutine
- `exec_test_extr_frac` — subroutine
- `exec_test_multinomal_test` — subroutine
- `exec_test_pca_all` — subroutine
- `exec_test_pca_imgvar` — subroutine
- `exec_test_sp_project` — subroutine
- `make_cs` — subroutine

---
## Module: simple_commanders_test_utils

Files:
- `main/commanders/test/simple_commanders_test_utils.f90`

Uses:
- `simple_ansi_ctrls`
- `simple_atoms`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_defs_fname`
- `simple_image`
- `simple_imghead`
- `simple_molecule_data`
- `simple_nice`
- `simple_ppca`
- `simple_segmentation`
- `simple_stackops`
- `simple_test_utils`
- `simple_testfuns`

Public symbols:
- `commander_test_ansi_colors` — type
- `commander_test_binoris_io_test` — type
- `commander_test_binoris_test` — type
- `commander_test_cif2mrc` — type
- `commander_test_cif2pdb` — type
- `commander_test_cmdline` — type
- `commander_test_install` — type
- `commander_test_nice` — type
- `commander_test_pdb2mrc` — type
- `commander_test_peak_thres_fdr` — type
- `commander_test_serialize` — type
- `commander_test_stringmatch` — type
- `exec_test_ansi_colors` — subroutine
- `exec_test_binoris_io_test` — subroutine
- `exec_test_binoris_test` — subroutine
- `exec_test_cif2mrc` — subroutine
- `exec_test_cif2pdb` — subroutine
- `exec_test_cmdline` — subroutine
- `exec_test_install` — subroutine
- `exec_test_nice` — subroutine
- `exec_test_pdb2mrc` — subroutine
- `exec_test_peak_thres_fdr` — subroutine
- `exec_test_serialize` — subroutine
- `exec_test_stringmatch` — subroutine

---
## Module: simple_commanders_validate

Files:
- `main/commanders/simple/simple_commanders_validate.f90`

Uses:
- `simple_commanders_abinitio2d`
- `simple_commanders_api`
- `simple_commanders_cavgs`
- `simple_commanders_pick`
- `simple_commanders_preprocess`
- `simple_commanders_project_core`
- `simple_commanders_project_mov`
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
- `simple_clustering_utils`
- `simple_commanders_api`
- `simple_commanders_reproject`
- `simple_dock_vols`
- `simple_image_msk`
- `simple_imgproc`
- `simple_kpca_svd`
- `simple_pca`
- `simple_pca_svd`
- `simple_ppca`
- `simple_procimgstk`
- `simple_projector`
- `simple_segmentation`
- `simple_simple_volinterp`
- `simple_strategy2d_utils`
- `simple_symanalyzer`
- `simple_volanalyzer`
- `simple_volcluster`
- `simple_volpft_symsrch`

Public symbols:
- `commander_centervol` — type
- `commander_dock_volpair` — type
- `commander_noisevol` — type
- `commander_postprocess` — type
- `commander_ppca_volvar` — type
- `commander_sharpvol` — type
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
- `exec_ppca_volvar` — subroutine
- `exec_sharpvol` — subroutine
- `exec_symaxis_search` — subroutine
- `exec_symmetrize_map` — subroutine
- `exec_symmetry_test` — subroutine
- `exec_volanalyze` — subroutine
- `exec_volcluster` — subroutine
- `exec_volops` — subroutine
- `postprocess_volume_from_files` — subroutine

---
## Module: simple_complex_ppca

Files:
- `main/pca/simple_complex_ppca.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `complex_ppca` — type

Private symbols:
- `denoise_complex_ppca` — subroutine
- `fit_covariance_complex_ppca` — subroutine
- `generate_complex_ppca` — subroutine
- `get_eigvals_complex_ppca` — function
- `get_feat_complex_ppca` — function
- `get_mean_complex_ppca` — function
- `hermitian_jacobi` — subroutine
- `kill_complex_ppca` — subroutine
- `master_complex_ppca` — subroutine
- `new_complex_ppca` — subroutine
- `outer_hermitian` — function
- `reconstruct_external_complex_ppca` — subroutine
- `sort_eigs_desc` — subroutine
- `stream_finalize_complex_ppca` — subroutine
- `stream_reset_complex_ppca` — subroutine
- `stream_update_complex_ppca` — subroutine

---
## Module: simple_convergence

Files:
- `main/simple_convergence.f90`

Uses:
- `cplot2d_wrapper_module`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_parameters`
- `simple_progress`

Private symbols:
- `append_stats` — subroutine
- `check_conv2D` — function
- `check_conv3D` — function
- `plot_projdirs` — subroutine
- `print_iteration` — subroutine
- `read` — subroutine

---
## Module: simple_core_module_api

Files:
- `main/apis/simple_core_api.f90`

Uses:
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
- `utils/simple_corrmat.f90`

Uses:
- `simple_builder`
- `simple_defs`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`

---
## Module: simple_ctf

Files:
- `main/ctf/simple_ctf.f90`

Uses:
- `simple_core_module_api`

---
## Module: simple_ctf_estimate_cost

Files:
- `main/ctf/simple_ctf_estimate_cost.f90`

Uses:
- `simple_core_module_api`
- `simple_ctf`
- `simple_image`
- `simple_opt_de`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`

Private symbols:
- `df_costfun` — subroutine
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
- `simple_starfile_wrappers`
- `simple_timer`

Private symbols:
- `calc_ctfres` — subroutine
- `calc_frc` — subroutine
- `calc_icefrac` — subroutine
- `calc_pspec_stats` — subroutine
- `calc_tilt` — subroutine
- `ctf2pspecimg` — subroutine
- `fit` — subroutine
- `fit_patches` — subroutine
- `fit_polynomial` — subroutine
- `ft2img` — subroutine
- `gen_profiles` — subroutine
- `gen_resmsk` — subroutine
- `gen_roavspec1d` — subroutine
- `gen_tiles` — subroutine
- `get_parms` — subroutine
- `get_pspec` — subroutine
- `kill` — subroutine
- `mask_central_disc` — subroutine
- `mic2spec` — subroutine
- `new` — subroutine
- `norm_pspec` — subroutine
- `pix2poly` — subroutine
- `pix2polyvals` — subroutine
- `plot_ctf` — subroutine
- `poly` — function
- `rank_spec` — subroutine
- `read_doc` — subroutine
- `refine` — subroutine
- `srch` — subroutine
- `subtr_backgr` — subroutine
- `write_diagnostic` — subroutine
- `write_doc` — subroutine
- `write_star` — subroutine

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
- `inv_cos_decay` — function
- `inv_nsampl_decay` — function
- `nsampl_decay` — function

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
- `set_cluster2D_defaults` — subroutine

---
## Module: simple_defs

Files:
- `defs/simple_defs.f90`

---
## Module: simple_defs_atoms

Files:
- `defs/simple_defs_atoms.f90`

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

---
## Module: simple_defs_fname

Files:
- `defs/simple_defs_fname.f90`

---
## Module: simple_defs_ori

Files:
- `defs/simple_defs_ori.f90`

Public symbols:
- `get_oriparam_flag` — function

---
## Module: simple_defs_stream

Files:
- `defs/simple_defs_stream.f90`

---
## Module: simple_defs_string

Files:
- `defs/simple_defs_string.f90`

---
## Module: simple_denoise_movies

Files:
- `main/motion/simple_denoise_movies.f90`

Uses:
- `simple_core_module_api`
- `simple_diff_map_denoise`
- `simple_diff_map_graphs`
- `simple_diffusion_maps`
- `simple_image`

Public symbols:
- `diffmap_denoise_image_stack` — subroutine

Private symbols:
- `validate_realspace_2d_stack` — subroutine

---
## Module: simple_denoise_project_strategy

Files:
- `main/strategies/parallelization/simple_denoise_project_strategy.f90`

Uses:
- `simple_builder`
- `simple_classaverager`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_default_clines`
- `simple_diff_map_graphs`
- `simple_diffusion_maps`
- `simple_eulspace_neigh_map`
- `simple_exec_helpers`
- `simple_image`
- `simple_image_msk`
- `simple_imgarr_utils`
- `simple_imgfile`
- `simple_imgproc`
- `simple_linalg`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_sp_project`
- `simple_srch_sort_loc`

Public symbols:
- `create_denoise_project_strategy` — function
- `denoise_project_master_strategy` — type
- `denoise_project_shmem_strategy` — type
- `denoise_project_strategy` — type
- `denoise_project_worker_strategy` — type

Private symbols:
- `apply_defaults` — subroutine
- `assign_global_so3_mixtures` — subroutine
- `build_diffmap_so3_mixture_graphs` — subroutine
- `calc_diffmap_reconstruction_error` — subroutine
- `calc_diffmap_residual_energy_ratio` — subroutine
- `cleanup_interface` — subroutine
- `count_data_lines` — subroutine
- `estimate_diffmap_denoise_icm_rank_stats` — subroutine
- `estimate_diffmap_denoise_rank` — subroutine
- `exec_interface` — subroutine
- `filter_classes_by_assignment` — subroutine
- `finalize_diffmap_worker_outputs` — subroutine
- `finalize_interface` — subroutine
- `graph_nystrom_residual_preimage` — subroutine
- `init_common` — subroutine
- `init_interface` — subroutine
- `master_cleanup` — subroutine
- `master_execute` — subroutine
- `master_finalize_run` — subroutine
- `master_initialize` — subroutine
- `pack_diffmap_so3_mixture_pinds` — subroutine
- `prepare_class_partitions` — subroutine
- `preserve_diffmap_frc2D` — subroutine
- `read_diffmap_worker_maps` — subroutine
- `read_int_file` — subroutine
- `run_denoise_project` — subroutine
- `run_diffmap_denoise_icm_rank_seed` — subroutine
- `run_diffmap_so3_mixture_graphs` — subroutine
- `select_diffmap_denoise_rank_icm` — subroutine
- `setup_diffmap_so3_mixture_map` — subroutine
- `shmem_cleanup` — subroutine
- `shmem_execute` — subroutine
- `shmem_finalize_run` — subroutine
- `shmem_initialize` — subroutine
- `sort_order_by_weight_desc` — subroutine
- `update_diffmap_denoise_icm_rank_site` — subroutine
- `validate_denoise_project` — subroutine
- `worker_cleanup` — subroutine
- `worker_execute` — subroutine
- `worker_finalize_run` — subroutine
- `worker_initialize` — subroutine
- `write_diffmap_denoised_class` — subroutine
- `write_diffmap_partial_project` — subroutine
- `write_diffmap_project` — subroutine
- `write_diffmap_stack_image` — subroutine

---
## Module: simple_diff_map_denoise

Files:
- `main/pca/simple_diff_map_denoise.f90`

Uses:
- `simple_core_module_api`
- `simple_diff_map_graphs`
- `simple_diffusion_maps`
- `simple_image`
- `simple_imgarr_utils`
- `simple_linalg`
- `simple_parameters`
- `simple_polarft_calc`

Public symbols:
- `calc_diffmap_reconstruction_error` — subroutine
- `calc_diffmap_residual_energy_ratio` — subroutine
- `estimate_diffmap_denoise_rank` — subroutine
- `graph_coeffproj_denoise` — subroutine
- `graph_nystrom_residual_preimage` — subroutine
- `select_spectral_rank_icm` — subroutine

Private symbols:
- `accumulate_pft_angular_mode` — subroutine
- `estimate_diffmap_denoise_icm_rank_stats` — subroutine
- `extract_pft_angular_mode` — subroutine
- `prepare_steerable_pft_params` — subroutine
- `project_pft_angular_mode` — subroutine
- `run_diffmap_denoise_icm_rank_seed` — subroutine
- `scatter_fourier_sample` — subroutine
- `so2_coeffproj_denoise` — subroutine
- `synthesize_pft_residual` — subroutine
- `update_diffmap_denoise_icm_rank_site` — subroutine

---
## Module: simple_diff_map_graphs

Files:
- `main/pca/simple_diff_map_graphs.f90`

Uses:
- `simple_core_module_api`
- `simple_ori`
- `simple_parameters`
- `simple_sp_project`
- `simple_srch_sort_loc`

Public symbols:
- `build_cls_split_graph` — subroutine
- `build_euclidean_knn_graph` — subroutine
- `build_gated_euclidean_graph_from_neighbors` — subroutine
- `build_gated_euclidean_knn_graph` — subroutine
- `build_orientation_knn_graph` — subroutine
- `diffmap_graph` — type
- `find_gated_euclidean_neighbors_rows` — subroutine
- `graph_matvec` — subroutine
- `projection_occupancy_weights` — subroutine

Private symbols:
- `diffmap_graph_degree` — subroutine
- `estimate_ferguson_bandwidth` — subroutine
- `find_euclidean_neighbors` — subroutine
- `find_orientation_neighbors` — subroutine
- `insert_neighbor` — subroutine
- `kill_diffmap_graph` — subroutine
- `make_singleton_graph` — subroutine
- `normalize_diffmap_graph` — subroutine
- `pack_scalar_knn_to_csr` — subroutine

---
## Module: simple_diffusion_maps

Files:
- `main/pca/simple_diffusion_maps.f90`

Uses:
- `simple_core_module_api`
- `simple_diff_map_graphs`
- `simple_linalg`

Public symbols:
- `diffusion_map_embedder` — type
- `embed_graph` — subroutine
- `normalize_coords` — subroutine

Private symbols:
- `embed` — subroutine
- `set_params` — subroutine

---
## Module: simple_discrete_stack_io

Files:
- `fileio/simple_discrete_stack_io.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_imgfile`

Private symbols:
- `cache_stack_info` — subroutine
- `clear_stack_cache` — subroutine
- `get_stack_info` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `open_1` — subroutine
- `read` — subroutine

---
## Module: simple_distr_comm

Files:
- `utils/comm/simple_distr_comm.f90`

Uses:
- `simple_socket_comm`
- `unix`

Public symbols:
- `distr_comm` — type

Private symbols:
- `init` — subroutine
- `server` — subroutine
- `thread_comm` — type

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
## Module: simple_edges_sqwins

Files:
- `main/interp/simple_edges_sqwins.f90`

Uses:
- `simple_defs`
- `simple_error`

---
## Module: simple_eer_factory

Files:
- `fileio/simple_eer_factory.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_tifflib`

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
- `kill` — subroutine
- `new` — subroutine
- `prep_gainref` — subroutine
- `read` — subroutine
- `set_dims` — subroutine

---
## Module: simple_error

Files:
- `fileio/simple_error.f90`

Uses:
- `simple_defs`

Public symbols:
- `simple_error_check` — subroutine
- `simple_exception` — subroutine

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
- `calc_dose_weights` — subroutine
- `fsc2optlp` — function
- `fsc2optlp_sub` — subroutine
- `fsc2ssnr` — function
- `gaussian_filter` — subroutine
- `get_resolution` — subroutine
- `get_resolution_at_fsc` — subroutine
- `lowpass_from_klim` — subroutine
- `lpstages` — subroutine
- `lpstages_fast` — subroutine
- `lpstages_setlims` — subroutine
- `mskdiam2lplimits` — subroutine
- `mskdiam2streamresthreshold` — subroutine
- `ssnr2fsc` — function
- `ssnr2optlp` — function

Private symbols:
- `calc_lpinfo` — subroutine
- `calc_scaleinfo` — subroutine
- `calc_scaleinfo` — subroutine
- `print_lpinfo` — subroutine
- `print_scaleinfo` — subroutine

---
## Module: simple_euclid_sigma2

Files:
- `main/simple_euclid_sigma2.f90`

Uses:
- `simple_core_module_api`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_sigma2_binfile`
- `simple_starfile_wrappers`

Public symbols:
- `average_sigma2_groups` — subroutine
- `consolidate_sigma2_groups` — subroutine
- `fill_sigma2_before_nyq` — subroutine
- `sigma2_star_from_iter` — function
- `split_sigma2_into_groups` — subroutine
- `test_unit` — subroutine
- `write_groups_starfile` — subroutine

Private symbols:
- `allocate_ptcls` — subroutine
- `calc_sigma2` — subroutine
- `consolidate_sigma2_history` — subroutine
- `get_kfromto` — function
- `init_from_group_header` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `parse_key_int_pair` — subroutine
- `parse_key_string_pair` — subroutine
- `read_groups` — subroutine
- `read_groups_starfile` — subroutine
- `read_part` — subroutine
- `read_sigma2_groups` — subroutine
- `write_info` — subroutine
- `write_sigma2` — subroutine

---
## Module: simple_eul_prob_tab

Files:
- `main/simple_eul_prob_tab.f90`

Uses:
- `simple_builder`
- `simple_eul_prob_tab_utils`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_type_defs`

Public symbols:
- `eul_prob_tab` — type

Private symbols:
- `advance_ref_head` — subroutine
- `advance_state_head` — subroutine
- `assign_candidate` — subroutine
- `assign_current_ref` — subroutine
- `assign_greedy_state_labels` — subroutine
- `assign_refs_for_state` — subroutine
- `begin_write` — subroutine
- `estimate_shift_seed` — subroutine
- `fill_tab` — subroutine
- `fill_tab_range` — subroutine
- `fill_tab_state_only` — subroutine
- `fill_tab_state_only_range` — subroutine
- `flush_candidate_buffers` — subroutine
- `gather_active_refs` — subroutine
- `get_particle_context` — subroutine
- `initialize_storage` — subroutine
- `kill` — subroutine
- `make_ref_candidate` — function
- `new_assignment` — subroutine
- `new_common` — subroutine
- `new_compact_global` — subroutine
- `new_global` — subroutine
- `new_state` — subroutine
- `new_worker` — subroutine
- `prepare_ref_score_vectors` — subroutine
- `process_particle_with_shift` — subroutine
- `process_particle_without_shift` — subroutine
- `read_assignment` — subroutine
- `read_state_tab` — subroutine
- `read_tab_to_glob` — subroutine
- `record_ref_eval` — subroutine
- `ref_assign` — subroutine
- `refine_best_refs` — subroutine
- `replace_ref_eval` — subroutine
- `reset_ref_frontier` — subroutine
- `score_ref` — subroutine
- `state_assign` — subroutine
- `write_assignment` — subroutine
- `write_state_tab` — subroutine
- `write_tab` — subroutine

---
## Module: simple_eul_prob_tab2D

Files:
- `main/simple_eul_prob_tab2D.f90`

Uses:
- `simple_builder`
- `simple_decay_funs`
- `simple_eul_prob_tab_utils`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`

Public symbols:
- `eul_prob_tab2D` — type

Private symbols:
- `advance_class_head` — subroutine
- `assign_frontier_ws` — type
- `assign_graph_ws` — type
- `assign_particles_from_frontier` — subroutine
- `assign_remaining_particles` — subroutine
- `begin_write` — subroutine
- `build_sparse_assignment_graph` — subroutine
- `commit_selected_assignment` — subroutine
- `dealloc_eval2D_sparse_ws` — subroutine
- `eval2D_sparse_ws` — type
- `fill_tab` — subroutine
- `fill_tab_prob_snhc` — subroutine
- `fill_tab_prob_snhc_range` — subroutine
- `fill_tab_range` — subroutine
- `finalize_assignment` — subroutine
- `finalize_candidate_assignment` — subroutine
- `flush_candidate_buffers` — subroutine
- `init_eval2D_sparse_ws` — subroutine
- `init_frontier` — subroutine
- `initialize_storage` — subroutine
- `kill` — subroutine
- `new_assignment` — subroutine
- `new_common` — subroutine
- `new_global` — subroutine
- `new_worker` — subroutine
- `process_particle` — subroutine
- `read_assignment` — subroutine
- `read_sparse_tab_to_glob` — subroutine
- `read_tab_to_glob` — subroutine
- `read_tabs_to_glob` — subroutine
- `record_sparse_eval` — subroutine
- `ref_assign_likelihood` — subroutine
- `ref_assign_sparse_likelihood` — subroutine
- `refine_best_classes` — subroutine
- `sample_frontier_likelihood` — subroutine
- `score_class` — subroutine
- `sync_frontier_class` — subroutine
- `update_frontier_after_assignment` — subroutine
- `write_assignment` — subroutine
- `write_tab` — subroutine

---
## Module: simple_eul_prob_tab_neigh

Files:
- `main/simple_eul_prob_tab_neigh.f90`

Uses:
- `simple_builder`
- `simple_decay_funs`
- `simple_eul_prob_tab`
- `simple_eul_prob_tab_utils`
- `simple_ori`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`

Public symbols:
- `eul_prob_tab_neigh` — type

Private symbols:
- `add_pooled_subspace` — subroutine
- `advance_ref_head` — subroutine
- `alloc_coarse_ws` — subroutine
- `alloc_eval_stoch_ws` — subroutine
- `alloc_eval_subspace_ws` — subroutine
- `assign_frontier_ws` — type
- `assign_graph_ws` — type
- `assign_greedy_state_labels` — subroutine
- `assign_particles_for_state` — subroutine
- `assign_particles_globally` — subroutine
- `assign_particles_single_state` — subroutine
- `assign_remaining_particles_from_best_touched_ref` — subroutine
- `build_geometric_neighborhood` — subroutine
- `build_pooled_neighborhood` — subroutine
- `build_sparse_assignment_graph` — subroutine
- `calc_previous_corr` — subroutine
- `coarse_search_ws` — type
- `commit_selected_assignment` — subroutine
- `consider_peak_subspace` — subroutine
- `dealloc_coarse_ws` — subroutine
- `dealloc_eval_ws` — subroutine
- `estimate_shift_seed` — subroutine
- `eval_ws` — type
- `evaluate_direct_stochastic_refs` — subroutine
- `evaluate_neighborhood` — subroutine
- `fill_tab_neigh` — subroutine
- `fill_tab_neigh_range` — subroutine
- `fill_tab_stoch_range` — subroutine
- `fill_tab_subspace_range` — subroutine
- `find_peak_subspaces` — subroutine
- `get_particle_context` — subroutine
- `init_eval_stoch_ws` — subroutine
- `init_eval_subspace_ws` — subroutine
- `init_frontier` — subroutine
- `kill_neigh` — subroutine
- `new_neigh` — subroutine
- `new_neigh_global` — subroutine
- `process_particle_stoch` — subroutine
- `process_particle_subspace` — subroutine
- `read_sparse_tab_to_glob` — subroutine
- `read_tabs_to_glob` — subroutine
- `record_sparse_eval` — subroutine
- `ref_assign_neigh` — subroutine
- `refine_best_neighbors` — subroutine
- `score_direct_ref` — subroutine
- `score_subspace_ref` — subroutine
- `seed_fallback_if_empty` — subroutine
- `sync_frontier_ref` — subroutine
- `update_frontier_after_assignment` — subroutine
- `validate_neigh_mode` — subroutine

---
## Module: simple_eul_prob_tab_utils

Files:
- `main/simple_eul_prob_tab_utils.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_math`
- `simple_oris`
- `simple_rnd`
- `simple_type_defs`

Public symbols:
- `prob_candidate` — type
- `prob_candidate_buffer` — type
- `prob_candidate_store` — type

---
## Module: simple_eulspace_neigh_map

Files:
- `utils/simple_eulspace_neigh_map.f90`

Uses:
- `simple_core_module_api`
- `simple_srchspace_map`

Public symbols:
- `eulspace_neigh_map` — type

Private symbols:
- `get_full2sub_map` — function
- `get_neighbors_list` — subroutine
- `get_neighbors_mask` — subroutine
- `get_neighbors_mask_pooled` — subroutine
- `kill` — subroutine
- `new_from_labels` — subroutine
- `new_from_spaces` — subroutine

---
## Module: simple_exec_abinitio3D

Files:
- `main/exec/simple_exec_abinitio3D.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_abinitio`
- `simple_commanders_resolest`
- `simple_commanders_volops`
- `simple_exec_helpers`
- `simple_string`

Public symbols:
- `exec_abinitio3D_commander` — subroutine

---
## Module: simple_exec_api

Files:
- `main/apis/simple_exec_api.f90`

Uses:
- `iso_fortran_env`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_exec_abinitio3d`
- `simple_exec_cavgproc`
- `simple_exec_cluster2d`
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
- `simple_exec_refine3d`
- `simple_exec_res`
- `simple_exec_sieve`
- `simple_exec_sim`
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
## Module: simple_exec_cluster2D

Files:
- `main/exec/simple_exec_cluster2D.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_abinitio2d`
- `simple_commanders_cavgs`
- `simple_commanders_cluster2d`
- `simple_commanders_mkcavgs`
- `simple_commanders_project_cls`
- `simple_exec_helpers`
- `simple_stream_abinitio2d_chunks`
- `simple_string`

Public symbols:
- `exec_cluster2D_commander` — subroutine

---
## Module: simple_exec_denoise

Files:
- `main/exec/simple_exec_denoise.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_cluster2d`
- `simple_commanders_denoise`
- `simple_commanders_flex_analysis`
- `simple_commanders_imgops`
- `simple_commanders_resolest`
- `simple_commanders_volops`

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
- `simple_commanders_imgops`
- `simple_commanders_resolest`

Public symbols:
- `exec_filter_commander` — subroutine

---
## Module: simple_exec_helpers

Files:
- `utils/simple_exec_helpers.f90`

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
- `simple_commanders_misc`
- `simple_commanders_project_ptcl`

Public symbols:
- `exec_other_commander` — subroutine

---
## Module: simple_exec_preproc

Files:
- `main/exec/simple_exec_preproc.f90`

Uses:
- `simple_cmdline`
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
- `simple_commanders_project_cls`
- `simple_commanders_project_core`
- `simple_commanders_project_mov`
- `simple_commanders_project_ptcl`
- `simple_commanders_relion`
- `simple_commanders_starproject`
- `single_commanders_trajectory`

Public symbols:
- `exec_project_commander` — subroutine

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
- `simple_commanders_pcg_recon`
- `simple_commanders_reproject`
- `simple_commanders_volops`

Public symbols:
- `exec_volume_commander` — subroutine

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
- `killimgbatch` — subroutine
- `prepimgbatch` — subroutine
- `set_extract_defaults` — subroutine

---
## Module: simple_fftw3

Files:
- `defs/simple_fftw3.f90`

Public symbols:
- `fftw_cleanup` — subroutine
- `fftw_cleanup_threads` — subroutine
- `fftw_destroy_plan` — subroutine
- `fftw_execute_dft` — subroutine
- `fftw_execute_dft_c2r` — subroutine
- `fftw_execute_dft_r2c` — subroutine
- `fftw_execute_r2r` — subroutine
- `fftw_execute_split_dft` — subroutine
- `fftw_execute_split_dft_c2r` — subroutine
- `fftw_execute_split_dft_r2c` — subroutine
- `fftw_export_wisdom` — subroutine
- `fftw_export_wisdom_to_file` — subroutine
- `fftw_flops` — subroutine
- `fftw_forget_wisdom` — subroutine
- `fftw_fprint_plan` — subroutine
- `fftw_free` — subroutine
- `fftw_free_str` — subroutine
- `fftw_iodim` — type
- `fftw_iodim64` — type
- `fftw_plan_with_nthreads` — subroutine
- `fftw_print_plan` — subroutine
- `fftw_set_timelimit` — subroutine
- `fftwf_cleanup` — subroutine
- `fftwf_cleanup_threads` — subroutine
- `fftwf_destroy_plan` — subroutine
- `fftwf_execute_dft` — subroutine
- `fftwf_execute_dft_c2r` — subroutine
- `fftwf_execute_dft_r2c` — subroutine
- `fftwf_execute_r2r` — subroutine
- `fftwf_execute_split_dft` — subroutine
- `fftwf_execute_split_dft_c2r` — subroutine
- `fftwf_execute_split_dft_r2c` — subroutine
- `fftwf_export_wisdom` — subroutine
- `fftwf_export_wisdom_to_file` — subroutine
- `fftwf_flops` — subroutine
- `fftwf_forget_wisdom` — subroutine
- `fftwf_fprint_plan` — subroutine
- `fftwf_free` — subroutine
- `fftwf_iodim` — type
- `fftwf_iodim64` — type
- `fftwf_plan_with_nthreads` — subroutine
- `fftwf_print_plan` — subroutine
- `fftwf_set_timelimit` — subroutine

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

---
## Module: simple_finch

Files:
- `utils/clustering/simple_finch.f90`

Uses:
- `simple_core_module_api`
- `simple_kd_tree`

Public symbols:
- `finch_hierarchy` — type
- `finch_partition_from_first_neighbors` — subroutine
- `finch_representatives` — subroutine
- `fit_finch` — subroutine
- `refine_finch_level` — subroutine
- `select_finch_level` — subroutine

Private symbols:
- `append_level` — subroutine
- `finch_get_coarsest_labels` — subroutine
- `finch_get_finest_labels` — subroutine
- `finch_get_labels` — subroutine
- `grow_ward_heap` — subroutine
- `initialize_refinement_nodes` — subroutine
- `kill_finch_hierarchy` — subroutine
- `kill_ward_heap` — subroutine
- `reserve_ward_heap` — subroutine
- `swap_ward_edges` — subroutine
- `trim_hierarchy` — subroutine
- `union_components` — subroutine
- `ward_edge_heap` — type
- `ward_heap_pop` — subroutine
- `ward_heap_push` — subroutine
- `weighted_cluster_centroids` — subroutine

---
## Module: simple_flex_analysis_strategy

Files:
- `main/flex/simple_flex_analysis_strategy.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_diff_map_denoise`
- `simple_diff_map_graphs`
- `simple_diffusion_maps`
- `simple_flex_diffmap_features`
- `simple_flex_diffmap_preimage`
- `simple_flex_diffmap_rec3d`
- `simple_image`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_sigma2_files`
- `simple_sp_project`

Public symbols:
- `create_flex_analysis_strategy` — function
- `fit_flex_analysis_embedding` — subroutine
- `flex_analysis_master_strategy` — type
- `flex_analysis_shmem_strategy` — type
- `flex_analysis_strategy` — type
- `flex_analysis_worker_strategy` — type
- `flex_embedding_result` — type
- `run_flex_preimage_basis_ab_test` — subroutine

Private symbols:
- `apply_defaults` — subroutine
- `build_flex_diffmap_graph` — subroutine
- `capture_result` — subroutine
- `cleanup_distributed_analysis_parts` — subroutine
- `cleanup_interface` — subroutine
- `collect_target_coefficients` — subroutine
- `embed_flex_graph` — subroutine
- `ensure_flex_sigma_support` — subroutine
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `finish_analysis_outputs` — subroutine
- `init_common` — subroutine
- `init_interface` — subroutine
- `init_model_context` — subroutine
- `kill_flex_embedding_result` — subroutine
- `log_flex_coordinate_scales` — subroutine
- `master_cleanup` — subroutine
- `master_execute` — subroutine
- `master_finalize_run` — subroutine
- `master_initialize` — subroutine
- `preimage_ab_state_fname` — function
- `prepare_particle_partitions` — subroutine
- `prepare_reconstruction_partitions` — subroutine
- `read_flex_reconstruction_assignment` — subroutine
- `read_graph_parts` — subroutine
- `read_int_file` — subroutine
- `run_flex_analysis` — subroutine
- `select_flex_spectral_rank` — subroutine
- `select_particles` — subroutine
- `shmem_cleanup` — subroutine
- `shmem_execute` — subroutine
- `shmem_finalize_run` — subroutine
- `shmem_initialize` — subroutine
- `validate_inputs` — subroutine
- `validate_selected_particles` — subroutine
- `worker_cleanup` — subroutine
- `worker_execute` — subroutine
- `worker_finalize_run` — subroutine
- `worker_initialize` — subroutine
- `write_coordinates` — subroutine
- `write_discrete_state_project` — subroutine
- `write_graph_part` — subroutine
- `write_graph_summary` — subroutine
- `write_preimage_basis_ab_metrics` — subroutine
- `write_preimage_particle_map` — subroutine
- `write_spectrum` — subroutine

---
## Module: simple_flex_diffmap_features

Files:
- `main/flex/simple_flex_diffmap_features.f90`

Uses:
- `simple_builder`
- `simple_classaverager`
- `simple_core_module_api`
- `simple_ctf`
- `simple_image`
- `simple_imgarr_utils`
- `simple_math_ft`
- `simple_memoize_ft_maps`
- `simple_ori`
- `simple_parameters`
- `simple_simple_volinterp`
- `simple_sp_project`
- `simple_stack_io`

Public symbols:
- `assemble_flex_diffmap_feature_parts` — subroutine
- `flex_projection_directions` — subroutine
- `map_flex_registered_to_native_project` — subroutine
- `prepare_flex_diffmap_feature_part` — subroutine
- `prepare_flex_diffmap_features` — subroutine
- `read_flex_diffmap_feature_parts` — subroutine
- `write_flex_mean_projection_stack` — subroutine

Private symbols:
- `apply_sigma_shell_whitening` — subroutine
- `read_mean_projection_stack` — subroutine
- `write_registered_project` — subroutine

---
## Module: simple_flex_diffmap_preimage

Files:
- `main/flex/simple_flex_diffmap_preimage.f90`

Uses:
- `simple_clustering_utils`
- `simple_core_module_api`
- `simple_srch_sort_loc`

Public symbols:
- `build_flex_preimage_kernel_weights` — subroutine
- `select_flex_diffmap_preimages` — subroutine
- `test_flex_preimage_bandwidth_decoupling` — subroutine

Private symbols:
- `dense_bandwidth` — function
- `median_default_real` — function
- `median_positive_d2` — function
- `raw_neff_state` — function

---
## Module: simple_flex_diffmap_rec3D

Files:
- `main/flex/simple_flex_diffmap_rec3D.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_flex_projected_latent_model`
- `simple_flex_reconstructor_latent_ops`
- `simple_gridding`
- `simple_image`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_parameters`
- `simple_reconstructor`
- `simple_sp_project`

Public symbols:
- `canonicalize_flex_preimage_coordinates` — subroutine
- `cleanup_flex_diffmap_rec_parts` — subroutine
- `reconstruct_flex_diffmap_local_linear_states` — subroutine
- `reconstruct_flex_diffmap_states` — subroutine
- `reconstruct_flex_diffmap_weighted_states` — subroutine
- `reduce_flex_diffmap_rec_parts` — subroutine
- `test_fake_preimage_against_reconstruct3D` — subroutine
- `test_flex_local_linear_preimage` — subroutine
- `write_flex_diffmap_rec_parts` — subroutine

Private symbols:
- `cleanup_reconstructors` — subroutine
- `flex_part_fname` — function
- `identity_matrix` — function
- `init_basis_reconstructor` — subroutine
- `init_mean_reconstructor` — subroutine
- `ll_solve` — subroutine
- `ll_solve_intercept` — function
- `ll_solve_slope` — function
- `prepare_project_fsc_lowpass_filters` — subroutine
- `prepare_unfiltered_model_params` — subroutine
- `reconstruct3D_reference` — subroutine
- `reconstruct3D_standard_residual_control` — subroutine
- `solve_upper_triangular` — subroutine
- `validate_model_tables` — subroutine
- `write_preimage_states` — subroutine
- `write_state` — subroutine

---
## Module: simple_flex_projected_latent_model

Files:
- `main/flex/simple_flex_projected_latent_model.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_flex_reconstructor_latent_ops`
- `simple_gridding`
- `simple_image`
- `simple_imgarr_utils`
- `simple_linalg`
- `simple_map_reduce`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_io`
- `simple_math`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_reconstructor`
- `simple_rnd`

Public symbols:
- `basis_fourier_energy` — function
- `canonicalize_projected_latent_basis` — subroutine
- `cleanup_planes` — subroutine
- `infer_latents_from_basis` — subroutine
- `initialize_latents` — subroutine
- `latent_covariance` — subroutine
- `latent_sdev` — function
- `orthonormalize_latents` — subroutine
- `prep_imgs4projected_model` — subroutine
- `projected_model_kfromto` — function
- `reduce_estep_latent_part_files` — subroutine
- `solve_coupled_basis_exp` — subroutine
- `test_projected_latent_canonicalization` — subroutine
- `test_projected_latent_mstep_stats_io` — subroutine
- `test_projected_model_plane_preparation` — subroutine
- `update_basis_from_latents` — subroutine
- `update_basis_from_mstep_stats_part_files` — subroutine
- `write_estep_latent_part_file` — subroutine
- `write_mstep_stats_part_file` — subroutine

Private symbols:
- `accumulate_projected_latent_mstep_2d_block` — subroutine
- `cap_fplane_for_projected_model` — subroutine
- `cleanup_plane` — subroutine
- `cleanup_runtime_batch` — subroutine
- `compute_canonical_transform` — subroutine
- `finalize_basis_for_projection` — subroutine
- `hermitian_plane_inner_product` — function
- `init_projected_latent_estep_part` — subroutine
- `init_projected_latent_mstep_2d_block` — subroutine
- `init_projected_latent_mstep_stats` — subroutine
- `insert_projected_latent_mstep_2d_block` — subroutine
- `kill_projected_latent_estep_part` — subroutine
- `kill_projected_latent_mstep_2d_block` — subroutine
- `kill_projected_latent_mstep_stats` — subroutine
- `log_comp_seconds` — subroutine
- `log_seconds` — subroutine
- `mix_projected_latent_basis` — subroutine
- `plane_energy` — function
- `prepare_projected_latent_estep_part` — subroutine
- `prepare_projected_latent_mstep_2d_block` — subroutine
- `projected_latent_estep_part` — type
- `projected_latent_mstep_2d_block` — type
- `projected_latent_mstep_stats` — type
- `projected_model_log_prefix` — function
- `read_particles` — subroutine
- `read_projected_latent_estep_part` — subroutine
- `reduce_projected_latent_estep_part` — subroutine
- `reduce_projected_latent_mstep_stats_file` — subroutine
- `regularize_basis_volume` — subroutine
- `reseed_latent_column` — subroutine
- `reset_projected_latent_estep_part` — subroutine
- `reset_projected_latent_mstep_2d_block` — subroutine
- `solve_latent_least_squares` — subroutine
- `solve_real_spd_complex` — subroutine
- `subtract_plane` — subroutine
- `subtract_scaled_plane` — subroutine
- `write_projected_latent_estep_part` — subroutine
- `write_projected_latent_mstep_stats` — subroutine

---
## Module: simple_flex_reconstructor_latent_ops

Files:
- `main/flex/simple_flex_reconstructor_latent_ops.f90`

Uses:
- `simple_core_module_api`
- `simple_ftiter`
- `simple_math`
- `simple_parameters`
- `simple_reconstructor`
- `simple_sp_project`

Public symbols:
- `accumulate_plane_oversamp_coupled_stats` — subroutine
- `accumulate_planes_oversamp_coupled_stats_batch` — subroutine
- `insert_plane_oversamp_coupled_scaled` — subroutine
- `insert_plane_oversamp_multi_scaled` — subroutine
- `insert_planes_oversamp_coupled_batch_scaled` — subroutine
- `project_fplane_mean` — subroutine
- `project_fplanes_mean_basis` — subroutine
- `test_cartesian_projection_contract` — subroutine
- `test_coupled_batch_accumulation` — subroutine

Private symbols:
- `cleanup_test_plane` — subroutine
- `ensure_latent_projection_plane` — subroutine
- `kb_apod_vecs_3d_fast` — subroutine
- `kb_apod_vecs_3d_fast` — subroutine
- `kb_apod_vecs_3d_fast` — subroutine
- `kb_apod_vecs_3d_fast` — subroutine
- `kb_apod_vecs_3d_fast` — subroutine
- `latent_projection_weights` — subroutine
- `weighted_expanded_cmat` — function

---
## Module: simple_forked_process

Files:
- `utils/simple_forked_process.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_error`
- `simple_fileio`
- `simple_memory_monitor`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `unix`

Public symbols:
- `forked_process` — type

Private symbols:
- `await_final_status` — subroutine
- `destroy` — subroutine
- `execute_test` — subroutine
- `get_failtime` — function
- `get_nrestarts` — function
- `get_pid` — function
- `get_queuetime` — function
- `get_starttime` — function
- `get_stoptime` — function
- `kill` — subroutine
- `sigterm_handler` — subroutine
- `skip` — subroutine
- `start` — subroutine
- `status` — function
- `terminate` — subroutine

---
## Module: simple_forked_process_tester

Files:
- `utils/simple_forked_process_tester.f90`

Uses:
- `simple_forked_process`
- `simple_string`
- `simple_syslib`
- `simple_test_utils`
- `unix`

Public symbols:
- `run_all_forked_process_tests` — subroutine

Private symbols:
- `test_destroy` — subroutine
- `test_fail_timestamps` — subroutine
- `test_kill` — subroutine
- `test_logfile_redirection` — subroutine
- `test_restart` — subroutine
- `test_start` — subroutine
- `test_terminate` — subroutine
- `test_timestamps` — subroutine

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
- `add_crossing_bars` — subroutine
- `add_curve` — subroutine
- `add_curve_segment` — subroutine
- `add_filled_band` — subroutine
- `add_filled_band_segment` — subroutine
- `add_filled_rect` — subroutine
- `add_invisible_bounds` — subroutine
- `add_segment` — subroutine
- `bar_edges` — subroutine
- `calc_fsc_area_score` — subroutine
- `fibonacci_sphere_dirs` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `score_cfsc_curve` — subroutine

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
- `simple_image`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `ftexp_shsrch` — type
- `test_ftexp_shsrch` — subroutine
- `test_ftexp_shsrch2` — subroutine

Private symbols:
- `ftexp_shsrch_cost_8` — function
- `ftexp_shsrch_fdfcost_8` — subroutine
- `ftexp_shsrch_gcost_8` — subroutine
- `ftexp_shsrch_kill` — subroutine
- `ftexp_shsrch_minimize` — function
- `ftexp_shsrch_new` — subroutine
- `profile_corrs` — subroutine
- `set_factr_pgtol` — subroutine
- `set_shsrch_tol` — subroutine
- `test_shifted_correlator` — subroutine

---
## Module: simple_ftiter

Files:
- `main/image/simple_ftiter.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `ftiter` — type

---
## Module: simple_gauss2Dfit

Files:
- `utils/simple_gauss2Dfit.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

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
## Module: simple_gpu_utils

Files:
- `utils/simple_gpu_utils.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `cuda_device_prop` — type

---
## Module: simple_gridding

Files:
- `main/interp/simple_gridding.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_projector`

Public symbols:
- `prep2D_inv_instrfun4mul` — function
- `prep3D_inv_instrfun4mul` — function

---
## Module: simple_gui_assembler

Files:
- `utils/gui/simple_gui_assembler.f90`

Uses:
- `simple_forked_process`
- `simple_gui_metadata_api`
- `simple_string`
- `unix`

Public symbols:
- `gui_assembler` — type

Private symbols:
- `assemble_stream_heartbeat` — subroutine
- `assemble_stream_initial_picking` — subroutine
- `assemble_stream_opening2D` — subroutine
- `assemble_stream_optics_assignment` — subroutine
- `assemble_stream_particle_sieving` — subroutine
- `assemble_stream_pool2D` — subroutine
- `assemble_stream_preprocess` — subroutine
- `assemble_stream_reference_picking` — subroutine
- `clear_hashes` — subroutine
- `forked_process_status` — subroutine
- `is_associated` — function
- `kill` — subroutine
- `new` — subroutine
- `set_stoptime` — subroutine
- `to_string` — function

---
## Module: simple_gui_assembler_tester

Files:
- `utils/gui/simple_gui_assembler_tester.f90`

Uses:
- `simple_gui_assembler`
- `simple_gui_metadata_api`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_gui_assembler_tests` — subroutine

Private symbols:
- `test_clear_hashes` — subroutine
- `test_initial_picking` — subroutine
- `test_new_kill` — subroutine
- `test_new_kill_reuse` — subroutine
- `test_opening2D` — subroutine
- `test_optics_assignment` — subroutine
- `test_particle_sieving` — subroutine
- `test_pool2D` — subroutine
- `test_preprocess` — subroutine
- `test_reference_picking` — subroutine
- `test_set_stoptime` — subroutine

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
- `simple_gui_metadata_stream_opening2d`
- `simple_gui_metadata_stream_optics_assignment`
- `simple_gui_metadata_stream_particle_sieving`
- `simple_gui_metadata_stream_picking`
- `simple_gui_metadata_stream_pool2d`
- `simple_gui_metadata_stream_pool2d_snapshot`
- `simple_gui_metadata_stream_preprocess`
- `simple_gui_metadata_stream_update`
- `simple_gui_metadata_timeplot`
- `simple_gui_metadata_types`

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
- `simple_gui_metadata_types`
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
- `get_max_points` — function
- `jsonise_override` — function
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_opening2D

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_opening2D.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_opening2D` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `set` — subroutine
- `set_user_input` — subroutine

---
## Module: simple_gui_metadata_stream_optics_assignment

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_optics_assignment.f90`

Uses:
- `json_kinds`
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
- `set` — subroutine
- `set_selection` — subroutine
- `set_user_input` — subroutine

---
## Module: simple_gui_metadata_stream_picking

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_picking.f90`

Uses:
- `json_kinds`
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
- `set` — subroutine
- `set_initial_ref_selection` — subroutine
- `set_user_input` — subroutine

---
## Module: simple_gui_metadata_stream_pool2D_snapshot

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_pool2D_snapshot.f90`

Uses:
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`
- `unix`

Public symbols:
- `gui_metadata_stream_pool2D_snapshot` — type

Private symbols:
- `get` — function
- `jsonise_override` — function
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_preprocess

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_preprocess.f90`

Uses:
- `json_kinds`
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
- `set` — subroutine

---
## Module: simple_gui_metadata_stream_update

Files:
- `utils/gui/metadata/stream/simple_gui_metadata_stream_update.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_defs`
- `simple_error`
- `simple_gui_metadata_base`
- `simple_string`

Public symbols:
- `gui_metadata_stream_update` — type

Private symbols:
- `get_astigmatism_update` — function
- `get_ctfres_update` — function
- `get_icescore_update` — function
- `get_increase_nmics` — function
- `get_mskdiam2D_update` — function
- `get_pickrefs_clusters` — function
- `get_pickrefs_selection` — function
- `get_pickrefs_selection_length` — function
- `get_sieverefs_selection` — function
- `get_sieverefs_selection_length` — function
- `get_snapshot2D_update` — subroutine
- `has_snapshot2D_update` — function
- `set_astigmatism_update` — subroutine
- `set_ctfres_update` — subroutine
- `set_icescore_update` — subroutine
- `set_increase_nmics` — subroutine
- `set_mskdiam2D_update` — subroutine
- `set_pickrefs_clusters` — subroutine
- `set_pickrefs_selection` — subroutine
- `set_pickrefs_selection_length` — subroutine
- `set_sieverefs_selection` — subroutine
- `set_sieverefs_selection_length` — subroutine
- `set_snapshot2D_update` — subroutine

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
- `test_jsonise_cavg2D` — subroutine
- `test_jsonise_histogram` — subroutine
- `test_jsonise_micrograph` — subroutine
- `test_jsonise_optics_group` — subroutine
- `test_jsonise_stream_initial_picking` — subroutine
- `test_jsonise_stream_opening2D` — subroutine
- `test_jsonise_stream_optics_assignment` — subroutine
- `test_jsonise_stream_particle_sieving` — subroutine
- `test_jsonise_stream_pool2D` — subroutine
- `test_jsonise_stream_pool2D_snapshot` — subroutine
- `test_jsonise_stream_preprocess` — subroutine
- `test_jsonise_stream_reference_picking` — subroutine
- `test_jsonise_stream_update` — subroutine
- `test_jsonise_timeplot` — subroutine
- `test_new_kill` — subroutine
- `test_serialise` — subroutine
- `test_serialise_cavg2D` — subroutine
- `test_serialise_histogram` — subroutine
- `test_serialise_micrograph` — subroutine
- `test_serialise_optics_group` — subroutine
- `test_serialise_stream_initial_picking` — subroutine
- `test_serialise_stream_opening2D` — subroutine
- `test_serialise_stream_optics_assignment` — subroutine
- `test_serialise_stream_particle_sieving` — subroutine
- `test_serialise_stream_pool2D` — subroutine
- `test_serialise_stream_pool2D_snapshot` — subroutine
- `test_serialise_stream_preprocess` — subroutine
- `test_serialise_stream_reference_picking` — subroutine
- `test_serialise_stream_update` — subroutine
- `test_serialise_timeplot` — subroutine
- `test_set_get_cavg2D` — subroutine
- `test_set_get_histogram` — subroutine
- `test_set_get_micrograph` — subroutine
- `test_set_get_optics_group` — subroutine
- `test_set_get_stream_initial_picking` — subroutine
- `test_set_get_stream_opening2D` — subroutine
- `test_set_get_stream_optics_assignment` — subroutine
- `test_set_get_stream_particle_sieving` — subroutine
- `test_set_get_stream_pool2D` — subroutine
- `test_set_get_stream_pool2D_snapshot` — subroutine
- `test_set_get_stream_preprocess` — subroutine
- `test_set_get_stream_reference_picking` — subroutine
- `test_set_get_stream_update` — subroutine
- `test_set_get_timeplot` — subroutine

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
- `simple_gui_metadata_api`

Public symbols:
- `max_metadata_size` — function

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
## Module: simple_guistats

Files:
- `utils/gui/simple_guistats.f90`

Uses:
- `json_module`
- `simple_core_module_api`
- `simple_image`

Private symbols:
- `create_section_array` — subroutine
- `deactivate_section` — subroutine
- `delete` — subroutine
- `ensure_section` — subroutine
- `generate_2D_jpeg` — subroutine
- `generate_2D_thumbnail` — subroutine
- `get_1` — subroutine
- `get_keyline` — function
- `hide` — subroutine
- `init` — subroutine
- `kill` — subroutine
- `merge` — subroutine
- `new_section` — subroutine
- `set_1` — subroutine
- `set_2` — subroutine
- `set_3` — subroutine
- `set_now` — subroutine
- `write` — subroutine
- `write_json` — subroutine

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
## Module: simple_hclust

Files:
- `utils/clustering/simple_hclust.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `test_hclust` — subroutine

Private symbols:
- `cluster` — subroutine
- `cut_tree` — subroutine
- `get_medoids` — subroutine
- `kill` — subroutine
- `new` — subroutine

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
- `kill` — subroutine
- `new_1` — subroutine
- `new_2` — subroutine
- `new_3` — subroutine
- `new_4` — subroutine
- `new_5` — subroutine
- `plot` — subroutine
- `quantize` — subroutine
- `reset` — subroutine
- `smooth` — subroutine
- `topdf` — subroutine
- `update` — subroutine
- `zero` — subroutine

---
## Module: simple_http_post

Files:
- `utils/comm/simple_http_post.f90`

Uses:
- `curl`
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
- `simple_http_post`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_http_post_tests` — subroutine

Private symbols:
- `test_create_and_kill` — subroutine
- `test_request_no_body` — subroutine
- `test_request_with_body` — subroutine
- `test_response_reset` — subroutine

---
## Module: simple_image

Files:
- `main/image/simple_image.f90`

Uses:
- `gnufor2`
- `simple_core_module_api`
- `simple_ctf`
- `simple_fftw3`
- `simple_ftiter`
- `simple_imgfile`
- `simple_neighs`
- `simple_winfuns`

Public symbols:
- `image` — type

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
- `dequeue` — subroutine
- `diameter_bin` — subroutine
- `diameter_cc` — subroutine
- `dilate` — subroutine
- `elim_cc` — subroutine
- `elim_ccs` — subroutine
- `elim_largestcc` — subroutine
- `empty_queue` — subroutine
- `enqueue` — subroutine
- `erode` — subroutine
- `feret_minmax` — subroutine
- `find_ccs` — subroutine
- `flood_fill` — subroutine
- `get_imat` — subroutine
- `get_nccs` — subroutine
- `grow_bins` — subroutine
- `inv_bimg` — subroutine
- `is_empty` — function
- `kill_bimg` — subroutine
- `masscen_cc` — subroutine
- `max_dist` — subroutine
- `new_bimg` — subroutine
- `order_ccs` — subroutine
- `polish_ccs` — subroutine
- `read_bimg` — subroutine
- `set_edgecc2background` — subroutine
- `set_imat` — subroutine
- `set_largestcc2background` — subroutine
- `size_ccs` — function
- `transfer2bimg` — subroutine
- `update_img_rmat` — subroutine
- `update_mask_2d` — subroutine
- `update_mask_3d` — subroutine
- `write_bimg` — subroutine

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
- `automask2D` — subroutine
- `automask2D_support_pix` — subroutine
- `density_inoutside_mask` — subroutine
- `image_msk` — type

Private symbols:
- `automask2D_binary_one` — subroutine
- `automask3D_1` — subroutine
- `automask3D_2` — subroutine
- `automask3D_binarize` — subroutine
- `automask3D_filter` — subroutine
- `automask3D_keep_largest_cc` — subroutine
- `binary_imat_to_pix` — subroutine
- `env_rproject` — subroutine
- `estimate_spher_mask_diam` — subroutine

---
## Module: simple_imgarr_utils

Files:
- `utils/simple_imgarr_utils.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_sp_project`
- `simple_stack_io`

---
## Module: simple_imgfile

Files:
- `fileio/simple_imgfile.f90`

Uses:
- `gnufor2`
- `iso_fortran_env`
- `simple_core_module_api`
- `simple_string_utils`
- `simple_tifflib`

Private symbols:
- `close` — subroutine
- `getDim` — function
- `getDims` — function
- `open` — subroutine
- `open_local` — subroutine
- `rSlices` — subroutine
- `rTiffSlices` — subroutine
- `setDims` — subroutine
- `setIform` — subroutine
- `setMean` — subroutine
- `setMinMax` — subroutine
- `setMode` — subroutine
- `setPixSz` — subroutine
- `setRMSD` — subroutine
- `slice2bytepos` — subroutine
- `slice2recpos` — subroutine
- `update_MRC_stats` — subroutine
- `wmrcSlices` — subroutine
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
- `find_ldim_nptcls` — subroutine
- `get_mrcfile_info` — subroutine
- `MrcImgHead` — type
- `SpiImgHead` — type
- `test_imghead` — subroutine
- `TiffImgHead` — type
- `update_stack_nimgs` — subroutine

Private symbols:
- `CloseTiff` — subroutine
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
- `getMode` — function
- `is_even` — function
- `kill` — subroutine
- `new` — subroutine
- `print_imghead` — subroutine
- `print_spihed` — subroutine
- `read` — subroutine
- `read_spihed` — subroutine
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

---
## Module: simple_imgproc

Files:
- `utils/simple_imgproc.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `make_pcavecs` — subroutine
- `make_pcavol` — subroutine

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
- `get_server_ips` — function
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

---
## Module: simple_jiffys

Files:
- `utils/simple_jiffys.f90`

Uses:
- `simple_defs`
- `simple_defs_fname`
- `simple_timer`

---
## Module: simple_jpg

Files:
- `extlibs/jpg/simple_jpg.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `test_jpg_export` — subroutine

Private symbols:
- `constructor` — subroutine
- `load_jpeg_i4` — function
- `load_jpeg_r4` — function
- `montage` — function
- `read_jpeg` — subroutine
- `save_jpeg_i4` — function
- `save_jpeg_r4` — function
- `save_jpeg_r4_3D` — function
- `setColorspace` — subroutine
- `setGamma` — subroutine
- `setHeight` — subroutine
- `setNormalised` — subroutine
- `setQuality` — subroutine
- `setup_jpeg` — subroutine
- `setWidth` — subroutine
- `write_jpeg` — subroutine

---
## Module: simple_kbinterpol

Files:
- `main/interp/simple_kbinterpol.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_edges_sqwins`

Public symbols:
- `kbinterpol` — type

---
## Module: simple_kd_tree

Files:
- `utils/structs/simple_kd_tree.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `kd_tree` — type
- `knn_table` — type

Private symbols:
- `build_kd_tree` — subroutine
- `build_node` — subroutine
- `grow_nodes` — subroutine
- `heap_consider` — subroutine
- `heap_sift_down` — subroutine
- `kill_kd_tree` — subroutine
- `kill_knn_table` — subroutine
- `query_all_kd_tree` — subroutine
- `query_kd_tree` — subroutine
- `query_node` — subroutine
- `reserve_nodes` — subroutine
- `select_kth` — subroutine
- `sort_neighbor_pairs` — subroutine
- `swap_int` — subroutine
- `swap_pair` — subroutine
- `validate_query` — subroutine

---
## Module: simple_kmedoids

Files:
- `utils/clustering/simple_kmedoids.f90`

Uses:
- `simple_core_module_api`

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
- `new_1` — subroutine
- `new_2` — subroutine

---
## Module: simple_kpca_svd

Files:
- `main/pca/simple_kpca_svd.f90`

Uses:
- `simple_core_module_api`
- `simple_pca`

Public symbols:
- `kpca_svd` — type

Private symbols:
- `center_columns` — subroutine
- `choose_nystrom_inds_from_data` — subroutine
- `compute_eigvecs` — subroutine
- `cosine_kernel` — subroutine
- `dense_mm` — subroutine
- `dense_tmm` — subroutine
- `generate_kpca` — subroutine
- `get_eigvals_kpca` — function
- `get_feat_kpca` — function
- `gram_symmetric` — subroutine
- `kernel_center` — subroutine
- `kill_kpca` — subroutine
- `master_kpca` — subroutine
- `master_nystrom` — subroutine
- `new_kpca` — subroutine
- `orthonormalize_cols` — subroutine
- `partial_eigh_sym` — subroutine
- `projected_kernel_col` — subroutine
- `rbf_kernel` — subroutine
- `select_local_support_inds` — subroutine
- `select_nystrom_inds` — subroutine
- `set_params_kpca` — subroutine

---
## Module: simple_linalg

Files:
- `utils/math/simple_linalg.f90`

Uses:
- `simple_defs`
- `simple_is_check_assert`

Private symbols:
- `sparse_matvec_sp_proc` — subroutine

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
- `front` — subroutine
- `getter` — subroutine
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
## Module: simple_magic_boxes

Files:
- `utils/simple_magic_boxes.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_math`

---
## Module: simple_make_cavgs_strategy

Files:
- `main/strategies/parallelization/simple_make_cavgs_strategy.f90`

Uses:
- `simple_builder`
- `simple_classaverager`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_exec_helpers`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`

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
- `split_nobjs_even` — function
- `split_pairs_in_parts` — subroutine

---
## Module: simple_matcher_2Dprep

Files:
- `main/strategies/search/simple_matcher_2Dprep.f90`

Uses:
- `simple_builder`
- `simple_butterworth`
- `simple_matcher_ptcl_io`
- `simple_matcher_smpl_and_lplims`
- `simple_pftc_srch_api`
- `simple_projector`
- `simple_strategy2d_utils`
- `simple_timer`

Public symbols:
- `calc_2Dref_offset` — subroutine
- `prep2Dref` — subroutine
- `prepimg4align` — subroutine

---
## Module: simple_matcher_3Drec

Files:
- `main/strategies/search/simple_matcher_3Drec.f90`

Uses:
- `simple_builder`
- `simple_classaverager`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_image`
- `simple_imgarr_utils`
- `simple_matcher_ptcl_io`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_reconstructor`
- `simple_refine3d_fnames`
- `simple_timer`

Public symbols:
- `calc_3Drec` — subroutine
- `calc_projdir3Drec` — subroutine
- `cleanup_rec_buffers` — subroutine
- `init_rec` — subroutine
- `prep_imgs4rec` — subroutine
- `write_state_partial` — subroutine

Private symbols:
- `group_pinds_by_state_eo` — subroutine
- `init_state_half_rec` — subroutine
- `kill_state_half_rec` — subroutine
- `mark_empty_state` — subroutine
- `set_state_vol_output` — subroutine
- `update_state_half_rec` — subroutine
- `write_state_half_partial` — subroutine

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
- `simple_parameters`
- `simple_strategy2d_alloc`

Public symbols:
- `prep_pftc4align2D` — subroutine

---
## Module: simple_matcher_ptcl_batch

Files:
- `main/strategies/search/simple_matcher_ptcl_batch.f90`

Uses:
- `simple_builder`
- `simple_euclid_sigma2`
- `simple_imgarr_utils`
- `simple_matcher_2dprep`
- `simple_matcher_ptcl_io`
- `simple_pftc_srch_api`

Public symbols:
- `alloc_ptcl_imgs` — subroutine
- `build_batch_particles2D` — subroutine
- `build_batch_particles3D` — subroutine
- `clean_batch_particles2D` — subroutine
- `clean_batch_particles3D` — subroutine
- `prep_sigmas_objfun` — subroutine

Private symbols:
- `polarize_batch_particles3D` — subroutine
- `polarize_batch_particles3D_den` — subroutine

---
## Module: simple_matcher_ptcl_io

Files:
- `main/strategies/search/simple_matcher_ptcl_io.f90`

Uses:
- `simple_builder`
- `simple_discrete_stack_io`
- `simple_imghead`
- `simple_pftc_srch_api`

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

Public symbols:
- `adopt_reprojection_model_range` — subroutine
- `estimate_lp_from_refs` — subroutine
- `materialize_reprojection_model_from_volumes` — subroutine
- `pick_lp_est_state` — subroutine
- `read_mask_filter_reproject_refvols` — subroutine
- `read_reprojection_model` — subroutine
- `remove_ref_section_files` — subroutine

Private symbols:
- `blend_lowres_eo_for_registration` — subroutine
- `calcrefvolshift_and_mapshifts2ptcls` — subroutine
- `read_mask_filter_refvols` — subroutine
- `read_reprojection_model_file_header` — subroutine
- `read_reprojection_model_header` — subroutine
- `regularize_ref_with_noise` — subroutine

---
## Module: simple_matcher_smpl_and_lplims

Files:
- `main/strategies/search/simple_matcher_smpl_and_lplims.f90`

Uses:
- `simple_builder`
- `simple_pftc_srch_api`
- `simple_refine3d_fnames`

Public symbols:
- `sample_ptcls4fillin` — subroutine
- `sample_ptcls4missing3D` — subroutine
- `sample_ptcls4update2D` — subroutine
- `sample_ptcls4update3D` — subroutine
- `set_bp_range2D` — subroutine
- `set_bp_range3D` — subroutine

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

---
## Module: simple_math_ctf

Files:
- `utils/math/simple_math_ctf.f90`

Uses:
- `simple_defs`
- `simple_memoize_ft_maps`

---
## Module: simple_math_ft

Files:
- `utils/math/simple_math_ft.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_is_check_assert`
- `simple_srch_sort_loc`

---
## Module: simple_mem_estimator

Files:
- `utils/simple_mem_estimator.f90`

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
- `utils/simple_memoize_ft_maps.f90`

Uses:
- `simple_core_module_api`
- `simple_ftiter`

Public symbols:
- `forget_ft_maps` — subroutine
- `ft_map_get_farray_shape` — function
- `memoize_ft_maps` — subroutine

---
## Module: simple_memory_monitor

Files:
- `utils/simple_memory_monitor.f90`

Uses:
- `simple_cmdline`
- `simple_defs`
- `simple_string`
- `simple_syslib`

Public symbols:
- `mem_monitor_finish` — subroutine
- `mem_monitor_init` — subroutine
- `mem_report` — subroutine

Private symbols:
- `log_sample` — subroutine
- `start_monitor` — subroutine

---
## Module: simple_micproc

Files:
- `utils/simple_micproc.f90`

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
- `cluster` — subroutine
- `clustering_rejection` — subroutine
- `flag_amorphous_carbon` — subroutine
- `flag_ice` — subroutine
- `read_mic` — subroutine
- `read_mic_subtr_backgr` — subroutine
- `read_mic_subtr_backgr_shrink` — subroutine
- `sample_filetab` — function

---
## Module: simple_micrograph_generator

Files:
- `main/simple_micrograph_generator.f90`

Uses:
- `simple_core_module_api`
- `simple_eer_factory`
- `simple_image`
- `simple_motion_correct_utils`
- `simple_ori`
- `simple_starfile_wrappers`
- `simple_string`

Public symbols:
- `mic_generator` — type

Private symbols:
- `cure_outliers_1` — subroutine
- `cure_outliers_2` — subroutine
- `display` — subroutine
- `generate_micrographs` — subroutine
- `get_moviename` — function
- `kill` — subroutine
- `new` — subroutine
- `parse_movie_metadata` — subroutine
- `parse_string` — subroutine
- `write_star` — subroutine

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
- `simple_sp_project`

Public symbols:
- `ensure_real_capacity` — subroutine
- `print_diam_stats` — subroutine
- `segdiampick_mics` — subroutine
- `segdiampick_mics_multi` — subroutine
- `segdiampick_preprocess` — subroutine

---
## Module: simple_molecule_data

Files:
- `main/nano/simple_molecule_data.f90`

Public symbols:
- `betagal_1jyx` — function
- `molecule_data` — type
- `sars_cov2_spkgp_6vxx` — function

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
- `calc_shifts` — subroutine
- `dealloc_ftexps` — subroutine
- `dealloc_images` — subroutine
- `fit_polynomial` — subroutine
- `gen_weights` — subroutine
- `get_coords` — subroutine
- `get_corrs` — subroutine
- `get_opt_shifts` — subroutine
- `get_shifts_toplot` — subroutine
- `get_weights` — subroutine
- `init_ftexps` — subroutine
- `init_images` — subroutine
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
- `set_maxits` — subroutine
- `set_rand_init_shifts` — subroutine
- `set_reslims` — subroutine
- `set_shsrch_tol` — subroutine
- `set_trs` — subroutine
- `set_weights` — subroutine

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
- `simple_motion_patched`
- `simple_opt_image_weights`
- `simple_parameters`
- `simple_starfile_wrappers`

Public symbols:
- `motion_correct_calc_bid` — subroutine
- `motion_correct_calc_opt_weights` — subroutine
- `motion_correct_iso` — subroutine
- `motion_correct_iso_calc_sums` — subroutine
- `motion_correct_iso_kill` — subroutine
- `motion_correct_iso_shift_frames` — subroutine
- `motion_correct_kill_common` — subroutine
- `motion_correct_mic2spec` — subroutine
- `motion_correct_patched` — subroutine
- `motion_correct_patched_calc_sums` — subroutine
- `motion_correct_patched_kill` — subroutine
- `motion_correct_write2star` — subroutine
- `motion_correct_write_poly` — subroutine

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
- `calc_eer_fraction` — subroutine
- `correct_gain` — subroutine
- `flip_gain` — subroutine

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

---
## Module: simple_motion_gain_helpers

Files:
- `main/motion/simple_motion_gain_helpers.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `normalized_inverse_average_intensity` — subroutine
- `read_movies_and_sum_frames` — subroutine

---
## Module: simple_motion_gain_tester

Files:
- `main/motion/simple_motion_gain_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_motion_gain_analysis`
- `simple_motion_gain_helpers`
- `simple_test_utils`

Public symbols:
- `run_all_motion_gain_tests` — subroutine

Private symbols:
- `create_constant_image` — subroutine
- `create_gainref` — subroutine
- `create_movie_stack` — subroutine
- `test_gain_flip_analyzer_batch_updates` — subroutine
- `test_read_movies_and_sum_frames_counts` — subroutine

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
- `apply_patch_poly` — function
- `correct` — subroutine
- `deallocate_patches` — subroutine
- `det_shifts` — subroutine
- `det_shifts_refine` — subroutine
- `fit_polynomial` — subroutine
- `gen_micrograph` — subroutine
- `gen_patch` — subroutine
- `get_local_shift` — subroutine
- `get_poly4star` — subroutine
- `get_poly_coeffs` — subroutine
- `get_polyfit_rmsd` — function
- `kill` — subroutine
- `new` — subroutine
- `patch_poly` — function
- `pix2polycoords` — subroutine
- `plot_shifts` — subroutine
- `robust_fit_polynomial` — subroutine
- `set_bfactor` — subroutine
- `set_fitshifts` — subroutine
- `set_fixed_frame` — subroutine
- `set_frameweights` — subroutine
- `set_interp_fixed_frame` — subroutine
- `set_nframes` — subroutine
- `set_poly_coeffs` — subroutine
- `set_size_frames_ref` — subroutine

---
## Module: simple_multi_dendro

Files:
- `utils/structs/simple_multi_dendro.f90`

Uses:
- `simple_binary_tree`
- `simple_core_module_api`
- `simple_hclust`
- `simple_srch_sort_loc`

Public symbols:
- `deserialize_multi_dendro` — subroutine
- `multi_dendro` — type
- `serialize_multi_dendro` — subroutine

Private symbols:
- `build_tree_balanced` — subroutine
- `build_tree_from_subdistmat` — subroutine
- `get_full2sub_map` — function
- `get_node` — function
- `get_root_node` — function
- `get_sub2full_map` — function
- `get_tree_refs` — function
- `get_tree_refs_static` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `new_multi_balanced` — subroutine

---
## Module: simple_multi_dendro_tester

Files:
- `utils/structs/simple_multi_dendro_tester.f90`

Uses:
- `simple_binary_tree`
- `simple_multi_dendro`
- `simple_test_utils`

Public symbols:
- `run_all_multi_dendro_tests` — subroutine

Private symbols:
- `assert_int_in_set` — subroutine
- `build_tree_from_full_dist` — subroutine
- `make_distmat_two_pairs` — subroutine
- `make_distmat_two_triples` — subroutine
- `test_build_tree_balanced_two_triples` — subroutine
- `test_build_tree_pairs_basic_structure` — subroutine
- `test_build_tree_singletons` — subroutine
- `test_build_tree_two_triples_basic_structure` — subroutine
- `test_get_tree_refs` — subroutine
- `test_initial_state` — subroutine
- `test_kill_and_reuse` — subroutine
- `test_new_counts_and_getters` — subroutine
- `test_new_multi_balanced_single_tree` — subroutine

---
## Module: simple_nanoparticle

Files:
- `main/nano/simple_nanoparticle.f90`

Uses:
- `simple_atoms`
- `simple_core_module_api`
- `simple_defs_atoms`
- `simple_image`
- `simple_image_bin`
- `simple_image_msk`
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
- `calc_cn_stats` — subroutine
- `calc_isotropic_disp` — subroutine
- `calc_longest_atm_dist` — subroutine
- `calc_zscore` — subroutine
- `center_on_atom` — subroutine
- `check_neighbors_cn` — subroutine
- `conv_denoise` — subroutine
- `discard_atoms` — subroutine
- `discard_small_ccs` — subroutine
- `fillin_atominfo` — subroutine
- `find_centers` — subroutine
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
- `pack_instance4stats` — subroutine
- `remove_lowly_contacted` — subroutine
- `remove_small_and_lowly_contacted` — subroutine
- `set_atm_info` — subroutine
- `set_atom` — subroutine
- `set_atomic_coords_from_pdb` — subroutine
- `set_atomic_coords_from_xyz` — subroutine
- `set_coords4stats` — subroutine
- `set_img` — subroutine
- `set_ncc` — subroutine
- `simulate_atoms` — subroutine
- `split_atom` — subroutine
- `split_atoms` — subroutine
- `update_ncc` — subroutine
- `validate_atoms` — subroutine
- `write_2D_slice` — subroutine
- `write_atominfo` — subroutine
- `write_centers_1` — subroutine
- `write_centers_2` — subroutine
- `write_centers_aniso` — subroutine
- `write_cn_atoms` — subroutine
- `write_cn_stats` — subroutine
- `write_csv_files` — subroutine
- `write_individual_atoms` — subroutine
- `write_np_stats` — subroutine

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
- `calc_contact_scores` — subroutine
- `dists_btw_common` — function
- `find_atoms_subset` — subroutine
- `find_couples` — subroutine
- `find_rMax` — function
- `fit_lattice` — subroutine
- `Kabsch_algo` — subroutine
- `phasecorr_one_atom` — subroutine
- `read_pdb2matrix` — subroutine
- `remove_atoms` — subroutine
- `run_cn_analysis` — subroutine
- `strain_analysis` — subroutine
- `thres_detect_conv_atom_denoised` — subroutine
- `write_matrix2pdb` — subroutine

Private symbols:
- `gaussian_kernel` — subroutine
- `nnbdl_lat` — subroutine
- `qr` — subroutine
- `write_ideal_lattice_pdb` — subroutine

---
## Module: simple_neighs

Files:
- `utils/math/simple_neighs.f90`

Uses:
- `simple_core_module_api`

---
## Module: simple_nice

Files:
- `utils/gui/simple_nice.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_core_module_api`
- `simple_histogram`
- `simple_socket_comm`
- `simple_sp_project`
- `unix`

Public symbols:
- `add_headline_cls2D` — subroutine
- `add_headline_micrographs` — subroutine
- `add_histograms_cls2D` — subroutine
- `add_histograms_micrographs` — subroutine
- `add_plot_optics` — subroutine
- `add_view_cls2D` — subroutine
- `add_view_ini3D` — subroutine
- `add_view_micrographs` — subroutine
- `add_view_optics` — subroutine
- `add_view_pick` — subroutine
- `add_view_vols` — subroutine
- `bar_plot_object` — subroutine
- `calculate_checksum` — function
- `create_interactive_plot` — subroutine
- `cycle` — subroutine
- `datestr` — function
- `doughnut_plot_object` — subroutine
- `fsc_object` — subroutine
- `generate_stat_json` — subroutine
- `get_real_keys_json` — subroutine
- `get_stat_json` — subroutine
- `image_thumbnail_object` — subroutine
- `init_1` — subroutine
- `init_2` — subroutine
- `nice_plot_bar` — type
- `nice_plot_doughnut` — type
- `nice_stat_root` — type
- `nice_stat_thumb_image` — type
- `nice_stat_thumb_image_meta` — type
- `nice_thread_comm_message` — type
- `nice_view_cls2D` — type
- `nice_view_ini3D` — type
- `nice_view_micrographs` — type
- `nice_view_optics` — type
- `nice_view_pick` — type
- `nice_view_vols` — type
- `remote_heartbeat` — subroutine
- `simple_nice_comm` — type
- `start_comm_thread` — subroutine
- `terminate` — subroutine
- `terminate_comm_thread` — subroutine
- `text_data_object_1` — subroutine
- `text_data_object_2` — subroutine
- `text_data_object_3` — subroutine
- `update_cls2D` — subroutine
- `update_from_project` — subroutine
- `update_ini3D` — subroutine
- `update_micrographs` — subroutine
- `update_optics` — subroutine
- `update_pick` — subroutine
- `update_vols` — subroutine

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

Uses:
- `simple_butterworth`
- `simple_core_module_api`
- `simple_image`
- `simple_neighs`
- `simple_tent_smooth`

Private symbols:
- `nu_highres_extension_stats` — type

---
## Module: simple_online_var

Files:
- `utils/math/simple_online_var.f90`

Uses:
- `simple_defs`
- `simple_rnd`
- `simple_stat`

Public symbols:
- `online_var` — type
- `test_online_var` — subroutine

Private symbols:
- `add_1` — subroutine
- `add_2` — subroutine
- `get_mean` — function
- `get_var` — function
- `reset_mean` — subroutine
- `serialize` — function
- `unserialize` — subroutine

---
## Module: simple_opt_bfgs

Files:
- `main/opt/simple_opt_bfgs.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`
- `simple_opt_subs`
- `simple_optimizer`

Public symbols:
- `opt_bfgs` — type

Private symbols:
- `bfgs_minimize` — subroutine
- `bfgsmin` — subroutine
- `kill_opt_bfgs` — subroutine
- `new_opt_bfgs` — subroutine

---
## Module: simple_opt_bfgs2

Files:
- `main/opt/simple_opt_bfgs2.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_helpers`
- `simple_opt_spec`
- `simple_opt_subs`
- `simple_optimizer`

Public symbols:
- `opt_bfgs2` — type

Private symbols:
- `bfgs2_iterate` — function
- `bfgs2_minimize` — subroutine
- `bfgs2_set` — subroutine
- `bfgs2_wrapper` — type
- `bfgs2min` — subroutine
- `change_direction` — subroutine
- `check_extremum` — subroutine
- `cubic` — function
- `interp_cubic` — function
- `interp_quad` — function
- `interpolate` — function
- `kill_opt_bfgs2` — subroutine
- `linear_minimize` — function
- `moveto` — subroutine
- `new_opt_bfgs2` — subroutine
- `prepare_wrapper` — subroutine
- `slope` — function
- `solve_quadratic` — function
- `update_position` — subroutine
- `wrap_df` — function
- `wrap_f` — function
- `wrap_fdf` — subroutine

---
## Module: simple_opt_bforce

Files:
- `main/opt/simple_opt_bforce.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `opt_bforce` — type

Private symbols:
- `bforce_minimize` — subroutine
- `kill_opt_bforce` — subroutine
- `new_opt_bforce` — subroutine
- `srch_not_done` — function

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
- `simple_opt_bfgs2`
- `simple_opt_bforce`
- `simple_opt_de`
- `simple_opt_lbfgsb`
- `simple_opt_simplex`
- `simple_opt_spec`
- `simple_opt_stde`
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

---
## Module: simple_opt_fr_cg

Files:
- `main/opt/simple_opt_fr_cg.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_helpers`
- `simple_opt_spec`
- `simple_opt_subs`
- `simple_optimizer`

Public symbols:
- `opt_fr_cg` — type

Private symbols:
- `fr_cg_iterate` — function
- `fr_cg_minimize` — subroutine
- `fr_cg_set` — subroutine
- `fr_cgmin` — subroutine
- `kill_opt_fr_cg` — subroutine
- `new_opt_fr_cg` — subroutine

---
## Module: simple_opt_helpers

Files:
- `main/opt/simple_opt_helpers.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`

Public symbols:
- `intermediate_point` — subroutine
- `minimize` — subroutine
- `take_step` — subroutine
- `test_gradient` — function

---
## Module: simple_opt_image_weights

Files:
- `main/motion/simple_opt_image_weights.f90`

Uses:
- `simple_core_module_api`
- `simple_ft_expanded`
- `simple_ftexp_shsrch`
- `simple_image`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `opt_image_weights` — type

Private symbols:
- `calc_Dmat` — subroutine
- `calc_opt_weights` — subroutine
- `create_ftexp_objs` — subroutine
- `dealloc_ftexp_objs` — subroutine
- `get_weights` — function
- `kill` — subroutine
- `new` — subroutine
- `opt_image_weights_cost_Dmat` — function
- `opt_image_weights_cost_Ref` — function
- `opt_image_weights_cost_sp_wrapper` — function
- `opt_image_weights_cost_wrapper` — function
- `opt_image_weights_fdf_Dmat` — subroutine
- `opt_image_weights_fdf_Ref` — subroutine
- `opt_image_weights_fdf_wrapper` — subroutine
- `opt_image_weights_gcost` — subroutine
- `opt_image_weights_gcost_Dmat` — subroutine
- `opt_image_weights_gcost_Ref` — subroutine
- `opt_image_weights_gcost_wrapper` — subroutine
- `set_Nrestarts` — subroutine

---
## Module: simple_opt_lbfgsb

Files:
- `main/opt/simple_opt_lbfgsb.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_helpers`
- `simple_opt_spec`
- `simple_opt_subs`
- `simple_optimizer`

Public symbols:
- `opt_lbfgsb` — type

Private symbols:
- `active` — subroutine
- `bmv` — subroutine
- `cauchy` — subroutine
- `cmprlb` — subroutine
- `daxpy` — subroutine
- `dcopy` — subroutine
- `dcsrch` — subroutine
- `dcstep` — subroutine
- `ddot` — function
- `dpofa` — subroutine
- `dscal` — subroutine
- `dtrsl` — subroutine
- `errclb` — subroutine
- `formk` — subroutine
- `formt` — subroutine
- `freev` — subroutine
- `hpsolb` — subroutine
- `kill_opt_lbfgsb` — subroutine
- `lbfgsb_minimize` — subroutine
- `lbfgsbmin` — subroutine
- `lnsrlb` — subroutine
- `mainlb` — subroutine
- `matupd` — subroutine
- `new_opt_lbfgsb` — subroutine
- `projgr` — subroutine
- `setulb` — subroutine
- `subsm` — subroutine

---
## Module: simple_opt_mask

Files:
- `utils/simple_opt_mask.f90`

Uses:
- `simple_core_module_api`
- `simple_image`

Public symbols:
- `estimate_spher_mask` — subroutine

---
## Module: simple_opt_particle_swarm

Files:
- `main/opt/simple_opt_particle_swarm.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `opt_particle_swarm` — type

Private symbols:
- `init` — subroutine
- `kill_particle_swarm` — subroutine
- `new_particle_swarm` — subroutine
- `particle_swarm_minimize` — subroutine
- `update_particle` — subroutine

---
## Module: simple_opt_powell

Files:
- `main/opt/simple_opt_powell.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`
- `simple_opt_subs`
- `simple_optimizer`

Public symbols:
- `opt_powell` — type

Private symbols:
- `kill_opt_powell` — subroutine
- `new_opt_powell` — subroutine
- `powell` — subroutine
- `powell_minimize` — subroutine

---
## Module: simple_opt_pr_cg

Files:
- `main/opt/simple_opt_pr_cg.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_helpers`
- `simple_opt_spec`
- `simple_opt_subs`
- `simple_optimizer`

Public symbols:
- `opt_pr_cg` — type

Private symbols:
- `kill_opt_pr_cg` — subroutine
- `new_opt_pr_cg` — subroutine
- `pr_cg_iterate` — function
- `pr_cg_minimize` — subroutine
- `pr_cg_set` — subroutine
- `pr_cgmin` — subroutine

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
- `init` — subroutine
- `kill_opt_simplex` — subroutine
- `new_opt_simplex` — subroutine
- `simplex_minimize` — subroutine

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
- `fun` — function
- `fun` — subroutine
- `fun` — subroutine
- `fun` — function
- `fun` — subroutine
- `fun` — subroutine
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

---
## Module: simple_opt_stde

Files:
- `main/opt/simple_opt_stde.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_helpers`
- `simple_opt_spec`
- `simple_opt_subs`
- `simple_optimizer`

Public symbols:
- `opt_stde` — type

Private symbols:
- `kill_opt_stde` — subroutine
- `new_opt_stde` — subroutine
- `stde_iterate` — function
- `stde_minimize` — subroutine
- `stde_set` — subroutine
- `stdemin` — subroutine

---
## Module: simple_opt_subs

Files:
- `main/opt/simple_opt_subs.f90`

Uses:
- `simple_core_module_api`
- `simple_opt_spec`

Public symbols:
- `amoeba` — subroutine
- `check_and_correct_vec` — subroutine
- `linmin` — subroutine
- `lnsrch` — subroutine
- `shc_selector` — function

Private symbols:
- `amoeba_private` — subroutine
- `amotry` — function
- `brent` — function
- `eval_move` — function
- `mnbrak` — subroutine

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

Public symbols:
- `ori` — type

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
- `simple_defs`
- `simple_ori`
- `simple_ori_utils`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_ori_tests` — subroutine

Private symbols:
- `test_2Dshift_set_get_and_rounding` — subroutine
- `test_copy_nonptcl_to_ptcl` — subroutine
- `test_copy_ptcl_to_nonptcl` — subroutine
- `test_copy_ptcl_to_ptcl` — subroutine
- `test_ctor_and_flags` — subroutine
- `test_delete_2Dclustering` — subroutine
- `test_delete_3Dalignment` — subroutine
- `test_euler_basic_set_get` — subroutine
- `test_euler_compose_vs_compeuler` — subroutine
- `test_euler_matrix_roundtrip` — subroutine
- `test_geodesic_metrics` — subroutine
- `test_get_dfx_dfy_and_setters` — subroutine
- `test_isthere_and_ischar` — subroutine
- `test_mirror2d_involution` — subroutine
- `test_mirror3d_involution` — subroutine
- `test_normal_and_mat_consistency` — subroutine
- `test_ori2str_and_str2ori_roundtrip` — subroutine
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
- `euler2m` — function
- `euler_compose` — subroutine
- `euler_mirror` — subroutine
- `m2euler` — function
- `make_transfmat` — function
- `rnd_romat` — subroutine

Private symbols:
- `drotmat` — function
- `euler2dm` — subroutine
- `euler_normal` — function
- `rotmat` — function

---
## Module: simple_oris

Files:
- `main/ori/simple_oris.f90`

Uses:
- `simple_ori`
- `simple_ori_api`

Public symbols:
- `oris` — type

---
## Module: simple_oris_tester

Files:
- `main/ori/simple_oris_tester.f90`

Uses:
- `simple_core_module_api`
- `simple_test_utils`

Public symbols:
- `run_all_oris_tests` — subroutine

Private symbols:
- `test_compress_and_masks` — subroutine
- `test_constructors_and_basic_props` — subroutine
- `test_extract_and_copy` — subroutine
- `test_getters_setters` — subroutine
- `test_misc_flags` — subroutine
- `test_proj_space_and_remap` — subroutine
- `test_randomization_and_symmetry` — subroutine
- `test_rotations_and_errors` — subroutine
- `test_sample4update_cnt_large` — subroutine
- `test_sampling_and_updatecnt` — subroutine
- `test_stats_and_ordering` — subroutine

---
## Module: simple_parameters

Files:
- `main/params/simple_parameters.f90`

Uses:
- `simple_atoms`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_decay_funs`
- `simple_parameters_registry`
- `simple_ui`
- `simple_ui_program`

Public symbols:
- `parameters` — type

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
## Module: simple_particle_extractor

Files:
- `main/simple_particle_extractor.f90`

Uses:
- `simple_core_module_api`
- `simple_eer_factory`
- `simple_image`
- `simple_motion_correct_utils`
- `simple_ori`
- `simple_starfile_wrappers`

Public symbols:
- `ptcl_extractor` — type

Private symbols:
- `cure_outliers` — subroutine
- `display` — subroutine
- `extract_particles` — subroutine
- `extract_particles_from_mic` — subroutine
- `extract_ptcl` — subroutine
- `get_local_shift` — subroutine
- `init_mask` — subroutine
- `init_mic` — subroutine
- `init_mov` — subroutine
- `kill` — subroutine
- `parse_movie_metadata` — subroutine
- `parse_string` — subroutine
- `pix2polycoords` — subroutine
- `post_process` — subroutine

---
## Module: simple_pca

Files:
- `main/pca/simple_pca.f90`

Public symbols:
- `pca` — type

Private symbols:
- `generic_generate` — subroutine
- `generic_get_feat` — function
- `generic_kill` — subroutine
- `generic_master` — subroutine
- `generic_new` — subroutine

---
## Module: simple_pca_svd

Files:
- `main/pca/simple_pca_svd.f90`

Uses:
- `simple_core_module_api`
- `simple_pca`

Public symbols:
- `pca_svd` — type

Private symbols:
- `generate_svd` — subroutine
- `get_feat_svd` — function
- `kill_svd` — subroutine
- `master_ori` — subroutine
- `master_svd` — subroutine
- `master_T` — subroutine
- `new_svd` — subroutine

---
## Module: simple_pcg_reconstruction

Files:
- `main/volume/simple_pcg_reconstruction.f90`

Uses:
- `simple_core_module_api`
- `simple_ctf`
- `simple_image`

Public symbols:
- `pcg_reconstruction` — type

Private symbols:
- `absT2_plane` — subroutine
- `adjoint_plane_add` — subroutine
- `apply_adjoint_all` — function
- `apply_normal` — function
- `apply_normal_kernel` — function
- `apply_normal_matrixfree` — function
- `apply_precond` — function
- `build_env` — subroutine
- `build_kernel` — subroutine
- `build_precond` — subroutine
- `build_transfer` — function
- `calibrate_kernel` — subroutine
- `crop_vol` — function
- `deapod_mul` — subroutine
- `dot_real_volume` — function
- `ensure_wimg` — subroutine
- `extract_native_plane` — function
- `fold_and_ifft` — function
- `forward_plane` — subroutine
- `fourier_dot` — function
- `get_env` — function
- `get_invenv` — function
- `get_lims2` — function
- `get_lims3` — function
- `kill` — subroutine
- `new` — subroutine
- `pad_vol` — function
- `prep_particles` — subroutine
- `scatter_plane` — subroutine
- `scatter_window` — subroutine
- `scatter_window_lims` — subroutine
- `set_deapod` — subroutine
- `set_op_mode` — subroutine
- `set_volume` — subroutine
- `solve` — subroutine
- `transfer_plane_cmplx` — subroutine

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

Private symbols:
- `claim_warmup_worker_ids` — subroutine
- `clear_registry_entry_by_fd` — subroutine
- `clear_stale_registry_entry` — subroutine
- `get_host_ips` — function
- `get_port` — function
- `get_queue_pressure_workers_required` — function
- `handle_heartbeat_msg` — subroutine
- `handle_new_task_msg` — subroutine
- `handle_terminate_msg` — subroutine
- `handle_unknown_msg` — subroutine
- `is_running` — function
- `kill` — subroutine
- `mark_worker_slots_launch_pending` — subroutine
- `new` — subroutine
- `pack_task_queue` — subroutine
- `queue_task` — function
- `send_idle_reply` — subroutine
- `set_warmup_cooldown_enabled` — subroutine
- `set_warmup_request` — subroutine
- `signal_ready` — subroutine
- `update_queue_pressure` — subroutine
- `update_worker_registry` — subroutine
- `worker_listener_thread` — subroutine

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
## Module: simple_pftc_shsrch_grad

Files:
- `main/pftc/simple_pftc_shsrch_grad.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_opt_factory`
- `simple_opt_spec`
- `simple_optimizer`

Public symbols:
- `pftc_shsrch_grad` — type

Private symbols:
- `coarse_search` — subroutine
- `coarse_search_opt_angle` — subroutine
- `grad_shsrch_costfun` — function
- `grad_shsrch_fdfcostfun` — subroutine
- `grad_shsrch_gcostfun` — subroutine
- `grad_shsrch_kill` — subroutine
- `grad_shsrch_minimize` — function
- `grad_shsrch_new` — subroutine
- `grad_shsrch_optimize_angle` — subroutine
- `grad_shsrch_optimize_angle_wrapper` — subroutine
- `grad_shsrch_set_indices` — subroutine
- `set_limits` — subroutine

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
- `simple_sigma2_binfile`
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
- `simple_pickseg`
- `simple_picksegdiam`
- `simple_string_utils`

Public symbols:
- `exec_gaupick` — subroutine
- `exec_refpick` — subroutine
- `exec_segdiampick` — subroutine
- `exec_segpick` — subroutine

---
## Module: simple_pickref

Files:
- `main/pick/simple_pickref.f90`

Uses:
- `simple_core_module_api`
- `simple_gui_utils`
- `simple_image`
- `simple_micproc`
- `simple_segmentation`

Public symbols:
- `multiref_merge` — subroutine
- `read_mic_raw_pickref` — subroutine

Private symbols:
- `detect_peaks` — subroutine
- `distance_filter` — subroutine
- `get_loc_sdevs` — subroutine
- `get_nboxes` — function
- `get_positions` — subroutine
- `get_scores` — subroutine
- `kill` — subroutine
- `match_boximgs` — subroutine
- `new` — subroutine
- `peak_vs_nonpeak_stats` — subroutine
- `refine_upscaled` — subroutine
- `refpick` — subroutine
- `remove_outliers` — subroutine
- `report_boxfile` — subroutine
- `report_thumb_den` — subroutine
- `setup_iterators` — subroutine
- `write_boxfile` — subroutine

---
## Module: simple_pickseg

Files:
- `main/pick/simple_pickseg.f90`

Uses:
- `simple_bspline_smoother`
- `simple_core_module_api`
- `simple_image`
- `simple_image_bin`
- `simple_parameters`
- `simple_segmentation`

Private symbols:
- `center_mass_cc` — function
- `get_boxsize` — function
- `get_ldim` — function
- `get_nboxes` — function
- `get_positions` — subroutine
- `get_smpd_shrink` — function
- `pick` — subroutine
- `report_boxfile` — subroutine

---
## Module: simple_picksegdiam

Files:
- `main/pick/simple_picksegdiam.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_image_bin`
- `simple_linked_list`
- `simple_micproc`
- `simple_nrtxtfile`

Private symbols:
- `dealloc_objs` — subroutine
- `dealloc_objs` — subroutine
- `get_diameters` — subroutine
- `kill` — subroutine
- `pick_1` — subroutine
- `pick_2` — subroutine
- `pick_3` — subroutine
- `write_diameters` — subroutine
- `write_pos_and_diams` — subroutine

---
## Module: simple_polarft_calc

Files:
- `main/pftc/simple_polarft_calc.f90`

Uses:
- `simple_pftc_api`

Public symbols:
- `polarft_calc` — type

---
## Module: simple_ppca

Files:
- `main/pca/simple_ppca.f90`

Uses:
- `simple_core_module_api`
- `simple_math`
- `simple_pca`
- `simple_rnd`

Public symbols:
- `ppca` — type

Private symbols:
- `calc_eigvals` — subroutine
- `em_opt` — subroutine
- `generate_ppca` — subroutine
- `get_eigvals_ppca` — function
- `get_feat_ppca` — function
- `get_signal_eigvals_ppca` — function
- `init` — subroutine
- `kill_ppca` — subroutine
- `master_ppca` — subroutine
- `new_ppca` — subroutine
- `reconstruct_external_ppca` — subroutine
- `set_verbose_ppca` — subroutine
- `slim_ppca` — subroutine

---
## Module: simple_preprocess_strategy

Files:
- `main/strategies/parallelization/simple_preprocess_strategy.f90`

Uses:
- `fox_dom`
- `simple_binoris_io`
- `simple_cmdline`
- `simple_commanders_api`
- `simple_ctf_estimate_iter`
- `simple_mini_stream_utils`
- `simple_motion_correct_iter`
- `simple_motion_correct_utils`
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
- `simple_commanders_cluster2d`
- `simple_commanders_denoise`
- `simple_commanders_distr`
- `simple_commanders_euclid`
- `simple_commanders_euclid_distr`
- `simple_commanders_flex_analysis`
- `simple_commanders_imgops`
- `simple_commanders_mask`
- `simple_commanders_misc`
- `simple_commanders_ori`
- `simple_commanders_pick`
- `simple_commanders_preprocess`
- `simple_commanders_prob`
- `simple_commanders_project_cls`
- `simple_commanders_project_core`
- `simple_commanders_project_ptcl`
- `simple_commanders_rec`
- `simple_commanders_rec_distr`
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
- `utils/simple_private_prgs.f90`

Uses:
- `simple_core_module_api`

Public symbols:
- `get_n_private_keys_required` — function
- `get_private_keys_required` — function
- `make_private_ui` — subroutine
- `print_cmdline_oldschool` — subroutine
- `print_private_cmdline` — subroutine

Private symbols:
- `get_n_private_keys_optional` — function
- `get_private_keys_optional` — function
- `init_cmd_dict` — subroutine
- `new_private_prgs` — subroutine
- `push_opt_key` — subroutine
- `push_req_key` — subroutine
- `set_name` — subroutine

---
## Module: simple_procimgstk

Files:
- `utils/simple_procimgstk.f90`

Uses:
- `simple_bspline_smoother`
- `simple_core_module_api`
- `simple_ctf`
- `simple_denoise_movies`
- `simple_image`
- `simple_sp_project`
- `simple_stack_io`

Public symbols:
- `add_noise_imgfile` — subroutine
- `apply_ctf_imgfile` — subroutine
- `bp_imgfile` — subroutine
- `bs_smoother_imgfile` — subroutine
- `clip_imgfile` — subroutine
- `copy_imgfile` — subroutine
- `diffmap_denoise_imgfile` — subroutine
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
- `fill_class_row` — subroutine
- `fill_particle_row` — subroutine
- `make_project` — subroutine
- `test_merge_pruned_stack_indexing_heterogeneous_ctf` — subroutine
- `test_validate_projfile_normalizes_legacy_stack_indices` — subroutine

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
- `kill_expanded` — subroutine
- `reset_expanded` — subroutine

---
## Module: simple_projector_pft

Files:
- `main/image/simple_projector_pft.f90`

Uses:
- `simple_core_module_api`
- `simple_polarft_calc`
- `simple_projector`

Public symbols:
- `fproject_polar` — subroutine
- `fproject_polar_oversamp` — subroutine

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

Private symbols:
- `fproject_polar_batch_opt_kernel` — subroutine

---
## Module: simple_projfile_utils

Files:
- `fileio/simple_projfile_utils.f90`

Uses:
- `simple_class_frcs`
- `simple_core_module_api`
- `simple_euclid_sigma2`
- `simple_image`
- `simple_oris`
- `simple_sp_project`

Public symbols:
- `absolutize_project_stack_paths` — subroutine
- `copy_particle_row` — subroutine
- `err` — subroutine
- `make_prefix_offsets` — subroutine
- `merge_chunk_projfiles` — subroutine
- `merge_selected_project_files` — subroutine
- `remap_row_ogid` — subroutine
- `repair` — subroutine
- `repair_particle_segment` — subroutine
- `repair_stack_counts_without_particles` — subroutine
- `repair_stack_nptcls_stk` — subroutine
- `repair_stack_ranges` — subroutine
- `require_row_field` — subroutine
- `set_stack_int_if_changed` — subroutine
- `validate_and_repair_project_file` — subroutine
- `validate_field_presence` — subroutine
- `validate_mic_sampling` — subroutine
- `validate_particle_field` — subroutine
- `validate_source_project` — subroutine
- `validate_stack_dimensions` — subroutine
- `validate_stack_particle_ranges` — subroutine
- `warn` — subroutine

---
## Module: simple_pspec_thumb_iter

Files:
- `main/simple_pspec_thumb_iter.f90`

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
- `main/simple_pspecs.f90`

Uses:
- `simple_core_module_api`
- `simple_fsc`
- `simple_image`
- `simple_oris`

Private symbols:
- `calc_distmat_1` — function
- `calc_distmat_2` — function
- `get_nspecs` — function
- `kill` — subroutine
- `new` — subroutine

---
## Module: simple_ptcl_sieve

Files:
- `main/sieve/simple_ptcl_sieve.f90`

Uses:
- `simple_cavg_quality_analysis`
- `simple_cavg_quality_feats`
- `simple_cavg_quality_helpers`
- `simple_cavg_quality_model`
- `simple_cavg_quality_types`
- `simple_class_compatibility`
- `simple_cmdline`
- `simple_commanders_cavgs`
- `simple_defs`
- `simple_defs_fname`
- `simple_error`
- `simple_fileio`
- `simple_gui_utils`
- `simple_image`
- `simple_image_bin`
- `simple_imgarr_utils`
- `simple_parameters`
- `simple_projfile_utils`
- `simple_qsys_env`
- `simple_rec_list`
- `simple_segmentation`
- `simple_sp_project`
- `simple_string`
- `simple_string_utils`
- `simple_syslib`
- `simple_timer`
- `unix`

Public symbols:
- `ptcl_sieve` — type

Private symbols:
- `append_chunk_coarse` — subroutine
- `append_chunk_fine` — subroutine
- `chunk2D_coarse_defaults` — type
- `chunk2D_fine_defaults` — type
- `chunk2D_state` — type
- `cleanup_chunk` — subroutine
- `collect_and_reject` — subroutine
- `combine_completed_chunks` — subroutine
- `cycle` — subroutine
- `draw_reason_key_swatches` — subroutine
- `generate_chunk_coarse_cline` — subroutine
- `generate_chunk_fine_cline` — subroutine
- `generate_chunks_coarse` — subroutine
- `generate_chunks_fine` — subroutine
- `import_existing_chunks_coarse` — subroutine
- `import_existing_chunks_fine` — subroutine
- `kill` — subroutine
- `merge_and_clear` — subroutine
- `new` — subroutine
- `paint_mask_outline` — subroutine
- `paint_reason_border` — subroutine
- `reject_cavgs` — subroutine
- `remove_duplicates` — function
- `set_final_ingestion` — subroutine
- `submit` — subroutine
- `timer_stop` — subroutine
- `unset_final_ingestion` — subroutine
- `write_rejection_reason_overlay_jpg` — subroutine

---
## Module: simple_ptcl_sieve_tester

Files:
- `main/sieve/simple_ptcl_sieve_tester.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_fname`
- `simple_parameters`
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
- `setup_workspace` — subroutine
- `teardown_workspace` — subroutine
- `test_cycle_empty_project_list` — subroutine
- `test_finished_semantics` — subroutine
- `test_import_existing_chunks_and_counts` — subroutine
- `test_new_accepts_tuning_overrides` — subroutine
- `test_new_kill_and_empty_queries` — subroutine
- `test_single_pass_ignores_incomplete_fine` — subroutine

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
- `simple_qsys_lsf`
- `simple_qsys_persistent_worker`
- `simple_qsys_slurm`
- `simple_syslib`

Private symbols:
- `worker_warmup_callback` — subroutine

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
- `qsys_env` — type

Private symbols:
- `diagnose_parallelism` — subroutine
- `exec_simple_prg_in_queue` — subroutine
- `exec_simple_prg_in_queue_async` — subroutine
- `exec_simple_prgs_in_queue_async` — subroutine
- `exists` — function
- `gen_script` — subroutine
- `gen_scripts_and_schedule_jobs` — subroutine
- `gen_subproject_scripts_and_schedule` — subroutine
- `get_exec_bin` — function
- `get_persistent_worker_server_address` — function
- `get_qsys` — function
- `init_qdescr_from_runtime` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `persistent_worker_warmup_callback` — subroutine
- `schedule_subproject_jobs` — subroutine
- `service_persistent_worker_warmup` — subroutine
- `standard_exec_path` — subroutine
- `start_persistent_workers` — subroutine
- `validate_requested_qsys` — subroutine

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
- `write_formatted` — subroutine
- `write_pbs_array_header` — subroutine
- `write_pbs_header` — subroutine

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
## Module: simple_r8lib

Files:
- `utils/math/simple_r8lib.f90`

Public symbols:
- `r8mat_cholesky_factor` — subroutine
- `r8mat_cholesky_solve` — subroutine

Private symbols:
- `r8mat_l_solve` — subroutine
- `r8mat_lt_solve` — subroutine

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

---
## Module: simple_rec3D_strategy

Files:
- `main/strategies/parallelization/simple_rec3D_strategy.f90`

Uses:
- `simple_builder`
- `simple_cmdline`
- `simple_commanders_rec_distr`
- `simple_commanders_volops`
- `simple_core_module_api`
- `simple_exec_helpers`
- `simple_matcher_3drec`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_refine3d_fnames`

Public symbols:
- `create_rec3D_strategy` — function
- `rec3D_distr_strategy` — type
- `rec3D_inmem_strategy` — type
- `rec3D_strategy` — type

Private symbols:
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_execute` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `exec_interface` — subroutine
- `finalize_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_execute` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `maybe_postprocess_reconstruct3D` — subroutine
- `sync_resolved_rec_params` — subroutine

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
- `get_processed_flags` — function
- `get_projfiles` — function
- `get_projnames` — function
- `iter_advance` — subroutine
- `iter_get_chunk_rec` — subroutine
- `iter_get_process_rec` — subroutine
- `iter_get_project_rec` — subroutine
- `iter_get_rec` — subroutine
- `iter_next` — subroutine
- `iter_replace` — subroutine
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
- `simple_fsc`
- `simple_image`
- `simple_kbinterpol`
- `simple_math`
- `simple_parameters`
- `simple_sp_project`

Public symbols:
- `reconstructor` — type

Private symbols:
- `add_conical_invtausq2rho` — subroutine
- `add_invtausq2rho` — subroutine
- `alloc_rho` — subroutine
- `apply_weight` — subroutine
- `compress_exp` — subroutine
- `dealloc_exp` — subroutine
- `dealloc_rho` — subroutine
- `expand_exp` — subroutine
- `get_kbwin` — function
- `insert_plane_oversamp` — subroutine
- `insert_plane_oversamp_opt` — subroutine
- `interp_cmat_exp` — function
- `kb_apod_vecs_3d_fast` — subroutine
- `kernel` — subroutine
- `pad_with_zeros` — subroutine
- `project_fplane` — subroutine
- `project_polar` — subroutine
- `read_raw_rho` — subroutine
- `read_rho` — subroutine
- `reset` — subroutine
- `reset_exp` — subroutine
- `sampl_dens_correct` — subroutine
- `set_sh_lim` — subroutine
- `sum_reduce` — subroutine
- `write_absfc_as_mrc` — subroutine
- `write_rho` — subroutine
- `write_rho_as_mrc` — subroutine

---
## Module: simple_reconstructor_eo

Files:
- `main/volume/simple_reconstructor_eo.f90`

Uses:
- `simple_core_module_api`
- `simple_fsc`
- `simple_image`
- `simple_image_msk`
- `simple_imgfile`
- `simple_parameters`
- `simple_reconstructor`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_vol_pproc_policy`

Public symbols:
- `reconstructor_eo` — type

Private symbols:
- `apply_weight` — subroutine
- `calc_fsc4sampl_dens_correct` — subroutine
- `compress_exp` — subroutine
- `expand_exp` — subroutine
- `get_kbwin` — function
- `get_res` — subroutine
- `grid_plane` — subroutine
- `grid_plane_compact` — subroutine
- `kill` — subroutine
- `kill_exp` — subroutine
- `load_state_mask_or_fallback` — subroutine
- `new` — subroutine
- `project_polar` — subroutine
- `read_eos` — subroutine
- `read_eos_parallel_io` — subroutine
- `read_even` — subroutine
- `read_odd` — subroutine
- `reset_all` — subroutine
- `reset_eoexp` — subroutine
- `reset_eos` — subroutine
- `reset_even` — subroutine
- `reset_odd` — subroutine
- `reset_sum` — subroutine
- `sampl_dens_correct_eos` — subroutine
- `sampl_dens_correct_sum` — subroutine
- `set_sh_lim` — subroutine
- `sum_eos` — subroutine
- `sum_reduce` — subroutine
- `write_eos` — subroutine
- `write_even` — subroutine
- `write_fsc2txt` — subroutine
- `write_odd` — subroutine

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

Public symbols:
- `calc_3Drec_gpu` — subroutine

Private symbols:
- `insert_all_slices` — subroutine
- `insert_slices_batch` — subroutine
- `prep_batch` — subroutine
- `update_ptcls_arrays` — subroutine

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
- `killimgbatch_local` — subroutine
- `prepimgbatch_local` — subroutine
- `set_reextract_defaults` — subroutine
- `validate_reextract_cline` — subroutine

---
## Module: simple_refine3D_fnames

Files:
- `defs/simple_refine3D_fnames.f90`

Uses:
- `simple_defs_fname`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `refine3D_partial_rec_glob` — function

---
## Module: simple_refine3D_strategy

Files:
- `main/strategies/parallelization/simple_refine3D_strategy.f90`

Uses:
- `simple_builder`
- `simple_cluster_seed`
- `simple_cmdline`
- `simple_commanders_euclid`
- `simple_commanders_prob`
- `simple_commanders_rec`
- `simple_commanders_rec_distr`
- `simple_commanders_volops`
- `simple_convergence`
- `simple_core_module_api`
- `simple_decay_funs`
- `simple_euclid_sigma2`
- `simple_exec_helpers`
- `simple_fsc`
- `simple_image`
- `simple_image_msk`
- `simple_matcher_refvol_utils`
- `simple_matcher_smpl_and_lplims`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_refine3d_fnames`
- `simple_sp_project`
- `simple_strategy3d_matcher`

Public symbols:
- `create_refine3D_strategy` — function
- `refine3D_distr_strategy` — type
- `refine3D_inmem_strategy` — type
- `refine3D_strategy` — type

Private symbols:
- `activate_ptcl3D_states_from_selection` — subroutine
- `assert_multistate_populations` — subroutine
- `cleanup_interface` — subroutine
- `distr_cleanup` — subroutine
- `distr_execute_iteration` — subroutine
- `distr_finalize_iteration` — subroutine
- `distr_finalize_run` — subroutine
- `distr_initialize` — subroutine
- `exec_iter_interface` — subroutine
- `finalize_iter_interface` — subroutine
- `finalize_run_interface` — subroutine
- `init_interface` — subroutine
- `inmem_cleanup` — subroutine
- `inmem_execute_iteration` — subroutine
- `inmem_finalize_iteration` — subroutine
- `inmem_finalize_run` — subroutine
- `inmem_initialize` — subroutine
- `invalidate_fresh_start_refs_from_volumes` — subroutine
- `materialize_reprojection_model` — subroutine
- `prepare_assembly_cline` — subroutine
- `refine3D_bench_state` — type
- `refine3D_stage_bench_state` — type
- `refresh_matching_lp_from_project` — subroutine
- `refresh_resolution_fields_from_fsc` — subroutine
- `remove_partial_rec_files` — subroutine
- `reset_refine3D_bench` — subroutine
- `seed_multistate_startup_labels` — subroutine
- `write_stage_bench_report` — subroutine
- `write_strategy_bench_report` — subroutine

---
## Module: simple_relion

Files:
- `main/star/simple_relion.f90`

Uses:
- `cplot2d_wrapper_module`
- `fox_dom`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_sp_project`
- `simple_starfile_wrappers`

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

---
## Module: simple_rnd

Files:
- `utils/math/simple_rnd.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_syslib`

---
## Module: simple_sauron

Files:
- `utils/simple_sauron.f90`

Uses:
- `simple_chash`
- `simple_defs`
- `simple_defs_ori`
- `simple_hash`
- `simple_string_utils`

---
## Module: simple_segmentation

Files:
- `utils/simple_segmentation.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_image_bin`
- `simple_neighs`

---
## Module: simple_sigma2_binfile

Files:
- `fileio/simple_sigma2_binfile.f90`

Uses:
- `simple_core_module_api`

Private symbols:
- `create_empty` — subroutine
- `get_resrange` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `new_from_file` — subroutine
- `open_and_check_header` — function
- `read` — subroutine
- `read_header` — subroutine
- `write` — subroutine
- `write_info` — subroutine

---
## Module: simple_sigma2_files

Files:
- `fileio/simple_sigma2_files.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_euclid_sigma2`
- `simple_parameters`
- `simple_polarft_calc`
- `simple_sp_project`

Public symbols:
- `carry_over_prior_sigma_files` — subroutine
- `load_sigma2_groups` — subroutine
- `pick_sigma_group_file` — subroutine

Private symbols:
- `copy_sigma_files_from_dir` — subroutine

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

Public symbols:
- `simimg` — subroutine

---
## Module: simple_socket_comm

Files:
- `utils/comm/simple_socket_comm.f90`

Public symbols:
- `simple_socket` — type

Private symbols:
- `accept` — subroutine
- `bind_any` — subroutine
- `c_accept` — function
- `c_bind` — function
- `c_close` — function
- `c_connect` — function
- `c_errno` — function
- `c_htons` — function
- `c_inet_pton` — function
- `c_listen` — function
- `c_read` — function
- `c_send` — function
- `c_setsockopt` — function
- `c_socket` — function
- `close_1` — subroutine
- `close_2` — subroutine
- `in_addr` — type
- `listen` — subroutine
- `open` — subroutine
- `read_1` — subroutine
- `read_2` — subroutine
- `send_1` — subroutine
- `set_options` — subroutine
- `sockaddr_in` — type
- `timeval` — type

---
## Module: simple_sp_project

Files:
- `main/project/simple_sp_project.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_ansi_ctrls`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_ori`
- `simple_discrete_stack_io`
- `simple_gui_utils`
- `simple_histogram`
- `simple_image`
- `simple_map_reduce`
- `simple_rec_list`
- `simple_stack_io`
- `simple_starfile`

---
## Module: simple_srch_sort_loc

Files:
- `utils/math/simple_srch_sort_loc.f90`

Uses:
- `simple_defs`
- `simple_error`

---
## Module: simple_srchspace_map

Files:
- `utils/simple_srchspace_map.f90`

Uses:
- `simple_error`
- `simple_srch_sort_loc`

Public symbols:
- `srchspace_map` — type

Private symbols:
- `get_full2sub_map` — function
- `get_inds_in_full` — function
- `get_sub2full_map` — function
- `kill` — subroutine
- `new_from_distmat` — subroutine
- `new_from_maps` — subroutine

---
## Module: simple_srchspace_map2D_io

Files:
- `fileio/simple_srchspace_map2D_io.f90`

Uses:
- `simple_error`
- `simple_fileio`
- `simple_srch_sort_loc`
- `simple_string`

Public symbols:
- `read_srchspace_map2D` — subroutine
- `test_srchspace_map2D_io` — subroutine
- `write_srchspace_map2D` — subroutine

---
## Module: simple_stack_io

Files:
- `fileio/simple_stack_io.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_imgfile`

Private symbols:
- `close` — subroutine
- `get_image` — subroutine
- `get_ldim` — function
- `get_nptcls` — function
- `open_1` — subroutine
- `open_2` — subroutine
- `read` — subroutine
- `read_whole` — subroutine
- `same_stk` — function
- `write` — subroutine
- `write_buffer` — subroutine

---
## Module: simple_stackops

Files:
- `utils/simple_stackops.f90`

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

Private symbols:
- `calc_part_boundaries` — subroutine
- `complete` — subroutine
- `get_relative_path_here` — function
- `get_relative_path_here` — function
- `get_relative_path_here` — function
- `init` — subroutine
- `write_mics_table` — subroutine
- `write_optics_table` — subroutine
- `write_ptcl2D_table` — subroutine
- `write_ptcl2D_table_parallel` — subroutine

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
- `fill_ptcl_stk` — subroutine
- `make_oris_n` — subroutine
- `setup_tmpdir` — subroutine
- `test_complete_noop_without_tmp` — subroutine
- `test_complete_overwrites_existing` — subroutine
- `test_complete_renames_tmp` — subroutine
- `test_init_creates_no_files` — subroutine
- `test_init_removes_stale_tmp` — subroutine
- `test_verbose_flag_no_content_change` — subroutine
- `test_write_mics_table_all_fields` — subroutine
- `test_write_mics_table_basic` — subroutine
- `test_write_mics_table_empty` — subroutine
- `test_write_mics_table_state_filter` — subroutine
- `test_write_optics_table_all_fields` — subroutine
- `test_write_optics_table_basic` — subroutine
- `test_write_optics_table_empty` — subroutine
- `test_write_optics_table_state_filter` — subroutine
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
- `starfile_table_type` — type

Private symbols:
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

---
## Module: simple_starproject

Files:
- `main/star/simple_starproject.f90`

Uses:
- `cplot2d_wrapper_module`
- `fox_dom`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_refine3d_fnames`
- `simple_rnd`
- `simple_sp_project`
- `simple_starproject_utils`

Private symbols:
- `apply_optics_offset` — subroutine
- `assign_initial_tiltgroups` — subroutine
- `assign_optics` — subroutine
- `assign_xml_tiltinfo` — subroutine
- `check_stk_params` — subroutine
- `cluster_tiltinfo` — subroutine
- `export_cls2D` — subroutine
- `export_iter3D` — subroutine
- `export_manifoldem_ptcls3D` — subroutine
- `export_mics` — subroutine
- `export_opticsgroups` — subroutine
- `export_ptcls2D` — subroutine
- `export_ptcls3D` — subroutine
- `export_stardata` — subroutine
- `export_stream2D` — subroutine
- `get_image_basename` — subroutine
- `import_cls2D` — subroutine
- `import_mics` — subroutine
- `import_ptcls` — subroutine
- `import_ptcls2D` — subroutine
- `import_ptcls3D` — subroutine
- `import_stardata` — subroutine
- `initialise` — subroutine
- `kill` — subroutine
- `plot_opticsgroups` — subroutine
- `populate_opticsmap` — subroutine
- `populate_stkmap` — subroutine
- `propagate_optics2D` — subroutine
- `propagate_optics3D` — subroutine
- `propagate_optics_box` — subroutine
- `read_starheaders` — subroutine
- `set_verbose` — subroutine
- `sort_optics_maxpop` — subroutine

---
## Module: simple_starproject_stream

Files:
- `main/star/simple_starproject_stream.f90`

Uses:
- `cplot2d_wrapper_module`
- `fox_dom`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_histogram`
- `simple_parameters`
- `simple_sp_project`
- `simple_starfile_wrappers`

Private symbols:
- `assign_optics` — subroutine
- `assign_optics_single` — subroutine
- `assign_shiftgroups` — subroutine
- `assign_tiltgroups` — subroutine
- `copy_micrographs_optics` — subroutine
- `copy_optics` — subroutine
- `get_relative_path_here` — function
- `get_relative_path_here` — function
- `get_relative_path_here` — function
- `get_relative_path_here` — function
- `h_clust` — subroutine
- `starfile_deinit` — subroutine
- `starfile_init` — subroutine
- `starfile_set_clusters2D_table` — subroutine
- `starfile_set_micrographs_table` — subroutine
- `starfile_set_optics_group_table` — subroutine
- `starfile_set_optics_table` — subroutine
- `starfile_set_particles2D_subtable` — subroutine
- `starfile_set_particles2D_table` — subroutine
- `starfile_set_pick_diameters_table` — subroutine
- `starfile_write_table` — subroutine
- `starpart` — type
- `stream_export_micrographs` — subroutine
- `stream_export_optics` — subroutine
- `stream_export_particles_2D` — subroutine
- `stream_export_pick_diameters` — subroutine
- `stream_export_picking_references` — subroutine

---
## Module: simple_starproject_tester

Files:
- `main/star/simple_starproject_tester.f90`

Uses:
- `omp_lib`
- `simple_cmdline`
- `simple_defs`
- `simple_parameters`
- `simple_relion`
- `simple_sp_project`
- `simple_starfile`
- `simple_starfile_wrappers`
- `simple_starproject`
- `simple_starproject_stream`
- `simple_starproject_utils`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_starproject_tests` — subroutine

Private symbols:
- `create_directory` — subroutine
- `delete_directory` — subroutine
- `force_openmp_threads` — subroutine
- `make_tiltinfo` — function
- `setup_tmpdir` — subroutine
- `test_optics_clustering_basic` — subroutine
- `test_parallel_export_readback` — subroutine
- `test_parallel_particle_export` — subroutine
- `test_relion_allocate_opticsgroups` — subroutine
- `test_relion_epu_tiltgroups` — subroutine
- `test_relion_find_movienames` — subroutine
- `test_relion_single_tiltgroup` — subroutine
- `test_relion_write_corrected_micrographs_star` — subroutine
- `test_relion_write_micrographs_star` — subroutine
- `test_relion_write_particles2D_star` — subroutine
- `test_relion_writer_micrographs` — subroutine
- `test_roundtrip_micrographs` — subroutine
- `test_star_export_micrographs` — subroutine
- `test_star_import_mics` — subroutine
- `test_starfile_append_mode` — subroutine
- `test_starfile_basic_init` — subroutine
- `test_starfile_multiple_tables` — subroutine
- `test_starfile_write_and_readback` — subroutine
- `test_write_omem_parallel_stability` — subroutine
- `test_xml_tiltinfo` — subroutine
- `tmpfile` — function
- `write_textfile` — subroutine

---
## Module: simple_starproject_utils

Files:
- `main/star/simple_starproject_utils.f90`

Uses:
- `simple_core_module_api`
- `simple_sp_project`

Public symbols:
- `center_boxes` — subroutine
- `enable_rlnflag` — subroutine
- `enable_splflag` — subroutine
- `enable_splflags` — subroutine
- `find_separators` — subroutine
- `get_rlnflagindex` — subroutine
- `h_clust` — subroutine
- `split_dataline` — subroutine

---
## Module: simple_stat

Files:
- `utils/math/simple_stat.f90`

Uses:
- `simple_defs`
- `simple_error`
- `simple_is_check_assert`
- `simple_srch_sort_loc`
- `simple_type_defs`

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
- `prep_strategy2D_batch` — subroutine
- `prep_strategy2D_glob` — subroutine
- `prep_strategy2D_thread` — subroutine
- `set_strategy2D_stoch_bound` — subroutine

---
## Module: simple_strategy2D_eval

Files:
- `main/strategies/search/simple_strategy2D_eval.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_oris`
- `simple_parameters`
- `simple_strategy2d`
- `simple_strategy2d_srch`

Public symbols:
- `strategy2D_eval` — type

Private symbols:
- `kill_eval` — subroutine
- `new_eval` — subroutine
- `srch_eval` — subroutine

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
- `simple_binoris_io`
- `simple_builder`
- `simple_classaverager`
- `simple_convergence`
- `simple_decay_funs`
- `simple_eul_prob_tab2d`
- `simple_imgarr_utils`
- `simple_matcher_pftc_prep`
- `simple_matcher_ptcl_batch`
- `simple_matcher_smpl_and_lplims`
- `simple_pftc_srch_api`
- `simple_progress`
- `simple_qsys_funs`
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

Public symbols:
- `cluster2D_exec` — subroutine
- `set_b_p_ptrs2D` — subroutine

Private symbols:
- `allocate_strategy_for_particle` — subroutine
- `build_batch_particles_local` — subroutine
- `cleanup_search_state` — subroutine
- `cluster2D_ctrl` — type
- `compute_neigh_frac` — subroutine
- `display` — subroutine
- `ensure_even_odd_partition` — subroutine
- `finalize_restoration_and_convergence` — subroutine
- `init_ctrl` — subroutine
- `maybe_write_bench` — subroutine
- `prepare_alignment_references` — subroutine
- `prepare_batches` — subroutine
- `prepare_class_averages_and_restoration` — subroutine
- `restore_class_averages_for_batch` — subroutine
- `sample_particles_for_update` — subroutine
- `write_orientations` — subroutine

---
## Module: simple_strategy2D_prob

Files:
- `main/strategies/search/simple_strategy2D_prob.f90`

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
- `simple_eul_prob_tab2d`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_strategy2d_alloc`

Private symbols:
- `assign_ori` — subroutine
- `inpl_srch` — subroutine
- `inpl_srch_first` — subroutine
- `inpl_srch_peaks` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `prep4srch` — subroutine
- `store_solution` — subroutine

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
- `per_ref_srch` — subroutine
- `srch_tseries` — subroutine

---
## Module: simple_strategy2D_utils

Files:
- `main/strategies/search/simple_strategy2D_utils.f90`

Uses:
- `simple_atoms`
- `simple_builder`
- `simple_clustering_utils`
- `simple_cmdline`
- `simple_histogram`
- `simple_image_bin`
- `simple_image_msk`
- `simple_molecule_data`
- `simple_parameters`
- `simple_pftc_shsrch_grad`
- `simple_pftc_srch_api`
- `simple_pspecs`
- `simple_segmentation`
- `simple_simple_volinterp`
- `simple_sp_project`

Public symbols:
- `align_and_score_cavg_clusters` — function
- `calc_cavg_offset` — subroutine
- `calc_cavg_pairwise_algninfo` — function
- `calc_cavg_sigstats_components` — subroutine
- `calc_cluster_cavgs_dmat` — function
- `calc_match_cavgs_dmat` — function
- `flag_non_junk_cavgs` — subroutine
- `id_junk` — subroutine
- `id_junk_and_prep_cavgs4clust` — subroutine
- `prep_cavgs4clust` — subroutine
- `test_strategy2D_utils` — subroutine
- `write_aligned_cavgs` — subroutine

Private symbols:
- `align_clusters2medoids` — function
- `calc_cc_and_res_dmats` — subroutine
- `calc_cc_and_res_dmats_ref` — subroutine
- `calc_scores` — subroutine
- `calc_sigstats_dmats` — subroutine
- `calc_sigstats_dmats_ref` — subroutine
- `cavg_sigstats_matrices` — type
- `identify_good_bad_clusters` — subroutine
- `kill_cavg_sigstats_matrices` — subroutine
- `match_imgs` — function
- `match_imgs2ref` — function
- `rank_clusters` — subroutine
- `rtsq_imgs` — subroutine

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

---
## Module: simple_strategy3D_eval

Files:
- `main/strategies/search/simple_strategy3D_eval.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d`
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
- `simple_strategy3d`
- `simple_strategy3d_alloc`
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
- `simple_strategy3d`
- `simple_strategy3d_alloc`
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
## Module: simple_strategy3D_greedy_smpl

Files:
- `main/strategies/search/simple_strategy3D_greedy_smpl.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d`
- `simple_strategy3d_alloc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`
- `simple_type_defs`

Public symbols:
- `strategy3D_greedy_smpl` — type

Private symbols:
- `kill_greedy_smpl` — subroutine
- `new_greedy_smpl` — subroutine
- `oris_assign_greedy_smpl` — subroutine
- `srch_greedy_smpl` — subroutine

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
- `simple_strategy3d`
- `simple_strategy3d_alloc`
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
## Module: simple_strategy3D_matcher

Files:
- `main/strategies/search/simple_strategy3D_matcher.f90`

Uses:
- `simple_binoris_io`
- `simple_builder`
- `simple_euclid_sigma2`
- `simple_eul_prob_tab`
- `simple_matcher_2dprep`
- `simple_matcher_3drec`
- `simple_matcher_ptcl_batch`
- `simple_matcher_refvol_utils`
- `simple_matcher_smpl_and_lplims`
- `simple_pftc_srch_api`
- `simple_qsys_funs`
- `simple_refine3d_fnames`
- `simple_strategy3d`
- `simple_strategy3d_alloc`
- `simple_strategy3d_eval`
- `simple_strategy3d_greedy`
- `simple_strategy3d_greedy_inpl`
- `simple_strategy3d_greedy_smpl`
- `simple_strategy3d_greedy_sub`
- `simple_strategy3d_prob`
- `simple_strategy3d_shc`
- `simple_strategy3d_shc_smpl`
- `simple_strategy3d_snhc_smpl`
- `simple_strategy3d_srch`
- `simple_syslib`

Public symbols:
- `refine3D_exec` — subroutine

Private symbols:
- `build_batch_particles_local` — subroutine
- `choose_and_run_strategy` — subroutine
- `ensure_even_odd_partition` — subroutine
- `init_ctrl` — subroutine
- `maybe_write_orientations` — subroutine
- `prepare_particles_batches` — subroutine
- `prepare_refs_sigmas_and_pftc` — subroutine
- `print_flags` — subroutine
- `refine3D_ctrl` — type
- `sample_particles_for_update` — subroutine
- `strategy3D_per_ptcl` — type

---
## Module: simple_strategy3D_prob

Files:
- `main/strategies/search/simple_strategy3D_prob.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d`
- `simple_strategy3d_alloc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`

Public symbols:
- `strategy3D_prob` — type

Private symbols:
- `assign_state_fixed` — subroutine
- `kill_prob` — subroutine
- `new_prob` — subroutine
- `oris_assign_prob` — subroutine
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
- `simple_strategy3d`
- `simple_strategy3d_alloc`
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
## Module: simple_strategy3D_shc_smpl

Files:
- `main/strategies/search/simple_strategy3D_shc_smpl.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_strategy3d`
- `simple_strategy3d_alloc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`
- `simple_type_defs`

Public symbols:
- `strategy3D_shc_smpl` — type

Private symbols:
- `kill_shc_smpl` — subroutine
- `new_shc_smpl` — subroutine
- `oris_assign_shc_smpl` — subroutine
- `srch_shc_smpl` — subroutine

---
## Module: simple_strategy3D_snhc_smpl

Files:
- `main/strategies/search/simple_strategy3D_snhc_smpl.f90`

Uses:
- `simple_builder`
- `simple_decay_funs`
- `simple_eul_prob_tab_utils`
- `simple_oris`
- `simple_parameters`
- `simple_pftc_srch_api`
- `simple_strategy3d`
- `simple_strategy3d_alloc`
- `simple_strategy3d_srch`
- `simple_strategy3d_utils`
- `simple_type_defs`

Public symbols:
- `strategy3D_snhc_smpl` — type

Private symbols:
- `kill_snhc_smpl` — subroutine
- `new_snhc_smpl` — subroutine
- `oris_assign_snhc_smpl` — subroutine
- `srch_snhc_smpl` — subroutine

---
## Module: simple_strategy3D_srch

Files:
- `main/strategies/search/simple_strategy3D_srch.f90`

Uses:
- `simple_builder`
- `simple_core_module_api`
- `simple_eul_prob_tab`
- `simple_parameters`
- `simple_pftc_shsrch_grad`
- `simple_strategy3d_alloc`

Private symbols:
- `inpl_srch` — subroutine
- `inpl_srch_first` — subroutine
- `inpl_srch_peaks` — subroutine
- `kill` — subroutine
- `new` — subroutine
- `prep4prob` — subroutine
- `prep4srch` — subroutine
- `store_solution` — subroutine

---
## Module: simple_strategy3D_utils

Files:
- `main/strategies/search/simple_strategy3D_utils.f90`

Uses:
- `simple_core_module_api`
- `simple_strategy3d_alloc`
- `simple_strategy3d_srch`

Public symbols:
- `assign_ori` — subroutine
- `extract_peak_ori` — subroutine
- `extract_peak_oris` — subroutine

---
## Module: simple_stream2D_state

Files:
- `main/stream/simple_stream2D_state.f90`

Uses:
- `json_module`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_qsys_env`
- `simple_sp_project`
- `simple_starproject`
- `simple_stream_chunk`

---
## Module: simple_stream_abinitio2D_chunks

Files:
- `main/stream/simple_stream_abinitio2D_chunks.f90`

Uses:
- `simple_stream_api`

Public symbols:
- `exec_stream_abinitio2D_chunks` — subroutine
- `generate_chunk_projects` — subroutine
- `init_one_chunk` — subroutine
- `stream_abinitio2D_chunks` — type
- `submit_one_chunk` — subroutine
- `wait_one_chunk` — subroutine

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
- `simple_euclid_sigma2`
- `simple_gui_utils`
- `simple_guistats`
- `simple_image`
- `simple_nice`
- `simple_parameters`
- `simple_progress`
- `simple_projfile_utils`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_rec_list`
- `simple_sp_project`
- `simple_stack_io`
- `simple_starproject_stream`
- `simple_stream2d_state`
- `simple_stream_chunk`
- `simple_stream_chunk2d_utils`
- `simple_stream_cluster2d_utils`
- `simple_stream_communicator`
- `simple_stream_utils`
- `simple_stream_watcher`

---
## Module: simple_stream_chunk

Files:
- `main/stream/simple_stream_chunk.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_environment`
- `simple_euclid_sigma2`
- `simple_image`
- `simple_parameters`
- `simple_qsys_env`
- `simple_qsys_funs`
- `simple_rec_list`
- `simple_sp_project`

Private symbols:
- `analyze2D` — subroutine
- `assign` — subroutine
- `average_into` — subroutine
- `calc_sigma2` — subroutine
- `copy` — subroutine
- `debug_print` — subroutine
- `display_iter` — subroutine
- `gen_final_cavgs` — subroutine
- `generate` — subroutine
- `get_projfile_fname` — function
- `init_chunk` — subroutine
- `is_available` — function
- `kill` — subroutine
- `print_info` — subroutine
- `read` — subroutine
- `remove_folder` — subroutine
- `split_sigmas_into` — subroutine
- `terminate_chunk` — subroutine
- `to_analyze2D` — function

---
## Module: simple_stream_chunk2D_utils

Files:
- `main/stream/simple_stream_chunk2D_utils.f90`

Uses:
- `simple_cmdline`
- `simple_core_module_api`
- `simple_defs_environment`
- `simple_gui_utils`
- `simple_parameters`
- `simple_rec_list`
- `simple_sp_project`
- `simple_stream2d_state`
- `simple_stream_chunk`
- `simple_stream_cluster2d_utils`

Public symbols:
- `analyze2D_new_chunks` — subroutine
- `init_chunk_clustering` — subroutine
- `memoize_chunks` — subroutine
- `update_chunks` — subroutine

Private symbols:
- `set_chunk_dimensions` — subroutine

---
## Module: simple_stream_cluster2D_utils

Files:
- `main/stream/simple_stream_cluster2D_utils.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_class_frcs`
- `simple_cmdline`
- `simple_commanders_cluster2d`
- `simple_core_module_api`
- `simple_euclid_sigma2`
- `simple_image`
- `simple_parameters`
- `simple_qsys_funs`
- `simple_rec_list`
- `simple_sp_project`
- `simple_stack_io`
- `simple_starproject`
- `simple_starproject_stream`
- `simple_stream2d_state`

Public symbols:
- `cleanup_root_folder` — subroutine
- `consolidate_sigmas` — subroutine
- `setup_downscaling` — subroutine
- `terminate_chunks` — subroutine
- `terminate_stream2D` — subroutine
- `tidy_2Dstream_iter` — subroutine
- `update_user_params2D` — subroutine
- `write_project_stream2D` — subroutine
- `write_repick_refs` — subroutine

Private symbols:
- `apply_snapshot_selection` — subroutine
- `debug_print` — subroutine
- `get_latest_optics_map` — subroutine
- `get_latest_optics_map` — subroutine
- `rank_cavgs` — subroutine
- `rescale_cavgs` — subroutine
- `rescale_refs` — subroutine
- `set_dimensions` — subroutine
- `set_resolution_limits` — subroutine
- `set_snapshot_time` — subroutine
- `write_raw_project` — subroutine

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
- `add_to_json_1` — subroutine
- `add_to_json_2` — subroutine
- `add_to_json_3` — subroutine
- `add_to_json_4` — subroutine
- `add_to_json_5` — subroutine
- `background_heartbeat` — subroutine
- `create` — subroutine
- `curl_request` — subroutine
- `destroy_arg` — subroutine
- `evaluate_checksum` — function
- `get_json_arg_1` — subroutine
- `get_json_arg_2` — subroutine
- `get_json_arg_3` — subroutine
- `get_json_arg_4` — subroutine
- `get_json_arg_5` — subroutine
- `join_background_heartbeat` — subroutine
- `remove_from_json_if_present` — subroutine
- `send_heartbeat` — subroutine
- `send_jobstats` — subroutine
- `term` — subroutine
- `update_json_1` — subroutine
- `update_json_2` — subroutine
- `update_json_3` — subroutine
- `update_json_4` — subroutine

---
## Module: simple_stream_p00_master

Files:
- `main/stream/simple_stream_p00_master.f90`

Uses:
- `simple_forked_process`
- `simple_gui_assembler`
- `simple_gui_metadata_api`
- `simple_gui_metadata_utils`
- `simple_http_post`
- `simple_memory_monitor`
- `simple_stream_api`
- `simple_stream_p01_preprocess_new`
- `simple_stream_p02_assign_optics_new`
- `simple_stream_p03_initial_analysis`
- `simple_stream_p04_refpick_extract_new`
- `simple_stream_p05_sieve_cavgs_new`
- `simple_stream_p06_pool2d_new`
- `simple_stream_state`
- `simple_syslib`
- `unix`

Public symbols:
- `stream_p00_master` — type

Private symbols:
- `append_pending` — subroutine
- `assign_optics_fork` — type
- `close_child_pipe_fds` — subroutine
- `close_fd_silent` — subroutine
- `close_pipe_except_fd` — subroutine
- `exec_stream_p00_master` — subroutine
- `init_cline_assign_optics` — subroutine
- `init_cline_opening2D` — subroutine
- `init_cline_particle_sieving` — subroutine
- `init_cline_pool2D` — subroutine
- `init_cline_preprocess` — subroutine
- `init_cline_reference_picking` — subroutine
- `init_ipc_pipe` — subroutine
- `init_metadata_assign_optics` — subroutine
- `init_metadata_opening2D` — subroutine
- `init_metadata_particle_sieving` — subroutine
- `init_metadata_pool2D` — subroutine
- `init_metadata_preprocess` — subroutine
- `init_metadata_reference_picking` — subroutine
- `initial_analysis_fork` — type
- `kill_ipc_pipe` — subroutine
- `metadata_listener` — subroutine
- `particle_sieving_fork` — type
- `pool2D_fork` — type
- `preprocess_fork` — type
- `reference_picking_fork` — type
- `send_framed_to_pipe` — subroutine
- `send_update_to_stage_pipes` — subroutine
- `sigint_handler` — subroutine
- `sigterm_handler` — subroutine
- `try_extract_framed_message` — subroutine
- `try_read_from_fd` — subroutine
- `xassign_optics` — subroutine
- `xinitial_analysis` — subroutine
- `xparticle_sieving` — subroutine
- `xpool2D` — subroutine
- `xpreprocess` — subroutine
- `xreference_picking` — subroutine

---
## Module: simple_stream_p01_preprocess_new

Files:
- `main/stream/simple_stream_p01_preprocess_new.f90`

Uses:
- `simple_gui_metadata_api`
- `simple_gui_metadata_utils`
- `simple_histogram`
- `simple_motion_correct_utils`
- `simple_stream_api`
- `simple_stream_state`
- `unix`

Public symbols:
- `stream_p01_preprocess` — type

Private symbols:
- `create_movies_set_project` — subroutine
- `exec_stream_p01_preprocess` — subroutine
- `import_previous_projects` — subroutine
- `mics_window_stats` — subroutine
- `send_to_preprocess_in_pipe` — subroutine
- `set_stat_thresholds` — subroutine
- `sigterm_handler` — subroutine
- `test_stat_thresholds` — subroutine
- `update_projects_list` — subroutine
- `write_mic_star_and_field` — subroutine

---
## Module: simple_stream_p02_assign_optics_new

Files:
- `main/stream/simple_stream_p02_assign_optics_new.f90`

Uses:
- `simple_gui_metadata_api`
- `simple_stream_api`
- `simple_stream_state`
- `unix`

Public symbols:
- `stream_p02_assign_optics` — type

Private symbols:
- `exec_stream_p02_assign_optics` — subroutine
- `send_to_assign_optics_in_pipe` — subroutine
- `sigterm_handler` — subroutine
- `wait_for_stream_folder` — subroutine

---
## Module: simple_stream_p03_initial_analysis

Files:
- `main/stream/simple_stream_p03_initial_analysis.f90`

Uses:
- `simple_abinitio_utils`
- `simple_cavg_quality_analysis`
- `simple_cavg_quality_model`
- `simple_cavg_quality_types`
- `simple_commanders_abinitio`
- `simple_commanders_abinitio2d`
- `simple_commanders_cavgs`
- `simple_commanders_denoise`
- `simple_commanders_mkcavgs`
- `simple_commanders_pick`
- `simple_commanders_reproject`
- `simple_commanders_sieve`
- `simple_defs`
- `simple_fileio`
- `simple_gui_metadata_api`
- `simple_image_msk`
- `simple_imgarr_utils`
- `simple_mini_stream_utils`
- `simple_procimgstk`
- `simple_projfile_utils`
- `simple_qsys_env`
- `simple_stream_api`
- `simple_stream_state`
- `unix`

Public symbols:
- `stream_p03_initial_analysis` — type

Private symbols:
- `balance_classes` — subroutine
- `cleanup_previous` — subroutine
- `duplicate_balanced_stack` — subroutine
- `exec_stream_p03_initial_analysis` — subroutine
- `micimporter` — subroutine
- `reestimate_box_size_from_selected_cavgs` — subroutine
- `restart_cleanup_allocatables` — subroutine
- `run_abinitio2D` — subroutine
- `run_abinitio3D_and_reproject` — subroutine
- `run_cavg_quality_selection` — subroutine
- `run_cavg_size_selection` — subroutine
- `run_cls_split` — subroutine
- `run_extract` — subroutine
- `run_make_cavgs` — subroutine
- `run_reextract` — subroutine
- `run_sieve_ptcls` — subroutine
- `send_available_cavgs2D` — subroutine
- `send_cavg2D_meta` — subroutine
- `send_meta` — subroutine
- `send_meta2D` — subroutine
- `send_micrograph_meta` — subroutine
- `send_pickref_meta` — subroutine
- `send_selected_pickrefs` — subroutine
- `send_to_initial_analysis_in_pipe` — subroutine
- `sigterm_handler` — subroutine
- `update_os_out_stk` — subroutine
- `write_quality_stack` — subroutine

---
## Module: simple_stream_p04_refpick_extract_new

Files:
- `main/stream/simple_stream_p04_refpick_extract_new.f90`

Uses:
- `simple_commanders_pick`
- `simple_gui_metadata_api`
- `simple_stream_api`
- `simple_stream_state`
- `simple_timer`
- `unix`

Public symbols:
- `stream_p04_refpick_extract` — type

Private symbols:
- `create_individual_project` — subroutine
- `exec_stream_pick_extract` — subroutine
- `import_previous_mics` — subroutine
- `send_cavg2D_meta` — subroutine
- `send_meta` — subroutine
- `send_micrograph_meta` — subroutine
- `send_pickrefs` — subroutine
- `send_to_refpick_in_pipe` — subroutine
- `sigterm_handler` — subroutine
- `update_projects_list` — subroutine
- `validate_ptcl2D_star_inputs` — subroutine
- `write_micrographs_starfile` — subroutine
- `write_project` — subroutine

---
## Module: simple_stream_p05_sieve_cavgs_new

Files:
- `main/stream/simple_stream_p05_sieve_cavgs_new.f90`

Uses:
- `simple_fileio`
- `simple_gui_metadata_api`
- `simple_ptcl_sieve`
- `simple_stream_api`
- `simple_stream_pool2d_utils`
- `simple_stream_state`
- `unix`

Public symbols:
- `stream_p05_sieve_cavgs` — type

Private symbols:
- `exec_stream_p05_sieve_cavgs` — subroutine
- `send_cavg2D_meta` — subroutine
- `send_cavgs2D_batch` — subroutine
- `send_meta` — subroutine
- `send_to_sieve_cavgs_in_pipe` — subroutine
- `sigterm_handler` — subroutine

---
## Module: simple_stream_p06_pool2D_new

Files:
- `main/stream/simple_stream_p06_pool2D_new.f90`

Uses:
- `simple_gui_metadata_api`
- `simple_gui_metadata_utils`
- `simple_stream2d_state`
- `simple_stream_api`
- `simple_stream_pool2d_utils`
- `simple_stream_state`
- `unix`

Public symbols:
- `stream_p06_pool2D` — type

Private symbols:
- `cleanup4restart` — subroutine
- `exec_stream_p06_pool2D` — subroutine
- `import_sets_into_pool` — subroutine
- `send_cavg2D_meta` — subroutine
- `send_cavgs2D` — subroutine
- `send_meta` — subroutine
- `send_meta_snapshot2D` — subroutine
- `send_to_pool2D_in_pipe` — subroutine
- `sigterm_handler` — subroutine
- `unpause_pool` — subroutine

---
## Module: simple_stream_pool2D_utils

Files:
- `main/stream/simple_stream_pool2D_utils.f90`

Uses:
- `simple_classaverager`
- `simple_euclid_sigma2`
- `simple_procimgstk`
- `simple_ran_tabu`
- `simple_stream_api`

Public symbols:
- `generate_pool_stats` — subroutine
- `get_pool_ptr` — subroutine
- `init_pool_clustering` — subroutine
- `iterate_pool` — subroutine
- `set_lpthres_type` — subroutine
- `set_pool_resolution_limits` — subroutine
- `update_match_class_states` — subroutine
- `update_mskdiam` — subroutine
- `update_pool` — subroutine
- `update_pool_aln_params` — subroutine
- `update_pool_status` — subroutine

Private symbols:
- `biased_stack_sampling` — subroutine
- `generate_pool_jpeg` — subroutine
- `set_iteration_time` — subroutine
- `set_pool_dimensions` — subroutine
- `uniform_stack_sampling` — subroutine
- `update_match_class_states_in_pool` — subroutine
- `update_pool_dims` — subroutine
- `update_pool_for_gui` — subroutine

---
## Module: simple_stream_state

Files:
- `main/stream/simple_stream_state.f90`

---
## Module: simple_stream_utils

Files:
- `main/stream/simple_stream_utils.f90`

Uses:
- `json_kinds`
- `json_module`
- `simple_cmdline`
- `simple_core_module_api`
- `simple_default_clines`
- `simple_gui_utils`
- `simple_image`
- `simple_image_msk`
- `simple_parameters`
- `simple_qsys_env`
- `simple_rec_list`
- `simple_sp_project`
- `simple_stack_io`
- `simple_stream_communicator`

Public symbols:
- `create_stream_project` — subroutine
- `get_latest_optics_map_id` — function
- `import_new_projects` — subroutine
- `init_stream_qenv` — subroutine
- `process_selected_refs` — subroutine
- `process_selected_refs_2` — subroutine
- `stream_datestr` — function
- `terminate_stream` — subroutine
- `update_user_params` — subroutine
- `wait_for_folder` — subroutine
- `wait_for_folder2` — subroutine

---
## Module: simple_stream_watcher

Files:
- `main/stream/simple_stream_watcher.f90`

Uses:
- `simple_core_module_api`
- `simple_progress`

---
## Module: simple_string

Files:
- `utils/simple_string.f90`

Uses:
- `simple_defs`
- `simple_defs_string`
- `simple_error`

---
## Module: simple_string_tester

Files:
- `utils/simple_string_tester.f90`

Uses:
- `simple_defs`
- `simple_string`
- `simple_test_utils`

Public symbols:
- `run_all_string_tests` — subroutine

Private symbols:
- `test_append_operator` — subroutine
- `test_assign_and_ctor` — subroutine
- `test_default_and_kill` — subroutine
- `test_ends_with` — subroutine
- `test_equality_operators` — subroutine
- `test_has_substr` — subroutine
- `test_is_blank_and_is_allocated` — subroutine
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
- `utils/simple_string_utils.f90`

Uses:
- `simple_defs`
- `simple_defs_string`
- `simple_error`
- `simple_string`

---
## Module: simple_sym

Files:
- `main/simple_sym.f90`

Uses:
- `simple_ori`
- `simple_ori_api`
- `simple_oris`

---
## Module: simple_symanalyzer

Files:
- `main/simple_symanalyzer.f90`

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
- `find_symaxis` — subroutine
- `symaverage` — subroutine

---
## Module: simple_syslib

Files:
- `fileio/simple_syslib.f90`

Uses:
- `iso_c_binding`
- `simple_defs`
- `simple_defs_fname`
- `simple_error`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `del_file` — subroutine
- `exec_cmdline` — subroutine
- `find_next_int_dir_prefix` — function
- `fsync` — function
- `get_absolute_pathname` — function
- `get_current_rss_bytes` — function
- `get_peak_rss_bytes` — function
- `get_sysinfo` — function
- `isdir` — function
- `list_dirs` — function
- `makedir` — function
- `memory_monitor_report_c` — function
- `memory_monitor_start_c` — function
- `memory_monitor_stop_c` — function
- `mkdir` — function
- `print_slurm_env` — subroutine
- `raise_sys_error` — subroutine
- `redirect_stdout_stderr` — subroutine
- `regexp_match` — function
- `removedir` — function
- `restore_stdout_stderr` — subroutine
- `rmdir` — function
- `safe_flush` — subroutine
- `simple_abspath` — function
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
- `simple_touch` — subroutine
- `sscanf` — function
- `symlink` — function
- `syslib_c2fortran_string` — subroutine
- `syslib_symlink` — subroutine
- `touch` — function
- `unlink` — function
- `wait_pid` — function

---
## Module: simple_syslib_tester

Files:
- `fileio/simple_syslib_tester.f90`

Uses:
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
- `test_simple_abspath` — subroutine
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
- `tent_smooth_3d` — subroutine

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
- `simple_test_exec_fft`
- `simple_test_exec_geometry`
- `simple_test_exec_highlevel`
- `simple_test_exec_io`
- `simple_test_exec_masks`
- `simple_test_exec_network`
- `simple_test_exec_numerics`
- `simple_test_exec_optimize`
- `simple_test_exec_parallel`
- `simple_test_exec_single`
- `simple_test_exec_stats`
- `simple_test_exec_utils`
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
## Module: simple_test_exec_fft

Files:
- `main/exec/simple_test_exec_fft.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_fft`

Public symbols:
- `exec_test_fft_commander` — subroutine

---
## Module: simple_test_exec_geometry

Files:
- `main/exec/simple_test_exec_geometry.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_geometry`

Public symbols:
- `exec_test_geometry_commander` — subroutine

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
## Module: simple_test_exec_io

Files:
- `main/exec/simple_test_exec_io.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_io`

Public symbols:
- `exec_test_io_commander` — subroutine

---
## Module: simple_test_exec_masks

Files:
- `main/exec/simple_test_exec_masks.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_masks`

Public symbols:
- `exec_test_masks_commander` — subroutine

---
## Module: simple_test_exec_network

Files:
- `main/exec/simple_test_exec_network.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_network`

Public symbols:
- `exec_test_network_commander` — subroutine

---
## Module: simple_test_exec_numerics

Files:
- `main/exec/simple_test_exec_numerics.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_numerics`

Public symbols:
- `exec_test_numerics_commander` — subroutine

---
## Module: simple_test_exec_optimize

Files:
- `main/exec/simple_test_exec_optimize.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_optimize`

Public symbols:
- `exec_test_optimize_commander` — subroutine

---
## Module: simple_test_exec_parallel

Files:
- `main/exec/simple_test_exec_parallel.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_parallel`

Public symbols:
- `exec_test_parallel_commander` — subroutine

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
## Module: simple_test_exec_stats

Files:
- `main/exec/simple_test_exec_stats.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_stats`

Public symbols:
- `exec_test_stats_commander` — subroutine

---
## Module: simple_test_exec_utils

Files:
- `main/exec/simple_test_exec_utils.f90`

Uses:
- `simple_cmdline`
- `simple_commanders_test_utils`

Public symbols:
- `exec_test_utils_commander` — subroutine

---
## Module: simple_test_ui_class

Files:
- `main/ui/simple_test/simple_test_ui_class.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_class_programs` — subroutine
- `new_strategy2D` — subroutine
- `new_ui_hash_test` — subroutine
- `new_units` — subroutine
- `print_test_class_programs` — subroutine

---
## Module: simple_test_ui_fft

Files:
- `main/ui/simple_test/simple_test_ui_fft.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_fft_programs` — subroutine
- `new_corrs2weights_test` — subroutine
- `new_eval_polarftcc` — subroutine
- `new_ft_expanded` — subroutine
- `new_gencorrs_fft` — subroutine
- `new_order_corr` — subroutine
- `new_phasecorr` — subroutine
- `new_rank_weights` — subroutine
- `new_rotate_ref` — subroutine
- `print_test_fft_programs` — subroutine

---
## Module: simple_test_ui_geometry

Files:
- `main/ui/simple_test/simple_test_ui_geometry.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_geometry_programs` — subroutine
- `new_angres` — subroutine
- `new_ori_test` — subroutine
- `new_oris_test` — subroutine
- `new_sym_test` — subroutine
- `new_uniform_euler` — subroutine
- `new_uniform_rot` — subroutine
- `print_test_geometry_programs` — subroutine

---
## Module: simple_test_ui_highlevel

Files:
- `main/ui/simple_test/simple_test_ui_highlevel.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_highlevel_programs` — subroutine
- `new_flex_preimage_basis_ab` — subroutine
- `new_flex_preimage_identity` — subroutine
- `new_mini_stream` — subroutine
- `new_pcg_recon_ctf_free` — subroutine
- `new_pcg_recon_ctf_hetero` — subroutine
- `new_pcg_recon_deapod` — subroutine
- `new_pcg_recon_kernel` — subroutine
- `new_ptcls_ppca_subproject_distr` — subroutine
- `new_reproject` — subroutine
- `new_simulate_particles` — subroutine
- `new_simulated_workflow` — subroutine
- `new_subproject_distr` — subroutine
- `print_test_highlevel_programs` — subroutine

---
## Module: simple_test_ui_io

Files:
- `main/ui/simple_test/simple_test_ui_io.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_io_programs` — subroutine
- `new_imgfile` — subroutine
- `new_inside_write` — subroutine
- `new_io` — subroutine
- `new_io_parallel` — subroutine
- `new_mrc2jpeg` — subroutine
- `new_mrc_validation` — subroutine
- `new_stack_io` — subroutine
- `new_star_export` — subroutine
- `new_starfile_test` — subroutine
- `print_test_io_programs` — subroutine

---
## Module: simple_test_ui_masks

Files:
- `main/ui/simple_test/simple_test_ui_masks.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_masks_programs` — subroutine
- `new_bounds_from_mask3D_test` — subroutine
- `new_graphene_mask` — subroutine
- `new_image_bin` — subroutine
- `new_mask` — subroutine
- `new_msk_routines` — subroutine
- `new_nano_mask` — subroutine
- `new_otsu_test` — subroutine
- `new_ptcl_center` — subroutine
- `print_test_masks_programs` — subroutine

---
## Module: simple_test_ui_network

Files:
- `main/ui/simple_test/simple_test_ui_network.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_network_programs` — subroutine
- `new_socket_client` — subroutine
- `new_socket_comm_distr` — subroutine
- `new_socket_io` — subroutine
- `new_socket_server` — subroutine
- `print_test_network_programs` — subroutine

---
## Module: simple_test_ui_numerics

Files:
- `main/ui/simple_test/simple_test_ui_numerics.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_numerics_programs` — subroutine
- `new_eigh_test` — subroutine
- `new_kbinterpol_fast` — subroutine
- `new_maxnloc_test` — subroutine
- `new_neigh` — subroutine
- `print_test_numerics_programs` — subroutine

---
## Module: simple_test_ui_optimize

Files:
- `main/ui/simple_test/simple_test_ui_optimize.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_optimize_programs` — subroutine
- `new_lbfgsb` — subroutine
- `new_lbfgsb_cosine` — subroutine
- `new_lplims` — subroutine
- `new_lpstages_test` — subroutine
- `new_opt_lp` — subroutine
- `print_test_optimize_programs` — subroutine

---
## Module: simple_test_ui_parallel

Files:
- `main/ui/simple_test/simple_test_ui_parallel.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_parallel_programs` — subroutine
- `new_coarrays` — subroutine
- `new_openacc` — subroutine
- `new_openmp` — subroutine
- `new_simd` — subroutine
- `print_test_parallel_programs` — subroutine

---
## Module: simple_test_ui_single

Files:
- `main/ui/simple_test/simple_test_ui_single.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_single_programs` — subroutine
- `new_atoms_stats` — subroutine
- `new_detect_atoms` — subroutine
- `new_simulate_nanoparticle` — subroutine
- `new_single_workflow` — subroutine
- `print_test_single_programs` — subroutine

---
## Module: simple_test_ui_stats

Files:
- `main/ui/simple_test/simple_test_ui_stats.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_stats_programs` — subroutine
- `new_class_sample_test` — subroutine
- `new_clustering` — subroutine
- `new_ctf_test` — subroutine
- `new_eo_diff` — subroutine
- `new_extr_frac` — subroutine
- `new_multinomal_test` — subroutine
- `new_pca_all` — subroutine
- `new_pca_imgvar` — subroutine
- `new_sp_project` — subroutine
- `print_test_stats_programs` — subroutine

---
## Module: simple_test_ui_utils

Files:
- `main/ui/simple_test/simple_test_ui_utils.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_test_utils_programs` — subroutine
- `new_ansi_colors` — subroutine
- `new_binoris_io_test` — subroutine
- `new_binoris_test` — subroutine
- `new_cif2mrc` — subroutine
- `new_cif2pdb` — subroutine
- `new_cmdline` — subroutine
- `new_install` — subroutine
- `new_nice` — subroutine
- `new_pdb2mrc` — subroutine
- `new_peak_thres_fdr` — subroutine
- `new_serialize` — subroutine
- `new_stringmatch` — subroutine
- `print_test_utils_programs` — subroutine

---
## Module: simple_test_utils

Files:
- `utils/simple_test_utils.f90`

Uses:
- `simple_defs`
- `simple_string`

Public symbols:
- `assert_char` — subroutine
- `assert_double` — subroutine
- `assert_false` — subroutine
- `assert_int` — subroutine
- `assert_real` — subroutine
- `assert_string_eq` — subroutine
- `assert_true` — subroutine
- `report_summary` — subroutine

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
- `TIFFGetField` — function
- `TIFFGetVersion` — function
- `TIFFIsTiled` — function
- `TIFFLastDirectory` — function
- `TIFFmalloc` — function
- `TIFFMuteWarnings` — subroutine
- `TIFFNumberOfStrips` — function
- `TIFFOpen` — function
- `TIFFPrintInfo` — subroutine
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
- `now` — subroutine
- `reset_timer` — subroutine
- `timer_loop_end` — subroutine
- `timer_loop_start` — subroutine
- `timer_profile_break` — subroutine
- `timer_profile_report` — subroutine
- `timer_profile_setup` — subroutine
- `timer_profile_start` — subroutine
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

---
## Module: simple_trajectory_chunker

Files:
- `main/nano/simple_trajectory_chunker.f90`

Uses:
- `simple_core_module_api`
- `simple_flex_analysis_strategy`
- `simple_srch_sort_loc`

Public symbols:
- `make_trajectory_chunk_plan` — subroutine
- `select_trajectory_chunk_plan` — subroutine
- `trajectory_chunk` — type
- `trajectory_chunk_plan` — type
- `trajectory_chunks_to_parts` — subroutine
- `write_trajectory_chunks_csv` — subroutine

Private symbols:
- `kill_trajectory_chunk_plan` — subroutine
- `make_boundary_bands` — subroutine
- `prepare_weighted_features` — subroutine

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

Private symbols:
- `calc_3bands_mask` — function
- `interp_peak` — subroutine
- `obscure_peak` — subroutine
- `range_convention` — subroutine
- `remove_lattice` — subroutine

---
## Module: simple_type_defs

Files:
- `defs/simple_type_defs.f90`

Uses:
- `simple_defs`
- `simple_string`

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
- `simple_ui_params_common`
- `simple_ui_program`
- `simple_ui_simple_group`
- `simple_ui_single_group`
- `simple_ui_stream_group`
- `simple_ui_test_group`

Public symbols:
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
- `write_ui_json` — subroutine

Private symbols:
- `create_program_entry` — subroutine
- `create_program_entry` — subroutine
- `create_section_list` — subroutine
- `create_section_list` — subroutine

---
## Module: simple_ui_abinitio3D

Files:
- `main/ui/simple/simple_ui_abinitio3D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_abinitio3D_programs` — subroutine
- `new_abinitio3D` — subroutine
- `new_abinitio3D_cavgs` — subroutine
- `new_abinitio3D_cavgs_reject` — subroutine
- `new_estimate_lpstages` — subroutine
- `new_noisevol` — subroutine
- `print_abinitio3D_programs` — subroutine

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
- `print_cavgproc_programs` — subroutine

---
## Module: simple_ui_cluster2D

Files:
- `main/ui/simple/simple_ui_cluster2D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_cluster2D_programs` — subroutine
- `new_abinitio2D` — subroutine
- `new_abinitio2D_chunks` — subroutine
- `new_bootstrap_cavgs` — subroutine
- `new_make_cavgs` — subroutine
- `new_map_cavgs_selection` — subroutine
- `new_sample_classes` — subroutine
- `new_unbootstrap_cavgs` — subroutine
- `new_write_classes` — subroutine
- `print_cluster2D_programs` — subroutine

---
## Module: simple_ui_denoise

Files:
- `main/ui/simple/simple_ui_denoise.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_denoise_programs` — subroutine
- `new_cls_split` — subroutine
- `new_denoise_project` — subroutine
- `new_flex_analysis` — subroutine
- `new_icm2D` — subroutine
- `new_icm3D` — subroutine
- `new_map_params_from_den` — subroutine
- `new_ppca_denoise` — subroutine
- `new_ppca_denoise_classes` — subroutine
- `new_ppca_volvar` — subroutine
- `print_denoise_programs` — subroutine

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
- `print_dock_programs` — subroutine

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
- `print_filter_programs` — subroutine

---
## Module: simple_ui_hash

Files:
- `main/ui/simple_ui_hash.f90`

Uses:
- `simple_string`
- `simple_ui_program`
- `simple_vrefhash`

Public symbols:
- `test_ui_hash` — subroutine
- `ui_hash` — type

Private symbols:
- `get_ref_ui_param_char` — subroutine
- `get_ref_ui_param_str` — subroutine
- `get_ref_ui_program_char` — subroutine
- `get_ref_ui_program_str` — subroutine
- `set_ref_ui_param_char` — subroutine
- `set_ref_ui_param_str` — subroutine
- `set_ref_ui_program_char` — subroutine
- `set_ref_ui_program_str` — subroutine

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
- `new_ctf_phaseflip` — subroutine
- `new_ctfops` — subroutine
- `new_normalize` — subroutine
- `new_scale` — subroutine
- `new_stack` — subroutine
- `new_stackops` — subroutine
- `print_image_programs` — subroutine

---
## Module: simple_ui_mask

Files:
- `main/ui/simple/simple_ui_mask.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_mask_programs` — subroutine
- `new_auto_spher_mask` — subroutine
- `new_automask2D` — subroutine
- `new_mask` — subroutine
- `print_mask_programs` — subroutine

---
## Module: simple_ui_modules

Files:
- `main/ui/simple_ui_modules.f90`

Uses:
- `simple_ansi_ctrls`
- `simple_core_module_api`
- `simple_ui_hash`
- `simple_ui_params_common`
- `simple_ui_program`
- `simple_ui_utils`

---
## Module: simple_ui_ori

Files:
- `main/ui/simple/simple_ui_ori.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_ori_programs` — subroutine
- `new_make_oris` — subroutine
- `new_oriops` — subroutine
- `new_oristats` — subroutine
- `new_vizoris` — subroutine
- `print_ori_programs` — subroutine

---
## Module: simple_ui_other

Files:
- `main/ui/simple/simple_ui_other.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_other_programs` — subroutine
- `new_cif2pdb` — subroutine
- `new_fractionate_movies` — subroutine
- `new_split_` — subroutine
- `new_split_stack` — subroutine
- `print_other_programs` — subroutine

---
## Module: simple_ui_param

Files:
- `main/ui/simple_ui_param.f90`

Uses:
- `simple_string`

Public symbols:
- `apply_gui_overrides` — subroutine
- `finalize` — subroutine
- `set_param_1` — subroutine
- `set_param_2` — subroutine

---
## Module: simple_ui_params_common

Files:
- `main/ui/simple_ui_params_common.f90`

Uses:
- `simple_core_module_api`
- `simple_ui_param`

Public symbols:
- `set_ui_params` — subroutine

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
- `print_preproc_programs` — subroutine

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
- `print_print_programs` — subroutine

---
## Module: simple_ui_program

Files:
- `main/ui/simple_ui_program.f90`

Uses:
- `json_module`
- `simple_ansi_ctrls`
- `simple_core_module_api`
- `simple_linked_list`
- `simple_ui_param`

Public symbols:
- `add_input_num` — subroutine
- `add_input_param` — subroutine
- `add_input_str` — subroutine
- `append_required_keys_from_list` — subroutine
- `create_section_from_list` — subroutine
- `get_executable` — function
- `get_name` — function
- `get_required_keys` — function
- `kill` — subroutine
- `new` — subroutine
- `print_cmdline` — subroutine
- `print_param_hash` — subroutine
- `print_param_list` — subroutine
- `print_prg_descr_long` — subroutine
- `print_ui` — subroutine
- `ui_program` — type
- `write2json` — subroutine

---
## Module: simple_ui_project

Files:
- `main/ui/simple/simple_ui_project.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_project_programs` — subroutine
- `new_export_manifoldem_starproject` — subroutine
- `new_export_relion` — subroutine
- `new_export_starproject` — subroutine
- `new_extract_subproj` — subroutine
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
- `new_ptcl3D_state_consensus` — subroutine
- `new_reimport_particles` — subroutine
- `new_replace_project_field` — subroutine
- `new_selection` — subroutine
- `new_update_project` — subroutine
- `new_validate_projfile` — subroutine
- `new_write_mic_filetab` — subroutine
- `new_zero_project_shifts` — subroutine
- `print_project_programs` — subroutine

---
## Module: simple_ui_refine3D

Files:
- `main/ui/simple/simple_ui_refine3D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_refine3D_programs` — subroutine
- `new_automask` — subroutine
- `new_bootstrap_rec3D` — subroutine
- `new_postprocess` — subroutine
- `new_reconstruct3D` — subroutine
- `new_refine3D` — subroutine
- `new_refine3D_auto` — subroutine
- `new_refine3D_multi` — subroutine
- `print_refine3D_programs` — subroutine

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
- `print_res_programs` — subroutine

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
- `print_sim_programs` — subroutine

---
## Module: simple_ui_simple_group

Files:
- `main/ui/simple_ui_simple_group.f90`

Uses:
- `simple_ui_abinitio3d`
- `simple_ui_cavgproc`
- `simple_ui_cluster2d`
- `simple_ui_denoise`
- `simple_ui_dock`
- `simple_ui_filter`
- `simple_ui_hash`
- `simple_ui_image`
- `simple_ui_mask`
- `simple_ui_ori`
- `simple_ui_other`
- `simple_ui_preproc`
- `simple_ui_print`
- `simple_ui_project`
- `simple_ui_refine3d`
- `simple_ui_res`
- `simple_ui_sim`
- `simple_ui_sym`
- `simple_ui_validate`
- `simple_ui_volume`

Public symbols:
- `add_simple_programs` — subroutine
- `print_simple_programs` — subroutine

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
- `print_single_programs` — subroutine

---
## Module: simple_ui_stream

Files:
- `main/ui/simple/simple_ui_stream.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_stream_programs` — subroutine
- `new_abinitio2D_stream` — subroutine
- `new_assign_optics` — subroutine
- `new_gen_pickrefs` — subroutine
- `new_master` — subroutine
- `new_pick_extract` — subroutine
- `new_preproc` — subroutine
- `new_sieve_cavgs` — subroutine
- `print_stream_programs` — subroutine

---
## Module: simple_ui_stream_group

Files:
- `main/ui/simple_ui_stream_group.f90`

Uses:
- `simple_ui_hash`
- `simple_ui_stream`

Public symbols:
- `add_stream_programs` — subroutine
- `print_stream_programs_group` — subroutine

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
- `print_symmetry_programs` — subroutine

---
## Module: simple_ui_test_group

Files:
- `main/ui/simple_ui_test_group.f90`

Uses:
- `simple_test_ui_class`
- `simple_test_ui_fft`
- `simple_test_ui_geometry`
- `simple_test_ui_highlevel`
- `simple_test_ui_io`
- `simple_test_ui_masks`
- `simple_test_ui_network`
- `simple_test_ui_numerics`
- `simple_test_ui_optimize`
- `simple_test_ui_parallel`
- `simple_test_ui_single`
- `simple_test_ui_stats`
- `simple_test_ui_utils`
- `simple_ui_hash`

Public symbols:
- `add_test_programs` — subroutine
- `print_test_programs` — subroutine

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
- `print_validate_programs` — subroutine

---
## Module: simple_ui_volume

Files:
- `main/ui/simple/simple_ui_volume.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_volume_programs` — subroutine
- `new_center` — subroutine
- `new_reconstruct3D_pcg` — subroutine
- `new_reproject` — subroutine
- `new_volops` — subroutine
- `print_volume_programs` — subroutine

---
## Module: simple_vol_pproc_policy

Files:
- `main/volume/simple_vol_pproc_policy.f90`

Uses:
- `simple_core_module_api`
- `simple_parameters`

Public symbols:
- `plan_state_postprocess` — subroutine
- `state_mask_is_compatible` — subroutine
- `vol_pproc_plan` — type

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
- `volpft_symsrch_costfun` — function
- `volpft_symsrch_kill` — subroutine
- `volpft_symsrch_scorefun` — function

---
## Module: simple_vrefhash

Files:
- `utils/structs/simple_vrefhash.f90`

Uses:
- `simple_error`
- `simple_string`
- `simple_string_utils`

Public symbols:
- `vrefhash` — type
- `vrefhash_bucket` — type
- `vrefhash_node` — type

Private symbols:
- `clear` — subroutine
- `del_char` — subroutine
- `del_str` — subroutine
- `destroy` — subroutine
- `finalize_node` — subroutine
- `find_node` — subroutine
- `get_ref_char` — subroutine
- `get_ref_str` — subroutine
- `init` — subroutine
- `keys` — function
- `keys_sorted` — function
- `set_ref_char` — subroutine
- `set_ref_str` — subroutine

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
## Module: simple_winfuns

Files:
- `main/interp/simple_winfuns.f90`

Uses:
- `simple_defs`

Public symbols:
- `winfuns` — type

Private symbols:
- `ifun` — function
- `wfun` — function

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
- `simple_commanders_cluster2d`
- `simple_commanders_imgproc`
- `simple_commanders_mkcavgs`
- `simple_commanders_sim`

Public symbols:
- `commander_analysis2D_nano` — type
- `commander_center2D_nano` — type
- `commander_cluster2D_nano` — type
- `exec_analysis2D_nano` — subroutine
- `exec_center2D_nano` — subroutine
- `exec_cluster2D_nano` — subroutine

---
## Module: single_commanders_nano3D

Files:
- `main/commanders/single/single_commanders_nano3D.f90`

Uses:
- `simple_abinitio_utils`
- `simple_commanders_abinitio`
- `simple_commanders_api`
- `simple_commanders_atoms`
- `simple_commanders_cluster2d`
- `simple_commanders_ori`
- `simple_commanders_rec`
- `simple_commanders_refine3d`
- `simple_commanders_reproject`
- `simple_flex_analysis_strategy`
- `simple_nanoparticle`
- `simple_refine3d_fnames`
- `simple_trajectory_chunker`

Public symbols:
- `commander_abinitio3D_nano` — type
- `commander_autorefine3D_nano` — type
- `commander_refine3D_nano` — type
- `commander_trajectory_reconstruct3D_distr` — type
- `exec_abinitio3D_nano` — subroutine
- `exec_autorefine3D_nano` — subroutine
- `exec_commander_trajectory_reconstruct3D_distr` — subroutine
- `exec_refine3D_nano` — subroutine

---
## Module: single_commanders_trajectory

Files:
- `main/commanders/single/single_commanders_trajectory.f90`

Uses:
- `simple_commanders_api`
- `simple_commanders_imgops`
- `simple_tseries_graphene_subtr`
- `single_tseries_tracker`

Public symbols:
- `commander_extract_substk` — type
- `commander_graphene_subtr` — type
- `commander_import_trajectory` — type
- `commander_track_particles` — type
- `commander_track_particles_distr` — type
- `commander_trajectory_backgr_subtr` — type
- `commander_trajectory_denoise` — type
- `commander_trajectory_swap_stack` — type
- `exec_extract_substk` — subroutine
- `exec_graphene_subtr` — subroutine
- `exec_import_trajectory` — subroutine
- `exec_track_particles` — subroutine
- `exec_track_particles_distr` — subroutine
- `exec_trajectory_backgr_subtr` — subroutine
- `exec_trajectory_denoise` — subroutine
- `exec_trajectory_swap_stack` — subroutine

---
## Module: single_commanders_tseries

Files:
- `main/commanders/single/single_commanders_tseries.f90`

Uses:
- `simple_bspline_smoother`
- `simple_commanders_api`
- `simple_denoise_movies`
- `simple_imgarr_utils`
- `simple_motion_correct_iter`
- `single_tseries_extractor`

Public symbols:
- `commander_tseries_extractor` — type
- `commander_tseries_import` — type
- `commander_tseries_make_pickavg` — type
- `commander_tseries_motion_correct` — type
- `commander_tseries_motion_correct_distr` — type
- `commander_tseries_prep4tracking` — type
- `downscale_frames` — subroutine
- `exec_tseries_extractor` — subroutine
- `exec_tseries_import` — subroutine
- `exec_tseries_make_pickavg` — subroutine
- `exec_tseries_motion_correct` — subroutine
- `exec_tseries_motion_correct_distr` — subroutine
- `exec_tseries_prep4tracking` — subroutine
- `find_outliers` — subroutine

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
## Module: single_tseries_extractor

Files:
- `main/nano/single_tseries_extractor.f90`

Uses:
- `simple_core_module_api`
- `simple_image`
- `simple_parameters`
- `simple_sp_project`
- `simple_stack_io`

Public symbols:
- `extract_trajectory` — subroutine
- `init_trajectory_extractor` — subroutine
- `kill_trajectory_extractor` — subroutine

Private symbols:
- `build_mapping` — subroutine
- `identify_neighbours` — subroutine
- `update_background_pspec` — subroutine

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

Private symbols:
- `center_reference` — function
- `cleanup_centering_images` — subroutine
- `identify_neighbours` — subroutine
- `update_background_pspec` — subroutine
- `write_tester_dump` — subroutine
- `write_trajectory` — subroutine

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
- `new_detect_atoms` — subroutine
- `new_simulate_nanoparticle` — subroutine
- `print_single_atom_programs` — subroutine

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
- `print_single_map_programs` — subroutine

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
- `new_cluster2D_nano` — subroutine
- `new_estimate_diam` — subroutine
- `print_single_nano2D_programs` — subroutine

---
## Module: single_ui_nano3D

Files:
- `main/ui/single/single_ui_nano3D.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_nano3D_programs` — subroutine
- `new_abinitio3D_nano` — subroutine
- `new_autorefine3D_nano` — subroutine
- `new_refine3D_nano` — subroutine
- `print_single_nano3D_programs` — subroutine

---
## Module: single_ui_trajectory

Files:
- `main/ui/single/single_ui_trajectory.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_trajectory_programs` — subroutine
- `new_extract_substk` — subroutine
- `new_graphene_subtr` — subroutine
- `new_import_trajectory` — subroutine
- `new_trajectory_denoise` — subroutine
- `new_trajectory_make_projavgs` — subroutine
- `new_trajectory_reconstruct3D` — subroutine
- `new_trajectory_swap_stack` — subroutine
- `print_single_trajectory_programs` — subroutine

---
## Module: single_ui_tseries

Files:
- `main/ui/single/single_ui_tseries.f90`

Uses:
- `simple_ui_modules`

Public symbols:
- `construct_single_tseries_programs` — subroutine
- `new_track_particles` — subroutine
- `new_tseries_extractor` — subroutine
- `new_tseries_import` — subroutine
- `new_tseries_make_pickavg` — subroutine
- `new_tseries_motion_correct` — subroutine
- `new_tseries_prep4tracking` — subroutine
- `print_single_tseries_programs` — subroutine

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
- `print_single_validate_programs` — subroutine

---
## Module: subroutine

Files:
- `main/class/simple_classaverager.f90`
- `main/class/simple_classaverager_core.f90`
- `main/class/simple_classaverager_restore.f90`
- `main/image/simple_image.f90`
- `main/image/simple_image_access.f90`
- `main/image/simple_image_arith.f90`
- `main/image/simple_image_calc.f90`
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
- `main/nu_filt/simple_nu_filter.f90`
- `main/nu_filt/simple_nu_filter_apply.f90`
- `main/nu_filt/simple_nu_filter_bank.f90`
- `main/nu_filt/simple_nu_filter_extend.f90`
- `main/nu_filt/simple_nu_filter_potts.f90`
- `main/nu_filt/simple_nu_filter_state.f90`
- `main/nu_filt/simple_nu_filter_stats.f90`
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
- `main/params/simple_parameters.f90`
- `main/params/simple_parameters_core.f90`
- `main/params/simple_parameters_parse.f90`
- `main/params/simple_parameters_phases.f90`
- `main/pftc/simple_polarft_access.f90`
- `main/pftc/simple_polarft_calc.f90`
- `main/pftc/simple_polarft_core.f90`
- `main/pftc/simple_polarft_corr.f90`
- `main/pftc/simple_polarft_ctf.f90`
- `main/pftc/simple_polarft_geom.f90`
- `main/pftc/simple_polarft_memo.f90`
- `main/pftc/simple_polarft_ops_io.f90`
- `main/project/simple_sp_project.f90`
- `main/project/simple_sp_project_cls.f90`
- `main/project/simple_sp_project_core.f90`
- `main/project/simple_sp_project_io.f90`
- `main/project/simple_sp_project_mic.f90`
- `main/project/simple_sp_project_optics.f90`
- `main/project/simple_sp_project_out.f90`
- `main/project/simple_sp_project_ptcl.f90`
- `main/project/simple_sp_project_stk.f90`
- `main/simple_abinitio_controller.f90`
- `main/simple_abinitio_utils.f90`

Uses:
- `simple_commanders_euclid`
- `simple_ctf`
- `simple_euclid_sigma2`
- `simple_gpu_utils`
- `simple_gridding`
- `simple_imgarr_utils`
- `simple_jpg`
- `simple_linalg`
- `simple_map_reduce`
- `simple_matcher_ptcl_io`
- `simple_math`
- `simple_math_ctf`
- `simple_math_ft`
- `simple_memoize_ft_maps`
- `simple_parameters`
- `simple_projector`
- `simple_projector_pft_batch`
- `simple_ran_tabu`
- `simple_sp_project`
- `simple_tent_smooth`

Public symbols:
- `add_edge` — subroutine
- `add_fsc` — subroutine
- `add_oriplot` — subroutine
- `append_and_thin_nu_highres_candidate` — subroutine
- `apply_dose_weighing` — subroutine
- `assert_polarize_pdim` — subroutine
- `cache_nu_extension_frontier_dmats` — subroutine
- `cache_nu_highres_extension_frontier_after_selection` — subroutine
- `calc_class_center_shift` — subroutine
- `calc_final_rec` — subroutine
- `calc_lplim_final_stage` — function
- `calc_rec` — subroutine
- `calculate_histogram` — subroutine
- `calculate_indices` — subroutine
- `calculate_optics_plot` — subroutine
- `calculate_plot` — subroutine
- `cavger_update_sums` — subroutine
- `check_file_formats` — subroutine
- `check_vol` — subroutine
- `clear_polar_memo` — subroutine
- `collect_nu_base_label_voxels` — subroutine
- `compact_nu_highres_dmat_bank_for_capacity` — subroutine
- `dealloc_cavgs` — subroutine
- `double_check_file_formats` — subroutine
- `dt_1d` — subroutine
- `ensure_phase_shift_fields` — subroutine
- `exec_refine3D` — subroutine
- `fill_nu_frontier_dmat_from_bank` — subroutine
- `gen_c1` — subroutine
- `gen_euclid_crvec` — subroutine
- `gen_hybrid_grad_for_rot_8_local` — subroutine
- `gen_many2many_euclids_cufft_kernel` — subroutine
- `gen_ortho_reprojs4viz` — subroutine
- `generate_random_volumes` — subroutine
- `get_projections` — subroutine
- `init_gridcorr_mats` — subroutine
- `inject_refine3D_volume` — subroutine
- `loc_var_masked` — subroutine
- `log_nu_aux_replacement_margin_stats` — subroutine
- `log_nu_objective_smoothing_bank` — subroutine
- `memo_ptcls_alloc_error` — subroutine
- `memoize_polar_point` — subroutine
- `mkfnames` — subroutine
- `normalize_input_volumes` — subroutine
- `normalized_masked_value` — subroutine
- `open_pft_or_ctf2_array_for_write` — subroutine
- `pack_nu_dmat_to_mask_vector` — subroutine
- `prep_class_command_lines` — subroutine
- `prep_final_rec_cline` — subroutine
- `prepare_nu_smooth_norm` — subroutine
- `print_states` — subroutine
- `prune_nu_highres_extension_bank` — subroutine
- `randomize_states` — subroutine
- `read_cavgs` — subroutine
- `read_masks` — subroutine
- `register_stage_volume` — subroutine
- `require_valid_stats_mask` — subroutine
- `reset_mats` — subroutine
- `set_action` — subroutine
- `set_boxcoords` — subroutine
- `set_face_indices` — subroutine
- `set_indstk_from_range` — subroutine
- `set_ldim_box_from_stk` — subroutine
- `set_lplims_from_frcs` — subroutine
- `set_lplims_from_input` — subroutine
- `set_symmetry_class_vars` — subroutine
- `shift_stack_slice2D` — subroutine
- `sort_oris` — function
- `strip_refine3D_planning_keys` — subroutine
- `symmetrize` — subroutine
- `write_abinitio_lowpass_snapshot` — subroutine
- `write_final_rec_outputs` — subroutine

---
## Module: unix

Files:
- `extlibs/unix/src/unix.f90`

---
## Module: unix_dirent

Files:
- `extlibs/unix/src/unix_dirent.F90`

Public symbols:
- `c_closedir` — function
- `c_opendir` — function
- `c_readdir` — function

Private symbols:
- `c_dirent` — type
- `c_dirent` — type

---
## Module: unix_errno

Files:
- `extlibs/unix/src/unix_errno.F90`

Public symbols:
- `c_errno` — function

---
## Module: unix_fcntl

Files:
- `extlibs/unix/src/unix_fcntl.F90`

Public symbols:
- `c_fcntl` — function
- `c_open` — function

---
## Module: unix_inet

Files:
- `extlibs/unix/src/unix_inet.F90`

Public symbols:
- `c_htonl` — function
- `c_htons` — function
- `c_inet_addr` — function

---
## Module: unix_ioctl

Files:
- `extlibs/unix/src/unix_ioctl.F90`

Public symbols:
- `c_ioctl` — function

---
## Module: unix_ipc

Files:
- `extlibs/unix/src/unix_ipc.F90`

Public symbols:
- `c_ftok` — function

---
## Module: unix_mqueue

Files:
- `extlibs/unix/src/unix_mqueue.F90`

Public symbols:
- `c_mq_close` — function
- `c_mq_getattr` — function
- `c_mq_open` — function
- `c_mq_receive` — function
- `c_mq_send` — function
- `c_mq_setattr` — function
- `c_mq_timedreceive` — function
- `c_mq_timedsend` — function
- `c_mq_unlink` — function

Private symbols:
- `c_mq_attr` — type

---
## Module: unix_msg

Files:
- `extlibs/unix/src/unix_msg.F90`

Public symbols:
- `c_msgctl` — function
- `c_msgget` — function
- `c_msgrcv` — function
- `c_msgsnd` — function

---
## Module: unix_netdb

Files:
- `extlibs/unix/src/unix_netdb.F90`

Public symbols:
- `c_freeaddrinfo` — subroutine
- `c_gai_strerror` — function
- `c_getaddrinfo` — function

Private symbols:
- `c_addrinfo` — type
- `c_addrinfo` — type
- `c_in_addr` — type
- `c_sockaddr` — type
- `c_sockaddr` — type
- `c_sockaddr_in` — type
- `c_sockaddr_in` — type

---
## Module: unix_pthread

Files:
- `extlibs/unix/src/unix_pthread.F90`

Public symbols:
- `c_pthread_cancel` — function
- `c_pthread_create` — function
- `c_pthread_detach` — function
- `c_pthread_exit` — subroutine
- `c_pthread_join` — function
- `c_pthread_mutex_destroy` — function
- `c_pthread_mutex_init` — function
- `c_pthread_mutex_lock` — function
- `c_pthread_mutex_trylock` — function
- `c_pthread_mutex_unlock` — function
- `c_pthread_setcancelstate` — function
- `c_pthread_setcanceltype` — function
- `c_pthread_t` — type

Private symbols:
- `c_pthread_mutex_t` — type

---
## Module: unix_regex

Files:
- `extlibs/unix/src/unix_regex.F90`

Public symbols:
- `c_regcomp` — function
- `c_regerror` — function
- `c_regexec` — function
- `c_regfree` — subroutine

Private symbols:
- `c_regex_t` — type
- `c_regex_t` — type
- `c_regmatch_t` — type

---
## Module: unix_semaphore

Files:
- `extlibs/unix/src/unix_semaphore.F90`

Public symbols:
- `c_sem_close` — function
- `c_sem_destroy` — function
- `c_sem_getvalue` — function
- `c_sem_init` — function
- `c_sem_open` — function
- `c_sem_post` — function
- `c_sem_timedwait` — function
- `c_sem_trywait` — function
- `c_sem_unlink` — function
- `c_sem_wait` — function

Private symbols:
- `c_sem_t` — type

---
## Module: unix_signal

Files:
- `extlibs/unix/src/unix_signal.F90`

Public symbols:
- `c_kill` — function
- `c_signal` — function

---
## Module: unix_socket

Files:
- `extlibs/unix/src/unix_socket.F90`

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

Public symbols:
- `c_fstat` — function
- `c_lstat` — function
- `c_mkdir` — function
- `c_mkfifo` — function
- `c_stat` — function
- `c_umask` — function

Private symbols:
- `c_stat_type` — type
- `c_stat_type` — type
- `c_stat_type` — type

---
## Module: unix_stdio

Files:
- `extlibs/unix/src/unix_stdio.F90`

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

Public symbols:
- `c_exit` — subroutine
- `c_free` — subroutine

---
## Module: unix_string

Files:
- `extlibs/unix/src/unix_string.F90`

Public symbols:
- `c_memcpy` — function
- `c_memset` — function
- `c_strerror` — function
- `c_strlen` — function

---
## Module: unix_syslog

Files:
- `extlibs/unix/src/unix_syslog.F90`

Public symbols:
- `c_closelog` — subroutine
- `c_openlog` — subroutine
- `c_syslog` — subroutine

---
## Module: unix_termios

Files:
- `extlibs/unix/src/unix_termios.F90`

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

Private symbols:
- `c_termios` — type
- `c_termios` — type

---
## Module: unix_time

Files:
- `extlibs/unix/src/unix_time.F90`

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

Private symbols:
- `c_timespec` — type
- `c_timeval` — type
- `c_timezone` — type
- `c_tm` — type

---
## Module: unix_types

Files:
- `extlibs/unix/src/unix_types.F90`

---
## Module: unix_unistd

Files:
- `extlibs/unix/src/unix_unistd.F90`

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

Public symbols:
- `c_uname` — function

Private symbols:
- `c_utsname` — type
- `c_utsname` — type

---
## Module: unix_wait

Files:
- `extlibs/unix/src/unix_wait.F90`

Public symbols:
- `c_wait` — function
- `c_waitpid` — function

---
