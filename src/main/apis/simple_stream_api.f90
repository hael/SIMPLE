!@descr: Aggregated public API for stream_exec.
module simple_stream_api
use simple_core_module_api
use simple_defs_environment
use json_kinds
use json_module
use simple_class_frcs,             only: class_frcs
use simple_cmdline,                only: cmdline
use simple_commander_base,         only: commander_base
use simple_gui_utils,              only: mic2thumb, mrc2jpeg_tiled
use simple_image,                  only: image
use simple_parameters,             only: parameters
use simple_progress,               only: progressfile_init, progressfile_update, progress_estimate_preprocess_stream
use simple_projfile_utils,         only: merge_chunk_projfiles
use simple_qsys_env,               only: qsys_env
use simple_qsys_funs,              only: qsys_watcher, qsys_cleanup, qsys_declare_part_finished
use simple_rec_list,               only: rec, project_rec, process_rec, chunk_rec, rec_list, rec_iterator
use simple_sp_project,             only: sp_project
use simple_stack_io,               only: stack_io
use simple_starproject_stream,     only: starproject_stream
use simple_stream_refine2D_utils, only: cleanup_root_folder, tidy_2Dstream_iter
use simple_stream_utils,           only: get_latest_optics_map_id, create_stream_project, init_stream_qenv, import_new_projects
use simple_stream_watcher,         only: stream_watcher, sniff_folders_SJ, workout_directory_structure 
end module simple_stream_api
