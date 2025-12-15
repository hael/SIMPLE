module simple_stream_module_api

! DEFINE API
use simple_class_frcs,             only: class_frcs
use simple_cmdline,                only: cmdline
use simple_commander_base,         only: commander_base
use simple_euclid_sigma2,          only: average_sigma2_groups, sigma2_star_from_iter
use simple_gui_utils,              only: mic2thumb
use simple_guistats,               only: guistats
use simple_image,                  only: image
use simple_nice,                   only: simple_nice_communicator
use simple_parameters,             only: parameters, params_glob
use simple_progress,               only: progressfile_init, progressfile_update
use simple_projfile_utils,         only: merge_chunk_projfiles
use simple_qsys_env,               only: qsys_env
use simple_qsys_funs,              only: qsys_watcher, qsys_cleanup, qsys_job_finished
use simple_rec_list,               only: rec, project_rec, process_rec, chunk_rec, rec_list, rec_iterator
use simple_sp_project,             only: sp_project
use simple_stack_io,               only: stack_io
use simple_starproject_stream,     only: starproject_stream
use simple_stream_chunk,           only: stream_chunk
use simple_stream_cluster2D_utils, only: cleanup_root_folder, consolidate_sigmas, setup_downscaling, terminate_chunks, terminate_stream2D,&
                                  &test_repick, tidy_2Dstream_iter, update_user_params2D, write_project_stream2D, write_repick_refs
use simple_stream_communicator,    only: stream_http_communicator 
use simple_stream_utils,           only: update_user_params, wait_for_folder, class_rejection, stream_datestr, process_selected_refs,&
                                  &get_latest_optics_map_id
use simple_stream_watcher,         only: stream_watcher

! PUBLICIZE API
public :: class_frcs
public :: cmdline
public :: commander_base
public :: average_sigma2_groups, sigma2_star_from_iter
public :: mic2thumb
public :: guistats
public :: image
public :: simple_nice_communicator
public :: parameters, params_glob
public :: progressfile_init, progressfile_update
public :: merge_chunk_projfiles
public :: qsys_env
public :: qsys_watcher, qsys_cleanup, qsys_job_finished
public :: rec, project_rec, process_rec, chunk_rec, rec_list, rec_iterator
public :: sp_project
public :: stack_io
public :: starproject_stream
public :: stream_chunk
public :: cleanup_root_folder, consolidate_sigmas, setup_downscaling, terminate_chunks, terminate_stream2D,&
         &test_repick, tidy_2Dstream_iter, update_user_params2D, write_project_stream2D, write_repick_refs
public :: stream_http_communicator 
public :: update_user_params, wait_for_folder, class_rejection, stream_datestr, process_selected_refs,&
         &get_latest_optics_map_id
public :: stream_watcher
end module simple_stream_module_api
