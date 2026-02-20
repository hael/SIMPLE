!@descr: API for the commanders
module simple_commanders_api
use simple_core_module_api
use simple_default_clines
use simple_nice
use simple_qsys_funs
use simple_binoris_io,       only: binread_nlines, binread_oritab, binwrite_oritab, binread_ctfparams_state_eo
use simple_builder,          only: builder
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_euclid_sigma2,    only: write_groups_starfile, euclid_sigma2, sigma2_star_from_iter, eucl_sigma2_glob
use simple_exec_helpers,     only: set_shmem_flag, set_master_num_threads
use simple_image,            only: image
use simple_image_bin,        only: image_bin
use simple_image_msk,        only: image_msk, automask2D
use simple_parameters,       only: parameters, params_glob
use simple_qsys_env,         only: qsys_env
use simple_sp_project,       only: sp_project, oritype2segment
use simple_stack_io,         only: stack_io
use simple_stream_utils,     only: terminate_stream
end module simple_commanders_api
