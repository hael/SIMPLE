!@descr: Facade pattern API to avoid circular dependencies involving polarft_calc
module simple_pftc_srch_api
use simple_core_module_api
use simple_class_frcs,      only: class_frcs
use simple_cmdline,         only: cmdline
use simple_ctf,             only: ctf
use simple_image,           only: image
use simple_imgarr_utils,    only: read_cavgs_into_imgarr,  read_stk_into_imgarr, dealloc_imgarr, alloc_imgarr, pack_imgarr, write_imgarr
use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps
use simple_parameters,      only: parameters
use simple_polarft_calc,    only: polarft_calc, pftc_glob, polaft_dims_from_file_header
use simple_sigma2_binfile,  only: sigma2_binfile
use simple_sp_project,      only: sp_project
use simple_stack_io,        only: stack_io
use simple_starproject,     only: starproject
end module simple_pftc_srch_api
