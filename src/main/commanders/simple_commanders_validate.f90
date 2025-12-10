
! concrete commander: cluster2D_stream for streaming 2D alignment and clustering of single-particle images
module simple_commanders_validate
include 'simple_lib.f08'
use simple_sp_project,     only: sp_project
use simple_parameters,     only: parameters
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_commanders_abinitio2D
use simple_commanders_preprocess
use simple_commanders_project
use simple_micproc  
use simple_mini_stream_utils
use simple_picksegdiam
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_mini_stream
  contains
    procedure :: execute      => exec_mini_stream
end type commander_mini_stream

type, extends(commander_base) :: commander_check_refpick
  contains
    procedure :: execute      => exec_check_refpick
end type commander_check_refpick

contains

    subroutine exec_mini_stream( self, cline )
        class(commander_mini_stream), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        character(len=*), parameter        :: PROJFILE_MINI_STREAM = 'mini_stream.simple'
        integer,          parameter        :: NCLS_MIN = 10, NCLS_MAX = 100
        real,             parameter        :: LPSTOP = 8.
        type(string),     allocatable      :: micnames(:)
        type(string)                       :: output_dir
        type(parameters)                   :: params
        type(sp_project)                   :: spproj
        type(cmdline)                      :: cline_new_proj, cline_import_movies, cline_ctf_estimate
        type(cmdline)                      :: cline_extract, cline_abinitio2D, cline_shape_rank
        type(commander_new_project)        :: xnew_project
        type(commander_import_movies)      :: ximport_movies
        type(commander_ctf_estimate_distr) :: xctf_estimate
        type(commander_extract_distr)      :: xextract
        type(commander_abinitio2D)         :: xabinitio2D
        type(commander_shape_rank_cavgs)   :: xshape_rank
        integer :: ncls, nmics, nptcls, box_in_pix
        real    :: mskdiam_estimate
        if( .not. cline%defined('mkdir')            ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('kv')               ) call cline%set('kv',              300.)
        if( .not. cline%defined('cs')               ) call cline%set('cs',               2.7)
        if( .not. cline%defined('fraca')            ) call cline%set('fraca',            0.1)
        if( .not. cline%defined('pspecsz')          ) call cline%set('pspecsz',          512)
        if( .not. cline%defined('hp')               ) call cline%set('hp',               30.)
        if( .not. cline%defined('lp')               ) call cline%set('lp',                5.)
        if( .not. cline%defined('dfmin')            ) call cline%set('dfmin',  DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')            ) call cline%set('dfmax',  DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')         ) call cline%set('ctfpatch',        'no')
        if( .not. cline%defined('nptcls_per_class') ) call cline%set('nptcls_per_cls',   200)
        if( .not. cline%defined('pick_roi')         ) call cline%set('pick_roi',       'yes')
        call params%new(cline)
        call read_filetable(params%filetab, micnames)
        nmics = size(micnames)
        ! output directory
        output_dir = PATH_HERE
        ! project creation
        call cline_new_proj%set('dir',                      PATH_HERE)
        call cline_new_proj%set('projname',     string('mini_stream'))
        call xnew_project%execute_safe(cline_new_proj)
        ! movie import
        call cline_import_movies%set('prg',           'import_movies')
        call cline_import_movies%set('mkdir',                    'no')
        call cline_import_movies%set('cs',                  params%cs)
        call cline_import_movies%set('fraca',            params%fraca)
        call cline_import_movies%set('kv',                  params%kv)
        call cline_import_movies%set('smpd',              params%smpd)
        call cline_import_movies%set('filetab',        params%filetab)
        call cline_import_movies%set('ctf',                     'yes')
        call cline_import_movies%set('projfile', PROJFILE_MINI_STREAM)
        call ximport_movies%execute_safe(cline_import_movies)
        ! CTF estimation
        call cline_ctf_estimate%set('prg',            'ctf_estimate')
        call cline_ctf_estimate%set('mkdir',                    'no')
        call cline_ctf_estimate%set('ctfpatch',      params%ctfpatch)
        call cline_ctf_estimate%set('dfmax',            params%dfmax)
        call cline_ctf_estimate%set('dfmin',            params%dfmin)
        call cline_ctf_estimate%set('hp',                  params%hp)
        call cline_ctf_estimate%set('lp',                  params%lp)
        call cline_ctf_estimate%set('nparts',                      1)
        call cline_ctf_estimate%set('nthr',              params%nthr)
        call cline_ctf_estimate%set('projfile', PROJFILE_MINI_STREAM)
        call xctf_estimate%execute_safe(cline_ctf_estimate)
        ! this is the actual test
        call spproj%read(string(PROJFILE_MINI_STREAM))
        call segdiampick_mics(spproj, params%pcontrast, nmics, params%moldiam_max, box_in_pix, mskdiam_estimate)
        ! segdiampick_mics updates the project file on disk
        ! extract
        call cline_extract%set('prg',                      'extract')
        call cline_extract%set('mkdir',                         'no')
        call cline_extract%set('nparts',                           1)
        call cline_extract%set('nthr',                   params%nthr)
        call cline_extract%set('projfile',      PROJFILE_MINI_STREAM)
        call xextract%execute_safe(cline_extract)
        ! 2D analysis
        call spproj%read(string(PROJFILE_MINI_STREAM))
        nptcls = spproj%os_ptcl2D%get_noris()
        ncls   = min(NCLS_MAX,max(NCLS_MIN,nptcls/params%nptcls_per_cls))
        call cline_abinitio2D%set('prg',                'abinitio2D')
        call cline_abinitio2D%set('mkdir',                      'no')
        call cline_abinitio2D%set('ncls',                       ncls)
        call cline_abinitio2D%set('sigma_est',              'global')
        call cline_abinitio2D%set('center',                    'yes')
        call cline_abinitio2D%set('autoscale',                 'yes')
        call cline_abinitio2D%set('lpstop',                   LPSTOP)
        call cline_abinitio2D%set('mskdiam',        mskdiam_estimate)
        call cline_abinitio2D%set('nthr',                params%nthr)
        call cline_abinitio2D%set('projfile',   PROJFILE_MINI_STREAM)
        call xabinitio2D%execute_safe(cline_abinitio2D)
        ! shape rank cavgs
        call cline_shape_rank%set('nthr',                params%nthr)
        call cline_shape_rank%set('projfile',   PROJFILE_MINI_STREAM)
        call xshape_rank%execute_safe(cline_shape_rank)
        ! destruct
        call spproj%kill
    end subroutine exec_mini_stream

    subroutine exec_check_refpick( self, cline )
        class(commander_check_refpick), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        character(len=*), parameter        :: DIR_THUMBS             = 'thumbnails/'
        character(len=*), parameter        :: PROJFILE_CHECK_REFPICK = 'check_refpick.simple'
        integer,          parameter        :: NCLS_MIN = 10, NCLS_MAX = 100
        real,             parameter        :: LPSTOP = 8.
        type(string),     allocatable      :: micnames(:), fnames(:)
        type(string)                       :: output_dir
        type(parameters)                   :: params
        type(sp_project)                   :: spproj
        type(cmdline)                      :: cline_new_proj, cline_import_movies, cline_ctf_estimate
        type(cmdline)                      :: cline_make_pickrefs, cline_pick_extract, cline_abinitio2D, cline_shape_rank
        type(commander_new_project)        :: xnew_project
        type(commander_make_pickrefs)      :: xmake_pickrefs
        type(commander_import_movies)      :: ximport_movies
        type(commander_ctf_estimate_distr) :: xctf_estimate
        type(commander_pick_extract)       :: xpickextract
        type(commander_abinitio2D)         :: xabinitio2D
        type(commander_shape_rank_cavgs)   :: xshape_rank
        integer :: ncls, nmics, nptcls 
        real    :: mskdiam_estimate
        if( .not. cline%defined('mkdir')            ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('kv')               ) call cline%set('kv',              300.)
        if( .not. cline%defined('cs')               ) call cline%set('cs',               2.7)
        if( .not. cline%defined('fraca')            ) call cline%set('fraca',            0.1)
        if( .not. cline%defined('pspecsz')          ) call cline%set('pspecsz',          512)
        if( .not. cline%defined('hp')               ) call cline%set('hp',               30.)
        if( .not. cline%defined('lp')               ) call cline%set('lp',                5.)
        if( .not. cline%defined('dfmin')            ) call cline%set('dfmin',  DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')            ) call cline%set('dfmax',  DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')         ) call cline%set('ctfpatch',        'no')
        if( .not. cline%defined('nptcls_per_class') ) call cline%set('nptcls_per_cls',   200)
        if( .not. cline%defined('pick_roi')         ) call cline%set('pick_roi',       'yes')
        call params%new(cline)
        call read_filetable(params%filetab, micnames)
        nmics = size(micnames)
        ! output directory
        output_dir = PATH_HERE
        ! project creation
        call cline_new_proj%set('dir',                         PATH_HERE)
        call cline_new_proj%set('projname',      string('check_refpick'))
        call xnew_project%execute_safe(cline_new_proj)
        ! movie import
        call cline_import_movies%set('prg',              'import_movies')
        call cline_import_movies%set('mkdir',                       'no')
        call cline_import_movies%set('cs',                     params%cs)
        call cline_import_movies%set('fraca',               params%fraca)
        call cline_import_movies%set('kv',                     params%kv)
        call cline_import_movies%set('smpd',                 params%smpd)
        call cline_import_movies%set('filetab',           params%filetab)
        call cline_import_movies%set('ctf',                        'yes')
        call cline_import_movies%set('projfile',  PROJFILE_CHECK_REFPICK)
        call ximport_movies%execute_safe(cline_import_movies)
        ! CTF estimation
        call cline_ctf_estimate%set('prg',                'ctf_estimate')
        call cline_ctf_estimate%set('mkdir',                        'no')
        call cline_ctf_estimate%set('ctfpatch',          params%ctfpatch)
        call cline_ctf_estimate%set('dfmax',                params%dfmax)
        call cline_ctf_estimate%set('dfmin',                params%dfmin)
        call cline_ctf_estimate%set('hp',                      params%hp)
        call cline_ctf_estimate%set('lp',                      params%lp)
        call cline_ctf_estimate%set('nparts',                          1)
        call cline_ctf_estimate%set('nthr',                  params%nthr)
        call cline_ctf_estimate%set('projfile',   PROJFILE_CHECK_REFPICK)
        call xctf_estimate%execute_safe(cline_ctf_estimate)
        ! make pickrefs
        call cline_make_pickrefs%set('mkdir',                       'no')
        call cline_make_pickrefs%set('pickrefs',         params%pickrefs)
        call cline_make_pickrefs%set('smpd',                 params%smpd)
        call cline_make_pickrefs%set('nthr',                 params%nthr)
        call xmake_pickrefs%execute_safe(cline_make_pickrefs)
        mskdiam_estimate = cline_make_pickrefs%get_rarg('mskdiam')
        ! pick and extract
        cline_pick_extract = cline
        call cline_pick_extract%set('mkdir',                        'no')
        call cline_pick_extract%set('pickrefs', PICKREFS_FBODY//params%ext%to_char())               
        call cline_pick_extract%set('stream',                      'yes')
        call cline_pick_extract%set('extract',                     'yes')
        call cline_pick_extract%set('dir',                    output_dir)
        call cline_pick_extract%set('fromp',                           1)
        call cline_pick_extract%set('top',                         nmics)
        call cline_pick_extract%set('projfile',   PROJFILE_CHECK_REFPICK)
        call xpickextract%execute_safe(cline_pick_extract)
        ! 2D analysis
        call spproj%read(string(PROJFILE_CHECK_REFPICK))
        nptcls = spproj%os_ptcl2D%get_noris()
        ncls   = min(NCLS_MAX,max(NCLS_MIN,nptcls/params%nptcls_per_cls))
        call cline_abinitio2D%set('prg',                    'abinitio2D')
        call cline_abinitio2D%set('mkdir',                          'no')
        call cline_abinitio2D%set('ncls',                           ncls)
        call cline_abinitio2D%set('sigma_est',                  'global')
        call cline_abinitio2D%set('center',                        'yes')
        call cline_abinitio2D%set('autoscale',                     'yes')
        call cline_abinitio2D%set('lpstop',                       LPSTOP)
        call cline_abinitio2D%set('mskdiam',            mskdiam_estimate)
        call cline_abinitio2D%set('nthr',                    params%nthr)
        call cline_abinitio2D%set('projfile',     PROJFILE_CHECK_REFPICK)
        call xabinitio2D%execute_safe(cline_abinitio2D)
        ! shape rank cavgs
        call cline_shape_rank%set('nthr',                    params%nthr)
        call cline_shape_rank%set('projfile',     PROJFILE_CHECK_REFPICK)
        call xshape_rank%execute_safe(cline_shape_rank)
        ! organize output in folders
        call simple_list_files('*jpg', fnames)
        call simple_mkdir(DIR_THUMBS)
        call move_files2dir(string(DIR_THUMBS), fnames)
        ! destruct
        call spproj%kill
        call fnames%kill
    end subroutine exec_check_refpick

end module simple_commanders_validate