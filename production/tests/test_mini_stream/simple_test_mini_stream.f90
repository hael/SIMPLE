program simple_test_mini_stream
include 'simple_lib.f08'
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_commanders_project,  only: commander_new_project, commander_import_movies
use simple_commanders_validate, only: commander_mini_stream
use simple_commanders_project,  only: commander_selection
use simple_stream_watcher,        only: stream_watcher
use simple_commanders_preprocess 
use simple_sp_project
implicit none
real,         parameter          :: CTFRES_THRES = 8.0, ICE_THRES = 1.0, OVERSHOOT = 1.2
type(string), allocatable        :: dataset_cmds(:)
type(string), allocatable        :: micstab(:), filetab(:), movfnames(:)
type(string)                     :: output_dir, imgkind
integer,      allocatable        :: orimap(:)
type(cmdline)                    :: cline, cline_dataset, cline_new_project, cline_import_movies, cline_preprocess
type(cmdline)                    :: cline_select, cline_mini_stream
type(parameters)                 :: params
type(commander_new_project)      :: xnew_project
type(commander_preprocess_distr) :: xpreprocess
type(commander_import_movies)    :: ximport_movies
type(commander_mini_stream)      :: xmini_stream
type(commander_selection)        :: xsel
type(sp_project)                 :: spproj
type(stream_watcher)               :: movie_buff
integer                          :: i, ndata_sets, n_nonzero, nmovf
type(string)                     :: abspath, projfile
character(len=*), parameter      :: filetab_file='filetab.txt'
! Parsing
if( command_argument_count() < 1 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_mini_stream fname=filetab.txt'
    call exit(-1)
else 
    call cline%parse_oldschool
endif
call cline%checkvar('fname',        1)
call cline%check()
call params%new(cline)
call read_filetable(params%fname, dataset_cmds)
ndata_sets=size(dataset_cmds)
! projname=system_name smpd=1.3 cs=2.7 kv=300 fraca=0.1 total_dose=53 dir_movies=/usr/local/data/movies gainref=gainref.mrc nparts=4 nthr=16 moldiam_max=200 nram=100
call simple_getcwd(abspath)
output_dir=abspath
do i = 1, ndata_sets
    call cline_dataset%read(dataset_cmds(i)%to_char())
    call cline_dataset%checkvar('projname',        1)
    call cline_dataset%checkvar('smpd',            2)
    call cline_dataset%checkvar('cs',              3)
    call cline_dataset%checkvar('kv',              4)
    call cline_dataset%checkvar('fraca',           5)
    call cline_dataset%checkvar('total_dose',      6)
    call cline_dataset%checkvar('dir_movies',      7)
    call cline_dataset%checkvar('gainref',         8)
    call cline_dataset%checkvar('nparts',          9)
    call cline_dataset%checkvar('nthr',           10)
    call cline_dataset%checkvar('moldiam_max',    11)
    call cline_dataset%checkvar('nran',           12)
    call cline_dataset%check()
    call params%new(cline_dataset)
    call cline_dataset%kill()
    ! project creation
    call cline_new_project%set('projname',  params%projname)
    call xnew_project%execute_safe(cline_new_project)
    call cline_new_project%kill()
    projfile = params%projname//'.simple'
    ! create filetab with a subset of overshoot randomly selected movies
    movie_buff = stream_watcher(1,params%dir_movies)
    call movie_buff%watch(nmovf, movfnames)
    call movie_buff%kill
    filetab = sample_filetab(movfnames, ceiling(real(params%nran)*OVERSHOOT))
    call write_filetable(string(filetab_file), filetab)
    ! movie import
    call cline_import_movies%set('prg',                'import_movies')
    call cline_import_movies%set('mkdir',                        'yes')
    call cline_import_movies%set('cs',                       params%cs)
    call cline_import_movies%set('fraca',                 params%fraca)
    call cline_import_movies%set('kv',                       params%kv)
    call cline_import_movies%set('smpd',                   params%smpd)
    call cline_import_movies%set('filetab',               filetab_file)
    call cline_import_movies%set('ctf',                          'yes')
    call ximport_movies%execute_safe(cline_import_movies)
    call cline_import_movies%kill()
    ! check either movies or micrographs
    call spproj%read(projfile)
    imgkind = spproj%get_mic_kind(1)
    if( imgkind.eq.'intg' )then
        ! nothing to do
    else
        ! preprocess
        call simple_chdir(output_dir//'/'//params%projname%to_char())
        call cline_preprocess%set('prg',                  'preprocess')
        call cline_preprocess%set('mkdir',                       'yes')
        call cline_preprocess%set('gainref',            params%gainref)
        call cline_preprocess%set('total_dose',      params%total_dose)
        call cline_preprocess%set('dfmin',               DFMIN_DEFAULT)
        call cline_preprocess%set('dfmax',               DFMAX_DEFAULT)
        call cline_preprocess%set('hp',                            30.)
        call cline_preprocess%set('lp',                             2.)
        call cline_preprocess%set('mcpatch',                      'no')
        call cline_preprocess%set('nparts',              params%nparts)
        call cline_preprocess%set('nthr',                  params%nthr)
        call cline_preprocess%check()
        call xpreprocess%execute_safe(cline_preprocess)
        call cline_preprocess%kill()
    endif
    ! reject based on CTF resolution and ice score
    call simple_chdir(output_dir//'/'//params%projname%to_char())
    call cline_select%delete('nran')
    call cline_select%set('prg',                           'selection')
    call cline_select%set('mkdir',                               'yes')
    call cline_select%set('oritype',                             'mic')
    call cline_select%set('ctfresthreshold',              CTFRES_THRES)
    call cline_select%set('icefracthreshold',                ICE_THRES)
    call xsel%execute_safe(cline_select)
    call cline_select%kill()
    ! state=0/1 should now be in project file on disk
    ! re-run for random selection
    call spproj%read(projfile)
    n_nonzero = spproj%get_n_insegment_state('mic', 1)
    if( n_nonzero > params%nran )then
        ! make random selection
        call simple_chdir(output_dir//'/'//params%projname%to_char())
        call cline_select%delete('ctfresthreshold')
        call cline_select%delete('icefracthreshold')
        call cline_select%set('prg',                       'selection')
        call cline_select%set('mkdir',                           'yes') 
        call cline_select%set('oritype',                         'mic')
        call cline_select%set('nran',                      params%nran)
        call xsel%execute_safe(cline_select)
        call cline_select%kill()
    endif
    call spproj%read(projfile)
    call spproj%get_mics_table(micstab, orimap)
    call simple_chdir(output_dir//'/'//params%projname%to_char())
    call write_filetable(string('intgs.txt'),micstab)
    ! mini stream 
    call cline_mini_stream%set('prg',                    'mini_stream')
    call cline_mini_stream%set('mkdir',                          'yes')
    call cline_mini_stream%set('filetab',                  'intgs.txt')
    call cline_mini_stream%set('smpd',                     params%smpd)
    call cline_mini_stream%set('fraca',                   params%fraca)
    call cline_mini_stream%set('kv',                         params%kv)
    call cline_mini_stream%set('cs',                         params%cs)
    call cline_mini_stream%set('moldiam_max',       params%moldiam_max)
    call cline_mini_stream%set('nparts',                 params%nparts)
    call cline_mini_stream%set('nthr',                     params%nthr)
    call xmini_stream%execute_safe(cline_mini_stream)
    call cline_dataset%kill()
    call cline_new_project%kill()
    call cline_import_movies%kill()
    call cline_preprocess%kill()
    call cline_mini_stream%kill()
    call cline_dataset%kill()
    call simple_chdir(output_dir)
enddo
end program simple_test_mini_stream
