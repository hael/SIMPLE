program simple_test_mini_stream
include 'simple_lib.f08'
use simple_cmdline,               only: cmdline
use simple_parameters,            only: parameters
use simple_commanders_project,    only: commander_new_project, commander_import_movies
use simple_commanders_validate,   only: commander_mini_stream
use simple_commanders_preprocess
use simple_sp_project
implicit none
character(len=LONGSTRLEN), allocatable :: dataset_cmds(:)
type(cmdline)                    :: cline, cline_dataset, cline_new_project, cline_import_movies, cline_preprocess, cline_mini_stream
type(parameters)                 :: params
type(commander_new_project)      :: xnew_project
type(commander_preprocess_distr) :: xpreprocess
type(commander_import_movies)    :: ximport_movies
type(commander_mini_stream)      :: xmini_stream
type(sp_project)                 :: spproj
integer,          allocatable :: orimap(:)
integer                       :: i, ndata_sets, status
character(len=LONGSTRLEN)     :: abspath, projfile
character(len=LONGSTRLEN), allocatable :: micstab(:)
character(len=:), allocatable :: output_dir
! Parsing
if( command_argument_count() < 1 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_mini_stream fname=filetab.txt'
    call exit(-1)
else 
    call cline%parse_oldschool
endif
call cline%checkvar('fname',        1)
call cline%check()
!call cline%printline
call params%new(cline)
call read_filetable(params%fname, dataset_cmds)
ndata_sets=size(dataset_cmds)
! projname=name_system smpd=1.3 cs=2.7 fraca=0.1 total_dose=53 dir_movies=/usr/local/data/movies gainref=gainref.mrc nparts=4 nthr=16 moldiam_max=200
call getcwd(abspath)
output_dir=trim(adjustl(abspath))
do i = 1, ndata_sets
    call cline_dataset%read(dataset_cmds(i))
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
    call cline_dataset%check()
    call params%new(cline_dataset)
    call cline_dataset%kill()
    ! project creation
    call cline_new_project%set('projname',  params%projname)
    call xnew_project%execute_safe(cline_new_project)
    call cline_new_project%kill()
    ! movie import
    call cline_import_movies%set('prg',          'import_movies')
    call cline_import_movies%set('mkdir',                  'yes')
    call cline_import_movies%set('cs',                 params%cs)
    call cline_import_movies%set('fraca',           params%fraca)
    call cline_import_movies%set('kv',                 params%kv)
    call cline_import_movies%set('smpd',             params%smpd)
    call cline_import_movies%set('dir_movies', params%dir_movies)
    call cline_import_movies%set('ctf',                    'yes')
    call ximport_movies%execute_safe(cline_import_movies)
    call cline_import_movies%kill()
    ! preprocess
    call simple_chdir(trim(trim(adjustl(output_dir))//'/'//params%projname))
    call cline_preprocess%set('prg',                'preprocess')
    call cline_preprocess%set('mkdir',                     'yes')
    call cline_preprocess%set('gainref',          params%gainref)
    call cline_preprocess%set('total_dose',    params%total_dose)
    call cline_preprocess%set('dfmin',             DFMIN_DEFAULT)
    call cline_preprocess%set('dfmax',             DFMAX_DEFAULT)
    call cline_preprocess%set('hp',                          30.)
    call cline_preprocess%set('lp',                           2.)
    call cline_preprocess%set('mcpatch',                    'no')
    call cline_preprocess%set('nparts',            params%nparts)
    call cline_preprocess%set('nthr',                params%nthr)
    call cline_preprocess%check()
    call xpreprocess%execute_safe(cline_preprocess)
    call cline_preprocess%kill()
    projfile = trim(params%projname)//'.simple'
    call spproj%read(projfile)
    call spproj%get_mics_table(micstab, orimap)
    call simple_chdir(trim(trim(adjustl(output_dir))//'/'//params%projname))
    call write_filetable('intgs.txt',micstab)
    ! mini stream 
    call cline_mini_stream%set('prg',               'mini_stream')
    call cline_mini_stream%set('mkdir',                     'yes')
    call cline_mini_stream%set('filetab',             'intgs.txt')
    call cline_mini_stream%set('smpd',                params%smpd)
    call cline_mini_stream%set('fraca',              params%fraca)
    call cline_mini_stream%set('kv',                    params%kv)
    call cline_mini_stream%set('cs',                    params%cs)
    call cline_mini_stream%set('moldiam_max',  params%moldiam_max)
    call cline_mini_stream%set('nparts',            params%nparts)
    call cline_mini_stream%set('nthr',                params%nthr)
    call xmini_stream%execute_safe(cline_mini_stream)
    call cline_dataset%kill()
    call cline_new_project%kill()
    call cline_import_movies%kill()
    call cline_preprocess%kill()
    call cline_mini_stream%kill()
    call cline_dataset%kill()
    call simple_chdir(trim(adjustl(output_dir)))
enddo
end program simple_test_mini_stream
