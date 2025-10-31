program simple_test_mini_stream
include 'simple_lib.f08'
use simple_cmdline,               only: cmdline
use simple_parameters,            only: parameters
use simple_commanders_project,    only: commander_new_project
use simple_commanders_validate,   only: commander_mini_stream
use simple_commanders_preprocess, only: commander_preprocess
implicit none
character(len=LONGSTRLEN), allocatable :: dataset_cmds(:)
type(cmdline)               :: cline, cline_dataset, cline_new_project, cline_preprocess, cline_mini_stream
type(parameters)            :: params
type(commander_new_project) :: xnew_project
type(commander_preprocess)  :: xpreprocess
type(commander_mini_stream)  :: xmini_stream
integer                     :: i, ndata_sets
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
! projname=name_system smpd=1.3 cs=2.7 fraca=0.1 total_dose=53 dir_movies=filetab.txt gainref=gainref.mrc nparts=4 nthr=16 moldiam_max=200
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
    call cline_dataset%checkvar('nthr',            9)
    call cline_dataset%checkvar('moldiam_max',    10)
    call cline_dataset%check
    call params%new(cline_dataset)
    ! create project
    call cline_new_project%set('projname',  params%projname)
    call cline_new_project%set('dir',       PATH_HERE//'/'//params%projname)
    call xnew_project%execute(cline_new_project)
    call cline_dataset%kill
    ! preprocess
    call cline_preprocess%set('smpd',       params%smpd)
    call cline_preprocess%set('gainref',    params%gainref)
    call cline_preprocess%set('cs',         params%cs)
    call cline_preprocess%set('kv',         params%kv)
    call cline_preprocess%set('fraca',      params%fraca)
    call cline_preprocess%set('total_dose', params%total_dose)
    call xpreprocess%execute(cline_preprocess)
    call cline_dataset%kill
    ! mini_stream 
    call cline_mini_stream%set('filetab', params%filetab)
    call cline_mini_stream%set('moldiam_max', params%moldiam_max)
    call xmini_stream%execute(cline_mini_stream)
    call cline_dataset%kill
    call cline_new_project%kill
    call cline_preprocess%kill
    call cline_dataset%kill
enddo
! copy movies to local storage?
end program simple_test_mini_stream
