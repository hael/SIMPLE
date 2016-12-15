program simple_euler_consensus
use simple_cmdline           ! singleton
use simple_defs              ! singleton
use simple_orialgn_ensemble, only: orialgn_ensemble
use simple_params,           only: params
use simple_oris,             only: oris
use simple_math,             only: sortmeans
use simple_jiffys,           only: simple_end, alloc_err, fopen_err, get_fileunit, nlines
implicit none
type(params)           :: p
type(orialgn_ensemble) :: dockem
type(oris)             :: consensus
integer                :: i, nd, ns, fnr, file_stat, alloc_stat, good, ngood, cnt
character(len=STDLEN)  :: line
real                   :: means(2)
real, allocatable      :: scores(:)
integer, allocatable   :: labels(:)
if( command_argument_count() < 2 )then
    write(*,'(a)', advance='no') 'simple_euler_consensus doclist=<list of alignment docs>'
    write(*,'(a)') ' outfile=<output consensus alignment doc> [scorelist=<list of score values>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('doclist', 1)
call cmdcheckvar('outfile', 2)
call cmdcheck
p = params()
nd = nlines(p%doclist)
if( defined_cmd_arg('scorelist') )then
    ns = nlines(p%scorelist)
    if( nd /= ns ) stop 'nr scores .ne. nr docs'
    allocate( scores(ns), labels(ns), stat=alloc_stat )
    call alloc_err("In: simple_euler_consensus", alloc_stat)
    fnr = get_fileunit( )
    open(unit=fnr, FILE=p%scorelist, STATUS='OLD', action='READ', iostat=file_stat)
    call fopen_err( 'In: simple_euler_consensus', file_stat )
    do i=1,ns
        read(fnr,*) scores(i)
    end do
    close(fnr)
    call sortmeans(scores, means, labels)
    good = 0
    if( means(1) > means(2) )then
        good = 1
    else
        good = 2
    endif
    write(*,'(a)') '>>> SOLUTIONS SELECTED FOR CONSENSUS CALCULATION'
    ngood = 0
    do i=1,ns
        if( labels(i) == good )then
            write(*,'(1X,A,I4,1X,A,f8.4)') 'SELECTED:', i, 'SCORE:', scores(i)
            ngood = ngood+1
        endif
    end do
else
    ngood = nd
endif
call dockem%new(ngood)
fnr = get_fileunit( )
open(unit=fnr, FILE=p%doclist, STATUS='OLD', action='READ', iostat=file_stat)
call fopen_err( 'In: simple_euler_consensus', file_stat )
if( defined_cmd_arg('scorelist') )then
    cnt = 0
    do i=1,nd
        read(fnr,*) line
        if( labels(i) == good )then
            cnt = cnt+1
            call dockem%read(cnt, line)
        endif
    end do
else
    do i=1,nd
        read(fnr,*) line
        call dockem%read(i, line)
    end do
endif
close(fnr)
consensus = dockem%l1median()
call consensus%write(p%outfile)
call simple_end('**** SIMPLE_EULER_CONSENSUS NORMAL STOP ****')
end program