!==Program simple_tseries_split
!
! <simple\_tseries\_split/begin> is a program for splitting time series.
! <simple\_tseries\_split/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
program simple_tseries_split
use simple_jiffys,  ! singleton
use simple_cmdline, only: cmdline
use simple_params,  only: params
use simple_build,   only: build
use simple_oris,    only: oris
use simple_ori,     only: ori
implicit none
type(params)                  :: p
type(build)                   :: b
type(cmdline)                 :: cline
integer                       :: iptcl, istart, istop, cnt, chunkcnt, numlen
type(oris)                    :: oset
type(ori)                     :: o
character(len=:), allocatable :: fname, oname
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_TSERIES_SPLIT stk=<input particle stack> oritab=oritab=<alignment doc>'
    write(*,'(a)') ' chunksz=<size of contigous series> jumpsz=<time step size>'
    stop
endif
call cline%parse
call cline%checkvar('stk', 1)
call cline%check
p = params(cline)                   ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
istart   = 1
istop    = istart+p%chunksz-1
chunkcnt = 0
numlen   = len(int2str(p%nptcls/p%chunksz+1))
do while( istop <= p%nptcls )
    chunkcnt = chunkcnt+1
    allocate(fname, source='stack_chunk'//int2str_pad(chunkcnt,numlen)//p%ext)
    allocate(oname, source='oris_chunk'//int2str_pad(chunkcnt,numlen)//'.txt')
    cnt = 0
    oset = oris(istop-istart+1)
    do iptcl=istart,istop
        cnt = cnt+1
        call b%img%read(p%stk, iptcl)
        call b%img%write(fname, cnt)
        o = b%a%get_ori(iptcl)
        call oset%set_ori(cnt, o)
    end do
    if( cnt == oset%get_noris() )then
        ! all ok
    else
        stop 'wrong number of oris allocated'
    endif
    call oset%write(oname)
    call oset%kill
    istart = istart+p%jumpsz
    istop  = istop+p%jumpsz
    deallocate(fname, oname)
end do
call simple_end('**** SIMPLE_TSERIES_SPLIT NORMAL STOP ****')
end program simple_tseries_split