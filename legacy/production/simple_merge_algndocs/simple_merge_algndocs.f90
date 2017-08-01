!==Program simple_merge_algndocs
!
! <merge_algndocs/begin> is a program for merging alignment documents produced by PRIME2D/3D 
! when run in distributed mode using \texttt{distr\_simple.pl}. <merge_algndocs/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_merge_algndocs
use simple_defs     ! singleton
use simple_jiffys   ! singleton
use simple_cmdline, only: cmdline
use simple_oris,    only: oris
use simple_params,  only: params
use simple_timing
implicit none
type(params)          :: p
type(cmdline)         :: cline
type(oris)            :: o, o_read
integer               :: i, j, istart, istop, ptcls_per_part, nentries, leftover, cnt, nentries_all, numlen
character(len=STDLEN) :: fname
logical               :: here, useoritab
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_MERGE_ALGNDOCS fbody=<file_body of algndocs> nptcls=<nr of particle images>'
    write(*,'(a)') ' ndocs=<nr of docs> outfile=<merged alignment doc> [oritab=<previous oritab>]'
    stop
endif
call cline%parse
call cline%checkvar('fbody',  1)
call cline%checkvar('nptcls', 2)
call cline%checkvar('ndocs',  3)
call cline%check
call cline%set('prg', 'merge_algndocs')
p = params(cline) ! parameters generated
useoritab = .false.
if( cline%defined('oritab') )then
    if( file_exists(p%oritab) ) useoritab = .true.
endif
if( useoritab )then
    if( nlines(p%oritab) /= p%nptcls )then
        stop 'the inputted nptcls is not consistent with the nptcls in oritab!'
    endif
    ! create object for orientations
    o = oris(p%nptcls)
    ! read previous orientations
    call o%read(p%oritab)
endif
ptcls_per_part = p%nptcls/p%ndocs
leftover       = p%nptcls-ptcls_per_part*p%ndocs
istop          = 0
numlen         = len(int2str(p%ndocs))
do i=1,p%ndocs
    fname  = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//'.txt'
    if( i == p%ndocs )then
        istart = istop+1
        istop  = p%nptcls
    else
        if( leftover == 0 )then
            istart = istop+1;
            istop  = istart+ptcls_per_part-1;
        else
            istop  = i*(ptcls_per_part+1)
            istart = istop-(ptcls_per_part+1)+1
            leftover = leftover-1
        endif
    endif
    ! calculate the number of all entries
    nentries_all = istop-istart+1 
    ! calculate the actual number of entries
    inquire(FILE=fname, EXIST=here)
    if( here )then
        nentries = nlines(fname)
    else
        nentries = 0
    endif
    ! check if oritab is there to fill-in blanks
    if( nentries < nentries_all )then
        if( .not. useoritab )then
            stop 'need previous oritab to fill-in blanks; simple_merge_algndocs'
        endif
    endif
    ! print partition info
    write(*,'(a,1x,i3,1x,a,1x,i6,1x,i6)') 'partition:', i, 'from/to:', istart, istop
    if( nentries > 0 )then
        o_read = oris(nentries)
        call o_read%read(fname)
    endif
    ! read
    if( useoritab )then ! read and fill-in from oritab
        cnt = 0
        do j=istart,istop
            cnt = cnt+1
            if( cnt <= nentries )then
                call o%set_ori(j,o_read%get_ori(cnt))
            else
                exit
            endif
        end do
    else                                ! just merge (all ptcls is there)
        call o%merge(o_read)
    endif
end do
call o%write(p%outfile)
call simple_end('**** SIMPLE_MERGE_ALGNDOCS NORMAL STOP ****')
end program 
