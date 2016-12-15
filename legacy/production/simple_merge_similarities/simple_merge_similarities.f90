!==Program simple_merge_similarities
!
! <simple\_merge_similarities/begin> is a program for splitting calculations between pairs of objects 
! into balanced partitions. <simple\_merge_similarities/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016
!
program simple_merge_similarities
use simple_map_reduce, only: merge_similarities_from_parts
use simple_jiffys,     only: simple_end, get_fileunit
use simple_params,     only: params
use simple_cmdline,    only: cmdline
use simple_timing
implicit none
type(params)      :: p
type(cmdline)     :: cline
real, allocatable :: simmat(:,:)
integer           :: filnum, io_stat
if( command_argument_count() < 2 )then
    write(*,'(a)') 'SIMPLE_MERGE_SIMILARITIES nptcls=<nr particles> npart=<nr partitions>'
    stop
endif
call cline%parse
call cline%check
p = params(cline) ! parameters generated
call cline%checkvar('nptcls', 1)
call cline%checkvar('npart',  2)
simmat = merge_similarities_from_parts(p%nptcls, p%npart)
filnum = get_fileunit()
open(unit=filnum, status='REPLACE', action='WRITE', file='smat.bin', access='STREAM')
write(unit=filnum,pos=1,iostat=io_stat) simmat
if( io_stat .ne. 0 )then
    write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to smat.bin'
    stop 'I/O error; simple_merge_similarities'
endif
close(filnum)
call simple_end('**** SIMPLE_MERGE_SIMILARITIES NORMAL STOP ****')
end program
