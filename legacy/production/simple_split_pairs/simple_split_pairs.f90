!==Program simple_split_pairs
!
! <split_pairs/begin> is a program for splitting calculations between pairs of objects 
! into balanced partitions. <split_pairs/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016
!
program simple_split_pairs
use simple_map_reduce, only: split_pairs_in_parts
use simple_jiffys,     only: simple_end
use simple_params,     only: params
use simple_cmdline,    only: cmdline
use simple_timing
implicit none
type(params)  :: p
type(cmdline) :: cline
if( command_argument_count() < 2 )then
    write(*,'(a)') 'SIMPLE_SPLIT_PAIRS nptcls=<number of particles> npart=<number of partitions>'
    stop
endif
call cline%parse
call cline%check
p = params(cline) ! parameters generated
call cline%checkvar('nptcls', 1)
call cline%checkvar('npart',  2)
call split_pairs_in_parts(p%nptcls, p%npart)
call simple_end('**** SIMPLE_SPLIT_PAIRS NORMAL STOP ****')
end program
