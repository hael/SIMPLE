!==Program simple_bootstrap
!
! <bootstrap/begin> is a program for bootstrap sampling. <bootstrap/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_bootstrap
use simple_cmdline ! singleton
use simple_stat,   only: moment
use simple_rnd,    only: seed_rnd, irnd_uni
use simple_params, only: params
use simple_jiffys, only: simple_end, txtfile2rarr, alloc_err
implicit none
type(params)      :: p
real, allocatable :: arr(:), bootarr(:)
logical           :: err
real              :: ave, sdev, var
integer           :: sz, alloc_stat, i
if( command_argument_count() < 1 )then
    write(*,'(a)',advance='no') 'SIMPLE_BOOTSTRAP infile=<text file with numbers>'
    write(*,'(a)') ' [nboot=<nr of bootstrap samples{5000}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('infile', 1)
call cmdcheck
p = params()  ! parameters generated
call seed_rnd ! random number generator seeded
arr = txtfile2rarr(p%infile)
call moment(arr, ave, sdev, var, err)
write(*,'(A,1X,F7.4)') '>>> AVERAGE:', ave
write(*,'(A,1X,F7.4)') '>>> STANDARD DEVIATION:', sdev
write(*,'(A,1X,F7.4)') '>>> VARIANCE:', var
allocate( bootarr(p%nboot), stat=alloc_stat )
call alloc_err("In: simple_bootstrap", alloc_stat)
sz = size(arr)
do i=1,p%nboot
    bootarr(i) = arr(irnd_uni(sz))
end do
call moment(bootarr, ave, sdev, var, err)
write(*,'(A,1X,F7.4)') '>>> BOOTSTRAP AVERAGE:', ave
write(*,'(A,1X,F7.4)') '>>> BOOTSTRAP STANDARD DEVIATION:', sdev
write(*,'(A,1X,F7.4)') '>>> BOOTSTRAP VARIANCE:', var
call simple_end('**** SIMPLE_BOOTSTRAP NORMAL STOP ****')
end program