!==Program simple_print_fsc
!
! <print_fsc/begin> is a program for printing the binary FSC files produced by PRIME3D <print_fsc/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_print_fsc
use simple_cmdline, only: cmdline
use simple_jiffys,  only: file2rarr, simple_end
use simple_math,    only: get_resolution, get_lplim
use simple_params,  only: params
use simple_image,   only: image
use simple_timing
implicit none
type(params)      :: p
type(image)       :: img
type(cmdline)     :: cline
real, allocatable :: res(:), fsc(:)
integer           :: k
real              :: res0143, res05
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_PRINT_FSC smpd=<sampling distance(in A)> box=<image size(in pixels)>'
    write(*,'(a)') ' fsc=<fsc_state1.bin>'
    stop
endif
call cline%parse
call cline%checkvar('smpd', 1)
call cline%checkvar('box',  2)
call cline%checkvar('fsc',  3)
call cline%check
call cline%set('prg', 'print_fsc')
p = params(cline) ! parameters generated
call img%new([p%box,p%box,1], p%smpd)
res = img%get_res() 
fsc = file2rarr(p%fsc)
do k=1,size(fsc) 
    write(*,'(A,1X,F6.2,1X,A,1X,F15.3)') '>>> RESOLUTION:', res(k), '>>> FSC:', fsc(k)
end do
! get & print resolution
call get_resolution(fsc, res, res05, res0143)
write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
! END GRACEFULLY
call simple_end('**** SIMPLE_PRINT_FSC NORMAL STOP ****')
end program
