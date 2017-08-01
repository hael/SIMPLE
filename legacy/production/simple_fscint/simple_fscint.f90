program simple_fscint
use simple_cmdline ! singleton
use simple_jiffys, only: file2rarr, simple_end, alloc_err
use simple_math,   only: get_resolution, get_lplim, ratint
use simple_params, only: params
use simple_image,  only: image
implicit none
type(params)      :: p
type(image)       :: img
real, allocatable :: res1(:), res2(:), fsc1(:), fsc2(:)
integer           :: k, alloc_stat
real              :: res0143, res05, dy
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_FSCINT smpd=<sampling distance(in A)> box=<image size(in pixels)>'
    write(*,'(a)') ' fsc=<fsc_state1.bin>'
    stop
endif
call parse_cmdline
call cmdcheckvar('smpd', 1)
call cmdcheckvar('box',  2)
call cmdcheckvar('fsc',  3)
call cmdcheck
p = params() ! parameters generated
call img%new([p%boxmatch,p%boxmatch,1], p%smpd)
res1 = img%get_res()
fsc1 = file2rarr(p%fsc)
if( size(fsc1) /= size(res1) ) write(*,'(a)')'nonconformable fsc sampling; simple_fscint'
call img%new([p%box,p%box,1], p%smpd)
res2 = img%get_res()
allocate(fsc2(size(res2)), stat=alloc_stat)
call alloc_err('In: simple_fscint', alloc_stat)
do k=1,size(res1) 
    write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res1(k), '>>> FSC:', fsc1(k)
end do
! get & print resolution
call get_resolution(fsc1, res1, res05, res0143)
write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
! get & print low-pass limit
k = get_lplim(fsc1)
write(*,'(A,1X,F6.2)') '>>> LOW-PASS LIMIT:', res1(k)
! interpolate the FSC with box-sampling
do k=1,size(res1) 
    call ratint(res1, fsc1, res2(k), fsc2(k), dy)
end do
! END GRACEFULLY
call simple_end('**** SIMPLE_FSCINT NORMAL STOP ****')
end program
