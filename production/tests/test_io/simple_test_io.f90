program simple_test_io
include 'simple_lib.f08'
use simple_image,    only: image
use simple_stack_io, only: stack_io
implicit none
type(image),      allocatable :: imgs(:)
character(len=:), allocatable :: benchfname
type(image)                   :: vol
type(stack_io)                :: stkio_r, stkio_w
integer,          parameter   :: NVOLRWS   = 10
integer,          parameter   :: NSTKRWS   = 10
integer,          parameter   :: NPTCLS    = 1024 * 10
!integer,          parameter   :: BOX       = 512
integer,          parameter   :: BOX       = 48 
integer,          parameter   :: ONE_M     = 1024**2
integer(kind=8),  parameter   :: NSTKBYTES = NSTKRWS * NPTCLS * BOX * BOX * 4
integer(kind=8),  parameter   :: NVOLBYTES = NVOLRWS * BOX    * BOX * BOX * 4
real,             parameter   :: SMPD      = 1.0
integer(timer_int_kind)       ::  t_stk_w,  t_stk_r,  t_vol_w,  t_vol_r
real(timer_int_kind)          :: rt_stk_w, rt_stk_r, rt_vol_w, rt_vol_r, rt_tot
integer :: iptcl, i, fnr
real    :: mb_per_s_stk_w, mb_per_s_stk_r, mb_per_s_vol_w, mb_per_s_vol_r, mb_per_s_w, mb_per_s_r

print *, 'simulating a stack of '//int2str(NPTCLS)//' particles'
allocate(imgs(NPTCLS))
do iptcl = 1, NPTCLS
    call imgs(iptcl)%new([BOX,BOX,1], SMPD)
    call imgs(iptcl)%ran()
end do

print *, 'simulating a volume'
call vol%new([BOX,BOX,BOX], SMPD)
call vol%ran()
call vol%write('random_vol.mrc')

print *, 'writing the stack '//int2str(NSTKRWS)//' times'
rt_tot  = 0.
t_stk_w = tic()
do i = 1, NSTKRWS
    call stkio_w%open('stack_of_random_imgs.mrcs', SMPD, 'write', box=BOX)
    do iptcl = 1, NPTCLS
        call stkio_w%write(iptcl, imgs(iptcl))
    end do
    call stkio_w%close
end do
rt_stk_w = toc(t_stk_w)
rt_tot   = rt_tot + rt_stk_w

print *, 'reading the stack '//int2str(NSTKRWS)//' times'
t_stk_r = tic()
do i = 1, NSTKRWS
    call stkio_r%open('stack_of_random_imgs.mrcs', SMPD, 'read')
    do iptcl = 1, NPTCLS
        call stkio_r%read(iptcl, imgs(iptcl))
    end do
    call stkio_r%close
end do
rt_stk_r = toc(t_stk_r)
rt_tot   = rt_tot + rt_stk_r

print *, 'writing the volume '//int2str(NVOLRWS)//' times'
t_vol_w = tic()
do i = 1, NVOLRWS
    call vol%write('random_vol.mrc')
end do
rt_vol_w = toc(t_vol_w)
rt_tot   = rt_tot + rt_vol_w

print *, 'reading the volume '//int2str(NVOLRWS)//' times'
t_vol_r = tic()
do i = 1, NVOLRWS
    call vol%read('random_vol.mrc')
end do
rt_vol_r = toc(t_vol_r)
rt_tot   = rt_tot + rt_vol_r

! calc MB / s
mb_per_s_stk_w = real(real(NSTKBYTES,dp)             / real(ONE_M,dp) / rt_stk_w)
mb_per_s_stk_r = real(real(NSTKBYTES,dp)             / real(ONE_M,dp) / rt_stk_r)
mb_per_s_vol_w = real(real(NVOLBYTES,dp)             / real(ONE_M,dp) / rt_vol_w)
mb_per_s_vol_r = real(real(NVOLBYTES,dp)             / real(ONE_M,dp) / rt_vol_r)
mb_per_s_w     = real(real(NSTKBYTES + NVOLBYTES,dp) / real(ONE_M,dp) / (rt_stk_w + rt_vol_w))
mb_per_s_r     = real(real(NSTKBYTES + NVOLBYTES,dp) / real(ONE_M,dp) / (rt_stk_r + rt_vol_r))

! write benchmark stats
benchfname = 'SIMPLE_IO_BENCH.txt'
call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
write(fnr,'(a)') '*** TIMINGS (s) ***'
write(fnr,'(a,1x,f9.2)') 'stack  writes : ', rt_stk_w
write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', rt_stk_r
write(fnr,'(a,1x,f9.2)') 'volume writes : ', rt_vol_w
write(fnr,'(a,1x,f9.2)') 'volume reads  : ', rt_vol_r
write(fnr,'(a,1x,f9.2)') 'total  time   : ', rt_tot
write(fnr,'(a)') ''
write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
write(fnr,'(a,1x,f9.2)') 'stack  writes : ', (rt_stk_w / rt_tot) * 100.
write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', (rt_stk_r / rt_tot) * 100.
write(fnr,'(a,1x,f9.2)') 'volume writes : ', (rt_vol_w / rt_tot) * 100.
write(fnr,'(a,1x,f9.2)') 'volume reads  : ', (rt_vol_r / rt_tot) * 100.
write(fnr,'(a,1x,f9.2)') '% accounted for : ', ((rt_stk_w+rt_stk_r+rt_vol_w+rt_vol_r)/rt_tot) * 100.
write(fnr,'(a)') ''
write(fnr,'(a)') '*** READ/WRITE SPEEDS (MB/s) ***'
write(fnr,'(a,1x,f9.2)') 'stack  writes : ', mb_per_s_stk_w
write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', mb_per_s_stk_r
write(fnr,'(a,1x,f9.2)') 'volume writes : ', mb_per_s_vol_w
write(fnr,'(a,1x,f9.2)') 'volume reads  : ', mb_per_s_vol_r
write(fnr,'(a,1x,f9.2)') 'total  writes : ', mb_per_s_w
write(fnr,'(a,1x,f9.2)') 'total  reads  : ', mb_per_s_r
call fclose(fnr)

end program simple_test_io
