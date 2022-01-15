program simple_test_io_parallel
include 'simple_lib.f08'
use simple_image,    only: image
use simple_stack_io, only: stack_io
use simple_imgfile,  only: imgfile
implicit none
character(len=:), allocatable :: benchfname
integer,          parameter   :: NVOLRWS   = 10
integer,          parameter   :: NVOLS     = 4
integer,          parameter   :: BOX       = 512
integer,          parameter   :: ONE_M     = 1024**2
integer(kind=8),  parameter   :: NVOLBYTES = NVOLRWS * BOX * BOX * BOX * 4
real,             parameter   :: SMPD      = 1.0
real(kind=c_float),            pointer :: rmat_ptr(:,:,:)=>null()  !< image pixels/voxels (in data)
complex(kind=c_float_complex), pointer :: cmat_ptr(:,:,:)=>null()  !< Fourier components
type(image)                   :: vols(NVOLS)
type(imgfile)                 :: ioimg(NVOLS)
integer(timer_int_kind)       ::  t_vol_w,  t_vol_w_para,  t_vol_r,  t_vol_r_para
real(timer_int_kind)          :: rt_vol_w, rt_vol_w_para, rt_vol_r, rt_vol_r_para
integer :: i, j, fnr

print *, 'simulating volumes'
do i = 1, NVOLS
    call vols(i)%new([BOX,BOX,BOX], SMPD)
    call vols(i)%ran()
    call vols(i)%write('random_vol'//int2str(i)//'.mrc')
end do

print *, 'writing the volume '//int2str(NVOLRWS)//' times'
t_vol_w = tic()
do j = 1, NVOLRWS
    do i = 1, NVOLS
        call vols(i)%write('random_vol'//int2str(i)//'.mrc')
    end do
end do
rt_vol_w = toc(t_vol_w)

print *, 'writing the volume '//int2str(NVOLRWS)//' times in parallel'
t_vol_w_para = tic()
do i = 1, NVOLS
    call ioimg(i)%open('random_vol'//int2str(i)//'.mrc', [BOX,BOX,BOX], SMPD, del_if_exists=.true., formatchar='M', readhead=.false., rwaction='write')
end do
do j = 1, NVOLRWS
    !$omp parallel do default(shared) private(i,rmat_ptr) schedule(static) num_threads(NVOLS)
    do i = 1, NVOLS
        call vols(i)%get_rmat_ptr(rmat_ptr)
        call ioimg(i)%wmrcSlices(1, BOX, rmat_ptr, [BOX,BOX,BOX], .false.) ! .false. is FT status
    end do
    !$omp end parallel do
end do
do i = 1, NVOLS
    call ioimg(i)%close
end do
rt_vol_w_para = toc(t_vol_w_para)

print *, 'reading the volume '//int2str(NVOLRWS)//' times'
t_vol_r = tic()
do j = 1, NVOLRWS
    do i = 1, NVOLS
        call vols(i)%read('random_vol'//int2str(i)//'.mrc')
    end do
end do
rt_vol_r = toc(t_vol_r)

print *, 'reading the volume '//int2str(NVOLRWS)//' times, in parallel'
t_vol_r_para = tic()
do i = 1, NVOLS
    call ioimg(i)%open('random_vol'//int2str(i)//'.mrc', [BOX,BOX,BOX], SMPD, formatchar='M', readhead=.false., rwaction='read')
end do
do j = 1, NVOLRWS
    !$omp parallel do default(shared) private(i,rmat_ptr) schedule(static) num_threads(NVOLS)
    do i = 1, NVOLS
        call vols(i)%get_rmat_ptr(rmat_ptr)
        call ioimg(i)%rSlices(1,BOX,rmat_ptr,is_mrc=.true.)
    end do
    !$omp end parallel do
end do
do i = 1, NVOLS
    call ioimg(i)%close
end do
rt_vol_r_para = toc(t_vol_r_para)

! cleanup
do i = 1, NVOLS
    call del_file('random_vol'//int2str(i)//'.mrc')
end do

! write benchmark stats
benchfname = 'SIMPLE_PARALLEL_IO_BENCH.txt'
call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
write(fnr,'(a)') '*** TIMINGS (s) ***'
write(fnr,'(a,1x,f9.2)') 'volume writes           : ', rt_vol_w
write(fnr,'(a,1x,f9.2)') 'volume writes, parallel : ', rt_vol_w_para
write(fnr,'(a,1x,f9.2)') 'volume reads            : ', rt_vol_r
write(fnr,'(a,1x,f9.2)') 'volume reads, parallel  : ', rt_vol_r_para
write(fnr,'(a)') '*** SPEEDUPS ***'
write(fnr,'(a,1x,f9.2)') 'parallel write speedup  : ', rt_vol_w / rt_vol_w_para
write(fnr,'(a,1x,f9.2)') 'parallel read  speedup  : ', rt_vol_r / rt_vol_r_para
call fclose(fnr)

end program simple_test_io_parallel
