program simple_test_fm
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_cmdline,          only: cmdline
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_corrmat
implicit none
type(cmdline)            :: cline
type(parameters)         :: p
type(image)              :: img1, img2, img3, img4
type(image), allocatable :: imgs(:)
real,        allocatable :: a(:,:), R(:,:), X(:,:), Y(:,:), T(:,:)
logical,     allocatable :: YM(:,:)
real    :: rotmat(2,2), sh(2), ang
integer :: N,M, i,j
logical :: mirr
if( command_argument_count() < 1 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_fm stk=<particles.ext> mskdiam=<mask radius(in pixels)>'
    write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
    stop
endif

N=256 ! image size
M=14  ! number of images

call cline%set('ctf',   'no')
call cline%set('objfun','cc')
call cline%set('smpd',1.)
call cline%set('box', N)
call cline%set('trs', real(N)/4.)
call cline%set('lp', 6.)
call cline%set('hp', 100.)
call cline%set('sh_inv', 'yes')
call p%new(cline)
p%ldim = [N,N,1]

! some dummy image
call img1%soft_ring(p%ldim, 1., 16.)
call img2%soft_ring(p%ldim, 1., 24.)
call img3%soft_ring(p%ldim, 1., 32.)
call img4%soft_ring(p%ldim, 1., 64.)
call img1%fft
call img2%fft
call img3%fft
call img4%fft
call img1%shift2Dserial([ 16.,-32.])
call img2%shift2Dserial([ 64., 0.])
call img3%shift2Dserial([-32., 16.])
call img1%add(img2)
call img1%add(img3)
call img1%add(img4)

! some transformed dummy images, image 1 is the original
allocate(imgs(2*M), T(2*M,3))
do i = 1,M
    call imgs(i)%copy(img1)
    j = i-1
    T(i,:) = real([2*j, -2*j, j*24])
    call transform_img(imgs(i), T(i,1:2), T(i,3), .false.)
enddo
do i = M+1,2*M
    j = i - M - 1
    call imgs(i)%copy(img1)
    T(i,:) = real([-2*j, 2*j, j*24])
    call transform_img(imgs(i), T(i,1:2), T(i,3), .true.)
enddo
do i = 1,2*M
    call imgs(i)%mask(120.,'soft',backgr=0.)
    call imgs(i)%write('imgs.mrc',i)
enddo

! restore images
call calc_inplane_fast_dev( imgs, p%hp, p%lp, A, R, X, Y, YM )

! write out restored images
call img1%ifft
do i = 1,2*M
    ! imgs -> img1
    call transform_img(imgs(i),[X(1,i),Y(1,i)], R(1,i), YM(i,1))
    call imgs(i)%mask(120.,'soft',backgr=0.)
    call imgs(i)%write('restored_imgs.mrc',i)

    ! img1 -> imgs
    if( YM(1,i) )then
        ang = R(1,i)
        call rotmat2d(ang, rotmat)
        sh     = matmul([X(1,i),Y(1,i)], transpose(rotmat))
        sh(1)  = -sh(1)
    else
        ang = 360.-R(1,i)
        call rotmat2d(ang, rotmat)
        sh = -matmul([X(1,i),Y(1,i)], rotmat)
    endif
    mirr = i > M
    print *,i,'mirror,   truth vs calc: ',  mirr,    '|', YM(1,i)
    print *,i,'offset,   truth vs calc: ',  T(i,1:2),'|', sh
    print *,i,'rotation, truth vs calc: ',  T(i,3),  '|', ang
    call img2%copy_fast(img1)
    call transform_img(img2, sh, ang, YM(i,1))
    call img2%mask(120.,'soft',backgr=0.)
    call img2%write('reverse_restored_imgs.mrc',i)
enddo

contains

    subroutine transform_img( img, shift, angle, mirror )
        class(image), intent(inout) :: img
        real,         intent(in)    :: shift(2), angle
        logical,      intent(in)    :: mirror
        if( mirror )then
            call img%ifft
            call img%mirror('y',fourier=.true.)
        endif
        call img%fft
        call img%shift2Dserial(shift)
        call img%ifft
        call img%rtsq(angle, 0.,0.)
    end subroutine transform_img

end program simple_test_fm
