program simple_test_fm
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_cmdline,          only: cmdline
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_polarizer,        only: polarizer
use simple_pftcc_shsrch_fm
use simple_corrmat
implicit none
type(cmdline)          :: cline
type(cmdline), allocatable :: completed_jobs_clines(:)
type(parameters)       :: p
type(image):: img1, img2, img3, img4, imgf2, imgf3, img
type(image), allocatable :: imgs(:)
type(polarft_corrcalc)   :: pftcc
type(polarizer)          :: pol
type(pftcc_shsrch_fm)    :: fmt
character(len=:), allocatable :: fname1, exec
character(len=LONGSTRLEN)     :: fname
character(len=STDLEN)         :: old_exec
real, allocatable :: a(:,:), R(:,:), X(:,:), Y(:,:)
real, allocatable :: T(:,:)
logical, allocatable :: YM(:,:)
real                   :: scale, v, smpd, minmax(2), trs, shift(2), corr, ave, sdev, maxv, minv, dfx,dfy,angast
real                   :: ang,offset(2),cc,ang2
integer :: N,M
integer                :: ldim1(3), ldim2(3), ldim_box(3), nptcls, i,j,nn,ithr, ncls, ncls_rejected, irot, nptcls_glob
if( command_argument_count() < 1 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_fast_corrcalc stk=<particles.ext> mskdiam=<mask radius(in pixels)>'
    write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
    stop
endif

N=256
M=14

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
do i = 1,2*M
    if( YM(1,i) )then
        Y(1,i) = -Y(1,i) !!!
        call transform_img(imgs(i),[X(1,i),Y(1,i)], R(1,i), YM(i,1))
    else
        call transform_img(imgs(i),[X(1,i),Y(1,i)], R(1,i), YM(i,1))
    endif
    print *,i,A(1,i),X(1,i),Y(1,i),R(1,i),YM(1,i)
    call imgs(i)%mask(120.,'soft',backgr=0.)
    call imgs(i)%write('restored_imgs.mrc',i)
enddo


contains

    subroutine transform_img( img, shift, angle, mirror )
        class(image), intent(inout) :: img
        real,         intent(in)    :: shift(2), angle
        logical,      intent(in)    :: mirror
        if( mirror )then
            call img%ifft
            call img%mirror('y')
        endif
        call img%fft
        call img%shift2Dserial(shift)
        call img%ifft
        call img%rtsq(angle, 0.,0.)
    end subroutine transform_img

end program simple_test_fm
