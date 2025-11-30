program simple_test_cc_gradient
! include 'simple_lib.f08'
! use simple_cmdline,          only: cmdline
! use simple_image,            only: image
! use simple_parameters,       only: parameters
! implicit none
! type(cmdline)       :: cline
! type(parameters)    :: p
! type(image)         :: img_copy, img_ori, img_grad
! logical             :: be_verbose=.false.
! integer, parameter  :: N_ITERS = 1
! integer             :: i
! real                :: norm_copy, alpha, xd, xy, yd, dd, yy, norm_ori
! if( command_argument_count() < 3 )then
!     write(logfhandle,'(a)',advance='no') 'simple_test_shiftsrch stk=<particles.ext> mskdiam=<mask radius(in pixels)>'
!     write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
!     stop
! endif
! call cline%parse_oldschool
! call cline%checkvar('stk',      1)
! call cline%checkvar('mskdiam',  2)
! call cline%checkvar('smpd',     3)
! call cline%check
! be_verbose = .false.
! if( cline%defined('verbose') )then
!     if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
!         be_verbose = .true.
!     endif
! endif
! call p%new(cline)
! call img_ori%new(p%ldim, p%smpd)
! call img_ori%read(p%stk, 1)
! call img_copy%copy(img_ori)
! call img_grad%copy(img_ori)
! call img_ori%fft
! call img_copy%fft
! call img_grad%fft
! call img_copy%set_cmat(sum(img_ori%get_cmat())/sum(p%ldim))
! ! call img_copy%write('shifted.mrc', 1)
! print *, cc_nonorm(img_ori, img_copy)
! print *, cc_here(img_ori, img_ori)
! norm_ori = sqrt(sum(csq_fast(img_ori%get_cmat())))
! print *, cc_nonorm(img_ori, img_ori), norm_ori**2
! print *, cc_here(img_ori, img_copy)
! call img_copy%ifft
! call img_copy%write('before_grad.mrc', 1)
! call img_copy%fft
! do i = 1, N_ITERS
!     ! computing the gradient
!     norm_copy = sqrt(sum(csq_fast(img_copy%get_cmat())))
!     call img_grad%set_cmat(real(img_ori%get_cmat() - cc_nonorm(img_ori, img_copy) * img_copy%get_cmat()/norm_copy**2)/norm_copy/norm_ori + complex(0., 0.) )
!     xd = cc_nonorm(img_ori,  img_grad)
!     xy = cc_nonorm(img_ori,  img_copy)
!     yd = cc_nonorm(img_copy, img_grad)
!     dd = cc_nonorm(img_grad, img_grad)
!     yy = cc_nonorm(img_copy, img_copy)
!     alpha = - (xd * yy - xy * yd)/(xd * yd - xy * dd)
!     call img_copy%set_cmat(img_copy%get_cmat() + alpha * img_grad%get_cmat())
! enddo
! print *, cc_here(img_ori, img_copy)
! call img_copy%ifft
! call img_copy%write('grad.mrc', 1)

! contains
!     function cc_here(img1, img2) result(f)
!         type(image), intent(in) :: img1, img2
!         real :: norm1, norm2, f
!         norm1 = sum(csq_fast(img1%get_cmat()))
!         norm2 = sum(csq_fast(img2%get_cmat()))
!         f     = real(sum(img1%get_cmat() * conjg(img2%get_cmat()))) / sqrt( norm1 * norm2 )
!     end function cc_here

!     function cc_nonorm(img1, img2) result(f)
!         type(image), intent(in) :: img1, img2
!         real :: f
!         f = real(sum(img1%get_cmat() * conjg(img2%get_cmat())))
!     end function cc_nonorm
end program simple_test_cc_gradient
