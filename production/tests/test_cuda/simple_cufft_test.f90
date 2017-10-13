module simple_cufft_tester
include 'simple_lib.f08'
use simple_image,       only: image, test_image
implicit none

public :: exec_cufft2D_test
private

! module global constants
integer, parameter :: LDIM(3)=[240,240,1], SQRAD=60, NTST=50, NNOISY=20
real,    parameter :: SMPD=1.77, TRS=10., HP=100.0, LP=8., SNR=0.2

! module global variables

type(image)              :: img, img_shifted
type(image), allocatable :: noisy_imgs(:)
integer                  :: x, y

contains


  subroutine exec_cufft2d_test
#ifdef PGI
    use cufft
      implicit none
      integer, parameter :: n=450
      complex :: a(n,n),b(n,n)
      complex, device :: a_d(n,n), b_d(n,n)
      real :: ar(n,n),br(n,n),x
      real, device :: ar_d(n,n), br_d(n,n)
      integer :: plan, ierr
      logical passing

  a = 1; a_d = a
  ar = 1; ar_d = ar

  ierr = cufftPlan2D(plan,n,n,CUFFT_C2C)
  ierr = ierr + cufftExecC2C(plan,a_d,b_d,CUFFT_FORWARD)
  b = b_d
  write(*,*) maxval(real(b)),sum(b),450*450
  ierr = ierr + cufftExecC2C(plan,b_d,b_d,CUFFT_INVERSE)
  b = b_d
  x = maxval(abs(a-b/(n*n)))
  write(*,*) 'Max error C2C: ', x
  passing = x .le. 1.0e-5

  ierr = ierr + cufftPlan2D(plan,n,n,CUFFT_R2C)
  ierr = ierr + cufftExecR2C(plan,ar_d,b_d)
  ierr = ierr + cufftPlan2D(plan,n,n,CUFFT_C2R)
  ierr = ierr + cufftExecC2R(plan,b_d,br_d)
  br = br_d
  x = maxval(abs(ar-br/(n*n)))
  write(*,*) 'Max error R2C/C2R: ', x
  passing = passing .and. (x .le. 1.0e-5)

  ierr = ierr + cufftDestroy(plan)
  print *,ierr
  passing = passing .and. (ierr .eq. 0)
  if (passing) then
    print *,"Test PASSED"
  else
    print *,"Test FAILED"
 endif
#else
 print *, 'This test is only for PGI compiler'
#endif
 
end subroutine exec_cufft2d_test



end module simple_cufft_tester
