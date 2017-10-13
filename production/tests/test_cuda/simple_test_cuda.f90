!------------------------------------------------------------------------------!
! SIMPLE               Elmlund & Elmlund Lab         simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Test program for simple CUDA
!
! @author
!
!
! DESCRIPTION:
!> CUDA implementation -- for PGI

!
! REVISION HISTORY:
! 06 Oct 2017 - Initial Version

!------------------------------------------------------------------------------
program simple_test_cuda
include 'simple_lib.f08'
#ifdef PGI
use simple_cufft_tester
use simple_image,            only: test_image_pgi_cufft
use simple_timer_cuda
use cudafor
implicit none
type(timer_cuda) :: ctimer
type(cudaEvent)  :: res

call ctimer%nowCU
write (*,'(A)') 'SIMPLE_CUDA timer setup'
call ctimer%CUtimer_setup
res=ctimer%ticU
write (*,'(A)') 'SIMPLE_CUDA CUFFT 2D test '
call exec_cufft2D_test
call test_image_pgi_cufft(.true.)
write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test): ', ctimer%tocU(res)
#else
print *, '>> SKIPPING: This test is only for PGI compiler'
#endif

end program simple_test_cuda
