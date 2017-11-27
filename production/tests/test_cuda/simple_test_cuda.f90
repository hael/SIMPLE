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
    use simple_image,            only: image
    use gnufor2
    use simple_cuda_tests
#if defined(PGI)
    use simple_timer_cuda
    use cudafor
    implicit none
    type (cudaDeviceProp) :: prop
    type (timer_cuda) :: ctimer
    type (cudaEvent)  :: ev1,ev2
#endif
    integer :: i, ierr,istat, cuVer, cuMem, cuFree,n
    real, allocatable :: x(:),y(:),y1(:)
    logical :: errflag
    integer(timer_int_kind) :: t1
    write (*,'(A)') 'SIMPLE CUDA TEST'
    errflag=.true.
    n=10
    allocate(x(n),y(n),y1(n))

#if defined(PGI)
    write (*,'(A)') 'TESTING CUDAFOR INTERFACE'
    call cuda_query_version
    call cuda_query_devices
    call cuda_query_peak_bandwidth

    write (*,'(A)') 'TESTING CUDA multi-GPU capability'
    call cuda_query_p2pAccess(errflag)
    call test_minimal_P2P(errflag)
    call test_transposeP2P

    errflag=.true.
    call test_cuda_precision(errflag)
    call test_acc(errflag)

#endif



#if 0
    write (*,'(A)') 'TESTING CUDAFOR TIMING'
    call ctimer%nowCU()
    write (*,'(A)') 'SIMPLE_CUDA timer setup'

    t1=tic()
    ev1=ctimer%ticU()
    ev2=ctimer%ticU()
    write (*,'(A)') 'SIMPLE_CUDA timer CPU/CUDA', toc(t1), ctimer%tocU(ev1), ctimer%tocU(ev2)
    call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
#endif


end program simple_test_cuda
