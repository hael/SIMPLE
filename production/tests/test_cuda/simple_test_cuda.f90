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
#if defined(PGI)
    use simple_cuda_tests
    use simple_timer_cuda
    use simple_cufft_test
    use simple_cufft_test2

    use cudafor
    implicit none
    type (cudaDeviceProp) :: prop
    type (timer_cuda) :: ctimer
    type (cudaEvent)  :: init,res
#endif
    integer :: i, ierr,istat, cuVer, cuMem, cuFree,n
  !  real, allocatable :: x(:),y(:),y1(:)
    logical :: errflag
    integer(timer_int_kind) :: t1
    write (*,'(A)') 'SIMPLE CUDA TEST'
#if defined(PGI)
    write (*,'(A)') 'TESTING CUDAFOR INTERFACE'
    call cuda_query_version
    call cuda_query_devices
    call cuda_query_peak_bandwidth

 !   call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
    errflag=.true.
    call test_cuda_precision(errflag)
    n=10
   ! allocate(x(n),y(n),y1(n))

#ifdef USE_OPENACC
    call test_acc(errflag)
#endif
    ! write (*,'(A)') 'TESTING CUDAFOR TIMING'
    ! call ctimer%nowCU
    ! write (*,'(A)') 'SIMPLE_CUDA timer setup'
    ! call ctimer%CUtimer_setup
    ! t1=tic()
    ! init=ctimer%ticU
    ! res=ctimer%ticU
    !    call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
#ifdef USE_OPENACC_ONLY
    call exec_oacc_cufft2d
#endif
    call cufftTest
    call cufftTest3D
    ! write (*,'(A)') 'SIMPLE_CUDA CUFFT 2D test '
    ! do i=1,10
    !     istat = cudaDeviceSynchronize()
    !     istat = cudaThreadSynchronize()
    !     call exec_cufft2D_test
    ! enddo
    ! write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test):  Average 2D FFT (450x450) sec', ctimer%tocU(res)/10.0
    ! call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
    ! print *,'Standard timer', toc()
    ! res=ctimer%ticU
    ! write (*,'(A)') 'SIMPLE_CUDA CUFFT 3D test '

    ! do i=1,10
    !     istat = cudaDeviceSynchronize()
    !     istat = cudaThreadSynchronize()
    !     call exec_cufft3D_test
    ! end do
    ! write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test):  Average 3D FFT (150x150x150) sec', ctimer%tocU(res)/10.0

    ! call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)

    ! print *,'Standard timer', toc()

    !  t1=tic()
    !  res=ctimer%ticU
    ! write (*,'(A)') 'SIMPLE_CUDA image class fft with fftw '
    ! do i=1,10
    !     call test_fftw(200,20*i,10*i,.false.)
    !      print *,'fftw ',i,toc()
    ! end do
    ! write (*,'(A,F9.4)') 'FFTW  Average 3D FFT sec', ctimer%tocU(res)/10.0
    ! print *,'Standard timer', toc(t1)
    ! call test_pgi4(8,8,1,.true.)
    ! call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
    ! call test_pgi3(8,8,1,.true.)
    !   print *,'simple_test_cuda:  Standard timer', toc()
    !    !return
    !    t1=tic()
    !    write (*,'(A)') 'SIMPLE_CUDA image class fft with pgi '
    !    res=ctimer%ticU
    !        call test_pgi(20*i,20*i,10*i,.true.)
    !        print *,'pgi fft ',i,toc()
    !    end do
    print *," ready"
    read(*,*)
   ! call test_pgi_ifft(.true.)
    !   ! write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test):  Average Simple_CUFFT test (300x300) sec', ctimer%tocU(res)/10.
    !   ! print *,'Standard timer', toc(t1)
    !
    ! istat = cudaThreadSynchronize()
    !  call simple_cuda_stop("In simple_test_cuda post fft sync ",__FILENAME__,__LINE__)
    !  istat = cudaDeviceSynchronize()
    !   call simple_cuda_stop("In simple_test_cuda post fft sync ",__FILENAME__,__LINE__)
    !call test_pgi1(256,256,1,.true.)
    !   call simple_cuda_stop("In simple_test_cuda post test_pgi2 ",__FILENAME__,__LINE__)
     !call test_pgi3(256,256,1,.true.)
    !   call simple_cuda_stop("In simple_test_cuda post test_pgi3 ",__FILENAME__,__LINE__)
    !   !call test_image_pgi_cufft(.true.)
    !   call simple_cuda_stop("In simple_test_cuda post fft sync ",__FILENAME__,__LINE__)
    ! write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test): ', ctimer%tocU(init)


    t1=tic()

    write (*,'(A)') 'SIMPLE_CUDA image class fft with cufft '
    do i=1,3
        call test_pgi(2048,2048,1,.false.)
       ! print *,'pgi 2D iteration/time',i,toc()
    end do
    write (*,'(A,F9.4)') 'PGI  Average 2D (2048x2048) FFT/IFFT ', toc(t1)/3.0
    t1=tic()
    do i=1,3
        call test_pgi(256,256,128,.false.)
       ! print *,'pgi 3D iteration/time',i,toc()
    end do
    write (*,'(A,F9.4)') 'PGI  Average 3D (128x256x100) FFT/IFFT ', toc(t1)/3.0

    istat = cudaDeviceSynchronize()
    istat = cudaThreadSynchronize()
    if(istat /= 0)call simple_cuda_stop("In simple_test_cuda  ",__FILENAME__,__LINE__)

#else
    print *, '>> SKIPPING: The tests for CUDA are only for the PGI Fortran compiler'
#endif

    t1=tic()
    write (*,'(A)') 'SIMPLE_CUDA image class fft with fftw '
    do i=1,3
        call test_fftw(2048,2048,1,.false.)
     !   print *,'fftw 2D iteration/time',i,toc()
    end do
    write (*,'(A,F15.8)') 'FFTW  Average 2D (2048x2048) FFT/IFFT time (sec):', toc(t1)/3.0

    t1=tic()
    do i=1,3
        call test_fftw(256,256,128,.false.)
     !   print *,'fftw 3D iteration/time',i,toc()
    end do
    write (*,'(A,F15.8)') 'FFTW  Average 3D (128x256x100) FFT/IFFT time (sec):', toc(t1)/3.0

#ifdef PGI

    ! do i=1,10
    !     x(i)=2**(i+4)-1
    !     t1=tic()
    !     call test_fftw(2**(i+4)-1,2**(i+4)-1,1,.false.)
    !     y(i)=toc(t1)
    ! end do
    ! do  i=1,10
    !     x(i)=2**(i+4)-1
    !     t1=tic()
    !     call test_pgi(2**(i+4)-1,2**(i+4)-1,1,.false.)
    !     y1(i)=toc(t1)
    !     istat = cudaDeviceSynchronize()
    !     istat = cudaThreadSynchronize()
    ! end do
    ! call plot(x,y,x,y1,terminal='dumb')
#endif
!  deallocate(x,y,y1)

end program simple_test_cuda
