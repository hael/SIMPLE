!------------------------------------------------------------------------------!
! SIMPLE               Elmlund & Elmlund Lab         simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Test program for simple CUDA FFT
!
! @author
!
!
! DESCRIPTION:
!> CUDA FFT implementation -- for PGI

!
! REVISION HISTORY:
! 27 Nov 2017 - Initial Version

!------------------------------------------------------------------------------
program simple_test_cufft
    include 'simple_lib.f08'

    use simple_cufft_test
#if defined(PGI)
    use simple_cufft_test2
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
    write (*,'(A)') 'SIMPLE CUFFT TEST'
    errflag=.true.
    n=10
    allocate(x(n),y(n),y1(n))

    !! CUFFT testing
    write (*,'(A)') 'TESTING FFTW INTERFACE'
    t1=tic()
    call test_fftw(2048,2048,1,.false.,1)
    write (*,'(A,F15.8)') 'FFTW 2D (2048x2048) FFT/IFFT, First pass ', toc(t1)
    t1=tic()
    write (*,'(A)') 'SIMPLE_CUDA image class fft with fftw '

    call test_fftw(2048,2048,1,.false.,5)
    !   print *,'fftw 2D iteration/time',i,toc()

    write (*,'(A,F15.8)') 'FFTW 2D (2048x2048) FFT/IFFT, 5 iterations, ave time (sec):', toc(t1)/5.0
    t1=tic()
    call test_fftw(256,256,100,.false.,1)
    write (*,'(A,F15.8)') 'FFTW 3D (256x256x100) FFT/IFFT, First pass ', toc(t1)
    t1=tic()
    call test_fftw(256,256,100,.false.,5)
    !   print *,'fftw 3D iteration/time',i,toc()
    write (*,'(A,F15.8)') 'FFTW 3D (256x256x100) FFT/IFFT, 5 iterations, ave time (sec):', toc(t1)/5.0



#if defined(PGI)
    write (*,'(A)') 'TESTING CUFFT INTERFACE'
    t1=tic()
    call test_pgi1(2048,2048,1,.false.,1)
    write (*,'(A,F9.4)') 'PGI 2D (2048x2048) FFT/IFFT, First pass', toc(t1)
    istat = cudaDeviceSynchronize()
    istat = istat + cudaThreadSynchronize()
    if(istat /= 0)call simple_cuda_stop("In simple_test_cuda  ",__FILENAME__,__LINE__)
    t1=tic()
    write (*,'(A)') 'SIMPLE_CUDA image class fft with cufft '
!    do i=1,5
        call test_pgi1(2048,2048,1,.false.,5)
        ! print *,'pgi 2D iteration/time',i,toc()
!    end do
    write (*,'(A,F9.4)') 'PGI 2D (2048x2048) FFT/IFFT 5 iterations, ave time (sec)', toc(t1)/5.0
    istat = cudaDeviceSynchronize()
    istat = istat + cudaThreadSynchronize()
    if(istat /= 0)call simple_cuda_stop("In simple_test_cuda  ",__FILENAME__,__LINE__)
    t1=tic()
    call test_pgi1(256,256,100,.false.,1)
    write(*,'(A,F9.4)') 'PGI 3D (256x256x100) FFT/IFFT, First pass', toc(t1)
    istat = cudaDeviceSynchronize()
    istat = istat + cudaThreadSynchronize()
    if(istat /= 0)call simple_cuda_stop("In simple_test_cuda  ",__FILENAME__,__LINE__)
    t1=tic()
!    do i=1,5
        call test_pgi1(256,256,100,.false.,5)
        ! print *,'pgi 3D iteration/time',i,toc()
!    end do
    write (*,'(A,F9.4)') 'PGI 3D (256x256x100) FFT/IFFT, 5 iterations, ave time (sec) ', toc(t1)/5.0

    istat = cudaDeviceSynchronize()
    istat = istat + cudaThreadSynchronize()
    if(istat /= 0)call simple_cuda_stop("In simple_test_cuda  ",__FILENAME__,__LINE__)

#ifdef USE_OPENACC
    write(*,*) '  Testing cuFFT OpenACC, streaming,  device-managed memory, 2D'
    call exec_oacc_cufft2d()
#endif






    write (*,'(A)') 'TESTING CUFFT INTERFACE'
    call cufftTest
    call cufftTest3D
    ! write (*,'(A)') 'SIMPLE_CUDA CUFFT 2D test '
    ! do i=1,10
    !     istat = cudaDeviceSynchronize()
    !     istat = cudaThreadSynchronize()
    !     call exec_cufft2D_test
    ! enddo
    ! write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test):  Average 2D FFT (450x450) sec', ctimer%tocU(ev2)/10.0
    ! call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
    ! print *,'Standard timer', toc()
    ! ev2=ctimer%ticU
    ! write (*,'(A)') 'SIMPLE_CUDA CUFFT 3D test '

    ! do i=1,10
    !     istat = cudaDeviceSynchronize()
    !     istat = cudaThreadSynchronize()
    !     call exec_cufft3D_test
    ! end do
    ! write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test):  Average 3D FFT (150x150x150) sec', ctimer%tocU(ev2)/10.0

    ! call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)

    ! print *,'Standard timer', toc()

    !  t1=tic()
    !  ev2=ctimer%ticU
    ! write (*,'(A)') 'SIMPLE_CUDA image class fft with fftw '
    ! do i=1,10
    !     call test_fftw(200,20*i,10*i,.false.)
    !      print *,'fftw ',i,toc()
    ! end do
    ! write (*,'(A,F9.4)') 'FFTW  Average 3D FFT sec', ctimer%tocU(ev2)/10.0
    ! print *,'Standard timer', toc(t1)
    ! call test_pgi4(8,8,1,.true.)
    ! call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
    ! call test_pgi3(8,8,1,.true.)
    !   print *,'simple_test_cuda:  Standard timer', toc()
    !    !return
    !    t1=tic()
    !    write (*,'(A)') 'SIMPLE_CUDA image class fft with pgi '
    !    ev2=ctimer%ticU
    !        call test_pgi(20*i,20*i,10*i,.true.)
    !        print *,'pgi fft ',i,toc()
    !    end do
    print *," Press enter to test_pgi_ifft"
    read(*,*)
    call test_pgi_ifft(.true.)
    !   ! write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test):  Average Simple_CUFFT test (300x300) sec', ctimer%tocU(ev2)/10.
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
    ! write (*,'(A,F9.4)') 'SIMPLE_CUDA timer (exec_cufft2D_test): ', ctimer%tocU(ev1)

#ifdef CUDA_MANAGED_MEM
    call managed_cufft2dTest()
#endif

#ifdef CUDA_MULTIGPU
    call test_multigpu_C2C_3Dfft
#endif


#else
    print *, '>> SKIPPING: The tests for CUDA '
#endif



#ifdef PGI
    write (*,'(A)') 'TESTING FFTW/CUFFT Worst-Case '
    do i=1,10
        x(i)=2**(i+4)-1
        t1=tic()
        call test_fftw(2**(i+4)-1,2**(i+4)-1,1,.false.)
        y(i)=toc(t1)
    end do
    do  i=1,10
        x(i)=2**(i+4)-1
        t1=tic()
        call test_pgi(2**(i+4)-1,2**(i+4)-1,1,.false.)
        y1(i)=toc(t1)

    end do
    call plot(x,y,x,y1,terminal='dumb')
#endif
    deallocate(x,y,y1)


end program simple_test_cuda
