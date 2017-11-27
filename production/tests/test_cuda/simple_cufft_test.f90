module simple_cufft_test
    include 'simple_lib.f08'
    use simple_image,       only: image, test_image
    use simple_fftshifter
    use gnufor2
    implicit none

#ifdef PGI
    public :: exec_cufft2D_test, exec_cufft3D_test ,exec_oacc_cufft2d,cufftTest,cufftTest3D!, fft_pgi_stream_test

#ifdef CUDA_MANAGED_MEM
    public :: test_managed_cufft2dTest
#endif

#ifdef CUDA_MULTIGPU
    public :: test_multigpu_C2C_3Dfft
#endif

#endif
    public :: test_fftw
    private

    ! module global constants
    integer, parameter :: SQRAD=60, NTST=50, NNOISY=20
    real,    parameter :: SMPD=1.77, TRS=10., HP=100.0, LP=8., SNR=0.2

    real(kind=sp), parameter :: tolerance = 1.e-05_sp
    ! module global variables

    ! type(image)              :: img, img_shifted
    ! type(image), allocatable :: noisy_imgs(:)
    ! integer                  :: x, y
    ! USE_OPENACC
#include "simple_local_flags.inc"
contains

#ifdef PGI


    subroutine exec_cufft2d_test

        use cufft
        implicit none
        integer, parameter :: n=450
        complex, allocatable :: a(:,:),b(:,:)
        complex, allocatable, device :: a_d(:,:), b_d(:,:)
        real, allocatable :: ar(:,:),br(:,:)
        real, allocatable, device :: ar_d(:,:), br_d(:,:)
        integer :: plan, ierr
        real::x
        logical passing

        allocate( a(n,n),b(n,n),a_d(n,n),b_d(n,n))
        allocate( ar(n,n),br(n,n),ar_d(n,n),br_d(n,n) )
        a = 1; a_d = a
        ar = 1; ar_d = ar

        ierr = cufftPlan2D(plan,n,n,CUFFT_C2C)
        ierr = ierr + cufftExecC2C(plan,a_d,b_d,CUFFT_FORWARD)
        b = b_d
        if(verbose)write(*,*) maxval(real(b)),sum(b),450*450
        ierr = ierr + cufftExecC2C(plan,b_d,b_d,CUFFT_INVERSE)
        b = b_d
        x = maxval(abs(a-b/(n*n)))
        if (verbose)write(*,*) 'Max error C2C: ', x
        passing = x .le. 1.0e-5

        ierr = ierr + cufftPlan2D(plan,n,n,CUFFT_R2C)
        ierr = ierr + cufftExecR2C(plan,ar_d,b_d)
        ierr = ierr + cufftPlan2D(plan,n,n,CUFFT_C2R)
        ierr = ierr + cufftExecC2R(plan,b_d,br_d)
        br = br_d
        x = maxval(abs(ar-br/(n*n)))
        if(verbose)write(*,*) 'Max error R2C/C2R: ', x
        passing = passing .and. (x .le. 1.0e-5)

        ierr = ierr + cufftDestroy(plan)
        if (ierr /= 0) print *, " exec_cufft2d_test error: ", ierr
        passing = passing .and. (ierr .eq. 0)
        if (passing) then
            if(verbose) print *,"Test PASSED"
        else
            print *,"Test FAILED"
        endif
        deallocate(a,b,a_d,b_d,ar,br,ar_d,br_d)
        !        print *, 'This test is only for PGI compiler'
    end subroutine exec_cufft2d_test

    subroutine exec_cufft3d_test

        use cufft
        implicit none
        integer, parameter :: n=150,m=150,o=150
        complex, allocatable :: a(:,:,:),b(:,:,:)
        complex, allocatable, device :: a_d(:,:,:), b_d(:,:,:)
        real, allocatable :: ar(:,:,:),br(:,:,:)
        real, allocatable, device :: ar_d(:,:,:), br_d(:,:,:)
        integer :: plan, ierr
        real :: x
        logical passing

        allocate(a(n,m,o),b(n,m,o), a_d(n,m,o), b_d(n,m,o), ar(n,m,o),br(n,m,o),&
            ar_d(n,m,o), br_d(n,m,o))
        a = 1; a_d = a
        ar = 1; ar_d = ar

        ierr = cufftPlan3D(plan,n,m,o,CUFFT_C2C)
        ierr = ierr + cufftExecC2C(plan,a_d,b_d,CUFFT_FORWARD)
        b = b_d
        if(verbose)write(*,*) maxval(real(b)),sum(b),450*450
        ierr = ierr + cufftExecC2C(plan,b_d,b_d,CUFFT_INVERSE)
        b = b_d
        x = maxval(abs(a-b/(n*m*o)))
        if(verbose) write(*,*) 'Max error C2C: ', x
        passing = x .le. 1.0e-5

        ierr = ierr + cufftPlan3D(plan,n,m,o,CUFFT_R2C)
        ierr = ierr + cufftExecR2C(plan,ar_d,b_d)
        ierr = ierr + cufftPlan3D(plan,n,m,o,CUFFT_C2R)
        ierr = ierr + cufftExecC2R(plan,b_d,br_d)
        br = br_d
        x = maxval(abs(ar-br/(n*n*n)))
        if(verbose)write(*,*) 'Max error R2C/C2R: ', x
        passing = passing .and. (x .le. 1.0e-5)

        ierr = ierr + cufftDestroy(plan)
        if (ierr /= 0) print *, " exec_cufft2d_test error: ", ierr
        passing = passing .and. (ierr .eq. 0)
        if (passing) then
            if(verbose) print *,"Test PASSED"
        else
            print *,"Test FAILED"
        endif

        ! print *, 'This test is only for PGI compiler'
        deallocate(a,b,a_d,b_d,ar,br,ar_d,br_d)

    end subroutine exec_cufft3d_test

    subroutine exec_oacc_cufft2d
        use cufft
        use openacc
        integer, parameter :: m=768, n=512
        complex, allocatable  :: a(:,:),b(:,:),c(:,:)
        real, allocatable     :: r(:,:),q(:,:)
        real :: xmx
        integer :: iplan1, iplan2, iplan3, ierr

        allocate(a(m,n),b(m,n),c(m,n))
        allocate(r(m,n),q(m,n))
        write(*,*) '  Testing cuFFT OpenACC, streaming,  device-managed memory, 2D'
        a = 1; r = 1
        xmx = -99.0

        ierr = cufftPlan2D(iplan1,m,n,CUFFT_C2C)
        ierr = ierr + cufftSetStream(iplan1,acc_get_cuda_stream(acc_async_sync))
        !$acc host_data use_device(a,b,c)
        ierr = ierr + cufftExecC2C(iplan1,a,b,CUFFT_FORWARD)
        ierr = ierr + cufftExecC2C(iplan1,b,c,CUFFT_INVERSE)
        !$acc end host_data

        ! scale c
        !$acc kernels
        c = c / (m*n)
        !$acc end kernels

        ! Check forward answer
        write(*,*) 'Max error C2C FWD: ', cmplx(maxval(real(b)) - sum(real(b)), &
            maxval(imag(b)))
        ! Check inverse answer
        write(*,*) 'Max error C2C INV: ', maxval(abs(a-c))

        ! Real transform
        ierr = ierr + cufftPlan2D(iplan2,m,n,CUFFT_R2C)
        ierr = ierr + cufftPlan2D(iplan3,m,n,CUFFT_C2R)
        ierr = ierr + cufftSetStream(iplan2,acc_get_cuda_stream(acc_async_sync))
        ierr = ierr + cufftSetStream(iplan3,acc_get_cuda_stream(acc_async_sync))

        !$acc host_data use_device(r,b,q)
        ierr = ierr + cufftExecR2C(iplan2,r,b)
        ierr = ierr + cufftExecC2R(iplan3,b,q)
        !$acc end host_data

        !$acc kernels
        xmx = maxval(abs(r-q/(m*n)))
        !$acc end kernels

        ! Check R2C + C2R answer
        write(*,*) 'Max error R2C/C2R: ', xmx
        deallocate(a,b,c,r,q)
        ierr = ierr + cufftDestroy(iplan1)
        ierr = ierr + cufftDestroy(iplan2)
        ierr = ierr + cufftDestroy(iplan3)

        if (ierr.eq.0) then
            print *,"test PASSED"
        else
            print *,"test FAILED"
        endif


    end subroutine exec_oacc_cufft2d

    !     Copyright (c) 2017, NVIDIA CORPORATION.  All rights reserved.
    !
    ! NVIDIA CORPORATION and its licensors retain all intellectual property
    ! and proprietary rights in and to this software, related documentation
    ! and any modifications thereto.  Any use, reproduction, disclosure or
    ! distribution of this software and related documentation without an express
    ! license agreement from NVIDIA CORPORATION is strictly prohibited.
    !
    subroutine cufftTest
        use simple_cufft

        complex, allocatable :: a(:,:),b(:,:)
        complex, device, allocatable :: a_d(:,:), b_d(:,:)
        real, allocatable :: rtimes(:)
        integer i,j,it
        integer :: nerrors
        integer(8)::t1
        real :: minr, mini
        integer,parameter :: n=600, iters=30
        integer :: plan, planType
        allocate(rtimes(iters-1))
        do it=1,iters
            t1=tic()
            ! allocate arrays on the host
            allocate(a(n,n), b(n,n))

            ! allocate arrays on the device
            allocate(a_d(n,n), b_d(n,n))

            !initialize arrays on host
            do i = 1, n
                do j=1,n
                    a(i,j) = cmplx(cos((i-1) * atan2(0.0,-1.0) / n),cos((j-1) * atan2(0.0,-1.0) / n))
                end do
            end do

            !copy arrays to device
            a_d = a

            ! Print initial array
            !print *, "Array A:"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") a

            ! set planType to either single or double precision
            !if (sp == singlePrecision) then
            planType = CUFFT_C2C
            !else
            !   planType = CUFFT_Z2Z
            !endif

            ! initialize the plan and execute the FFTs.

            call cufftPlan2D(plan,n,n,planType)
            call cufftExec(plan,planType,a_d,b_d,CUFFT_INVERSE)

            ! Copy results back to host
            b = b_d
            !print *, "Forward B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

            call cufftExec(plan,planType,b_d,b_d,CUFFT_FORWARD)

            ! Copy results back to host
            b = b_d
            !print *, "Inverse B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

            ! Scale
            b = b / (n*n)
            !print *, "Scaled B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b
            if (it .ne. 1 ) rtimes(it-1) =  toc(t1)
            nerrors = 0;minr=1e20;mini=1e20
            do i = 1, n
                do j=1,n
                    if ( abs(real(a(i,j)) - real(b(i,j))) .gt. 2.0e-6 ) then
                        if ( abs(real(a(i,j)) - real(b(i,j))) < minr ) then
                            minr = abs(real(a(i,j)) - real(b(i,j)))
                        end if
                        nerrors = nerrors + 1
                    endif
                    if ( abs(imag(a(i,j)) - imag(b(i,j))) .gt. 2.0e-6 ) then
                        if ( abs(imag(a(i,j)) - imag(b(i,j))) < mini ) then
                            mini = abs(imag(a(i,j)) - imag(b(i,j)))
                        end if
                        nerrors = nerrors + 1
                    endif
                end do
            end do
            if (nerrors .ne. 0) then
                print *, "Test FAILED", minr, mini
            else
                if (it == 1 ) then
                    print *, "Test PASSED (first)",toc(t1)
                else
                    print *, "Test PASSED", rtimes(it-1)
                endif

            endif

            !release memory on the host and device
            deallocate(a, b, a_d, b_d)

            ! Destroy the plan
            call cufftDestroy(plan)
        end do
        print *," cufftTest  600x600 Cmplx FWD/BWD In-place, average runtime (excl first)", sum(rtimes)/real(iters-1)
        deallocate(rtimes)
    end subroutine cufftTest


    !     Copyright (c) 2017, NVIDIA CORPORATION.  All rights reserved.
    !
    ! NVIDIA CORPORATION and its licensors retain all intellectual property
    ! and proprietary rights in and to this software, related documentation
    ! and any modifications thereto.  Any use, reproduction, disclosure or
    ! distribution of this software and related documentation without an express
    ! license agreement from NVIDIA CORPORATION is strictly prohibited.
    !
    subroutine cufftTest3D
        use simple_cufft

        complex, allocatable :: a(:,:,:),b(:,:,:)
        complex, device, allocatable :: a_d(:,:,:), b_d(:,:,:)
        real, allocatable :: rtimes(:)
        integer i,j,k,it
        integer :: nerrors
        integer(8)::t1
        real :: minr, mini
        integer,parameter :: n=64, iters=10
        integer :: plan, planType
        allocate(rtimes(iters-1))
        do it=1,iters
            t1=tic()
            ! allocate arrays on the host
            allocate(a(n,n,n), b(n,n,n))

            ! allocate arrays on the device
            allocate(a_d(n,n,n), b_d(n,n,n))

            !initialize arrays on host
            do i = 1, n
                do j=1,n
                    do k=1,n
                        a(i,j,k) = cmplx(cos((i-1) * atan2(0.0,-1.0) / n)*cos((k-n/2)* atan2(0.0,-1.0) / n),&
                            cos((j-1) * atan2(0.0,-1.0) / n)*cos((k-n/2)* atan2(0.0,-1.0) / n))
                    end do
                end do
            end do


            !copy arrays to device
            a_d = a

            ! Print initial array
            !print *, "Array A:"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") a

            ! set planType to either single or double precision
            !if (sp == singlePrecision) then
            planType = CUFFT_C2C
            !else
            !   planType = CUFFT_Z2Z
            !endif

            ! initialize the plan and execute the FFTs.

            call cufftPlan3D(plan,n,n,n,planType)
            call cufftExec(plan,planType,a_d,b_d,CUFFT_INVERSE)

            ! Copy results back to host
            b = b_d
            !print *, "Forward B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

            call cufftExec(plan,planType,b_d,b_d,CUFFT_FORWARD)

            ! Copy results back to host
            b = b_d
            !print *, "Inverse B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

            ! Scale
            b = b / (n*n*n)
            !print *, "Scaled B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b
            if (it .ne. 1 ) rtimes(it-1) =  toc(t1)
            nerrors = 0;minr=1e20;mini=1e20
            do i = 1, n
                do j=1,n
                     do k=1,n
                         if ( abs(real(a(i,j,k)) - real(b(i,j,k))) .gt. 2.0e-6 ) then
                             if ( abs(real(a(i,j,k)) - real(b(i,j,k))) < minr ) then
                                 minr = abs(real(a(i,j,k)) - real(b(i,j,k)))
                        end if
                        nerrors = nerrors + 1
                    endif
                    if ( abs(imag(a(i,j,k)) - imag(b(i,j,k))) .gt. 2.0e-6 ) then
                        if ( abs(imag(a(i,j,k)) - imag(b(i,j,k))) < mini ) then
                            mini = abs(imag(a(i,j,k)) - imag(b(i,j,k)))
                        end if
                        nerrors = nerrors + 1
                    endif
                end do
            end do
            end do
            if (nerrors .ne. 0) then
                print *, "Test FAILED", minr, mini
            else
                if (it == 1 ) then
                    print *, "Test PASSED (first)",toc(t1)
                else
                    print *, "Test PASSED", rtimes(it-1)
                endif

            endif

            !release memory on the host and device
            deallocate(a, b, a_d, b_d)

            ! Destroy the plan
            call cufftDestroy(plan)
        end do
        print *," cufftTest  64^3 Cmplx FWD/BWD In-place, average runtime (excl first)", sum(rtimes)/real(iters-1)
        deallocate(rtimes)
    end subroutine cufftTest3D


    subroutine test_managed_cufft2dTest
        use cudafor
        use cufft
        implicit none
        integer, parameter :: m=768, n=512
#ifdef CUDA_MANAGED_MEM
        complex, managed :: a(m,n),b(m,n)
        real, managed :: ar(m,n),br(m,n)
        real    x
        integer plan, ierr
        logical passing

        a = 1; ar = 1

        ierr = cufftPlan2D(plan,m,n,CUFFT_C2C)
        ierr = ierr + cufftExecC2C(plan,a,b,CUFFT_FORWARD)
        ierr = ierr + cufftExecC2C(plan,b,b,CUFFT_INVERSE)
        ierr = ierr + cudaDeviceSynchronize()
        x = maxval(abs(a-b/(m*n)))
        write(*,*) 'Max error C2C: ', x
        passing = x .le. 1.0e-5

        ierr = ierr + cufftPlan2D(plan,m,n,CUFFT_R2C)
        ierr = ierr + cufftExecR2C(plan,ar,b)
        ierr = ierr + cufftPlan2D(plan,m,n,CUFFT_C2R)
        ierr = ierr + cufftExecC2R(plan,b,br)
        ierr = ierr + cudaDeviceSynchronize()
        x = maxval(abs(ar-br/(m*n)))
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
#endif
    end subroutine test_managed_cufft2dTest


#ifdef CUDA_MULTIGPU
    subroutine test_multigpu_C2C_3Dfft
        use cufft
        ! 3D Complex-to-Complex Transforms using Two GPUs
        ! In this example a three-dimensional complex-to-complex transform is
        ! applied to the input data using two GPUs.

        ! Demonstrate how to use CUFFT to perform 3-d FFTs using 2 GPUs
        complex, allocatable :: host_data_input(:,:,:), host_data_output(:,:,:)
        ! Create an empty plan
        type(cufftPlan) plan
        integer ::  res,nGPUs, whichGPUs(2),worksize(2), nx , ny , nz, size_of_data

        res = cufftCreate(&plan_input);
        if (res /= CUFFT_SUCCESS)then
            print *,"*Create failed"
            return
        end if


        ! cufftXtSetGPUs() - Define which GPUs to use
        nGPUs = 2
        whichGPUs(0) = 0; whichGPUs(1) = 1;
        res = cufftXtSetGPUs (plan_input, nGPUs, whichGPUs)
        if (res /= CUFFT_SUCCESS)then
            print *,"*XtSetGPUs failed"
            return
        end if

        ! Initialize FFT input data


        nx = 64, ny = 128, nz = 32;
        size_of_data = 8 * nx * ny * nz;
        allocate(host_data_input(nx,ny,nz),host_data_output(nx,ny,nz)

        initialize_3d_data (nx, ny, nz, host_data_input, host_data_output)

        ! cufftMakePlan3d() - Create the plan
        res = cufftMakePlan3d (plan_input, nz, ny, nx, CUFFT_C2C, worksize)
        if (res /= CUFFT_SUCCESS)then
            print *,"*MakePlan* failed"
            return
        end if

        ! cufftXtMalloc() - Malloc data on multiple GPUs
        cudaLibXtDesc *device_data_input;
        res = cufftXtMalloc (plan_input, &device_data_input,&
            CUFFT_XT_FORMAT_INPLACE)
        if (res /= CUFFT_SUCCESS) then
            print *,"*XtMalloc failed"
            return
        endif

        ! cufftXtMemcpy() - Copy data from host to multiple GPUs
        res = cufftXtMemcpy (plan_input, device_data_input,&
            host_data_input, CUFFT_COPY_HOST_TO_DEVICE)
        if (res /= CUFFT_SUCCESS) then
            print *,"*XtMemCpy failed"; return
        endif

        ! cufftXtExecDescriptorC2C() - Execute FFT on multiple GPUs
        res = cufftXtExecDescriptorC2C (plan_input, device_data_input,&
            device_data_input, CUFFT_FORWARD)
        if (res /= CUFFT_SUCCESS) then
            print *,"*XtExec failed"; return
        endif

        ! cufftXtMemcpy() - Copy data from multiple GPUs to host
        res = cufftXtMemcpy (plan_input, host_data_output,&
            device_data_input, CUFFT_COPY_DEVICE_TO_HOST)
        if (res /= CUFFT_SUCCESS) then
            print *,"*XtMemCpy failed"; return;
        endif

        ! Print output and check results
        res = output_3d_results (nx, ny, nz,&
            host_data_input, host_data_output)
        if (res /= 0) then
            print *,"*    output_3d_results failed"
            return
        endif


        ! cufftXtFree() - Free GPU memory
        result = cufftXtFree(device_data_input)
        if (res /= CUFFT_SUCCESS) then
            print *,"*XtFree failed"
            return
        endif

        !  cufftDestroy() - Destroy FFT plan
        res = cufftDestroy(plan_input)
        if (res /= CUFFT_SUCCESS) then
            print *,"*Destroy failed "
            return
        endif
        deallocate(host_data_input,host_data_output);


        !Read more at: http://docs.nvidia.com/cuda/cufft/index.html#ixzz4xED7Iv1J
        !Follow us: @GPUComputing on Twitter | NVIDIA on Facebook

    end subroutine test_multigpu_C2C_3Dfft

#endif

#endif

    subroutine test_fftw( ld1, ld2, ld3, doplot ,iter)
!        use simple_fftshifter
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        integer, intent(in),optional:: iter
        type(image)          :: img, img_2
        integer              :: ldim(3),iter_here,it
        real                 :: smpd=2.
        logical              :: passed
        iter_here=1
        if(present(iter))iter_here=iter
        verbose=.false.
        !debug=.true.
        write(*,*) '**info(unit_test, FFTW  test_fftw'
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)

        call img%gauimg(10)
        call img%gauimg2(5, floor(real(ld1)/5.),floor(real(ld2)/3.))
        call img%gauimg2(8, -floor(real(ld1)/5.),floor(real(ld2)/6.))
        call img%add_gauran(0.1)
        if(doplot)call img%vis(geomorsphr=.false.)
        img_2=img
        do it=1,iter_here
            call img_2%fwd_ft()
            call img_2%bwd_ft()
        end do
        write(*,'(a,1ES20.10)') '**info(unit_test, FFTW   L2 norm sum', img.lonesum.img_2
        if(doplot)call img_2%vis()
        call img_2%kill
        call img%kill

    end subroutine test_fftw

end module simple_cufft_test
