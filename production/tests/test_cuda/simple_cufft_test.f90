module simple_cufft_test
    include 'simple_lib.f08'
    use simple_image,       only: image, test_image
    use simple_fftshifter
    use gnufor2
    implicit none

#ifdef PGI
    public :: test_cuda, test_acc, exec_cufft2D_test, exec_cufft3D_test, test_pgi, test_pgi1, &
        test_pgi2, test_pgi3, test_pgi4 ,exec_oacc_cufft2d,cufftTest
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
    subroutine test_cuda(flag)
        use cudafor
        use simple_defs
        implicit none
        logical, intent(inout):: flag
        call test_precision(flag)
        call sum_accuracy(flag)

    contains
        !> Floating-point precision test
        subroutine test_precision(flag)
            logical, intent(inout):: flag

            real :: x, y, dist
            double precision:: x_dp, y_dp, dist_dp
            x=Z'3F1DC57A'
            y=Z'3F499AA3'
            dist= x**2 +y**2

            x_dp=real(x,8)
            y_dp=real(y,8)
            dist_dp= x_dp**2 +y_dp**2

            print *, 'Result with operands in single precision:'
            print '((2x,z8)) ', dist

            print *, 'Result in double precision with operands'
            print *, 'promoted to double precision:'
            print '((2x,z16))', dist_dp

            print *, 'Result in single precision with operands'
            print *, 'promoted to double precision:'
            print '((2x,z8))', real(dist_dp,4)
        end subroutine test_precision

        !>  Floating-point precision test
        subroutine sum_accuracy(flag)
            logical, intent(inout):: flag

            real, allocatable :: x(:)
            real :: sum_intrinsic,sum_cpu, sum_kahan, sum_pairwise, &
                comp, y, tmp
            double precision :: sum_cpu_dp
            integer :: i,inext,icurrent,  N=10000000

            allocate (x(N))
            x=7.

            ! Summation using intrinsic
            sum_intrinsic=sum(x)

            ! Recursive summation
            sum_cpu=0.
            sum_cpu_dp=0.d0
            do i=1,N
                ! accumulator in single precision
                sum_cpu=sum_cpu+x(i)
                ! accumulator in double precision
                sum_cpu_dp=sum_cpu_dp+x(i)
            end do

            ! Kahan summation
            sum_kahan=0.
            comp=0. ! running compensation to recover lost low-order bits

            do i=1,N
                y    = comp +x(i)
                tmp  = sum_kahan + y     ! low-order bits may be lost
                comp = (sum_kahan-tmp)+y ! (sum-tmp) recover low-order bits
                sum_kahan = tmp
            end do
            sum_kahan=sum_kahan +comp

            ! Pairwise summation
            icurrent=N
            inext=ceiling(real(N)/2)
            do while (inext >1)
                do i=1,inext
                    if ( 2*i <= icurrent) x(i)=x(i)+x(i+inext)
                end do
                icurrent=inext
                inext=ceiling(real(inext)/2)
            end do
            sum_pairwise=x(1)+x(2)

            write(*, "('Summming ',i10, &
                ' elements of magnitude ',f3.1)") N,7.
            write(*, "('Sum with intrinsic function       =',f12.1, &
                '   Error=', f12.1)")  &
                sum_intrinsic, 7.*N-sum_intrinsic
            write(*, "('Recursive sum with SP accumulator =',f12.1, &
                '   Error=', f12.1)")  sum_cpu, 7.*N-sum_cpu
            write(*, "('Recursive sum with DP accumulator =',f12.1, &
                '   Error=', f12.1)")  sum_cpu_dp, 7.*N-sum_cpu_dp
            write(*, "('Pairwise sum in SP                =',f12.1, &
                '   Error=', f12.1)")  sum_pairwise, 7.*N-sum_pairwise
            write(*, "('Compensated sum in SP             =',f12.1, &
                '   Error=', f12.1)")  sum_kahan, 7.*N-sum_kahan

            deallocate(x)
        end subroutine sum_accuracy

    end subroutine test_cuda

    ! Step3  Jacobi relaxation
    ! NVIDIA blog example https://github.com/parallel-forall/code-samples
    subroutine test_acc(flag)
#ifdef USE_OPENACC
        use openacc
        logical, intent(inout):: flag
        integer, parameter :: fp_kind=kind(1.0)
        integer, parameter :: n=4096, m=4096, iter_max=1000
        integer :: i, j, iter
        real(fp_kind), dimension (:,:), allocatable :: A, Anew
        real(fp_kind), dimension (:),   allocatable :: y0
        real(fp_kind) :: pi=2.0_fp_kind*asin(1.0_fp_kind), tol=1.0e-5_fp_kind, error=1.0_fp_kind
        real(fp_kind) :: start_time, stop_time

        allocate ( A(0:n-1,0:m-1), Anew(0:n-1,0:m-1) )
        allocate ( y0(0:m-1) )

        A = 0.0_fp_kind

        ! Set B.C.
        y0 = sin(pi* (/ (j,j=0,m-1) /) /(m-1))

        A(0,:)   = 0.0_fp_kind
        A(n-1,:) = 0.0_fp_kind
        A(:,0)   = y0
        A(:,m-1) = y0*exp(-pi)

        call acc_init(acc_device_nvidia)

        write(*,'(a,i5,a,i5,a)') 'Jacobi relaxation Calculation:', n, ' x', m, ' mesh'

        call cpu_time(start_time)

        iter=0

        !$omp parallel do shared(Anew)
        do i=1,m-1
            Anew(0,i)   = 0.0_fp_kind
            Anew(n-1,i) = 0.0_fp_kind
        end do
        !$omp end parallel do

        !$omp parallel do shared(Anew)
        do i=1,n-1
            Anew(i,0)   = y0(i)
            Anew(i,m-1) = y0(i)*exp(-pi)
        end do
        !$omp end parallel do

        !$acc data copy(A), create(Anew)
        do while ( error .gt. tol .and. iter .lt. iter_max )
            error=0.0_fp_kind
            !$omp parallel do shared(m, n, Anew, A) reduction( max:error )
            !$acc kernels loop gang(32), vector(16)
            do j=1,m-2
                !$acc loop gang(16), vector(32)
                do i=1,n-2
                    Anew(i,j) = 0.25_fp_kind * ( A(i+1,j  ) + A(i-1,j  ) + &
                        A(i  ,j-1) + A(i  ,j+1) )
                    error = max( error, abs(Anew(i,j)-A(i,j)) )
                end do
                !$acc end loop
            end do
            !$acc end kernels
            !$omp end parallel do

            if(mod(iter,100).eq.0 ) write(*,'(i5,f10.6)'), iter, error
            iter = iter +1


            !$omp parallel do shared(m, n, Anew, A)
            !$acc kernels loop
            do j=1,m-2
                !$acc loop gang(16), vector(32)
                do i=1,n-2
                    A(i,j) = Anew(i,j)
                end do
                !$acc end loop
            end do
            !$acc end kernels
            !$omp end parallel do

        end do
        !$acc end data

        call cpu_time(stop_time)
        write(*,'(a,f10.3,a)')  ' completed in ', stop_time-start_time, ' seconds'

        deallocate (A,Anew,y0)
#else
        logical, intent(inout):: flag
#endif

    end subroutine test_acc




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

    !
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

        implicit none

        complex(sp), allocatable :: a(:,:),b(:,:)
        complex(sp), device, allocatable :: a_d(:,:), b_d(:,:)
        integer i,j,it
        integer :: nerrors
        real :: minr, mini
        integer :: n=600
        integer :: plan, planType
        do it=1,300
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
            call cufftExec(plan,planType,a_d,b_d,CUFFT_FORWARD)

            ! Copy results back to host
            b = b_d
            !print *, "Forward B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

            call cufftExec(plan,planType,b_d,b_d,CUFFT_INVERSE)

            ! Copy results back to host
            b = b_d
            !print *, "Inverse B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

            ! Scale
            b = b / (n*n)
            !print *, "Scaled B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

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
                print *, "Test PASSED"
            endif

            !release memory on the host and device
            deallocate(a, b, a_d, b_d)

            ! Destroy the plan
            call cufftDestroy(plan)
        end do
    end subroutine cufftTest



    !> \brief fft  forward Fourier transform
    !!
    subroutine fft_pgi_cuda_test( self )
        ! use simple_fftshifter
        use openacc
        use cudafor
        !! IN-PLACE 3D FFT using CUDA's cufft
        use simple_cufft
        implicit none
        type(image), intent(inout)      :: self
        real(sp),    allocatable         :: rinput(:,:,:)
        !real(sp),    allocatable, device :: rinput_d(:,:,:)
        !complex(sp), allocatable, device :: coutput_d(:,:,:)
        complex(sp), allocatable, device :: inplace_d(:,:,:) !< inplace real and complex FFT data
        complex(sp), allocatable         :: coutput(:,:,:)
        integer      :: plan, planType, ierr
        integer      :: i,j,k,h,l,n,istat
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        integer(timer_int_kind) :: t1


        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        if(istat.ne.0)then
            call simple_cuda_stop("In simple_image::fft start ",__FILENAME__,__LINE__)
        endif


        verbose=.false.
        if(verbose) t1=tic()
        ldim = self%get_ldim()
        cdim = self%get_array_shape()
        nscale =  real(product(ldim))
        !   lims = self%loop_lims(3)
        VerbosePrint "In simple_cufft::fft_pgi  PGI test ldim ",ldim
        VerbosePrint "In simple_cufft::fft_pgi  PGI test cdim ", cdim
        VerbosePrint "In simple_cufft::fft_pgi  PGI test nc ", nscale
        !    VerbosePrint "In simple_cufft::fft_pgi  PGI test loop lims ", lims
        ierr=0
        ! allocate arrays on the host
        !  if(allocated(coutput))deallocate(coutput)
        allocate(coutput(cdim(1),cdim(2),cdim(3)),source=cmplx(0.0),stat=istat)
        if(istat /= 0)call allocchk("In simple_image::fft coutput",istat)
        !    if(allocated(rinput))deallocate(rinput)
        !  allocate( rinput(ldim(1), ldim(2), ldim(3)),source=0.,stat=istat)
        !    if(istat.ne.0)call allocchk("In simple_image::pgi_fft rinput",istat)

        allocate( rinput(ldim(1), ldim(2), ldim(3)),stat=istat)
        rinput = self%get_rmat()
        if(istat /= 0)call allocchk("In simple_image::fft rinput ",istat)

        ! call gnufor_image(rinput(:ldim(1),:ldim(2),1),palette='gray')
        !coutput(:ldim(1), :ldim(2), :ldim(3))=cmplx(rinput(:ldim(1), :ldim(2), :ldim(3)),0.)

        !call fftshift(rinput)
        !VerbosePrint "In simple_cufft::fft_pgi  PGI test ", shape(coutput), shape(rinput)
        ! istat=cudaDeviceSynchronize()
        ! istat=istat+cudaThreadSynchronize()
        ! if(istat.ne.0)then
        !     call simple_cuda_stop("In simple_image::fft sync ",__FILENAME__,__LINE__)
        ! endif
        !        call fftshift(rinput)
        VerbosePrint "In simple_cufft::fft_pgi  PGI start "
        ! if(.not.allocated(coutput)) call simple_stop(' coutput failed ')
        ! if(.not.allocated(rinput)) call simple_stop(' rinput failed ')
        ! allocate arrays on the device
        allocate(inplace_d(cdim(1),cdim(2),cdim(3)),source=cmplx(0.,0.),stat=istat)
        if(istat /= 0) call allocchk("In simple_image::fft inplace_d",stat)
        !if(allocated(rinput_d))deallocate( rinput_d)
        !allocate( rinput_d(ldim(1),ldim(2),ldim(3)),source=0.,stat=istat)
        !if(istat /= 0)call allocchk("In simple_image::fft rinput_d",istat)
        ! if(allocated(coutput_d))deallocate(coutput_d)
        !allocate(coutput_d(cdim(1),cdim(2),cdim(3)),stat=istat)
        !if(istat /= 0) call allocchk("In simple_image::fft rinput_d",stat)
        ! set planType to either single or double precision
        planType = CUFFT_C2C    !! (fp_kind == singlePrecision)
        ! VerbosePrint "In simple_cufft::fft_pgi  PGI device allocated"

        !copy arrays to device
        inplace_d = cmplx(rinput(1:ldim(1)-1:2,:ldim(2),:ldim(3)),rinput(2:ldim(1):2,:ldim(2),:ldim(3)))
        !rinput_d = rinput
        VerbosePrint "In simple_cufft::fft_pgi  PGI  arrays copied to device "


        ! Initialize the plan for real to complex transform
        call cufftPlan3D(plan,ldim(1),ldim(2),ldim(3), planType)
        VerbosePrint "In simple_cufft::fft_pgi  PGI cufftPlan completed"

        ! Execute  Forward transform in place
        !call cufftExec(plan,planType,rinput_d,coutput_d,CUFFT_FORWARD )
        call cufftExec(plan,planType,inplace_d,inplace_d,CUFFT_FORWARD )
        VerbosePrint "In simple_cufft::fft_pgi  PGI cufftexec completed "

        ! now scale the values so that a ifft of the output yields the
        ! original image back, rather than a scaled version
        !$acc kernels
        !        coutput_d = coutput_d / nscale
        inplace_d = inplace_d / nscale
        !$acc end kernels
        VerbosePrint "In simple_cufft::fft_pgi  PGI acc kernel completed"

        ! Copy results back to host
        coutput = inplace_d
        VerbosePrint "In simple_cufft::fft_pgi  PGI copied results back to host"
        ! istat=cudaThreadSynchronize()
        ! istat=istat+cudaDeviceSynchronize()
        ! if(istat.ne.0)then
        !     call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
        ! endif
        ! VerbosePrint "In simple_cufft::fft_pgi  PGI synchronized "
        ! call ifftshift(coutput,lims)
        ! VerbosePrint " COUTPUT ", real(coutput(1:3,1,1)), cabs(coutput(1:3,1,1))
        VerbosePrint "In simple_cufft::fft_pgi  set cmat ", shape(coutput), LBOUND(coutput), UBOUND(coutput)
        ! call self%set_cmat(coutput(:ldim(1),:ldim(2),:ldim(3)))
        ! VerbosePrint "In simple_cufft::fft_pgi  PGI destroyed "

        !      call self%set_cmat( coutput )

        call self%set_ft(.true.)
        call self%set_cmat(coutput)
        ! release memory on the host and on the device
        !  if(allocated(rinput_d))then
        !      deallocate(rinput_d)!,stat=istat)
        !if(istat.ne.0)call allocchk("In simple_image::fft dealloc device coutput_d ",istat)
        !  endif
        !  if(allocated(coutput_d))then
        istat=0
        ! deallocate(coutput_d, rinput_d,stat=istat)
        deallocate(inplace_d,stat=istat)
        VerbosePrint "In simple_cufft::fft_pgi Released memory on the device PGI  ",istat
        if(istat /= 0)call allocchk("In simple_image::fft dealloc device coutput_d ",istat)
        !   endif



        ! call self%set_ft(.true.)
        ! call self%set_cmat(coutput(:ldim(1), :ldim(2), :ldim(3)))
        ! ! omp parallel do collapse(3) schedule(static) default(shared) private(h,k,l,phys) proc_bind(close)
        ! do h=lims(1,1),lims(1,2)
        !     do l=lims(2,1),lims(2,2)
        !         do k=lims(3,1),lims(3,2)
        !             phys = self%comp_addr_phys([h,l,k])
        !             ! if (h > 0) then
        !             !     phys(1) = h + 1
        !             !     phys(2) = l + 1 + MERGE(ldim(2),0,l  < 0)
        !             !     phys(3) = k + 1 + MERGE(ldim(3),0,k  < 0)
        !             ! else
        !             !     phys(1) = -h + 1
        !             !     phys(2) = -l + 1 + MERGE(ldim(2),0,-l  < 0)
        !             !     phys(3) = -k + 1 + MERGE(ldim(3),0,-k  < 0)
        !             ! end if
        !             call self%set_cmat_at_ind(phys(1),phys(2),phys(3),coutput(h-lims(1,1)+1,l-lims(2,1)+1,k-lims(3,1)+1))
        !         end do
        !     end do
        ! end do
        ! ! omp end parallel do
        !   call gnufor_image(real(coutput(:,:,1)),text='gpo')
        !   call gnufor_image(aimag(coutput(:,:,1)),text='gpo')
        !   call gnufor_image(cabs(coutput(:,:,1)),text='gpo')
        !   call gnufor_image(atan2(real(coutput(:,:,1)),aimag(coutput(:,:,1))),text='gpo')

        VerbosePrint "In simple_cufft::fft_pgi  deallocating host arrays"
        ! Remove host storage
        ! if(allocated(rinput))then
        !     VerbosePrint "In simple_cufft::fft_pgi  deallocating rinput"
        !     deallocate(rinput,stat=istat)
        !     if(istat /= 0)call allocchk("In simple_image::fft dealloc host rinput",istat)
        ! end if
        ! if(allocated(rinput).and.allocated(coutput))then
        !     VerbosePrint "In simple_cufft::fft_pgi  deallocating coutput and rinput"
        !     deallocate(rinput, coutput,stat=istat)
        !     if(istat /= 0)call allocchk("In simple_image::fft dealloc host arrays",istat)
        ! end if
        ! if(allocated(coutput))then

        if(allocated(coutput))deallocate(coutput,stat=istat)
        VerbosePrint "In simple_cufft::fft_pgi  deallocated coutput", istat
        if(istat /= 0) call allocchk("In simple_image::fft dealloc host coutput",istat)
        ! end if
        ! if(allocated(rinput))then
        if(allocated(rinput))deallocate(rinput,stat=istat)
        VerbosePrint "In simple_cufft::fft_pgi  deallocated rinput", istat
        if(istat /= 0)call allocchk("In simple_image::fft dealloc host rinput",istat)
        ! end if
        VerbosePrint "In simple_cufft::fft_pgi deallocating done "

        ! Destroy the plans
        call cufftDestroy(plan)
        !if (ierr /= 0) print *, " fft_pgi_cuda_test error: ", ierr
        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        if(istat.ne.0)then
            call simple_cuda_stop("In simple_image::fft end ",__FILENAME__,__LINE__)
        endif
        VerbosePrint "In simple_cufft::fft_pgi finished "
        !return
    end subroutine fft_pgi_cuda_test

    subroutine fft_pgi_stream_test (self)

        use simple_cufft
        use cudafor
        implicit none
        type(image), intent(inout)  :: self
        complex(fp_kind), allocatable, dimension(:,:,:),pinned :: A,B,C
        complex(fp_kind), allocatable, dimension(:,:,:),device :: A_d,B_d
        integer, parameter :: num_streams=4
        integer:: nx, ny, nomega, ifr, i,j, istat,stream_index
        integer(8):: clock_start,clock_end,clock_rate
        integer(kind=cuda_stream_kind) :: stream(num_streams)
        type(c_ptr):: plan
        real:: elapsed_time
        real(fp_kind):: scale

        nx=512; ny=512;  nomega=196
        scale=1./real(nx*ny,fp_kind)

        ! Initialize FFT plan
        call cufftPlan2d(plan,nx,ny,CUFFT_C2C)

        ! Create streams
        do i = 1,num_streams
            istat= cudaStreamCreate(stream(i))
        end do

        call SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate

        ! Allocate arrays on CPU and GPU
        allocate(A(nx,ny,nomega),B(nx,ny,nomega),C(nx,ny,nomega))
        allocate(A_d(nx,ny,num_streams),B_d(nx,ny,num_streams))

        ! Initialize arrays on CPU
        A=cmplx(1.,1.,fp_kind); B=cmplx(1.,1.,fp_kind); C=cmplx(0.,0.,fp_kind)

        ! Measure only the transfer time
        istat=cudaThreadSynchronize()

        print *,"I/O only"
        call SYSTEM_CLOCK(COUNT=clock_start) ! Start timing

        do ifr=1,nomega
            istat= cudaMemcpy(A_d(1,1,1),A(1,1,ifr),nx*ny)
            istat= cudaMemcpy(B_d(1,1,1),B(1,1,ifr),nx*ny)
            istat= cudaMemcpy(C(1,1,ifr),A_d(1,1,1),nx*ny)
        end do

        istat=cudaThreadSynchronize()
        call SYSTEM_CLOCK(COUNT=clock_end) ! End timing
        elapsed_time=REAL(clock_end-clock_start)/REAL(clock_rate)
        print *,"Elapsed time :",elapsed_time

        ! Measure the transfer time H2D, FFT , IFFT and transfer time D2H

        print '(/a)',"Single stream  loop"
        istat=cudaThreadSynchronize()
        call SYSTEM_CLOCK(COUNT=clock_start) ! Start timing
        stream_index = 1
        call cufftSetStream(plan,stream(stream_index))
        do ifr=1,nomega
            istat= cudaMemcpy(A_d(1,1,stream_index),A(1,1,ifr),nx*ny)
            istat= cudaMemcpy(B_d(1,1,stream_index),B(1,1,ifr),nx*ny)
            call cufftExecC2C(plan ,A_d(1,1,stream_index),&
                A_d(1,1,stream_index),CUFFT_FORWARD)
            call cufftExecC2C(plan ,B_d(1,1,stream_index),&
                B_d(1,1,stream_index),CUFFT_FORWARD)

            ! Convolution and scaling of the  arrays
            !$cuf kernel do(2) <<<*,(16,16),stream=stream(stream_index)>>>
            do j=1,ny
                do i=1,nx
                    B_d(i,j,stream_index)= A_d(i,j,stream_index)*&
                        B_d(i,j,stream_index)*scale
                end do
            end do

            call cufftExecC2C(plan ,B_d(1,1,stream_index),&
                B_d(1,1,stream_index),CUFFT_INVERSE)
            istat=cudaMemcpy( C(1,1,ifr),B_d(1,1,stream_index),nx*ny)
        end do

        istat=cudaThreadSynchronize()
        call SYSTEM_CLOCK(COUNT=clock_end) ! End timing
        elapsed_time=REAL(clock_end-clock_start)/REAL(clock_rate)
        print *,"Elapsed time :",elapsed_time

        ! Overlap I/O and compute using multiple streams and async copies
        print '(/a)',"Do loop with multiple streams"
        call SYSTEM_CLOCK(COUNT=clock_start) ! Start timing

        do ifr=1,nomega

            ! assign a stream for the current plan
            stream_index = mod(ifr,num_streams)+1

            ! Set the stream used by CUFFT
            call cufftSetStream(plan,stream(stream_index))

            ! Send A to GPU
            istat= cudaMemcpyAsync(A_d(1,1,stream_index),A(1,1,ifr),&
                nx*ny, stream(stream_index))

            ! Execute forward FFTs on GPU
            call cufftExecC2C(plan ,A_d(1,1,stream_index),&
                A_d(1,1,stream_index),CUFFT_FORWARD)

            ! Send B to GPU
            istat= cudaMemcpyAsync(B_d(1,1,stream_index), &
                B(1,1,ifr),nx*ny, stream(stream_index))

            ! Execute forward FFTs on GPU
            call cufftExecC2C(plan ,B_d(1,1,stream_index),&
                B_d(1,1,stream_index),CUFFT_FORWARD)

            ! Convolution and scaling of the  arrays
            !$cuf kernel do(2) <<<*,(16,16),stream=stream(stream_index)>>>
            do j=1,ny
                do i=1,nx
                    B_d(i,j,stream_index)= A_d(i,j,stream_index)* &
                        B_d(i,j,stream_index)*scale
                end do
            end do

            ! Execute inverse FFTs on GPU
            call cufftExecC2C(plan ,B_d(1,1,stream_index), &
                B_d(1,1,stream_index),CUFFT_INVERSE)

            ! Copy results back
            istat=cudaMemcpyAsync( C(1,1,ifr),B_d(1,1,stream_index), &
                nx*ny, stream=stream(stream_index))

        end do

        istat=cudaThreadSynchronize()
        call SYSTEM_CLOCK(COUNT=clock_end) ! Start timing
        elapsed_time=REAL(clock_end-clock_start)/REAL(clock_rate)
        print *,"Elapsed time :",elapsed_time

        if (istat .eq. 0) then
            print *,"Test Passed"
        else
            print *,"Test Failed"
        endif

        deallocate(A,B,C); deallocate(A_d,B_d)
        call cufftDestroy(plan)

    end subroutine fft_pgi_stream_test

    subroutine ifft_pgi_cuda_test( self )
        !#ifndef PGI
        !        class(image), intent(inout) :: self
        !        print *,'CUDA ifft cannot be performed with current compiler'
        !#else
        !! IN-PLACE 3D FFT using CUDA's cufft
        use simple_cufft
        use simple_fftshifter
        use openacc
        use cudafor

        type(image), intent(inout)       :: self
        real(dp),    allocatable         :: routput(:,:,:),rtmp(:,:,:)
        real(dp),    allocatable, device :: routput_d(:,:,:)
        complex(dp), allocatable, device :: cinput_d(:,:,:)
        complex(dp), allocatable         :: cinput(:,:,:)
        complex(sp), allocatable         :: ctmp(:,:,:)
        ! type(image)  :: c1,c2,c3
        integer      :: plan, planType
        integer      :: i,j,k,h,l,n,istat,mk,mh,ml,fsize
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        complex(dp)      :: comp, cswap
        integer(timer_int_kind) :: t1

        verbose=.false.
        !   if(verbose) t1=tic()
        ldim = self%get_ldim()
        cdim = self%get_array_shape()
        fsize = self%get_filtsz()
        lims = self%loop_lims(2)
        VerbosePrint "In simple_cufft::ifft_pgi  PGI cmat lims", lims
        nscale =  real(product(ldim))
        VerbosePrint "In simple_cufft::ifft_pgi  PGI LDIM/CDIM", ldim, cdim

        ! allocate arrays on the host
        if(allocated(routput))deallocate(routput)
        allocate(routput(ldim(1), ldim(2), ldim(3)),source=0._dp,stat=istat)
        if(istat.ne.0)call allocchk("In simple_image::ifft routput",istat)
        if(allocated(cinput))deallocate( cinput )
        allocate(cinput((ldim(1)/2)+1, ldim(2), ldim(3)),source=dcmplx(0._dp,0._dp),stat=istat)
        if(istat.ne.0)call allocchk("In simple_image::ifft cinput",istat)
        allocate(ctmp(cdim(1),cdim(2),cdim(3)),source=cmplx(0.,0.),stat=istat)
        if(istat.ne.0)call allocchk("In simple_image::ifft ctmp",istat)

        ! allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        !rtmp = self%get_rmat()
        ctmp = self%get_cmatfull()
        !! VerbosePrint "In simple_cufft::ifft_pgi  PGI routput ", shape(routput), LBOUND(routput), UBOUND(routput)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI cinput ", shape(cinput), LBOUND(cinput), UBOUND(cinput)
        VerbosePrint "In simple_cufft::ifft_pgi  PGI ctmp ", shape(ctmp), LBOUND(ctmp), UBOUND(ctmp)
        !call ifftshift(ctmp)
        !call gnufor_image(rtmp(:,:,1), palette='gray')
        !call gnufor_image(real(ctmp(:,:,1)), palette='gray')
        !call gnufor_image(aimag(ctmp(:,:,1)), palette='gray')

        cinput=dcmplx(ctmp(:(ldim(1)/2)+1,:ldim(2),:ldim(3)))
        print *,"shape cinput", shape(cinput)
        !cinput(:cdim(1),:ldim(2),:ldim(3)) = cmplx(rtmp(1:ldim(1):2,:ldim(2),:ldim(3)),rtmp(2:ldim(1):2,:ldim(2),:ldim(3)))
        ! cinput(cdim(1)+2:ldim(1),1:ldim(2)-1,:ldim(3)) =  cmplx(rtmp(:cdim(1)-2,2:ldim(2),:ldim(3))
        !  cinput(cdim(1)+1:ldim(1),:cdim(2),:cdim(3)) = ctmp(cdim(1):-1:2, :cdim(2), :cdim(3))
        !conjg(ctmp(2:cdim(1), :cdim(2),:cdim(3)))
        ! !omp parallel do collapse(3) default(shared) private(cswap,i,j,k)&
        ! !omp schedule(static) proc_bind(close)
        !         do i=1,ldim(1)/2
        !             do j=1,ldim(2)/2
        !                 do k=1,ldim(3)/2
        !                     !(1)
        !                     cswap = cinput(i,j,k)
        !                     cinput(i,j,k) = cinput(ldim(1)/2+i,ldim(2)/2+j,ldim(3)/2+k)
        !                     cinput(ldim(1)/2+i,ldim(2)/2+j,ldim(3)/2+k) = cswap
        !                     !(2)
        !                     cswap = cinput(i,ldim(2)/2+j,ldim(3)/2+k)
        !                     cinput(i,ldim(2)/2+j,ldim(3)/2+k) = cinput(ldim(1)/2+i,j,k)
        !                     cinput(ldim(1)/2+i,j,k) = cswap
        !                     !(3)
        !                     cswap = cinput(ldim(1)/2+i,j,ldim(3)/2+k)
        !                     cinput(ldim(1)/2+i,j,ldim(3)/2+k) = cinput(i,ldim(2)/2+j,k)
        !                     cinput(i,ldim(2)/2+j,k) = cswap
        !                     !(4)
        !                     cswap = cinput(i,j,ldim(3)/2+k)
        !                     cinput(i,j,ldim(3)/2+k) = cinput(ldim(1)/2+i,ldim(2)/2+j,k)
        !                     cinput(ldim(1)/2+i,ldim(2)/2+j,k) = cswap
        !                 end do
        !             end do
        !         end do
        !         !omp end parallel do

        ! call ifftshift(cinput)
        !     call gnufor_image((real(cinput(:,:,1),sp)),palette='gray')
        !      call gnufor_image(real(aimag(cinput(:,:,1)),sp),palette='gray')

        ! !omp parallel do collapse(3) schedule(static) default(shared) private(h,k,l,phys,comp) proc_bind(close)
        ! do h=lims(1,1),lims(1,2)
        !     do k=lims(2,1),lims(2,2)
        !         do l=lims(3,1),lims(3,2)
        !             !                        phys = self%comp_addr_phys([h,k,l])
        !             if (h > 0) then
        !                 phys(1) = h + 1
        !                 phys(2) = k + 1 + MERGE(ldim(2),0,k  < 0)
        !                 phys(3) = l + 1 + MERGE(ldim(3),0,l  < 0)
        !             else
        !                 phys(1) = -h + 1
        !                 phys(2) = -k + 1 + MERGE(ldim(2),0,-k  < 0)
        !                 phys(3) = -l + 1 + MERGE(ldim(3),0,-l  < 0)
        !             end if
        !             comp = self%get_cmat_at_ind(phys(1),phys(2),phys(3))
        !             ! if(h<0)comp=conjg(comp)
        !             cinput(phys(1),phys(2),phys(3))=comp
        !         end do
        !     end do
        ! end do
        ! !omp end parallel do
        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        if(istat.ne.0)then
            call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        endif
        allocate(routput_d(ldim(1), ldim(2), ldim(3)),stat=istat)
        if(istat.ne.0) call allocchk("In simple_image::ifft routput_d",istat)

        allocate(cinput_d((ldim(1)/2)+1, ldim(2), ldim(3)),stat=istat)
        if(istat.ne.0) call allocchk("In simple_image::ifft cinput_d",istat)
        ! set planType to either single or double precision
        planType = CUFFT_Z2D    !! (fp_kind == singlePrecision)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI device allocated ", size(cinput_d), size(routput)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI host -> device ", shape(cinput), LBOUND(cinput), UBOUND(cinput)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI host -> device ", shape(cinput_d), LBOUND(cinput_d), UBOUND(cinput_d)
        ! copy input to device
        cinput_d = cinput
        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI "
        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI host -> device ", shape(cinput_d), lbound(cinput_d), ubound(cinput_d)


        ! Initialize the plan for complex to real transform
        call cufftPlan3d(plan,ldim(1),ldim(2),ldim(3),planType)
        VerbosePrint "In simple_cufft::ifft_pgi  PGI "

        ! Execute  Backward transform in place
        call cufftExec(plan,planType,cinput_d,routput_d,CUFFT_INVERSE )
        VerbosePrint "In simple_cufft::ifft_pgi  PGI exec completed "

        ! ! Copy results back to host
        ! istat=cudaThreadSynchronize()
        ! istat=istat+cudaDeviceSynchronize()
        ! if(istat.ne.0)then
        !     call simple_cuda_stop("In simple_image::ifft post fft sync ",__FILENAME__,__LINE__)
        ! endif
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI copy results back to host"
        ! Copy results back to host
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI device -> host ", size(cinput_d), size(routput)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI host shape ", shape(routput), LBOUND(routput), UBOUND(routput)
        routput = routput_d
        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI device -> host ", size(cinput_d), size(routput)
        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI deallocate arrays on device"

        ! release memory on the host and on the device
        ! if(allocated(routput_d))then
        deallocate(routput_d)!,stat=istat)
        !if(istat /= 0)call allocchk("In simple_image::ifft dealloc routput_d",istat)
        !  end if
        !  if(allocated(cinput_d))then
        deallocate(cinput_d)!,stat=istat)
        !if(istat /= 0)call allocchk("In simple_image::ifft dealloc cinput_d",istat)
        !  end if

        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI destroy cufft plans"
        ! Destroy the plans
        call cufftDestroy(plan)
        istat=cudaThreadSynchronize()
        istat=istat+cudaDeviceSynchronize()
        if(istat /= 0) then
            call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        endif
        VerbosePrint "In simple_cufft::ifft_pgi  PGIC UFFT vs FFTW3"
        !call fftshift(routput)
        call gnufor_image((real(routput(:,:,1))),text='go',palette='hot')
        call self%set_rmat( real(routput(:ldim(1),:ldim(2),:ldim(3))) )
        call self%set_ft(.false.)
        call self%shift_phorig
        ! Cleanup
        if(allocated(rtmp))then
            deallocate( rtmp, stat=istat )
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating rtmp",istat)
        end if
        if(allocated(ctmp))then
            deallocate( ctmp, stat=istat )
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating ctmp",istat)
        end if
        if(allocated(routput))then
            ! deallocate(routput,stat=istat)
            ! if(istat /= 0)call allocchk("In simple_image::ifft deallocatiing routput",istat)
        end if
        if(allocated(cinput))then
            deallocate(cinput,stat=istat)
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating cinput",istat)
        end if

        !  #endif
        return
    end subroutine ifft_pgi_cuda_test


    subroutine ifft_pgi_cuda_test2( self )

        use simple_cufft
        use simple_fftshifter
        use openacc
        use cudafor

        type(image), intent(inout)       :: self
        ! real(sp),    allocatable         :: routput(:,:,:)
        !complex(sp),    allocatable, device :: routput_d(:,:,:)
        !  complex(sp), allocatable, device :: cinput_d(:,:,:)
        complex(sp), allocatable, device :: inplace_d(:,:,:)
        ! complex(sp), allocatable         :: cinput(:,:,:), coutput(:,:,:)
        complex(sp), allocatable         :: cinout(:,:,:),ctmp(:,:,:)
        ! type(image)  :: c1,c2,c3
        integer      :: plan, planType, xdim,ydim
        integer      :: i,j,k,h,l,n,istat,mk,mh,ml,fsize
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        complex(sp)      :: comp, cswap
        integer(timer_int_kind) :: t1
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        verbose=.true.
        if(verbose)then
            t1=tic()
            print *,"In ifft_pgi_cuda_test2 "
        end if
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)

        ldim = self%get_ldim()
        cdim = self%get_array_shape()
        cdim(1)=cdim(1)
        fsize = self%get_filtsz()
        lims = self%loop_lims(3)
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI cmat lims, fsize", lims, fsize
        nscale =  real(product(ldim))
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI LDIM/nscale", ldim, nscale
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI CDIM       ", cdim
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)

        ! allocate arrays on the host
        ! if(allocated(coutput))deallocate(coutput)
        ! allocate(coutput(ldim(1), ldim(2), ldim(3)),source=cmplx(0.,0.),stat=istat)
        ! if(istat.ne.0)call allocchk("In simple_image::ifft routput",istat)
        ! if(allocated(cinout))deallocate( cinout )
        ! print *,"Half plus 1", ceiling(real(ldim(1))/2.)+1
        !allocate(cinout(ldim(1), ldim(2), ldim(3)),stat=istat)
        ! if(istat.ne.0)call allocchk("In simple_image::ifft cinout",istat)
        ! allocate(ctmp(lims(1,1):lims(1,2),lims(2,1):lims(2,2),lims(3,1):lims(3,2)),source=cmplx(0.,0.),stat=istat)
        allocate(ctmp(cdim(1), ldim(2), ldim(3)),source=cmplx(0.,0.),stat=istat)
        if(istat /= 0)call allocchk("In simple_image::ifft ctmp",istat)

        ! allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        !  rtmp = self%get_rmat()
        ctmp = self%get_cmatfull()
        !call ifftshift(ctmp)
        call gnufor_image(real(ctmp(:,:,1)), palette='gray')
        call gnufor_image(aimag(ctmp(:,:,1)), palette='gray')
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI ctmp ",size(ctmp),lbound(ctmp), ubound(ctmp)
        k=1

        if(is_even(ldim(1)))then
            xdim = ldim(1)/2
            ydim = ldim(2)/2
        else
            xdim = (ldim(1)-1)/2
            ydim = (ldim(2)-1)/2
        endif
        !allocate(cinout(-xdim:xdim,-ydim:ydim,1:1),stat=istat)
        allocate(cinout(ldim(1), ldim(2), ldim(3)),stat=istat)
        cinout = cmplx(0.,0.)
        !$omp parallel do collapse(2) default(shared) private(h,k,phys)
        do h=-xdim,xdim
            do k=-ydim,ydim
                phys = self%comp_addr_phys([h,k,0])
                if(h<0)then
                    cinout(h+xdim+1,k+ydim+1,1) = ctmp(phys(1),phys(2),1)
                else
                    cinout(h+xdim+1,k+ydim+1,1) = conjg(ctmp(phys(1),phys(2),1))
                endif

            end do
        end do
        !$omp end parallel do
        ! call fftshift(cinout)
        call gnufor_image(real(cinout(:,:,1)), palette='gray')
        call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        !    read(*,*)
        !         deallocate( cinout )
        !         allocate(cinout(ldim(1), ldim(2), ldim(3)),stat=istat)

        !         !$omp parallel do collapse(2) default(shared) private(i,j)
        !   do i=1,cdim(1)
        !          do j=1,floor(real(ldim(2))/2.)
        !              k=1
        ! !             do k=1,floor(real(ldim(3))/2.)+1
        !              cinout(i,j,k) =( ctmp(i,j,k))

        !                !if(i/=cdim(1))cinout(ldim(1)-i+1,j,k) =  conjg(ctmp(i,ldim(2)-j+1,k))

        !                !  cinout(i,ldim(2)-j+1,k) = ctmp(i,ldim(2)-j+1,k)
        !             !     if(i/=1)cinout(ldim(1)-i+2,ldim(2)-j+1,k) =  ( ctmp(i,j,k)) !ldim(2)-j+1,k) )
        !                  ! if(i/=cdim(1))cinout(ldim(1)-i+1,ldim(2)-j+1,k) =  ( ctmp(i,ldim(2)-j+1,k) )
        !         !             routput((i-1)*2+1,j,k) = real(cinout(i,j,k))
        !         !             routput((i-1)*2+2,j,k) = aimag(cinout(i,j,k))
        !  !                end do
        !              end do
        !              do j=floor(real(ldim(2))/2.)+1,ldim(2)
        !              k=1
        ! !             do k=1,floor(real(ldim(3))/2.)+1
        !             ! cinout(i,j,k) = ( ctmp(i,j,k))

        !                !if(i/=cdim(1))cinout(ldim(1)-i+1,j,k) =  conjg(ctmp(i,ldim(2)-j+1,k))

        !                !  cinout(i,ldim(2)-j+1,k) = ctmp(i,ldim(2)-j+1,k)
        !              !    if(i/=1)cinout(ldim(1)-i+2,ldim(2)-j+1,k) =  ( ctmp(i,j,k)) !ldim(2)-j+1,k) )

        !         !             routput((i-1)*2+1,j,k) = real(cinout(i,j,k))
        !         !             routput((i-1)*2+2,j,k) = aimag(cinout(i,j,k))
        !  !                end do
        !                     end do
        !          end do
        !          !$omp end parallel do
        !         ! cinout(cdim(1)-1:ldim(1), cdim(2)-1:ldim(2), cdim(3)-1:ldim(3))=ctmp(1:cdim(1), 1:ceiling(cdim(2)/2), 1:ceiling(cdim(3)/2)
        !         ! cinout(cdim(1)-1:ldim(1), cdim(2)-1:ldim(2), cdim(3)-1:ldim(3))=ctmp(1:cdim(1), ceiling(cdim(2)/2)+1:cdim(2), ceiling(cdim(3)/2)+1:cdim(2)
        !         ! cinout(cdim(1)-1:ldim(1), cdim(2)-1:ldim(2), cdim(3)-1:ldim(3))=ctmp(1:cdim(1), 1:ceiling(cdim(2)/2), 1:ceiling(cdim(3)/2)


        !         VerbosePrint "In simple_cufft::ifft_pgi2  PGI ctmp ",size(ctmp), lbound(ctmp), ubound(ctmp)
        !         !call ifftshift(cinout)

        !         !call gnufor_image(cabs(ctmp(:cdim(1),:ldim(2),1)), text='go',palette='hot')
        !         call gnufor_image(real(cinout(:,:,1)), palette='gray')
        !         call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        !         VerbosePrint "In simple_cufft::ifft_pgi2  ctmp shape", shape(ctmp)
        !        ! cinout(:cdim(1), :ldim(2), :ldim(3)) = cmplx(ctmp(:cdim(1), :ldim(2), :ldim(3)))
        !        ! VerbosePrint "In simple_cufft::ifft_pgi2  cinout shape", shape(cinout)

        !         call fftshift(cinout)
        !         call gnufor_image(real(cinout(:,:,1)), palette='gray')
        !         call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        ! Synchronize
        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Synchronize completed"

        !allocate on device
        ! if(allocated(routput_d)) deallocate( routput_d )
        ! VerbosePrint "In simple_cufft::ifft_pgi2 device allocation"
        ! call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        allocate(inplace_d(ldim(1), ldim(2), ldim(3)))!,stat=istat)
        call simple_cuda_stop("In ifft TEST2 routput_d",__FILENAME__,__LINE__)
        !allocate(routput_d(cdim(1), ldim(2), ldim(3)))!,stat=istat)
        !call simple_cuda_stop("In ifft TEST2 routput_d",__FILENAME__,__LINE__)
        !if(istat.ne.0) call allocchk("In simple_image::ifft routput_d",istat)
        VerbosePrint "In simple_cufft::ifft_pgi2 routput_d allocation completed"
        ! if(allocated(cinout_d)) deallocate( cinout_d )
        ! allocate(cinout_d(cdim(1), cdim(2), cdim(3)))!,stat=istat)
        ! call simple_cuda_stop("In ifft TEST2 cinout_d",__FILENAME__,__LINE__)
        !if(istat.ne.0) call allocchk("In simple_image::ifft cinout_d",istat)
        VerbosePrint "In simple_cufft::ifft_pgi2 device allocation completed"

        ! set planType to either single or double precision
        planType = CUFFT_C2C
        VerbosePrint "In simple_cufft::ifft_pgi2 set planType completed"

        ! copy input to device
        ! routput_d = cinout
        inplace_d = cinout
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 copy input to device completed"

        ! Initialize the plan for complex to real transform
        call cufftPlan3d(plan,ldim(1),ldim(2),ldim(3),planType)
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Initialize the plan completed"

        ! Execute  Backward transform in place
        call cufftExec(plan,planType,inplace_d,inplace_d,CUFFT_INVERSE )
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Execute  Backward transform completed"

        ! Copy results back to host
        cinout = inplace_d

        VerbosePrint "In simple_cufft::ifft_pgi2  PGI ctmp ",size(ctmp), ubound(ctmp),lbound(ctmp)
        ! !$omp parallel do default(shared) private(i,j,k)
        ! !$acc kernels loop
        ! do j=1,cdim(2)
        !     do k=1,cdim(3)
        !         !$acc loop gang(16), vector(32)
        !         do i=1,cdim(1)-1
        !             routput((i-1)*2+1,j,k) = real(cinout(i,j,k))
        !             routput((i-1)*2+2,j,k) = aimag(cinout(i,j,k))
        !         end do
        !         !$acc end loop
        !     end do
        ! end do
        ! !$acc end kernels
        ! !$omp end parallel do
        ! !$omp parallel do default(shared) private(i)
        !  do i=1,cdim(1)-1
        !      routput((i-1)*2+1, :ldim(2), :ldim(3)) = real(coutput(:cdim(1),:cdim(2),:cdim(3)))
        !      routput((i-1)*2+2, :ldim(2), :ldim(3)) = aimag(coutput(:cdim(1),:cdim(2),:cdim(3)))
        ! end do
        ! !$omp end parallel do
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Copy results back to host completed"

        ! release memory on the host and on the device
        !  if(allocated(routput_d))then
        !      call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        deallocate(inplace_d)!,stat=istat)
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        !if(istat /= 0)call allocchk("In simple_image::ifft dealloc routput_d",istat)
        !  end if
        ! if(allocated(cinout_d))then
        !     call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        !     deallocate(cinout_d)!,stat=istat)
        !     call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        !     !if(istat /= 0)call allocchk("In simple_image::ifft dealloc cinout_d",istat)
        ! end if
        !        call fftshift(ctmp)
        call gnufor_image(real(cinout(:,:,1)), palette='gray')
        call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        call gnufor_image(cabs(cinout(:,:,1)), palette='gray')
        !call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        !call fftshift(routput)
        !        call self%set_cmat( cinout )
        call self%set_ft(.false.)
        call self%set_rmat( real(cinout) )
        call self%shift_phorig
        ! Cleanup
        ! if(allocated(rtmp))then
        !     deallocate( rtmp, stat=istat )
        !     if(istat /= 0)call allocchk("In simple_image::ifft deallocating rtmp",istat)
        ! end if
        if(allocated(ctmp))then
            deallocate( ctmp, stat=istat )
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating ctmp",istat)
        end if
        ! if(allocated(coutput))then
        !     deallocate(coutput,stat=istat)
        !     if(istat /= 0)call allocchk("In simple_image::ifft deallocatiing coutput",istat)
        ! end if
        if(allocated(cinout))then
            deallocate(cinout,stat=istat)
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating cinout",istat)
        end if
        ! Destroy the plans
        call cufftDestroy(plan)
        !if (ierr /= 0) print *, " fft_pgi_cuda_test error: ", ierr
        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        if(istat.ne.0)then
            call simple_cuda_stop("In simple_image::ifft end ",__FILENAME__,__LINE__)
        endif
        call self%set_ft(.false.)
        VerbosePrint "In simple_cufft::ifft_pgi finished "
    end subroutine ifft_pgi_cuda_test2
    subroutine ifft_pgi_cuda_test3( self )

        use simple_cufft
        !  use simple_fftshifter
        !  use openacc
        !  use cudafor

        type(image), intent(inout)       :: self

        complex(sp), allocatable, device :: inplace_d(:,:,:)

        complex(sp), allocatable         :: cinout(:,:,:) !,ctmp(:,:,:)
        !   complex(sp), allocatable         :: coutput(:,:,:)
        integer      :: plan, planType, xdim,ydim
        integer      :: i,j,k,h,l,n,istat,mk,mh,ml,fsize
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        complex(sp)      :: comp, cswap
        integer(timer_int_kind) :: t1
        !   call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        verbose=.true.
        ! if(verbose)then
        !     t1=tic()
        !     print *,"In ifft_pgi_cuda_test2 "
        ! end if
        ! call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        call self%set_ft(.false.) ;print *,"In simple_image::ifft set ft"
        ldim = self%get_ldim()
        allocate(cinout(ldim(1),ldim(2),ldim(3)),inplace_d(ldim(1),ldim(2),ldim(3)))
        ! cdim = self%get_array_shape()
        ! cdim(1)=cdim(1)
        ! fsize = self%get_filtsz()
        ! lims = self%loop_lims(3)
        ! VerbosePrint "In simple_cufft::ifft_pgi2  PGI cmat lims, fsize", lims, fsize
        ! nscale =  real(product(ldim))
        ! VerbosePrint "In simple_cufft::ifft_pgi2  PGI LDIM/nscale", ldim, nscale
        ! VerbosePrint "In simple_cufft::ifft_pgi2  PGI CDIM       ", cdim
        !  call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)

        ! allocate(ctmp(cdim(1), ldim(2), ldim(3)),source=cmplx(0.,0.),stat=istat)
        ! if(istat /= 0)call allocchk("In simple_image::ifft ctmp",istat)

        ! allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        !  rtmp = self%get_rmat()
        ! ctmp = self%get_cmatfull()
        !call ifftshift(ctmp)
        !   call gnufor_image(real(ctmp(:,:,1)), palette='gray')
        !  call gnufor_image(aimag(ctmp(:,:,1)), palette='gray')
        ! VerbosePrint "In simple_cufft::ifft_pgi2  PGI ctmp ",size(ctmp),lbound(ctmp), ubound(ctmp)
        k=1

        if(is_even(ldim(1)))then
            xdim = ldim(1)/2
            ydim = ldim(2)/2
        else
            xdim = (ldim(1)-1)/2
            ydim = (ldim(2)-1)/2
        endif
        !allocate(cinout(-xdim:xdim,-ydim:ydim,1:1),stat=istat)

        cinout = cmplx(0.,0.)
        !$omp parallel do collapse(2) default(shared) private(h,k,phys,comp)
        do h=-xdim,xdim
            do k=-ydim,ydim
                phys = self%comp_addr_phys([h,k,0])
                comp= self%get_cmat_at_ind(phys(1),phys(2),1)
                if(h>=0)then
                    cinout(h+xdim+1,k+ydim+1,1) = comp
                else
                    cinout(h+xdim+1,k+ydim+1,1) = conjg(comp)
                endif

            end do
        end do
        !$omp end parallel do
        ! call fftshift(cinout)
        !call gnufor_image(real(cinout(:,:,1)), palette='gray')
        !call gnufor_image(aimag(cinout(:,:,1)), palette='gray')

        ! Synchronize
        !  istat=cudaDeviceSynchronize()
        !  istat=istat+cudaThreadSynchronize()
        !  call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Synchronize completed"


        ! allocate(inplace_d(ldim(1), ldim(2), ldim(3)),stat=istat)
        ! call simple_cuda_stop("In ifft TEST2 routput_d",__FILENAME__,__LINE__)



        ! set planType to either single or double precision
        planType = CUFFT_C2C
        !VerbosePrint "In simple_cufft::ifft_pgi2 set planType completed"

        ! copy input to device
        ! routput_d = cinout
        inplace_d = cinout
        VerbosePrint "In simple_cufft::ifft_pgi2 copy input to device completed"

        ! Initialize the plan for complex to real transform
        call cufftPlan3d(plan,ldim(1),ldim(2),ldim(3),planType)
        VerbosePrint "In simple_cufft::ifft_pgi2 Initialize the plan completed"

        ! Execute  Backward transform in place
        call cufftExec(plan,planType,inplace_d,inplace_d,CUFFT_INVERSE )
        VerbosePrint "In simple_cufft::ifft_pgi2 Execute  Backward transform completed"

        ! Copy results back to host
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI inplace_d ",size(inplace_d), ubound(inplace_d),lbound(inplace_d)
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI cinout ",size(cinout), ubound(cinout),lbound(cinout)
        cinout = inplace_d
        VerbosePrint "In simple_cufft::ifft_pgi2 Copied results back to host completed"
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI cinout ",size(cinout), ubound(cinout),lbound(cinout)
        ! Save to image obj

        !call self%set_rmat( cabs(cinout) ) ;print *,"In simple_image::ifft cabs inout"
        !$omp parallel do collapse(2) default(shared) private(i,j)
        do i=1,ldim(1)
            do j=1,ldim(2)
                call self%set_rmat_at_ind(i, j, 1, cabs( cinout(i,j,1)) )
            end do
        end do
        !$omp end parallel do
        !  call self%shift_phorig ;print *,"In simple_image::ifft shift"
        ! Destroy the plans
        call cufftDestroy(plan)

        ! release memory on the host and on the device
        print *,"In simple_image::ifft deallocating cinout"

        ! deallocate(cinout ) ; print *," deallocated cinout"
        !   if(istat /= 0)call allocchk("In simple_image::ifft deallocating cinout",istat)
        !   deallocate(ctmp); print *," deallocated ctmp"
        deallocate(cinout, inplace_d ) ; print *," deallocated cinout"
        !deallocate( inplace_d) ; print *," deallocated inplace_d"

        ! Synchronize
        !  istat=cudaDeviceSynchronize()
        !  istat=istat+cudaThreadSynchronize()
        !   call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Synchronize completed"

        VerbosePrint "In simple_cufft::ifft_pgi finished "
    end subroutine ifft_pgi_cuda_test3

    subroutine test_pgi( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2
        integer              :: i, j, k, cnt, lfny, ldim(3)
        real                 :: input, msk, ave, sdev, var, med, xyz(3), pow
        real                 :: imcorr, recorr, corr, corr_lp, maxv, minv
        real                 :: smpd=2.
        logical              :: passed, test(6)


        verbose=.false.

        write(*,'(a)') '**info(unit_test, PGI CUFFT: test_pgi  '
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)
        !  call img_3%new(ldim, smpd)
        ! call img_4%new(ldim, smpd)

        call img%gauimg(10)
        call img%add_gauran(1.)

        img_2 = img

        if(doplot)write(*,'(a)') '**info(unit_test PGI CUFFT: images created '
        call fft_pgi_cuda_test(img_2)
        if(doplot)write(*,'(a)') '**info(simple_image PGI test): forward FFT completed '
        call ifft_pgi_cuda_test2(img_2)
        if(doplot)write(*,'(a)') '**info(simple_image PGI test): inverse FFT completed '
        write(*,'(a,ES20.10)') "**info(simple_image PGI test): L2 norm sum ", img.lonesum.img_2




        ! call img_3%kill
        ! call img_4%kill
        ! call img3d%kill

        !img_4 = img_2
        !    call img_2%ft2img('amp',img_4)
        !    if(doplot)call img_4%vis(geomorsphr=.false.)
        !  write(*,'(a)') '**info(simple_image PGI test): testing fftshifter fwd'
        !  write(*,'(a)') 'Press enter to test fftshift '
        !    read(*,*)
        !     call img_3%fwd_fftshift()
        !    img_4 = img_2
        if(doplot)call img_2%vis(geomorsphr=.false.)
        if(doplot)call img%vis(geomorsphr=.false.)
        !     write(*,'(a)') 'Press enter to test fftw3 '
        !    read(*,*)

        !     write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '
        !     call img_3%fwd_ft()
        !     call img_3%ft2img('amp',img_4)
        !     if(doplot)call img_4%vis(geomorsphr=.false.)
        !     call img_3%bwd_ft()
        ! !    write(*,'(a)') '**info(simple_image PGI test): testing subtr '
        !   !  read *
        !     print *," L1 norm sum Orig v CUFFT", img.lonesum.img_2
        !     print *," L1 norm sum Orig v CUFFT", img.lonesum.img_2
        !     print *," L1 norm sum Orig v CUFFT", img.lonesum.img_2

        !  call img%subtr(img_2)
        !    call img%ft2img('amp',img_3)
        !    if(doplot)call img_3%vis(geomorsphr=.false.)
        !    write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '

        !   call img%add(img_2)
        ! read(*,*)
        !  write(*,'(a)') '**info(simple_image PGI test): comparing cufft and fftw3 '



        ! do i=1,ldim(1)
        !     do j=1,ldim(2)
        !         do k=1,ldim(3)
        !             if (abs(cabs(img%get_cmat_at_ind(i,j,k)) - cabs(img_2%get_cmat_at_ind(i,j,k))) >= TINY )then
        !                 print *,'wrong FFT i/j/k fftw cufft', i,j,k,real(img%get_cmat_at_ind(i,j,k)) ,&
        !                     &aimag(img%get_cmat_at_ind(i,j,k)) , real( img_2%get_cmat_at_ind(i,j,k)), &
        !                     & aimag( img_2%%get_cmat_at_ind(i,j,k))
        !                 call simple_stop("image FFT inverse PGI test  stopped")
        !             end if
        !         end do
        !     end do
        ! end do
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftw3 inverse fft'

        !   call img%ifft
        !  write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '
        ! call ifft_pgi_cuda_test(img_2)
        !  print *," L1 norm sum", img.lonesum.img_2

        ! do i=1,ldim(1)
        !     VerbosePrint "In simple_cufft::ifft_pgi testing ",i
        !     do j=1,ldim(2),2
        !         do k=1,ldim(3)
        !             if (abs( img%%get_rmat_at_ind(i,j,k) - img_2%get_rmat_at_ind(i,j,k)) >= 1.e-06 )then
        !                 print *,'wrong FFT i/j/k fftw /cufft', i,j,k,img%%get_rmat_at_ind(i,j,k) , img_2%%get_rmat_at_ind(i,j,k)
        !                 call simple_stop("image FFT PGI test stopped")
        !             end if
        !         end do

        !     end do
        ! end do
        VerbosePrint " test_pgi done"
        call img_2%kill
        call img%kill


    end subroutine test_pgi


    subroutine test_pgi1( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2
        integer              :: i, j, k, cnt, lfny, ldim(3)
        real                 :: input, msk, ave, sdev, var, med, xyz(3), pow
        real                 :: imcorr, recorr, corr, corr_lp, maxv, minv
        real                 :: smpd=2.
        logical              :: passed, test(6)


        verbose=.false.

        write(*,'(a)') '**info(unit_test, PGI CUFFT: test_pgi  '
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)
        !  call img_3%new(ldim, smpd)
        ! call img_4%new(ldim, smpd)

        call img%gauimg(10)
        call img%add_gauran(1.)
        if(doplot)call img%vis(geomorsphr=.false.)
        img_2 = img

        if(doplot)write(*,'(a)') '**info(unit_test PGI CUFFT: images created '
        call img_2%fft()
        if(doplot)call img_2%vis(geomorsphr=.false.)
        if(doplot)write(*,'(a)') '**info(simple_image PGI test): forward FFT completed '
        call img_2%ifft()
        if(doplot)write(*,'(a)') '**info(simple_image PGI test): inverse FFT completed '
        write(*,'(a,ES20.10)') "**info(simple_image PGI test): L1 norm sum ", img.lonesum.img_2
        if(doplot)call img_2%vis()


        !  write(*,'(a)') '**info(simple_image PGI test): comparing cufft and fftw3 '



        ! do i=1,ldim(1)
        !     do j=1,ldim(2)
        !         do k=1,ldim(3)
        !             if (abs(cabs(img%get_cmat_at_ind(i,j,k)) - cabs(img_2%get_cmat_at_ind(i,j,k))) >= TINY )then
        !                 print *,'wrong FFT i/j/k fftw cufft', i,j,k,real(img%get_cmat_at_ind(i,j,k)) ,&
        !                     &aimag(img%get_cmat_at_ind(i,j,k)) , real( img_2%get_cmat_at_ind(i,j,k)), &
        !                     & aimag( img_2%%get_cmat_at_ind(i,j,k))
        !                 call simple_stop("image FFT inverse PGI test  stopped")
        !             end if
        !         end do
        !     end do
        ! end do
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftw3 inverse fft'

        !   call img%ifft
        !  write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '
        ! call ifft_pgi_cuda_test(img_2)
        !  print *," L1 norm sum", img.lonesum.img_2

        ! do i=1,ldim(1)
        !     VerbosePrint "In simple_cufft::ifft_pgi testing ",i
        !     do j=1,ldim(2),2
        !         do k=1,ldim(3)
        !             if (abs( img%%get_rmat_at_ind(i,j,k) - img_2%get_rmat_at_ind(i,j,k)) >= 1.e-06 )then
        !                 print *,'wrong FFT i/j/k fftw /cufft', i,j,k,img%%get_rmat_at_ind(i,j,k) , img_2%%get_rmat_at_ind(i,j,k)
        !                 call simple_stop("image FFT PGI test stopped")
        !             end if
        !         end do

        !     end do
        ! end do
        VerbosePrint " test_pgi done"
        call img_2%kill
        call img%kill

    end subroutine test_pgi1


    subroutine test_pgi2( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2, img_3, img_4!, img_5, img3d
        !    type(image)          :: imgs(20)
        complex, allocatable :: ctmp(:,:,:), ctmp2(:,:,:)
        integer              :: i, j, k,h,l, cnt, lfny, ldim(3), lims(3,2), phys(3)
        real                 :: input, msk, ave, sdev, var, med, xyz(3), pow
        real                 :: imcorr, recorr, corr, corr_lp, maxv, minv
        !   real, allocatable    :: pcavec1(:), pcavec2(:), spec(:), res(:),rsh(:,:,:)
        real                 :: smpd=2.
        logical              :: passed, test(6)


        verbose=.false.
        !debug=.true.
        !write(*,'(a)') '**info(simple_image PGI unit_test, part 1): testing basal constructors'

        ! call img%new([ld1,ld2,1], 1.)

        ! call img_3%new([ld1,ld2,1], 1.)

        ! call img3d%new([ld1,ld2,1], 1.)
        ! call img_4%new([ld1,ld2,1], 1.)

        ! call img_4%set_rmat_at_ind(5,5,1, 1.)
        ! call img_4%set_rmat_at_ind(10,3,1, 1.)
        ! call img_4%set_rmat_at_ind(27,26,1, 1.)
        ! call img_4%set_rmat_at_ind(29,25,1, 1.)


        ! if( .not. img%exists() ) call simple_stop('ERROR, in constructor or in exists function, 1')
        ! if( .not. img3d%exists() ) call simple_stop('ERROR, in constructor or in exists function, 2')


        ! write(*,'(a)') '**info(simple_image_unit_test, PGI CUFFT 1): testing '
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)
        call img_3%new(ldim, smpd)
        !  call img_4%new(ldim, smpd)
        !     call img_5%new(ldim, smpd)
        call img%gauimg(10)
        call img%gauimg2(10, -floor(ld1/4.),-floor(ld2/3.))
        call img%gauimg2(10, floor(ld1/5.),-floor(ld2/8.))
        call img%gauimg2(10, floor(ld1/3.),floor(ld2/6.))
        call img%add_gauran(1.)
        !   call img%set_rmat_at_ind(floor(ld1/4.),floor(ld2/5.),1, 1.)
        !   call img%set_rmat_at_ind(floor(ld1/6.),floor(ld2/3.),1, 1.)
        if(doplot)call img%vis()
        img_2 = img
        img_3 = img
        !     if(doplot)call img_2%vis
        !  ave = (img.lonesum.img_2)
        !  print *," L1 norm sum", ave
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftshifter bwd'
        !     write(*,'(a)') 'Press enter to test ifftshift '
        !     read(*,*)
        !  call img_2%bwd_fftshift()
        !  img_4= img_2
        !  if(doplot)call img_4%vis(geomorsphr=.false.)


        !write(*,'(a)') 'Press enter to test cufft '
        !read(*,*)
        write(*,'(a)') '**info(simple_image PGI test): testing cufft '

        call fft_pgi_cuda_test(img_2)
        call simple_cuda_stop("In PGI TEST2 ",__FILENAME__,__LINE__)
        write(*,'(a)') '**info(simple_image PGI test): testing fftw '
        call img_3%fwd_ft()
        print *,"FFTW"
        call img_3%print_cmat()
        !call img_2%ft2img('amp',img_4)
        print *,"CUFFT"
        call img_2%print_cmat()


        !     print *,"Diff"
        !     call img_4%print_cmat()

        if(doplot)call img_3%vis(geomorsphr=.false.)
        if(doplot)call img_2%vis(geomorsphr=.false.)

        !   img_4 = img_2
        !   img_5 = img_3
        print *," L1 norm sum FFTW v CUFFT", img_3.lonesum.img_2
        !  lims=img_3%loop_lims(2)
        !   write(*,'(a)') '**info(simple_image PGI test): cufft returned '
        !   write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '
        !     call ifft_pgi_cuda_test(img_2)
        read(*,*)
        ! if(doplot)call img_2%vis
        ! print *," L1 norm sum", img.lonesum.img_2

        ! call img_3%kill
        ! call img_4%kill
        ! call img3d%kill

        !img_4 = img_2
        !    call img_2%ft2img('amp',img_4)
        !    if(doplot)call img_4%vis(geomorsphr=.false.)
        !  write(*,'(a)') '**info(simple_image PGI test): testing fftshifter fwd'
        !  write(*,'(a)') 'Press enter to test fftshift '
        !    read(*,*)
        !     call img_3%fwd_fftshift()
        !    img_4 = img_2
        !    if(doplot)call img_2%vis(geomorsphr=.false.)


        !     write(*,'(a)') 'Press enter to test fftw3 '
        !    read(*,*)


        !        call img_3%ft2img('amp',img_4)
        ! ! if(doplot)call img_3%vis(geomorsphr=.false.)
        !   call img_3%bwd_ft()
        !   if(doplot)call img_3%vis()
        !   if(doplot)call img_2%vis()
        ! ! !    write(*,'(a)') '**info(simple_image PGI test): testing subtr '
        ! ! !  read *
        !   print *," L1 norm sum Orig v CUFFT", img.lonesum.img_2
        !  print *," L1 norm sum Orig v FFTW", img.lonesum.img_3
        !   print *," L1 norm sum FFTW v CUFFT", img_3.lonesum.img_2

        !  call img%subtr(img_2)
        !    call img%ft2img('amp',img_3)
        !    if(doplot)call img_3%vis(geomorsphr=.false.)
        !    write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '

        !   call img%add(img_2)
        !read(*,*)
        write(*,'(a)') '**info(simple_image PGI test): comparing FWD cufft and fftw3 '
        passed=.true.
        if(allocated(ctmp)) deallocate(ctmp)
        if( allocated(ctmp2) ) deallocate(ctmp2)
        !$ allocate(ctmp(ldim(1)/2 +1,ldim(2),ldim(3)),ctmp2(ldim(1)/2 +1,ldim(2),ldim(3)) )
        ctmp2= img_2%get_cmatfull(); ctmp= img_3%get_cmatfull()
        VerbosePrint " Ldim/shapes ", ldim, shape(ctmp), shape(ctmp2)

        do i=1,ldim(1)/2 +1 !lims(1,1),lims(1,2)
            do j=1,ldim(2)!lims(2,1),lims(2,2)
                do k=1,ldim(3)!lims(3,1),lims(3,2)
                    ! phys = img_3%comp_addr_phys([h,k,l])
                    !if (abs(cabs(img_2%get_cmat_at(phys)) - cabs(img_3%get_cmat_at(phys))) >= 1.e-06 )then
                    !   print *,'wrong FFT i/j/k fftw cufft', h,k,l,img_3%get_cmat_at(phys) ,&
                    !      & img_2%get_cmat_at(phys)
                    ! passed = .false.
                    !                 call simple_stop("image FFT inverse PGI test  stopped")
                    if ((abs(real(ctmp2(i,j,k)) - real(ctmp(i,j,k))) >= 1.e-06 ) .or. &
                        &(abs(aimag(ctmp2(i,j,k)) - aimag(ctmp(i,j,k))) >= 1.e-06)) then
                        print *,'wrong IFFT i/j/k fftw cufft', i,j,k,ctmp2(i,j,k),ctmp(i,j,k)
                    end if
                end do
            end do
        end do
        write(*,'(a)') '**info(simple_image PGI test): compare fftw and cufft (images) '
        ctmp = ctmp - ctmp2
        call gnufor_image(real(ctmp(:ldim(1),:ldim(2),1)), palette='gray')
        call gnufor_image(aimag(ctmp(:ldim(1),:ldim(2),1)), palette='gray')
        call gnufor_image(cabs(ctmp(:ldim(1),:ldim(2),1)), palette='gray')
        call gnufor_image(atan2(real(ctmp(:ldim(1),:ldim(2),1)),aimag(ctmp(:ldim(1),:ldim(2),1))), palette='gray')
        read(*,*)
        img_4 = img_2 - img_3
        call img_4%set_ft(.true.)
        if(doplot)call img_4%vis(geomorsphr=.false.)
        if(doplot)call img_4%vis(geomorsphr=.true.)
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftw3 inverse fft'

        !   call img%ifft
        !  write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '
        ! call ifft_pgi_cuda_test(img_2)
        !  print *," L1 norm sum", img.lonesum.img_2

        !  do i=1,ldim(1)
        ! !     VerbosePrint "In simple_cufft::ifft_pgi testing ",i
        !      do j=1,ldim(2)
        !          do k=1,ldim(3)
        !              if (abs( img_2%get([i,j,k]) - img_3%get([i,j,k])) >= 1.e-06 )then
        !                 print *,'wrong IFFT i/j/k fftw /cufft', i,j,img_3%get([i,j,k]) , img_2%get([i,j,k])
        ! !                 call simple_stop("image FFT PGI test stopped")
        !              end if
        !          end do
        !      end do
        !  end do
        !     VerbosePrint "In simple_cufft::ifft_pgi testing done"
        read(*,*)
        call img_4%kill
        !      call img_5%kill
        call img_2%kill
        call img_3%kill
        call img%kill        ! write(*,'(a)') '**info(simple_image_unit_test, part 7): testing fftshift'

        if( allocated(ctmp) ) deallocate(ctmp)
        if( allocated(ctmp2) ) deallocate(ctmp2)
        ! passed = .false.
        ! msk=50
        ! call img%gauimg(10)
        ! write(*,'(a)') '**info(simple_image PGI unit_test, part 2):'
        ! if( doplot ) call img%vis
        ! write(*,'(a)') '**info(simple_image PGI unit_test, part 2):'
        ! call img%serialize(pcavec1, msk)
        ! call img%shift([-9.345,-5.786,0.])
        ! if( doplot ) call img%vis
        ! call img%shift([9.345,5.786,0.])
        ! call img%serialize(pcavec2, msk)
        ! if( doplot ) call img%vis
        ! if( pearsn(pcavec1, pcavec2) > 0.99 ) passed = .true.
        if( .not. passed )  call simple_stop('pgi fft test failed')

    end subroutine test_pgi2


    subroutine test_pgi3( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        use gnufor2
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2, img_3, img_4!, img_5
        complex, allocatable ::  ctmp(:,:,:)
        integer              :: i, j, k,h,l, cnt, lfny, ldim(3), lims(3,2), phys(3)
        real, allocatable    :: rtmp(:,:,:),rtmp2(:,:,:)
        real                 :: smpd=2.
        logical              :: passed, test(6)

        verbose=doplot

        VerbosePrint "In simple_cufft::ifft_pgi testing start"
        ! write(*,'(a)') '**info(simple_image_unit_test, PGI CUFFT 1): testing '
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)
        call img_3%new(ldim, smpd)
        call img_4%new(ldim, smpd)
        ! call img_5%new(ldim, smpd)
        img=2.; img_2=0.; img_3=0.; img_4=0.
        call img%gauimg(10)
        call img%add_gauran(.5)
        ! call img%set_rmat_at_ind(8,5,1, 1.5)
        ! call img%set_rmat_at_ind(6,3,1, 1.)
        !$ allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        rtmp = img%get_rmat()
        call fftshift(rtmp)
        VerbosePrint "In simple_cufft::ifft_pgi testing fftshift"
        !   rtmp(8,5,1) = 1.5
        !  rtmp(6,3,1)= 1.
        !  rtmp(6:8,:,1)=0.
        call img_2%set_rmat(rtmp); call img_2%set_ft(.false.)
        call gnufor_image(rtmp(:,:,1),palette='gray')
        !call gnufor_image(rtmp(ldim(1):-1:1,ldim(2):-1:1,1),palette='gray')
        allocate(ctmp(ldim(1),ldim(2),ldim(3)))
        ctmp = cmplx(rtmp,abs(1.0-rtmp))

        call img_2%set_cmat(ctmp); call img_2%set_ft(.true.)
        call img_3%set_cmat(ctmp); call img_3%set_ft(.true.)
        VerbosePrint "In simple_cufft::ifft_pgi testing set_cmat"
        !     if(doplot)call img_2%vis
        !  ave = (img.lonesum.img_2)
        !  print *," L1 norm sum", ave
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftshifter bwd'
        !     write(*,'(a)') 'Press enter to test ifftshift '
        !     read(*,*)
        !  call img_2%bwd_fftshift()
        !  img_4= img_2
        !  if(doplot)call img_4%vis(geomorsphr=.false.)
        !
        !  write(*,'(a)') 'Press enter to test cufft '
        ! read(*,*)
        write(*,'(a)') '**info(simple_image PGI test): testing cufft '

        call ifft_pgi_cuda_test(img_2)
        call simple_cuda_stop("In PGI TEST2 ",__FILENAME__,__LINE__)
        call img_3%bwd_ft()
        if(doplot)call img_3%vis()
        print *,"FFTW"
        call img_3%print_rmat
        !call img_2%ft2img('amp',img_4)
        if(doplot)call img_2%vis()
        print *,"CUFFT"
        call img_2%print_rmat()
        print *," test_pgi3 L1 norm sum FFTW v CUFFT", img_3.lonesum.img_2



        !   write(*,'(a)') '**info(simple_image PGI test): cufft returned '
        !   write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '


        ! if(doplot)call img_2%vis
        ! print *," L1 norm sum", img.lonesum.img_2

        ! call img_3%kill
        ! call img_4%kill
        ! call img3d%kill

        !img_4 = img_2
        !    call img_2%ft2img('amp',img_4)
        !    if(doplot)call img_4%vis(geomorsphr=.false.)
        !  write(*,'(a)') '**info(simple_image PGI test): testing fftshifter fwd'
        !  write(*,'(a)') 'Press enter to test fftshift '
        !    read(*,*)
        !     call img_3%fwd_fftshift()
        !    img_4 = img_2
        !    if(doplot)call img_2%vis(geomorsphr=.false.)


        !     write(*,'(a)') 'Press enter to test fftw3 '
        !    read(*,*)


        !        call img_3%ft2img('amp',img_4)
        ! if(doplot)call img_3%vis(geomorsphr=.false.)

        lims=img_3%loop_lims(2)
        VerbosePrint "In simple_cufft::ifft_pgi testing set_cmat"
        write(*,'(a)') 'Press enter to test fftw '
        VerbosePrint "In simple_cufft::ifft_pgi testing set_cmat"
        read(*,*)


        !  call img%subtr(img_2)
        !    call img%ft2img('amp',img_3)
        !    if(doplot)call img_3%vis(geomorsphr=.false.)
        !    write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '
        ldim = img_2%get_ldim()
        if(allocated(rtmp)) deallocate(rtmp)
        if( allocated(rtmp2) ) deallocate(rtmp2)
        !$ allocate(rtmp(ldim(1),ldim(2),ldim(3)),rtmp2(ldim(1),ldim(2),ldim(3)) )
        rtmp2= img_2%get_rmat(); rtmp= img_3%get_rmat()
        print *," test_pgi3 Ldim/shapes ", ldim, shape(rtmp), shape(rtmp2)
        do i=1,ldim(1)
            !     VerbosePrint "In simple_cufft::ifft_pgi testing ",i
            do j=1,ldim(2)
                do k=1,ldim(3)
                    if (abs( rtmp(i,j,k) - rtmp2(i,j,k)) >= 1.e-08 )then
                        write(*,'(A,3I3,1x,2F20.10)'), ' fftw /cufft', i,j,k,rtmp(i,j,k) , rtmp2(i,j,k)
                        !                 call simple_stop("image FFT PGI test stopped")
                    end if
                end do
            end do
        end do

        write(*,'(a)') 'Press enter to test difference '
        read(*,*)
        img_4= img_2 - img_3
        call img_4%print_rmat()
        if(doplot) call img_4%vis()




        if( allocated(rtmp) ) deallocate(rtmp)
        if( allocated(rtmp2) ) deallocate(rtmp2)
        if( allocated(ctmp) ) deallocate(ctmp)

        ! call img_5%kill
        call img_4%kill
        call img_3%kill
        call img_2%kill
        call img%kill
        VerbosePrint "In simple_cufft::ifft_pgi testing done"
    end subroutine test_pgi3


    subroutine test_pgi4( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        use gnufor2
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2, img_3
        complex, allocatable :: ctmp(:,:,:)
        integer              :: i, j, k,h,l, cnt, lfny, ldim(3), phys(3), cdim(3)
        !   complex(kind=c_float_complex), pointer :: cmat(:,:,:)=>null()
        real, allocatable    :: rtmp(:,:,:),rtmp2(:,:,:)
        real                 :: smpd=2.
        logical              :: passed, test(6)
        integer(timer_int_kind):: t4
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)

        verbose=doplot
        if(verbose)t4=tic()
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        VerbosePrint "In test_pgi4 start"
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        ! write(*,'(a)') '**info(simple_image_unit_test, PGI CUFFT 1): testing '
        passed = .false.
        ldim = [ld1,ld2,1]
        !cdim = [9,16,1]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)
        call img_3%new(ldim, smpd)
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        img=2.; img_2=0.; img_3=0.
        call img%gauimg(10)
        call img%gauimg2(5, floor(real(ld1)/5.),floor(real(ld1)/3.))
        call img%gauimg2(8, -floor(real(ld1)/5.),floor(real(ld1)/6.))

        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        rtmp = img%get_rmat()
        call ifftshift(rtmp)
        call img%set_rmat(rtmp)
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        if(verbose)print *, "testing fftshift",toc()
        !call fftshift(rtmp)
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        VerbosePrint "In test_pgi4 testing fftshift"
        !call img%set_rmat_at_ind(17,16,1,10.)
        !rtmp(6,3,1) = 1.
        if(doplot)call img%vis()

        call img%set_ft(.true.)
        ! call img%fwd_ft()
        img_2 = img
        img_3 = img
        cdim = img%get_array_shape()
        !call img_2%set_rmat(rtmp); call img_2%set_ft(.false.)
        !call gnufor_image(rtmp(:,:,1),palette='gray')
        ! if(verbose)read(*,*)
        !call gnufor_image(rtmp(ldim(1):-1:1,ldim(2):-1:1,1),palette='gray')
        !call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)

        if(doplot)call img%vis(geomorsphr=.true.)
        !  if(doplot)call img_2%vis(geomorsphr=.true.)
        !  if(doplot)call img_3%vis(geomorsphr=.true.)
        print *,' pgi4  ctmp array shape ', cdim
        allocate(ctmp(cdim(1),cdim(2),cdim(3)))
        ctmp = img%get_cmatfull()
        !call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        !        call gnufor_image(real(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
        !         call gnufor_image(aimag(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
        !         call ifftshift(ctmp)
        !    call gnufor_image(real(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
        !    call gnufor_image(aimag(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
        ! read(*,*)
        VerbosePrint "In test_pgi4 testing cufft ifft_pgi_cuda"
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        ! if(verbose)read(*,*)
        if(verbose)t4=tic()
        call ifft_pgi_cuda_test3(img_2)
        if(img_2%is_ft())print *, " cufft ifft test failed",toc(t4)
        if(doplot)call img_2%vis()
        if(verbose)print *, "test completed backward cufft ifft",toc(t4)
        read(*,*)
        if(doplot)call img_2%vis()
        if(verbose)print *,"CUFFT "
        !if(verbose)call img_2%print_rmat()
        !if(verbose)read(*,*)
        if(verbose)t4=tic()
        call img_3%bwd_ft()
        if(verbose)print *, "testing backward cufft ifft",toc(t4)
        if(doplot)call img_3%vis()
        if(verbose)print *,"FFTW"
        !f(verbose)call img_3%print_rmat

        print *," L1 norm sum FFTW v CUFFT", img_3.lonesum.img_2
        call img%set_ft(.false.)
        img = img_3-img_2
        if(doplot)call img%vis()
        if(verbose)read(*,*)

        if(allocated(rtmp)) deallocate(rtmp)
        if( allocated(rtmp2) ) deallocate(rtmp2)
        allocate(rtmp(ldim(1),ldim(2),ldim(3)),rtmp2(ldim(1),ldim(2),ldim(3)) )
        rtmp2= img_2%get_rmat(); rtmp= img_3%get_rmat()
        ! print *," Ldim/shapes ", ldim, shape(rtmp), shape(rtmp2)
        ! cnt=0
        ! do i=1,ldim(1)
        !     !     VerbosePrint "In test_pgi4 testing diff",i
        !     do j=1,ldim(2)
        !         do k=1,ldim(3)
        !             if (abs( rtmp(i,j,k) - rtmp2(i,j,k)) >= 1.e-06 )then
        !                 write(*,'(A,3I0,1x,2ES20.10)'), ' fftw /cufft', i,j,k,rtmp(i,j,k) , rtmp2(i,j,k)
        !                 !                 call simple_stop("image FFT PGI test stopped")
        !                 cnt=cnt+1
        !             end if
        !         end do
        !     end do
        ! end do
        ! write(*,'(A,3I0,1x,2ES20.10)'), ' fftw /cufft', i,j,k,rtmp(i,j,k) , rtmp2(i,j,k)
        ! VerbosePrint "In test_pgi4 testing done"

        if( allocated(rtmp) ) deallocate(rtmp)
        if( allocated(rtmp2) ) deallocate(rtmp2)
        if( allocated(ctmp) ) deallocate(ctmp)

        call img_3%kill
        call img_2%kill
        call img%kill
        return
    end subroutine test_pgi4


    !     subroutine test_multigpu_C2C_3Dfft
    !         use cufft
    !         ! 3D Complex-to-Complex Transforms using Two GPUs
    !         ! In this example a three-dimensional complex-to-complex transform is applied to the input data using two GPUs.

    !         ! Demonstrate how to use CUFFT to perform 3-d FFTs using 2 GPUs
    !         complex, allocatable :: host_data_input(:,:,:), host_data_output(:,:,:)
    !         ! Create an empty plan
    !         type(cufftPlan) plan
    !         integer ::  res,nGPUs, whichGPUs(2),worksize(2), nx , ny , nz, size_of_data

    !         res = cufftCreate(&plan_input);
    !         if (res /= CUFFT_SUCCESS) { printf ("*Create failed\n"); return; }

    !     ! cufftXtSetGPUs() - Define which GPUs to use
    !     nGPUs = 2
    !     whichGPUs(0) = 0; whichGPUs(1) = 1;
    !     res = cufftXtSetGPUs (plan_input, nGPUs, whichGPUs)
    !     if (res /= CUFFT_SUCCESS)then
    !         print *,"*XtSetGPUs failed"
    !         return
    !     end if

    !     ! Initialize FFT input data


    !      nx = 64, ny = 128, nz = 32;
    !     size_of_data = 8 * nx * ny * nz;
    !     allocate(host_data_input(nx,ny,nz),host_data_output(nx,ny,nz)

    !     initialize_3d_data (nx, ny, nz, host_data_input, host_data_output);

    ! ! cufftMakePlan3d() - Create the plan
    !     res = cufftMakePlan3d (plan_input, nz, ny, nx, CUFFT_C2C, worksize);
    !     if (res /= CUFFT_SUCCESS)then
    !         print *,"*MakePlan* failed"
    !         return
    !         end if

    ! ! cufftXtMalloc() - Malloc data on multiple GPUs
    !     cudaLibXtDesc *device_data_input;
    !     result = cufftXtMalloc (plan_input, &device_data_input,
    !         CUFFT_XT_FORMAT_INPLACE);
    !     if (result != CUFFT_SUCCESS) { printf ("*XtMalloc failed\n"); return; }
    ! //
    ! // cufftXtMemcpy() - Copy data from host to multiple GPUs
    !     result = cufftXtMemcpy (plan_input, device_data_input,
    !         host_data_input, CUFFT_COPY_HOST_TO_DEVICE);
    !     if (result != CUFFT_SUCCESS) { printf ("*XtMemcpy failed\n"); return; }
    ! //
    ! // cufftXtExecDescriptorC2C() - Execute FFT on multiple GPUs
    !     result = cufftXtExecDescriptorC2C (plan_input, device_data_input,
    !         device_data_input, CUFFT_FORWARD);
    !     if (result != CUFFT_SUCCESS) { printf ("*XtExec* failed\n"); return; }
    ! //
    ! // cufftXtMemcpy() - Copy data from multiple GPUs to host
    !     result = cufftXtMemcpy (plan_input, host_data_output,
    !         device_data_input, CUFFT_COPY_DEVICE_TO_HOST);
    !     if (result != CUFFT_SUCCESS) { printf ("*XtMemcpy failed\n"); return; }
    ! //
    ! // Print output and check results
    !     int output_return = output_3d_results (nx, ny, nz,
    !         host_data_input, host_data_output);
    !     if (output_return != 0) { return; }
    ! //
    ! // cufftXtFree() - Free GPU memory
    !     result = cufftXtFree(device_data_input);
    !     if (result != CUFFT_SUCCESS) { printf ("*XtFree failed\n"); return; }
    ! //
    ! // cufftDestroy() - Destroy FFT plan
    !     result = cufftDestroy(plan_input);
    !     if (result != CUFFT_SUCCESS) { printf ("*Destroy failed: code\n"); return; }
    !     free(host_data_input); free(host_data_output);


    ! Read more at: http://docs.nvidia.com/cuda/cufft/index.html#ixzz4xED7Iv1J
    ! Follow us: @GPUComputing on Twitter | NVIDIA on Facebook

    ! end



#endif



    subroutine test_fftw( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2
        integer              :: ldim(3)
        real                 :: smpd=2.
        logical              :: passed

        verbose=.false.
        !debug=.true.
        write(*,*) '**info(unit_test, FFTW  test_fftw'
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)

        call img%gauimg(10)
        call img%add_gauran(1.)

        img_2=img

        call img_2%fwd_ft()
        call img_2%bwd_ft()
        write(*,'(a,1ES20.10)') '**info(unit_test, FFTW   L2 norm sum', img.lonesum.img_2
        call img_2%kill
        call img%kill

    end subroutine test_fftw

end module simple_cufft_test
