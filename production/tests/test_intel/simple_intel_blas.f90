module simple_intel_blas
    include 'simple_lib.f08'
    use simple_image,       only: image, test_image
    use gnufor2
    implicit none


#ifdef INTEL
    include 'mkl.fi'
    public :: test_intel_blas

    public ::  test_mkl_fftw
    public ::  basic_sp_real_dft_2d
    public ::  config_thread_limit
#endif
    public :: test_omp
    private

    ! module global constants
    integer, parameter :: SQRAD=60, NTST=50, NNOISY=20
    real,    parameter :: SMPD=1.77, TRS=10., HP=100.0, LP=8., SNR=0.2

    real(kind=sp), parameter :: tolerance = 1.e-05_sp
    ! module global variables

    ! type(image)              :: img, img_shifted
    ! type(image), allocatable :: noisy_imgs(:)
    ! integer                  :: x, y


#include "simple_local_flags.inc"
contains
    ! Step3  Jacobi relaxation
    ! NVIDIA blog example https://github.com/parallel-forall/code-samples
    ! ACC removed for comparison
    subroutine test_omp(errflag)

        logical, intent(inout):: errflag
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

        do while ( error .gt. tol .and. iter .lt. iter_max )
            error=0.0_fp_kind
            !$omp parallel do collapse(2) default(shared) private(i,j) reduction( max:error )
            do j=1,m-2
                do i=1,n-2
                    Anew(i,j) = 0.25_fp_kind * ( A(i+1,j  ) + A(i-1,j  ) + &
                        A(i  ,j-1) + A(i  ,j+1) )
                    error = max( error, abs(Anew(i,j)-A(i,j)) )
                end do
            end do
            !$omp end parallel do

            if(mod(iter,100).eq.0 ) write(*,'(i5,f10.6)'), iter, error
            iter = iter +1

            !$omp parallel do  collapse(2) private(i,j) default(shared)
            do j=1,m-2
                do i=1,n-2
                    A(i,j) = Anew(i,j)
                end do
            end do
            !$omp end parallel do
        end do

        call cpu_time(stop_time)
        write(*,'(a,f10.3,a)')  ' completed in ', stop_time-start_time, ' seconds'

        deallocate (A,Anew,y0)

    end subroutine test_omp

    !>  https://software.intel.com/en-us/mkl-developer-reference-fortran-blas-code-examples
#ifdef INTEL
    subroutine test_intel_blas(flag)
        use simple_defs
        implicit none
        logical, intent(inout):: flag

        write(*,'(a)') '**info(simple_intel unit tests, part 1): testing blas sdot'
        call dot_main(flag)
        write(*,'(a)') '**info(simple_intel unit tests, part 2): testing blas scopy'
        call copy_main(flag)
        write(*,'(a)') '**info(simple_intel unit tests, part 3): testing blas sger'
        call ger_main(flag)
        write(*,'(a)') '**info(simple_intel unit tests, part 4): testing blas ssymm'
        call symm_main(flag)
        write(*,'(a)') '**info(simple_intel unit tests, test_intel_blas completed'
    contains


        subroutine dot_main(errflag)
            logical, intent(inout):: errflag
            real x(10), y(10), sdot, res
            integer n, incx, incy, i
            !external sdot
            n = 5
            incx = 2
            incy = 1
            do i = 1, 10
                x(i) = 2.0e0
                y(i) = 1.0e0
            end do
            write(*,'(a)') '**info(simple_intel unit tests, part 1): sdot'
            res = sdot (n, x, incx, y, incy)
            print*, 'SDOT = ', res
            if (res /= 10.0) then
                print*,' SDOT did not return expected 10.0, result: ', res
                errflag=.true.
            end if

        end subroutine dot_main

        subroutine copy_main(errflag)
            logical, intent(inout):: errflag
            real x(10), y(10)
            integer n, incx, incy, i
            n = 3
            incx = 3
            incy = 1
            do i = 1, 10
                x(i) = i
            end do
            call scopy (n, x, incx, y, incy)
            print*, 'Y = ', (y(i), i = 1, n)
            if (y(1) /= 1.0 .or. y(4) /= 4.0 .or. y(7) /= 7.0) then
                print*, ' SCOPY did not return expected (1,...,10), result: ', y
                errflag=.true.
            end if
        end subroutine copy_main

        ! The following example illustrates a call to the BLAS Level 2 routine sger. This
        ! routine performs a matrix-vector operation a := alpha*x*y' + a.
        subroutine ger_main(errflag)
            logical, intent(inout):: errflag
            real a(5,3), x(10), y(10), alpha
            integer m, n, incx, incy, i, j, lda
            m = 2
            n = 3
            lda = 5
            incx = 2
            incy = 1
            alpha = 0.5
            do i = 1, 10
                x(i) = 1.0
                y(i) = 1.0
            end do
            do i = 1, m
                do j = 1, n
                    a(i,j) = j
                end do
            end do
            call sger (m, n, alpha, x, incx, y, incy, a, lda)
            if (a(1,1) /= 1.5 .or. y(4) /= 4.0 .or. y(7) /= 7.0) then
                errflag=.true.
                print*, 'Matrix A: '
                do i = 1, m
                    print*, (a(i,j), j = 1, n)
                end do
            end if
        end subroutine ger_main

        !        The following example illustrates a call to the BLAS Level 3 routine ssymm. This
        !        routine performs a matrix-matrix operation
        !c :=  alpha*a*b' + beta*c.

        subroutine symm_main(errflag)
            logical, intent(inout):: errflag
            real a(3,3), b(3,2), c(3,3), alpha, beta
            integer m, n, lda, ldb, ldc, i, j
            character uplo, side
            uplo = 'u'
            side = 'l'
            m = 3
            n = 2
            lda = 3
            ldb = 3
            ldc = 3
            alpha = 0.5
            beta = 2.0
            do i = 1, m
                do j = 1, m
                    a(i,j) = 1.0
                end do
            end do
            do i = 1, m
                do j = 1, n
                    c(i,j) = 1.0
                    b(i,j) = 2.0
                end do
            end do
            call ssymm (side, uplo, m, n, alpha,  a, lda, b, ldb, beta, c, ldc)
            print*, 'Matrix C: '
            do i = 1, m
                print*, (c(i,j), j = 1, n)
            end do
            if (c(1,1) /= 5.0 .or. c(3,3) /= 5.0) then
                print*, ' SSYMM did not return expected (5), result: ', c
                errflag=.true.
            end if
        end subroutine symm_main

    end subroutine test_intel_blas






    subroutine test_mkl_fftw (errflag)
        use MKL_DFTI
        logical, intent(inout) :: errflag

        call test_1d(errflag)
        call test_2d(errflag)
        call test_2d_real_inplace(errflag)
        call test_mkl_fft_openmp(errflag)
        call test_mkl_fft2d_array_descr_main(errflag)
        call test_mkl_fft2d_shared_descr_main(errflag)

    contains

        subroutine test_1d (errflag)

            !One-dimensional In-place FFT
            ! Fortran example.
            ! 1D complex to complex, and real to conjugate-even

            logical, intent(inout) :: errflag
            complex :: X(32)
            real :: Y(34)
            type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
            integer :: Status
            !...put input data into X(1),...,X(32); Y(1),...,Y(32)

            ! Perform a complex to complex transform
            status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_SINGLE,&
                DFTI_COMPLEX, 1, 32 )
            status = DftiCommitDescriptor( My_Desc1_Handle )
            status = DftiComputeForward( My_Desc1_Handle, X )
            status = DftiFreeDescriptor(My_Desc1_Handle)
            ! result is given by {X(1),X(2),...,X(32)}

            ! Perform a real to complex conjugate-even transform
            status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,&
                DFTI_REAL, 1, 32)
            status = DftiCommitDescriptor(My_Desc2_Handle)
            status = DftiComputeForward(My_Desc2_Handle, Y)
            status = DftiFreeDescriptor(My_Desc2_Handle)
            ! result is given in CCS format.
            if (status .ne. 0) then
                errflag=.true.
                if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                    print *, 'Error: ', DftiErrorMessage(status)
                endif
            endif

        end subroutine test_1d


        subroutine test_2d (errflag)
            ! Fortran example.
            ! 2D complex to complex, and real to conjugate-even

            logical, intent(inout) :: errflag
            complex ::  X_2D(32,100)
            real :: Y_2D(34, 102)
            complex ::  X(3200)
            real :: Y(3468)
            equivalence (X_2D, X)
            equivalence (Y_2D, Y)
            type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
            integer :: status, L(2)
            !...put input data into X_2D(j,k), Y_2D(j,k), 1<=j=32,1<=k<=100
            !...set L(1) = 32, L(2) = 100
            !...the transform is a 32-by-100

            ! Perform a complex to complex transform
            status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_SINGLE,&
                DFTI_COMPLEX, 2, L)
            status = DftiCommitDescriptor( My_Desc1_Handle)
            status = DftiComputeForward( My_Desc1_Handle, X)
            status = DftiFreeDescriptor(My_Desc1_Handle)
            ! result is given by X_2D(j,k), 1<=j<=32, 1<=k<=100

            ! Perform a real to complex conjugate-even transform
            status = DftiCreateDescriptor( My_Desc2_Handle, DFTI_SINGLE,&
                DFTI_REAL, 2, L)
            status = DftiCommitDescriptor( My_Desc2_Handle)
            status = DftiComputeForward( My_Desc2_Handle, Y)
            status = DftiFreeDescriptor(My_Desc2_Handle)
            ! result is given by the complex value z(j,k) 1<=j<=32; 1<=k<=100
            ! and is stored in CCS format
            if (status .ne. 0) then
                errflag=.true.
                if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                    print *, 'Error: ', DftiErrorMessage(status)
                end if
            end if

        end subroutine test_2d

        ! The following code illustrates real multi-dimensional transforms with CCE format storage
        ! of conjugate-even complex matrix. Example “Two-Dimensional REAL In-place FFT (Fortran
        ! Interface)” is two-dimensional in-place transform and Example “Two-Dimensional REAL
        ! Out-of-place FFT (Fortran Interface)” is two-dimensional out-of-place transform in Fortran
        ! interface. Note that the data and result parameters in computation functions are all
        ! declared as assumed-size rank-1 array DIMENSION(0:*). Therefore two-dimensional array must
        ! be transformed to one-dimensional array by EQUIVALENCE statement or other facilities of
        ! Fortran.

        !Two-Dimensional REAL In-place FFT
        ! Fortran example.
        ! 2D and real to conjugate-even

        subroutine test_2d_real_inplace (errflag)


            logical, intent(inout) :: errflag
            real :: X_2D(34,100) ! 34  = (32/2 + 1)*2
            real :: X(3400)
            equivalence (X_2D, X)
            type(DFTI_DESCRIPTOR), POINTER :: My_Desc_Handle
            integer :: status, L(2)
            integer :: strides_in(3)
            integer :: strides_out(3)
            ! ...put input data into X_2D(j,k), 1<=j=32,1<=k<=100
            L(1) = 32; L(2) = 100
            strides_in(1) = 0; strides_in(2) = 1; strides_in(3) = 34
            strides_out(1) = 0; strides_out(2) = 1; strides_out(3) = 17
            ! ...the transform is a 32-by-100
            ! Perform a real to complex conjugate-even transform
            status = DftiCreateDescriptor( My_Desc_Handle, DFTI_SINGLE,&
                DFTI_REAL, 2, L )
            status = DftiSetValue(My_Desc_Handle, DFTI_CONJUGATE_EVEN_STORAGE,&
                DFTI_COMPLEX_COMPLEX)
            status = DftiSetValue(My_Desc_Handle, DFTI_INPUT_STRIDES, strides_in)
            status = DftiSetValue(My_Desc_Handle, DFTI_OUTPUT_STRIDES, strides_out)
            status = DftiCommitDescriptor( My_Desc_Handle)
            status = DftiComputeForward( My_Desc_Handle, X )
            status = DftiFreeDescriptor(My_Desc_Handle)
            ! result is given by the complex value z(j,k) 1<=j<=17; 1<=k<=100 and
            ! is stored in real matrix X_2D in CCE format.
            if (status .ne. 0) then
                errflag=.true.
                if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                    print *, 'Error: ', DftiErrorMessage(status)
                endif
            endif

        end subroutine test_2d_real_inplace

        ! Using Parallel Mode with Multiple Descriptors Initialized in a Parallel Region Note that
        ! in this example, the program can be transformed to become single-threaded at the customer
        ! level but using parallel mode within Intel MKL. To achieve this, you need to set the
        ! parameter DFTI_NUMBER_OF_TRANSFORMS = 4 and to set the corresponding parameter
        ! DFTI_INPUT_DISTANCE = 5000.

        ! Specify the number of OpenMP threads for Example “Using Parallel Mode with Multiple
        ! Descriptors Initialized in One Thread” like this:

        ! set MKL_NUM_THREADS = 1 for Intel MKL to work in the single-threaded mode (obligatory);

        ! set OMP_NUM_THREADS = 4 for the customer program to work in the multi-threaded mode.


        subroutine test_mkl_fft_openmp (errflag)

            logical, intent(inout) :: errflag
            integer nth, len(2)
            ! 4 OMP threads, each does 2D FFT 50x100 points
            parameter (nth = 4, len = (/50, 100/))
            complex x(len(2)*len(1), nth)

            type(dfti_descriptor), pointer :: myFFT
            integer th, status, mystatus

            ! assume x is initialized and do 2D FFTs
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(myFFT, mystatus) reduction(+:status)
            do th = 1, nth
                mystatus = mystatus + DftiCreateDescriptor (myFFT, DFTI_SINGLE, DFTI_COMPLEX, 2, len)
                mystatus = mystatus + DftiCommitDescriptor (myFFT)
                mystatus = mystatus + DftiComputeForward (myFFT, x(:, th))
                mystatus = mystatus + DftiFreeDescriptor (myFFT)
                status = status + mystatus
            end do
            !$OMP END PARALLEL DO

            if (status .ne. 0) then
                errflag=.true.
                if (.not. DftiErrorClass(status,DFTI_NO_ERROR))then
                    print *, 'Error: ', DftiErrorMessage(status)
                endif
            endif

        end subroutine test_mkl_fft_openmp

        subroutine test_mkl_fft2d_array_descr_main (errflag)
            logical, intent(inout) :: errflag


            integer nth, len(2)
            ! 4 OMP threads, each does 2D FFT 50x100 points
            parameter (nth = 4, len = (/50, 100/))
            complex x(len(2)*len(1), nth)

            type thread_data
                type(dfti_descriptor), pointer :: FFT
            end type thread_data
            type(thread_data) :: workload(nth)

            integer th, status, mystatus

            do th = 1, nth
                status = DftiCreateDescriptor (workload(th)%FFT, DFTI_SINGLE, DFTI_COMPLEX, 2, len)
                status = DftiCommitDescriptor (workload(th)%FFT)
            end do
            ! assume x is initialized and do 2D FFTs
            !$OMP PARALLEL DO SHARED(len, x, workload) PRIVATE(mystatus)
            do th = 1, nth
                mystatus = DftiComputeForward (workload(th)%FFT, x(:, th))
            end do
            !$OMP END PARALLEL DO
            do th = 1, nth
                status = DftiFreeDescriptor (workload(th)%FFT)
            end do

            if (status .ne. 0) then
                errflag=.true.
                if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                    print *, 'Error: ', DftiErrorMessage(status)
                endif
            endif
        end subroutine test_mkl_fft2d_array_descr_main

        ! Using Parallel Mode with a Common Descriptor The following Example “Using Parallel Mode
        ! with a Common Descriptor” illustrates a parallel customer program with a common descriptor
        ! used in several threads.

        subroutine  test_mkl_fft2d_shared_descr_main(errflag)
            logical, intent(inout) :: errflag


            integer nth, len(2)
            ! 4 OMP threads, each does 2D FFT 50x100 points
            parameter (nth = 4, len = (/50, 100/))
            complex x(len(2)*len(1), nth)
            type(dfti_descriptor), pointer :: FFT

            integer th, status, mystatus

            status = DftiCreateDescriptor (FFT, DFTI_SINGLE, DFTI_COMPLEX, 2, len)
            status = DftiCommitDescriptor (FFT)
            ! assume x is initialized and do 2D FFTs
            !$OMP PARALLEL DO SHARED(len, x, FFT) PRIVATE(mystatus)
            do th = 1, nth
                mystatus = DftiComputeForward (FFT, x(:, th))
            end do
            !$OMP END PARALLEL DO
            status = DftiFreeDescriptor (FFT)

            if (status .ne. 0) then
                errflag=.true.
                if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                    print *, 'Error: ', DftiErrorMessage(status)
                endif
            endif
        end subroutine test_mkl_fft2d_shared_descr_main

    end subroutine test_mkl_fftw


    ! basic_sp_real_dft_2d.f90
    !===============================================================================
    ! Copyright 2011-2017 Intel Corporation All Rights Reserved.
    !
    ! The source code,  information  and material  ("Material") contained  herein is
    ! owned by Intel Corporation or its  suppliers or licensors,  and  title to such
    ! Material remains with Intel  Corporation or its  suppliers or  licensors.  The
    ! Material  contains  proprietary  information  of  Intel or  its suppliers  and
    ! licensors.  The Material is protected by  worldwide copyright  laws and treaty
    ! provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
    ! modified, published,  uploaded, posted, transmitted,  distributed or disclosed
    ! in any way without Intel's prior express written permission.  No license under
    ! any patent,  copyright or other  intellectual property rights  in the Material
    ! is granted to  or  conferred  upon  you,  either   expressly,  by implication,
    ! inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
    ! property rights must be express and approved by Intel in writing.
    !
    ! Unless otherwise agreed by Intel in writing,  you may not remove or alter this
    ! notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
    ! suppliers or licensors in any way.
    !===============================================================================

    ! Content:
    ! A simple example of single-precision real-to-complex out-of-place 2D
    ! FFT using Intel(R) MKL DFTI
    !
    !*****************************************************************************
    subroutine basic_sp_real_dft_2d(errflag)


        use MKL_DFTI, forget => DFTI_SINGLE, DFTI_SINGLE => DFTI_SINGLE_R
        logical, intent(inout) :: errflag
        ! Sizes of 2D transform
        integer, parameter :: N1 = 7
        integer, parameter :: N2 = 13

        ! Arbitrary harmonic used to test FFT
        integer, parameter :: H1 = 3
        integer, parameter :: H2 = 1

        ! Need single precision
        integer, parameter :: WP = selected_real_kind(6,37)

        ! Execution status
        integer :: status = 0, ignored_status

        ! Strides define data layout for real and complex domains
        integer cstrides(3), rstrides(3)

        ! Data arrays
        real(WP), allocatable :: x_real (:,:)
        complex(WP), allocatable :: x_cmplx (:,:)

        ! DFTI descriptor handle
        type(DFTI_DESCRIPTOR), POINTER :: hand

        hand => null()

        print *,"Example basic_sp_real_dft_2d"
        print *,"Forward-Backward single-precision real out-of-place 2D transform"
        print *,"Configuration parameters:"
        print *,"DFTI_PRECISION              = DFTI_SINGLE"
        print *,"DFTI_FORWARD_DOMAIN         = DFTI_REAL"
        print *,"DFTI_DIMENSION              = 2"
        print '(" DFTI_LENGTHS                = /"I0","I0"/" )', N1, N2
        print *,"DFTI_PLCEMENT               =  DFTI_NOT_INPLACE"
        print *,"DFTI_CONJUGATE_EVEN_STORAGE = DFTI_COMPLEX_COMPLEX"

        print *,"Create DFTI descriptor for real transform"
        status = DftiCreateDescriptor(hand, DFTI_SINGLE, DFTI_REAL, 2, (/N1, N2/))
        if (0 /= status) goto 999

        print *,"Set out-of-place"
        status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        if (0 /= status) goto 999

        print *,"Set CCE storage"
        status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE,                   &
            &                    DFTI_COMPLEX_COMPLEX)
        if (0 /= status) goto 999

        rstrides = [0, 1, N1]
        cstrides = [0, 1, INT(N1/2.0)+1]

        print '(" Set input  strides = "3(I0:", "))', rstrides
        status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides)
        if (0 /= status) goto 999

        print '(" Set output strides = "3(I0:", "))', cstrides
        status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides)
        if (0 /= status) goto 999

        print *,"Commit DFTI descriptor"
        status = DftiCommitDescriptor(hand)
        if (0 /= status) goto 999

        print *,"Allocate data arrays"
        allocate ( x_real(N1, N2), STAT = status)
        if (0 /= status) goto 999
        allocate ( x_cmplx(INT(N1/2.0)+1, N2), STAT = status)
        if (0 /= status) goto 999

        print *,"Initialize data for real-to-complex FFT"
        call init_r(x_real, N1, N2, H1, H2)

        print *,"Compute forward transform"
        status = DftiComputeForward(hand, x_real(:,1), x_cmplx(:,1))
        if (0 /= status) goto 999

        print *,"Verify the result"
        status = verify_c(x_cmplx, N1, N2, H1, H2)
        if (0 /= status) goto 999


        print *,"Reconfigure DFTI descriptor for backward transform"

        print '(" Set input  strides = "3(I0:", "))', cstrides
        status = DftiSetValue(hand, DFTI_INPUT_STRIDES, cstrides)
        if (0 /= status) goto 999

        print '(" Set output strides = "3(I0:", "))', rstrides
        status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, rstrides)
        if (0 /= status) goto 999

        print *,"Recommit DFTI descriptor"
        status = DftiCommitDescriptor(hand)
        if (0 /= status) goto 999

        print *,"Initialize data for complex-to-real FFT"
        call init_c(x_cmplx, N1, N2, H1, H2)

        print *,"Compute forward transform"
        status = DftiComputeBackward(hand, x_cmplx(:,1), x_real(:,1))
        if (0 /= status) goto 999

        print *,"Verify the result"
        status = verify_r(x_real, N1, N2, H1, H2)
        if (0 /= status) goto 999

100     continue

        print *,"Release the DFTI descriptor"
        ignored_status = DftiFreeDescriptor(hand)

        if (allocated(x_real) .or. allocated(x_cmplx)) then
            print *,"Deallocate data arrays"
        endif
        if (allocated(x_real))     deallocate(x_real)
        if (allocated(x_cmplx))    deallocate(x_cmplx)

        if (status == 0) then
            print *, "TEST PASSED"
           return !  call exit(0)
        else
            print *, "TEST FAILED"
            errflag=.true.
            return ! call exit(1)
        end if

999     print '("  Error, status = ",I0)', status
        goto 100

    contains

        ! Compute mod(K*L,M) accurately
        pure real(WP) function moda(k,l,m)
            integer, intent(in) :: k,l,m
            integer*8 :: k8
            k8 = k
            moda = real(mod(k8*l,m),WP)
        end function moda

        ! Initialize x(:,:) to harmonic H
        subroutine init_r(x, N1, N2, H1, H2)
            integer N1, N2, H1, H2
            real(WP) :: x(:,:)

            integer k1, k2
            real(WP), parameter:: TWOPI = 6.2831853071795864769_WP

            forall (k1=1:N1, k2=1:N2)
                x(k1,k2) = 2.0_WP * cos( TWOPI * ( moda(H1,k1-1,N1) / real(N1,WP) &
                    &    +                           moda(H2,k2-1,N2) / real(N2,WP))) / real((N1*N2),WP)
            end forall
            if (mod(2*(N1-H1),N1)==0 .and. mod(2*(N2-H2),N2)==0) then
                x(1:N1,1:N2) = x(1:N1,1:N2) / 2
            end if
        end subroutine init_r

        ! Initialize x(:,:) to produce unit peak at x(H1,H2)
        subroutine init_c(x, N1, N2, H1, H2)
            integer N1, N2, H1, H2
            complex(WP) :: x(:,:)

            complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)
            integer k1,k2

            forall (k1=1:N1/2+1, k2=1:N2)
                x(k1,k2) = exp( -I_TWOPI * ( moda(H1,k1-1,N1) / real(N1,WP) &
                    &    +                          moda(H2,k2-1,N2) / real(N2,WP))) / cmplx(N1*N2)
            end forall
        end subroutine init_c


        ! Verify that x(:,:) has unit peak at (H1,H2)
        integer function verify_c(x, N1, N2, H1, H2)
            integer N1, N2, H1, H2
            complex(WP) :: x(:,:)

            integer k1, k2
            real(WP) err, errthr, maxerr
            complex(WP) res_exp, res_got

            ! Note, this simple error bound doesn't take into account error of
            ! input data
            errthr = 2.5 * log(real(N1*N2,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
            print '("  Check if err is below errthr " G10.3)', errthr

            maxerr = 0.0_WP
            do k2 = 1, N2
                do k1 = 1, N1/2+1
                    if (mod(k1-1-H1,N1)==0 .and. mod(k2-1-H2,N2)==0) then
                        res_exp = 1.0_WP
                    else if (mod(-k1+1-H1,N1)==0 .and. mod(-k2+1-H2,N2)==0) then
                        res_exp = 1.0_WP
                    else
                        res_exp = 0.0_WP
                    end if
                    res_got = x(k1, k2)
                    err = abs(res_got - res_exp)
                    maxerr = max(err,maxerr)
                    if (.not.(err < errthr)) then
                        print '("  x("I0","I0"):"$)', k1,k2
                        print '(" expected ("G14.7","G14.7"),"$)', res_exp
                        print '(" got ("G14.7","G14.7"),"$)', res_got
                        print '(" err "G10.3)', err
                        print *,"  Verification FAILED"
                        verify_c = 100
                        return
                    end if
                end do
            end do
            print '("  Verified,  maximum error was " G10.3)', maxerr
            verify_c = 0
        end function verify_c

        ! Verify that x(:,:) is unit peak at x(H1,H2)
        integer function verify_r(x, N1, N2, H1, H2)
            integer N1, N2, H1, H2
            real(WP) :: x(:,:)

            integer k1, k2
            real(WP) err, errthr, maxerr, res_exp, res_got

            ! Note, this simple error bound doesn't take into account error of
            ! input data
            errthr = 2.5 * log(real(N1*N2,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
            print '("  Check if err is below errthr " G10.3)', errthr

            maxerr = 0.0_WP
            do k2 = 1, N2
                do k1 = 1, N1
                    if (mod(k1-1-H1, N1)==0 .AND. mod(k2-1-H2,N2)==0) then
                        res_exp = 1.0_WP
                    else
                        res_exp = 0.0_WP
                    end if
                    res_got = x(k1, k2)
                    err = abs(res_got - res_exp)
                    maxerr = max(err,maxerr)
                    if (.not.(err < errthr)) then
                        print '("  x("I0","I0"):"$)', k1, k2
                        print '(" expected "G14.7","$)', res_exp
                        print '(" got "G14.7","$)', res_got
                        print '(" err "G10.3)', err
                        print *,"  Verification FAILED"
                        verify_r = 100
                        return
                    end if
                end do
            end do
            print '("  Verified,  maximum error was " G10.3)', maxerr
            verify_r = 0
        end function verify_r
    end subroutine basic_sp_real_dft_2d

    ! config_thread_limit.f90
    !===============================================================================
    ! Copyright 2012-2017 Intel Corporation All Rights Reserved.
    !
    ! The source code,  information  and material  ("Material") contained  herein is
    ! owned by Intel Corporation or its  suppliers or licensors,  and  title to such
    ! Material remains with Intel  Corporation or its  suppliers or  licensors.  The
    ! Material  contains  proprietary  information  of  Intel or  its suppliers  and
    ! licensors.  The Material is protected by  worldwide copyright  laws and treaty
    ! provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
    ! modified, published,  uploaded, posted, transmitted,  distributed or disclosed
    ! in any way without Intel's prior express written permission.  No license under
    ! any patent,  copyright or other  intellectual property rights  in the Material
    ! is granted to  or  conferred  upon  you,  either   expressly,  by implication,
    ! inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
    ! property rights must be express and approved by Intel in writing.
    !
    ! Unless otherwise agreed by Intel in writing,  you may not remove or alter this
    ! notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
    ! suppliers or licensors in any way.
    !===============================================================================

    ! Content:
    ! An example of using DFTI_THREAD_LIMIT configuration parameter.
    ! The parameter specifies maximum number of OpenMP threads FFT can use.
    !
    ! Values:
    !   0 (default) = use number of threads specified by
    !                 mkl_[domain_]set_num_threads()
    !   Any positive integer N = use not more than N threads
    !
    !*****************************************************************************
    subroutine config_thread_limit(errflag)
#if defined(_OPENMP)
        use omp_lib
#endif
        logical, intent(inout):: errflag

        ! Need double precision
        integer, parameter :: WP = selected_real_kind(15,307)

        ! Number of parallel user threads
#if defined(_OPENMP)
        integer, parameter :: NUT = 2
#endif

        integer :: failed = 0
        integer :: thr_failed, me, team

        print *,"Example config_thread_limit"

        ! Enable nested parallel OpenMP sections (maybe oversubscribed)
        ! Use kind(4) literals because some compilers fail on this in i8 mode
#if defined(_OPENMP)
        call omp_set_nested(.true._4)
        call omp_set_dynamic(.false._4)
#endif

        ! Enable threading of Intel(R) MKL called from OpenMP parallel sections
        call mkl_set_dynamic(0)

#if defined(_OPENMP)
        print '(" Run parallel FFTs on "I0" parallel threads")', NUT
        !$omp parallel num_threads(NUT) default(shared) private(thr_failed, me)
#else
        print '(" Run parallel FFT on a single thread")'
#endif

#if defined(_OPENMP)
        me = omp_get_thread_num()
        team = omp_get_num_threads()
#else
        me = 0
        team = 1
#endif

        if (me == 0) then
            print '(" Thread "I0": parallel team is "I0" threads")', me, team
        endif

        if (me == 0) then
            thr_failed = run_dft(me, 2, 100,200,300, -1,-2,-3)
        else
            thr_failed = run_dft(me, 3, 200,300,100, -1,-2,-3)
        endif
        if (0 /= thr_failed) failed = thr_failed

#if defined(_OPENMP)
        !$omp end parallel
#endif

        if (failed == 0) then
            print *,"TEST PASSED"
            return ! call exit(0)
        else
            print *,"TEST FAILED"
            errflag=.true.
            return ! call exit(1)
        endif

    contains

        integer function run_dft(tid,tlimit,n1,n2,n3,h1,h2,h3)

            use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

            integer :: tid              ! Id of this thread
            integer :: tlimit           ! Thread limit
            integer :: N1,N2,N3         ! Sizes of 3D transform
            integer :: H1,H2,H3         ! Arbitrary harmonic used to test the FFT


            ! Execution status
            integer :: status = 0, ignored_status

            ! DFTI descriptor handle
            type(DFTI_DESCRIPTOR), POINTER :: hand

            ! Data array
            complex(WP), allocatable :: x (:,:,:)

            ! Temporary variable
            integer :: tl

            print '(" Thread "I0": Create DFTI descriptor for "I0"x"I0"x"I0" FFT")', &
                & tid, N1, N2, N3
            status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 3, [N1,N2,N3])
            if (0 /= status) goto 999

            print '(" Thread "I0": Set number of user threads "I0)', tid, tlimit
            status = DftiSetValue(hand, DFTI_THREAD_LIMIT, tlimit)
            if (0 /= status) goto 999

            ! If tlimit > 1 check if we linked with sequential MKL
            if (tlimit > 1) then
                ! Get thread limit of uncommitted descriptor
                status = DftiGetValue(hand, DFTI_THREAD_LIMIT, tl)
                if (0 /= status) goto 999
                print '(" Thread "I0": uncommitted descriptor thread limit "I0)', tid, tl
            endif

            print '(" Thread "I0": Commit DFTI descriptor")', tid
            status = DftiCommitDescriptor(hand)
            if (0 /= status) goto 999

            ! Get thread limit of committed descriptor
            status = DftiGetValue(hand, DFTI_THREAD_LIMIT, tl)
            if (0 /= status) goto 999
            print '(" Thread "I0": committed descriptor thread limit "I0)', tid, tl

            print '(" Thread "I0": Allocate data array")', tid
            allocate ( x(N1, N2, N3), STAT = status)
            if (0 /= status) goto 999

            print '(" Thread "I0": Initialize input for forward transform")', tid
            call init(x, N1, N2, N3, H1, H2, H3)

            print '(" Thread "I0": Compute forward transform")', tid
            status = DftiComputeForward(hand, x(:, 1, 1))
            if (0 /= status) goto 999

            print '(" Thread "I0": Verify the result")', tid
            status = verificate(x, N1, N2, N3, H1, H2, H3)
            if (0 /= status) goto 999

100         continue

            print '(" Thread "I0": Release the DFTI descriptor")', tid
            ignored_status = DftiFreeDescriptor(hand)

            if (allocated(x)) then
                print '(" Thread "I0": Deallocate data array")', tid
                deallocate(x)
            endif

            if (status == 0) then
                print '(" Thread "I0": Subtest Passed")', tid
            else
                print '(" Thread "I0": Subtest Failed")', tid
            endif

            run_dft = status
            return

999         print '(" Thread "I0":  Error, status = ",I0)', tid, status
            goto 100

        end function run_dft

        ! Compute mod(K*L,M) accurately
        pure real(WP) function moda(k,l,m)
            integer, intent(in) :: k,l,m
            integer*8 :: k8
            k8 = k
            moda = real(mod(k8*l,m),WP)
        end function moda

        ! Initialize arrays with harmonic /H1, H2, H3/
        subroutine init(x, N1, N2, N3, H1, H2, H3)
            integer N1, N2, N3, H1, H2, H3
            complex(WP) :: x(:,:,:)

            integer k1, k2, k3
            complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)

            forall (k1=1:N1, k2=1:N2, k3=1:N3)
                x(k1,k2,k3) = exp( I_TWOPI * ( &
                    &    moda(  k1-1,H1, N1)/N1 &
                    &    + moda(k2-1,H2, N2)/N2 &
                    &    + moda(k3-1,H3, N3)/N3 )) / (N1*N2*N3)
            end forall
        end subroutine init

        ! Verify that x(:,:,:) are unit peaks at /H1, H2, H3/
        integer function verificate(x, N1, N2, N3, H1, H2, H3)
            integer N1, N2, N3, H1, H2, H3
            complex(WP) :: x(:,:,:)

            integer k1, k2, k3
            real(WP) err, errthr, maxerr
            complex(WP) res_exp, res_got

            ! Note, this simple error bound doesn't take into account error of
            ! input data
            errthr = 5.0 * log(real(N1*N2*N3,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
            print '("  Check if err is below errthr " G10.3)', errthr

            maxerr = 0.0_WP
            do k3 = 1, N3
                do k2 = 1, N2
                    do k1 = 1, N1
                        if (mod(k1-1-H1, N1)==0 .AND.                                &
                            &        mod  (k2-1-H2, N2)==0 .AND.                     &
                            &        mod  (k3-1-H3, N3)==0) then
                            res_exp = 1.0_WP
                        else
                            res_exp = 0.0_WP
                        end if
                        res_got = x(k1,k2,k3)
                        err = abs(res_got - res_exp)
                        maxerr = max(err,maxerr)
                        if (.not.(err <= errthr)) then
                            print '("  x("I0","I0","I0"):"$)', k1, k2, k3
                            print '(" expected ("G24.17","G24.17"),"$)', res_exp
                            print '(" got ("G24.17","G24.17"),"$)', res_got
                            print '(" err "G10.3)', err
                            print *," Verification FAILED"
                            verificate = 100
                            return
                        end if
                    end do
                end do
            end do
            print '("  Verified,  maximum error was " G10.3)', maxerr
            verificate = 0
        end function verificate

    end subroutine config_thread_limit

#endif
end module simple_intel_blas
