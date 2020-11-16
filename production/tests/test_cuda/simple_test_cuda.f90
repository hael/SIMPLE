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
! 06 Oct 2017 - Initial Version -- PGI
! 07 May 2018 - FortCUDA implementaiton
!------------------------------------------------------------------------------
program simple_test_cuda
include 'simple_lib.f08'
use simple_image,            only: image
use gnufor2
use CUDA
use simple_cuda
use simple_timer_cuda
use simple_cuda_kernels
use, intrinsic :: ISO_C_BINDING
implicit none

type (timer_cuda) :: ctimer
type (cudaEvent_t)  :: ev1,ev2
integer(timer_int_kind) :: t1
integer(c_int):: runtimeVersion,driverVersion,pValue, deviceCount,free,total
integer (KIND(cudaLimitStackSize))::limitStackSize= cudaLimitStackSize
integer (KIND(cudaSuccess)) :: err
logical :: error_found

    error_found=.false.

    print *," CUDA Runtime functions "
    call cuda_query_version
    call cuda_query_driver_version
    call cuda_query_device_count
    call cuda_query_thread_limit
    call cuda_thread_synchronize
    call cuda_print_mem_info


    call check_cuda_device
    !  call set_cuda_device(0)

    write (*,'(A)') 'TESTING CUDA FORTRAN PRECISION '
    call test_cuda_precision( error_found )

    write (*,'(A)') 'SIMPLE_CUDA timer test'

    write (*,*)""
    write (*,'(A)') 'TESTING CUDA FORTRAN KERNELS'
    !call test_FortCUDA_kernels(0.)
    !call test_fortran_mul1dComplex_kernels
    call test_fortran_squaremul2dComplex_kernels
    call test_fortran_squaremuladd2dComplex_kernels
    !call test_fortran_mul2dComplex_kernels


contains

    subroutine test_cuda_precision(flag)
        implicit none
        logical, intent(inout):: flag
        write(logfhandle,"(a)") '  CUDA Query: Test precision'
        call test_precision(flag)
        write(logfhandle,"(a)") '  CUDA Query: Sum accuracy'
        call sum_accuracy(flag)
    end subroutine test_cuda_precision

        !> Floating-point precision test
        subroutine test_precision(flag)
            logical, intent(inout):: flag

            real, parameter :: x=REAL(Z'3F1DC57A'), y=REAL(Z'3F499AA3')
            real :: dist
            double precision:: x_dp, y_dp, dist_dp


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

            write(logfhandle, "('Summming ',i10,' elements of magnitude ',f3.1)") N,7.
            write(logfhandle, "('Sum with intrinsic function       =',f12.1,'   Error=', f12.1)")  &
                sum_intrinsic, 7.*N-sum_intrinsic
            write(logfhandle, "('Recursive sum with SP accumulator =',f12.1,'   Error=', f12.1)")  sum_cpu, 7.*N-sum_cpu
            write(logfhandle, "('Recursive sum with DP accumulator =',f12.1,'   Error=', f12.1)")  sum_cpu_dp, 7.*N-sum_cpu_dp
            write(logfhandle, "('Pairwise sum in SP                =',f12.1,'   Error=', f12.1)")  sum_pairwise, 7.*N-sum_pairwise
            write(logfhandle, "('Compensated sum in SP             =',f12.1,'   Error=', f12.1)")  sum_kahan, 7.*N-sum_kahan

            deallocate(x)
    end subroutine sum_accuracy

    subroutine simple_test_cuda_timer

    ctimer = timer_cuda()

    write (*,'(A)') 'TESTING CUDAFOR TIMING'
    ! call ctimer%nowCU()
    t1=tic()
    ev1=ctimer%ticU()
    ev2=ctimer%ticU()
    write (*,'(A)') 'SIMPLE_CUDA timer CPU/CUDA'
    print *, " Simple_timer ", toc(t1)
    print *, " CUDA Event timer 1", ctimer%tocU(ev1)
    print *, " CUDA Event timer 2",ctimer%tocU(ev2)
    call ctimer%kill_()

    end subroutine simple_test_cuda_timer


end program simple_test_cuda
