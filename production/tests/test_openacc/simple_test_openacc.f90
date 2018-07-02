program simple_test_openacc
    include 'simple_lib.f08'
    !    use omp_lib
    use openacc
    use simple_oacc_vecadd
    use simple_oacc_omp
    use simple_rbc
    implicit none

!#ifndef OPENACC
!    write(*,*) " _OPENACC not defined use -fopenacc in FFLAGS or enable USE_OPENACC in cmake build"
!#else
    !!  warning "C Preprocessor got here!"  _OPENACC
!    write(*,'(a,1x,i0)')"OpenACC preprocessor version:",  _OPENACC
!#endif


#ifdef OPENACC
    print *,' simple_test_openacc OpenACC is enabled '
    print *,' OpenACC version : ', openacc_version
    
     call test_oacc_basics
     call test_oacc_vecadd
     call test_oacc_vecadd_nocopy

     call test_oacc_omp_matrixmul
     call test_oacc_omp_matrixmul2
     call test_gang_static_addition
     call test_nested1
     call test_nested2
     !call test_oacc_omp
!     call test_oacc_reduction
     
     call run_rbc_serial
     call run_rbc_omp
     call run_rbc_oacc

     ! call test_matrix_mul_omp
!      call test_matrix_mul_openacc
!     call  saxpy_test_omp
!     call  saxpy_test_acc

!     call test_pi_omp
!     call test_pi_openacc
    

#else
    print *,' simple_test_openacc OpenACC is disabled '
#endif
contains

    subroutine test_openacc_runtime
        integer :: i, n,acc_device_nvidia
        real :: a
        real, allocatable :: x(:), y(:)
        n=1000
        acc_device_nvidia=1
        a=atan2(1.0,1.0)
        ! System setup routines

!        call acc_init(acc_device_nvidia)
!        call acc_set_device_type(acc_device_nvidia)
!        call acc_set_device_num(acc_device_nvidia)
        ! Synchronization routines
        ! call acc_async_wait(i)
        allocate(x(n),y(n))
        !$acc kernels
        do i = 1, n
            x(i) = 1.0
            y(i) = 2.0
        end do
        call acc_async_wait(i)
        do i = 1, n
            y(i) = y(i) + a * x(i)
        end do
        !$acc end kernels
        call acc_async_wait_all()

    end subroutine test_openacc_runtime



    subroutine test_oacc_basics

        integer, parameter :: N=10, M=500, P=30
        integer :: i, j, k
        real    :: particles(N,M,P), total
#if defined(PGI)|| (defined(GNU) && (__GNUC__ >= 6 ))
        particles = 1.0
        if (openacc_version >= 201306) then
            !$acc kernels
            do i=1,N
                do j=1,M
                    do k=1,P
                        particles(i,j,k) = particles(i,j,k) * 2.0
                    end do
                end do
            end do
            !$acc end kernels

            print *,'test_oacc_basics:kernels: '
            print *,'                sum should be: ', N*M*P*2, 'sum is: ', sum(particles)

            !$acc parallel loop
            do i=1,N
                do j=1,M
                    do k=1,P
                        particles(i,j,k) = particles(i,j,k) * 2.0
                    end do
                end do
            end do
            !$acc end parallel loop

            print *, 'test_oacc_basics:parallel loop: '
            print *,'                 sum should be: ', N*M*P*4, 'sum is: ', sum(particles)

            !$acc parallel loop num_gangs(32) vector_length(128)
            do i=1,N
                do j=1,M
                    do k=1,P
                        particles(i,j,k) = particles(i,j,k) * 2.0
                    end do
                end do
            end do
            !$acc end parallel loop

            print *, 'test_oacc_basics:parallel loop num_gangs(32) vector_length(128): '
            print *,'                 sum should be: ', N*M*P*8, 'sum is: ', sum(particles)

            !$acc parallel loop collapse(3)
            do i=1,N
                do j=1,M
                    do k=1,P
                        particles(i,j,k) = particles(i,j,k) * 2.0
                    end do
                end do
            end do
            !$acc end parallel loop

            print *, 'test_oacc_basics:parallel loop collapse(3): '
            print *,'                 sum should be: ', N*M*P*16, 'sum is: ', sum(particles)

            particles = 1.0

            total = 0.
            !$acc parallel loop collapse(3) reduction(+:total)
            do i=1,N
                do j=1,M
                    do k=1,P
                        particles(i,j,k) = particles(i,j,k) * 2.0
                        total = total + particles(i,j,k)
                    end do
                end do
            end do
            !$acc end parallel loop

            print *, 'test_oacc_basics:parallel loop collapse(3) reduction(+:total): '
            print *,'                 sum should be: ', N*M*P*2, 'sum is: ', total
        else
            print *, 'test_oacc_basics:parallel loop: sorry, unimplemented: directive not yet implemented in your gfortran version'
        endif
#endif

    end subroutine test_oacc_basics


    subroutine saxpy_test_omp
        !$use omp_lib
        implicit none
!        integer, parameter :: dp = selected_real_kind(15)
!        integer, parameter :: ip = selected_int_kind(15)
        integer(8) :: i,n
        real(8),dimension(:),allocatable :: x, y
        real(8) :: a,start_time, end_time
        n=500000000
        allocate(x(n),y(n))
        !$omp parallel sections
        !$omp section
        x = 1.0
        !$omp section
        y = 1.0
        !$omp end parallel sections
        a = 2.0
         start_time = omp_get_wtime()
        !$omp parallel do default(shared) private(i)
        do i = 1, n
            y(i) = y(i) + a * x(i)
        end do
        !$omp end parallel do
         end_time = omp_get_wtime()
        deallocate(x,y)
        print *, 'SAXPY Time (OpenMP): ', end_time - start_time, 'in secs'
    end subroutine saxpy_test_omp

    subroutine saxpy_test_acc
        use omp_lib
        implicit none
        integer(8) :: i,n
        real(8),dimension(:),allocatable :: x, y
        real(8) :: a,start_time, end_time
        n=500000000
        allocate(x(n),y(n))
        a = 2.0
        !$acc data create(x,y) copyin(a)
        !$acc parallel
        x(:) = 1.0
        !$acc end parallel
        !$acc parallel
        y(:) = 1.0
        !$acc end parallel
        start_time = omp_get_wtime()
        !$acc parallel 
        do i = 1, n
            y(i) = y(i) + a * x(i)
        end do
      !  y=a*x+y
        !$acc end parallel
        end_time = omp_get_wtime()
        !$acc end data
        deallocate(x,y)
        print *,'SAXPY Time (OpenACC): ', end_time - start_time, 'in secs'
    end subroutine saxpy_test_acc

    subroutine test_pi_omp
        !$use omp_lib
        implicit none
        integer, parameter :: dp=selected_real_kind(14)
        integer, parameter :: ip=selected_int_kind(15)
        integer :: i
        integer, parameter :: n=10000000000
        real(dp) :: x,pi,sum,step,start_time,end_time
        integer, dimension(8) :: value
        sum = 0d0
        step = 1.d0/float(n)

        ! call date_and_time(VALUES=value)
        ! start_time = float(value(6)*60) + float(value(7)) + float(value(8))/1000d0
        start_time = omp_get_wtime()
        !$omp parallel do private(i,x) reduction(+:sum)
        do i = 0, n
            x = (i + 0.5d0) * step
            sum = sum + 4.d0 / (1.d0 + x ** 2)
        end do
        !$omp end parallel do
        pi = step * sum
        ! call date_and_time(VALUES=value)
        ! end_time = float(value(6)*60) + float(value(7)) + float(value(8))/1000d0
        ! if ( start_time > end_time ) end_time = end_time + 3600d0
        end_time = omp_get_wtime()

        print '(a,f17.15)', "pi = ", pi
        print '(a,f9.3,a)', "OpenMP time to compute =",end_time - start_time, " seconds"
    end subroutine test_pi_omp

    subroutine test_pi_openacc
        !$use openacc
        implicit none
        integer, parameter :: dp=selected_real_kind(14)
        integer, parameter :: ip=selected_int_kind(15)
        integer :: i
        integer, parameter :: n=10000000000
        real(dp) :: x,pi,sumpi,step,start_time,end_time
        integer, dimension(8) :: value
        sumpi = 0d0
        step = 1.d0/float(n)
        ! call date_and_time(VALUES=value)
        ! start_time = float(value(6)*60) + float(value(7)) + float(value(8))/1000d0
        start_time = omp_get_wtime()
        !$acc data copyin(step) copyout(sumpi)
        !$acc parallel loop private(x) reduction(+:sumpi)
        do i = 0, n
            x = (i + 0.5d0) * step
            sumpi = sumpi + 4.d0 / (1.d0 + x ** 2)
        end do
        !$acc end parallel loop
        pi = step * sumpi
        !$acc end data

        ! call date_and_time(VALUES=value)
        ! end_time = float(value(6)*60) + float(value(7)) + float(value(8))/1000d0
        ! if ( start_time > end_time ) end_time = end_time + 3600d0
        end_time = omp_get_wtime()

        print '(a,f17.15)', "pi = ", pi
        print '(a,f9.3,a)', "OpenACC time to compute =",end_time - start_time, " seconds"
    end subroutine test_pi_openacc



    subroutine test_matrix_mul_omp
        !$use omp_lib
      integer :: i,j,k
        integer, parameter :: nra=1500, nca=2000, ncb=1000
        real(8) :: a(nra,nca) , b(nca,ncb) , c(nra,ncb)
        real(8) :: flops, tmp, mean
        real(8) :: init_time, start_time, end_time
        integer :: c1, c2, c3, cr
        integer, dimension(8) :: value
        flops = 2.d0 * real(nra) * real(nca) * real(ncb)
        init_time = omp_get_wtime()
  
        c = 0.d0
        do i = 1,nra
            do j = 1,nca
                a(i,j) = i + j
            end do
        end do
        do i = 1,nca
            do j = 1,ncb
                b(i,j) = i * j
            end do
        end do
  
        start_time = omp_get_wtime()
        !$omp parallel do private(i,j,k,tmp)
        do j = 1, nca
            do k = 1, ncb
                tmp = 0.d0
                !$omp parallel do reduction(+:tmp)
                do i = 1, nra
                    tmp = tmp + a(i,j) * b(j,k)
                end do
                !$omp end parallel do
                c(i,k) = tmp
            end do
        end do
        !$omp end parallel do
  
        end_time = omp_get_wtime()
mean = sum(sum(c,2),1)/ real(nra*nca)
        print '(a,f6.3,a,f6.3,a,f7.3,a,f7.3)', 'OpenMP Init Time: ', start_time - init_time, &
            ' Calc Time: ', end_time - start_time, &
            ' GFlops: ', 1d-9 * flops/(end_time - start_time), ' sum ', mean
    end subroutine test_matrix_mul_omp

    subroutine test_matrix_mul_openacc
        implicit none
     !   integer, parameter :: dp = selected_real_kind(8)
        integer :: i,j,k
        integer, parameter :: nra=1500, nca=2000, ncb=1000
        real(8) :: a(nra,nca) , b(nca,ncb) , c(nra,ncb)
        real(8) :: flops, tmp, mean
        real(8) :: init_time, start_time, end_time
        integer :: c1, c2, c3, cr
        integer, dimension(8) :: value
        flops = 2d0 * real(nra) * real(nca) * real(ncb)
        init_time = omp_get_wtime()
        !$acc data create(a,b,c)

        c = 0.d0
        do i = 1,nra
            do j = 1,nca
                a(i,j) = i + j
            end do
        end do
        do i = 1,nca
            do j = 1,ncb
                b(i,j) = i * j
            end do
        end do
        start_time = omp_get_wtime()
        !$acc parallel loop private(i,j,k,tmp)
        do j = 1, nca
            do k = 1, ncb
                tmp = 0.d0
                !$acc loop reduction(+:tmp)
                do i = 1, nra
                    tmp = tmp + a(i,j) * b(j,k)
                end do
                c(i,k) = tmp
            end do
        end do
        !$acc end parallel loop
        !$acc end data
         end_time = omp_get_wtime()
         mean = sum(sum(c,2),1)/ real(nra*nca)
        print '(a,f6.3,a,f6.3,a,f7.3,a,f7.3)', 'OpenACC Init Time: ', start_time - init_time, &
            ' Calc Time: ', end_time - start_time, &
            ' GFlops: ', 1d-9 * flops/(end_time - start_time), ' sum ', mean
    end subroutine test_matrix_mul_openacc

end program simple_test_openacc
