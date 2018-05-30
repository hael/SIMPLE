program simple_test_openacc
! #include "openacc_lib.h"
    use openacc
    use simple_oacc_vecadd
    use simple_oacc_omp
    implicit none

#ifndef _OPENACC
    write(*,*) " _OPENACC not defined use -fopenacc in FFLAGS or enable USE_OPENACC in cmake build"
#else
!!  warning "C Preprocessor got here!"  _OPENACC
    write(*,'(a,1x,i0)')"OpenACC preprocessor version:",  _OPENACC
#endif


#ifdef  _OPENACC
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
    call test_oacc_omp
    call test_oacc_reduction
#else
    print *,' simple_test_openacc OpenACC is disabled '
#endif
contains

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

end program simple_test_openacc
