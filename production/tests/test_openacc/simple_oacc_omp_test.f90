! 2.10.2. Multi-Threaded Program Utilizing Multiple Devices
! This simple example shows how to run a multi-threaded host program that utilizes
! multiple devices.

! The program starts by having each thread call acc_set_device_num so each thread
! will use a different GPU. Within the computational OpenMP parallel region, each thread
! copies the data it needs to its GPU and proceeds.
module simple_oacc_omp
    use openacc
   use omp_lib
! #include "openacc_lib.h"

    implicit none

    public:: test_oacc_omp, test_oacc_omp_matrixmul, test_oacc_omp_matrixmul2, test_gang_static_addition, test_oacc_reduction, test_nested1, test_nested2

contains
    subroutine test_oacc_omp
        integer, parameter :: N = 10000
        real(8) x(N), y(N), z
        integer, allocatable :: offs(:)
        real(8), allocatable :: zs(:)
        real(8) ddot

        integer nthr,i,nsec,j
        if (openacc_version >= 201306) then
            ! Max at 4 threads for now
            nthr = omp_get_max_threads()
            if (nthr .gt. 4) nthr = 4
            call omp_set_num_threads(nthr)
            ! Run on host
            call random_number(x)
            call random_number(y)
            z = dot_product(x,y)  !! BLAS ddot(N,x,1,y,1)
            print *,"test_oacc_omp: Host Serial",z
            ! Attach each thread to a device
            !$omp PARALLEL private(i)
            i = omp_get_thread_num()
            call acc_set_device_num(i, acc_device_nvidia)
            !$omp end parallel
            ! Break up the array into sections
            nsec = N / nthr
            allocate(offs(nthr),zs(nthr))
            offs = (/ (i*nsec,i=0,nthr-1) /)
            zs = 0.0d0
            ! Decompose the problem across devices
#ifdef PGI
            !$omp parallel private(i,j,z)
#else
            !acc parallel private(i,j,z)
#endif
            do
                i = omp_get_thread_num() + 1
                z = 0.0d0

               !$acc kernels loop copyin(x(offs(i)+1:offs(i)+nsec),y(offs(i)+1:offs(i)+nsec))
                do j = offs(i)+1, offs(i)+nsec
                    z = z + x(j) * y(j)
                end do

                zs(i) = z
            enddo
#ifdef PGI
            !$omp end  parallel
#else
            !acc  end parallel
#endif

            z = sum(zs)
            print *,"Multi-Device Parallel",z

        else
            print *,"test_oacc_omp: OPENACC not supported. ", openacc_version
        endif
    end subroutine test_oacc_omp


    subroutine test_oacc_omp_matrixmul ()
        integer, parameter :: n_size=1000
        real, dimension(:,:) :: a(n_size,n_size)
        real, dimension(:,:) :: b(n_size,n_size)
        real, dimension(:,:) :: c(n_size,n_size)
        integer i,j,k
        print *,"test_oacc_omp_matrixmul: ACC matmul"
#if defined(_OPENACC)
        if(openacc_version >= 201306)then
            !     Initialize matrices (values differ from C version)
            do i=1, n_size
                do j=1, n_size
                    a(i,j) = i + j;
                    b(i,j) = i - j;
                    c(i,j) = 0.;
                enddo
            enddo

            !$acc data copyin(a,b) copy(c)
            !$acc kernels loop
            !     Compute matrix multiplication.
            do i=1, n_size
                do j=1, n_size
                    do k = 1, n_size
                        c(i,j) = c(i,j) + a(i,k) * b(k,j)
                    enddo
                enddo
            enddo
            !$acc end data
            print *,"test_oacc_omp_matrixmul: Completed successfully"
        else
            print *,"test_oacc_omp_matrixmul: OPENACC version not supported ", openacc_version
        end if
#else
        print *,"test_oacc_omp_matrixmul: OPENACC not supported."
#endif

    end subroutine test_oacc_omp_matrixmul

    subroutine test_oacc_omp_matrixmul2
        integer :: i, j, k, m, n, p
        integer :: t1, t2, dt, count_rate, count_max
        integer :: td(4)=(/16,64,256,1024/)
        real, allocatable, dimension(:,:) :: a, b, c
        real :: tmp, secs
 print *,"test_oacc_omp_matrixmul2: "
        call system_clock(count_max=count_max, count_rate=count_rate)

        n = 100
        allocate( a(n,n), b(n,n), c(n,n) )

        !$acc data create(a,b) copyout(c)
        do m = 1,4

            call system_clock(t1)
            do p = 1,td(m)           ! td={16,64,256,1024}

                ! Initialize matrices
                !$acc kernels
                do j=1,n
                    do i=1,n
                        a(i,j) = real(i + j)
                        b(i,j) = real(i - j)
                    enddo
                enddo

                ! Compute matrix multiplication.
                do j=1,n
                    do i=1,n
                        tmp = 0.0
                        do k=1,n
                            tmp = tmp + a(i,k) * b(k,j)
                        enddo
                        c(i,j) = tmp
                    enddo
                enddo
                !$acc end kernels

            enddo

            call system_clock(t2)
            dt = t2-t1
            secs = real(dt)/real(count_rate)
            write(*,*)"For domain_num=",td(m)," wall clock time is ", secs

        enddo
        !$acc end data
        print *, C(1,1), C(n,n)
        deallocate(a, b, c)
 print *,"test_oacc_omp_matrixmul2: Completed"
    end subroutine test_oacc_omp_matrixmul2


    !> libgomp.oacc-fortran/parallel-reduction.f90
    subroutine test_oacc_reduction
        integer, parameter :: n = 10
        integer s1, s2
        include "openacc_lib.h"
        print *,"test_oacc_reduction: "
        s1 = 0
        s2 = 0

        !$acc parallel reduction(+:s1,s2) num_gangs (n)
        !! copy(s1)
        s1 = s1 + 1
        s2 = s2 + 1
        !$acc end parallel

        if (acc_get_device_type () .ne. acc_device_host) then
            if (s1 .ne. n) call abort
            if (s2 .ne. n) call abort
        else
            if (s1 .ne. 1) call abort
            if (s2 .ne. 1) call abort
        end if

        ! Test reductions inside subroutines

        s1 = 0
        s2 = 0
        call redsub (s1, s2, n)

        if (acc_get_device_type () .ne. acc_device_host) then
            if (s1 .ne. n) call abort
        else
            if (s2 .ne. 1) call abort
        end if

        print *,"test_oacc_reduction Completed successfully "
    end subroutine test_oacc_reduction

    subroutine redsub(s1, s2, n)
        integer :: s1, s2, n

        !$acc parallel reduction(+:s1,s2) num_gangs (10)
        !! copy(s1)
        s1 = s1 + 1
        s2 = s2 + 1
        !$acc end parallel


    end subroutine redsub


    subroutine test_gang_static_addition
        integer, parameter :: n = 100
        integer i, a(n), b(n)
        print *,"test_gang_static_addition "
        do i = 1, n
            b(i) = i
        end do

#if defined(PGI) || (defined(GNU) &&  __GNUC__ >= 6)

        if(openacc_version >= 201306) then

            !$acc parallel loop gang (static:*) num_gangs (10)
            do i = 1, n
                a(i) = b(i) + 0
            end do
            !$acc end parallel loop

            call test_gang (a, b, 0, n)

            !$acc parallel loop gang (static:1) num_gangs (10)
            do i = 1, n
                a(i) = b(i) + 1
            end do
            !$acc end parallel loop

            call test_gang (a, b, 1, n)

            !$acc parallel loop gang (static:2) num_gangs (10)
            do i = 1, n
                a(i) = b(i) + 2
            end do
            !$acc end parallel loop

            call test_gang (a, b, 2, n)

            !$acc parallel loop gang (static:5) num_gangs (10)
            do i = 1, n
                a(i) = b(i) + 5
            end do
            !$acc end parallel loop

            call test_gang (a, b, 5, n)

            !$acc parallel loop gang (static:20) num_gangs (10)
            do i = 1, n
                a(i) = b(i) + 20
            end do
            !$acc end parallel loop

            !$acc kernels loop gang (num:5, static:*)
            do i = 1, n
                a(i) = b(i) + 20
            end do
            !$acc end kernels loop

            !$acc kernels loop gang (static:20, num:30)
            do i = 1, n
                a(i) = b(i) + 20
            end do
            !$acc end kernels loop

            call test_gang (a, b, 20, n)
            print *,"test_gang_static_addition Completed"
        end if

#else
print *,"OpenACC test_gang_static_addition Not completed"
#endif

    contains

        subroutine test_gang (a, b, sarg, n)
            integer n
            integer a (n), b(n), sarg
            integer i

            do i = 1, n
                if (a(i) .ne. b(i) + sarg) STOP 1
            end do
        end subroutine test_gang
    end subroutine test_gang_static_addition


 subroutine test_nested1
    integer :: i, j, k, a(1:3, 4:6, 5:7)
    logical :: l
    l = .false.
    a(:, :, :) = 0
    print *,"test_oacc nested1"
    !$acc parallel reduction (.or.:l)
    !$acc loop worker vector collapse(4 - 1)
      do 164 i = 1, 3
        do 164 j = 4, 6
          do 164 k = 5, 7
            a(i, j, k) = i + j + k
164      end do
    !$acc loop worker vector reduction(.or.:l) collapse(2)
firstdo: do i = 1, 3
        do j = 4, 6
          do k = 5, 7
            if (a(i, j, k) .ne. (i + j + k)) l = .true.
          end do
        end do
      end do firstdo
    !$acc end parallel
    if (l) STOP 1
    print *,"test_oacc nested1: Completed"

  end subroutine test_nested1

  subroutine test_nested2
    integer :: a(3,3,3), k, kk, kkk, l, ll, lll
    a = 0
    print *,"test_oacc nested2"

    !$acc parallel num_workers(8)
    ! Use "gang(static:1)" here and below to effectively turn gang-redundant
    ! execution mode into something like gang-single.
    !$acc loop gang(static:1) collapse(1)
      do 115 k=1,3
         !$acc loop collapse(2)
  dokk: do kk=1,3
          do kkk=1,3
            a(k,kk,kkk) = 1
          enddo
        enddo dokk
115   continue
    !$acc loop gang(static:1) collapse(1)
      do k=1,3
         if (any(a(k,1:3,1:3).ne.1)) STOP 2
      enddo
    ! Use "gang(static:1)" here and below to effectively turn gang-redundant
    ! execution mode into something like gang-single.
    !$acc loop gang(static:1) collapse(1)
 dol: do 120 l=1,3
    !$acc loop collapse(2)
  doll: do ll=1,3
          do lll=1,3
            a(l,ll,lll) = 2
          enddo
        enddo doll
120   end do dol
    !$acc loop gang(static:1) collapse(1)
     do l=1,3
        if (any(a(l,1:3,1:3).ne.2)) STOP 3
     enddo
     !$acc end parallel
     print *,"test_oacc nested2: Completed"

  end subroutine test_nested2


end module simple_oacc_omp
