module simple_test_omp_basics
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
implicit none
private
public :: test_omp_basics, test_internal_omp, test_parallel_omp, test_shared_race_condition
contains

    subroutine test_omp_basics (n)
        integer, intent(in) ::n
        integer i
        real b(n), a(n)
        !$omp parallel do
        do i=2,n
            b(i) = (a(i) + a(i-1)) / 2.0
        enddo
        !$omp end parallel do

    end subroutine test_omp_basics

    !> example 4.1f
    subroutine test_internal_omp
        call omp_set_nested(.true.)
        call omp_set_max_active_levels(8)
        call omp_set_dynamic(.false.)
        call omp_set_num_threads(2)
        !$omp parallel
        call omp_set_num_threads(3)
        !$omp parallel
        call omp_set_num_threads(4)
        !$omp single
        ! The following should print:
        ! Inner: max_act_lev= 8 , num_thds= 3 , max_thds= 4
        ! Inner: max_act_lev= 8 , num_thds= 3 , max_thds= 4
        print *, "Inner: max_act_lev=", omp_get_max_active_levels(), &
            & ", num_thds=", omp_get_num_threads(),&
            & ", max_thds=", omp_get_max_threads()
        !$omp end single
        !$omp end parallel
        !$omp barrier
        !$omp single
        ! The following should print:
        ! Outer: max_act_lev= 8 , num_thds= 2 , max_thds= 3
        print *, "Outer: max_act_lev=", omp_get_max_active_levels(),&
            & ", num_thds=", omp_get_num_threads(), &
            & ", max_thds=", omp_get_max_threads()
        !$omp end single
        !$omp end parallel


    end subroutine test_internal_omp

    !> example 5.1f
    subroutine test_parallel_omp
        real array(10000)
        call sub(array, 10000)
    contains
        subroutine subdomain(x, istart, ipoints)
            integer istart, ipoints
            real x(*)
            integer i

            do i=1,ipoints
                x(istart+i) = 123.456
            enddo

        end subroutine subdomain

        subroutine sub(x, npoints)
            real x(*)
            integer npoints
            integer iam, nt, ipoints, istart
            !$omp parallel default(private) shared(x,npoints)
            iam = omp_get_thread_num()
            nt = omp_get_num_threads()
            ipoints = npoints/nt
            istart = iam * ipoints
            if (iam .eq. nt-1) then
                ipoints = npoints - istart
            endif
            call subdomain(x,istart,ipoints)
            !$omp end parallel
        end subroutine sub
    end  subroutine test_parallel_omp

    !> Example 10.1f
    subroutine ploop_1(a,n)
        include "omp_lib.h" ! or use omp_lib
        real a(*)
        integer i, myoffset, n
        !$omp parallel private(myoffset)
        myoffset = omp_get_thread_num()*n
        do i = 1, n
            a(myoffset+i) = float(i)
        enddo
        !$omp end parallel
    end subroutine ploop_1
    ! in exceptional cases, loop iteration variables can be made shared, as in the following example:
    !>example 10.2f
    subroutine ploop_2(a,b,n,i1,i2)
        real a(*), b(*)
        integer i1, i2, n
        !$omp parallel shared(a,b,i1,i2)
        !$omp sections
        !$omp section
        do i1 = i1, n
            if (a(i1).ne.0.0) exit
        enddo
        !$omp section
        do i2 = i2, n
            if (b(i2).ne.0.0) exit
        enddo
        !$omp end sections
        !$omp single
        if (i1.le.n) print *, 'items in a up to ', i1, 'are all zero.'
        if (i2.le.n) print *, 'items in b up to ', i2, 'are all zero.'
        !$omp end single
        !$omp end parallel
    end subroutine ploop_2

    !> example 34.1f
    subroutine test_shared_race_condition
        use omp_lib

        real a(20)
        integer mythread
        !$omp parallel shared(a) private(mythread)
        mythread = omp_get_thread_num()
        if (mythread .eq. 0) then
            call sub_shared(a(1:10)) ! compiler may introduce writes to a(6:10)
        else
            a(6:10) = 12
        endif
        !$omp end parallel
    end subroutine test_shared_race_condition
    subroutine sub_shared(x)
        real x(*)
        x(1:5) = 4
    end subroutine sub_shared

end module simple_test_omp_basics
