
module list
implicit none
    type node
        integer :: payload
        type (node), pointer :: next
        type(node), pointer :: left, right
    end type node
contains
    subroutine processlist(p)
        type (node), pointer :: p
        ! do work here
    end subroutine processlist
    subroutine increment_list_items (head)
        type (node), pointer :: head
        type (node), pointer :: p
        !$omp parallel private(p)
        !$omp single
        p => head
        do
            !$omp task
            ! p is firstprivate by default
            call processlist(p)
            !$omp end task
            p => p%next
            if ( .not. associated (p) ) exit
        end do
        !$omp end single
        !$omp end parallel
    end subroutine increment_list_items

end module list
    module m
        ! interface
        !     function fetch_and_add(p)
        !         integer :: fetch_and_add
        !         integer, intent(inout) :: p
        !     end function fetch_and_add
        !     function atomic_read(p)
        !         integer :: atomic_read
        !         integer, intent(in) :: p
        !     end function atomic_read
        ! end interface
        type locktype
            integer ticketnumber
            integer turn
        end type locktype
    contains
  !> Example 24.2f
    function atomic_read(p)
        integer :: atomic_read
        integer, intent(in) :: p
        ! Guarantee that the entire value of p is read atomically. No part of
        ! p can change during the read operation.
        !$omp atomic read
        atomic_read = p
        return
    end function atomic_read
    subroutine atomic_write(p, value)
        integer, intent(out) :: p
        integer, intent(in) :: value
        ! Guarantee that value is stored atomically into p. No part of p can change
        ! until after the entire write operation is completed.
        !$omp atomic write
        p = value
    end subroutine atomic_write
    subroutine work
        integer i
        do i=1,1000
        end do
    end subroutine work

    function fetch_and_add(p)
        integer:: fetch_and_add
        integer, intent(inout) :: p
        ! Atomically read the value of p and then increment it. The previous value is
        ! returned. This can be used to implement a simple lock as shown below.
        !$omp atomic capture
        fetch_and_add = p
        p = p + 1
        !$omp end atomic
    end function fetch_and_add
    ! Use fetch_and_add to implement a lock

        subroutine do_locked_work(lock)
            type(locktype), intent(inout) :: lock
            integer myturn
            integer junk
            ! obtain the lock
            myturn = fetch_and_add(lock%ticketnumber)
            do while (atomic_read(lock%turn) .ne. myturn)

                continue
            enddo
            ! Do some work. The flush is needed to ensure visibility of variables
            ! not involved in atomic directives
            !$omp flush
            call workm
            !$omp flush
            ! Release the lock
            junk = fetch_and_add(lock%turn)
        end subroutine do_locked_work
        subroutine workm
            integer i
            do i=1,100
            end do
        end subroutine workm
    end module m


module simple_test_omp_basics
    include 'simple_lib.f08'
    !$ use omp_lib
    !$ use omp_lib_kinds
    implicit none
    private
    public :: test_omp_basics, test_internal_omp,&
        test_parallel_sections_omp, test_parallel_loops_omp,&
        test_shared_race_condition, test_collapse_order,&
        test_simd_example2, test_simd_example3, test_simd_example4,&
        test_simd_example5, test_first_private, test_first_private2,&
        test_standalone_ok, test_ordered_example
contains

    subroutine test_omp_basics (n)
        integer, intent(in) ::n
        integer i
       ! integer(kind=omp_sched_kind) kind
     !   integer modifier
        real b(n), a(n)
        !$omp parallel do
        do i=2,n
            b(i) = (a(i) + a(i-1)) / 2.0
        enddo
        !$omp end parallel do
        !$omp parallel
        if(  omp_get_thread_num() == 0) then
            print *, "Number of active parallel regions              (int)    :", omp_get_active_level()
            print *, "Ancestor thread ID                             (int)    :", omp_get_ancestor_thread_num(omp_get_level())
            print *, "Dynamic teams setting                          (logical):", omp_get_dynamic()
            print *, "Number of parallel regions                     (int)    :", omp_get_level()
            print *, "Maximal number of active regions               (int)    :", omp_get_max_active_levels()
            print *, "Maximal number of threads of parallel region   (int)    :", omp_get_max_threads ()
            print *, "Nested parallel regions                        (logical):", omp_get_nested()
            print *, "Number of processors online                    (int)    :", omp_get_num_procs()
            print *, "Size of the active team                        (int)    :", omp_get_num_threads()
            !   print *, "Obtain the runtime scheduling method         :", omp_get_schedule ()
            print *, "Number of threads in a team                    (int)    :", omp_get_team_size(omp_get_level())
            print *, "Maximal number of threads                      (int)    :", omp_get_thread_limit()
            print *, "Current thread ID                              (int)    :", omp_get_thread_num()
            print *, "Whether a parallel region is active            (true)   :", omp_in_parallel()
        end if
        !$omp end parallel
        print *, "Whether a parallel region is active            (false)  :", omp_in_parallel()

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
    subroutine test_parallel_sections_omp
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
    end  subroutine test_parallel_sections_omp

    subroutine test_parallel_loops_omp
        real array(10000), array2(10000)
        call ploop_1(array, 10)
        ! call ploop_2(array, array2, 10000, 345,678)
    contains
        !> Example 10.1f
        subroutine ploop_1(a,n)
            !$ use omp_lib
            real a(*)
            integer i, myoffset, n
            !$omp parallel private(myoffset)
            myoffset = omp_get_thread_num()*n
            do i = 1, n
                a(myoffset+i) = real(i)
            enddo
            !$omp end parallel
        end subroutine ploop_1
        ! in exceptional cases, loop iteration variables can be made shared, as in the following example:
        !>example 10.2f
        subroutine ploop_2(a,b,n,i1,i2)
            !$ use omp_lib
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
    end subroutine test_parallel_loops_omp
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
    subroutine workcollapse_ordered(a,i,j)
        real a(*)
        integer i,j
        a(1:5) = 4
    end subroutine
    !>Example 12.3f
    subroutine test_collapse_order
        ! include 'omp_lib.h'
        integer j,k
        real :: a(100)
        !$omp parallel num_threads(2)
        !$omp do collapse(2) ordered private(j,k) schedule(static,3)
        do k = 1,3
            do j = 1,2
#if OPENMP_VERSION >= 201511
                !$omp ordered
                print *, omp_get_thread_num(), k, j
                !$omp end ordered
#endif
                call workcollapse_ordered(a,j,k)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine test_collapse_order

    RECURSIVE INTEGER FUNCTION fib(n) RESULT(res)
        INTEGER n, i, j
        IF ( n .LT. 2) THEN
            res = n
        ELSE
            !$OMP TASK SHARED(i)
            i = fib( n-1 )
            !$OMP END TASK
            !$OMP TASK SHARED(j)
            j = fib( n-2 )
            !$OMP END TASK
            !$OMP TASKWAIT
            res = i+j
        END IF

    END FUNCTION fib

    !>Example 16.2f
    RECURSIVE SUBROUTINE traverse ( P )
        use list
        ! TYPE Node
        !     TYPE(Node), POINTER :: left, right
        ! END TYPE Node
        TYPE(Node) :: P
        IF (associated(P%left)) THEN
            !$OMP TASK ! P is firstprivate by default
            call traverse(P%left)
            !$OMP END TASK
        ENDIF
        IF (associated(P%right)) THEN
            !$OMP TASK ! P is firstprivate by default
            call traverse(P%right)
            !$OMP END TASK
        ENDIF
        !$OMP TASKWAIT
        CALL processlist ( P%next )
    END SUBROUTINE traverse

    subroutine first_private(a)
        real(8) a
    end subroutine first_private
    !> Example_16_5f
    subroutine test_first_private
        real(8) ::  item(10000)
        integer i
        !$omp parallel
        !$omp single ! loop iteration variable i is private
        do i=1,10000
            !$omp task
            ! i is firstprivate, item is shared
            call first_private(item(i))
            !$omp end task
        end do
        !$omp end single
        !$omp end parallel
    end SUBROUTINE test_first_private
    !> Example_16_6f
    subroutine  test_first_private2
        real(8) ::  item(10000)
        integer i
        !$omp parallel
        !$omp single
        !$omp task untied
        ! loop iteration variable i is private
        do i=1,10000
            !$omp task ! i is firstprivate, item is shared
            call first_private(item(i))
            !$omp end task
        end do
        !$omp end task
        !$omp end single
        !$omp end parallel
    end subroutine test_first_private2
    subroutine check_solution(state)
        character, pointer :: state(:)
    end subroutine check_solution
    recursive subroutine bin_search(pos, n, state)
        use omp_lib
        integer :: pos, n
        character, pointer :: state(:)
        character, target, dimension(n) :: new_state1, new_state2
        integer, parameter :: LIMIT = 3
        if (pos .eq. n) then
            call check_solution(state)
            return
        endif
        !$omp task final(pos > LIMIT) mergeable
        if (.not. omp_in_final()) then
            new_state1(1:pos) = state(1:pos)
            state => new_state1
        endif
        state(pos+1) = 'z'
        call bin_search(pos+1, n, state)
        !$omp end task
        !$omp task final(pos > LIMIT) mergeable
        if (.not. omp_in_final()) then
            new_state2(1:pos) = state(1:pos)
            state => new_state2
        endif
        state(pos+1) = 'y'
        call bin_search(pos+1, n, state)
        !$omp end task
        !$omp taskwait

    end subroutine bin_search



    !> Example 27.2f
    subroutine test_standalone_ok()
        integer a
        a = 1
        if (a .ne. 0) then
            !$omp flush(a)
        endif
        if (a .ne. 0) then
            !$omp barrier
        endif
        if (a .ne. 0) then
            !$omp taskwait
        endif
        if (a .ne. 0) then
            !$omp taskyield
        endif
        goto 100
100     continue
        !$omp flush(a)
        goto 200
200     continue
        !$omp barrier
        goto 300
300     continue
        !$omp taskwait
        goto 400
400     continue
        !$omp taskyield
    end subroutine test_standalone_ok

    !> Example 28.1f
    subroutine workordered(k)
        integer k
        !$omp ordered
        write(*,*) k
        !$omp end ordered
    end subroutine workordered
    subroutine ordered_good(n)
        integer n,i
        !$omp do ordered
        do i = 1,n
            if (i <= 10) then
                !$omp ordered
                call workordered(i)
                !$omp end ordered
            endif
            if (i > 10) then
                !$omp ordered
                call workordered(i+1)
                !$omp end ordered
            endif
        enddo
    end subroutine ordered_good
    subroutine subordered(lb, ub, stride)
        integer lb, ub, stride
        integer i
        !$omp parallel do ordered schedule(dynamic)
        do i=lb,ub,stride
            call workordered(i)
        end do
        !$omp end parallel do
    end subroutine subordered
    subroutine test_ordered_example
        call subordered(1,100,5)
    end subroutine test_ordered_example

    subroutine simd_tar(a,b,c,n,ioff_ptr)
        implicit none
        double precision :: a(*),b(*),c(*)
        integer :: n, i
        integer, pointer :: ioff_ptr

        !$omp simd
        do i = 1,n
            a(i) = a(i) * b(i) * c(i+ioff_ptr)
        end do
    end subroutine simd_tar

    subroutine test_simd_example2
        implicit none
        integer, parameter :: N=32
        integer :: i
        double precision :: a(N), b(N)
        do i = 1,N
            a(i) = i-1
            b(i) = N-(i-1)
        end do
        call worksimd2(a, b, N )
        do i = 1,N
            print*, i,a(i)
        end do
        print *,"End of simd example2"
    end subroutine test_simd_example2

    function add1(a,b,fact) result(c)
        !$omp declare simd(add1) uniform(fact)
        double precision :: a,b,fact, c
        c = a + b + fact
    end function add1

    function add2(a,b,i, fact,n) result(c)
        !$omp declare simd(add2) uniform(a,b,fact) linear(i:1)
        integer,value :: i,n
        double precision :: a(n),b(n),fact, c
        c = a(i) + b(i) + fact
    end function add2

    subroutine worksimd2(a, b, n )
        implicit none
        integer :: n, i
        double precision :: a(n),b(n)
        real(8) :: tmp
!!        double precision, external :: add1, add2
        !$omp simd private(tmp)
        do i = 1,n
            tmp = add1(a(i), b(i), 1.0d0)
            a(i) = add2(a, b, i, 1.0d0,n) + tmp
            a(i) = a(i) + b(i) + 1.0d0
        end do
    end subroutine worksimd2
    subroutine test_simd_example3
        implicit none
        integer, parameter :: N=32
        integer :: i
        double precision :: a(N), b(N), cuml
        do i = 1,N
            a(i) = i-1
            b(i) = N-(i-1)
        end do
        call worksimd3(a, b, N, cuml )
        do i = 1,N
            print*, i,a(i)
        end do
        print *,"End of simd example3"
    end subroutine test_simd_example3
    subroutine worksimd3( a, b, n, sum )
        implicit none
        integer :: i, n
        double precision :: a(n), b(n), sum, tmp

        sum = 0.0d0
        !$omp simd private(tmp) reduction(+:sum)
        do i = 1,n
            tmp = a(i) + b(i)
            sum = sum + tmp
        end do
    end subroutine worksimd3
    subroutine test_simd_example4
        implicit none
        integer, parameter :: N=32, m=12
        integer :: i
        real :: a(N)
        do i = 1,N
            a(i) = i-1
        end do
        call worksimd4(a, N , m)
        do i = 1,N
            print*, i,a(i)
        end do
        !$omp flush
        print *,"End of simd example4"
    end subroutine test_simd_example4
    subroutine worksimd4( b, n, m )
        implicit none
        integer :: i,n,m
        real :: b(n)
        !$omp simd safelen(16)
        do i = m+1, n
            b(i) = b(i-m) - 1.0
        end do
    end subroutine worksimd4
    subroutine test_simd_example5
        implicit none
        integer, parameter :: N=32, m=12
        integer :: i
        double precision ::a(N),  b(N),c(N)
        do i = 1,N
            a(i) = i-1
            b(i) = N-(i-1)
        end do
        call worksimd5(a,b,c,N)
        do i = 1,N
            print*, i,c(i)
        end do
        !$omp flush
        print *,"End of simd example5"
    end subroutine test_simd_example5
    subroutine worksimd5( a, b, c, n )
        implicit none
        integer :: i,j,n
        double precision :: a(n,n), b(n,n), c(n,n), tmp
        !$omp do simd collapse(2) private(tmp)
        do j = 1,n
            do i = 1,n
                tmp = a(i,j) + b(i,j)
                c(i,j) = tmp
            end do
        end do
    end subroutine worksimd5
end module simple_test_omp_basics
