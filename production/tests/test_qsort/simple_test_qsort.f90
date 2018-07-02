
program simple_test_sort
    include 'simple_lib.f08'
    use simple_qsort_mt
    use simple_qsort
    use simple_qsort_c, only: sortp_1r4
    implicit none

    integer (8), parameter :: nmax =  10000000
    integer, parameter :: it_max=10
    real,    allocatable, dimension(:) :: A
    integer, allocatable, dimension(:) :: Aind
    real,    allocatable, dimension(:,:) :: B
    integer, allocatable, dimension(:) :: Bind
    integer(8) :: count1, count2, rate
    integer ::  i
    integer, dimension(12) :: seedx
    integer, dimension(33) :: seed8
    integer (8) :: n,nA
    real :: t1, t2, t3, t4
    seedx = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)
    seed8 = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,&
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,1, 2, 3, 4, 5, 6, 7, 8, 9 /)
    call system_clock(count_rate=rate)
!     write(*,*) __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__
!     write(*,*) FC_COMPILER_VERSION
!#if(__GNUC__ >= 7) || defined(INTEL)
     call random_seed(put = seed8)
!#else
!     call random_seed(put = seedx)
!#endif
     t1=0.;t2=0.;t3=0.;t4=0.
    do n=31,8,-1
        nA = 2_8**n
        if(nA > nmax) exit
        allocate (A(nA))
        call make_data(A,nA)
        call system_clock(count1)
        call hpsort(A)
        call system_clock(count2)
        t1 = real(count2-count1)/(real(rate))
        ! write (*,*) "First and last in sorted list"
        ! write (*,*) A(1), A(nA)
        ! write (*,*) "Execution time in seconds:"
        ! write (*,*) real(count2-count1)/real(rate)

         deallocate(A)
        allocate (A(nA))
        call make_data(A,nA)
        !    write (*,*) "Qsort"
        call system_clock(count1)
        call qsort(A)
        call system_clock(count2)
        t2 = real(count2-count1)/(real(rate))

        ! write (*,*) "First and last in sorted list"
        ! write (*,*) A(1), A(nA)
        ! write (*,*) "Execution time in seconds:"
        ! write (*,*) real(count2-count1)/real(rate)
        !  t2=real(count2-count1)/real(rate)
        deallocate(A)
        allocate (A(nA))
        allocate(Aind(nA))
        call make_data(A,nA)
        Aind = (/(i, i=1, nA )/)
        !    write (*,*) "Qsort C"
        call system_clock(count1)
        call sortp_1r4(INT(nA,4),Aind,A)
        call system_clock(count2)
        t3 = real(count2-count1)/real(rate)
        deallocate(Aind)
        deallocate(A)
        allocate (A(nA))
        call make_data(A,nA)
        !    write (*,*) "Multi-threaded Qsort & Merge"
        call system_clock(count1)
        call MTSort(A,nA,"Ascending")   ! Missing optional 4th argument means use all available threads.  To specify, add 4th argument.
        call system_clock(count2)
        t4 = real(count2-count1)/(real(rate))
        deallocate(A)
        print *,nA,t1,t2,t3,t4
    end do

!     !! ITERATIONS
!     n=12
!     nA = 2_8**n

!     allocate (A(nA))
!     write (*,*) "Heap sort"
!     t1=0.
!     do i=1,it_max
!         call make_data(A,nA)

!         call system_clock(count1)
!         call hpsort(A)
!         call system_clock(count2)
!         t1 = t1+real(count2-count1)
!     end do
!     t1 = t1/(real(it_max*rate))
!     write (*,*) "First and last in sorted list"
!     write (*,*) A(1), A(nA)
!     write (*,*) "Execution time in seconds:"
!     write (*,*) t1
!     deallocate(A)

!     allocate (A(nA))
!     t2=0.
!     write (*,*) "Qsort"
!     do i=1,it_max
!         call make_data(A,nA)
!         call system_clock(count1)
!         call qsort(A)
!         call system_clock(count2)
!         t2 = t2+real(count2-count1)
!     end do
!     t2 = t2/(real(it_max*rate))
!     write (*,*) "First and last in sorted list"
!     write (*,*) A(1), A(nA)
!     write (*,*) "Execution time in seconds:"
!     write (*,*) t2
!     deallocate(A)


!     allocate (A(nA))
!     allocate(Aind(nA))
!     t3=0.
!     write (*,*) "Qsort C"
!     do i=1,it_max
!         call make_data(A,nA)
!         Aind = (/(i, i=1, nA )/)
!         call system_clock(count1)
!         call sortp_1r4(INT(nA,4),Aind,A)
!         call system_clock(count2)
!         t3 = t3+real(count2-count1)
!     end do
!     t3 = t3/(real(it_max*rate))
!     write (*,*) "First and last in sorted list"
!     write (*,*) A(Aind(1)), A(Aind(nA))
!     write (*,*) "Execution time in seconds:"
!     write (*,*) t3
!     deallocate(Aind)
!     deallocate(A)


!     allocate (A(nA))
!     t4=0.
!     write (*,*) "Multi-threaded Qsort & Merge"
!     do i=1,it_max
!         call make_data(A,nA)
!         call system_clock(count1)
!         call MTSort(A,nA,"Ascending")   ! Missing optional 4th argument means use all available threads.  To specify, add 4th argument.
!         call system_clock(count2)
!         t4 = t4+real(count2-count1)
!     end do
!     t4 = t4/(real(it_max*rate))

!     write (*,*) "First and last in sorted list"
!     write (*,*) A(1), A(nA)
!     write (*,*) "Execution time in seconds:"
!     write (*,*) t4
!     deallocate(A)

    n=9
    nA = 2_8**n

    allocate (B(nA,nA))
    write (*,*) "Heap sort"
    t1=0.
    call random_number(B)
    do i=1,nA
        call system_clock(count1)
        call hpsort(B(:,i))
        call system_clock(count2)
        t1 = t1+real(count2-count1)
    end do
    t1 = t1/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(1,1), B(nA,nA)
    write (*,*) "Execution time in seconds:"
    write (*,*) t1
    deallocate(B)

    allocate (B(nA,nA))
    t2=0.
    write (*,*) "Qsort"
    call random_number(B)
    do i=1,nA

        call system_clock(count1)
        call qsort(B(:,i))
        call system_clock(count2)
        t2 = t2+real(count2-count1)
    end do
    t2 = t2/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(1,1), B(nA,nA)
    write (*,*) "Execution time in seconds:"
    write (*,*) t2
    deallocate(B)


    allocate (B(nA,nA))
    allocate(Aind(nA))
    t3=0.
    write (*,*) "Qsort C"
    call random_number(B)
    do i=1,nA
        Aind = (/(i, i=1, nA )/)
        call system_clock(count1)
        call sortp_1r4(INT(nA,4),Aind,B(:,i))
        call system_clock(count2)
        t3 = t3+real(count2-count1)
    end do
    t3 = t3/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(nA,Aind(1)), B(nA,Aind(nA))
    write (*,*) "Execution time in seconds:"
    write (*,*) t3
    deallocate(Aind)
    deallocate(B)


    allocate (B(nA,nA))
    t4=0.
    write (*,*) "Multi-threaded Qsort & Merge 2D ",nA
    call random_number(B)
    do i=1,nA
        call system_clock(count1)
        call MTSort(B(:,i),nA,"Ascending")   ! Missing optional 4th argument means use all available threads.  To specify, add 4th argument.
        call system_clock(count2)
        t4 = t4+real(count2-count1)
    end do
    t4 = t4/(real(nA*rate))

    write (*,*) "First and last in sorted list"
    write (*,*) B(1,1), B(nA, nA)
    write (*,*) "Execution time in seconds:"
    write (*,*) t4
    deallocate(B)

contains

    subroutine make_data(A,nA)

        ! DUMMY ARGUMENTS
        integer (8), intent(in) :: nA
        real (4), dimension(nA), intent(out) :: A

        ! LOCAL VARIABLES
        integer (8) :: i
        real :: random

        do i = 1, nA
            call random_number(random)
            A(i) = random
        end do

    end subroutine make_data

end program simple_test_sort
