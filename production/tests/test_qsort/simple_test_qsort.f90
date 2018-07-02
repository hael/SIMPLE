
program simple_test_sort
    include 'simple_lib.f08'
    use simple_ansi_ctrls
    use simple_qsort_mt
    use simple_qsort
    use simple_qsort_c, only: sortp_1r4
    implicit none

    integer (8), parameter :: nmax =  1000000
    integer, parameter :: it_max=10
    real,    allocatable, dimension(:) :: A
    integer, allocatable, dimension(:) :: Aind
    real,    allocatable, dimension(:,:) :: B
    !  integer, allocatable, dimension(:) :: Bind
    integer(8) :: count1, count2, count3, count4, rate
    integer(8) ::  i, thr
    integer, dimension(12) :: seedx
    integer, dimension(33) :: seed8
    integer (8) :: n,nA
    real :: t1, t2, t3, t4,  t1min, t2min, t3min, t4min, t2m
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
    t1=0.;t2=0.;t3=0.;t4=0.;t2m=0.
    write(*,*) "                   N       HPSORT            QSORT     RecursiveQsort      C-QSORT         MT-QSORT"
    do n=31,6,-1
        nA = (2_8**n)
        if(nA > nmax) cycle
        allocate (A(nA));allocate ( Aind(nA))
        Aind = (/(i, i=1, nA )/)
        call make_data(A,nA)
        call system_clock(count1)
        call hpsort(A, Aind)
        call system_clock(count2)
        t1 = real(count2-count1)/(real(rate))
        ! write (*,*) "First and last in sorted list"
        ! write (*,*) A(1), A(nA)
        ! write (*,*) "Execution time in seconds:"
        ! write (*,*) real(count2-count1)/real(rate)
        deallocate(A); deallocate(Aind)

        allocate (A(nA))
        call make_data(A,nA)
        !    write (*,*) "Qsort"
        call system_clock(count1)
        call qsort(A)
        call system_clock(count2)
        t2 = real(count2-count1)/(real(rate))
        deallocate(A)

        allocate (A(nA))
        call make_data(A,nA)
        !    write (*,*) "Qsort"
        call system_clock(count1)
        call quicksort_m(A,1_8,nA)
        call system_clock(count2)
        t2m = real(count2-count1)/(real(rate))

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

        if (t1 < t2 .and. t1< t2m .and. t1<t3 .and. t1<t4 )then
            write(*,*) nA, format_str(trim(real2str(t1)),C_RED), t2, t2m, t3, t4
        else  if (t2 < t1 .and. t2< t2m .and. t2<t3 .and. t2<t4 )then
            write(*,*)nA,t1, format_str(trim(real2str(t2)),C_RED), t2m, t3, t4
        else if (t2m < t2 .and. t2m< t1 .and. t2m<t3 .and. t2m<t4 )then
            write(*,*)nA, t1, t2, format_str(trim(real2str(t2m)),C_RED), t3, t4
        else if (t3 < t2 .and. t3< t2m .and. t3<t1 .and. t1<t4 )then
            write(*,*)nA, t1, t2, t2m, format_str(trim(real2str(t3)),C_RED), t4
        else
            write(*,*)nA, t1, t2, t2m, t3, format_str(trim(real2str(t4)),C_RED)
        end if
    end do

    ! !! ITERATIONS
    ! n=12
    ! nA = 2_8**n

    ! allocate (A(nA))
    ! write (*,*) "Heap sort"
    ! t1=0.
    ! do i=1,it_max
    !     call make_data(A,nA)

    !     call system_clock(count1)
    !     call hpsort(A)
    !     call system_clock(count2)
    !     t1 = t1+real(count2-count1)
    ! end do
    ! t1 = t1/(real(it_max*rate))
    ! write (*,*) "First and last in sorted list"
    ! write (*,*) A(1), A(nA)
    ! write (*,*) "Execution time in seconds:"
    ! write (*,*) t1
    ! deallocate(A)

    ! allocate (A(nA))
    ! t2=0.
    ! write (*,*) "Qsort"
    ! do i=1,it_max
    !     call make_data(A,nA)
    !     call system_clock(count1)
    !     call qsort(A)
    !     call system_clock(count2)
    !     t2 = t2+real(count2-count1)
    ! end do
    ! t2 = t2/(real(it_max*rate))
    ! write (*,*) "First and last in sorted list"
    ! write (*,*) A(1), A(nA)
    ! write (*,*) "Execution time in seconds:"
    ! write (*,*) t2
    ! deallocate(A)


    ! allocate (A(nA))
    ! allocate(Aind(nA))
    ! t3=0.
    ! write (*,*) "Qsort C"
    ! do i=1,it_max
    !     call make_data(A,nA)
    !     Aind = (/(i, i=1, nA )/)
    !     call system_clock(count1)
    !     call sortp_1r4(INT(nA,4),Aind,A)
    !     call system_clock(count2)
    !     t3 = t3+real(count2-count1)
    ! end do
    ! t3 = t3/(real(it_max*rate))
    ! write (*,*) "First and last in sorted list"
    ! write (*,*) A(Aind(1)), A(Aind(nA))
    ! write (*,*) "Execution time in seconds:"
    ! write (*,*) t3
    ! deallocate(Aind)
    ! deallocate(A)


    ! allocate (A(nA))
    ! t4=0.
    ! write (*,*) "Multi-threaded Qsort & Merge"
    ! do i=1,it_max
    !     call make_data(A,nA)
    !     call system_clock(count1)
    !     call MTSort(A,nA,"Ascending")   ! Missing optional 4th argument means use all available threads.  To specify, add 4th argument.
    !     call system_clock(count2)
    !     t4 = t4+real(count2-count1)
    ! end do
    ! t4 = t4/(real(it_max*rate))

    ! write (*,*) "First and last in sorted list"
    ! write (*,*) A(1), A(nA)
    ! write (*,*) "Execution time in seconds:"
    ! write (*,*) t4
    ! deallocate(A)


    ! ITERATIVE  2D test
    write (*,*) "2D test for 2500x2500 real matrix"
    n=11
    nA = 2500 ! 2_8**n
    t1=0.;t2=0.;t3=0.;t4=0.
    t1min=huge(0.);t2min=t1min;t3min=t1min;t4min=t1min
    allocate (B(nA,nA));allocate (Aind(nA))
    write (*,*) "Heap sort"
    t1=0.
    call random_number(B)
    call system_clock(count3)
    do i=1,nA
         Aind = (/(i, i=1, nA )/)
        call system_clock(count1)
        call hpsort(B(:,i),Aind)
        call system_clock(count2)
        if(real(count2-count1) < t1min) t1min = (real(count2-count1))
        t1 = t1+real(count2-count1)
    end do
    call system_clock(count4)
    t1 = t1/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(1,1), B(nA,nA)
    write (*,*) "Execution time in seconds:  mean ", t1, " min ", t1min/real(rate), " total ", real(count4-count3)/real(rate)
    deallocate(B, Aind)

    allocate (B(nA,nA))
    t2=0.
    write (*,*) "Qsort"
    call random_number(B)
    call system_clock(count3)
    do i=1,nA
        call system_clock(count1)
        call qsort(B(:,i))
        call system_clock(count2)
        if(real(count2-count1) < t2min) t2min = (real(count2-count1))
        t2 = t2+real(count2-count1)
    end do
    call system_clock(count4)

    t2 = t2/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(1,1), B(nA,nA)
    write (*,*) "Execution time in seconds: mean ", t2, " min ", t2min/real(rate), " total ", real(count4-count3)/real(rate)

    deallocate(B)

    allocate (B(nA,nA))
    t2=0.
    write (*,*) "Quicksort (microbenchmark)"
    call random_number(B)
    t2=0.;t2min=huge(0.)
    call system_clock(count3)
    do i=1, nA
        call system_clock(count1)
        call quicksort_m(B(:,i), 1_8, nA)
        call system_clock(count2)
        t2 = t2+real(count2-count1)
        if(real(count2-count1) < t2min) t2min = (real(count2-count1))
    end do
    call system_clock(count4)

    t2 = t2/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(1,1), B(nA,nA)
    write (*,*) "Execution time in seconds: mean ", t2, " min ", t2min/real(rate), " total ", real(count4-count3)/real(rate)
    deallocate(B)


    allocate (B(nA,nA))
    allocate(Aind(nA))
    t3=0.
    write (*,*) "Qsort C"
    call random_number(B)
    call system_clock(count3)

    do i=1,nA
        Aind = (/(i, i=1, nA )/)
        call system_clock(count1)
        call sortp_1r4(INT(nA,4),Aind,B(:,i))
        call system_clock(count2)
        if(real(count2-count1) < t3min) t3min = (real(count2-count1))
        t3 = t3+real(count2-count1)
    end do
    call system_clock(count4)

    t3 = t3/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(nA,Aind(1)), B(nA,Aind(nA))
    write (*,*) "Execution time in seconds: mean ", t3, " min ", t3min/real(rate), " total ", real(count4-count3)/real(rate)
    deallocate(Aind)
    deallocate(B)


!     allocate (B(nA,nA))
!     t4=0.
!     write (*,*) "Multi-threaded Qsort & Merge 2D  (8 threads)"
!     call random_number(B)
!     call system_clock(count3)
!     do i=1,nA
!         call system_clock(count1)
!         call MTSort(B(:,i),nA,"Ascending",8)   ! Missing optional 4th argument means use all available threads.  To specify, add 4th argument.
!         call system_clock(count2)
!         t4 = t4+real(count2-count1)
!         if(real(count2-count1) < t4min) t4min = (real(count2-count1))
!     end do
!     call system_clock(count4)

!     t4 = t4/(real(nA*rate))

!     write (*,*) "First and last in sorted list"
!     write (*,*) B(1,1), B(nA, nA)
!     write (*,*) "Execution time in seconds: mean ", t4, " min ", t4min/real(rate), " total ", real(count4-count3)/real(rate)
!     deallocate(B)

!     allocate (B(nA,nA))
!     t4=0.
!     write (*,*) "Multi-threaded Qsort & Merge 2D  (2 threads)"
!     call random_number(B)
!     call system_clock(count3)

!     do i=1,nA
!         call system_clock(count1)
!         call MTSort(B(:,i),nA,"Ascending",2)   ! Missing optional 4th argument means use all available threads.  To specify, add 4th argument.
!         call system_clock(count2)
!         t4 = t4+real(count2-count1)
!         if(real(count2-count1) < t4min) t4min = (real(count2-count1))
!     end do
!     call system_clock(count4)

!     t4 = t4/(real(nA*rate))

!     write (*,*) "First and last in sorted list"
!     write (*,*) B(1,1), B(nA, nA)
!     write (*,*) "Execution time in seconds: mean ", t4, " min ", t4min/real(rate), " total ", real(count4-count3)/real(rate)
!     deallocate(B)
!     allocate (B(nA,nA))
!     t4=0.
!     write (*,*) "Multi-threaded Qsort & Merge 2D  (4 threads)"
!     call random_number(B)
!     call system_clock(count3)
!     do i=1,nA
!         call system_clock(count1)
!         call MTSort(B(:,i),nA,"Ascending",4)   ! Missing optional 4th argument means use all available threads.  To specify, add 4th argument.
!         call system_clock(count2)
!         t4 = t4+real(count2-count1)
!         if(real(count2-count1) < t4min) t4min = (real(count2-count1))
!     end do
!     call system_clock(count4)

!     t4 = t4/(real(nA*rate))

!     write (*,*) "First and last in sorted list"
!     write (*,*) B(1,1), B(nA, nA)
!     write (*,*) "Execution time in seconds: mean ", t4, " min ", t4min/real(rate), " total ", real(count4-count3)/real(rate)
!     deallocate(B)

!     allocate (B(nA,nA))
!     t4=0.; t4min=huge(0.)
!     write (*,*) "Multi-threaded Qsort & Merge 2D  (OMP loop + 2 threads per sort)"
!     call random_number(B)
!     call system_clock(count3)
!     !$omp parallel do private(i)
!     ! reduction(+:t4) reduction(min:t4min)
!     do i=1,nA
! !            call system_clock(count1)
!         call MTSort(B(:,i),nA,"Ascending",2)   ! Missing optional 4th argument means use all available threads.  To specify, add 4th argument.
! !            call system_clock(count2)
! !            t4 = t4+real(count2-count1)
! !            if(real(count2-count1) < t4min) t4min = (real(count2-count1))
!     end do
!     !$omp end parallel do
!     call system_clock(count4)
! !    t4 = t4/(real(nA*rate))
!     write (*,*) "First and last in sorted list"
!     write (*,*) B(1,1), B(nA, nA)
!     write (*,*) "Execution time in seconds: total ", real(count4-count3)/real(rate), " eff mean ", real(count4-count3)/real(nA*rate)
!     deallocate(B)

    allocate (B(nA,nA))
!    t2=0.
    write (*,*) "Quicksort (microbenchmark)"
    call random_number(B)
!    t2=0.;t2min=huge(0.)
    call system_clock(count3)
    !$omp parallel do private(i)
    do i=1, nA
 !       call system_clock(count1)
        call quicksort_m(B(:,i), 1_8, nA)
 !       call system_clock(count2)
 !       t2 = t2+real(count2-count1)
 !       if(real(count2-count1) < t2min) t2min = (real(count2-count1))
    end do
    !$omp end parallel do
    call system_clock(count4)
    !    t2 = t2/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(1,1), B(nA,nA)
    write (*,*) "Execution time in seconds: total ", real(count4-count3)/real(rate), " eff mean ", real(count4-count3)/real(nA*rate)
    deallocate(B)

  allocate (B(nA,nA))
!    t2=0.
    write (*,*) "Quicksort (microbenchmark)"
    call random_number(B)
!    t2=0.;t2min=huge(0.)
    call system_clock(count3)
    !$omp parallel do private(i)
    do i=1, nA
 !       call system_clock(count1)
        call quicksort_m(B(i,:), 1_8, nA)
 !       call system_clock(count2)
 !       t2 = t2+real(count2-count1)
 !       if(real(count2-count1) < t2min) t2min = (real(count2-count1))
    end do
    !$omp end parallel do
    call system_clock(count4)
    !    t2 = t2/(real(nA*rate))
    write (*,*) "First and last in sorted list"
    write (*,*) B(1,1), B(nA,nA)
    write (*,*) "Execution time in seconds: total ", real(count4-count3)/real(rate), " eff mean ", real(count4-count3)/real(nA*rate)
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
