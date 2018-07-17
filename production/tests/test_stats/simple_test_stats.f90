

program simple_test_stats
    !$ use omp_lib
  include 'simple_lib.f08'
  use simple_ansi_ctrls
    implicit none
    real(8) :: startt, stopt
    integer (8), parameter :: nmax =  1000000
    integer, parameter :: it_max=100
    real,    allocatable, dimension(:) :: A
    integer, allocatable, dimension(:) :: Aind
    real,    allocatable, dimension(:,:) :: B
    real:: val(5)
    !  integer, allocatable, dimension(:) :: Bind
    integer(8) :: count1, count2, count3, count4, rate
    integer(8) ::  i, thr
    integer, dimension(12) :: seedx
    integer, dimension(33) :: seed8
    integer (8) :: n,nA
    integer (4) :: m,nAsp,it, pos1,pos2
    real :: t1, t2, t3, t4,  t1min, t2min, t3min, t4min, t2m, t5
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
    write(*,*) "   N          MEDIAN    MEDIAN_NOCOPY     BAPPROX      BMEDIAN       SELEC"
    do m=18,4,-1
        nAsp = INT(2**m,4)-1
        if(nAsp > nmax) cycle
   ! do nAsp=10,1000


        allocate (A(nAsp));allocate ( Aind(nAsp))
        Aind = (/(i, i=1, nAsp )/)
        do it=1,it_max
            call make_data4(A,nAsp)

        call system_clock(count1)
        val(1) = median(A)
        call system_clock(count2)
        t1 = real(count2-count1)/(real(rate))
        enddo
        deallocate(A); deallocate(Aind)

        allocate (A(nAsp))
        do it=1,it_max
           call make_data4(A,nAsp)
        !    write (*,*) "Qsort"
        call system_clock(count1)
        val(2) = median_nocopy(A)
        call system_clock(count2)
        t2 = real(count2-count1)/(real(rate))
        enddo
        deallocate(A)

        allocate (A(nAsp))
        do it=1,it_max
            call make_data4(A,nAsp)
            !    write (*,*) "Qsort"
            call system_clock(count1)
            val(3)= bapprox(nAsp,A)
            call system_clock(count2)
            t3 = real(count2-count1)/(real(rate))
        enddo
        ! write (*,*) "First and last in sorted list"
        ! write (*,*) A(1), A(nAsp)
        ! write (*,*) "Execution time in seconds:"
        ! write (*,*) real(count2-count1)/real(rate)
        !  t2=real(count2-count1)/real(rate)
        deallocate(A)



        allocate (A(nAsp))
        allocate(Aind(nAsp))
        do it=1,it_max
           call make_data4(A,nAsp)
           call system_clock(count1)
           if( is_even(nAsp) )then
               val(4)= bmedian(nAsp-1,A(1:nAsp-1))
           else
               val(4)= bmedian(nAsp,A)
           end if
        call system_clock(count2)
        t4 = real(count2-count1)/real(rate)
        enddo
        deallocate(Aind)
        deallocate(A)

        allocate (A(nAsp))
        allocate(Aind(nAsp))
        do it=1,it_max
            call make_data4(A,nAsp)
            call system_clock(count1)
            if( is_even(nAsp) )then
                pos1 = nAsp/2
                pos2 = pos1+1
                val(5)= quickselect(pos1,nAsp,A)
                val(5)= val(5)+ quickselect(pos2,nAsp,A)
                val(5)= val(5)/2.
            else
                pos1 = nint(real(nAsp)/2.)
                val(5)= quickselect(pos1,nAsp,A)
            endif


            call system_clock(count2)
            t5 = real(count2-count1)/real(rate)
        enddo
        deallocate(Aind)
        deallocate(A)



        call print_table(INT(nAsp,8), t1, t2, t3, t4,t5)

    end do
    t1=0.;t2=0.;t3=0.;t4=0.;t2m=0.

    write(*,*) ""
    write(*,*) "EVEN "

    write(*,*) "   N          MEDIAN    MEDIAN_NOCOPY     BAPPROX      BMEDIAN       SELEC"
    do m=18,4,-1
        nAsp = INT(2**m,4)
        if(nAsp > nmax) cycle
    !do nAsp=10,1000
        allocate (A(nAsp));allocate ( Aind(nAsp))
        Aind = (/(i, i=1, nAsp )/)
        do it=1,it_max
        call make_data4(A,nAsp)
        call system_clock(count1)
        val(1) = median(A)
        call system_clock(count2)
        t1 = real(count2-count1)/(real(rate))
        enddo
        deallocate(A); deallocate(Aind)

        allocate (A(nAsp))
        do it=1,it_max
           call make_data4(A,nAsp)
        !    write (*,*) "Qsort"
        call system_clock(count1)
        val(2) = median_nocopy(A)
        call system_clock(count2)
        t2 = real(count2-count1)/(real(rate))
        enddo
        deallocate(A)

        allocate (A(nAsp))
        do it=1,it_max
            call make_data4(A,nAsp)
            !    write (*,*) "Qsort"
            call system_clock(count1)
            val(3)= bapprox(nAsp,A)
            call system_clock(count2)
            t3 = real(count2-count1)/(real(rate))
        enddo
        ! write (*,*) "First and last in sorted list"
        ! write (*,*) A(1), A(nAsp)
        ! write (*,*) "Execution time in seconds:"
        ! write (*,*) real(count2-count1)/real(rate)
        !  t2=real(count2-count1)/real(rate)
        deallocate(A)



        allocate (A(nAsp))
        allocate(Aind(nAsp))
        do it=1,it_max
           call make_data4(A,nAsp)
           call system_clock(count1)
           val(4)= bapproxmp(nAsp,A)

!           if( is_even(nAsp) )then
!               val(4)= bmedian(nAsp-1,A(1:nAsp-1))
!           else
!               val(4)= bmedian(nAsp,A)
!           end if
        call system_clock(count2)
        t4 = real(count2-count1)/real(rate)
        enddo
        deallocate(Aind)
        deallocate(A)

        allocate (A(nAsp))
        allocate(Aind(nAsp))
        do it=1,it_max
            call make_data4(A,nAsp)
            call system_clock(count1)
            if( is_even(nAsp) )then
                pos1 = nAsp/2
                pos2 = pos1+1
                val(5)= quickselect(pos1,nAsp,A)
                val(5)= val(5)+ quickselect(pos2,nAsp,A)
                val(5)= val(5)/2.
            else
                pos1 = nint(real(nAsp)/2.)
                val(5)= quickselect(pos1,nAsp,A)
            endif


            call system_clock(count2)
            t5 = real(count2-count1)/real(rate)
        enddo
        deallocate(Aind)
        deallocate(A)



        call print_table(INT(nAsp,8), t1, t2, t3, t4,t5)

    end do


contains

    subroutine make_data8(A,nA)

        ! DUMMY ARGUMENTS
        integer (8), intent(in) :: nA
        real (4), dimension(nA), intent(out) :: A

        ! LOCAL VARIABLES
        integer (8) :: i
        real :: random

        do i = 1, nA
            call random_number(random)
            A(i) = 25.0*random
        end do

    end subroutine make_data8

    subroutine make_data4(A,nAsp)

        ! DUMMY ARGUMENTS
        integer (4), intent(in) :: nAsp
        real (4), dimension(nAsp), intent(out) :: A

        ! LOCAL VARIABLES
        integer (4) :: i
        real :: random

        do i = 1, nAsp
            call random_number(random)
            A(i) = 25.0*random
        end do

    end subroutine make_data4

    subroutine print_table(nA, t1, t2, t3, t4, t5)
      integer(8) :: nA
      real :: t1, t2, t3, t4
      real, optional :: t5
      character(len=:), allocatable :: redc
      character(len=:), allocatable :: rede
      redc=achar(27)//'[31m'
      rede=achar(27)//'[0m'



      write(*,'(I6,1x)', advance='no') nA
      if (t1 < t2  .and. t1<t3 .and. t1<t4 )then
!         write(*,'(a,1x)', advance='no')  adjustr(format_str(trim(real2str(t1)),C_RED)//" ")
          write(*,'(a,ES15.8,a,1x)', advance='no') redc,t1,rede
      else
         write(*,'(ES15.8,2x)', advance='no') t1
      endif
      write(*,'(3x)', advance='no')
      if (t2 < t1  .and. t2<t3 .and. t2<t4 )then
         !write(*,'(a,1x)', advance='no')  adjustr(format_str(trim(real2str(t2)),C_RED))
          write(*,'(a,ES15.8,a,1x)', advance='no') redc,t2,rede
      else
         write(*,'(ES15.8,2x)', advance='no') t2
      endif
      write(*,'(3x)', advance='no')
      if (t3 < t2  .and. t3<t1 .and. t3<t4 )then
!          write(*,'(a,1x)', advance='no') adjustr(format_str(trim(real2str(t3)),C_RED))
          write(*,'(a,ES15.8,a,1x)', advance='no') redc,t3,rede
      else
         write(*,'(ES15.8,2x)', advance='no') t3
      endif
      write(*,'(3x)', advance='no')
      if (t4 < t1 .and. t4<t2 .and. t4<t3 )then
!          write(*,'(a,1x)', advance='no') adjustr(format_str(trim(real2str(t4)),C_RED))
          write(*,'(a,ES15.8,a,1x)', advance='no') redc,t4,rede
      else
         write(*,'(ES15.8,2x)', advance='no') t4
      endif
      write(*,'(3x)', advance='no')

      if(present(t5))then
         if (t5 < t1 .and. t5< t2  .and. t5<t3 .and. t5<t4 )then
             !write(*,'(a,1x)', advance='no')  adjustr(format_str(trim(real2str(t5)),C_RED))
             write(*,'(a,ES15.8,a,1x)', advance='no') redc,t5,rede

         else
            write(*,'(ES20.8,1x)', advance='no') t5
         endif
      endif
      write(*,*) ""

    end subroutine print_table



    ! http://www.stat.cmu.edu/~ryantibs/median/
    ! Some sample Fortran code for the binapprox algorithm.
    real function bapproxmp(n,x)
        !$ use omp_lib
      integer n
      real x(n)
      integer i,bottom,counts(0:1000),bin,k,j,count,medbin
      real mu,sigma,scalefac,leftend,rghtend,xsq

      !     Compute the mean and standard deviation
      if(n<1000)then
          mu=sum(x)/n
          sigma=sqrt(dot_product(x-mu,x-mu)/n)
      else
      !$omp parallel do reduction(+:mu) private(i)
      do i=1,n
          mu = mu+x(i)
      end do
      !$omp end parallel do
      mu=mu/n
      xsq = dotprod(x - mu, x-mu,n)
      sigma = sqrt(xsq/n)
      endif
      !     Bin x across the interval [mu-sigma, mu+sigma]
      bottom = 0
      counts = 0
      scalefac = 1000/(2*sigma)
      leftend = mu-sigma
      rghtend = mu+sigma
      !$omp parallel do reduction(+:bottom) private(i,bin)
      do  i = 1,n
         if (x(i).lt.leftend) then
            bottom = bottom+1
         else if (x(i).lt.rghtend) then
            bin = int((x(i)-leftend) * scalefac)
            counts(bin) = counts(bin)+1
         endif
     end do
     !$omp end parallel do
      !     If n is odd
      if (mod(n,2).eq.1) then
         !        Find the bin that contains the median
         k = (n+1)/2
         count = bottom
         do  i = 0,1000
            count = count+counts(i)
            if (count.ge.k) then
               bapproxmp = (i+0.5)/scalefac + leftend
               return
            endif
         end do
         !     If n is even
      else
         !        Find the bins that contain the medians
         k = n/2
         count = bottom
         do  i = 0,1000
            count = count+counts(i)
            if (count.ge.k) then
               j = i
               do while  (count.eq.k)
                  j = j+1
                  count = count+counts(j)
               end do
               bapproxmp = (i+j+1)/(2*scalefac) + leftend
               return
            endif
         end do
      endif
    end function bapproxmp


    function dotprod (B, C, n) result(sum)
        real :: B(N), C(N), sum
        integer :: N, i
        sum = 0.0e0
        !$omp target map(to: B, C) map(tofrom: sum)
        !$omp teams num_teams(8) thread_limit(16) reduction(+:sum)
        !$omp distribute parallel do reduction(+:sum) &
        !$omp& dist_schedule(static, 1024) schedule(static, 64)
        do i = 1, N
            sum = sum + B(i) * C(i)
        end do
        !$omp end teams
        !$omp end target
    end function dotprod


  end program simple_test_stats
