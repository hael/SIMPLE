!! build with:
!! FC=mpifort cmake -DUSE_MPI=ON <simple_src>
!! make -j install
!! source add2.bashrc
!! mpiexec -np 4 simple_test_mpi

program simple_test_mpi
include 'simple_lib.f08'
#if defined( USE_MPIF08_MODULE )
#if USE_MPIF08_MODULE == 1
    use mpi_f08
#else
    use mpi
#endif
#else
    use mpi
#endif

implicit none
#include "simple_local_flags.inc"
    integer ierr, rank

    call MPI_INIT ( ierr )
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call hello
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call hello_usempif08
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call ring
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call block_distribution
    call MPI_FINALIZE ( ierr )

contains

    subroutine hello
        integer :: rank,  ierr
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank,ierr)
        print *, "Hello world, ", rank

    end subroutine hello

    subroutine hello_usempif08
        implicit none
        integer :: rank, size, len, ierr
        character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: version

        !        call MPI_INIT()
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_GET_LIBRARY_VERSION(version, len,ierr)

        write(*, '("Hello, world, I am ", i2, " of ", i2, ": ", a)') &
            rank, size, version

        !       call MPI_FINALIZE()
    end subroutine hello_usempif08


    subroutine ring
        implicit none
        integer :: rank, size, tag, next, from, i, message, ierr
        include 'mpif.h'

        ! Start up MPI

        !  call MPI_INIT()
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size,ierr)

        ! Calculate the rank of the next process in the ring.  Use the modulus
        ! operator so that the last process "wraps around" to rank zero.

        tag = 201
        next = mod((rank + 1), size)
        from = mod((rank + size - 1), size)

        ! If we are the "master" process (i.e., MPI_COMM_WORLD rank 0), put
        ! the number of times to go around the ring in the message.

        if (rank .eq. 0) then
            message = 10

            write(*, '("Process 0 sending ", i2, " to ", i2, " tag ", i3, " (", i2, " processes in ring)")') &
                &message, next, tag, size
            call MPI_SEND(message, 1, MPI_INTEGER, next, tag, MPI_COMM_WORLD,ierr)
            write(*, '("Process 0 sent to ", i2)') next
        endif

        ! Pass the message around the ring.  The exit mechanism works as
        ! follows: the message (a positive integer) is passed around the ring.
        ! Each time it passes rank 0, it is decremented.  When each processes
        ! receives a message containing a 0 value, it passes the message on to
        ! the next process and then quits.  By passing the 0 message first,
        ! every process gets the 0 message and can quit normally.

        i = 1
10      call MPI_Recv(message, i, MPI_INTEGER, from, tag, MPI_COMM_WORLD, &
            MPI_STATUS_IGNORE,ierr)

        if (rank .eq. 0) then
            message = message - 1
            write(*, '("Process 0 decremented value: ", i2)') message
        endif

        call MPI_SEND(message, 1, MPI_INTEGER, next, tag, MPI_COMM_WORLD,ierr)

        if (message .eq. 0) then
            write(*, '("Process ", i2, " exiting")') rank
            goto 20
        endif
        goto 10

        ! The last process does one extra send to process 0, which needs to be
        ! received before the program can exit

20      if (rank .eq. 0) then
            call MPI_RECV(message, 1, MPI_INTEGER, from, tag, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE,ierr)
        endif

        ! All done

        !call MPI_FINALIZE()
    end subroutine ring


    subroutine block_distribution

        implicit none
        integer, parameter                            :: dp = selected_real_kind(15,307)
        integer(8), parameter :: array_size = 1000000
        integer                                       :: ierr, num_procs, my_id, ista, iend, i, j
        integer, dimension(:), allocatable            :: ista_idx, iend_idx
        real(kind = dp)                               :: time1, time2
        real(kind = dp), dimension(:), allocatable    :: a

        !   call MPI_INIT ( ierr )
        call MPI_COMM_RANK ( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE ( MPI_COMM_WORLD, num_procs, ierr )

        if( my_id == num_procs - 1 )then
            print *, " MPI block distribution test "
            time1 = MPI_Wtime()
        end if

        !Distribute loop with block distribution
        call para_range ( 1, array_size, num_procs, my_id, ista, iend )
        allocate ( a( ista : iend ), ista_idx( num_procs ), iend_idx( num_procs ) )

        !Initialisation and saving ista and iend
        do i = ista, iend
            a(i) = sqrt( dble(i) / 3.0d+0 )
            ista_idx( my_id + 1 ) = ista
            iend_idx( my_id + 1 ) = iend
        end do
        call MPI_BARRIER ( MPI_COMM_WORLD, ierr )
        if( my_id == num_procs - 1 )then
            time2 = MPI_Wtime()
            print *, 'Initialisation real time = ', time2 - time1, 'second(s)'
            time1=MPI_Wtime()
        end if

        !Performing main calculation for all processors (including master and slaves)
        do j = 1, 1000
            do i = ista_idx( my_id + 1 ), iend_idx( my_id + 1 )
                a(i) = a(i) + sqrt( dble(i) )
            end do
        end do

        call MPI_BARRIER ( MPI_COMM_WORLD, ierr )

        if( my_id == num_procs - 1 ) then
            time2 = MPI_Wtime()
            print *, a( array_size )
            print *, 'Elapsed real time = ', time2 - time1, 'second(s)'
            print *, 'Sum ', sum(a)

            call block_distribution_openmp
            call block_distribution_serial
        end if
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        ! call MPI_FINALIZE ( ierr )
        deallocate ( a )
    end subroutine block_distribution

    !-----------------------------------------------------------------------------------------
    subroutine para_range ( n1, n2, num_procs, my_id, ista, iend )
        implicit none
        integer(8) :: n2
        integer                                       :: n1, num_procs, my_id, ista, iend, &
            iwork1, iwork2

        iwork1 = int( n2 - n1 + 1 ) / num_procs
        iwork2 = mod( int(n2 - n1 + 1), num_procs )
        ista = my_id * iwork1 + n1 + min( my_id, iwork2 )
        iend = ista + iwork1 - 1
        if( iwork2 > my_id ) then
            iend = iend + 1
        end if

    end subroutine para_range


    subroutine block_distribution_openmp
        !$ use omp_lib
        implicit none
        integer, parameter                         :: dp = selected_real_kind(15,307)
        integer(8), parameter :: nA =1000000
        integer                                    :: i, j
        ! real(kind = 8)                            :: time1, time2
        double precision :: time1,time2
        real(kind = dp), dimension(nA)        :: a
        real(kind = dp) :: atmp
        print *," OpenMP block distribution test"
        !Initialisation
        !$omp parallel do private(i)
        do i = 1, nA
            a(i) = sqrt( dble(i) / 3.0d+0 )
        end do
        !$omp end parallel do
         time1 = omp_get_wtime()
          !$omp parallel do private(j,i) &
          !$omp schedule( runtime )
             do j = 1, 1000
                  do i = 1, nA
                     a(i) = a(i) + sqrt( dble(i) )
                  end do
             end do
          !$omp end parallel do
         time2 =  omp_get_wtime()
         print *, a(1000000)
        if(abs(1000577.3502691896_dp -  a(1000000)) > 1e-14) print*, 'Error in calculation'
         print *, 'Elapsed real time = ', time2 - time1, 'second(s)'

        !  !$omp parallel do private(i)
        !  do i = 1, nA
        !      a(i) = sqrt( dble(i) / 3.0d+0 )
        !  end do
        !  !$omp end parallel do

        ! time1 = omp_get_wtime()
        ! !$omp parallel do private(i)
        ! !schedule( runtime )
        ! do i = 1, nA
        !     atmp=0.
        !     !$omp parallel do reduction(+:atmp) private(j)
        !     do j = 1, 1000
        !         atmp = atmp + sqrt( dble(i) )
        !     end do
        !     !$omp end parallel do
        !     a(i) = a(i) + atmp
        ! end do
        ! !$omp end parallel do
        ! time2 =  omp_get_wtime()
        ! print *, a(1000000)
        ! if(abs(1000577.3502691896_dp -  a(1000000)) > 1e-6) print*, 'Error in calculation'

        ! print *, 'Elapsed real time = ', time2 - time1, 'second(s)'

        print *, ' sum(a)', sum(a)

    end subroutine block_distribution_openmp

    subroutine block_distribution_serial
        implicit none
        integer, parameter                         :: dp = selected_real_kind(15,307)
        integer                                    :: i, j
        integer(kind = timer_int_kind)             :: time1, time2, crate
        real(8) :: elapsed
        real(kind = dp), dimension(1000000)        :: a
        print *, "Serial block distribution test"
        call system_clock(count=time1, count_rate=crate)
        !Initialisation
        do i = 1, 1000000
            a(i) = sqrt( dble(i) / 3.0d+0 )
        end do
        time1 = tic()
        do j = 1, 1000
            do i = 1, 1000000
                a(i) = a(i) + sqrt( dble(i) )
            end do
        end do
        call system_clock(count=time2)
        elapsed = real(time2-time1,8)/real(crate,8)
        print *, a(1000000)
        if(abs(1000577.3502691896_dp -  a(1000000)) > 1e-6 )&
            THROW_HARD('Error in Serial block distribution')

        print *, 'Elapsed real time = ', elapsed, 'second(s)'
          print *, ' sum(a)', sum(a)
    end subroutine block_distribution_serial


end program simple_test_mpi
