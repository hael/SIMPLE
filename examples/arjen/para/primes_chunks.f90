! primes_chunks.f90
!     Determine the first 1000 primes using "chunks"
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program primes_chunks

    implicit none

    integer, dimension(2)    :: range
    integer                  :: number_tasks
    integer, dimension(1000) :: prime
    integer                  :: number_primes
    logical                  :: new_results
    logical                  :: new_task
    logical                  :: ready

    ready         = .false.
    new_results   = .false.
    new_task      = .false.
    number_tasks  = 0
    number_primes = 0

    range(2)      = 0

    !
    ! Determine the primes: iterate over small ranges of
    ! integers and gather the results.
    !
    do while ( .not. ready )
        range(1) = range(2) + 1
        range(2) = range(2) + 100

        call find_primes
    enddo

    write(*,'(10i5)') prime

contains

!
! Subroutine to get a task and search for
! primes inside the new range
!
subroutine find_primes

    integer, dimension(100) :: new_prime
    integer                 :: lower
    integer                 :: upper
    logical                 :: isprime
    logical                 :: got_task
    integer                 :: np
    integer                 :: i
    integer                 :: j
    integer                 :: residue
    integer                 :: maxindex

    np       = 0

    lower = range(1)
    upper = range(2)

    do i = lower,upper
        isprime = .true.
        do j = 2,int(sqrt(real(i)))
            residue = mod(i,j)
            if ( residue == 0 ) then
                isprime = .false.
                exit
            endif
        enddo
        if ( isprime ) then
            np = np + 1
            new_prime(np) = i
        endif
    enddo

     maxindex = min( size(prime) - number_primes, np )
     prime(number_primes+1:number_primes+maxindex) = new_prime(1:maxindex)
     number_primes = number_primes + maxindex

     ready = number_primes >= size(prime)

end subroutine find_primes

end program primes_chunks
