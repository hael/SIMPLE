! primes_plain.f90 --
!     Determine the first 1000 primes - no parallel programming involved
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program primes_plain

    implicit none

    integer, dimension(1000) :: prime
    integer                  :: number_primes
    integer                  :: candidate
    integer                  :: residue
    integer                  :: j
    logical                  :: isprime

    number_primes = 0
    candidate     = 2

    do while ( number_primes < size(prime) )
        isprime = .true.
        do j = 2,int(sqrt(real(candidate)))
            residue = mod(candidate,j)
            if ( residue == 0 ) then
                isprime = .false.
                exit
            endif
        enddo
        if ( isprime ) then
            number_primes = number_primes + 1
            prime(number_primes) = candidate
        endif

        candidate = candidate + 1
    enddo

    write(*,'(10i5)' ) prime

end program primes_plain
