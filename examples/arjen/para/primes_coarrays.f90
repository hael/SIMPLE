! primes_coarrays.f90 --
!     Determine the first 1000 primes - using coarrays
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program primes_coarrays
!!    use, intrinsic :: iso_fortran_env
    implicit none

!!    type(lock_type), codimension[*]          :: primes_lock

    integer, dimension(2)                    :: range_priv
    integer, dimension(2), codimension[*]    :: range
    integer, dimension(1000), codimension[*] :: prime
    integer, codimension[*]                  :: number_primes
    logical, codimension[*]                  :: new_results
    logical, codimension[*]                  :: new_task
    logical, codimension[*]                  :: ready

    ready         = .false.
    new_results   = .false.
    new_task      = .false.
    number_primes = 0

    range_priv(2) = 0

    !
    ! Create tasks in image 1 and handle them in all images
    !
    do while ( .not. ready )
        if ( this_image() == 1 ) then
            call add_task
        endif

        sync all

        call get_task
    enddo

    if ( this_image() == 1 ) then
        write(*,'(10i5)') prime
    endif
    write(*,*) 'Image done:', this_image()

    sync all
    stop

contains

!
! Subroutine to post a new task (consisting of a
! range of integers in which to look for primes)
!
! Loop over the images to see which one wants a
! new task
!
subroutine add_task
    integer :: i

    do i = 1,num_images()
        if ( .not. new_task[i] ) then
            range_priv(1) = range_priv(2) + 1
            range_priv(2) = range_priv(2) + 100

            range(:)[i]   = range_priv(:)
            new_task[i]   = .true.

            write(*,*) 'Assign', i, range_priv
        endif
    enddo

end subroutine add_task

!
! Subroutine to get a task and search for
! primes inside the new range
!
subroutine get_task

    integer, dimension(100) :: new_prime
    integer                 :: lower
    integer                 :: upper
    logical                 :: isprime
    integer                 :: np
    integer                 :: i
    integer                 :: j
    integer                 :: residue
    integer                 :: maxindex

    np       = 0

    if ( new_task ) then
        lower = range(1)
        upper = range(2)

        !
        ! Determine what primes we have in this range
        !
        write(*,*) 'Range: ',lower, upper, this_image()
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
        write(*,*) 'Found:', np, ' primes'
        write(*,'(5x,10i5)' ) new_prime(1:np)

        !
        ! Now put the results in image 1
        !
     !!   lock( primes_lock )

        critical

        number_primes = number_primes[1]
        write(*,*) 'Storing primes:', number_primes, np, this_image()

        maxindex = min( size(prime) - number_primes, np )
        prime(number_primes+1:number_primes+maxindex)[1] = &
            new_prime(1:maxindex)

        number_primes[1] = number_primes + maxindex
        ready            = number_primes[1] >= size(prime)
        if ( ready ) then
            do i = 1,num_images()
                ready[i] = .true.
            enddo
        endif


        !
        ! We are done - get a new task, if any
        !
        new_task = .false.
        write(*,*) 'Done - new task:', this_image(), ready
     !!   unlock( primes_lock )

        end critical
    endif

    write(*,*) 'Done?', this_image(), ready

end subroutine get_task

end program primes_coarrays
