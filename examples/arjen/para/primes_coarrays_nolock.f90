! primes_coarrays_nolock.f90 --
!     Determine the first 1000 primes - using coarrays
!
!     Alternative version: image 1 collects all the data
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program primes_coarrays
    implicit none

    integer, dimension(2)                    :: range_priv
    integer, dimension(2), codimension[*]    :: range
    integer, dimension(1000)                 :: prime
    integer                                  :: total_primes
    integer, dimension(100), codimension[*]  :: primes_in_task
    integer, codimension[*]                  :: number_in_task
    logical, codimension[*]                  :: new_results
    logical, codimension[*]                  :: new_task
    logical, codimension[*]                  :: ready

    ready          = .false.
    new_results    = .true.   ! Indicates the image has results available
    new_task       = .false.  ! Indicates a new task has been issued for the image
    number_in_task = 0        ! Number of primes found in task

    range_priv(2) = 0
    total_primes  = 0

    write(*,*) 'Image', this_image(), new_results

    sync all
    write(*,*) 'Image2', this_image(), new_results

    !
    ! Collect the found primes in image 1, create new tasks
    ! for all images
    !
    do while ( .not. ready )
        if ( this_image() == 1 ) then
            call collect_results
            call create_tasks
            sync images(*)
        else
            sync images(1)
        endif

        call get_task
    enddo

    if ( this_image() == 1 ) then
        write(*,*) 'Primes:'
        write(*,'(20i5)') prime
    endif

contains

!
! Subroutine to collect the results from all
! images (run by image 1)
!
subroutine collect_results
    integer :: i
    integer :: np
    integer :: maxindex

    do i = 1,num_images()
        write(*,*) 'Examine', i, new_results[i]
        sync images(i)
        if ( new_results[i] ) then
            np       = number_in_task[i]

            maxindex = min( size(prime) - total_primes, np )
            prime(total_primes+1:total_primes+maxindex) = &
                primes_in_task(1:maxindex)[i]

            write(*,*) 'Got primes:', i, np
            total_primes = total_primes + maxindex
            write(*,*) 'Total:', total_primes
        endif
    enddo

    ready = total_primes >= size(prime)

    if ( ready ) then
        do i = 1,num_images()
            ready[i] = .true.
        enddo
    endif
end subroutine collect_results

!
! Subroutine to post new tasks (consisting of a
! range of integers in which to look for primes)
!
! Loop over the images to see which one wants a
! new task
!
subroutine create_tasks
    integer :: i

    do i = 1,num_images()
        if ( new_results[i] ) then
            new_results[i] = .false.
            range_priv(1)  = range_priv(2) + 1
            range_priv(2)  = range_priv(2) + 100

            range(:)[i]    = range_priv(:)
            new_task[i]    = .true.

            write(*,*) 'Assign', i, range_priv
        endif
    enddo

end subroutine create_tasks

!
! Subroutine to get a task and search for
! primes inside the new range
!
subroutine get_task

    integer                 :: lower
    integer                 :: upper
    logical                 :: isprime
    integer                 :: np
    integer                 :: i
    integer                 :: j
    integer                 :: residue

    if ( new_task ) then
        np = 0
        new_task = .false.
        lower    = range(1)
        upper    = range(2)

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
                primes_in_task(np) = i
            endif
        enddo
        write(*,*) this_image(), ' found:', np, ' primes'
        write(*,'(5x,10i5)' ) primes_in_task(1:np)

        number_in_task = np
        new_results    = .true.
    endif

end subroutine get_task

end program primes_coarrays
