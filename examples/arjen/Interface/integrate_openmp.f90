! integrate_openmp.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Integrate using OpenMP
!
module functions
    implicit none
contains
real function f(x,a)

    real :: x, a

    f = x + a
end function f
end module functions

module integration
    implicit none

    type :: integration_status
        logical :: done  = .false.
        logical :: ready = .false.
        logical :: next  = .false.
        real    :: xvalue
        real    :: fvalue
    end type

contains
subroutine start_integration( status )
    type(integration_status) :: status

    status%done  = .false.
    status%next  = .false.
    status%ready = .false.
end subroutine start_integration

subroutine set_value( status, value )
    type(integration_status) :: status
    real                     :: value

    !
    ! This is one of three places where the threads may
    ! get in each others' ways
    !
    !$omp critical
    status%fvalue = value
    status%ready  = .true.
    !$omp end critical
end subroutine set_value


subroutine get_next( status, x, next )
    type(integration_status) :: status
    real                     :: x
    logical                  :: next

    !$omp flush
    if ( status%done ) then
        next = .false.
    else
        do while ( .not. status%next )
            !$omp flush
        enddo

        !
        ! This is one of three places where the threads may
        ! get in each others' ways
        !
        !$omp critical
        x            = status%xvalue
        status%next  = .false.
        status%ready = .false.
        !$omp end critical

        next = .true.
    endif
end subroutine get_next

subroutine integrate( status, xbegin, xend, steps, result )
    type(integration_status) :: status
    real                     :: xbegin, xend, result
    integer                  :: steps

    real                     :: deltx, x
    integer                  :: i

    result = 0.0
    status%xvalue = xbegin

    deltx  = (xend - xbegin) / steps

    do i = 0,steps
        x = xbegin + deltx * i

        !
        ! This is one of three places where the threads may
        ! get in each others' ways
        !
        !$omp critical
        status%xvalue = x
        status%next   = .true.
        status%ready  = .false.
        !$omp flush
        !$omp end critical

        do while (.not. status%ready)
            !$omp flush
        enddo

        if ( i == 0 .or. i == steps ) then
            result = result + 0.5 * status%fvalue
        else
            result = result + status%fvalue
        endif
    enddo

    ! Done
    status%done = .true.
end subroutine integrate

end module integration

! Program illustrating the use
program test_integrate
    use integration
    use functions

    implicit none

    real    :: x, a, xbegin, xend, result
    integer :: steps
    logical :: next
    type(integration_status) :: status

    a      = 1.0
    xbegin =  0.0
    xend   = 10.0
    steps  = 10

    ! All can be shared!

    call start_integration( status )
    !$omp parallel sections
    !$omp section
        call integrate( status, xbegin, xend, steps, result )
    !$omp section
        do
            call get_next( status, x, next )
            if ( .not. next ) exit
            call set_value( status, f(x,a) )
        enddo
    !$omp end parallel sections

    write(*,*) 'Result: ', result

end program test_integrate
