! test_kinds.f90
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Small test for an idea of generically defining the kind
!
module single
    implicit none

    integer, parameter :: wp = kind(1.0)

! -- begin includeable code

    interface print_kind
        module procedure print_kind_wp
    end interface

contains
subroutine print_kind_wp( x )
    real(kind=wp), intent(in) :: x

    write(*,*) 'Kind:', kind(x)
end subroutine print_kind_wp

! -- end includeable code

end module single

module double
    implicit none

    integer, parameter :: wp = kind(1.0d0)

! -- begin includeable code

    interface print_kind
        module procedure print_kind_wp
    end interface

contains
subroutine print_kind_wp( x )
    real(kind=wp), intent(in) :: x

    write(*,*) 'Kind:', kind(x)
end subroutine print_kind_wp

! -- end includeable code

end module double

module unite
    use single
    use double
end module unite

program test_kinds

    use unite

    call print_kind( 1.0   )
    call print_kind( 1.0d0 )

end program test_kinds
