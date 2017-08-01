! abstract_sort.f90 --
!     Show how to define a generic sortable type
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
! Module for sortable objects
!
module sortable_types

    type, abstract :: sortable
        ! No particular data
    contains
        procedure(islower), deferred :: islower
        procedure(assignx), deferred :: assign_data
        procedure(copy), deferred    :: copy_data    ! See below
        generic :: operator(<)       => islower
        generic :: assignment(=)     => assign_data
    end type

    abstract interface
        logical function islower( item1, item2 )
            import sortable
            class(sortable), intent(in) :: item1, item2
        end function islower
    end interface

    abstract interface
        subroutine assignx( item1, item2 )
            import sortable
            class(sortable), intent(inout) :: item1
            class(sortable), intent(in)    :: item2
        end subroutine assignx
    end interface

    ! Workaround for compilers that do not support "source ="
    abstract interface
        subroutine copy( item1, item2 )
            import sortable
            class(sortable), intent(in)               :: item1
            class(sortable), intent(out), allocatable :: item2
        end subroutine copy
    end interface

contains

subroutine sort( array )
    class(sortable), dimension(:), intent(inout), target :: array

    class(sortable), allocatable :: tmp
    class(sortable), pointer     :: first_element

    integer                       :: i
    integer                       :: j

    !
    ! Workaround for compilers that do not support "source =" yet:
    !
    call array(1)%copy_data( tmp )

    !
    ! Otherwise:
    !
    ! allocate( tmp, source = array(1) )

    do i = 1,size(array)
        do j = i+1,size(array)
            if ( array(j) < array(i) ) then
                tmp      = array(i)
                array(i) = array(j)
                array(j) = tmp
            endif
        enddo
    enddo
end subroutine sort

end module sortable_types

!
! Module for printable and sortable objects
!
module printable_sortable_types
    use sortable_types

    implicit none

    type, abstract, extends(sortable) :: printable_sortable
        ! No particular data
    contains
        procedure(print_item), deferred :: print
    end type printable_sortable

    abstract interface
        subroutine print_item( item )
            import printable_sortable
            class(printable_sortable), intent(in) :: item
        end subroutine print_item
    end interface

contains

subroutine print_array( array )

    class(printable_sortable), dimension(:) :: array

    integer :: i

    do i = 1,size(array)
        call array(i)%print
    enddo

end subroutine print_array

end module printable_sortable_types

!
! Module and program to test the above
!
module addresses
    use printable_sortable_types

    implicit none

    type, extends(printable_sortable) :: address_type
         character(len=20) :: name
         character(len=20) :: city
    contains
         procedure :: copy_data   => copy_address
         procedure :: assign_data => assign_address
         procedure :: islower     => islower_address
         procedure :: print       => print_address
    end type address_type

contains

subroutine assign_address( item1, item2 )
    class(address_type), intent(inout) :: item1
    class(sortable), intent(in)        :: item2

    select type (item2)
        type is (address_type)
            item1%name = item2%name
            item1%city = item2%city
    end select
end subroutine assign_address

subroutine copy_address( item1, item2 )
    class(address_type), intent(in)               :: item1
    class(sortable),     intent(out), allocatable :: item2

    allocate( address_type :: item2 )
    select type (item2)
        type is (address_type)
            item2%name = item1%name
            item2%city = item1%city
    end select

end subroutine copy_address

logical function islower_address( item1, item2 )
    class(address_type), intent(in) :: item1
    class(sortable),     intent(in) :: item2

    select type (item2)
        type is (address_type)
            if ( item1%name /= item2%name ) then
                islower_address = item1%name < item2%name
            else
                islower_address = item1%city < item2%city
            endif
    end select

end function islower_address

subroutine print_address( item )
    class(address_type), intent(in) :: item

    write(*,'(a,4x,a)') item%name, item%city

end subroutine print_address

end module addresses


program test_addresses
    use addresses

    implicit none

    type(address_type), dimension(6) :: address

    address = (/ address_type( "John", "London" ),   &
                 address_type( "Jan", "Amsterdam" ), &
                 address_type( "Jan", "Warsaw" ),    &
                 address_type( "Jean", "Paris" ),    &
                 address_type( "Giovanni", "Rome" ), &
                 address_type( "Juan", "Madrid" )    /)

    call sort( address )
    call print( address )

end program test_addresses
