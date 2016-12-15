! fishes.f90 --
!     Ilustration of non-type bound procedure pointers:
!     fish with three different life stages
!
!     Note:
!     initialising the procedure pointers like below is
!     a Fortran 2008 feature
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!

! flow_fields --
!     Routines and types for dealing with the flow field
!     in the sea
!
!     Note: dummy only
!
module flow_fields

    implicit none

    type flow_field_type
        ! Nothing
        integer :: dummy
    end type flow_field_type

end module flow_fields


! food_situations --
!     Routines and types for dealing with the distribution of food
!     The type of food that is of interest depends on the age
!     of the fish
!
!     Note: dummy only
!
module food_situations

    implicit none

    type food_situation_type
        ! Nothing
        integer :: dummy
    end type food_situation_type

end module food_situations


! fishes --
!     Routines and types to simulate the behaviour of fish in
!     various stages of their lives.
!
!     Note: dummy only
!
module fishes

    use flow_fields
    use food_situations

    implicit none

    type fish_type
        real                         :: x, y
        real                         :: age
        procedure(behaving), pointer :: behave                 => behave_juvenile
        procedure(update),   pointer :: update_position        => update_position_flow
        procedure(tofood),   pointer :: swim_to_food           => swim_to_food_mature
        procedure(tomating), pointer :: swim_to_mating_grounds => swim_to_mating_grounds_mating
    end type fish_type

    real, parameter :: age_reach_adulthood = 0.5   ! Grown up at age 1/2
    real, parameter :: age_reach_mating    = 5.0   ! Mating start at age 5

    abstract interface
        subroutine behaving( this, deltt, flow_field, food_situation )
            import :: fish_type, flow_field_type, food_situation_type
            class(fish_type)          :: this
            real                      :: deltt
            type(flow_field_type)     :: flow_field
            type(food_situation_type) :: food_situation
        end subroutine behaving

        subroutine update( this, deltt, flow_field )
            import :: fish_type, flow_field_type
            class(fish_type)          :: this
            real                      :: deltt
            type(flow_field_type)     :: flow_field
        end subroutine update

        subroutine tofood( this, deltt, food_situation )
            import :: fish_type, food_situation_type
            class(fish_type)          :: this
            real                      :: deltt
            type(food_situation_type) :: food_situation
        end subroutine tofood

        subroutine tomating( this, deltt )
            import :: fish_type
            class(fish_type)          :: this
            real                      :: deltt
        end subroutine tomating
    end interface

contains

subroutine behave_juvenile( this, deltt, flow_field, food_situation )
    class(fish_type)          :: this
    real                      :: deltt
    type(flow_field_type)     :: flow_field
    type(food_situation_type) :: food_situation

    !
    ! No motion of their own
    !
    call this%update_position( deltt, flow_field )

    !
    ! Update age
    !
    this%age = this%age + deltt

    if ( this%age >= age_reach_adulthood ) then
        this%behave => behave_adult
    endif
end subroutine behave_juvenile

subroutine behave_adult( this, deltt, flow_field, food_situation )
    class(fish_type)          :: this
    real                      :: deltt
    type(flow_field_type)     :: flow_field
    type(food_situation_type) :: food_situation

    !
    ! New position: where is the food?
    !
    call this%update_position( deltt, flow_field )

    call this%swim_to_food( deltt, food_situation )

    !
    ! Update age - time to mate?
    !
    this%age = this%age + deltt

    if ( this%age >= age_reach_mating ) then
        this%behave => behave_migrate
    endif
end subroutine behave_adult

subroutine behave_migrate( this, deltt, flow_field, food_situation )
    class(fish_type)          :: this
    real                      :: deltt
    type(flow_field_type)     :: flow_field
    type(food_situation_type) :: food_situation

    !
    ! New position: to the mating grounds
    !
    call this%update_position( deltt, flow_field )

    call this%swim_to_mating_grounds( deltt )

    !
    ! Update age - time to mate?
    !
    this%age = this%age + deltt
end subroutine behave_migrate

subroutine update_position_flow( this, deltt, flow_field )
    class(fish_type)          :: this
    real                      :: deltt
    type(flow_field_type)     :: flow_field

    write(*,*) 'Go with the flow'
end subroutine update_position_flow

subroutine swim_to_food_mature( this, deltt, food_situation )
    class(fish_type)          :: this
    real                      :: deltt
    type(food_situation_type) :: food_situation

    write(*,*) 'Swim to food'
end subroutine swim_to_food_mature

subroutine swim_to_mating_grounds_mating( this, deltt )
    class(fish_type)          :: this
    real                      :: deltt

    write(*,*) 'Swim to mating grounds'
end subroutine swim_to_mating_grounds_mating

end module fishes


! test_fishes --
!     Program to illustrate the simulation loop
!
program test_fishes

    use fishes

    implicit none

    integer, parameter                 :: number = 1
    type(fish_type), dimension(number) :: fish
    type(flow_field_type)              :: flow
    type(food_situation_type)          :: food

    real                               :: deltt
    integer                            :: i
    integer                            :: time

    !
    ! Initialise the information on fish, flow and food:
    ! fish in a square of 100x100 km
    !
    call random_number( fish%x )
    call random_number( fish%y )

    fish%x   =  100000.0 * fish%x
    fish%y   =  100000.0 * fish%y
    fish%age =  0.0

    deltt    =  0.1

    do time = 1,100
        do i = 1,size(fish)
            call fish(i)%behave( deltt, flow, food )
        enddo
    enddo
end program test_fishes
