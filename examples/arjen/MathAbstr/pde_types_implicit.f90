! pde_types_implicit.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!
!     Program to demonstrate how derived types and user-defined and
!     overloaded operations can be used to solve a mathematical problem
!     in a way that closely resembles the formulation.
!
!     The problem in question is this relatively simple ordinary differential equation:
!
!        dC
!        --  = r/H (Cs - C) - k(C)
!        dt
!
!     where:
!        k(C) = 1.0 if C > 0 and 0.0 otherwise
!
!        r is the reaeration coefficient and H is the depth
!
!        Cs is the saturation concentration
!
!     Note:
!     This version use an _implicit_ method to solve the equation
!
!     Time discretisation allows us to rewrite the equation as:
!
!        C(t+Dt) - C(t)
!        --------------  = r/H (Cs - C(t+Dt)) - k(C(t))
!              Dt
!
!     We take the linear terms at the _next_ time level and the non-linear
!     term at the current time level, because otherwise we can not evaluate
!     it (it does not allow linearisation).
!
!     Rearranging this equation gives:
!
!        C(t+Dt) = C(t) + Dt * ( r/H (Cs - C(t+Dt)) - k(C(t)) )
!
!     By using a clever combination of derived types and user-defined
!     operations we can actually solve the equation _in this form_!
!
!     Note 1:
!     For simplicity we use an ordinary differential equation,
!     so there is no matrix to be inverted.
!
!     Note 2:
!     We leave out any grid definition. It would only
!     confuse matters.
!
module pde_types_zero_dimension
    implicit none

    type field_centre_0d
        real :: value
    end type field_centre_0d

    type field_centre_0d_new
        real :: value
        real :: coeff                ! Would be a matrix if dimension > 0
    end type field_centre_0d_new


    interface operator(>)
        module procedure greater_than_field_centre_0d
    end interface

    interface operator(+)
        module procedure add_field_centre_0d_new
        module procedure add_field_centre_0d_old_new
    end interface

    interface operator(-)
        module procedure subtract_field_scalar_centre_0d_new
        module procedure subtract_field_centre_scalar_0d_new
        module procedure subtract_field_centre_0d_new
    end interface

    interface operator(*)
        module procedure multiply_field_centre_scalar_0d_new
    end interface

    interface assignment(=)
        module procedure initialize_field_centre_0d      ! Assign a uniform value
        module procedure initialize_field_centre_0d_new  ! Assign a uniform value
        module procedure assign_field_centre_0d_new
    end interface

    interface replace
        module procedure assign_field_centre_0d_from_new
    end interface

contains

subroutine initialize_field_centre_0d( field, uniform_value )
    type(field_centre_0d), intent(inout) :: field
    real, intent(in)                     :: uniform_value

    field%value = uniform_value

end subroutine initialize_field_centre_0d

subroutine initialize_field_centre_0d_new( field, uniform_value )
    type(field_centre_0d_new), intent(inout) :: field
    real, intent(in)                         :: uniform_value

    field%value = uniform_value
    field%coeff = 1.0

end subroutine initialize_field_centre_0d_new

!
! Not an assignment pur sang - the value of the field_new variable must
! be set to zero
!
subroutine assign_field_centre_0d_from_new( field, field_new )
    type(field_centre_0d),     intent(inout) :: field
    type(field_centre_0d_new), intent(inout) :: field_new

    field%value     = field_new%value
    field_new%value = 0.0

end subroutine assign_field_centre_0d_from_new

!
! The central routine to implement the implicit method!
!
subroutine assign_field_centre_0d_new( field, rhs )
    type(field_centre_0d_new), intent(inout) :: field
    type(field_centre_0d_new), intent(in)    :: rhs

    field%value = rhs%value / (1.0 - rhs%coeff)

end subroutine assign_field_centre_0d_new

subroutine print_field_centre_0d( time, field )
    integer, intent(in)               :: time
    type(field_centre_0d), intent(in) :: field

    write(*,'(i10,f10.4)') time, field%value

end subroutine print_field_centre_0d

function greater_than_field_centre_0d( field, value ) result(condition)
    type(field_centre_0d), intent(in)     :: field
    real, intent(in)                      :: value
    logical                               :: condition

    condition = field%value > value
end function greater_than_field_centre_0d

function binary_minus_field_centre_scalar_0d_new( scalar, field_new )
    real, intent(in)                      :: scalar
    type(field_centre_0d_new), intent(in) :: field_new
    type(field_centre_0d_new)             :: binary_minus_field_centre_scalar_0d_new

    binary_minus_field_centre_scalar_0d_new%value = scalar - field_new%value
    binary_minus_field_centre_scalar_0d_new%coeff = -field_new%coeff

end function binary_minus_field_centre_scalar_0d_new

function multiply_field_centre_scalar_0d_new( scalar, field_new )
    real, intent(in)                      :: scalar
    type(field_centre_0d_new), intent(in) :: field_new
    type(field_centre_0d_new)             :: multiply_field_centre_scalar_0d_new

    multiply_field_centre_scalar_0d_new%value = scalar * field_new%value
    multiply_field_centre_scalar_0d_new%coeff = scalar * field_new%coeff

end function multiply_field_centre_scalar_0d_new

function add_field_centre_0d_new( field1, field2 )
    type(field_centre_0d_new), intent(in) :: field1
    type(field_centre_0d_new), intent(in) :: field2
    type(field_centre_0d_new)             :: add_field_centre_0d_new

    add_field_centre_0d_new%value = field1%value + field2%value
    add_field_centre_0d_new%coeff = field1%coeff + field2%coeff

end function add_field_centre_0d_new

function add_field_centre_0d_old_new( field1, field2 )
    type(field_centre_0d), intent(in)     :: field1
    type(field_centre_0d_new), intent(in) :: field2
    type(field_centre_0d_new)             :: add_field_centre_0d_old_new

    add_field_centre_0d_old_new%value = field1%value + field2%value
    add_field_centre_0d_old_new%coeff = field2%coeff

end function add_field_centre_0d_old_new

function subtract_field_centre_0d_new( field1, field2 )
    type(field_centre_0d_new), intent(in) :: field1
    type(field_centre_0d_new), intent(in) :: field2
    type(field_centre_0d_new)             :: subtract_field_centre_0d_new

    subtract_field_centre_0d_new%value = field1%value - field2%value
    subtract_field_centre_0d_new%coeff = field1%coeff - field2%coeff

end function subtract_field_centre_0d_new

function subtract_field_scalar_centre_0d_new( scalar, field_new )
    real, intent(in)                      :: scalar
    type(field_centre_0d_new), intent(in) :: field_new
    type(field_centre_0d_new)             :: subtract_field_scalar_centre_0d_new

    subtract_field_scalar_centre_0d_new%value = scalar - field_new%value
    subtract_field_scalar_centre_0d_new%coeff = - field_new%coeff

end function subtract_field_scalar_centre_0d_new

function subtract_field_centre_scalar_0d_new( field_new, scalar )
    type(field_centre_0d_new), intent(in) :: field_new
    real, intent(in)                      :: scalar
    type(field_centre_0d_new)             :: subtract_field_centre_scalar_0d_new

    subtract_field_centre_scalar_0d_new%value = field_new%value - scalar
    subtract_field_centre_scalar_0d_new%coeff = field_new%coeff

end function subtract_field_centre_scalar_0d_new

end module pde_types_zero_dimension

! pde_types_impl --
!     Intermediate module for providing a uniform set of type names
!
module pde_types_impl
    use pde_types_zero_dimension, field_centre      => field_centre_0d, &
                                  field_centre_new  => field_centre_0d_new
end module pde_types_impl

program test_pde_impl
    use pde_types_impl
    implicit none

    type(field_centre)      :: conc
    type(field_centre_new)  :: cnew

    integer                 :: time
    integer                 :: notimes

    real                    :: deltt
    real                    :: decay0
    real                    :: decay         ! We have a zero-dimensional model, so scalar
    real                    :: reaer
    real                    :: depth
    real                    :: csaturation

    !
    ! Set the coefficients
    !
    decay0      = 1.0    ! 1/d
    reaer       = 1.1    ! m/d  -- with 0.1 we get negative concentrations
    depth       = 2.5    ! m
    csaturation = 8.0    ! g/m3

    !
    ! Initial condition
    !
    conc = 0.0
    cnew = 0.0

    !
    ! Integrate the equation
    !
    deltt   =   1.0 !0.1
    notimes = 100

    do time = 1,notimes

        decay = merge( decay0, 0.0, conc > 0.0 )

        !
        ! Use the implicit form of the discretised equation
        !
        cnew  = conc + deltt &
                       * ( reaer/depth * (csaturation - cnew) - decay )

        !
        ! Transfer the new concentration
        !
        call replace( conc, cnew )

        if ( mod(time,20) == 0 ) then
            call print( time, conc )
        endif
    enddo
contains
subroutine print( time, conc )
    integer, intent(in)            :: time
    type(field_centre), intent(in) :: conc

    write(*,'(i10,f10.4)') time, conc%value
end subroutine print

end program test_pde_impl
