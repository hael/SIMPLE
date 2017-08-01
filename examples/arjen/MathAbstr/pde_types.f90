! pde_types.f90 --
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
!     The problem in question is this partial differential equation:
!
!        dC      dC    d    dC
!        --  + u --  = -- D --  - k(C)
!        dt      dx    dx   dx
!
!     where k(C) = 1.0 if C > 0 and 0.0 otherwise
!
!     This one-dimensional equation can be extended to two or more dimensions
!     and it can be formulated on different coordinate systems. In vector
!     notation:
!
!        dC
!        --  + u . grad C = div (D grad C) - k(C)
!        dt
!
!     One way to solve it numerically is by rewriting it in
!     conservative form (using the fact that in a conservative
!     flow field the divergence of the flow velocity is zero) and integrating
!     it over a finite volume. By further more using the flow rate through
!     the faces of the volume we get a formulation that is independent of
!     the chosen coordinate system - the geometrical factors are incorporated
!     in the area of the faces.
!
!     The above equation can now be written as:
!
!        d VC                  dC
!        ----  = sum (-QC + DA --)  - V k(C)
!         dt                   dx
!
!     where V is the volume of the finite cell, Q the flow rate through each
!     face and A the area of the face. The summation is over all faces.
!
!     To implement a simple solution we use the following types:
!
!     - type(field_centre):
!       A variable whose values are defined in the centre of the volume
!
!     - type(field_face):
!       A variable whose values are defined on the faces of the volume and
!       therefore constitute a vector quantity
!
!     The following operations have to be implemented:
!
!     - Multiplication of type(field_face) and type(field_centre):
!       This can in fact be done in various ways, as we need to estimate
!       the value of the centrally defined quantity at each face. A common
!       method is to the use the value in the upstream (upwind) volume
!       adjacent to the face.
!
!     - Multiplication of type(field_face) and type(field_face):
!       The dispersion coefficient D is actually a 2D tensor, but usually
!       it is either isotropic or at most a diagonal tensor. Hence we
!       can model it as a vector (type(field_face)). The result is another
!       type(field_face) variable.
!
!     - Gradient of type(field_centre):
!       This is a vector quantity, that is, a type(field_face) result.
!       The actual computation can be straightforward.
!
!     - Divergence of a type(field_face) variable:
!       Here we sum the values at the faces over the whole volume to
!       get the representative value in the centre. This gives us a
!       type(field_centre) result.
!
!     - Assignment of a type(field_centre) result to a type(field_centre)
!       variable:
!       Straightforward, but note that this step could be optimised by
!       moving the allocated arrays (see also note 3).
!
!     The solution presented here is - in principle - independent of
!     the type of coordinate system or the dimension of the problem.
!     However, the problem of how to specify boundary conditions and
!     possibly initial conditions makes the geometry surface again.
!
!     To minimize the changes we require we use an intermediate module
!     that renames the various types. This way, if we change the
!     dimensionality, we only need to adjust the specification of
!     boundary and initial conditions, not the types or the integration
!     itself.
!
!     Note 1:
!     For convenience similar more operations are defined.
!
!     Note 2:
!     The method presented below is an explicit method. We can formulate
!     implicit solution methods in roughly the same way, though the
!     implementation will be radically different.
!
!     Note 3:
!     We could optimise the memory allocation and deallocation by using
!     a pool of arrays, rather than allocate and (automatically) deallocate
!     the "value" array in each operation.
!
module pde_types_one_dimension
    implicit none

    type grid_type_1d
        integer                         :: size
        real                            :: deltx
        real, dimension(:), allocatable :: volume
    end type grid_type_1d

    type field_centre_1d
        real                            :: deltx
        real, dimension(:), allocatable :: value
    end type field_centre_1d

    type field_face_1d
        real                              :: deltx
        real, dimension(:,:), allocatable :: value
    end type field_face_1d

    interface define_field
        module procedure define_field_centre_1d
        module procedure define_field_face_1d
    end interface

    interface initialize
        module procedure initialize_field_centre_1d
        module procedure initialize_field_face_1d
    end interface

    interface boundary_condition
        module procedure boundary_condition_field_centre_1d
    end interface

    interface print
        module procedure print_field_centre_1d
    end interface

    interface operator(>)
        module procedure greater_than_field_centre_1d
    end interface

    interface operator(+)
        module procedure add_field_centre_1d
        module procedure add_field_face_1d
    end interface

    interface operator(-)
        module procedure binary_minus_field_centre_1d
        module procedure unary_minus_field_centre_1d
        module procedure unary_minus_field_face_1d
    end interface

    interface operator(*)
        module procedure multiply_field_centre_1d
        module procedure multiply_field_face_1d
        module procedure multiply_field_face_centre_1d
    end interface

    interface operator(/)
        module procedure divide_scalar_by_field_centre_1d
    end interface

    interface assignment(=)
        module procedure initialize_field_centre_1d    ! Assign a uniform value
        module procedure initialize_field_face_1d
        module procedure assign_field_centre_1d
        module procedure assign_field_centre_1d_array
    end interface

    interface operator(.grad.)
        module procedure gradient_field_centre_1d
    end interface

    interface operator(.div.)
        module procedure divergence_field_face_1d
    end interface

contains
subroutine define_grid( grid, length, number_cells )
    type(grid_type_1d), intent(out) :: grid
    real, intent(in)                :: length
    integer, intent(in)             :: number_cells

    allocate( grid%volume(number_cells+2) ) ! Outer two cells do not really participate
                                            ! but it is easier if all arrays are the
                                            ! same size

    grid%size   = number_cells + 2
    grid%deltx  = length / number_cells
    grid%volume = grid%deltx

end subroutine define_grid

subroutine define_field_centre_1d( field, grid )
    type(field_centre_1d), intent(out) :: field
    type(grid_type_1d), intent(in)     :: grid

    allocate( field%value(grid%size) )
    field%deltx = grid%deltx
    field%value = 0.0

end subroutine define_field_centre_1d

subroutine define_field_face_1d( field, grid )
    type(field_face_1d), intent(out)   :: field
    type(grid_type_1d), intent(in)     :: grid

    allocate( field%value(grid%size,1) ) ! It is a one-dimensional grid after all
    field%deltx = grid%deltx
    field%value = 0.0

end subroutine define_field_face_1d

subroutine initialize_field_centre_1d( field, uniform_value )
    type(field_centre_1d), intent(inout) :: field
    real, intent(in)                     :: uniform_value

    field%value = uniform_value

end subroutine initialize_field_centre_1d

subroutine initialize_field_face_1d( field, uniform_value )
    type(field_face_1d), intent(inout) :: field
    real, intent(in)                   :: uniform_value

    field%value = uniform_value

end subroutine initialize_field_face_1d

subroutine boundary_condition_field_centre_1d( field, left_value, right_value )
    type(field_centre_1d), intent(inout) :: field
    real, intent(in)                     :: left_value
    real, intent(in)                     :: right_value

    field%value(1)                 = left_value
    field%value(size(field%value)) = right_value

end subroutine boundary_condition_field_centre_1d

subroutine print_field_centre_1d( time, field )
    integer, intent(in)               :: time
    type(field_centre_1d), intent(in) :: field

    write(*,'(i10,5f10.4,/,(10x,5f10.4))') time, field%value(2::10)

end subroutine print_field_centre_1d

function greater_than_field_centre_1d( field, value ) result(condition)
    type(field_centre_1d), intent(in)     :: field
    real, intent(in)                      :: value
    logical, dimension(size(field%value)) :: condition

    condition = field%value > value
end function greater_than_field_centre_1d

subroutine assign_field_centre_1d_array( field, array )
    type(field_centre_1d), intent(inout)  :: field
    real, dimension(:), intent(in)        :: array

    field%value = array

end subroutine assign_field_centre_1d_array

subroutine assign_field_centre_1d( field1, field2 )
    type(field_centre_1d), intent(inout)  :: field1
    type(field_centre_1d), intent(in)     :: field2

    field1%value = field2%value

end subroutine assign_field_centre_1d

function add_field_centre_1d( field1, field2 )
    type(field_centre_1d), intent(in)     :: field1
    type(field_centre_1d), intent(in)     :: field2
    type(field_centre_1d)                 :: add_field_centre_1d

    allocate( add_field_centre_1d%value(size(field1%value)) )

    add_field_centre_1d%value = field1%value + field2%value

end function add_field_centre_1d

function add_field_face_1d( field1, field2 )
    type(field_face_1d), intent(in)     :: field1
    type(field_face_1d), intent(in)     :: field2
    type(field_face_1d)                 :: add_field_face_1d

    allocate( add_field_face_1d%value(size(field1%value,1),1) )

    add_field_face_1d%value = field1%value + field2%value

end function add_field_face_1d

function binary_minus_field_centre_1d( field1, field2 )
    type(field_centre_1d), intent(in)     :: field1
    type(field_centre_1d), intent(in)     :: field2
    type(field_centre_1d)                 :: binary_minus_field_centre_1d

    allocate( binary_minus_field_centre_1d%value(size(field1%value)) )

    binary_minus_field_centre_1d%value = field1%value - field2%value

end function binary_minus_field_centre_1d

function unary_minus_field_centre_1d( field )
    type(field_centre_1d), intent(in)    :: field
    type(field_centre_1d)                :: unary_minus_field_centre_1d

    allocate( unary_minus_field_centre_1d%value(size(field%value)) )

    unary_minus_field_centre_1d%value = -field%value

end function unary_minus_field_centre_1d

function unary_minus_field_face_1d( field )
    type(field_face_1d), intent(in)     :: field
    type(field_face_1d)                 :: unary_minus_field_face_1d

    allocate( unary_minus_field_face_1d%value(size(field%value,1),1) )

    unary_minus_field_face_1d%value = -field%value

end function unary_minus_field_face_1d

function multiply_field_centre_1d( field1, field2 )
    type(field_centre_1d), intent(in)     :: field1
    type(field_centre_1d), intent(in)     :: field2
    type(field_centre_1d)                 :: multiply_field_centre_1d

    allocate( multiply_field_centre_1d%value(size(field1%value)) )

    multiply_field_centre_1d%value = field1%value * field2%value

end function multiply_field_centre_1d

function multiply_field_face_1d( field1, field2 )
    type(field_face_1d), intent(in)     :: field1
    type(field_face_1d), intent(in)     :: field2
    type(field_face_1d)                 :: multiply_field_face_1d

    allocate( multiply_field_face_1d%value(size(field1%value,1),1) )

    multiply_field_face_1d%value = field1%value * field2%value

end function multiply_field_face_1d

! multiply_field_face_centre_1d --
!     Use the upstream approach to multiply the flow field by
!     the concentration
!
function multiply_field_face_centre_1d( field1, field2 )
    type(field_face_1d),   intent(in)     :: field1
    type(field_centre_1d), intent(in)     :: field2
    type(field_face_1d)                   :: multiply_field_face_centre_1d

    integer                               :: sz

    sz = size(field1%value,1)
    allocate( multiply_field_face_centre_1d%value(size(field1%value,1),1) )

    multiply_field_face_centre_1d%value(sz,1)     = 0.0
    multiply_field_face_centre_1d%value(1:sz-1,1) = &
        field1%value(1:sz-1,1) * merge( field2%value(1:sz-1), field2%value(2:sz), field1%value(1:sz-1,1) >= 0.0 )

end function multiply_field_face_centre_1d

function divide_scalar_by_field_centre_1d( value, field )
    real, intent(in)                      :: value
    type(field_centre_1d), intent(in)     :: field
    type(field_centre_1d)                 :: divide_scalar_by_field_centre_1d

    allocate( divide_scalar_by_field_centre_1d%value(size(field%value)) )

    divide_scalar_by_field_centre_1d = value / field%value

end function divide_scalar_by_field_centre_1d

function gradient_field_centre_1d( field )
    type(field_centre_1d), intent(in)     :: field
    type(field_face_1d)                   :: gradient_field_centre_1d

    integer                               :: sz

    sz = size(field%value)
    allocate( gradient_field_centre_1d%value(sz,1) )

    gradient_field_centre_1d%value(sz,1)     = 0.0
    gradient_field_centre_1d%value(1:sz-1,1) = (field%value(2:sz) - field%value(1:sz-1)) / field%deltx

end function gradient_field_centre_1d

function divergence_field_face_1d( field )
    type(field_face_1d), intent(in)       :: field
    type(field_centre_1d)                 :: divergence_field_face_1d

    integer                               :: sz

    sz = size(field%value)
    allocate( divergence_field_face_1d%value(sz) )

    divergence_field_face_1d%value(1)    = 0.0
    divergence_field_face_1d%value(2:sz) = (field%value(2:sz,1) - field%value(1:sz-1,1)) ! No need to include the grid cell size

end function divergence_field_face_1d

end module pde_types_one_dimension

! pde_types --
!     Intermediate module for providing a uniform set of type names
!
module pde_types
    use pde_types_one_dimension, field_centre => field_centre_1d, &
                                 field_face   => field_face_1d,   &
                                 grid_type    => grid_type_1d
end module pde_types

program test_pde
    use pde_types
    implicit none

    type(grid_type)    :: grid

    type(field_centre) :: conc
    type(field_centre) :: decay
    type(field_centre) :: volume

    type(field_face)   :: flow
    type(field_face)   :: disparea

    integer            :: time
    integer            :: notimes

    real               :: deltt
    real               :: decay0

    !
    ! Set up the grid:
    ! 1 m wide, 50 grid cells
    !
    call define_grid( grid, 1.0, 50 )

    call define_field( volume,   grid )
    call define_field( conc,     grid )
    call define_field( decay,    grid )
    call define_field( flow,     grid )
    call define_field( disparea, grid )

    volume%value = grid%volume ! Should be hidden

    !
    ! Initial and boundary conditions
    ! (boundary conditions are taken to be constant)
    !
    conc = 0.0
    call boundary_condition( conc, 1.0, 0.0 )

    !
    ! Flow field and dispersion coefficient
    ! Both are constant over time and space
    ! (for simplicity only)
    !
    ! Note: both include the area of each face
    !
    flow     = 0.1
    disparea = 0.1

    !
    ! Decay rate
    !
    decay0 = 1.0

    write(88,'(a,/,52e12.4)') 'Concentration', conc%value

    !
    ! Integrate the equation
    !
    deltt   =   0.001
    notimes = 1000

    write(88,*) 'Flow: ', grid%volume(1) / (flow%value(1,1)*deltt)
    write(88,*) 'Disp: ', grid%volume(1) / (disparea%value(1,1)*deltt / grid%deltx)

    do time = 1,notimes
         !
         ! We need to renew the boundary condition - it will drift off due to
         ! indiscriminate decay process implementation otherwise
         ! (This is a consequence of treating all grid cells in the same
         ! way. An improvement would be to update only the _active_ cells)
         !
         call boundary_condition( conc, 1.0, 0.0 )

         decay = merge( decay0, 0.0, conc > 0.0 )

         conc  = conc + deltt / volume &
                        * ( .div. (-flow * conc + disparea * .grad. conc) &
                            - decay * volume )

        if ( mod(time,20) == 0 ) then
            call print( time, conc )
            write(88,'(a,/,52e12.4)') 'Concentration', conc%value
        endif
    enddo

contains
function printable_centre( v )
    type(field_centre)             :: v
    real, dimension(size(v%value)) :: printable_centre

    printable_centre = v%value
end function printable_centre
function printable_face( v )
    type(field_face)                   :: v
    real, dimension(size(v%value,1),1) :: printable_face

    printable_face = v%value
end function printable_face

end program test_pde
