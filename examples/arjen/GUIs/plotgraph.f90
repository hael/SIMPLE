! plotgraph.f90 --
!     Use PLplot to draw the graph of a function
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program plotgraph
    use plplot

    implicit none

    real(kind=plflt)                 :: param
    real(kind=plflt), dimension(201) :: x, y
    real(kind=plflt)                 :: xmin, ymin, xmax, ymax
    integer                          :: justify, axis
    integer                          :: i

    ! Ask for the parameter "a"
    write(*,*) 'Value for the parameter "a":'
    read(*,*) param


    ! Parse the command-line arguments
    call plparseopts(PL_PARSE_FULL)

    ! Initialise the library
    call plinit

    ! Set up the viewport with default axes
    xmin =  0.0_plflt
    xmax = 10.0_plflt
    ymin = -1.0_plflt
    ymax =  1.0_plflt

    justify = 0
    axis    = 0
    call plenv( xmin, xmax, ymin, ymax, justify, axis )
    call pllab( '(x)', '(y)', 'Function: f(x) = exp(-x) * cos(ax)' )

    ! Compute the values and draw the graph
    do i = 1,size(x)
        x(i) = (i-1) * 0.05_plflt
        y(i) = func( param, x(i) )
    enddo

    call plline( x, y )

    call plend
contains

real function func( param, x )
    real(kind=plflt) :: param, x

    func = exp(-x) * cos(param*x)

end function func

end program plotgraph
