!> Main Program for Video 4 of VSCode for Fortran
!!
!! @author Lukas Lamm
program main
    implicit none (type, external)
    external :: sgesv
    real     :: a(2, 2)  ! Matrix A.
    real     :: b(2)     ! Vector b/x.
    real     :: pivot(2) ! Pivot indices (list of swap operations).
    integer  :: rc       ! Return code.

    a = reshape([ 2., 3., 1., 1. ], [ 2, 2 ])
    b = [ 5., 6. ]

    call sgesv(2, 1, a, 2, pivot, b, 2, rc)

    if (rc /= 0) then
        print '(a, i0)', 'Error: ', rc
        stop
    end if

    print '("Solution (x1, x2): ", f0.4, ", ", f0.4)', b
end program main