program simple_test_matinv
    include 'simple_lib.f08'
    implicit none

    real        :: a(3,3), b(3,3)
    integer     :: errflg

    print *, "First Test"
    a = 0.
    a(1,1) = 1.
    a(2,2) = 1.
    a(3,3) = 1.
    b = 0.
    call matinv(a, b, 3, errflg)
    print *, errflg
    print *, a
    print *, b

    print *, "Second Test"
    a = 0.
    b = 0.
    a(1,2) = 1.
    a(2,1) = 1.
    a(3,3) = 1.
    call matinv(a, b, 3, errflg)
    print *, errflg
    print *, a
    print *, b

    print *, "Third Test"
    a = 0.
    b = 0.
    !a(1,1) = 0.0001
    a(1,3) = 1.
    a(2,2) = 1.
    a(3,1) = 1.
    call matinv(a, b, 3, errflg)
    print *, errflg
    print *, a
    print *, b

end program simple_test_matinv