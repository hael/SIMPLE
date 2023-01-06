program simple_test_matinv
    include 'simple_lib.f08'
    implicit none

    real        :: a(3,3), b(3,3)
    integer     :: errflg

    print *, "First Test. Produces correct result"
    a = 0.
    a(1,1) = 1.
    a(2,2) = 1.
    a(3,3) = 1.
    b = 0.
    call matinv(a, b, 3, errflg)
    print *, "erflg", errflg
    print *, "a", a
    print *, "b", b

    print *, ""
    print *, "Second Test.  Produces correct result"
    a = 0.
    b = 0.
    a(1,2) = 1.
    a(2,1) = 1.
    a(3,3) = 1.
    call matinv(a, b, 3, errflg)
    print *, "erflg", errflg
    print *, "a", a
    print *, "b", b

    print *, ""
    print *, "Third Test. Produces incorrect result"
    a = 0.
    b = 0.
    ! Adding the following lines will produce a correct result
    !a(1,1) = 0.0001
    a(1,3) = 1.
    a(2,2) = 1.
    a(3,1) = 1.
    call matinv(a, b, 3, errflg)
    print *, "erflg", errflg
    print *, "a", a
    print *, "b", b

end program simple_test_matinv