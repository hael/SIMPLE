program simple_test_math
    include 'simple_lib.f08'
    use simple_math, only: butterworth
    implicit none
    
    integer :: n = 8
    real    :: s = 1., fc = 1.
    real    :: val
    
    val = butterworth(s, n, fc)
    if (abs(val - 1./sqrt(2.)) <= 10*epsilon(val)) then
        write(*, *) 'Test (s = 1, fc = 1, n = 8) passed!'
    else
        write(*, *) epsilon(val), abs(val - 1./sqrt(2.)), 'value of Butterworth poly at freq = 1 (cut-off frequency = 1) should be 1/sqrt(2)!'
    end if

    fc  = 2.
    s   = 2.
    val = butterworth(s, n, fc)
    if (abs(val - 1./sqrt(2.)) <= 10*epsilon(val)) then
        write(*, *) 'Test (s = 2, fc = 2, n = 8) passed!'
    else
        write(*, *) epsilon(val), abs(val - 1./sqrt(2.)), 'value of Butterworth poly at freq = 2 (cut-off frequency = 2) should be 1/sqrt(2)!'
    end if
end program simple_test_math
    