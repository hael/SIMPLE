program simple_test_nlopt
include 'simple_lib.f08'
use simple_nlopt_tester
implicit none
call seed_rnd
call exec_nlopt_test()
end program simple_test_nlopt
