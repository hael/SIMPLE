program simple_test_ft_expanded
use simple_ft_expanded_tester
use simple_defs
implicit none
print *, 'SIMPLE ', trim(adjustl(SIMPLE_LIB_VERSION)), " (", trim(adjustl(build_descr)), ")"
call exec_ft_expanded_test
end program simple_test_ft_expanded
