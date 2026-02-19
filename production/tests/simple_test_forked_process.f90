program simple_test_forked_process
use simple_forked_process_tester
use simple_ipc_mq_tester
implicit none

call run_all_ipc_mq_tests()
call run_all_forked_process_tests()

end program simple_test_forked_process

