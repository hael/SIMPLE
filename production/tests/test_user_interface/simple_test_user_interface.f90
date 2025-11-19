program simple_test_user_interface
include 'simple_lib.f08'
use simple_user_interface, only: validate_ui_json
write(*,*)'VALIDATING UI JSON FILE:'
call validate_ui_json
write(*,*)'PASSED UI JSON FILE TEST'
end program simple_test_user_interface
