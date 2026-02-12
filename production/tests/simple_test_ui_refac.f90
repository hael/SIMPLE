program simple_test_ui_refac
use simple_core_module_api
use simple_user_interface, only: validate_ui_json, test_ui_refactoring_func
implicit none
write(*,*)'VALIDATING UI JSON FILE:'
call validate_ui_json
write(*,*)'PASSED UI JSON FILE TEST'
call test_ui_refactoring_func
end program simple_test_ui_refac
