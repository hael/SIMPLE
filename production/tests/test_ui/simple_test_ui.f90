program simple_test_ui
use simple_user_interface ! use all in there
implicit none
call new_cluster2D
call cluster2D%print_ui
! call cluster2D%write2json()
end program simple_test_ui
