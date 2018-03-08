program simple_test_ui
use simple_user_interface ! use all in there
implicit none
call make_user_interface
call cluster2D%print_ui()
call cluster2D%write2json()
call refine3D%print_ui()
call refine3D%write2json()
call initial_3Dmodel%print_ui()
call initial_3Dmodel%write2json()
call postprocess%print_ui()
call postprocess%write2json()
call extract%print_ui()
call extract%write2json()
call scale%print_ui()
call scale%write2json()
call map2ptcls%print_ui()
call map2ptcls%write2json()
call cavgassemble%print_ui()
call cavgassemble%write2json()
call make_cavgs%print_ui()
call make_cavgs%write2json()
end program simple_test_ui
