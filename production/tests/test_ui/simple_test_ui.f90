program simple_test_ui
use simple_user_interface ! use all in there
implicit none
class(simple_program), pointer :: ptr2prg

call make_user_interface

call get_prg_ptr('cluster2D', ptr2prg)
call ptr2prg%print_cmdline()
call ptr2prg%write2json()

! call get_prg_ptr('cluster2D_stream', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('refine3D', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('initial_3Dmodel', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('preprocess', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('extract', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('motion_correct', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()
!
! call get_prg_ptr('ctf_estimate', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('pick', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('postprocess', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('make_cavgs', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('cluster3D', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('scale', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()
!
! call get_prg_ptr('mask', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

! call get_prg_ptr('volops', ptr2prg)
! call ptr2prg%print_cmdline()
! call ptr2prg%write2json()

end program simple_test_ui
