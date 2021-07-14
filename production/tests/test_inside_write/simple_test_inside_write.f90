program simple_test_inside_write
use simple_sp_project, only: sp_project
implicit none
type(sp_project) :: spproj
integer          :: i
call spproj%os_mic%new(10, is_ptcl=.false.)
do i=1,10
    call spproj%os_mic%set(i, 'intg', 'intg')
end do
call spproj%os_ptcl2D%new(100, is_ptcl=.true.)
call spproj%os_ptcl2D%rnd_oris
call spproj%write('original_proj.simple')
call spproj%write('updated_proj.simple')

! insert a field into file
call spproj%os_stk%new(10, is_ptcl=.false.)
do i=1,10
    call spproj%os_stk%set(i, 'stk', 'stk')
end do
call spproj%write_segment_inside('stk', 'updated_proj.simple')
end program simple_test_inside_write
