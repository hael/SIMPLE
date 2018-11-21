program simple_test_sp_project
include 'simple_lib.f08'
use simple_oris,       only: oris
use simple_sp_project, only: sp_project
use simple_binoris,    only: binoris
implicit none

type(sp_project)     :: myproject
type(binoris)        :: bos
integer, allocatable :: strlens(:)

! prepare stack oris in project
call myproject%os_stk%new(2)
! motion_correct
call myproject%os_stk%set(1, 'movie',       'movie1.mrc')
call myproject%os_stk%set(1, 'intg',        'movie1_intg.mrc')
call myproject%os_stk%set(1, 'forctf',      'movie1_forctf.mrc')
call myproject%os_stk%set(1, 'pspec',       'movie1_pspec.mrc')
call myproject%os_stk%set(1, 'thumb',       'movie1_thumb.mrc')
call myproject%os_stk%set(1, 'intg_frames', 'movie1_intg_frames1-14.mrc')
call myproject%os_stk%set(1, 'smpd',        1.3)
call myproject%os_stk%set(2, 'movie',       'movie2.mrc')
call myproject%os_stk%set(2, 'intg',        'movie2_intg.mrc')
call myproject%os_stk%set(2, 'forctf',      'movie2_forctf.mrc')
call myproject%os_stk%set(2, 'pspec',       'movie2_pspec.mrc')
call myproject%os_stk%set(2, 'thumb',       'movie2_thumb.mrc')
call myproject%os_stk%set(2, 'intg_frames', 'movie2_intg_frames1-14.mrc')
call myproject%os_stk%set(2, 'smpd',        1.3)
! write/read
call myproject%write('myproject.simple')
call myproject%read('myproject.simple')
! call myproject%print_header
call myproject%write_segment2txt('stk', 'myproject_os_stk_1.txt')
! ctf_estimate
call myproject%os_stk%set(1, 'kv',       300.)
call myproject%os_stk%set(1, 'cs',       2.7 )
call myproject%os_stk%set(1, 'fraca',    0.1 )
call myproject%os_stk%set(1, 'dfx',      1.2 )
call myproject%os_stk%set(1, 'dfy',      1.3 )
call myproject%os_stk%set(1, 'angast',   30. )
call myproject%os_stk%set(1, 'phshift',  0.  )
call myproject%os_stk%set(1, 'ctf_estimatecc', 0.8 )
call myproject%os_stk%set(1, 'ctfres',   5.2 )
call myproject%os_stk%set(2, 'kv',       300.)
call myproject%os_stk%set(2, 'cs',       2.7 )
call myproject%os_stk%set(2, 'fraca',    0.1 )
call myproject%os_stk%set(2, 'dfx',      1.5 )
call myproject%os_stk%set(2, 'dfy',      1.7 )
call myproject%os_stk%set(2, 'angast',   60. )
call myproject%os_stk%set(2, 'phshift',  0.  )
call myproject%os_stk%set(2, 'ctf_estimatecc', 0.75)
call myproject%os_stk%set(2, 'ctfres',   3.8 )
! write/read
call myproject%write('myproject.simple')
call myproject%read('myproject.simple')
call myproject%print_info('myproject.simple')
call myproject%write_segment2txt('stk', 'myproject_os_stk_2.txt')

! generate 3 algndocs for testing merging
call myproject%os_ptcl3D%new(9)
call myproject%os_ptcl3D%set_euler(1, [1.,1.,1.])
call myproject%os_ptcl3D%set_euler(2, [1.,1.,1.])
call myproject%os_ptcl3D%set_euler(3, [1.,1.,1.])
print *, 'writing doc1'
call myproject%write('doc1.simple', [1,3])
call myproject%os_ptcl3D%new(9)
call myproject%os_ptcl3D%set_euler(4, [2.,2.,2.])
call myproject%os_ptcl3D%set_euler(5, [2.,2.,2.])
call myproject%os_ptcl3D%set_euler(6, [2.,2.,2.])
print *, 'writing doc2'
call myproject%write('doc2.simple', [4,6])
call myproject%os_ptcl3D%new(9)
call myproject%os_ptcl3D%set_euler(7, [3.,3.,3.])
call myproject%os_ptcl3D%set_euler(8, [3.,3.,3.])
call myproject%os_ptcl3D%set_euler(9, [3.,3.,3.])
print *, 'writing doc3'
call myproject%write('doc3.simple', [7,9])

end program simple_test_sp_project
