program simple_test_sp_project
include 'simple_lib.f08'
use simple_sp_project, only: sp_project
use simple_binoris_io, only: binwrite_oritab
implicit none
integer, parameter :: NMICS = 87
integer, parameter :: NPTCLS = 5646
type(sp_project)   :: project1, project2, project3
character(len=:), allocatable :: template
real    :: eullims(3,2)
integer :: i
call seed_rnd
eullims(1,:) = [0.,360.]
eullims(2,:) = [0.,180.]
eullims(3,:) = [0.,360.]
! prepare dummy project
call project1%os_mic%new(NMICS, is_ptcl=.false.)
do i = 1,NMICS
    template = 'FoilHole_'//int2str_pad(i,16)
    call project1%os_mic%set(i, 'movie',    template//'_fractions.tiff')
    call project1%os_mic%set(i, 'intg',     template//'_intg.mrc')
    call project1%os_mic%set(i, 'forctf',   template//'_forctf.mrc')
    call project1%os_mic%set(i, 'thumb',    template//'_thumb.jpg')
    call project1%os_mic%set(i, 'smpd',     1.34)
    call project1%os_mic%set(i, 'kv',       300.)
    call project1%os_mic%set(i, 'cs',       2.7 )
    call project1%os_mic%set(i, 'fraca',    0.1 )
    call project1%os_mic%set(i, 'dfx',      1.+ran3())
    call project1%os_mic%set(i, 'dfy',      1.+ran3())
    call project1%os_mic%set(i, 'angast',   360.*ran3())
    call project1%os_mic%set(i, 'phshift',  0.)
enddo
! prepare dummy project
call project1%os_ptcl2D%new(NPTCLS, .true.)
do i = 1,NPTCLS
    call project1%os_ptcl2D%set_dfx(i, 1.+ran3())
    call project1%os_ptcl2D%set_dfy(i, 1.+ran3())
    call project1%os_ptcl2D%set(i, 'corr',    ran3())
    call project1%os_ptcl2D%set_class(i, nint(ran3()*100.))
    call project1%os_ptcl2D%set_state(i, nint(ran3()))
enddo
call project1%os_ptcl2D%set_all2single('w', 1.0)
project1%os_ptcl3D = project1%os_ptcl2D
call project1%os_ptcl2D%rnd_oris( trs=5.0, eullims=eullims)
call project1%os_ptcl2D%delete_entry('z')
call project1%os_ptcl3D%rnd_oris( trs=5.0, eullims=eullims)
call project1%os_ptcl3D%delete_entry('class')
call project1%os_ptcl3D%delete_entry('z')
do i = 1,NPTCLS
    call project1%os_ptcl3D%set(i, 'proj', nint(1000.*ran3()))
enddo
call project1%update_projinfo('myproject.simple')
! write/read
call project1%write('myproject.simple')
call project2%read('myproject.simple')
! compare
do i = 1,NMICS
    if( trim(project1%os_mic%ori2str(i)) /= trim(project2%os_mic%ori2str(i)) )then
        write(*,*)'1 TEST FAILED COMPARING ',trim(project1%os_mic%ori2str(i)),' AND ',trim(project2%os_mic%ori2str(i))
        stop
    endif
enddo
write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_mic'
do i = 1,NPTCLS
    if( trim(project1%os_ptcl2D%ori2str(i)) /= trim(project2%os_ptcl2D%ori2str(i)) )then
        write(*,*)'2 TEST FAILED COMPARING ',trim(project1%os_ptcl2D%ori2str(i)),' AND ',trim(project2%os_ptcl2D%ori2str(i))
        stop
    endif
enddo
write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_ptcl2D'
do i = 1,NPTCLS
    if( trim(project1%os_ptcl3D%ori2str(i)) /= trim(project2%os_ptcl3D%ori2str(i)) )then
        write(*,*)'3 TEST FAILED COMPARING ',trim(project1%os_ptcl3D%ori2str(i)),' AND ',trim(project2%os_ptcl3D%ori2str(i))
        stop
    endif
enddo
write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_ptcl3D'
call project2%kill
! merging
call project3%os_ptcl2D%new(NPTCLS, .true.)
call project3%update_projinfo('project3.simple')
call binwrite_oritab('doc_1.simple', project1, project1%os_ptcl2D, [  1,   1882],  isegment=PTCL2D_SEG)
call binwrite_oritab('doc_2.simple', project1, project1%os_ptcl2D, [1883,  3764],  isegment=PTCL2D_SEG)
call binwrite_oritab('doc_3.simple', project1, project1%os_ptcl2D, [3765, NPTCLS], isegment=PTCL2D_SEG)
call project3%merge_algndocs(NPTCLS, 3, 'ptcl2D', 'doc_', 1 )
call del_file('doc_1.simple')
call del_file('doc_2.simple')
call del_file('doc_3.simple')
do i = 1,NPTCLS
    if( trim(project1%os_ptcl2D%ori2str(i)) /= trim(project3%os_ptcl2D%ori2str(i)) )then
        write(*,*)'TEST FAILED COMPARING ',trim(project1%os_ptcl2D%ori2str(i)),' AND ',trim(project3%os_ptcl2D%ori2str(i))
        stop
    endif
enddo
write(*,*)'TEST SUCCES WRITE PARTS/MERGE/COMPARE os_ptcl2D'
call project3%kill
! some i/o functionalities
call project1%print_segment_json( 'ptcl3D', 'myproject.simple', [NPTCLS-1,NPTCLS])
call project1%write_segment2txt('mic', 'some_mics.txt', [2,5])
write(*,*)'TEST SUCCES WRITE/PRINT'
call del_file('myproject.simple')
call del_file('project3.simple')
call del_file('some_mics.txt')
end program simple_test_sp_project
