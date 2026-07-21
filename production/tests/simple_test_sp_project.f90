program simple_test_sp_project
use simple_core_module_api
use simple_sp_project, only: sp_project
use simple_binoris_io, only: binwrite_oritab
use simple_image,       only: image
implicit none
integer, parameter :: NMICS = 87
integer, parameter :: NPTCLS = 5646
type(sp_project)   :: project1, project2, project3, zero_phase_project, mapped_phase_project
type(sp_project)   :: missing_phase_project, normalized_phase_project, read_phase_project
type(sp_project)   :: heterogeneous_phase_project
type(ctfparams)    :: zero_phase_ctf, mapped_phase_ctf, mapped_phase_result
type(image)        :: zero_phase_img
type(oris)         :: heterogeneous_phase_oris
type(string)       :: template, str1, str2, cwd_orig, probe_dir, proj_abs, projinfo_before, projinfo_after
type(string)       :: heterogeneous_stack_names(1)
real    :: eullims(3,2)
integer :: i, status
call seed_rnd
! New CTF-bearing records always carry an explicit numerical phase. Zero is
! the conventional CTF value, not an absent or provenance-dependent value.
zero_phase_ctf%smpd     = 1.
zero_phase_ctf%kv       = 300.
zero_phase_ctf%cs       = 2.7
zero_phase_ctf%fraca    = 0.1
zero_phase_ctf%dfx      = 1.5
zero_phase_ctf%dfy      = 1.6
zero_phase_ctf%angast   = 5.
zero_phase_ctf%ctfflag  = CTFFLAG_YES
zero_phase_ctf%phshift  = 0.
call zero_phase_img%new([8,8,1], zero_phase_ctf%smpd)
call zero_phase_img%write(string('zero_phase_record.mrc'), 1)
call zero_phase_project%add_single_movie(string('zero_phase_record.mrc'), zero_phase_ctf)
call zero_phase_project%add_stk(string('zero_phase_record.mrc'), zero_phase_ctf)
call assert_zero_phase_record(zero_phase_project%os_mic,    'os_mic')
call assert_zero_phase_record(zero_phase_project%os_stk,    'os_stk')
call assert_zero_phase_record(zero_phase_project%os_ptcl2D, 'os_ptcl2D')
call assert_zero_phase_record(zero_phase_project%os_ptcl3D, 'os_ptcl3D')
! Per-particle phases remain authoritative when a stack contains mixed phases.
! Prepending an existing stack also exercises local input versus global project
! particle indexing in add_stktab_2.
call zero_phase_img%write(string('heterogeneous_phase_stack.mrc'), 1)
call zero_phase_img%write(string('heterogeneous_phase_stack.mrc'), 2)
call heterogeneous_phase_oris%new(2, is_ptcl=.true.)
call heterogeneous_phase_oris%set_dfx(1, zero_phase_ctf%dfx)
call heterogeneous_phase_oris%set_dfy(1, zero_phase_ctf%dfy)
call heterogeneous_phase_oris%set(1, 'phshift', PI/4.)
call heterogeneous_phase_oris%set_dfx(2, zero_phase_ctf%dfx)
call heterogeneous_phase_oris%set_dfy(2, zero_phase_ctf%dfy)
call heterogeneous_phase_oris%set(2, 'phshift', PIO2)
heterogeneous_stack_names(1) = 'heterogeneous_phase_stack.mrc'
call heterogeneous_phase_project%add_stk(string('zero_phase_record.mrc'), zero_phase_ctf)
call heterogeneous_phase_project%add_stktab(heterogeneous_stack_names, zero_phase_ctf, &
    &heterogeneous_phase_oris)
call assert_phase_at(heterogeneous_phase_project%os_stk, 2, 'heterogeneous os_stk', 0.)
call assert_phase_at(heterogeneous_phase_project%os_ptcl2D, 2, 'heterogeneous os_ptcl2D first', PI/4.)
call assert_phase_at(heterogeneous_phase_project%os_ptcl2D, 3, 'heterogeneous os_ptcl2D second', PIO2)
call heterogeneous_phase_oris%kill
call heterogeneous_phase_project%kill
call zero_phase_img%kill
call zero_phase_project%kill
mapped_phase_ctf = zero_phase_ctf
mapped_phase_ctf%phshift = PI + PI/4.
call mapped_phase_project%add_single_movie(string('zero_phase_record.mrc'), mapped_phase_ctf)
call mapped_phase_project%add_stk(string('zero_phase_record.mrc'), mapped_phase_ctf)
call assert_phase_record(mapped_phase_project%os_mic,    'os_mic',    PI/4.)
call assert_phase_record(mapped_phase_project%os_stk,    'os_stk',    PI/4.)
call assert_phase_record(mapped_phase_project%os_ptcl2D, 'os_ptcl2D', PI/4.)
call assert_phase_record(mapped_phase_project%os_ptcl3D, 'os_ptcl3D', PI/4.)
mapped_phase_result = mapped_phase_project%get_micparams(1)
call assert_phase_value(mapped_phase_result%phshift, 'os_mic -> ctfparams', PI/4.)
mapped_phase_result = mapped_phase_project%get_ctfparams('stk', 1)
call assert_phase_value(mapped_phase_result%phshift, 'os_stk -> ctfparams', PI/4.)
mapped_phase_result = mapped_phase_project%get_ctfparams('ptcl2D', 1)
call assert_phase_value(mapped_phase_result%phshift, 'os_ptcl2D -> ctfparams', PI/4.)
mapped_phase_result = mapped_phase_project%get_ctfparams('ptcl3D', 1)
call assert_phase_value(mapped_phase_result%phshift, 'os_ptcl3D -> ctfparams', PI/4.)
call mapped_phase_project%kill
! Serialization is the final schema boundary: even partially assembled project
! records materialize the identity phase rather than persisting an absent field.
call missing_phase_project%os_mic%new(1, is_ptcl=.false.)
call missing_phase_project%os_stk%new(1, is_ptcl=.false.)
call missing_phase_project%os_ptcl2D%new(1, is_ptcl=.true.)
call missing_phase_project%os_ptcl3D%new(1, is_ptcl=.true.)
call missing_phase_project%os_mic%set(1, 'smpd', 1.)
call missing_phase_project%os_mic%write(string('missing_phase_mic.txt'))
call read_phase_project%os_mic%new(1, is_ptcl=.false.)
call read_phase_project%read_segment('mic', string('missing_phase_mic.txt'))
call assert_zero_phase_record(read_phase_project%os_mic, 'read-path os_mic')
call read_phase_project%kill
call missing_phase_project%update_projinfo(string('normalized_phase.simple'))
call missing_phase_project%write(string('normalized_phase.simple'))
call normalized_phase_project%read(string('normalized_phase.simple'))
call assert_zero_phase_record(normalized_phase_project%os_mic,    'serialized os_mic')
call assert_zero_phase_record(normalized_phase_project%os_stk,    'serialized os_stk')
call assert_zero_phase_record(normalized_phase_project%os_ptcl2D, 'serialized os_ptcl2D')
call assert_zero_phase_record(normalized_phase_project%os_ptcl3D, 'serialized os_ptcl3D')
call missing_phase_project%kill
call normalized_phase_project%kill
call del_file('normalized_phase.simple')
call del_file('zero_phase_record.mrc')
call del_file('heterogeneous_phase_stack.mrc')
call del_file('missing_phase_mic.txt')
write(*,*)'TEST SUCCESS PHASE PRESENCE, CANONICALIZATION, AND CTF PARAMETER MAPPING'
eullims(1,:) = [0.,360.]
eullims(2,:) = [0.,180.]
eullims(3,:) = [0.,360.]
! prepare dummy project
call project1%os_mic%new(NMICS, is_ptcl=.false.)
do i = 1,NMICS
    template = 'FoilHole_'//int2str_pad(i,16)
    call project1%os_mic%set(i, 'movie',    template%to_char()//'_fractions.tiff')
    call project1%os_mic%set(i, 'intg',     template%to_char()//'_intg.mrc')
    call project1%os_mic%set(i, 'forctf',   template%to_char()//'_forctf.mrc')
    call project1%os_mic%set(i, 'thumb',    template%to_char()//'_thumb.jpg')
    call project1%os_mic%set(i, 'smpd',     1.34)
    call project1%os_mic%set(i, 'kv',       300.)
    call project1%os_mic%set(i, 'cs',       2.7 )
    call project1%os_mic%set(i, 'fraca',    0.1 )
    call project1%os_mic%set(i, 'dfx',      1.+ran3())
    call project1%os_mic%set(i, 'dfy',      1.+ran3())
    call project1%os_mic%set(i, 'angast',   360.*ran3())
    call project1%os_mic%set(i, 'phshift',  PI*real(mod(i,7))/14.)
enddo
! prepare dummy project
call project1%os_ptcl2D%new(NPTCLS, .true.)
do i = 1,NPTCLS
    call project1%os_ptcl2D%set_dfx(i, 1.+ran3())
    call project1%os_ptcl2D%set_dfy(i, 1.+ran3())
    call project1%os_ptcl2D%set(i, 'phshift', PI*real(mod(i,11))/22.)
    call project1%os_ptcl2D%set(i, 'corr',    ran3())
    call project1%os_ptcl2D%set_class(i, nint(ran3()*100.))
    call project1%os_ptcl2D%set_state(i, nint(ran3()))
enddo
project1%os_ptcl3D = project1%os_ptcl2D
call project1%os_ptcl2D%rnd_oris( trs=5.0, eullims=eullims)
call project1%os_ptcl2D%delete_entry('z')
call project1%os_ptcl3D%rnd_oris( trs=5.0, eullims=eullims)
call project1%os_ptcl3D%delete_entry('class')
call project1%os_ptcl3D%delete_entry('z')
do i = 1,NPTCLS
    call project1%os_ptcl3D%set(i, 'proj', nint(1000.*ran3()))
enddo
call project1%update_projinfo(string('myproject.simple'))
! write/read
call project1%write(string('myproject.simple'))
call project2%read_segment('projinfo', string('myproject.simple'))
projinfo_before = project2%projinfo%ori2str(1)
call project2%kill
call simple_getcwd(cwd_orig)
proj_abs  = filepath(cwd_orig, string('myproject.simple'))
probe_dir = filepath(cwd_orig, string('sp_project_read_probe_'//int2str(get_process_id())))
call simple_mkdir(probe_dir)
call simple_chdir(probe_dir, status)
if( status /= 0 )then
    write(*,*)'TEST FAILED CHDIR INTO READ PROBE'
    stop
endif
call project2%read(proj_abs)
call project2%kill
call simple_chdir(cwd_orig, status)
if( status /= 0 )then
    write(*,*)'TEST FAILED CHDIR OUT OF READ PROBE'
    stop
endif
call project2%read_segment('projinfo', string('myproject.simple'))
projinfo_after = project2%projinfo%ori2str(1)
if( projinfo_before /= projinfo_after )then
    write(*,*)'TEST FAILED: sp_project%read mutated on-disk projinfo'
    write(*,*)'BEFORE: ', projinfo_before%to_char()
    write(*,*)'AFTER : ', projinfo_after%to_char()
    stop
endif
call project2%kill
call simple_rmdir(probe_dir, status)
write(*,*)'TEST SUCCES READ DOES NOT MUTATE projinfo'
call project2%read(string('myproject.simple'))
! compare
do i = 1,NMICS
    str1 = project1%os_mic%ori2str(i)
    str2 = project2%os_mic%ori2str(i)
    if( project1%os_mic%ori2str(i) /= project2%os_mic%ori2str(i) )then
        write(*,*)'1 TEST FAILED COMPARING ', str1%to_char(),' AND ', str2%to_char()
        stop
    endif
    call str1%kill
    call str2%kill
enddo
write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_mic'
do i = 1,NPTCLS
    str1 = project1%os_ptcl2D%ori2str(i)
    str2 = project2%os_ptcl2D%ori2str(i)
    if( project1%os_ptcl2D%ori2str(i) /= project2%os_ptcl2D%ori2str(i) )then
        write(*,*)'2 TEST FAILED COMPARING ', str1%to_char(),' AND ', str2%to_char()
        stop
    endif
    call str1%kill
    call str2%kill
enddo
write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_ptcl2D'
do i = 1,NPTCLS
    str1 = project1%os_ptcl3D%ori2str(i)
    str2 = project2%os_ptcl3D%ori2str(i)
    if( project1%os_ptcl3D%ori2str(i) /= project2%os_ptcl3D%ori2str(i) )then
        write(*,*)'3 TEST FAILED COMPARING ', str1%to_char(),' AND ', str2%to_char()
        stop
    endif
    call str1%kill
    call str2%kill
enddo
write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_ptcl3D'
call project2%kill
! merging
call project3%os_ptcl2D%new(NPTCLS, .true.)
call project3%update_projinfo(string('project3.simple'))
call binwrite_oritab(string('doc_1.simple'), project1, project1%os_ptcl2D, [  1,   1882],  isegment=PTCL2D_SEG)
call binwrite_oritab(string('doc_2.simple'), project1, project1%os_ptcl2D, [1883,  3764],  isegment=PTCL2D_SEG)
call binwrite_oritab(string('doc_3.simple'), project1, project1%os_ptcl2D, [3765, NPTCLS], isegment=PTCL2D_SEG)
call project3%merge_algndocs(NPTCLS, 3, 'ptcl2D', 'doc_', 1 )
call del_file('doc_1.simple')
call del_file('doc_2.simple')
call del_file('doc_3.simple')
do i = 1,NPTCLS
    str1 = project1%os_ptcl2D%ori2str(i)
    str2 = project3%os_ptcl2D%ori2str(i)
    if( project1%os_ptcl2D%ori2str(i) /= project3%os_ptcl2D%ori2str(i) )then
        write(*,*)'TEST FAILED COMPARING ',str1%to_char(),' AND ', str2%to_char()
        stop
    endif
enddo
write(*,*)'TEST SUCCES WRITE PARTS/MERGE/COMPARE os_ptcl2D'
call project3%kill
! some i/o functionalities
call project1%print_segment_json( 'ptcl3D', string('myproject.simple'), [NPTCLS-1,NPTCLS])
call project1%write_segment2txt('mic', string('some_mics.txt'), [2,5])
write(*,*)'TEST SUCCES WRITE/PRINT'
call del_file('myproject.simple')
call del_file('project3.simple')
call del_file('some_mics.txt')

contains

    subroutine assert_zero_phase_record(os, segment_name)
        class(oris),      intent(in) :: os
        character(len=*), intent(in) :: segment_name
        if( .not. os%isthere(1, 'phshift') )then
            write(*,*)'TEST FAILED: phshift missing from ', trim(segment_name)
            error stop 1
        endif
        if( abs(os%get(1, 'phshift')) > 1.e-6 )then
            write(*,*)'TEST FAILED: phshift is not zero in ', trim(segment_name)
            error stop 1
        endif
    end subroutine assert_zero_phase_record

    subroutine assert_phase_record(os, segment_name, expected)
        class(oris),      intent(in) :: os
        character(len=*), intent(in) :: segment_name
        real,             intent(in) :: expected
        if( .not. os%isthere(1, 'phshift') )then
            write(*,*)'TEST FAILED: phshift missing from ', trim(segment_name)
            error stop 1
        endif
        call assert_phase_value(os%get(1, 'phshift'), segment_name, expected)
    end subroutine assert_phase_record

    subroutine assert_phase_at(os, i, segment_name, expected)
        class(oris),      intent(in) :: os
        integer,          intent(in) :: i
        character(len=*), intent(in) :: segment_name
        real,             intent(in) :: expected
        if( .not. os%isthere(i, 'phshift') )then
            write(*,*)'TEST FAILED: phshift missing from ', trim(segment_name)
            error stop 1
        endif
        call assert_phase_value(os%get(i, 'phshift'), segment_name, expected)
    end subroutine assert_phase_at

    subroutine assert_phase_value(actual, mapping_name, expected)
        real,             intent(in) :: actual, expected
        character(len=*), intent(in) :: mapping_name
        if( abs(actual-expected) > 1.e-6 )then
            write(*,*)'TEST FAILED: phshift mapping ', trim(mapping_name), actual, expected
            error stop 1
        endif
    end subroutine assert_phase_value

end program simple_test_sp_project
