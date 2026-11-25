!@descr: unit tests of the select2D save (cls2D_selection%save): no folder before a save, <n>_selection/, prg=selection's segments
! A six-class fixture (class 4 with pop 0, hidden) of 40 particles: one without a class and one at state 0 in a
! class that stays selected (select2D has no append, so it takes its class's state). A session saves a selection,
! then prg=selection append=no runs in process on the same input with the same 0/1 states; the two projects must
! agree but for the folder and the jobproc row of the job. The input sits in a numbered run folder (the selection
! goes beside it) or in a plain folder (the selection goes inside it)
! (doc/implementation_notes/planned/select2D_viewer.md, decisions 12, 13 and 14).
module simple_select2D_tester
use simple_test_utils
use simple_error,           only: simple_exception
use simple_string,          only: string
use simple_chash,           only: chash
use simple_oris,            only: oris
use simple_image,           only: image
use simple_sp_project,      only: sp_project
use simple_cls2D_selection, only: cls2D_selection
use simple_fileio,          only: filepath, stemname, simple_copy_file, files_identical, fopen, fclose
use simple_syslib,          only: file_exists, dir_exists, simple_list_dirs, simple_chdir, simple_mkdir
implicit none
private
public :: run_all_select2D_tests
#include "simple_local_flags.inc"

integer,          parameter :: NCLS        = 6
integer,          parameter :: NPTCLS      = 40
integer,          parameter :: HIDDEN      = 4    ! pop 0: not shown, always written 0
integer,          parameter :: DESELECTED  = 2    ! the class the session deselects first
integer,          parameter :: DESELECTED2 = 5    ! and the one the second save deselects as well
integer,          parameter :: ZERO_PTCL   = 3    ! state 0 in class 3, which stays selected
integer,          parameter :: BOX         = 32   ! parameters refuses a box under 26 (prg=selection reads it from stk)
real,             parameter :: SMPD        = 2.0
integer,          parameter :: CLASSES(5)  = [1, 2, 3, 5, 6]
character(len=*), parameter :: INPUT_NAME  = 'select2D_input.simple'
character(len=*), parameter :: SNAPSHOT    = 'select2D_input_snapshot.simple'
character(len=*), parameter :: STATES_FILE = 'select2D_states.txt'
character(len=*), parameter :: STK_NAME    = 'select2D_cavgs.mrc'

contains

    subroutine run_all_select2D_tests()
        write(*,'(A)') '**** running all select2D save tests ****'
        call test_save(run_layout=.true.)
        call test_save(run_layout=.false.)
    end subroutine run_all_select2D_tests

    !> one session, the input in the run folder 7_solve2D/ (run_layout) or in the plain folder ws/: no folder
    !! before the first save, the first save makes <n>_selection/ in the workspace root, the input stays
    !! byte-identical, the saved project equals prg=selection append=no's, a second save rewrites it
    subroutine test_save( run_layout )
        use simple_ui,                       only: make_ui
        use simple_cmdline,                  only: cmdline
        use simple_commanders_project_core, only: commander_selection
        logical, intent(in) :: run_layout
        class(sp_project),     allocatable :: disp, saved_proj, ref_proj
        class(cmdline),        allocatable :: cline
        type(cls2D_selection)              :: sel
        type(commander_selection)          :: xselection
        type(chash)                        :: job_descr
        type(string)                       :: cwd_saved, root, projdir, workroot, input, saved, ref, ref_name
        character(len=:),      allocatable :: tag, first_dir, ref_dir, next_dir
        integer                            :: nfail0, states(NCLS), njobs, ndirs0
        write(*,'(A)') 'test_save, run folder layout='//merge('yes', 'no ', run_layout)
        call make_ui
        tag = 'select2D_save_'//trim(merge('run_folder  ', 'plain_folder', run_layout))
        nfail0 = tests_failed
        call enter_fixture(tag, cwd_saved, root)
        allocate(disp, saved_proj, ref_proj, cline)
        ! the workspace root (workroot) holds the selection folders; numbers continue its numbered folders
        if( run_layout )then
            projdir   = filepath(root, '7_solve2D')
            workroot  = root
            first_dir = '8_selection'
            ref_dir   = '9_selection'
            next_dir  = '10_selection'
            ndirs0    = 1
        else
            projdir   = filepath(root, 'ws')
            workroot  = projdir
            first_dir = '1_selection'
            ref_dir   = '2_selection'
            next_dir  = '3_selection'
            ndirs0    = 0
        endif
        call simple_mkdir(projdir, verbose=.false.)
        input = filepath(projdir, INPUT_NAME)
        call make_fixture(projdir, input)
        call simple_copy_file(input, string(SNAPSHOT))
        ! the session as the commander opens it: the classes and their stack only
        call disp%read_segment('cls2D', input)
        call disp%read_segment('out',   input)
        call job_descr%new()
        call job_descr%push('prg',      'select2D')
        call job_descr%push('projfile', input%to_char())
        call sel%new(disp, input, job_descr)
        call assert_int(NCLS - 1, sel%get_nvisible(), 'a class with pop 0 is not shown')
        call assert_int(ndirs0, count_dirs(workroot), 'opening a session makes no directory')
        call assert_int(0, get_saved_len(sel), 'no saved project before the first save')
        ! first save: class DESELECTED off
        call sel%click(DESELECTED, .false.)
        call sel%save()
        saved = sel%get_saved_projfile()
        call assert_string_eq(workroot%to_char()//'/'//first_dir//'/'//INPUT_NAME, saved, &
            &'the first save writes <workspace root>/'//first_dir//'/<input name>')
        call assert_true(file_exists(saved), 'the saved project exists')
        call assert_int(ndirs0 + 1, count_dirs(workroot), 'the first save makes one directory')
        call assert_true(files_identical(input, string(SNAPSHOT)), 'the input project is not written')
        call assert_false(sel%is_dirty(), 'a save clears the unsaved mark')
        states = expected_cls_states([DESELECTED])
        call saved_proj%read(saved)
        call assert_states(saved_proj, states, 'first save')
        call assert_string_eq(saved%to_char(), saved_proj%projinfo%get_str(1, 'projfile'), 'projinfo names the saved project')
        call assert_string_eq(stemname_of(saved), saved_proj%projinfo%get_str(1, 'cwd'), 'projinfo cwd is the selection folder')
        njobs = saved_proj%jobproc%get_noris()
        call assert_int(1, njobs, 'the save records one job')
        if( njobs == 1 )then
            call assert_string_eq('select2D', saved_proj%jobproc%get_str(1, 'prg'), 'the job is select2D')
            call assert_string_eq(saved%to_char(), saved_proj%jobproc%get_str(1, 'projfile'), 'the job names the saved project')
        endif
        ! prg=selection append=no (select2D has no append) on the same input with the same states, in process
        ! from the workspace root
        call write_states(filepath(root, STATES_FILE), states)
        call cline%set('prg',      'selection')
        call cline%set('projfile', input)
        call cline%set('oritype',  'cls2D')
        call cline%set('infile',   filepath(root, STATES_FILE))
        call cline%set('append',   'no')
        call cline%set('mkdir',    'yes')
        call cline%set('nthr',     1)
        call simple_chdir(workroot)
        call xselection%execute(cline)
        call simple_chdir(root)
        ref = filepath(filepath(workroot, ref_dir), INPUT_NAME)
        call assert_true(file_exists(ref), 'prg=selection writes its project to '//ref_dir//'/')
        if( file_exists(ref) )then
            call ref_proj%read(ref)
            call assert_same_records(saved_proj%os_mic,    ref_proj%os_mic,    'mic')
            call assert_same_records(saved_proj%os_stk,    ref_proj%os_stk,    'stk')
            call assert_same_records(saved_proj%os_ptcl2D, ref_proj%os_ptcl2D, 'ptcl2D')
            call assert_same_records(saved_proj%os_ptcl3D, ref_proj%os_ptcl3D, 'ptcl3D')
            call assert_same_records(saved_proj%os_cls2D,  ref_proj%os_cls2D,  'cls2D')
            call assert_same_records(saved_proj%os_cls3D,  ref_proj%os_cls3D,  'cls3D')
            call assert_same_records(saved_proj%os_out,    ref_proj%os_out,    'out')
            call assert_same_records(saved_proj%os_optics, ref_proj%os_optics, 'optics')
            call assert_same_records(saved_proj%compenv,   ref_proj%compenv,   'compenv')
            ref_name = ref_proj%projinfo%get_str(1, 'projname')
            call assert_string_eq(ref_name%to_char(), saved_proj%projinfo%get_str(1, 'projname'), &
                &'projinfo: the same project name as prg=selection')
            call assert_int(ref_proj%jobproc%get_noris() + 1, njobs, 'jobproc: prg=selection rows plus the select2D job')
            call ref_proj%kill
        endif
        call saved_proj%kill
        ! second save: DESELECTED2 off as well, into the same folder
        call sel%click(DESELECTED2, .false.)
        call sel%save()
        call assert_string_eq(saved%to_char(), sel%get_saved_projfile(), 'a later save writes the same project')
        call assert_int(ndirs0 + 2, count_dirs(workroot), 'a later save makes no new directory')
        call assert_false(dir_exists(filepath(workroot, next_dir)), 'no '//next_dir//' after the second save')
        call assert_true(files_identical(input, string(SNAPSHOT)), 'the input project is still not written')
        call saved_proj%read(saved)
        call assert_states(saved_proj, expected_cls_states([DESELECTED, DESELECTED2]), 'second save')
        call assert_int(1, saved_proj%jobproc%get_noris(), 'a later save starts again from the input: one job')
        call saved_proj%kill
        call sel%kill
        call disp%kill
        call job_descr%kill
        call cline%kill
        deallocate(disp, saved_proj, ref_proj, cline)
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_save

    ! ---- fixture ------------------------------------------------------------

    !> six class averages in dir (class HIDDEN with pop 0) and 40 particles: particle NPTCLS without a class,
    !! particle ZERO_PTCL at state 0; every class and every other particle at state 1
    subroutine make_fixture( dir, input )
        class(string), intent(in) :: dir, input
        class(sp_project), allocatable :: proj
        type(image)  :: img
        type(string) :: stk
        integer      :: i, icls, pops(NCLS)
        allocate(proj)
        stk = filepath(dir, STK_NAME)
        call img%new([BOX,BOX,1], SMPD)
        do icls = 1,NCLS
            call img%gauran(0., 1.)
            call img%write(stk, icls)
        enddo
        call img%kill
        call proj%os_stk%new(1, is_ptcl=.false.)
        call proj%os_stk%set(1, 'stk',    'select2D_particles.mrc')
        call proj%os_stk%set(1, 'smpd',   SMPD)
        call proj%os_stk%set(1, 'box',    BOX)
        call proj%os_stk%set(1, 'nptcls', NPTCLS)
        call proj%os_stk%set(1, 'fromp',  1)
        call proj%os_stk%set(1, 'top',    NPTCLS)
        call proj%os_stk%set(1, 'ctf',    'no')
        call proj%os_stk%set_state(1, 1)
        call proj%os_ptcl2D%new(NPTCLS, is_ptcl=.true.)
        call proj%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        pops = 0
        do i = 1,NPTCLS
            call proj%os_ptcl2D%set(i, 'stkind', 1)
            call proj%os_ptcl3D%set(i, 'stkind', 1)
            if( i < NPTCLS )then
                icls = ptcl_class(i)
                call proj%os_ptcl2D%set(i, 'class', icls)
                call proj%os_ptcl3D%set(i, 'class', icls)
                pops(icls) = pops(icls) + 1
            endif
            call proj%os_ptcl2D%set_state(i, merge(0, 1, i == ZERO_PTCL))
            call proj%os_ptcl3D%set_state(i, merge(0, 1, i == ZERO_PTCL))
        enddo
        ! makes cls2D (state 1) and cls3D as a copy of it
        call proj%add_cavgs2os_out(stk, SMPD, imgkind='cavg')
        do icls = 1,NCLS
            call proj%os_cls2D%set(icls, 'pop', pops(icls))
            call proj%os_cls2D%set(icls, 'res', 4. + real(icls))
        enddo
        if( pops(HIDDEN) /= 0 ) THROW_HARD('broken fixture: the hidden class has particles')
        call proj%update_projinfo(input)
        call proj%write(input)
        call proj%kill
        deallocate(proj)
    end subroutine make_fixture

    integer function ptcl_class( i )
        integer, intent(in) :: i
        ptcl_class = CLASSES(mod(i - 1, size(CLASSES)) + 1)
    end function ptcl_class

    ! ---- expectations, from the fixture's definition --------------------------

    !> 1 for a shown class the session left selected, 0 for the deselected ones and the hidden one
    function expected_cls_states( deselected ) result( states )
        integer, intent(in) :: deselected(:)
        integer :: states(NCLS), icls
        do icls = 1,NCLS
            states(icls) = 1
            if( icls == HIDDEN .or. any(deselected == icls) ) states(icls) = 0
        enddo
    end function expected_cls_states

    !> cls2D and cls3D carry the class states; a particle takes its class's state whatever its state in the
    !! input (select2D has no append), none without a class
    subroutine assert_states( proj, cls_states, label )
        class(sp_project), intent(inout) :: proj
        integer,           intent(in)    :: cls_states(NCLS)
        character(len=*),  intent(in)    :: label
        integer :: i, icls, expected, nbad2D, nbad3D, nbadcls
        nbadcls = 0
        do icls = 1,NCLS
            if( proj%os_cls2D%get_state(icls) /= cls_states(icls) ) nbadcls = nbadcls + 1
            if( proj%os_cls3D%get_state(icls) /= cls_states(icls) ) nbadcls = nbadcls + 1
        enddo
        call assert_int(0, nbadcls, label//': cls2D and cls3D carry the selection')
        nbad2D = 0
        nbad3D = 0
        do i = 1,NPTCLS
            if( i == NPTCLS )then
                expected = 0
            else
                expected = cls_states(ptcl_class(i))
            endif
            if( proj%os_ptcl2D%get_state(i) /= expected ) nbad2D = nbad2D + 1
            if( proj%os_ptcl3D%get_state(i) /= expected ) nbad3D = nbad3D + 1
        enddo
        call assert_int(0, nbad2D, label//': ptcl2D states follow the classes')
        call assert_int(0, nbad3D, label//': ptcl3D states follow the classes')
        call assert_int(1, proj%os_ptcl2D%get_state(ZERO_PTCL), label//': a particle at state 0 in a kept class takes its class''s state')
    end subroutine assert_states

    ! ---- helpers ------------------------------------------------------------

    subroutine assert_same_records( os1, os2, label )
        class(oris),      intent(inout) :: os1, os2
        character(len=*), intent(in)    :: label
        type(string) :: str1, str2
        integer      :: i, nbad
        call assert_int(os2%get_noris(), os1%get_noris(), label//': as many records as prg=selection')
        if( os1%get_noris() /= os2%get_noris() ) return
        nbad = 0
        do i = 1,os1%get_noris()
            str1 = os1%ori2str(i)
            str2 = os2%ori2str(i)
            if( str1 /= str2 ) nbad = nbad + 1
        enddo
        call assert_int(0, nbad, label//': every record equals prg=selection''s')
    end subroutine assert_same_records

    integer function count_dirs( root )
        class(string), intent(in) :: root
        type(string), allocatable :: dirs(:)
        dirs       = simple_list_dirs(root)
        count_dirs = size(dirs)
    end function count_dirs

    integer function get_saved_len( sel )
        class(cls2D_selection), intent(in) :: sel
        type(string) :: fname
        fname         = sel%get_saved_projfile()
        get_saved_len = fname%strlen_trim()
    end function get_saved_len

    function stemname_of( fname ) result( stem )
        class(string), intent(in) :: fname
        character(len=:), allocatable :: stem
        type(string) :: s
        s    = stemname(fname)
        stem = s%to_char()
    end function stemname_of

    subroutine write_states( fname, states )
        class(string), intent(in) :: fname
        integer,       intent(in) :: states(:)
        integer :: funit, icls
        call fopen(funit, file=fname, status='replace', action='write')
        do icls = 1,size(states)
            write(funit,'(i0)') states(icls)
        enddo
        call fclose(funit)
    end subroutine write_states

end module simple_select2D_tester
