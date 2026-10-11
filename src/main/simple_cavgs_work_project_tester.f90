!@descr: unit tests of the class-average work project (simple_cavgs_work_project)
! Six generated class averages with _even/_odd companions, the third deselected in cls2D: the row layout,
! the selection, the project info, the transient sigma2 state, the nparts cap, the even/odd diagnostics
! and the file removal.
module simple_cavgs_work_project_tester
use simple_cavgs_work_project, only: cavgs_work_project
use simple_cmdline,            only: cmdline
use simple_fileio,             only: basename
use simple_image,              only: image
use simple_oris,               only: oris
use simple_parameters,         only: parameters
use simple_sp_project,         only: sp_project
use simple_string,             only: string
use simple_syslib,             only: file_exists, simple_touch
use simple_test_utils
implicit none
private
public :: run_all_cavgs_work_project_tests

integer,          parameter :: NCAVGS        = 6
integer,          parameter :: DESELECTED    = 3
integer,          parameter :: BOX           = 8
real,             parameter :: SMPD          = 2.0
character(len=*), parameter :: STK_FBODY     = 'tmp_cavgs_work_project_cavgs'
character(len=*), parameter :: WORK_PROJNAME = 'tmp_cavgs_work_project'
character(len=*), parameter :: SIGMA2_FNAME  = 'tmp_cavgs_work_project_sigma2.bin'

contains

    subroutine run_all_cavgs_work_project_tests()
        write(*,'(A)') '**** running all cavgs work project tests ****'
        call test_build()
        call test_cap_nparts()
        call test_eo_diagnostics()
    end subroutine run_all_cavgs_work_project_tests

    ! ---- fixtures -------------------------------------------------------------

    !> the merged stack and its _even/_odd companions, registered as the cavg stack, class DESELECTED
    !! deselected in cls2D; the input carries its own sigma2 state, which the work project must not inherit
    subroutine make_input( spproj )
        class(sp_project), intent(inout) :: spproj
        character(len=5), parameter :: SUFFIXES(3) = [character(len=5) :: '', '_even', '_odd']
        type(image) :: img
        integer     :: i, isfx
        call img%new([BOX,BOX,1], SMPD)
        do isfx = 1,3
            do i = 1,NCAVGS
                call img%write(string(STK_FBODY//trim(SUFFIXES(isfx))//'.mrc'), i)
            enddo
        enddo
        call img%kill
        call spproj%projinfo%new(1, is_ptcl=.false.)
        call spproj%projinfo%set(1, 'projname',     'tmp_cavgs_work_project_input')
        call spproj%projinfo%set(1, 'projfile',     'tmp_cavgs_work_project_input.simple')
        call spproj%projinfo%set(1, 'sigma2_state', 'tmp_cavgs_work_project_input_sigma2.bin')
        call spproj%add_cavgs2os_out(string(STK_FBODY//'.mrc'), SMPD, imgkind='cavg')
        call spproj%os_cls2D%new(NCAVGS, is_ptcl=.false.)
        do i = 1,NCAVGS
            call spproj%os_cls2D%set_state(i, merge(0, 1, i == DESELECTED))
        enddo
    end subroutine make_input

    ! ---- tests ----------------------------------------------------------------

    subroutine test_build()
        class(sp_project), allocatable :: spproj, work, info
        type(cavgs_work_project) :: wp
        type(string)             :: cwd_saved, fixture_root, sigma2_path
        integer :: icls, irow, nfail_before, istate
        logical :: l_found
        write(*,'(A)') 'test_build'
        nfail_before = tests_failed
        call enter_fixture('tmp_cavgs_work_project_build', cwd_saved, fixture_root)
        allocate(spproj, work, info)
        call make_input(spproj)
        ! a stale sigma2 state of an earlier run is removed, never inherited
        call simple_touch(string(SIGMA2_FNAME))
        call wp%new(spproj, string(WORK_PROJNAME//'.simple'), string(SIGMA2_FNAME))
        call assert_false(file_exists(string(SIGMA2_FNAME)), 'a stale sigma2 state is removed')
        call assert_int(NCAVGS, wp%get_ncavgs(), 'the registered class averages are counted')
        call assert_true(file_exists(string(WORK_PROJNAME//'.simple')), 'the work project is written')
        call work%read(string(WORK_PROJNAME//'.simple'))
        call assert_int(2*NCAVGS, work%os_ptcl3D%get_noris(), 'one ptcl3D row per even and per odd average')
        call assert_int(2*NCAVGS, work%os_ptcl2D%get_noris(), 'as many ptcl2D rows')
        call assert_int(2,        work%os_stk%get_noris(),    'the even and the odd stack')
        call assert_string_eq(STK_FBODY//'_even.mrc', basename(work%os_stk%get_str(1, 'stk')), 'stack 1 is the even stack')
        call assert_string_eq(STK_FBODY//'_odd.mrc',  basename(work%os_stk%get_str(2, 'stk')), 'stack 2 is the odd stack')
        call assert_string_eq('no', work%os_stk%get_str(1, 'ctf'), 'the even averages carry no CTF')
        call assert_string_eq('no', work%os_stk%get_str(2, 'ctf'), 'the odd averages carry no CTF')
        do icls = 1,NCAVGS
            istate = merge(0, 1, icls == DESELECTED)
            irow   = wp%even_row(icls)
            call assert_int(icls,   irow,                                    'even rows come first')
            call assert_int(0,      work%os_ptcl3D%get_int(irow, 'eo'),      'even row: eo 0')
            call assert_int(1,      work%os_ptcl3D%get_int(irow, 'stkind'),  'even row: the even stack')
            call assert_int(icls,   work%os_ptcl3D%get_int(irow, 'indstk'),  'even row: the class index in the stack')
            call assert_int(icls,   work%os_ptcl3D%get_int(irow, 'class'),   'even row: its class')
            call assert_int(istate, work%os_ptcl3D%get_state(irow),          'even row: the cls2D selection')
            call assert_int(istate, work%os_ptcl2D%get_state(irow),          'even row: the selection in ptcl2D')
            irow   = wp%odd_row(icls)
            call assert_int(NCAVGS + icls, irow,                             'odd rows follow')
            call assert_int(1,      work%os_ptcl3D%get_int(irow, 'eo'),      'odd row: eo 1')
            call assert_int(2,      work%os_ptcl3D%get_int(irow, 'stkind'),  'odd row: the odd stack')
            call assert_int(icls,   work%os_ptcl3D%get_int(irow, 'indstk'),  'odd row: the class index in the stack')
            call assert_int(icls,   work%os_ptcl3D%get_int(irow, 'class'),   'odd row: its class')
            call assert_int(istate, work%os_ptcl3D%get_state(irow),          'odd row: the cls2D selection')
            call assert_int(istate, work%os_ptcl2D%get_state(irow),          'odd row: the selection in ptcl2D')
        enddo
        ! read_segment leaves the project info as written; read would rewrite it from the file name
        call info%read_segment('projinfo', string(WORK_PROJNAME//'.simple'))
        call assert_string_eq(WORK_PROJNAME, info%projinfo%get_str(1, 'projname'), 'the work project is named after its file')
        call assert_string_eq(WORK_PROJNAME//'.simple', basename(info%projinfo%get_str(1, 'projfile')), 'it records its file')
        call work%get_sigma2_state_path(sigma2_path, l_found)
        call assert_true(l_found, 'the work project registers a sigma2 state')
        if( l_found ) call assert_string_eq(SIGMA2_FNAME, basename(sigma2_path), 'its own transient sigma2 state')
        call simple_touch(string(SIGMA2_FNAME))
        call wp%delete_files
        call assert_false(file_exists(string(WORK_PROJNAME//'.simple')), 'delete_files removes the work project')
        call assert_false(file_exists(string(SIGMA2_FNAME)), 'delete_files removes the sigma2 state')
        call wp%kill
        call work%kill
        call info%kill
        call spproj%kill
        call sigma2_path%kill
        deallocate(spproj, work, info)
        call leave_fixture(cwd_saved, fixture_root, nfail_before)
    end subroutine test_build

    subroutine test_cap_nparts()
        class(parameters), allocatable :: params
        class(sp_project), allocatable :: spproj
        class(cmdline),    allocatable :: cline, cline_shmem
        type(cavgs_work_project) :: wp
        type(string)             :: cwd_saved, fixture_root
        integer :: nfail_before
        write(*,'(A)') 'test_cap_nparts'
        nfail_before = tests_failed
        call enter_fixture('tmp_cavgs_work_project_nparts', cwd_saved, fixture_root)
        allocate(params, spproj, cline, cline_shmem)
        call make_input(spproj)
        call wp%new(spproj, string(WORK_PROJNAME//'.simple'), string(SIGMA2_FNAME))
        ! more parts than rows: one part per row
        params%nparts = 64
        call cline%set('nparts', 64)
        call wp%cap_nparts(params, cline)
        call assert_int(2*NCAVGS, params%nparts,            'nparts is capped at the row count')
        call assert_int(2*NCAVGS, cline%get_iarg('nparts'), 'the capped nparts is on the command line')
        call assert_int(2,        params%numlen,            'numlen follows the capped nparts')
        ! fewer parts than rows: kept
        params%nparts = 4
        call cline%set('nparts', 4)
        call wp%cap_nparts(params, cline)
        call assert_int(4, params%nparts, 'fewer parts than rows are kept')
        call assert_int(1, params%numlen, 'numlen of four parts')
        ! shared memory: nothing to cap
        params%nparts = 64
        call wp%cap_nparts(params, cline_shmem)
        call assert_int(64, params%nparts, 'a shared-memory run is left alone')
        call assert_false(cline_shmem%defined('nparts'), 'no nparts is added to a shared-memory command line')
        call wp%delete_files
        call wp%kill
        call spproj%kill
        call cline%kill
        call cline_shmem%kill
        deallocate(params, spproj, cline, cline_shmem)
        call leave_fixture(cwd_saved, fixture_root, nfail_before)
    end subroutine test_cap_nparts

    !> even and odd rows agree everywhere but in the odd row of class 1 (30 degrees off in theta) and the
    !! odd state of class 2; the deselected class counts in neither diagnostic
    subroutine test_eo_diagnostics()
        class(sp_project), allocatable :: spproj
        class(oris),       allocatable :: os
        type(cavgs_work_project) :: wp
        type(string)             :: cwd_saved, fixture_root
        integer :: icls, istate, nfail_before
        real    :: theta
        write(*,'(A)') 'test_eo_diagnostics'
        nfail_before = tests_failed
        call enter_fixture('tmp_cavgs_work_project_eo', cwd_saved, fixture_root)
        allocate(spproj, os)
        call make_input(spproj)
        call wp%new(spproj, string(WORK_PROJNAME//'.simple'), string(SIGMA2_FNAME))
        call os%new(2*NCAVGS, is_ptcl=.true.)
        do icls = 1,NCAVGS
            istate = merge(0, 1, icls == DESELECTED)
            call os%set_euler(wp%even_row(icls), [real(20*icls), 60., 0.])
            call os%set_state(wp%even_row(icls), istate)
            theta = 60.
            if( icls == 1 ) theta = 90.
            call os%set_euler(wp%odd_row(icls), [real(20*icls), theta, 0.])
            call os%set_state(wp%odd_row(icls), merge(2, istate, icls == 2))
        enddo
        call assert_real(30. / real(NCAVGS - 1), wp%eo_angular_distance(os, 'c1'), 1.e-3, &
            &'the average angular distance over the active classes')
        call assert_real(100. * real(NCAVGS - 2) / real(NCAVGS - 1), wp%eo_state_overlap(os), 1.e-3, &
            &'the state overlap over the classes with an active row')
        call wp%delete_files
        call wp%kill
        call os%kill
        call spproj%kill
        deallocate(spproj, os)
        call leave_fixture(cwd_saved, fixture_root, nfail_before)
    end subroutine test_eo_diagnostics

end module simple_cavgs_work_project_tester
