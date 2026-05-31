!@descr: unit tests for SIMPLE project merging
module simple_project_merge_tester
use simple_projfile_utils, only: merge_selected_project_files
use simple_sp_project,     only: sp_project
use simple_string,         only: string
use simple_syslib,         only: del_file, file_exists
use simple_test_utils
use simple_type_defs,      only: CTFFLAG_YES, ctfparams
implicit none

private
public :: run_all_project_merge_tests

contains

    subroutine run_all_project_merge_tests()
        write(*,'(A)') '**** running all project merge tests ****'
        call test_merge_2d_selection_heterogeneous_ctf()
    end subroutine run_all_project_merge_tests

    subroutine test_merge_2d_selection_heterogeneous_ctf()
        type(sp_project) :: proj1, proj2, merged, reread
        type(string), allocatable :: project_files(:)
        type(string) :: projfile1, projfile2, merged_file
        type(ctfparams) :: ctfvars
        integer, parameter :: NPTCLS1 = 3, NPTCLS2 = 2, NCLS = 2
        integer :: i
        write(*,'(A)') 'test_merge_2d_selection_heterogeneous_ctf'
        projfile1   = 'merge_project_src1.simple'
        projfile2   = 'merge_project_src2.simple'
        merged_file = 'merge_project_merged.simple'
        call cleanup_files(projfile1, projfile2, merged_file)
        call make_project(proj1, projfile1, NPTCLS1, NCLS, 200.0, 1.0, 0.07, &
            [1, 0, 1], [1, 2, 1], 1)
        call make_project(proj2, projfile2, NPTCLS2, NCLS, 300.0, 2.7, 0.10, &
            [0, 1], [1, 2], 1)
        allocate(project_files(2))
        project_files(1) = projfile1
        project_files(2) = projfile2
        call merge_selected_project_files(project_files, merged_file, merged, write_proj=.true.)
        call assert_true(file_exists(merged_file), 'merge_projects creates missing output project')
        call reread%read(merged_file)
        call assert_int(2, reread%os_stk%get_noris(), 'merged stack count')
        call assert_int(NPTCLS1 + NPTCLS2, reread%os_ptcl2D%get_noris(), 'merged ptcl2D count')
        call assert_int(NPTCLS1 + NPTCLS2, reread%os_ptcl3D%get_noris(), 'merged ptcl3D count')
        call assert_int(0, reread%os_cls2D%get_noris(), 'merged cls2D is intentionally empty')
        call assert_int(0, reread%os_cls3D%get_noris(), 'merged cls3D is intentionally empty')
        call assert_int(0, reread%os_out%get_noris(), 'merged out is intentionally empty')
        call assert_real(200.0, reread%os_stk%get(1, 'kv'), 1.0e-6, 'project 1 kv preserved')
        call assert_real(1.0, reread%os_stk%get(1, 'cs'), 1.0e-6, 'project 1 cs preserved')
        call assert_real(0.07, reread%os_stk%get(1, 'fraca'), 1.0e-6, 'project 1 fraca preserved')
        call assert_real(300.0, reread%os_stk%get(2, 'kv'), 1.0e-6, 'project 2 kv preserved')
        call assert_real(2.7, reread%os_stk%get(2, 'cs'), 1.0e-6, 'project 2 cs preserved')
        call assert_real(0.10, reread%os_stk%get(2, 'fraca'), 1.0e-6, 'project 2 fraca preserved')
        call assert_int(1, reread%os_stk%get_fromp(1), 'project 1 stack fromp')
        call assert_int(NPTCLS1, reread%os_stk%get_top(1), 'project 1 stack top')
        call assert_int(NPTCLS1 + 1, reread%os_stk%get_fromp(2), 'project 2 stack fromp remapped')
        call assert_int(NPTCLS1 + NPTCLS2, reread%os_stk%get_top(2), 'project 2 stack top remapped')
        call assert_int(1, reread%os_ptcl2D%get_int(1, 'stkind'), 'project 1 particle stkind')
        call assert_int(2, reread%os_ptcl2D%get_int(NPTCLS1 + 1, 'stkind'), 'project 2 particle stkind remapped')
        do i = 1,NPTCLS1
            call assert_int(1, reread%os_ptcl2D%get_int(i, 'stkind'), 'project 1 particle stkind range')
            call assert_int(i, reread%os_ptcl2D%get_int(i, 'indstk'), 'project 1 particle indstk preserved')
        enddo
        do i = 1,NPTCLS2
            call assert_int(2, reread%os_ptcl2D%get_int(NPTCLS1 + i, 'stkind'), &
                &'project 2 particle stkind range remapped')
            call assert_int(i, reread%os_ptcl2D%get_int(NPTCLS1 + i, 'indstk'), &
                &'project 2 particle indstk preserved')
        enddo
        call assert_int(0, reread%os_ptcl2D%get_class(1), 'project 1 particle class reset')
        call assert_int(0, reread%os_ptcl2D%get_class(NPTCLS1 + 1), 'project 2 particle class reset')
        call assert_int(1, reread%os_ptcl2D%get_state(1), 'project 1 selected ptcl2D state preserved')
        call assert_int(0, reread%os_ptcl2D%get_state(2), 'project 1 rejected ptcl2D state preserved')
        call assert_int(0, reread%os_ptcl2D%get_state(NPTCLS1 + 1), 'project 2 rejected ptcl2D state preserved')
        call assert_int(1, reread%os_ptcl2D%get_state(NPTCLS1 + 2), 'project 2 selected ptcl2D state preserved')
        call assert_int(1, reread%os_ptcl3D%get_state(1), 'project 1 ptcl3D state preserved')
        call assert_int(0, reread%os_ptcl3D%get_state(NPTCLS1 + 1), 'project 2 ptcl3D state preserved')
        call assert_int(1, reread%os_ptcl2D%get_int(1, 'ogid'), 'project 1 particle ogid')
        call assert_int(2, reread%os_ptcl2D%get_int(NPTCLS1 + 1, 'ogid'), 'project 2 particle ogid remapped')
        call assert_int(2, reread%os_stk%get_int(2, 'ogid'), 'project 2 stack ogid remapped')
        call assert_int(0, reread%os_optics%get_noris(), 'merge does not require os_optics')
        ctfvars = reread%get_ctfparams('ptcl2D', 1)
        call assert_true(ctfvars%ctfflag == CTFFLAG_YES, 'project 1 CTF flag resolved')
        call assert_real(200.0, ctfvars%kv, 1.0e-6, 'project 1 particle CTF kv resolved')
        call assert_real(1.0, ctfvars%cs, 1.0e-6, 'project 1 particle CTF cs resolved')
        call assert_real(0.07, ctfvars%fraca, 1.0e-6, 'project 1 particle CTF fraca resolved')
        ctfvars = reread%get_ctfparams('ptcl2D', NPTCLS1 + 1)
        call assert_true(ctfvars%ctfflag == CTFFLAG_YES, 'project 2 CTF flag resolved')
        call assert_real(300.0, ctfvars%kv, 1.0e-6, 'project 2 particle CTF kv resolved')
        call assert_real(2.7, ctfvars%cs, 1.0e-6, 'project 2 particle CTF cs resolved')
        call assert_real(0.10, ctfvars%fraca, 1.0e-6, 'project 2 particle CTF fraca resolved')
        ctfvars = reread%get_ctfparams('ptcl3D', NPTCLS1 + 1)
        call assert_real(300.0, ctfvars%kv, 1.0e-6, 'project 2 ptcl3D CTF kv resolved')
        call assert_real(2.7, ctfvars%cs, 1.0e-6, 'project 2 ptcl3D CTF cs resolved')
        call assert_real(0.10, ctfvars%fraca, 1.0e-6, 'project 2 ptcl3D CTF fraca resolved')
        call cleanup_files(projfile1, projfile2, merged_file)
        call proj1%kill
        call proj2%kill
        call merged%kill
        call reread%kill
        if( allocated(project_files) ) deallocate(project_files)
    end subroutine test_merge_2d_selection_heterogeneous_ctf

    subroutine make_project(proj, projfile, nptcls, ncls, kv, cs, fraca, states, classes, ogid)
        type(sp_project), intent(inout) :: proj
        type(string),     intent(in)    :: projfile
        integer,          intent(in)    :: nptcls, ncls, states(:), classes(:), ogid
        real,             intent(in)    :: kv, cs, fraca
        integer :: i
        call proj%kill
        call proj%os_stk%new(1, is_ptcl=.false.)
        call proj%os_stk%set(1, 'stk',        projfile%to_char()//'.mrcs')
        call proj%os_stk%set(1, 'ctf',        'yes')
        call proj%os_stk%set(1, 'phaseplate', 'no')
        call proj%os_stk%set(1, 'smpd',       1.25)
        call proj%os_stk%set(1, 'kv',         kv)
        call proj%os_stk%set(1, 'cs',         cs)
        call proj%os_stk%set(1, 'fraca',      fraca)
        call proj%os_stk%set(1, 'box',        128)
        call proj%os_stk%set(1, 'nptcls',     nptcls)
        call proj%os_stk%set(1, 'fromp',      1)
        call proj%os_stk%set(1, 'top',        nptcls)
        call proj%os_stk%set_ogid(1, ogid)
        call proj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        call proj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        do i = 1,nptcls
            call fill_particle_row(proj%os_ptcl2D, i, states(i), classes(i), ogid)
            call fill_particle_row(proj%os_ptcl3D, i, states(i), 0,          ogid)
        enddo
        call proj%os_cls2D%new(ncls, is_ptcl=.false.)
        call proj%os_cls3D%new(ncls, is_ptcl=.false.)
        do i = 1,ncls
            call fill_class_row(proj%os_cls2D, i, i, modulo(i, 2), ogid)
            call fill_class_row(proj%os_cls3D, i, i, modulo(i, 2), ogid)
        enddo
        call proj%os_out%new(1, is_ptcl=.false.)
        call proj%os_out%set(1, 'imgkind', 'cavg')
        call proj%os_out%set(1, 'stk',     projfile%to_char()//'_cavgs.mrcs')
        call proj%os_out%set(1, 'nptcls',  ncls)
        call proj%os_out%set(1, 'fromp',   1)
        call proj%os_out%set(1, 'top',     ncls)
        call proj%os_out%set_ogid(1, ogid)
        call proj%update_projinfo(projfile)
        call proj%write(projfile)
    end subroutine make_project

    subroutine fill_particle_row(os, i, state, cls, ogid)
        use simple_oris, only: oris
        type(oris), intent(inout) :: os
        integer,    intent(in)    :: i, state, cls, ogid
        call os%set_stkind(i, 1)
        call os%set(i, 'indstk', i)
        call os%set_state(i, state)
        call os%set_class(i, cls)
        call os%set_dfx(i, 1.50 + 0.01 * real(i))
        call os%set_dfy(i, 1.60 + 0.01 * real(i))
        call os%set(i, 'angast', 10.0 * real(i))
        call os%set(i, 'phshift', 0.0)
        call os%set(i, 'corr', 0.5 + 0.01 * real(i))
        call os%set(i, 'sampled', i)
        call os%set(i, 'updatecnt', i + 10)
        call os%set_ogid(i, ogid)
    end subroutine fill_particle_row

    subroutine fill_class_row(os, i, cls, state, ogid)
        use simple_oris, only: oris
        type(oris), intent(inout) :: os
        integer,    intent(in)    :: i, cls, state, ogid
        call os%set_class(i, cls)
        call os%set_state(i, state)
        call os%set(i, 'pop', 100 + i)
        call os%set(i, 'corr', 0.25 + 0.01 * real(i))
        call os%set_ogid(i, ogid)
    end subroutine fill_class_row

    subroutine cleanup_files(projfile1, projfile2, merged_file)
        type(string), intent(in) :: projfile1, projfile2, merged_file
        call del_file(projfile1)
        call del_file(projfile2)
        call del_file(merged_file)
    end subroutine cleanup_files

end module simple_project_merge_tester
