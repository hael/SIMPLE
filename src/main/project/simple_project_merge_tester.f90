!@descr: unit tests for SIMPLE project merging
module simple_project_merge_tester
use simple_projfile_utils, only: merge_selected_project_files, fix_project_file, remap_project_paths
use simple_sp_project,     only: sp_project
use simple_fileio,         only: filepath
use simple_string,         only: string
use simple_syslib,         only: del_file, file_exists, simple_getcwd, simple_mkdir, simple_rmdir
use simple_test_utils
implicit none

private
public :: run_all_project_merge_tests

contains

    subroutine run_all_project_merge_tests()
        write(*,'(A)') '**** running all project merge tests ****'
        call test_merge_pruned_stack_indexing_heterogeneous_ctf()
        call test_fix_projfile_repairs_unpruned_legacy_stack()
        call test_fix_projfile_repairs_wrong_segment_stkind()
        call test_fix_projfile_refuses_particles_without_stacks()
        call test_remap_project_paths()
        call test_remap_project_paths_allows_unmatched_scope()
        call test_remap_project_paths_scoped_roots()
        call test_sigma2_state_path_registration()
    end subroutine run_all_project_merge_tests

    subroutine test_merge_pruned_stack_indexing_heterogeneous_ctf()
        type(sp_project) :: proj1, proj2, merged, reread
        type(string), allocatable :: project_files(:)
        type(string) :: projfile1, projfile2, merged_file
        integer, parameter :: NPTCLS1 = 3, NPTCLS2 = 2, NCLS = 2
        integer, parameter :: NPTCLS_STK1 = 6, NPTCLS_STK2 = 4
        integer, parameter :: INDSTKS1(NPTCLS1) = [1, 4, 6]
        integer, parameter :: INDSTKS2(NPTCLS2) = [4, 2]
        integer :: i, stkind, ind_in_stk
        write(*,'(A)') 'test_merge_pruned_stack_indexing_heterogeneous_ctf'
        projfile1   = 'merge_project_src1.simple'
        projfile2   = 'merge_project_src2.simple'
        merged_file = 'merge_project_merged.simple'
        call cleanup_files(projfile1, projfile2, merged_file)
        call make_project(proj1, projfile1, NPTCLS1, NPTCLS_STK1, NCLS, 200.0, 1.0, 0.07, &
            [1, 1, 1], [1, 2, 1], INDSTKS1, 1)
        call make_project(proj2, projfile2, NPTCLS2, NPTCLS_STK2, NCLS, 300.0, 2.7, 0.10, &
            [1, 1], [1, 2], INDSTKS2, 1)
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
        call assert_string_eq('yes', reread%os_stk%get_str(1, 'ctf'), 'project 1 ctf flag preserved')
        call assert_real(0.0, reread%os_ptcl2D%get(1, 'phshift'), 1.0e-6, &
            &'project 1 numerical phase preserved')
        call assert_real(200.0, reread%os_stk%get(1, 'kv'), 1.0e-6, 'project 1 kv preserved')
        call assert_real(1.0, reread%os_stk%get(1, 'cs'), 1.0e-6, 'project 1 cs preserved')
        call assert_real(0.07, reread%os_stk%get(1, 'fraca'), 1.0e-6, 'project 1 fraca preserved')
        call assert_string_eq('yes', reread%os_stk%get_str(2, 'ctf'), 'project 2 ctf flag preserved')
        call assert_real(300.0, reread%os_stk%get(2, 'kv'), 1.0e-6, 'project 2 kv preserved')
        call assert_real(2.7, reread%os_stk%get(2, 'cs'), 1.0e-6, 'project 2 cs preserved')
        call assert_real(0.10, reread%os_stk%get(2, 'fraca'), 1.0e-6, 'project 2 fraca preserved')
        call assert_int(NPTCLS1, reread%os_stk%get_int(1, 'nptcls'), 'project 1 project particle count')
        call assert_int(NPTCLS_STK1, reread%os_stk%get_int(1, 'nptcls_stk'), 'project 1 physical stack count')
        call assert_int(NPTCLS2, reread%os_stk%get_int(2, 'nptcls'), 'project 2 project particle count')
        call assert_int(NPTCLS_STK2, reread%os_stk%get_int(2, 'nptcls_stk'), 'project 2 physical stack count')
        call assert_int(1, reread%os_stk%get_fromp(1), 'project 1 stack fromp')
        call assert_int(NPTCLS1, reread%os_stk%get_top(1), 'project 1 stack top')
        call assert_int(NPTCLS1 + 1, reread%os_stk%get_fromp(2), 'project 2 stack fromp remapped')
        call assert_int(NPTCLS1 + NPTCLS2, reread%os_stk%get_top(2), 'project 2 stack top remapped')
        call assert_int(1, reread%os_ptcl2D%get_int(1, 'stkind'), 'project 1 particle stkind')
        call assert_int(2, reread%os_ptcl2D%get_int(NPTCLS1 + 1, 'stkind'), 'project 2 particle stkind remapped')
        do i = 1,NPTCLS1
            call assert_int(1, reread%os_ptcl2D%get_int(i, 'stkind'), 'project 1 particle stkind range')
            call assert_int(INDSTKS1(i), reread%os_ptcl2D%get_int(i, 'indstk'), &
                &'project 1 physical indstk preserved')
            call reread%map_ptcl_ind2stk_ind('ptcl2D', i, stkind, ind_in_stk)
            call assert_int(1, stkind, 'project 1 mapped stkind')
            call assert_int(INDSTKS1(i), ind_in_stk, 'project 1 mapped physical indstk')
        enddo
        do i = 1,NPTCLS2
            call assert_int(2, reread%os_ptcl2D%get_int(NPTCLS1 + i, 'stkind'), &
                &'project 2 particle stkind range remapped')
            call assert_int(INDSTKS2(i), reread%os_ptcl2D%get_int(NPTCLS1 + i, 'indstk'), &
                &'project 2 physical indstk preserved')
            call reread%map_ptcl_ind2stk_ind('ptcl2D', NPTCLS1 + i, stkind, ind_in_stk)
            call assert_int(2, stkind, 'project 2 mapped stkind')
            call assert_int(INDSTKS2(i), ind_in_stk, 'project 2 mapped physical indstk')
            call reread%map_ptcl_ind2stk_ind('ptcl3D', NPTCLS1 + i, stkind, ind_in_stk)
            call assert_int(2, stkind, 'project 2 ptcl3D mapped stkind')
            call assert_int(INDSTKS2(i), ind_in_stk, 'project 2 ptcl3D mapped physical indstk')
        enddo
        call assert_int(0, reread%os_ptcl2D%get_class(1), 'project 1 particle class reset')
        call assert_int(0, reread%os_ptcl2D%get_class(NPTCLS1 + 1), 'project 2 particle class reset')
        call assert_int(1, reread%os_ptcl2D%get_state(1), 'project 1 ptcl2D state preserved')
        call assert_int(1, reread%os_ptcl2D%get_state(NPTCLS1), 'project 1 final ptcl2D state preserved')
        call assert_int(1, reread%os_ptcl2D%get_state(NPTCLS1 + 1), 'project 2 ptcl2D state preserved')
        call assert_int(1, reread%os_ptcl2D%get_state(NPTCLS1 + NPTCLS2), 'project 2 final ptcl2D state preserved')
        call assert_int(1, reread%os_ptcl3D%get_state(1), 'project 1 ptcl3D state preserved')
        call assert_int(1, reread%os_ptcl3D%get_state(NPTCLS1 + 1), 'project 2 ptcl3D state preserved')
        call assert_int(1, reread%os_ptcl2D%get_int(1, 'ogid'), 'project 1 particle ogid')
        call assert_int(2, reread%os_ptcl2D%get_int(NPTCLS1 + 1, 'ogid'), 'project 2 particle ogid remapped')
        call assert_int(2, reread%os_stk%get_int(2, 'ogid'), 'project 2 stack ogid remapped')
        call assert_int(0, reread%os_optics%get_noris(), 'merge does not require os_optics')
        call cleanup_files(projfile1, projfile2, merged_file)
        call proj1%kill
        call proj2%kill
        call merged%kill
        call reread%kill
        if( allocated(project_files) ) deallocate(project_files)
    end subroutine test_merge_pruned_stack_indexing_heterogeneous_ctf

    subroutine test_fix_projfile_repairs_unpruned_legacy_stack()
        ! A project written before the stack-index fix: the stack row has no
        ! nptcls_stk, one particle has no indstk and one an out-of-range value,
        ! and the stack range is wrong. The stack file holds one image per project
        ! row, so fix_projfile can prove the physical indices: nptcls_stk comes from
        ! the header and indstk from the project rows.
        use simple_image, only: image
        type(sp_project) :: proj, fixed
        type(image)      :: img
        type(string)     :: projfile, fixed_file, stkfile, staged_file
        integer, parameter :: NPTCLS = 3, NCLS = 2
        integer, parameter :: LEGACY_INDSTKS(NPTCLS) = [99, 0, 3]
        integer :: i, stkind, ind_in_stk
        write(*,'(A)') 'test_fix_projfile_repairs_unpruned_legacy_stack'
        projfile   = 'fix_project_src.simple'
        fixed_file = 'fix_project_src_fixed.simple'
        staged_file = 'fix_project_src.tmp'
        stkfile    = 'fix_project_src_stack.mrc'
        call del_file(projfile)
        call del_file(fixed_file)
        call del_file(staged_file)
        call del_file(stkfile)
        call img%new([8,8,1], 1.25)
        do i = 1,NPTCLS
            call img%write(stkfile, i)
        enddo
        call img%kill
        call make_project(proj, projfile, NPTCLS, 0, NCLS, 200.0, 1.0, 0.07, &
            [1, 0, 1], [1, 2, 1], LEGACY_INDSTKS, 1)
        call proj%os_stk%set(1, 'stk', stkfile)
        call proj%os_ptcl2D%delete_entry(2, 'indstk')
        call proj%os_ptcl3D%delete_entry(2, 'indstk')
        call proj%os_stk%set(1, 'top', NPTCLS + 5)   ! a stack range the tool repairs
        call proj%write(projfile)
        call assert_false(proj%os_stk%isthere(1, 'nptcls_stk'), 'legacy source lacks nptcls_stk')
        call fix_project_file(projfile, fixed_file)
        call assert_true(fixed_file == projfile, 'fix_projfile returns the replaced project path')
        call assert_false(file_exists('fix_project_src_fixed.simple'), 'fix_projfile leaves no fixed copy')
        call assert_false(file_exists(staged_file), 'fix_projfile removes its staged project')
        call fixed%read(projfile)
        call assert_int(1, fixed%os_stk%get_noris(), 'fixed stack count')
        call assert_int(NPTCLS, fixed%os_ptcl2D%get_noris(), 'fixed ptcl2D count includes state 0')
        call assert_int(NPTCLS, fixed%os_ptcl3D%get_noris(), 'fixed ptcl3D count includes state 0')
        call assert_int(1, fixed%os_stk%get_fromp(1), 'fixed fromp')
        call assert_int(NPTCLS, fixed%os_stk%get_top(1), 'fixed top repaired')
        call assert_int(NPTCLS, fixed%os_stk%get_int(1, 'nptcls'), 'fixed project nptcls')
        call assert_int(NPTCLS, fixed%os_stk%get_int(1, 'nptcls_stk'), 'nptcls_stk taken from the stack header')
        do i = 1,NPTCLS
            call assert_int(1, fixed%os_ptcl2D%get_int(i, 'stkind'), 'fixed ptcl2D stkind')
            call assert_int(i, fixed%os_ptcl2D%get_int(i, 'indstk'), 'fixed ptcl2D indstk from the project row')
            call assert_int(i, fixed%os_ptcl3D%get_int(i, 'indstk'), 'fixed ptcl3D indstk from the project row')
            call fixed%map_ptcl_ind2stk_ind('ptcl2D', i, stkind, ind_in_stk)
            call assert_int(1, stkind, 'fixed project maps stkind')
            call assert_int(i, ind_in_stk, 'fixed project maps the physical index')
        enddo
        call assert_int(0, fixed%os_ptcl2D%get_state(2), 'state 0 row remains present after the fix')
        call del_file(projfile)
        call del_file(stkfile)
        call proj%kill
        call fixed%kill
    end subroutine test_fix_projfile_repairs_unpruned_legacy_stack

    subroutine test_fix_projfile_repairs_wrong_segment_stkind()
        ! Two unpruned legacy stacks of two images each. ptcl2D is consistent, but
        ! ptcl3D row 1 claims stack 2: a valid stack index that does not own the
        ! row. fix_projfile must replace it by the owning stack before deriving
        ! indstk, so that both segments agree and both map.
        use simple_image, only: image
        type(sp_project) :: proj, fixed
        type(image)      :: img
        type(string)     :: projfile, fixed_file, staged_file, stkfiles(2)
        integer, parameter :: NSTK = 2, NPER = 2, NPTCLS = NSTK*NPER
        integer :: i, istk, stkind, ind_in_stk, nerrors
        write(*,'(A)') 'test_fix_projfile_repairs_wrong_segment_stkind'
        projfile    = 'fix_project_two_stacks.simple'
        fixed_file  = 'fix_project_two_stacks_fixed.simple'
        staged_file = 'fix_project_two_stacks.tmp'
        stkfiles(1) = 'fix_project_two_stacks_1.mrc'
        stkfiles(2) = 'fix_project_two_stacks_2.mrc'
        call del_file(projfile)
        call del_file(fixed_file)
        call del_file(staged_file)
        call img%new([8,8,1], 1.25)
        do istk = 1,NSTK
            call del_file(stkfiles(istk))
            do i = 1,NPER
                call img%write(stkfiles(istk), i)
            enddo
        enddo
        call img%kill
        call proj%os_stk%new(NSTK, is_ptcl=.false.)
        do istk = 1,NSTK
            call proj%os_stk%set(istk, 'stk',    stkfiles(istk))
            call proj%os_stk%set(istk, 'ctf',    'no')
            call proj%os_stk%set(istk, 'smpd',   1.25)
            call proj%os_stk%set(istk, 'box',    8)
            call proj%os_stk%set(istk, 'nptcls', NPER)
            call proj%os_stk%set(istk, 'fromp',  (istk-1)*NPER + 1)
            call proj%os_stk%set(istk, 'top',    istk*NPER)
        enddo
        call proj%os_ptcl2D%new(NPTCLS, is_ptcl=.true.)
        call proj%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        do i = 1,NPTCLS
            istk = (i-1)/NPER + 1
            call proj%os_ptcl2D%set_stkind(i, istk)
            call proj%os_ptcl3D%set_stkind(i, istk)
            call proj%os_ptcl2D%set_state(i, 1)
            call proj%os_ptcl3D%set_state(i, 1)
        enddo
        call proj%os_ptcl3D%set_stkind(1, 2)   ! in range, but stack 2 does not own row 1
        call proj%update_projinfo(projfile)
        call proj%write(projfile)
        call fix_project_file(projfile, fixed_file, nerrors)
        call assert_int(0, nerrors, 'two-stack legacy project is fixable')
        call assert_true(fixed_file == projfile, 'fix_projfile returns the replaced two-stack project path')
        call assert_false(file_exists('fix_project_two_stacks_fixed.simple'), 'fix_projfile leaves no two-stack copy')
        call assert_false(file_exists(staged_file), 'fix_projfile removes its two-stack staged project')
        call fixed%read(projfile)
        do i = 1,NPTCLS
            istk = (i-1)/NPER + 1
            call assert_int(istk, fixed%os_ptcl2D%get_int(i, 'stkind'), 'fixed ptcl2D stkind')
            call assert_int(istk, fixed%os_ptcl3D%get_int(i, 'stkind'), 'fixed ptcl3D stkind is the owning stack')
            call assert_int(fixed%os_ptcl2D%get_int(i, 'indstk'), fixed%os_ptcl3D%get_int(i, 'indstk'), &
                &'fixed segments agree on indstk')
            call fixed%map_ptcl_ind2stk_ind('ptcl2D', i, stkind, ind_in_stk)
            call assert_int(istk, stkind, 'fixed ptcl2D maps its stack')
            call assert_int(i - (istk-1)*NPER, ind_in_stk, 'fixed ptcl2D maps its physical index')
            call fixed%map_ptcl_ind2stk_ind('ptcl3D', i, stkind, ind_in_stk)
            call assert_int(istk, stkind, 'fixed ptcl3D maps its stack')
            call assert_int(i - (istk-1)*NPER, ind_in_stk, 'fixed ptcl3D maps its physical index')
        enddo
        call del_file(projfile)
        do istk = 1,NSTK
            call del_file(stkfiles(istk))
        enddo
        call proj%kill
        call fixed%kill
    end subroutine test_fix_projfile_repairs_wrong_segment_stkind

    subroutine test_fix_projfile_refuses_particles_without_stacks()
        ! Particle rows without any stack row cannot be mapped to images:
        ! fix_projfile reports the error and leaves the input unchanged.
        type(sp_project) :: proj
        type(string)     :: projfile, fixed_file, staged_file
        integer :: i, nerrors
        write(*,'(A)') 'test_fix_projfile_refuses_particles_without_stacks'
        projfile   = 'fix_project_no_stacks.simple'
        fixed_file = 'fix_project_no_stacks_fixed.simple'
        staged_file = 'fix_project_no_stacks.tmp'
        call del_file(projfile)
        call del_file(fixed_file)
        call del_file(staged_file)
        call proj%os_ptcl2D%new(3, is_ptcl=.true.)
        call proj%os_ptcl3D%new(3, is_ptcl=.true.)
        do i = 1,3
            call proj%os_ptcl2D%set_state(i, 1)
            call proj%os_ptcl3D%set_state(i, 1)
        enddo
        call proj%update_projinfo(projfile)
        call proj%write(projfile)
        call fix_project_file(projfile, fixed_file, nerrors)
        call assert_true(nerrors > 0, 'particles without stack rows are an error')
        call assert_true(file_exists(projfile), 'unfixable input project remains in place')
        call assert_false(file_exists('fix_project_no_stacks_fixed.simple'), 'unfixable project leaves no fixed copy')
        call assert_false(file_exists(staged_file), 'unfixable project leaves no staged project')
        call del_file(projfile)
        call proj%kill
    end subroutine test_fix_projfile_refuses_particles_without_stacks

    subroutine test_remap_project_paths()
        type(sp_project) :: proj
        type(string) :: old_root, new_root, movie_file, intg_file, stk_file, box_file
        integer :: nremapped, status
        write(*,'(A)') 'test_remap_project_paths'
        old_root = '/retired/dataset/'
        new_root = 'remap_project_paths_target/'
        call simple_rmdir(new_root, status)
        call simple_mkdir(new_root)
        call simple_mkdir(new_root//'movies')
        call simple_mkdir(new_root//'micrographs')
        call simple_mkdir(new_root//'particles')
        call simple_mkdir(new_root//'boxes')
        movie_file = new_root//'movies/movie_001.eer'
        intg_file  = new_root//'micrographs/micrograph_001.mrc'
        stk_file   = new_root//'particles/particles_001.mrcs'
        box_file   = new_root//'boxes/micrograph_001.box'
        call create_empty_file(movie_file)
        call create_empty_file(intg_file)
        call create_empty_file(stk_file)
        call create_empty_file(box_file)
        call proj%os_mic%new(1, is_ptcl=.false.)
        call proj%os_mic%set(1, 'movie', '/retired/dataset/movies/movie_001.eer')
        call proj%os_mic%set(1, 'intg',  '/retired/dataset/micrographs/micrograph_001.mrc')
        call proj%os_stk%new(2, is_ptcl=.false.)
        call proj%os_stk%set(1, 'stk',     '/retired/dataset/particles/particles_001.mrcs')
        call proj%os_stk%set(1, 'boxfile', '/retired/dataset/boxes/micrograph_001.box')
        ! A textual prefix without a path boundary must remain unchanged.
        call proj%os_stk%set(2, 'stk', '/retired/dataset_backup/leave_unchanged.mrcs')
        call remap_project_paths(proj, old_root, new_root, nremapped)
        call assert_int(4, nremapped, 'remap_project_paths updates supported path fields')
        call assert_string_eq(movie_file%to_char(), proj%os_mic%get_str(1, 'movie'), 'movie path remapped')
        call assert_string_eq(intg_file%to_char(), proj%os_mic%get_str(1, 'intg'), 'intg path remapped')
        call assert_string_eq(stk_file%to_char(), proj%os_stk%get_str(1, 'stk'), 'stack path remapped')
        call assert_string_eq(box_file%to_char(), proj%os_stk%get_str(1, 'boxfile'), 'boxfile path remapped')
        call assert_string_eq('/retired/dataset_backup/leave_unchanged.mrcs', &
            &proj%os_stk%get_str(2, 'stk'), 'path-component boundary is respected')
        call proj%kill
        call simple_rmdir(new_root, status)
        call assert_int(0, status, 'remap_project_paths test directory cleanup')
    end subroutine test_remap_project_paths

    subroutine test_remap_project_paths_allows_unmatched_scope()
        type(sp_project) :: proj
        type(string) :: new_root
        integer :: nremapped, status
        write(*,'(A)') 'test_remap_project_paths_allows_unmatched_scope'
        new_root = 'remap_project_paths_unmatched_target/'
        call simple_rmdir(new_root, status)
        call simple_mkdir(new_root)
        call proj%os_mic%new(1, is_ptcl=.false.)
        call proj%os_mic%set(1, 'movie', '/different/root/movie.eer')
        call remap_project_paths(proj, string('/not/present'), new_root, nremapped, &
            &scope='mic', require_match=.false.)
        call assert_int(0, nremapped, 'unmatched global fallback scope is skipped')
        call assert_string_eq('/different/root/movie.eer', proj%os_mic%get_str(1, 'movie'), &
            &'unmatched scope leaves paths unchanged')
        call proj%kill
        call simple_rmdir(new_root, status)
        call assert_int(0, status, 'unmatched scope test directory cleanup')
    end subroutine test_remap_project_paths_allows_unmatched_scope

    subroutine test_remap_project_paths_scoped_roots()
        type(sp_project) :: proj
        type(string) :: root, mic_root, ptcl_root, cavg_root, cavg_path, vol_root
        type(string) :: movie_file, intg_file, mic_box_file
        type(string) :: stk_file, stk_box_file
        type(string) :: cavg_file, frcs2D_file, sigma2_file, vol_file, fsc_file, frcs3D_file
        integer :: nremapped, status
        write(*,'(A)') 'test_remap_project_paths_scoped_roots'
        root      = 'remap_project_paths_scoped_target/'
        mic_root  = root//'mic/'
        ptcl_root = root//'ptcl/'
        cavg_path = root//'cavg'
        cavg_root = cavg_path//'/'
        vol_root  = root//'vol/'
        call simple_rmdir(root, status)
        call simple_mkdir(root)
        call simple_mkdir(mic_root)
        call simple_mkdir(ptcl_root)
        call simple_mkdir(cavg_root)
        call simple_mkdir(vol_root)
        movie_file   = mic_root//'movie.eer'
        intg_file    = mic_root//'intg.mrc'
        mic_box_file = mic_root//'mic.box'
        stk_file     = ptcl_root//'raw.mrcs'
        stk_box_file = ptcl_root//'ptcl.box'
        cavg_file    = cavg_root//'classes.mrcs'
        frcs2D_file  = cavg_root//'frcs2D.bin'
        sigma2_file  = cavg_root//'sigma2.bin'
        vol_file     = vol_root//'map.mrc'
        fsc_file     = vol_root//'fsc.bin'
        frcs3D_file  = vol_root//'frcs3D.bin'
        call create_empty_file(movie_file)
        call create_empty_file(intg_file)
        call create_empty_file(mic_box_file)
        call create_empty_file(stk_file)
        call create_empty_file(stk_box_file)
        call create_empty_file(cavg_file)
        call create_empty_file(frcs2D_file)
        call create_empty_file(sigma2_file)
        call create_empty_file(vol_file)
        call create_empty_file(fsc_file)
        call create_empty_file(frcs3D_file)
        call proj%os_mic%new(1, is_ptcl=.false.)
        call proj%os_mic%set(1, 'movie',   '/legacy/mic/movie.eer')
        call proj%os_mic%set(1, 'intg',    '/legacy/mic/intg.mrc')
        call proj%os_mic%set(1, 'boxfile', '/legacy/mic/mic.box')
        call proj%os_stk%new(1, is_ptcl=.false.)
        call proj%os_stk%set(1, 'stk',     '/legacy/ptcl/raw.mrcs')
        call proj%os_stk%set(1, 'boxfile', '/legacy/ptcl/ptcl.box')
        call proj%os_ptcl3D%new(1, is_ptcl=.true.)
        call proj%os_ptcl3D%set_stkind(1, 1)
        call proj%os_out%new(3, is_ptcl=.false.)
        call proj%os_out%set(1, 'stk',     '/legacy/cavg/classes.mrcs')
        call proj%os_out%set(1, 'stkpath', '/legacy/cavg')
        call proj%os_out%set(1, 'imgkind', 'cavg')
        call proj%os_out%set(1, 'sigma2',  '/legacy/cavg/sigma2.bin')
        call proj%os_out%set(2, 'imgkind', 'frc2D')
        call proj%os_out%set(2, 'frcs',    '/legacy/cavg/frcs2D.bin')
        call proj%os_out%set(3, 'imgkind', 'frc3D')
        call proj%os_out%set(3, 'frcs',    '/legacy/vol/frcs3D.bin')
        call proj%os_out%set(3, 'vol',     '/legacy/vol/map.mrc')
        call proj%os_out%set(3, 'fsc',     '/legacy/vol/fsc.bin')
        call remap_project_paths(proj, string('/legacy/mic'), mic_root, nremapped, scope='mic')
        call assert_int(3, nremapped, 'scoped mic path count')
        call remap_project_paths(proj, string('/legacy/ptcl'), ptcl_root, nremapped, scope='ptcl')
        call assert_int(2, nremapped, 'scoped particle path count')
        call remap_project_paths(proj, string('/legacy/cavg'), cavg_root, nremapped, scope='cavg')
        call assert_int(4, nremapped, 'scoped class-average path count')
        call remap_project_paths(proj, string('/legacy/vol'), vol_root, nremapped, scope='vol')
        call assert_int(3, nremapped, 'scoped volume path count')
        call assert_string_eq(movie_file%to_char(), proj%os_mic%get_str(1, 'movie'), 'scoped movie path')
        call assert_string_eq(intg_file%to_char(), proj%os_mic%get_str(1, 'intg'), 'scoped intg path')
        call assert_string_eq(mic_box_file%to_char(), proj%os_mic%get_str(1, 'boxfile'), 'scoped mic box path')
        call assert_string_eq(stk_file%to_char(), proj%os_stk%get_str(1, 'stk'), 'scoped raw stack path')
        call assert_string_eq(stk_box_file%to_char(), proj%os_stk%get_str(1, 'boxfile'), 'scoped stack box path')
        call assert_string_eq(cavg_file%to_char(), proj%os_out%get_str(1, 'stk'), 'scoped class-average path')
        call assert_string_eq(cavg_path%to_char(), proj%os_out%get_str(1, 'stkpath'), 'scoped class path')
        call assert_string_eq(frcs2D_file%to_char(), proj%os_out%get_str(2, 'frcs'), 'scoped 2D FRC path')
        call assert_string_eq(sigma2_file%to_char(), proj%os_out%get_str(1, 'sigma2'), 'scoped sigma2 path')
        call assert_string_eq(vol_file%to_char(), proj%os_out%get_str(3, 'vol'), 'scoped volume path')
        call assert_string_eq(fsc_file%to_char(), proj%os_out%get_str(3, 'fsc'), 'scoped FSC path')
        call assert_string_eq(frcs3D_file%to_char(), proj%os_out%get_str(3, 'frcs'), 'scoped 3D FRC path')
        call proj%kill
        call simple_rmdir(root, status)
        call assert_int(0, status, 'scoped remap test directory cleanup')
    end subroutine test_remap_project_paths_scoped_roots

    subroutine test_sigma2_state_path_registration()
        type(sp_project) :: proj
        type(string) :: requested, registered, cwd, expected
        logical :: found
        write(*,'(A)') 'test_sigma2_state_path_registration'
        requested = 'sigma2_state.bin'
        call simple_getcwd(cwd)
        expected = filepath(cwd, requested)
        call proj%get_sigma2_state_path(registered, found)
        call assert_false(found, 'canonical sigma2 path initially absent')
        call proj%set_sigma2_state_path(requested)
        call proj%get_sigma2_state_path(registered, found)
        call assert_true(found, 'canonical sigma2 path registered')
        call assert_string_eq(expected%to_char(), registered, 'canonical sigma2 path retained')
        call proj%kill
    end subroutine test_sigma2_state_path_registration

    subroutine create_empty_file( fname )
        class(string), intent(in) :: fname
        integer :: funit, io_stat
        open(newunit=funit, file=fname%to_char(), status='replace', action='write', iostat=io_stat)
        call assert_int(0, io_stat, 'create remap target file')
        if( io_stat == 0 ) close(funit)
    end subroutine create_empty_file

    subroutine make_project(proj, projfile, nptcls, nptcls_stk, ncls, kv, cs, fraca, states, classes, &
        indstks, ogid)
        type(sp_project), intent(inout) :: proj
        type(string),     intent(in)    :: projfile
        integer,          intent(in)    :: nptcls, nptcls_stk, ncls, states(:), classes(:), indstks(:), ogid
        real,             intent(in)    :: kv, cs, fraca
        integer :: i
        call proj%kill
        call proj%os_stk%new(1, is_ptcl=.false.)
        call proj%os_stk%set(1, 'stk',        projfile%to_char()//'.mrcs')
        call proj%os_stk%set(1, 'ctf',        'yes')
        call proj%os_stk%set(1, 'smpd',       1.25)
        call proj%os_stk%set(1, 'kv',         kv)
        call proj%os_stk%set(1, 'cs',         cs)
        call proj%os_stk%set(1, 'fraca',      fraca)
        call proj%os_stk%set(1, 'phshift',    0.)
        call proj%os_stk%set(1, 'box',        128)
        call proj%os_stk%set(1, 'nptcls',     nptcls)
        if( nptcls_stk > 0 ) call proj%os_stk%set(1, 'nptcls_stk', nptcls_stk)
        call proj%os_stk%set(1, 'fromp',      1)
        call proj%os_stk%set(1, 'top',        nptcls)
        call proj%os_stk%set_ogid(1, ogid)
        call proj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        call proj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        do i = 1,nptcls
            call fill_particle_row(proj%os_ptcl2D, i, states(i), classes(i), indstks(i), ogid)
            call fill_particle_row(proj%os_ptcl3D, i, states(i), 0,          indstks(i), ogid)
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

    subroutine fill_particle_row(os, i, state, cls, indstk, ogid)
        use simple_oris, only: oris
        type(oris), intent(inout) :: os
        integer,    intent(in)    :: i, state, cls, indstk, ogid
        call os%set_stkind(i, 1)
        call os%set(i, 'indstk', indstk)
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
