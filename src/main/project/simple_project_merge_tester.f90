!@descr: unit tests for SIMPLE project merging
module simple_project_merge_tester
use simple_projfile_utils, only: merge_selected_project_files, validate_and_repair_project_file, remap_project_paths
use simple_sp_project,     only: sp_project
use simple_string,         only: string
use simple_syslib,         only: del_file, file_exists, simple_mkdir, simple_rmdir
use simple_test_utils
implicit none

private
public :: run_all_project_merge_tests

contains

    subroutine run_all_project_merge_tests()
        write(*,'(A)') '**** running all project merge tests ****'
        call test_merge_pruned_stack_indexing_heterogeneous_ctf()
        call test_validate_projfile_normalizes_legacy_stack_indices()
        call test_remap_project_paths()
        call test_remap_project_paths_allows_unmatched_scope()
        call test_remap_project_paths_scoped_roots()
    end subroutine run_all_project_merge_tests

    subroutine test_merge_pruned_stack_indexing_heterogeneous_ctf()
        type(sp_project) :: proj1, proj2, merged, reread
        type(string), allocatable :: project_files(:)
        type(string) :: projfile1, projfile2, merged_file
        integer, parameter :: NPTCLS1 = 3, NPTCLS2 = 2, NCLS = 2
        integer, parameter :: NPTCLS_STK1 = 6, NPTCLS_STK2 = 0
        integer, parameter :: INDSTKS1(NPTCLS1) = [1, 4, 6]
        integer, parameter :: INDSTKS2(NPTCLS2) = [99, 2]
        integer, parameter :: MAPPED_INDSTKS2(NPTCLS2) = [1, 2]
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
        call assert_false(proj2%os_stk%isthere(1, 'nptcls_stk'), 'project 2 source lacks nptcls_stk')
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
            call assert_int(MAPPED_INDSTKS2(i), reread%os_ptcl2D%get_int(NPTCLS1 + i, 'indstk'), &
                &'project 2 legacy indstk normalized')
            call reread%map_ptcl_ind2stk_ind('ptcl2D', NPTCLS1 + i, stkind, ind_in_stk)
            call assert_int(2, stkind, 'project 2 mapped stkind')
            call assert_int(MAPPED_INDSTKS2(i), ind_in_stk, 'project 2 mapped physical indstk')
            call reread%map_ptcl_ind2stk_ind('ptcl3D', NPTCLS1 + i, stkind, ind_in_stk)
            call assert_int(2, stkind, 'project 2 ptcl3D mapped stkind')
            call assert_int(MAPPED_INDSTKS2(i), ind_in_stk, 'project 2 ptcl3D mapped physical indstk')
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

    subroutine test_validate_projfile_normalizes_legacy_stack_indices()
        type(sp_project) :: proj, validated
        type(string) :: projfile, validated_file
        integer, parameter :: NPTCLS = 3, NCLS = 2
        integer, parameter :: BAD_INDSTKS(NPTCLS) = [99, 0, 3]
        integer :: i, stkind, ind_in_stk
        write(*,'(A)') 'test_validate_projfile_normalizes_legacy_stack_indices'
        projfile       = 'validate_project_src.simple'
        validated_file = 'validate_project_src_validated.simple'
        call del_file(projfile)
        call del_file(validated_file)
        call make_project(proj, projfile, NPTCLS, 0, NCLS, 200.0, 1.0, 0.07, &
            [1, 0, 1], [1, 2, 1], BAD_INDSTKS, 1)
        call proj%os_ptcl2D%delete_entry(2, 'indstk')
        call proj%os_ptcl3D%delete_entry(2, 'indstk')
        call proj%write(projfile)
        call assert_false(proj%os_stk%isthere(1, 'nptcls_stk'), 'legacy source lacks nptcls_stk')
        call validate_and_repair_project_file(projfile, validated_file)
        call assert_true(file_exists(validated_file), 'validate_projfile creates validated project')
        call validated%read(validated_file)
        call assert_int(1, validated%os_stk%get_noris(), 'validated stack count')
        call assert_int(NPTCLS, validated%os_ptcl2D%get_noris(), 'validated ptcl2D count includes state 0')
        call assert_int(NPTCLS, validated%os_ptcl3D%get_noris(), 'validated ptcl3D count includes state 0')
        call assert_int(1, validated%os_stk%get_fromp(1), 'validated fromp')
        call assert_int(NPTCLS, validated%os_stk%get_top(1), 'validated top')
        call assert_int(NPTCLS, validated%os_stk%get_int(1, 'nptcls'), 'validated project nptcls')
        call assert_int(NPTCLS, validated%os_stk%get_int(1, 'nptcls_stk'), 'validated fallback nptcls_stk')
        do i = 1,NPTCLS
            call assert_int(1, validated%os_ptcl2D%get_int(i, 'stkind'), 'validated ptcl2D stkind')
            call assert_int(i, validated%os_ptcl2D%get_int(i, 'indstk'), 'validated ptcl2D legacy indstk')
            call validated%map_ptcl_ind2stk_ind('ptcl2D', i, stkind, ind_in_stk)
            call assert_int(1, stkind, 'validated ptcl2D mapped stkind')
            call assert_int(i, ind_in_stk, 'validated ptcl2D mapped indstk')
            call assert_int(1, validated%os_ptcl3D%get_int(i, 'stkind'), 'validated ptcl3D stkind')
            call assert_int(i, validated%os_ptcl3D%get_int(i, 'indstk'), 'validated ptcl3D legacy indstk')
            call validated%map_ptcl_ind2stk_ind('ptcl3D', i, stkind, ind_in_stk)
            call assert_int(1, stkind, 'validated ptcl3D mapped stkind')
            call assert_int(i, ind_in_stk, 'validated ptcl3D mapped indstk')
        enddo
        call assert_int(0, validated%os_ptcl2D%get_state(2), 'state 0 row remains present after validation')
        call del_file(projfile)
        call del_file(validated_file)
        call proj%kill
        call validated%kill
    end subroutine test_validate_projfile_normalizes_legacy_stack_indices

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
        type(string) :: stk_file, stk_den_file, stk_box_file
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
        stk_den_file = ptcl_root//'den.mrcs'
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
        call create_empty_file(stk_den_file)
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
        call proj%os_stk%set(1, 'stk_den', '/legacy/ptcl/den.mrcs')
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
        call assert_int(3, nremapped, 'scoped particle path count')
        call remap_project_paths(proj, string('/legacy/cavg'), cavg_root, nremapped, scope='cavg')
        call assert_int(4, nremapped, 'scoped class-average path count')
        call remap_project_paths(proj, string('/legacy/vol'), vol_root, nremapped, scope='vol')
        call assert_int(3, nremapped, 'scoped volume path count')
        call assert_string_eq(movie_file%to_char(), proj%os_mic%get_str(1, 'movie'), 'scoped movie path')
        call assert_string_eq(intg_file%to_char(), proj%os_mic%get_str(1, 'intg'), 'scoped intg path')
        call assert_string_eq(mic_box_file%to_char(), proj%os_mic%get_str(1, 'boxfile'), 'scoped mic box path')
        call assert_string_eq(stk_file%to_char(), proj%os_stk%get_str(1, 'stk'), 'scoped raw stack path')
        call assert_string_eq(stk_den_file%to_char(), proj%os_stk%get_str(1, 'stk_den'), 'scoped denoised stack path')
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
