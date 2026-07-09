!@descr: unit test routines for the ptcl_sieve orchestrator
module simple_ptcl_sieve_tester
use simple_test_utils
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_ptcl_sieve, only: ptcl_sieve
use simple_sp_project, only: sp_project
use simple_rec_list,   only: rec_list
use simple_string,     only: string
use simple_defs_fname, only: METADATA_EXT, ABINITIO2D_FINISHED
implicit none
private
public :: run_all_ptcl_sieve_tests

contains

    subroutine run_all_ptcl_sieve_tests()
        write(*,'(A)') '**** running all ptcl_sieve tests ****'
        call test_new_kill_and_empty_queries()
        call test_import_existing_chunks_and_counts()
        call test_finished_semantics()
        call test_cycle_empty_project_list()
    end subroutine run_all_ptcl_sieve_tests

    subroutine test_new_kill_and_empty_queries()
        type(ptcl_sieve)          :: sieve
        type(parameters)          :: params
        type(string)              :: ws_dir, cwd_saved, jpeg, stk
        integer, allocatable      :: inds(:), pops(:), sel(:)
        real, allocatable         :: res(:)
        integer                   :: xtiles, ytiles
        logical                   :: ok

        write(*,'(A)') 'test_new_kill_and_empty_queries'

        call setup_workspace(string('new_kill'), ws_dir, cwd_saved)
        call init_test_params(params)
        call sieve%new(params, string('completed'))

        call assert_int(0, sieve%get_n_chunks_coarse(), 'new() initializes zero coarse chunks')
        call assert_int(0, sieve%get_n_chunks_fine(),   'new() initializes zero fine chunks')
        call assert_int(0, sieve%get_n_chunks_running(),'new() initializes zero running chunks')
        call assert_false(sieve%get_finished(),         'fresh sieve is not finished')

        ok = sieve%get_latest(inds, pops, res, jpeg, stk, xtiles, ytiles, sel)
        call assert_false(ok, 'get_latest() returns false when no latest product exists')

        call sieve%set_final_ingestion()
        call sieve%kill()
        call sieve%kill() ! idempotence check

        call teardown_workspace(ws_dir, cwd_saved)
    end subroutine test_new_kill_and_empty_queries

    subroutine test_import_existing_chunks_and_counts()
        type(ptcl_sieve) :: sieve
        type(parameters) :: params
        type(string)     :: ws_dir, cwd_saved

        write(*,'(A)') 'test_import_existing_chunks_and_counts'

        call setup_workspace(string('import_counts'), ws_dir, cwd_saved)

        ! Coarse #1 contributes to pass-1 non-rejected count.
        call make_chunk_project('coarse', 1, 10, 6, 2)
        call simple_touch(string('chunks_coarse/chunk_coarse_1/' // ABINITIO2D_FINISHED))
        call simple_touch(string('chunks_coarse/chunk_coarse_1/REJECTION_FINISHED'))

        ! Coarse #2 is already complete, so excluded from pass-1 non-rejected count.
        call make_chunk_project('coarse', 2, 8, 3, 1)
        call simple_touch(string('chunks_coarse/chunk_coarse_2/' // ABINITIO2D_FINISHED))
        call simple_touch(string('chunks_coarse/chunk_coarse_2/REJECTION_FINISHED'))
        call simple_touch(string('chunks_coarse/chunk_coarse_2/COMPLETE'))

        ! Coarse #3 is failed and excluded from non-rejected counts.
        call make_chunk_project('coarse', 3, 7, 4, 1)
        call simple_touch(string('chunks_coarse/chunk_coarse_3/REJECTION_FAILED'))

        ! Fine #1 contributes to pass-2 non-rejected count.
        call make_chunk_project('fine', 1, 9, 5, 2)
        call simple_touch(string('chunks_fine/chunk_fine_1/' // ABINITIO2D_FINISHED))
        call simple_touch(string('chunks_fine/chunk_fine_1/REJECTION_FINISHED'))

        call init_test_params(params)
        call sieve%new(params, string('completed'))

        call assert_int(3, sieve%get_n_chunks_coarse(), 'imported three coarse chunks')
        call assert_int(1, sieve%get_n_chunks_fine(),   'imported one fine chunk')
        call assert_int(0, sieve%get_n_chunks_running(),'imported chunks start as non-running')
        call assert_int(6, sieve%get_n_pass_1_non_rejected_ptcls(), 'pass-1 selected count excludes complete/failed coarse chunks')
        call assert_int(5, sieve%get_n_pass_2_non_rejected_ptcls(), 'pass-2 selected count includes rejection-complete fine chunks')
        call assert_false(sieve%get_finished(), 'not finished while coarse/fine chunks remain incomplete')

        call sieve%kill()
        call teardown_workspace(ws_dir, cwd_saved)
    end subroutine test_import_existing_chunks_and_counts

    subroutine test_finished_semantics()
        type(ptcl_sieve) :: sieve
        type(parameters) :: params
        type(string)     :: ws_dir, cwd_saved

        write(*,'(A)') 'test_finished_semantics'

        call setup_workspace(string('finished'), ws_dir, cwd_saved)

        ! Coarse complete case.
        call make_chunk_project('coarse', 1, 10, 6, 2)
        call simple_touch(string('chunks_coarse/chunk_coarse_1/' // ABINITIO2D_FINISHED))
        call simple_touch(string('chunks_coarse/chunk_coarse_1/REJECTION_FINISHED'))
        call simple_touch(string('chunks_coarse/chunk_coarse_1/COMPLETE'))

        ! Failed coarse chunk still satisfies completion criteria.
        call make_chunk_project('coarse', 2, 6, 0, 1)
        call simple_touch(string('chunks_coarse/chunk_coarse_2/REJECTION_FAILED'))

        call init_test_params(params)

        ! coarse_only: all coarse chunks complete/failed => finished.
        call sieve%new(params, string('completed'), coarse_only=.true.)
        call assert_true(sieve%get_finished(), 'coarse_only run finishes when all coarse chunks are complete/failed')
        call sieve%kill()

        ! two-tier mode: no fine chunks => coarse completion is terminal.
        call sieve%new(params, string('completed'))
        call assert_true(sieve%get_finished(), 'two-tier run with no fine chunks finishes at coarse completion')
        call sieve%kill()

        ! Add incomplete fine chunk; now two-tier mode is not finished.
        call make_chunk_project('fine', 1, 9, 4, 1)
        call simple_touch(string('chunks_fine/chunk_fine_1/' // ABINITIO2D_FINISHED))
        call simple_touch(string('chunks_fine/chunk_fine_1/REJECTION_FINISHED'))

        call sieve%new(params, string('completed'))
        call assert_false(sieve%get_finished(), 'two-tier run not finished while any fine chunk is incomplete')
        call sieve%kill()

        ! Mark fine complete; now finished.
        call simple_touch(string('chunks_fine/chunk_fine_1/COMPLETE'))
        call sieve%new(params, string('completed'))
        call assert_true(sieve%get_finished(), 'two-tier run finishes once all fine chunks are complete/failed')
        call sieve%kill()

        call teardown_workspace(ws_dir, cwd_saved)
    end subroutine test_finished_semantics

    subroutine test_cycle_empty_project_list()
        type(ptcl_sieve) :: sieve
        type(parameters) :: params
        type(rec_list)   :: project_list
        type(string)     :: ws_dir, cwd_saved

        write(*,'(A)') 'test_cycle_empty_project_list'

        call setup_workspace(string('cycle_empty'), ws_dir, cwd_saved)
        call init_test_params(params)
        call sieve%new(params, string('completed'))

        call sieve%cycle(project_list)
        call assert_int(0, sieve%get_n_chunks_coarse(), 'cycle on empty list creates no coarse chunks')
        call assert_int(0, sieve%get_n_chunks_fine(),   'cycle on empty list creates no fine chunks')

        call sieve%kill()
        call teardown_workspace(ws_dir, cwd_saved)
    end subroutine test_cycle_empty_project_list

    subroutine init_test_params(params)
        type(parameters), intent(inout) :: params
        type(cmdline)                   :: cline

        call cline%set('prg',        'abinitio2D')
        call cline%set('mkdir',      'yes')
        call cline%set('split_mode', 'even')
        call cline%set('nchunks',    1)
        call cline%set('nthr',       1)
        call cline%set('nparts',     1)
        call cline%set('nptcls',     16)
        call cline%set('mskdiam',    120.0)
        call cline%set('walltime',   60)
        call cline%set('qsys_name',  'local')
        call params%new(cline)
    end subroutine init_test_params

    subroutine setup_workspace(tag, ws_dir, cwd_saved)
        type(string), intent(in)    :: tag
        type(string), intent(inout) :: ws_dir, cwd_saved

        call simple_getcwd(cwd_saved)
        ws_dir = filepath(cwd_saved, string('PTCL_SIEVE_TEST_' // tag%to_char() // '_' // int2str(get_process_id())))
        call exec_cmdline('rm -rf ' // ws_dir%to_char())
        call simple_mkdir(ws_dir)
        call simple_chdir(ws_dir)
        call simple_mkdir('completed')
    end subroutine setup_workspace

    subroutine teardown_workspace(ws_dir, cwd_saved)
        type(string), intent(in) :: ws_dir, cwd_saved

        call simple_chdir(cwd_saved)
        call exec_cmdline('rm -rf ' // ws_dir%to_char())
    end subroutine teardown_workspace

    subroutine make_chunk_project(tier, id, nptcls, nsel, nmics)
        character(len=*), intent(in) :: tier
        integer,          intent(in) :: id, nptcls, nsel, nmics
        type(sp_project)             :: proj
        type(string)                 :: chunk_dir, projfile, stem
        integer                      :: i

        if( tier == 'coarse' ) then
            stem = string('chunk_coarse_' // int2str(id))
            chunk_dir = string('chunks_coarse/' // stem%to_char())
        else
            stem = string('chunk_fine_' // int2str(id))
            chunk_dir = string('chunks_fine/' // stem%to_char())
        end if

        call simple_mkdir(string('chunks_' // tier))
        call simple_mkdir(chunk_dir)
        projfile = string(chunk_dir%to_char() // '/' // stem%to_char() // METADATA_EXT)

        call proj%kill()

        call proj%os_mic%new(max(1, nmics), is_ptcl=.false.)
        do i = 1, max(1, nmics)
            call proj%os_mic%set_state(i, 1)
            call proj%os_mic%set(i, 'nptcls', max(1, nptcls / max(1, nmics)))
        end do

        call proj%os_stk%new(1, is_ptcl=.false.)
        call proj%os_stk%set_state(1, 1)
        call proj%os_stk%set(1, 'nptcls', nptcls)
        call proj%os_stk%set(1, 'fromp',  1)
        call proj%os_stk%set(1, 'top',    nptcls)
        call proj%os_stk%set(1, 'stk',    stem%to_char() // '.mrcs')

        call proj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        do i = 1, nptcls
            call proj%os_ptcl2D%set_class(i, 1)
            call proj%os_ptcl2D%set_stkind(i, 1)
            if( i <= nsel ) then
                call proj%os_ptcl2D%set_state(i, 1)
            else
                call proj%os_ptcl2D%set_state(i, 0)
            end if
        end do
        proj%os_ptcl3D = proj%os_ptcl2D

        call proj%os_cls2D%new(1, is_ptcl=.false.)
        call proj%os_cls2D%set_class(1, 1)
        call proj%os_cls2D%set_state(1, merge(1, 0, nsel > 0))
        call proj%os_cls2D%set(1, 'pop', nsel)
        call proj%os_cls2D%set(1, 'res', 10.0)

        call proj%update_projinfo(projfile)
        call proj%write(projfile)
        call proj%kill()
    end subroutine make_chunk_project

end module simple_ptcl_sieve_tester
