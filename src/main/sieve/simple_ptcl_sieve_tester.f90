!@descr: unit test routines for the ptcl_sieve orchestrator
module simple_ptcl_sieve_tester
use simple_test_utils
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_ptcl_sieve, only: ptcl_sieve
use simple_sp_project, only: sp_project
use simple_image,      only: image
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
        call test_single_pass_ignores_incomplete_fine()
        call test_new_accepts_tuning_overrides()
        call test_cycle_empty_project_list()
        call test_collect_and_reject_hard_gates()
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
        call assert_int(0, sieve%get_n_total_particles(),'new() initializes zero total particles')
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

    subroutine test_single_pass_ignores_incomplete_fine()
        type(ptcl_sieve) :: sieve
        type(parameters) :: params
        type(string)     :: ws_dir, cwd_saved

        write(*,'(A)') 'test_single_pass_ignores_incomplete_fine'

        call setup_workspace(string('single_pass'), ws_dir, cwd_saved)

        ! Coarse chunk is complete.
        call make_chunk_project('coarse', 1, 10, 6, 2)
        call simple_touch(string('chunks_coarse/chunk_coarse_1/' // ABINITIO2D_FINISHED))
        call simple_touch(string('chunks_coarse/chunk_coarse_1/REJECTION_FINISHED'))
        call simple_touch(string('chunks_coarse/chunk_coarse_1/COMPLETE'))

        ! Fine chunk exists but is incomplete.
        call make_chunk_project('fine', 1, 9, 4, 1)
        call simple_touch(string('chunks_fine/chunk_fine_1/' // ABINITIO2D_FINISHED))
        call simple_touch(string('chunks_fine/chunk_fine_1/REJECTION_FINISHED'))

        ! Baseline (two-tier): incomplete fine chunk prevents finished state.
        call init_test_params(params)
        call sieve%new(params, string('completed'))
        call assert_false(sieve%get_finished(), 'two-tier run is not finished when fine chunk is incomplete')
        call sieve%kill()

        ! single_pass=yes: coarse-only terminal semantics apply.
        call init_test_params(params, single_pass='yes')
        call sieve%new(params, string('completed'))
        call assert_true(sieve%get_finished(), 'single_pass run ignores incomplete fine tier for finished state')
        call sieve%kill()

        call teardown_workspace(ws_dir, cwd_saved)
    end subroutine test_single_pass_ignores_incomplete_fine

    subroutine test_new_accepts_tuning_overrides()
        type(ptcl_sieve) :: sieve
        type(parameters) :: params
        type(string)     :: ws_dir, cwd_saved

        write(*,'(A)') 'test_new_accepts_tuning_overrides'

        call setup_workspace(string('override_init'), ws_dir, cwd_saved)
        call init_test_params(params, lpstart=12.0, lpstop_coarse=18.0, lpstop_fine=9.0, &
                              box_coarse=96, box_fine=80, nsample_coarse=500, nsample_fine=250, &
                              ncls_coarse=64, ncls_fine=48)
        call sieve%new(params, string('completed'))

        call assert_int(0, sieve%get_n_chunks_coarse(), 'override init keeps empty coarse chunk list')
        call assert_int(0, sieve%get_n_chunks_fine(),   'override init keeps empty fine chunk list')

        call sieve%kill()
        call teardown_workspace(ws_dir, cwd_saved)
    end subroutine test_new_accepts_tuning_overrides

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

    !> collect_and_reject on one completed coarse chunk, single pass, without the learned model:
    !! a class average with a strong centred component passes the hard gates, a blank one is
    !! rejected, and the selection reaches the particles, the sentinels, the exported project,
    !! the previews and the latest-product metadata (Ruben's stream test sieve_cavgs, moved here
    !! by the stream review, plan section 9.7)
    subroutine test_collect_and_reject_hard_gates()
        character(len=*), parameter :: CHUNK_STEM = 'chunk_coarse_1'
        character(len=*), parameter :: CAVG_STACK = 'cavgs_iter001'//MRC_EXT
        real,             parameter :: SMPD       = 2.0
        integer,          parameter :: CAVG_BOX   = 64
        integer,          parameter :: NCLASSES   = 2
        integer,          parameter :: NPARTICLES = 8
        integer,          parameter :: NPER_CLASS = NPARTICLES / NCLASSES
        type(ptcl_sieve)     :: sieve
        type(parameters)     :: params_sieve
        type(cmdline)        :: cline_sieve
        type(image)          :: cavg_good, cavg_bad, feature
        type(sp_project)     :: chunk_project, result
        type(string)         :: ws_dir, cwd_saved, chunk_dir, completed_path, chunk_projfile, cavg_path
        type(string)         :: completed_projfile, rejection_reason, latest_jpeg, latest_stk
        type(string)         :: selected_jpeg, rejected_jpeg, reasons_jpeg, reasons_key
        integer, allocatable :: latest_inds(:), latest_pops(:), latest_selection(:)
        real,    allocatable :: latest_res(:)
        integer              :: i, icls, ldim(3), nimages, xtiles, ytiles
        logical              :: has_latest

        write(*,'(A)') 'test_collect_and_reject_hard_gates'

        call setup_workspace(string('collect_reject'), ws_dir, cwd_saved)
        completed_path = filepath(ws_dir, 'completed')
        chunk_dir      = filepath(ws_dir, 'chunks_coarse')
        call simple_mkdir(chunk_dir)
        chunk_dir      = filepath(chunk_dir, CHUNK_STEM)
        call simple_mkdir(chunk_dir)
        chunk_projfile = filepath(chunk_dir, CHUNK_STEM//METADATA_EXT)
        cavg_path      = filepath(chunk_dir, CAVG_STACK)

        ! class 1: a strong centred square on weak noise; class 2: blank (no component)
        call cavg_good%new([CAVG_BOX, CAVG_BOX, 1], SMPD, wthreads=.false.)
        call cavg_good%gauran(0., 0.02)
        call feature%new([CAVG_BOX, CAVG_BOX, 1], SMPD, wthreads=.false.)
        call feature%square(8)
        call feature%mul(5.)
        call cavg_good%add(feature)
        call cavg_good%write(cavg_path, 1, del_if_exists=.true.)
        call cavg_bad%new([CAVG_BOX, CAVG_BOX, 1], SMPD, wthreads=.false.)
        call cavg_bad%zero()
        call cavg_bad%write(cavg_path, 2)
        call feature%kill()
        call cavg_good%kill()
        call cavg_bad%kill()

        ! a completed coarse abinitio2D chunk awaiting rejection
        call chunk_project%os_mic%new(1, is_ptcl=.false.)
        call chunk_project%os_mic%set_state(1, 1)
        call chunk_project%os_mic%set(1, 'imgkind', 'mic')
        call chunk_project%os_mic%set(1, 'nptcls',  NPARTICLES)
        call chunk_project%os_mic%set(1, 'smpd',    SMPD)
        call chunk_project%os_ptcl2D%new(NPARTICLES, is_ptcl=.true.)
        do i = 1, NPARTICLES
            icls = 1 + (i - 1) / NPER_CLASS
            call chunk_project%os_ptcl2D%set_class(i, icls)
            call chunk_project%os_ptcl2D%set_state(i, 1)
            call chunk_project%os_ptcl2D%set_stkind(i, 1)
            call chunk_project%os_ptcl2D%set(i, 'indstk', i)
        enddo
        chunk_project%os_ptcl3D = chunk_project%os_ptcl2D
        call chunk_project%add_cavgs2os_out(cavg_path, SMPD, imgkind='cavg')
        do icls = 1, NCLASSES
            call chunk_project%os_cls2D%set_class(icls, icls)
            call chunk_project%os_cls2D%set_state(icls, 1)
            call chunk_project%os_cls2D%set(icls, 'pop',  NPER_CLASS)
            call chunk_project%os_cls2D%set(icls, 'res',  10.0)
            call chunk_project%os_cls2D%set(icls, 'corr', 0.9)
        enddo
        chunk_project%os_cls3D = chunk_project%os_cls2D
        call chunk_project%update_projinfo(chunk_projfile)

        ! the production collector in coarse-only mode; the learned model is off, so only the
        ! deterministic hard gates decide
        call cline_sieve%set('prg',            'sieve_cavgs')
        call cline_sieve%set('projfile',       chunk_projfile)
        call cline_sieve%set('dir_target',     ws_dir)
        call cline_sieve%set('ncls',           NCLASSES)
        call cline_sieve%set('nptcls_per_cls', NPER_CLASS)
        call cline_sieve%set('nchunksperset',  1)
        call cline_sieve%set('nchunks',        1)
        call cline_sieve%set('nparts',         1)
        call cline_sieve%set('nthr',           1)
        call cline_sieve%set('nptcls_coarse',  NPARTICLES)
        call cline_sieve%set('ncls_coarse',    NCLASSES)
        call cline_sieve%set('box_coarse',     CAVG_BOX)
        call cline_sieve%set('nsample_coarse', NPARTICLES)
        call cline_sieve%set('lpstart',        20.0)
        call cline_sieve%set('lpstop_coarse',  15.0)
        call cline_sieve%set('mskdiam',        80.0)
        call cline_sieve%set('single_pass',    'yes')
        call cline_sieve%set('use_model',      'no')
        call cline_sieve%set('qsys_name',      'local')
        call cline_sieve%set('walltime',       60)
        call cline_sieve%set('mkdir',          'no')
        call chunk_project%update_compenv(cline_sieve)
        call chunk_project%write(chunk_projfile)
        call chunk_project%kill()
        call simple_touch(filepath(chunk_dir, ABINITIO2D_FINISHED))
        call params_sieve%new(cline_sieve)
        call sieve%new(params_sieve, completed_path)
        call sieve%collect_and_reject()

        ! the collector's counters
        call assert_int(1,          sieve%get_n_chunks_coarse(),         'the completed coarse chunk is imported')
        call assert_int(NPER_CLASS, sieve%get_n_coarse_accepted_ptcls(), 'coarse accepted particles = the class that passes')
        call assert_int(NPER_CLASS, sieve%get_n_coarse_rejected_ptcls(), 'coarse rejected particles = the blank class')
        call assert_int(NPER_CLASS, sieve%get_n_accepted_ptcls(),        'final accepted particles')
        call assert_int(NPER_CLASS, sieve%get_n_rejected_ptcls(),        'final rejected particles')
        call assert_int(NPARTICLES, sieve%get_n_total_particles(),       'total particles')
        call assert_int(1,          sieve%get_n_accepted_micrographs(),  'accepted micrographs')
        call assert_true(sieve%get_finished(),                           'the coarse-only sieve finishes')
        call assert_true(file_exists(filepath(chunk_dir, 'REJECTION_FINISHED')), 'the rejection sentinel exists')
        call assert_true(file_exists(filepath(chunk_dir, 'COMPLETE')),           'the completion sentinel exists')

        ! the selection in the chunk project, mapped to the particles
        call result%read(chunk_projfile)
        call assert_true(all(result%os_cls2D%get_all_asint('state') == [1, 0]), 'class 1 selected, class 2 rejected')
        call assert_int(NPER_CLASS, result%os_ptcl2D%count_state_gt_zero(), 'the selection reaches the 2D particles')
        call assert_int(NPER_CLASS, result%os_ptcl3D%count_state_gt_zero(), 'the selection reaches the 3D particles')
        ! the orientation reader splits character values at blanks, so the round-tripped part of
        ! the reason is its tier prefix
        rejection_reason = result%os_cls2D%get_str(2, 'rejection_reason')
        call assert_true(rejection_reason%has_substr('coarse_reject'), 'the rejected class records the coarse rejection')
        call result%kill()

        ! the exported project and the previews
        completed_projfile = filepath(completed_path, CHUNK_STEM//METADATA_EXT)
        selected_jpeg      = filepath(chunk_dir, CHUNK_STEM//'_selected'//JPG_EXT)
        rejected_jpeg      = filepath(chunk_dir, CHUNK_STEM//'_rejected'//JPG_EXT)
        reasons_jpeg       = filepath(chunk_dir, CHUNK_STEM//'_all_reasons'//JPG_EXT)
        reasons_key        = reasons_jpeg//'.key.txt'
        call assert_true(file_exists(completed_projfile), 'the chunk project is exported to the completed directory')
        call assert_true(file_exists(selected_jpeg),      'the selected-classes preview exists')
        call assert_true(file_exists(rejected_jpeg),      'the rejected-classes preview exists')
        call assert_true(file_exists(reasons_jpeg),       'the rejection-reason report exists')
        call assert_true(file_exists(reasons_key),        'the rejection-reason key exists')

        ! the latest product the GUI shows
        has_latest = sieve%get_latest(latest_inds, latest_pops, latest_res, latest_jpeg, latest_stk, &
            &xtiles, ytiles, latest_selection)
        call assert_true(has_latest, 'the latest class-average product is available')
        if( has_latest )then
            call assert_int(NCLASSES, size(latest_inds), 'latest product: two classes')
            if( size(latest_inds) == NCLASSES )then
                call assert_true(all(latest_inds      == [1, 2]),                 'latest product: class indices')
                call assert_true(all(latest_pops      == [NPER_CLASS, NPER_CLASS]), 'latest product: populations')
                call assert_true(all(latest_selection == [1, 0]),                 'latest product: selection')
                call assert_true(all(abs(latest_res - 10.0) <= 0.01),             'latest product: resolutions')
            endif
            call assert_true(xtiles * ytiles >= NCLASSES, 'latest product: the preview tiles hold every class')
            call assert_true(file_exists(latest_jpeg), 'latest product: the preview exists')
            call assert_true(file_exists(latest_stk),  'latest product: the retained stack exists')
            if( file_exists(latest_stk) )then
                call find_ldim_nptcls(latest_stk, ldim, nimages)
                call assert_int(NCLASSES, nimages, 'latest product: the retained stack holds every class')
                call assert_true(all(ldim(1:2) == [CAVG_BOX, CAVG_BOX]), 'latest product: the retained stack keeps the box')
            endif
        endif

        if( allocated(latest_inds)      ) deallocate(latest_inds)
        if( allocated(latest_pops)      ) deallocate(latest_pops)
        if( allocated(latest_res)       ) deallocate(latest_res)
        if( allocated(latest_selection) ) deallocate(latest_selection)
        call sieve%kill()
        call cline_sieve%kill()
        call teardown_workspace(ws_dir, cwd_saved)
    end subroutine test_collect_and_reject_hard_gates

    subroutine init_test_params(params, single_pass, lpstart, lpstop_coarse, lpstop_fine, box_coarse, box_fine, &
                                nsample_coarse, nsample_fine, ncls_coarse, ncls_fine)
        type(parameters), intent(inout) :: params
        type(cmdline)                   :: cline
        character(len=*), optional, intent(in) :: single_pass
        real, optional, intent(in) :: lpstart, lpstop_coarse, lpstop_fine
        integer, optional, intent(in) :: box_coarse, box_fine, nsample_coarse, nsample_fine, ncls_coarse, ncls_fine

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
        if( present(single_pass)   ) call cline%set('single_pass',    trim(single_pass))
        if( present(lpstart)       ) call cline%set('lpstart',        lpstart)
        if( present(lpstop_coarse) ) call cline%set('lpstop_coarse',  lpstop_coarse)
        if( present(lpstop_fine)   ) call cline%set('lpstop_fine',    lpstop_fine)
        if( present(box_coarse)    ) call cline%set('box_coarse',     box_coarse)
        if( present(box_fine)      ) call cline%set('box_fine',       box_fine)
        if( present(nsample_coarse)) call cline%set('nsample_coarse', nsample_coarse)
        if( present(nsample_fine)  ) call cline%set('nsample_fine',   nsample_fine)
        if( present(ncls_coarse)   ) call cline%set('ncls_coarse',    ncls_coarse)
        if( present(ncls_fine)     ) call cline%set('ncls_fine',      ncls_fine)
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
