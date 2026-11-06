!@descr: unit tests of solve3D stage rules, symmetry randomization, and solve3D_cavgs state reseeding
! In-memory labels and scores of the even (1..ncavgs) and odd (ncavgs+1..2*ncavgs) halves of the
! class averages; no volume or job. The symmetry randomization test writes one small project file.
module simple_solve3D_reseed_tester
use simple_solve3D_utils, only: reseed_state_labels, solve3D_symsrch_stage, solve3D_stage_skipping, &
    &symmetry_randomize_poses
use simple_test_utils
implicit none
private
public :: run_all_solve3D_reseed_tests

contains

    subroutine run_all_solve3D_reseed_tests()
        write(*,'(A)') '**** running all solve3D_cavgs stage rule tests ****'
        call test_symsrch_stage()
        call test_stage_skipping()
        call test_gold_standard_stages()
        call test_symmetry_randomize()
        call test_no_weak_state()
        call test_one_empty_state()
        call test_all_in_one_state()
        call test_weak_by_fraction()
        call test_too_few_to_reseed()
    end subroutine run_all_solve3D_reseed_tests

    !> the symmetry axis search: after stage 4 in a run that searches an axis, stage 3 otherwise
    subroutine test_symsrch_stage()
        write(*,'(A)') 'test_symsrch_stage'
        call assert_int(3, solve3D_symsrch_stage(.false.), 'no axis search: stage 3')
        call assert_int(4, solve3D_symsrch_stage(.true.),  'axis search after stage 4')
    end subroutine test_symsrch_stage

    !> stage omission and the mid-stage resolution exit: single-state runs only
    subroutine test_stage_skipping()
        write(*,'(A)') 'test_stage_skipping'
        call assert_true(solve3D_stage_skipping(1),        'single state: stages may be skipped')
        call assert_false(solve3D_stage_skipping(2),       'two states: no stage skipped')
        call assert_false(solve3D_stage_skipping(4),       'four states: no stage skipped')
    end subroutine test_stage_skipping

    !> Exercise emitted stage commands, including a skipped transition and stale matching limits.
    subroutine test_gold_standard_stages()
        use simple_solve3D_utils, only: lpinfo, cline_refine3D, set_cline_refine3D, solve3D_addon_ctx, &
            &l_nonuniform, l_srch4symaxis, l_cavg_ini, l_state_continue_mode, &
            &l_refine3D_mode_override, l_refine3D_lp_override, l_refine3D_lpstop_override, &
            &l_refine3D_trs_override, l_prob_beta_end_override, update_frac, nptcls_eff, nstages_refine3D
        use simple_parameters, only: parameters
        use simple_cmdline,    only: cmdline
        use simple_type_defs,  only: lp_crop_inf
        use simple_sp_project, only: sp_project
        use simple_string,     only: string
        class(parameters), allocatable :: params
        class(cmdline), allocatable    :: saved_cline
        class(sp_project), allocatable :: proj
        type(lp_crop_inf), allocatable :: saved_lpinfo(:)
        type(solve3D_addon_ctx) :: addon
        type(string) :: cwd, root
        logical :: saved_flags(9)
        real    :: saved_update_frac
        integer :: saved_nptcls_eff, saved_nstages, nfail0, backend
        write(*,'(A)') 'test_gold_standard_stages'
        nfail0 = tests_failed
        call enter_fixture('solve3D_gold_standard', cwd, root)
        allocate(params, saved_cline, proj)
        saved_cline = cline_refine3D
        saved_flags = [l_nonuniform, l_srch4symaxis, l_cavg_ini, l_state_continue_mode, &
            &l_refine3D_mode_override, l_refine3D_lp_override, l_refine3D_lpstop_override, &
            &l_refine3D_trs_override, l_prob_beta_end_override]
        saved_update_frac = update_frac
        saved_nptcls_eff  = nptcls_eff
        saved_nstages     = nstages_refine3D
        call move_alloc(lpinfo, saved_lpinfo)
        allocate(lpinfo(8))
        lpinfo%lp = [20., 16., 12., 10., 8., 6., 5., 4.5]
        lpinfo%box_crop = 64
        lpinfo%trslim   = 2.
        l_nonuniform = .true.
        l_srch4symaxis = .false.
        l_cavg_ini = .false.
        l_state_continue_mode = .false.
        l_refine3D_mode_override = .false.
        l_refine3D_lp_override = .false.
        l_refine3D_lpstop_override = .false.
        l_refine3D_trs_override = .false.
        l_prob_beta_end_override = .false.
        update_frac = 0.25
        nptcls_eff = 10000
        nstages_refine3D = 8
        params%nsample = 1000
        params%nstates = 1
        params%box = 64
        params%smpd = 2.
        params%gold_standard = 'yes'
        params%automsk = 'yes'
        params%envfsc = 'yes'
        params%projfile = 'stages.simple'
        call proj%os_ptcl3D%new(2, is_ptcl=.true.)
        call proj%os_ptcl3D%set_all2single('state', 1.)
        call proj%os_ptcl3D%set_all2single('res05', 5.3)
        call proj%update_projinfo(params%projfile)
        call proj%write(params%projfile)
        do backend = 1, 2
            params%rec_backend = 'gridding'
            if( backend == 2 ) params%rec_backend = 'pcg'
            params%filt_mode = 'nonuniform_lpset'
            call set_cline_refine3D(params, 5, 1, .false.)
            call assert_string_eq('nonuniform_lpset', cline_refine3D%get_carg('filt_mode'), 'stage 5 uses merged NU references')
            call assert_real(6., cline_refine3D%get_rarg('lp'), 1.e-6, 'FSC=0.5 still advances the planned 8 A rung')
            call assert_false(cline_refine3D%defined('lpstop'), 'NU stage 5 has no added stage ceiling')
            ! Stage 6 is omitted: stage 7 must clear the preceding LP-set command.
            call set_cline_refine3D(params, 7, 1, .false.)
            call assert_string_eq('nonuniform', cline_refine3D%get_carg('filt_mode'), 'stage 7 uses half-specific NU references')
            call assert_false(cline_refine3D%defined('lp'), 'transition clears the preceding scheduled lp')
            call assert_false(cline_refine3D%defined('lpstop'), 'gold-standard stage has no scheduled ceiling')
            call assert_real(0.25, cline_refine3D%get_rarg('update_frac'), 1.e-6, 'the transition keeps fractional updates')
            call assert_string_eq('yes', cline_refine3D%get_carg('trail_rec'), 'the transition keeps trailing reconstruction')
            call assert_string_eq(params%rec_backend, cline_refine3D%get_carg('rec_backend'), 'backend selection is preserved')
            l_refine3D_lpstop_override = .true.
            params%lpstop = 7.5
            call set_cline_refine3D(params, 8, 1, .false.)
            call assert_false(cline_refine3D%defined('lp'), 'stage 8 still selects its matching band adaptively')
            call assert_real(7.5, cline_refine3D%get_rarg('lpstop'), 1.e-6, 'explicit user ceiling survives')
            call assert_string_eq('no', cline_refine3D%get_carg('automsk'), 'every stage keeps spherical support')
            call assert_string_eq('no', cline_refine3D%get_carg('envfsc'), 'FSC uses no generated envelope')
            l_refine3D_lpstop_override = .false.
        enddo
        params%filt_mode = 'nonuniform'
        call set_cline_refine3D(params, 6, 1, .false.)
        call assert_false(cline_refine3D%defined('lp'), 'stage 6 removes scheduled lp with plain NU input')
        call assert_string_eq('nonuniform', cline_refine3D%get_carg('filt_mode'), 'stage 6 keeps plain NU half references')
        addon%active = .true.
        addon%l_register = .true.
        addon%frozen_rec = 'frozen_context.txt'
        call set_cline_refine3D(params, 7, 1, .false., addon)
        call assert_true(cline_refine3D%defined('lp'), 'add-on registration keeps its fixed matching band')
        call assert_string_eq('none', cline_refine3D%get_carg('filt_mode'), 'add-on registration remains unfiltered')
        addon%l_register = .false.
        call set_cline_refine3D(params, 7, 1, .false., addon)
        call assert_false(cline_refine3D%defined('lp'), 'the add-on stage proper restores half-specific matching')
        call assert_string_eq('nonuniform', cline_refine3D%get_carg('filt_mode'), 'the add-on stage uses plain NU')
        call addon%frozen_rec%kill
        params%gold_standard = 'no'
        call set_cline_refine3D(params, 7, 1, .false.)
        call assert_true(cline_refine3D%defined('lp'), 'opt-out retains scheduled LP-set matching')
        call assert_string_eq('nonuniform_lpset', cline_refine3D%get_carg('filt_mode'), 'opt-out retains merged NU references')
        params%gold_standard = 'yes'
        params%nstates = 2
        call set_cline_refine3D(params, 7, 1, .false.)
        call assert_true(cline_refine3D%defined('lp'), 'extended multi-state schedule retains LP-set matching')
        call assert_string_eq('nonuniform_lpset', cline_refine3D%get_carg('filt_mode'), 'state count guards the switch')
        params%nstates = 1
        call set_cline_refine3D(params, 7, 1, .true.)
        call assert_true(cline_refine3D%defined('lp'), 'class-average route retains its scheduled matching band')
        call assert_string_eq('none', cline_refine3D%get_carg('filt_mode'), 'class-average route keeps NU off')
        params%filt_mode = 'none'
        l_nonuniform = .false.
        call set_cline_refine3D(params, 7, 1, .false.)
        call assert_false(cline_refine3D%defined('lp'), 'without NU, gold-standard matching uses the half-map FSC')
        call assert_false(cline_refine3D%defined('lpstop'), 'without NU, gold-standard matching still has no stage ceiling')
        call cline_refine3D%kill
        cline_refine3D = saved_cline
        deallocate(lpinfo)
        call move_alloc(saved_lpinfo, lpinfo)
        l_nonuniform = saved_flags(1)
        l_srch4symaxis = saved_flags(2)
        l_cavg_ini = saved_flags(3)
        l_state_continue_mode = saved_flags(4)
        l_refine3D_mode_override = saved_flags(5)
        l_refine3D_lp_override = saved_flags(6)
        l_refine3D_lpstop_override = saved_flags(7)
        l_refine3D_trs_override = saved_flags(8)
        l_prob_beta_end_override = saved_flags(9)
        update_frac = saved_update_frac
        nptcls_eff = saved_nptcls_eff
        nstages_refine3D = saved_nstages
        call saved_cline%kill
        call proj%kill
        deallocate(params, saved_cline, proj)
        call leave_fixture(cwd, root, nfail0)
    end subroutine test_gold_standard_stages

    !> the randomization acts on the poses the stage wrote into the project file, not on the in-memory
    !! copy of the previous stage: every pose in the file ends as itself or its c2 mate
    subroutine test_symmetry_randomize()
        use simple_sp_project, only: sp_project
        use simple_sym,        only: sym
        use simple_ori,        only: ori
        use simple_string,     only: string
        integer,          parameter :: NPTCLS = 12
        character(len=*), parameter :: PROJFILE = 'symrand.simple'
        class(sp_project), allocatable :: proj_file, proj_stale, proj_back
        type(sym)    :: symobj
        type(ori)    :: o_ref, o_mate, o_back, o_stale
        type(string) :: cwd, root
        integer      :: i, nfail0, nkept, nmoved, nstale
        real         :: mat_back(3,3), mat_stale(3,3,NPTCLS)
        write(*,'(A)') 'test_symmetry_randomize'
        nfail0 = tests_failed
        call enter_fixture('symrand', cwd, root)
        allocate(proj_file, proj_stale, proj_back)
        call proj_file%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            call proj_file%os_ptcl3D%set_euler(i, [real(mod(37*i,360)), 20.+6.*real(i), real(mod(53*i,360))])
        enddo
        call proj_file%update_projinfo(string(PROJFILE))
        call proj_file%write(string(PROJFILE))
        ! the in-memory project is the previous stage's: other poses than the file's
        call proj_stale%read(string(PROJFILE))
        do i = 1, NPTCLS
            call proj_stale%os_ptcl3D%set_euler(i, [real(mod(91*i+7,360)), 100.-5.*real(i), real(mod(29*i+3,360))])
        enddo
        do i = 1, NPTCLS
            call proj_stale%os_ptcl3D%get_ori(i, o_stale)
            mat_stale(:,:,i) = o_stale%get_mat()
        enddo
        call symobj%new('c2')
        call symmetry_randomize_poses(symobj, proj_stale, string(PROJFILE))
        call proj_back%read(string(PROJFILE))
        nkept = 0; nmoved = 0; nstale = 0
        do i = 1, NPTCLS
            call proj_file%os_ptcl3D%get_ori(i, o_ref)
            call symobj%apply(o_ref, 2, o_mate)
            call proj_back%os_ptcl3D%get_ori(i, o_back)
            mat_back = o_back%get_mat()
            if( maxval(abs(mat_back - o_ref%get_mat())) < 1.e-4 )  nkept  = nkept  + 1
            if( maxval(abs(mat_back - o_mate%get_mat())) < 1.e-4 ) nmoved = nmoved + 1
            if( maxval(abs(mat_back - mat_stale(:,:,i))) < 1.e-4 ) nstale = nstale + 1
        enddo
        call assert_int(NPTCLS, nkept + nmoved, 'each pose of the file ends as itself or its mate')
        call assert_int(0, nstale, 'the in-memory poses of the previous stage do not reach the file')
        call proj_file%kill; call proj_stale%kill; call proj_back%kill
        deallocate(proj_file, proj_stale, proj_back)
        call symobj%kill
        call leave_fixture(cwd, root, nfail0)
    end subroutine test_symmetry_randomize

    !> populated states, and a single state, are left as they are
    subroutine test_no_weak_state()
        integer, parameter :: NCAVGS = 12
        integer :: states(2*NCAVGS), nmoved
        real    :: scores(2*NCAVGS)
        integer :: i
        write(*,'(A)') 'test_no_weak_state'
        do i = 1, NCAVGS
            states(i)        = mod(i - 1, 3) + 1
            states(NCAVGS+i) = states(i)
            scores(i)        = real(i) / 100.
            scores(NCAVGS+i) = scores(i)
        enddo
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(0, nmoved, 'three populated states: nothing moves')
        states = 1
        call reseed_state_labels(states, scores, NCAVGS, 1, nmoved)
        call assert_int(0, nmoved, 'one state: nothing to reseed')
    end subroutine test_no_weak_state

    !> an empty state takes the donor's worst-fitting classes, both halves together, up to half
    !! the donor; deselected classes and the other state stay as they are
    subroutine test_one_empty_state()
        integer, parameter :: NCAVGS = 12
        real,    parameter :: SCORES1(6) = [0.9, 0.2, 0.8, 0.1, 0.7, 0.3] ! classes 4, 2 and 6 fit worst
        integer :: states(2*NCAVGS), expected(2*NCAVGS), nmoved, i
        real    :: scores(2*NCAVGS)
        write(*,'(A)') 'test_one_empty_state'
        ! classes 1-6 in state 1, 7-11 in state 2, 12 deselected, state 3 empty
        states(1:6)   = 1
        states(7:11)  = 2
        states(12)    = 0
        states(NCAVGS+1:2*NCAVGS) = states(1:NCAVGS)
        scores        = 0.5
        scores(1:6)   = SCORES1
        scores(NCAVGS+1:NCAVGS+6) = SCORES1
        expected = states
        do i = 1, NCAVGS
            if( any(i == [2, 4, 6]) )then
                expected(i)        = 3
                expected(NCAVGS+i) = 3
            endif
        enddo
        ! 22 selected halves: an equal share is 7, half the donor's 12 is 6
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(6, nmoved, 'half the donor moves')
        call assert_true(all(states == expected), 'the three worst-fitting classes of state 1, both halves')
        call assert_true(states(12) == 0 .and. states(2*NCAVGS) == 0, 'the deselected class stays deselected')
    end subroutine test_one_empty_state

    !> every class in one state: the two empty states each take an equal share, the worst-fitting first
    subroutine test_all_in_one_state()
        integer, parameter :: NCAVGS = 12
        integer :: states(2*NCAVGS), nmoved, i
        real    :: scores(2*NCAVGS)
        write(*,'(A)') 'test_all_in_one_state'
        states = 1
        do i = 1, NCAVGS
            scores(i)        = real(i) / 100.
            scores(NCAVGS+i) = scores(i)
        enddo
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(16, nmoved, 'two equal shares of 8 move')
        call assert_true(all(states(1:4) == 2) .and. all(states(NCAVGS+1:NCAVGS+4) == 2),&
            &'state 2 takes the four worst-fitting classes')
        call assert_true(all(states(5:8) == 3) .and. all(states(NCAVGS+5:NCAVGS+8) == 3),&
            &'state 3 the next four')
        call assert_true(all(states(9:12) == 1) .and. all(states(NCAVGS+9:NCAVGS+12) == 1),&
            &'state 1 keeps the best-fitting four')
    end subroutine test_all_in_one_state

    !> a state above the search's floor but under 2% of the selected halves is reseeded
    subroutine test_weak_by_fraction()
        integer, parameter :: NCAVGS = 200
        integer :: states(2*NCAVGS), nmoved, i
        real    :: scores(2*NCAVGS)
        write(*,'(A)') 'test_weak_by_fraction'
        states(1:100)   = 1
        states(101:197) = 2
        states(198:200) = 3 ! 6 halves of 400
        states(NCAVGS+1:2*NCAVGS) = states(1:NCAVGS)
        do i = 1, NCAVGS
            scores(i)        = 0.5 + real(i) / 1000.
            scores(NCAVGS+i) = scores(i)
        enddo
        ! an equal share is 133, half the donor's 200 is 100
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(100, nmoved, 'half the donor moves')
        call assert_int(106, count(states == 3), 'state 3 holds its own and the moved halves')
        call assert_true(all(states(1:50) == 3) .and. all(states(51:100) == 1),&
            &'the worst-fitting 50 classes of state 1 move')
    end subroutine test_weak_by_fraction

    !> too few classes to leave a reseeded state above the search's floor: nothing moves, so no
    !! state is left empty while another is relabelled
    subroutine test_too_few_to_reseed()
        integer, parameter :: NCAVGS = 4
        integer :: states(2*NCAVGS), states_in(2*NCAVGS), nmoved, i
        real    :: scores(2*NCAVGS)
        write(*,'(A)') 'test_too_few_to_reseed'
        states = 1
        do i = 1, NCAVGS
            scores(i)        = real(i) / 100.
            scores(NCAVGS+i) = scores(i)
        enddo
        states_in = states
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(0, nmoved, 'nothing moves')
        call assert_true(all(states == states_in), 'the labels are unchanged')
    end subroutine test_too_few_to_reseed

end module simple_solve3D_reseed_tester
