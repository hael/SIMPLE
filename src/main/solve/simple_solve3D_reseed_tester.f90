!@descr: unit tests of solve3D stage rules, symmetry randomization, and solve3D_cavgs state reseeding
! In-memory labels and scores of the even (1..ncavgs) and odd (ncavgs+1..2*ncavgs) halves of the
! class averages; no volume or job. The symmetry randomization and sticky-release tests each write one
! small project file.
module simple_solve3D_reseed_tester
use simple_solve3D_utils, only: reseed_state_labels, solve3D_symsrch_stage, solve3D_stage_skipping, &
    &symmetry_randomize_poses, solve3D_size_lplims
use simple_test_utils
implicit none
private
public :: run_all_solve3D_reseed_tests

contains

    subroutine run_all_solve3D_reseed_tests()
        write(*,'(A)') '**** running all solve3D_cavgs stage rule tests ****'
        call test_symsrch_stage()
        call test_stage_skipping()
        call test_size_lplims()
        call test_ptcl_lplims()
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

    !> the class-average route's size limits: mskdiam/8 within 10-30 A, mskdiam/26.7 within 4.5-6 A
    subroutine test_size_lplims()
        real :: lims(2)
        write(*,'(A)') 'test_size_lplims'
        lims = solve3D_size_lplims(160.)
        call assert_real(20.,    lims(1), 1.e-4, '160 A mask (PfCRT): starts at 20 A')
        call assert_real(6.,     lims(2), 1.e-4, '160 A mask: ends at 6 A')
        lims = solve3D_size_lplims(80.)
        call assert_real(10.,    lims(1), 1.e-4, '80 A mask (streptavidin): starts at 10 A')
        call assert_real(4.5,    lims(2), 1.e-4, '80 A mask: the 4.5 A bound')
        lims = solve3D_size_lplims(100.)
        call assert_real(12.5,   lims(1), 1.e-4, '100 A mask: mskdiam/8')
        lims = solve3D_size_lplims(50.)
        call assert_real(10.,    lims(1), 1.e-4, '50 A mask: the 10 A floor')
        lims = solve3D_size_lplims(300.)
        call assert_real(30.,    lims(1), 1.e-4, '300 A mask: the 30 A cap')
        call assert_real(6.,     lims(2), 1.e-4, '300 A mask: the 6 A bound')
    end subroutine test_size_lplims

    !> the particle ladder: mskdiam/10 within 10-30 A to 4.5 A, each stage halving the distance to the stop;
    !! an explicit lpstart or lpstop sets its end as given, coarser or finer
    subroutine test_ptcl_lplims()
        use simple_solve3D_utils, only: lpinfo, nstages_refine3D, l_cavgs_mode, l_refine3D_trs_override, &
            &set_lplims_ptcls, solve3D_ptcl_lplims
        use simple_parameters,    only: parameters
        use simple_type_defs,     only: lp_crop_inf
        use simple_string_utils,  only: int2str
        real, parameter :: LPS_EXPORT_GATE(8) = [15.5, 10., 7.25, 5.875, 5.1875, 4.84375, 4.671875, 4.5]
        real, parameter :: LPS_FINE_STOP(8)   = [15.5, 9.65, 6.725, 5.2625, 4.53125, 4.165625, 3.982813, 3.8]
        class(parameters), allocatable :: params
        type(lp_crop_inf), allocatable :: saved_lpinfo(:)
        real    :: lims(2), stage1_lp
        integer :: saved_nstages, i
        logical :: saved_cavgs_mode, saved_trs_override
        write(*,'(A)') 'test_ptcl_lplims'
        lims = solve3D_ptcl_lplims(155.)
        call assert_real(15.5,  lims(1), 1.e-4, '155 A mask (export gate): starts at mskdiam/10')
        call assert_real(4.5,   lims(2), 1.e-4, '155 A mask: ends at 4.5 A')
        lims = solve3D_ptcl_lplims(160.)
        call assert_real(16.,   lims(1), 1.e-4, '160 A mask (PfCRT): starts at 16 A')
        lims = solve3D_ptcl_lplims(80.)
        call assert_real(10.,   lims(1), 1.e-4, '80 A mask (streptavidin): the 10 A floor')
        call assert_real(4.5,   lims(2), 1.e-4, '80 A mask: ends at 4.5 A')
        lims = solve3D_ptcl_lplims(400.)
        call assert_real(30.,   lims(1), 1.e-4, '400 A mask: the 30 A cap')
        allocate(params)
        saved_nstages      = nstages_refine3D
        saved_cavgs_mode   = l_cavgs_mode
        saved_trs_override = l_refine3D_trs_override
        if( allocated(lpinfo) ) call move_alloc(lpinfo, saved_lpinfo)
        l_refine3D_trs_override = .false.
        params%box     = 288
        params%smpd    = 0.822
        params%mskdiam = 155.
        nstages_refine3D = 8
        call set_lplims_ptcls(params, 0., 0., stage1_lp)
        call assert_real(19.375, stage1_lp, 1.e-4, 'export gate: stage 1 separately starts at mskdiam/8')
        call assert_false(l_cavgs_mode, 'the particle ladder clears the class-average mode')
        call assert_int(8, size(lpinfo), 'one rung per stage')
        do i = 1,8
            call assert_real(LPS_EXPORT_GATE(i), lpinfo(i)%lp, 1.e-3, 'export gate, stage '//int2str(i)//': halving')
        end do
        ! an explicit lpstop finer than 4.5 A sets the end as given
        call set_lplims_ptcls(params, 0., 3.8, stage1_lp)
        call assert_real(19.375, stage1_lp, 1.e-4, 'a finer stop leaves the stage-1 size rule unchanged')
        do i = 1,8
            call assert_real(LPS_FINE_STOP(i), lpinfo(i)%lp, 1.e-3, 'lpstop 3.8, stage '//int2str(i)//': as given')
        end do
        ! an explicit lpstart coarser than the size rule sets the start as given
        call set_lplims_ptcls(params, 25., 0., stage1_lp)
        call assert_real(25., stage1_lp, 1.e-4, 'an explicit lpstart also controls stage 1')
        call assert_real(25., lpinfo(1)%lp, 1.e-4, 'lpstart 25: the first rung as given')
        call assert_real(4.5, lpinfo(8)%lp, 1.e-4, 'lpstart 25: the last rung from the size rule')
        ! fewer stages: the same two ends over the stages the run executes
        nstages_refine3D = 5
        call set_lplims_ptcls(params, 0., 0., stage1_lp)
        call assert_int(5, size(lpinfo), 'five stages: five rungs')
        call assert_real(4.5, lpinfo(5)%lp, 1.e-4, 'five stages: ends at 4.5 A')
        ! a start finer than the stop holds the ladder at the stop
        nstages_refine3D = 4
        call set_lplims_ptcls(params, 4., 5., stage1_lp)
        call assert_real(5., stage1_lp, 1.e-4, 'an explicit start finer than the stop respects the stop')
        call assert_true(all(abs(lpinfo(:)%lp - 5.) < 1.e-4), 'a start finer than the stop: flat at the stop')
        call set_lplims_ptcls(params, 0., 25., stage1_lp)
        call assert_real(25., stage1_lp, 1.e-4, 'a stop coarser than the stage-1 size rule remains a ceiling')
        deallocate(lpinfo)
        if( allocated(saved_lpinfo) ) call move_alloc(saved_lpinfo, lpinfo)
        nstages_refine3D        = saved_nstages
        l_cavgs_mode            = saved_cavgs_mode
        l_refine3D_trs_override = saved_trs_override
        deallocate(params)
    end subroutine test_ptcl_lplims

    !> stage omission and the mid-stage resolution exit: single-state runs only
    subroutine test_stage_skipping()
        write(*,'(A)') 'test_stage_skipping'
        call assert_true(solve3D_stage_skipping(1),        'single state: stages may be skipped')
        call assert_false(solve3D_stage_skipping(2),       'two states: no stage skipped')
        call assert_false(solve3D_stage_skipping(4),       'four states: no stage skipped')
    end subroutine test_stage_skipping

    !> Exercise emitted stage commands, including a skipped transition, stale matching limits and the sampling.
    subroutine test_gold_standard_stages()
        use simple_solve3D_utils, only: lpinfo, cline_refine3D, set_cline_refine3D, solve3D_addon_ctx, &
            &l_nonuniform, l_srch4symaxis, l_cavg_ini, l_state_continue_mode, &
            &l_refine3D_mode_override, l_refine3D_lp_override, l_refine3D_lpstop_override, &
            &l_refine3D_trs_override, l_prob_beta_end_override, update_frac, nptcls_eff, nstages_refine3D, &
            &set_lplims_ptcls, l_cavgs_mode
        use simple_parameters, only: parameters
        use simple_cmdline,    only: cmdline
        use simple_type_defs,  only: lp_crop_inf
        use simple_sp_project, only: sp_project
        use simple_string,     only: string
        use simple_defs_conv,  only: OVERLAP_3D, OVERLAP_3D_STAGE
        use simple_defs,       only: logfhandle
        use simple_final_rec,  only: final_rec_mask_settings
        class(parameters), allocatable :: params
        class(cmdline), allocatable    :: saved_cline
        class(sp_project), allocatable :: proj
        type(lp_crop_inf), allocatable :: saved_lpinfo(:), base_lpinfo(:)
        type(solve3D_addon_ctx) :: addon
        type(string) :: cwd, root, final_automsk
        logical :: saved_flags(10), l_final_envfsc, l_final_nu_mask
        real    :: saved_update_frac, stage1_lp, planned_lp(8)
        real, parameter :: START_DIAMS(6) = [50., 80., 155., 160., 170., 400.]
        real, parameter :: START_LPS(6)   = [10., 10., 19.375, 20., 21.25, 30.]
        real, parameter :: LADDER_LPS(6)  = [10., 10., 15.5, 16., 17., 30.]
        integer :: saved_nptcls_eff, saved_nstages, nfail0, backend, istage, imode, i, stage1_box
        integer :: saved_logfhandle, log_unit
        character(len=7), parameter :: mask_modes(4) = [character(len=7) :: 'no', 'density', 'nu', 'order']
        write(*,'(A)') 'test_gold_standard_stages'
        nfail0 = tests_failed
        call enter_fixture('solve3D_gold_standard', cwd, root)
        allocate(params, saved_cline, proj)
        saved_cline = cline_refine3D
        saved_flags = [l_nonuniform, l_srch4symaxis, l_cavg_ini, l_state_continue_mode, &
            &l_refine3D_mode_override, l_refine3D_lp_override, l_refine3D_lpstop_override, &
            &l_refine3D_trs_override, l_prob_beta_end_override, l_cavgs_mode]
        saved_update_frac = update_frac
        saved_nptcls_eff  = nptcls_eff
        saved_nstages     = nstages_refine3D
        call move_alloc(lpinfo, saved_lpinfo)
        allocate(lpinfo(8))
        lpinfo%lp = [20., 16., 12., 10., 8., 6., 5., 4.5]
        lpinfo%box_crop = 64
        lpinfo%trslim   = 2.
        base_lpinfo = lpinfo
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
        params%gold_std_stage = 0
        params%automsk = 'density'
        params%envfsc = 'yes'
        params%projfile = 'stages.simple'
        call proj%os_ptcl3D%new(2, is_ptcl=.true.)
        call proj%os_ptcl3D%set_all2single('state', 1.)
        call proj%os_ptcl3D%set_all2single('res05', 5.3)
        call proj%update_projinfo(params%projfile)
        call proj%write(params%projfile)
        ! The first-stage matching and Gaussian limits are separate from the ladder and its FSC bounds.
        saved_logfhandle = logfhandle
        open(newunit=log_unit, status='scratch', action='write')
        logfhandle = log_unit
        params%box  = 240
        params%smpd = 1.3
        call write_res05(40.)
        do i = 1,size(START_DIAMS)
            params%mskdiam = START_DIAMS(i)
            call set_lplims_ptcls(params, 0., 0., stage1_lp)
            planned_lp = lpinfo%lp
            stage1_box = lpinfo(1)%box_crop
            call set_cline_refine3D(params, 1, 1, .false., stage1_lp=stage1_lp)
            call assert_real(START_LPS(i), cline_refine3D%get_rarg('lp'), 1.e-4, 'stage 1 matches at the size-rule limit')
            call assert_real(START_LPS(i), cline_refine3D%get_rarg('lpstop'), 1.e-4, 'stage 1 ceiling follows its limit')
            call assert_real(START_LPS(i), cline_refine3D%get_rarg('gaufreq'), 1.e-4, 'stage 1 Gaussian follows its limit')
            call assert_string_eq('yes', cline_refine3D%get_carg('gauref'), 'stage 1 retains Gaussian filtering')
            call assert_int(stage1_box, cline_refine3D%get_iarg('box_crop'), 'stage 1 retains the ladder crop')
            call assert_true(all(abs(lpinfo%lp-planned_lp) < 1.e-6), 'stage 1 does not mutate the ladder')
            do istage = 2,5
                call set_cline_refine3D(params, istage, 1, .false., stage1_lp=stage1_lp)
                call assert_real(LADDER_LPS(i), cline_refine3D%get_rarg('lp'), 1.e-4, &
                    &'later FSC hold-back remains bounded by the original ladder start')
                if( istage == 2 ) call assert_real(LADDER_LPS(i), cline_refine3D%get_rarg('gaufreq'), 1.e-4, &
                    &'stage 2 Gaussian retains its FSC-derived limit')
            enddo
        enddo
        params%mskdiam = 170.
        call set_lplims_ptcls(params, 25., 0., stage1_lp)
        call set_cline_refine3D(params, 1, 1, .false., stage1_lp=stage1_lp)
        call assert_real(25., cline_refine3D%get_rarg('lp'), 1.e-4, 'explicit lpstart controls stage-1 matching')
        call assert_real(25., cline_refine3D%get_rarg('gaufreq'), 1.e-4, 'explicit lpstart controls stage-1 Gaussian')
        call set_lplims_ptcls(params, 0., 0., stage1_lp)
        call set_cline_refine3D(params, 1, 1, .true., stage1_lp=stage1_lp)
        call assert_real(17., cline_refine3D%get_rarg('lp'), 1.e-4, 'class-average stages retain their own ladder')
        call assert_real(17., cline_refine3D%get_rarg('gaufreq'), 1.e-4, 'class-average Gaussian retains its own ladder')
        logfhandle = saved_logfhandle
        close(log_unit)
        lpinfo = base_lpinfo
        deallocate(base_lpinfo)
        params%box     = 64
        params%smpd    = 2.
        params%mskdiam = 0.
        call write_res05(5.3)
        ! no stage is skipped: every stage runs 3 iterations, then its exits decide
        call set_cline_refine3D(params, 4, 1, .false.)
        call assert_int(15, cline_refine3D%get_iarg('maxits'), 'stage 4 keeps its budget behind FSC=0.5')
        call assert_int(3,  cline_refine3D%get_iarg('minits'), 'stage 4 runs at least 3 iterations')
        ! the resolution exit from stage 3 on aims at the next stage's limit
        call set_cline_refine3D(params, 3, 1, .false.)
        call assert_real(10., cline_refine3D%get_rarg('lpexit'), 1.e-6, 'stage 3 exits past stage 4''s limit')
        call set_cline_refine3D(params, 5, 1, .false.)
        call assert_real(6., cline_refine3D%get_rarg('lpexit'), 1.e-6, 'stage 5 exits past stage 6''s limit')
        ! single state, from stage 2 until gold standard: FSC=0.5 at stage entry sets the limit, both ways
        call set_cline_refine3D(params, 1, 1, .false.)
        call assert_real(20., cline_refine3D%get_rarg('lp'), 1.e-6, 'stage 1 keeps its planned limit')
        call set_cline_refine3D(params, 2, 1, .false.)
        call assert_real(5.3, cline_refine3D%get_rarg('lp'),      1.e-6, 'stage 2 matches at FSC=0.5')
        call assert_real(5.3, cline_refine3D%get_rarg('gaufreq'), 1.e-6, 'stage 2''s Gaussian follows its limit')
        call set_cline_refine3D(params, 3, 1, .false.)
        call assert_real(5.3, cline_refine3D%get_rarg('lp'),     1.e-6, 'stage 3 matches at FSC=0.5')
        call assert_real(5.3, cline_refine3D%get_rarg('lpstop'), 1.e-6, 'stage 3''s ceiling is its limit')
        call write_res05(15.)
        call set_cline_refine3D(params, 4, 1, .false.)
        call assert_real(15., cline_refine3D%get_rarg('lp'), 1.e-6, 'a coarser FSC=0.5 than planned holds stage 4 back')
        call write_res05(30.)
        call set_cline_refine3D(params, 4, 1, .false.)
        call assert_real(20., cline_refine3D%get_rarg('lp'), 1.e-6, 'the hold-back stops at the ladder start')
        call write_res05(9.)
        call set_cline_refine3D(params, 4, 1, .false.)
        call assert_real(9., cline_refine3D%get_rarg('lp'), 1.e-6, 'a finer FSC=0.5 than planned sets stage 4''s limit')
        call write_res05(3.)
        call set_cline_refine3D(params, 4, 1, .false.)
        call assert_real(4.5, cline_refine3D%get_rarg('lp'), 1.e-6, 'the ladder cap bounds the FSC=0.5 limit')
        l_refine3D_lpstop_override = .true.
        params%lpstop = 3.8
        call set_cline_refine3D(params, 4, 1, .false.)
        call assert_real(3.8, cline_refine3D%get_rarg('lp'), 1.e-6, 'a finer command-line lpstop replaces the cap')
        l_refine3D_lpstop_override = .false.
        call write_res05(5.3)
        params%nstates = 2
        call set_cline_refine3D(params, 4, 1, .false.)
        call assert_real(8., cline_refine3D%get_rarg('lp'), 1.e-6, 'multi-state stages keep the one-rung advance')
        params%nstates = 1
        l_nonuniform = .false.
        call set_cline_refine3D(params, 5, 1, .false.)
        call assert_real(5.3, cline_refine3D%get_rarg('lp'), 1.e-6, 'without NU, stage 5 matches at FSC=0.5')
        l_nonuniform = .true.
        do backend = 1, 2
            params%rec_backend = 'gridding'
            if( backend == 2 ) params%rec_backend = 'pcg'
            params%filt_mode = 'nonuniform_lpset'
            call set_cline_refine3D(params, 5, 1, .false.)
            call assert_string_eq('nonuniform_lpset', cline_refine3D%get_carg('filt_mode'), 'stage 5 uses merged NU references')
            call assert_real(5.3, cline_refine3D%get_rarg('lp'), 1.e-6, 'FSC=0.5 sets NU stage 5''s limit')
            call assert_false(cline_refine3D%defined('lpstop'), 'NU stage 5 has no added stage ceiling')
            ! Stage 6 is passed over here: stage 7 must clear the preceding LP-set command.
            call set_cline_refine3D(params, 7, 1, .false.)
            call assert_string_eq('nonuniform', cline_refine3D%get_carg('filt_mode'), 'stage 7 uses half-specific NU references')
            call assert_false(cline_refine3D%defined('lp'), 'transition clears the preceding scheduled lp')
            call assert_false(cline_refine3D%defined('lpstop'), 'gold-standard stage has no scheduled ceiling')
            call assert_real(0.25, cline_refine3D%get_rarg('update_frac'), 1.e-6, 'the transition keeps fractional updates')
            call assert_string_eq('yes', cline_refine3D%get_carg('trail_rec'), 'the transition keeps trailing reconstruction')
            call assert_string_eq(params%rec_backend, cline_refine3D%get_carg('rec_backend'), 'backend selection is preserved')
            call assert_string_eq('density', cline_refine3D%get_carg('automsk'), 'automsk=density: stage 7 masks by the density envelope')
            call assert_string_eq('yes', cline_refine3D%get_carg('envfsc'), 'automsk=density: stage 7 FSC uses the envelope')
            l_refine3D_lpstop_override = .true.
            params%lpstop = 7.5
            call set_cline_refine3D(params, 8, 1, .false.)
            call assert_false(cline_refine3D%defined('lp'), 'stage 8 still selects its matching band adaptively')
            call assert_real(7.5, cline_refine3D%get_rarg('lpstop'), 1.e-6, 'explicit user ceiling survives')
            call assert_string_eq('density', cline_refine3D%get_carg('automsk'), 'automsk=density: stage 8 masks by the density envelope')
            call assert_string_eq('yes', cline_refine3D%get_carg('envfsc'), 'automsk=density: stage 8 FSC uses the envelope')
            params%automsk = 'no'
            call set_cline_refine3D(params, 8, 1, .false.)
            call assert_string_eq('no', cline_refine3D%get_carg('automsk'), 'automsk=no: stage 8 keeps spherical support')
            call assert_string_eq('no', cline_refine3D%get_carg('envfsc'), 'automsk=no: FSC uses no generated envelope')
            params%automsk = 'density'
            l_refine3D_lpstop_override = .false.
        enddo
        do imode = 1, size(mask_modes)
            params%automsk = mask_modes(imode)
            do istage = 1, 8
                call set_cline_refine3D(params, istage, 1, .false.)
                if( istage < 7 .or. trim(mask_modes(imode)) == 'no' )then
                    call assert_string_eq('no', cline_refine3D%get_carg('automsk'), 'early or default stage has no mask')
                    call assert_string_eq('no', cline_refine3D%get_carg('envfsc'), 'early or default stage has no envelope FSC')
                else
                    call assert_string_eq(trim(mask_modes(imode)), cline_refine3D%get_carg('automsk'), &
                        &'late stage forwards the requested mask')
                    call assert_string_eq('yes', cline_refine3D%get_carg('envfsc'), 'late stage enables envelope FSC')
                endif
            enddo
        enddo
        ! the final reconstruction takes its mask settings from the last stage's command line only: a run
        ! that ends before stage 7 ships unmasked although automsk was requested, one that reaches it masked
        params%automsk   = 'order'
        params%filt_mode = 'nonuniform'
        call set_cline_refine3D(params, 5, 1, .false.)
        call final_rec_mask_settings(cline_refine3D, final_automsk, l_final_envfsc, l_final_nu_mask)
        call assert_string_eq('no', final_automsk, 'a run ending at stage 5 ships without the requested mask')
        call assert_false(l_final_envfsc,  'a run ending at stage 5 ships without envelope FSC')
        call assert_false(l_final_nu_mask, 'a run ending at stage 5 ships without the NU mask filter mode')
        call set_cline_refine3D(params, 8, 1, .false.)
        call final_rec_mask_settings(cline_refine3D, final_automsk, l_final_envfsc, l_final_nu_mask)
        call assert_string_eq('order', final_automsk, 'a run ending at stage 8 ships with the order mask')
        call assert_true(l_final_envfsc,  'a run ending at stage 8 ships with envelope FSC')
        call assert_true(l_final_nu_mask, 'a run ending at stage 8 keeps the NU filter mode for the order mask')
        params%automsk = 'density'
        call set_cline_refine3D(params, 8, 1, .false.)
        call final_rec_mask_settings(cline_refine3D, final_automsk, l_final_envfsc, l_final_nu_mask)
        call assert_string_eq('density', final_automsk, 'a run ending at stage 8 ships with the density mask')
        call assert_true(l_final_envfsc,  'the density mask implies envelope FSC')
        call assert_false(l_final_nu_mask, 'the density mask ships a classical map')
        params%automsk = 'density'
        params%filt_mode = 'nonuniform'
        call set_cline_refine3D(params, 6, 1, .false.)
        call assert_false(cline_refine3D%defined('lp'), 'stage 6 removes scheduled lp with plain NU input')
        call assert_string_eq('nonuniform', cline_refine3D%get_carg('filt_mode'), 'stage 6 keeps plain NU half references')
        call assert_string_eq('gridding', cline_refine3D%get_carg('rec_backend'), 'PCG starts at stage 7')
        call assert_string_eq('no', cline_refine3D%get_carg('automsk'), 'automsk=density: stage 6 keeps spherical support')
        call assert_string_eq('no', cline_refine3D%get_carg('envfsc'), 'automsk=density: stage 6 FSC uses no envelope')
        call set_cline_refine3D(params, 7, 1, .false.)
        call assert_string_eq('pcg', cline_refine3D%get_carg('rec_backend'), 'stage 7 runs the requested PCG backend')
        addon%active = .true.
        addon%l_register = .true.
        addon%frozen_rec = 'frozen_context.txt'
        call set_cline_refine3D(params, 7, 1, .false., addon)
        call assert_true(cline_refine3D%defined('lp'), 'add-on registration keeps its fixed matching band')
        call assert_string_eq('none', cline_refine3D%get_carg('filt_mode'), 'add-on registration remains unfiltered')
        call assert_string_eq('no', cline_refine3D%get_carg('automsk'), 'add-on registration keeps spherical support')
        addon%l_register = .false.
        call set_cline_refine3D(params, 7, 1, .false., addon)
        call assert_false(cline_refine3D%defined('lp'), 'the add-on stage proper restores half-specific matching')
        call assert_string_eq('nonuniform', cline_refine3D%get_carg('filt_mode'), 'the add-on stage uses plain NU')
        call assert_string_eq('density', cline_refine3D%get_carg('automsk'), 'the add-on stage replays automsk=density')
        call addon%frozen_rec%kill
        params%gold_std_stage = 9
        call set_cline_refine3D(params, 7, 1, .false.)
        call assert_true(cline_refine3D%defined('lp'), 'a later switch retains scheduled LP-set matching')
        call assert_string_eq('nonuniform_lpset', cline_refine3D%get_carg('filt_mode'), 'a later switch retains merged NU references')
        params%gold_std_stage = 0
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
        ! no stage passes a score-based sampling control: one state and several draw from the whole class;
        ! one state refines one sticky subset through stage 4: drawn in stage 1 by the 2D objective, best first,
        ! kept in stages 2-4, free from stage 5, where trailing reconstruction starts
        do istage = 1, 8
            call set_cline_refine3D(params, istage, 1, .false.)
            call assert_int(min(3, cline_refine3D%get_iarg('maxits')), cline_refine3D%get_iarg('minits'), &
                &'ordinary stages run the three-iteration minimum, then early-stop')
            call assert_false(cline_refine3D%defined('greedy_sampling'), 'one state: no greedy sampling control')
            call assert_false(cline_refine3D%defined('frac_best'),       'one state: no best-fraction control')
            select case(istage)
                case(1)
                    call assert_string_eq('best', cline_refine3D%get_carg('sticky_sampling'), 'one state: stage 1 draws the sticky subset best first')
                case(2:4)
                    call assert_string_eq('keep', cline_refine3D%get_carg('sticky_sampling'), 'one state: stages 2-4 keep it')
                case DEFAULT
                    call assert_false(cline_refine3D%defined('sticky_sampling'), 'one state: free sampling from stage 5')
            end select
            ! an intermediate stage from stage 3 moves on at overlap 0.95; stages 1-2 and the last keep 0.99
            if( istage >= 3 .and. istage < 8 )then
                call assert_real(OVERLAP_3D_STAGE, cline_refine3D%get_rarg('overlap'), 1.e-6, 'one state: an intermediate stage moves on at 0.95')
            else
                call assert_real(OVERLAP_3D, cline_refine3D%get_rarg('overlap'), 1.e-6, 'one state: stages 1-2 and the last keep 0.99')
            endif
            ! trailing reconstruction starts with the free sampling
            if( istage <= 4 )then
                call assert_string_eq('no',  cline_refine3D%get_carg('trail_rec'), 'one state: no trailing over the sticky subset')
            else
                call assert_string_eq('yes', cline_refine3D%get_carg('trail_rec'), 'one state: trailing from stage 5')
            endif
        end do
        l_srch4symaxis = .true.
        call set_cline_refine3D(params, 4, 1, .false.)
        call assert_int(cline_refine3D%get_iarg('maxits'), cline_refine3D%get_iarg('minits'), 'the axis-search stage runs its full budget')
        call assert_real(OVERLAP_3D, cline_refine3D%get_rarg('overlap'), 1.e-6, 'the axis-search stage keeps the strict overlap')
        call assert_false(cline_refine3D%defined('lpexit'), 'the axis-search stage has no resolution exit')
        call set_cline_refine3D(params, 5, 1, .false.)
        call assert_int(3, cline_refine3D%get_iarg('minits'), 'the stage after the axis search returns to the three-iteration minimum')
        nstages_refine3D = 4
        call set_cline_refine3D(params, 4, 1, .true.)
        call assert_int(cline_refine3D%get_iarg('maxits'), cline_refine3D%get_iarg('minits'), 'the class-average axis-search stage keeps its actual budget')
        nstages_refine3D = 8
        l_srch4symaxis = .false.
        l_cavg_ini = .true.
        call set_cline_refine3D(params, 4, 4, .false.)
        call assert_int(cline_refine3D%get_iarg('maxits'), cline_refine3D%get_iarg('minits'), 'the particle entry after class averages runs its full budget')
        call assert_real(OVERLAP_3D, cline_refine3D%get_rarg('overlap'), 1.e-6, 'the particle entry keeps the strict overlap')
        call assert_string_eq('best', cline_refine3D%get_carg('sticky_sampling'), 'the particle entry after class averages holds a sticky subset')
        call set_cline_refine3D(params, 5, 4, .false.)
        call assert_int(3, cline_refine3D%get_iarg('minits'), 'the stage after particle entry returns to the three-iteration minimum')
        addon%active = .true.
        addon%l_register = .false.
        addon%frozen_rec = 'frozen_context_for_stage_rule'
        call set_cline_refine3D(params, 4, 4, .false., addon=addon)
        call assert_int(3, cline_refine3D%get_iarg('minits'), 'an add-on has no full-budget stage, only the three-iteration minimum')
        call assert_false(cline_refine3D%defined('sticky_sampling'), 'an add-on keeps its cohort''s sampling: no sticky subset')
        addon%active = .false.
        call addon%frozen_rec%kill
        l_cavg_ini = .false.
        ! a run that starts at stage 2 has no sticky subset to keep: it draws one and holds it through stage 4
        call set_cline_refine3D(params, 2, 2, .false.)
        call assert_string_eq('best', cline_refine3D%get_carg('sticky_sampling'), 'one state: a stage-2 start draws its sticky subset')
        call set_cline_refine3D(params, 4, 2, .false.)
        call assert_string_eq('keep', cline_refine3D%get_carg('sticky_sampling'), 'one state: a stage-2 start keeps it in stage 4')
        ! a run that enters at stage 4 (supplied volumes, cavg_ini_ext) holds a sticky subset for that stage and
        ! trails from stage 5
        call set_cline_refine3D(params, 4, 4, .false.)
        call assert_string_eq('best', cline_refine3D%get_carg('sticky_sampling'), 'one state: a stage-4 entry draws its sticky subset')
        call assert_string_eq('no', cline_refine3D%get_carg('trail_rec'), 'one state: a stage-4 entry does not trail over it')
        call set_cline_refine3D(params, 5, 4, .false.)
        call assert_false(cline_refine3D%defined('sticky_sampling'), 'one state: a stage-4 entry samples freely from stage 5')
        call assert_string_eq('yes', cline_refine3D%get_carg('trail_rec'), 'one state: a stage-4 entry trails from stage 5')
        ! a run that enters after the sticky stages samples freely and trails from its first stage
        call set_cline_refine3D(params, 5, 5, .false.)
        call assert_false(cline_refine3D%defined('sticky_sampling'), 'one state: a stage-5 entry has no sticky subset')
        call assert_string_eq('yes', cline_refine3D%get_carg('trail_rec'), 'one state: a stage-5 entry trails from its first stage')
        ! full sampling and the class-average route hold none
        params%nsample = 9500
        call set_cline_refine3D(params, 1, 1, .false.)
        call assert_false(cline_refine3D%defined('sticky_sampling'), 'full sampling: no sticky subset')
        call set_cline_refine3D(params, 5, 1, .false.)
        call assert_string_eq('no', cline_refine3D%get_carg('trail_rec'), 'full sampling: no trailing')
        params%nsample = 1000
        call set_cline_refine3D(params, 1, 1, .true.)
        call assert_false(cline_refine3D%defined('sticky_sampling'), 'class-average route: no sticky subset')
        call check_early_stage_sticky(params, 1, 1)
        call check_early_stage_sticky(params, 2, 1)
        call check_early_stage_sticky(params, 2, 2)
        call check_early_stage_sticky(params, 4, 1)
        call check_early_stage_sticky(params, 4, 2)
        call check_early_stage_sticky(params, 4, 4)
        call check_sticky_release(params, 'none',  .true.)
        call check_sticky_release(params, 'none',  .false.)
        call check_sticky_release(params, 'class', .true.)
        call check_sticky_release(params, 'class', .false.)
        call check_sticky_release(params, 'cavg',  .true.)
        call check_sticky_release(params, 'cavg',  .false.)
        nstages_refine3D = 8
        params%nstates = 3
        do istage = 1, 8
            call set_cline_refine3D(params, istage, 1, .false.)
            call assert_false(cline_refine3D%defined('greedy_sampling'), 'several states: no greedy sampling control')
            call assert_false(cline_refine3D%defined('frac_best'),       'several states: no best-fraction control')
            call assert_false(cline_refine3D%defined('sticky_sampling'), 'several states: no sticky subset')
            if( istage < 8 )then
                call assert_string_eq('no',  cline_refine3D%get_carg('trail_rec'), 'several states: no trailing before stage 8')
            else
                call assert_string_eq('yes', cline_refine3D%get_carg('trail_rec'), 'several states: trailing from stage 8')
            endif
        end do
        call check_sigma2_schedules(params)
        params%nstates = 1
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
        l_cavgs_mode = saved_flags(10)
        update_frac = saved_update_frac
        nptcls_eff = saved_nptcls_eff
        nstages_refine3D = saved_nstages
        call saved_cline%kill
        call proj%kill
        deallocate(params, saved_cline, proj)
        call leave_fixture(cwd, root, nfail0)

    contains

        subroutine write_res05( res05 )
            real, intent(in) :: res05
            call proj%os_ptcl3D%set_all2single('res05', res05)
            call proj%write(params%projfile)
        end subroutine write_res05

    end subroutine test_gold_standard_stages

    !> Hand-planned boundaries: short runs reserve two measured stages; full sampling and states share the schedule.
    subroutine check_sigma2_schedules( params )
        use simple_solve3D_utils, only: cline_refine3D, set_cline_refine3D, nstages_refine3D
        use simple_parameters, only: parameters
        class(parameters), intent(inout) :: params
        real :: bounds(0:8)
        integer :: nstages, stage, mode
        write(*,'(A)') 'check_sigma2_schedules'
        params%sigma2_anneal = 'yes'
        do mode=1,3
            params%nstates = 1
            params%nsample = 1000
            if( mode == 2 ) params%nsample = 100000
            if( mode == 3 ) params%nstates = 2
            do nstages=1,8
                nstages_refine3D = nstages
                bounds = 1.
                select case(nstages)
                    case(3)
                        bounds(0:1) = [-1.,1.]
                    case(4)
                        bounds(0:2) = [-1.,0.,1.]
                    case(5)
                        bounds(0:3) = [-1.,-0.5,0.,1.]
                    case(6)
                        bounds(0:4) = [-1.,-2./3.,-1./3.,0.,1.]
                    case(7)
                        bounds(0:5) = [-1.,-0.75,-0.5,-0.25,0.,1.]
                    case(8)
                        bounds(0:6) = [-1.,-0.75,-0.5,-0.25,0.,0.5,1.]
                end select
                do stage=1,nstages
                    call set_cline_refine3D(params, stage, 1, .false.)
                    call assert_string_eq('yes',cline_refine3D%get_carg('sigma2_anneal'),'stage enables scoring annealing')
                    call assert_real(bounds(stage-1),cline_refine3D%get_rarg('sigma2_weight_start'),1.e-6, &
                        &'stage starts at its planned annealing boundary')
                    call assert_real(bounds(stage),cline_refine3D%get_rarg('sigma2_weight_end'),1.e-6, &
                        &'stage ends at its planned annealing boundary')
                enddo
            enddo
        enddo
        params%nstates = 1
        params%nsample = 1000
        nstages_refine3D = 8
        bounds(3:8) = [-1.,0.,0.5,1.,1.,1.]
        do stage=4,8
            call set_cline_refine3D(params,stage,4,.false.)
            call assert_real(bounds(stage-1),cline_refine3D%get_rarg('sigma2_weight_start'),1.e-6,'late entry ramp start')
            call assert_real(bounds(stage),cline_refine3D%get_rarg('sigma2_weight_end'),1.e-6,'late entry ramp end')
        enddo
        call set_cline_refine3D(params,6,6,.false.)
        call assert_real(-1.,cline_refine3D%get_rarg('sigma2_weight_start'),0.,'one pre-final stage starts white')
        call assert_real(1.,cline_refine3D%get_rarg('sigma2_weight_end'),0.,'one pre-final stage completes both ramps')
        call set_cline_refine3D(params,7,7,.false.)
        call assert_real(1.,cline_refine3D%get_rarg('sigma2_weight_start'),0.,'two-stage entry is model based')
        call set_cline_refine3D(params,4,1,.true.)
        call assert_false(cline_refine3D%defined('sigma2_anneal'),'class-average initialization keeps its scoring')
        params%sigma2_anneal = 'no'
        call cline_refine3D%set('sigma2_anneal','yes')
        call cline_refine3D%set('sigma2_weight',-1.)
        call set_cline_refine3D(params,1,1,.false.)
        call assert_false(cline_refine3D%defined('sigma2_anneal'),'disabled schedule removes an inherited switch')
        call assert_false(cline_refine3D%defined('sigma2_weight'),'disabled schedule removes a stale worker weight')
        call assert_false(cline_refine3D%defined('sigma2_weight_start'),'disabled schedule adds no stage bounds')
    end subroutine check_sigma2_schedules

    !> Run emitted early-stop settings through the owner's sampler, over several iterations and the stage handoff.
    subroutine check_early_stage_sticky( params, last_stage, first_stage )
        use simple_solve3D_utils,           only: cline_refine3D, set_cline_refine3D, nstages_refine3D
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4prob3D
        use simple_parameters,              only: parameters
        use simple_builder,                 only: builder
        use simple_string,                  only: string
        class(parameters), intent(in) :: params
        integer,           intent(in) :: last_stage, first_stage
        integer, parameter :: NPTCLS = 16, NITER = 6
        class(builder),    allocatable, target :: build
        class(parameters), allocatable :: sampling
        type(string) :: balance
        integer, allocatable :: inds(:), expected(:)
        integer :: istage, iter, nsamp, nround, i
        write(*,'(A,2I3)') 'check_early_stage_sticky', first_stage, last_stage
        allocate(build, sampling)
        call build%spproj%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        build%spproj_field => build%spproj%os_ptcl3D
        call build%spproj_field%set_all2single('state', 1.)
        sampling%nptcls = NPTCLS
        sampling%l_prob_align_mode = .true.
        sampling%l_update_frac = .true.
        sampling%l_centring = .false.
        nstages_refine3D = last_stage
        nround = 0
        call cline_refine3D%kill
        do istage = first_stage, last_stage
            call set_cline_refine3D(params, istage, first_stage, .false.)
            sampling%startit = cline_refine3D%get_iarg('startit')
            sampling%update_frac = cline_refine3D%get_rarg('update_frac')
            balance = cline_refine3D%get_carg('balance')
            sampling%balance = balance%to_char()
            sampling%l_sticky_sampling = cline_refine3D%defined('sticky_sampling')
            sampling%l_sticky_keep = cline_refine3D%get_carg('sticky_sampling') == 'keep'
            ! no ptcl2D field: the best-first draw falls back to a random one
            sampling%l_sticky_best = cline_refine3D%get_carg('sticky_sampling') == 'best'
            do iter = sampling%startit, sampling%startit + NITER - 1
                sampling%which_iter = iter
                call sample_ptcls4prob3D(sampling, build, nsamp, inds)
                nround = nround + 1
                call assert_int(4, nsamp, 'an early-stop sticky subset keeps its requested size')
                if( .not. allocated(expected) ) expected = inds
                if( size(inds) == size(expected) )then
                    call assert_true(all(inds == expected), 'an early-stop sticky subset survives every iteration and the stage handoff')
                endif
                call assert_int(4, count(build%spproj_field%get_all_asint('updatecnt') > 0), &
                    &'only the original sticky subset acquires update history')
                call assert_int(nround, maxval(build%spproj_field%get_all_asint('sampled')), 'each sticky iteration stamps a new round')
                call assert_true(all([(build%spproj_field%get_updatecnt(expected(i)) == nround, i=1,size(expected))]), &
                    &'each sticky iteration increments every selected particle once')
            end do
            call cline_refine3D%set('endit', sampling%startit + NITER - 1)
        end do
        call build%spproj%os_ptcl3D%kill
        nullify(build%spproj_field)
    end subroutine check_early_stage_sticky

    !> Stages 1-5 through the owner's sampler under one balance mode (none, class or cavg), with the commander's
    !! boundary call before every stage on a project file. Stage 1 draws the sticky subset by the 2D objective,
    !! best first within every unit's quota; stages 2-4 keep it and nothing is released; at stage 5 every update
    !! count restarts and the first free sweep reaches the sticky subset. Without the release (the negative
    !! control) that sweep skips the subset.
    subroutine check_sticky_release( params, mode, l_release )
        use simple_solve3D_utils,           only: cline_refine3D, set_cline_refine3D, nstages_refine3D, &
            &release_sticky_subset
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4prob3D
        use simple_parameters,              only: parameters
        use simple_builder,                 only: builder
        use simple_string,                  only: string
        use simple_type_defs,               only: class_sample
        use simple_class_sample_io,         only: write_class_samples, deallocate_class_samples
        use simple_defs_fname,              only: CLASS_SAMPLING_FILE
        use simple_oris,                    only: class_sample_sweep
        class(parameters), intent(in) :: params
        character(len=*),  intent(in) :: mode
        logical,           intent(in) :: l_release
        integer, parameter :: NPTCLS = 16, NITER = 3, NSTICKY = 4
        class(builder),    allocatable, target :: build
        class(parameters), allocatable :: sampling
        type(class_sample), allocatable :: units(:)
        type(string) :: projfile
        integer, allocatable :: inds(:), sticky(:), best(:)
        logical :: drawn(NPTCLS), l_sticky(NPTCLS), l_sticky_prev
        integer :: istage, iter, nsamp, nits, nupdated, i
        write(*,'(A,A,L2)') 'check_sticky_release ', mode, l_release
        allocate(build, sampling)
        call build%spproj%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        build%spproj_field => build%spproj%os_ptcl3D
        call build%spproj_field%set_all2single('state', 1.)
        ! the 2D objective of rows 1-16: 5 10 15 3 8 13 1 6 | 31 36 24 29 34 22 27 32
        call build%spproj%os_ptcl2D%new(NPTCLS, is_ptcl=.true.)
        call build%spproj%os_ptcl2D%set_all2single('state', 1.)
        do i = 1, NPTCLS
            call build%spproj%os_ptcl2D%set(i, 'corr', real(mod(5*i, 17)) + merge(20., 0., i > 8))
        end do
        ! the units and the best rows of every unit's quota (4 of 16 drawn)
        select case(mode)
            case('none')
                best = [9,10,13,16]                         ! the four best rows
            case('class')
                call make_units([8,8], [0,0])               ! 2 per unit
                best = [3,6,10,13]
            case('cavg')
                call make_units([4,4,8], [1,2,2])           ! group quotas 2 and 2: 2 from rows 1-4, 1 from each other unit
                best = [2,3,6,10]
        end select
        projfile = 'sticky_release.simple'
        call build%spproj%update_projinfo(projfile)
        sampling%nptcls = NPTCLS
        sampling%l_prob_align_mode = .true.
        sampling%l_update_frac = .true.
        sampling%l_centring = .false.
        nstages_refine3D = NSTICKY + 1
        drawn = .false.
        l_sticky_prev = .false.
        call cline_refine3D%kill
        do istage = 1, NSTICKY + 1
            call set_cline_refine3D(params, istage, 1, .false.)
            sampling%startit = cline_refine3D%get_iarg('startit')
            sampling%update_frac = cline_refine3D%get_rarg('update_frac')
            sampling%balance = mode
            sampling%l_sticky_sampling = cline_refine3D%defined('sticky_sampling')
            sampling%l_sticky_keep = .false.
            sampling%l_sticky_best = .false.
            if( sampling%l_sticky_sampling )then
                sampling%l_sticky_keep = cline_refine3D%get_carg('sticky_sampling') == 'keep'
                sampling%l_sticky_best = cline_refine3D%get_carg('sticky_sampling') == 'best'
            endif
            if( l_release )then
                ! the commander's stage boundary, on the project file the stage starts from
                call build%spproj%write(projfile)
                nupdated = count(build%spproj_field%get_all_asint('updatecnt') > 0)
                call release_sticky_subset(build%spproj, projfile, l_sticky_prev)
                if( istage <= NSTICKY )then
                    call assert_int(nupdated, count(build%spproj_field%get_all_asint('updatecnt') > 0), &
                        &'no release while the stages are sticky')
                else
                    call assert_int(0, count(build%spproj_field%get_all_asint('updatecnt') > 0), &
                        &'the first free stage restarts every update count')
                    call build%spproj%read_segment('ptcl3D', projfile)
                    call assert_int(0, count(build%spproj_field%get_all_asint('updatecnt') > 0), &
                        &'the released counts are in the project file')
                endif
            endif
            nits = NITER
            if( istage == 1 ) call assert_true(sampling%l_sticky_best, 'stage 1 draws the sticky subset best first')
            if( istage > NSTICKY )then
                call assert_false(sampling%l_sticky_sampling, 'stage 5 samples freely')
                ! one free sweep: every particle once, if the draws start from a common count
                if( allocated(units) )then
                    nits = class_sample_sweep(units, nint(sampling%update_frac * real(NPTCLS)))
                else
                    nits = NPTCLS / nint(sampling%update_frac * real(NPTCLS))
                endif
            endif
            do iter = sampling%startit, sampling%startit + nits - 1
                sampling%which_iter = iter
                call sample_ptcls4prob3D(sampling, build, nsamp, inds)
                if( istage <= NSTICKY )then
                    sticky = inds
                    call assert_int(size(best), size(sticky), 'the sticky subset keeps its size')
                    if( size(sticky) == size(best) )then
                        call assert_true(all(sticky == best), 'the sticky subset is every unit''s best rows by the 2D objective')
                    endif
                else
                    drawn(inds) = .true.
                endif
            end do
            call cline_refine3D%set('endit', sampling%startit + nits - 1)
        end do
        l_sticky = [(any(sticky == i), i = 1, NPTCLS)]
        if( l_release )then
            call assert_true(all(drawn), 'after the release the first free sweep reaches every particle')
        else
            call assert_false(any(drawn .and. l_sticky), 'without the release the first free sweep skips the sticky subset')
            call assert_true(all(drawn .neqv. l_sticky), 'without the release the first free sweep takes every other particle')
        endif
        if( allocated(units) ) call deallocate_class_samples(units)
        call build%spproj%os_ptcl2D%kill
        call build%spproj%os_ptcl3D%kill
        nullify(build%spproj_field)

    contains

        !> consecutive rows per unit, written as the run's sampling units
        subroutine make_units( pops, groups )
            integer, intent(in) :: pops(:), groups(:)
            integer :: k, j, first
            allocate(units(size(pops)))
            first = 1
            do k = 1, size(pops)
                units(k)%clsind  = k
                units(k)%pop     = pops(k)
                units(k)%group   = groups(k)
                units(k)%nsample = 0
                allocate(units(k)%pinds(pops(k)), source=[(j, j=first,first+pops(k)-1)])
                first = first + pops(k)
            end do
            call write_class_samples(units, string(CLASS_SAMPLING_FILE))
        end subroutine make_units

    end subroutine check_sticky_release

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
