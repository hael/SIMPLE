!@descr: startup handoff and first-update sampling contract
module simple_refine3D_startup_tester
implicit none
private
public :: run_all_refine3D_startup_tests

contains

    subroutine run_all_refine3D_startup_tests()
        use simple_refine3D_startup, only: refine3D_startup, STARTUP_ORDINARY, STARTUP_REUSE_CHAIN, STARTUP_SEED_CHAIN
        use simple_cmdline,          only: cmdline
        use simple_parameters,       only: parameters
        use simple_test_utils,       only: assert_true, assert_false, assert_real, assert_int
        class(cmdline), allocatable    :: cline
        class(parameters), allocatable :: params
        type(refine3D_startup)         :: producer, consumer
        integer                        :: mode
        logical                        :: full_update
        allocate(cline, params)
        params%startit = 3
        params%which_iter = 3
        params%maxits = 12
        do mode = STARTUP_ORDINARY, STARTUP_SEED_CHAIN
            params%update_frac = 0.25
            params%l_update_frac = .true.
            params%l_trail_rec = .true.
            call cline%set('update_frac', 0.25)
            call producer%new(mode)
            call producer%write_cline(cline)
            call consumer%read_cline(cline)
            call assert_false(cline%defined('trail_seed'), 'startup handoff is consumed at command entry')
            call assert_true(consumer%preserves_history() .eqv. (mode == STARTUP_REUSE_CHAIN), &
                &'only a reused chain preserves the previous update history')
            full_update = consumer%begin_iteration(params, cline)
            call assert_true(full_update .eqv. (mode == STARTUP_SEED_CHAIN), &
                &'only a missing chain requires a full first update')
            if( mode == STARTUP_SEED_CHAIN )then
                call assert_real(1., params%update_frac, 0., 'seed iteration updates every active row')
                call assert_real(1., cline%get_rarg('update_frac'), 0., 'seed command uses the full sample')
                call assert_false(params%l_update_frac, 'seed iteration disables fractional matching')
                call assert_false(params%l_trail_rec, 'seed iteration does not blend a previous chain')
                call assert_true(cline%defined('trail_seed'), 'seed iteration asks reconstruction to publish its chain')
            else
                call assert_real(0.25, params%update_frac, 0., 'ordinary/reuse iteration keeps the requested sample')
                call assert_true(params%l_update_frac, 'ordinary/reuse iteration remains fractional')
                call assert_true(params%l_trail_rec, 'ordinary/reuse iteration retains trailing reconstruction')
            endif
            call consumer%finish_iteration(params, cline)
            call assert_real(0.25, params%update_frac, 0., 'completed iteration restores the requested fraction')
            call assert_real(0.25, cline%get_rarg('update_frac'), 0., 'next command uses the requested fraction')
            call assert_true(params%l_update_frac, 'completed iteration restores fractional matching')
            call assert_true(params%l_trail_rec, 'completed iteration restores trailing reconstruction')
            call assert_false(cline%defined('trail_seed'), 'completed iteration leaves no reconstruction handshake')
            call assert_int(12, params%maxits, 'startup transition preserves the iteration budget')
            params%which_iter = 4
            full_update = consumer%begin_iteration(params, cline)
            call assert_false(full_update, 'the next iteration uses the requested sample')
            params%which_iter = 3
        enddo
        ! Full-update workflows need no fractional seed override.
        call consumer%new(STARTUP_SEED_CHAIN)
        params%l_update_frac = .false.
        params%update_frac = 1.
        full_update = consumer%begin_iteration(params, cline)
        call assert_false(full_update, 'a full-update workflow needs no temporary sampling override')
        call consumer%kill
        call producer%kill
        call cline%kill
    end subroutine run_all_refine3D_startup_tests

end module simple_refine3D_startup_tester
