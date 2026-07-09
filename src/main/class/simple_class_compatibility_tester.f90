!@descr: unit test routines for the class_compatibility support model
module simple_class_compatibility_tester
use simple_test_utils
use simple_core_module_api
use simple_class_compatibility, only: class_compatibility, support_model_metrics
use simple_sp_project,          only: sp_project
use simple_imgarr_utils,        only: write_imgarr, dealloc_imgarr
use simple_image,               only: image
use simple_string,              only: string
implicit none
private
public :: run_all_class_compatibility_tests

real, parameter :: AXIS_TOL = 1.0e-4

contains

    subroutine run_all_class_compatibility_tests()
        write(*,'(A)') '**** running all class_compatibility tests ****'
        call test_defaults_and_lifecycle()
        call test_train_1_too_small_set_no_fit()
        call test_train_1_from_project_sets_valid_fit()
        call test_train_2_retrain_sets_delta_and_convergence()
        call test_infer_without_valid_model_preserves_selection()
        call test_infer_rejects_size_outliers()
        call test_kill_support_model_resets_after_fit()
    end subroutine run_all_class_compatibility_tests

    subroutine test_defaults_and_lifecycle()
        type(class_compatibility) :: comp
        type(support_model_metrics) :: metrics

        write(*,'(A)') 'test_defaults_and_lifecycle'
        call comp%new()
        call comp%get_support_model_metrics(metrics)
        call assert_false(metrics%valid,       'default model invalid before training')
        call assert_false(metrics%delta_valid, 'default delta_valid is false')
        call assert_false(metrics%converged,   'default converged is false')
        call assert_real(0.0, metrics%axis_c,  AXIS_TOL, 'default axis_c')
        call assert_real(0.0, metrics%axis_b,  AXIS_TOL, 'default axis_b')
        call assert_real(0.0, metrics%axis_a,  AXIS_TOL, 'default axis_a')
        call assert_false(comp%converged(),    'converged() reports false before training')

        call comp%kill_support_model()
        call comp%kill()
    end subroutine test_defaults_and_lifecycle

    subroutine test_train_1_from_project_sets_valid_fit()
        type(class_compatibility) :: comp
        type(sp_project)          :: proj
        type(support_model_metrics) :: metrics
        type(string)              :: stkname
        integer, allocatable      :: states(:)

        write(*,'(A)') 'test_train_1_from_project_sets_valid_fit'

        stkname = string('class_compat_train1.mrcs')
        call make_blob_stack(stkname, [8, 9, 8, 9, 10, 9])

        allocate(states(6), source=1)
        call init_project_for_cavgs(proj, stkname, states)

        call comp%new()
        call comp%train(proj)
        call comp%get_support_model_metrics(metrics)

        call assert_true(metrics%valid,            'train_1 produces a valid fit')
        call assert_true(metrics%axis_a >= metrics%axis_b, 'axis_a >= axis_b after train_1')
        call assert_true(metrics%axis_b >= metrics%axis_c, 'axis_b >= axis_c after train_1')
        call assert_false(metrics%delta_valid,     'first fit has no previous delta')
        call assert_false(metrics%converged,       'first fit is not converged')

        call comp%kill()
        call proj%kill()
        if( allocated(states) ) deallocate(states)
        call del_file(stkname)
    end subroutine test_train_1_from_project_sets_valid_fit

    subroutine test_train_1_too_small_set_no_fit()
        type(class_compatibility) :: comp
        type(sp_project)          :: proj
        type(support_model_metrics) :: metrics
        type(string)              :: stkname
        integer, allocatable      :: states(:)

        write(*,'(A)') 'test_train_1_too_small_set_no_fit'

        stkname = string('class_compat_train1_small.mrcs')
        call make_blob_stack(stkname, [8, 9, 10])

        allocate(states(3), source=0)
        states(1:2) = 1
        call init_project_for_cavgs(proj, stkname, states)

        call comp%new()
        call comp%train(proj)
        call comp%get_support_model_metrics(metrics)

        call assert_false(metrics%valid,      'train_1 with <3 selected classes stays invalid')
        call assert_false(metrics%delta_valid,'train_1 small set has no valid deltas')
        call assert_false(metrics%converged,  'train_1 small set is not converged')

        call comp%kill()
        call proj%kill()
        if( allocated(states) ) deallocate(states)
        call del_file(stkname)
    end subroutine test_train_1_too_small_set_no_fit

    subroutine test_train_2_retrain_sets_delta_and_convergence()
        type(class_compatibility) :: comp
        type(support_model_metrics) :: metrics
        type(string)              :: stkname

        write(*,'(A)') 'test_train_2_retrain_sets_delta_and_convergence'

        stkname = string('class_compat_train2.mrcs')
        call make_blob_stack(stkname, [9, 9, 10, 9, 10, 9])

        call comp%new()
        call comp%train(stkname)
        call comp%get_support_model_metrics(metrics)
        call assert_true(metrics%valid,        'first train_2 fit is valid')
        call assert_false(metrics%delta_valid, 'first train_2 has no valid delta')
        call assert_false(metrics%converged,   'first train_2 is not converged')

        call comp%train(stkname)
        call comp%get_support_model_metrics(metrics)
        call assert_true(metrics%valid,        'second train_2 fit remains valid')
        call assert_true(metrics%delta_valid,  'second train_2 has valid deltas')
        call assert_true(metrics%converged,    'second identical train_2 converges')
        call assert_real(0.0, metrics%delta_c, AXIS_TOL, 'delta_c near zero on identical retrain')
        call assert_real(0.0, metrics%delta_b, AXIS_TOL, 'delta_b near zero on identical retrain')
        call assert_real(0.0, metrics%delta_a, AXIS_TOL, 'delta_a near zero on identical retrain')

        call comp%kill()
        call del_file(stkname)
    end subroutine test_train_2_retrain_sets_delta_and_convergence

    subroutine test_infer_rejects_size_outliers()
        type(class_compatibility) :: comp
        type(sp_project)          :: proj
        type(string)              :: train_stk, infer_stk
        integer, allocatable      :: states(:), out_states(:)

        write(*,'(A)') 'test_infer_rejects_size_outliers'

        train_stk = string('class_compat_infer_train.mrcs')
        infer_stk = string('class_compat_infer_eval.mrcs')

        ! Train on compact, consistent objects.
        call make_blob_stack(train_stk, [9, 9, 10, 9, 10, 9])
        ! Evaluate one in-family object plus two obvious outliers.
        call make_blob_stack(infer_stk, [9, 3, 20])

        allocate(states(3), source=1)
        call init_project_for_cavgs(proj, infer_stk, states)

        call comp%new()
        call comp%train(train_stk)
        call comp%infer(proj)

        out_states = proj%os_cls2D%get_all_asint('state')
        call assert_true(out_states(1) > 0,          'in-family cavg remains selected')
        call assert_true(count(out_states <= 0) >= 1,'at least one outlier is rejected')

        call comp%kill()
        call proj%kill()
        if( allocated(states) ) deallocate(states)
        if( allocated(out_states) ) deallocate(out_states)
        call del_file(train_stk)
        call del_file(infer_stk)
    end subroutine test_infer_rejects_size_outliers

    subroutine test_infer_without_valid_model_preserves_selection()
        type(class_compatibility) :: comp
        type(sp_project)          :: proj
        type(string)              :: infer_stk
        integer, allocatable      :: states(:), out_states(:)

        write(*,'(A)') 'test_infer_without_valid_model_preserves_selection'

        infer_stk = string('class_compat_infer_untrained.mrcs')
        call make_blob_stack(infer_stk, [9, 3, 20])

        allocate(states(3), source=1)
        call init_project_for_cavgs(proj, infer_stk, states)

        call comp%new()
        call comp%infer(proj)

        out_states = proj%os_cls2D%get_all_asint('state')
        call assert_int(3, count(out_states > 0), 'infer without valid model preserves selection')

        call comp%kill()
        call proj%kill()
        if( allocated(states) ) deallocate(states)
        if( allocated(out_states) ) deallocate(out_states)
        call del_file(infer_stk)
    end subroutine test_infer_without_valid_model_preserves_selection

    subroutine test_kill_support_model_resets_after_fit()
        type(class_compatibility) :: comp
        type(support_model_metrics) :: metrics
        type(string)              :: stkname

        write(*,'(A)') 'test_kill_support_model_resets_after_fit'

        stkname = string('class_compat_reset.mrcs')
        call make_blob_stack(stkname, [9, 9, 10, 9, 10, 9])

        call comp%new()
        call comp%train(stkname)
        call comp%get_support_model_metrics(metrics)
        call assert_true(metrics%valid, 'fit is valid before kill_support_model')

        call comp%kill_support_model()
        call comp%get_support_model_metrics(metrics)
        call assert_false(metrics%valid,       'kill_support_model resets valid')
        call assert_false(metrics%delta_valid, 'kill_support_model resets delta_valid')
        call assert_false(metrics%converged,   'kill_support_model resets converged')
        call assert_real(0.0, metrics%axis_c,  AXIS_TOL, 'kill_support_model resets axis_c')
        call assert_real(0.0, metrics%axis_b,  AXIS_TOL, 'kill_support_model resets axis_b')
        call assert_real(0.0, metrics%axis_a,  AXIS_TOL, 'kill_support_model resets axis_a')

        call comp%kill()
        call del_file(stkname)
    end subroutine test_kill_support_model_resets_after_fit

    subroutine init_project_for_cavgs(proj, stkname, states)
        type(sp_project), intent(inout) :: proj
        type(string),     intent(in)    :: stkname
        integer,          intent(in)    :: states(:)
        integer :: i, ncls

        ncls = size(states)
        call proj%kill()
        call proj%os_cls2D%new(ncls, is_ptcl=.false.)
        do i = 1, ncls
            call proj%os_cls2D%set_class(i, i)
            call proj%os_cls2D%set_state(i, states(i))
        end do

        call proj%os_out%new(1, is_ptcl=.false.)
        call proj%os_out%set(1, 'imgkind', 'cavg')
        call proj%os_out%set(1, 'stk',     stkname)
        call proj%os_out%set(1, 'nptcls',  ncls)
        call proj%os_out%set(1, 'fromp',   1)
        call proj%os_out%set(1, 'top',     ncls)
        call proj%os_out%set(1, 'smpd',    1.0)
        call proj%os_out%set(1, 'box',     64)
    end subroutine init_project_for_cavgs

    subroutine make_blob_stack(stkname, half_widths)
        type(string), intent(in) :: stkname
        integer,      intent(in) :: half_widths(:)
        type(image), allocatable :: imgs(:)
        integer :: i, n

        n = size(half_widths)
        allocate(imgs(n))
        do i = 1, n
            call imgs(i)%new([64, 64, 1], 1.0)
            call imgs(i)%zero()
            call paint_square_blob(imgs(i), half_widths(i), 1.0)
        end do
        call write_imgarr(imgs, stkname)
        call dealloc_imgarr(imgs)
    end subroutine make_blob_stack

    subroutine paint_square_blob(img, half_width, val)
        type(image), intent(inout) :: img
        integer,     intent(in)    :: half_width
        real,        intent(in)    :: val
        integer :: i, j, cx, cy
        integer :: ldim(3)

        ldim = img%get_ldim()
        cx = (ldim(1) + 1) / 2
        cy = (ldim(2) + 1) / 2
        do j = max(1, cy - half_width), min(ldim(2), cy + half_width)
            do i = max(1, cx - half_width), min(ldim(1), cx + half_width)
                call img%set([i, j, 1], val)
            end do
        end do
    end subroutine paint_square_blob

end module simple_class_compatibility_tester
