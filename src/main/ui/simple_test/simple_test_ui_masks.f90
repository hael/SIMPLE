!@descr: module defining the user interfaces for masks test programs in the simple_test_exec suite
module simple_test_ui_masks
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('masks', 'Masks', 60)
type(ui_program), target :: nano_mask
type(ui_program), target :: score_volume_shape

contains

    subroutine construct_test_masks_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_nano_mask(tsttab)
        call new_score_volume_shape(tsttab)
    end subroutine construct_test_masks_programs

    subroutine new_nano_mask( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call nano_mask%new(&
        &'nano_mask',&                         ! name
        &'nano_mask ',&                        ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call nano_mask%add_input(UI_IMG, stk, required_override=.false.)
        ! parameter input/output
        call nano_mask%add_input(UI_PARM, smpd, required_override=.false.)
        ! alternative inputs
        ! <no additional inputs>
        !call nano_mask%add_input(UI_PARM, )
        ! search controls
        !call nano_mask%add_input(UI_SRCH, )
        ! filter controls
        !call nano_mask%add_input(UI_FILT, )
        ! mask controls
        call nano_mask%add_input(UI_MASK, mskdiam, required_override=.false.)
        ! computer controls
        !call nano_mask%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('nano_mask', nano_mask, tsttab, UI_CATEGORY)
    end subroutine new_nano_mask

    subroutine new_score_volume_shape( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call score_volume_shape%new(&
        &'score_volume_shape',&
        &'Score volume shape descriptors',&
        &'is a test program for evaluating the volume shape descriptors',&
        &'simple_test_exec',&
        &.false.)
        call score_volume_shape%add_input(UI_IMG, 'vol1', 'file', 'Volume', &
        &'Volume to score', 'input volume e.g. vol.mrc', .true., '')
        call score_volume_shape%add_input(UI_PARM, smpd, required_override=.false.)
        call add_ui_program('score_volume_shape', score_volume_shape, tsttab, UI_CATEGORY)
    end subroutine new_score_volume_shape

end module simple_test_ui_masks
