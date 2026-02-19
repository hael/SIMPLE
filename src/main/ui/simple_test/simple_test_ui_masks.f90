!@descr: module defining the user interfaces for masks test programs in the simple_test_exec suite
module simple_test_ui_masks
use simple_ui_modules
implicit none

type(ui_program), target :: bounds_from_mask3D_test
type(ui_program), target :: graphene_mask
type(ui_program), target :: image_bin
type(ui_program), target :: mask
type(ui_program), target :: msk_routines
type(ui_program), target :: nano_mask
type(ui_program), target :: otsu_test
type(ui_program), target :: ptcl_center

contains

    subroutine construct_test_masks_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_bounds_from_mask3D_test(tsttab)
        call new_graphene_mask(tsttab)
        call new_image_bin(tsttab)
        call new_mask(tsttab)
        call new_msk_routines(tsttab)
        call new_nano_mask(tsttab)
        call new_otsu_test(tsttab)
        call new_ptcl_center(tsttab)
    end subroutine construct_test_masks_programs

    subroutine print_test_masks_programs( logfhandle )
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('MASKS:', C_UNDERLINED)
        write(logfhandle,'(A)') bounds_from_mask3D_test%name%to_char()
        write(logfhandle,'(A)') graphene_mask%name%to_char()
        write(logfhandle,'(A)') image_bin%name%to_char()
        write(logfhandle,'(A)') mask%name%to_char()
        write(logfhandle,'(A)') msk_routines%name%to_char()
        write(logfhandle,'(A)') nano_mask%name%to_char()
        write(logfhandle,'(A)') otsu_test%name%to_char()
        write(logfhandle,'(A)') ptcl_center%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_masks_programs

    subroutine new_bounds_from_mask3D_test( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call bounds_from_mask3D_test%new(&
        &'bounds_from_mask3D_test',&           ! name
        &'bounds_from_mask3D_test ',&          ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call bounds_from_mask3D_test%add_input(UI_IO, )
        ! parameter input/output
        !call bounds_from_mask3D_test%add_input(UI_IMG, )
        ! alternative inputs
        !call bounds_from_mask3D_test%add_input(UI_PARM, )
        ! search controls
        !call bounds_from_mask3D_test%add_input(UI_SRCH, )
        ! filter controls
        !call bounds_from_mask3D_test%add_input(UI_FILT, )
        ! mask controls
        !call bounds_from_mask3D_test%add_input(UI_MASK, )
        ! computer controls
        !call bounds_from_mask3D_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('bounds_from_mask3D_test', bounds_from_mask3D_test, tsttab)
    end subroutine new_bounds_from_mask3D_test

    subroutine new_graphene_mask( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call graphene_mask%new(&
        &'graphene_mask',&                     ! name
        &'graphene_mask ',&                    ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call graphene_mask%add_input(UI_IO, )
        ! parameter input/output
        !call graphene_mask%add_input(UI_IMG, )
        ! alternative inputs
        !call graphene_mask%add_input(UI_PARM, )
        ! search controls
        !call graphene_mask%add_input(UI_SRCH, )
        ! filter controls
        !call graphene_mask%add_input(UI_FILT, )
        ! mask controls
        !call graphene_mask%add_input(UI_MASK, )
        ! computer controls
        !call graphene_mask%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('graphene_mask', graphene_mask, tsttab)
    end subroutine new_graphene_mask

    subroutine new_image_bin( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call image_bin%new(&
        &'image_bin',&                         ! name
        &'image_bin ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call image_bin%add_input(UI_IO, )
        ! parameter input/output
        !call image_bin%add_input(UI_IMG, )
        ! alternative inputs
        !call image_bin%add_input(UI_PARM, )
        ! search controls
        !call image_bin%add_input(UI_SRCH, )
        ! filter controls
        !call image_bin%add_input(UI_FILT, )
        ! mask controls
        !call image_bin%add_input(UI_MASK, )
        ! computer controls
        !call image_bin%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('image_bin', image_bin, tsttab)
    end subroutine new_image_bin

    subroutine new_mask( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call mask%new(&
        &'mask',&                         ! name
        &'mask ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call mask%add_input(UI_IO, )
        ! parameter input/output
        !call mask%add_input(UI_IMG, )
        ! alternative inputs
        !call mask%add_input(UI_PARM, )
        ! search controls
        !call mask%add_input(UI_SRCH, )
        ! filter controls
        !call mask%add_input(UI_FILT, )
        ! mask controls
        !call mask%add_input(UI_MASK, )
        ! computer controls
        !call mask%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('mask', mask, tsttab)
    end subroutine new_mask

    subroutine new_msk_routines( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call msk_routines%new(&
        &'msk_routines',&                      ! name
        &'msk_routines ',&                     ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call msk_routines%add_input(UI_IO, )
        ! parameter input/output
        !call msk_routines%add_input(UI_IMG, )
        ! alternative inputs
        !call msk_routines%add_input(UI_PARM, )
        ! search controls
        !call msk_routines%add_input(UI_SRCH, )
        ! filter controls
        !call msk_routines%add_input(UI_FILT, )
        ! mask controls
        !call msk_routines%add_input(UI_MASK, )
        ! computer controls
        !call msk_routines%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('msk_routines', msk_routines, tsttab)
    end subroutine new_msk_routines

    subroutine new_nano_mask( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call nano_mask%new(&
        &'nano_mask',&                         ! name
        &'nano_mask ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call nano_mask%add_input(UI_IO, )
        ! parameter input/output
        !call nano_mask%add_input(UI_IMG, )
        ! alternative inputs
        !call nano_mask%add_input(UI_PARM, )
        ! search controls
        !call nano_mask%add_input(UI_SRCH, )
        ! filter controls
        !call nano_mask%add_input(UI_FILT, )
        ! mask controls
        !call nano_mask%add_input(UI_MASK, )
        ! computer controls
        !call nano_mask%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('nano_mask', nano_mask, tsttab)
    end subroutine new_nano_mask

    subroutine new_otsu_test( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call otsu_test%new(&
        &'otsu_test',&                         ! name
        &'otsu_test ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call otsu_test%add_input(UI_IO, )
        ! parameter input/output
        !call otsu_test%add_input(UI_IMG, )
        ! alternative inputs
        !call otsu_test%add_input(UI_PARM, )
        ! search controls
        !call otsu_test%add_input(UI_SRCH, )
        ! filter controls
        !call otsu_test%add_input(UI_FILT, )
        ! mask controls
        !call otsu_test%add_input(UI_MASK, )
        ! computer controls
        !call otsu_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('otsu_test', otsu_test, tsttab)
    end subroutine new_otsu_test

    subroutine new_ptcl_center( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call ptcl_center%new(&
        &'ptcl_center',&                       ! name
        &'ptcl_center ',&                      ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call ptcl_center%add_input(UI_IO, )
        ! parameter input/output
        !call ptcl_center%add_input(UI_IMG, )
        ! alternative inputs
        !call ptcl_center%add_input(UI_PARM, )
        ! search controls
        !call ptcl_center%add_input(UI_SRCH, )
        ! filter controls
        !call ptcl_center%add_input(UI_FILT, )
        ! mask controls
        !call ptcl_center%add_input(UI_MASK, )
        ! computer controls
        !call ptcl_center%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('ptcl_center', ptcl_center, tsttab)
    end subroutine new_ptcl_center

end module simple_test_ui_masks
