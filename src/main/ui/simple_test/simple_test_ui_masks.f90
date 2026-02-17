!@ descr: module defining the user interfaces for masks programs in the simple_test_exec suite
module simple_test_ui_masks
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_mask
type(ui_program), target :: simple_test_msk_routines
type(ui_program), target :: simple_test_otsu
type(ui_program), target :: simple_test_bounds_from_mask3D
type(ui_program), target :: simple_test_graphene_mask
type(ui_program), target :: simple_test_nano_mask
type(ui_program), target :: simple_test_ptcl_center
type(ui_program), target :: simple_test_image_bin

contains

    subroutine construct_masks_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_mask(prgtab)
        call new_simple_test_msk_routines(prgtab)
        call new_simple_test_otsu(prgtab)
        call new_simple_test_bounds_from_mask3D(prgtab)
        call new_simple_test_graphene_mask(prgtab)
        call new_simple_test_nano_mask(prgtab)
        call new_simple_test_ptcl_center(prgtab)
        call new_simple_test_image_bin(prgtab)
    end subroutine construct_masks_programs

    subroutine print_masks_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('MASKS:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_mask%name%to_char()
        write(logfhandle,'(A)') simple_test_msk_routines%name%to_char()
        write(logfhandle,'(A)') simple_test_otsu%name%to_char()
        write(logfhandle,'(A)') simple_test_bounds_from_mask3D%name%to_char()
        write(logfhandle,'(A)') simple_test_graphene_mask%name%to_char()
        write(logfhandle,'(A)') simple_test_nano_mask%name%to_char()
        write(logfhandle,'(A)') simple_test_ptcl_center%name%to_char()
        write(logfhandle,'(A)') simple_test_image_bin%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_masks_programs

    subroutine new_simple_test_mask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_mask', simple_test_mask, prgtab)
    end subroutine new_simple_test_mask

    subroutine new_simple_test_msk_routines( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_msk_routines', simple_test_msk_routines, prgtab)
    end subroutine new_simple_test_msk_routines

    subroutine new_simple_test_otsu( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_otsu', simple_test_otsu, prgtab)
    end subroutine new_simple_test_otsu

    subroutine new_simple_test_bounds_from_mask3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_bounds_from_mask3D', simple_test_bounds_from_mask3D, prgtab)
    end subroutine new_simple_test_bounds_from_mask3D

    subroutine new_simple_test_graphene_mask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_graphene_mask', simple_test_graphene_mask, prgtab)
    end subroutine new_simple_test_graphene_mask

    subroutine new_simple_test_nano_mask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_nano_mask', simple_test_nano_mask, prgtab)
    end subroutine new_simple_test_nano_mask

    subroutine new_simple_test_ptcl_center( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_ptcl_center', simple_test_ptcl_center, prgtab)
    end subroutine new_simple_test_ptcl_center

    subroutine new_simple_test_image_bin( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_image_bin', simple_test_image_bin, prgtab)
    end subroutine new_simple_test_image_bin

end module simple_test_ui_masks
