!@descr: execution of test masks processing commanders
module simple_test_exec_masks
use simple_cmdline,               only: cmdline
use simple_commanders_test_masks, only: commander_test_bounds_from_mask3D_test, &
                                        commander_test_graphene_mask, commander_test_image_bin, &
                                        commander_test_mask, commander_test_msk_routines, &
                                        commander_test_nano_mask, commander_test_otsu_test, &
                                        commander_test_ptcl_center
implicit none

public :: exec_masks_commander
private

type(commander_test_bounds_from_mask3D_test) :: xbounds_from_mask3D_test
type(commander_test_graphene_mask)           :: xgraphene_mask
type(commander_test_image_bin)               :: ximage_bin
type(commander_test_mask)                    :: xmask
type(commander_test_msk_routines)            :: xmsk_routines
type(commander_test_nano_mask)               :: xnano_mask
type(commander_test_otsu_test)               :: xotsu_test
type(commander_test_ptcl_center)             :: xptcl_center

contains

    subroutine exec_masks_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'bounds_from_mask3D_test' )
                call xbounds_from_mask3D_test%execute(cline)
            case( 'graphene_mask' )
                call xgraphene_mask%execute(cline)
            case( 'image_bin' )
                call ximage_bin%execute(cline)
            case( 'mask' )
                call xmask%execute(cline)
            case( 'msk_routines' )
                call xmsk_routines%execute(cline)
            case( 'nano_mask' )
                call xnano_mask%execute(cline)
            case( 'otsu_test' )
                call xotsu_test%execute(cline)
            case( 'ptcl_center' )
                call xptcl_center%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_masks_commander

end module simple_test_exec_masks
