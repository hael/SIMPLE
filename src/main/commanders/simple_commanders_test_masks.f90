!@descr: for all masks tests
module simple_commanders_test_masks
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_bounds_from_mask3D_test
  contains
    procedure :: execute      => exec_test_bounds_from_mask3D_test
end type commander_test_bounds_from_mask3D_test

type, extends(commander_base) :: commander_test_graphene_mask
  contains
    procedure :: execute      => exec_test_graphene_mask
end type commander_test_graphene_mask

type, extends(commander_base) :: commander_test_image_bin
  contains
    procedure :: execute      => exec_test_image_bin
end type commander_test_image_bin

type, extends(commander_base) :: commander_test_mask
  contains
    procedure :: execute      => exec_test_mask
end type commander_test_mask

type, extends(commander_base) :: commander_test_msk_routines
  contains
    procedure :: execute      => exec_test_msk_routines
end type commander_test_msk_routines

type, extends(commander_base) :: commander_test_nano_mask
  contains
    procedure :: execute      => exec_test_nano_mask
end type commander_test_nano_mask

type, extends(commander_base) :: commander_test_otsu_test
  contains
    procedure :: execute      => exec_test_otsu_test
end type commander_test_otsu_test

type, extends(commander_base) :: commander_test_ptcl_center
  contains
    procedure :: execute      => exec_test_ptcl_center
end type commander_test_ptcl_center

contains

subroutine exec_test_bounds_from_mask3D_test( self, cline )
    class(commander_test_bounds_from_mask3D_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_BOUNDS_FROM_MASK3D_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_bounds_from_mask3D_test

subroutine exec_test_graphene_mask( self, cline )
    class(commander_test_graphene_mask),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_GRAPHENE_MASK_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_graphene_mask

subroutine exec_test_image_bin( self, cline )
    class(commander_test_image_bin),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_IMAGE_BIN_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_image_bin

subroutine exec_test_mask( self, cline )
    class(commander_test_mask),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_MASK_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mask

subroutine exec_test_msk_routines( self, cline )
    class(commander_test_msk_routines),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_MSK_ROUTINES_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_msk_routines

subroutine exec_test_nano_mask( self, cline )
    class(commander_test_nano_mask),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_NANO_MASK_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_nano_mask

subroutine exec_test_otsu_test( self, cline )
    class(commander_test_otsu_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_OTSU_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_otsu_test

subroutine exec_test_ptcl_center( self, cline )
    class(commander_test_ptcl_center),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_PTCL_CENTER_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ptcl_center

end module simple_commanders_test_masks
