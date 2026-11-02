!@descr: unit tests for image serialisation (serialize), the pixel vectors of stack operations and the nanoparticle tools
! Moved from the serialize test program by the utils review (plan, section 9.7); the program
! displayed a masked square and checked nothing. serialize packs the pixels inside a logical mask
! in column-major order; stack operations and the nanoparticle tools use it.
module simple_image_serialize_tester
use simple_core_module_api
use simple_image, only: image
use simple_test_utils
implicit none
private
public :: run_all_image_serialize_tests

integer, parameter :: BOX  = 64
real,    parameter :: SMPD = 1.0

contains

    subroutine run_all_image_serialize_tests()
        write(*,'(A)') '**** running all image serialisation tests ****'
        call test_masked_serialize()
    end subroutine run_all_image_serialize_tests

    !> a 64^2 image with a different value in every pixel and a disc mask of radius 20: the vector
    !! holds one value per masked pixel, the masked pixels in column-major order
    subroutine test_masked_serialize()
        type(image)          :: img, img_msk
        real,    allocatable :: pcavec(:), orig(:,:,:)
        logical, allocatable :: l_mask(:,:,:)
        write(*,'(A)') 'test_masked_serialize'
        call make_distinct_image(img, orig)
        call img_msk%disc([BOX,BOX,1], SMPD, 20., l_mask)
        call assert_true(count(l_mask) > 0 .and. count(l_mask) < BOX*BOX, 'the disc mask selects part of the image')
        pcavec = img%serialize(l_mask)
        call assert_int(count(l_mask), size(pcavec), 'serialize: one value per masked pixel')
        if( size(pcavec) == count(l_mask) )then
            call assert_true(all(pcavec == pack(orig, l_mask)), 'serialize: the masked pixels in column-major order')
        endif
        call img%kill
        call img_msk%kill
    end subroutine test_masked_serialize

    !> pixel (i,j) holds i + 100 j, so no two pixels are equal
    subroutine make_distinct_image( img, orig )
        type(image),          intent(inout) :: img
        real,    allocatable, intent(out)   :: orig(:,:,:)
        integer :: i, j
        allocate(orig(BOX,BOX,1))
        do j = 1, BOX
            do i = 1, BOX
                orig(i,j,1) = real(i) + 100. * real(j)
            enddo
        enddo
        call img%new([BOX,BOX,1], SMPD)
        call img%set_rmat(orig, .false.)
    end subroutine make_distinct_image

end module simple_image_serialize_tester
