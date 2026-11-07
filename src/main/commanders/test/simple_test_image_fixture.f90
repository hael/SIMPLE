!@descr: image test fixture shared by testers: a centred square (a cube in a volume) of ones on zero
module simple_test_image_fixture
use simple_image, only: image
implicit none
private
public :: draw_square

contains

    !> the image becomes, in real space, a centred square of side 2*sqrad (a cube in a volume) of ones on zero
    subroutine draw_square( img, sqrad )
        class(image), intent(inout) :: img
        integer,      intent(in)    :: sqrad
        real, allocatable :: rmat(:,:,:)
        integer :: ldim(3), lo(3), hi(3)
        ldim = img%get_ldim()
        if( any(ldim(1:2) <= sqrad) .or. (ldim(3) > 1 .and. ldim(3) <= sqrad) )then
            error stop 'draw_square: the image is too small for the square'
        endif
        allocate(rmat(ldim(1),ldim(2),ldim(3)), source=0.)
        lo = ldim/2 - sqrad + 1
        hi = ldim/2 + sqrad
        if( ldim(3) == 1 )then
            lo(3) = 1
            hi(3) = 1
        endif
        rmat(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 1.
        call img%set_rmat(rmat, .false.)
    end subroutine draw_square

end module simple_test_image_fixture
