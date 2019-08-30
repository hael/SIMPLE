module simple_nano_tseries_subs
use simple_image
include 'simple_lib.f08'
implicit none

public :: subtr_backgr
private

integer, parameter :: NNN=8

type neigh_ptr
    real(kind=c_float), pointer :: p(:,:,:)
end type neigh_ptr

contains

    subroutine subtr_backgr( particle_img, neigh_imgs )
        class(image),    intent(in) :: particle_img, neigh_imgs(NNN)
        real(kind=c_float), pointer :: ptcl_ptr(:,:,:)
        type(neigh_ptr) :: nnptrs(NNN)
        integer :: i
        ! get the data
        call particle_img%get_rmat_ptr(ptcl_ptr)
        do i = 1, NNN
            call neigh_imgs(i)%get_rmat_ptr(ptcl_ptr)
            call neigh_imgs(i)%get_rmat_ptr(nnptrs(i)%p)
        end do
    end subroutine subtr_backgr

end module simple_nano_tseries_subs
