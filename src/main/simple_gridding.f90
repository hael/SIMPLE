! utilities for convolution interpolation (gridding)
module simple_gridding
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_projector,  only: projector
use simple_kbinterpol, only: kbinterpol
implicit none

public :: mul_w_instr, gen_instrfun_img
private
#include "simple_local_flags.inc"

contains

    !>  \brief generates instrument function image for division of real-space images
    subroutine gen_instrfun_img( instrfun_img, kbwin )
        class(image),      intent(inout) :: instrfun_img
        class(kbinterpol), intent(in)    :: kbwin
        real, allocatable :: w(:)
        real    :: arg
        integer :: ldim(3), center(3), i,j,k
        ldim = instrfun_img%get_ldim()
        if( any(ldim==0) .or. instrfun_img%is_ft() .or. .not.instrfun_img%square_dims() )then
            THROW_HARD('Erroneous image in gen_instrfun_img')
        endif
        allocate(w(ldim(1)),source=1.)
        center = ldim/2+1
        ! kaiser-besel window
        do i = 1,ldim(1)
            arg  = real(i-center(1))/real(ldim(1))
            w(i) = kbwin%instr(arg)
        end do
        if( instrfun_img%is_2d() )then
            !$omp parallel do collapse(2) private(i,j) default(shared) proc_bind(close) schedule(static)
            do i = 1,ldim(1)
                do j = 1,ldim(2)
                    call instrfun_img%set([i,j,1], w(i)*w(j))
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do collapse(3) private(i,j,k) default(shared) proc_bind(close) schedule(static)
            do i = 1,ldim(1)
                do j = 1,ldim(2)
                    do k = 1,ldim(3)
                        call instrfun_img%set([i,j,k], w(i)*w(j)*w(k))
                    enddo
                enddo
            enddo
            !$omp end parallel do
        endif
        deallocate(w)
    end subroutine gen_instrfun_img

    !> \brief  for dividing a real image with the instrument function
    subroutine mul_w_instr( vol, kbwin )
        class(image),      intent(inout) :: vol
        class(kbinterpol), intent(in)    :: kbwin
        real, allocatable :: w1(:), w2(:), w3(:)
        integer :: ldim(3), i, j, k, lims(3,2)
        real    :: arg
        ! get the limits
        ldim = vol%get_ldim()
        if( vol%is_ft() ) THROW_HARD('Volume must not be FTed')
        lims(:,1) = 1
        lims(:,2) = ldim
        ! make the window
        allocate( w1(lims(1,1):lims(1,2)), source=1., stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk("In: divide_w_instr; simple_gridding")
        allocate( w2(lims(2,1):lims(2,2)), source=1., stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk("In: divide_w_instr; simple_gridding")
        allocate( w3(lims(3,1):lims(3,2)), source=1., stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk("In: divide_w_instr; simple_gridding")
        ! calculate the values
        call calc_w(lims(1,1), lims(1,2), w1)
        if( vol%square_dims() )then
            w2 = w1
            if( vol%is_3d() ) w3 = w1
        else
            call calc_w(lims(2,1), lims(2,2), w2)
            if( vol%is_3d() )call calc_w(lims(3,1), lims(3,2), w3)
        endif
        ! divide the image
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
        do i=lims(1,1),lims(1,2)
            do j=lims(2,1),lims(2,2)
                do k=lims(3,1),lims(3,2)
                    call vol%mul_rmat_at([i,j,k], w1(i)*w2(j)*w3(k))
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(w1,w2,w3)

        contains

            subroutine calc_w( from, to, w )
                integer, intent(in)    :: from, to
                real,    intent(inout) :: w(from:to)
                real    :: ci, w_zero
                integer :: i, length
                length = to - from + 1
                w_zero = kbwin%instr(0.)
                ci     = -real(length)/2.
                do i = from,to
                    arg  = ci/real(length)
                    w(i) = kbwin%instr(arg) / w_zero
                    ci   = ci+1.
                end do
            end subroutine calc_w

    end subroutine mul_w_instr

end module simple_gridding
