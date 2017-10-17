! utilities for convolution interpolation (gridding)
module simple_gridding
#include "simple_lib.f08"
    
use simple_image,      only: image
use simple_kbinterpol, only: kbinterpol

implicit none

contains

    !>  \brief  prepare image for gridding interpolation in Fourier space
    subroutine prep4cgrid( img, img4grid, msk, kbwin )
        class(image),       intent(inout) :: img, img4grid
        real,               intent(in)    :: msk
        class(kbinterpol), intent(in)     :: kbwin
        call img%bwd_ft                      ! make sure not FTed
        call img%zero_background(msk)        ! remove median backround
        call img%pad(img4grid, backgr=0.)    ! padding in real space
        call divide_w_instr(img4grid, kbwin) ! division w instr in real space
        call img4grid%fwd_ft                 ! return the Fourier transform
    end subroutine prep4cgrid

    !> \brief  for dividing a real or complex image with the instrument function
    subroutine divide_w_instr( img, kbwin )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(image),      intent(inout) :: img
        class(kbinterpol), intent(in)    :: kbwin
        real, allocatable :: w1(:), w2(:), w3(:)
        integer :: ldim(3), i, j, k, lims(3,2)
        real    :: arg
        ! get the limits
        ldim = img%get_ldim()
        if( img%is_ft() )then
            lims = img%loop_lims(2)
        else
            lims(:,1) = 1
            lims(:,2) = ldim
        endif
        ! make the window
        allocate( w1(lims(1,1):lims(1,2)), w2(lims(2,1):lims(2,2)),&
        w3(lims(3,1):lims(3,2)), stat=alloc_stat )
        allocchk("In: divide_w_instr; simple_gridding")
        ! calculate the values
        call calc_w(lims(1,:), ldim(1), w1)
        if( img%square_dims() )then
            w2 = w1
            if(img%is_3d()) w3 = w1
        else
            call calc_w(lims(2,:), ldim(2), w2)
            call calc_w(lims(3,:), ldim(3), w3)
        endif
        if( img%is_2d() ) w3 = 1.
        ! divide the image
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
        do i=lims(1,1),lims(1,2)
            do j=lims(2,1),lims(2,2)
                do k=lims(3,1),lims(3,2)
                    call img%div_rmat_at([i,j,k], w1(i)*w2(j)*w3(k))
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(w1,w2,w3)

        contains

            subroutine calc_w( lims_here, ldim_here, w )
                integer, intent(in)    :: lims_here(2), ldim_here
                real,    intent(inout) :: w(lims_here(1):lims_here(2))
                real    :: ci
                integer :: i
                ci = -real(ldim_here)/2.
                do i=lims_here(1),lims_here(2)
                    if( img%is_ft() )then
                        arg = real(i)/real(ldim_here)
                    else
                        arg = ci/real(ldim_here)
                    endif
                    w(i) = kbwin%instr(arg)
                    ci = ci+1.
                end do
            end subroutine calc_w

    end subroutine divide_w_instr

    !> \brief  for dividing a real image with the instrument function
    subroutine mul_w_instr( vol, kbwin )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(image),      intent(inout) :: vol
        class(kbinterpol), intent(in)    :: kbwin
        real, allocatable :: w1(:), w2(:), w3(:)
        integer :: ldim(3), i, j, k, lims(3,2)
        real    :: arg
        ! get the limits
        ldim = vol%get_ldim()
        if( vol%is_ft() )stop 'Volume must not be FTed'
        lims(:,1) = 1
        lims(:,2) = ldim
        ! make the window
        allocate( w1(lims(1,1):lims(1,2)), w2(lims(2,1):lims(2,2)),&
        w3(lims(3,1):lims(3,2)), stat=alloc_stat )
        allocchk("In: divide_w_instr; simple_gridding w1 w2")
        ! calculate the values
        call calc_w(lims(1,:), ldim(1), w1)
        if( vol%square_dims() )then
            w2 = w1
            w3 = w1
        else
            call calc_w(lims(2,:), ldim(2), w2)
            call calc_w(lims(3,:), ldim(3), w3)
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

            subroutine calc_w( lims_here, ldim_here, w )
                integer, intent(in)    :: lims_here(2), ldim_here
                real,    intent(inout) :: w(lims_here(1):lims_here(2))
                real    :: ci, w_zero
                integer :: i
                w_zero = kbwin%instr(0.)
                ci = -real(ldim_here)/2.
                do i=lims_here(1),lims_here(2)
                    arg  = ci/real(ldim_here)
                    w(i) = kbwin%instr(arg) / w_zero
                    ci   = ci+1.
                end do
            end subroutine calc_w

    end subroutine mul_w_instr

end module simple_gridding
